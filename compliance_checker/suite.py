"""
Compliance Checker suite runner
"""

import codecs
import inspect
import itertools
import os
import re
import subprocess
import sys
import textwrap
import warnings
from collections import defaultdict
from datetime import datetime, timezone
from importlib.metadata import entry_points
from operator import itemgetter
from pathlib import Path
from urllib.parse import urlparse

import requests
from lxml import etree as ET
from netCDF4 import Dataset
from owslib.sos import SensorObservationService
from owslib.swe.sensor.sml import SensorML
from packaging.version import parse

from compliance_checker import __version__, tempnc
from compliance_checker.base import BaseCheck, GenericFile, Result, fix_return_value
from compliance_checker.protocols import cdl, netcdf, opendap, zarr

# Ensure output is encoded as Unicode when checker output is redirected or piped
if sys.stdout.encoding is None:
    sys.stdout = codecs.getwriter("utf8")(sys.stdout)
if sys.stderr.encoding is None:
    sys.stderr = codecs.getwriter("utf8")(sys.stderr)


def extract_docstring_summary(docstring):
    """
    Returns a dedented docstring without parameter information
    :param docstring: A docstring
    :type docstring: str
    :returns: str
    """
    # return a dedented, then indented two spaces docstring with leading and
    # trailing whitespace removed.
    return re.sub(
        r"^(?=.)",
        "  ",
        textwrap.dedent(
            re.split(r"\n\s*:\w", docstring, flags=re.MULTILINE)[0],
        ).strip(),
        flags=re.MULTILINE,
    )


class CheckSuite:
    checkers = (
        {}
    )  # Base dict of checker names to BaseCheck derived types, override this in your CheckSuite implementation
    templates_root = "compliance_checker"  # modify to load alternative Jinja2 templates

    def __init__(self, options=None):
        self.col_width = 40
        self.options = options or {}

    @classmethod
    def _get_generator_plugins(cls):
        """
        Return a list of classes from external plugins that are used to
        generate checker classes
        """

        if not hasattr(cls, "suite_generators"):
            gens = entry_points(
                group="compliance_checker.generators",
            )
            cls.suite_generators = [x.load() for x in gens]

        return cls.suite_generators

    def _print_suites(self, verbose=0):
        """
        Prints out available check suites.  If the verbose argument is True,
        includes the internal module version number of the check and also displays
        "latest" meta-versions.
        :param check_suite: Check suite object
        :param verbose: Integer indicating whether to print verbose output
        :type verbose: int
        """
        for checker in sorted(self.checkers.keys()):
            version = getattr(self.checkers[checker], "_cc_checker_version", "???")
            if verbose > 0:
                print(f" - {checker} (v{version})")
            elif ":" in checker and not checker.endswith(
                ":latest",
            ):  # Skip the "latest" output
                print(f" - {checker}")

    def _print_checker(self, checker_obj):
        """
        Prints each available check and a description with an abridged
        docstring for a given checker object
        :param checker_obj: Checker object on which to operate
        :type checker_obj: subclass of compliance_checker.base.BaseChecker
        """

        check_functions = self._get_checks(checker_obj, {}, defaultdict(lambda: None))
        for c, _ in check_functions:
            print(f"- {c.__name__}")
            if c.__doc__ is not None:
                u_doc = c.__doc__
                print(f"\n{extract_docstring_summary(u_doc)}\n")

    @classmethod
    def add_plugin_args(cls, parser):
        """
        Add command line arguments for external plugins that generate checker
        classes
        """
        for gen in cls._get_generator_plugins():
            gen.add_arguments(parser)

    @classmethod
    def load_generated_checkers(cls, args):
        """
        Load checker classes from generator plugins
        """

        for gen in cls._get_generator_plugins():
            checkers = gen.get_checkers(args)
            cls.checkers.update(checkers)

    @classmethod
    def load_all_available_checkers(cls):
        """
        Helper method to retrieve all sub checker classes derived from various
        base classes.
        """
        cls._load_checkers(
            entry_points(group="compliance_checker.suites"),
        )

    @classmethod
    def _load_checkers(cls, checkers):
        """
        Loads up checkers in an iterable into the class checkers dict
        :param checkers: An iterable containing the checker objects
        """

        for c in checkers:
            try:
                check_obj = c.load()
                if hasattr(check_obj, "_cc_spec") and hasattr(
                    check_obj,
                    "_cc_spec_version",
                ):
                    check_version_str = ":".join(
                        (check_obj._cc_spec, check_obj._cc_spec_version),
                    )
                    cls.checkers[check_version_str] = check_obj
                # TODO: remove this once all checkers move over to the new
                #       _cc_spec, _cc_spec_version
                else:
                    # if _cc_spec and _cc_spec_version attributes aren't
                    # present, fall back to using name attribute
                    checker_name = getattr(check_obj, "name", None) or getattr(
                        check_obj,
                        "_cc_spec",
                        None,
                    )
                    warnings.warn(
                        "Checker for {} should implement both "
                        '"_cc_spec" and "_cc_spec_version" '
                        'attributes. "name" attribute is deprecated. '
                        "Assuming checker is latest version.",
                        DeprecationWarning,
                        stacklevel=2,
                    )
                    # append "unknown" to version string since no versioning
                    # info was provided
                    cls.checkers[f"{checker_name}:unknown"] = check_obj

            except Exception as e:
                print("Could not load", c, ":", e, file=sys.stderr)
        # find the latest version of versioned checkers and set that as the
        # default checker for compliance checker if no version is specified
        ver_checkers = sorted([c.split(":", 1) for c in cls.checkers if ":" in c])
        for spec, versions in itertools.groupby(ver_checkers, itemgetter(0)):
            version_nums = [v[-1] for v in versions]
            try:
                latest_version = str(max(parse(v) for v in version_nums))
            # if the version can't be parsed, do it according to character collation
            except ValueError:
                latest_version = max(version_nums)
            cls.checkers[spec] = cls.checkers[spec + ":latest"] = cls.checkers[
                ":".join((spec, latest_version))
            ]

    def _get_checks(self, checkclass, include_checks, skip_checks):
        """
        Helper method to retrieve check methods from a Checker class.  Excludes
        any checks in `skip_checks`.

        The name of the methods in the Checker class should start with "check_"
        for this method to find them.
        :param checkclass BaseCheck: The checker class being considered
        :param skip_checks list: A list of strings with the names of the check
                                 methods to skip or include, depending on the
                                 value of `skip_flag`.
        :param skip_flag bool: A boolean parameter to determine whether to
                               skip over checks specified (True) or only
                               include the checks specified (False).
        """
        meths = inspect.getmembers(checkclass, inspect.isroutine)
        # return all check methods not among the skipped checks
        returned_checks = []
        if include_checks:
            for fn_name, fn_obj in meths:
                if fn_name in include_checks:
                    returned_checks.append((fn_obj, skip_checks[fn_name]))
        else:
            for fn_name, fn_obj in meths:
                if (
                    fn_name.startswith("check_")
                    and skip_checks[fn_name] != BaseCheck.HIGH
                ):
                    returned_checks.append((fn_obj, skip_checks[fn_name]))

        return returned_checks

    def _run_check(self, check_method, ds, max_level):
        """
        Runs a check and appends a result to the values list.
        @param bound method check_method: a given check method
        @param netCDF4 dataset ds
        @param int max_level: check level
        @return list: list of Result objects
        """
        val = check_method(ds)
        if hasattr(val, "__iter__"):
            # Handle OrderedDict when we need to modify results in a superclass
            # i.e. some checks in CF 1.7 which extend CF 1.6 behaviors
            if isinstance(val, dict):
                val_iter = val.values()
            else:
                val_iter = val
            check_val = []
            for v in val_iter:
                res = fix_return_value(
                    v,
                    check_method.__func__.__name__,
                    check_method,
                    check_method.__self__,
                )
                if max_level is None or res.weight > max_level:
                    check_val.append(res)

            return check_val
        else:
            check_val = fix_return_value(
                val,
                check_method.__func__.__name__,
                check_method,
                check_method.__self__,
            )
            if max_level is None or check_val.weight > max_level:
                return [check_val]
            else:
                return []

    def _get_check_versioned_name(self, check_name):
        """
        The compliance checker allows the user to specify a
        check without a version number but we want the report
        to specify the version number.

        Returns the check name with the version number it checked
        """
        if ":" not in check_name or ":latest" in check_name:
            check_name = ":".join(
                (check_name.split(":")[0], self.checkers[check_name]._cc_spec_version),
            )
        return check_name

    def _get_check_url(self, check_name):
        """
        Return the check's reference URL if it exists. If not, return empty str.
        @param check_name str: name of the check being run returned by
                               _get_check_versioned_name()
        """
        return getattr(self.checkers[check_name], "_cc_url", "")

    def _get_valid_checkers(self, ds, checker_names):
        """
        Returns a filtered list of 2-tuples: (name, valid checker) based on the ds object's type and
        the user selected names.
        """

        if len(self.checkers) <= 0:
            raise ValueError("No checkers could be found.")

        if len(checker_names) == 0:
            checker_names = list(self.checkers.keys())

        args = [
            (name, self.checkers[name])
            for name in checker_names
            if name in self.checkers
        ]
        valid = []

        all_checked = {a[1] for a in args}  # only class types
        checker_queue = set(args)
        while len(checker_queue):
            name, a = checker_queue.pop()
            # is the current dataset type in the supported filetypes
            # for the checker class?
            if type(ds) in a().supported_ds:
                valid.append((name, a))

            # add subclasses of SOS checks
            if "ioos_sos" in name:
                for subc in a.__subclasses__():
                    if subc not in all_checked:
                        all_checked.add(subc)
                        checker_queue.add((name, subc))

        return valid

    @classmethod
    def _process_skip_checks(cls, skip_checks):
        """
        Processes an iterable of skip_checks with strings and returns a dict
        with <check_name>: <max_skip_level> pairs
        """

        check_dict = defaultdict(lambda: None)
        # A is for "all", "M" is for medium, "L" is for low
        check_lookup = {"A": BaseCheck.HIGH, "M": BaseCheck.MEDIUM, "L": BaseCheck.LOW}

        for skip_check_spec in skip_checks:
            split_check_spec = skip_check_spec.split(":")
            check_name = split_check_spec[0]
            if len(split_check_spec) < 2:
                check_max_level = BaseCheck.HIGH
            else:
                try:
                    check_max_level = check_lookup[split_check_spec[1]]
                except KeyError:
                    warnings.warn(
                        f"Skip specifier '{split_check_spec[1]}' on check '{check_name}' not found,"
                        " defaulting to skip entire check",
                        stacklevel=2,
                    )
                    check_max_level = BaseCheck.HIGH

            check_dict[check_name] = check_max_level

        return check_dict

    def run(self, ds, skip_checks, *checker_names):
        warnings.warn(
            "suite.run is deprecated, use suite.run_all in calls instead",
            stacklevel=2,
        )
        return self.run_all(ds, checker_names, skip_checks=skip_checks)

    def run_all(self, ds, checker_names, include_checks=None, skip_checks=None):
        """
        Runs this CheckSuite on the dataset with all the passed Checker instances.

        Returns a dictionary mapping checker names to a 2-tuple of their grouped scores and errors/exceptions while running checks.
        """

        ret_val = {}
        checkers = self._get_valid_checkers(ds, checker_names)

        if skip_checks is not None:
            skip_check_dict = CheckSuite._process_skip_checks(skip_checks)
        else:
            skip_check_dict = defaultdict(lambda: None)

        if include_checks:
            include_dict = dict.fromkeys(include_checks, 0)
        else:
            include_dict = {}

        if len(checkers) == 0:
            print(
                "No valid checkers found for tests '{}'".format(
                    ",".join(checker_names),
                ),
            )

        for checker_name, checker_class in checkers:
            # TODO: maybe this a little more reliable than depending on
            #       a string to determine the type of the checker -- perhaps
            #       use some kind of checker object with checker type and
            #       version baked in
            checker_type_name = checker_name.split(":")[0]
            checker_opts = self.options.get(checker_type_name, {})

            # instantiate a Checker object
            try:
                checker = checker_class(options=checker_opts)
            # hacky fix for no options in constructor
            except TypeError:
                checker = checker_class()
            # TODO? : Why is setup(ds) called at all instead of just moving the
            #         checker setup into the constructor?
            # setup method to prep
            checker.setup(ds)

            checks = self._get_checks(checker, include_dict, skip_check_dict)
            vals = []
            errs = {}  # check method name -> (exc, traceback)

            for c, max_level in checks:
                try:
                    vals.extend(self._run_check(c, ds, max_level))
                except Exception as e:
                    errs[c.__func__.__name__] = (e, sys.exc_info()[2])

            # score the results we got back
            groups = self.scores(vals)

            # invoke finalizer explicitly
            del checker

            ret_val[checker_name] = groups, errs

        return ret_val

    @classmethod
    def passtree(cls, groups, limit):
        for r in groups:
            if r.children:
                x = cls.passtree(r.children, limit)
                if r.weight >= limit and x is False:
                    return False

            if r.weight >= limit and r.value[0] != r.value[1]:
                return False

        return True

    def build_structure(self, check_name, groups, source_name, limit=1):
        """
        Compiles the checks, results and scores into an aggregate structure which looks like:

            {
              "scored_points": 396,
              "low_count": 0,
              "possible_points": 400,
              "testname": "gliderdac",
              "medium_count": 2,
              "source_name": ".//rutgers/ru01-20140120T1444/ru01-20140120T1649.nc",
              "high_count": 0,
              "all_priorities" : [...],
              "high_priorities": [...],
              "medium_priorities" : [...],
              "low_priorities" : [...]
            }

        @param check_name  The test which was run
        @param groups      List of results from compliance checker
        @param source_name Source of the dataset, used for title
        """
        aggregates = {}

        aggregates["scored_points"] = 0
        aggregates["possible_points"] = 0
        high_priorities = []
        medium_priorities = []
        low_priorities = []
        all_priorities = []

        aggregates["high_count"] = 0
        aggregates["medium_count"] = 0
        aggregates["low_count"] = 0

        def named_function(result):
            for child in result.children:
                all_priorities.append(child)
                named_function(child)

        # For each result, bin them into the appropriate category, put them all
        # into the all_priorities category and add up the point values
        for res in groups:
            if res.weight < limit:
                continue
            # If the result has 0 possible points, then it was not valid for
            # this dataset and contains no meaningful information
            if res.value[1] == 0:
                continue
            aggregates["scored_points"] += res.value[0]
            aggregates["possible_points"] += res.value[1]
            if res.weight == 3:
                high_priorities.append(res)
                if res.value[0] < res.value[1]:
                    aggregates["high_count"] += 1
            elif res.weight == 2:
                medium_priorities.append(res)
                if res.value[0] < res.value[1]:
                    aggregates["medium_count"] += 1
            else:
                low_priorities.append(res)
                if res.value[0] < res.value[1]:
                    aggregates["low_count"] += 1
            all_priorities.append(res)
            # Some results have children
            # We don't render children inline with the top three tables, but we
            # do total the points and display the messages
            named_function(res)

        aggregates["high_priorities"] = high_priorities
        aggregates["medium_priorities"] = medium_priorities
        aggregates["low_priorities"] = low_priorities
        aggregates["all_priorities"] = all_priorities
        aggregates["testname"] = self._get_check_versioned_name(check_name)
        aggregates["source_name"] = source_name
        aggregates["scoreheader"] = self.checkers[check_name]._cc_display_headers
        aggregates["cc_spec_version"] = self.checkers[check_name]._cc_spec_version
        aggregates["cc_url"] = self._get_check_url(aggregates["testname"])
        aggregates["report_timestamp"] = datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ",
        )
        aggregates["cc_version"] = __version__
        return aggregates

    def dict_output(self, check_name, groups, source_name, limit):
        """
        Builds the results into a JSON structure and writes it to the file buffer.

        @param check_name      The test which was run
        @param groups          List of results from compliance checker
        @param output_filename Path to file to save output
        @param source_name     Source of the dataset, used for title
        @param limit           Integer value for limiting output
        """
        aggregates = self.build_structure(check_name, groups, source_name, limit)
        return self.serialize(aggregates)

    def serialize(self, o):
        """
        Returns a safe serializable object that can be serialized into JSON.

        @param o Python object to serialize
        """
        if isinstance(o, (list, tuple)):
            return [self.serialize(i) for i in o]
        if isinstance(o, dict):
            return {k: self.serialize(v) for k, v in o.items()}
        if isinstance(o, datetime):
            return o.isoformat()
        if isinstance(o, Result):
            return self.serialize(o.serialize())
        return o

    def checker_html_output(self, check_name, groups, source_name, limit):
        """
        Renders the HTML output for a single test using Jinja2 and returns it
        as a string.

        @param check_name      The test which was run
        @param groups          List of results from compliance checker
        @param source_name     Source of the dataset, used for title
        @param limit           Integer value for limiting output
        """
        from jinja2 import Environment, PackageLoader

        self.j2 = Environment(
            loader=PackageLoader(self.templates_root, "data/templates"),
        )
        template = self.j2.get_template("ccheck.html.j2")

        template_vars = self.build_structure(check_name, groups, source_name, limit)
        return template.render(**template_vars)

    def html_output(self, checkers_html):
        """
        Renders the HTML output for multiple tests and returns it as a string.

        @param checkers_html     List of HTML for single tests as returned by
                                 checker_html_output
        """
        # Note: This relies on checker_html_output having been called so that
        # self.j2 is initialised
        template = self.j2.get_template("ccheck_wrapper.html.j2")
        return template.render(checkers=checkers_html)

    def get_points(self, groups, limit):
        score_list = []
        score_only_list = []

        for g in groups:
            if g.weight >= limit:
                score_only_list.append(g.value)

        # checks where all pertinent sections passed
        all_passed = sum(x[0] == x[1] for x in score_only_list)
        out_of = len(score_only_list)

        # sorts lists into high/medium/low order
        score_list.sort(key=lambda x: x.weight, reverse=True)

        return score_list, all_passed, out_of

    def standard_output(self, ds, limit, check_name, groups):
        """
        Generates the Terminal Output for Standard cases

        Returns the dataset needed for the verbose output, as well as the failure flags.
        """

        score_list, points, out_of = self.get_points(groups, limit)
        issue_count = out_of - points

        # Let's add the version number to the check name if it's missing
        check_name = self._get_check_versioned_name(check_name)
        check_url = self._get_check_url(check_name)
        width = 2 * self.col_width

        # NOTE: printing and use of .center()
        # Nested .format() calls should be avoided when possible.
        # As a future enhancement, a string.Template string might work best here
        # but for the time being individual lines are printed and centered with
        # .center()

        print("\n")
        print("-" * width)
        print("IOOS Compliance Checker Report".center(width))
        print(f"Version {__version__}".center(width))
        print(
            "Report generated {}".format(
                datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            ).center(width),
        )
        print(f"{check_name}".center(width))
        print(f"{check_url}".center(width))
        print("-" * width)
        if issue_count > 0:
            print("Corrective Actions".center(width))
            plural = "" if issue_count == 1 else "s"
            print(
                f"{os.path.basename(ds)} has {issue_count} potential issue{plural}",
            )

        return [groups, points, out_of]

    def standard_output_generation(self, groups, limit, points, out_of, check):
        """
        Generates the Terminal Output
        """
        if points < out_of:
            self.reasoning_routine(groups, check, priority_flag=limit)
        else:
            print("All tests passed!")

    def reasoning_routine(self, groups, check, priority_flag=3, _top_level=True):
        """
        print routine performed
        @param list groups: the Result groups
        @param str check: checker name
        @param int priority_flag: indicates the weight of the groups
        @param bool _top_level: indicates the level of the group so as to
                                print out the appropriate header string
        """

        def weight_sort(result):
            return result.weight

        groups_sorted = sorted(groups, key=weight_sort, reverse=True)

        # create dict of the groups -> {level: [reasons]}
        result = {
            key: [v for v in valuesiter if v.value[0] != v.value[1]]
            for key, valuesiter in itertools.groupby(groups_sorted, key=weight_sort)
        }
        priorities = self.checkers[check]._cc_display_headers

        def process_table(res, check):
            """Recursively calls reasoning_routine to parse out child reasons
            from the parent reasons.
            @param Result res: Result object
            @param str check: checker name"""

            issue = res.name
            if not res.children:
                reasons = res.msgs
            else:
                child_reasons = self.reasoning_routine(
                    res.children,
                    check,
                    _top_level=False,
                )
                # there shouldn't be messages if there are children
                # is this a valid assumption?
                reasons = child_reasons

            return issue, reasons

        # iterate in reverse to the min priority requested;
        # the higher the limit, the more lenient the output
        proc_strs = ""
        for level in range(3, priority_flag - 1, -1):
            level_name = priorities.get(level, level)
            # print headers
            proc_strs = []

            # skip any levels that aren't in the result
            if level not in result:
                continue

            # skip any empty result levels
            if len(result[level]) > 0:
                # only print priority headers at top level, i.e. non-child
                # datasets
                if _top_level:
                    width = 2 * self.col_width
                    print("\n")
                    print("{:^{width}}".format(level_name, width=width))
                    print("-" * width)

                data_issues = [process_table(res, check) for res in result[level]]

                has_printed = False
                for issue, reasons in data_issues:
                    # if this isn't the first printed issue, add a newline
                    # separating this and the previous level
                    if has_printed:
                        print("")
                    # join alphabetized reasons together
                    reason_str = "\n".join(
                        f"* {r}" for r in sorted(reasons, key=lambda x: x[0])
                    )
                    proc_str = f"{issue}\n{reason_str}"
                    print(proc_str)
                    proc_strs.append(proc_str)
                    has_printed = True
        return "\n".join(proc_strs)

    def process_doc(self, doc):
        """
        Attempt to parse an xml string conforming to either an SOS or SensorML
        dataset and return the results
        """
        xml_doc = ET.fromstring(doc)
        if xml_doc.tag == "{http://www.opengis.net/sos/1.0}Capabilities":
            ds = SensorObservationService(None, xml=doc)
            # SensorObservationService does not store the etree doc root,
            # so maybe use monkey patching here for now?
            ds._root = xml_doc

        elif xml_doc.tag == "{http://www.opengis.net/sensorML/1.0.1}SensorML":
            ds = SensorML(xml_doc)
        else:
            raise ValueError(f"Unrecognized XML root element: {xml_doc.tag}")
        return ds

    def generate_dataset(self, cdl_path):
        """
        Use ncgen to generate a netCDF file from a .cdl file
        Returns the path to the generated netcdf file. If ncgen fails, uses
        sys.exit(1) to terminate program so a long stack trace is not reported
        to the user.

        :param str cdl_path: Absolute path to cdl file that is used to generate netCDF file
        """
        if isinstance(cdl_path, str):
            cdl_path = Path(cdl_path)
        ds_str = cdl_path.with_suffix(".nc")

        # generate netCDF-4 file
        iostat = subprocess.run(
            ["ncgen", "-k", "nc4", "-o", ds_str, cdl_path],
            stderr=subprocess.PIPE,
        )
        if iostat.returncode != 0:
            # if not successful, create netCDF classic file
            print(
                "netCDF-4 file could not be generated from cdl file with " + "message:",
            )
            print(iostat.stderr.decode())
            print("Trying to create netCDF Classic file instead.")
            iostat = subprocess.run(
                ["ncgen", "-k", "nc3", "-o", ds_str, cdl_path],
                stderr=subprocess.PIPE,
            )
            if iostat.returncode != 0:
                # Exit program if neither a netCDF Classic nor a netCDF-4 file
                # could be created.
                print(
                    "netCDF Classic file could not be generated from cdl file"
                    + "with message:",
                )
                print(iostat.stderr.decode())
                sys.exit(1)
        return ds_str

    def load_dataset(self, ds_str):
        """
        Returns an instantiated instance of either a netCDF file or an SOS
        mapped DS object.

        :param str ds_str: URL of the resource to load
        """
        if isinstance(ds_str, Path):
            ds_str = str(ds_str)

        # If it's a remote URL load it as a remote resource, otherwise treat it
        # as a local resource.
        pr = urlparse(ds_str)
        if pr.netloc:
            return self.load_remote_dataset(ds_str)
        else:
            return self.load_local_dataset(ds_str)

    def check_remote_netcdf(self, ds_str):
        if netcdf.is_remote_netcdf(ds_str):
            response = requests.get(ds_str, allow_redirects=True, timeout=60)
            try:
                return Dataset(
                    urlparse(response.url).path,
                    memory=response.content,
                )
            except OSError:
                # handle case when netCDF C libs weren't compiled with
                # in-memory support by using tempfile
                with tempnc(response.content) as _nc:
                    return Dataset(_nc)

    def load_remote_dataset(self, ds_str):
        """
        Returns a dataset instance for the remote resource, either OPeNDAP or SOS

        :param str ds_str: URL to the remote resource
        """

        url_parsed = urlparse(ds_str)
        # ERDDAP TableDAP request

        nc_remote_result = self.check_remote_netcdf(ds_str)
        if nc_remote_result:
            return nc_remote_result

        # if application/x-netcdf wasn't detected in the Content-Type headers
        # and this is some kind of erddap tabledap form, then try to get the
        # .ncCF file from ERDDAP
        elif "tabledap" in ds_str and not url_parsed.query:
            # modify ds_str to contain the full variable request
            variables_str = opendap.create_DAP_variable_str(ds_str)

            # join to create a URL to an .ncCF resource
            ds_str = f"{ds_str}.ncCF?{variables_str}"

        nc_remote_result = self.check_remote_netcdf(ds_str)
        if nc_remote_result:
            return nc_remote_result

        # if it's just an OPeNDAP endpoint, use that
        elif opendap.is_opendap(ds_str):
            return Dataset(ds_str)

        # Check if the HTTP response is XML, if it is, it's likely SOS so
        # we'll attempt to parse the response as SOS.
        # Some SOS servers don't seem to support HEAD requests.
        # Issue GET instead if we reach here and can't get the response
        response = requests.get(ds_str, allow_redirects=True, timeout=60)
        content_type = response.headers.get("content-type")
        if content_type.split(";")[0] == "text/xml":
            return self.process_doc(response.content)
        elif content_type.split(";")[0] == "application/x-netcdf":
            return Dataset(
                urlparse(response.url).path,
                memory=response.content,
            )
        else:
            raise ValueError(
                f"Unknown service with content-type: {content_type}",
            )

    def load_local_dataset(self, ds_str):
        """
        Returns a dataset instance for the local resource

        :param ds_str: Path to the resource
        """
        if cdl.is_cdl(ds_str):
            ds_str = self.generate_dataset(ds_str)

        if zarr.is_zarr(ds_str):
            return Dataset(zarr.as_zarr(ds_str))

        if netcdf.is_netcdf(ds_str):
            return Dataset(ds_str)

        # Assume this is just a Generic File if it exists
        if os.path.isfile(ds_str):
            return GenericFile(ds_str)

        raise ValueError("File is an unknown format")

    def scores(self, raw_scores):
        """
        Transforms raw scores from a single checker into a fully tallied and grouped scoreline.
        """
        grouped = self._group_raw(raw_scores)

        return grouped

    def _group_raw(self, raw_scores, cur=None, level=1):
        """
        Internal recursive method to group raw scores into a cascading score summary.
        Only top level items are tallied for scores.
        @param list raw_scores: list of raw scores (Result objects)
        """

        def trim_groups(r):
            if isinstance(r.name, tuple) or isinstance(r.name, list):
                new_name = r.name[1:]
            else:
                new_name = []

            return Result(r.weight, r.value, new_name, r.msgs)

        # CHECK FOR TERMINAL CONDITION: all raw_scores.name are single length
        # @TODO could have a problem here with scalar name, but probably still works
        terminal = [len(x.name) for x in raw_scores]
        if terminal == [0] * len(raw_scores):
            return []

        def group_func(r):
            """
            Takes a Result object and slices off the first element of its name
            if its's a tuple. Otherwise, does nothing to the name. Returns the
            Result's name and weight in a tuple to be used for sorting in that
            order in a groupby function.
            @param Result r
            @return tuple (str, int)
            """
            if isinstance(r.name, tuple) or isinstance(r.name, list):
                if len(r.name) == 0:
                    retval = ""
                else:
                    retval = r.name[0:1][0]
            else:
                retval = r.name
            return retval, r.weight

        # END INTERNAL FUNCS ##########################################

        # NOTE until this point, *ALL* Results in raw_scores are
        # individual Result objects.

        # sort then group by name, then by priority weighting
        grouped = itertools.groupby(sorted(raw_scores, key=group_func), key=group_func)

        # NOTE: post-grouping, grouped looks something like
        # [(('Global Attributes', 1), <itertools._grouper at 0x7f10982b5390>),
        #  (('Global Attributes', 3), <itertools._grouper at 0x7f10982b5438>),
        #  (('Not a Global Attr', 1), <itertools._grouper at 0x7f10982b5470>)]
        #  (('Some Variable', 2),     <itertools._grouper at 0x7f10982b5400>),

        ret_val = []

        for k, v in grouped:  # iterate through the grouped tuples
            k = k[0]  # slice ("name", weight_val) --> "name"
            v = list(v)  # from itertools._grouper to list

            cv = self._group_raw(list(map(trim_groups, v)), k, level + 1)
            if len(cv):
                # if this node has children, max weight of children + sum of all the scores
                max_weight = max([x.weight for x in cv])
                sum_scores = tuple(map(sum, list(zip(*([x.value for x in cv])))))
                msgs = []

            else:
                max_weight = max([x.weight for x in v])
                sum_scores = tuple(
                    map(sum, list(zip(*([self._translate_value(x.value) for x in v])))),
                )
                msgs = sum([x.msgs for x in v], [])

            ret_val.append(
                Result(
                    name=k,
                    weight=max_weight,
                    value=sum_scores,
                    children=cv,
                    msgs=msgs,
                ),
            )

        return ret_val

    def _translate_value(self, val):
        """
        Turns shorthand True/False/None checks into full scores (1, 1)/(0, 1)/(0, 0).
        Leaves full scores alone.
        """
        if val is True:
            return (1, 1)
        elif val is False:
            return (0, 1)
        elif val is None:
            return (0, 0)

        return val
