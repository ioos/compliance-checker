"""
Compliance Checker suite runner
"""
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys
import subprocess
import inspect
import itertools
from operator import itemgetter
from netCDF4 import Dataset
from lxml import etree as ET
from distutils.version import StrictVersion
from compliance_checker.base import fix_return_value, Result, GenericFile
from owslib.sos import SensorObservationService
from owslib.swe.sensor.sml import SensorML
from compliance_checker.protocols import opendap, netcdf, cdl
from compliance_checker import MemoizedDataset
try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse
from datetime import datetime
import requests
import textwrap
import codecs
from pkg_resources import working_set
import tabulate


tabulate.PRESERVE_WHITESPACE = True
# Ensure output is encoded as Unicode when checker output is redirected or piped
if sys.stdout.encoding is None:
    sys.stdout = codecs.getwriter('utf8')(sys.stdout)
if sys.stderr.encoding is None:
    sys.stderr = codecs.getwriter('utf8')(sys.stderr)

class CheckSuite(object):

    checkers = {}       # Base dict of checker names to BaseCheck derived types, override this in your CheckSuite implementation

    def __init__(self):
        self.col_width = 60

    @classmethod
    def _get_generator_plugins(cls):
        """
        Return a list of classes from external plugins that are used to
        generate checker classes
        """
        if not hasattr(cls, 'suite_generators'):
            gens = working_set.iter_entry_points('compliance_checker.generators')
            cls.suite_generators = [x.load() for x in gens]

        return cls.suite_generators

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
        for x in working_set.iter_entry_points('compliance_checker.suites'):
            try:
                xl = x.load()
                cls.checkers[':'.join((xl._cc_spec, xl._cc_spec_version))] = xl
            # TODO: remove this once all checkers move over to the new
            #       _cc_spec, _cc_spec_version
            except AttributeError:
                # if there are versioned classes, it will get overwritten by the
                # latest version later.  If there are not, it will be assigned
                # the checker as the main class
                # TODO: nix name attribute in plugins.  Keeping in for now
                #       to provide backwards compatibility
                cls.checkers[getattr(xl, 'name', None) or xl._cc_spec] = xl

            except Exception as e:
                print("Could not load", x, ":", e, file=sys.stderr)
        # find the latest version of versioned checkers and set that as the
        # default checker for compliance checker if no version is specified
        ver_checkers = sorted([c.split(':', 1) for c
                               in cls.checkers if ':' in c])
        for spec, versions in itertools.groupby(ver_checkers, itemgetter(0)):
            version_nums = [v[-1] for v in versions]
            try:
                latest_version = str(max(StrictVersion(v) for v
                                         in version_nums))
            # if the version can't be parsed as a StrictVersion, parse
            # according to character collation
            except ValueError:
                latest_version = max(version_nums)
            cls.checkers[spec] = cls.checkers[spec + ':latest'] = \
                cls.checkers[':'.join((spec, latest_version))]

    def _get_checks(self, checkclass, skip_checks):
        """
        Helper method to retreive check methods from a Checker class.  Excludes
        any checks in `skip_checks`.

        The name of the methods in the Checker class should start with "check_" for this
        method to find them.
        """
        meths = inspect.getmembers(checkclass, inspect.ismethod)
        # return all check methods not among the skipped checks
        return [x[1] for x in meths if x[0].startswith("check_") and
                x[0] not in skip_checks]

    def _run_check(self, check_method, ds):
        val = check_method(ds)

        if isinstance(val, list):
            return [fix_return_value(v, check_method.__func__.__name__, check_method, check_method.__self__) for v in val]

        return [fix_return_value(val, check_method.__func__.__name__, check_method, check_method.__self__)]

    def _get_check_versioned_name(self, check_name):
        """
        The compliance checker allows the user to specify a
        check without a version number but we want the report
        to specify the version number.

        Returns the check name with the version number it checked
        """
        if ':' not in check_name or ':latest' in check_name:
            check_name = ':'.join((check_name.split(':')[0],
                                   self.checkers[check_name]._cc_spec_version))
        return check_name

    def _get_check_url(self, check_name):
        """
        Return the check's reference URL if it exists. If not, return emtpy str.
        @param check_name str: name of the check being run returned by
                               _get_check_versioned_name()
        """
        return getattr(self.checkers[check_name], '_cc_url', '')

    def _get_valid_checkers(self, ds, checker_names):
        """
        Returns a filtered list of 2-tuples: (name, valid checker) based on the ds object's type and
        the user selected names.
        """
        assert len(self.checkers) > 0, "No checkers could be found."

        if len(checker_names) == 0:
            checker_names = list(self.checkers.keys())

        args = [(name, self.checkers[name]) for name in checker_names if name in self.checkers]
        valid = []

        all_checked = set([a[1] for a in args])  # only class types
        checker_queue = set(args)
        while len(checker_queue):
            name, a = checker_queue.pop()
            # is the current dataset type in the supported filetypes
            # for the checker class?
            if type(ds) in a().supported_ds:
                valid.append((name, a))

            # add any subclasses of the checker class
            for subc in a.__subclasses__():
                if subc not in all_checked:
                    all_checked.add(subc)
                    checker_queue.add((name, subc))

        return valid

    def run(self, ds, skip_checks, *checker_names):
        """
        Runs this CheckSuite on the dataset with all the passed Checker instances.

        Returns a dictionary mapping checker names to a 2-tuple of their grouped scores and errors/exceptions while running checks.
        """

        ret_val = {}
        checkers = self._get_valid_checkers(ds, checker_names)

        if len(checkers) == 0:
            print("No valid checkers found for tests '{}'".format(",".join(checker_names)))

        for checker_name, checker_class in checkers:

            checker = checker_class()
            checker.setup(ds)

            checks = self._get_checks(checker, skip_checks)
            vals = []
            errs = {}   # check method name -> (exc, traceback)

            for c in checks:
                try:
                    vals.extend(self._run_check(c, ds))
                except Exception as e:
                    errs[c.__func__.__name__] = (e, sys.exc_info()[2])

            # score the results we got back
            groups = self.scores(vals)

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
        '''
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
        '''
        aggregates = {}

        aggregates['scored_points'] = 0
        aggregates['possible_points'] = 0
        high_priorities = []
        medium_priorities = []
        low_priorities = []
        all_priorities = []

        aggregates['high_count'] = 0
        aggregates['medium_count'] = 0
        aggregates['low_count'] = 0

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
            aggregates['scored_points'] += res.value[0]
            aggregates['possible_points'] += res.value[1]
            if res.weight == 3:
                high_priorities.append(res)
                if res.value[0] < res.value[1]:
                    aggregates['high_count'] += 1
            elif res.weight == 2:
                medium_priorities.append(res)
                if res.value[0] < res.value[1]:
                    aggregates['medium_count'] += 1
            else:
                low_priorities.append(res)
                if res.value[0] < res.value[1]:
                    aggregates['low_count'] += 1
            all_priorities.append(res)
            # Some results have children
            # We don't render children inline with the top three tables, but we
            # do total the points and display the messages
            named_function(res)

        aggregates['high_priorities'] = high_priorities
        aggregates['medium_priorities'] = medium_priorities
        aggregates['low_priorities'] = low_priorities
        aggregates['all_priorities'] = all_priorities
        aggregates['testname'] = self._get_check_versioned_name(check_name)
        aggregates['source_name'] = source_name
        aggregates['scoreheader'] = self.checkers[check_name]._cc_display_headers
        aggregates['cc_spec_version'] = self.checkers[check_name]._cc_spec_version
        aggregates['cc_url'] = self._get_check_url(aggregates['testname'])
        return aggregates

    def dict_output(self, check_name, groups, source_name, limit):
        '''
        Builds the results into a JSON structure and writes it to the file buffer.

        @param check_name      The test which was run
        @param groups          List of results from compliance checker
        @param output_filename Path to file to save output
        @param source_name     Source of the dataset, used for title
        @param limit           Integer value for limiting output
        '''
        aggregates = self.build_structure(check_name, groups, source_name, limit)
        return self.serialize(aggregates)

    def serialize(self, o):
        '''
        Returns a safe serializable object that can be serialized into JSON.

        @param o Python object to serialize
        '''
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
        '''
        Renders the HTML output for a single test using Jinja2 and returns it
        as a string.

        @param check_name      The test which was run
        @param groups          List of results from compliance checker
        @param source_name     Source of the dataset, used for title
        @param limit           Integer value for limiting output
        '''
        from jinja2 import Environment, PackageLoader
        self.j2 = Environment(loader=PackageLoader('compliance_checker', 'data/templates'))
        template = self.j2.get_template('ccheck.html.j2')

        template_vars = self.build_structure(check_name, groups, source_name, limit)
        return template.render(**template_vars)

    def html_output(self, checkers_html):
        '''
        Renders the HTML output for multiple tests and returns it as a string.

        @param checkers_html     List of HTML for single tests as returned by
                                 checker_html_output
        '''
        # Note: This relies on checker_html_output having been called so that
        # self.j2 is initialised
        template = self.j2.get_template('ccheck_wrapper.html.j2')
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
        width = 2*self.col_width
        print('\n')
        print("-" * width)
        print('{:^{width}}'.format("IOOS Compliance Checker Report", width=width))
        print('{:^{width}}'.format(check_name, width=width))
        print('{:^{width}}'.format(check_url, width=width))
        print("-" * width)
        if issue_count > 0:
            print('{:^{width}}'.format("Corrective Actions", width=width))
            plural = '' if issue_count == 1 else 's'
            print("{} has {} potential issue{}".format(os.path.basename(ds), issue_count, plural))

        return [groups, points, out_of]

    def standard_output_generation(self, groups, limit, points, out_of, check):
        '''
        Generates the Terminal Output
        '''
        if points < out_of:
            self.reasoning_routine(groups, 0, check, priority_flag=limit)
        else:
            print("All tests passed!")

    def reasoning_routine(self, groups, indent, check, line=True, priority_flag=3,
                          _top_level=True):
        """
        print routine performed
        """

        sort_fn = lambda x: x.weight
        groups_sorted = sorted(groups, key=sort_fn, reverse=True)
        result = {key: [v for v in valuesiter if v.value[0] != v.value[1]]
                    for key, valuesiter in itertools.groupby(groups_sorted,
                                                             key=sort_fn)}
        wrapper = textwrap.TextWrapper(initial_indent='',
                                       width=max(int(self.col_width / 2**indent), self.col_width / 2))

        priorities = self.checkers[check]._cc_display_headers
        def process_table(res, check):
            issue = wrapper.fill("{}:".format(res.name))
            if not res.children:
                reason = wrapper.fill('\n'.join(res.msgs))
            else:
                child_reasons = self.reasoning_routine(res.children,
                                                       indent + 1,
                                                       check,
                                                        _top_level=False)
                # there shouldn't be messages if there are children
                # is this a valid assumption?
                reason = "\n{}".format(child_reasons)

            return issue, reason

        # iterate up to the min priority requested
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
                    print('{:^{width}}'.format(level_name, width=width))
                    print("-" * width)

                data_issues = [process_table(res, check) for res in result[level]]

                if _top_level:
                    proc_str = tabulate.tabulate(data_issues, ('Name', 'Reasoning'),
                                        'plain')
                    print(proc_str)
                else:
                    proc_str = tabulate.tabulate(data_issues, tablefmt='plain')
                proc_strs.append(proc_str)
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
            raise ValueError("Unrecognized XML root element: {}".format(xml_doc.tag))
        return ds

    def generate_dataset(self, cdl_path):
        '''
        Use ncgen to generate a netCDF file from a .cdl file
        Returns the path to the generated netcdf file

        :param str cdl_path: Absolute path to cdl file that is used to generate netCDF file
        '''
        if '.cdl' in cdl_path:  # it's possible the filename doesn't have the .cdl extension
            ds_str = cdl_path.replace('.cdl', '.nc')
        else:
            ds_str = cdl_path + '.nc'
        subprocess.call(['ncgen', '-o', ds_str, cdl_path])
        return ds_str

    def load_dataset(self, ds_str):
        """
        Returns an instantiated instance of either a netCDF file or an SOS
        mapped DS object.

        :param str ds_str: URL of the resource to load
        """
        # If it's a remote URL load it as a remote resource, otherwise treat it
        # as a local resource.
        pr = urlparse(ds_str)
        if pr.netloc:
            return self.load_remote_dataset(ds_str)
        return self.load_local_dataset(ds_str)

    def load_remote_dataset(self, ds_str):
        '''
        Returns a dataset instance for the remote resource, either OPeNDAP or SOS

        :param str ds_str: URL to the remote resource
        '''

        if opendap.is_opendap(ds_str):
            return Dataset(ds_str)
        else:
            # Check if the HTTP response is XML, if it is, it's likely SOS so
            # we'll attempt to parse the response as SOS
            response = requests.get(ds_str, allow_redirects=True)
            if 'text/xml' in response.headers['content-type']:
                return self.process_doc(response.content)

            raise ValueError("Unknown service with content-type: {}".format(response.headers['content-type']))

    def load_local_dataset(self, ds_str):
        '''
        Returns a dataset instance for the local resource

        :param ds_str: Path to the resource
        '''
        if cdl.is_cdl(ds_str):
            ds_str = self.generate_dataset(ds_str)

        if netcdf.is_netcdf(ds_str):
            return MemoizedDataset(ds_str)

        # Assume this is just a Generic File if it exists
        if os.path.isfile(ds_str):
            return GenericFile(ds_str)

        raise ValueError("File is an unknown format")

    def scores(self, raw_scores):
        """
        Transforms raw scores from a single checker into a fully tallied and grouped scoreline.
        """
        grouped = self._group_raw(raw_scores)

        return (grouped)

    def _group_raw(self, raw_scores, cur=None, level=1):
        """
        Internal recursive method to group raw scores into a cascading score summary.

        Only top level items are tallied for scores.
        """

        def build_group(label=None, weight=None, value=None, sub=None):
            label = label
            weight = weight
            value = self._translate_value(value)
            sub = sub or []

            return Result(weight=weight,
                          value=value,
                          name=label,
                          children=sub)

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
            Slices off first element (if list/tuple) of classification or just returns it if scalar.
            """
            if isinstance(r.name, tuple) or isinstance(r.name, list):
                if len(r.name) == 0:
                    retval = ''
                else:
                    retval = r.name[0:1][0]
            else:
                retval = r.name
            return retval

        grouped = itertools.groupby(sorted(raw_scores, key=group_func),
                                    key=group_func)

        ret_val = []

        for k, v in grouped:

            v = list(v)

            cv = self._group_raw(list(map(trim_groups, v)), k, level + 1)
            if len(cv):
                # if this node has children, max weight of children + sum of all the scores
                max_weight = max([x.weight for x in cv])
                sum_scores = tuple(map(sum, list(zip(*([x.value for x in cv])))))
                msgs = []
            else:
                max_weight = max([x.weight for x in v])
                sum_scores = tuple(map(sum, list(zip(*([self._translate_value(x.value) for x in v])))))
                msgs = sum([x.msgs for x in v], [])

            ret_val.append(Result(name=k, weight=max_weight, value=sum_scores, children=cv, msgs=msgs))

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
