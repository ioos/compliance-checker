import json
import os
import sys
import traceback
from collections import OrderedDict
from contextlib import contextmanager

from compliance_checker.suite import CheckSuite


# Py 3.4+ has contextlib.redirect_stdout to redirect stdout to a different
# stream, but use this decorated function in order to redirect output in
# previous versions
@contextmanager
def stdout_redirector(stream):
    old_stdout = sys.stdout
    sys.stdout = stream
    try:
        yield
    finally:
        sys.stdout = old_stdout


class ComplianceChecker:
    """
    Compliance Checker runner class.

    Ties together the entire compliance checker framework, is used from
    the command line or can be used via import.
    """

    # Consider using __init__ instead of so many classmethods
    @classmethod
    def run_checker(
        cls,
        ds_loc,
        checker_names,
        verbose,
        criteria,
        skip_checks=None,
        include_checks=None,
        output_filename="-",
        output_format="text",
        options=None,
    ):
        """
        Static check runner.

        @param  ds_loc          Dataset location (url or file)
        @param  checker_names   List of string names to run, should match keys of checkers dict (empty list means run all)
        @param  verbose         Verbosity of the output (0, 1, 2)
        @param  criteria        Determines failure (lenient, normal, strict)
        @param  output_filename Path to the file for output
        @param  skip_checks     Names of checks to skip
        @param  include_checks  Names of checks to include
        @param  output_format   Format of the output(s)

        @returns                If the tests failed (based on the criteria)
        """
        all_groups = []
        cs = CheckSuite(options=options or {})
        # using OrderedDict is important here to preserve the order
        # of multiple datasets which may be passed in
        score_dict = OrderedDict()
        if not isinstance(ds_loc, str):
            locs = ds_loc
        # if single dataset, put in list
        else:
            locs = [ds_loc]

        # Make sure output format is a list
        if isinstance(output_format, str):
            output_format = [output_format]

        for loc in locs:  # loop through each dataset and run specified checks
            ds = cs.load_dataset(loc)

            score_groups = cs.run_all(ds, checker_names, include_checks, skip_checks)
            for group in score_groups.values():
                all_groups.append(group[0])
            # TODO: consider wrapping in a proper context manager instead
            if hasattr(ds, "close"):
                ds.close()

            if not score_groups:
                raise ValueError(
                    "No checks found, please check the name of the checker(s) and that they are installed",
                )
            else:
                score_dict[loc] = score_groups

        # define a score limit to truncate the output to the strictness level
        # specified by the user
        if criteria == "normal":
            limit = 2
        elif criteria == "strict":
            limit = 1
        elif criteria == "lenient":
            limit = 3

        for out_fmt in output_format:
            if out_fmt == "text":
                if output_filename == "-":
                    cls.stdout_output(cs, score_dict, verbose, limit)
                # need to redirect output from stdout since print functions are
                # presently used to generate the standard report output
                else:
                    if len(output_format) > 1:
                        # Update file name if needed
                        output_filename = f"{os.path.splitext(output_filename)[0]}.txt"
                    with open(output_filename, "w", encoding="utf-8") as f:
                        with stdout_redirector(f):
                            cls.stdout_output(cs, score_dict, verbose, limit)

            elif out_fmt == "html":
                # Update file name if needed
                if len(output_format) > 1 and output_filename != "-":
                    output_filename = f"{os.path.splitext(output_filename)[0]}.html"
                cls.html_output(cs, score_dict, output_filename, ds_loc, limit)

            elif out_fmt in {"json", "json_new"}:
                # Update file name if needed
                if len(output_format) > 1 and output_filename != "-":
                    output_filename = f"{os.path.splitext(output_filename)[0]}.json"
                cls.json_output(cs, score_dict, output_filename, ds_loc, limit, out_fmt)

            else:
                raise TypeError(f"Invalid format {out_fmt}")

            errors_occurred = cls.check_errors(score_groups, verbose)

        return (
            all(cs.passtree(groups, limit) for groups in all_groups),
            errors_occurred,
        )

    @classmethod
    def stdout_output(cls, cs, score_dict, verbose, limit):
        """
        Calls output routine to display results in terminal, including scoring.
        Goes to verbose function if called by user.

        @param cs           Compliance Checker Suite
        @param score_dict   Dict with dataset name as key, list of results as
                            value
        @param verbose      Integer value for verbosity level
        @param limit        The degree of strictness, 1 being the strictest, and going up from there.
        """

        for ds, score_groups in score_dict.items():
            for checker, (groups, _errors) in score_groups.items():
                score_list, points, out_of = cs.standard_output(
                    ds,
                    limit,
                    checker,
                    groups,
                )
                # send list of grouped result objects to stdout & reasoning_routine
                cs.standard_output_generation(
                    groups,
                    limit,
                    points,
                    out_of,
                    check=checker,
                )
        return groups

    @classmethod
    def html_output(cls, cs, score_dict, output_filename, ds_loc, limit):
        """
        Generates rendered HTML output for the compliance score(s)
        @param cs              Compliance Checker Suite
        @param score_groups    List of results
        @param output_filename The file path to output to
        @param ds_loc          List of source datasets
        @param limit           The degree of strictness, 1 being the strictest, and going up from there.
        """
        checkers_html = []
        for ds, score_groups in score_dict.items():
            for checker, (groups, _errors) in score_groups.items():
                checkers_html.append(cs.checker_html_output(checker, groups, ds, limit))

        html = cs.html_output(checkers_html)
        if output_filename == "-":
            print(html)
        else:
            with open(output_filename, "w", encoding="utf8") as f:
                f.write(html)

        return groups

    @classmethod
    def json_output(
        cls,
        cs,
        score_dict,
        output_filename,
        ds_loc,
        limit,
        output_type="json",
    ):
        """
        Generates JSON output for the ocmpliance score(s)
        @param cs              Compliance Checker Suite
        @param score_groups    List of results
        @param output_filename The file path to output to
        @param ds_loc          List of source datasets
        @param limit           The degree of strictness, 1 being the strictest,
                               and going up from there.
        @param output_type     Either 'json' or 'json_new'. json_new is the new
                               json output format that supports multiple datasets
        """
        results = {}
        # json output keys out at the top level by
        if len(score_dict) > 1 and output_type != "json_new":
            raise ValueError(
                "output_type must be set to 'json_new' if outputting multiple datasets to a single json file or stdout",
            )

        if output_type == "json":
            for ds, score_groups in score_dict.items():
                for checker, rpair in score_groups.items():
                    groups, errors = rpair
                    results[checker] = cs.dict_output(
                        checker,
                        groups,
                        ds,
                        limit,
                    )
        elif output_type == "json_new":
            for ds, score_groups in score_dict.items():
                for checker, rpair in score_groups.items():
                    groups, errors = rpair
                    results[ds] = {}
                    results[ds][checker] = cs.dict_output(checker, groups, ds, limit)
        json_results = json.dumps(results, indent=2, ensure_ascii=False)

        if output_filename == "-":
            print(json_results)
        else:
            with open(output_filename, "w", encoding="utf8") as f:
                f.write(json_results)

        return groups

    @classmethod
    def check_errors(cls, score_groups, verbose):
        """
        Reports any errors (exceptions) that occurred during checking to stderr.
        Goes to verbose function if called by user.

        @param score_groups List of results
        @param verbose      Integer value for verbosity level
        """
        errors_occurred = False
        for checker, rpair in score_groups.items():
            errors = rpair[-1]
            if len(errors):
                errors_occurred = True
                print(
                    f"WARNING: The following exceptions occurred during the {checker} checker (possibly indicate compliance checker issues):",
                    file=sys.stderr,
                )
                for check_name, epair in errors.items():
                    print(
                        f"{checker}.{check_name}: {epair[0]}",
                        file=sys.stderr,
                    )

                    if verbose > 0:
                        traceback.print_tb(
                            epair[1].tb_next.tb_next,
                        )  # skip first two as they are noise from the running itself @TODO search for check_name
                        print(file=sys.stderr)

        return errors_occurred
