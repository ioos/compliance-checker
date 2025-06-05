#!/usr/bin/env python

import argparse
import sys
import warnings
from collections import defaultdict
from textwrap import dedent

from compliance_checker import __version__
from compliance_checker.cf.util import download_cf_standard_name_table
from compliance_checker.runner import CheckSuite, ComplianceChecker


def _print_checker_name_header(checker_str):
    """
    Helper function to prints a checker name surrounded by a border of "="
    :param checker_suite: A check suite string name
    :type checker: str
    """
    print("{0}\n {1} \n{0}".format("=" * (len(checker_str) + 2), checker_str))


def parse_options(opts):
    """
    Helper function to parse possible options. Splits option into key/value
    pairs and optionally a value for the checker option. The separator
    is a colon.

    :param opts: Iterable of strings with options
    :rtype: dict
    :return: Dictionary with keys as checker type (i.e. "cf", "acdd").
             Each value is a dictionary where keys are checker options and values
             are checker option values or None if not provided.
    """
    options_dict = defaultdict(dict)
    for opt_str in opts:
        try:
            checker_type, checker_opt, *checker_val = opt_str.split(":", 2)
            checker_val = checker_val[0] if checker_val else None
        except ValueError:
            warnings.warn(f"Could not split option {opt_str}, ignoring", stacklevel=2)
        else:
            options_dict[checker_type][checker_opt] = checker_val
    return options_dict


def main():
    # Load all available checker classes
    check_suite = CheckSuite()
    check_suite.load_all_available_checkers()

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test",
        "-t",
        "--test=",
        "-t=",
        default=[],
        action="append",
        help=(
            "Select the Checks you want to perform. Defaults to 'acdd'"
            " if unspecified.  Versions of standards can be specified via "
            "`-t <test_standard>:<version>`.  If `<version>` is omitted, or "
            'is "latest", the latest version of the test standard is used.'
        ),
    )

    parser.add_argument(
        "--criteria",
        "-c",
        help=(
            "Define the criteria for the checks. "
            "Either Strict, Normal, or Lenient.  Defaults to Normal."
        ),
        default="normal",
        choices=["lenient", "normal", "strict"],
    )

    parser.add_argument(
        "--verbose",
        "-v",
        help="Increase output. May be specified up to three times.",
        action="count",
        default=0,
    )

    parser.add_argument(
        "--describe-checks",
        "-D",
        help=(
            "Describes checks for checkers specified using "
            "`-t`. If `-t` is not specified, lists checks "
            "from all available checkers."
        ),
        action="store_true",
    )

    include_exclude = parser.add_mutually_exclusive_group()

    include_exclude.add_argument(
        "--skip-checks",
        "-s",
        help=dedent(
            """
                                    Specifies tests to skip. Can take the form
                                    of either `<check_name>` or
                                    `<check_name>:<skip_level>`.  The first
                                    form skips any checks matching the name.
                                    In the second form <skip_level> may be
                                    specified as "A", "M", or "L".  "A" skips
                                    all checks and is equivalent to calling
                                    the first form. "M" will only show high
                                    priority output from the given check and
                                    will skip medium and low.  "L" will show
                                    both high and medium priority issues, while
                                    skipping low priority issues.  Cannot be
                                    used with `-i`/`--include-checks` option.
                                    """,
        ),
        action="append",
    )

    include_exclude.add_argument(
        "--include-checks",
        "-i",
        help=dedent(
            """
                                    Specifies checks to include. Can only take the form
                                    of `<check_name>`.  Cannot be specified along with
                                    `-s`/`skip_checks`.
                                    """,
        ),
        action="append",
    )

    parser.add_argument(
        "-f",
        "--format",
        default=[],
        action="append",
        help=(
            "Output format(s). Options are 'text', 'html', 'json', 'json_new'."
            " The difference between the 'json' and the 'json_new'"
            " formats is that the 'json' format has the check as the top level"
            " key, whereas the 'json_new' format has the dataset name(s) as the"
            " main key in the output follow by any checks as subkeys.  Also, "
            "'json' format can be only be run against one input file, whereas "
            "'json_new' can be run against multiple files."
        ),
        choices=["text", "html", "json", "json_new"],
    )

    parser.add_argument(
        "-o",
        "--output",
        default=[],
        action="append",
        help=(
            "Output filename(s).  If '-' is supplied, output to stdout."
            " Can either be one or many files.  If one file is supplied,"
            " but the checker is run against many files, all the output"
            " from the checks goes to that file (does not presently work "
            "with 'json' format).  If more than one output file is "
            "supplied, the number of input datasets supplied must match "
            "the number of output files."
        ),
    )

    parser.add_argument(
        "-O",
        "--option",
        default=[],
        action="append",
        help=dedent(
            """
                                    Additional options to be passed to the
                                    checkers.  Multiple options can be specified
                                    via multiple invocations of this switch.
                                    Options should be prefixed with a the
                                    checker name followed by the option,
                                    potentially followed by a value, e.g.
                                    '<checker>:<option_name>[:<option_value>]'

                                    Available options:
                                    'cf:enable_appendix_a_checks' - Allow check
                                    results against CF Appendix A for attribute
                                    location and data types.
                                    """,
        ),
    )

    parser.add_argument(
        "-V",
        "--version",
        action="store_true",
        help="Display the IOOS Compliance Checker version information.",
    )

    parser.add_argument(
        "dataset_location",
        nargs="*",
        help=(
            "Defines the location of the dataset to be checked. The location "
            "can be a local netCDF file, a remote OPeNDAP endpoint, a remote "
            "netCDF file which returns content-type header of "
            "'application/x-netcdf', or an ERDDAP TableDAP endpoint. "
            "Note that the ERDDAP TableDAP endpoint will currently attempt "
            "to fetch the entire TableDAP dataset."
        ),
    )

    parser.add_argument(
        "-l",
        "--list-tests",
        action="store_true",
        help="List the available tests",
    )

    parser.add_argument(
        "-d",
        "--download-standard-names",
        help=(
            "Specify a version of the cf standard name table"
            " to download as packaged version. Either specify"
            ' a version number (e.g. "72") to fetch a '
            'specific version or "latest" to get the '
            "latest CF standard name table."
        ),
    )

    # Add command line args from generator plugins
    check_suite.add_plugin_args(parser)

    args = parser.parse_args()

    check_suite.load_generated_checkers(args)

    if args.version:
        print(f"IOOS compliance checker version {__version__}")
        sys.exit(0)

    options_dict = parse_options(args.option) if args.option else defaultdict(dict)

    if args.describe_checks:
        error_stat = 0
        if args.test:
            checker_names = set(args.test)
        else:
            # skip "latest" meta-versions (":latest" or no explicit version
            # specifier)
            checker_names = [
                c
                for c in check_suite.checkers
                if ":" in c and not c.endswith(":latest")
            ]

        for checker_name in sorted(checker_names):
            if checker_name not in check_suite.checkers:
                print(
                    f"Cannot find checker '{checker_name}' with which to "
                    "describe checks",
                    file=sys.stderr,
                )
                error_stat = 1
            else:
                _print_checker_name_header(checker_name)
                check_suite._print_checker(check_suite.checkers[checker_name])
        sys.exit(error_stat)

    if args.list_tests:
        print("IOOS compliance checker available checker suites:")
        check_suite._print_suites(args.verbose)
        return 0

    if args.download_standard_names:
        download_cf_standard_name_table(args.download_standard_names)

    if len(args.dataset_location) == 0:
        parser.print_help()
        sys.exit(1)

    # Check the number of output files
    if not args.output:
        args.output = "-"
    output_len = len(args.output)
    if not (output_len == 1 or output_len == len(args.dataset_location)):
        print(
            "The number of output files must either be one or the same as the number of datasets",
            file=sys.stderr,
        )
        sys.exit(2)

    # Run the compliance checker
    # 2 modes, concatenated output file or multiple output files
    return_values = []
    had_errors = []
    if output_len == 1:
        if args.format != "json":
            print(
                f"Running Compliance Checker on the datasets from: {args.dataset_location}",
                file=sys.stderr,
            )
        return_value, errors = ComplianceChecker.run_checker(
            args.dataset_location,
            args.test or ["acdd"],
            args.verbose,
            args.criteria,
            args.skip_checks,
            args.include_checks,
            args.output[0],
            args.format or ["text"],
            options=options_dict,
        )
        return_values.append(return_value)
        had_errors.append(errors)
    else:
        for output, dataset in zip(args.output, args.dataset_location):
            if args.format != "json":
                print(
                    f"Running Compliance Checker on the dataset from: {dataset}",
                    file=sys.stderr,
                )
            return_value, errors = ComplianceChecker.run_checker(
                [dataset],
                args.test or ["acdd"],
                args.verbose,
                args.criteria,
                args.skip_checks,
                args.include_checks,
                output,
                args.format or ["text"],
                options=options_dict,
            )
            return_values.append(return_value)
            had_errors.append(errors)

    if any(had_errors):
        sys.exit(2)
    if all(return_values):
        sys.exit(0)
    sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
