#!/usr/bin/env python
from __future__ import print_function

import argparse
import sys
from compliance_checker.runner import ComplianceChecker, CheckSuite
from compliance_checker.cf.util import download_cf_standard_name_table
from compliance_checker import __version__


def main():
    # Load all available checker classes
    check_suite = CheckSuite()
    check_suite.load_all_available_checkers()

    parser = argparse.ArgumentParser()
    parser.add_argument('--test', '-t', '--test=', '-t=', default=[],
                        action='append',
                        help=("Select the Checks you want to perform. Defaults to 'acdd'"
                              " if unspecified.  Versions of standards can be specified via "
                              "`-t <test_standard>:<version>`.  If `<version>` is omitted, or "
                              "is \"latest\", the latest version of the test standard is used."))

    parser.add_argument('--criteria', '-c',
                        help=("Define the criteria for the checks. "
                              "Either Strict, Normal, or Lenient.  Defaults to Normal."),
                        nargs='?', default='normal',
                        choices=['lenient', 'normal', 'strict'])

    parser.add_argument('--verbose', '-v',
                        help="Increase output. May be specified up to three times.",
                        action="count",
                        default=0)

    parser.add_argument('--skip-checks', '-s',
                        help="Specifies tests to skip",
                        action='append')

    parser.add_argument('-f', '--format', default=[], action='append',
                        help=("Output format(s). Options are 'text', 'html', 'json', 'json_new'."
                              " The difference between the 'json' and the 'json_new'"
                              " formats is that the 'json' format has the check as the top level"
                              " key, whereas the 'json_new' format has the dataset name(s) as the"
                              " main key in the output follow by any checks as subkeys.  Also, "
                              "'json' format can be only be run against one input file, whereas "
                              "'json_new' can be run against multiple files."))

    parser.add_argument('-o', '--output', default=[], action='append',
                        help=("Output filename(s).  If '-' is supplied, output to stdout."
                              " Can either be one or many files.  If one file is supplied,"
                              " but the checker is run against many files, all the output"
                              " from the checks goes to that file (does not presently work "
                              "with 'json' format).  If more than one output file is "
                              "supplied, the number of input datasets supplied must match "
                              "the number of output files."))

    parser.add_argument('-V', '--version', action='store_true',
                        help='Display the IOOS Compliance Checker version information.')

    parser.add_argument('dataset_location', nargs='*',
                        help="Defines the location of the dataset to be checked.")

    parser.add_argument('-l', '--list-tests', action='store_true',
                        help='List the available tests')

    parser.add_argument('-d', '--download-standard-names',
                        help=("Specify a version of the cf standard name table"
                              " to download as packaged version"))

    # Add command line args from generator plugins
    check_suite.add_plugin_args(parser)

    args = parser.parse_args()

    check_suite.load_generated_checkers(args)

    if args.version:
        print("IOOS compliance checker version %s" % __version__)
        return 0

    if args.list_tests:
        print("IOOS compliance checker available checker suites:")
        for checker in sorted(check_suite.checkers.keys()):
            version = getattr(check_suite.checkers[checker],
                              '_cc_checker_version', "???")
            if args.verbose:
                print(" - {} (v{})".format(checker, version))
            elif ':' in checker and not checker.endswith(':latest'):  # Skip the "latest" output
                print(" - {}".format(checker))
        return 0

    if args.download_standard_names:
        download_cf_standard_name_table(args.download_standard_names)

    # Check the number of output files
    if not args.output:
        args.output = '-'
    output_len = len(args.output)
    if not (output_len == 1 or output_len == len(args.dataset_location)):
        print('The number of output files must either be one or the same as the number of datasets', file=sys.stderr)
        sys.exit(2)

    # Check the output formats
    format_choices = ['text', 'html', 'json', 'json_new']
    for out_format in args.format:
        if out_format not in format_choices:
            print(("Error: argument -f/--format: invalid choice: '{}'"
                   " (choose from 'text', 'html', 'json', 'json_new')".format(out_format)))
            sys.exit(2)

    # Run the compliance checker
    # 2 modes, concatenated output file or multiple output files
    return_values = []
    had_errors = []
    if output_len == 1:
        if args.format != 'json':
            print("Running Compliance Checker on the datasets from: {}".format(args.dataset_location), file=sys.stderr)
        return_value, errors = ComplianceChecker.run_checker(args.dataset_location,
                                                             args.test or ['acdd'],
                                                             args.verbose,
                                                             args.criteria,
                                                             args.skip_checks,
                                                             args.output[0],
                                                             args.format or ['text'])
        return_values.append(return_value)
        had_errors.append(errors)
    else:
        for output, dataset in zip(args.output, args.dataset_location):
            if args.format != 'json':
                print("Running Compliance Checker on the dataset from: {}".format(dataset), file=sys.stderr)
            return_value, errors = ComplianceChecker.run_checker([dataset],
                                                                args.test or ['acdd'],
                                                                args.verbose,
                                                                args.criteria,
                                                                args.skip_checks,
                                                                output,
                                                                args.format or ['text'])
            return_values.append(return_value)
            had_errors.append(errors)

    if any(had_errors):
        return 2
    if all(return_values):
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
