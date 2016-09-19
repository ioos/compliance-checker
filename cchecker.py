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
                        help="Select the Checks you want to perform.  Defaults to 'acdd' if unspecified")

    parser.add_argument('--criteria', '-c',
                        help="Define the criteria for the checks.  Either Strict, Normal, or Lenient.  Defaults to Normal.",
                        nargs='?', default='normal',
                        choices=['lenient', 'normal', 'strict'])

    parser.add_argument('--verbose', '-v',
                        help="Increase output. May be specified up to three times.",
                        action="count",
                        default=0)

    parser.add_argument('--skip-checks', '-s',
                        help="Specifies tests to skip",
                        action='append')

    parser.add_argument('-f', '--format', default='text',
                        choices=['text', 'html', 'json'], help='Output format')
    parser.add_argument('-o', '--output', default='-', action='store',
                        help='Output filename')
    parser.add_argument('-V', '--version', action='store_true',
                        help='Display the IOOS Compliance Checker version information.')
    parser.add_argument('dataset_location', nargs='*',
                        help="Defines the location of the dataset to be checked.")
    parser.add_argument('-l', '--list-tests', action='store_true', help='List the available tests')
    parser.add_argument('-d', '--download-standard-names', help='Specify a version of the cf standard name table to download as packaged version')

    args = parser.parse_args()

    if args.version:
        print("IOOS compliance checker version %s" % __version__)
        return 0

    if args.list_tests:
        print("IOOS compliance checker available checker suites (code version):")
        for checker in sorted(check_suite.checkers.keys()):
            version = getattr(check_suite.checkers[checker], '_cc_checker_version', "???")
            print(" - {} ({})".format(checker, version))
        return 0

    if args.download_standard_names:
        download_cf_standard_name_table(args.download_standard_names)

    return_values = []
    had_errors = []
    for dataset in args.dataset_location:
        if args.format != 'json':
            print("Running Compliance Checker on the dataset from: {}".format(dataset), file=sys.stderr)
        return_value, errors = ComplianceChecker.run_checker(args.dataset_location[0],
                                                             args.test or ['acdd'],
                                                             args.verbose,
                                                             args.criteria,
                                                             args.skip_checks,
                                                             args.output,
                                                             args.format)
        return_values.append(return_value)
        had_errors.append(errors)

    if any(had_errors):
        return 2
    if all(return_values):
        return 0
    return 1

if __name__ == "__main__":
    sys.exit(main())
