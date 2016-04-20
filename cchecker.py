#!/usr/bin/env python
from __future__ import print_function

import argparse
import sys
from compliance_checker.runner import ComplianceChecker, CheckSuite
from compliance_checker import __version__
from ordered_set import OrderedSet


def main():
    # Load all available checker classes
    check_suite = CheckSuite()
    check_suite.load_all_available_checkers()

    class SupportedTests(argparse.Action):
        """
A custom argparse Action to ensure only the supported tests are passed in
        """
        def __call__(self, parser, namespace, values, option_string=None):
            choices_sorted = sorted(check_suite.checkers.keys())
            choices = set(check_suite.checkers.keys())
            specified = OrderedSet(values.split(','))
            arg_diff = specified - choices
            if len(specified - choices) != 0:
                err_msg = "Test arguments {} not in valid values of {}".format(list(arg_diff), choices_sorted)
                raise argparse.ArgumentError(self, err_msg)
            else:
                setattr(namespace, self.dest, specified)


    parser = argparse.ArgumentParser()
    parser.add_argument('--test', '-t', '--test=', '-t=', required=True,
                        #nargs='+',
                        action=SupportedTests,
                        help="Select the Checks you want to perform.")
                        #choices=sorted(check_suite.checkers.keys()))

    parser.add_argument('--criteria', '-c',
                        help="Define the criteria for the checks.  Either Strict, Normal, or Lenient.  Defaults to Normal.",
                        nargs='?', default='normal',
                        choices = ['lenient', 'normal', 'strict'])

    parser.add_argument('--verbose', '-v',
                        help="Increase output. May be specified up to three times.",
                        action="count")

    parser.add_argument('-f', '--format', default='text',
                        choices=['text', 'html', 'json'], help='Output format')
    parser.add_argument('-o', '--output', default='-', action='store',
                        help='Output filename')
    parser.add_argument('-V', '--version', action='store_true',
                        help='Display the IOOS Compliance Checker version information.')
    parser.add_argument('dataset_location', nargs='*',
                        help="Defines the location of the dataset to be checked.")

    args = parser.parse_args()
    args.test = args.test or ['acdd']

    if args.version:
        print("IOOS compliance checker version %s" % __version__)
        return 0

    return_values = []
    had_errors = []
    for dataset in args.dataset_location:
        if args.format != 'json':
            print("Running Compliance Checker on the dataset from: {}".format(dataset), file=sys.stderr)
        return_value, errors = ComplianceChecker.run_checker(args.dataset_location[0],
                                                             args.test,
                                                             args.verbose,
                                                             args.criteria,
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
