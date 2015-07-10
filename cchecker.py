#!/usr/bin/env python
import pprint
import argparse
import sys
from compliance_checker.runner import ComplianceChecker, ComplianceCheckerCheckSuite

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', '-t', '--test=', '-t=', action='append', help= "Select the Checks you want to perform.",  choices=ComplianceCheckerCheckSuite.checkers.keys())
    parser.add_argument('--criteria', '-c', help="Define the criteria for the checks.  Either Strict, Normal, or Lenient.  Defaults to Normal.", nargs='?', default='normal', choices = ['lenient', 'normal', 'strict'])
    parser.add_argument('--verbose' , '-v', help="Increase output. May be specified up to three times.", action="count")
    parser.add_argument('-f', '--format', default='text', choices=['text', 'html'], help='Output format')
    parser.add_argument('-o', '--output', default='-', action='store', help='Output filename')
    parser.add_argument('dataset_location', nargs='+', help= "Defines the location of the dataset to be checked.")

    args = parser.parse_args()
    args.test = args.test or ['acdd']

    return_values = []
    for dataset in args.dataset_location:
        print "Running Compliance Checker on the dataset from: %s" % dataset
        return_value = ComplianceChecker.run_checker(args.dataset_location[0],
                                      args.test,
                                      args.verbose,
                                      args.criteria,
                                      args.output,
                                      args.format)
        return_values.append(return_value)


    if all(return_values):
        return 0
    return 1

if __name__ == "__main__":
    sys.exit(main())
