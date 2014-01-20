#!/usr/bin/env python
import pprint
import argparse
import json


from compliance_checker.suite import CheckSuite
from compliance_checker.acdd import ACDDCheck
from compliance_checker.cf import CFCheck

def main(ds_loc, tests_to_run, verbose, dict_of_checks, json_flag, criteria):
        
    cs = CheckSuite()
    #if statement to determine if we are running all of the checks
    if tests_to_run == 'all':
        tests_to_run = dict_of_checks.keys()
    
    #print "Running Compliance Checker on the dataset from: %s" % ds_loc
    tests_sent = [dict_of_checks[t]() for t in tests_to_run]
    fail_flag, json_dict = cs.run(ds_loc, criteria, tests_to_run, verbose, json_flag, *tests_sent) #Add more arguements if more checks are added
    if json_flag:
        print(json.dumps(json_dict))

if __name__ == "__main__":
    
    #Dictionary of checks
    check_dict = {
                    #CFChecker
                    'cf' : CFCheck,

                    #ACDDChecker
                    'acdd' : ACDDCheck,

                    #Additional Checkers can be added here

                 }

    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_location', nargs = '?', help = "Defines the location of the dataset to be checked.")
    parser.add_argument('--test', '-t', '--test=', '-t=', help = "Select the Checks you want to perform.  Either All, CF, or ACDD.  Defaults to All.", nargs='+', default='all', choices=check_dict.keys())
    parser.add_argument('--criteria', '-c', help = "Define the criteria for the checks.  Either Strict, Normal, or Lenient.  Defaults to Normal.", nargs='?', default='normal', choices = ['lenient', 'normal', 'strict'])
    parser.add_argument('--verbose' , '-v', help = "Increase Output Verbosity", action="count")
    parser.add_argument('--output' , '-o', help = 'Provide Standard or JSON Output.  Note: turning on JSON overrides verbosity', nargs = '+', default = [], choices = ['json'])

    args = parser.parse_args()


    main(args.dataset_location, args.test, args.verbose, check_dict, args.output, args.criteria)

