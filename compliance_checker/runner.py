from compliance_checker.acdd import ACDDBaseCheck
from compliance_checker.cf import CFBaseCheck
from compliance_checker.ioos import IOOSBaseCheck
from compliance_checker.suite import CheckSuite

class ComplianceChecker(object):
    """
    Compliance Checker runner class.

    Ties together the entire compliance checker framework, is used from
    the command line or can be used via import.
    """
    checkers = {
        'cf' : CFBaseCheck,
        'acdd' : ACDDBaseCheck,
        'ioos' : IOOSBaseCheck,
    }

    @classmethod
    def run_checker(cls, ds_loc, tests_to_run, verbose, criteria):
        """
        Static check runner.

        @param  ds_loc          Dataset location (url or file)
        @param  tests_to_run    List of string names to run, should match keys of checkers dict (empty list means run all)
        @param  verbose         Verbosity of the output (0, 1, 2)
        @param  criteria        Determines failure (lenient, normal, strict)

        @returns                If the tests failed (based on the criteria)
        """
        cs = CheckSuite()
        #if statement to determine if we are running all of the checks
        if not tests_to_run:
            tests_to_run = cls.checkers.keys()

        tests_sent = [v for k,v in cls.checkers.iteritems() if k in tests_to_run]
        fail_flag = cs.run(ds_loc, criteria, tests_to_run, verbose, *tests_sent)
        return fail_flag


