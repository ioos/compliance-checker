from compliance_checker.acdd import ACDDBaseCheck
from compliance_checker.cf import CFBaseCheck
from compliance_checker.ioos import IOOSBaseCheck
from compliance_checker.suite import CheckSuite

class ComplianceCheckerCheckSuite(CheckSuite):
    """
    CheckSuite that defines all the possible Checker classes for the application.
    """
    checkers = {
        'cf' : CFBaseCheck,
        'acdd' : ACDDBaseCheck,
        'ioos' : IOOSBaseCheck,
    }

class ComplianceChecker(object):
    """
    Compliance Checker runner class.

    Ties together the entire compliance checker framework, is used from
    the command line or can be used via import.
    """
    @classmethod
    def run_checker(cls, ds_loc, checker_names, verbose, criteria):
        """
        Static check runner.

        @param  ds_loc          Dataset location (url or file)
        @param  checker_names    List of string names to run, should match keys of checkers dict (empty list means run all)
        @param  verbose         Verbosity of the output (0, 1, 2)
        @param  criteria        Determines failure (lenient, normal, strict)

        @returns                If the tests failed (based on the criteria)
        """
        retval = True

        cs = ComplianceCheckerCheckSuite()
        ds = cs.load_dataset(ds_loc)
        score_groups = cs.run(ds, *checker_names)

        if criteria == 'normal':
            limit = 2
        elif criteria == 'strict':
            limit = 1
        elif criteria == 'lenient':
            limit = 3

        #Calls output routine to display results in terminal, including scoring.  Goes to verbose function if called by user.
        # @TODO cleanup
        for check_name, groups in score_groups.iteritems():
            score_list, check_number, points, out_of = cs.standard_output(limit, check_name, groups)
            if not verbose:
                cs.non_verbose_output_generation(score_list, limit, points, out_of)
            else:
                cs.verbose_output_generation(groups, limit, points, out_of)

        return cs.passtree(groups, limit)

