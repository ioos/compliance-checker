import traceback

from compliance_checker.acdd import ACDDBaseCheck
from compliance_checker.cf import CFBaseCheck
from compliance_checker.ioos import IOOSBaseCheck
from compliance_checker.suite import CheckSuite
from compliance_checker.glider_dac import GliderCheck

class ComplianceCheckerCheckSuite(CheckSuite):
    """
    CheckSuite that defines all the possible Checker classes for the application.
    """
    checkers = {
        'cf'        : CFBaseCheck,
        'acdd'      : ACDDBaseCheck,
        'ioos'      : IOOSBaseCheck,
        'gliderdac' : GliderCheck
    }

class ComplianceChecker(object):
    """
    Compliance Checker runner class.

    Ties together the entire compliance checker framework, is used from
    the command line or can be used via import.
    """
    @classmethod
    def run_checker(cls, ds_loc, checker_names, verbose, criteria, output_filename='stdout', output_format='stdout'):
        """
        Static check runner.

        @param  ds_loc          Dataset location (url or file)
        @param  checker_names    List of string names to run, should match keys of checkers dict (empty list means run all)
        @param  verbose         Verbosity of the output (0, 1, 2)
        @param  criteria        Determines failure (lenient, normal, strict)
        @param  output_filename Path to the file for output
        @param  output_format   Format of the output

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

        if output_filename == 'stdout' and output_format == 'stdout':
            groups = cls.stdout_output(cs, score_groups, verbose, limit)

        elif output_filename != 'stdout' and output_format == 'html':
            groups = cls.html_output(cs, score_groups, output_filename, verbose, limit, ds_loc)

        else:
            raise TypeError('Invalid format %s' % output_format)

        return cs.passtree(groups, limit)

    @classmethod
    def stdout_output(cls, cs, score_groups, verbose, limit):
        '''
        Calls output routine to display results in terminal, including scoring.
        Goes to verbose function if called by user.

        :param list score_groups: Score as returned by the compliance checker suite
        '''
        for checker, rpair in score_groups.iteritems():
            groups, errors = rpair

            if len(errors):
                print "The following exceptions occured during the %s checker (possibly indicate compliance checker issues):" % checker

                for check_name, epair in errors.iteritems():
                    print "%s.%s: %s" % (checker, check_name, epair[0].message)
                    if verbose > 0:
                        traceback.print_tb(epair[1].tb_next.tb_next)    # skip first two as they are noise from the running itself @TODO search for check_name
                        print

            score_list, points, out_of = cs.standard_output(limit, checker, groups)
            if not verbose:
                cs.non_verbose_output_generation(score_list, groups, limit, points, out_of)
            else:
                cs.verbose_output_generation(groups, limit, points, out_of)
        return groups

    @classmethod
    def html_output(cls, cs, score_groups, output_filename, verbose, limit, ds_loc):
        '''
        '''
        for checker, rpair in score_groups.iteritems():
            groups, errors = rpair
            cs.html_output(limit, checker, groups, output_filename, ds_loc)
        return groups

