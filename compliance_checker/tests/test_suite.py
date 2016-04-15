from compliance_checker import acdd
from pkg_resources import resource_filename
from compliance_checker.suite import CheckSuite
import unittest

static_files = {
    '2dim'                       : resource_filename('compliance_checker', 'tests/data/2dim-grid.nc'),
    'bad_region'                 : resource_filename('compliance_checker', 'tests/data/bad_region.nc'),
    'bad_data_type'              : resource_filename('compliance_checker', 'tests/data/bad_data_type.nc'),
}


class TestSuite(unittest.TestCase):
    # @see
    # http://www.saltycrane.com/blog/2012/07/how-prevent-nose-unittest-using-docstring-when-verbosity-2/

    def shortDescription(self):
        return None

    # override __str__ and __repr__ behavior to show a copy-pastable nosetest name for ion tests
    #  ion.module:TestClassName.test_function_name
    def __repr__(self):
        name = self.id()
        name = name.split('.')
        if name[0] not in ["ion", "pyon"]:
            return "%s (%s)" % (name[-1], '.'.join(name[:-1]))
        else:
            return "%s ( %s )" % (name[-1], '.'.join(name[:-2]) + ":" + '.'.join(name[-2:]))
    __str__ = __repr__

    def test_suite(self):
        # BWA: what's the purpose of this test?  Just to see if the suite
        # runs without errors?
        cs = CheckSuite()
        cs.load_all_available_checkers()
        ds = cs.load_dataset(static_files['2dim'])
        vals = cs.run(ds, 'acdd')

        # run no longer returns the summed score, so this test.. just runs
        # assert vals['acdd'][0] == (43.5, 78)

    def test_unicode_formatting(self):
        cs = CheckSuite()
        cs.load_all_available_checkers()
        ds = cs.load_dataset(static_files['bad_region'])
        score_groups = cs.run(ds, 'cf')

        limit = 2
        for checker, rpair in score_groups.items():
            groups, errors = rpair
            score_list, points, out_of = cs.standard_output(limit, checker, groups)
            # This asserts that print is able to generate all of the unicode output
            cs.non_verbose_output_generation(score_list, groups, limit, points, out_of)

    def test_group_func(self):
        # This is checking for issue #183, where group_func results in
        # IndexError: list index out of range
        cs = CheckSuite()
        cs.load_all_available_checkers()
        ds = cs.load_dataset(static_files['bad_data_type'])
        score_groups = cs.run(ds, 'cf')

        limit = 2
        for checker, rpair in score_groups.items():
            groups, errors = rpair
            score_list, points, out_of = cs.standard_output(limit, checker, groups)
            # This asserts that print is able to generate all of the unicode output
            cs.non_verbose_output_generation(score_list, groups, limit, points, out_of)
