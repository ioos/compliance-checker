from compliance_checker import acdd
from pkg_resources import resource_filename
from compliance_checker.suite import CheckSuite
from compliance_checker.base import Result, BaseCheck
import numpy as np
import unittest
import os

static_files = {
    '2dim'                       : resource_filename('compliance_checker', 'tests/data/2dim-grid.nc'),
    'bad_region'                 : resource_filename('compliance_checker', 'tests/data/bad_region.nc'),
    'bad_data_type'              : resource_filename('compliance_checker', 'tests/data/bad_data_type.nc'),
    'test_cdl'                   : resource_filename('compliance_checker', 'tests/data/test_cdl.cdl'),
    'test_cdl_nc'                : resource_filename('compliance_checker', 'tests/data/test_cdl_nc_file.nc'),
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

    def test_skip_checks(self):
        """Tests that checks are properly skipped when specified"""
        cs = CheckSuite()
        cs.load_all_available_checkers()
        ds = cs.load_dataset(static_files['2dim'])
        # exclude title from the check attributes
        score_groups = cs.run(ds, ['check_high'], 'acdd')
        assert all(sg.name not in {'Conventions', 'title', 'keywords',
                                   'summary'} for sg in score_groups['acdd'][0])

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

    def test_score_grouping(self):
        # Testing the grouping of results for output, which can fail
        # if some assumptions are not met, e.g. if a Result object has
        # a value attribute of unexpected type
        cs = CheckSuite()
        res = [Result(BaseCheck.MEDIUM, True, 'one'),
               Result(BaseCheck.MEDIUM, (1, 3), 'one'),
               Result(BaseCheck.MEDIUM, None, 'one'),
               Result(BaseCheck.MEDIUM, True, 'two'),
               Result(BaseCheck.MEDIUM, np.isnan(1), 'two')  # value is type numpy.bool_
        ]
        score = cs.scores(res)
        self.assertEqual(score[0].name, 'one')
        self.assertEqual(score[0].value, (2, 4))
        self.assertEqual(score[1].name, 'two')
        self.assertEqual(score[1].value, (1, 2))

    def test_cdl_file(self):
        # Testing whether you can run compliance checker on a .cdl file
        cs = CheckSuite()
        cs.load_all_available_checkers()

        # Load the cdl file
        ds = cs.load_dataset(static_files['test_cdl'])
        vals = cs.run(ds, 'cf')

        limit = 2
        for checker, rpair in vals.items():
            groups, errors = rpair
            score_list, cdl_points, cdl_out_of = cs.standard_output(limit, checker, groups)
            # This asserts that print is able to generate all of the unicode output
            cs.non_verbose_output_generation(score_list, groups, limit, cdl_points, cdl_out_of)

        # Ok now load the nc file that it came from
        ds = cs.load_dataset(static_files['test_cdl_nc'])
        vals = cs.run(ds, 'cf')

        limit = 2
        for checker, rpair in vals.items():
            groups, errors = rpair
            score_list, nc_points, nc_out_of = cs.standard_output(limit, checker, groups)
            # This asserts that print is able to generate all of the unicode output
            cs.non_verbose_output_generation(score_list, groups, limit, nc_points, nc_out_of)

        nc_file_path = static_files['test_cdl'].replace('.cdl', '.nc')
        self.addCleanup(os.remove, nc_file_path)

        # Ok the scores should be equal!
        self.assertEqual(nc_points, cdl_points)
        self.assertEqual(nc_out_of, cdl_out_of)
