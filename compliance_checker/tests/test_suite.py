from pkg_resources import resource_filename
from compliance_checker.suite import CheckSuite
from compliance_checker.base import Result, BaseCheck, GenericFile
import numpy as np
import unittest
import os

static_files = {
    '2dim'                       : resource_filename('compliance_checker', 'tests/data/2dim-grid.nc'),
    'bad_region'                 : resource_filename('compliance_checker', 'tests/data/bad_region.nc'),
    'bad_data_type'              : resource_filename('compliance_checker', 'tests/data/bad_data_type.nc'),
    'test_cdl'                   : resource_filename('compliance_checker', 'tests/data/test_cdl.cdl'),
    'test_cdl_nc'                : resource_filename('compliance_checker', 'tests/data/test_cdl_nc_file.nc'),
    'empty'                    : resource_filename('compliance_checker', 'tests/data/non-comp/empty.file'),
}


class TestSuite(unittest.TestCase):
    # @see
    # http://www.saltycrane.com/blog/2012/07/how-prevent-nose-unittest-using-docstring-when-verbosity-2/

    def setUp(self):
        self.cs = CheckSuite()
        self.cs.load_all_available_checkers()

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
            return "%s ( %s )" % (name[-1], '.'.join(name[:-2]) + ":" +
                                  '.'.join(name[-2:]))
    __str__ = __repr__

    def test_suite(self):
        # BWA: what's the purpose of this test?  Just to see if the suite
        # runs without errors?
        ds = self.cs.load_dataset(static_files['2dim'])
        self.cs.run(ds, 'acdd')

    def test_unicode_formatting(self):
        ds = self.cs.load_dataset(static_files['bad_region'])
        score_groups = self.cs.run(ds, 'cf')

        limit = 2
        for checker, rpair in score_groups.items():
            groups, errors = rpair
            score_list, points, out_of = self.cs.standard_output(ds.filepath(),
                                                            limit, checker,
                                                            groups)
            # This asserts that print is able to generate all of the unicode output
            self.cs.standard_output_generation(groups, limit, points, out_of, checker)

    def test_skip_checks(self):
        """Tests that checks are properly skipped when specified"""
        ds = self.cs.load_dataset(static_files['2dim'])
        # exclude title from the check attributes
        score_groups = self.cs.run(ds, ['check_high'], 'acdd')
        assert all(sg.name not in {'Conventions', 'title', 'keywords',
                                   'summary'} for sg in score_groups['acdd'][0])

    def test_group_func(self):
        # This is checking for issue #183, where group_func results in
        # IndexError: list index out of range
        ds = self.cs.load_dataset(static_files['bad_data_type'])
        score_groups = self.cs.run(ds, 'cf')

        limit = 2
        for checker, rpair in score_groups.items():
            groups, errors = rpair
            score_list, points, out_of = self.cs.standard_output(ds.filepath(),
                                                            limit, checker,
                                                            groups)
            # This asserts that print is able to generate all of the unicode output
            self.cs.standard_output_generation(groups, limit, points, out_of, checker)

    def test_score_grouping(self):
        # Testing the grouping of results for output, which can fail
        # if some assumptions are not met, e.g. if a Result object has
        # a value attribute of unexpected type
        res = [
            Result(BaseCheck.MEDIUM, True, 'one'),
            Result(BaseCheck.MEDIUM, (1, 3), 'one'),
            Result(BaseCheck.MEDIUM, None, 'one'),
            Result(BaseCheck.MEDIUM, True, 'two'),
            Result(BaseCheck.MEDIUM, np.isnan(1), 'two')  # value is type numpy.bool_
        ]
        score = self.cs.scores(res)
        self.assertEqual(score[0].name, 'one')
        self.assertEqual(score[0].value, (2, 4))
        self.assertEqual(score[1].name, 'two')
        self.assertEqual(score[1].value, (1, 2))

    def test_cdl_file(self):
        # Testing whether you can run compliance checker on a .cdl file
        # Load the cdl file
        ds = self.cs.load_dataset(static_files['test_cdl'])
        vals = self.cs.run(ds, 'cf')

        limit = 2
        for checker, rpair in vals.items():
            groups, errors = rpair
            score_list, cdl_points, cdl_out_of = self.cs.standard_output(ds.filepath(),
                                                                    limit,
                                                                    checker,
                                                                    groups)
            # This asserts that print is able to generate all of the unicode output
            self.cs.standard_output_generation(groups, limit, cdl_points, cdl_out_of, checker)
        ds.close()

        # Ok now load the nc file that it came from
        ds = self.cs.load_dataset(static_files['test_cdl_nc'])
        vals = self.cs.run(ds, 'cf')

        limit = 2
        for checker, rpair in vals.items():
            groups, errors = rpair
            score_list, nc_points, nc_out_of = self.cs.standard_output(ds.filepath(),
                                                                  limit,
                                                                  checker,
                                                                  groups)
            # This asserts that print is able to generate all of the unicode output
            self.cs.standard_output_generation(groups, limit, nc_points, nc_out_of, checker)
        ds.close()

        nc_file_path = static_files['test_cdl'].replace('.cdl', '.nc')
        self.addCleanup(os.remove, nc_file_path)

        # Ok the scores should be equal!
        self.assertEqual(nc_points, cdl_points)
        self.assertEqual(nc_out_of, cdl_out_of)

    def test_load_local_dataset_GenericFile(self):
        resp = self.cs.load_local_dataset(static_files['empty'])
        assert isinstance(resp, GenericFile) ==  True

    def test_standard_output_score_header(self):
        """
        Check that the output score header only checks the number of
        of potential issues, rather than the weighted score
        """
        ds = self.cs.load_dataset(static_files['bad_region'])
        score_groups = self.cs.run(ds, [], 'cf')
        limit = 2
        groups, errors = score_groups['cf']
        score_list, all_passed, out_of = self.cs.standard_output(
                                                        ds.filepath(),
                                                        limit, 'cf',
                                                        groups)

        self.assertEqual((all_passed, out_of), (30, 47))

