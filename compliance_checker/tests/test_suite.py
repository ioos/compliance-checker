# coding=utf-8
import os
import unittest
from pathlib import Path

import numpy as np
from pkg_resources import resource_filename

from compliance_checker.base import BaseCheck, GenericFile, Result
from compliance_checker.suite import CheckSuite

static_files = {
    "2dim": resource_filename("compliance_checker", "tests/data/2dim-grid.nc"),
    "bad_region": resource_filename("compliance_checker", "tests/data/bad_region.nc"),
    "bad_data_type": resource_filename(
        "compliance_checker", "tests/data/bad_data_type.nc"
    ),
    "test_cdl": resource_filename("compliance_checker", "tests/data/test_cdl.cdl"),
    "test_cdl_nc": resource_filename(
        "compliance_checker", "tests/data/test_cdl_nc_file.nc"
    ),
    "empty": resource_filename("compliance_checker", "tests/data/non-comp/empty.file"),
    "ru07": resource_filename(
        "compliance_checker", "tests/data/ru07-20130824T170228_rt0.nc"
    ),
    "netCDF4": resource_filename(
        "compliance_checker", "tests/data/test_cdl_nc4_file.cdl"
    ),
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
        name = name.split(".")
        if name[0] not in ["ion", "pyon"]:
            return "%s (%s)" % (name[-1], ".".join(name[:-1]))
        else:
            return "%s ( %s )" % (
                name[-1],
                ".".join(name[:-2]) + ":" + ".".join(name[-2:]),
            )

    __str__ = __repr__

    def test_suite(self):
        # BWA: what's the purpose of this test?  Just to see if the suite
        # runs without errors?
        ds = self.cs.load_dataset(static_files["2dim"])
        self.cs.run(ds, [], "acdd")

    def test_suite_pathlib(self):
        path_obj = Path(static_files["2dim"])
        ds = self.cs.load_dataset(path_obj)
        self.cs.run(ds, [], "acdd")

    def test_unicode_formatting(self):
        ds = self.cs.load_dataset(static_files["bad_region"])
        score_groups = self.cs.run(ds, [], "cf")

        limit = 2
        for checker, rpair in score_groups.items():
            groups, errors = rpair
            score_list, points, out_of = self.cs.standard_output(
                ds.filepath(), limit, checker, groups
            )
            # This asserts that print is able to generate all of the unicode
            # output
            self.cs.standard_output_generation(groups, limit, points, out_of, checker)

    def test_generate_dataset_netCDF4(self):
        """
        Tests that suite.generate_dataset works with cdl file with netCDF4
        features.
        """
        # create netCDF4 file
        ds_name = self.cs.generate_dataset(static_files["netCDF4"])
        # check if correct name is return
        assert ds_name == static_files["netCDF4"].replace(".cdl", ".nc")
        # check if netCDF4 file was created
        assert os.path.isfile(static_files["netCDF4"].replace(".cdl", ".nc"))

    def test_include_checks(self):
        ds = self.cs.load_dataset(static_files["bad_data_type"])
        score_groups = self.cs.run_all(ds, ["cf:1.7"], ["check_standard_name"])
        checks_run = score_groups["cf:1.7"][0]
        assert len(checks_run) == 1
        first_check = checks_run[0]
        assert first_check.name == "§3.3 Standard Name"
        assert first_check.value[0] < first_check.value[1]

    def test_skip_checks(self):
        """Tests that checks are properly skipped when specified"""
        ds = self.cs.load_dataset(static_files["2dim"])
        # exclude title from the check attributes
        score_groups = self.cs.run_all(ds, ["acdd"], skip_checks=["check_high"])
        assert all(
            sg.name not in {"Conventions", "title", "keywords", "summary"}
            for sg in score_groups["acdd"][0]
        )

    def test_skip_check_level(self):
        """Checks level limited skip checks"""
        ds = self.cs.load_dataset(static_files["ru07"])
        score_groups = self.cs.run_all(
            ds,
            ["cf"],
            skip_checks=[
                "check_flags:A",
                "check_convention_possibly_var_attrs:M",
                "check_standard_name:L",
            ],
        )

        name_set = {sg.name for sg in score_groups["cf"][0]}
        # flattened set of messages
        msg_set = {msg for sg in score_groups["cf"][0] for msg in sg.msgs}

        expected_excluded_names = {
            "§3.5 flag_meanings for lat",
            "§3.5 flag_meanings for lon",
            "§3.5 lat is a valid flags variable",
            "§3.5 lat is a valid flags variable",
            "§3.5 lon is a valid flags variable",
        }

        self.assertTrue(len(expected_excluded_names & name_set) == 0)

        # should skip references
        ref_msg = "references global attribute should be a non-empty string"
        self.assertTrue(ref_msg not in msg_set)
        # check_standard_name is high priority, but we requested only low,
        # so the standard_name check should still exist
        standard_name_hdr = "§3.3 Standard Name"
        self.assertTrue(standard_name_hdr in name_set)

    def test_group_func(self):
        # This is checking for issue #183, where group_func results in
        # IndexError: list index out of range
        ds = self.cs.load_dataset(static_files["bad_data_type"])
        score_groups = self.cs.run(ds, [], "cf")

        limit = 2
        for checker, rpair in score_groups.items():
            groups, errors = rpair
            score_list, points, out_of = self.cs.standard_output(
                ds.filepath(), limit, checker, groups
            )
            # This asserts that print is able to generate all of the unicode output
            self.cs.standard_output_generation(groups, limit, points, out_of, checker)

    def test_score_grouping(self):
        # Testing the grouping of results for output, which can fail
        # if some assumptions are not met, e.g. if a Result object has
        # a value attribute of unexpected type
        res = [
            Result(BaseCheck.MEDIUM, True, "one"),
            Result(BaseCheck.MEDIUM, (1, 3), "one"),
            Result(BaseCheck.MEDIUM, None, "one"),
            Result(BaseCheck.MEDIUM, True, "two"),
            Result(BaseCheck.MEDIUM, np.isnan(1), "two"),  # value is type numpy.bool_
        ]
        score = self.cs.scores(res)
        self.assertEqual(score[0].name, "one")
        self.assertEqual(score[0].value, (2, 4))
        self.assertEqual(score[1].name, "two")
        self.assertEqual(score[1].value, (1, 2))

    def test_cdl_file(self):
        # Testing whether you can run compliance checker on a .cdl file
        # Load the cdl file
        ds = self.cs.load_dataset(static_files["test_cdl"])
        vals = self.cs.run(ds, [], "cf")

        limit = 2
        for checker, rpair in vals.items():
            groups, errors = rpair
            score_list, cdl_points, cdl_out_of = self.cs.standard_output(
                ds.filepath(), limit, checker, groups
            )
            # This asserts that print is able to generate all of the unicode output
            self.cs.standard_output_generation(
                groups, limit, cdl_points, cdl_out_of, checker
            )
        ds.close()

        # Ok now load the nc file that it came from
        ds = self.cs.load_dataset(static_files["test_cdl_nc"])
        vals = self.cs.run(ds, [], "cf")

        limit = 2
        for checker, rpair in vals.items():
            groups, errors = rpair
            score_list, nc_points, nc_out_of = self.cs.standard_output(
                ds.filepath(), limit, checker, groups
            )
            # This asserts that print is able to generate all of the unicode output
            self.cs.standard_output_generation(
                groups, limit, nc_points, nc_out_of, checker
            )
        ds.close()

        nc_file_path = static_files["test_cdl"].replace(".cdl", ".nc")
        self.addCleanup(os.remove, nc_file_path)

        # Ok the scores should be equal!
        self.assertEqual(nc_points, cdl_points)
        self.assertEqual(nc_out_of, cdl_out_of)

    def test_load_local_dataset_GenericFile(self):
        resp = self.cs.load_local_dataset(static_files["empty"])
        assert isinstance(resp, GenericFile)

    def test_standard_output_score_header(self):
        """
        Check that the output score header only checks the number of
        of potential issues, rather than the weighted score
        """
        ds = self.cs.load_dataset(static_files["bad_region"])
        score_groups = self.cs.run(ds, [], "cf")
        limit = 2
        groups, errors = score_groups["cf"]
        score_list, all_passed, out_of = self.cs.standard_output(
            ds.filepath(), limit, "cf", groups
        )
        assert all_passed < out_of

    def test_netCDF4_features(self):
        """
        Check if a proper netCDF4 file with netCDF4-datatypes is created.
        """
        # create and open dataset
        ds = self.cs.load_dataset(static_files["netCDF4"])
        # check if netCDF type of global attributes is correct
        assert isinstance(ds.global_att_of_type_int, np.int32)
        # check if netCDF4 type of global attributes is correct
        assert isinstance(ds.global_att_of_type_int64, np.int64)
        # check if netCDF type of variable is correct
        assert ds["tas"].dtype is np.dtype("float32")
        # check if netCDF4 type of variable is correct
        assert ds["mask"].dtype is np.dtype("int64")
