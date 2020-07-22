#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import os
import sqlite3

from itertools import chain
from tempfile import gettempdir

import numpy as np
import pytest

from netCDF4 import Dataset

from compliance_checker import cfutil
from compliance_checker.cf import (
    CF1_6Check,
    CF1_7Check,
    dimless_vertical_coordinates_1_6,
    dimless_vertical_coordinates_1_7,
)
from compliance_checker.cf.appendix_d import no_missing_terms
from compliance_checker.cf.util import (
    StandardNameTable,
    create_cached_data_dir,
    download_cf_standard_name_table,
    is_time_variable,
    is_vertical_coordinate,
    units_convertible,
    units_temporal,
)
from compliance_checker.suite import CheckSuite
from compliance_checker.tests import BaseTestCase
from compliance_checker.tests.helpers import MockTimeSeries, MockVariable
from compliance_checker.tests.resources import STATIC_FILES


def get_results(results):
    """
    Returns a tuple of the value scored, possible, and a list of messages
    in the result set.
    """
    out_of = 0
    scored = 0
    if isinstance(results, dict):
        results_list = results.values()
    else:
        results_list = results
    for r in results_list:
        if isinstance(r.value, tuple):
            out_of += r.value[1]
            scored += r.value[0]
        else:
            out_of += 1
            scored += int(r.value)

    # Store the messages
    messages = []
    for r in results_list:
        messages.extend(r.msgs)

    return scored, out_of, messages


class TestCF1_6(BaseTestCase):
    def setUp(self):
        """Initialize a CF1_6Check object."""

        self.cf = CF1_6Check()

    # --------------------------------------------------------------------------------
    # Helper Methods
    # --------------------------------------------------------------------------------

    def new_nc_file(self):
        """
        Make a new temporary netCDF file for the scope of the test
        """
        nc_file_path = os.path.join(gettempdir(), "example.nc")
        if os.path.exists(nc_file_path):
            raise IOError("File Exists: %s" % nc_file_path)
        nc = Dataset(nc_file_path, "w")
        self.addCleanup(os.remove, nc_file_path)
        self.addCleanup(nc.close)
        return nc

    def test_coord_data_vars(self):
        """Check that coordinate data variables are properly handled"""
        ds = MockTimeSeries()
        ds.createDimension("siglev", 20)

        temp = ds.createVariable(
            "temp",
            np.float64,
            dimensions=("time",),
            fill_value=np.float(99999999999999999999.0),
        )
        temp.coordinates = "sigma noexist"
        ds.createVariable("sigma", np.float64, dimensions=("siglev",))
        self.cf.setup(ds)
        # time is a NUG coordinate variable, sigma is not, but is referred to in
        # variables, so both should show up in cf_coord_data_vars.
        # noexist does not exist in the dataset's variables, so it is not
        # present in coord_data_vars
        self.assertEqual(self.cf.coord_data_vars, {"time", "sigma"})

    def load_dataset(self, nc_dataset):
        """
        Return a loaded NC Dataset for the given path
        """
        if not isinstance(nc_dataset, str):
            raise ValueError("nc_dataset should be a string")

        nc_dataset = Dataset(nc_dataset, "r")
        self.addCleanup(nc_dataset.close)
        return nc_dataset

    # --------------------------------------------------------------------------------
    # Compliance Tests
    # --------------------------------------------------------------------------------

    def test_check_data_types(self):
        """
        Invoke check_data_types() and loop through all variables to check data
        types. Pertains to 2.2. The netCDF data types char, byte, short, int,
        float or real, and double are all acceptable. NetCDF4 allows string as
        data type, which is also acceptable.
        """

        # check default netCDF data types
        dataset = self.load_dataset(STATIC_FILES["rutgers"])
        result = self.cf.check_data_types(dataset)
        assert result.value[0] == result.value[1]

        # check if variables of type `string` is properly processed
        dataset = self.load_dataset(STATIC_FILES["string"])
        if dataset.file_format != 'NETCDF4':
            raise RuntimeError(
                "netCDF file of wrong format (not netCDF4) was created for checking"
            )
        result = self.cf.check_data_types(dataset)
        assert result.value[0] == result.value[1]

        # check bad data types
        dataset = self.load_dataset(STATIC_FILES["bad_data_type"])
        result = self.cf.check_data_types(dataset)

        # TODO
        # the acdd_reformat_rebase branch has a new .nc file
        # which constructs the temp variable with an int64 dtype --
        # upon rebasing, this should work as expected
        # assert result.msgs[0] == u'The variable temp failed because the datatype is int64'
        # assert result.value == (6, 7)

    def test_check_child_attr_data_types(self):
        """
        Tests check_child_attr_data_types() to ensure the attributes specified in Section 2.5.1
        have a matching data type to their parent variables."""

        # create dataset using MockDataset (default constructor gives it time dimension)
        ds = MockTimeSeries()
        ds.createVariable(
            "temp", np.float64, dimensions=("time")
        )  # add variable "temp" with dimension "time"

        # check where no special data attrs are present, should result good
        result = self.cf.check_child_attr_data_types(
            ds
        )  # checks all special attrs for all variables
        self.assert_result_is_good(result)

        # delete the dataset and start over to create the variable with _FillValue at time of creation
        del ds
        ds = MockTimeSeries()
        ds.createVariable(
            "temp",
            np.float64,
            dimensions=("time",),
            fill_value=np.float(99999999999999999999.0),
        )

        # give temp _FillValue as a float, expect good result
        result = self.cf.check_child_attr_data_types(ds)
        self.assert_result_is_good(result)

        # give temp valid_range as an array of floats, all should check out
        ds.variables["temp"].setncattr("valid_range", np.array([35.0, 38.0]))
        result = self.cf.check_child_attr_data_types(ds)
        self.assert_result_is_good(result)

        # dimensions would probably not be time for platform,
        # but this makes for an easy sanity check against string-like
        # variables and attributes
        var = ds.createVariable("platform", "S1", dimensions=("time",), fill_value="")

        # this probably doesn't make much sense -- more for _FillValue,
        # but _FillVaue data type checks are done at variable creation time?
        # Can't set manually
        var.setncattr("valid_max", -999)
        result = self.cf.check_child_attr_data_types(ds)
        self.assert_result_is_bad(result)
        # str or bytes should work
        var.setncattr("valid_max", "@")
        result = self.cf.check_child_attr_data_types(ds)
        self.assert_result_is_good(result)
        var.setncattr("valid_max", b"@")
        result = self.cf.check_child_attr_data_types(ds)
        self.assert_result_is_good(result)

        # now give invalid integer for valid_min; above two should still check out, this one should fail
        ds.variables["temp"].setncattr("valid_min", 45)
        result = self.cf.check_child_attr_data_types(ds)
        self.assert_result_is_bad(result)

        # now give invalid string for valid_max
        ds.variables["temp"].setncattr("valid_max", "eighty")
        result = self.cf.check_child_attr_data_types(ds)
        self.assert_result_is_bad(result)

        # TODO for CF-1.7: actual_range, actual_min/max

    def test_appendix_a(self):
        dataset = self.load_dataset(STATIC_FILES["bad_data_type"])
        # Ordinarily, options would be specified in the checker constructor, but
        # we set them manually here so we don't have to monkey patch `setUp`
        self.cf.options = {"enable_appendix_a_checks"}
        new_check = copy.deepcopy(self.cf)
        self.cf.setup(dataset)
        aa_results = self.cf.check_appendix_a(dataset)
        flat_messages = {msg for res in aa_results for msg in res.msgs}
        self.assertIn(
            '[Appendix A] Attribute "compress" should not be present in non-coordinate data (D) variable "temp". This attribute may only appear in coordinate data (C).',
            flat_messages,
        )
        self.assertIn("add_offset must be a numeric type", flat_messages)
        nc_obj = MockTimeSeries()
        nc_obj._FillValue = "-9999.00"
        new_check.setup(nc_obj)
        res2 = new_check.check_appendix_a(nc_obj)
        flat_messages = {msg for res in res2 for msg in res.msgs}
        self.assertIn(
            '[Appendix A] Attribute "_FillValue" should not be present in global (G) attributes. This attribute may only appear in coordinate data (C) and non-coordinate data (D).',
            flat_messages,
        )

    def test_naming_conventions(self):
        """
        Section 2.3 Naming Conventions

        Variable, dimension and attr names should begin with a letter and be composed of letters, digits, and underscores.
        """

        # compliant dataset
        dataset = self.load_dataset(STATIC_FILES["rutgers"])
        results = self.cf.check_naming_conventions(dataset)
        scored, out_of, messages = get_results(results)
        assert scored == out_of

        # non-compliant dataset
        dataset = self.load_dataset(STATIC_FILES["bad"])
        results = self.cf.check_naming_conventions(dataset)
        scored, out_of, messages = get_results(results)
        assert len(results) == 3
        assert scored < out_of
        assert len([r for r in results if r.value[0] < r.value[1]]) == 2
        assert all(r.name == u"§2.3 Naming Conventions" for r in results)

        # another non-compliant dataset
        dataset = self.load_dataset(STATIC_FILES["chap2"])
        results = self.cf.check_naming_conventions(dataset)
        scored, out_of, messages = get_results(results)
        assert len(results) == 3
        assert scored < out_of
        assert len([r for r in results if r.value[0] < r.value[1]]) == 2
        assert all(r.name == u"§2.3 Naming Conventions" for r in results)

    def test_check_names_unique(self):
        """
        2.3 names should not be distinguished purely by case, i.e., if case is disregarded, no two names should be the same.
        """
        dataset = self.load_dataset(STATIC_FILES["rutgers"])
        result = self.cf.check_names_unique(dataset)

        num_var = len(dataset.variables)
        expected = (num_var,) * 2

        self.assertEqual(result.value, expected)

        dataset = self.load_dataset(STATIC_FILES["chap2"])
        result = self.cf.check_names_unique(dataset)
        assert result.value == (6, 7)
        assert (
            result.msgs[0]
            == u"Variables are not case sensitive. Duplicate variables named: not_unique"
        )

    def test_check_dimension_names(self):
        """
        2.4 A variable may have any number of dimensions, including zero, and the dimensions must all have different names.
        """

        dataset = self.load_dataset(STATIC_FILES["bad_data_type"])
        result = self.cf.check_dimension_names(dataset)
        assert result.value == (6, 7)

        dataset = self.load_dataset(STATIC_FILES["chap2"])
        result = self.cf.check_dimension_names(dataset)
        assert result.msgs[0] == u"no_reason has two or more dimensions named time"

    def test_check_dimension_order(self):
        """
        2.4 If any or all of the dimensions of a variable have the interpretations of "date or time" (T), "height or depth" (Z),
        "latitude" (Y), or "longitude" (X) then we recommend, those dimensions to appear in the relative order T, then Z, then Y,
        then X in the CDL definition corresponding to the file. All other dimensions should, whenever possible, be placed to the
        left of the spatiotemporal dimensions.
        """
        dataset = self.load_dataset(STATIC_FILES["bad_data_type"])
        result = self.cf.check_dimension_order(dataset)
        assert result.value == (5, 6)
        assert result.msgs[0] == (
            u"really_bad's dimensions are not in the recommended order "
            "T, Z, Y, X. They are latitude, power"
        )

        dataset = self.load_dataset(STATIC_FILES["dimension_order"])
        result = self.cf.check_dimension_order(dataset)
        self.assertEqual((3, 3), result.value)
        self.assertEqual([], result.msgs)

    def test_check_fill_value_outside_valid_range(self):
        """
        2.5.1 The _FillValue should be outside the range specified by valid_range (if used) for a variable.
        """

        dataset = self.load_dataset(STATIC_FILES["bad_data_type"])
        result = self.cf.check_fill_value_outside_valid_range(dataset)
        assert result.msgs[0] == (
            u"salinity:_FillValue (1.0) should be outside the "
            "range specified by valid_min/valid_max (-10, 10)"
        )

        dataset = self.load_dataset(STATIC_FILES["chap2"])
        result = self.cf.check_fill_value_outside_valid_range(dataset)
        assert result.value == (1, 2)
        assert result.msgs[0] == (
            u"wind_speed:_FillValue (12.0) should be outside the "
            "range specified by valid_min/valid_max (0.0, 20.0)"
        )

    def test_check_conventions_are_cf_16(self):
        """
        §2.6.1 the NUG defined global attribute Conventions to the string value
        "CF-1.6"
        """
        # :Conventions = "CF-1.6"
        dataset = self.load_dataset(STATIC_FILES["rutgers"])
        result = self.cf.check_conventions_version(dataset)
        self.assertTrue(result.value)

        # :Conventions = "CF-1.6 ,ACDD" ;
        dataset = self.load_dataset(STATIC_FILES["conv_multi"])
        result = self.cf.check_conventions_version(dataset)
        self.assertTrue(result.value)

        # :Conventions = "NoConvention"
        dataset = self.load_dataset(STATIC_FILES["conv_bad"])
        result = self.cf.check_conventions_version(dataset)
        self.assertFalse(result.value)
        assert result.msgs[0] == (
            u"§2.6.1 Conventions global attribute does not contain " '"CF-1.6"'
        )

    def test_check_convention_globals(self):
        """
        Load up a dataset and ensure title and history global attrs are checked
        properly (§2.6.2).
        """

        # check for pass
        dataset = self.load_dataset(STATIC_FILES["rutgers"])
        result = self.cf.check_convention_globals(dataset)
        assert result.value[0] == result.value[1]

        # check if it doesn't exist that we pass
        dataset = self.load_dataset(STATIC_FILES["bad_data_type"])
        result = self.cf.check_convention_globals(dataset)
        assert result.value[0] != result.value[1]
        assert (
            result.msgs[0]
            == u"§2.6.2 global attribute title should exist and be a non-empty string"
        )

    def test_check_convention_possibly_var_attrs(self):
        """
        §2.6.2 The units attribute is required for all variables that represent dimensional quantities
        (except for boundary variables defined in Section 7.1, "Cell Boundaries" and climatology variables
        defined in Section 7.4, "Climatological Statistics").

        Units are not required for dimensionless quantities. A variable with no units attribute is assumed
        to be dimensionless. However, a units attribute specifying a dimensionless unit may optionally be
        included.

        - units required
        - type must be recognized by udunits
        - if std name specified, must be consistent with standard name table, must also be consistent with a
          specified cell_methods attribute if present
        """

        dataset = self.load_dataset(STATIC_FILES["rutgers"])
        result = self.cf.check_convention_possibly_var_attrs(dataset)
        # 10x comment attrs
        # 1x institution
        # 1x source
        # 1x EMPTY references
        assert result.value[0] != result.value[1]
        assert (
            result.msgs[0]
            == u"§2.6.2 references global attribute should be a non-empty string"
        )

        # load bad_data_type.nc
        dataset = self.load_dataset(STATIC_FILES["bad_data_type"])
        result = self.cf.check_convention_possibly_var_attrs(dataset)
        # no references
        # institution is a 10L
        # no source

        # comments don't matter unless they're empty

        assert result.value[0] != result.value[1]
        assert (
            result.msgs[0]
            == u"§2.6.2 salinity:institution should be a non-empty string"
        )

    def test_check_standard_name(self):
        """
        3.3 A standard name is associated with a variable via the attribute standard_name which takes a
        string value comprised of a standard name optionally followed by one or more blanks and a
        standard name modifier
        """
        dataset = self.load_dataset(STATIC_FILES["2dim"])
        results = self.cf.check_standard_name(dataset)
        for each in results:
            self.assertTrue(each.value)

        # load failing ds
        dataset = self.load_dataset(STATIC_FILES["bad_data_type"])
        results = self.cf.check_standard_name(dataset)
        score, out_of, messages = get_results(results)

        # 9 vars checked, 8 fail
        assert len(results) == 9
        assert score < out_of
        assert all(r.name == u"§3.3 Standard Name" for r in results)

        # load different ds --  ll vars pass this check
        dataset = self.load_dataset(STATIC_FILES["reduced_horizontal_grid"])
        results = self.cf.check_standard_name(dataset)
        score, out_of, messages = get_results(results)
        assert score == out_of

    def test_cell_bounds(self):
        dataset = self.load_dataset(STATIC_FILES["grid-boundaries"])
        results = self.cf.check_cell_boundaries(dataset)
        score, out_of, messages = get_results(results)
        assert (score, out_of) == (2, 2)

        dataset = self.load_dataset(STATIC_FILES["cf_example_cell_measures"])
        results = self.cf.check_cell_boundaries(dataset)

        dataset = self.load_dataset(STATIC_FILES["bad_data_type"])
        results = self.cf.check_cell_boundaries(dataset)

        dataset = self.load_dataset(STATIC_FILES["bounds_bad_order"])
        results = self.cf.check_cell_boundaries(dataset)
        score, out_of, messages = get_results(results)
        # Make sure that the rgrid coordinate variable isn't checked for standard_name
        assert (score, out_of) == (0, 2)

        dataset = self.load_dataset(STATIC_FILES["bounds_bad_num_coords"])
        results = self.cf.check_cell_boundaries(dataset)
        score, out_of, messages = get_results(results)
        assert (score, out_of) == (0, 2)

        dataset = self.load_dataset(STATIC_FILES["1d_bound_bad"])
        results = self.cf.check_cell_boundaries(dataset)
        score, out_of, messages = get_results(results)
        assert (score, out_of) == (0, 2)

    def test_cell_measures(self):
        dataset = self.load_dataset(STATIC_FILES["cell_measure"])
        results = self.cf.check_cell_measures(dataset)
        score, out_of, messages = get_results(results)
        assert score == out_of
        assert score > 0

        dataset = self.load_dataset(STATIC_FILES["bad_cell_measure1"])
        results = self.cf.check_cell_measures(dataset)
        score, out_of, messages = get_results(results)
        message = (
            "The cell_measures attribute for variable PS is formatted incorrectly.  "
            "It should take the form of either 'area: cell_var' or 'volume: cell_var' "
            "where cell_var is the variable describing the cell measures"
        )
        assert message in messages

        dataset = self.load_dataset(STATIC_FILES["bad_cell_measure2"])
        results = self.cf.check_cell_measures(dataset)
        score, out_of, messages = get_results(results)
        message = u"Cell measure variable box_area referred to by PS is not present in dataset variables"
        assert message in messages

    def test_climatology_cell_methods(self):
        """
        Checks that climatology cell_methods strings are properly validated
        """
        dataset = self.load_dataset(STATIC_FILES["climatology"])
        results = self.cf.check_climatological_statistics(dataset)
        # cell methods in this file is
        # "time: mean within days time: mean over days"
        score, out_of, messages = get_results(results)
        self.assertEqual(score, out_of)
        temp_var = dataset.variables["temperature"] = MockVariable(
            dataset.variables["temperature"]
        )
        temp_var.cell_methods = "INVALID"
        results = self.cf.check_climatological_statistics(dataset)
        score, out_of, messages = get_results(results)
        self.assertNotEqual(score, out_of)
        # incorrect time units
        temp_var.cell_methods = "time: mean within years time: mean over days"
        results = self.cf.check_climatological_statistics(dataset)
        score, out_of, messages = get_results(results)
        self.assertNotEqual(score, out_of)
        # can only have third method over years if first two are within and
        # over days, respectively
        temp_var.cell_methods = (
            "time: mean within years time: mean over years time: sum over years"
        )
        results = self.cf.check_climatological_statistics(dataset)
        score, out_of, messages = get_results(results)
        self.assertNotEqual(score, out_of)
        # this, on the other hand, should work.
        temp_var.cell_methods = (
            "time: mean within days time: mean over days time: sum over years"
        )
        results = self.cf.check_climatological_statistics(dataset)
        score, out_of, messages = get_results(results)
        self.assertEqual(score, out_of)
        # parenthesized comment to describe climatology
        temp_var.cell_methods = (
            "time: sum within days time: maximum over days (ENSO years)"
        )
        results = self.cf.check_climatological_statistics(dataset)
        score, out_of, messages = get_results(results)
        self.assertEqual(score, out_of)

    def test_check_ancillary_variables(self):
        """
        Test to ensure that ancillary variables are properly checked
        """

        dataset = self.load_dataset(STATIC_FILES["rutgers"])
        results = self.cf.check_ancillary_variables(dataset)
        result_dict = {result.name: result for result in results}
        result = result_dict[u"§3.4 Ancillary Data"]
        assert result.value == (2, 2)

        dataset = self.load_dataset(STATIC_FILES["bad_reference"])
        results = self.cf.check_ancillary_variables(dataset)
        result_dict = {result.name: result for result in results}
        result = result_dict[u"§3.4 Ancillary Data"]
        assert result.value == (1, 2)
        assert u"temp_qc is not a variable in this dataset" == result.msgs[0]

    def test_download_standard_name_table(self):
        """
        Test that a user can download a specific standard name table
        """
        version = "35"

        data_directory = create_cached_data_dir()
        location = os.path.join(
            data_directory, "cf-standard-name-table-test-{0}.xml".format(version)
        )
        download_cf_standard_name_table(version, location)

        # Test that the file now exists in location and is the right version
        self.assertTrue(os.path.isfile(location))
        std_names = StandardNameTable(location)
        self.assertEqual(std_names._version, version)
        self.addCleanup(os.remove, location)

    def test_bad_standard_name_table(self):
        """
        Test that failure in case a bad standard name table is passed.
        """
        # would this ever actually be reached by the code?
        with pytest.raises(IOError):
            StandardNameTable("dummy_non_existent_file.ext")

        nc_obj = MockTimeSeries()
        nc_obj.standard_name_table = "dummy_non_existent_file.ext"
        self.assertFalse(self.cf._find_cf_standard_name_table(nc_obj))

        nc_obj.standard_name_table = np.array([], np.float64)
        self.assertFalse(self.cf._find_cf_standard_name_table(nc_obj))

        nc_obj.standard_name_vocabulary = "CF Standard Name Table vNN???"
        with pytest.warns(
            UserWarning,
            match="Cannot extract CF standard name version "
            "number from standard_name_vocabulary string",
        ):
            self.assertFalse(self.cf._find_cf_standard_name_table(nc_obj))

    def test_check_flags(self):
        """Test that the check for flags works as expected."""

        dataset = self.load_dataset(STATIC_FILES["rutgers"])
        results = self.cf.check_flags(dataset)
        scored, out_of, messages = get_results(results)

        # only 4 variables in this dataset do not have perfect scores
        imperfect = [r.value for r in results if r.value[0] < r.value[1]]
        assert len(imperfect) == 4

    def test_check_flag_masks(self):
        dataset = self.load_dataset(STATIC_FILES["ghrsst"])
        results = self.cf.check_flags(dataset)
        scored, out_of, messages = get_results(results)
        # This is an example of a perfect dataset for flags
        assert scored > 0
        assert scored == out_of

    def test_check_bad_units(self):
        """Load a dataset with units that are expected to fail (bad_units.nc).
        There are 6 variables in this dataset, three of which should give
        an error:
            - time, with units "s" (should be <units> since <epoch>)
            - lat, with units "degrees_E" (should be degrees)
            - lev, with units "level" (deprecated)"""

        dataset = self.load_dataset(STATIC_FILES["2dim"])
        results = self.cf.check_units(dataset)
        for result in results:
            self.assert_result_is_good(result)

        # Not sure why bad_data_type was being used, we have a dataset specifically for bad units
        # dataset = self.load_dataset(STATIC_FILES['bad_data_type'])

        dataset = self.load_dataset(STATIC_FILES["bad_units"])
        all_results = self.cf.check_units(dataset)

        # use itertools.chain() to unpack the lists of messages
        results_list = list(chain(*(r.msgs for r in all_results if r.msgs)))

        # check the results only have '§3.1 Units' as the header
        assert all(r.name == u"§3.1 Units" for r in all_results)

        # check that all the expected variables have been hit
        assert all(
            any(s in msg for msg in results_list) for s in ["time", "lat", "lev"]
        )

    def test_latitude(self):
        """
        Section 4.1 Latitude Coordinate
        """
        # Check compliance
        dataset = self.load_dataset(STATIC_FILES["example-grid"])
        results = self.cf.check_latitude(dataset)
        score, out_of, messages = get_results(results)
        assert score == out_of

        # Verify non-compliance -- 9/12 pass
        dataset = self.load_dataset(STATIC_FILES["bad"])
        results = self.cf.check_latitude(dataset)
        scored, out_of, messages = get_results(results)
        assert len(results) == 12
        assert scored < out_of
        assert len([r for r in results if r.value[0] < r.value[1]]) == 3
        assert (r.name == u"§4.1 Latitude Coordinate" for r in results)

        # check with another ds -- all 6 vars checked pass
        dataset = self.load_dataset(STATIC_FILES["rotated_pole_grid"])
        results = self.cf.check_latitude(dataset)
        scored, out_of, messages = get_results(results)
        assert len(results) == 6
        assert scored == out_of
        assert (r.name == u"§4.1 Latitude Coordinate" for r in results)

        # hack to avoid writing to read-only file
        dataset.variables["rlat"] = MockVariable(dataset.variables["rlat"])
        rlat = dataset.variables["rlat"]
        rlat.name = "rlat"
        # test with a bad value
        rlat.units = "degrees_north"
        results = self.cf.check_latitude(dataset)
        scored, out_of, messages = get_results(results)
        wrong_format = u"Grid latitude variable '{}' should use degree equivalent units without east or north components. Current units are {}"
        self.assertTrue(wrong_format.format(rlat.name, rlat.units) in messages)
        rlat.units = "radians"
        results = self.cf.check_latitude(dataset)
        scored, out_of, messages = get_results(results)
        self.assertTrue(wrong_format.format(rlat.name, rlat.units) in messages)

    def test_longitude(self):
        """
        Section 4.2 Longitude Coordinate
        """
        # Check compliance
        dataset = self.load_dataset(STATIC_FILES["example-grid"])
        results = self.cf.check_longitude(dataset)
        score, out_of, messages = get_results(results)
        assert score == out_of

        # Verify non-compliance -- 12 checked, 3 fail
        dataset = self.load_dataset(STATIC_FILES["bad"])
        results = self.cf.check_longitude(dataset)
        scored, out_of, messages = get_results(results)
        assert len(results) == 12
        assert scored < out_of
        assert len([r for r in results if r.value[0] < r.value[1]]) == 3
        assert all(r.name == u"§4.2 Longitude Coordinate" for r in results)

        # check different dataset # TODO can be improved for check_latitude too
        dataset = self.load_dataset(STATIC_FILES["rotated_pole_grid"])
        results = self.cf.check_latitude(dataset)
        scored, out_of, messages = get_results(results)
        assert (scored, out_of) == (6, 6)
        # hack to avoid writing to read-only file
        dataset.variables["rlon"] = MockVariable(dataset.variables["rlon"])
        rlon = dataset.variables["rlon"]
        rlon.name = "rlon"
        # test with a bad value
        rlon.units = "degrees_east"
        results = self.cf.check_longitude(dataset)
        scored, out_of, messages = get_results(results)
        wrong_format = u"Grid longitude variable '{}' should use degree equivalent units without east or north components. Current units are {}"
        self.assertTrue(wrong_format.format(rlon.name, rlon.units) in messages)
        rlon.units = "radians"
        results = self.cf.check_longitude(dataset)
        scored, out_of, messages = get_results(results)
        self.assertTrue(wrong_format.format(rlon.name, rlon.units) in messages)

    def test_is_vertical_coordinate(self):
        """
        Section 4.3 Qualifiers for Vertical Coordinate

        NOTE: The standard doesn't explicitly say that vertical coordinates must be a
        coordinate type.
        """
        # Make something that I can attach attrs to
        mock_variable = MockVariable

        # Proper name/standard_name
        known_name = mock_variable()
        known_name.standard_name = "depth"
        self.assertTrue(is_vertical_coordinate("not_known", known_name))

        # Proper Axis
        axis_set = mock_variable()
        axis_set.axis = "Z"
        self.assertTrue(is_vertical_coordinate("not_known", axis_set))

        # Proper units
        units_set = mock_variable()
        units_set.units = "dbar"
        self.assertTrue(is_vertical_coordinate("not_known", units_set))

        # Proper units/positive
        positive = mock_variable()
        positive.units = "m"
        positive.positive = "up"
        self.assertTrue(is_vertical_coordinate("not_known", positive))

    def test_vertical_dimension(self):
        """
        Section 4.3.1 Dimensional Vertical Coordinate
        """
        # Check for compliance
        dataset = self.load_dataset(STATIC_FILES["example-grid"])
        results = self.cf.check_dimensional_vertical_coordinate(dataset)
        assert len(results) == 1
        assert all(r.name == u"§4.3 Vertical Coordinate" for r in results)

        # non-compliance -- one check fails
        dataset = self.load_dataset(STATIC_FILES["illegal-vertical"])
        results = self.cf.check_dimensional_vertical_coordinate(dataset)
        scored, out_of, messages = get_results(results)
        assert len(results) == 1
        assert all(r.name == u"§4.3 Vertical Coordinate" for r in results)
        assert scored < out_of

    def test_appendix_d(self):
        """
        CF 1.6
        Appendix D
        The definitions given here allow an application to compute dimensional
        coordinate values from the dimensionless ones and associated variables.
        The formulas are expressed for a gridpoint (n,k,j,i) where i and j are
        the horizontal indices, k is the vertical index and n is the time index.
        A coordinate variable is associated with its definition by the value of
        the standard_name attribute. The terms in the definition are associated
        with file variables by the formula_terms attribute. The formula_terms
        attribute takes a string value, the string being comprised of
        blank-separated elements of the form "term: variable", where term is a
        keyword that represents one of the terms in the definition, and variable
        is the name of the variable in a netCDF file that contains the values
        for that term. The order of elements is not significant.
        """

        # For each of the listed dimensionless vertical coordinates,
        # verify that the formula_terms match the provided set of terms
        self.assertTrue(
            no_missing_terms(
                "atmosphere_ln_pressure_coordinate",
                {"p0", "lev"},
                dimless_vertical_coordinates_1_6,
            )
        )
        self.assertTrue(
            no_missing_terms(
                "atmosphere_sigma_coordinate",
                {"sigma", "ps", "ptop"},
                dimless_vertical_coordinates_1_6,
            )
        )
        self.assertTrue(
            no_missing_terms(
                "atmosphere_hybrid_sigma_pressure_coordinate",
                {"a", "b", "ps"},
                dimless_vertical_coordinates_1_6,
            )
        )
        # test alternative terms for
        # 'atmosphere_hybrid_sigma_pressure_coordinate'
        self.assertTrue(
            no_missing_terms(
                "atmosphere_hybrid_sigma_pressure_coordinate",
                {"ap", "b", "ps"},
                dimless_vertical_coordinates_1_6,
            )
        )
        # check that an invalid set of terms fails
        self.assertFalse(
            no_missing_terms(
                "atmosphere_hybrid_sigma_pressure_coordinate",
                {"a", "b", "p"},
                dimless_vertical_coordinates_1_6,
            )
        )
        self.assertTrue(
            no_missing_terms(
                "atmosphere_hybrid_height_coordinate",
                {"a", "b", "orog"},
                dimless_vertical_coordinates_1_6,
            )
        )
        # missing terms should cause failure
        self.assertFalse(
            no_missing_terms(
                "atmosphere_hybrid_height_coordinate",
                {"a", "b"},
                dimless_vertical_coordinates_1_6,
            )
        )
        # excess terms should cause failure
        self.assertFalse(
            no_missing_terms(
                "atmosphere_hybrid_height_coordinate",
                {"a", "b", "c", "orog"},
                dimless_vertical_coordinates_1_6,
            )
        )
        self.assertTrue(
            no_missing_terms(
                "atmosphere_sleve_coordinate",
                {"a", "b1", "b2", "ztop", "zsurf1", "zsurf2"},
                dimless_vertical_coordinates_1_6,
            )
        )
        self.assertTrue(
            no_missing_terms(
                "ocean_sigma_coordinate",
                {"sigma", "eta", "depth"},
                dimless_vertical_coordinates_1_6,
            )
        )
        self.assertTrue(
            no_missing_terms(
                "ocean_s_coordinate",
                {"s", "eta", "depth", "a", "b", "depth_c"},
                dimless_vertical_coordinates_1_6,
            )
        )
        self.assertTrue(
            no_missing_terms(
                "ocean_sigma_z_coordinate",
                {"sigma", "eta", "depth", "depth_c", "nsigma", "zlev"},
                dimless_vertical_coordinates_1_6,
            )
        )
        self.assertTrue(
            no_missing_terms(
                "ocean_double_sigma_coordinate",
                {"sigma", "depth", "z1", "z2", "a", "href", "k_c"},
                dimless_vertical_coordinates_1_6,
            )
        )

    def test_dimensionless_vertical(self):
        """
        Section 4.3.2
        """
        # Check affirmative compliance
        dataset = self.load_dataset(STATIC_FILES["dimensionless"])
        results = self.cf.check_dimensionless_vertical_coordinates(dataset)
        scored, out_of, messages = get_results(results)

        # all variables checked (2) pass
        assert len(results) == 2
        assert scored == out_of
        assert all(r.name == u"§4.3 Vertical Coordinate" for r in results)

        # Check negative compliance -- 3 out of 4 pass

        dataset = self.load_dataset(STATIC_FILES["bad"])
        results = self.cf.check_dimensionless_vertical_coordinates(dataset)
        scored, out_of, messages = get_results(results)
        assert len(results) == 4
        assert scored <= out_of
        assert len([r for r in results if r.value[0] < r.value[1]]) == 2
        assert all(r.name == u"§4.3 Vertical Coordinate" for r in results)

        # test with an invalid formula_terms
        dataset.variables["lev2"] = MockVariable(dataset.variables["lev2"])
        lev2 = dataset.variables["lev2"]
        lev2.formula_terms = "a: var1 b:var2 orog:"

        # create a malformed formula_terms attribute and check that it fails
        # 2/4 still pass
        results = self.cf.check_dimensionless_vertical_coordinates(dataset)
        scored, out_of, messages = get_results(results)

        assert len(results) == 4
        assert scored <= out_of
        assert len([r for r in results if r.value[0] < r.value[1]]) == 2
        assert all(r.name == u"§4.3 Vertical Coordinate" for r in results)

    def test_is_time_variable(self):
        var1 = MockVariable()
        var1.standard_name = "time"
        self.assertTrue(is_time_variable("not_time", var1))

        var2 = MockVariable()
        self.assertTrue(is_time_variable("time", var2))

        self.assertFalse(is_time_variable("not_time", var2))

        var3 = MockVariable()
        var3.axis = "T"
        self.assertTrue(is_time_variable("maybe_time", var3))

        var4 = MockVariable()
        var4.units = "seconds since 1900-01-01"
        self.assertTrue(is_time_variable("maybe_time", var4))

    def test_dimensionless_standard_names(self):
        """Check that dimensionless standard names are properly detected"""
        std_names_xml_root = self.cf._std_names._root
        # canonical_units are K, should be False
        self.assertFalse(
            cfutil.is_dimensionless_standard_name(
                std_names_xml_root, "sea_water_temperature"
            )
        )
        # canonical_units are 1, should be True
        self.assertTrue(
            cfutil.is_dimensionless_standard_name(
                std_names_xml_root, "sea_water_practical_salinity"
            )
        )
        # canonical_units are 1e-3, should be True
        self.assertTrue(
            cfutil.is_dimensionless_standard_name(
                std_names_xml_root, "sea_water_salinity"
            )
        )

    def test_check_time_coordinate(self):
        dataset = self.load_dataset(STATIC_FILES["example-grid"])
        results = self.cf.check_time_coordinate(dataset)
        for r in results:
            self.assertTrue(r.value)

        dataset = self.load_dataset(STATIC_FILES["bad"])
        results = self.cf.check_time_coordinate(dataset)

        scored, out_of, messages = get_results(results)

        assert u"time does not have correct time units" in messages
        assert (scored, out_of) == (1, 2)

    def test_check_calendar(self):
        """Load a dataset with an invalid calendar attribute (non-comp/bad.nc).
        This dataset has a variable, "time" with  calendar attribute "nope"."""

        dataset = self.load_dataset(STATIC_FILES["example-grid"])
        results = self.cf.check_calendar(dataset)
        for r in results:
            self.assertTrue(r.value)

        dataset = self.load_dataset(STATIC_FILES["bad"])
        results = self.cf.check_calendar(dataset)
        scored, out_of, messages = get_results(results)

        assert (
            u"§4.4.1 Variable time should have a valid calendar: 'nope' is not a valid calendar"
            in messages
        )

    def test_check_aux_coordinates(self):
        dataset = self.load_dataset(STATIC_FILES["illegal-aux-coords"])
        results = self.cf.check_aux_coordinates(dataset)
        result_dict = {result.name: result for result in results}
        result = result_dict[u"§5 Coordinate Systems"]
        assert result.msgs == []  # shouldn't have any messages
        assert result.value == (4, 4)

    def test_check_grid_coordinates(self):
        dataset = self.load_dataset(STATIC_FILES["2dim"])
        results = self.cf.check_grid_coordinates(dataset)
        scored, out_of, messages = get_results(results)

        result_dict = {result.name: result for result in results}
        result = result_dict[
            u"§5.6 Horizontal Coorindate Reference Systems, Grid Mappings, Projections"
        ]
        assert result.value == (2, 2)
        assert (scored, out_of) == (2, 2)

    def test_check_two_dimensional(self):
        dataset = self.load_dataset(STATIC_FILES["2dim"])
        results = self.cf.check_grid_coordinates(dataset)
        for r in results:
            self.assertTrue(r.value)
        # Need the bad testing
        dataset = self.load_dataset(STATIC_FILES["bad2dim"])
        results = self.cf.check_grid_coordinates(dataset)
        scored, out_of, messages = get_results(results)

        # all variables checked fail (2)
        assert len(results) == 2
        assert scored < out_of
        assert all(
            r.name
            == u"§5.6 Horizontal Coorindate Reference Systems, Grid Mappings, Projections"
            for r in results
        )

    def test_check_reduced_horizontal_grid(self):
        dataset = self.load_dataset(STATIC_FILES["rhgrid"])
        results = self.cf.check_reduced_horizontal_grid(dataset)
        scored, out_of, messages = get_results(results)
        assert scored == out_of
        assert len(results) == 1
        assert all(r.name == u"§5.3 Reduced Horizontal Grid" for r in results)

        # load failing ds -- one variable has failing check
        dataset = self.load_dataset(STATIC_FILES["bad-rhgrid"])
        results = self.cf.check_reduced_horizontal_grid(dataset)
        scored, out_of, messages = get_results(results)
        assert scored != out_of
        assert len(results) == 2
        assert len([r for r in results if r.value[0] < r.value[1]]) == 1
        assert all(r.name == u"§5.3 Reduced Horizontal Grid" for r in results)

    def test_check_grid_mapping(self):
        dataset = self.load_dataset(STATIC_FILES["mapping"])
        results = self.cf.check_grid_mapping(dataset)

        assert len(results) == 6
        assert len([r.value for r in results.values() if r.value[0] < r.value[1]]) == 1
        expected_name = (
            u"§5.6 Horizontal Coorindate Reference Systems, Grid Mappings, Projections"
        )
        assert all(r.name == expected_name for r in results.values())

    def test_is_geophysical(self):

        # check whether string type variable, which are not `cf_role`, are
        # properly processed
        dataset = self.load_dataset(STATIC_FILES["string"])
        if dataset.file_format != 'NETCDF4':
            raise RuntimeError(
                "netCDF file of wrong format (not netCDF4) was created for checking"
            )
        try:
            result = cfutil.is_geophysical(dataset, "j")
        except AttributeError:
            pytest.fail(
                "Test probably fails because var.dtype.kind or var.dtype.char "
                "was tested on string-type variable. Consider checking for "
                "`var.dtype is str`"
            )
        assert not result
        # assert False

    # TODO: overhaul to use netCDF global attributes or mocks and variable
    #       attributes
    def test_check_attr_type(self):
        """
        Check that the check_attr_type method checks
        grid_mapping attribute types correctly.
        """

        # test good
        att_name = "test_att"
        att = np.int64(45)
        att_type = "N"  # numeric
        res = self.cf._check_attr_type(att_name, att_type, att)
        self.assertTrue(res[0])
        self.assertEqual(res[1], None)

        # create a temporary variable and test this only
        nc_obj = MockTimeSeries()
        nc_obj.createVariable("temperature", "d", ("time",))
        nc_obj.variables["temperature"].setncattr("test_att", np.float64(45))
        att_name = "test_att"
        _var = nc_obj.variables["temperature"]
        att = np.float64(45)
        att_type = "D"  # numeric, types should match
        res = self.cf._check_attr_type(att_name, att_type, att, _var)
        self.assertTrue(res[0])
        self.assertEqual(res[1], None)

        att_name = "test_att"
        att = "yo"
        att_type = "S"  # string
        res = self.cf._check_attr_type(att_name, att_type, att)
        self.assertTrue(res[0])
        self.assertEqual(res[1], None)

        # test bad
        att_name = "test_att"
        att = np.int64(45)
        att_type = "S"  # string, but att type is numeric
        res = self.cf._check_attr_type(att_name, att_type, att)
        self.assertFalse(res[0])
        self.assertEqual(res[1], "test_att must be a string")

        # test bad
        att_name = "test_att"
        att = "bad"
        att_type = "N"  # numeric, but att type is string
        res = self.cf._check_attr_type(att_name, att_type, att)
        self.assertFalse(res[0])
        self.assertEqual(res[1], "test_att must be a numeric type")

        # create a temporary variable and test this only
        nc_obj = MockTimeSeries()
        nc_obj.createVariable("temperature", "d", ("time",))
        nc_obj.variables["temperature"].setncattr("test_att", np.int32(45))
        _var = nc_obj.variables["temperature"]
        att_name = "test_att"
        att = np.int32(2)
        att_type = "D"  # should be same datatypes
        res = self.cf._check_attr_type(att_name, att_type, att, _var)
        self.assertFalse(res[0])
        self.assertEqual(
            res[1], "test_att must be numeric and must be equivalent to float64 dtype"
        )

    def test_check_grid_mapping_attr_condition(self):
        """
        Ensure the check_grid_mapping_attr_condition() method works as expected.
        """

        # test passes
        attr_name = "latitude_of_projection_origin"
        val = 0
        res = self.cf._check_grid_mapping_attr_condition(val, attr_name)
        self.assertTrue(res[0])

        attr_name = "longitude_of_projection_origin"
        val = 0
        res = self.cf._check_grid_mapping_attr_condition(val, attr_name)
        self.assertTrue(res[0])

        attr_name = "longitude_of_prime_meridian"
        val = 0
        res = self.cf._check_grid_mapping_attr_condition(val, attr_name)
        self.assertTrue(res[0])

        attr_name = "scale_factor_at_central_meridian"
        val = 1
        res = self.cf._check_grid_mapping_attr_condition(val, attr_name)
        self.assertTrue(res[0])

        attr_name = "scale_factor_at_projection_origin"
        val = 1
        res = self.cf._check_grid_mapping_attr_condition(val, attr_name)
        self.assertTrue(res[0])

        attr_name = "standard_parallel"
        val = 0
        res = self.cf._check_grid_mapping_attr_condition(val, attr_name)
        self.assertTrue(res[0])

        attr_name = "straight_vertical_longitude_from_pole"
        val = 0
        res = self.cf._check_grid_mapping_attr_condition(val, attr_name)
        self.assertTrue(res[0])

    def test_check_geographic_region(self):
        dataset = self.load_dataset(STATIC_FILES["bad_region"])
        results = self.cf.check_geographic_region(dataset)
        scored, out_of, messages = get_results(results)

        # only one variable failed this check in this ds out of 2
        assert len(results) == 2
        assert scored < out_of
        assert (
            u"6.1.1 'Neverland' specified by 'neverland' is not a valid region"
            in messages
        )

    def test_check_packed_data(self):
        dataset = self.load_dataset(STATIC_FILES["bad_data_type"])
        results = self.cf.check_packed_data(dataset)
        score, out_of, messages = get_results(results)

        msgs = [
            u"Attributes add_offset and scale_factor have different data type.",
            u"Type of tempvalid_min attribute (int32) does not match variable type (int64)",
            u"Type of temp:valid_max attribute (int32) does not match variable type (int64)",
            u"Type of salinityvalid_min attribute (int32) does not match variable type (float64)",
            u"Type of salinity:valid_max attribute (int32) does not match variable type (float64)",
        ]

        self.assertEqual(len(results), 4)
        self.assertTrue(score < out_of)
        self.assertTrue(all(m in messages for m in msgs))

    def test_compress_packed(self):
        """Tests compressed indexed coordinates"""
        dataset = self.load_dataset(STATIC_FILES["reduced_horizontal_grid"])
        results = self.cf.check_compression_gathering(dataset)
        self.assertTrue(results[0].value)

        dataset = self.load_dataset(STATIC_FILES["bad_data_type"])
        results = self.cf.check_compression_gathering(dataset)
        self.assertFalse(results[0].value)
        self.assertFalse(results[1].value)

    def test_check_all_features_are_same_type(self):
        dataset = self.load_dataset(STATIC_FILES["rutgers"])
        result = self.cf.check_all_features_are_same_type(dataset)
        assert result

        dataset = self.load_dataset(STATIC_FILES["featureType"])
        result = self.cf.check_all_features_are_same_type(dataset)
        assert result

    def test_featureType_is_case_insensitive(self):
        """
        Tests that the featureType attribute is case insensitive
        """
        nc = self.new_nc_file()
        nc.featureType = "timeseriesprofile"
        result = self.cf.check_feature_type(nc)
        self.assertTrue(result.value == (1, 1))

        nc.featureType = "timeSeriesProfile"
        result = self.cf.check_feature_type(nc)
        self.assertTrue(result.value == (1, 1))

        nc.featureType = "traJectorYpRofile"
        result = self.cf.check_feature_type(nc)
        self.assertTrue(result.value == (1, 1))

        # This one should fail
        nc.featureType = "timeseriesprofilebad"
        result = self.cf.check_feature_type(nc)
        self.assertTrue(result.value == (0, 1))

    def test_check_units(self):
        """
        Ensure that container variables are not checked for units but geophysical variables are
        """
        dataset = self.load_dataset(STATIC_FILES["units_check"])
        results = self.cf.check_units(dataset)

        # We don't keep track of the variables names for checks that passed, so
        # we can make a strict assertion about how many checks were performed
        # and if there were errors, which there shouldn't be.
        # FIXME (badams): find a better way of grouping together results by
        #                 variable checked instead of checking the number of
        #                 points scored, which should be deprecated, and
        #                 furthermore is fragile and breaks tests when check
        #                 definitions change
        scored, out_of, messages = get_results(results)
        assert scored == 24
        assert out_of == 24
        assert messages == []

    def test_check_duplicates(self):
        """
        Test to verify that the check identifies duplicate axes. Load the
        duplicate_axis.nc dataset and verify the duplicate axes are accounted
        for.
        """

        dataset = self.load_dataset(STATIC_FILES["duplicate_axis"])
        results = self.cf.check_duplicate_axis(dataset)
        scored, out_of, messages = get_results(results)

        # only one check run here, so we can directly compare all the values
        assert scored != out_of
        assert messages[0] == u"'temp' has duplicate axis X defined by [lon_rho, lon_u]"

    def test_check_multi_dimensional_coords(self):
        """
        Test to verify that multi dimensional coordinates are checked for
        sharing names with dimensions
        """
        dataset = self.load_dataset(STATIC_FILES["multi-dim-coordinates"])
        results = self.cf.check_multi_dimensional_coords(dataset)
        scored, out_of, messages = get_results(results)

        # 4 variables were checked in this ds, 2 of which passed
        assert len(results) == 4
        assert len([r for r in results if r.value[0] < r.value[1]]) == 2
        assert all(r.name == u"§5 Coordinate Systems" for r in results)

    def test_64bit(self):
        dataset = self.load_dataset(STATIC_FILES["ints64"])
        suite = CheckSuite()
        suite.checkers = {"cf": CF1_6Check}
        suite.run(dataset, "cf")

    def test_variable_feature_check(self):

        # non-compliant dataset -- 1/1 fail
        dataset = self.load_dataset(STATIC_FILES["bad-trajectory"])
        results = self.cf.check_variable_features(dataset)
        scored, out_of, messages = get_results(results)
        assert len(results) == 1
        assert scored < out_of
        assert len([r for r in results if r.value[0] < r.value[1]]) == 1
        assert all(r.name == u"§9.1 Features and feature types" for r in results)

        # compliant dataset
        dataset = self.load_dataset(STATIC_FILES["trajectory-complete"])
        results = self.cf.check_variable_features(dataset)
        scored, out_of, messages = get_results(results)
        assert scored == out_of

        # compliant(?) dataset
        dataset = self.load_dataset(STATIC_FILES["trajectory-implied"])
        results = self.cf.check_variable_features(dataset)
        scored, out_of, messages = get_results(results)
        assert scored == out_of

    def test_check_cell_methods(self):
        """Load a dataset (climatology.nc) and check the cell methods.
        This dataset has variable "temperature" which has valid cell_methods
        format, cell_methods attribute, and valid names within the
        cell_methods attribute."""

        dataset = self.load_dataset(STATIC_FILES["climatology"])
        results = self.cf.check_cell_methods(dataset)
        scored, out_of, messages = get_results(results)

        # use itertools.chain() to unpack the lists of messages
        results_list = list(chain(*(r.msgs for r in results if r.msgs)))

        # check the results only have expected headers
        assert set([r.name for r in results]).issubset(
            set([u"§7.1 Cell Boundaries", u"§7.3 Cell Methods"])
        )

        # check that all the expected variables have been hit
        assert all("temperature" in msg for msg in results_list)

        # check that all the results have come back passing
        assert all(r.value[0] == r.value[1] for r in results)

        # create a temporary variable and test this only
        nc_obj = MockTimeSeries()
        nc_obj.createVariable("temperature", "d", ("time",))

        temp = nc_obj.variables["temperature"]
        temp.cell_methods = "lat: lon: mean depth: mean (interval: 20 meters)"
        results = self.cf.check_cell_methods(nc_obj)
        # invalid components lat, lon, and depth -- expect score == (6, 9)
        scored, out_of, messages = get_results(results)
        assert scored != out_of

        temp.cell_methods = "lat: lon: mean depth: mean (interval: x whizbangs)"
        results = self.cf.check_cell_methods(nc_obj)
        scored, out_of, messages = get_results(results)

        # check non-standard comments are gauged correctly
        temp.cell_methods = (
            "lat: lon: mean depth: mean (comment: should not go here interval: 2.5 m)"
        )
        results = self.cf.check_cell_methods(nc_obj)
        scored, out_of, messages = get_results(results)

        self.assertTrue(
            u'§7.3.3 The non-standard "comment:" element must come after any standard elements in cell_methods for variable temperature'
            in messages
        )

        # standalone comments require no keyword
        temp.cell_methods = "lon: mean (This is a standalone comment)"
        results = self.cf.check_cell_methods(nc_obj)
        scored, out_of, messages = get_results(results)
        assert "standalone" not in messages

        # check that invalid keywords dealt with
        temp.cell_methods = (
            "lat: lon: mean depth: mean (invalid_keyword: this is invalid)"
        )
        results = self.cf.check_cell_methods(nc_obj)
        scored, out_of, messages = get_results(results)
        self.assertTrue(
            u'§7.3.3 Invalid cell_methods keyword "invalid_keyword:" for variable temperature. Must be one of [interval, comment]'
            in messages
        )

        # check that "parenthetical elements" are well-formed (they should not be)
        temp.cell_methods = (
            "lat: lon: mean depth: mean (interval 0.2 m interval: 0.01 degrees)"
        )
        results = self.cf.check_cell_methods(nc_obj)
        scored, out_of, messages = get_results(results)
        assert (
            u"§7.3.3 Parenthetical content inside temperature:cell_methods is not well formed: interval 0.2 m interval: 0.01 degrees"
            in messages
        )

    # --------------------------------------------------------------------------------
    # Utility Method Tests
    # --------------------------------------------------------------------------------

    def test_temporal_unit_conversion(self):
        self.assertTrue(units_convertible("hours", "seconds"))
        self.assertFalse(units_convertible("hours", "hours since 2000-01-01"))

    def test_units_temporal(self):
        self.assertTrue(units_temporal("hours since 2000-01-01"))
        self.assertFalse(units_temporal("hours"))
        self.assertFalse(units_temporal("days since the big bang"))


class TestCF1_7(BaseTestCase):
    """Extends the CF 1.6 tests. Most of the tests remain the same."""

    def setUp(self):
        """Initialize a CF1_7Check object."""

        self.cf = CF1_7Check()

    def test_check_actual_range(self):
        """Test the check_actual_range method works as expected"""

        # using a with block closes the ds; for checks operating on the data, we need
        # to initialize and then manually close

        dataset = MockTimeSeries()
        dataset.createVariable("a", "d", ("time",))  # dtype=double, dims=time
        # test that if the variable doesn't have an actual_range attr, no score
        result = self.cf.check_actual_range(dataset)
        assert result == []
        dataset.close()

        # NOTE this is a data check
        # if variable values are equal, actual_range should not exist
        dataset = MockTimeSeries()
        dataset.createVariable("a", "d", ("time",))  # dtype=double, dims=time
        dataset.variables["a"][0:500] = 0  # set all 500 vals to 0
        dataset.variables["a"].setncattr("actual_range", [1])
        result = self.cf.check_actual_range(dataset)
        score, out_of, messages = get_results(result)
        assert score < out_of
        assert len(messages) == 1
        assert messages[0] == u"actual_range of 'a' must be 2 elements"
        dataset.close()

        dataset = MockTimeSeries()
        dataset.createVariable("a", "d", ("time",))  # dtype=double, dims=time
        dataset.variables["a"][0] = 0  # set some arbitrary val so not all equal
        dataset.variables["a"].setncattr("actual_range", [1])
        result = self.cf.check_actual_range(dataset)
        score, out_of, messages = get_results(result)
        assert score < out_of
        assert len(messages) == 1
        assert messages[0] == "actual_range of 'a' must be 2 elements"
        dataset.close()

        # NOTE this is a data check
        # check equality to min and max values
        dataset = MockTimeSeries()
        dataset.createVariable("a", "d", ("time",))
        dataset.variables["a"][0] = -299  # set some arbitrary minimum
        dataset.variables["a"][1] = 10e36  # set some arbitrary max > _FillValue default
        dataset.variables["a"].setncattr("actual_range", [0, 0])  # should fail
        result = self.cf.check_actual_range(dataset)
        score, out_of, messages = get_results(result)
        assert score < out_of
        assert len(messages) == 1
        assert (
            messages[0]
            == "actual_range elements of 'a' inconsistent with its min/max values"
        )
        dataset.close()

        # check equality to valid_range attr
        dataset = MockTimeSeries()
        dataset.createVariable("a", "d", ("time",))
        dataset.variables["a"][0] = -299  # set some arbitrary val to not all equal
        dataset.variables["a"][1] = 10e36  # set some arbitrary max > _FillValue default
        dataset.variables["a"].setncattr("valid_range", [1, 3])  # should conflict
        dataset.variables["a"].setncattr("actual_range", [-299, 10e36])
        result = self.cf.check_actual_range(dataset)
        score, out_of, messages = get_results(result)
        assert score < out_of
        assert len(messages) == 1
        assert messages[0] == '"a"\'s actual_range must be within valid_range'
        dataset.close()

        # check equality to valid_min and valid_max values
        dataset = MockTimeSeries()
        dataset.createVariable("a", "d", ("time",))
        dataset.variables["a"][0] = -299  # set some arbitrary minimum
        dataset.variables["a"][1] = 10e36  # set some arbitrary max > _FillValue default
        dataset.variables["a"].setncattr("valid_min", 42)  # conflicting valid_min/max
        dataset.variables["a"].setncattr("valid_max", 45)
        dataset.variables["a"].setncattr("actual_range", [-299, 10e36])
        result = self.cf.check_actual_range(dataset)
        score, out_of, messages = get_results(result)
        assert score < out_of
        assert len(messages) == 2
        assert (
            messages[0] == '"a"\'s actual_range first element must be >= valid_min (42)'
        )
        assert (
            messages[1]
            == '"a"\'s actual_range second element must be <= valid_max (45)'
        )
        dataset.close()

    def test_check_cell_boundaries(self):
        """Check our over-ridden check_cell_boundaries emthod behaves as expected"""

        dataset = self.load_dataset(STATIC_FILES["grid-boundaries"])
        results = self.cf.check_cell_boundaries(dataset)
        score, out_of, messages = get_results(results)
        assert (score, out_of) == (2, 2)

        dataset = self.load_dataset(STATIC_FILES["cf_example_cell_measures"])
        results = self.cf.check_cell_boundaries(dataset)

        dataset = self.load_dataset(STATIC_FILES["bad_data_type"])
        results = self.cf.check_cell_boundaries(dataset)

        dataset = self.load_dataset(STATIC_FILES["bounds_bad_order"])
        results = self.cf.check_cell_boundaries(dataset)
        score, out_of, messages = get_results(results)
        # Make sure that the rgrid coordinate variable isn't checked for standard_name
        assert (score, out_of) == (0, 2)

        dataset = self.load_dataset(STATIC_FILES["bounds_bad_num_coords"])
        results = self.cf.check_cell_boundaries(dataset)
        score, out_of, messages = get_results(results)
        assert (score, out_of) == (0, 2)

        dataset = self.load_dataset(STATIC_FILES["1d_bound_bad"])
        results = self.cf.check_cell_boundaries(dataset)
        score, out_of, messages = get_results(results)
        assert (score, out_of) == (0, 2)

        # if the variable has formula_terms, the bounds var must also
        with MockTimeSeries() as dataset:
            dataset.createVariable("a", "d", ("time",))
            dataset.createVariable("b", "d", ("time",))
            dataset.variables["a"].setncattr("bounds", "b")  # set bounds variable
            dataset.variables["a"].setncattr("formula_terms", "test")
            results = self.cf.check_cell_boundaries(dataset)
            score, out_of, messages = get_results(results)
            assert score < out_of
            assert (
                "'a' has 'formula_terms' attr, bounds variable 'b' must also have 'formula_terms'"
                in messages
            )

    def test_cell_measures(self):
        """Over-ride the test_cell_measures from CF1_6"""

        # create a temporary variable and test this only
        with MockTimeSeries() as dataset:
            dataset.createVariable("PS", "d", ("time",))  # dtype=double, dims=time
            dataset.variables["PS"].setncattr("cell_measures", "area: cell_area")
            # ensure the cell_measures var is in the dataset
            dataset.createVariable("cell_area", "d", ("time",))
            dataset.variables["cell_area"].setncattr("units", "m2")

            # run the check
            results = self.cf.check_cell_measures(dataset)
            score, out_of, messages = get_results(results)
            assert (score == out_of) and (score > 0)

        # same thing, but test that the cell_area variable is in
        # the global attr "external_variables"

        with MockTimeSeries() as dataset:
            dataset.createVariable("PS", "d", ("time",))  # dtype=double, dims=time
            dataset.variables["PS"].setncattr("cell_measures", "area: cell_area")
            dataset.setncattr("external_variables", ["cell_area"])

            # run the check
            results = self.cf.check_cell_measures(dataset)
            score, out_of, messages = get_results(results)
            assert score > 0
            assert score == out_of

        # now test a dataset with a poorly formatted cell_measure attr
        dataset = self.load_dataset(STATIC_FILES["bad_cell_measure1"])
        results = self.cf.check_cell_measures(dataset)
        score, out_of, messages = get_results(results)
        message = (
            "The cell_measures attribute for variable PS is formatted incorrectly.  "
            "It should take the form of either 'area: cell_var' or 'volume: cell_var' "
            "where cell_var is the variable describing the cell measures"
        )
        assert message in messages

        # test a dataset where the cell_measure attr is not in the dataset or external_variables
        # check for the variable should fail
        dataset = self.load_dataset(STATIC_FILES["bad_cell_measure2"])
        results = self.cf.check_cell_measures(dataset)
        score, out_of, messages = get_results(results)
        message = u"Cell measure variable box_area referred to by PS is not present in dataset variables"
        assert message in messages

    def test_process_vdatum(self):
        # first, we set up a mock SQLite database
        conn_str = ":memory:"
        conn = sqlite3.connect(conn_str)
        cur = conn.cursor()
        # create alias and vertical datum tables without
        # triggers
        cur.execute(
            """
        CREATE TABLE alias_name(
            table_name TEXT NOT NULL CHECK (table_name IN (
                    'unit_of_measure', 'celestial_body', 'ellipsoid',
                    'area', 'prime_meridian', 'geodetic_datum', 'vertical_datum', 'geodetic_crs',
                    'projected_crs', 'vertical_crs', 'compound_crs', 'conversion', 'grid_transformation',
                    'helmert_transformation', 'other_transformation', 'concatenated_operation')),
            auth_name TEXT NOT NULL CHECK (length(auth_name) >= 1),
            code TEXT NOT NULL CHECK (length(code) >= 1),
            alt_name TEXT NOT NULL CHECK (length(alt_name) >= 2),
            source TEXT
        );
        """
        )
        cur.execute(
            """
        CREATE TABLE vertical_datum (
            auth_name TEXT NOT NULL CHECK (length(auth_name) >= 1),
            code TEXT NOT NULL CHECK (length(code) >= 1),
            name TEXT NOT NULL CHECK (length(name) >= 2),
            description TEXT,
            scope TEXT,
            area_of_use_auth_name TEXT NOT NULL,
            area_of_use_code TEXT NOT NULL,
            deprecated BOOLEAN NOT NULL CHECK (deprecated IN (0, 1)),
            CONSTRAINT pk_vertical_datum PRIMARY KEY (auth_name, code)
        );
        """
        )
        cur.execute(
            """INSERT INTO alias_name VALUES
                       ('vertical_datum', 'EPSG', '5103', 'NAVD88', 'EPSG');
                    """
        )

        cur.execute(
            """INSERT INTO vertical_datum VALUES
                    ('EPSG', '5101', 'Ordnance Datum Newlyn', NULL, NULL,
                     'EPSG', '2792', '0')"""
        )

        cur.close()

        self.assertTrue(self.cf._process_v_datum_str("NAVD88", conn))
        self.assertTrue(self.cf._process_v_datum_str("Ordnance Datum Newlyn", conn))
        # NAD83 isn't a vertical datum to begin with, expect failure
        self.assertFalse(self.cf._process_v_datum_str("NAD83", conn))

    def test_check_grid_mapping_crs_wkt(self):
        dataset = self.load_dataset(STATIC_FILES["mapping"])
        valid_crs_check = copy.deepcopy(self.cf)
        dataset.variables["wgs84"] = MockVariable(dataset.variables["wgs84"])
        dataset.variables["wgs84"].crs_wkt = 1
        results = self.cf.check_grid_mapping(dataset)
        score, out_of, messages = get_results(results)
        self.assertIn("crs_wkt attribute must be a string", messages)
        # test with an invalid OGC CRS WKT string
        dataset.variables["wgs84"].crs_wkt = "EPSG:3785"
        results = self.cf.check_grid_mapping(dataset)
        # reuses and appends to old messages, but this is OK since we only need
        # to check that the invalid CRS string message was added
        score, out_of, messages = get_results(results)
        begin_crs_err_msg = "Cannot parse crs_wkt attribute to CRS using Proj4"
        invalid_crs_str = any(s.startswith(begin_crs_err_msg) for s in messages)
        self.assertTrue(invalid_crs_str)

        self.assertIn("crs_wkt attribute must be a string", messages)
        score, out_of, messages = get_results(results)

        valid_crs_wkt = """PROJCS ["OSGB 1936 / British National Grid",
      GEOGCS ["OSGB 1936",
        DATUM ["OSGB 1936", SPHEROID ["Airy 1830", 6377563.396, 299.3249646]],
        PRIMEM ["Greenwich", 0],
        UNIT ["degree", 0.0174532925199433]],
      PROJECTION ["Transverse Mercator"],
      PARAMETER ["False easting", 400000],
      PARAMETER ["False northing", -100000],
      PARAMETER ["Longitude of natural origin", -2.0],
      PARAMETER ["Latitude of natural origin", 49.0],
      PARAMETER ["Scale factor at natural origin", 0.9996012717],
      UNIT ["metre", 1.0]]"""

        dataset.variables["wgs84"].crs_wkt = valid_crs_wkt
        results = valid_crs_check.check_grid_mapping(dataset)
        score, out_of, messages = get_results(results)
        # without false_easting warning in current file
        msg_len = len(
            [
                m
                for m in messages
                if m
                != "false_easting is a required attribute for grid mapping stereographic"
            ]
        )
        self.assertEqual(msg_len, 0)

    def test_check_grid_mapping_coordinates(self):
        """
        Checks that coordinates variables referred to by a grid mapping
        are well-formed and exist.
        """
        dataset = self.load_dataset(STATIC_FILES["grid_mapping_coordinates"])
        valid_grid_mapping = copy.deepcopy(self.cf)
        valid_grid_mapping_2 = copy.deepcopy(self.cf)
        dataset.variables["temp"] = MockVariable(dataset.variables["temp"])
        results = self.cf.check_grid_mapping(dataset)
        self.assertEqual(results["temp"].value[0], results["temp"].value[1])
        malformed_sep = "crsOSGB: x y : lat lon"
        dataset.variables["temp"].grid_mapping = malformed_sep
        results = valid_grid_mapping.check_grid_mapping(dataset)
        self.assertIn(
            "Could not consume entire grid_mapping expression, please check for well-formedness",
            results["temp"].msgs,
        )
        self.assertLess(*results["temp"].value)
        malformed_var = "crsOSGB: x y_null z_null"
        dataset.variables["temp"].grid_mapping = malformed_var
        results = valid_grid_mapping_2.check_grid_mapping(dataset)
        self.assertEqual(
            [
                "Coordinate-related variable y_null referenced by grid_mapping variable crsOSGB must exist in this dataset",
                "Coordinate-related variable z_null referenced by grid_mapping variable crsOSGB must exist in this dataset",
            ],
            results["temp"].msgs,
        )
        self.assertLess(*results["temp"].value)

    def test_check_grid_mapping_vert_datum_geoid_name(self):
        """Checks that geoid_name works proerly"""
        dataset = self.load_dataset(STATIC_FILES["mapping"])
        dataset.variables["wgs84"] = MockVariable(dataset.variables["wgs84"])
        dataset.variables["wgs84"].geoid_name = "NAVD88"
        dataset.variables["wgs84"].geopotential_datum_name = "WGS84"
        geoid_name_good = copy.deepcopy(self.cf)
        geopotential_datum_name_bad = copy.deepcopy(self.cf)
        results = self.cf.check_grid_mapping(dataset)
        score, out_of, messages = get_results(results)
        self.assertIn(
            "Cannot have both 'geoid_name' and 'geopotential_datum_name' attributes in grid mapping variable 'wgs84'",
            messages,
        )
        del dataset.variables["wgs84"].geopotential_datum_name
        results = geoid_name_good.check_grid_mapping(dataset)
        self.assertEqual(*results["wgs84"].value)
        # WGS84 isn't a valid vertical datum name, of course
        dataset.variables["wgs84"].geopotential_datum_name = "WGS84"
        del dataset.variables["wgs84"].geoid_name
        results = geopotential_datum_name_bad.check_grid_mapping(dataset)
        self.assertLess(*results["wgs84"].value)
        self.assertIn(
            "Vertical datum value 'WGS84' for attribute 'geopotential_datum_name' in grid mapping variable 'wgs84' is not valid",
            results["wgs84"].msgs,
        )

    def test_check_conventions_are_cf_1_7(self):
        """Ensure the check_conventions_are_cf_1_7() check works as expected"""

        # create a temporary variable and test this only
        with MockTimeSeries() as dataset:
            # no Conventions attribute
            result = self.cf.check_conventions_version(dataset)
            self.assertFalse(result.value)

        with MockTimeSeries() as dataset:
            # incorrect Conventions attribute
            dataset.setncattr("Conventions", "CF-1.9999")
            result = self.cf.check_conventions_version(dataset)
            self.assertFalse(result.value)

        with MockTimeSeries() as dataset:
            # correct Conventions attribute
            dataset.setncattr("Conventions", "CF-1.7, ACDD-1.3")
            result = self.cf.check_conventions_version(dataset)
            self.assertTrue(result.value)

    def test_appendix_d(self):
        """
        CF 1.7
        Appendix D

        As the CF-1.7 dimensionless vertical coordinates dict extends the 1.6 version,
        this test only examines the extensions made there.
        """

        # For each of the listed dimensionless vertical coordinates,
        # verify that the formula_terms match the provided set of terms
        self.assertTrue(
            no_missing_terms(
                "ocean_s_coordinate_g1",
                {"s", "C", "eta", "depth", "depth_c"},
                dimless_vertical_coordinates_1_7,
            )
        )
        self.assertTrue(
            no_missing_terms(
                "ocean_s_coordinate_g2",
                {"s", "C", "eta", "depth", "depth_c"},
                dimless_vertical_coordinates_1_7,
            )
        )

    def test_check_dimensionless_vertical_coordinate_1_7(self):
        """
        Unit test for _check_dimensionless_vertical_coordinate_1_7 method.
        """
        deprecated_units = ["level", "layer", "sigma_level"]

        ret_val = []

        # create mock dataset for test; create three variables, one as dimensionless
        with MockTimeSeries() as dataset:
            dataset.createVariable("lev", "d")  # dtype=double, dims=1
            dataset.variables["lev"].setncattr(
                "standard_name", "atmosphere_sigma_coordinate"
            )
            dataset.variables["lev"].setncattr(
                "formula_terms", "sigma: lev ps: PS ptop: PTOP"
            )

            dataset.createVariable("PS", "d", ("time",))  # dtype=double, dims=time
            dataset.createVariable("PTOP", "d", ("time",))  # dtype=double, dims=time

            # run the check
            self.cf._check_dimensionless_vertical_coordinate_1_7(
                dataset,
                "lev",
                deprecated_units,
                ret_val,
                dimless_vertical_coordinates_1_7,
            )

            # one should have failed, as no computed_standard_name is assigned
            score, out_of, messages = get_results(ret_val)
            assert score == 0
            assert out_of == 1

            # this time, assign compufted_standard_name
            ret_val = []
            dataset.variables["lev"].setncattr("computed_standard_name", "air_pressure")

            # run the check
            self.cf._check_dimensionless_vertical_coordinate_1_7(
                dataset,
                "lev",
                deprecated_units,
                ret_val,
                dimless_vertical_coordinates_1_7,
            )

            # computed_standard_name is assigned, should pass
            score, out_of, messages = get_results(ret_val)
            assert score == out_of

    def test_dimensionless_vertical(self):
        """
        Section 4.3.2 check, but for CF-1.7 implementation. With the refactor in
        place, these are more of integration tests, but kept here for simplicity.
        """
        # Check affirmative compliance
        dataset = self.load_dataset(STATIC_FILES["dimensionless"])
        dataset.variables["lev"] = MockVariable(dataset.variables["lev"])
        dataset.variables["lev"].computed_standard_name = "air_pressure"
        results = self.cf.check_dimensionless_vertical_coordinates(dataset)
        scored, out_of, messages = get_results(results)

        # all variables checked (2) pass
        assert len(results) == 3
        assert scored == out_of
        assert all(r.name == u"§4.3 Vertical Coordinate" for r in results)

        # make one variable's computed_standard_name incorrect, one should fail
        dataset.variables["lev"].computed_standard_name = "definitely_not_right"
        results = self.cf.check_dimensionless_vertical_coordinates(dataset)
        scored, out_of, messages = get_results(results)

        assert len(results) == 3
        assert scored < out_of
        assert all(r.name == u"§4.3 Vertical Coordinate" for r in results)

    def test_check_attr_type(self):
        """
         Ensure the _check_attr_type method works as expected.
        """

        # create a temporary variable and test this only
        nc_obj = MockTimeSeries()
        nc_obj.createVariable("temperature", "d", ("time",))
        nc_obj.variables["temperature"].setncattr("test_att", np.float64(45))
        att_name = "test_att"
        _var = nc_obj.variables["temperature"]

        # first, test all valid checks show that it's valid
        attr = "my_attr_value"  # string
        attr_type = "S"
        result = self.cf._check_attr_type(att_name, attr_type, attr)
        self.assertTrue(result[0])

        attr = np.int64(1)
        attr_type = "N"
        self.assertTrue(self.cf._check_attr_type(att_name, attr_type, attr)[0])

        attr = np.float64(45)
        attr_type = "D"
        self.assertTrue(self.cf._check_attr_type(att_name, attr_type, attr, _var)[0])

        # check failures
        attr = "my_attr_value"
        attr_type = "N"  # should be numeric
        self.assertFalse(self.cf._check_attr_type(att_name, attr_type, attr)[0])

        attr = np.int(64)
        attr_type = "S"  # should be string
        self.assertFalse(self.cf._check_attr_type(att_name, attr_type, attr)[0])

        nc_obj = MockTimeSeries()
        nc_obj.createVariable("temperature", "d", ("time",))
        nc_obj.variables["temperature"].setncattr("test_att", np.int32(45))
        _var = nc_obj.variables["temperature"]
        attr = np.int32(45)
        attr_type = "D"  # should match
        self.assertFalse(self.cf._check_attr_type(att_name, attr_type, attr, _var)[0])

    def test_check_grid_mapping_attr_condition(self):
        """
        Ensure the CF-1.7 implementation of _check_grid_mapping_attr_condition()
        works as expected.
        """

        # test good

        att_name = "horizontal_datum_name"
        att = "Monte Mario (Rome)"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertTrue(res[0])

        att_name = "prime_meridian_name"
        att = "Athens"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertTrue(res[0])

        att_name = "reference_ellipsoid_name"
        att = "Airy 1830"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertTrue(res[0])

        att_name = "towgs84"
        att = np.array([0, 0, 0], dtype=np.float64)  # len 3
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertTrue(res[0])

        att_name = "towgs84"
        att = np.array([0, 0, 0, 0, 0, 0], dtype=np.float64)  # len 6
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertTrue(res[0])

        att_name = "towgs84"
        att = np.array([0, 0, 0, 0, 0, 0, 0], dtype=np.float64)  # len 7
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertTrue(res[0])

        att_name = "geographic_crs_name"
        att = "NAD83(CSRS98)"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertTrue(res[0])

        att_name = "geoid_name"
        att = "Mayotte 1950"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertTrue(res[0])

        att_name = "geopotential_datum_name"
        att = "NAVD88"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertTrue(res[0])

        att_name = "projected_crs_name"
        att = "Anguilla 1957 / British West Indies Grid"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertTrue(res[0])

        # test bad

        att_name = "horizontal_datum_name"
        att = "bad"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertFalse(res[0])

        att_name = "prime_meridian_name"
        att = "bad"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertFalse(res[0])

        att_name = "reference_ellipsoid_name"
        att = "goofy goober"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertFalse(res[0])

        att_name = "towgs84"
        att = np.array([0, 0, 0], dtype=np.int64)  # len 3, wrong dtype
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertFalse(res[0])

        att_name = "towgs84"
        att = np.array([0, 0, 0, 0], dtype=np.int64)  # len 4
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertFalse(res[0])

        att_name = "towgs84"
        att = np.float64(0)  # single value, right dtype
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertFalse(res[0])

        att_name = "geographic_crs_name"
        att = "badbadbadbadbadnotinhere"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertFalse(res[0])

        att_name = "geoid_name"
        att = "yooooooo"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertFalse(res[0])

        att_name = "geopotential_datum_name"
        att = "NAVBAD BAD"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertFalse(res[0])

        att_name = "projected_crs_name"
        att = "Teddy Bruschi"
        res = self.cf._check_grid_mapping_attr_condition(att, att_name)
        self.assertFalse(res[0])

    def test_check_gmattr_existence_condition_geoid_name_geoptl_datum_name(self):
        # create mock dataset for test; create three variables, one as dimensionless

        # test good (either-or)
        dataset = MockTimeSeries()
        dataset.createVariable("lev", "d")  # dtype=double, dims=1
        dataset.variables["lev"].setncattr("geoid_name", "blah")
        res = self.cf._check_gmattr_existence_condition_geoid_name_geoptl_datum_name(
            dataset.variables["lev"]
        )
        self.assertTrue(res[0])
        dataset.close()

        dataset = MockTimeSeries()
        dataset.createVariable("lev", "d")  # dtype=double, dims=1
        dataset.variables["lev"].setncattr("geopotential_datum_name", "blah")
        res = self.cf._check_gmattr_existence_condition_geoid_name_geoptl_datum_name(
            dataset.variables["lev"]
        )
        self.assertTrue(res[0])
        dataset.close()

        # bad
        dataset = MockTimeSeries()
        dataset.createVariable("lev", "d")  # dtype=double, dims=1
        dataset.variables["lev"].setncattr("geopotential_datum_name", "blah")
        dataset.variables["lev"].setncattr("geoid_name", "blah")
        res = self.cf._check_gmattr_existence_condition_geoid_name_geoptl_datum_name(
            dataset.variables["lev"]
        )
        self.assertFalse(res[0])
        dataset.close()

    def test_check_gmattr_existence_condition_ell_pmerid_hdatum(self):

        # test good (all)
        dataset = MockTimeSeries()
        dataset.createVariable("lev", "d")  # dtype=double, dims=1
        dataset.variables["lev"].setncattr("reference_ellipsoid_name", "blah")
        dataset.variables["lev"].setncattr("prime_meridian_name", "blah")
        dataset.variables["lev"].setncattr("horizontal_datum_name", "blah")
        res = self.cf._check_gmattr_existence_condition_ell_pmerid_hdatum(
            dataset.variables["lev"]
        )
        self.assertTrue(res[0])
        dataset.close()

        # test bad (not all)
        dataset = MockTimeSeries()
        dataset.createVariable("lev", "d")  # dtype=double, dims=1
        dataset.variables["lev"].setncattr("reference_ellipsoid_name", "blah")
        res = self.cf._check_gmattr_existence_condition_ell_pmerid_hdatum(
            dataset.variables["lev"]
        )
        self.assertFalse(res[0])
        dataset.close()

        # test bad (not all)
        dataset = MockTimeSeries()
        dataset.createVariable("lev", "d")  # dtype=double, dims=1
        dataset.variables["lev"].setncattr("reference_ellipsoid_name", "blah")
        dataset.variables["lev"].setncattr("prime_meridian_name", "blah")
        res = self.cf._check_gmattr_existence_condition_ell_pmerid_hdatum(
            dataset.variables["lev"]
        )
        self.assertFalse(res[0])
        dataset.close()
