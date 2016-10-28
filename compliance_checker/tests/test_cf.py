#!/usr/bin/env python
#-*- coding: utf-8 -*-

from compliance_checker.suite import CheckSuite
from compliance_checker.cf import CFBaseCheck, dimless_vertical_coordinates
from compliance_checker.cf.util import is_vertical_coordinate, is_time_variable, units_convertible, units_temporal, StandardNameTable, create_cached_data_dir, download_cf_standard_name_table
from netCDF4 import Dataset
from tempfile import gettempdir
from compliance_checker.tests.resources import STATIC_FILES
from compliance_checker.tests import BaseTestCase

import unittest
import os
import re


class MockVariable(object):
    '''
    For mocking a dataset variable
    '''
    pass


class TestCF(BaseTestCase):

    def setUp(self):
        '''
        Initialize the dataset
        '''
        self.cf = CFBaseCheck()

    # --------------------------------------------------------------------------------
    # Helper Methods
    # --------------------------------------------------------------------------------

    def new_nc_file(self):
        '''
        Make a new temporary netCDF file for the scope of the test
        '''
        nc_file_path = os.path.join(gettempdir(), 'example.nc')
        if os.path.exists(nc_file_path):
            raise IOError('File Exists: %s' % nc_file_path)
        nc = Dataset(nc_file_path, 'w')
        self.addCleanup(os.remove, nc_file_path)
        self.addCleanup(nc.close)
        return nc

    def load_dataset(self, nc_dataset):
        '''
        Return a loaded NC Dataset for the given path
        '''
        if not isinstance(nc_dataset, str):
            raise ValueError("nc_dataset should be a string")

        nc_dataset = Dataset(nc_dataset, 'r')
        self.addCleanup(nc_dataset.close)
        return nc_dataset

    def get_results(self, results):
        '''
        Returns a tuple of the value scored, possible, and a list of messages
        in the result set.
        '''
        out_of = 0
        scored = 0
        for r in results:
            if isinstance(r.value, tuple):
                out_of += r.value[1]
                scored += r.value[0]
            else:
                out_of += 1
                scored += int(r.value)

        # Store the messages
        messages = []
        for r in results:
            messages.extend(r.msgs)

        return scored, out_of, messages

    # --------------------------------------------------------------------------------
    # Compliance Tests
    # --------------------------------------------------------------------------------

    def test_check_data_types(self):
        """
        2.2 The netCDF data types char, byte, short, int, float or real, and double are all acceptable
        """
        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        result = self.cf.check_data_types(dataset)
        assert result.value[0] == result.value[1]

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_data_types(dataset)
        assert result.value == (6, 7)

        dataset = self.load_dataset(STATIC_FILES['chap2'])
        result = self.cf.check_data_types(dataset)
        assert result.value == (6, 7)
        assert 'The variable bad_dtype failed because the datatype is uint16' == result.msgs[0]

    def test_naming_conventions(self):
        '''
        Section 2.3 Naming Conventions

        Variable, dimension and attribute names should begin with a letter and be composed of letters, digits, and underscores.
        '''
        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        results = self.cf.check_naming_conventions(dataset)
        num_var = len(dataset.variables)
        result_dict = {result.name: result for result in results}
        result = result_dict[u'§2.3 Naming Conventions for variables']
        assert result.value == (num_var, num_var)

        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_naming_conventions(dataset)
        result_dict = {result.name: result for result in results}
        result = result_dict[u'§2.3 Naming Conventions for variables']
        assert result.value == (13, 14)
        assert u'variable _poor_dim should begin with a letter and be composed of letters, digits, and underscores' == result.msgs[0]
        score, out_of, messages = self.get_results(results)
        assert (score, out_of) == (49, 51)

        dataset = self.load_dataset(STATIC_FILES['chap2'])
        results = self.cf.check_naming_conventions(dataset)
        result_dict = {result.name: result for result in results}
        result = result_dict[u'§2.3 Naming Conventions for variables']
        assert result.value == (6, 7)
        assert u'variable bad name should begin with a letter and be composed of letters, digits, and underscores' == result.msgs[0]

        result = result_dict[u'§2.3 Naming Conventions for attributes']
        assert result.msgs[0] == ('attribute no_reason:_bad_attr should begin with a letter and be '
                                  'composed of letters, digits, and underscores')
        assert result.msgs[1] == ('global attribute bad global should begin with a letter and be '
                                  'composed of letters, digits, and underscores')

        score, out_of, messages = self.get_results(results)
        assert (score, out_of) == (16, 19)

    def test_check_names_unique(self):
        """
        2.3 names should not be distinguished purely by case, i.e., if case is disregarded, no two names should be the same.
        """
        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        result = self.cf.check_names_unique(dataset)

        num_var = len(dataset.variables)
        expected = (num_var,) * 2

        self.assertEqual(result.value, expected)

        dataset = self.load_dataset(STATIC_FILES['chap2'])
        result = self.cf.check_names_unique(dataset)
        assert result.value == (6, 7)
        assert result.msgs[0] == 'Variables are not case sensitive. Duplicate variables named: not_unique'

    def test_check_dimension_names(self):
        """
        2.4 A variable may have any number of dimensions, including zero, and the dimensions must all have different names.
        """

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_dimension_names(dataset)
        assert result.value == (6, 7)

        dataset = self.load_dataset(STATIC_FILES['chap2'])
        result = self.cf.check_dimension_names(dataset)
        assert result.msgs[0] == u'no_reason has two or more dimensions named time'

    def test_check_dimension_order(self):
        """
        2.4 If any or all of the dimensions of a variable have the interpretations of "date or time" (T), "height or depth" (Z),
        "latitude" (Y), or "longitude" (X) then we recommend, those dimensions to appear in the relative order T, then Z, then Y,
        then X in the CDL definition corresponding to the file. All other dimensions should, whenever possible, be placed to the
        left of the spatiotemporal dimensions.
        """
        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_dimension_order(dataset)
        assert result.value == (5, 7)
        assert result.msgs[0] == ("salinity's dimensions are not in the recommended order "
                                  "T, X, Y, Z. They are Y, U, U")
        assert result.msgs[1] == ("really_bad's dimensions are not in the recommended order "
                                  "T, X, Y, Z. They are Y, U")

    def test_check_fill_value_outside_valid_range(self):
        """
        2.5.1 The _FillValue should be outside the range specified by valid_range (if used) for a variable.
        """

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_fill_value_outside_valid_range(dataset)
        assert result.msgs[0] == ('salinity:_FillValue (1.0) should be outside the '
                                  'range specified by valid_min/valid_max (-10, 10)')

        dataset = self.load_dataset(STATIC_FILES['chap2'])
        result = self.cf.check_fill_value_outside_valid_range(dataset)
        assert result.value == (1, 2)
        assert result.msgs[0] == ('wind_speed:_FillValue (12.0) should be outside the '
                                  'range specified by valid_min/valid_max (0.0, 20.0)')

    def test_check_conventions_are_cf_16(self):
        """
        2.6.1 the NUG defined global attribute Conventions to the string value "CF-1.6"
        """
        # :Conventions = "CF-1.6"
        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        result = self.cf.check_conventions_are_cf_16(dataset)
        self.assertTrue(result.value)

        # :Conventions = "CF-1.6 ,ACDD" ;
        dataset = self.load_dataset(STATIC_FILES['conv_multi'])
        result = self.cf.check_conventions_are_cf_16(dataset)
        self.assertTrue(result.value)

        # :Conventions = "NoConvention"
        dataset = self.load_dataset(STATIC_FILES['conv_bad'])
        result = self.cf.check_conventions_are_cf_16(dataset)
        self.assertFalse(result.value)
        assert result.msgs[0] == 'Conventions global attribute does not contain "CF-1.6"'

    def test_check_convention_globals(self):
        """
        2.6.2 title/history global attributes, must be strings. Do not need to exist.
        """
        # check for pass
        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        result = self.cf.check_convention_globals(dataset)
        assert result.value == (2, 2)
        # check if it doesn't exist that we pass
        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_convention_globals(dataset)
        assert result.value == (0, 2)
        assert result.msgs[0] == 'global attribute title should exist and be a non-empty string'

    def test_check_convention_possibly_var_attrs(self):
        """
        3.1 The units attribute is required for all variables that represent dimensional quantities
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
        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        result = self.cf.check_convention_possibly_var_attrs(dataset)
        # 10x comment attrs
        # 1x institution
        # 1x source
        # 1x EMPTY references
        assert result.value == (15, 16)
        assert result.msgs[0] == "references global attribute should be a non-empty string"

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_convention_possibly_var_attrs(dataset)
        # no references
        # institution is a 10L
        # no source
        # comments doment matter unless they're empty
        assert result.value == (1, 4)
        assert result.msgs[0] == 'salinity:institution should be a non-empty string'
        assert result.msgs[1] == 'source should be defined'
        assert result.msgs[2] == 'references should be defined'

    def test_check_standard_name(self):
        """
        3.3 A standard name is associated with a variable via the attribute standard_name which takes a
        string value comprised of a standard name optionally followed by one or more blanks and a
        standard name modifier
        """
        dataset = self.load_dataset(STATIC_FILES['2dim'])
        result = self.cf.check_standard_name(dataset)
        for each in result:
            self.assertTrue(each.value)

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        results = self.cf.check_standard_name(dataset)
        result_dict = {result.name: result for result in results}
        result = result_dict[u'§3.3 Variable time has valid standard_name attribute']
        assert result.value == (0, 2)
        assert 'standard_name undefined is not defined' in result.msgs[1]

        result = result_dict[u'§3.3 Variable latitude has valid standard_name attribute']
        assert result.value == (0, 2)
        assert 'standard_name undefined is not defined' in result.msgs[1]

        result = result_dict[u'§3.3 Variable salinity has valid standard_name attribute']
        assert result.value == (1, 2)
        assert 'standard_name Chadwick is not defined in Standard Name Table' in result.msgs[0]

        result = result_dict[u'§3.3 standard_name modifier for salinity is valid']
        assert result.value == (0, 1)

        assert len(result_dict) == 8

        dataset = self.load_dataset(STATIC_FILES['reduced_horizontal_grid'])
        results = self.cf.check_standard_name(dataset)
        score, out_of, messages = self.get_results(results)
        # Make sure that the rgrid coordinate variable isn't checked for standard_name
        # time, lat, lon exist with three checks each
        assert (score, out_of) == (6, 6)

    def test_check_ancillary_variables(self):
        '''
        Test to ensure that ancillary variables are properly checked
        '''

        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        results = self.cf.check_ancillary_variables(dataset)
        result_dict = {result.name: result for result in results}
        result = result_dict[u'§3.4 Ancillary Variables defined by temperature']
        assert result.value == (2, 2)

        dataset = self.load_dataset(STATIC_FILES['bad_reference'])
        results = self.cf.check_ancillary_variables(dataset)
        result_dict = {result.name: result for result in results}
        result = result_dict[u'§3.4 Ancillary Variables defined by temp']
        assert result.value == (1, 2)
        assert "temp_qc is not a variable in this dataset" == result.msgs[0]

    def test_download_standard_name_table(self):
        """
        Test that a user can download a specific standard name table
        """
        version = '35'

        data_directory = create_cached_data_dir()
        location = os.path.join(data_directory, 'cf-standard-name-table-test-{0}.xml'.format(version))
        download_cf_standard_name_table(version, location)

        # Test that the file now exists in location and is the right version
        self.assertTrue(os.path.isfile(location))
        std_names = StandardNameTable(location)
        self.assertEqual(std_names._version, version)
        self.addCleanup(os.remove, location)

    def test_check_flags(self):
        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        results = self.cf.check_flags(dataset)
        result_dict = {result.name: result for result in results}
        result = result_dict[u'§3.5 time_qc is a valid flags variable']
        assert result.value == (1, 1)
        result = result_dict[u'§3.5 flag_meanings for time_qc']
        assert result.value == (3, 3)
        result = result_dict[u'§3.5 flag_values for time_qc']
        assert result.value == (4, 4)
        # lat(time);
        #   lat:flag_meanings = "";
        result = result_dict[u'§3.5 lat is a valid flags variable']
        assert result.value == (0, 1)
        result = result_dict[u'§3.5 flag_meanings for lat']
        assert result.value == (2, 3)
        assert "flag_meanings can't be empty" == result.msgs[0]

    def test_check_flag_masks(self):
        dataset = self.load_dataset(STATIC_FILES['ghrsst'])
        results = self.cf.check_flags(dataset)
        scored, out_of, messages = self.get_results(results)
        # This is an example of a perfect dataset for flags
        assert scored > 0
        assert scored == out_of

    def test_check_bad_units(self):

        dataset = self.load_dataset(STATIC_FILES['2dim'])
        results = self.cf.check_units(dataset)
        for result in results:
            self.assert_result_is_good(result)

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        results = self.cf.check_units(dataset)
        result_dict = {result.name: result for result in results}

        result = result_dict[u'§3.1 Variable time contains valid CF units']
        # Since no units are specified, they can't be deprecated
        assert result.value == (1, 3)
        assert 'units attribute is required for time' in result.msgs
        assert 'units attribute for time needs to be a string' in result.msgs

        # it's Degrees_E which is a valid udunits. The preferred units are
        # degrees_east and they are checked in the check_longitude check
        result = result_dict[u'§3.1 Variable longitude\'s units are contained in UDUnits']
        assert result.value == (1, 1)

        result = result_dict[u'§3.1 Variable temp contains valid CF units']
        assert result.value == (3, 3)

        result = result_dict[u'§3.1 Variable temp\'s units are contained in UDUnits']
        assert result.value == (1, 1)

        dataset = self.load_dataset(STATIC_FILES['bad_units'])
        results = self.cf.check_units(dataset)
        result_dict = {result.name: result for result in results}

        # time(time)
        #   time:units = "s"
        result = result_dict[u'§3.1 Variable time contains valid CF units']
        # They are valid and even valid UDUnits
        assert result.value == (3, 3)
        result = result_dict[u"§3.1 Variable time's units are contained in UDUnits"]
        assert result.value == (1, 1)

        # But they are not appropriate for time
        result = result_dict[u"§3.1 Variable time's units are appropriate for the standard_name time"]
        assert result.value == (0, 1)

        # lat;
        #   lat:units = "degrees_E";
        # Should all be good
        result = result_dict[u"§3.1 Variable lat's units are appropriate for the standard_name latitude"]
        assert result.value == (0, 1)

        # lev;
        #   lev:units = "level";
        # level is deprecated
        result = result_dict[u"§3.1 Variable lev contains valid CF units"]
        assert result.value == (2, 3)
        assert 'units for lev, "level" are deprecated by CF 1.6' in result.msgs

        # temp_count(time);
        #   temp_count:standard_name = "atmospheric_temperature number_of_observations";
        #   temp_count:units = "1";
        result = result_dict[u"§3.1 Variable temp_count's units are appropriate for "
                             u"the standard_name atmospheric_temperature number_of_observations"]
        assert result.value == (1, 1)

    def test_coordinate_types(self):
        '''
        Section 4 Coordinate Types

        We strongly recommend that coordinate variables be used for all coordinate types whenever they are applicable.
        '''
        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_coordinate_vars_for_all_coordinate_types(dataset)
        for each in result:
            self.assertTrue(each.value)

    def test_check_coordinate_types(self):

        dataset = self.load_dataset(STATIC_FILES['coordinate_types'])
        results = self.cf.check_coordinate_types(dataset)
        result_dict = {result.name: result for result in results}
        result = result_dict[u'§4 time is a valid coordinate type']
        assert result.value == (1, 1)

        result = result_dict[u'§4 time has suggested standard_name for coordinate type']
        assert result.value == (2, 2)

        result = result_dict[u'§4 lat has suggested mapping from axis to standard_name']
        assert result.value == (2, 3)
        assert result.msgs[0] == 'standard_name for axis X is suggested to be longitude. Is currently latitude'

        result = result_dict[u'§4 temperature has suggested standard_name for coordinate type']
        assert result.value == (1, 2)

        scored, out_of, messages = self.get_results(results)
        assert (scored, out_of) == (21, 25)

        dataset = self.load_dataset(STATIC_FILES['reduced_horizontal_grid'])
        results = self.cf.check_coordinate_types(dataset)
        scored, out_of, messages = self.get_results(results)
        # CF recommends coordinate types defining lat/lon/z be coordinate variables
        assert (scored, out_of) == (10, 12)

    def test_latitude(self):
        '''
        Section 4.1 Latitude Coordinate
        '''
        # Check compliance
        dataset = self.load_dataset(STATIC_FILES['example-grid'])
        results = self.cf.check_latitude(dataset)
        # Four checks per latitude variable
        assert len(results) == 4
        score, out_of, messages = self.get_results(results)
        assert (score, out_of) == (4, 4)

        # Verify non-compliance
        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_latitude(dataset)

        result_dict = {result.name: result for result in results}
        scored, out_of, messages = self.get_results(results)

        result = result_dict[u'§4.1 Latitude variable lat has required units attribute']
        assert result.value == (0, 1)
        assert result.msgs[0] == "latitude variable 'lat' must define units"

        result = result_dict[u'§4.1 Latitude variable lat uses recommended units']
        assert result.value == (0, 1)

        result = result_dict[u'§4.1 Latitude variable lat defines units using degrees_north']
        assert result.value == (0, 1)

        result = result_dict[u'§4.1 Latitude variable lat defines either standard_name or axis']
        assert result.value == (1, 1)

        result = result_dict[u'§4.1 Latitude variable lat_uv has required units attribute']
        assert result.value == (1, 1)

        result = result_dict[u'§4.1 Latitude variable lat_uv uses recommended units']
        assert result.value == (1, 1)

        result = result_dict[u'§4.1 Latitude variable lat_uv defines units using degrees_north']
        assert result.value == (0, 1)

        result = result_dict[u'§4.1 Latitude variable lat_uv defines either standard_name or axis']
        assert result.value == (1, 1)

        assert (scored, out_of) == (6, 12)

        dataset = self.load_dataset(STATIC_FILES['rotated_pole_grid'])
        results = self.cf.check_latitude(dataset)
        scored, out_of, messages = self.get_results(results)
        assert (scored, out_of) == (7, 7)

    def test_longitude(self):
        '''
        Section 4.2 Longitude Coordinate
        '''
        # Check compliance
        dataset = self.load_dataset(STATIC_FILES['example-grid'])
        results = self.cf.check_longitude(dataset)
        score, out_of, messages = self.get_results(results)
        assert (score, out_of) == (4, 4)

        # Verify non-compliance
        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_longitude(dataset)

        result_dict = {result.name: result for result in results}
        scored, out_of, messages = self.get_results(results)

        result = result_dict[u'§4.1 Longitude variable lon has required units attribute']
        assert result.value == (0, 1)
        assert result.msgs[0] == "longitude variable 'lon' must define units"

        result = result_dict[u'§4.1 Longitude variable lon uses recommended units']
        assert result.value == (0, 1)

        result = result_dict[u'§4.1 Longitude variable lon defines units using degrees_east']
        assert result.value == (0, 1)

        result = result_dict[u'§4.1 Longitude variable lon defines either standard_name or axis']
        assert result.value == (1, 1)

        result = result_dict[u'§4.1 Longitude variable lon_uv has required units attribute']
        assert result.value == (1, 1)

        result = result_dict[u'§4.1 Longitude variable lon_uv uses recommended units']
        assert result.value == (1, 1)

        result = result_dict[u'§4.1 Longitude variable lon_uv defines units using degrees_east']
        assert result.value == (0, 1)

        result = result_dict[u'§4.1 Longitude variable lon_uv defines either standard_name or axis']
        assert result.value == (1, 1)

        assert (scored, out_of) == (6, 12)

        dataset = self.load_dataset(STATIC_FILES['rotated_pole_grid'])
        results = self.cf.check_latitude(dataset)
        scored, out_of, messages = self.get_results(results)
        assert (scored, out_of) == (7, 7)

    def test_is_vertical_coordinate(self):
        '''
        Section 4.3 Qualifiers for Vertical Coordinate

        NOTE: The standard doesn't explicitly say that vertical coordinates must be a
        coordinate type.
        '''
        # Make something that I can attach attrs to
        mock_variable = MockVariable

        # Proper name/standard_name
        known_name = mock_variable()
        known_name.standard_name = 'depth'
        self.assertTrue(is_vertical_coordinate('not_known', known_name))

        # Proper Axis
        axis_set = mock_variable()
        axis_set.axis = 'Z'
        self.assertTrue(is_vertical_coordinate('not_known', axis_set))

        # Proper units
        units_set = mock_variable()
        units_set.units = 'dbar'
        self.assertTrue(is_vertical_coordinate('not_known', units_set))

        # Proper units/positive
        positive = mock_variable()
        positive.units = 'm'
        positive.positive = 'up'
        self.assertTrue(is_vertical_coordinate('not_known', positive))

    def test_vertical_dimension(self):
        '''
        Section 4.3.1 Dimensional Vertical Coordinate
        '''
        # Check for compliance
        dataset = self.load_dataset(STATIC_FILES['example-grid'])
        results = self.cf.check_dimensional_vertical_coordinate(dataset)
        assert len(results) == 1
        assert results[0].name == u'§4.3.1 height is a valid vertical coordinate'
        assert results[0].value == (2, 2)

        dataset = self.load_dataset(STATIC_FILES['illegal-vertical'])
        results = self.cf.check_dimensional_vertical_coordinate(dataset)
        assert len(results) == 1
        assert results[0].name == u'§4.3.1 z is a valid vertical coordinate'
        assert results[0].value == (0, 2)
        assert results[0].msgs[0] == 'units must be defined for vertical coordinates, there is no default'
        assert results[0].msgs[1] == ("vertical coordinates not defining pressure must include a positive attribute that "
                                      "is either 'up' or 'down'")

    def test_appendix_d(self):
        '''
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
        '''

        dimless = dict(dimless_vertical_coordinates)

        def verify(std_name, test_str):
            regex_matches = re.match(dimless[std_name], test_str)
            self.assertIsNotNone(regex_matches)

        # For each of the listed dimensionless vertical coordinates,
        # verify that the formula_terms match the provided regex
        verify('atmosphere_ln_pressure_coordinate',
               "p0: var1 lev: var2")
        verify('atmosphere_sigma_coordinate',
               "sigma: var1 ps: var2 ptop: var3")
        verify('atmosphere_hybrid_sigma_pressure_coordinate',
               "a: var1 b: var2 ps: var3 p0: var4")
        verify('atmosphere_hybrid_height_coordinate',
               "a: var1 b: var2 orog: var3")
        verify('atmosphere_sleve_coordinate',
               "a: var1 b1: var2 b2: var3 ztop: var4 zsurf1: var5 zsurf2: var6")
        verify('ocean_sigma_coordinate',
               "sigma: var1 eta: var2 depth: var3")
        verify('ocean_s_coordinate',
               "s: var1 eta: var2 depth: var3 a: var4 b: var5 depth_c: var6")
        verify('ocean_sigma_z_coordinate',
               "sigma: var1 eta: var2 depth: var3 depth_c: var4 nsigma: var5 zlev: var6")
        verify('ocean_double_sigma_coordinate',
               "sigma: var1 depth: var2 z1: var3 z2: var4 a: var5 href: var6 k_c: var7")

    def test_dimensionless_vertical(self):
        '''
        Section 4.3.2
        '''
        # Check affirmative compliance
        dataset = self.load_dataset(STATIC_FILES['dimensionless'])
        results = self.cf.check_dimensionless_vertical_coordinate(dataset)

        result_dict = {result.name: result for result in results}
        result = result_dict[u'§4.3.2 lev does not contain deprecated units']
        assert result.value == (1, 1)
        result = result_dict[u'§4.3.2 lev has valid formula_terms']
        assert result.value == (6, 6)

        # Check negative compliance
        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_dimensionless_vertical_coordinate(dataset)

        result_dict = {result.name: result for result in results}
        result = result_dict[u'§4.3.2 lev1 does not contain deprecated units']
        assert result.value == (1, 1)
        result = result_dict[u'§4.3.2 lev1 has valid formula_terms']
        assert result.value == (0, 1)
        assert result.msgs[0] == u'formula_terms is a required attribute and must be a non-empty string'

        result = result_dict[u'§4.3.2 lev2 has valid formula_terms']
        assert result.value == (3, 6)
        assert result.msgs[0] == 'variable var1 referenced by formula_terms does not exist'
        assert result.msgs[1] == 'variable var2 referenced by formula_terms does not exist'
        assert result.msgs[2] == 'variable var3 referenced by formula_terms does not exist'

    def test_is_time_variable(self):
        var1 = MockVariable()
        var1.standard_name = 'time'
        self.assertTrue(is_time_variable('not_time', var1))

        var2 = MockVariable()
        self.assertTrue(is_time_variable('time', var2))

        self.assertFalse(is_time_variable('not_time', var2))

        var3 = MockVariable()
        var3.axis = 'T'
        self.assertTrue(is_time_variable('maybe_time', var3))

        var4 = MockVariable()
        var4.units = 'seconds since 1900-01-01'
        self.assertTrue(is_time_variable('maybe_time', var4))

    def test_check_time_coordinate(self):
        dataset = self.load_dataset(STATIC_FILES['example-grid'])
        results = self.cf.check_time_coordinate(dataset)
        for r in results:
            self.assertTrue(r.value)

        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_time_coordinate(dataset)

        scored, out_of, messages = self.get_results(results)

        assert 'time does not have correct time units' in messages
        assert (scored, out_of) == (1, 2)

    def test_check_calendar(self):
        dataset = self.load_dataset(STATIC_FILES['example-grid'])
        results = self.cf.check_calendar(dataset)
        for r in results:
            self.assertTrue(r.value)

        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_calendar(dataset)
        scored, out_of, messages = self.get_results(results)

        assert "Variable time should have a valid calendar: 'nope' is not a valid calendar" in messages

    def test_check_aux_coordinates(self):
        dataset = self.load_dataset(STATIC_FILES['illegal-aux-coords'])
        results = self.cf.check_aux_coordinates(dataset)
        result_dict = {result.name: result for result in results}
        result = result_dict[u"§5.0 Auxiliary Coordinates of h_temp must have a subset of h_temp's dimensions"]
        assert result.value == (2, 4)
        regx = (r"dimensions for auxiliary coordinate variable lat \([xy]c, [xy]c\) are not a subset of dimensions for variable "
                r"h_temp \(xc\)")
        assert re.match(regx, result.msgs[0]) is not None

        result = result_dict[u"§5.0 Auxiliary Coordinates of sal must have a subset of sal's dimensions"]
        assert result.value == (4, 4)

    def test_check_grid_coordinates(self):
        dataset = self.load_dataset(STATIC_FILES['2dim'])
        results = self.cf.check_grid_coordinates(dataset)
        scored, out_of, messages = self.get_results(results)

        result_dict = {result.name: result for result in results}
        result = result_dict[u'§5.6 Grid Feature T is associated with true latitude and true longitude']
        assert result.value == (2, 2)
        assert (scored, out_of) == (2, 2)

    def test_check_two_dimensional(self):
        dataset = self.load_dataset(STATIC_FILES['2dim'])
        results = self.cf.check_grid_coordinates(dataset)
        for r in results:
            self.assertTrue(r.value)
        # Need the bad testing
        dataset = self.load_dataset(STATIC_FILES['bad2dim'])
        results = self.cf.check_grid_coordinates(dataset)
        scored, out_of, messages = self.get_results(results)

        result_dict = {result.name: result for result in results}
        # Missing association
        result = result_dict[u'§5.6 Grid Feature T is associated with true latitude and true longitude']
        assert result.msgs[0] == 'T is not associated with a coordinate defining true latitude and sharing a subset of dimensions'
        # Dimensions aren't a subet of the variables'
        result = result_dict[u'§5.6 Grid Feature C is associated with true latitude and true longitude']
        assert result.msgs[0] == 'C is not associated with a coordinate defining true latitude and sharing a subset of dimensions'

    def test_check_reduced_horizontal_grid(self):
        dataset = self.load_dataset(STATIC_FILES['rhgrid'])
        results = self.cf.check_reduced_horizontal_grid(dataset)

        result_dict = {result.name: result for result in results}
        result = result_dict[u'§5.3 PS is a valid reduced horizontal grid']
        assert result.value == (7, 7)

        dataset = self.load_dataset(STATIC_FILES['bad-rhgrid'])
        results = self.cf.check_reduced_horizontal_grid(dataset)

        result_dict = {result.name: result for result in results}
        result = result_dict[u'§5.3 PSa is a valid reduced horizontal grid']
        assert result.value == (6, 7)
        assert result.msgs[0] == "PSa must be associated with a valid longitude coordinate"
        result = result_dict[u'§5.3 PSb is a valid reduced horizontal grid']
        # The dimensions don't line up but another §5.0 check catches it.
        assert result.value == (7, 7)

    def test_check_grid_mapping(self):
        dataset = self.load_dataset(STATIC_FILES['mapping'])
        results = self.cf.check_grid_mapping(dataset)

        result_dict = {result.name: result for result in results}
        result = result_dict[u'§5.6 Grid Mapping Variable epsg must define a valid grid mapping']
        assert result.value == (7, 8)
        assert result.msgs[0] == 'false_easting is a required attribute for grid mapping stereographic'

        result = result_dict[u'§5.6 Grid Mapping Variable wgs84 must define a valid grid mapping']
        assert result.value == (3, 3)

        result = result_dict[u'§5.6 Variable lat defining a grid mapping has valid grid_mapping attribute']
        assert result.value == (2, 2)

    def test_check_scalar_coordinate_system(self):
        dataset = self.load_dataset(STATIC_FILES['scalar_coordinate_variable'])
        results = self.cf.check_scalar_coordinate_system(dataset)
        self.assertEqual(len(results), 2)
        for r in results:
            if r.name[1] == 'HEIGHT':
                self.assertEqual(r.value, (0, 1))
            elif r.name[1] == 'DEPTH':
                self.assertEqual(r.value, (2, 2))
            else:
                self.assertTrue(False, 'Unexpected variable in results of check_scalar_coordinate_system')

    def test_check_geographic_region(self):
        dataset = self.load_dataset(STATIC_FILES['bad_region'])
        results = self.cf.check_geographic_region(dataset)

        self.assertFalse(results[0].value)
        self.assertTrue(results[1].value)

    # def test_check_cell_boundaries(self):
    #    dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
    #    results = self.cf.check_cell_boundaries(dataset)
    #    print results
    #    self.assertTrue(results[0].value)

    def test_check_packed_data(self):
        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        results = self.cf.check_packed_data(dataset)
        self.assertEqual(len(results), 4)
        self.assertFalse(results[0].value)
        self.assertTrue(results[1].value)
        self.assertTrue(results[2].value)
        self.assertFalse(results[3].value)

    def test_check_instance_variable(self):
        dataset = self.load_dataset(STATIC_FILES['bad-instance'])
        results = self.cf.check_instance_variables(dataset)

        result_dict = {result.name: result for result in results}
        result = result_dict[u'§9.5 Instance Variable station_name is not referenced as a coordinate variable']
        assert result.value == (0, 1)

    def test_check_compression(self):
        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        results = self.cf.check_compression(dataset)
        assert results[0].value == (2, 2)
        assert results[1].value == (0, 2)

    def test_check_all_features_are_same_type(self):
        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        results = self.cf.check_all_features_are_same_type(dataset)
        assert results is None

        dataset = self.load_dataset(STATIC_FILES['featureType'])
        results = self.cf.check_all_features_are_same_type(dataset)
        self.assertTrue(results.value)

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        results = self.cf.check_all_features_are_same_type(dataset)
        self.assertFalse(results.value)

    def test_check_orthogonal_multidim_array(self):
        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        results = self.cf.check_orthogonal_multidim_array(dataset)
        for each in results:
            self.assertTrue(each.value)

    def test_check_incomplete_multidim_array(self):
        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        results = self.cf.check_incomplete_multidim_array(dataset)
        for each in results:
            self.assertTrue(each.value)

    def test_check_contiguous_ragged_array(self):
        dataset = self.load_dataset(STATIC_FILES['cont_ragged'])
        results = self.cf.check_contiguous_ragged_array(dataset)
        for each in results:
            self.assertTrue(each.value)

    def test_check_indexed_ragged_array(self):
        dataset = self.load_dataset(STATIC_FILES['index_ragged'])
        results = self.cf.check_indexed_ragged_array(dataset)
        for each in results:
            self.assertTrue(each.value)

    def test_check_feature_type(self):
        dataset = self.load_dataset(STATIC_FILES['index_ragged'])
        results = self.cf.check_feature_type(dataset)
        self.assertTrue(results.value)

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        results = self.cf.check_feature_type(dataset)
        self.assertFalse(results.value)

    def test_check_coordinates_and_metadata(self):
        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        results = self.cf.check_coordinates_and_metadata(dataset)
        self.assertFalse(results[0].value)
        self.assertTrue(results[1].value)
        self.assertFalse(results[2].value)

        dataset = self.load_dataset(STATIC_FILES['index_ragged'])
        results = self.cf.check_coordinates_and_metadata(dataset)
        self.assertTrue(results[-1].value)

        dataset = self.load_dataset(STATIC_FILES['coordinates_and_metadata'])
        results = self.cf.check_coordinates_and_metadata(dataset)
        self.assertTrue(len(results) == 2)
        self.assertFalse(results[0].value)
        self.assertFalse(results[1].value)

    def test_check_missing_data(self):
        dataset = self.load_dataset(STATIC_FILES['index_ragged'])
        results = self.cf.check_missing_data(dataset)
        for each in results:
            self.assertTrue(each.value)

        dataset = self.load_dataset(STATIC_FILES['bad_missing_data'])
        results = self.cf.check_missing_data(dataset)
        for each in results:
            self.assertFalse(each.value)

    def test_check_units(self):
        '''
        Ensure that container variables are not checked for units but geophysical variables are
        '''
        dataset = self.load_dataset(STATIC_FILES['units_check'])
        results = self.cf.check_units(dataset)

        # We don't keep track of the variables names for checks that passed, so
        # we can make a strict assertion about how many checks were performed
        # and if there were errors, which there shouldn't be.
        scored, out_of, messages = self.get_results(results)
        assert scored == 20
        assert out_of == 20
        assert messages == []

    def test_check_duplicates(self):
        '''
        Test to verify that the check identifies duplicate axes
        '''
        dataset = self.load_dataset(STATIC_FILES['duplicate_axis'])
        results = self.cf.check_duplicate_axis(dataset)
        result_dict = {result.name: result for result in results}
        result = result_dict[u'§5.0 Variable temp does not contain duplicate coordinates']
        assert result.msgs[0] == 'duplicate axis X defined by lon_rho'

    def test_check_multi_dimensional_coords(self):
        '''
        Test to verify that multi dimensional coordinates are checked for
        sharing names with dimensions
        '''
        dataset = self.load_dataset(STATIC_FILES['multi-dim-coordinates'])
        results = self.cf.check_multi_dimensional_coords(dataset)
        scored, out_of, messages = self.get_results(results)
        result_dict = {result.name: result for result in results}
        result = result_dict[u'§5.0 multidimensional coordinate xlon should not have the same name as dimension']
        assert result.msgs[0] == 'xlon shares the same name as one of its dimensions'
        result = result_dict[u'§5.0 multidimensional coordinate xlat should not have the same name as dimension']
        assert result.msgs[0] == 'xlat shares the same name as one of its dimensions'

        assert (scored, out_of) == (2, 4)

    def test_64bit(self):
        dataset = self.load_dataset(STATIC_FILES['ints64'])
        suite = CheckSuite()
        suite.checkers = {
            'cf'        : CFBaseCheck
        }
        suite.run(dataset, 'cf')

    # --------------------------------------------------------------------------------
    # Utility Method Tests
    # --------------------------------------------------------------------------------

    def test_temporal_unit_conversion(self):
        self.assertTrue(units_convertible('hours', 'seconds'))
        self.assertFalse(units_convertible('hours', 'hours since 2000-01-01'))

    def test_units_temporal(self):
        self.assertTrue(units_temporal('hours since 2000-01-01'))
        self.assertFalse(units_temporal('hours'))
        self.assertFalse(units_temporal('days since the big bang'))




def breakpoint(scope=None, global_scope=None):
    import traceback
    from IPython.config.loader import Config
    ipy_config = Config()
    ipy_config.PromptManager.in_template = '><> '
    ipy_config.PromptManager.in2_template = '... '
    ipy_config.PromptManager.out_template = '--> '
    ipy_config.InteractiveShellEmbed.confirm_exit = False

    # First import the embeddable shell class
    from IPython.frontend.terminal.embed import InteractiveShellEmbed
    from mock import patch
    if scope is not None:
        locals().update(scope)
    if global_scope is not None:
        globals().update(global_scope)

    # Update namespace of interactive shell
    # TODO: Cleanup namespace even further
    # Now create an instance of the embeddable shell. The first argument is a
    # string with options exactly as you would type them if you were starting
    # IPython at the system command line. Any parameters you want to define for
    # configuration can thus be specified here.
    with patch("IPython.core.interactiveshell.InteractiveShell.init_virtualenv"):
        ipshell = InteractiveShellEmbed(config=ipy_config,
                                        banner1="Entering Breakpoint Shell",
                                        exit_msg = 'Returning...')

        stack = traceback.extract_stack(limit=2)
        message = 'File %s, line %s, in %s' % stack[0][:-1]

        try:
            import growl
            growl.growl('breakpoint', 'Ready')
        except:
            pass
        ipshell('(%s) Breakpoint @ %s' % ('breakpoint', message))
