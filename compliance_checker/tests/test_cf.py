#!/usr/bin/env python

from compliance_checker.suite import CheckSuite
from compliance_checker.cf import CFBaseCheck, dimless_vertical_coordinates
from compliance_checker.cf.util import is_vertical_coordinate, is_time_variable, units_convertible, units_temporal, StandardNameTable, create_cached_data_dir, download_cf_standard_name_table
from netCDF4 import Dataset
from tempfile import gettempdir
from compliance_checker.tests.resources import STATIC_FILES

import unittest
import os
import re


class MockVariable(object):
    '''
    For mocking a dataset variable
    '''
    pass


class TestCF(unittest.TestCase):
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
        self.assertTrue(result.value)

        dpair = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_data_types(dpair)
        assert result.value == (5, 6)

    def test_naming_conventions(self):
        '''
        Section 2.3 Naming Conventions

        Variable, dimension and attribute names should begin with a letter and be composed of letters, digits, and underscores.
        '''
        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        result = self.cf.check_naming_conventions(dataset)
        num_var = len(dataset.variables)

        expected = (num_var,) * 2
        self.assertEqual(result.value, expected)

        dataset = self.load_dataset(STATIC_FILES['bad'])
        result = self.cf.check_naming_conventions(dataset)
        num_var = len(dataset.variables)
        expected = (num_var - 1, num_var)
        self.assertEqual(result.value, expected)
        assert '_poor_dim' in result.msgs[0]

    def test_check_names_unique(self):
        """
        2.3 names should not be distinguished purely by case, i.e., if case is disregarded, no two names should be the same.
        """
        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        result = self.cf.check_names_unique(dataset)

        num_var = len(dataset.variables)
        expected = (num_var,) * 2

        self.assertEqual(result.value, expected)

        # TODO: Add bad unique names to bad.nc

    def test_check_dimension_names(self):
        """
        2.4 A variable may have any number of dimensions, including zero, and the dimensions must all have different names.
        """

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_dimension_names(dataset)
        assert result.value == (5, 6)

    def test_check_dimension_order(self):
        """
        2.4 If any or all of the dimensions of a variable have the interpretations of "date or time" (T), "height or depth" (Z),
        "latitude" (Y), or "longitude" (X) then we recommend, those dimensions to appear in the relative order T, then Z, then Y,
        then X in the CDL definition corresponding to the file. All other dimensions should, whenever possible, be placed to the
        left of the spatiotemporal dimensions.
        """
        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_dimension_order(dataset)
        assert result.value == (11, 12)

    def test_check_fill_value_outside_valid_range(self):
        """
        2.5.1 The _FillValue should be outside the range specified by valid_range (if used) for a variable.
        """

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        results = self.cf.check_fill_value_outside_valid_range(dataset)
        assert sum((result.value for result in results)) == 1
        assert len(results) == 2

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

    def test_check_convention_globals(self):
        """
        2.6.2 title/history global attributes, must be strings. Do not need to exist.
        """
        # check for pass
        dataset = self.load_dataset(STATIC_FILES['rutgers'])
        result = self.cf.check_convention_globals(dataset)
        for each in result:
            self.assertTrue(each.value)
        # check if it doesn't exist that we pass
        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_convention_globals(dataset)
        for each in result:
            self.assertTrue(each.value)

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
        for each in result:
            self.assertTrue(each.value)

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_convention_possibly_var_attrs(dataset)
        for each in result:
            self.assertFalse(each.value)

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
        result = self.cf.check_standard_name(dataset)
        for each in result:
            self.assertFalse(each.value)

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
        dataset = self.load_dataset(STATIC_FILES['self_referencing'])
        results = self.cf.check_flags(dataset)
        scored, out_of, messages = self.get_results(results)

        self.assertEqual(scored, 46)
        self.assertEqual(out_of, 59)
        self.assertEqual(messages.count('flag_values must be a list'), 6)
        m_str = r"'flag_values' attribute for variable '\w+' does not have same type \(fv: [<>]?\w+, v: [<>]?\w+\)"
        # make sure flag_values attribute where not equal to variable type
        # has the proper message
        self.assertEqual(sum(bool(re.match(m_str, msg)) for msg in messages), 7)

    def test_check_bad_units(self):

        dataset = self.load_dataset(STATIC_FILES['2dim'])
        result = self.cf.check_units(dataset)
        for each in result:
            self.assertTrue(each.value)

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_units(dataset)
        for each in result:
            self.assertFalse(each.value)

    def test_coordinate_types(self):
        '''
        Section 4 Coordinate Types

        We strongly recommend that coordinate variables be used for all coordinate types whenever they are applicable.
        '''
        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_coordinate_vars_for_all_coordinate_types(dataset)
        for each in result:
            self.assertTrue(each.value)

    def test_check_coordinate_axis_attr(self):

        dataset = self.load_dataset(STATIC_FILES['2dim'])
        result = self.cf.check_coordinate_axis_attr(dataset)
        for each in result:
            self.assertTrue(each.value)

        dataset = self.load_dataset(STATIC_FILES['bad_data_type'])
        result = self.cf.check_coordinate_axis_attr(dataset)
        for each in result:
            if each.name[1] in ['time', 'latitude']:
                self.assertTrue(each.value)
            if each.name[1] in ['salinity']:
                if each.name[2] not in ['does_not_depend_on_mult_coord_vars']:
                    self.assertFalse(each.value)

    def test_latitude(self):
        '''
        Section 4.1 Latitude Coordinate
        '''
        # Check compliance
        dataset = self.load_dataset(STATIC_FILES['example-grid'])
        results = self.cf.check_latitude(dataset)
        for r in results:
            if isinstance(r.value, tuple):
                self.assertEqual(r.value[0], r.value[1])
            else:
                self.assertTrue(r.value)

        # Verify non-compliance
        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_latitude(dataset)

        scored, out_of, messages = self.get_results(results)

        assert 'lat does not have units attribute' in messages
        assert 'lat_uv units are acceptable, but not recommended' in messages
        assert 'lat_like does not have units attribute' in messages

        assert scored == 5
        assert out_of == 12

    def test_longitude(self):
        '''
        Section 4.2 Longitude Coordinate
        '''
        # Check compliance
        dataset = self.load_dataset(STATIC_FILES['example-grid'])
        results = self.cf.check_longitude(dataset)
        for r in results:
            if isinstance(r.value, tuple):
                self.assertEqual(r.value[0], r.value[1])
            else:
                self.assertTrue(r.value)

        # Verify non-compliance
        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_longitude(dataset)

        scored, out_of, messages = self.get_results(results)

        assert 'lon does not have units attribute' in messages
        assert 'lon_uv units are acceptable, but not recommended' in messages
        assert 'lon_like does not have units attribute' in messages

        assert scored == 5
        assert out_of == 12

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

    def test_vertical_coordinate(self):
        '''
        Section 4.3 Vertical (Height or Depth) coordinate
        '''
        # Check compliance

        dataset = self.load_dataset(STATIC_FILES['example-grid'])
        results = self.cf.check_vertical_coordinate(dataset)
        for r in results:
            self.assertTrue(r.value)

        # Check non-compliance
        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_vertical_coordinate(dataset)

        scored, out_of, messages = self.get_results(results)

        assert 'height does not have units' in messages
        assert 'vertical variable depth needs to define positive attribute'
        assert 'vertical variable depth2 needs to define positive attribute'

    def test_vertical_dimension(self):
        '''
        Section 4.3.1 Dimensional Vertical Coordinate
        '''
        # Check for compliance
        dataset = self.load_dataset(STATIC_FILES['example-grid'])
        results = self.cf.check_dimensional_vertical_coordinate(dataset)
        for r in results:
            self.assertTrue(r.value)

        # Check for non-compliance
        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_dimensional_vertical_coordinate(dataset)
        for r in results:
            self.assertFalse(r.value)

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
        for r in results:
            self.assertTrue(r.value)

        # Check negative compliance
        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_dimensionless_vertical_coordinate(dataset)

        scored, out_of, messages = self.get_results(results)

        assert 'formula_terms missing from dimensionless coordinate lev1' in messages
        assert 'formula_terms not defined for dimensionless coordinate lev1' in messages
        assert 'var1 missing for dimensionless coordinate lev2' in messages
        assert 'var2 missing for dimensionless coordinate lev2' in messages
        assert 'var3 missing for dimensionless coordinate lev2' in messages
        assert scored == 1
        assert out_of == 4

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

        assert 'bad_time_1 does not have units' in messages
        assert 'bad_time_2 doesn not have correct time units' in messages
        assert scored == 1
        assert out_of == 3

    def test_check_calendar(self):
        dataset = self.load_dataset(STATIC_FILES['example-grid'])
        results = self.cf.check_calendar(dataset)
        for r in results:
            self.assertTrue(r.value)

        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_calendar(dataset)
        scored, out_of, messages = self.get_results(results)

        assert 'Variable bad_time_1 should have a calendar attribute' in messages
        assert "Variable bad_time_2 should have a valid calendar: 'nope' is not a valid calendar" in messages

    def test_self_referencing(self):
        '''
        This test captures a check where a coordinate has circular references
        '''
        dataset = self.load_dataset(STATIC_FILES['self_referencing'])
        results = self.cf.check_two_dimensional(dataset)

        scored, out_of, messages = self.get_results(results)
        assert "Variable LATITUDE's coordinate references itself" in messages
        assert scored == 1
        assert out_of == 3

        dataset = self.load_dataset(STATIC_FILES['valid_coordinates'])
        results = self.cf.check_two_dimensional(dataset)
        scored, out_of, messages = self.get_results(results)
        assert out_of == 4
        assert scored == 4
        assert "Variable CD_310's coordinate references itself" not in messages

    def test_check_independent_axis_dimensions(self):
        dataset = self.load_dataset(STATIC_FILES['example-grid'])
        results = self.cf.check_independent_axis_dimensions(dataset)
        for r in results:
            self.assertTrue(r.value)

        dataset = self.load_dataset(STATIC_FILES['bad'])
        results = self.cf.check_independent_axis_dimensions(dataset)

        scored, out_of, messages = self.get_results(results)
        assert 'The lev dimension for the variable lev1 does not have an associated coordinate variable, but is a Lat/Lon/Time/Height dimension.' \
            in messages
        assert 'The lev dimension for the variable lev2 does not have an associated coordinate variable, but is a Lat/Lon/Time/Height dimension.' \
            in messages
        assert 'The time dimension for the variable bad_time_1 does not have an associated coordinate variable, but is a Lat/Lon/Time/Height dimension.' \
            in messages
        assert 'The time dimension for the variable bad_time_2 does not have an associated coordinate variable, but is a Lat/Lon/Time/Height dimension.' \
            in messages
        assert 'The time dimension for the variable column_temp does not have an associated coordinate variable, but is a Lat/Lon/Time/Height dimension.' \
            in messages
        assert scored == 6
        assert out_of == 11

    def test_check_two_dimensional(self):
        dataset = self.load_dataset(STATIC_FILES['2dim'])
        results = self.cf.check_two_dimensional(dataset)
        for r in results:
            self.assertTrue(r.value)
        # Need the bad testing
        dataset = self.load_dataset(STATIC_FILES['bad2dim'])
        results = self.cf.check_two_dimensional(dataset)

        scored, out_of, messages = self.get_results(results)

        assert "Variable T's coordinate, lat, is not a coordinate or auxiliary variable" in messages
        assert "coordinate lat is not a correct lat/lon variable" in messages
        assert "Variable C's coordinate, lat_p, does not share dimension x with the variable" in messages

    def test_check_reduced_horizontal_grid(self):
        dataset = self.load_dataset(STATIC_FILES['rhgrid'])
        results = self.cf.check_reduced_horizontal_grid(dataset)
        rd = { r.name[1] : r.value for r in results }
        self.assertTrue(rd['PS'])

        dataset = self.load_dataset(STATIC_FILES['bad-rhgrid'])
        results = self.cf.check_reduced_horizontal_grid(dataset)
        rd = { r.name[1] : (r.value, r.msgs) for r in results }

        for name, (value, msg) in rd.items():
            self.assertFalse(value)

        self.assertIn('Coordinate longitude is not a proper variable', rd['PSa'][1])
        self.assertIn("Coordinate latitude's dimension, latdim, is not a dimension of PSb", rd['PSb'][1])
        assert 'PSc' not in list(rd.keys())

    def test_check_horz_crs_grid_mappings_projections(self):
        dataset = self.load_dataset(STATIC_FILES['mapping'])
        results = self.cf.check_horz_crs_grid_mappings_projections(dataset)
        rd = { r.name[1] : r.value for r in results }
        assert rd['wgs84'] == (3, 3)
        assert rd['epsg']  == (7, 8)

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
        assert scored == 4
        assert out_of == 4
        assert messages == []

    def test_64bit(self):
        dataset = self.load_dataset(STATIC_FILES['ints64'])
        suite = CheckSuite()
        suite.checkers = {
            'cf'        : CFBaseCheck
        }
        suite.run(dataset, 'cf')

    def test_time_units(self):
        dataset = self.load_dataset(STATIC_FILES['time_units'])
        results = self.cf.check_units(dataset)
        scored, out_of, messages = self.get_results(results)
        assert 'units are days since 1970-01-01, standard_name units should be K' in messages
        assert scored == 1
        assert out_of == 2

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
