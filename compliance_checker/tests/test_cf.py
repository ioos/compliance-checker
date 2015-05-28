#!/usr/bin/env python

from compliance_checker.cf import CFBaseCheck, BaseCheck, dimless_vertical_coordinates
from compliance_checker.cf.util import is_vertical_coordinate, is_time_variable, units_convertible
from compliance_checker.base import DSPair
from wicken.netcdf_dogma import NetCDFDogma
from netCDF4 import Dataset
from tempfile import gettempdir
from pkg_resources import resource_filename

import unittest
import os
import re


static_files = {
        'rutgers'              : resource_filename('compliance_checker', 'tests/data/ru07-20130824T170228_rt0.nc'),
        'conv_multi'           : resource_filename('compliance_checker', 'tests/data/conv_multi.nc'),
        'conv_bad'             : resource_filename('compliance_checker', 'tests/data/conv_bad.nc'),
        'example-grid'         : resource_filename('compliance_checker', 'tests/data/example-grid.nc'),
        'badname'              : resource_filename('compliance_checker', 'tests/data/non-comp/badname.netcdf'),
        'bad'                  : resource_filename('compliance_checker', 'tests/data/non-comp/bad.nc'),
        'dimensionless'        : resource_filename('compliance_checker', 'tests/data/dimensionless.nc'),
        '2dim'                 : resource_filename('compliance_checker', 'tests/data/2dim-grid.nc'),
        'bad2dim'              : resource_filename('compliance_checker', 'tests/data/non-comp/bad2dim.nc'),
        'rhgrid'               : resource_filename('compliance_checker', 'tests/data/rhgrid.nc'),
        'bad-rhgrid'           : resource_filename('compliance_checker', 'tests/data/non-comp/bad-rhgrid.nc'),
        'bad_data_type'        : resource_filename('compliance_checker', 'tests/data/bad_data_type.nc'),
        'mapping'              : resource_filename('compliance_checker', 'tests/data/mapping.nc'),
        'bad_region'           : resource_filename('compliance_checker', 'tests/data/bad_region.nc'),
        'featureType'          : resource_filename('compliance_checker', 'tests/data/example-grid.nc'),
        'cont_ragged'          : resource_filename('compliance_checker', 'tests/data/cont_ragged.nc'),
        'index_ragged'         : resource_filename('compliance_checker', 'tests/data/index_ragged.nc'),
        'bad_missing_data'     : resource_filename('compliance_checker', 'tests/data/bad_missing_data.nc'),
        'self-referencing-var' : resource_filename('compliance_checker', 'tests/data/self-referencing-var.nc')
        }

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

    #--------------------------------------------------------------------------------
    # Helper Methods
    #--------------------------------------------------------------------------------

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

    def get_pair(self, nc_dataset):
        '''
        Return a pairwise object for the dataset
        '''
        if isinstance(nc_dataset, basestring):
            nc_dataset = Dataset(nc_dataset, 'r')
            self.addCleanup(nc_dataset.close)
        dogma = NetCDFDogma('nc', self.cf.beliefs(), nc_dataset)
        pair = DSPair(nc_dataset, dogma)
        return pair

    #--------------------------------------------------------------------------------
    # Compliance Tests
    #--------------------------------------------------------------------------------

    def test_check_data_types(self):
        """
        2.2 The netCDF data types char, byte, short, int, float or real, and double are all acceptable
        """
        dataset = self.get_pair(static_files['rutgers'])
        result = self.cf.check_data_types(dataset)
        self.assertTrue(result.value)


        dpair = self.get_pair(static_files['bad_data_type'])
        result = self.cf.check_data_types(dpair)
        assert result.value == (5, 6)


    def test_naming_conventions(self):
        '''
        Section 2.3 Naming Conventions

        Variable, dimension and attribute names should begin with a letter and be composed of letters, digits, and underscores.
        '''
        dataset = self.get_pair(static_files['rutgers'])
        result = self.cf.check_naming_conventions(dataset)
        num_var = len(dataset.dataset.variables)
        
        expected = (num_var,) * 2
        self.assertEquals(result.value, expected)

        dataset = self.get_pair(static_files['bad'])
        result = self.cf.check_naming_conventions(dataset)
        num_var = len(dataset.dataset.variables)
        expected = (num_var-1, num_var)
        self.assertEquals(result.value, expected)
        assert '_poor_dim' in result.msgs [0]

    def test_check_names_unique(self):
        """
        2.3 names should not be distinguished purely by case, i.e., if case is disregarded, no two names should be the same.
        """
        dataset = self.get_pair(static_files['rutgers'])
        result = self.cf.check_names_unique(dataset)

        num_var = len(dataset.dataset.variables)
        expected = (num_var,) * 2

        self.assertEquals(result.value, expected)

        #TODO: Add bad unique names to bad.nc

    def test_check_dimension_names(self):
        """
        2.4 A variable may have any number of dimensions, including zero, and the dimensions must all have different names.
        """

        dataset = self.get_pair(static_files['bad_data_type'])
        result = self.cf.check_dimension_names(dataset)
        assert result.value == (5, 6)

    def test_check_dimension_order(self):
        """
        2.4 If any or all of the dimensions of a variable have the interpretations of "date or time" (T), "height or depth" (Z),
        "latitude" (Y), or "longitude" (X) then we recommend, those dimensions to appear in the relative order T, then Z, then Y,
        then X in the CDL definition corresponding to the file. All other dimensions should, whenever possible, be placed to the
        left of the spatiotemporal dimensions.
        """
        dataset = self.get_pair(static_files['bad_data_type'])
        result = self.cf.check_dimension_order(dataset)
        assert result.value == (11, 12)

    def test_check_fill_value_outside_valid_range(self):
        """
        2.5.1 The _FillValue should be outside the range specified by valid_range (if used) for a variable.
        """

        dataset = self.get_pair(static_files['bad_data_type'])
        result = self.cf.check_fill_value_outside_valid_range(dataset)
        assert result.value == (1, 2)

    def test_check_conventions_are_cf_16(self):
        """
        2.6.1 the NUG defined global attribute Conventions to the string value "CF-1.6"
        """
        # :Conventions = "CF-1.6"
        dataset = self.get_pair(static_files['rutgers'])
        result = self.cf.check_conventions_are_cf_16(dataset)
        self.assertTrue(result.value)

        # :Conventions = "CF-1.6 ,ACDD" ;
        dataset = self.get_pair(static_files['conv_multi'])
        result = self.cf.check_conventions_are_cf_16(dataset)
        self.assertTrue(result.value)

        # :Conventions = "NoConvention"
        dataset = self.get_pair(static_files['conv_bad'])
        result = self.cf.check_conventions_are_cf_16(dataset)
        self.assertFalse(result.value)

    def test_check_convention_globals(self):
        """
        2.6.2 title/history global attributes, must be strings. Do not need to exist.
        """
        #check for pass
        dataset = self.get_pair(static_files['rutgers'])
        result = self.cf.check_convention_globals(dataset)
        for each in result:
            self.assertTrue(each.value)
        #check if it doesn't exist that we pass
        dataset = self.get_pair(static_files['bad_data_type'])
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
        dataset = self.get_pair(static_files['rutgers'])
        result = self.cf.check_convention_possibly_var_attrs(dataset)
        for each in result:
            self.assertTrue(each.value)

        dataset = self.get_pair(static_files['bad_data_type'])
        result = self.cf.check_convention_possibly_var_attrs(dataset)
        for each in result:
            self.assertFalse(each.value)  

    def test_check_standard_name(self):
        """
        3.3 A standard name is associated with a variable via the attribute standard_name which takes a
        string value comprised of a standard name optionally followed by one or more blanks and a
        standard name modifier
        """
        dataset = self.get_pair(static_files['2dim'])
        result = self.cf.check_standard_name(dataset)
        for each in result:
            self.assertTrue(each.value)

        dataset = self.get_pair(static_files['bad_data_type'])
        result = self.cf.check_standard_name(dataset)
        for each in result:
            self.assertFalse(each.value)  




    def test_check_units(self):

        dataset = self.get_pair(static_files['2dim'])
        result = self.cf.check_units(dataset)
        for each in result:
            self.assertTrue(each.value)

        dataset = self.get_pair(static_files['bad_data_type'])
        result = self.cf.check_units(dataset)
        for each in result:
            self.assertFalse(each.value)  



    def test_coordinate_types(self):
        '''
        Section 4 Coordinate Types

        We strongly recommend that coordinate variables be used for all coordinate types whenever they are applicable.
        '''
        dataset = self.get_pair(static_files['bad_data_type'])
        result = self.cf.check_coordinate_vars_for_all_coordinate_types(dataset)
        for each in result:
            print each
            self.assertTrue(each.value)

    def test_check_coordinate_axis_attr(self):

        dataset = self.get_pair(static_files['2dim'])
        result = self.cf.check_coordinate_axis_attr(dataset)
        for each in result:
            self.assertTrue(each.value)

        dataset = self.get_pair(static_files['bad_data_type'])
        result = self.cf.check_coordinate_axis_attr(dataset)
        for each in result:
            print each
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
        dataset = self.get_pair(static_files['example-grid'])
        results = self.cf.check_latitude(dataset)
        for r in results:
            if isinstance(r.value, tuple):
                self.assertEquals(r.value[0], r.value[1])
            else:
                self.assertTrue(r.value)
        
        # Verify non-compliance
        dataset = self.get_pair(static_files['bad'])
        results = self.cf.check_latitude(dataset)
        # Store the results in a dict
        rd = {}
        for r in results:
            rd[r.name[1:]] = r.value
        # ('lat', 'has_units') should be False
        self.assertFalse(rd[('lat', 'has_units')])
        # ('lat', 'correct_units') should be (0,3)
        self.assertEquals(rd[('lat', 'correct_units')], (0,3))
        # ('lat_uv', 'has_units') should be True
        self.assertTrue(rd[('lat_uv', 'has_units')])
        # ('lat_uv', 'correct_units') should be (2,3)
        self.assertEquals(rd[('lat_uv', 'correct_units')], (2,3))
        # ('lat_like', 'has_units') should be True
        self.assertTrue(rd[('lat_like', 'has_units')])
        # ('lat_like', 'correct_units') should be (1,3)
        self.assertEquals(rd[('lat_like', 'correct_units')], (1,3))
        

    def test_longitude(self):
        '''
        Section 4.2 Longitude Coordinate
        '''
        # Check compliance
        dataset = self.get_pair(static_files['example-grid'])
        results = self.cf.check_longitude(dataset)
        for r in results:
            if isinstance(r.value, tuple):
                self.assertEquals(r.value[0], r.value[1])
            else:
                self.assertTrue(r.value)
        
        # Verify non-compliance
        dataset = self.get_pair(static_files['bad'])
        results = self.cf.check_longitude(dataset)
        # Store the results in a dict
        rd = {}
        for r in results:
            rd[r.name[1:]] = r.value
        # ('lon', 'has_units') should be False
        self.assertFalse(rd[('lon', 'has_units')])
        # ('lon', 'correct_units') should be (0,3)
        self.assertEquals(rd[('lon', 'correct_units')], (0,3))
        # ('lon_uv', 'has_units') should be True
        self.assertTrue(rd[('lon_uv', 'has_units')])
        # ('lon_uv', 'correct_units') should be (2,3)
        self.assertEquals(rd[('lon_uv', 'correct_units')], (2,3))
        # ('lon_like', 'has_units') should be True
        self.assertTrue(rd[('lon_like', 'has_units')])
        # ('lon_like', 'correct_units') should be (1,3)
        self.assertEquals(rd[('lon_like', 'correct_units')], (1,3))

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

        dataset = self.get_pair(static_files['example-grid'])
        results = self.cf.check_vertical_coordinate(dataset)
        for r in results:
            self.assertTrue(r.value)

        # Check non-compliance
        dataset = self.get_pair(static_files['bad'])
        results = self.cf.check_vertical_coordinate(dataset)
        
        # Store the results by the tuple
        rd = { r.name[1:] : r.value for r in results }
        # ('height', 'has_units') should be False
        self.assertFalse(rd[('height', 'has_units')])
        # ('height', 'correct_units') should be False
        self.assertFalse(rd[('height', 'correct_units')])
        # ('depth', 'has_units') should be True
        self.assertTrue(rd[('depth', 'has_units')])
        # ('depth', 'correct_units') should be False
        self.assertFalse(rd[('depth', 'correct_units')])
        # ('depth2', 'has_units') should be False
        self.assertTrue(rd[('depth2', 'has_units')])
        # ('depth2', 'correct_units') should be False
        self.assertFalse(rd[('depth2', 'correct_units')])
        

    def test_vertical_dimension(self):
        '''
        Section 4.3.1 Dimensional Vertical Coordinate
        '''
        # Check for compliance
        dataset = self.get_pair(static_files['example-grid'])
        results = self.cf.check_dimensional_vertical_coordinate(dataset)
        for r in results:
            self.assertTrue(r.value)

        # Check for non-compliance
        dataset = self.get_pair(static_files['bad'])
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
        dataset = self.get_pair(static_files['dimensionless'])
        results = self.cf.check_dimensionless_vertical_coordinate(dataset)
        for r in results:
            self.assertTrue(r.value)

        # Check negative compliance
        dataset = self.get_pair(static_files['bad'])
        results = self.cf.check_dimensionless_vertical_coordinate(dataset)
        rd = { r.name[1:] : r.value for r in results }
        
        # ('lev1', 'formula_terms') should be False
        self.assertFalse(rd[('lev1', 'formula_terms')])
        
        # ('lev2', 'formula_terms') should be True
        self.assertTrue(rd[('lev2', 'formula_terms')])
        # ('lev2', 'terms_exist') should be False
        self.assertFalse(rd[('lev2', 'terms_exist')])
            
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
        dataset = self.get_pair(static_files['example-grid'])
        results = self.cf.check_time_coordinate(dataset)
        for r in results:
            self.assertTrue(r.value)


        dataset = self.get_pair(static_files['bad'])
        results = self.cf.check_time_coordinate(dataset)
        rd = {r.name[1:] : r.value for r in results }
        self.assertFalse(rd[('bad_time_1', 'has_units')])
        self.assertTrue(rd[('bad_time_2', 'has_units')])
        self.assertFalse(rd[('bad_time_2', 'correct_units')])

    def test_check_calendar(self):
        dataset = self.get_pair(static_files['example-grid'])
        results = self.cf.check_calendar(dataset)
        for r in results:
            self.assertTrue(r.value)

        dataset = self.get_pair(static_files['bad'])
        results = self.cf.check_calendar(dataset)
        rd = {r.name[1:] : r.value for r in results }
        self.assertFalse(rd[('bad_time_1', 'has_calendar')])
        self.assertFalse(rd[('bad_time_1', 'valid_calendar')])
        self.assertTrue(rd[('bad_time_2', 'has_calendar')])
        self.assertFalse(rd[('bad_time_2', 'valid_calendar')])

    def test_check_independent_axis_dimensions(self):
        dataset = self.get_pair(static_files['example-grid'])
        results = self.cf.check_independent_axis_dimensions(dataset)
        for r in results:
            self.assertTrue(r.value)

        dataset = self.get_pair(static_files['bad'])
        results = self.cf.check_independent_axis_dimensions(dataset)
        
        false_variable_name_list = ['lev1', 'lev2', 'bad_time_1', 'bad_time_2', 'column_temp']
        
        for each in results:
            if each.name[1] in false_variable_name_list:
                self.assertFalse(each.value)
            else:
                self.assertTrue(each.value)

    def test_check_two_dimensional(self):
        dataset = self.get_pair(static_files['2dim'])
        results = self.cf.check_two_dimensional(dataset)
        for r in results:
            self.assertTrue(r.value)


        # Need the bad testing
        dataset = self.get_pair(static_files['bad2dim'])
        results = self.cf.check_two_dimensional(dataset)
        self.assertTrue(results[0].value)
        self.assertFalse(results[1].value)
        self.assertFalse(results[2].value)
        self.assertTrue(results[3].value)
        self.assertFalse(results[4].value)
        self.assertTrue(results[5].value)

        # Test the self referencing variables
        dataset = self.get_pair(static_files['self-referencing-var'])
        try:
            results = self.cf.check_two_dimensional(dataset)
            self.assertFalse(results[0].value)
        except:
            self.assertTrue(False)

    def test_check_reduced_horizontal_grid(self):
        dataset = self.get_pair(static_files['rhgrid'])
        results = self.cf.check_reduced_horizontal_grid(dataset)
        rd = { r.name[1] : r.value for r in results }
        self.assertTrue(rd['PS'])

        dataset = self.get_pair(static_files['bad-rhgrid'])
        results = self.cf.check_reduced_horizontal_grid(dataset)
        rd = { r.name[1] : (r.value, r.msgs) for r in results }

        for name, (value, msg) in rd.iteritems():
            self.assertFalse(value)

        self.assertIn('Coordinate longitude is not a proper variable', rd['PSa'][1])
        self.assertIn("Coordinate latitude's dimension, latdim, is not a dimension of PSb", rd['PSb'][1])
        assert 'PSc' not in rd.keys()


    def test_check_horz_crs_grid_mappings_projections(self):
        dataset = self.get_pair(static_files['mapping'])
        results = self.cf.check_horz_crs_grid_mappings_projections(dataset)
        rd = { r.name[1] : r.value for r in results }
        assert rd['wgs84'] == (3, 3)
        assert rd['epsg']  == (7, 8)


    def test_check_scalar_coordinate_system(self):
        dataset = self.get_pair(static_files['bad_data_type'])
        results = self.cf.check_scalar_coordinate_system(dataset)
        assert results[0].value == (1, 2)

    def test_check_geographic_region(self):
        dataset = self.get_pair(static_files['bad_region'])
        results = self.cf.check_geographic_region(dataset)

        self.assertFalse(results[0].value)
        self.assertTrue(results[1].value)

    def test_check_alternative_coordinates(self):
        dataset = self.get_pair(static_files['bad_data_type'])
        results = self.cf.check_alternative_coordinates(dataset)
        self.assertTrue(results[0].value)


    #def test_check_cell_boundaries(self):
    #    dataset = self.get_pair(static_files['bad_data_type'])
    #    results = self.cf.check_cell_boundaries(dataset)
    #    print results
    #    self.assertTrue(results[0].value)


    def test_check_packed_data(self):
        dataset = self.get_pair(static_files['bad_data_type'])
        results = self.cf.check_packed_data(dataset)
        self.assertFalse(results[0].value)
        self.assertTrue(results[1].value)


    def test_check_compression(self):
        dataset = self.get_pair(static_files['bad_data_type'])
        results = self.cf.check_compression(dataset)
        assert results[0].value == (2,2)
        assert results[1].value == (0,2)

    def test_check_all_features_are_same_type(self):
        dataset = self.get_pair(static_files['rutgers'])
        results = self.cf.check_all_features_are_same_type(dataset)
        assert results == None

        dataset = self.get_pair(static_files['featureType'])
        results = self.cf.check_all_features_are_same_type(dataset)
        self.assertTrue(results.value)     

        dataset = self.get_pair(static_files['bad_data_type'])
        results = self.cf.check_all_features_are_same_type(dataset)
        self.assertFalse(results.value)   

    def test_check_orthogonal_multidim_array(self):
        dataset = self.get_pair(static_files['rutgers'])
        results = self.cf.check_orthogonal_multidim_array(dataset)
        for each in results:
            self.assertTrue(each.value)

    def test_check_incomplete_multidim_array(self):
        dataset = self.get_pair(static_files['bad_data_type'])
        results = self.cf.check_incomplete_multidim_array(dataset)
        for each in results:
            self.assertTrue(each.value)

    def test_check_contiguous_ragged_array(self):
        dataset = self.get_pair(static_files['cont_ragged'])
        results = self.cf.check_contiguous_ragged_array(dataset)
        for each in results:
            self.assertTrue(each.value)

    def test_check_indexed_ragged_array(self):
        dataset = self.get_pair(static_files['index_ragged'])
        results = self.cf.check_indexed_ragged_array(dataset)
        for each in results:
            self.assertTrue(each.value)


    def test_check_feature_type(self):
        dataset = self.get_pair(static_files['index_ragged'])
        results = self.cf.check_feature_type(dataset)
        self.assertTrue(results.value)

        dataset = self.get_pair(static_files['bad_data_type'])
        results = self.cf.check_feature_type(dataset)
        self.assertFalse(results.value)



    def test_check_coordinates_and_metadata(self):
        dataset = self.get_pair(static_files['bad_data_type'])
        results = self.cf.check_coordinates_and_metadata(dataset)
        self.assertFalse(results[0].value)
        self.assertTrue(results[1].value)
        self.assertFalse(results[2].value)

        dataset = self.get_pair(static_files['index_ragged'])
        results = self.cf.check_coordinates_and_metadata(dataset)
        self.assertTrue(results[-1].value)

    def test_check_missing_data(self):
        dataset = self.get_pair(static_files['index_ragged'])
        results = self.cf.check_missing_data(dataset)
        for each in results:
            self.assertTrue(each.value)

        dataset = self.get_pair(static_files['bad_missing_data'])
        results = self.cf.check_missing_data(dataset)
        for each in results:
            self.assertFalse(each.value)

    #--------------------------------------------------------------------------------
    # Utility Method Tests
    #--------------------------------------------------------------------------------

    def test_temporal_unit_conversion(self):
        self.assertTrue(units_convertible('hours', 'seconds'))
        self.assertFalse(units_convertible('hours', 'hours since 2000-01-01'))



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
