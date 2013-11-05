#!/usr/bin/env python

from compliance_checker.cf import CFCheck, BaseCheck
from compliance_checker.suite import DSPair, NetCDFDogma
from netCDF4 import Dataset
from tempfile import gettempdir

import unittest
import os


static_files = {
        'rutgers'      : 'test-data/ru07-20130824T170228_rt0.nc',
        'example-grid' : 'test-data/example-grid.nc',
        'badname'      : 'test-data/non-comp/badname.netcdf',
        'bad'          : 'test-data/non-comp/bad.nc',
        }


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
        self.cf = CFCheck()

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

    def test_filename(self):
        '''
        Section 2.1 Filenames

        NetCDF files should have the file name extension
        '''
        dataset = self.get_pair(static_files['rutgers'])
        result = self.cf.check_filename_extension(dataset)
        self.assertTrue(result.value)


        dpair = self.get_pair(static_files['badname'])
        result = self.cf.check_filename_extension(dpair)
        # Verify that the non-compliant file returns a negative result
        self.assertFalse(result.value)

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
        self.assertEquals(result.msgs, ['_poor_dim'])


        
    def test_coordinate_types(self):
        '''
        Section 4 Coordinate Types

        We strongly recommend that coordinate variables be used for all coordinate types whenever they are applicable.
        '''

        pass

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
        mock_variable = type('MockVariable', (object,), {})

        # Proper name/standard_name
        known_name = mock_variable()
        known_name.standard_name = 'depth'
        self.assertTrue(self.cf._is_vertical_coordinate('not_known', known_name))

        # Proper Axis
        axis_set = mock_variable()
        axis_set.axis = 'Z'
        self.assertTrue(self.cf._is_vertical_coordinate('not_known', axis_set))

        # Proper units
        units_set = mock_variable()
        units_set.units = 'dbar'
        self.assertTrue(self.cf._is_vertical_coordinate('not_known', units_set))

        # Proper units/positive
        positive = mock_variable()
        positive.units = 'm'
        positive.positive = 'up'
        self.assertTrue(self.cf._is_vertical_coordinate('not_known', positive))

    def test_vertical_coordinate(self):
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
        

