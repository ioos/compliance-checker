from compliance_checker.glider_dac import GliderCheck
from pkg_resources import resource_filename
from netCDF4 import Dataset
from compliance_checker.base import DSPair
from wicken.netcdf_dogma import NetCDFDogma
import unittest

static_files = {
    'glider_std': resource_filename('compliance_checker', 'tests/data/gliders/IOOS_Glider_NetCDF_v2.0.nc'),
    'bad_location': resource_filename('compliance_checker', 'tests/data/gliders/bad_location.nc'),
    'bad_qc': resource_filename('compliance_checker', 'tests/data/gliders/bad_qc.nc'),
    'bad_metadata': resource_filename('compliance_checker', 'tests/data/gliders/bad_metadata.nc'),
}

class TestGliderCheck(unittest.TestCase):
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
    
    def get_pair(self, nc_dataset):
        '''
        Return a pairwise object for the dataset
        '''
        print nc_dataset
        if isinstance(nc_dataset, basestring):
            nc_dataset = Dataset(nc_dataset, 'r')
            self.addCleanup(nc_dataset.close)
        dogma = NetCDFDogma('nc', self.check.beliefs(), nc_dataset)
        pair = DSPair(nc_dataset, dogma)
        return pair
    
    def setUp(self):
        self.check = GliderCheck()
    
    def test_location(self):
        '''
        Checks that a file with the proper lat and lon do work
        '''
        dataset = self.get_pair(static_files['glider_std'])
        result = self.check.check_locations(dataset)
        self.assertTrue(result.value)

    def test_location_fail(self):
        '''
        Ensures the checks fail for location
        '''
        dataset = self.get_pair(static_files['bad_location'])
        result = self.check.check_locations(dataset)
        self.assertEquals(result.value, (0,1))

    def test_ctd_fail(self):
        '''
        Ensures the ctd checks fail for temperature
        '''
        dataset = self.get_pair(static_files['bad_qc'])
        result = self.check.check_ctd_variables(dataset)
        self.assertEquals(result.value, (55,56))
    
    def test_ctd_vars(self):
        '''
        Ensures the ctd checks for the correct file
        '''
        dataset = self.get_pair(static_files['glider_std'])
        result = self.check.check_ctd_variables(dataset)
        self.assertEquals(result.value, (56,56))

    def test_global_fail(self):
        '''
        Tests that the global checks fail where appropriate
        '''
        dataset = self.get_pair(static_files['bad_qc'])
        result = self.check.check_global_attributes(dataset)
        self.assertEquals(result.value, (41,64))

    def test_global(self):
        '''
        Tests that the global checks work
        '''
        dataset = self.get_pair(static_files['glider_std'])
        result = self.check.check_global_attributes(dataset)
        self.assertEquals(result.value, (64,64))


    def test_metadata(self):
        '''
        Tests that empty attributes fail
        '''
        dataset = self.get_pair(static_files['bad_metadata'])
        result = self.check.check_global_attributes(dataset)
        self.assertEquals(result.value, (41,64))

    def test_standard_names(self):
        '''
        Tests that the standard names succeed/fail
        '''

        dataset = self.get_pair(static_files['bad_metadata'])
        result = self.check.check_standard_names(dataset)
        self.assertEquals(result.value, (0, 1))
        
        dataset = self.get_pair(static_files['glider_std'])
        result = self.check.check_standard_names(dataset)
        self.assertEquals(result.value, (1, 1))
