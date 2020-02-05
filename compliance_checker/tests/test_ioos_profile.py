from compliance_checker.ioos import IOOS0_1Check, IOOS1_1Check, IOOS1_2Check
from compliance_checker.tests.resources import STATIC_FILES
from compliance_checker.tests import BaseTestCase
from compliance_checker.tests.helpers import MockTimeSeries, MockVariable
from compliance_checker.tests.test_cf import get_results
from netCDF4 import Dataset
import numpy as np
import os


class TestIOOS0_1(BaseTestCase):
    '''
    Tests for the IOOS Inventory Metadata v0.1
    '''

    def setUp(self):
        # Use the NCEI Gold Standard Point dataset for IOOS checks
        self.ds = self.load_dataset(STATIC_FILES['ncei_gold_point_1'])

        self.ioos = IOOS0_1Check()

    def test_cc_meta(self):
        assert self.ioos._cc_spec == 'ioos'
        assert self.ioos._cc_spec_version == '0.1'

    def test_global_attributes(self):
        '''
        Tests that all global attributes checks are working
        '''

        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        nc_obj = Dataset(os.devnull, 'w', diskless=True)
        self.addCleanup(nc_obj.close)

        results = self.ioos.check_global_attributes(nc_obj)
        for result in results:
            self.assert_result_is_bad(result)

        attrs = [
            'acknowledgement',
            'publisher_email',
            'institution',
            'publisher_name',
            'Conventions'
        ]
        for attr in attrs:
            setattr(nc_obj, attr, 'test')

        results = self.ioos.check_global_attributes(nc_obj)
        for result in results:
            self.assert_result_is_good(result)

    def test_variable_attributes(self):
        '''
        Tests that the platform variable attributes check is working
        '''

        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        nc_obj = Dataset(os.devnull, 'w', diskless=True)
        self.addCleanup(nc_obj.close)

        # The dataset needs at least one variable to check that it's missing
        # all the required attributes.
        nc_obj.createDimension('time', 1)
        nc_obj.createVariable('platform', 'S1', ())

        platform = nc_obj.variables['platform']

        results = self.ioos.check_variable_attributes(nc_obj)
        for result in results:
            self.assert_result_is_bad(result)

        platform.long_name = 'platform'
        platform.short_name = 'platform'
        platform.source = 'glider'
        platform.ioos_name = 'urn:ioos:station:glos:leorgn'
        platform.wmo_id = '1234'
        platform.comment = 'test'

        results = self.ioos.check_variable_attributes(nc_obj)
        for result in results:
            self.assert_result_is_good(result)

    def test_variable_units(self):
        '''
        Tests that the variable units test is working
        '''

        # this check tests that units attribute is present on EVERY variable

        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        nc_obj = Dataset(os.devnull, 'w', diskless=True)
        self.addCleanup(nc_obj.close)

        # The dataset needs at least one variable to check that it's missing
        # all the required attributes.
        nc_obj.createDimension('time', 1)
        nc_obj.createVariable('sample_var', 'd', ('time',))

        sample_var = nc_obj.variables['sample_var']

        results = self.ioos.check_variable_units(nc_obj)
        self.assert_result_is_bad(results)

        sample_var.units = 'm'
        sample_var.short_name = 'sample_var'

        results = self.ioos.check_variable_units(nc_obj)
        self.assert_result_is_good(results)

    def test_altitude_units(self):
        '''
        Tests that the altitude variable units test is working
        '''

        results = self.ioos.check_altitude_units(self.ds)
        self.assert_result_is_good(results)

        # Now test an nc file with a 'z' variable without units
        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        nc_obj = Dataset(os.devnull, 'w', diskless=True)
        self.addCleanup(nc_obj.close)

        # The dataset needs at least one variable to check that it's missing
        # all the required attributes.
        nc_obj.createDimension('time', 1)
        nc_obj.createVariable('z', 'd', ('time',))
        z = nc_obj.variables['z']
        z.short_name = 'sample_var'

        results = self.ioos.check_variable_units(nc_obj)
        self.assert_result_is_bad(results)


class TestIOOS1_1(BaseTestCase):
    '''
    Tests for the compliance checker implementation of IOOS Metadata Profile
    for NetCDF, Version 1.1
    '''

    def setUp(self):
        # Use the IOOS 1_1 dataset for testing
        self.ds = self.load_dataset(STATIC_FILES['ioos_gold_1_1'])

        self.ioos = IOOS1_1Check()

    def test_cc_meta(self):
        assert self.ioos._cc_spec == 'ioos'
        assert self.ioos._cc_spec_version == '1.1'

    def test_required_attributes(self):
        '''
        Tests that required attributes test is working properly
        '''

        results = self.ioos.check_high(self.ds)
        for result in results:
            self.assert_result_is_good(result)

    def test_recomended_attributes(self):
        '''
        Tests that recommended attributes test is working properly
        '''

        results = self.ioos.check_recommended(self.ds)
        for result in results:
            self.assert_result_is_good(result)

    def test_bad_platform_variables(self):
        '''
        Tests that the platform variable attributes check is working
        '''

        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        nc_obj = Dataset(os.devnull, 'w', diskless=True)
        self.addCleanup(nc_obj.close)

        # The dataset needs at least one variable to check that it's missing
        # all the required attributes.
        nc_obj.createDimension('time', 1)
        nc_obj.platform = 'platform'
        # global attribute 'platform' points to variable that does not exist in dataset

        results = self.ioos.check_platform_variables(nc_obj)
        for result in results:
            self.assert_result_is_bad(result)

    def test_good_platform_variables(self):
        '''
        Tests that the platform variable attributes check is working
        '''

        results = self.ioos.check_platform_variables(self.ds)
        for result in results:
            self.assert_result_is_good(result)

    def test_bad_geophysical_vars_fill_value(self):
        '''
        Tests that the geophysical variable _FillValue check is working
        '''

        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        nc_obj = Dataset(os.devnull, 'w', diskless=True)
        self.addCleanup(nc_obj.close)

        # The dataset needs at least one variable to check that it's missing
        # all the required attributes.
        nc_obj.createDimension('time', 1)
        nc_obj.createVariable('sample_var', 'd', ('time',))
        # Define some variable attributes but don't specify _FillValue
        sample_var = nc_obj.variables['sample_var']
        sample_var.units = 'm'
        sample_var.short_name = 'temp'
        # global attribute 'platform' points to variable that does not exist in dataset

        results = self.ioos.check_geophysical_vars_fill_value(nc_obj)
        for result in results:
            self.assert_result_is_bad(result)

    def test_good_geophysical_vars_fill_value(self):
        '''
        Tests that the geophysical variable _FillValue check is working
        '''
        results = self.ioos.check_geophysical_vars_fill_value(self.ds)
        for result in results:
            self.assert_result_is_good(result)

    def test_bad_geophysical_vars_standard_name(self):
        '''
        Tests that the platform variable attributes check is working
        '''

        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        nc_obj = Dataset(os.devnull, 'w', diskless=True)
        self.addCleanup(nc_obj.close)

        # The dataset needs at least one variable to check that it's missing
        # all the required attributes.
        nc_obj.createDimension('time', 1)
        nc_obj.createVariable('sample_var', 'd', ('time',))
        # Define some variable attributes but don't specify _FillValue
        sample_var = nc_obj.variables['sample_var']
        sample_var.units = 'm'
        sample_var.short_name = 'temp'
        # global attribute 'platform' points to variable that does not exist in dataset

        results = self.ioos.check_geophysical_vars_standard_name(nc_obj)
        for result in results:
            self.assert_result_is_bad(result)

    def test_good_geophysical_vars_standard_name(self):
        '''
        Tests that the geophysical variable _FillValue check is working
        '''
        results = self.ioos.check_geophysical_vars_standard_name(self.ds)
        for result in results:
            self.assert_result_is_good(result)

    def test_bad_units(self):
        '''
        Tests that the valid units check is working
        '''

        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        nc_obj = Dataset(os.devnull, 'w', diskless=True)
        self.addCleanup(nc_obj.close)

        # The dataset needs at least one variable to check that it's missing
        # all the required attributes.
        nc_obj.createDimension('time', 1)
        nc_obj.createVariable('temperature', 'd', ('time',))
        # Define some variable attributes but don't specify _FillValue
        sample_var = nc_obj.variables['temperature']
        sample_var.units = 'degC'   # Not valid units
        sample_var.short_name = 'temp'
        # global attribute 'platform' points to variable that does not exist in dataset

        results = self.ioos.check_geophysical_vars_standard_name(nc_obj)
        for result in results:
            self.assert_result_is_bad(result)

    def test_good_units(self):
        '''
        Tests that the valid units check is working
        '''
        results = self.ioos.check_units(self.ds)
        for result in results:
            self.assert_result_is_good(result)

class TestIOOS1_2(BaseTestCase):
    '''
    Tests for the compliance checker implementation of IOOS Metadata Profile
    for NetCDF, Version 1.1
    '''

    # for reference
        #ds.createDimension('siglev', 20)

        #temp = ds.createVariable("temp", np.float64, dimensions=("time",),
        #                         fill_value=np.float(99999999999999999999.))
        #temp.coordinates = "sigma noexist"
        #ds.createVariable("sigma", np.float64, dimensions=('siglev',))

    def setUp(self):
        self.ioos = IOOS1_2Check() 

    def test_check_vars_have_attrs(self):

        # create geophysical variable
        ds = MockTimeSeries() # time, lat, lon, depth 
        temp = ds.createVariable("temp", np.float64, dimensions=("time",))

        # should fail here
        results = self.ioos.check_vars_have_attrs(ds)
        scored, out_of, messages = get_results(results)
        self.assertLess(scored, out_of)

        # should pass
        ds = MockTimeSeries() # time, lat, lon, depth 
        temp = ds.createVariable("temp", np.float64, fill_value=9999999999.) # _FillValue
        temp.setncattr("missing_value", 9999999999.)
        temp.setncattr("standard_name", "sea_surface_temperature")
        temp.setncattr("standard_name_uri", "http://cfconventions.org/Data/cf-standard-names/64/build/cf-standard-name-table.html")
        temp.setncattr("units", "degree_C")
        temp.setncattr("platform", "myPlatform")

        results = self.ioos.check_vars_have_attrs(ds)
        scored, out_of, messages = get_results(results)
        self.assertEqual(scored, out_of)

    def test_check_single_platform(self):

        ds = MockTimeSeries()

        # should fail as no platform attrs exist
        result = self.ioos.check_single_platform(ds)
        self.assertFalse(result.value)

        # should pass with single platform
        ds.setncattr("platform", "myPlatform")
        result = self.ioos.check_single_platform(ds)
        self.assertTrue(result.value)

        # add a different platform
        ds.variables["time"].setncattr("platform", "myPlatform2")
        result = self.ioos.check_single_platform(ds)
        self.assertFalse(result.value)

        # make all the same
        ds.variables["time"].setncattr("platform", "myPlatform")
        result = self.ioos.check_single_platform(ds)
        self.assertTrue(result.value)

    def test_check_gts_ingest(self):
        ds = MockTimeSeries() # time, lat, lon, depth 

        # no gts_ingest, should pass
        results = self.ioos.check_gts_ingest(ds)
        scored, out_of, messages = get_results(results)
        self.assertEqual(scored, out_of)

        # global
        ds.setncattr("gts_ingest", "true")
        results = self.ioos.check_gts_ingest(ds)
        scored, out_of, messages = get_results(results)
        self.assertEqual(scored, out_of)

        # give one variable the gts_ingest attribute
        ds.variables["time"].setncattr("gts_ingest", "true")
        results = self.ioos.check_gts_ingest(ds)
        scored, out_of, messages = get_results(results)
        self.assertEqual(scored, out_of)

        # give a poor value
        ds.variables["time"].setncattr("gts_ingest", "blah")
        results = self.ioos.check_gts_ingest(ds)
        scored, out_of, messages = get_results(results)
        self.assertLess(scored, out_of)

    def test_check_instrument_variables(self):

        ds = MockTimeSeries() # time, lat, lon, depth 

        # no instrument variable, should pass
        results = self.ioos.check_instrument_variables(ds)
        scored, out_of, messages = get_results(results)
        self.assertEqual(scored, out_of)

        temp = ds.createVariable("temp", np.float64, dimensions=("time",))
        temp.setncattr("cf_role", "timeseries")
        temp.setncattr("standard_name", "sea_surface_temperature")
        temp.setncattr("units", "degree_C")
        temp.setncattr("axis", "Y")
        temp.setncattr("instrument", "myInstrument")
        temp[:] = 45.
        instr = ds.createVariable("myInstrument", np.float64, dimensions=("time",))

        # give instrument variable with component
        instr.setncattr("component", "someComponent")
        results = self.ioos.check_instrument_variables(ds)
        scored, out_of, messages = get_results(results)
        self.assertEqual(scored, out_of)

        # give discriminant
        instr.setncattr("discriminant", "someDiscriminant")
        results = self.ioos.check_instrument_variables(ds)
        scored, out_of, messages = get_results(results)
        self.assertEqual(scored, out_of)

        # bad component
        instr.setncattr("component", 45)
        results = self.ioos.check_instrument_variables(ds)
        scored, out_of, messages = get_results(results)
        self.assertLess(scored, out_of)
        
