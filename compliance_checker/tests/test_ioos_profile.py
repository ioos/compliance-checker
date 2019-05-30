from compliance_checker.ioos import IOOS0_1Check, IOOS1_1Check, IOOS1_2Check
from compliance_checker.tests.resources import STATIC_FILES
from compliance_checker.tests import BaseTestCase
from compliance_checker.tests.helpers import MockTimeSeries
import numpy as np
from netCDF4 import Dataset
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
        sample_var.units = 'invalidUnits'   # Not valid units
        sample_var.short_name = 'temp'
        # global attribute 'platform' points to variable that does not exist in dataset

        results = self.ioos.check_units(nc_obj)
        self.assert_result_is_good(results[0]) # should be good
        self.assert_result_is_bad(results[1])  # 'temp' variable, should fail
        self.assert_result_is_good(results[2])  # should be good

    def test_good_units(self):
        '''
        Tests that the valid units check is working
        '''
        results = self.ioos.check_units(self.ds)
        for result in results:
            self.assert_result_is_good(result)

    def test_check_instrument_vars(self):
        """Test the intrument_variable check works appropriately"""

        ds = MockTimeSeries()

        # no instrument variable, yet, should fail
        results = self.ioos.check_instrument_variables(ds)
        assert len(results) == 0 # no results given for skipped variable

        # give a fake geophysical variable
        temp = ds.createVariable('temp', 'float32', ('time',))
        temp[:] = np.linspace(0., 30., num=500)

        # give it an instrument attr but no discriminant, should still fail
        ds.variables["temp"].setncattr("instrument", "instrument")
        instrument = ds.createVariable('instrument', 'c')
        results = self.ioos.check_instrument_variables(ds)
        self.assert_result_is_bad(results[0])

        # give it discriminant attr
        instrument.setncattr("discriminant", "instrumentDiscriminant")
        results = self.ioos.check_instrument_variables(ds)
        self.assert_result_is_good(results[0])


class TestIOOS1_2(BaseTestCase):
    '''
    Tests for the compliance checker implementation of IOOS Metadata Profile
    for NetCDF, Version 1.1
    '''
    def setUp(self):
        self.ds = MockTimeSeries()
        # give a fake geophysical variable
        temp = self.ds.createVariable('temp', 'float32', ('time',))
        temp[:] = np.linspace(0., 30., num=500)
        self.ioos = IOOS1_2Check()

    def test_cc_meta(self):
        assert self.ioos._cc_spec == 'ioos'
        assert self.ioos._cc_spec_version == '1.2'

    def test_required_attributes(self):
        """Check the dataset has all required attrs""" 
        # set the required attrs to the in-memory ds
        self.ds.setncattr('creator_country', "USA")
        self.ds.setncattr('creator_email', "noreply@oceansmap.com")
        self.ds.setncattr('creator_institution', "RPS Ocean Science")
        self.ds.setncattr('creator_sector', "industry")
        self.ds.setncattr('creator_url', "www.rpsgroup.com")
        self.ds.setncattr('featureType', "timeSeries")
        self.ds.setncattr('id', "myID")
        self.ds.setncattr('info_url', 'www.rpsgroup.com')
        self.ds.setncattr('naming_authority', "IOOS")
        self.ds.setncattr('platform_name', "Dr. Floats-a-Lot")
        self.ds.setncattr('platform_vocabulary', "https://mmisw.org/ont/ioos/platform")
        self.ds.setncattr('publisher_country', "USA")
        self.ds.setncattr('publisher_email', "noreply@oceansmap.com")
        self.ds.setncattr('publisher_name', "RPS Ocean Science")
        self.ds.setncattr('publisher_url', "rpsgroup.com")
        self.ds.setncattr('standard_name_vocabulary', "NetCDF Climate and Forecast (CF) Metadata Convention Standard Name Table v26")
        self.ds.setncattr('title', 'my Dataset Title')

        results = self.ioos.check_high(self.ds)
        for result in results:
            self.assert_result_is_good(result)

    def test_recomended_attributes(self):
        '''
        Tests that recommended attributes test is working properly
        '''

        self.ds.setncattr('contributor_email', "noreply@oceansmap.com")
        self.ds.setncattr('contributor_name', "RPS Ocean Science")
        self.ds.setncattr('contributor_role', "DMAC")
        self.ds.setncattr('contributor_role_vocabulary', "http://vocab.nerc.ac.uk/collection/G04/current/")
        self.ds.setncattr('contributor_url', "www.rspgroup.com")
        self.ds.setncattr('creator_address', "test rd")
        self.ds.setncattr('creator_city', "test")
        self.ds.setncattr('creator_name', "RPS Ocean Science")
        self.ds.setncattr('creator_phone', "1-800-pls-call")
        self.ds.setncattr('creator_state', "Lil Rhody")
        self.ds.setncattr('creator_postalcode',"38890")
        self.ds.setncattr('institution', "RPS")
        self.ds.setncattr("instrument", "instrument")
        self.ds.setncattr('keywords', "test")
        self.ds.setncattr('license', "test")
        self.ds.setncattr('platform_id', "48-b25")
        self.ds.setncattr('publisher_address', "123 Yellow Brick Rd")
        self.ds.setncattr('publisher_city', "Emerald City")
        self.ds.setncattr('publisher_phone', "1-800-pls-call")
        self.ds.setncattr('publisher_state', "Lil Rhody")
        self.ds.setncattr('publisher_postalcode', "02879")
        self.ds.setncattr('summary', "Nah.")

        results = self.ioos.check_recommended(self.ds)
        for result in results:
            self.assert_result_is_good(result)
        
    def test_check_instrument_vars(self):
        """Test the check_instrument_vars behaves as expected"""

        # no instrument variable, yet, should fail
        results = self.ioos.check_instrument_variables(self.ds)
        assert len(results) == 0 # no results given for skipped variable

        temp = self.ds.variables["temp"]

        # give it an instrument attr but no discriminant, should still fail
        temp.setncattr("instrument", "instrument")
        instrument = self.ds.createVariable('instrument', 'c')
        results = self.ioos.check_instrument_variables(self.ds)
        self.assert_result_is_bad(results[0])

        # give it component attr, one should be good, other should fail
        instrument.setncattr("component", "instrumentDiscriminant")
        results = self.ioos.check_instrument_variables(self.ds)
        self.assert_result_is_good(results[0])
        self.assert_result_is_bad(results[1])

        # give it discriminant attr, one should be good, other should fail
        instrument.setncattr("discriminant", "instrumentDiscriminant")
        results = self.ioos.check_instrument_variables(self.ds)
        self.assert_result_is_good(results[0])
        self.assert_result_is_good(results[1])

    def test_check_geophysical_vars_standard_name(self):
        """Test checking the geophysical vars behaves as expected"""

        temp = self.ds.variables["temp"]

        # results will be a list of length 2, since we only have 1 variable
        results = self.ioos.check_geophysical_vars_standard_name(self.ds)
        self.assert_result_is_bad(results[0])
        self.assert_result_is_bad(results[1])

        # give the temp variable a standard name, should still fail
        temp.setncattr("standard_name", "sea_water_temperature")
        results = self.ioos.check_geophysical_vars_standard_name(self.ds)
        self.assert_result_is_good(results[0])
        self.assert_result_is_bad(results[1]) # only first should be good

        # finally, give standard_name_uri
        temp.setncattr("standard_name_uri", "http://vocab.nerc.ac.uk/collection/P07/current/CFSN0335/")
        results = self.ioos.check_geophysical_vars_standard_name(self.ds)
        self.assert_result_is_good(results[0])
        self.assert_result_is_good(results[1])

    def test_check_units(self):
        '''
        Tests that the valid units check is working
        '''
        # first run, should have fails as 'temp' has no units
        results = self.ioos.check_units(self.ds)
        assert(not all([r.value[0] for r in results]))

        # add units to 'temp', all should pass
        self.ds.variables["temp"].units = "degC"
        results = self.ioos.check_units(self.ds)
        for result in results:
            self.assert_result_is_good(result)

    def test_check_platform_variable(self):
        "Test the check_platform_variable check is behaving as expected"""

        # no platform variable, yet, should fail
        result = self.ioos.check_platform_variable(self.ds)
        self.assert_result_is_bad(result)

        # give it a platform variable but no cf_role, should fail
        self.ds.variables["temp"].setncattr("platform", "platform_var")
        self.ds.setncattr("platform", "platform_var")
        platform = self.ds.createVariable('platform_var', 'c')
        result = self.ioos.check_platform_variable(self.ds)
        self.assert_result_is_bad(result)

        # give a cf_role
        platform.setncattr("cf_role", "timeSeries")
        result = self.ioos.check_platform_variable(self.ds)
        self.assert_result_is_good(result)

    def test_check_wmo_platform_id(self):
        """Test the check_wmo_platform_id functions as expected"""

        # no wmo_platform_id, should be fine
        result = self.ioos.check_wmo_platform_id(self.ds)
        self.assert_result_is_good(result)

        # compliant IDs
        self.ds.setncattr("wmo_platform_id", "77777")
        result = self.ioos.check_wmo_platform_id(self.ds)
        self.assert_result_is_good(result)
        
        self.ds.setncattr("wmo_platform_id", "5555555")
        result = self.ioos.check_wmo_platform_id(self.ds)
        self.assert_result_is_good(result)

        self.ds.setncattr("wmo_platform_id", "aaaabcd")
        result = self.ioos.check_wmo_platform_id(self.ds)
        self.assert_result_is_good(result)

        self.ds.setncattr("wmo_platform_id", "712a4hi")
        result = self.ioos.check_wmo_platform_id(self.ds)
        self.assert_result_is_good(result)

        self.ds.setncattr("wmo_platform_id", "2a4hi")
        result = self.ioos.check_wmo_platform_id(self.ds)
        self.assert_result_is_good(result)

        self.ds.setncattr("wmo_platform_id", "zzzzz")
        result = self.ioos.check_wmo_platform_id(self.ds)
        self.assert_result_is_good(result)

        # non-compliant IDs
        self.ds.setncattr("wmo_platform_id", "a4h")
        result = self.ioos.check_wmo_platform_id(self.ds)
        self.assert_result_is_bad(result)

        self.ds.setncattr("wmo_platform_id", "5555")
        result = self.ioos.check_wmo_platform_id(self.ds)
        self.assert_result_is_bad(result)

        self.ds.setncattr("wmo_platform_id", "qwerty")
        result = self.ioos.check_wmo_platform_id(self.ds)
        self.assert_result_is_bad(result)

        # non-compliant, integer and float codes
        self.ds.setncattr("wmo_platform_id", 77777)
        result = self.ioos.check_wmo_platform_id(self.ds)
        self.assert_result_is_bad(result)
        
        self.ds.setncattr("wmo_platform_id", 5.8091)
        result = self.ioos.check_wmo_platform_id(self.ds)
        self.assert_result_is_bad(result)
