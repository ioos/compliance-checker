from compliance_checker.ioos import (IOOS0_1Check, IOOS1_1Check, IOOS1_2Check,
                                     NamingAuthorityValidator)
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
        ds = MockTimeSeries(default_fill_value=9999999999.) # time, lat, lon, depth 
        temp = ds.createVariable("temp", np.float64, fill_value=9999999999.) # _FillValue
        temp.setncattr("missing_value", 9999999999.)
        temp.setncattr("standard_name", "sea_surface_temperature")
        temp.setncattr("standard_name_uri", "http://cfconventions.org/Data/cf-standard-names/64/build/cf-standard-name-table.html")
        temp.setncattr("units", "degree_C")
        temp.setncattr("platform", "myPlatform")

        ds.variables["time"].setncattr("platform", "myPlatform")
        ds.variables["time"].setncattr("standard_name", "time")
        ds.variables["time"].setncattr("standard_name_uri", "time")
        ds.variables["time"].setncattr("units", "hours since 1970-01-01T00:00:00")
        ds.variables["time"].setncattr("missing_value", 9999999999.)

        results = self.ioos.check_vars_have_attrs(ds)
        scored, out_of, messages = get_results(results)
        self.assertEqual(scored, out_of)

    def test_check_contributor_role_and_vocabulary(self):
        ds = MockTimeSeries() # time, lat, lon, depth 

        # no contributor_role or vocab, fail both
        results = self.ioos.check_contributor_role_and_vocabulary(ds)
        self.assertFalse(all(r.value for r in results))

        # bad contributor_role and vocab
        ds.setncattr("contributor_role", "bad")
        ds.setncattr("contributor_role_vocabulary", "bad")
        results = self.ioos.check_contributor_role_and_vocabulary(ds)
        self.assertFalse(all(r.value for r in results))

        # good role, bad vocab
        ds.setncattr("contributor_role", "contributor")
        results = self.ioos.check_contributor_role_and_vocabulary(ds)
        self.assertTrue(results[0].value)
        self.assertFalse(results[1].value)

        # bad role, good vocab
        ds.setncattr("contributor_role", "bad")
        ds.setncattr("contributor_role_vocabulary", "http://vocab.nerc.ac.uk/collection/G04/current/")
        results = self.ioos.check_contributor_role_and_vocabulary(ds)
        self.assertFalse(results[0].value)
        self.assertTrue(results[1].value)
        
        # good role, good vocab
        ds.setncattr("contributor_role", "contributor")
        ds.setncattr("contributor_role_vocabulary", "http://vocab.nerc.ac.uk/collection/G04/current/")
        results = self.ioos.check_contributor_role_and_vocabulary(ds)
        self.assertTrue(results[0].value)
        self.assertTrue(results[1].value)

    def test_check_creator_and_publisher_type(self):
        """
        Checks the creator_type and publisher_type global attributes with
        the following values:
        Empty: Valid, defaults to "person" when not specified, which is
               contained in the list of valid values.
        Bad values: Invalid, not contained in list of valid values.
        Good values: Valid, contained in list.
        """
        ds = MockTimeSeries()
        # values which are not set/specified default to person, which is valid
        result_list = self.ioos.check_creator_and_publisher_type(ds)
        self.assertTrue(all(res.value for res in result_list))
        # create invalid values for attribute
        ds.setncattr('creator_type', 'PI')
        ds.setncattr('publisher_type', 'Funder')
        result_list = self.ioos.check_creator_and_publisher_type(ds)
        err_regex = (r"^If specified, \w+_type must be in value list "
                     r"\(\['group', 'institution', 'person', 'position'\]\)$")
        for res in result_list:
            self.assertFalse(res.value)
            self.assertRegex(res.msgs[0], err_regex)
        # good values
        ds.setncattr('creator_type', 'person')
        ds.setncattr('publisher_type', 'institution')
        result_list = self.ioos.check_creator_and_publisher_type(ds)
        self.assertTrue(all(res.value for res in result_list))

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
        # no ancillary vars, should fail
        ds.variables["time"].setncattr("gts_ingest", "true")
        results = self.ioos.check_gts_ingest(ds)
        scored, out_of, messages = get_results(results)
        self.assertLess(scored, out_of)

        # set ancillary var with bad standard name
        tmp = ds.createVariable("tmp", np.byte, ("time",))
        tmp.setncattr("standard_name", "bad")
        ds.variables["time"].setncattr("ancillary_variables", "tmp")
        results = self.ioos.check_gts_ingest(ds)
        scored, out_of, messages = get_results(results)
        self.assertLess(scored, out_of)

         # good ancillary var standard name
        tmp.setncattr("standard_name", "aggregate_quality_flag")
        results = self.ioos.check_gts_ingest(ds)
        scored, out_of, messages = get_results(results)
        self.assertEqual(scored, out_of)

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

    def test_check_wmo_platform_code(self):
        ds = MockTimeSeries() # time, lat, lon, depth

        # no wmo_platform_code, pass
        result = self.ioos.check_wmo_platform_code(ds)
        self.assertTrue(result.value)

        # valid code
        ds.setncattr("wmo_platform_code", "12345")
        result = self.ioos.check_wmo_platform_code(ds)
        self.assertTrue(result.value)

        # valid code
        ds.setncattr("wmo_platform_code", "7654321")
        result = self.ioos.check_wmo_platform_code(ds)
        self.assertTrue(result.value)

        # non-numeric, fail
        ds.setncattr("wmo_platform_code", "abcd1")
        result = self.ioos.check_wmo_platform_code(ds)
        self.assertFalse(result.value)

        # invalid length, fail
        ds.setncattr("wmo_platform_code", "123")
        result = self.ioos.check_wmo_platform_code(ds)
        self.assertFalse(result.value)

    def test_check_standard_name(self):
        ds = MockTimeSeries() # time, lat, lon, depth

        # no standard names
        results = self.ioos.check_standard_name(ds)
        scored, out_of, messages = get_results(results)
        self.assertLess(scored, out_of)

        # give standard names to all variables
        ds.variables["time"].setncattr("standard_name", "time")
        ds.variables["lon"].setncattr("standard_name", "longitude")
        ds.variables["lat"].setncattr("standard_name", "latitude")
        ds.variables["depth"].setncattr("standard_name", "depth")
        results = self.ioos.check_standard_name(ds)
        scored, out_of, messages = get_results(results)
        self.assertEqual(scored, out_of)

        # add a QARTOD variable, no standard name - should fail
        qr = ds.createVariable("depth_qc", np.byte)
        qr.setncattr("flag_meanings", "blah")
        results = self.ioos.check_standard_name(ds)
        scored, out_of, messages = get_results(results)
        self.assertLess(scored, out_of)

        # bad standard name
        qr.setncattr("standard_name", "blah")
        results = self.ioos.check_standard_name(ds)
        scored, out_of, messages = get_results(results)
        self.assertLess(scored, out_of)

        # good standard name
        qr.setncattr("standard_name", "spike_test_quality_flag")
        results = self.ioos.check_standard_name(ds)
        scored, out_of, messages = get_results(results)
        self.assertEqual(scored, out_of)

    def test_naming_authority_validation(self):
        test_attr_name = "naming_authority"
        validator = NamingAuthorityValidator()
        # check URL - should pass
        self.assertTrue(validator.validate(test_attr_name,
                                           "https://ioos.us")[0])
        # check reverse DNS - should pass
        self.assertTrue(validator.validate(test_attr_name,
                                           "edu.ucar.unidata")[0])
        # email address is neither of the above, so should fail
        bad_result = validator.validate(test_attr_name,
                                        "webmaster.ioos.us@noaa.gov")
        self.assertFalse(bad_result[0])
        self.assertEqual(bad_result[1],
                         "naming_authority should either be a URL or a "
                         "reversed DNS name (e.g \"edu.ucar.unidata\")")

    def test_check_platform_cf_role(self):
        """
        Check that cf_role inside platform variables only allows certain
        values, namely "profile_id", "timeseries_id", or "trajectory_id"
        """
        ds = MockTimeSeries()
        plat_var = ds.createVariable("platform", np.int8, ())
        ds.variables['depth'].platform = "platform"
        self.ioos.setup(ds)
        results = self.ioos.check_platform_variable_cf_role(ds)
        # don't set attribute, should raise error about attribute not
        # existing
        self.assertEqual(len(results), 1)
        score, out_of = results[0].value
        self.assertLess(score, out_of)
        # set to invalid value
        plat_var.setncattr("cf_role", "bad_value")
        results = self.ioos.check_platform_variable_cf_role(ds)
        self.assertLess(score, out_of)
        expected_vals = ["profile_id", "timeseries_id", "trajectory_id"]
        expect_msg = ("attribute cf_role in variable platform present, but not "
                      "in expected value list ({})".format(expected_vals))
        self.assertEqual(results[0].msgs, [expect_msg])
        # set to valid value
        plat_var.setncattr("cf_role", "timeseries_id")
        results = self.ioos.check_platform_variable_cf_role(ds)
        score, out_of = results[0].value
        self.assertEqual(score, out_of)

    def test_check_single_platform(self):

        ds = MockTimeSeries() # time, lat, lon, depth 

        # no global attr but also no platform variables, should pass
        results = self.ioos.check_single_platform(ds)
        self.assertTrue(results[0].value)

        # give platform global, no variables, fail
        ds.setncattr("platform", "buoy")
        results = self.ioos.check_single_platform(ds)
        self.assertFalse(results[0].value)

        # global attribute, one platform variable, correct cf_role & featureType, pass
        ds.setncattr("featureType", "profile")
        ds.createDimension("profile", 1)
        temp = ds.createVariable("temp", "d", ("time"))
        temp.setncattr("platform", "platform_var")
        plat = ds.createVariable("platform_var", np.byte)
        cf_role_var = ds.createVariable("cf_role_var", np.byte, ("profile",))
        cf_role_var.setncattr("cf_role", "timeseries_id")
        results = self.ioos.check_single_platform(ds)
        self.assertTrue(all(r.value for r in results))

        # global attr, multiple platform variables, correct cf_role & featureType, fail
        plat2 = ds.createVariable("platform_var_2", np.byte)
        temp2 = ds.createVariable("temp2", "d", ("time"))
        temp2.setncattr("platform", "platform_var2")
        results = self.ioos.check_single_platform(ds)
        self.assertFalse(results[0].value)

        # no global attr, one platform var, correct cf_role & featureType, fail
        ds.delncattr("platform")
        self.assertFalse(results[0].value)

        # global attr, one platform var, correct featureType, incorrect cf_role var dimension
        ds = MockTimeSeries() # time, lat, lon, depth 
        ds.setncattr("featureType", "trajectoryprofile")
        ds.createDimension("trajectory", 2) # should only be 1
        temp = ds.createVariable("temp", "d", ("time"))
        temp.setncattr("platform", "platform_var")
        plat = ds.createVariable("platform_var", np.byte)
        cf_role_var = ds.createVariable("cf_role_var", np.byte, ("trajectory",))
        cf_role_var.setncattr("cf_role", "trajectory_id")
        results = self.ioos.check_single_platform(ds)
        self.assertFalse(results[0].value)

    def test_check_platform_vocabulary(self):
        ds = MockTimeSeries() # time, lat, lon, depth 
        ds.setncattr("platform_vocabulary", "http://google.com")
        self.assertTrue(self.ioos.check_platform_vocabulary(ds).value)

        ds.setncattr("platform_vocabulary", "bad")
        self.assertFalse(self.ioos.check_platform_vocabulary(ds).value)

    def test_check_qartod_variables_references(self):
        ds = MockTimeSeries() # time, lat, lon, depth

        # no QARTOD variables
        results = self.ioos.check_qartod_variables_references(ds)
        scored, out_of, messages = get_results(results)
        self.assertEqual(scored, out_of)

        # QARTOD variable without references (fail)
        qr = ds.createVariable("depth_qc", np.byte)
        qr.setncattr("flag_meanings", "blah")
        qr.setncattr("standard_name", "spike_test_quality_flag")
        results = self.ioos.check_qartod_variables_references(ds)
        self.assertFalse(all(r.value for r in results))

        # QARTOD variable with references (pass)
        qr.setncattr("references", "http://services.cormp.org/quality.php")
        results = self.ioos.check_qartod_variables_references(ds)
        self.assertTrue(all(r.value for r in results))

        # QARTOD variable with bad references (fail)
        qr.setncattr("references", r"p9q384ht09q38@@####???????////??//\/\/\/\//\/\74ht")
        results = self.ioos.check_qartod_variables_references(ds)
        self.assertFalse(all(r.value for r in results))
