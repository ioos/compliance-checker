import unittest
from compliance_checker.acdd import ACDD1_1Check, ACDD1_3Check, ACDDBaseCheck
from netCDF4 import Dataset
import os

def to_singleton_var(l):
    """
    Get the first value of a list if this implements iterator protocol
    and is not a string
    """
    return [x[0] if hasattr(x, '__iter__') and not isinstance(x, str) else x
            for x in l]

def check_varset_nonintersect(true_set, name_list):
    """
    Helper function intended to check that the series we scraped from the HTML
    table is equal to the list of values we have
    """
    return len(true_set ^ set(to_singleton_var(name_list))) == 0


# TODO: move common atts to Base ACDD check test

class TestACDDCommon(unittest.TestCase):
    def setUp(self):
        self.acdd = ACDDBaseCheck()
        self.ds = Dataset(filename=os.devnull, mode='w', diskless=True)

    def tearDown(self):
        self.ds.close()

    def test_verify_geospatial_bounds(self):
        """Tests the geospatial_bounds function"""
        self.ds.geospatial_bounds = 'POLYGON ((40.26 -111.29, 41.26 -111.29, 41.26 -110.29, 40.26 -110.29, 40.26 -111.29))'
        # give arbitrary rating and check that things passed
        result = self.acdd.verify_geospatial_bounds(self.ds)(1)
        self.assertTrue(result.value)

    def test_verify_valid_title(self):
        self.ds.title = '@@@@@@@invalid@@@@@'

    def check_valid_date(self):
        self.ds.date_created = '2011-04-27T00:00:00Z'
        self.ds.date_created = '2011-04-27T00:00:00Z'


class TestACDD1_1(unittest.TestCase):

    # Adapted using `pandas.read_html` from URL
    # http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery_1-1
    expected = {'Suggested': {'contributor_name', 'contributor_role',
                                'publisher_name', 'publisher_url',
                                'publisher_email', 'date_modified',
                                'date_issued', 'geospatial_lat_units',
                                'geospatial_lat_resolution',
                    'geospatial_lon_units', 'geospatial_lon_resolution',
                    'geospatial_vertical_units',
                    'geospatial_vertical_resolution',
                    'geospatial_vertical_positive'},
                'Highly Recommended': {'title', 'summary', 'keywords'},
                    'Recommended': {'id', 'naming_authority',
                                    'keywords_vocabulary', 'cdm_data_type',
                                    'history', 'comment', 'date_created',
                                    'creator_name', 'creator_url',
                                    'creator_email', 'institution',
                                    'project', 'processing_level',
                                    # Tested separately
                                    #'acknowledgement',
                                    'geospatial_bounds',
                                    'geospatial_lat_min',
                                    'geospatial_lat_max',
                                    'geospatial_lon_min',
                                    'geospatial_lon_max',
                                    'geospatial_vertical_min',
                                    'geospatial_vertical_max',
                                    'time_coverage_start',
                                    'time_coverage_end',
                                    'time_coverage_duration',
                                    'time_coverage_resolution',
                                    'standard_name_vocabulary', 'license'},
                'Highly Recommended Variable Attributes': {'long_name',
                                                'standard_name', 'units',
                                                'coverage_content_type'}}

    def setUp(self):
        # TODO: Find or make a canonical ACDD 1.1 reference file
        self.acdd = ACDD1_1Check()

    def test_cc_meta(self):
        assert self.acdd._cc_spec == 'acdd'
        assert self.acdd._cc_spec_version == '1.1'

    def test_high_rec_present(self):
        """Checks that all highly recommended attributes are present"""
        assert check_varset_nonintersect(self.expected['Highly Recommended'],
                                        self.acdd.high_rec_atts)

    def test_rec_present(self):
        """Checks that all recommended attributes are present"""
        # 'geospatial_bounds' attribute currently has its own separate check
        # from the list of required atts
        assert check_varset_nonintersect(self.expected['Recommended'],
                                         self.acdd.rec_atts)

    def test_sug_present(self):
        """Checks that all suggested attributes are present"""
        assert check_varset_nonintersect(self.expected['Suggested'],
                                         self.acdd.sug_atts)


class TestACDD1_3(unittest.TestCase):
    # Adapted using `pandas.read_html` from URL
    # http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery_1-3
    expected = {'Suggested': {'creator_type', 'creator_institution',
                              'publisher_type', 'publisher_institution',
                              'program', 'contributor_name', 'contributor_role',
                              'geospatial_lat_units',
                              'geospatial_lat_resolution',
                              'geospatial_lon_units',
                              'geospatial_lon_resolution',
                              'geospatial_vertical_units',
                              'geospatial_vertical_resolution',
                              'date_modified', 'date_issued',
                              'date_metadata_modified', 'product_version',
                              'keywords_vocabulary', 'platform',
                              'platform_vocabulary', 'instrument',
                              'instrument_vocabulary', 'cdm_data_type',
                              'metadata_link', 'references'},
            'Highly Recommended': {'title', 'summary', 'keywords',
                                   'Conventions'},
            'Recommended': {'id', 'naming_authority', 'history', 'source',
                            'processing_level', 'comment',
                            # Tested separately
                            #'acknowledgement',
                            'license',
                            'standard_name_vocabulary', 'date_created',
                            'creator_name', 'creator_email', 'creator_url',
                            'institution', 'project', 'publisher_name',
                            'publisher_email', 'publisher_url',
                            'geospatial_bounds', 'geospatial_bounds_crs',
                            'geospatial_bounds_vertical_crs',
                            'geospatial_lat_min', 'geospatial_lat_max',
                            'geospatial_lon_min', 'geospatial_lon_max',
                            'geospatial_vertical_min',
                            'geospatial_vertical_max',
                            'geospatial_vertical_positive',
                            'time_coverage_start',
                            'time_coverage_end',
                            'time_coverage_duration',
                            'time_coverage_resolution'},
                'Highly Recommended Variable Attributes': {'long_name',
                                                           'standard_name',
                                                           'units',
                                                           'coverage_content_type'}
            }

    def setUp(self):
        # TODO: Find or make a canonical ACDD 1.3 reference file
        #ds = Dataset("/Users/asadeveloper/Downloads/hycomglobalnavy_2012120300.nc")
        # data originally obtained from
        # http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery_1-3
        self.acdd = ACDD1_3Check()
        # create in memory netCDF file for testing purposes
        self.ds = Dataset(filename=os.devnull, mode='w', diskless=True)

    def tearDown(self):
        self.ds.close()

    def test_cc_meta(self):
        assert self.acdd._cc_spec == 'acdd'
        assert self.acdd._cc_spec_version == '1.3'

    def test_high_rec_present(self):
        """Checks that all highly recommended attributes are present"""
        assert check_varset_nonintersect(self.expected['Highly Recommended'],
                                        self.acdd.high_rec_atts)

    def test_high_rec(self):
        self.acdd.check_high(self.ds)

    def test_rec_present(self):
        """Checks that all recommended attributes are present"""
        # 'geospatial_bounds' attribute currently has its own separate check
        # from the list of required atts
        assert check_varset_nonintersect(self.expected['Recommended'],
                                         self.acdd.rec_atts)

    def test_sug_present(self):
        """Checks that all suggested attributes are present"""
        assert check_varset_nonintersect(self.expected['Suggested'],
                                         self.acdd.sug_atts)

    def test_var_coverage_content_type(self):
        var = self.ds.createVariable('foo', 'i4')
        var.coverage_content_type = 'modelResult'
        result = self.acdd.check_var_coverage_content_type(self.ds)
        assert len(result) == 0
