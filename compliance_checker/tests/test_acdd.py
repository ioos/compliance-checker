import unittest
from compliance_checker.acdd import ACDDBaseCheck
# from netCDF4 import Dataset
import os
import pandas as pd

# not updated

def to_singleton_var(l):
    """Get the first value of a list if this is a list"""
    return [x[0] if hasattr(x, '__iter__') else x for x in l]

def check_varset_nonintersect(true_set, name_list):
    """
    Helper function intended to check that the series we scraped from the HTML
    table is equal to the list of values we have
    """
    return len(true_set ^ set(to_singleton_var(name_list))) == 0



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
                                        'acknowledgement', 'geospatial_bounds',
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
        #ds = Dataset("/Users/asadeveloper/Downloads/hycomglobalnavy_2012120300.nc")
        self.acdd = ACDDBaseCheck(version='1.1')

    # TODO: Break into multiple tests
    def test_high_rec_present(self):
        """Checks that all highly recommended attributes are present"""
        assert check_varset_nonintersect(self.expected['Highly Recommended'],
                                        self.acdd.high_rec_atts)

    def test_rec_present(self):
        """Checks that all recommended attributes are present"""
        # 'geospatial_bounds' attribute currently has its own separate check
        # from the list of required atts
        rec_atts = self.acdd.rec_atts + ['geospatial_bounds']
        assert check_varset_nonintersect(self.expected['Recommended'],
                                         rec_atts)

    def test_sug_present(self):
        """Checks that all suggested attributes are present"""
        assert check_varset_nonintersect(self.expected['Suggested'],
                                         self.acdd.sug_atts)

    #assert self.acdd.check_high(ds) is True


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
                            'acknowledgement', 'license',
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
        self.acdd = ACDDBaseCheck(version='1.3')

    def test_high_rec_present(self):
        """Checks that all highly recommended attributes are present"""
        assert check_varset_nonintersect(self.expected['Highly Recommended'],
                                        self.acdd.high_rec_atts)

    def test_rec_present(self):
        """Checks that all recommended attributes are present"""
        # 'geospatial_bounds' attribute currently has its own separate check
        # from the list of required atts
        rec_atts = self.acdd.rec_atts + ['geospatial_bounds']
        assert check_varset_nonintersect(self.expected['Recommended'],
                                         rec_atts)

    def test_sug_present(self):
        """Checks that all suggested attributes are present"""
        assert check_varset_nonintersect(self.expected['Suggested'],
                                         self.acdd.sug_atts)

    #assert self.acdd.check_high(ds) is True
