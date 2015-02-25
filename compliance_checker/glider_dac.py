#!/usr/bin/env python
'''
compliance_checker.glider_dac

Compliance Test Suite for the IOOS National Glider Data Assembly Center
https://github.com/ioos/ioosngdac/wiki
'''

from compliance_checker.base import BaseCheck, BaseNCCheck, Result

class GliderCheck(BaseNCCheck):
    @classmethod
    def beliefs(cls): 
        '''
        Not applicable for gliders
        '''
        return {}

    @classmethod
    def make_result(cls, level, score, out_of, name, messages):
        return Result(level, (score, out_of), name, messages)

    def setup(self, ds):
        pass

    def check_locations(self, ds):
        '''
        Validates that lat and lon are indeed variables
        '''
        level = BaseCheck.HIGH
        out_of = 1
        score = 0
        messages = []
        test = ('lat' in ds.dataset.variables and 'lon' in ds.dataset.variables)
        if test:
            score += 1
        else:
            messages.append("lat and lon don't exist")
        return self.make_result(level, score, out_of, 'Verifies lat and lon are variables', messages)

    def check_location_dimensions(self, ds):
        '''
        Validates that lat and lon are valid timeseries variables
        '''
        level = BaseCheck.MEDIUM
        out_of = 1
        score = 0
        messages = []
        if ('lat' in ds.dataset.variables and 'lon' in ds.dataset.variables):
            test = ds.dataset.variables['lat'].dimensions == ('time',)
            test &= ds.dataset.variables['lon'].dimensions == ('time',)
            score = int(test)
            if not test:
                messages.append('Latitude and Longitude do not use time as their only dimension') 
        else:
            messages.append('Latitude and Longitude were not found in the dataset')
        return self.make_result(level, score, out_of, 'Lat and Lon are Time Series', messages)

    def check_variables(self, ds):
        '''
        Verifies the dataset has the required variables
        '''
        required_variables = [
            'trajectory',
            'time',
            'lat',
            'lon',
            'pressure',
            'depth',
            'temperature',
            'conductivity',
            'density',
            'profile_id',
            'profile_time',
            'profile_lat',
            'profile_lon',
            'time_uv',
            'lat_uv',
            'lon_uv',
            'u',
            'v',
            'platform',
            'instrument_ctd'
        ]
        level = BaseCheck.HIGH
        out_of = len(required_variables)
        score = 0
        messages = []
        for variable in required_variables:
            test = variable in ds.dataset.variables
            score += int(test)
            if not test:
                messages.append("%s is a required variable" % variable)
        return self.make_result(level, score, out_of, 'Required Variables', messages)

    def check_qc_variables(self, ds):
        '''
        Verifies the dataset has all the required QC variables
        '''
        required_variables = [
            'time_qc',
            'lat_qc',
            'lon_qc',
            'pressure_qc',
            'depth_qc',
            'temperature_qc',
            'conductivity_qc',
            'density_qc',
            'profile_time_qc',
            'profile_lat_qc',
            'profile_lon_qc',
            'time_uv_qc',
            'lat_uv_qc',
            'lon_uv_qc',
            'u_qc',
            'v_qc'
        ]
        level = BaseCheck.MEDIUM
        out_of = len(required_variables)
        score = 0
        messages = []
        for variable in required_variables:
            test = variable in ds.dataset.variables
            score += int(test)
            if not test:
                messages.append("%s is a required qc variable" % variable)
        return self.make_result(level, score, out_of, 'Required QC Variables', messages)

    def check_global_attributes(self, ds):
        '''
        Verifies the base metadata in the global attributes
        '''
        attribute_fields = [
            'Conventions',
            'Metadata_Conventions',
            'comment',
            'contributor_name',
            'contributor_role',
            'creator_email',
            'creator_name',
            'creator_url',
            'date_created',
            'date_issued',
            'date_modified',
            'format_version',
            'history',
            'id',
            'institution',
            'keywords',
            'keywords_vocabulary',
            'license',
            'metadata_link',
            'naming_authority',
            'platform_type',
            'processing_level',
            'project',
            'publisher_email',
            'publisher_name',
            'publisher_url',
            'references',
            'sea_name',
            'source',
            'standard_name_vocabulary',
            'summary',
            'title',
            'wmo_id'
        ]
        level = BaseCheck.MEDIUM
        out_of = len(attribute_fields)
        score = 0
        messages = []
        for field in attribute_fields:
            test = hasattr(ds.dataset, field)
            score += int(test)
            if not test:
                messages.append('%s global attribute is missing' % field)
        
        return self.make_result(level, score, out_of, 'Required Global Attributes', messages)

    def check_primary_variable_attributes(self, ds):
        '''
        Verifies that each primary variable has the necessary metadata
        '''
        level = BaseCheck.MEDIUM
        out_of = 0
        score = 0
        messages = []
        primary_variables = [
            'time',
            'lat',
            'lon',
            'pressure',
            'depth',
            'temperature',
            'conductivity',
            'salinity',
            'density',
            'profile_time',
            'profile_lat',
            'profile_lon',
            'time_uv',
            'lat_uv',
            'lon_uv',
            'u',
            'v'
        ]
        required_attributes = [
            '_FillValue',
            'units',
            'standard_name',
            'observation_type'
        ]
        for var in primary_variables:
            for attribute in required_attributes:
                test = hasattr(ds.dataset.variables[var], attribute)
                out_of += 1
                score += int(test)
                if not test:
                    messages.append('%s variable is missing attribute %s' % (var, attribute))

        return self.make_result(level, score, out_of, 'Required Variable Attributes', messages)
