#!/usr/bin/env python
'''
compliance_checker.glider_dac

Compliance Test Suite for the IOOS National Glider Data Assembly Center
https://github.com/ioos/ioosngdac/wiki
'''

from compliance_checker.base import BaseCheck, BaseNCCheck, Result
import numpy as np

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
        out_of = 26
        score = 0
        messages = []

        test = 'lat' in ds.dataset.variables
        score += int(test)
        if not test:
            messages.append('lat is a required variable')
            return self.make_result(level, score, out_of, 'Lat and Lon are Time Series', messages)

        test = 'lon' in ds.dataset.variables
        score += int(test)
        if not test:
            messages.append('lon is a required variable')
            return self.make_result(level, score, out_of, 'Lat and Lon are Time Series', messages)

        required_coordinate_attributes = [
            '_FillValue',
            'ancillary_variables',
            'comment',
            'coordinate_reference_frame',
            'long_name',
            'observation_type',
            'platform',
            'reference',
            'standard_name',
            'units',
            'valid_max',
            'valid_min'
        ]
        for attribute in required_coordinate_attributes:
            test = hasattr(ds.dataset.variables['lat'], attribute)
            score += int(test)
            if not test:
                messages.append('%s attribute is required for lat' % attribute)
            
            test = hasattr(ds.dataset.variables['lon'], attribute)
            score += int(test)
            if not test:
                messages.append('%s attribute is required for lat' % attribute)

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

        required_attributes = [
            'flag_meanings',
            'flag_values',
            'long_name',
            'standard_name',
            'valid_max',
            'valid_min'
        ]
        level = BaseCheck.MEDIUM
        out_of = len(required_variables)
        score = 0
        messages = []
        for variable in required_variables:
            test = variable in ds.dataset.variables
            if not test:
                messages.append("%s is a required qc variable" % variable)
                continue
            for field in required_attributes:
                if not hasattr(ds.dataset.variables[variable], field):
                    messages.append('%s is missing attribute %s' % (variable, field))
                    test = False
                    break
            score += int(test)
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
            'title'
        ]
        level = BaseCheck.MEDIUM
        out_of = 0
        score = 0
        messages = []
        for field in attribute_fields:
            v = getattr(ds.dataset, field, '')
            test = v != ''
            score += int(test)
            out_of += 1
            if not test:
                messages.append('%s global attribute is missing' % field)

            if isinstance(v, basestring):
                test = len(v.strip()) > 0
            else:
                test = True
            score += int(test)
            out_of += 1
            if not test:
                messages.append('%s global attribute can not be empty' % field)
        
        return self.make_result(level, score, out_of, 'Required Global Attributes', messages)

    def check_wmo(self, ds):
        '''
        Verifies that the data has a WMO ID but not necessarily filled out
        '''
        level = BaseCheck.MEDIUM
        score = 0
        out_of = 1
        messages = []
        test = hasattr(ds.dataset, 'wmo_id')
        score += int(test)
        if not test:
            messages.append("WMO ID is a required attribute but can be empty if the dataset doesn't have a WMO ID")

        return self.make_result(level, score, out_of, 'WMO ID', messages)

    def check_summary(self, ds):
        level = BaseCheck.MEDIUM
        out_of = 1
        score = 0
        messages = []
        if hasattr(ds.dataset, 'summary') and ds.dataset.summary:
            score += 1
        else:
            messages.append('Dataset must define summary')
        return self.make_result(level, score, out_of, 'Summary defined', messages)



    def check_primary_variable_attributes(self, ds):
        '''
        Verifies that each primary variable has the necessary metadata
        '''
        level = BaseCheck.MEDIUM
        out_of = 0
        score = 0
        messages = []
        primary_variables = [
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

    def check_dimensions(self, ds):
        '''
        NetCDF files submitted by the individual glider operators contain 2
        dimension variables:
         - time
         - traj
        '''
        level = BaseCheck.HIGH
        score = 0
        messages = []

        required_dimensions = [
            'time',
            'traj_strlen'
        ]
        out_of = len(required_dimensions)

        for dimension in ds.dataset.dimensions:
            test =  dimension in required_dimensions
            score += int(test)
            if not test:
                messages.append('%s is not a valid dimension' % dimension)
        return self.make_result(level, score, out_of, 'Required Dimensions', messages)

    def check_trajectory_variables(self, ds):
        '''
        The trajectory variable stores a character array that identifies the
        deployment during which the data was gathered. This variable is used by
        the DAC to aggregate all individual NetCDF profiles containing the same
        trajectory value into a single trajectory profile data set. This value
        should be a character array that uniquely identifies the deployment and
        each individual NetCDF file from the deployment data set should have
        the same value.
        '''
        level = BaseCheck.MEDIUM
        out_of = 5
        score = 0
        messages = []

        test = 'trajectory' in ds.dataset.variables
        score += int(test)
        if not test:
            messages.append('trajectory variable not found')
            return self.make_result(level, score, out_of, 'Trajectory Variable', messages)
        test = ds.dataset.variables['trajectory'].dimensions == ('traj_strlen',)
        score += int(test)
        if not test:
            messages.append('trajectory has an invalid dimension')
        test = hasattr(ds.dataset.variables['trajectory'], 'cf_role')
        score += int(test)
        if not test:
            messages.append('trajectory is missing cf_role')
        test = hasattr(ds.dataset.variables['trajectory'], 'comment')
        score += int(test)
        if not test:
            messages.append('trajectory is missing comment')
        test = hasattr(ds.dataset.variables['trajectory'], 'long_name')
        score += int(test)
        if not test:
            messages.append('trajectory is missing long_name')
        return self.make_result(level, score, out_of, 'Trajectory Variable', messages)

    def check_time_series_variables(self, ds):
        '''
        Verifies that the time coordinate variable is correct
        '''

        level = BaseCheck.HIGH
        out_of = 16
        score = 0
        messages = []

        test = 'time' in ds.dataset.variables
        score += int(test)
        if not test:
            messages.append('Required coordinate variable time is missing')
            return self.make_result(level, score, out_of, 'Time Series Variable', messages)

        test = ds.dataset.variables['time'].dtype.str == '<f8'
        score += int(test)
        if not test:
            messages.append('Invalid variable type for time, it should be float64')

        test = ds.dataset.variables['time'].ancillary_variables == 'time_qc'
        score += int(test)
        if not test:
            messages.append('Invalid ancillary_variables attribute for time, should be "time_qc"')

        test = ds.dataset.variables['time'].calendar == 'gregorian'
        score += int(test)
        if not test:
            messages.append('Invalid calendar for time, should be "gregorian"')

        test = ds.dataset.variables['time'].long_name == 'Time'
        score += int(test)
        if not test:
            messages.append('Invalid long_name for time, should be "Time"')

        test = ds.dataset.variables['time'].observation_type == 'measured'
        score += int(test)
        if not test:
            messages.append('Invalid observation_type for time, should be "measured"')

        test = ds.dataset.variables['time'].standard_name == 'time'
        score += int(test)
        if not test:
            messages.append('Invalid standard name for time, should be "time"')

        test = hasattr(ds.dataset.variables['time'], 'units')
        score += int(test)
        if not test:
            messages.append('No units defined for time')

        test = 'time_qc' in ds.dataset.variables
        score += int(test)
        if not test:
            messages.append('time_qc is not defined')
            return self.make_result(level, score, out_of, 'Time Series Variable', messages)

        required_time_qc_attributes = [
            '_FillValue',
            'flag_meanings',
            'flag_values',
            'long_name',
            'standard_name',
            'valid_max',
            'valid_min'
        ]
        for attribute in required_time_qc_attributes:
            test = hasattr(ds.dataset.variables['time_qc'], attribute)
            score += int(test)
            if not test:
                messages.append('%s attribute is required for time_qc' % attribute)

        return self.make_result(level, score, out_of, 'Time Series Variable', messages)

    def check_depth_coordinates(self, ds):
        '''
        Verifies that the pressure coordinate/data variable is correct
        '''

        level = BaseCheck.HIGH
        out_of = 34
        score = 0
        messages = []

        coordinates = ['pressure', 'depth']

        for coordinate in coordinates:
            test = coordinate in ds.dataset.variables
            score += int(test)

            if not test:
                messages.append('Required coordinate variable %s is missing' % coordinate)
                return self.make_result(level, score, out_of, 'Depth/Pressure Variables', messages)

        required_values = {
        }

        data_vars = { i: ds.dataset.variables[i] for i in coordinates}
        for key, value in required_values.iteritems():
            for data_var in data_vars:
                if not hasattr(data_vars[data_var], key):
                    messages.append('%s variable requires attribute %s with a value of %s' % (data_var, key, value))
                    continue
                if not test:
                    messages.append('%s.%s != %s' % (data_var, key, value))
                    continue
                score += 1

        required_attributes = [
            '_FillValue',
            'ancillary_variables',
            'accuracy',
            'ancillary_variables',
            'comment',
            'instrument',
            'long_name',
            'platform',
            'positive',
            'precision',
            'reference_datum',
            'resolution',
            'standard_name',
            'units',
            'valid_max',
            'valid_min'
        ]

        for field in required_attributes:
            for data_var in data_vars:
                if not hasattr(data_vars[data_var], field):
                    messages.append('%s variable requires attribute %s' % (data_var, field))
                    continue
                score += 1

        return self.make_result(level, score, out_of, 'Depth/Pressure Variables', messages)

    def check_ctd_variables(self, ds):
        '''
        Verifies that the CTD Variables are the correct data type and contain
        the correct metadata
        '''

        level = BaseCheck.HIGH
        out_of = 56
        score = 0
        messages = []

        variables = [
            'temperature',
            'conductivity',
            'salinity',
            'density'
        ]

        required_fields = [
            'accuracy',
            'ancillary_variables',
            'instrument',
            'long_name',
            'observation_type',
            'platform',
            'precision',
            'resolution',
            'standard_name',
            'units',
            'valid_max',
            'valid_min'
        ]

        for var in variables:
            if var not in ds.dataset.variables:
                messages.append('%s variable missing' % var)
                continue
            nc_var = ds.dataset.variables[var]

            test = nc_var.dtype.str == '<f8'
            score += int(test)
            if not test:
                messages.append('%s variable is incorrect data type' % var)

            for field in required_fields:
                if not hasattr(nc_var, field):
                    messages.append('%s variable is missing required attribute %s' % (var, field))
                    continue
                score += 1

            score += 1
        return self.make_result(level, score, out_of, 'CTD Variables', messages)

    def check_standard_names(self, ds):
        '''
        Verifies that the standard names are correct.
        '''
        level = BaseCheck.MEDIUM
        out_of = 1
        score = 0
        messages = []
        std_names = {
            'salinity' : 'sea_water_practical_salinity'
        }

        for var in std_names:
            if var not in ds.dataset.variables:
                messages.append("Can't verify standard name for %s: %s is missing." % (var, var))
                continue
            nc_var = ds.dataset.variables[var]
            test = getattr(nc_var, 'standard_name', None) == std_names[var]
            score += int(test)
            if not test:
                messages.append("Invalid standard name for %s: %s" % (var, std_names[var]))

        return self.make_result(level, score, out_of, 'Standard Names', messages)


    def check_profile_vars(self, ds):
        '''
        Verifies that the profile variables are the of the correct data type
        and contain the correct metadata
        '''

        level = BaseCheck.MEDIUM
        out_of = 74
        score = 0
        messages = []

        data_struct = {
            'profile_id' : {
                'dtype' : '<i4',
                'fields' : [
                    '_FillValue',
                    'comment',
                    'long_name',
                    'valid_min',
                    'valid_max'
                ]
            },
            'profile_time' : {
                'dtype' : '<f8',
                'fields' : [
                    '_FillValue',
                    'comment',
                    'long_name',
                    'observation_type',
                    'platform',
                    'standard_name',
                    'units'
                ]
            },
            'profile_lat' : {
                'dtype' : '<f8',
                'fields' : [
                    '_FillValue',
                    'comment',
                    'long_name',
                    'observation_type',
                    'platform',
                    'standard_name',
                    'units',
                    'valid_min',
                    'valid_max'
                ]
            },
            'profile_lon' : {
                'dtype' : '<f8',
                'fields' : [
                    '_FillValue',
                    'comment',
                    'long_name',
                    'observation_type',
                    'platform',
                    'standard_name',
                    'units',
                    'valid_min',
                    'valid_max'
                ]
            },
            'lat_uv' : {
                'dtype' : '<f8',
                'fields' : [
                    '_FillValue',
                    'comment',
                    'long_name',
                    'observation_type',
                    'platform',
                    'standard_name',
                    'units',
                    'valid_min',
                    'valid_max'
                ]
            },
            'lon_uv' : {
                'dtype' : '<f8',
                'fields' : [
                    '_FillValue',
                    'comment',
                    'long_name',
                    'observation_type',
                    'platform',
                    'standard_name',
                    'units',
                    'valid_min',
                    'valid_max'
                ]
            },
            'u' : {
                'dtype' : '<f8',
                'fields' : [
                    '_FillValue',
                    'comment',
                    'long_name',
                    'observation_type',
                    'platform',
                    'standard_name',
                    'units',
                    'valid_min',
                    'valid_max'
                ]
            },
            'v' : {
                'dtype' : '<f8',
                'fields' : [
                    '_FillValue',
                    'comment',
                    'long_name',
                    'observation_type',
                    'platform',
                    'standard_name',
                    'units',
                    'valid_min',
                    'valid_max'
                ]
            }
        }

        for profile_var in data_struct:
            test = profile_var in ds.dataset.variables
            if not test:
                messages.append("Required Variable %s is missing" % profile_var)
                continue

            nc_var = ds.dataset.variables[profile_var]
            dtype = np.dtype(data_struct[profile_var]['dtype'])
            test = nc_var.dtype.str == dtype.str
            score += int(test)
            if not test:
                messages.append('%s variable has incorrect dtype, should be %s' % (profile_var, dtype.name))

            for field in data_struct[profile_var]['fields']:
                test = hasattr(nc_var, field)
                if not test:
                    messages.append('%s variable is missing required attribute %s' % (profile_var, field))
                    continue
                score += 1

        return self.make_result(level, score, out_of, 'Profile Variables', messages)

    def check_container_variables(self, ds):
        '''
        Verifies that the dimensionless container variables are the correct
        data type and contain the required metadata
        '''
        
        level = BaseCheck.MEDIUM
        out_of = 19
        score = 0
        messages = []
        
        data_struct = {
            'platform' : {
                'dtype' : '<i4',
                'fields' : [
                    '_FillValue',
                    'comment',
                    'id',
                    'instrument',
                    'long_name',
                    'type',
                    'wmo_id'
                ]
            },
            'instrument_ctd' : {
                'dtype' : '<i4',
                'fields' : [
                    '_FillValue',
                    'calibration_date',
                    'calibration_report',
                    'comment',
                    'factory_calibrated',
                    'long_name',
                    'make_model',
                    'platform',
                    'serial_number',
                    'type'
                ]
            }
        }

        for container_var in data_struct:
            test = container_var in ds.dataset.variables
            if not test:
                messages.append("Required Variable %s is missing" % container_var)
                continue

            nc_var = ds.dataset.variables[container_var]
            dtype = np.dtype(data_struct[container_var]['dtype'])
            test = nc_var.dtype.str == dtype.str
            score += int(test)
            if not test:
                messages.append('%s variable has incorrect dtype, should be %s' % (container_var, dtype.name))

            for field in data_struct[container_var]['fields']:
                test = hasattr(nc_var, field)
                if not test:
                    messages.append('%s variable is missing required attribute %s' % (container_var, field))
                    continue
                score += 1

        return self.make_result(level, score, out_of, 'Container Variables', messages)
