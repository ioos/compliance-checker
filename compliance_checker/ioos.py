import json
import itertools
from compliance_checker.base import BaseCheck, BaseNCCheck, BaseSOSGCCheck, BaseSOSDSCheck, check_has, score_group, Result
from pkgutil import get_data
from datetime import datetime
import time
from lxml import etree as ET


class IOOSBaseCheck(BaseCheck):
    pass

class IOOSNCCheck(BaseNCCheck, IOOSBaseCheck):
    # belefs
    @classmethod
    def beliefs(cls):
        f = get_data("compliance_checker", "data/ioos-metamap-ncml.json")
        beliefs = json.loads(f)

        # strip out metadata
        return {k:v for k,v in beliefs.iteritems() if not k.startswith("__")}

    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        return [
            'platform_sponsor',
            'platform_type',
            'station_publisher_name',
            'station_publisher_email',
            'station_id',
            'station_long_name',
            'station_short_name',
            'station_wmo_id',
            'time_period',
            'station_description',
            'station_location_lat',
            'station_location_lon',

        ]

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        return [
            # @TODO: only in sensors

            'service_contact_email',
            'service_contact_name',
            'service_provider_name',

            'data_format_template_version',

            'variable_names',
            'variable_units'
        ]

    @check_has(BaseCheck.LOW)
    def check_suggested(self, ds):
        return [
            'altitude_units'
        ]


class IOOSSOSGCCheck(BaseSOSGCCheck, IOOSBaseCheck):
    # beliefs
    @classmethod
    def beliefs(cls):
        f = get_data("compliance_checker", "data/ioos-metamap-sos-gc.json")
        beliefs = json.loads(f)

        # strip out metadata
        return {k:v for k,v in beliefs.iteritems() if not k.startswith("__")}

    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        return [

        ]

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        return [
            'service_contact_email',
            'service_contact_name',
            'service_provider_name',

            'service_title',
            'service_type_name',
            'service_type_version',

            'data_format_template_version',
            'variable_names',
        ]

    @check_has(BaseCheck.LOW)
    def check_suggested(self, ds):
        return [
            'altitude_units'
        ]

class IOOSSOSDSCheck(BaseSOSDSCheck, IOOSBaseCheck):
    # beliefs
    @classmethod
    def beliefs(cls):
        f = get_data("compliance_checker", "data/ioos-metamap-sos-ds.json")
        beliefs = json.loads(f)

        # strip out metadata
        return {k:v for k,v in beliefs.iteritems() if not k.startswith("__")}

    def check_valid_network_id(self,ds):
        if hasattr(ds.dogma,'network_id'):
            asset_type = ['network']
            authority = ['wmo', 'noaa.nos.co-ops','usace','fiu','test']
            label = []
            component = []
    
    
            network_id = getattr(ds.dogma,'network_id', '')
            network_id_sep = network_id.split(':')
        
            valid = False
            fails = ['The URN convention for network_id does not match IOOS Specification.']
            if str(network_id_sep[0]) == "urn" and str(network_id_sep[1]) == "ioos" and str(network_id_sep[2]).lower() in asset_type and str(network_id_sep[3]).lower() in authority:
                valid = True
                fails = []
    
            return Result(BaseCheck.HIGH, valid, 'network_id_convention', fails)
        else:
            pass

    def check_short_name(self,ds):
        short_name = getattr(ds.dogma,'short_name', '')
        valid = False
        fails = ['The shortName value is incorrect.']
        if type(short_name) == type(str()) and len(short_name) > 0:
            fails = []
            valid = True
        return Result(BaseCheck.HIGH, valid, 'short_name_convention', fails)

    def check_long_name(self,ds):
        long_name = getattr(ds.dogma,'long_name', '')
        valid = False
        fails = ['The longName value is incorrect.']
        if type(long_name) == type(str()) and len(long_name) > 0:
            fails = []
            valid = True
        return Result(BaseCheck.HIGH, valid, 'long_name_convention', fails)

    def check_lat_lon_values(self,ds):
        ret_val = []
        if hasattr(ds.dogma, 'network_station_locations'):
            station_locs = getattr(ds.dogma,'network_station_locations', '')
            if hasattr(ds.dogma, 'lower_corner') and hasattr(ds.dogma, 'upper_corner'):
                lower_corner = getattr(ds.dogma,'lower_corner', '').split(' ')
                upper_corner = getattr(ds.dogma,'upper_corner', '').split(' ')

            valid = 0
            checks = 0
            fails = []
            for each in station_locs:
                checks = checks+1
                each = each.split(' ')

                if (float(each[0]) >= float(lower_corner[0]) and float(each[1]) >= float(lower_corner[1])):
                    if (float(each[0]) <= float(upper_corner[0]) and float(each[1]) <= float(upper_corner[1])):
                        valid = valid + 1
                    else:
                        fails.append('The station locations are above the bounding box.')
                else:
                    fails.append('The station locations are below the bounding box.')


            ret_val.append(Result(BaseCheck.HIGH, (valid,checks), 'legal_station_locations', fails))
            ret_val.append(Result(BaseCheck.HIGH, (valid,checks), 'station_locations', fails))
            return ret_val

        else:
            pass

    def check_time_values(self,ds):
        if hasattr(ds.dogma, 'network_station_times'):
            station_times = getattr(ds.dogma,'network_station_times', '')
    
            valid = 0
            checks = 0
            fails = []
            for each in station_times:
                failed = 0
                valid_every = 0
                checks = checks+1
                each = each.split(' ')
                for every in each:
                    try:
                        datetime.strptime(every, '%Y-%m-%dT%H:%M:%SZ')
                        valid_every = valid_every + 1
                    except:
                        try: 
                            datetime.strptime(every, '%Y-%m-%dT%H:%M:%S.%fZ')
                            valid_every = valid_every + 1
                        except:
                            failed = 1

                    if valid_every == 2:
                        valid = valid+1
                if failed == 1:
                    fails.append('The station time is not in proper format (YYYY-MM-DDTHH:MM:SS)')


            return Result(BaseCheck.HIGH, (valid,checks), 'network_station_times', fails)
        else:
            pass
    def check_network_station_ids(self,ds):
        if hasattr(ds.dogma, 'network_id'):
            fails = []
            ret_val = []
    
            valid_stationID = 0
            checks = 0
            valid = 0
            failed = 0
    
            for each in  ds.dataset._root.xpath("//*[local-name() = 'component']"):
    
                every_component = 1
                every_stationID = 0
                checks = checks + 1
    
                for every in each.iterdescendants():
                    if every.attrib == {'name':'stationID'}:
                        every_stationID = every_stationID +1


                if every_component <= every_stationID:
                    valid_stationID = valid_stationID+1
                else:
                    fails.append('Not every station has an ID')
    
    
            return Result(BaseCheck.HIGH, (valid_stationID,checks) , 'network_station_ids', fails)
        else:
            pass

    def check_oberved_properties(self,ds):
        fails = []
        ret_val = []

        valid_output = 0
        valid_quantity = 0
        checks = 0
        valid = 0
        failed = 0

        for each in  ds.dataset._root.xpath("//*[local-name() = 'component']"):

            every_component = 1
            every_output = 0
            every_quantity = 0
            checks = checks + 1

            for every in each.iterdescendants():
                if every.tag == '{http://www.opengis.net/sensorML/1.0.1}output':
                    every_output = every_output +1
                if every.tag == '{http://www.opengis.net/swe/1.0.1}Quantity':
                    every_quantity = every_quantity +1


            if every_component <= every_output:
                valid_output = valid_output+1
            else:
                fails.append('Not every component has an output.')


            valid_quantity = valid_quantity+every_output
            if every_output == every_quantity:
                pass
            else:
                fails.append('Not every output has a quantity field.')
                failed = failed+1

        ret_val = [Result(BaseCheck.HIGH, (valid_output,checks) , 'observed_output_count', fails),Result(BaseCheck.HIGH, (valid_quantity-failed,valid_quantity) , 'observed_quantities_count', fails)]
        return ret_val

    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        if hasattr(ds.dogma,'network_id'):
            return [
                'ioosservicemetadata',
                'epsg',
                'short_name',
                'long_name',
                'parent_network',
                'station_publisher_name',
                'station_publisher_email',
                'station_publisher_country',
                'station_long_name',
                'station_short_name',
                'station_wmo_id',
                'time_period',
                'operator_email',
                'operator_country',
                'operator_name',
                'station_description',
                
            ]
        else:
            return [
                'ioosservicemetadata',
                'short_name',
                'long_name',
                'parent_network',
                'platform_type',
                'station_publisher_name',
                'station_publisher_email',
                'station_publisher_country',
                'station_id',
                'station_long_name',
                'station_short_name',
                'station_wmo_id',
                'time_period',
                'operator_email',
                'operator_country',
                'operator_name',
                'station_description',
                'station_location'
                
            ]

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        if hasattr(ds.dogma,'network_id'):
            return [
            'sensor_descriptions',
            'sensor_ids',
            'sensor_names',
            'data_format_template_version',
            'variable_names',
            'variable_units',
            'network_id'
        ]
        else:
            return [
            'sensor_descriptions',
            'sensor_ids',
            'sensor_names',
            'data_format_template_version',
            'variable_names',
            'variable_units',
            'operator_sector'
            ]
    @check_has(BaseCheck.LOW)
    def check_suggested(self, ds):
        return [
        ]
    
