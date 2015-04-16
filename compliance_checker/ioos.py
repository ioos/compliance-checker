import json
import itertools
from compliance_checker.base import BaseCheck, BaseNCCheck, BaseSOSGCCheck, BaseSOSDSCheck, check_has, score_group, Result
from pkgutil import get_data

class IOOSBaseCheck(BaseCheck):
    register_checker = True
    name = 'ioos'

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
            'operator_email',
            'operator_name',
            'station_description',
            'station_location_lat',
            'station_location_lon',

        ]

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        return [
            'sensor_descriptions',
            'sensor_ids',
            'sensor_names',

            'data_format_template_version',

            'variable_names',
            'variable_units',
            'network_id',
            'operator_sector',
        ]

    @check_has(BaseCheck.LOW)
    def check_suggested(self, ds):
        return [
        ]
