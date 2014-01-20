import json
import itertools
from compliance_checker.base import BaseCheck, check_has, score_group, Result

class ioosCheck(BaseCheck):

    @classmethod
    def beliefs(cls):
        with open("ioos-ncml.json") as f:
            beliefs = json.load(f)

        # strip out metadata
        return {k:v for k,v in beliefs.iteritems() if not k.startswith("__")}

    ###############################################################################
    #
    # HIGHLY RECOMMENDED
    # 
    ###############################################################################

    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        return [
            'operator_email',
            'operator_name',
            'platform_sponsor',
            'platform_type',
            'publisher',
            'publisher_email',
            'station_id',
            'station_long_name',
            'station_name',
            'station_short_name',
            'station_wmo_id',
            'time_period_begin',
            'time_period_end',
            'time_period_interval',
            'time_period_interval_units',
                        ]

    ###############################################################################
    #
    # RECOMMENDED
    #
    ###############################################################################

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        return [
            'location_units',
            'observed_property',
            'observed_property_time_last',
            'operator_address',
            'operator_sector',
            'sensor_descriptions',
            'sensor_ids',
            'sensor_names',
            'service_contact_address',
            'service_contact_email',
            'service_contact_name',
            'service_keyword_1',
            'service_provider_name',
            'service_title',
            'service_type_name',
            'service_type_version',
            'sos_template_version',
            'station_deployment_end',
            'station_deployment_start',
            'variable_altitudes',
            'variable_names',
            'variable_units'
        ]

    ###############################################################################
    #
    # SUGGESTED
    #
    ###############################################################################

    @check_has(BaseCheck.LOW)
    def check_suggested(self, ds):
        return [
            'altitude_units'
                ]

    ###############################################################################
    #
    # HIGHLY RECOMMENDED VARIABLE ATTRS
    #
    ###############################################################################

    def _get_vars(self, ds, attr_filter=None):
        vars = ds.dogma._eval_xpath('//ncml:variable')

        if attr_filter is not None:
            attrs = itertools.chain.from_iterable((v.xpath('ncml:attribute[@name="%s"]/@value' % attr_filter, namespaces=ds.dogma._namespaces) or [None] for v in vars))
            names = (v.attrib.get('name', 'unknown') for v in vars)

            attrs = zip(attrs, names)

            return attrs

        return vars

    @score_group('varattr')
    def check_var_long_name(self, ds):
        vars = self._get_vars(ds, #Variables])

        retval = [Result(BaseCheck.HIGH, v[0] is not None, (v[1], #Variables])) for v in vars]
        return retval

    @score_group('varattr')
    def check_var_standard_name(self, ds):
        vars = self._get_vars(ds, #Variables) 

        retval = [Result(BaseCheck.HIGH, v[0] is not None, (v[1], #Variables])) for v in vars]
        return retval

    @score_group('varattr')
    def check_var_units(self, ds):
        vars = self._get_vars(ds, #Variables])

        retval = [Result(BaseCheck.HIGH, v[0] is not None, (v[1], #Variables])) for v in vars]
        return retval

    @score_group('varattr')
    def check_var_coverage_content_type(self, ds):
        vars = self._get_vars(ds, #Variables)
        allowed = [#Variables]

        retval = [Result(BaseCheck.HIGH, v[0] is not None and v[0] in allowed, (v[1], #Variables])) for v in vars]
        return retval

