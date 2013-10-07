import json
import itertools
from compliance_checker.base import BaseCheck, check_has, score_group

class ACDDCheck(BaseCheck):

    @classmethod
    def beliefs(cls):
        with open("acdd-ncml.json") as f:
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
        return ['title', 'summary', 'keywords']

    ###############################################################################
    #
    # RECOMMENDED
    #
    ###############################################################################

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        return [
            'id',
            'naming_authority',
            'keywords_vocabulary',
            ('cdm_data_type', ['vector', 'grid', 'textTable', 'tin', 'stereoModel', 'video']),
            'history',
            'comment',
            'date_created',
            'creator_name',
            'creator_url',
            'creator_email',
            'institution',
            'project',
            'processing_level',
            'acknowledgement',
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
            'standard_name_vocabulary',
            'license'
        ]

    ###############################################################################
    #
    # SUGGESTED
    #
    ###############################################################################

    @check_has(BaseCheck.LOW)
    def check_suggested(self, ds):
        return [
            'contributor_name',
            ('contributor_role', ['principalInvestigator', 'author']),
            'publisher_name',       # publisher,dataCenter
            'publisher_url',        # publisher
            'publisher_email',      # publisher
            'date_modified',
            'date_issued',
            'geospatial_lat_units',
            'geospatial_lat_resolution',
            'geospatial_lon_units',
            'geospatial_lon_resolution',
            'geospatial_vertical_units',
            'geospatial_vertical_resolution',
            'geospatial_vertical_positive'
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
        vars = self._get_vars(ds, 'long_name')

        retval = [(BaseCheck.HIGH, v[0] is not None, (v[1], "var_long_name")) for v in vars]
        return retval

    @score_group('varattr')
    def check_var_standard_name(self, ds):
        vars = self._get_vars(ds, 'standard_name')

        retval = [(BaseCheck.HIGH, v[0] is not None, (v[1], "var_std_name")) for v in vars]
        return retval

    @score_group('varattr')
    def check_var_units(self, ds):
        vars = self._get_vars(ds, 'units')

        retval = [(BaseCheck.HIGH, v[0] is not None, (v[1], "var_units")) for v in vars]
        return retval

    @score_group('varattr')
    def check_var_coverage_content_type(self, ds):
        vars = self._get_vars(ds, 'coverage_content_type')
        allowed = ['image','thematicClassification','physicalMeasurement','auxiliaryInformation','qualityInformation','referenceInformation','modelResult','coordinate']

        retval = [(BaseCheck.HIGH, v[0] is not None and v[0] in allowed, (v[1], "var_coverage_content_type")) for v in vars]
        return retval

