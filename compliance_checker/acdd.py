import json
import itertools
import numpy as np
from pkgutil import get_data

from dateutil.parser import parse as parse_dt
from cf_units import Unit

from compliance_checker.base import BaseCheck, BaseNCCheck, check_has, score_group, Result
from compliance_checker.cf.cf import _possiblexunits, _possibleyunits
from compliance_checker.cf.util import is_time_variable, is_vertical_coordinate

class ACDDBaseCheck(BaseCheck):

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
            names = (v.get('name', 'unknown') for v in vars)

            attrs = zip(attrs, names)

            return attrs

        return vars

    def _get_msg(self, vpair, attr):
        vval, vname = vpair
        if vval is None:
            return ["Var %s missing attr %s" % (vname, attr)]

        return []

    @score_group('varattr')
    def check_var_long_name(self, ds):
        vars = self._get_vars(ds, 'long_name')

        retval = [Result(BaseCheck.HIGH, v[0] is not None, (v[1], "var_long_name"), self._get_msg(v, 'long_name')) for v in vars]
        return retval

    @score_group('varattr')
    def check_var_standard_name(self, ds):
        vars = self._get_vars(ds, 'standard_name')

        retval = [Result(BaseCheck.HIGH, v[0] is not None, (v[1], "var_std_name"), self._get_msg(v, 'standard_name')) for v in vars]
        return retval

    @score_group('varattr')
    def check_var_units(self, ds):
        vars = self._get_vars(ds, 'units')

        retval = [Result(BaseCheck.HIGH, v[0] is not None, (v[1], "var_units"), self._get_msg(v, 'units')) for v in vars]
        return retval

    @score_group('varattr')
    def check_var_coverage_content_type(self, ds):
        vars = self._get_vars(ds, 'coverage_content_type')
        allowed = ['image','thematicClassification','physicalMeasurement','auxiliaryInformation','qualityInformation','referenceInformation','modelResult','coordinate']

        ret_val = []
        for v in vars:
            vval, vname = v
            msgs = []

            count = 0
            if vval is not None:
                count += 1

                if vval in allowed:
                    count += 1
                else:
                    msgs.append("coverage_content_type present but value (%s) not in allowed values (%s)" % (vval, allowed))
            else:
                msgs.append("Var %s missing attr %s" % (vname, 'coverage_content_type'))

            ret_val.append(Result(BaseCheck.HIGH,
                                  (count, 2),
                                  (v[1], "var_coverage_content_type"),
                                  msgs))

        return ret_val

    ###############################################################################
    #
    # DATA EXTENTS MATCHING ATTRIBUTES
    #
    ###############################################################################

    def check_lat_extents(self, ds):
        """
        Check that the values of geospatial_lat_min/geospatial_lat_max approximately match the data.
        """
        if not (hasattr(ds.dataset, 'geospatial_lat_min') and hasattr(ds.dataset, 'geospatial_lat_max')):
            return

        lat_min = ds.dataset.geospatial_lat_min
        lat_max = ds.dataset.geospatial_lat_max

        # identify lat var(s) as per CF 4.1
        lat_vars = {}       # var -> number of criteria passed
        for name, var in ds.dataset.variables.iteritems():

            # must have units
            if not hasattr(var, 'units'):
                continue

            lat_vars[var] = 0

            # units in this set
            if var.units in _possibleyunits:
                lat_vars[var] += 1

            # standard name of "latitude"
            if hasattr(var, 'standard_name') and var.standard_name == 'latitude':
                lat_vars[var] += 1

            # axis of "Y"
            if hasattr(var, 'axis') and var.axis == 'Y':
                lat_vars[var] += 1

        # trim out any zeros
        lat_vars = {k:v for k, v in lat_vars.iteritems() if v > 0}

        if len(lat_vars) == 0:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'geospatial_lat_extents_match',
                          'Could not find lat variable to test extent of geospatial_lat_min/max, see CF-1.6 spec chapter 4.1')

        # sort by criteria passed
        final_lats = sorted(lat_vars, key=lambda x: lat_vars[x], reverse=True)

        obs_mins = {var._name:np.nanmin(var) for var in final_lats if not np.isnan(var).all()}
        obs_maxs = {var._name:np.nanmax(var) for var in final_lats if not np.isnan(var).all()}

        min_pass = any((np.isclose(lat_min, min_val) for min_val in obs_mins.itervalues()))
        max_pass = any((np.isclose(lat_max, max_val) for max_val in obs_maxs.itervalues()))

        allpass = sum((min_pass, max_pass))

        msgs = []
        if not min_pass:
            msgs.append("Data for possible latitude variables (%s) did not match geospatial_lat_min value (%s)" % (obs_mins, lat_min))
        if not max_pass:
            msgs.append("Data for possible latitude variables (%s) did not match geospatial_lat_max value (%s)" % (obs_maxs, lat_max))

        return Result(BaseCheck.MEDIUM,
                      (allpass, 2),
                      'geospatial_lat_extents_match',
                      msgs)

    def check_lon_extents(self, ds):
        """
        Check that the values of geospatial_lon_min/geospatial_lon_max approximately match the data.
        """
        if not (hasattr(ds.dataset, 'geospatial_lon_min') and hasattr(ds.dataset, 'geospatial_lon_max')):
            return

        lon_min = ds.dataset.geospatial_lon_min
        lon_max = ds.dataset.geospatial_lon_max

        # identify lon var(s) as per CF 4.2
        lon_vars = {}       # var -> number of criteria passed
        for name, var in ds.dataset.variables.iteritems():

            # must have units
            if not hasattr(var, 'units'):
                continue

            lon_vars[var] = 0

            # units in this set
            if var.units in _possiblexunits:
                lon_vars[var] += 1

            # standard name of "longitude"
            if hasattr(var, 'standard_name') and var.standard_name == 'longitude':
                lon_vars[var] += 1

            # axis of "Y"
            if hasattr(var, 'axis') and var.axis == 'X':
                lon_vars[var] += 1

        # trim out any zeros
        lon_vars = {k:v for k, v in lon_vars.iteritems() if v > 0}

        if len(lon_vars) == 0:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'geospatial_lon_extents_match',
                          'Could not find lon variable to test extent of geospatial_lon_min/max, see CF-1.6 spec chapter 4.2')

        # sort by criteria passed
        final_lons = sorted(lon_vars, key=lambda x: lon_vars[x], reverse=True)

        obs_mins = {var._name:np.nanmin(var) for var in final_lons if not np.isnan(var).all()}
        obs_maxs = {var._name:np.nanmax(var) for var in final_lons if not np.isnan(var).all()}

        min_pass = any((np.isclose(lon_min, min_val) for min_val in obs_mins.itervalues()))
        max_pass = any((np.isclose(lon_max, max_val) for max_val in obs_maxs.itervalues()))

        allpass = sum((min_pass, max_pass))

        msgs = []
        if not min_pass:
            msgs.append("Data for possible longitude variables (%s) did not match geospatial_lon_min value (%s)" % (obs_mins, lon_min))
        if not max_pass:
            msgs.append("Data for possible longitude variables (%s) did not match geospatial_lon_max value (%s)" % (obs_maxs, lon_max))

        return Result(BaseCheck.MEDIUM,
                      (allpass, 2),
                      'geospatial_lon_extents_match',
                      msgs)

    def check_vertical_extents(self, ds):
        """
        Check that the values of geospatial_vertical_min/geospatial_vertical_max approximately match the data.
        """
        if not (hasattr(ds.dataset, 'geospatial_vertical_min') and hasattr(ds.dataset, 'geospatial_vertical_max')):
            return

        vert_min = ds.dataset.geospatial_vertical_min
        vert_max = ds.dataset.geospatial_vertical_max

        # identify vertical vars as per CF 4.3
        v_vars = [var for name, var in ds.dataset.variables.iteritems() if is_vertical_coordinate(name, var)]

        if len(v_vars) == 0:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'geospatial_vertical_extents_match',
                          'Could not find vertical variable to test extent of geospatial_vertical_min/geospatial_vertical_max, see CF-1.6 spec chapter 4.3')

        obs_mins = {var._name:np.nanmin(var) for var in v_vars if not np.isnan(var).all()}
        obs_maxs = {var._name:np.nanmax(var) for var in v_vars if not np.isnan(var).all()}

        min_pass = any((np.isclose(vert_min, min_val) for min_val in obs_mins.itervalues()))
        max_pass = any((np.isclose(vert_max, max_val) for max_val in obs_maxs.itervalues()))

        allpass = sum((min_pass, max_pass))

        msgs = []
        if not min_pass:
            msgs.append("Data for possible vertical variables (%s) did not match geospatial_vertical_min value (%s)" % (obs_mins, vert_min))
        if not max_pass:
            msgs.append("Data for possible vertical variables (%s) did not match geospatial_vertical_max value (%s)" % (obs_maxs, vert_max))

        return Result(BaseCheck.MEDIUM,
                      (allpass, 2),
                      'geospatial_vertical_extents_match',
                      msgs)


    def check_time_extents(self, ds):
        """
        Check that the values of time_coverage_start/time_coverage_end approximately match the data.
        """
        if not (hasattr(ds.dataset, 'time_coverage_start') and hasattr(ds.dataset, 'time_coverage_end')):
            return

        epoch = parse_dt("1970-01-01 00:00:00 UTC")
        t_min = (parse_dt(ds.dataset.time_coverage_start) - epoch).total_seconds()
        t_max = (parse_dt(ds.dataset.time_coverage_end) - epoch).total_seconds()

        # identify t vars as per CF 4.4
        t_vars = [var for name, var in ds.dataset.variables.iteritems() if is_time_variable(name, var)]

        if len(t_vars) == 0:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'time_coverage_extents_match',
                          'Could not find time variable to test extent of time_coverage_start/time_coverage_end, see CF-1.6 spec chapter 4.4')

        obs_mins = {var._name:Unit(str(var.units)).get_converter("seconds since 1970-01-01").evaluate(np.nanmin(var)) for var in t_vars}
        obs_maxs = {var._name:Unit(str(var.units)).get_converter("seconds since 1970-01-01").evaluate(np.nanmax(var)) for var in t_vars}

        min_pass = any((np.isclose(t_min, min_val) for min_val in obs_mins.itervalues()))
        max_pass = any((np.isclose(t_max, max_val) for max_val in obs_maxs.itervalues()))

        allpass = sum((min_pass, max_pass))

        msgs = []
        if not min_pass:
            msgs.append("Data for possible time variables (%s) did not match time_coverage_start value (%s)" % (obs_mins, t_min))
        if not max_pass:
            msgs.append("Data for possible time variables (%s) did not match time_coverage_end value (%s)" % (obs_maxs, t_max))

        return Result(BaseCheck.MEDIUM,
                      (allpass, 2),
                      'time_coverage_extents_match',
                      msgs)

class ACDDNCCheck(BaseNCCheck, ACDDBaseCheck):
    @classmethod
    def beliefs(cls):
        f = get_data("compliance_checker", "data/acdd-ncml.json")
        beliefs = json.loads(f)

        # strip out metadata
        return {k:v for k,v in beliefs.iteritems() if not k.startswith("__")}

