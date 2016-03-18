import itertools
import numpy as np

from dateutil.parser import parse as parse_dt
from cf_units import Unit

from compliance_checker.base import BaseCheck, BaseNCCheck, check_has, score_group, Result
from compliance_checker.cf.util import is_time_variable, is_vertical_coordinate, _possiblexunits, _possibleyunits


from pygeoif import from_wkt

class ACDDBaseCheck(BaseCheck):

    register_checker = True
    name = 'acdd'

    _supported_versions = {'1.1', '1.3'}

    def __init__(self, version='1.1'):
        if version in self._supported_versions:
            self._cc_spec_version = version
        else:
            raise NotImplementedError("Version {} not found in valid versions".format(version))

        common_rec_atts = [
            'id',
            'naming_authority',
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
            'standard_name_vocabulary',
            'license'
        ]

        common_sug_atts = [
            'contributor_name',
            'contributor_role',
            'date_modified',
            'date_issued',
            'geospatial_lat_units',
            'geospatial_lat_resolution',
            'geospatial_lon_units',
            'geospatial_lon_resolution',
            'geospatial_vertical_units',
            'geospatial_vertical_resolution'
        ]

        if self._cc_spec_version == '1.1':
            self.high_rec_atts = ['title',
                    'summary',
                    'keywords']

            common_rec_atts.extend([
                            'keywords_vocabulary',
                            ('cdm_data_type', ['Grid', 'Image', 'Point',
                                               'Radial', 'Station', 'Swath',
                                               'Trajectory'])])

            common_sug_atts.extend([
                                'publisher_name',       # publisher,dataCenter
                                'publisher_url',        # publisher
                                'publisher_email',      # publisher
                                'geospatial_vertical_positive'
                              ])

        elif self._cc_spec_version == '1.3':
            self.high_rec_atts = ['title',
                    'summary',
                    'keywords',
                    # TODO: Requires at least 'ACDD-1.3' to be present
                    # in attribute
                    ('Conventions', ['ACDD-1.3'])]

            common_rec_atts.extend(
                ['geospatial_vertical_positive',
                 'geospatial_bounds_crs',
                 'geospatial_bounds_vertical_crs',
                 'publisher_name',       # publisher,dataCenter
                 'publisher_url',        # publisher
                 'publisher_email',      # publisher
                 'source'])

            common_sug_atts.extend([
                # 1.3.1, technically
                ('creator_type', ['person', 'group', 'institution',
                                  'position']),
                'creator_institution',
                ('cdm_data_type', ['Grid', 'Image', 'Point', 'Radial',
                                   'Station', 'Swath', 'Trajectory']),
                'platform',
                # TODO: make dependent on platform
                'platform_vocabulary',
                'keywords_vocabulary',
                'instrument',
                'metadata_link',
                'product_version',
                'references',
                ('publisher_type', ['person', 'group', 'institution',
                                    'position']),
                'instrument_vocabulary',
                'date_metadata_modified',
                'program',
                'publisher_institution',
            ])

        self.rec_atts = common_rec_atts
        self.sug_atts = common_sug_atts
    ###############################################################################
    #
    # HIGHLY RECOMMENDED
    #
    ###############################################################################
        # set up attributes accoriding to version
    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        return self.high_rec_atts
    ###############################################################################
    #
    # RECOMMENDED
    #
    ###############################################################################

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        return self.rec_atts
    ###############################################################################
    #
    # SUGGESTED
    #
    ###############################################################################

    @check_has(BaseCheck.LOW)
    def check_suggested(self, ds):
        return self.sug_atts

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

            attrs = list(zip(attrs, names))

            return attrs

        return vars

    def _get_msg(self, vpair, attr):
        vval, vname = vpair
        if vval is None:
            return ["Var %s missing attr %s" % (vname, attr)]

        return []

    @score_group('varattr')
    def check_var_long_name(self, ds):
        results = []
        # We don't check certain container variables for units
        platform_variable_name = getattr(ds, 'platform', None)
        for variable in ds.variables:
            msgs = []
            if variable in ('crs', platform_variable_name):
                continue
            # If the variable is a QC flag, we don't need units
            long_name = getattr(ds.variables[variable], 'long_name', None)
            check = long_name is not None
            if not check:
                msgs.append("Var %s missing attr long_name" % variable)
            results.append(Result(BaseCheck.HIGH, check, (variable, "var_std_name"), msgs))

        return results

    @score_group('varattr')
    def check_var_standard_name(self, ds):
        results = []
        # We don't check certain container variables for units
        platform_variable_name = getattr(ds, 'platform', None)
        for variable in ds.variables:
            msgs = []
            if variable in ('crs', platform_variable_name):
                continue
            # If the variable is a QC flag, we don't need units
            std_name = getattr(ds.variables[variable], 'standard_name', None)
            check = std_name is not None
            if not check:
                msgs.append("Var %s missing attr standard_name" % variable)
            results.append(Result(BaseCheck.HIGH, check, (variable, "var_std_name"), msgs))

        return results

    @score_group('varattr')
    def check_var_coverage_content_type(self, ds):
        results = []
        platform_variable_name = getattr(ds, 'platform', None)
        for variable in ds.variables:
            msgs = []
            if variable in {'crs', platform_variable_name}:
                continue
            ctype = getattr(ds.variables[variable],
                            'coverage_content_type', None)
            check = ctype is not None
            if not check:
                msgs.append("Var %s missing attr coverage_content_type" %
                            variable)
                results.append(Result(BaseCheck.HIGH, check,
                                    (variable, "coverage_content_type"),
                                    msgs))
                return results
            # ISO 19115-1 codes
            valid_ctypes = {'image', 'thematicClassification', 'physicalMeasurement',
                            'auxiliaryInformation', 'qualityInformation',
                            'referenceInformation', 'modelResult', 'coordinate'}
            if not ctype in valid_ctypes:
                msgs.append("Var %s does not have a coverage_content_type in %s"
                            % (variable, sorted(valid_ctypes)))

        return results

    @score_group('varattr')
    def check_var_units(self, ds):
        results = []
        # We don't check certain container variables for units
        platform_variable_name = getattr(ds, 'platform', None)
        for variable in ds.variables:
            msgs = []
            if variable in ('crs', platform_variable_name):
                continue
            # If the variable is a QC flag, we don't need units
            std_name = getattr(ds.variables[variable], 'standard_name', None)
            if std_name is not None:
                if 'status_flag' in std_name:
                    continue
            # Check units and dims for variable
            unit_check = hasattr(ds.variables[variable], 'units')
            no_dim_check = (getattr(ds.variables[variable], 'dimensions') == tuple())
            # Check if we have no dimensions.  If no dims, skip test
            if no_dim_check:
                continue
            # Check if we have no units
            if not unit_check:
                msgs.append("Var %s missing attr units" % variable)
            results.append(Result(BaseCheck.HIGH, unit_check, (variable, "var_units"), msgs))

        return results

    ###############################################################################
    #
    # DATA EXTENTS MATCHING ATTRIBUTES
    #
    ###############################################################################

    def check_lat_extents(self, ds):
        """
        Check that the values of geospatial_lat_min/geospatial_lat_max approximately match the data.
        """
        if not (hasattr(ds, 'geospatial_lat_min') and hasattr(ds, 'geospatial_lat_max')):
            return

        lat_min = ds.geospatial_lat_min
        lat_max = ds.geospatial_lat_max

        # identify lat var(s) as per CF 4.1
        lat_vars = {}       # var -> number of criteria passed
        for name, var in ds.variables.items():

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
        lat_vars = {k: v for k, v in lat_vars.items() if v > 0}

        if len(lat_vars) == 0:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'geospatial_lat_extents_match',
                          ['Could not find lat variable to test extent of geospatial_lat_min/max, see CF-1.6 spec chapter 4.1'])

        # sort by criteria passed
        final_lats = sorted(lat_vars, key=lambda x: lat_vars[x], reverse=True)

        obs_mins = {var._name: np.nanmin(var) for var in final_lats if not np.isnan(var).all()}
        obs_maxs = {var._name: np.nanmax(var) for var in final_lats if not np.isnan(var).all()}

        min_pass = any((np.isclose(lat_min, min_val) for min_val in obs_mins.values()))
        max_pass = any((np.isclose(lat_max, max_val) for max_val in obs_maxs.values()))

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
        if not (hasattr(ds, 'geospatial_lon_min') and hasattr(ds, 'geospatial_lon_max')):
            return

        lon_min = ds.geospatial_lon_min
        lon_max = ds.geospatial_lon_max

        # identify lon var(s) as per CF 4.2
        lon_vars = {}       # var -> number of criteria passed
        for name, var in ds.variables.items():

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
        lon_vars = {k: v for k, v in lon_vars.items() if v > 0}

        if len(lon_vars) == 0:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'geospatial_lon_extents_match',
                          ['Could not find lon variable to test extent of geospatial_lon_min/max, see CF-1.6 spec chapter 4.2'])

        # sort by criteria passed
        final_lons = sorted(lon_vars, key=lambda x: lon_vars[x], reverse=True)

        obs_mins = {var._name: np.nanmin(var) for var in final_lons if not np.isnan(var).all()}
        obs_maxs = {var._name: np.nanmax(var) for var in final_lons if not np.isnan(var).all()}

        min_pass = any((np.isclose(lon_min, min_val) for min_val in obs_mins.values()))
        max_pass = any((np.isclose(lon_max, max_val) for max_val in obs_maxs.values()))

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

    def check_geospatial_bounds(self, ds):
        """Checks that the geospatial bounds is well formed OGC WKT"""
        var = getattr(ds, 'geospatial_bounds', None)
        check = var is not None
        if not check:
            return Result(BaseCheck.MEDIUM, False,
                          'geospatial_bounds_valid_wkt',
                          ["Attr geospatial_bounds not present"])

        try:
            from_wkt(ds.geospatial_bounds)
        except AttributeError:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'geospatial_bounds_valid_wkt',
                          ['Could not parse WKT, possible bad value for WKT'])
        # parsed OK
        else:
            return Result(BaseCheck.MEDIUM, True, 'geospatial_bounds_valid_wkt',
                          ())


    def check_vertical_extents(self, ds):
        """
        Check that the values of geospatial_vertical_min/geospatial_vertical_max approximately match the data.
        """
        if not (hasattr(ds, 'geospatial_vertical_min') and hasattr(ds, 'geospatial_vertical_max')):
            return

        vert_min = ds.geospatial_vertical_min
        vert_max = ds.geospatial_vertical_max

        # identify vertical vars as per CF 4.3
        v_vars = [var for name, var in ds.variables.items() if is_vertical_coordinate(name, var)]

        if len(v_vars) == 0:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'geospatial_vertical_extents_match',
                          ['Could not find vertical variable to test extent of geospatial_vertical_min/geospatial_vertical_max, see CF-1.6 spec chapter 4.3'])

        obs_mins = {var._name: np.nanmin(var) for var in v_vars if not np.isnan(var).all()}
        obs_maxs = {var._name: np.nanmax(var) for var in v_vars if not np.isnan(var).all()}

        min_pass = any((np.isclose(vert_min, min_val) for min_val in obs_mins.values()))
        max_pass = any((np.isclose(vert_max, max_val) for max_val in obs_maxs.values()))

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
        if not (hasattr(ds, 'time_coverage_start') and hasattr(ds, 'time_coverage_end')):
            return

        epoch = parse_dt("1970-01-01 00:00:00 UTC")
        t_min = (parse_dt(ds.time_coverage_start) - epoch).total_seconds()
        t_max = (parse_dt(ds.time_coverage_end) - epoch).total_seconds()

        # identify t vars as per CF 4.4
        t_vars = [var for name, var in ds.variables.items() if is_time_variable(name, var)]

        if len(t_vars) == 0:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'time_coverage_extents_match',
                          ['Could not find time variable to test extent of time_coverage_start/time_coverage_end, see CF-1.6 spec chapter 4.4'])

        obs_mins = {var._name: Unit(str(var.units)).convert(np.nanmin(var), "seconds since 1970-01-01") for var in t_vars}
        obs_maxs = {var._name: Unit(str(var.units)).convert(np.nanmax(var), "seconds since 1970-01-01") for var in t_vars}

        min_pass = any((np.isclose(t_min, min_val) for min_val in obs_mins.values()))
        max_pass = any((np.isclose(t_max, max_val) for max_val in obs_maxs.values()))

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
    pass
