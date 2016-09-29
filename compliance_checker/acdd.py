from __future__ import unicode_literals
import itertools
import numpy as np
import numpy.ma as ma

from dateutil.parser import parse as parse_dt
from cf_units import Unit

from compliance_checker.base import (BaseCheck, BaseNCCheck, check_has,
                                     score_group, Result, ratable_result)
from compliance_checker.cf.util import (is_time_variable,
                                        is_vertical_coordinate,
                                       _possiblexunits, _possibleyunits)
from compliance_checker.util import datetime_is_iso
from pygeoif import from_wkt

class ACDDBaseCheck(BaseCheck):

    _cc_spec = 'acdd'
    _cc_description = 'Attribute Conventions for Dataset Discovery (ACDD)'
    _cc_url = 'http://wiki.esipfed.org/index.php?title=Category:Attribute_Conventions_Dataset_Discovery'

    def __init__(self):
        self.high_rec_atts = ['title',
                              'keywords',
                              'summary']

        self.rec_atts = [
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
            # mutliple spellings get tested, move to check fn instead
            #'acknowledgement',
            ('geospatial_bounds', self.verify_geospatial_bounds),
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

        self.sug_atts = [
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

    ###############################################################################
    #
    # HIGHLY RECOMMENDED ATTRIBUTES
    #
    ###############################################################################
    # set up attributes according to version
    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        return self.high_rec_atts
    ###############################################################################
    #
    # RECOMMENDED ATTRIBUTES
    #
    ###############################################################################

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        return self.rec_atts
    ###############################################################################
    #
    # SUGGESTED ATTRIBUTES
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

    ################################################################################
    #
    # HIGHLY RECOMMENDED CHECKS
    #
    ###############################################################################


    def check_acknowledgment(self, ds):
        """Check if acknowledgment/acknowledgment attr is present"""
        if not (hasattr(ds, 'acknowledgment') or
                hasattr(ds, 'acknowledgement')):
            return Result(BaseCheck.MEDIUM, False,
                          'acknowledgment/acknowledgement', msgs = [])
        else:
            return Result(BaseCheck.HIGH, True,
                          'acknowledgment/acknowledgement',
                          msgs=["Neither 'acknowledgment' nor 'acknowledgement' attributes present"])


    ###############################################################################
    #
    # RECOMMENDED CHECKS
    #
    ###############################################################################


    ###############################################################################
    #
    # SUGGESTED
    #
    ###############################################################################


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

    def verify_geospatial_bounds(self, ds):
        """Checks that the geospatial bounds is well formed OGC WKT"""
        var = getattr(ds, 'geospatial_bounds', None)
        check = var is not None
        if not check:
            return ratable_result(False,
                          'geospatial_bounds',
                          ["Attr geospatial_bounds not present"])

        try:
            # TODO: verify that WKT is valid given CRS (defaults to EPSG:4326
            #       in ACDD.
            from_wkt(ds.geospatial_bounds)
        except AttributeError:
            return ratable_result(False,
                          'geospatial_bounds',
                          ['Could not parse WKT, possible bad value for WKT'])
        # parsed OK
        else:
            return ratable_result(True, 'geospatial_bounds',
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
        v_vars = [(var._name, ma.array(var)) for name, var in
                  ds.variables.items() if is_vertical_coordinate(name, var)]

        if len(v_vars) == 0:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'geospatial_vertical_extents_match',
                          ['Could not find vertical variable to test extent of geospatial_vertical_min/geospatial_vertical_max, see CF-1.6 spec chapter 4.3'])

        obs_mins = {var[0]: np.nanmin(var[1]) for var in v_vars if not np.isnan(var[1]).all()}
        obs_maxs = {var[0]: np.nanmax(var[1]) for var in v_vars if not np.isnan(var[1]).all()}

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

        # allows non-ISO 8601 formatted dates
        epoch = parse_dt("1970-01-01 00:00:00 UTC")
        try:
            t_min = (parse_dt(ds.time_coverage_start) - epoch).total_seconds()
            t_max = (parse_dt(ds.time_coverage_end) - epoch).total_seconds()
        except:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'time_coverage_extents_match',
                          ['time_coverage variables are not formatted properly'])
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

    def verify_convention_version(self, ds):
        """
        Verify that the version in the Conventions field is correct
        """
        for convention in ds.Conventions.replace(' ','').split(','):
            if convention == 'ACDD-' + self._cc_spec_version:
                return ratable_result(
                        (2,2),
                        'Conventions',
                        [])
        # Conventions attribute is present, but does not include
        # proper ACDD version
        return ratable_result(
                (1,2),
                'Conventions',
                ["Attr Conventions does not contain 'ACDD-{}'".format(
                            self._cc_spec_version)])


class ACDDNCCheck(BaseNCCheck, ACDDBaseCheck):
    pass

class ACDD1_1Check(ACDDNCCheck):

    _cc_spec_version = '1.1'
    _cc_description = 'Attribute Conventions for Dataset Discovery (ACDD) 1.1'
    register_checker = True

    def __init__(self):
        super(ACDD1_1Check, self).__init__()
        self.rec_atts.extend(['keywords_vocabulary',
                        ('cdm_data_type', ['Grid', 'Image', 'Point',
                                           'Radial', 'Station', 'Swath',
                                           'Trajectory'])])

        self.sug_atts.extend(['publisher_name',       # publisher,dataCenter
                                'publisher_url',        # publisher
                                'publisher_email',      # publisher
                                'geospatial_vertical_positive'
                              ])


class ACDD1_3Check(ACDDNCCheck):

    _cc_spec_version = '1.3'
    _cc_description = 'Attribute Conventions for Dataset Discovery (ACDD) 1.3'
    register_checker = True

    def __init__(self):
        super(ACDD1_3Check, self).__init__()
        self.high_rec_atts.extend([('Conventions',
                                    self.verify_convention_version)])

        self.rec_atts.extend(['geospatial_vertical_positive',
                              'geospatial_bounds_crs',
                              'geospatial_bounds_vertical_crs',
                              'publisher_name',       # publisher,dataCenter
                              'publisher_url',        # publisher
                              'publisher_email',      # publisher
                              'source'])

        self.sug_atts.extend([  # 1.3.1, technically
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

    def check_platform_uses_vocab(self, ds):
        #Checks if platform vocab is in vocab list
        if not hasattr(ds, u'platform'):
            return
        if not hasattr(ds, u'platform_vocabulary'):
            return
        platform_names_good = [platform for platform in getattr(ds, u'platform') if platform in getattr(ds, u'platform_vocabulary')]
        return Result(BaseCheck.LOW, (len(platform_names_good),len(getattr(ds,'platform_uses_vocabulary'))), 'platforms_valid', msgs = [u'Platform_vocabulary not present'])

    def check_instrument_uses_vocab(self, ds):
        #Checks if instrument vocab is in vocab list
        if not hasattr(ds, u'instrument'):
            return
        if not hasattr(ds, u'instrument_vocabulary'):
            return
        instrument_names_good = [instrument for instrument in
                                 getattr(ds, u'instrument') if instrument
                                 in getattr(ds, u'instrument_vocabulary')]
        return Result(BaseCheck.LOW, (len(instrument_names_good),
                      len(getattr(ds,'instrument_uses_vocabulary'))),
                      'instruments_valid',
                      msgs=[u'Instrument_vocabulary not present'])

    def check_metadata_link(self, ds):
        #Checks if metadata link is formed in a rational manner
        if not hasattr(ds, u'metadata_link'):
            return
        msgs = []
        meta_link = getattr(ds, 'metadata_link')
        if not 'http' in meta_link:
            msgs.append('Metadata URL should include http:// or https://')
        if not '.' in meta_link:
            msgs.append('Metadata URL is malformed')
        valid_link = len(msgs) == 0
        return Result(BaseCheck.LOW, valid_link,  'metadata_link_valid', msgs)

    def check_date_modified_is_iso(self, ds):
        #Checks if date modified field is ISO compliant
        if not hasattr(ds, u'date_modified'):
            return
        date_modified_check, msgs = datetime_is_iso(getattr(ds, u'date_modified'))
        return Result(BaseCheck.MEDIUM, date_modified_check, 'date_modified_is_iso', msgs)

    def check_date_issued_is_iso(self, ds):
        #Checks if date issued field is ISO compliant
        if not hasattr(ds, u'date_issued'):
            return
        date_issued_check, msgs = datetime_is_iso(getattr(ds, u'date_issued'))
        return Result(BaseCheck.MEDIUM, date_issued_check, 'date_issued_is_iso', msgs)

    def check_date_metadata_modified_is_iso(self, ds):
        #Checks if date metadata modified field is ISO compliant
        if not hasattr(ds, u'date_metadata_modified'):
            return
        date_metadata_modified_check, msgs = datetime_is_iso(getattr(ds, u'date_metadata_modified'))
        return Result(BaseCheck.MEDIUM, date_metadata_modified_check, 'date_metadata_modified_is_iso', msgs)

    def check_id_has_no_blanks(self, ds):
        #Check if there are blanks in the id field
        if not hasattr(ds, u'id'):
            return
        if ' ' in getattr(ds, u'id'):
            return Result(BaseCheck.MEDIUM, False, 'no_blanks_in_id', msgs = [u'There should be no blanks in the id field'])
        else:
            return Result(BaseCheck.MEDIUM, True, 'no_blanks_in_id', msgs = [])

    # BWA: can describe license in plain text
    #def check_license(self, ds):
    #    #Checks if license is from accepted list
    #    if not hasattr(ds, u'license'):
    #        return
    #    license_list = {'none', 'freely distributed'}
    #    if getattr(ds, u'license').lower() in license_list:
    #        return Result(BaseCheck.MEDIUM, True, 'valid_license', msgs=['The license is valid'])
    #    elif '.' in getattr(ds, u'license'):
    #        return Result(BaseCheck.MEDIUM, True, 'valid_license', msgs=['The license is a url'])
    #    else:
    #        return Result(BaseCheck.MEDIUM, False, 'valid_license', msgs=['The license is not a url or in the accepted list'])


    def check_date_created(self, ds):
        #Check if date created is ISO
        if not hasattr(ds, u'date_created'):
            return
        date_created_check, msgs = datetime_is_iso(getattr(ds, u'date_created'))
        return Result(BaseCheck.MEDIUM, date_created_check,
                      'date_created_is_iso', msgs)

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

