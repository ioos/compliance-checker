'''
Checks for the Attribute Conventions for Dataset Discovery (ACDD)

This module contains classes defined as checks part of the compliance checker
project for the verification and scoring of attributes for datasets.
'''

from __future__ import unicode_literals

import numpy as np
import numpy.ma as ma
from netCDF4 import num2date
from datetime import timedelta

from compliance_checker.base import (BaseCheck, BaseNCCheck, check_has,
                                     score_group, Result, ratable_result)
from compliance_checker.cf.util import is_vertical_coordinate, _possiblexunits, _possibleyunits
from compliance_checker.util import datetime_is_iso, dateparse
from compliance_checker import cfutil
from pygeoif import from_wkt


class ACDDBaseCheck(BaseCheck):

    _cc_spec = 'acdd'
    _cc_description = 'Attribute Conventions for Dataset Discovery (ACDD)'
    _cc_url = 'http://wiki.esipfed.org/index.php?title=Category:Attribute_Conventions_Dataset_Discovery'

    def __init__(self):

        self.high_rec_atts = [
            'title',
            'keywords',
            'summary'
        ]

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

        # This variable is used to cache the results of applicable variables so
        # the method isn't executed repeatedly.
        self._applicable_variables = None

    # set up attributes according to version
    @check_has(BaseCheck.HIGH)
    def check_high(self, ds):
        '''
        Performs a check on each highly recommended attributes' existence in the dataset

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        return self.high_rec_atts

    @check_has(BaseCheck.MEDIUM)
    def check_recommended(self, ds):
        '''
        Performs a check on each recommended attributes' existence in the dataset

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        return self.rec_atts

    @check_has(BaseCheck.LOW)
    def check_suggested(self, ds):
        return self.sug_atts

    def get_applicable_variables(self, ds):
        '''
        Returns a list of variable names that are applicable to ACDD Metadata
        Checks for variables. This includes geophysical and coordinate
        variables only.

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        if self._applicable_variables is None:
            self.applicable_variables = cfutil.get_geophysical_variables(ds)
            varname = cfutil.get_time_variable(ds)
            if varname:
                self.applicable_variables.append(varname)
            varname = cfutil.get_lon_variable(ds)
            if varname:
                self.applicable_variables.append(varname)
            varname = cfutil.get_lat_variable(ds)
            if varname:
                self.applicable_variables.append(varname)
            varname = cfutil.get_z_variable(ds)
            if varname:
                self.applicable_variables.append(varname)
        return self.applicable_variables

    @score_group('varattr')
    def check_var_long_name(self, ds):
        '''
        Checks each applicable variable for the long_name attribute

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        results = []

        # ACDD Variable Metadata applies to all coordinate variables and
        # geophysical variables only.

        for variable in self.get_applicable_variables(ds):
            msgs = []
            long_name = getattr(ds.variables[variable], 'long_name', None)
            check = long_name is not None
            if not check:
                msgs.append("Var %s missing attr long_name" % variable)
            results.append(Result(BaseCheck.HIGH, check, (variable, "var_std_name"), msgs))

        return results

    @score_group('varattr')
    def check_var_standard_name(self, ds):
        '''
        Checks each applicable variable for the standard_name attribute

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        results = []
        for variable in self.get_applicable_variables(ds):
            msgs = []
            std_name = getattr(ds.variables[variable], 'standard_name', None)
            check = std_name is not None
            if not check:
                msgs.append("Var %s missing attr standard_name" % variable)
            results.append(Result(BaseCheck.HIGH, check, (variable, "var_std_name"), msgs))

        return results

    @score_group('varattr')
    def check_var_units(self, ds):
        '''
        Checks each applicable variable for the units attribute

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        results = []
        for variable in self.get_applicable_variables(ds):
            msgs = []
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

    def check_acknowledgment(self, ds):
        '''
        Check if acknowledgment/acknowledgment attr is present.

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        check = False
        messages = []
        if hasattr(ds, 'acknowledgment') or hasattr(ds, 'acknowledgement'):
            check = True
        else:
            messages.append("acknowledgement global attribute is recommended")

        return Result(BaseCheck.MEDIUM, check, 'acknowledgment/acknowledgement', msgs=messages)

    def check_lat_extents(self, ds):
        '''
        Check that the values of geospatial_lat_min/geospatial_lat_max
        approximately match the data.

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
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
        '''
        Check that the values of geospatial_lon_min/geospatial_lon_max
        approximately match the data.

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
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
            return ratable_result(True, 'geospatial_bounds', tuple())

    def _check_total_z_extents(self, ds, z_variable):
        '''
        Check the entire array of Z for minimum and maximum and compare that to
        the vertical extents defined in the global attributes

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str z_variable: Name of the variable representing the Z-Axis
        '''
        vert_min = ds.geospatial_vertical_min
        vert_max = ds.geospatial_vertical_max
        msgs = []
        total = 2

        zvalue = ds.variables[z_variable][:]
        # If the array has fill values, which is allowed in the case of point
        # features
        if hasattr(zvalue, 'mask'):
            zvalue = zvalue[~zvalue.mask]
        zmin = zvalue.min()
        zmax = zvalue.max()
        if not np.isclose(vert_min, zmin):
            msgs.append("geospatial_vertical_min != min(%s) values, %s != %s" % (
                z_variable,
                vert_min,
                zmin
            ))
        if not np.isclose(vert_max, zmax):
            msgs.append("geospatial_vertical_max != max(%s) values, %s != %s" % (
                z_variable,
                vert_min,
                zmax
            ))

        return Result(BaseCheck.MEDIUM,
                      (total - len(msgs), total),
                      'geospatial_vertical_extents_match',
                      msgs)

    def _check_scalar_vertical_extents(self, ds, z_variable):
        '''
        Check the scalar value of Z compared to the vertical extents which
        should also be equivalent

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str z_variable: Name of the variable representing the Z-Axis
        '''
        vert_min = ds.geospatial_vertical_min
        vert_max = ds.geospatial_vertical_max
        msgs = []
        total = 2

        zvalue = ds.variables[z_variable][:].item()
        if not np.isclose(vert_min, vert_max):
            msgs.append("geospatial_vertical_min != geospatial_vertical_max for scalar depth values, %s != %s" % (
                vert_min,
                vert_max
            ))

        if not np.isclose(vert_max, zvalue):
            msgs.append("geospatial_vertical_max != %s values, %s != %s" % (
                z_variable,
                vert_max,
                zvalue
            ))

        return Result(BaseCheck.MEDIUM,
                      (total - len(msgs), total),
                      'geospatial_vertical_extents_match',
                      msgs)

    def check_vertical_extents(self, ds):
        """
        Check that the values of geospatial_vertical_min/geospatial_vertical_max approximately match the data.

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
        if not (hasattr(ds, 'geospatial_vertical_min') and hasattr(ds, 'geospatial_vertical_max')):
            return

        z_variable = cfutil.get_z_variable(ds)
        if not z_variable:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'geospatial_vertical_extents_match',
                          ['Could not find vertical variable to test extent of geospatial_vertical_min/geospatial_vertical_max, see CF-1.6 spec chapter 4.3'])
        if ds.variables[z_variable].dimensions == tuple():
            return self._check_scalar_vertical_extents(ds, z_variable)

        return self._check_total_z_extents(ds, z_variable)

    def check_time_extents(self, ds):
        """
        Check that the values of time_coverage_start/time_coverage_end approximately match the data.
        """
        if not (hasattr(ds, 'time_coverage_start') and hasattr(ds, 'time_coverage_end')):
            return

        # allows non-ISO 8601 formatted dates
        try:
            t_min = dateparse(ds.time_coverage_start)
            t_max = dateparse(ds.time_coverage_end)
        except:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'time_coverage_extents_match',
                          ['time_coverage variables are not formatted properly. Please ensure they are valid ISO-8601 time strings'])

        timevar = cfutil.get_time_variable(ds)

        if not timevar:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'time_coverage_extents_match',
                          ['Could not find time variable to test extent of time_coverage_start/time_coverage_end, see CF-1.6 spec chapter 4.4'])

        # Time should be monotonically increasing, so we make that assumption here so we don't have to download THE ENTIRE ARRAY
        try:
            time0 = num2date(ds.variables[timevar][0], ds.variables[timevar].units)
            time1 = num2date(ds.variables[timevar][-1], ds.variables[timevar].units)
        except:
            return Result(BaseCheck.MEDIUM,
                          False,
                          'time_coverage_extents_match',
                          ['Failed to retrieve and convert times for variables %s.' % timevar])

        start_dt = abs(time0 - t_min)
        end_dt = abs(time1 - t_max)

        score = 2
        msgs = []
        if start_dt > timedelta(hours=1):
            msgs.append("Date time mismatch between time_coverage_start and actual "
                        "time values %s (time_coverage_start) != %s (time[0])" % (t_min.isoformat(), time0.isoformat()))
            score -= 1
        if end_dt > timedelta(hours=1):
            msgs.append("Date time mismatch between time_coverage_end and actual "
                        "time values %s (time_coverage_end) != %s (time[N])" % (t_max.isoformat(), time1.isoformat()))
            score -= 1

        return Result(BaseCheck.MEDIUM,
                      (score, 2),
                      'time_coverage_extents_match',
                      msgs)

    def verify_convention_version(self, ds):
        """
        Verify that the version in the Conventions field is correct
        """
        for convention in ds.Conventions.replace(' ', '').split(','):
            if convention == 'ACDD-' + self._cc_spec_version:
                return ratable_result((2, 2),
                                      'Conventions',
                                      [])
        # Conventions attribute is present, but does not include
        # proper ACDD version
        messages = [
            "Attr Conventions does not contain 'ACDD-{}'".format(self._cc_spec_version)
        ]
        return ratable_result((1, 2),
                              'Conventions',
                              messages)


class ACDDNCCheck(BaseNCCheck, ACDDBaseCheck):
    pass


class ACDD1_1Check(ACDDNCCheck):

    _cc_spec_version = '1.1'
    _cc_description = 'Attribute Conventions for Dataset Discovery (ACDD) 1.1'
    register_checker = True

    def __init__(self):
        super(ACDD1_1Check, self).__init__()
        self.rec_atts.extend([
            'keywords_vocabulary'
        ])

        self.sug_atts.extend([
            'publisher_name',       # publisher,dataCenter
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
        self.high_rec_atts.extend([
            ('Conventions', self.verify_convention_version)
        ])

        self.rec_atts.extend(['geospatial_vertical_positive',
                              'geospatial_bounds_crs',
                              'geospatial_bounds_vertical_crs',
                              'publisher_name',       # publisher,dataCenter
                              'publisher_url',        # publisher
                              'publisher_email',      # publisher
                              'source'])

        self.sug_atts.extend([
            ('creator_type', [
                'person',
                'group',
                'institution',
                'position'
            ]),
            'creator_institution',
            'platform',
            'platform_vocabulary',
            'keywords_vocabulary',
            'instrument',
            'metadata_link',
            'product_version',
            'references',
            ('publisher_type', [
                'person',
                'group',
                'institution',
                'position'
            ]),
            'instrument_vocabulary',
            'date_metadata_modified',
            'program',
            'publisher_institution',
        ])

    def check_metadata_link(self, ds):
        '''
        Checks if metadata link is formed in a rational manner

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        if not hasattr(ds, u'metadata_link'):
            return
        msgs = []
        meta_link = getattr(ds, 'metadata_link')
        if 'http' not in meta_link:
            msgs.append('Metadata URL should include http:// or https://')
        valid_link = (len(msgs) == 0)
        return Result(BaseCheck.LOW, valid_link, 'metadata_link_valid', msgs)

    def check_date_modified_is_iso(self, ds):
        '''
        Checks if date modified field is ISO compliant

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        if not hasattr(ds, u'date_modified'):
            return
        date_modified_check, msgs = datetime_is_iso(getattr(ds, u'date_modified'))
        return Result(BaseCheck.MEDIUM, date_modified_check, 'date_modified_is_iso', msgs)

    def check_date_issued_is_iso(self, ds):
        '''
        Checks if date issued field is ISO compliant

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        if not hasattr(ds, u'date_issued'):
            return
        date_issued_check, msgs = datetime_is_iso(getattr(ds, u'date_issued'))
        return Result(BaseCheck.MEDIUM, date_issued_check, 'date_issued_is_iso', msgs)

    def check_date_metadata_modified_is_iso(self, ds):
        '''
        Checks if date metadata modified field is ISO compliant

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        if not hasattr(ds, u'date_metadata_modified'):
            return
        date_metadata_modified_check, msgs = datetime_is_iso(getattr(ds, u'date_metadata_modified'))
        return Result(BaseCheck.MEDIUM, date_metadata_modified_check, 'date_metadata_modified_is_iso', msgs)

    def check_id_has_no_blanks(self, ds):
        '''
        Check if there are blanks in the id field

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        if not hasattr(ds, u'id'):
            return
        if ' ' in getattr(ds, u'id'):
            return Result(BaseCheck.MEDIUM, False, 'no_blanks_in_id',
                          msgs=[u'There should be no blanks in the id field'])
        else:
            return Result(BaseCheck.MEDIUM, True, 'no_blanks_in_id', msgs=[])

    def check_date_created(self, ds):
        '''
        Check if date created is ISO-8601

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        if not hasattr(ds, u'date_created'):
            return
        date_created_check, msgs = datetime_is_iso(getattr(ds, u'date_created'))
        return Result(BaseCheck.MEDIUM, date_created_check,
                      'date_created_is_iso', msgs)

    @score_group('varattr')
    def check_var_coverage_content_type(self, ds):
        '''
        Check coverage content type against valid ISO-19115-1 codes

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        results = []
        for variable in cfutil.get_geophysical_variables(ds):
            msgs = []
            ctype = getattr(ds.variables[variable],
                            'coverage_content_type', None)
            check = ctype is not None
            if not check:
                msgs.append("Var %s missing attr coverage_content_type" %
                            variable)
                results.append(Result(BaseCheck.HIGH, check,
                                      (variable, "coverage_content_type"),
                                      msgs))
                continue

            # ISO 19115-1 codes
            valid_ctypes = {
                'image',
                'thematicClassification',
                'physicalMeasurement',
                'auxiliaryInformation',
                'qualityInformation',
                'referenceInformation',
                'modelResult',
                'coordinate'
            }
            if ctype not in valid_ctypes:
                msgs.append("Var %s does not have a coverage_content_type in %s"
                            % (variable, sorted(valid_ctypes)))

        return results

