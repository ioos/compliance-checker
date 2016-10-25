#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import re
from functools import wraps
from collections import defaultdict
import numpy as np
import os
import six

from compliance_checker.base import BaseCheck, BaseNCCheck, score_group, Result, TestCtx
from compliance_checker.cf.appendix_d import dimless_vertical_coordinates
from compliance_checker.cf.util import NCGraph, StandardNameTable, units_known, units_convertible, units_temporal, map_axes, find_coord_vars, is_time_variable, is_vertical_coordinate, create_cached_data_dir, download_cf_standard_name_table, _possiblet, _possiblez, _possiblex, _possibley, _possibleaxis, _possibleaxisunits
from compliance_checker import cfutil
from compliance_checker.cf.util import _possibleyunits, _possiblexunits

try:
    basestring
except NameError:
    basestring = str


def print_exceptions(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception as e:
            from traceback import print_exc
            print_exc()
    return wrapper

__stdname_table__ = "v29"


def guess_dim_type(dimension):
    """
    Guesses the type of dimension of a variable X/Y/Z/T

    If can't figure it out, None is returned.
    """

    dimclasses = {'T': _possiblet,
                  'Z': _possiblez,
                  'Y': _possibley,
                  'X': _possiblex}

    for dcname, dcvals in dimclasses.items():
        if dimension in dcvals:
            return dcname

    return None


def is_variable(name, var):
    dims = var.dimensions
    if (name,) == dims:
        # Coordinate Type
        return False
    # Probably a variable
    return True


# helper to see if we should do DSG tests
def is_likely_dsg(func):
    @wraps(func)
    def _dec(s, ds):
        if hasattr(ds, 'featureType'):
            return func(s, ds)

        # @TODO: skips if we have formalized skips
        return None

    return _dec


class CFBaseCheck(BaseCheck):
    register_checker = True
    _cc_spec = 'cf'
    # TODO: break out into subclasses once CF-1.7 is a working standard
    _cc_spec_version = '1.6'
    _cc_description = 'Climate and Forecast Conventions (CF)'
    _cc_url = 'http://cfconventions.org'

    @classmethod
    def beliefs(cls):  # @TODO
        return {}
    """
    CF Convention Checker (1.6)

    These checks are translated documents:
        http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.html
        http://cf-pcmdi.llnl.gov/conformance/requirements-and-recommendations/1.6/
    """

    def __init__(self):
        # The compliance checker can be run on multiple datasets in a single
        # instantiation, so caching values has be done by the unique identifier
        # for each dataset loaded.

        # Each default dict is a key, value mapping from the dataset object to
        # a list of variables
        self._coord_vars       = defaultdict(list)
        self._ancillary_vars   = defaultdict(list)
        self._clim_vars        = defaultdict(list)
        self._metadata_vars    = defaultdict(list)
        self._boundary_vars    = defaultdict(list)
        self._geophysical_vars = defaultdict(list)
        self._aux_coords       = defaultdict(list)

        self._std_names        = StandardNameTable()

    ################################################################################
    #
    # Helper Methods - var classifications, etc
    #
    ################################################################################

    def setup(self, ds):
        self._find_coord_vars(ds)
        self._find_aux_coord_vars(ds)
        self._find_ancillary_vars(ds)
        self._find_clim_vars(ds)
        self._find_boundary_vars(ds)
        self._find_metadata_vars(ds)
        self._find_cf_standard_name_table(ds)
        self._find_geophysical_vars(ds)

    def _find_cf_standard_name_table(self, ds):
        '''
        Parse out the `standard_name_vocabulary` attribute and download that
        version of the cf standard name table

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        # Get the standard name vocab
        standard_name_vocabulary = getattr(ds, 'standard_name_vocabulary', '')

        # Try to parse this attribute to get version
        version = None
        if 'cf standard name table' in standard_name_vocabulary.lower():
            version = standard_name_vocabulary.split()[-1]
        else:
            # Can't parse the attribute, use the packaged version
            return 0

        if version.startswith('v'):  # i.e 'v34' -> '34' drop the v
            version = version[1:]

        # If the packaged version is what we're after, then we're good
        if version == self._std_names._version:
            print("Using packaged standard name table v{0}".format(version))
            return 0

        # Try to download the version specified
        try:
            data_directory = create_cached_data_dir()
            location = os.path.join(data_directory, 'cf-standard-name-table-test-{0}.xml'.format(version))
            # Did we already download this before?
            if not os.path.isfile(location):
                download_cf_standard_name_table(version, location)
                print("Using downloaded standard name table v{0}".format(version))
            else:
                print("Using cached standard name table v{0} from {1}".format(version, location))

            self._std_names = StandardNameTable(location)
            return 1
        except Exception:
            # There was an error downloading the CF table. That's ok, we'll just use the packaged version
            print("Error fetching standard name table. Using packaged v{0}".format(self._std_names._version))
            return 0

    def _find_coord_vars(self, ds, refresh=False):
        '''
        Returns a list of variable names that identify as coordinate variables.

        The result is cached by the passed in dataset object inside of this
        checker. Pass refresh=True to redo the cached value.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        :rtype: list
        :return: A list of variables names (str) that are defined as coordinate
                 variables in the dataset ds.
        '''
        if ds in self._coord_vars and refresh is False:
            return self._coord_vars[ds]

        self._coord_vars[ds] = cfutil.get_coordinate_variables(ds)

        return self._coord_vars[ds]

    def _find_aux_coord_vars(self, ds, refresh=False):
        '''
        Returns a list of auxiliary coordinate variables

        An auxiliary coordinate variable is any netCDF variable that contains
        coordinate data, but is not a coordinate variable (in the sense of the term
        defined by CF).

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        :rtype: list
        :return: List of variable names (str) that are defined to be auxiliary
                 coordinate variables.
        '''
        if self._aux_coords.get(ds, None) and refresh is False:
            return self._aux_coords[ds]

        self._aux_coords[ds] = cfutil.get_auxiliary_coordinate_variables(ds)
        return self._aux_coords[ds]

    def _find_ancillary_vars(self, ds, refresh=False):
        '''
        Returns a list of variable names that are defined as ancillary
        variables in the dataset ds.

        An ancillary variable generally is a metadata container and referenced
        from other variables via a string reference in an attribute.

        - via ancillary_variables (3.4)
        - "grid mapping var" (5.6)
        - TODO: more?

        The result is cached by the passed in dataset object inside of this
        checker. Pass refresh=True to redo the cached value.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        :rtype: list
        :return: List of variable names (str) that are defined as ancillary
                 variables in the dataset ds.
        '''

        # Used the cached version if it exists and is not empty
        if self._ancillary_vars.get(ds, None) and refresh is False:
            return self._ancillary_vars[ds]

        # Invalidate the cache at all costs
        self._ancillary_vars[ds] = []

        for name, var in ds.variables.items():
            if hasattr(var, 'ancillary_variables'):
                for anc_name in var.ancillary_variables.split(" "):
                    if anc_name in ds.variables:
                        self._ancillary_vars[ds].append(anc_name)

            if hasattr(var, 'grid_mapping'):
                gm_name = var.grid_mapping
                if gm_name in ds.variables:
                    self._ancillary_vars[ds].append(gm_name)

        return self._ancillary_vars[ds]

    def _find_metadata_vars(self, ds, refresh=False):
        '''
        Returns a list of netCDF variable instances for those that are likely metadata variables

        :param netCDF4.Dataset ds: An open netCDF dataset
        '''
        if self._metadata_vars.get(ds, None) and refresh is False:
            return self._metadata_vars[ds]

        self._metadata_vars[ds] = []
        for name, var in ds.variables.items():

            if name in self._find_ancillary_vars(ds) or name in self._find_coord_vars(ds):
                continue

            if name in ('platform_name', 'station_name', 'instrument_name', 'station_id', 'platform_id', 'surface_altitude'):
                self._metadata_vars[ds].append(name)

            elif getattr(var, 'cf_role', '') != '':
                self._metadata_vars[ds].append(name)

            elif getattr(var, 'standard_name', None) is None and len(var.dimensions) == 0:
                self._metadata_vars[ds].append(name)

        return self._metadata_vars[ds]

    def _find_geophysical_vars(self, ds, refresh=False):
        '''
        Returns a list of geophysical variables

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: dict
        '''
        if self._geophysical_vars.get(ds, None) and refresh is False:
            return self._geophysical_vars[ds]

        self._geophysical_vars[ds] = cfutil.get_geophysical_variables(ds)

        return self._geophysical_vars[ds]

    def _find_clim_vars(self, ds, refresh=False):
        '''
        Returns a list of variables that are likely to be climatology variables based on CF §7.4

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        '''

        if self._clim_vars.get(ds, None) and refresh is False:
            return self._clim_vars[ds]

        climatology_variable = cfutil.get_climatology_variable(ds)
        if climatology_variable:
            self._clim_vars[ds].append(climatology_variable)

        return self._clim_vars[ds]

    def _find_boundary_vars(self, ds, refresh=False):
        '''
        Returns dictionary of boundary variables mapping the variable instance
        to the name of the variable acting as a boundary variable.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        '''
        if self._boundary_vars.get(ds, None) and refresh is False:
            return self._boundary_vars[ds]

        self._boundary_vars[ds] = cfutil.get_cell_boundary_variables(ds)

        return self._boundary_vars[ds]

    ###############################################################################
    #
    # Chapter 2: NetCDF Files and Components
    #
    ###############################################################################

    def check_data_types(self, ds):
        '''
        Checks the data type of all netCDF variables to ensure they are valid
        data types under CF.

        CF §2.2 The netCDF data types char, byte, short, int, float or real, and
        double are all acceptable

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: Result
        '''
        fails = []
        total = len(ds.variables)

        for k, v in ds.variables.items():
            if v.dtype not in [np.character,
                               np.dtype('c'),
                               np.dtype('b'),
                               np.dtype('i4'),
                               np.int32,
                               np.float32,
                               np.double,
                               'int16',
                               'float32'
                               ]:

                fails.append('The variable {} failed because the datatype is {}'.format(k, v.datatype))
        return Result(BaseCheck.HIGH, (total - len(fails), total), '§2.2 Valid netCDF data types', msgs=fails)

    def check_naming_conventions(self, ds):
        '''
        Checks the variable names to ensure they are valid CF variable names under CF.

        CF §2.3 Variable, dimension and attribute names should begin with a letter
        and be composed of letters, digits, and underscores.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: Result
        '''
        ret_val = []
        variable_naming = TestCtx(BaseCheck.MEDIUM, '§2.3 Naming Conventions for variables')
        dimension_naming = TestCtx(BaseCheck.MEDIUM, '§2.3 Naming Conventions for dimensions')
        attribute_naming = TestCtx(BaseCheck.MEDIUM, '§2.3 Naming Conventions for attributes')

        rname = re.compile("^[A-Za-z][A-Za-z0-9_]*$")

        for name, variable in ds.variables.items():
            variable_naming.assert_true(rname.match(name) is not None,
                                        "variable {} should begin with a letter and be composed of "
                                        "letters, digits, and underscores".format(name))

            # Keep track of all the attributes, we'll need to check them
            for attr in variable.ncattrs():
                # ignore reserved attrs
                if attr == '_FillValue':
                    continue
                attribute_naming.assert_true(rname.match(attr) is not None,
                                             "attribute {}:{} should begin with a letter and be composed of "
                                             "letters, digits, and underscores".format(name, attr))

        ret_val.append(variable_naming.to_result())

        for dimension in ds.dimensions:
            dimension_naming.assert_true(rname.match(dimension) is not None,
                                         "dimension {} should begin with a latter and be composed of "
                                         "letters, digits, and underscores".format(dimension))
        ret_val.append(dimension_naming.to_result())

        for global_attr in ds.ncattrs():
            attribute_naming.assert_true(rname.match(global_attr) is not None,
                                         "global attribute {} should begin with a letter and be composed of "
                                         "letters, digits, and underscores".format(global_attr))
        ret_val.append(attribute_naming.to_result())

        return ret_val

    def check_names_unique(self, ds):
        '''
        Checks the variable names for uniqueness regardless of case.

        CF §2.3 names should not be distinguished purely by case, i.e., if case
        is disregarded, no two names should be the same.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: Result
        '''
        fails = []
        total = len(ds.variables)
        names = defaultdict(int)

        for k in ds.variables:
            names[k.lower()] += 1

        fails = ['Variables are not case sensitive. Duplicate variables named: %s' % k for k, v in names.items() if v > 1]
        return Result(BaseCheck.MEDIUM, (total - len(fails), total), '§2.3 Unique variable names', msgs=fails)

    def check_dimension_names(self, ds):
        '''
        Checks variables contain no duplicate dimension names.

        CF §2.4 A variable may have any number of dimensions, including zero,
        and the dimensions must all have different names.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: Result
        '''
        fails = []
        total = len(ds.variables)

        for k, v in ds.variables.items():
            dims = defaultdict(int)
            for d in v.dimensions:
                dims[d] += 1

            for dimension, count in dims.items():
                if count > 1:
                    fails.append("%s has two or more dimensions named %s" % (k, dimension))

        return Result(BaseCheck.HIGH, (total - len(fails), total), '§2.4 Unique dimensions', msgs=fails)

    def check_dimension_order(self, ds):
        '''
        Checks each variable's dimension order to ensure that the order is
        consistent and in order under CF §2.4

        CF §2.4 If any or all of the dimensions of a variable have the
        interpretations of "date or time" (T), "height or depth" (Z),
        "latitude" (Y), or "longitude" (X) then we recommend, those dimensions
        to appear in the relative order T, then Z, then Y, then X in the CDL
        definition corresponding to the file. All other dimensions should,
        whenever possible, be placed to the left of the spatiotemporal
        dimensions.
        '''
        valid_dimension_order = TestCtx(BaseCheck.MEDIUM, '§2.4 Dimension Order')
        expected = ['T', 'Z', 'Y', 'X']
        coord_vars = self._find_coord_vars(ds)
        # Build a map from coordinate variable to axis
        coord_axis_map = {}

        for coord_name in coord_vars:
            coord_var = ds.variables[coord_name]
            axis = getattr(coord_var, 'axis', None)
            standard_name = getattr(coord_var, 'standard_name', None)
            # axis takes precedence over standard_name
            if axis in expected:
                coord_axis_map[coord_name] = axis
            elif standard_name == 'time':
                coord_axis_map[coord_name] = 'T'
            elif standard_name == 'longitude':
                coord_axis_map[coord_name] = 'X'
            elif standard_name == 'latitude':
                coord_axis_map[coord_name] = 'Y'
            elif standard_name in ['height', 'depth', 'altitude']:
                coord_axis_map[coord_name] = 'Z'
            elif cfutil.is_compression_coordinate(ds, coord_name):
                coord_axis_map[coord_name] = 'C'
            else:
                # mark the coordinate variable as unknown
                coord_axis_map[coord_name] = 'U'

        # If a dimension does not have a coordinate variable mark it as unknown
        # 'U'
        for dimension in ds.dimensions:
            if dimension not in coord_axis_map:
                coord_axis_map[dimension] = 'U'

        # Check each variable's dimension order
        for name, variable in ds.variables.items():
            if variable.dimensions:
                dimension_order = self._get_dimension_order(ds, name, coord_axis_map)
                valid_dimension_order.assert_true(self._dims_in_order(dimension_order),
                                                  "{}'s dimensions are not in the recommended order "
                                                  "T, X, Y, Z. They are {}"
                                                  "".format(name,
                                                            ", ".join(dimension_order)))

        return valid_dimension_order.to_result()

    def _get_dimension_order(self, ds, name, coord_axis_map):
        '''
        Returns a list of strings corresponding to the named axis of the dimensions for a variable.

        Example::
            self._get_dimension_order(ds, 'temperature', coord_axis_map)
            --> ['T', 'Y', 'X']

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Name of the variable
        :param dict coord_axis_map: A dictionary mapping each coordinate variable and dimension to a named axis
        :rtype: list
        '''

        retval = []
        variable = ds.variables[name]
        for dim in variable.dimensions:
            retval.append(coord_axis_map[dim])
        return retval

    def _dims_in_order(self, dimension_order):
        '''
        Returns True if the dimensions are in order U*, T, Z, Y, X

        :param list dimension_order: A list of axes
        '''
        regx = re.compile(r'^U*T?Z?(?:(?:Y?X?)|(?:C?))$')
        dimension_string = ''.join(dimension_order)
        return regx.match(dimension_string) is not None

    def check_fill_value_outside_valid_range(self, ds):
        '''
        Checks each variable's _FillValue to ensure that it's in valid_range or
        between valid_min and valid_max according to CF §2.5.1

        CF §2.5.1 The _FillValue should be outside the range specified by
        valid_range (if used) for a variable.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        valid_fill_range = TestCtx(BaseCheck.MEDIUM,
                                   '§2.5.1 Fill Values should be outside the range specified by valid_range')

        for name, variable in ds.variables.items():
            # If the variable doesn't have a defined _FillValue don't check it.

            if not hasattr(variable, '_FillValue'):
                continue

            fill_value = variable._FillValue

            attrs = variable.ncattrs()

            if 'valid_range' in attrs:
                rmin, rmax = variable.valid_range
                spec_by = 'valid_range'

            elif 'valid_min' in attrs and 'valid_max' in attrs:
                rmin = variable.valid_min
                rmax = variable.valid_max
                spec_by = 'valid_min/valid_max'
            else:
                continue

            valid = (fill_value < rmin or fill_value > rmax)

            valid_fill_range.assert_true(valid,
                                         "{}:_FillValue ({}) should be outside the range specified by {} ({}, {})"
                                         "".format(name, fill_value, spec_by, rmin, rmax))

        return valid_fill_range.to_result()

    def check_conventions_are_cf_16(self, ds):
        '''
        Check the global attribute conventions to contain CF-1.6

        CF §2.6.1 the NUG defined global attribute Conventions to the string
        value "CF-1.6"

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: Result
        '''

        valid_conventions = ['CF-1.6']
        if hasattr(ds, 'Conventions'):
            conventions = re.split(',|\s+', getattr(ds, 'Conventions', ''))
            if any((c.strip() in valid_conventions for c in conventions)):
                valid = True
                reasoning = []
            else:
                valid = False
                reasoning = ['Conventions global attribute does not contain "CF-1.6"']
        else:
            valid = False
            reasoning = ['Conventions field is not present']
        return Result(BaseCheck.MEDIUM, valid, '§2.6.1 Global Attribute Conventions includes CF-1.6', msgs=reasoning)

    def check_convention_globals(self, ds):
        '''
        Check the common global attributes are strings if they exist.

        CF §2.6.2 title/history global attributes, must be strings. Do not need
        to exist.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        attrs = ['title', 'history']

        valid_globals = TestCtx(BaseCheck.MEDIUM, '§2.6.2 Recommended Global Attributes')

        for attr in attrs:
            dataset_attr = getattr(ds, attr, None)
            is_string = isinstance(dataset_attr, basestring)
            valid_globals.assert_true(is_string and len(dataset_attr),
                                      "global attribute {} should exist and be a non-empty string"
                                      "".format(attr))

        return valid_globals.to_result()

    def check_convention_possibly_var_attrs(self, ds):
        """
        Check variable and global attributes are strings for recommended attributes under CF §2.6.2

        CF §2.6.2 institution, source, references, and comment, either global
        or assigned to individual variables.  When an attribute appears both
        globally and as a variable attribute, the variable's version has
        precedence.  Must be strings.
        """
        attrs = ['institution', 'source', 'references', 'comment']

        valid_attributes = TestCtx(BaseCheck.MEDIUM, '§2.6.2 Recommended Attributes')

        attr_bin = set()
        # If the attribute is defined for any variable, check it and mark in
        # the set that we've seen it at least once.
        for name, variable in ds.variables.items():
            for attribute in variable.ncattrs():
                varattr = getattr(variable, attribute)
                if attribute in attrs:
                    is_string = isinstance(varattr, basestring)
                    valid_attributes.assert_true(is_string and len(varattr) > 0,
                                                 "{}:{} should be a non-empty string"
                                                 "".format(name, attribute))
                    attr_bin.add(attribute)

        # Check all the global attributes too and mark if we've seen them
        for attribute in ds.ncattrs():
            dsattr = getattr(ds, attribute)
            if attribute in attrs:
                is_string = isinstance(dsattr, basestring)
                valid_attributes.assert_true(is_string and len(dsattr) > 0,
                                             "{} global attribute should be a non-empty string"
                                             "".format(attribute))
                attr_bin.add(attribute)
        # Make sure we've seen each attribute at least once.

        valid_attributes.assert_true('institution' in attr_bin,
                                     "institution should be defined")
        valid_attributes.assert_true('source' in attr_bin,
                                     "source should be defined")
        valid_attributes.assert_true('references' in attr_bin,
                                     "references should be defined")

        # comment is optional and only needs to be a string and non-empty if it
        # exists.

        return valid_attributes.to_result()

    ###############################################################################
    #
    # Chapter 3: Description of the Data
    #
    ###############################################################################

    def check_units(self, ds):
        '''
        Check the units attribute for all variables to ensure they are CF
        compliant under CF §3.1

        CF §3.1 The units attribute is required for all variables that represent dimensional quantities
        (except for boundary variables defined in Section 7.1, "Cell Boundaries" and climatology variables
        defined in Section 7.4, "Climatological Statistics").

        Units are not required for dimensionless quantities. A variable with no units attribute is assumed
        to be dimensionless. However, a units attribute specifying a dimensionless unit may optionally be
        included.

        - units required
        - type must be recognized by udunits
        - if std name specified, must be consistent with standard name table, must also be consistent with a
          specified cell_methods attribute if present

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []

        coordinate_variables = self._find_coord_vars(ds)
        auxiliary_coordinates = self._find_aux_coord_vars(ds)
        geophysical_variables = self._find_geophysical_vars(ds)
        unit_required_variables = coordinate_variables + auxiliary_coordinates + geophysical_variables

        for name in unit_required_variables:
            # For reduced horizontal grids, the compression index variable does
            # not require units.
            if cfutil.is_compression_coordinate(ds, name):
                continue

            variable = ds.variables[name]

            standard_name = getattr(variable, 'standard_name', None)
            standard_name, standard_name_modifier = self._split_standard_name(standard_name)

            valid_units = self._check_valid_cf_units(ds, name)
            ret_val.append(valid_units)

            valid_udunits = self._check_valid_udunits(ds, name)
            ret_val.append(valid_udunits)

            if standard_name is None:
                continue

            valid_standard_units = self._check_valid_standard_units(ds, name)
            ret_val.append(valid_standard_units)

        return ret_val

    def _split_standard_name(self, standard_name):
        '''
        Returns a tuple of the standard_name and standard_name modifier

        Nones are used to represent the absence of a modifier or standard_name

        :rtype: tuple
        '''
        standard_name_modifier = None
        if not isinstance(standard_name, basestring):
            return (None, None)

        if ' ' in standard_name:
            standard_name, standard_name_modifier = standard_name.split(' ', 1)

        return (standard_name, standard_name_modifier)

    def _check_valid_cf_units(self, ds, variable_name):
        '''
        Checks that the variable contains units attribute, the attribute is a
        string and the value is not deprecated by CF

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str variable_name: Name of the variable to be checked
        '''
        # This list is straight from section 3
        deprecated = ['level', 'layer', 'sigma_level']
        variable = ds.variables[variable_name]

        units = getattr(variable, 'units', None)
        standard_name = getattr(variable, 'standard_name', None)
        standard_name, standard_name_modifier = self._split_standard_name(standard_name)
        unitless_standard_names = cfutil.get_unitless_standard_names()

        should_be_unitless = standard_name in unitless_standard_names

        # 1) Units must exist
        valid_units = TestCtx(BaseCheck.HIGH, '§3.1 Variable {} contains valid CF units'.format(variable_name))
        valid_units.assert_true(should_be_unitless or units is not None,
                                'units attribute is required for {}'.format(variable_name))

        # 2) units attribute must be a string
        valid_units.assert_true(should_be_unitless or isinstance(units, basestring),
                                'units attribute for {} needs to be a string'.format(variable_name))

        # 3) units are not deprecated
        valid_units.assert_true(units not in deprecated,
                                'units for {}, "{}" are deprecated by CF 1.6'.format(variable_name, units))
        return valid_units.to_result()

    def _check_valid_udunits(self, ds, variable_name):
        '''
        Checks that the variable's units are contained in UDUnits

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str variable_name: Name of the variable to be checked
        '''
        variable = ds.variables[variable_name]

        units = getattr(variable, 'units', None)
        standard_name = getattr(variable, 'standard_name', None)
        standard_name, standard_name_modifier = self._split_standard_name(standard_name)
        unitless_standard_names = cfutil.get_unitless_standard_names()

        # If the variable is supposed to be unitless, it automatically passes
        should_be_unitless = standard_name in unitless_standard_names

        valid_udunits = TestCtx(BaseCheck.LOW,
                                "§3.1 Variable {}'s units are contained in UDUnits".format(variable_name))
        are_udunits = (units is not None and units_known(units))
        valid_udunits.assert_true(should_be_unitless or are_udunits,
                                  'units for {}, "{}" are not recognized by udunits'.format(variable_name, units))
        return valid_udunits.to_result()

    def _check_valid_standard_units(self, ds, variable_name):
        '''
        Checks that the variable's units are appropriate for the standard name
        according to the CF standard name table and coordinate sections in CF
        1.6

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str variable_name: Name of the variable to be checked
        '''
        variable = ds.variables[variable_name]
        units = getattr(variable, 'units', None)
        standard_name = getattr(variable, 'standard_name', None)

        unitless_standard_names = cfutil.get_unitless_standard_names()

        # If the variable is supposed to be unitless, it automatically passes
        should_be_unitless = standard_name in unitless_standard_names

        valid_standard_units = TestCtx(BaseCheck.HIGH,
                                       "§3.1 Variable {}'s units are appropriate for "
                                       "the standard_name {}".format(variable_name,
                                                                     standard_name or "unspecified"))

        standard_name, standard_name_modifier = self._split_standard_name(standard_name)

        standard_entry = self._std_names.get(standard_name, None)
        if standard_entry is not None:
            canonical_units = standard_entry.canonical_units
        else:
            canonical_units = None

        # Other standard_name modifiers have the same units as the
        # unmodified standard name or are not checked for units.

        if standard_name_modifier == 'number_of_observations':
            canonical_units = '1'

        # This section represents the different cases where simple udunits
        # comparison isn't comprehensive enough to determine if the units are
        # appropriate under CF

        # UDUnits accepts "s" as a unit of time but it should be <unit> since <epoch>
        if standard_name == 'time':
            valid_standard_units.assert_true(units_convertible(units, 'seconds since 1970-01-01'),
                                             'time must be in a valid units format <unit> since <epoch> '
                                             'not {}'.format(units))

        # UDunits can't tell the difference between east and north facing coordinates
        elif standard_name == 'latitude':
            # degrees is allowed if using a transformed grid
            allowed_units = [i.lower() for i in _possibleyunits] + ['degrees']
            valid_standard_units.assert_true(units.lower() in allowed_units,
                                             'variables defining latitude must use degrees_north '
                                             'or degrees if defining a transformed grid. Currently '
                                             '{}'.format(units))
        # UDunits can't tell the difference between east and north facing coordinates
        elif standard_name == 'longitude':
            # degrees is allowed if using a transformed grid
            allowed_units = [i.lower() for i in _possiblexunits] + ['degrees']
            valid_standard_units.assert_true(units.lower() in allowed_units,
                                             'variables defining longitude must use degrees_east '
                                             'or degrees if defining a transformed grid. Currently '
                                             '{}'.format(units))
        # Standard Name table agrees the unit should be unitless
        elif should_be_unitless:
            valid_standard_units.assert_true(True, '')

        else:
            valid_standard_units.assert_true(units_convertible(canonical_units, units),
                                             'units for variable {} must be convertible to {} '
                                             'currently they are {}'.format(variable_name, canonical_units, units))

        return valid_standard_units.to_result()

    def check_standard_name(self, ds):
        '''
        Check a variables's standard_name attribute to ensure that it meets CF
        compliance.

        CF §3.3 A standard name is associated with a variable via the attribute
        standard_name which takes a string value comprised of a standard name
        optionally followed by one or more blanks and a standard name modifier

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []

        coord_vars = self._find_coord_vars(ds)
        aux_coord_vars = self._find_aux_coord_vars(ds)
        axis_vars = cfutil.get_axis_variables(ds)
        flag_vars = cfutil.get_flag_variables(ds)
        geophysical_vars = self._find_geophysical_vars(ds)

        variables_requiring_standard_names = coord_vars + aux_coord_vars + axis_vars + flag_vars + geophysical_vars
        for name in variables_requiring_standard_names:
            # Compression indices used in reduced horizontal grids or
            # compression schemes do not require attributes other than compress
            if cfutil.is_compression_coordinate(ds, name):
                continue

            ncvar = ds.variables[name]
            standard_name = getattr(ncvar, 'standard_name', None)

            standard_name, standard_name_modifier = self._split_standard_name(standard_name)
            # §1.3 The long_name and standard_name attributes are used to
            # describe the content of each variable. For backwards
            # compatibility with COARDS neither is required, but use of at
            # least one of them is strongly recommended.

            # If standard_name is not defined but long_name is, don't continue
            # the check for this variable
            if standard_name is None:
                long_name = getattr(ncvar, 'long_name', None)
                if long_name is not None:
                    continue

            valid_std_name = TestCtx(BaseCheck.HIGH, '§3.3 Variable {} has valid standard_name attribute'.format(name))

            valid_std_name.assert_true(isinstance(standard_name, basestring),
                                       "variable {}'s attribute standard_name must be a string".format(name))

            valid_std_name.assert_true(standard_name in self._std_names,
                                       "standard_name {} is not defined in Standard Name Table v{}".format(
                                           standard_name or 'undefined',
                                           self._std_names._version))

            ret_val.append(valid_std_name.to_result())

            # 2) optional - if modifiers, should be in table
            if standard_name_modifier is not None:
                valid_modifier = TestCtx(BaseCheck.HIGH, "§3.3 standard_name modifier for {} is valid".format(name))
                allowed = ['detection_minimum',
                           'number_of_observations',
                           'standard_error',
                           'status_flag']
                valid_modifier.assert_true(standard_name_modifier in allowed,
                                           "standard_name modifier {} is not a valid modifier "
                                           "according to appendix C".format(standard_name_modifier))

                ret_val.append(valid_modifier.to_result())

        return ret_val

    def check_ancillary_variables(self, ds):
        '''
        Checks the ancillary_variable attribute for all variables to ensure
        they are CF compliant.

        CF §3.4 It is a string attribute whose value is a blank separated list
        of variable names.  The nature of the relationship between variables
        associated via ancillary_variables must be determined by other
        attributes. The variables listed by the ancillary_variables attribute
        will often have the standard name of the variable which points to them
        including a modifier (Appendix C, Standard Name Modifiers) to indicate
        the relationship.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []

        for ncvar in ds.get_variables_by_attributes(ancillary_variables=lambda x: x is not None):
            name = ncvar.name
            valid_ancillary = TestCtx(BaseCheck.HIGH, "§3.4 Ancillary Variables defined by {}".format(name))
            ancillary_variables = ncvar.ancillary_variables

            valid_ancillary.assert_true(isinstance(ancillary_variables, basestring),
                                        "ancillary_variables attribute defined by {} "
                                        "should be string".format(name))

            # Can't perform the second check if it's not a string
            if not isinstance(ancillary_variables, basestring):
                ret_val.append(valid_ancillary.to_result())
                continue

            for ancillary_variable in ancillary_variables.split():
                valid_ancillary.assert_true(ancillary_variables in ds.variables,
                                            "{} is not a variable in this dataset".format(ancillary_variable))

            ret_val.append(valid_ancillary.to_result())

        return ret_val

    def check_flags(self, ds):
        '''
        Check the flag_values, flag_masks and flag_meanings attributes for
        variables to ensure they are CF compliant.

        CF §3.5 The attributes flag_values, flag_masks and flag_meanings are
        intended to make variables that contain flag values self describing.
        Status codes and Boolean (binary) condition flags may be expressed with
        different combinations of flag_values and flag_masks attribute
        definitions.

        The flag_values and flag_meanings attributes describe a status flag
        consisting of mutually exclusive coded values.

        The flag_meanings attribute is a string whose value is a blank
        separated list of descriptive words or phrases, one for each flag
        value. Each word or phrase should consist of characters from the
        alphanumeric set and the following five: '_', '-', '.', '+', '@'.

        The flag_masks and flag_meanings attributes describe a number of
        independent Boolean conditions using bit field notation by setting
        unique bits in each flag_masks value.

        The flag_masks, flag_values and flag_meanings attributes, used
        together, describe a blend of independent Boolean conditions and
        enumerated status codes. A flagged condition is identified by a bitwise
        AND of the variable value and each flag_masks value; a result that
        matches the flag_values value indicates a true condition.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []

        for name in cfutil.get_flag_variables(ds):
            variable = ds.variables[name]
            flag_values = getattr(variable, "flag_values", None)
            flag_masks = getattr(variable, "flag_masks", None)

            valid_flags_var = TestCtx(BaseCheck.HIGH, '§3.5 {} is a valid flags variable'.format(name))
            # Check that the variable defines mask or values
            valid_flags_var.assert_true(flag_values is not None or flag_masks is not None,
                                        "{} does not define either flag_masks or flag_values".format(name))
            ret_val.append(valid_flags_var.to_result())

            valid_meanings = self._check_flag_meanings(ds, name)
            ret_val.append(valid_meanings)

            # check flag_values
            if flag_values is not None:
                valid_values = self._check_flag_values(ds, name)
                ret_val.append(valid_values)

            # check flag_masks
            if flag_masks is not None:
                valid_masks = self._check_flag_masks(ds, name)
                ret_val.append(valid_masks)

            if flag_values is not None and flag_masks is not None:
                allv = list(map(lambda a, b: a & b == a, list(zip(flag_values, flag_masks))))

                allvr = Result(BaseCheck.MEDIUM, all(allv), '§3.5 flags for {}'.format(name))
                if not allvr.value:
                    allvr.msgs = ["flag masks and flag values combined don't equal flag value"]

                ret_val.append(allvr)

        return ret_val

    def _check_flag_values(self, ds, name):
        '''
        Checks a variable's flag_values attribute for compliance under CF

        - flag_values exists as an array
        - unique elements in flag_values
        - flag_values si the same dtype as the variable
        - flag_values is the same length as flag_meanings

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Name of variable to check
        '''
        variable = ds.variables[name]

        flag_values = variable.flag_values
        flag_meanings = getattr(variable, 'flag_meanings', None)
        valid_values = TestCtx(BaseCheck.HIGH, '§3.5 flag_values for {}'.format(name))

        # flag_values must be a list of values, not a string or anything else
        valid_values.assert_true(isinstance(flag_values, np.ndarray),
                                 "flag_values must be an array of values not {}".format(type(flag_values)))

        # We can't perform any more checks
        if not isinstance(flag_values, np.ndarray):
            return valid_values.to_result()

        # the flag values must be independent, no repeating values
        flag_set = set(flag_values)
        valid_values.assert_true(len(flag_set) == len(flag_values),
                                 "flag_values must be independent and can not be repeated")

        # the data type for flag_values should be the same as the variable
        valid_values.assert_true(variable.dtype == flag_values.dtype,
                                 "flag_values ({}) must be the same data type as {} ({})"
                                 "".format(flag_values.dtype, name, variable.dtype))

        if isinstance(flag_meanings, basestring):
            flag_meanings = flag_meanings.split()
            valid_values.assert_true(len(flag_meanings) == len(flag_values),
                                     "flag_meanings and flag_values should have the same number "
                                     "of elements.")

        return valid_values.to_result()

    def _check_flag_masks(self, ds, name):
        '''
        Check a variable's flag_masks attribute for compliance under CF

        - flag_masks exists as an array
        - flag_masks is the same dtype as the variable
        - variable's dtype can support bit-field
        - flag_masks is the same length as flag_meanings

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Variable name
        '''
        variable = ds.variables[name]

        flag_masks = variable.flag_masks
        flag_meanings = getattr(ds, 'flag_meanings', None)

        valid_masks = TestCtx(BaseCheck.HIGH, '§3.5 flag_masks for {}'.format(name))

        valid_masks.assert_true(isinstance(flag_masks, np.ndarray),
                                "flag_masks must be an array of values not {}".format(type(flag_masks)))

        if not isinstance(flag_masks, np.ndarray):
            return valid_masks.to_result()

        valid_masks.assert_true(variable.dtype == flag_masks.dtype,
                                "flag_masks ({}) mustbe the same data type as {} ({})"
                                "".format(flag_masks.dtype, name, variable.dtype))

        type_ok = (np.issubdtype(variable.dtype, int) or
                   np.issubdtype(variable.dtype, 'S') or
                   np.issubdtype(variable.dtype, 'b'))

        valid_masks.assert_true(type_ok, "{}'s data type must be capable of bit-field expression")

        if isinstance(flag_meanings, basestring):
            flag_meanings = flag_meanings.split()
            valid_masks.assert_true(len(flag_meanings) == len(flag_masks),
                                    "flag_meanings and flag_masks should have the same number "
                                    "of elements.")

        return valid_masks.to_result()

    def _check_flag_meanings(self, ds, name):
        '''
        Check a variable's flag_meanings attribute for compliance under CF

        - flag_meanings exists
        - flag_meanings is a string
        - flag_meanings elements are valid strings

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Variable name
        '''
        variable = ds.variables[name]
        flag_meanings = getattr(variable, 'flag_meanings', None)
        valid_meanings = TestCtx(BaseCheck.HIGH, '§3.5 flag_meanings for {}'.format(name))

        valid_meanings.assert_true(flag_meanings is not None,
                                   "flag_meanings attribute is required for flag variables")

        valid_meanings.assert_true(isinstance(flag_meanings, basestring),
                                   "flag_meanings attribute must be a string")

        # We can't perform any additional checks if it's not a string
        if not isinstance(flag_meanings, basestring):
            return valid_meanings.to_result()

        valid_meanings.assert_true(len(flag_meanings) > 0,
                                   "flag_meanings can't be empty")

        flag_regx = re.compile("^[0-9A-Za-z_\-.+@]+$")
        meanings = flag_meanings.split()
        for meaning in meanings:
            if flag_regx.match(meaning) is None:
                valid_meanings.assert_true(False,
                                           "flag_meanings attribute defined an illegal flag meaning "
                                           "{}".format(meaning))
        return valid_meanings.to_result()

    ###############################################################################
    #
    # Chapter 4: Coordinate Types
    #
    ###############################################################################

    def check_coordinate_types(self, ds):
        '''
        Check the axis attribute of coordinate variables

        CF §4 The attribute axis may be attached to a coordinate variable and
        given one of the values X, Y, Z or T which stand for a longitude,
        latitude, vertical, or time axis respectively. Alternatively the
        standard_name attribute may be used for direct identification.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []
        coord_types = self._find_coord_vars(ds) + self._find_aux_coord_vars(ds)

        for name in coord_types:
            # Coordinate compressions should not be checked as a valid
            # coordinate, which they are not. They are a mechanism to project
            # an array of indices onto a 2-d grid containing valid coordinates.
            if cfutil.is_compression_coordinate(ds, name):
                continue
            variable = ds.variables[name]

            valid_coord = TestCtx(BaseCheck.MEDIUM, '§4 {} is a valid coordinate type'.format(name))

            axis = getattr(variable, 'axis', None)
            standard_name = getattr(variable, 'standard_name', None)

            valid_coord.assert_true(axis is not None or standard_name is not None,
                                    "coordinate types are should define either axis or standard_name attributes")
            ret_val.append(valid_coord.to_result())

            if axis is not None:
                valid_axis = self._check_axis(ds, name)
                ret_val.append(valid_axis)

            if standard_name is not None:
                valid_standard_name = self._check_coord_standard_name(ds, name)
                ret_val.append(valid_standard_name)

            if axis is not None and standard_name is not None:
                valid_mapping = self._check_coord_mapping(ds, name)
                ret_val.append(valid_mapping)

            valid_coordinate_type = self._check_coordinate_type(ds, name)
            ret_val.append(valid_coordinate_type)

        return ret_val

    def _check_axis(self, ds, name):
        '''
        Checks that the axis attribute is a string and an allowed value

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Name of the variable
        '''
        allowed_axis = ['T', 'X', 'Y', 'Z']
        variable = ds.variables[name]
        axis = variable.axis

        valid_axis = TestCtx(BaseCheck.HIGH, '§4 {} contains a valid axis'.format(name))
        axis_is_string = isinstance(axis, basestring),
        valid_axis.assert_true(axis_is_string and len(axis) > 0,
                               "axis attribute must be a non-empty string")

        # If axis isn't a string we can't continue any checks
        if not axis_is_string or len(axis) == 0:
            return valid_axis.to_result()

        valid_axis.assert_true(axis in allowed_axis,
                               "axis attribute must be T, X, Y, or Z, "
                               "currently {}".format(axis))

        return valid_axis.to_result()

    def _check_coord_standard_name(self, ds, name):
        '''
        Checks that the standard_name attribute for a coordinate type is a suggested value.

        :param netCDF4.Dataset ds: An open netCDF Dataset
        :param str name: Name of the variable
        '''
        allowed_standard_names = [
            'time',
            'longitude',
            'latitude',
            'height',
            'depth',
            'altitude'
        ]

        variable = ds.variables[name]
        standard_name = variable.standard_name

        # §4.5 Discrete Axis states that it is only recommended that the
        # coordinate types map to coordinate positions time, lat, lon etc.
        # Discrete axes are also ok.
        valid_standard_name = TestCtx(BaseCheck.LOW, '§4 {} has suggested standard_name for coordinate type'.format(name))
        valid_standard_name.assert_true(isinstance(standard_name, basestring),
                                        "standard_name is not a string")

        if not isinstance(standard_name, basestring):
            return valid_standard_name.to_result()

        valid_standard_name.assert_true(standard_name in allowed_standard_names,
                                        "standard_name attribute for coordinate types is suggested to be "
                                        "time, longitude, latitude, height, depth or altitude")

        return valid_standard_name.to_result()

    def _check_coord_mapping(self, ds, name):
        '''
        Checks that the axis maps to a suggested coordinate

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Name of the variable
        '''
        variable = ds.variables[name]
        axis = variable.axis
        standard_name = variable.standard_name

        allowed_map = {
            'T': ['time'],
            'X': ['longitude'],
            'Y': ['latitude'],
            'Z': ['height', 'depth', 'altitude']
        }

        valid_coord_mapping = TestCtx(BaseCheck.LOW, '§4 {} has suggested mapping from axis to standard_name'.format(name))
        axis_is_string = isinstance(axis, basestring)
        standard_name_is_string = isinstance(standard_name, basestring)

        valid_coord_mapping.assert_true(axis_is_string and standard_name_is_string,
                                        "axis and standard_name must be strings")

        if not standard_name_is_string or not axis_is_string:
            return valid_coord_mapping.to_result()

        valid_coord_mapping.assert_true(axis in allowed_map,
                                        "axis must be T, X, Y, or Z")
        if axis not in allowed_map:
            return valid_coord_mapping.to_result()

        valid_coord_mapping.assert_true(standard_name in allowed_map[axis],
                                        "standard_name for axis {} is suggested to be "
                                        "{}. Is currently {}"
                                        "".format(axis,
                                                  ', '.join(allowed_map[axis]),
                                                  standard_name))

        return valid_coord_mapping.to_result()

    def _check_coordinate_type(self, ds, name):
        '''
        Checks that the coordinate type is a coordinate variable

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Name of the variable
        '''
        variable = ds.variables[name]
        valid_coordinate_type = TestCtx(BaseCheck.LOW, '§4 recommends that coordinate types be coordinate variables')
        is_coordinate_variable = variable.dimensions == (name,)
        is_dimensionless = variable.ndim == 0
        # dimensionless coordinate types are loosely considered coordinate
        # variables of dimension size 1
        valid_coordinate_type.assert_true(is_coordinate_variable or is_dimensionless,
                                          "{} is not a coordinate variable".format(name))
        return valid_coordinate_type.to_result()

    def check_coordinate_vars_for_all_coordinate_types(self, ds):
        '''
        Check that coordinate variables exist for X, Y, Z, and T axes of the
        physical world.

        CF §4 We strongly recommend that coordinate variables be used for all
        coordinate types whenever they are applicable.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []
        # 1. Verify that for any known or common coordinate name as a dmension
        #    there is a coordinate variable for that dimension.
        known_coordinate_names = ('longitude', 'lon'   , 'x',
                                  'latitude' , 'lat'   , 'y',
                                  'vertical' , 'height', 'z',
                                  'time'               , 't')
        for k, v in ds.dimensions.items():
            if k.lower() in known_coordinate_names:
                valid = k in ds.variables
                result = Result(BaseCheck.MEDIUM, valid, '§4 Coordinate Variables')
                if not valid:
                    result.msgs = ['No coordinate variable for coordinate type %s' % k]

                ret_val.append(result)

        # @TODO: Additional verifiable requirements

        return ret_val

    def check_latitude(self, ds):
        '''
        Check variable(s) that define latitude and are defined correctly according to CF.

        CF §4.1 Variables representing latitude must always explicitly include
        the units attribute; there is no default value.  The recommended unit
        of latitude is degrees_north. Also acceptable are degree_north,
        degree_N, degrees_N, degreeN, and degreesN.

        Optionally, the latitude type may be indicated additionally by
        providing the standard_name attribute with the value latitude, and/or
        the axis attribute with the value Y.

        - Four checks per latitude variable
        - (H) latitude has units attribute
        - (M) latitude has an allowed units attribute
        - (L) latitude uses degrees_north (if not in rotated pole)
        - (M) latitude defines either standard_name or axis

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []

        allowed_lat_units = [
            'degrees_north',
            'degree_north',
            'degree_n',
            'degrees_n',
            'degreen',
            'degreesn'
        ]

        # Determine the grid mappings in this dataset
        grid_mapping = []
        grid_mapping_variables = cfutil.get_grid_mapping_variables(ds)
        for name in grid_mapping_variables:
            variable = ds.variables[name]
            grid_mapping_name = getattr(variable, 'grid_mapping_name', None)
            if grid_mapping_name:
                grid_mapping.append(grid_mapping_name)

        latitude_variables = cfutil.get_latitude_variables(ds)
        for latitude in latitude_variables:
            variable = ds.variables[latitude]
            units = getattr(variable, 'units', None)
            units_is_string = isinstance(units, basestring)
            standard_name = getattr(variable, 'standard_name', None)
            axis = getattr(variable, 'axis', None)

            # Check that latitude defines units
            valid_latitude = TestCtx(BaseCheck.HIGH, '§4.1 Latitude variable {} has required units attribute'.format(latitude))
            valid_latitude.assert_true(units is not None,
                                       "latitude variable '{}' must define units".format(latitude))
            ret_val.append(valid_latitude.to_result())

            # Check that latitude uses allowed units
            allowed_units = TestCtx(BaseCheck.MEDIUM, '§4.1 Latitude variable {} uses recommended units'.format(latitude))
            if 'rotated_latitude_longitude' in grid_mapping and standard_name == 'grid_latitude':
                allowed_units.assert_true(units == 'degrees',
                                          "latitude variable '{}' should use degrees for units in rotated pole grid"
                                          "".format(latitude))
            else:
                allowed_units.assert_true(units_is_string and units.lower() in allowed_lat_units,
                                          "latitude variable '{}' should define valid units for latitude"
                                          "".format(latitude))
            ret_val.append(allowed_units.to_result())

            # Check that latitude uses degrees_north
            recommended_units = TestCtx(BaseCheck.LOW, '§4.1 Latitude variable {} defines units using degrees_north'.format(latitude))
            if standard_name == 'latitude':
                recommended_units.assert_true(units == 'degrees_north',
                                              "CF recommends latitude variable '{}' to use units degrees_north"
                                              "".format(latitude))
                ret_val.append(recommended_units.to_result())

            # Check that latitude defines either standard_name or axis
            definition = TestCtx(BaseCheck.MEDIUM, '§4.1 Latitude variable {} defines either standard_name or axis'.format(latitude))
            definition.assert_true(standard_name == 'latitude' or axis == 'Y',
                                   "latitude variable '{}' should define standard_name='latitude' or axis='Y'"
                                   "".format(latitude))
            ret_val.append(definition.to_result())

        return ret_val

    def check_longitude(self, ds):
        '''
        Check variable(s) that define longitude and are defined correctly according to CF.

        CF §4.2 Variables representing longitude must always explicitly include
        the units attribute; there is no default value.  The recommended unit
        of longitude is degrees_east. Also acceptable are degree_east,
        degree_E, degrees_E, degreeE, and degreesE.

        Optionally, the longitude type may be indicated additionally by
        providing the standard_name attribute with the value longitude, and/or
        the axis attribute with the value X.

        - Four checks per longitude variable
        - (H) longitude has units attribute
        - (M) longitude has an allowed units attribute
        - (L) longitude uses degrees_east (if not in rotated pole)
        - (M) longitude defines either standard_name or axis

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []
        allowed_lon_units = [
            'degrees_east',
            'degree_east',
            'degree_e',
            'degrees_e',
            'degreee',
            'degreese'
        ]

        # Determine the grid mappings in this dataset
        grid_mapping = []
        grid_mapping_variables = cfutil.get_grid_mapping_variables(ds)
        for name in grid_mapping_variables:
            variable = ds.variables[name]
            grid_mapping_name = getattr(variable, 'grid_mapping_name', None)
            if grid_mapping_name:
                grid_mapping.append(grid_mapping_name)

        longitude_variables = cfutil.get_longitude_variables(ds)
        for longitude in longitude_variables:
            variable = ds.variables[longitude]
            units = getattr(variable, 'units', None)
            units_is_string = isinstance(units, basestring)
            standard_name = getattr(variable, 'standard_name', None)
            axis = getattr(variable, 'axis', None)

            # Check that longitude defines units
            valid_longitude = TestCtx(BaseCheck.HIGH, '§4.1 Longitude variable {} has required units attribute'.format(longitude))
            valid_longitude.assert_true(units is not None,
                                        "longitude variable '{}' must define units".format(longitude))
            ret_val.append(valid_longitude.to_result())

            # Check that longitude uses allowed units
            allowed_units = TestCtx(BaseCheck.MEDIUM, '§4.1 Longitude variable {} uses recommended units'.format(longitude))
            if 'rotated_latitude_longitude' in grid_mapping and standard_name == 'grid_longitude':
                allowed_units.assert_true(units == 'degrees',
                                          "longitude variable '{}' should use degrees for units in rotated pole grid"
                                          "".format(longitude))
            else:
                allowed_units.assert_true(units_is_string and units.lower() in allowed_lon_units,
                                          "longitude variable '{}' should define valid units for longitude"
                                          "".format(longitude))
            ret_val.append(allowed_units.to_result())

            # Check that longitude uses degrees_east
            recommended_units = TestCtx(BaseCheck.LOW, '§4.1 Longitude variable {} defines units using degrees_east'.format(longitude))
            if standard_name == 'longitude':
                recommended_units.assert_true(units == 'degrees_east',
                                              "CF recommends longitude variable '{}' to use units degrees_east"
                                              "".format(longitude))
                ret_val.append(recommended_units.to_result())

            # Check that longitude defines either standard_name or axis
            definition = TestCtx(BaseCheck.MEDIUM, '§4.1 Longitude variable {} defines either standard_name or axis'.format(longitude))
            definition.assert_true(standard_name == 'longitude' or axis == 'Y',
                                   "longitude variable '{}' should define standard_name='longitude' or axis='X'"
                                   "".format(longitude))
            ret_val.append(definition.to_result())

        return ret_val

    def check_vertical_coordinate(self, ds):
        '''
        Check variables defining height or depth are defined correctly
        according to CF.

        CF §4.3 Variables representing dimensional height or depth axes must
        always explicitly include the units attribute; there is no default
        value.

        The attribute positive is required if the vertical axis units are not a
        valid unit of pressure. The positive attribute may have the value up or
        down (case insensitive). This attribute may be applied to either
        coordinate variables or auxillary coordinate variables that contain
        vertical coordinate data.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''

        ret_val = []
        for k, v in ds.variables.items():
            if is_vertical_coordinate(k, v):
                # Vertical variables MUST have units
                has_units = hasattr(v, 'units')
                result = Result(BaseCheck.HIGH,
                                has_units,
                                '§4.3 Vertical coordinates contain valid attributes')
                ret_val.append(result)

                # If it's not pressure then it must have positive defined
                if not has_units:
                    result = Result(BaseCheck.HIGH,
                                    False,
                                    '§4.3 Vertical coordinates contain valid attributes', ['%s does not have units' % k])
                    ret_val.append(result)
                    continue

                # Do we have pressure?
                is_pressure = units_convertible('dbar', v.units)
                if is_pressure:
                    result = Result(BaseCheck.HIGH,
                                    True,
                                    '§4.3 Vertical coordinates contain valid attributes')
                # What about positive?
                elif getattr(v, 'positive', '').lower() in ('up', 'down'):
                    result = Result(BaseCheck.HIGH,
                                    True,
                                    '§4.3 Vertical coordinates contain valid attributes')
                # Not-compliant
                else:
                    result = Result(BaseCheck.HIGH,
                                    False,
                                    '§4.3 Vertical coordinates contain valid attributes',
                                    ['vertical variable %s needs to define positive attribute' % k])
                ret_val.append(result)
        return ret_val

    def check_dimensional_vertical_coordinate(self, ds):
        '''
        Check units for variables defining vertical position are valid under
        CF.

        CF §4.3.1 The units attribute for dimensional coordinates will be a string
        formatted as per the udunits.dat file.

        The acceptable units for vertical (depth or height) coordinate variables
        are:
        - units of pressure as listed in the file udunits.dat. For vertical axes
          the most commonly used of these include include bar, millibar,
          decibar, atmosphere (atm), pascal (Pa), and hPa.
        - units of length as listed in the file udunits.dat. For vertical axes
          the most commonly used of these include meter (metre, m), and
          kilometer (km).
        - other units listed in the file udunits.dat that may under certain
          circumstances reference vertical position such as units of density or
          temperature.

        Plural forms are also acceptable.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []
        for k, v in ds.variables.items():
            # If this is not a vertical coordinate
            if not is_vertical_coordinate(k, v):
                continue

            # If this is not height or depth
            vertical_coordinates = ('height', 'depth')
            if k not in vertical_coordinates and \
                    getattr(v, 'standard_name', '') not in vertical_coordinates:
                continue

            # Satisfies 4.3.1
            # Pressure or length is okay
            is_pressure = units_convertible(getattr(v, 'units', '1'), 'dbar')
            is_length   = units_convertible(getattr(v, 'units', '1'), 'm')
            is_temp     = units_convertible(getattr(v, 'units', '1'), 'degrees_C')
            is_density  = units_convertible(getattr(v, 'units', '1'), 'kg m-3')

            if is_pressure or is_length:
                result = Result(BaseCheck.HIGH, True,
                                '§4.3.1 Vertical dimension coordinates contain valid attributes',
                                [])

            # Temperature or Density are okay as well
            elif is_temp or is_density:
                result = Result(BaseCheck.HIGH, True,
                                '§4.3.1 Vertical dimension coordinates contain valid attributes',
                                [])
            else:
                result = Result(BaseCheck.HIGH, False,
                                '§4.3.1 Vertical dimension coordinates contain valid attributes',
                                ['incorrect vertical units'])
            ret_val.append(result)

        return ret_val

    def check_dimensionless_vertical_coordinate(self, ds):
        '''
        Check the validity of dimensionless coordinates under CF

        CF §4.3.2 The units attribute is not required for dimensionless
        coordinates.

        The standard_name attribute associates a coordinate with its definition
        from Appendix D, Dimensionless Vertical Coordinates. The definition
        provides a mapping between the dimensionless coordinate values and
        dimensional values that can positively and uniquely indicate the
        location of the data.

        A new attribute, formula_terms, is used to associate terms in the
        definitions with variables in a netCDF file.  To maintain backwards
        compatibility with COARDS the use of these attributes is not required,
        but is strongly recommended.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []

        dimless = dict(dimless_vertical_coordinates)
        for k, v in ds.variables.items():
            std_name = getattr(v, 'standard_name', '')
            if std_name not in dimless:
                continue
            # Determine if the regex matches for formula_terms
            valid_formula = re.match(dimless[std_name],
                                     getattr(v, 'formula_terms', ''))

            if valid_formula is not None:
                result = Result(BaseCheck.MEDIUM,
                                True,
                                '§4.3.2 Dimensionless Coordinates and formula_terms')
            else:
                result = Result(BaseCheck.MEDIUM,
                                False,
                                '§4.3.2 Dimensionless Coordinates and formula_terms',
                                ['formula_terms missing from dimensionless coordinate %s' % k])
            ret_val.append(result)

            # Determine that each of the terms actually exists
            # If formula_terms wasn't defined then this fails
            if not valid_formula:
                result = Result(BaseCheck.MEDIUM,
                                False,
                                '§4.3.2 Dimensionless Coordinates and formula_terms',
                                ['formula_terms not defined for dimensionless coordinate %s' % k])
                ret_val.append(result)
                continue

            # Check the terms
            missing_terms = []
            groups = valid_formula.groups()
            for i in range(1, len(groups), 2):
                varname = groups[i]
                if varname not in ds.variables:
                    missing_terms.append(varname)
            # Report the missing terms
            result = Result(BaseCheck.MEDIUM,
                            not missing_terms,
                            '§4.3.2 Dimensionless Coordinates and formula_terms',
                            ['%s missing for dimensionless coordinate %s' % (i, k) for i in missing_terms])

            ret_val.append(result)

        return ret_val

    def check_time_coordinate(self, ds):
        '''
        Check variables defining time are valid under CF

        CF §4.4 Variables representing time must always explicitly include the
        units attribute; there is no default value.

        The units attribute takes a string value formatted as per the
        recommendations in the Udunits package.

        The acceptable units for time are listed in the udunits.dat file. The
        most commonly used of these strings (and their abbreviations) includes
        day (d), hour (hr, h), minute (min) and second (sec, s). Plural forms
        are also acceptable. The reference time string (appearing after the
        identifier since) may include date alone; date and time; or date, time,
        and time zone. The reference time is required. A reference time in year
        0 has a special meaning (see Section 7.4, "Climatological Statistics").

        Recommend that the unit year be used with caution. It is not a calendar
        year.  For similar reasons the unit month should also be used with
        caution.

        A time coordinate is identifiable from its units string alone.
        Optionally, the time coordinate may be indicated additionally by
        providing the standard_name attribute with an appropriate value, and/or
        the axis attribute with the value T.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''

        ret_val = []
        for k, v in ds.variables.items():
            if not is_time_variable(k, v):
                continue
            # Has units
            has_units = hasattr(v, 'units')
            if not has_units:
                result = Result(BaseCheck.HIGH,
                                False,
                                '§4.4 Time coordinate variable and attributes', ['%s does not have units' % k])
                ret_val.append(result)
                continue
            # Correct and identifiable units
            result = Result(BaseCheck.HIGH,
                            True,
                            '§4.4 Time coordinate variable and attributes')
            ret_val.append(result)
            correct_units = units_temporal(v.units)
            reasoning = None
            if not correct_units:
                reasoning = ['%s doesn not have correct time units' % k]
            result = Result(BaseCheck.HIGH,
                            correct_units,
                            '§4.4 Time coordinate variable and attributes',
                            reasoning)
            ret_val.append(result)

        return ret_val

    def check_calendar(self, ds):
        '''
        Check the calendar attribute for variables defining time and ensure it
        is a valid calendar prescribed by CF.

        CF §4.4.1 In order to calculate a new date and time given a base date, base
        time and a time increment one must know what calendar to use.

        The values currently defined for calendar are:
        - gregorian or standard
        - proleptic_gregorian
        - noleap or 365_day
        - all_leap or 366_day
        - 360_day
        - julian
        - none

        The calendar attribute may be set to none in climate experiments that
        simulate a fixed time of year.
        The time of year is indicated by the date in the reference time of the
        units attribute.

        If none of the calendars defined above applies, a non-standard calendar
        can be defined. The lengths of each month are explicitly defined with
        the month_lengths attribute of the time axis.

        If leap years are included, then two other attributes of the time axis
        should also be defined:

        leap_year, leap_month

        The calendar attribute is not required when a non-standard calendar is
        being used. It is sufficient to define the calendar using the
        month_lengths attribute, along with leap_year, and leap_month as
        appropriate. However, the calendar attribute is allowed to take
        non-standard values and in that case defining the non-standard calendar
        using the appropriate attributes is required.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        valid_calendars = [
            'gregorian',
            'standard',
            'proleptic_gregorian',
            'noleap',
            '365_day',
            'all_leap',
            '366_day',
            '360_day',
            'julian',
            'none'
        ]

        ret_val = []

        for k, v in ds.variables.items():
            if not is_time_variable(k, v):
                continue
            reasoning = None
            has_calendar = hasattr(v, 'calendar')
            if not has_calendar:
                reasoning = ['Variable %s should have a calendar attribute' % k]
            result = Result(BaseCheck.LOW,
                            has_calendar,
                            '§4.4.1 Time and calendar',
                            reasoning)
            ret_val.append(result)

            if has_calendar:
                reasoning = None
                valid_calendar = v.calendar in valid_calendars
                if not valid_calendar:
                    reasoning = ["Variable %s should have a valid calendar: '%s' is not a valid calendar" % (k, v.calendar)]
                    result = Result(BaseCheck.LOW,
                                    valid_calendar,
                                    '§4.4.1 Time and calendar',
                                    reasoning)
                ret_val.append(result)

        return ret_val

    ###############################################################################
    #
    # Chapter 5: Coordinate Systems
    #
    ###############################################################################

    def _is_station_var(self, var):
        if getattr(var, 'standard_name', None) in ('platform_name', 'station_name', 'instrument_name'):
            return True
        return False

    def _get_coord_vars(self, ds):
        coord_vars = []
        for name, var in ds.variables.items():
            if (name,) == var.dimensions:
                coord_vars.append(name)
        return coord_vars

    def check_two_dimensional(self, ds):
        """
        5.2 The latitude and longitude coordinates of a horizontal grid that was
        not defined as a Cartesian product of latitude and longitude axes, can
        sometimes be represented using two-dimensional coordinate variables.
        These variables are identified as coordinates by use of the coordinates
        attribute.

        For each variable, if the variable has a coordinates attribute:
          for each coordinate defined, verify that the coordinate:
            is either a coordinate variable OR comprises coordinate variables

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        reported_reference_variables = []

        for name, var in ds.variables.items():
            self_reference_variables = set()
            g = NCGraph(ds, name, var, self_reference_variables)

            reasoning = []

            for coord in g.coords:
                if coord in reported_reference_variables:
                    continue
                if coord in self_reference_variables:
                    reasoning.append("Variable %s's coordinate references itself" % (coord))

                    result = Result(BaseCheck.HIGH,
                                    False,
                                    ('§5.2 Latitude and longitude coordinates of a horizontal grid', coord, 'coordinates_reference_itself'),
                                    reasoning)
                    ret_val.append(result)
                else:
                    result = Result(BaseCheck.HIGH,
                                    True,
                                    ('§5.2 Latitude and longitude coordinates of a horizontal grid', coord, 'coordinates_reference_itself'),
                                    [])
                    ret_val.append(result)
                reported_reference_variables.append(coord)

            # Determine if 2-D coordinate variables (Lat and Lon are of shape (i,j)
            valid = True
            for _, graph in g.coords.items():
                try:
                    assert graph.ndim == 2
                except AssertionError:
                    valid = False
                except AttributeError:
                    # When graph is None
                    pass

            if len(g.coords) == 2 and valid is True:
                # ------------------------------------------------------------
                # Check all the dims are coordinate variables
                # ------------------------------------------------------------
                valid_dims = True
                reasoning = []
                for dim in g.dims.keys():
                    if dim not in ds.variables:
                        valid_dims = False
                        reasoning.append("Variable %s's dimension %s is not a coordinate variable" % (name, dim))

                result = Result(BaseCheck.HIGH,
                                valid_dims,
                                ('§5.2 Latitude and longitude coordinates of a horizontal grid', name, '2d_hgrid_valid_dimensions'),
                                reasoning)
                ret_val.append(result)

                # ------------------------------------------------------------
                # Check that the coordinates are correct
                # ------------------------------------------------------------
                valid_2d = True
                reasoning = []
                for cname, coord in g.coords.items():
                    if coord is None:
                        valid_2d = False
                        reasoning.append("Variable %s's coordinate, %s, is not a coordinate or auxiliary variable" % (name, cname))
                        continue
                    for dim in coord.dims.keys():
                        if dim not in g.dims:
                            valid_2d = False
                            reasoning.append("Variable %s's coordinate, %s, does not share dimension %s with the variable" % (name, cname, dim))
                result = Result(BaseCheck.MEDIUM,
                                valid_2d,
                                ('§5.2 Latitude and longitude coordinates of a horizontal grid', name, 'valid_coordinates'),
                                reasoning)
                ret_val.append(result)
                # ------------------------------------------------------------
                # Can make lat/lon?
                # ------------------------------------------------------------

                lat_check = False
                lon_check = False
                reasoning = []
                for cname, coord in g.coords.items():

                    if coord is not None and coord.units in ['degrees_north', 'degree_north', 'degrees_N', 'degree_N', 'degreesN', 'degreeN']:
                        lat_check = True
                    elif coord is not None and coord.units in ['degrees_east', 'degree_east', 'degrees_E', 'degree_E', 'degreesE', 'degreeE']:
                        lon_check = True
                    else:
                        reasoning.append("coordinate {} is not a correct lat/lon variable".format(cname))

                result = Result(BaseCheck.HIGH,
                                lat_check and lon_check,
                                ('§5.2 Latitude and longitude coordinates of a horizontal grid', name, 'lat_lon_correct'),
                                reasoning)
                ret_val.append(result)
            else:
                continue  # Not a 2d horizontal grid

        return ret_val

    def check_reduced_horizontal_grid(self, ds):
        """
        5.3 A "reduced" longitude-latitude grid is one in which the points are
        arranged along constant latitude lines with the number of points on a
        latitude line decreasing toward the poles.

        Recommend that this type of gridded data be stored using the compression
        scheme described in Section 8.2, "Compression by Gathering". The
        compressed latitude and longitude auxiliary coordinate variables are
        identified by the coordinates attribute.

        """
        ret_val = []
        coord_vars = self._get_coord_vars(ds)

        for name, var in ds.variables.items():
            if name in coord_vars:
                continue
            if not hasattr(var, 'coordinates'):
                continue

            reasoning = []
            valid_in_variables = True
            valid_dim = True
            valid_coord = True
            valid_cdim = True
            result = None

            coords = var.coordinates.split(' ')
            for coord in coords:
                is_reduced_horizontal_grid = True
                if coord not in ds.variables:
                    valid_in_variables = False
                    reasoning.append("Coordinate %s is not a proper variable" % coord)
                    continue

                for dim_name in ds.variables[coord].dimensions:

                    if dim_name not in var.dimensions:
                        valid_dim = False
                        reasoning.append("Coordinate %s's dimension, %s, is not a dimension of %s" % (coord, dim_name, name))
                        continue

                    if dim_name not in coord_vars:
                        valid_coord = False
                        reasoning.append("Coordinate %s's dimension, %s, is not a coordinate variable" % (coord, dim_name))
                        continue

                    dim = ds.variables[dim_name]
                    if not hasattr(dim, 'compress'):
                        is_reduced_horizontal_grid = False
                        continue

                    compress_dims = dim.compress.split(' ')
                    for cdim in compress_dims:
                        if cdim not in ds.dimensions:
                            valid_cdim = False
                            reasoning.append("Dimension %s compresses non-existent dimension, %s" % (dim_name, cdim))
                            continue
            if is_reduced_horizontal_grid is True:
                result = Result(BaseCheck.MEDIUM,
                                (valid_in_variables and valid_dim and valid_coord and valid_cdim),
                                ('§5.3 Is reduced horizontal grid', name, 'is_reduced_horizontal_grid'),
                                reasoning)

            if result:
                ret_val.append(result)

        return ret_val

    # grid mapping dictionary, appendix F
    grid_mapping_dict = {
        'albers_conical_equal_area': [('longitude_of_central_meridian', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'azimuthal_equidistant': [('longitude_of_projection_origin', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'lambert_cylindrical_equal_area': [('longitude_of_central_meridian', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate'), ('standard_parallel', 'scale_factor_at_projection_origin')],
        'lambert_azimuthal_equal_area': [('longitude_of_projection_origin', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'lambert_conformal_conic': [('standard_parallel', 'longitude_of_central_meridian', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'latitude_longitude': [(), (), ('longitude', 'latitude')],
        'mercator': [('longitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate'), ('standard_parallel', 'scale_factor_at_projection_origin')],
        'orthographic': [('longitude_of_projection_origin', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'polar_stereographic': [('straight_vertical_longitude_from_pole', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate'), ('standard_parallel', 'scale_factor_at_projection_origin')],
        'rotated_latitude_longitude': [('grid_north_pole_latitude', 'grid_north_pole_longitude'), ('north_pole_grid_longitude'), ('grid_latitude', 'grid_longitude')],
        'stereographic': [('longitude_of_projection_origin', 'latitude_of_projection_origin', 'scale_factor_at_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'transverse_mercator': [('scale_factor_at_central_meridian', 'longitude_of_central_meridian', 'latitude_of_projection_origin', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')],
        'vertical_perspective': [('longitude_of_projection_origin', 'latitude_of_projection_origin', 'perspective_point_height', 'false_easting', 'false_northing'), (), ('projection_x_coordinate', 'projection_y_coordinate')]
    }

    def check_horz_crs_grid_mappings_projections(self, ds):
        """
        5.6 When the coordinate variables for a horizontal grid are not
        longitude and latitude, it is required that the true latitude and
        longitude coordinates be supplied via the coordinates attribute. If in
        addition it is desired to describe the mapping between the given
        coordinate variables and the true latitude and longitude coordinates,
        the attribute grid_mapping may be used to supply this description.

        This attribute is attached to data variables so that variables with
        different mappings may be present in a single file. The attribute takes
        a string value which is the name of another variable in the file that
        provides the description of the mapping via a collection of attached
        attributes. This variable is called a grid mapping variable and is of
        arbitrary type since it contains no data. Its purpose is to act as a
        container for the attributes that define the mapping.

        The one attribute that all grid mapping variables must have is
        grid_mapping_name which takes a string value that contains the mapping's
        name. The other attributes that define a specific mapping depend on the
        value of grid_mapping_name. The valid values of grid_mapping_name along
        with the attributes that provide specific map parameter values are
        described in Appendix F, Grid Mappings.

        When the coordinate variables for a horizontal grid are longitude and
        latitude, a grid mapping variable with grid_mapping_name of
        latitude_longitude may be used to specify the ellipsoid and prime
        meridian.


        In order to make use of a grid mapping to directly calculate latitude
        and longitude values it is necessary to associate the coordinate
        variables with the independent variables of the mapping. This is done by
        assigning a standard_name to the coordinate variable. The appropriate
        values of the standard_name depend on the grid mapping and are given in
        Appendix F, Grid Mappings.
        """

        ret_val = []
        reasoning = []

        for name, var in ds.variables.items():
            valid_mapping_count = 0
            total_mapping_count = 0
            if hasattr(var, 'grid_mapping_name'):
                total_mapping_count = 1

                mapping = getattr(var, 'grid_mapping_name', '')
                if mapping in iter(self.grid_mapping_dict.keys()):
                    valid_mapping_count = valid_mapping_count + 1
                else:
                    reasoning.append('The grid_mapping_name attribute is not an accepted value.  See Appendix F.')

                for each in self.grid_mapping_dict[mapping][0]:
                    total_mapping_count = total_mapping_count + 1
                    if each in dir(var):
                        valid_mapping_count = valid_mapping_count + 1
                    else:
                        reasoning.append('The map parameters are not accepted values.  See Appendix F.')

                if len(self.grid_mapping_dict[mapping]) >= 4:
                    for each in self.grid_mapping_dict[mapping][3:]:
                        every_flag = 0
                        total_mapping_count = total_mapping_count + 1
                        for every in each:
                            if every in dir(var):
                                valid_mapping_count = valid_mapping_count + 1
                                every_flag = every_flag + 1

                        if every_flag == 0:
                            reasoning.append('Neither of the "either/or" parameters are present')
                        if every_flag == 2:
                            valid_mapping_count = valid_mapping_count - 2

                total_mapping_count = total_mapping_count + len(self.grid_mapping_dict[mapping][2])
                for name_again, var_again in ds.variables.items():
                    if hasattr(var_again, 'standard_name'):
                        if var_again.standard_name in self.grid_mapping_dict[mapping][2]:
                            valid_mapping_count = valid_mapping_count + 1

                result = Result(BaseCheck.MEDIUM,
                                (valid_mapping_count, total_mapping_count),
                                ('§5.6 Grid mapping projection present', name, 'horz_crs_grid_mappings_projections'),
                                reasoning)

                ret_val.append(result)

        return ret_val

    def check_scalar_coordinate_system(self, ds):
        """
        5.7 When a variable has an associated coordinate which is single-valued, that coordinate may be represented as a
        scalar variable. Since there is no associated dimension these scalar coordinate variables should be attached to a
        data variable via the coordinates attribute.

        Once a name is used for a scalar coordinate variable it can not be used for a 1D coordinate variable. For this
        reason we strongly recommend against using a name for a scalar coordinate variable that matches the name of any
        dimension in the file.
        """
        ret_val = []

        for name, var in ds.variables.items():
            valid_scalar_coordinate_var = 0
            total_scalar_coordinate_var = 0
            reasoning = []

            if not hasattr(var, 'coordinates'):
                continue

            for coordinate in getattr(var, 'coordinates', '').split(" "):
                if coordinate in ds.variables:
                    if ds.variables[coordinate].shape == ():
                        total_scalar_coordinate_var += 1
                        if coordinate not in list(ds.dimensions.keys()):
                            valid_scalar_coordinate_var += 1
                        else:
                            reasoning.append('Scalar coordinate var (%s) of var (%s) is correct size but is present in the dimensions list, which is not allowed.' % (coordinate, name))

            if total_scalar_coordinate_var > 0:
                result = Result(BaseCheck.MEDIUM,
                                (valid_scalar_coordinate_var, total_scalar_coordinate_var),
                                ('§5.7 Scalar coordinate variables', name, 'scalar_coordinates'),
                                reasoning)
                ret_val.append(result)

        return ret_val

    ###############################################################################
    #
    # Chapter 6: Labels and Alternative Coordinates
    #
    ###############################################################################

    def check_geographic_region(self, ds):
        """
        6.1.1 When data is representative of geographic regions which can be identified by names but which have complex
        boundaries that cannot practically be specified using longitude and latitude boundary coordinates, a labeled
        axis should be used to identify the regions.

        Recommend that the names be chosen from the list of standardized region names whenever possible. To indicate
        that the label values are standardized the variable that contains the labels must be given the standard_name
        attribute with the value region.
        """
        ret_val = []
        reasoning = []
        region_list = [
            'africa',
            'antarctica',
            'arabian_sea',
            'aral_sea',
            'arctic_ocean',
            'asia',
            'atlantic_ocean',
            'australia',
            'baltic_sea',
            'barents_opening',
            'barents_sea',
            'beaufort_sea',
            'bellingshausen_sea',
            'bering_sea',
            'bering_strait',
            'black_sea',
            'canadian_archipelago',
            'caribbean_sea',
            'caspian_sea',
            'central_america',
            'chukchi_sea',
            'contiguous_united_states',
            'denmark_strait',
            'drake_passage',
            'east_china_sea',
            'english_channel',
            'eurasia',
            'europe',
            'faroe_scotland_channel',
            'florida_bahamas_strait',
            'fram_strait',
            'global',
            'global_land',
            'global_ocean',
            'great_lakes',
            'greenland',
            'gulf_of_alaska',
            'gulf_of_mexico',
            'hudson_bay',
            'iceland_faroe_channel',
            'indian_ocean',
            'indonesian_throughflow',
            'indo_pacific_ocean',
            'irish_sea',
            'lake_baykal',
            'lake_chad',
            'lake_malawi',
            'lake_tanganyika',
            'lake_victoria',
            'mediterranean_sea',
            'mozambique_channel',
            'north_america',
            'north_sea',
            'norwegian_sea',
            'pacific_equatorial_undercurrent',
            'pacific_ocean',
            'persian_gulf',
            'red_sea',
            'ross_sea',
            'sea_of_japan',
            'sea_of_okhotsk',
            'south_america',
            'south_china_sea',
            'southern_ocean',
            'taiwan_luzon_straits',
            'weddell_sea',
            'windward_passage',
            'yellow_sea'
        ]

        for name, var in ds.variables.items():
            if getattr(var, 'standard_name', '') == 'region':
                if ''.join(var[:].astype(str)).lower() in region_list:
                    result = Result(BaseCheck.LOW,
                                    True,
                                    ('§6.1.1 Geographic region specified', name, 'geographic_region'),
                                    reasoning)
                else:
                    reasoning.append('The Region Value is not from the allowable list.')
                    result = Result(BaseCheck.LOW,
                                    False,
                                    ('§6.1.1 Geographic region specified', name, 'geographic_region'),
                                    reasoning)
                ret_val.append(result)
        return ret_val

    ###############################################################################
    #
    # Chapter 7: Data Representative of Cells
    #
    ###############################################################################

    def check_cell_boundaries(self, ds):
        """
        Checks the dimensions of cell boundary variables to ensure they are CF compliant.

        7.1 To represent cells we add the attribute bounds to the appropriate coordinate variable(s). The value of bounds
        is the name of the variable that contains the vertices of the cell boundaries. We refer to this type of variable as
        a "boundary variable." A boundary variable will have one more dimension than its associated coordinate or auxiliary
        coordinate variable. The additional dimension should be the most rapidly varying one, and its size is the maximum
        number of cell vertices.

        Applications that process cell boundary data often times need to determine whether or not adjacent cells share an
        edge. In order to facilitate this type of processing the following restrictions are placed on the data in boundary
        variables:

        Bounds for 1-D coordinate variables

            For a coordinate variable such as lat(lat) with associated boundary variable latbnd(x,2), the interval endpoints
            must be ordered consistently with the associated coordinate, e.g., for an increasing coordinate, lat(1) > lat(0)
            implies latbnd(i,1) >= latbnd(i,0) for all i

            If adjacent intervals are contiguous, the shared endpoint must be represented indentically in each instance where
            it occurs in the boundary variable. For example, if the intervals that contain grid points lat(i) and lat(i+1) are
            contiguous, then latbnd(i+1,0) = latbnd(i,1).

        Bounds for 2-D coordinate variables with 4-sided cells

            In the case where the horizontal grid is described by two-dimensional auxiliary coordinate variables in latitude
            lat(n,m) and longitude lon(n,m), and the associated cells are four-sided, then the boundary variables are given
            in the form latbnd(n,m,4) and lonbnd(n,m,4), where the trailing index runs over the four vertices of the cells.

        Bounds for multi-dimensional coordinate variables with p-sided cells

            In all other cases, the bounds should be dimensioned (...,n,p), where (...,n) are the dimensions of the auxiliary
            coordinate variables, and p the number of vertices of the cells. The vertices must be traversed anticlockwise in the
            lon-lat plane as viewed from above. The starting vertex is not specified.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        """
        ret_val = []
        reasoning = []

        # Here variable is usually a coordinate variable and is mapped to a boundary variable
        for variable_name, boundary_variable_name in cfutil.get_cell_boundary_map(ds):
            variable = ds.variables[variable_name]
            boundary_variable = ds.variables[boundary_variable_name]
            valid = True

            if boundary_variable.ndim != variable.ndim + 1:
                valid = False
                reasoning.append('The number of dimensions of the Coordinate Variable is %s, but the '
                                 'number of dimensions of the Boundary Variable is %s.' %
                                 (variable_name.ndim, boundary_variable_name.ndim))

            result = Result(BaseCheck.MEDIUM,
                            valid,
                            ('§7.1 Cell boundaries', variable_name, 'cell_boundaries'),
                            reasoning)
            ret_val.append(result)
            reasoning = []

        return ret_val

    def check_cell_measures(self, ds):
        """
        7.2 To indicate extra information about the spatial properties of a variable's grid cells, a cell_measures attribute may
        be defined for a variable. This is a string attribute comprising a list of blank-separated pairs of words of the form
        "measure: name". "area" and "volume" are the only defined measures.

        The "name" is the name of the variable containing the measure values, which we refer to as a "measure variable". The
        dimensions of the measure variable should be the same as or a subset of the dimensions of the variable to which they are
        related, but their order is not restricted.

        The variable must have a units attribute and may have other attributes such as a standard_name.
        """
        ret_val = []
        reasoning = []
        for name, var in ds.variables.items():
            for dim in var.dimensions:
                if getattr(var, 'cell_measures', ''):
                    measures = getattr(var, 'coordinates', '')
                    measures = measures.split(': ')
                    if measures[0] not in ['area', 'volume']:
                        reasoning.append("The 'measures' field is not equal to 'area' or 'volume'.")
                        return Result(BaseCheck.MEDIUM,
                                      False,
                                      ('§7.2 Cell measures', name, 'cell_measures'),
                                      reasoning)
                    for every, attri in ds.variables.items():
                        if every == measures[1]:
                            for dimi in attri.dimensions:
                                if dimi in var.dimensions:
                                    valid = True
                                else:
                                    reasoning.append('The measure variable dimensions are not a set or subset of the cell_measure variable.')
                                    valid = False

                    result = Result(BaseCheck.MEDIUM,
                                    valid,
                                    ('§7.2 Cell measures', name, 'cell_measures'),
                                    reasoning)
                    ret_val.append(result)

        return ret_val

    def check_cell_methods(self, ds):
        """
        7.3 To describe the characteristic of a field that is represented by cell values, we define the cell_methods attribute
        of the variable. This is a string attribute comprising a list of blank-separated words of the form "name: method". Each
        "name: method" pair indicates that for an axis identified by name, the cell values representing the field have been
        determined or derived by the specified method.

        name can be a dimension of the variable, a scalar coordinate variable, a valid standard name, or the word "area"

        values of method should be selected from the list in Appendix E, Cell Methods, which includes point, sum, mean, maximum,
        minimum, mid_range, standard_deviation, variance, mode, and median. Case is not significant in the method name. Some
        methods (e.g., variance) imply a change of units of the variable, as is indicated in Appendix E, Cell Methods.

        Because the default interpretation for an intensive quantity differs from that of an extensive quantity and because this
        distinction may not be understood by some users of the data, it is recommended that every data variable include for each
        of its dimensions and each of its scalar coordinate variables the cell_methods information of interest (unless this
        information would not be meaningful). It is especially recommended that cell_methods be explicitly specified for each
        spatio-temporal dimension and each spatio-temporal scalar coordinate variable.
        """

        _areatype_names     = ["bare_ground",
                               "all_area_types",
                               "burnt_vegetation",
                               "c3_plant_functional_types",
                               "c4_plant_functional_types",
                               "clear_sky",
                               "cloud",
                               "crops",
                               "floating_ice",
                               "ice_free_land",
                               "ice_free_sea",
                               "lake_ice_or_sea_ice",
                               "land",
                               "land_ice",
                               "natural_grasses",
                               "pastures",
                               "primary_deciduous_trees",
                               "primary_evergreen_trees",
                               "sea",
                               "sea_ice",
                               "secondary_deciduous_trees",
                               "secondary_evergreen_trees",
                               "shrubs"
                               "snow",
                               "trees"
                               "vegetation"]
        methods = [ 'point',
                    'sum',
                    'mean',
                    'maximum',
                    'minimum',
                    'mid_range',
                    'standard_deviation',
                    'variance',
                    'mode',
                    'median']

        ret_val = []
        reasoning = []
        paragraph = ''
        pvars = re.compile('\(.*?\)|(\w*?):')
        psep = re.compile('((?P<var>\w+): (?P<method>\w+) ?(?P<where>where (?P<wtypevar>\w+) ?(?P<over>over (?P<otypevar>\w+))?| ?)(?P<brace>\(((?P<brace_wunit>\w+): (\d+) (?P<unit>\w+)|(?P<brace_opt>\w+): (\w+))\))*)')

        for name, var in ds.variables.items():
            if getattr(var, 'cell_methods', '') :
                method = getattr(var, 'cell_methods', '')

                total_name_count = 0
                cell_dims = []
                for match in re.finditer(pvars, method):
                    if (match.groups()[0] is not None):
                        cell_dims.append(match.groups()[0])
                        total_name_count = total_name_count + 1

                # print "cell_methods_check: number DIMs", total_name_count

                # check that the name is valid
                valid_name_count = 0
                for match in re.finditer(psep, method):
                    # print 'dict ', match.groupdict()
                    if match.group('var') in ds.variables[name].dimensions:
                        valid_name_count = valid_name_count + 1
                    elif match.group('var') == 'area':
                        valid_name_count = valid_name_count + 1
                    elif match.group('var') in getattr(var, "coordinates", ""):
                        valid_name_count = valid_name_count + 1
                    else:
                        reasoning.append('The name field does not match a dimension, area or coordinate.')

                result = Result(BaseCheck.MEDIUM,
                                (valid_name_count, total_name_count),
                                ('§7.3 Cell Methods', name, 'cell_methods_name'),
                                reasoning)
                ret_val.append(result)

                # Checks if the method value of the 'name: method' pair is acceptable
                reasoning = []
                methods = ['point', 'sum', 'mean', 'maximum', 'minimum', 'mid_range', 'standard_deviation', 'variance', 'mode', 'median']

                valid_method_count = 0
                for match in re.finditer(psep, method):
                    # print 'dict ', match.groupdict()
                    if match.group('method') in methods:
                        valid_method_count = valid_method_count + 1
                    else:
                        reasoning.append('The method field does not match a valid method value.')

                total_method_count = total_name_count  # all dims must have a valid method

                result = Result(BaseCheck.MEDIUM,
                                (valid_method_count, total_method_count),
                                ('§7.3 Cell Methods', name, 'cell_methods_method'),
                                reasoning)
                ret_val.append(result)

                # check the method modifier 'name: method (modifier)'
                reasoning = []
                valid_brace_count = 0
                total_brace_count = 0

                for match in re.finditer(psep, method):
                    if match.group('brace') is not None:
                        total_brace_count = total_brace_count + 1
                        if match.group('brace_wunit') == 'interval':
                            valid_brace_count = valid_brace_count + 1
                        elif match.group('brace_wunit') in ['comment', 'area']:
                            valid_brace_count = valid_brace_count + 1
                        else:
                            reasoning.append('The method modifier not valid.')

                result = Result(BaseCheck.MEDIUM,
                                (valid_brace_count, total_brace_count),
                                ('§7.3 Cell Methods', name, 'cell_methods_method_modifier'),
                                reasoning)
                ret_val.append(result)

                # Checks the 'method where' formats
                reasoning = []
                valid_area_count = 0
                total_area_count = 0

                for match in re.finditer(psep, method):
                    if len(match.group('where')) != 0:
                        if match.group('wtypevar') in _areatype_names:
                            total_area_count = total_area_count + 1
                            if match.group('otypevar') is not None:
                                if match.group('otypevar') in _areatype_names:
                                    valid_area_count = valid_area_count + 1
                                else:
                                    reasoning.append('The "name: method where type over _areatype_names" (' + match.group('otypevar') + ') format is not correct.')
                            else:
                                valid_area_count = valid_area_count + 1
                        else:
                            reasoning.append('The "name: method where _areatype_names" (' + match.group('wvartype') + ') format is not correct.')

                result = Result(BaseCheck.MEDIUM,
                                (valid_area_count, total_area_count),
                                ('§7.3 Cell Methods', name, 'cell_methods_area'),
                                reasoning)
                ret_val.append(result)

        # Checks the Climatology Variables - 7.4
        reasoning = []
        paragraph = []
        total_climate_count = 0
        valid_climate_count = 0
        for name, var in ds.variables.items():
            if getattr(var, 'climatology', ''):
                climate_dim = ds.variables[name].dimensions
                clim_method = getattr(var, 'climatology', '')

                for each in climate.split(" "):
                    paragraph.append(each)

                total_climate_count = total_climate_count + 1
                for name_again, var_again in ds.variables.items():
                    if getattr(var_again, "cell_methods", ""):
                        climate = getattr(var, 'cell_methods', '')
                        name_dim = ds.variables[name_again].dimensions
                        if len(climate_dim) > 0:
                            if climate_dim[0] in name_dim:
                                case1 = re.search(r"time: \w* within years time: \w* over years", climate)
                                case2 = re.search(r"time: \w* within days time: \w* over days$", climate)
                                case3 = re.search(r"time: \w* within days time: \w* over days time: \w* over years", climate)

                        if (case1 or case2 or case3) and len(ds.variables[clim_method].shape) == 2 and ds.variables[clim_method].shape[1] == 2 and ds.variables[clim_method].shape[0] == ds.variables[name_again].shape[0] :

                            valid_climate_count = 1
                        if not (case1 or case2 or case3):
                            reasoning.append('The "time: method within years/days over years/days" format is not correct.')

                        if not (len(ds.variables[clim_method].shape) == 2 and ds.variables[clim_method].shape[1] == 2 and ds.variables[clim_method].shape[0] == ds.variables[name_again].shape[0]):
                            reasoning.append('The dimensions of the climatology varaible is incorrect.')

                result = Result(BaseCheck.MEDIUM,
                                (valid_climate_count, total_climate_count),
                                ('§7.3 Cell Methods', name, 'cell_methods_climatology'),
                                reasoning)
                ret_val.append(result)

        return ret_val

    # def check_cell_methods_for_multi_axes(self, ds):
        """
        7.3.1 If a data value is representative of variation over a combination of axes, a single method should be prefixed by the
        names of all the dimensions involved (listed in any order, since in this case the order must be immaterial).

        There is no way to check this.  A warning should be posted explaining this method to the user!"

        """

    # def check_spacing_and_extra_info(self, ds):
        """
        7.3.2 To indicate more precisely how the cell method was applied, extra information may be included in parentheses ()
        after the identification of the method. This information includes standardized and non-standardized parts.

        The only standardized information is to provide the typical interval between the original data values to which the method
        was applied, in the situation where the present data values are statistically representative of original data values which
        had a finer spacing.

        The syntax is (interval: value unit), where value is a numerical value and unit is a string that can be recognized by
        UNIDATA's Udunits package.

        If the cell method applies to a combination of axes, they may have a common original interval. Alternatively, they may have
        separate intervals, which are matched to the names of axes by position.

        If there is both standardized and non-standardized information, the non-standardized follows the standardized information
        and the keyword comment:. If there is no standardized information, the keyword comment: should be omitted.

        A dimension of size one may be the result of "collapsing" an axis by some statistical operation, for instance by
        calculating a variance from time series data. We strongly recommend that dimensions of size one be retained (or scalar
        coordinate variables be defined) to enable documentation of the method (through the cell_methods attribute) and its
        domain (through the cell_bounds attribute).
        """

    # def check_stats_applying_to_portions_of_cells(self, ds):
        """
        7.3.3 By default, the statistical method indicated by cell_methods is assumed to have been evaluated over the entire
        horizontal area of the cell. Sometimes, however, it is useful to limit consideration to only a portion of a cell.

        One of two conventions may be used.

        The first convention is a method that can be used for the common case of a single area-type. In this case, the
        cell_methods attribute may include a string of the form "name: method where type".

        The second convention is the more general. In this case, the cell_methods entry is of the form "name: method where
        _areatype_names". Here _areatype_names is a string-valued auxiliary coordinate variable or string-valued scalar coordinate variable
        with a standard_name of area_type. The variable _areatype_names contains the name(s) of the selected portion(s) of the grid
        cell to which the method is applied.

        If the method is mean, various ways of calculating the mean can be distinguished in the cell_methods attribute with
        a string of the form "mean where type1 [over type2]". Here, type1 can be any of the possibilities allowed for _areatype_names
        or type (as specified in the two paragraphs preceding above Example). The same options apply to type2, except it is
        not allowed to be the name of an auxiliary coordinate variable with a dimension greater than one (ignoring the
        dimension accommodating the maximum string length)
        """

    # def check_cell_methods_with_no_coords(self, ds):
        """
        7.3.4 To provide an indication that a particular cell method is relevant to the data without having to provide a
        precise description of the corresponding cell, the "name" that appears in a "name: method" pair may be an
        appropriate standard_name (which identifies the dimension) or the string, "area" (rather than the name of a scalar
        coordinate variable or a dimension with a coordinate variable). This convention cannot be used, however, if the name
        of a dimension or scalar coordinate variable is identical to name.

        Recommend that whenever possible, cell bounds should be supplied by giving the variable a dimension of size one
        and attaching bounds to the associated coordinate variable.
        """

    # def check_climatological_statistics(self, ds):
        """
        7.4 A climatological time coordinate variable does not have a bounds attribute. Instead, it has a climatology
        attribute, which names a variable with dimensions (n,2), n being the dimension of the climatological time axis.
        Using the units and calendar of the time coordinate variable, element (i,0) of the climatology variable specifies
        the beginning of the first subinterval and element (i,1) the end of the last subinterval used to evaluate the
        climatological statistics with index i in the time dimension. The time coordinates should be values that are
        representative of the climatological time intervals, such that an application which does not recognise climatological
        time will nonetheless be able to make a reasonable interpretation.

        Valid values of the cell_methods attribute must be in one of the forms from the following list.

        - time: method1 within years   time: method2 over years
        - time: method1 within days   time: method2 over days
        - time: method1 within days   time: method2 over days   time: method3 over years

        The methods which can be specified are those listed in Appendix E, Cell Methods and each entry in the cell_methods
        attribute may also, contain non-standardised information in parentheses after the method.
        """

    ###############################################################################
    #
    # Chapter 8: Reduction of Dataset Size
    #
    ###############################################################################

    def check_packed_data(self, ds):
        """
        8.1 Simple packing may be achieved through the use of the optional NUG defined attributes scale_factor and
        add_offset. After the data values of a variable have been read, they are to be multiplied by the scale_factor,
        and have add_offset added to them.

        The units of a variable should be representative of the unpacked data.

        If the scale_factor and add_offset attributes are of the same data type as the associated variable, the unpacked
        data is assumed to be of the same data type as the packed data. However, if the scale_factor and add_offset
        attributes are of a different data type from the variable (containing the packed data) then the unpacked data
        should match the type of these attributes, which must both be of type float or both be of type double. An additional
        restriction in this case is that the variable containing the packed data must be of type byte, short or int. It is
        not advised to unpack an int into a float as there is a potential precision loss.

        When data to be packed contains missing values the attributes that indicate missing values (_FillValue, valid_min,
        valid_max, valid_range) must be of the same data type as the packed data.
        """
        ret_val = []
        for name, var in ds.variables.items():

            add_offset = getattr(var, 'add_offset', None)
            scale_factor = getattr(var, 'scale_factor', None)
            if not (add_offset or scale_factor):
                continue

            valid = True
            reasoning = []

            # if only one of these attributes is defined, assume they
            # are the same type (value doesn't matter here)
            if not add_offset:
                add_offset = scale_factor
            if not scale_factor:
                scale_factor = add_offset

            if type(add_offset) != type(scale_factor):
                valid = False
                reasoning.append("Attributes add_offset and scale_factor have different data type.")
            elif type(scale_factor) != var.dtype:
                # Check both attributes are type float or double
                if not type(scale_factor) in [np.float, np.float16, np.float32, np.float64, np.float128]:
                    valid = False
                    reasoning.append("Attributes add_offset and scale_factor are not of type float or double.")
                else:
                    # Check variable type is byte, short or int
                    if var.dtype not in [np.int, np.int8, np.int16, np.int32, np.int64]:
                        valid = False
                        reasoning.append("Variable is not of type byte, short, or int.")

            result = Result(BaseCheck.MEDIUM, valid, ('§8.1 Packed Data', name, 'packed_data'), reasoning)
            ret_val.append(result)
            reasoning = []

            valid = True
            # test further with  _FillValue , valid_min , valid_max , valid_range
            if hasattr(var, "_FillValue"):
                if var._FillValue.dtype != var.dtype:
                    valid = False
                    reasoning.append("Type of _FillValue attribute (%s) does not match variable type (%s)" %
                                     (var._FillValue.dtype, var.dtype))
            if hasattr(var, "valid_min"):
                if var.valid_min.dtype != var.dtype:
                    valid = False
                    reasoning.append("Type of valid_min attribute (%s) does not match variable type (%s)" %
                                     (var.valid_min.dtype, var.dtype))
            if hasattr(var, "valid_max"):
                if var.valid_max.dtype != var.dtype:
                    valid = False
                    reasoning.append("Type of valid_max attribute (%s) does not match variable type (%s)" %
                                     (var.valid_max.dtype, var.dtype))
            if hasattr(var, "valid_range"):
                if var.valid_range.dtype != var.dtype:
                    valid = False
                    reasoning.append("Type of valid_range attribute (%s) does not match variable type (%s)" %
                                     (var.valid_range.dtype, var.dtype))

            result = Result(BaseCheck.MEDIUM, valid, ('§8.1 Packed Data', name, 'fillvalue_valid_range_attributes'), reasoning)
            ret_val.append(result)

        return ret_val

    def check_compression(self, ds):
        """
        8.2 To save space in the netCDF file, it may be desirable to eliminate points from data arrays that are invariably
        missing. Such a compression can operate over one or more adjacent axes, and is accomplished with reference to a list
        of the points to be stored.

        The list is stored as the coordinate variable for the compressed axis of the data array. Thus, the list variable and
        its dimension have the same name. The list variable has a string attribute compress, containing a blank-separated
        list of the dimensions which were affected by the compression in the order of the CDL declaration of the uncompressed
        array.
        """
        ret_val = []

        for name, var in ds.variables.items():
            valid_dim = 0
            valid_form = 0
            reasoning = []
            if hasattr(var, 'compress'):
                totals = 2
                if name in var.dimensions and var.ndim == 1:
                    valid_dim = 1
                else:
                    reasoning.append("The 'compress' attribute is not assigned to a coordinate variable.")
                if all([each in list(ds.dimensions.keys()) for each in getattr(var, 'compress', '').split(" ")]):
                    valid_form = 1
                else:
                    reasoning.append("The 'compress' attribute is not in the form of a coordinate.")

                result = Result(BaseCheck.MEDIUM,
                                (valid_form + valid_dim, totals),
                                ('§8.2 Dataset Compression', name, 'compressed_data'),
                                reasoning)
                ret_val.append(result)

        return ret_val
    ###############################################################################
    #
    # Chapter 9: Discrete Sampling Geometries
    #
    ###############################################################################

    @is_likely_dsg
    def check_all_features_are_same_type(self, ds):
        """
        9.1 The features contained within a collection must always be of the same type; and all the collections in a CF file
        must be of the same feature type.

        point, timeSeries, trajectory, profile, timeSeriesProfile, trajectoryProfile.

        The space-time coordinates that are indicated for each feature are mandatory.  However a featureType may also include
        other space-time coordinates which are not mandatory (notably the z coordinate).
        """
        flag = 0
        x = ''
        y = ''
        z = ''
        t = ''

        flag = 0
        for coord_name in self._find_coord_vars(ds):
            coord_var = ds.variables[coord_name]
            if getattr(coord_var, "grid_mapping_name", ""):
                # DO GRIDMAPPING CHECKS FOR X,Y,Z,T
                flag = 1
                for name_again, var_again in ds.variables.items():
                    if getattr(var_again, "standard_name", "") == self.grid_mapping_dict[getattr(coord_var, "grid_mapping_name", "")][2][0]:
                        x = name_again
                    if getattr(var_again, "standard_name", "") == self.grid_mapping_dict[getattr(coord_var, "grid_mapping_name", "")][2][1]:
                        y = name_again

        for coord_name in self._find_coord_vars(ds):
            coord_var = ds.variables[coord_name]
            # DO STANDARD SEARCH
            if getattr(coord_name, 'units', '').lower() in ['pa', 'kpa', 'mbar', 'bar', 'atm', 'hpa', 'dbar'] or getattr(coord_name, 'positive', '') or getattr(coord_name, 'standard_name', '') == 'z' or getattr(coord_var, 'axis', '') == 'z':
                z = coord_name
            if coord_name.lower() in ['lon', 'longitude'] and flag == 0:
                x = coord_name
            elif coord_name.lower()in ['lat', 'latitude'] and flag == 0:
                y = coord_name
            elif coord_name.lower() == 'time':
                t = coord_name

            if getattr(coord_var, '_CoordinateAxisType', ''):
                axis_type = getattr(coord_var, '_CoordinateAxisType', '')
                if axis_type.lower() in ['lon', 'longitude'] and flag == 0:
                    x = coord_name
                elif axis_type.lower()in ['lat', 'latitude'] and flag == 0:
                    y = coord_name
                elif axis_type.lower() == 'time':
                    t = coord_name

        valid = False
        feature_tuple_list = []

        # create shape size tuple
        if x == '' or y == '' or t == '':
            return
        elif z == '':
            feature_tuple = (ds.variables[x].ndim, ds.variables[y].ndim, ds.variables[t].ndim)
        else:
            feature_tuple = (ds.variables[x].ndim, ds.variables[y].ndim, ds.variables[t].ndim, ds.variables[z].ndim)

        feature_tuple_list.append(feature_tuple)

        data_vars = self._find_geophysical_vars(ds)

        feature_map = {}
        for var_name in data_vars:
            variable = ds.variables[var_name]
            feature = variable.dimensions
            feature_map[var_name] = feature

        features = list(feature_map.values())
        valid = all((features[0] == feature for feature in features))
        reasoning = []
        if not valid:
            reasoning.append("At least one of the variables has a different feature type than the rest of the variables.")
            feature_mess = []
            for var_name, feature in feature_map.items():
                feature_mess.append("%s(%s)" % (var_name, ', '.join(feature) ))
            reasoning.append(' '.join(feature_mess))

        return Result(BaseCheck.HIGH, valid, '§9.1 Feature Types are all the same', reasoning)

    @is_likely_dsg
    def check_orthogonal_multidim_array(self, ds):
        """
        9.3.1 The orthogonal multidimensional array representation, the simplest representation, can be used if each feature
        instance in the collection has identical coordinates along the element axis of the features.

        Data variables have both an instance dimension and an element dimension.  The dimensions may be given in any order.
        If there is a need for either the instance or an element dimension to be the netCDF unlimited dimension (so that more
        features or more elements can be appended), then that dimension must be the outer dimension of the data variable
        i.e. the leading dimension in CDL.
        """
        ret_val = []
        reasoning = []

        for name, var in ds.variables.items():
            reasoning = []
            if not hasattr(var, 'count_variable') and not hasattr(var, 'index_variable'):
                if hasattr(var, '_FillValue'):
                    pass
                else:
                    result = Result(BaseCheck.MEDIUM,
                                    True,
                                    ('§9.3.1 Orthogonal Multidimensional Array', name, 'orthogonal_multidimensional'),
                                    reasoning)
                    ret_val.append(result)
        return ret_val

    @is_likely_dsg
    def check_incomplete_multidim_array(self, ds):
        """
        9.3.2 The incomplete multidimensional array representation can used if the features within a collection do not all have
        the same number of elements, but sufficient storage space is available to allocate the number of elements required by
        the longest feature to all features.  That is, features that are shorter than the longest feature must be padded with
        missing values to bring all instances to the same storage size.

        Data variables have both an instance dimension and an element dimension.  The dimensions may be given in any order.
        If there is a need for either the instance or an element dimension to be the netCDF unlimited dimension (so that more
        features or more elements can be appended), thlen that dimension must be the outer dimension of the data variable
        i.e. the leading dimension in CDL.
        """

        ret_val = []
        for name, var in ds.variables.items():
            reasoning = []
            if not hasattr(var, 'count_variable') and not hasattr(var, 'index_variable'):
                if hasattr(var, '_FillValue'):
                    result = Result(BaseCheck.MEDIUM,
                                    True,
                                    ('§9.3.2 Incomplete Multidimensional Array', name, 'ragged_multidimensional'),
                                    reasoning)
                    ret_val.append(result)
                else:
                    pass

        return ret_val

    @is_likely_dsg
    def check_contiguous_ragged_array(self, ds):
        """
        9.3.3 The contiguous ragged array representation can be used only if the size of each feature is known at the time
        that it is created.

        In this representation, the file contains a count variable, which must be of type integer and must have the instance
        dimension as its sole dimension.  The count variable contains the number of elements that each feature has. This
        representation and its count variable are identifiable by the presence of an attribute, sample_dimension, found on
        the count variable, which names the sample dimension being counted. For indices that correspond to features, whose
        data have not yet been written, the count variable should  have a value of zero or a missing value.

        In the ragged array representations, the instance dimension (i), which sequences the individual features within the
        collection, and the element dimension, which sequences the data elements of each feature (o and p), both occupy the
        same dimension (the sample dimension).   If the sample dimension is the netCDF unlimited dimension, new data can be
        appended to the file.
        """
        ret_val = []
        reasoning = []
        for name, var in ds.variables.items():
            if getattr(var, 'sample_dimension', ''):
                result = Result(BaseCheck.MEDIUM,
                                True,
                                ('§9.3.3 Continuous ragged array', name, 'continuous_ragged'),
                                reasoning)
                ret_val.append(result)
            else:
                return []

        return ret_val

    @is_likely_dsg
    def check_indexed_ragged_array(self, ds):
        """
        9.3.4 The indexed ragged array representation stores the features interleaved along the sample dimension in the data
        variable.

        In this representation, the file contains an index variable, which must be of type integer, and must have the sample
        dimension as its single dimension. The index variable contains the zero-based index of the feature to which each
        element belongs. This representation is identifiable by the presence of an attribute, instance_dimension, on the index
        variable, which names the dimension of the instance variables. For those indices of the sample dimension, into which
        data have not yet been written, the index variable should be pre-filled with missing values.

        In the ragged array representations, the instance dimension (i), which sequences the individual features within the
        collection, and the element dimension, which sequences the data elements of each feature (o and p), both occupy the
        same dimension (the sample dimension).   If the sample dimension is the netCDF unlimited dimension, new data can be
        appended to the file.
        """
        ret_val = []
        reasoning = []
        for name, var in ds.variables.items():
            if getattr(var, 'instance_dimension', ''):
                result = Result(BaseCheck.MEDIUM,
                                True,
                                ('§9.3.4 Indexed ragged array', name, 'continuous_ragged'),
                                reasoning)
                ret_val.append(result)
            else:
                return []

        return ret_val

    @is_likely_dsg
    def check_feature_type(self, ds):
        """
        9.4 A global attribute, featureType, is required for all Discrete Geometry representations except the orthogonal
        multidimensional array representation, for which it is highly recommended.

        The value assigned to the featureType attribute is case-insensitive.
        """
        reasoning = []
        feature_list = ['point', 'timeseries', 'trajectory', 'profile', 'timeseriesprofile', 'trajectoryprofile']

        if getattr(ds, 'featureType', '').lower() in feature_list:
            return Result(BaseCheck.MEDIUM,
                          True, '§9.4 featureType attribute',
                          reasoning)

        elif getattr(ds, 'featureType', ''):
            reasoning.append('The featureType is provided and is not from the featureType list.')
            return Result(BaseCheck.MEDIUM,
                          False, '§9.4 featureType attribute',
                          reasoning)

    @is_likely_dsg
    def check_coordinates_and_metadata(self, ds):
        """
        9.5 Every feature within a Discrete Geometry CF file must be unambiguously associated with an extensible collection
        of instance variables that identify the feature and provide other metadata as needed to describe it.  Every element of
        every feature must be unambiguously associated with its space and time coordinates and with the feature that contains
        it.

        The coordinates attribute must be attached to every data variable to indicate the spatiotemporal coordinate variables
        that are needed to geo-locate the data.


        Where feasible a variable with the attribute cf_role should be included.  The only acceptable values of cf_role for
        Discrete Geometry CF data sets are timeseries_id, profile_id, and trajectory_id.   The variable carrying the cf_role
        attribute may have any data type.  When a variable is assigned this attribute, it must provide a unique identifier
        for each feature instance.

        CF files that contain timeSeries, profile or trajectory featureTypes, should include only a single occurrence of a
        cf_role attribute;  CF files that contain timeSeriesProfile or trajectoryProfile may contain two occurrences,
        corresponding to the two levels of structure in these feature types.

        CF Discrete Geometries provides a mechanism to encode both the nominal and the precise positions, while retaining the
        semantics of the idealized feature type. Only the set of coordinates which are regarded as the nominal (default or
        preferred) positions should be indicated by the attribute axis, which should be assigned string values to indicate
        the orientations of the axes (X, Y, Z, or T).

        Auxiliary coordinate variables containing the nominal and the precise positions should be listed in the relevant
        coordinates attributes of data variables. In orthogonal representations the nominal positions could be  coordinate
        variables, which do not need to be listed in the coordinates attribute, rather than auxiliary coordinate variables.

        Coordinate bounds may optionally be associated with coordinate variables and auxiliary coordinate variables using
        the bounds attribute.

        If there is a vertical coordinate variable or auxiliary coordinate variable, it must be identified by the means
        specified in section 4.3.   The use of the attribute axis=Z is recommended for clarity.  A standard_name attribute
        that identifies the vertical coordinate is recommended.
        """
        ret_val = []
        reasoning = []

        non_data_list = []

        # Build a list of variables that are not geophysical variables
        for name, var in ds.variables.items():
            if name in self._find_coord_vars(ds):
                non_data_list.append(name)
            elif name in self._find_ancillary_vars(ds):
                non_data_list.append(name)

            # If there's a variable that is referenced in another variable's
            # coordinates attribute, put that on the non-data list as well.
            elif hasattr(var, 'coordinates'):
                for coord in getattr(var, 'coordinates', '').split(' '):
                    if coord in ds.variables:
                        non_data_list.append(coord)

            # If the variable has a dimension and is not a station identifier
            elif var.dimensions == tuple() or var.dimensions == (name,):
                non_data_list.append(name)

            elif self._is_station_var(var):
                non_data_list.append(name)

            elif hasattr(var, 'cf_role'):
                non_data_list.append(name)

        for data_entry in ds.variables:
            # If the variable is not a geophysical variable, skip checking it
            if data_entry in non_data_list:
                continue

            # If the variable has a coordinates attribute, then it passes
            if getattr(ds.variables[data_entry], 'coordinates', ''):
                result = Result(BaseCheck.MEDIUM,
                                True,
                                ('§9.5 Discrete Geometry', data_entry, 'check_coordinates'),
                                reasoning)
                ret_val.append(result)
                reasoning = []
            else:
                reasoning.append('The variable {} does not have associated coordinates'.format(data_entry))
                result = Result(BaseCheck.MEDIUM,
                                False,
                                ('§9.5 Discrete Geometry', data_entry, 'check_coordinates'),
                                reasoning)
                ret_val.append(result)
                reasoning = []

        role_list = [getattr(var, 'cf_role', '').split(' ') for name, var in ds.variables.items() if hasattr(var, 'cf_role')]
        single_role = ['timeseries', 'profile', 'trajectory']
        dual_role = ['timeseries', 'profile', 'trajectory', 'timeseriesprofile', 'trajectoryprofile']
        if getattr(ds, 'featureType', '').lower() in single_role and len(np.ravel(role_list)) == 1:
            reasoning = []
            valid = True
        elif getattr(ds, 'featureType', '').lower() in dual_role and len(np.ravel(role_list)) in [1, 2]:
            reasoning = []
            valid = True
        else:
            valid = False
            reasoning = []
            reasoning.append('The cf_role featureType is not properly defined.')

        result = Result(BaseCheck.MEDIUM,
                        valid,
                        ('§9.5 Discrete Geometry', 'dataset'),
                        reasoning)
        ret_val.append(result)

        return ret_val

    @is_likely_dsg
    def check_missing_data(self, ds):
        """
        9.6 Auxiliary coordinate variables (spatial and time) must contain missing values to indicate a void in data storage
        in the file but must not have missing data for any other reason.

        It is not permitted for auxiliary coordinate variables to have missing values for elements where there is non-missing
        data. Where any auxiliary coordinate variable contains a missing value, all other coordinate, auxiliary coordinate
        and data values corresponding to that element should also contain missing values. Data variables should (as usual)
        also contain missing values to indicate when there is no valid data available for the element, although the
        coordinates are valid.

        Similarly, for indices where the instance variable identified by cf_role contains a missing value indicator, all other
        instance variable should also contain missing values.
        """

        # Data intensive check; consider flagging as optional
        ret_val = []

        name_list = set(ds.variables.keys())
        dim_list = set(ds.dimensions.keys())

        for name, var in ds.variables.items():
            if hasattr(var, 'coordinates'):
                aux_index_dict = {}
                dim_index_dict = {}
                reasoning = []
                valid = False
                aux_valid = False

                # check if _FillValue attribute present attribute
                if hasattr(var, '_FillValue'):
                    dim_index_dict[name] = dict()
                    aux_index_dict[name] = dict()
                    for coordinate in getattr(var, 'coordinates', '').split(" "):
                        indices = np.where(ds.variables[coordinate][:].mask)[0]

                        # if the coordinate name is in the list of variables,
                        # but has no associated dimension name
                        if coordinate in name_list and \
                           coordinate not in dim_list:

                            dim_index_dict[name][coordinate] = indices
                            aux_index_dict[name][coordinate] = indices

                        # if this is a coordinate variable, i.e. a variable
                        # with the same variable and coordinate name
                        elif coordinate in name_list and coordinate in dim_list:
                            dim_index_dict[name][coordinate] = indices
                        else:
                            dim_index_dict[name][coordinate] = np.array()
                # Check to see that all coordinate variable mising data
                # locations are the same

                # could be refactored to make much more efficient
                def check_index_dicts(idx_dict):
                    """
                    Checks that all indices align and are of the same value
                    """
                    # Base case: In the even there are no arrays to check then
                    # automatically valid
                    valid = True
                    idx_gen = (d[k] for d in six.itervalues(dim_index_dict)
                               for k in d if d[k].size > 0)
                    first_idx = next(idx_gen, None)
                    if first_idx is not None:
                        # if there is only one element, all will return True,
                        # which is also correct
                        valid = all(np.array_equal(val, first_idx)
                                    for val in idx_gen)
                    return valid

                aux_valid = check_index_dicts(aux_index_dict)

                valid = check_index_dicts(dim_index_dict)

                if not aux_valid:
                    reasoning.append('The auxillary coordinates do not have the same missing data locations')
                # dimensions not valid
                if not valid:
                    reasoning.append('The coordinate variables do not have the same missing data locations as the auxillary coordinates')

                # Check to see that all coordinate variable mising data is reflceted in the dataset
                valid_missing = True

                if hasattr(var, '_FillValue'):
                    x_indices = np.where(var[:].mask)[0]

                    for coordinate in var.coordinates.split(" "):
                        coordinate_ind_list = dim_index_dict[name][coordinate]
                        valid_missing = np.array_equal(x_indices,
                                                       coordinate_ind_list)
                    if valid_missing is False:
                        reasoning.append('The data does not have the same missing data locations as the coordinates')

                result = Result(BaseCheck.MEDIUM,
                                valid and aux_valid and valid_missing,
                                ('§9.6 Missing Data', name, 'missing_data'),
                                reasoning)
                ret_val.append(result)
        return ret_val

    def _find_container_variables(self, ds):
        container_vars = []
        platform_name = getattr(ds, 'platform', None)
        if platform_name is not None:
            container_vars.append(platform_name)
        for k, v in ds.variables.items():
            if k in ('crs', 'instrument', 'station'):
                if not v.shape:  # Empty dimension
                    container_vars.append(k)
            platform_name = getattr(v, 'platform', None)
            if platform_name is not None:
                container_vars.append(platform_name)
            instrument_name = getattr(v, 'instrument_name', None)
            if instrument_name is not None:
                container_vars.append(instrument_name)
        return list(set(container_vars))


class CFNCCheck(BaseNCCheck, CFBaseCheck):

    @classmethod
    def beliefs(cls):  # @TODO
        return {}
