#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals, division
from compliance_checker.base import BaseCheck, BaseNCCheck, Result, TestCtx
from compliance_checker.cf.appendix_d import dimless_vertical_coordinates
from compliance_checker.cf.appendix_f import grid_mapping_dict
from compliance_checker.cf import util
from compliance_checker import cfutil
from functools import wraps
from collections import defaultdict
import numpy as np
import os
import re

import logging

logger = logging.getLogger(__name__)


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

        self._std_names        = util.StandardNameTable()

    ################################################################################
    #
    # Helper Methods - var classifications, etc
    #
    ################################################################################

    def setup(self, ds):
        """
        Initialize various special variable types within the class.
        Mutates a number of instance variables.

        :param netCDF4.Dataset ds: An open netCDF dataset
        """
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
        version of the cf standard name table.  If the standard name table has
        already been downloaded, use the cached version.  Modifies `_std_names`
        attribute to store standard names.  Returns True if the file exists and
        False if it fails to download.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: bool
        '''
        # Get the standard name vocab
        standard_name_vocabulary = getattr(ds, 'standard_name_vocabulary', '')

        # Try to parse this attribute to get version
        version = None
        if 'cf standard name table' in standard_name_vocabulary.lower():
            version = standard_name_vocabulary.split()[-1]
        else:
            # Can't parse the attribute, use the packaged version
            return False

        if version.startswith('v'):  # i.e 'v34' -> '34' drop the v
            version = version[1:]

        # If the packaged version is what we're after, then we're good
        if version == self._std_names._version:
            print("Using packaged standard name table v{0}".format(version))
            return False

        # Try to download the version specified
        try:
            data_directory = util.create_cached_data_dir()
            location = os.path.join(data_directory, 'cf-standard-name-table-test-{0}.xml'.format(version))
            # Did we already download this before?
            if not os.path.isfile(location):
                util.download_cf_standard_name_table(version, location)
                print("Using downloaded standard name table v{0}".format(version))
            else:
                print("Using cached standard name table v{0} from {1}".format(version, location))

            self._std_names = util.StandardNameTable(location)
            return True
        except Exception:
            # There was an error downloading the CF table. That's ok, we'll just use the packaged version
            print("Error fetching standard name table. Using packaged v{0}".format(self._std_names._version))
            return False

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
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        :rtype: list
        :return:   List of variable names (str) that are likely metadata
                   variable candidates.

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
        Returns a list of geophysical variables.  Modifies
        `self._geophysical_vars`

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param bool refresh: if refresh is set to True, the cache is
                             invalidated.
        :rtype: list
        :return: A list containing strings with geophysical variable
                 names.
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
        :rtype: list
        :return: A list containing strings with geophysical variable
                 names.
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
        :rtype: list
        :return: A list containing strings with boundary variable names.
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
        :rtype: compliance_checker.base.Result
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
        :rtype: compliance_checker.base.Result
        '''
        ret_val = []
        variable_naming = TestCtx(BaseCheck.MEDIUM, '§2.3 Naming Conventions for variables')
        dimension_naming = TestCtx(BaseCheck.MEDIUM, '§2.3 Naming Conventions for dimensions')
        attribute_naming = TestCtx(BaseCheck.MEDIUM, '§2.3 Naming Conventions for attributes')

        ignore_attributes = [
            '_FillValue',
            'DODS',
            '_ChunkSizes',
            '_Coordinate',
            '_Unsigned'
        ]

        rname = re.compile("^[A-Za-z][A-Za-z0-9_]*$")

        for name, variable in ds.variables.items():
            variable_naming.assert_true(rname.match(name) is not None,
                                        "variable {} should begin with a letter and be composed of "
                                        "letters, digits, and underscores".format(name))

            # Keep track of all the attributes, we'll need to check them
            for attr in variable.ncattrs():
                if attr in ignore_attributes:
                    continue
                # Special attributes made by THREDDS
                if attr.startswith('DODS'):
                    continue
                # Ignore model produced attributes
                if attr.startswith('_Coordinate'):
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
            if global_attr.startswith('DODS'):
                continue
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
        :rtype: compliance_checker.base.Result
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
        :rtype: compliance_checker.base.Result
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

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        '''
        valid_dimension_order = TestCtx(BaseCheck.MEDIUM, '§2.4 Dimension Order')
        # Build a map from coordinate variable to axis
        coord_axis_map = self._get_coord_axis_map(ds)

        # Check each variable's dimension order, excluding climatology and
        # bounds variables
        any_clim = cfutil.get_climatology_variable(ds)
        any_bounds = cfutil.get_cell_boundary_variables(ds)
        for name, variable in ds.variables.items():
            # Skip bounds/climatology variables, as they should implicitly
            # have the same order except for the bounds specific dimension.
            # This is tested later in the respective checks
            if name in any_bounds or name == any_clim:
                continue

            # Skip strings/labels
            if hasattr(variable.dtype, 'char') and variable.dtype.char == 'S':
                continue
            elif variable.dtype == str:
                continue

            if variable.dimensions:
                dimension_order = self._get_dimension_order(ds, name, coord_axis_map)
                valid_dimension_order.assert_true(self._dims_in_order(dimension_order),
                                                  "{}'s dimensions are not in the recommended order "
                                                  "T, Z, Y, X. They are {}"
                                                  "".format(name, self._get_pretty_dimension_order(ds, name)))

        return valid_dimension_order.to_result()

    def _get_coord_axis_map(self, ds):
        '''
        Returns a dictionary mapping each coordinate to a letter identifier
        describing the _kind_ of coordinate.

        :param netCDF4.Dataset ds: An open netCDF dataset

        :rtype: dict
        :return: A dictionary with variable names mapped to axis abbreviations,
                 i.e. {'longitude': 'X', ... 'pressure': 'Z'}
        '''
        expected = ['T', 'Z', 'Y', 'X']
        coord_vars = self._find_coord_vars(ds)
        coord_axis_map = {}

        # L - Unlimited Coordinates
        # T - Time coordinates
        # Z - Depth/Altitude Coordinate
        # Y - Y-Coordinate (latitude)
        # X - X-Coordinate (longitude)
        # A - Auxiliary Coordinate
        # I - Instance Coordinate

        time_variables = cfutil.get_time_variables(ds)
        lat_variables = cfutil.get_latitude_variables(ds)
        lon_variables = cfutil.get_longitude_variables(ds)
        z_variables = cfutil.get_z_variables(ds)

        for coord_name in coord_vars:
            coord_var = ds.variables[coord_name]
            axis = getattr(coord_var, 'axis', None)
            standard_name = getattr(coord_var, 'standard_name', None)

            # Unlimited dimensions must come first
            if ds.dimensions[coord_name].isunlimited():
                coord_axis_map[coord_name] = 'L'
            # axis takes precedence over standard_name
            elif axis in expected:
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
            elif coord_name in time_variables:
                coord_axis_map[coord_name] = 'T'
            elif coord_name in z_variables:
                coord_axis_map[coord_name] = 'Z'
            elif coord_name in lat_variables:
                coord_axis_map[coord_name] = 'Y'
            elif coord_name in lon_variables:
                coord_axis_map[coord_name] = 'X'
            else:
                # mark the coordinate variable as unknown
                coord_axis_map[coord_name] = 'U'

        for dimension in self._get_instance_dimensions(ds):
            if dimension not in coord_axis_map:
                coord_axis_map[dimension] = 'I'

        # Dimensions of auxiliary coordinate variables will be marked with A.
        # This is useful to help determine if the dimensions are used like a
        # mapping from grid coordinates to physical lat/lon
        for coord_name in self._find_aux_coord_vars(ds):
            coord_var = ds.variables[coord_name]
            # Skip label auxiliary coordinates
            if coord_var.dtype.char == 'S':
                continue
            for dimension in coord_var.dimensions:
                if dimension not in coord_axis_map:
                    coord_axis_map[dimension] = 'A'

        # If a dimension does not have a coordinate variable mark it as unknown
        # 'U'
        for dimension in ds.dimensions:
            if dimension not in coord_axis_map:
                coord_axis_map[dimension] = 'U'

        return coord_axis_map

    def _get_instance_dimensions(self, ds):
        '''
        Returns a list of dimensions marked as instance dimensions

        :param netCDF4.Dataset ds: An open netCDF dataset

        :rtype: list
        :returns: A list of variable dimensions
        '''
        ret_val = []
        for variable in ds.get_variables_by_attributes(cf_role=lambda x: isinstance(x, basestring)):
            if variable.ndim > 0:
                ret_val.append(variable.dimensions[0])
        return ret_val

    def _get_pretty_dimension_order(self, ds, name):
        '''
        Returns a comma seperated string of the dimensions for a specified
        variable

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: A string with a valid NetCDF variable name for the
                         dataset
        :rtype: str
        :return: A comma separated string of the variable's dimensions
        '''
        dim_names = []
        for dim in ds.variables[name].dimensions:
            dim_name = dim
            if ds.dimensions[dim].isunlimited():
                dim_name += ' (Unlimited)'
            dim_names.append(dim_name)
        return ', '.join(dim_names)

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
        :return: A list of strings corresponding to the named axis of the dimensions for a variable
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
        :rtype: bool
        :return: Returns True if the dimensions are in order U*, T, Z, Y, X,
                 False otherwise
        '''
        regx = re.compile(r'^L?I?U*T?Z?(?:(?:Y?X?)|(?:C?)|(?:A+))$')
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
        :return: List of Results
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
                if isinstance(variable.valid_range, basestring):
                    valid_fill_range.assert_true(False, '{}:valid_range must be a numeric type not a string'.format(name))
                    continue
                rmin, rmax = variable.valid_range
                spec_by = 'valid_range'

            elif 'valid_min' in attrs and 'valid_max' in attrs:
                if isinstance(variable.valid_min, basestring):
                    valid_fill_range.assert_true(False, '{}:valid_min must be a numeric type not a string'.format(name))
                if isinstance(variable.valid_max, basestring):
                    valid_fill_range.assert_true(False, '{}:valid_max must be a numeric type not a string'.format(name))
                if isinstance(variable.valid_min, basestring) or \
                   isinstance(variable.valid_max, basestring):
                    continue
                rmin = variable.valid_min
                rmax = variable.valid_max
                spec_by = 'valid_min/valid_max'
            else:
                continue

            if np.isnan(fill_value):
                valid = True
            else:
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
        :rtype: compliance_checker.base.Result
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
        :return: List of Results
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

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of Results
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
        - if standard name specified, must be consistent with standard name table, must also be consistent with a
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

        for name in set(unit_required_variables):
            # For reduced horizontal grids, the compression index variable does
            # not require units.
            if cfutil.is_compression_coordinate(ds, name):
                continue

            variable = ds.variables[name]

            # Skip instance coordinate variables
            if getattr(variable, 'cf_role', None) is not None:
                continue

            # Skip labels
            if variable.dtype.char == 'S':
                continue

            standard_name = getattr(variable, 'standard_name', None)
            standard_name, standard_name_modifier = self._split_standard_name(standard_name)

            units = getattr(variable, 'units', None)

            valid_units = self._check_valid_cf_units(ds, name)
            ret_val.append(valid_units)

            if isinstance(units, basestring):
                valid_udunits = self._check_valid_udunits(ds, name)
                ret_val.append(valid_udunits)

            if isinstance(standard_name, basestring):
                valid_standard_units = self._check_valid_standard_units(ds, name)
                ret_val.append(valid_standard_units)

        return ret_val

    def _split_standard_name(self, standard_name):
        '''
        Returns a tuple of the standard_name and standard_name modifier

        Nones are used to represent the absence of a modifier or standard_name

        :rtype: tuple
        :return: 2-tuple of standard_name and modifier as strings
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
        :rtype:
        :return: List of results
        '''
        # This list is straight from section 3
        deprecated = ['level', 'layer', 'sigma_level']
        variable = ds.variables[variable_name]

        units = getattr(variable, 'units', None)
        standard_name_full = getattr(variable, 'standard_name', None)
        standard_name, standard_name_modifier = self._split_standard_name(standard_name_full)
        std_name_unitless = cfutil.get_unitless_standard_names(self._std_names._root,
                                                               standard_name)
        # Is this even in the database? also, if there is no standard_name,
        # there's no way to know if it is unitless.
        should_be_unitless = (variable.ndim == 0 or variable.dtype.char == 'S' or
                              std_name_unitless or standard_name is None)

        # 1) Units must exist
        valid_units = TestCtx(BaseCheck.HIGH, '§3.1 Variable {} contains valid CF units'.format(variable_name))
        valid_units.assert_true(should_be_unitless or units is not None,
                                'units attribute is required for {}'.format(variable_name))

        # Don't bother checking the rest
        if units is None and not should_be_unitless:
            return valid_units.to_result()
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
        std_name_unitless = cfutil.get_unitless_standard_names(self._std_names._root,
                                                               standard_name)

        # If the variable is supposed to be unitless, it automatically passes
        should_be_unitless = (variable.ndim == 0 or variable.dtype.char == 'S'
                              or std_name_unitless)

        valid_udunits = TestCtx(BaseCheck.LOW,
                                "§3.1 Variable {}'s units are contained in UDUnits".format(variable_name))
        are_udunits = (units is not None and util.units_known(units))
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

        valid_standard_units = TestCtx(BaseCheck.HIGH,
                                       "§3.1 Variable {}'s units are appropriate for "
                                       "the standard_name {}".format(variable_name,
                                                                     standard_name or "unspecified"))

        # If the variable is supposed to be unitless, it automatically passes
        std_name_unitless = cfutil.get_unitless_standard_names(self._std_names._root,
                                                               standard_name)

        standard_name, standard_name_modifier = self._split_standard_name(standard_name)

        standard_entry = self._std_names.get(standard_name, None)
        if standard_entry is not None:
            canonical_units = standard_entry.canonical_units
        else:
            # Any unit comparisons with None returns False
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
            valid_standard_units.assert_true(util.units_convertible(units, 'seconds since 1970-01-01'),
                                             'time must be in a valid units format <unit> since <epoch> '
                                             'not {}'.format(units))

        # UDunits can't tell the difference between east and north facing coordinates
        elif standard_name == 'latitude':
            # degrees is allowed if using a transformed grid
            allowed_units = cfutil.VALID_LAT_UNITS + ['degrees']
            valid_standard_units.assert_true(units.lower() in allowed_units,
                                             'variables defining latitude must use degrees_north '
                                             'or degrees if defining a transformed grid. Currently '
                                             '{}'.format(units))
        # UDunits can't tell the difference between east and north facing coordinates
        elif standard_name == 'longitude':
            # degrees is allowed if using a transformed grid
            allowed_units = cfutil.VALID_LON_UNITS + ['degrees']
            valid_standard_units.assert_true(units.lower() in allowed_units,
                                             'variables defining longitude must use degrees_east '
                                             'or degrees if defining a transformed grid. Currently '
                                             '{}'.format(units))
        # Standard Name table agrees the unit should be unitless
        elif std_name_unitless:
            valid_standard_units.assert_true(True, '')

        elif canonical_units is not None:
            valid_standard_units.assert_true(util.units_convertible(canonical_units, units),
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
        for name in set(variables_requiring_standard_names):
            # Compression indices used in reduced horizontal grids or
            # compression schemes do not require attributes other than compress
            if cfutil.is_compression_coordinate(ds, name):
                continue

            ncvar = ds.variables[name]

            # §9 doesn't explicitly allow instance variables as coordinates but
            # it's loosely implied. Just in case, skip it.
            if hasattr(ncvar, 'cf_role'):
                continue

            # Unfortunately, §6.1 allows for string types to be listed as
            # coordinates.
            if ncvar.dtype.char == 'S':
                continue

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
                                       "variable {}'s attribute standard_name must be a non-empty string "
                                       "or it should define a long_name attribute.".format(name))

            if isinstance(standard_name, basestring):
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
        :rtype: compliance_checker.base.Result
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
        :rtype: compliance_checker.base.Result
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
        :rtype: compliance_checker.base.Result
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

        for variable in ds.get_variables_by_attributes(axis=lambda x: x is not None):
            name = variable.name
            # Coordinate compressions should not be checked as a valid
            # coordinate, which they are not. They are a mechanism to project
            # an array of indices onto a 2-d grid containing valid coordinates.
            if cfutil.is_compression_coordinate(ds, name):
                continue

            variable = ds.variables[name]
            # Even though it's not allowed in CF 1.6, it is allowed in CF 1.7
            # and we see people do it, often.
            if hasattr(variable, 'cf_role'):
                continue

            # §6.1 allows for labels to be referenced as auxiliary coordinate
            # variables, which should not be checked like the rest of the
            # coordinates.
            if variable.dtype.char == 'S':
                continue

            axis = getattr(variable, 'axis', None)

            if axis is not None:
                valid_axis = self._check_axis(ds, name)
                ret_val.append(valid_axis)

        return ret_val

    def _check_axis(self, ds, name):
        '''
        Checks that the axis attribute is a string and an allowed value, namely
        one of 'T', 'X', 'Y', or 'Z'.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Name of the variable
        :rtype: compliance_checker.base.Result
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
            if standard_name == 'latitude' and units != 'degrees_north':
                # This is only a recommendation and we won't penalize but we
                # will include a recommended action.
                msg = ("CF recommends latitude variable '{}' to use units degrees_north"
                       "".format(latitude))
                recommended_units = Result(BaseCheck.LOW,
                                           True,
                                           '§4.1 Latitude variable {} defines units using degrees_north'.format(latitude),
                                           [msg])
                ret_val.append(recommended_units)

            y_variables = ds.get_variables_by_attributes(axis='Y')
            # Check that latitude defines either standard_name or axis
            definition = TestCtx(BaseCheck.MEDIUM, '§4.1 Latitude variable {} defines either standard_name or axis'.format(latitude))
            definition.assert_true(standard_name == 'latitude' or axis == 'Y' or y_variables != [],
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
            if standard_name == 'longitude' and units != 'degrees_east':
                # This is only a recommendation and we won't penalize but we
                # will include a recommended action.
                msg = ("CF recommends longitude variable '{}' to use units degrees_east"
                       "".format(longitude))
                recommended_units = Result(BaseCheck.LOW,
                                           True,
                                           '§4.1 Longitude variable {} defines units using degrees_east'.format(longitude),
                                           [msg])
                ret_val.append(recommended_units)

            x_variables = ds.get_variables_by_attributes(axis='X')
            # Check that longitude defines either standard_name or axis
            definition = TestCtx(BaseCheck.MEDIUM, '§4.1 Longitude variable {} defines either standard_name or axis'.format(longitude))
            definition.assert_true(standard_name == 'longitude' or axis == 'Y' or x_variables != [],
                                   "longitude variable '{}' should define standard_name='longitude' or axis='X'"
                                   "".format(longitude))
            ret_val.append(definition.to_result())

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
        z_variables = cfutil.get_z_variables(ds)
        dimless_standard_names = [name for name, regx in dimless_vertical_coordinates]
        for name in z_variables:
            variable = ds.variables[name]
            standard_name = getattr(variable, 'standard_name', None)
            units = getattr(variable, 'units', None)
            positive = getattr(variable, 'positive', None)
            # Skip the variable if it's dimensionless
            if hasattr(variable, 'formula_terms'):
                continue
            if standard_name in dimless_standard_names:
                continue

            valid_vertical_coord = TestCtx(BaseCheck.HIGH,
                                           "§4.3.1 {} is a valid vertical coordinate"
                                           "".format(name))
            valid_vertical_coord.assert_true(isinstance(units, basestring) and units,
                                             "units must be defined for vertical coordinates, there is no default")

            if not util.units_convertible('bar', units):
                valid_vertical_coord.assert_true(positive in ('up', 'down'),
                                                 "vertical coordinates not defining pressure must include "
                                                 "a positive attribute that is either 'up' or 'down'")

            # _check_valid_standard_units, part of the Chapter 3 checks,
            # already verifies that this coordinate has valid units

            ret_val.append(valid_vertical_coord.to_result())

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
        z_variables = cfutil.get_z_variables(ds)
        deprecated_units = [
            'level',
            'layer',
            'sigma_level'
        ]
        for name in z_variables:
            variable = ds.variables[name]
            standard_name = getattr(variable, 'standard_name', None)
            units = getattr(variable, 'units', None)
            formula_terms = getattr(variable, 'formula_terms', None)
            # Skip the variable if it's dimensional
            if formula_terms is None and standard_name not in dimless:
                continue

            is_not_deprecated = TestCtx(BaseCheck.LOW,
                                        "§4.3.2 {} does not contain deprecated units"
                                        "".format(name))

            is_not_deprecated.assert_true(units not in deprecated_units,
                                          "units are deprecated by CF in variable {}: {}"
                                          "".format(name, units))
            ret_val.append(is_not_deprecated.to_result())
            ret_val.append(self._check_formula_terms(ds, name))

        return ret_val

    def _check_formula_terms(self, ds, coord):
        '''
        Checks a dimensionless vertical coordinate contains valid formula_terms

        - formula_terms is a non-empty string
        - formula_terms matches regx
        - every variable defined in formula_terms exists

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        '''
        variable = ds.variables[coord]
        dimless = dict(dimless_vertical_coordinates)
        standard_name = getattr(variable, 'standard_name', None)
        formula_terms = getattr(variable, 'formula_terms', None)
        valid_formula_terms = TestCtx(BaseCheck.HIGH,
                                      '§4.3.2 {} has valid formula_terms'
                                      ''.format(coord))

        valid_formula_terms.assert_true(isinstance(formula_terms, basestring) and formula_terms,
                                        'formula_terms is a required attribute and must be a non-empty string')
        # We can't check any more
        if not formula_terms:
            return valid_formula_terms.to_result()

        valid_formula_terms.assert_true(standard_name in dimless,
                                        "unknown standard_name for dimensionless vertical coordinate: {}"
                                        "".format(standard_name))
        if standard_name not in dimless:
            return valid_formula_terms.to_result()

        regx_match = re.match(dimless[standard_name], formula_terms)
        valid_formula_terms.assert_true(regx_match is not None,
                                        "formula_terms are invalid for {}, please see appendix D of CF 1.6"
                                        "".format(standard_name))

        if regx_match is None:
            return valid_formula_terms.to_result()

        # The pattern for formula terms is always component: variable_name
        # the regex grouping always has component names in even positions and
        # the corresponding variable name in even positions.
        match_groups = regx_match.groups()

        for i in range(int(len(match_groups) / 2)):
            variable_name = match_groups[i * 2 + 1]
            valid_formula_terms.assert_true(variable_name in ds.variables,
                                            "variable {} referenced by formula_terms does not exist"
                                            "".format(variable_name))

        return valid_formula_terms.to_result()

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
        for name in cfutil.get_time_variables(ds):
            variable = ds.variables[name]
            # Has units
            has_units = hasattr(variable, 'units')
            if not has_units:
                result = Result(BaseCheck.HIGH,
                                False,
                                '§4.4 Time coordinate variable and attributes',
                                ['%s does not have units' % name])
                ret_val.append(result)
                continue
            # Correct and identifiable units
            result = Result(BaseCheck.HIGH,
                            True,
                            '§4.4 Time coordinate variable and attributes')
            ret_val.append(result)
            correct_units = util.units_temporal(variable.units)
            reasoning = None
            if not correct_units:
                reasoning = ['%s does not have correct time units' % name]
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

        # if has a calendar, check that it is within the valid values
        # otherwise no calendar is valid
        for time_var in \
            ds.get_variables_by_attributes(calendar=lambda c: c is not None):
            reasoning=None
            valid_calendar = time_var.calendar in valid_calendars

            if not valid_calendar:
                reasoning = ["Variable %s should have a valid calendar: '%s' is not a valid calendar" % (time_var.name, time_var.calendar)]

            # passes if the calendar is valid, otherwise notify of invalid
            # calendar

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
        """
        Returns True if the NetCDF variable is associated with a station, False
        otherwise.

        :param netCDF4.Variable var: a variable in an existing NetCDF dataset
        :rtype: bool
        :return: Status of whether variable appears to be associated with a
                 station
        """

        if getattr(var, 'standard_name', None) in ('platform_name', 'station_name', 'instrument_name'):
            return True
        return False

    def _get_coord_vars(self, ds):
        coord_vars = []
        for name, var in ds.variables.items():
            if (name,) == var.dimensions:
                coord_vars.append(name)
        return coord_vars

    def check_aux_coordinates(self, ds):
        '''
        Chapter 5 paragraph 3

        The dimensions of an auxiliary coordinate variable must be a subset of
        the dimensions of the variable with which the coordinate is associated,
        with two exceptions. First, string-valued coordinates (Section 6.1,
        "Labels") have a dimension for maximum string length. Second, in the
        ragged array representations of data (Chapter 9, Discrete Sampling
        Geometries), special methods are needed to connect the data and
        coordinates.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''

        ret_val = []
        geophysical_variables = self._find_geophysical_vars(ds)
        for name in geophysical_variables:
            variable = ds.variables[name]
            coordinates = getattr(variable, 'coordinates', None)
            # We use a set so we can assert
            dim_set = set(variable.dimensions)
            # No auxiliary coordinates, no check
            if not isinstance(coordinates, basestring) or coordinates == '':
                continue

            valid_aux_coords = TestCtx(BaseCheck.HIGH,
                                       "§5.0 Auxiliary Coordinates of {} must have a subset of {}'s dimensions"
                                       "".format(name, name))

            for aux_coord in coordinates.split():
                valid_aux_coords.assert_true(aux_coord in ds.variables,
                                             "auxiliary coordinate specified by the coordinates attribute, {}, "
                                             "is not a variable in this dataset"
                                             "".format(aux_coord))
                if aux_coord not in ds.variables:
                    continue

                # §6.1 Allows for "labels" to be referenced as coordinates
                if ds.variables[aux_coord].dtype.char == 'S':
                    continue

                aux_coord_dims = set(ds.variables[aux_coord].dimensions)
                valid_aux_coords.assert_true(aux_coord_dims.issubset(dim_set),
                                             "dimensions for auxiliary coordinate variable {} ({}) "
                                             "are not a subset of dimensions for variable {} ({})"
                                             "".format(aux_coord,
                                                       ', '.join(aux_coord_dims),
                                                       name,
                                                       ', '.join(dim_set)))
            ret_val.append(valid_aux_coords.to_result())
        return ret_val

    def check_duplicate_axis(self, ds):
        '''
        Checks that no variable contains two coordinates defining the same
        axis.

        Chapter 5 paragraph 6

        If an axis attribute is attached to an auxiliary coordinate variable,
        it can be used by applications in the same way the `axis` attribute
        attached to a coordinate variable is used. However, it is not
        permissible for a [geophysical variable] to have both a coordinate
        variable and an auxiliary coordinate variable, or more than one of
        either type of variable, having an `axis` attribute with any given
        value e.g. there must be no more than one axis attribute for X for any
        [geophysical variable].

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        :return: List of results
        '''

        ret_val = []
        geophysical_variables = self._find_geophysical_vars(ds)
        for name in geophysical_variables:
            no_duplicates = TestCtx(BaseCheck.HIGH, '§5.0 Variable {} does not contain duplicate coordinates'.format(name))
            axis_map = cfutil.get_axis_map(ds, name)
            axes = []
            # For every coordinate associated with this variable, keep track of
            # which coordinates define an axis and assert that there are no
            # duplicate axis attributes defined in the set of associated
            # coordinates.
            for axis, coordinates in axis_map.items():
                for coordinate in coordinates:
                    axis_attr = getattr(ds.variables[coordinate], 'axis', None)
                    no_duplicates.assert_true(axis_attr is None or axis_attr not in axes,
                                              "duplicate axis {} defined by {}".format(axis_attr, coordinate))

                    if axis_attr and axis_attr not in axes:
                        axes.append(axis_attr)

            ret_val.append(no_duplicates.to_result())

        return ret_val

    def check_multi_dimensional_coords(self, ds):
        '''
        Checks that no multidimensional coordinate shares a name with its
        dimensions.

        Chapter 5 paragraph 4

        We recommend that the name of a [multidimensional coordinate] should
        not match the name of any of its dimensions.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []

        # This can only apply to auxiliary coordinate variables
        for coord in self._find_aux_coord_vars(ds):
            variable = ds.variables[coord]
            if variable.ndim < 2:
                continue
            not_matching = TestCtx(BaseCheck.MEDIUM,
                                   '§5.0 multidimensional coordinate {} should not have the same '
                                   'name as dimension'.format(coord))

            not_matching.assert_true(coord not in variable.dimensions,
                                     '{} shares the same name as one of its dimensions'
                                     ''.format(coord))
            ret_val.append(not_matching.to_result())

        return ret_val

    def check_grid_coordinates(self, ds):
        """
        5.6 When the coordinate variables for a horizontal grid are not
        longitude and latitude, it is required that the true latitude and
        longitude coordinates be supplied via the coordinates attribute.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        latitudes = cfutil.get_true_latitude_variables(ds)
        longitudes = cfutil.get_true_longitude_variables(ds)

        check_featues = [
            '2d-regular-grid',
            '2d-static-grid',
            '3d-regular-grid',
            '3d-static-grid',
            'mapped-grid',
            'reduced-grid'
        ]

        # This one is tricky because there's a very subtle difference between
        # latitude as defined in Chapter 4 and "true" latitude as defined in
        # chapter 5.

        # For each geophysical variable that defines a grid, assert it is
        # associated with a true latitude or longitude coordinate.

        for variable in self._find_geophysical_vars(ds):
            # We use a set so we can do set-wise comparisons with coordinate
            # dimensions
            dimensions = set(ds.variables[variable].dimensions)
            # If it's not a grid, skip it
            if cfutil.guess_feature_type(ds, variable) not in check_featues:
                continue
            has_coords = TestCtx(BaseCheck.HIGH,
                                 '§5.6 Grid Feature {} is associated with true latitude and true longitude'
                                 ''.format(variable))

            # axis_map is a defaultdict(list) mapping the axis to a list of
            # coordinate names. For example:
            # {'X': ['lon'], 'Y':['lat'], 'Z':['lev']}
            # The mapping comes from the dimensions of the variable and the
            # contents of the `coordinates` attribute only.
            axis_map = cfutil.get_axis_map(ds, variable)

            # Make sure we can find latitude and it's dimensions are a subset
            found_lat = False
            for lat in axis_map['Y']:
                is_subset_dims = set(ds.variables[lat].dimensions).issubset(dimensions)

                if is_subset_dims and lat in latitudes:
                    found_lat = True
                    break
            has_coords.assert_true(found_lat,
                                   '{} is not associated with a coordinate defining true latitude '
                                   'and sharing a subset of dimensions'.format(variable))

            # Make sure we can find longitude and it's dimensions are a subset
            found_lon = False
            for lon in axis_map['X']:
                is_subset_dims = set(ds.variables[lon].dimensions).issubset(dimensions)

                if is_subset_dims and lon in longitudes:
                    found_lon = True
                    break
            has_coords.assert_true(found_lon,
                                   '{} is not associated with a coordinate defining true longitude '
                                   'and sharing a subset of dimensions'.format(variable))

            ret_val.append(has_coords.to_result())
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

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        # Create a set of coordinate varaibles defining `compress`
        lats = set(cfutil.get_latitude_variables(ds))
        lons = set(cfutil.get_longitude_variables(ds))

        for name in self._find_geophysical_vars(ds):
            coords = getattr(ds.variables[name], 'coordinates', None)
            axis_map = cfutil.get_axis_map(ds, name)
            # If this variable has no coordinate that defines compression
            if 'C' not in axis_map:
                continue

            valid_rgrid = TestCtx(BaseCheck.HIGH, '§5.3 {} is a valid reduced horizontal grid'.format(name))
            # Make sure reduced grid features define coordinates
            valid_rgrid.assert_true(isinstance(coords, basestring) and coords,
                                    "reduced grid feature {} must define coordinates attribute"
                                    "".format(name))
            # We can't check anything else if there are no defined coordinates
            if not isinstance(coords, basestring) and coords:
                continue

            coord_set = set(coords.split())

            # Make sure it's associated with valid lat and valid lon
            valid_rgrid.assert_true(len(coord_set.intersection(lons)) > 0,
                                    '{} must be associated with a valid longitude coordinate'.format(name))
            valid_rgrid.assert_true(len(coord_set.intersection(lats)) > 0,
                                    '{} must be associated with a valid latitude coordinate'.format(name))
            valid_rgrid.assert_true(len(axis_map['C']) == 1,
                                    '{} can not be associated with more than one compressed coordinates: '
                                    '({})'.format(name, ', '.join(axis_map['C'])))

            for compressed_coord in axis_map['C']:
                coord = ds.variables[compressed_coord]
                compress = getattr(coord, 'compress', None)
                valid_rgrid.assert_true(isinstance(compress, basestring) and compress,
                                        "compress attribute for compression coordinate {} must be a non-empty string"
                                        "".format(compressed_coord))
                if not isinstance(compress, basestring):
                    continue
                for dim in compress.split():
                    valid_rgrid.assert_true(dim in ds.dimensions,
                                            "dimension {} referenced by {}:compress must exist"
                                            "".format(dim, compressed_coord))
            ret_val.append(valid_rgrid.to_result())

        return ret_val

    # grid mapping dictionary, appendix F

    def check_grid_mapping(self, ds):
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

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """

        ret_val = []
        grid_mapping_variables = cfutil.get_grid_mapping_variables(ds)

        # Check the grid_mapping attribute to be a non-empty string and that it's reference exists
        for variable in ds.get_variables_by_attributes(grid_mapping=lambda x: x is not None):
            grid_mapping = getattr(variable, 'grid_mapping', None)
            defines_grid_mapping = TestCtx(BaseCheck.HIGH,
                                           "§5.6 Variable {} defining a grid mapping has valid grid_mapping attribute"
                                           "".format(variable.name))
            defines_grid_mapping.assert_true(isinstance(grid_mapping, basestring) and grid_mapping,
                                             "grid_mapping attribute must be a space-separated non-empty string")

            if isinstance(grid_mapping, basestring):
                for grid_var_name in grid_mapping.split():
                    defines_grid_mapping.assert_true(grid_var_name in ds.variables,
                                                     "grid mapping variable {} must exist in this dataset"
                                                     "".format(grid_var_name))
            ret_val.append(defines_grid_mapping.to_result())

        # Check the grid mapping variables themselves
        for grid_var_name in grid_mapping_variables:
            valid_grid_mapping = TestCtx(BaseCheck.HIGH,
                                         "§5.6 Grid Mapping Variable {} must define a valid grid mapping"
                                         "".format(grid_var_name))
            grid_var = ds.variables[grid_var_name]

            grid_mapping_name = getattr(grid_var, 'grid_mapping_name', None)

            # Grid mapping name must be in appendix F
            valid_grid_mapping.assert_true(grid_mapping_name in grid_mapping_dict,
                                           "{} is not a valid grid_mapping_name. See Appendix F for valid grid mappings"
                                           "".format(grid_mapping_name))

            # The grid_mapping_dict has a values of:
            # - required attributes
            # - optional attributes (can't check)
            # - required standard_names defined
            # - at least one of these attributes must be defined

            # We can't do any of the other grid mapping checks if it's not a valid grid mapping name
            if grid_mapping_name not in grid_mapping_dict:
                ret_val.append(valid_grid_mapping.to_result())
                continue

            grid_mapping = grid_mapping_dict[grid_mapping_name]
            required_attrs = grid_mapping[0]
            # Make sure all the required attributes are defined
            for req in required_attrs:
                valid_grid_mapping.assert_true(hasattr(grid_var, req),
                                               "{} is a required attribute for grid mapping {}"
                                               "".format(req, grid_mapping_name))

            # Make sure that exactly one of the exclusive attributes exist
            if len(grid_mapping_dict) == 4:
                at_least_attr = grid_mapping_dict[3]
                number_found = 0
                for attr in at_least_attr:
                    if hasattr(grid_var, attr):
                        number_found += 1
                valid_grid_mapping.assert_true(number_found == 1,
                                               "grid mapping {} must define exactly one of these attributes: "
                                               "{}".format(grid_mapping_name, ' or '.join(at_least_attr)))

            # Make sure that exactly one variable is defined for each of the required standard_names
            expected_std_names = grid_mapping[2]
            for expected_std_name in expected_std_names:
                found_vars = ds.get_variables_by_attributes(standard_name=expected_std_name)
                valid_grid_mapping.assert_true(len(found_vars) == 1,
                                               "grid mapping {} requires exactly one variable with standard_name "
                                               "{} to be defined".format(grid_mapping_name, expected_std_name))

            ret_val.append(valid_grid_mapping.to_result())

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

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
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

        for var in ds.get_variables_by_attributes(standard_name='region'):
            valid_region = TestCtx(BaseCheck.MEDIUM,
                                   "§6.1.1 Geographic region specified by {} is valid"
                                   "".format(var.name))
            valid_region.assert_true(''.join(var[:].astype(str)).lower() in region_list,
                                     "{} is not a valid region"
                                     "".format(''.join(var[:].astype(str))))
            ret_val.append(valid_region.to_result())
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
        :return: List of results
        """

        # Note that test does not check monotonicity
        ret_val = []
        reasoning = []

        for variable_name, boundary_variable_name in cfutil.get_cell_boundary_map(ds).items():
            variable = ds.variables[variable_name]
            valid = True
            reasoning = []
            if boundary_variable_name not in ds.variables:
                valid = False
                reasoning.append("Boundary variable {} referenced by {} not "
                                 "found in dataset variables".format(boundary_variable.name,
                                                                     variable.name))
            else:
                boundary_variable = ds.variables[boundary_variable_name]
            # The number of dimensions in the bounds variable should always be
            # the number of dimensions in the referring variable + 1
            if (boundary_variable.ndim < 2):
                valid = False
                reasoning.append('Boundary variable {} should have at least two'
                                 'dimensions to enclose the base case of a one dimensionsal variable'.format(boundary_variable.name))
            if (boundary_variable.ndim != variable.ndim + 1):
                valid = False
                reasoning.append('The number of dimensions of the variable %s is %s, but the '
                                 'number of dimensions of the boundary variable %s is %s. The boundary variable '
                                 'should have %s dimensions' %
                                 (variable.name, variable.ndim,
                                  boundary_variable.name,
                                  boundary_variable.ndim,
                                  variable.ndim + 1))
            if (variable.dimensions[:] !=
                  boundary_variable.dimensions[:variable.ndim]):
                valid = False
                reasoning.append(u"Boundary variable coordinates are in improper order: {}. Bounds-specific dimensions should be last".format(
                                boundary_variable.dimensions))

            # ensure p vertices form a valid simplex given previous a...n
            # previous auxiliary coordinates
            if (ds.dimensions[boundary_variable.dimensions[-1]].size <
                len(boundary_variable.dimensions[:-1]) + 1):
                valid = False
                reasoning.append("Boundary variable dimension {} must have at least {} elements to form a simplex/closed cell with previous dimensions {}.".format(boundary_variable.name,
                                                                                                                                                          len(variable.dimensions) + 1,
                                                                                                                                                          boundary_variable.dimensions[:-1]))
            result = Result(BaseCheck.MEDIUM, valid,
                            "§7.1 Cell boundaries are valid for variable {}".format(variable_name),
                            reasoning)
            ret_val.append(result)

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

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        reasoning = []
        var_names = ds.get_variables_by_attributes(cell_measures=lambda c:
                                                   c is not None)
        for var_name in var_names:
            var = ds.variables[var_name]
            search_str = '^(?:area|volume): (\w+)$'
            search_res = re.search(search_str, var.cell_measures)
            if not search_res:
                valid = False
                reasoning.append("The cell_measures attribute for variable {} "
                                 "is formatted incorrectly.  It should take the"
                                 " form of either 'area: cell_var' or "
                                 "'volume: cell_var' where cell_var is the "
                                 "variable describing the cell measures".format(
                                     var_name))
            else:
                valid = True
                cell_meas_var_name = search_res.groups[0]
                # TODO: cache previous results
                if not cell_meas_var_name in ds.variables:
                    valid = False
                    reasoning.append("Cell measure variable {} referred to by "
                                     "{} is not present in dataset variables".format(
                                                var_name, cell_meas_var_name))
                else:
                    cell_meas_var = ds.variables[cell_meas_var_name]
                    if not hasattr(cell_meas_var, 'units'):
                        valid = False
                        reasoning.append("Cell measure variable {} is required "
                                         "to have units attribute defined.".format(
                                                        cell_meas_var_name))
                    if not set(cell_meas_var.dimensions).issubset(
                                               var.dimensions):
                        valid = False
                        reasoning.append("Cell measure variable {} must have "
                                         "dimensions which are a subset of "
                                         "those defined in variable {}.".format(
                                                  cell_meas_var_name, var_name))

            result = Result(BaseCheck.MEDIUM,
                            valid,
                            ('§7.2 Cell measures', var_name, 'cell_measures'),
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

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """

        methods = [
            "point",
            "sum",
            "mean",
            "maximum",
            "minimum",
            "mid_range",
            "standard_deviation",
            "variance",
            "mode",
            "median"
        ]

        ret_val = []
        # The basic format is `name: method (
        psep = re.compile(r'((?P<var>\w+): (?P<method>\w+) ?(?P<where>where (?P<wtypevar>\w+) '
                          '?(?P<over>over (?P<otypevar>\w+))?| ?)(?P<brace>\(((?P<brace_wunit>\w+): '
                          '(\d+) (?P<unit>\w+)|(?P<brace_opt>\w+): (\w+))\))*)')

        for var in ds.get_variables_by_attributes(cell_methods=lambda x: x is not None):
            if not getattr(var, 'cell_methods', ''):
                continue

            method = getattr(var, 'cell_methods', '')

            valid_attribute = TestCtx(BaseCheck.HIGH,
                                      '§7.1 {} has a valid cell_methods attribute format'.format(var.name))
            valid_attribute.assert_true(re.match(psep, method) is not None,
                                        '"{}" is not a valid format for cell_methods attribute'
                                        ''.format(method))
            ret_val.append(valid_attribute.to_result())

            valid_cell_names = TestCtx(BaseCheck.MEDIUM,
                                       '§7.3 {} has valid names in cell_methods attribute'.format(var.name))

            # check that the name is valid
            for match in re.finditer(psep, method):
                valid = False
                if match.group('var') in var.dimensions:
                    valid = True
                elif match.group('var') == 'area':
                    valid = True
                elif match.group('var') in getattr(var, "coordinates", ""):
                    valid = True

                valid_cell_names.assert_true(valid,
                                             'cell_methods name component {} does not match a dimension, area or auxiliary coordinate'
                                             ''.format(match.group('var')))

            ret_val.append(valid_cell_names.to_result())

            # Checks if the method value of the 'name: method' pair is acceptable
            valid_cell_methods = TestCtx(BaseCheck.MEDIUM,
                                         '§7.3 {} has valid methods in cell_methods attribute'.format(var.name))

            for match in re.finditer(psep, method):
                valid_cell_methods.assert_true(match.group('method') in methods,
                                               '{}:cell_methods contains an invalid method: {}'
                                               ''.format(var.name, match.group('method')))

            ret_val.append(valid_cell_methods.to_result())

            valid_modifier = TestCtx(BaseCheck.MEDIUM,
                                     '§7.3.3 {} has valid cell_methods modifiers'.format(var.name))

            for match in re.finditer(psep, method):
                if match.group('brace') is not None:
                    valid_modifier.assert_true(match.group('brace_wunit') in ('interval', 'comment', 'area'),
                                               '{}:cell_methods contains an invalid modifier: {}. It should be one '
                                               'of interval, comment or area.'
                                               ''.format(var.name, match.group('brace_wunit')))

            if valid_modifier.out_of > 0:
                ret_val.append(valid_modifier.to_result())

        return ret_val

    def check_climatological_statistics(self, ds):
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

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        reasoning = []
        ret_val = []
        total_climate_count = 0
        valid_climate_count = 0

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

        # first, to determine whether or not we have a climatological time
        # variable, we need to make sure it has the attribute "climatology",
        # but not the attribute "bounds"
        meth_regex = "(?:{})".format("|".join(methods))
        clim_containing_vars = ds.get_variables_by_attributes(
                                        climatology=lambda s: s is not None)
        clim_var = clim_containing_vars[0] if clim_containing_vars else None
        if clim_var:
            if hasattr(clim_var, 'bounds'):
                reasoning.append('Variable {} has a climatology attribute and cannot also have a bounds attribute.'.format(clim_var.name))
                result = Result(BaseCheck.MEDIUM,
                                False,
                                ('§7.3 Cell Methods', clim_var, 'cell_methods_climatology'),
                                reasoning)
                ret_val.append(result)
                return ret_val
            # make sure the climatology variable referenced actually exists
            elif clim_var.climatology not in ds.variables:
                reasoning.append("Variable {} referenced in time's climatology attribute does not exist".format(ds.variables['time'].climatology))
                result = Result(BaseCheck.MEDIUM,
                                False,
                                ('§7.3 Cell Methods', clim_var, 'cell_methods_climatology'),
                                reasoning)
                ret_val.append(result)
                return ret_val
            # handle 1-d and 2d coordinate bounds
            if (clim_var.ndim + 1 != ds.variables[clim_var.climatology].ndim):
                valid = False
                # Probably realistically need two dimensions in majority of
                # practical cases.
                reasoning.append('The number of dimensions of the climatology variable %s is %s, but the '
                                 'number of dimensions of the referencing variable %s is %s. The climatology variable '
                                 'should have %s dimensions' %
                                 (ds.variables[clim_var.climatology].name,
                                  ds.variables[clim_var.climatology].ndim,
                                  clim_var.name,
                                  clim_var.ndim,
                                  clim_var.ndim + 1))
                return ret_val
            # check that coordinate bounds are in the proper order.
            # make sure last elements are boundary variable specific dimensions
            elif (clim_var.dimensions[:] !=
                  ds.variables[clim_var.climatology].dimensions[:clim_var.ndim]):
                valid = False
                reasoning.append(u"Climatology variable coordinates are in improper order: {}. Bounds-specific dimensions should be last".format(
                                ds.variables[clim_var.climatology].dimensions))
                return ret_val
            elif ds.dimensions[ds.variables[clim_var.climatology].dimensions[-1]].size != 2:
                valid = False
                reasoning.append(u"Climatology dimension {} should only contain two elements".format(
                                boundary_variable.dimensions))
        # catchall
        return ret_val


        # otherwise match the following values with for variable with
        # `cell_methods` attributes
        # time: method1 within years time: method2 over years
        # time: method1 within days time: method2 over days
        # time: method1 within days time: method2 over days time: method3 over years
        # optionally followed by parentheses for explaining additional
        # info, e.g.
        # "time: method1 within years time: method2 over years (sidereal years)"

        meth_regex = "(?:{})".format("|".join(methods))
        re_string = (r"^time: {0} within (?:years time: {0} over years|"
                     r"days time: {0} over days"
                     r"(?: time: {0} over years)?)(?: \([^)]+\))?$".format(meth_regex))
        # find any variables with a valid climatological cell_methods
        for cell_method_var in ds.get_variables_by_attributes(cell_methods=lambda s: s is not
                                                              None):
            total_climate_count += 1
            if not re.search(re_string, cell_method_var.cell_methods):
                reasoning.append('The "time: method within years/days over years/days" format is not correct in variable {}.'.format(cell_method_var.name))
            else:
                valid_climate_count += 1

            result = Result(BaseCheck.MEDIUM,
                            (valid_climate_count, total_climate_count),
                            ('§7.4 Climatological Statistics', clim_var, 'cell_methods_climatology'),
                            reasoning)
            ret_val.append(result)

        return ret_val


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

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
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
                if not isinstance(scale_factor, (float, np.floating)):
                    valid = False
                    reasoning.append("Attributes add_offset and scale_factor are not of type float or double.")
                else:
                    # Check variable type is byte, short or int
                    if var.dtype not in [np.int, np.int8, np.int16, np.int32, np.int64]:
                        valid = False
                        reasoning.append("Variable is not of type byte, short, or int.")

            result = Result(BaseCheck.MEDIUM,
                            valid,
                            '§8.1 Packed Data defined by {} contains valid packing'
                            ''.format(name),
                            reasoning)
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

            result = Result(BaseCheck.MEDIUM,
                            valid,
                            '§8.1 Packed Data defined by {} contains valid data types'.format(name),
                            reasoning)
            ret_val.append(result)

        return ret_val


    def check_compression_gathering(self, ds):
        """
        At the current time the netCDF interface does not provide for packing
        data. However a simple packing may be achieved through the use of the
        optional NUG defined attributes scale_factor and add_offset . After the
        data values of a variable have been read, they are to be multiplied by
        the scale_factor , and have add_offset added to them. If both
        attributes are present, the data are scaled before the offset is added.
        When scaled data are written, the application should first subtract the
        offset and then divide by the scale factor. The units of a variable
        should be representative of the unpacked data.

        This standard is more restrictive than the NUG with respect to the use
        of the scale_factor and add_offset attributes; ambiguities and
        precision problems related to data type conversions are resolved by
        these restrictions. If the scale_factor and add_offset attributes are
        of the same data type as the associated variable, the unpacked data is
        assumed to be of the same data type as the packed data. However, if the
        scale_factor and add_offset attributes are of a different data type
        from the variable (containing the packed data) then the unpacked data
        should match the type of these attributes, which must both be of type
        float or both be of type double . An additional restriction in this
        case is that the variable containing the packed data must be of type
        byte , short or int . It is not advised to unpack an int into a float
        as there is a potential precision loss.

        When data to be packed contains missing values the attributes that
        indicate missing values ( _FillValue , valid_min , valid_max ,
                                 valid_range ) must be of the same data type as
        the packed data. See Section 2.5.1, “Missing Data” for a discussion of
        how applications should treat variables that have attributes indicating
        both missing values and transformations defined by a scale and/or
        offset.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        for compress_var in ds.get_variables_by_attributes(compress=lambda s: s is not None):
            valid = True
            reasoning = []
            # puts the referenced variable being compressed into a set
            compress_set = set(compress_var.compress.split(' '))
            if compress_var.ndim != 1:
                valid = False
                reasoning.append("Compression variable {} may only have one dimension".format(compress_var.name))
            # ensure compression variable is a proper index, and thus is an
            # signed or unsigned integer type of some sort
            if compress_var.dtype.kind not in {'i', 'u'}:
                valid = False
                reasoning.append("Compression variable {} must be an integer type to form a proper array index".format(compress_var.name))
            # make sure all the variables referred to are contained by the
            # variables.
            if not compress_set.issubset(ds.dimensions):
                not_in_dims = sorted(compress_set.difference(ds.dimensions))
                valid = False
                reasoning.append("The following dimensions referenced by the compress attribute of variable {} do not exist: {}".format(compress_var.name, not_in_dims))

            result = Result(BaseCheck.MEDIUM,
                            valid,
                            '§8.2 Compression by gathering for variable {}'.format(compress_var.name),
                            reasoning)
            ret_val.append(result)

        return ret_val

    ###############################################################################
    #
    # Chapter 9: Discrete Sampling Geometries
    #
    ###############################################################################

    def check_all_features_are_same_type(self, ds):
        """
        Check that the feature types in a dataset are all the same.

        9.1 The features contained within a collection must always be of the same type; and all the collections in a CF file
        must be of the same feature type.

        point, timeSeries, trajectory, profile, timeSeriesProfile, trajectoryProfile.

        The space-time coordinates that are indicated for each feature are mandatory.  However a featureType may also include
        other space-time coordinates which are not mandatory (notably the z coordinate).

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """
        all_the_same = TestCtx(BaseCheck.HIGH,
                               '§9.1 Feature Types are all the same')
        feature_types_found = defaultdict(list)
        for name in self._find_geophysical_vars(ds):
            feature = cfutil.guess_feature_type(ds, name)
            # If we can't figure out the feature type, don't penalize, just
            # make a note of it in the messages
            if feature is not None:
                feature_types_found[feature].append(name)
            else:
                all_the_same.messages.append("Unidentifiable feature for variable {}"
                                             "".format(name))
        feature_description = ', '.join(['{} ({})'.format(ftr, ', '.join(vrs)) for ftr, vrs in feature_types_found.items()])

        all_the_same.assert_true(len(feature_types_found) < 2,
                                 "Different feature types discovered in this dataset: {}"
                                 "".format(feature_description))

        return all_the_same.to_result()

    def check_feature_type(self, ds):
        """
        Check the global attribute featureType for valid CF featureTypes

        9.4 A global attribute, featureType, is required for all Discrete Geometry representations except the orthogonal
        multidimensional array representation, for which it is highly recommended.

        The value assigned to the featureType attribute is case-insensitive.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """
        feature_list = ['point', 'timeSeries', 'trajectory', 'profile', 'timeSeriesProfile', 'trajectoryProfile']

        feature_type = getattr(ds, 'featureType', None)

        valid_feature_type = TestCtx(BaseCheck.HIGH, '§9.1 Dataset contains a valid featureType')
        valid_feature_type.assert_true(feature_type is None or feature_type in feature_list,
                                       "{} is not a valid CF featureType. It must be one of {}"
                                       "".format(feature_type, ', '.join(feature_list)))

        return valid_feature_type.to_result()

    def check_cf_role(self, ds):
        """
        Check variables defining cf_role for legal cf_role values.

        §9.5 The only acceptable values of cf_role for Discrete Geometry CF
        data sets are timeseries_id, profile_id, and trajectory_id

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """
        ret_val = []
        valid_roles = ['timeseries_id', 'profile_id', 'trajectory_id']
        for variable in ds.get_variables_by_attributes(cf_role=lambda x: x is not None):
            name = variable.name
            valid_cf_role = TestCtx(BaseCheck.HIGH, '§9.5 {} contains a valid cf_role attribute'.format(name))
            cf_role = variable.cf_role
            valid_cf_role.assert_true(cf_role in valid_roles,
                                      "{} is not a valid cf_role value. It must be one of {}"
                                      "".format(name, ', '.join(valid_roles)))
        return ret_val

    def check_variable_features(self, ds):
        '''
        Checks the variable feature types match the dataset featureType attribute

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []
        feature_list = ['point', 'timeSeries', 'trajectory', 'profile', 'timeSeriesProfile', 'trajectoryProfile']
        # Don't bother checking if it's not a legal featureType
        feature_type = getattr(ds, 'featureType', None)
        if feature_type not in feature_list:
            return []

        feature_type_map = {
            'point': [
                'point'
            ],
            'timeSeries': [
                'timeseries',
                'multi-timeseries-orthogonal',
                'multi-timeseries-incomplete',
            ],
            'trajectory': [
                'cf-trajectory',
                'single-trajectory',
            ],
            'profile': [
                'profile-orthogonal',
                'profile-incomplete'
            ],
            'timeSeriesProfile': [
                'timeseries-profile-single-station',
                'timeseries-profile-multi-station',
                'timeseries-profile-single-ortho-time',
                'timeseries-profile-multi-ortho-time',
                'timeseries-profile-ortho-depth',
                'timeseries-profile-incomplete'
            ],
            'trajectoryProfile': [
                'trajectory-profile-orthogonal',
                'trajectory-profile-incomplete'
            ]
        }
        for name in self._find_geophysical_vars(ds):
            variable_feature = cfutil.guess_feature_type(ds, name)
            # If we can't figure it out, don't check it.
            if variable_feature is None:
                continue
            matching_feature = TestCtx(BaseCheck.MEDIUM,
                                       '§9.1 Feature Type for {} is valid {}'
                                       ''.format(name, feature_type))
            matching_feature.assert_true(variable_feature in feature_type_map[feature_type],
                                         '{} is not a {}, it is detected as a {}'
                                         ''.format(name, feature_type, variable_feature))
            ret_val.append(matching_feature.to_result())

        return ret_val

    def check_hints(self, ds):
        '''
        Checks for potentially mislabeled metadata and makes suggestions for how to correct

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []

        ret_val.extend(self._check_hint_bounds(ds))

        return ret_val

    def _check_hint_bounds(self, ds):
        '''
        Checks for variables ending with _bounds, if they are not cell methods,
        make the recommendation

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        '''
        ret_val = []
        boundary_variables = cfutil.get_cell_boundary_variables(ds)
        for name in ds.variables:
            if name.endswith('_bounds') and name not in boundary_variables:
                msg = ('{} might be a cell boundary variable but there are no variables that define it '
                       'as a boundary using the `bounds` attribute.'.format(name))
                result = Result(BaseCheck.LOW,
                                True,
                                '§7.1 {} is a potential cell boundary variable'.format(name),
                                [msg])
                ret_val.append(result)

        return ret_val


class CFNCCheck(BaseNCCheck, CFBaseCheck):

    @classmethod
    def beliefs(cls):  # @TODO
        return {}
