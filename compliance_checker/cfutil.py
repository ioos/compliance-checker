#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
compliance_checker/cfutil.py
'''
from cf_units import Unit
from pkg_resources import resource_filename
from collections import defaultdict
import csv
import json
import re

# For python2/python3 support
try:
    basestring
except NameError:
    basestring = str


_UNITLESS_DB = None
_SEA_NAMES = None

VALID_LAT_UNITS = [
    'degrees_north',
    'degree_north',
    'degree_n',
    'degrees_n',
    'degreen',
    'degreesn'
]
VALID_LON_UNITS = [
    'degrees_east',
    'degree_east',
    'degree_e',
    'degrees_e',
    'degreee',
    'degreese'
]


# We can't import appendix d without getting circular imports
DIMENSIONLESS_VERTICAL_COORDINATES = [
    'ocean_s_coordinate',
    'ocean_s_coordinate_g1',
    'ocean_s_coordinate_g2',
    'atmosphere_hybrid_sigma_pressure_coordinate',
    'atmosphere_hybrid_height_coordinate',
    'ocean_double_sigma_coordinate',
    'ocean_sigma_z_coordinate',
    'ocean_sigma_coordinate',
    'atmosphere_sigma_coordinate',
    'atmosphere_ln_pressure_coordinate',
    'atmosphere_sleve_coordinate'
]


def get_unitless_standard_names(xml_tree, units):
    '''
    Returns True if the units are unitless. Unitless includes units that have
    no units and units that are defined as '1'.
    '''
    found_standard_name = xml_tree.find(".//entry[@id='{}']".format(units))
    if found_standard_name is not None:
        canonical_units = found_standard_name.find('canonical_units')
        return canonical_units is None or canonical_units.text == '1'
    # if the standard name is not found, assume we need units for the time being
    else:
        return False


def get_sea_names():
    '''
    Returns a list of NODC sea names

    source of list: https://www.nodc.noaa.gov/General/NODC-Archive/seanamelist.txt
    '''
    global _SEA_NAMES
    if _SEA_NAMES is None:
        buf = {}
        with open(resource_filename('compliance_checker', 'data/seanames.csv'), 'r') as f:
            reader = csv.reader(f)
            for code, sea_name in reader:
                buf[sea_name] = code
        _SEA_NAMES = buf
    return _SEA_NAMES


def is_unitless(ds, variable):
    '''
    Returns true if the variable is unitless

    Note units of '1' are considered whole numbers or parts but still represent
    physical units and not the absence of units.

    :param netCDF4.Dataset ds: An open netCDF dataset
    :param str variable: Name of the variable
    '''
    units = getattr(ds.variables[variable], 'units', None)
    if units is None or units == '':
        return True
    return False


def is_geophysical(ds, variable):
    '''
    Returns true if the dataset's variable is likely a geophysical variable

    :param netCDF4.Dataset ds: An open netCDF dataset
    :param str variable: Name of the variable
    '''
    ncvar = ds.variables[variable]

    if getattr(ncvar, 'cf_role', None):
        return False

    # Check for axis
    if getattr(ncvar, 'axis', None):
        return False

    standard_name = getattr(ncvar, 'standard_name', '')
    unitless = is_unitless(ds, variable)

    # Is the standard name associated with coordinates
    if standard_name in ('time', 'latitude', 'longitude', 'height', 'depth', 'altitude'):
        return False

    if variable in get_coordinate_variables(ds):
        return False

    if variable in get_auxiliary_coordinate_variables(ds):
        return False

    # Is it dimensionless and unitless?
    if len(ncvar.shape) == 0 and unitless:
        return False

    # Is it a QC Flag?
    if 'status_flag' in standard_name or hasattr(ncvar, 'flag_meanings'):
        return False

    # Is it a §7.1 Cell Boundaries variable
    if variable in get_cell_boundary_variables(ds):
        return False

    if variable == get_climatology_variable(ds):
        return False

    # Is it a string but with no defined units?
    if ncvar.dtype.char == 'S':
        return False

    # Is it an instrument descriptor?
    if variable in get_instrument_variables(ds):
        return False

    # What about a platform descriptor?
    if variable in get_platform_variables(ds):
        return False

    # Skip count variables too
    if hasattr(ncvar, 'sample_dimension'):
        return False

    return True


def get_coordinate_variables(ds):
    '''
    Returns a list of variable names that identify as coordinate variables.

    A coordinate variable is a netCDF variable with exactly one dimension. The
    name of this dimension must be equivalent to the variable name.

    From CF §1.2 Terminology

    It is a one-dimensional variable with the same name as its dimension [e.g.,
    time(time) ], and it is defined as a numeric data type with values that are
    ordered monotonically. Missing values are not allowed in coordinate
    variables.

    :param netCDF4.Dataset ds: An open netCDF dataset
    '''
    coord_vars = []
    for dimension in ds.dimensions:
        if dimension in ds.variables:
            if ds.variables[dimension].dimensions == (dimension,):
                coord_vars.append(dimension)
    return coord_vars


def get_auxiliary_coordinate_variables(ds):
    '''
    Returns a list of auxiliary coordinate variables

    An auxiliary coordinate variable is any netCDF variable that contains
    coordinate data, but is not a coordinate variable (in the sense of the term
    defined by CF).

    :param netCDf4.Dataset ds: An open netCDF dataset
    '''
    aux_vars = []
    # get any variables referecned by the coordinates attribute
    for ncvar in ds.get_variables_by_attributes(coordinates=lambda x: isinstance(x, basestring)):
        # split the coordinates into individual variable names
        referenced_variables = ncvar.coordinates.split(' ')
        # if the variable names exist, add them
        for referenced_variable in referenced_variables:
            if referenced_variable in ds.variables and referenced_variable not in aux_vars:
                aux_vars.append(referenced_variable)

    # axis variables are automatically in
    for variable in get_axis_variables(ds):
        if variable not in aux_vars:
            aux_vars.append(variable)

    # Last are any variables that define the common coordinate standard names
    coordinate_standard_names = ['time', 'longitude', 'latitude', 'height', 'depth', 'altitude']
    coordinate_standard_names += DIMENSIONLESS_VERTICAL_COORDINATES

    # Some datasets like ROMS use multiple variables to define coordinates
    for ncvar in ds.get_variables_by_attributes(standard_name=lambda x: x in coordinate_standard_names):
        if ncvar.name not in aux_vars:
            aux_vars.append(ncvar.name)

    # Remove any that are purely coordinate variables
    ret_val = []
    for aux_var in aux_vars:
        if ds.variables[aux_var].dimensions == (aux_var,):
            continue
        ret_val.append(aux_var)

    return ret_val


def get_cell_boundary_map(ds):
    '''
    Returns a dictionary mapping a variable to it's boundary variable. The
    returned dictionary maps a string variable name to the name of the boundary
    variable.

    :param netCDF4.Dataset nc: netCDF dataset
    '''
    boundary_map = {}
    for variable in ds.get_variables_by_attributes(bounds=lambda x: x is not None):
        if variable.bounds in ds.variables:
            boundary_map[variable.name] = variable.bounds
    return boundary_map


def get_cell_boundary_variables(ds):
    '''
    Returns a list of variable names for variables that represent cell
    boundaries through the `bounds` attribute

    :param netCDF4.Dataset nc: netCDF dataset
    '''
    boundary_variables = []
    has_bounds = ds.get_variables_by_attributes(bounds=lambda x: x is not None)
    for var in has_bounds:
        if var.bounds in ds.variables:
            boundary_variables.append(var.bounds)
    return boundary_variables


def get_geophysical_variables(ds):
    '''
    Returns a list of variable names for the variables detected as geophysical
    variables.

    :param netCDF4.Dataset nc: An open netCDF dataset
    '''

    parameters = []
    for variable in ds.variables:
        if is_geophysical(ds, variable):
            parameters.append(variable)
    return parameters


def get_z_variable(nc):
    '''
    Returns the name of the variable that defines the Z axis or height/depth

    :param netCDF4.Dataset nc: netCDF dataset
    '''
    z_variables = get_z_variables(nc)
    if not z_variables:
        return None

    # Priority is standard_name, units
    for var in z_variables:
        ncvar = nc.variables[var]
        if getattr(ncvar, 'standard_name', None) in ('depth', 'height', 'altitude'):
            return var

    for var in z_variables:
        ncvar = nc.variables[var]
        units = getattr(ncvar, 'units', None)
        if isinstance(units, basestring):
            if units_convertible(units, 'bar'):
                return var
            if units_convertible(units, 'm'):
                return var

    return z_variables[0]


def get_z_variables(nc):
    '''
    Returns a list of all variables matching definitions for Z

    :param netcdf4.dataset nc: an open netcdf dataset object
    '''
    z_variables = []
    # Vertical coordinates will be identifiable by units of pressure or the
    # presence of the positive attribute with a value of up/down

    # optionally, the vertical type may be indicated by providing the
    # standard_name attribute or axis='Z'

    total_coords = get_coordinate_variables(nc) + get_auxiliary_coordinate_variables(nc)
    for coord_name in total_coords:
        if coord_name in z_variables:
            continue
        coord_var = nc.variables[coord_name]
        units = getattr(coord_var, 'units', None)
        positive = getattr(coord_var, 'positive', None)
        standard_name = getattr(coord_var, 'standard_name', None)
        axis = getattr(coord_var, 'axis', None)
        # If there are no units, we can't identify it as a vertical coordinate
        # by checking pressure or positive
        if units is not None:
            if units_convertible(units, 'bar'):
                z_variables.append(coord_name)
            elif isinstance(positive, basestring):
                if positive.lower() in ['up', 'down']:
                    z_variables.append(coord_name)
        # if axis='Z' we're good
        if coord_name not in z_variables and axis == 'Z':
            z_variables.append(coord_name)
        if coord_name not in z_variables and standard_name in ('depth', 'height', 'altitude'):
            z_variables.append(coord_name)
        if coord_name not in z_variables and standard_name in DIMENSIONLESS_VERTICAL_COORDINATES:
            z_variables.append(coord_name)

    return z_variables


def get_lat_variable(nc):
    '''
    Returns the first variable matching latitude

    :param netcdf4.dataset nc: an open netcdf dataset object
    '''
    latitudes = get_latitude_variables(nc)
    if latitudes:
        return latitudes[0]
    return None


def get_latitude_variables(nc):
    '''
    Returns a list of all variables matching definitions for latitude

    :param netcdf4.dataset nc: an open netcdf dataset object
    '''
    latitude_variables = []
    # standard_name takes precedence
    for variable in nc.get_variables_by_attributes(standard_name="latitude"):
        latitude_variables.append(variable.name)

    # Then axis
    for variable in nc.get_variables_by_attributes(axis='Y'):
        if variable.name not in latitude_variables:
            latitude_variables.append(variable.name)

    for variable in nc.get_variables_by_attributes(units=lambda x: x is not None and x.lower() in VALID_LAT_UNITS):
        if variable.name not in latitude_variables:
            latitude_variables.append(variable.name)

    return latitude_variables


def get_true_latitude_variables(nc):
    '''
    Returns a list of variables defining true latitude.

    CF Chapter 4 refers to latitude as a coordinate variable that can also be
    used in non-standard coordinate systems like rotated pole and other
    projections. Chapter 5 refers to a concept of true latitude where the
    variabe defines latitude in a standard projection.

    True latitude, for lack of a better definition, is simply latitude where
    the standard_name is latitude or the units are degrees_north.

    :param netCDF4.Dataset nc: An open netCDF dataset
    '''
    lats = get_latitude_variables(nc)
    true_lats = []
    for lat in lats:
        standard_name = getattr(nc.variables[lat], "standard_name", None)
        units = getattr(nc.variables[lat], "units", None)
        if standard_name == 'latitude':
            true_lats.append(lat)
        elif isinstance(units, basestring) and units.lower() in VALID_LAT_UNITS:
            true_lats.append(lat)
    return true_lats


def get_lon_variable(nc):
    '''
    Returns the variable for longitude

    :param netCDF4.Dataset nc: netCDF dataset
    '''
    longitudes = get_longitude_variables(nc)
    if longitudes:
        return longitudes[0]
    return None


def get_longitude_variables(nc):
    '''
    Returns a list of all variables matching definitions for longitude

    :param netcdf4.dataset nc: an open netcdf dataset object
    '''
    longitude_variables = []
    # standard_name takes precedence
    for variable in nc.get_variables_by_attributes(standard_name="longitude"):
        longitude_variables.append(variable.name)

    # Then axis
    for variable in nc.get_variables_by_attributes(axis='X'):
        if variable.name not in longitude_variables:
            longitude_variables.append(variable.name)

    for variable in nc.get_variables_by_attributes(units=lambda x: x is not None and x.lower() in VALID_LON_UNITS):
        if variable.name not in longitude_variables:
            longitude_variables.append(variable.name)

    return longitude_variables


def get_true_longitude_variables(nc):
    '''
    Returns a list of variables defining true longitude.

    CF Chapter 4 refers to longitude as a coordinate variable that can also be
    used in non-standard coordinate systems like rotated pole and other
    projections. Chapter 5 refers to a concept of true longitude where the
    variabe defines longitude in a standard projection.

    True longitude, for lack of a better definition, is simply longitude where
    the standard_name is longitude or the units are degrees_north.

    :param netCDF4.Dataset nc: An open netCDF dataset
    '''
    lons = get_longitude_variables(nc)
    true_lons = []
    for lon in lons:
        standard_name = getattr(nc.variables[lon], "standard_name", None)
        units = getattr(nc.variables[lon], "units", None)
        if standard_name == 'longitude':
            true_lons.append(lon)
        elif isinstance(units, basestring) and units.lower() in VALID_LON_UNITS:
            true_lons.append(lon)
    return true_lons


def get_platform_variables(ds):
    '''
    Returns a list of platform variable NAMES

    :param netCDF4.Dataset ds: An open netCDF4 Dataset
    '''
    candidates = []
    for variable in ds.variables:
        platform = getattr(ds.variables[variable], 'platform', '')
        if platform and platform in ds.variables:
            if platform not in candidates:
                candidates.append(platform)

    platform = getattr(ds, 'platform', '')
    if platform and platform in ds.variables:
        if platform not in candidates:
            candidates.append(platform)
    return candidates


def get_instrument_variables(ds):
    '''
    Returns a list of instrument variables

    :param netCDF4.Dataset ds: An open netCDF4 Dataset
    '''
    candidates = []
    for variable in ds.variables:
        instrument = getattr(ds.variables[variable], 'instrument', '')
        if instrument and instrument in ds.variables:
            if instrument not in candidates:
                candidates.append(instrument)

    instrument = getattr(ds, 'instrument', '')
    if instrument and instrument in ds.variables:
        if instrument not in candidates:
            candidates.append(instrument)
    return candidates


def get_time_variable(ds):
    '''
    Returns the likeliest variable to be the time coordiante variable

    :param netCDF4.Dataset ds: An open netCDF4 Dataset
    '''
    for var in ds.variables:
        if getattr(ds.variables[var], 'axis', '') == 'T':
            return var
    else:
        candidates = ds.get_variables_by_attributes(standard_name='time')
        if len(candidates) == 1:
            return candidates[0].name
        else:  # Look for a coordinate variable time
            for candidate in candidates:
                if candidate.dimensions == (candidate.name,):
                    return candidate.name

    # If we still haven't found the candidate
    time_variables = set(get_time_variables(ds))
    coordinate_variables = set(get_coordinate_variables(ds))
    if len(time_variables.intersection(coordinate_variables)) == 1:
        return list(time_variables.intersection(coordinate_variables))[0]

    auxiliary_coordinates = set(get_auxiliary_coordinate_variables(ds))
    if len(time_variables.intersection(auxiliary_coordinates)) == 1:
        return list(time_variables.intersection(auxiliary_coordinates))[0]
    return None


def get_time_variables(ds):
    '''
    Returns a list of variables describing the time coordinate

    :param netCDF4.Dataset ds: An open netCDF4 Dataset
    '''
    time_variables = []
    for variable in ds.get_variables_by_attributes(standard_name='time'):
        time_variables.append(variable.name)

    for variable in ds.get_variables_by_attributes(axis='T'):
        if variable.name not in time_variables:
            time_variables.append(variable.name)

    regx = r'^(?:day|d|hour|hr|h|minute|min|second|s)s? since .*$'
    for variable in ds.get_variables_by_attributes(units=lambda x: isinstance(x, basestring)):
        if re.match(regx, variable.units) and variable.name not in time_variables:
            time_variables.append(variable.name)

    return time_variables


def get_axis_variables(ds):
    '''
    Returns a list of variables that define an axis of the dataset

    :param netCDF4.Dataset ds: An open netCDF4 Dataset
    '''
    axis_variables = []
    for ncvar in ds.get_variables_by_attributes(axis=lambda x: x is not None):
        axis_variables.append(ncvar.name)
    return axis_variables


def get_climatology_variable(ds):
    '''
    Returns the variable describing climatology bounds if it exists.

    Climatology variables are similar to cell boundary variables that describe
    the climatology bounds.

    See Example 7.8 in CF 1.6

    :param netCDF4.Dataset ds: An open netCDF4 Dataset
    :rtype: str or None
    '''
    time = get_time_variable(ds)
    # If there's no time dimension there's no climatology bounds
    if not time:
        return None
    # Climatology variable is simply whatever time points to under the
    # `climatology` attribute.
    if hasattr(ds.variables[time], 'climatology'):
        if ds.variables[time].climatology in ds.variables:
            return ds.variables[time].climatology
    return None


def get_flag_variables(ds):
    '''
    Returns a list of variables that are defined as flag variables

    :param netCDF4.Dataset ds: An open netCDF4 Dataset
    '''
    flag_variables = []
    for name, ncvar in ds.variables.items():
        standard_name = getattr(ncvar, 'standard_name', None)
        if isinstance(standard_name, basestring) and 'status_flag' in standard_name:
            flag_variables.append(name)
        elif hasattr(ncvar, 'flag_meanings'):
            flag_variables.append(name)
    return flag_variables


def get_grid_mapping_variables(ds):
    '''
    Returns a list of grid mapping variables

    :param netCDF4.Dataset ds: An open netCDF4 Dataset
    '''
    grid_mapping_variables = []
    for ncvar in ds.get_variables_by_attributes(grid_mapping=lambda x: x is not None):
        if ncvar.grid_mapping in ds.variables:
            grid_mapping_variables.append(ncvar.grid_mapping)
    return grid_mapping_variables


def get_axis_map(ds, variable):
    '''
    Returns an axis_map dictionary that contains an axis key and the coordinate
    names as values.

    For example::

        {'X': ['longitude'], 'Y': ['latitude'], 'T': ['time']}

    The axis C is for compressed coordinates like a reduced grid, and U is for
    unknown axis. This can sometimes be physical quantities representing a
    continuous discrete axis, like temperature or density.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: Variable name
    '''
    all_coords = get_coordinate_variables(ds) + get_auxiliary_coordinate_variables(ds)

    latitudes = get_latitude_variables(ds)
    longitudes = get_longitude_variables(ds)
    times = get_time_variables(ds)
    heights = get_z_variables(ds)

    coordinates = getattr(ds.variables[variable], "coordinates", None)
    if not isinstance(coordinates, basestring):
        coordinates = ''

    # For example
    # {'x': ['longitude'], 'y': ['latitude'], 't': ['time']}
    axis_map = defaultdict(list)
    for coord_name in all_coords:

        if is_compression_coordinate(ds, coord_name):
            axis = 'C'
        elif coord_name in times:
            axis = 'T'
        elif coord_name in longitudes:
            axis = 'X'
        elif coord_name in latitudes:
            axis = 'Y'
        elif coord_name in heights:
            axis = 'Z'
        else:
            axis = 'U'

        if coord_name in ds.variables[variable].dimensions:
            if coord_name not in axis_map[axis]:
                axis_map[axis].append(coord_name)

        elif coord_name in coordinates:
            if coord_name not in axis_map[axis]:
                axis_map[axis].append(coord_name)

    # Sometimes dimensionless coordinates are considered coordinate variables.
    # I know, it's confusing.
    if 'X' not in axis_map and longitudes:
        axis_map['X'] = longitudes

    if 'Y' not in axis_map and latitudes:
        axis_map['Y'] = latitudes

    if 'T' not in axis_map and times:
        axis_map['T'] = times

    if 'Z' not in axis_map and heights:
        axis_map['Z'] = heights

    return axis_map


def is_coordinate_variable(ds, variable):
    '''
    Returns True if the variable is a coordinate variable

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: Variable name
    '''
    if variable not in ds.variables:
        return False
    return ds.variables[variable].dimensions == (variable,)


def is_compression_coordinate(ds, variable):
    '''
    Returns True if the variable is a coordinate variable that defines a
    compression scheme.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: Variable name
    '''
    # Must be a coordinate variable
    if not is_coordinate_variable(ds, variable):
        return False
    # must have a string attribute compress
    compress = getattr(ds.variables[variable], 'compress', None)
    if not isinstance(compress, basestring):
        return False
    if not compress:
        return False
    # This should never happen or be allowed
    if variable in compress:
        return False
    # Must point to dimensions
    for dim in compress.split():
        if dim not in ds.dimensions:
            return False
    return True


def coordinate_dimension_matrix(nc):
    '''
    Returns a dictionary of coordinates mapped to their dimensions

    :param netCDF4.Dataset nc: An open netCDF dataset
    '''
    retval = {}
    x = get_lon_variable(nc)
    if x:
        retval['x'] = nc.variables[x].dimensions
    y = get_lat_variable(nc)
    if y:
        retval['y'] = nc.variables[y].dimensions

    z = get_z_variable(nc)
    if z:
        retval['z'] = nc.variables[z].dimensions

    t = get_time_variable(nc)
    if t:
        retval['t'] = nc.variables[t].dimensions
    return retval


def is_point(nc, variable):
    '''
    Returns true if the variable is a point feature type

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(o), y(o), z(o), t(o)
    # X(o)

    dims = nc.variables[variable].dimensions

    cmatrix = coordinate_dimension_matrix(nc)
    first_coord = None
    if 't' in cmatrix:
        first_coord = cmatrix['t']
        if len(cmatrix['t']) > 1:
            return False
    if 'x' in cmatrix:
        if first_coord is None:
            first_coord = cmatrix['x']
        if first_coord != cmatrix['x']:
            return False
        if len(cmatrix['x']) > 1:
            return False
    if 'y' in cmatrix:
        if first_coord is None:
            first_coord = cmatrix['y']
        if first_coord != cmatrix['y']:
            return False
        if len(cmatrix['y']) > 1:
            return False
    if 'z' in cmatrix:
        if first_coord is None:
            first_coord = cmatrix['z']
        if first_coord != cmatrix['z']:
            return False
        if len(cmatrix['z']) > 1:
            return False
    if first_coord and dims != first_coord:
        return False
    # Point is indistinguishable from trajectories where the instance dimension
    # is implied (scalar)
    traj_ids = nc.get_variables_by_attributes(cf_role="trajectory_id")
    if traj_ids:
        return False

    return True


def is_timeseries(nc, variable):
    '''
    Returns true if the variable is a time series feature type.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''

    # x, y, z, t(t)
    # X(t)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)
    time_variables = get_time_variables(nc)

    if len(dims) != 1:
        return False
    dim = dims[0]
    if dim not in time_variables:
        return False
    # No other coordinates can vary with time
    if 'x' in cmatrix:
        if len(cmatrix['x']) != 0:
            return False
    if 'y' in cmatrix:
        if len(cmatrix['y']) != 0:
            return False
    if 'z' in cmatrix:
        if len(cmatrix['z']) != 0:
            return False

    return True


def is_multi_timeseries_orthogonal(nc, variable):
    '''
    Returns true if the variable is a orthogonal multidimensional array
    representation of time series. For more information on what this means see
    CF 1.6 §H.2.1

    http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#_orthogonal_multidimensional_array_representation_of_time_series

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(i), y(i), z(i), t(o)
    # X(i, o)
    dims = nc.variables[variable].dimensions

    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 't'):
        if req not in cmatrix:
            return False
    if len(cmatrix['x']) != 1 or cmatrix['x'] != cmatrix['y']:
        return False
    if 'z' in cmatrix and cmatrix['x'] != cmatrix['z']:
        return False

    timevar = get_time_variable(nc)
    if cmatrix['t'] != (timevar,):
        return False

    i = cmatrix['x'][0]
    o = cmatrix['t'][0]
    if dims == (i, o):
        return True
    return False


def is_multi_timeseries_incomplete(nc, variable):
    '''
    Returns true if the variable is an incomplete multidimensional array
    representation of time series. For more information on what this means see
    CF 1.6 §H.2.2

    http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#_incomplete_multidimensional_array_representation_of_time_series

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''

    # x(i), y(i), z(i), t(i, o)
    # X(i, o)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 't'):
        if req not in cmatrix:
            return False
    if len(cmatrix['x']) != 1:
        return False
    if cmatrix['x'] != cmatrix['y']:
        return False
    if len(cmatrix['t']) != 2:
        return False
    if cmatrix['x'][0] != cmatrix['t'][0]:
        return False

    i = cmatrix['x'][0]
    o = cmatrix['t'][1]

    if dims == (i, o):
        return True
    return False


def is_cf_trajectory(nc, variable):
    '''
    Returns true if the variable is a CF trajectory feature type

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(i, o), y(i, o), z(i, o), t(i, o)
    # X(i, o)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 't'):
        if req not in cmatrix:
            return False
    if len(cmatrix['x']) != 2:
        return False
    if cmatrix['x'] != cmatrix['y']:
        return False
    if cmatrix['x'] != cmatrix['t']:
        return False
    if 'z' in cmatrix and cmatrix['x'] != cmatrix['z']:
        return False
    if dims == cmatrix['x']:
        return True
    return False


def is_single_trajectory(nc, variable):
    '''
    Returns true if the variable is a single trajectory feature

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(o), y(o), z(o), t(o)
    # X(o)
    # cf_role must be trajectory
    dims = nc.variables[variable].dimensions

    cmatrix = coordinate_dimension_matrix(nc)
    # Time is required for trajectories
    if 't' not in cmatrix:
        return False
    # Coordinates changing with time are optional
    if 'x' in cmatrix:
        if cmatrix['x'] != cmatrix['t']:
            return False
    if 'y' in cmatrix:
        if cmatrix['y'] != cmatrix['t']:
            return False
    if 'z' in cmatrix:
        if cmatrix['z'] != cmatrix['t']:
            return False
    if dims != cmatrix['t']:
        return False
    traj_ids = nc.get_variables_by_attributes(cf_role="trajectory_id")
    if len(traj_ids) != 1:
        return False
    return True


def is_profile_orthogonal(nc, variable):
    '''
    Returns true if the variable is a orthogonal profile feature type

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # Every profile has the exact same depths, think thermister or ADCP
    # x(i), y(i), z(j), t(i)
    # X(i, j)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 'z', 't'):
        if req not in cmatrix:
            return False
    if len(cmatrix['x']) != 1:
        return False
    if cmatrix['x'] != cmatrix['y']:
        return False
    if cmatrix['x'] != cmatrix['t']:
        return False
    if len(cmatrix['z']) != 1:
        return False

    i = cmatrix['x'][0]
    j = cmatrix['z'][0]

    if dims == (i, j):
        return True
    return False


def is_profile_incomplete(nc, variable):
    '''
    Returns true if the variable is a incomplete profile feature type

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # Every profile may have different depths
    # x(i), y(i), z(i, j), t(i)
    # X(i, j)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 'z', 't'):
        if req not in cmatrix:
            return False
    if len(cmatrix['x']) != 1:
        return False
    if cmatrix['x'] != cmatrix['y']:
        return False
    if cmatrix['x'] != cmatrix['t']:
        return False
    if len(cmatrix['z']) != 2:
        return False
    if cmatrix['z'][0] != cmatrix['x'][0]:
        return False

    i = cmatrix['x'][0]
    j = cmatrix['z'][1]

    if dims == (i, j):
        return True
    return False


def is_timeseries_profile_single_station(nc, variable):
    '''
    Returns true if the variable is a time-series profile that represents a
    single station and each profile is the same length.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''

    # x, y, z(z), t(t)
    # X(t, z)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 'z', 't'):
        if req not in cmatrix:
            return False
    if len(cmatrix['x']) != 0:
        return False
    if cmatrix['x'] != cmatrix['y']:
        return False

    z = get_z_variable(nc)
    if cmatrix['z'] != (z,):
        return False
    t = get_time_variable(nc)
    if cmatrix['t'] != (t,):
        return False

    if dims == (t, z):
        return True
    return False


def is_timeseries_profile_multi_station(nc, variable):
    '''
    Returns true if the variable is a time-series profile that represents multiple stations with orthogonal time and depth

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(i), y(i), z(z), t(t)
    # X(i, t, z)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 'z', 't'):
        if req not in cmatrix:
            return False
    if len(cmatrix['x']) != 1:
        return False
    if cmatrix['x'] != cmatrix['y']:
        return False
    i = cmatrix['x'][0]

    z = get_z_variable(nc)
    if cmatrix['z'] != (z,):
        return False
    t = get_time_variable(nc)
    if cmatrix['t'] != (t,):
        return False

    if dims == (i, t, z):
        return True
    return False


def is_timeseries_profile_single_ortho_time(nc, variable):
    '''
    Returns true if the variable is a time-series profile that represents a
    single station with orthogonal time only.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x, y, z(t, j), t(t)
    # X(t, j)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 'z', 't'):
        if req not in cmatrix:
            return False

    if len(cmatrix['x']) != 0:
        return False
    if cmatrix['x'] != cmatrix['y']:
        return False

    t = get_time_variable(nc)
    if cmatrix['t'] != (t,):
        return False

    if len(cmatrix['z']) != 2:
        return False

    if cmatrix['z'][0] != t:
        return False

    j = cmatrix['z'][1]

    if dims == (t, j):
        return True
    return False


def is_timeseries_profile_multi_ortho_time(nc, variable):
    '''
    Returns true if the variable is a time-series profile that represents a
    multi station with orthogonal time only.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(i), y(i), z(i, t, j), t(t)
    # X(i, t, j)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 'z', 't'):
        if req not in cmatrix:
            return False

    if len(cmatrix['x']) != 1:
        return False
    if cmatrix['x'] != cmatrix['y']:
        return False

    t = get_time_variable(nc)
    if cmatrix['t'] != (t,):
        return False

    if len(cmatrix['z']) != 3:
        return False

    if cmatrix['z'][1] != t:
        return False
    if cmatrix['z'][0] != cmatrix['x'][0]:
        return False

    i = cmatrix['x'][0]
    j = cmatrix['z'][2]

    if dims == (i, t, j):
        return True
    return False


def is_timeseries_profile_ortho_depth(nc, variable):
    '''
    Returns true if the variable is a time-series profile with orthogonal depth
    only.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(i), y(i), z(z), t(i, j)
    # X(i, j, z)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 'z', 't'):
        if req not in cmatrix:
            return False

    if len(cmatrix['x']) != 1:
        return False
    if cmatrix['x'] != cmatrix['y']:
        return False

    z = get_z_variable(nc)
    if cmatrix['z'] != (z,):
        return False

    i = cmatrix['x'][0]

    if len(cmatrix['t']) != 2:
        return False
    if cmatrix['t'][0] != i:
        return False

    j = cmatrix['t'][1]

    if dims == (i, j, z):
        return True
    return False


def is_timeseries_profile_incomplete(nc, variable):
    '''
    Returns true if the variable is a time-series profile incomplete depth and
    incomplete time.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(i), y(i), z(i, j, k), t(i, j)
    # X(i, j, k)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 'z', 't'):
        if req not in cmatrix:
            return False

    if len(cmatrix['x']) != 1:
        return False
    if cmatrix['x'] != cmatrix['y']:
        return False
    i = cmatrix['x'][0]

    if len(cmatrix['t']) != 2:
        return False
    if cmatrix['t'][0] != i:
        return False
    j = cmatrix['t'][1]

    if len(cmatrix['z']) != 3:
        return False
    if cmatrix['z'][0] != i:
        return False
    if cmatrix['z'][1] != j:
        return False
    k = cmatrix['z'][2]

    if dims == (i, j, k):
        return True
    return False


def is_trajectory_profile_orthogonal(nc, variable):
    '''
    Returns true if the variable is a trajectory profile with orthogonal
    depths.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(i, o), y(i, o), z(z), t(i, o)
    # X(i, o, z)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 'z', 't'):
        if req not in cmatrix:
            return False

    if len(cmatrix['x']) != 2:
        return False
    if cmatrix['x'] != cmatrix['y']:
        return False
    if cmatrix['x'] != cmatrix['t']:
        return False

    i, o = cmatrix['x']

    z = get_z_variable(nc)
    if cmatrix['z'] != (z,):
        return False

    if dims == (i, o, z):
        return True
    return False


def is_trajectory_profile_incomplete(nc, variable):
    '''
    Returns true if the variable is a trajectory profile with incomplete
    depths.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(i, o), y(i, o), z(i, o, j), t(i, o)
    # X(i, o, j)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 'z', 't'):
        if req not in cmatrix:
            return False

    if len(cmatrix['x']) != 2:
        return False
    if cmatrix['x'] != cmatrix['y']:
        return False
    if cmatrix['x'] != cmatrix['t']:
        return False

    i, o = cmatrix['x']

    if len(cmatrix['z']) != 3:
        return False

    if cmatrix['z'][0] != i:
        return False
    if cmatrix['z'][1] != o:
        return False

    j = cmatrix['z'][2]

    if dims == (i, o, j):
        return True
    return False


def is_2d_regular_grid(nc, variable):
    '''
    Returns True if the variable is a 2D Regular grid.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(x), y(y), t(t)
    # X(t, y, x)

    if is_mapped_grid(nc, variable):
        return False

    dims = nc.variables[variable].dimensions

    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 't'):
        if req not in cmatrix:
            return False

    x = get_lon_variable(nc)
    y = get_lat_variable(nc)
    t = get_time_variable(nc)

    if cmatrix['x'] != (x,):
        return False
    if cmatrix['y'] != (y,):
        return False
    if cmatrix['t'] != (t,):
        return False

    # Relaxed dimension ordering
    if len(dims) == 3 and x in dims and y in dims and t in dims:
        return True
    return False


def is_2d_static_grid(nc, variable):
    '''
    Returns True if the variable is a 2D Regular grid that does not vary with
    time.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(x), y(y)
    # X(y, x)

    if is_mapped_grid(nc, variable):
        return False

    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y'):
        if req not in cmatrix:
            return False

    x = get_lon_variable(nc)
    y = get_lat_variable(nc)

    if cmatrix['x'] != (x,):
        return False

    if cmatrix['y'] != (y,):
        return False

    if len(dims) != 2 or x not in dims or y not in dims:
        return False

    return True


def is_3d_regular_grid(nc, variable):
    '''
    Returns True if the variable is a 3D Regular grid.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(x), y(y), z(z), t(t)
    # X(t, z, y, x)

    if is_mapped_grid(nc, variable):
        return False

    dims = nc.variables[variable].dimensions

    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 'z', 't'):
        if req not in cmatrix:
            return False

    x = get_lon_variable(nc)
    y = get_lat_variable(nc)
    z = get_z_variable(nc)
    t = get_time_variable(nc)

    if cmatrix['x'] != (x,):
        return False
    if cmatrix['y'] != (y,):
        return False
    if cmatrix['z'] != (z,):
        return False
    if cmatrix['t'] != (t,):
        return False

    # Relaxed dimension ordering
    if len(dims) == 4 and x in dims and y in dims and t in dims and z in dims:
        return True
    return False


def is_3d_static_grid(nc, variable):
    '''
    Returns True if the variable is a 2D Regular grid that does not vary with
    time.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(x), y(y), z(z)
    # X(z, y, x)

    if is_mapped_grid(nc, variable):
        return False

    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 'z'):
        if req not in cmatrix:
            return False

    x = get_lon_variable(nc)
    y = get_lat_variable(nc)
    z = get_z_variable(nc)

    if cmatrix['x'] != (x,):
        return False

    if cmatrix['y'] != (y,):
        return False

    if cmatrix['z'] != (z,):
        return False

    if len(dims) != 3 or x not in dims or y not in dims or z not in dims:
        return False

    return True


def is_mapped_grid(nc, variable):
    '''
    Returns true if the feature-type of variable corresponds to a mapped grid
    type. Characterized by Appedix F of CF-1.6

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(j, i), y(j, i), z?, t?
    # F(t?, z?, j, i)

    # The important and really defining characteristic of mapped grids is that
    # the true latitude and longitude coordinates are functions of (j,i) and
    # that the geophysical variables are also functions of (j,i) in their
    # dimensions.
    dims = nc.variables[variable].dimensions
    # For cases like ROMS, the coordinates are mapped using the coordinates attribute
    variable_coordinates = getattr(nc.variables[variable], 'coordinates', '').split()

    lons = get_longitude_variables(nc)
    for lon in lons:
        if lon in variable_coordinates:
            break
    else:
        lon = get_lon_variable(nc)

    if lon is None:
        return False

    lats = get_latitude_variables(nc)
    for lat in lats:
        if lat in variable_coordinates:
            break
    else:
        lat = get_lat_variable(nc)

    if lat is None:
        return False

    x = nc.variables[lon].dimensions
    y = nc.variables[lat].dimensions

    if len(x) != 2:
        return False
    if x != y:
        return False

    comma_dimension = ','.join(dims)
    # Dimensions must be in the same order and the mapping coordinates i and j
    # must be in the same order
    if ','.join(x) not in comma_dimension:
        return False

    return True


def is_reduced_grid(nc, variable):
    '''
    Returns True if the feature-type of the variable corresponds to a reduced
    horizontal grid.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    axis_map = get_axis_map(nc, variable)

    if 'X' not in axis_map:
        return False
    if 'Y' not in axis_map:
        return False
    if 'C' not in axis_map:
        return False

    compressed_coordinates = axis_map['C']
    if len(compressed_coordinates) > 1:
        return False
    compressed_coordinate = axis_map['C'][0]
    for dim in nc.variables[compressed_coordinate].compress.split():
        if dim not in nc.dimensions:
            return False
    return True


def guess_feature_type(nc, variable):
    '''
    Returns a string describing the feature type for this variable

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    if is_point(nc, variable):
        return 'point'
    if is_timeseries(nc, variable):
        return 'timeseries'
    if is_multi_timeseries_orthogonal(nc, variable):
        return 'multi-timeseries-orthogonal'
    if is_multi_timeseries_incomplete(nc, variable):
        return 'multi-timeseries-incomplete'
    if is_cf_trajectory(nc, variable):
        return 'cf-trajectory'
    if is_single_trajectory(nc, variable):
        return 'single-trajectory'
    if is_profile_orthogonal(nc, variable):
        return 'profile-orthogonal'
    if is_profile_incomplete(nc, variable):
        return 'profile-incomplete'
    if is_timeseries_profile_single_station(nc, variable):
        return 'timeseries-profile-single-station'
    if is_timeseries_profile_multi_station(nc, variable):
        return 'timeseries-profile-multi-station'
    if is_timeseries_profile_single_ortho_time(nc, variable):
        return 'timeseries-profile-single-ortho-time'
    if is_timeseries_profile_multi_ortho_time(nc, variable):
        return 'timeseries-profile-multi-ortho-time'
    if is_timeseries_profile_ortho_depth(nc, variable):
        return 'timeseries-profile-ortho-depth'
    if is_timeseries_profile_incomplete(nc, variable):
        return 'timeseries-profile-incomplete'
    if is_trajectory_profile_orthogonal(nc, variable):
        return 'trajectory-profile-orthogonal'
    if is_trajectory_profile_incomplete(nc, variable):
        return 'trajectory-profile-incomplete'
    if is_2d_regular_grid(nc, variable):
        return '2d-regular-grid'
    if is_2d_static_grid(nc, variable):
        return '2d-static-grid'
    if is_3d_regular_grid(nc, variable):
        return '3d-regular-grid'
    if is_3d_static_grid(nc, variable):
        return '3d-static-grid'
    if is_mapped_grid(nc, variable):
        return 'mapped-grid'
    if is_reduced_grid(nc, variable):
        return 'reduced-grid'


def units_convertible(units1, units2, reftimeistime=True):
    """
    Return True if a Unit representing the string units1 can be converted
    to a Unit representing the string units2, else False.

    :param str units1: A string representing the units
    :param str units2: A string representing the units
    """
    try:
        u1 = Unit(units1)
        u2 = Unit(units2)
    except ValueError:
        return False
    return u1.is_convertible(u2)
