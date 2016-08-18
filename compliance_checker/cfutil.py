#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
compliance_checker/cfutil.py
'''
from pkg_resources import resource_filename
import csv
import json


_UNITLESS_DB = None
_SEA_NAMES = None


def get_unitless_standard_names():
    '''
    Returns a list of valid standard_names that are allowed to be unitless
    '''
    global _UNITLESS_DB
    if _UNITLESS_DB is None:
        with open(resource_filename('compliance_checker', 'data/unitless.json'), 'r') as f:
            _UNITLESS_DB = json.load(f)
    return _UNITLESS_DB


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


def is_geophysical(ds, variable):
    '''
    Returns true if the dataset's variable is likely a geophysical variable
    '''
    ncvar = ds.variables[variable]
    # Does it have a standard name and units?
    standard_name = getattr(ncvar, 'standard_name', '')
    if not standard_name:
        return False
    if standard_name and standard_name not in get_unitless_standard_names():
        units = getattr(ncvar, 'units', '')
        if units == '':
            return False
    if not hasattr(ncvar, 'standard_name') or not hasattr(ncvar, 'units'):
        return False
    if getattr(ncvar, 'standard_name') in ('time', 'latitude', 'longitude', 'height', 'depth', 'altitude'):
        return False
    # Is it dimensionless?
    if len(ncvar.shape) == 0:
        return False
    # Is it a QC Flag?
    if 'status_flag' in ncvar.standard_name:
        return False

    if getattr(ncvar, 'cf_role', None):
        return False

    if getattr(ncvar, 'axis', None):
        return False

    return True


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
    axis_z = nc.get_variables_by_attributes(axis='Z')
    if axis_z:
        return axis_z[0].name
    valid_standard_names = ('depth', 'height', 'altitude')
    z = nc.get_variables_by_attributes(standard_name=lambda x: x in valid_standard_names)
    if z:
        return z[0].name
    return


def get_lat_variable(nc):
    '''
    Returns the variable for latitude

    :param netcdf4.dataset nc: an open netcdf dataset object
    '''
    if 'latitude' in nc.variables:
        return 'latitude'
    latitudes = nc.get_variables_by_attributes(standard_name="latitude")
    if latitudes:
        return latitudes[0].name
    return None


def get_lon_variable(nc):
    '''
    Returns the variable for longitude

    :param netCDF4.Dataset nc: netCDF dataset
    '''
    if 'longitude' in nc.variables:
        return 'longitude'
    longitudes = nc.get_variables_by_attributes(standard_name="longitude")
    if longitudes:
        return longitudes[0].name
    return None


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

    return None


def get_crs_variable(ds):
    '''
    Returns the name of the variable identified by a grid_mapping attribute

    :param netCDF4.Dataset ds: An open netCDF4 Dataset
    '''
    for var in ds.variables:
        grid_mapping = getattr(ds.variables[var], 'grid_mapping', '')
        if grid_mapping and grid_mapping in ds.variables:
            return grid_mapping
    return None


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
    for req in ('x', 'y', 't'):
        if req not in cmatrix:
            return False
    t = get_time_variable(nc)
    if cmatrix['x'] != cmatrix['y'] or cmatrix['x'] != cmatrix['t']:
        return False
    # This is a trajectory
    if cmatrix['t'] == (t,):
        return False
    if len(cmatrix['x']) != 1:
        return False
    if 'z' in cmatrix and cmatrix['x'] != cmatrix['z']:
        return False
    if dims == cmatrix['x']:
        return True
    return False


def is_timeseries(nc, variable):
    '''
    Returns true if the variable is a time series feature type.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''

    # x, y, z, t(o)
    # X(o)
    dims = nc.variables[variable].dimensions

    cmatrix = coordinate_dimension_matrix(nc)
    for req in ('x', 'y', 't'):
        if req not in cmatrix:
            return False
    if len(cmatrix['x']) != 0:
        return False
    if len(cmatrix['y']) != 0:
        return False
    if 'z' in cmatrix and len(cmatrix['z']) != 0:
        return False
    timevar = get_time_variable(nc)

    # time has to be a coordinate variable in this case
    if cmatrix['t'] != (timevar,):
        return False
    if dims == cmatrix['t']:
        return True
    return False


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
    # x(t), y(t), z(t), t(t)
    # X(t)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ('x', 'y', 't'):
        if req not in cmatrix:
            return False
    t = get_time_variable(nc)
    if cmatrix['x'] != (t,):
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


def is_3d_regular_grid(nc, variable):
    '''
    Returns True if the variable is a 3D Regular grid.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    '''
    # x(x), y(y), z(z), t(t)
    # X(t, z, y, x)

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
