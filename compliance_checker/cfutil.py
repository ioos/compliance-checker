#!/usr/bin/env python
"""
compliance_checker/cfutil.py
"""
import csv
import re
import warnings
from collections import defaultdict
from functools import lru_cache, partial
from importlib.resources import files

from cf_units import Unit
from netCDF4 import Dataset

_UNITLESS_DB = None
_SEA_NAMES = None

VALID_LAT_UNITS = {
    "degrees_north",
    "degree_north",
    "degree_n",
    "degrees_n",
    "degreen",
    "degreesn",
}
VALID_LON_UNITS = {
    "degrees_east",
    "degree_east",
    "degree_e",
    "degrees_e",
    "degreee",
    "degreese",
}


# We can't import appendix d without getting circular imports
DIMENSIONLESS_VERTICAL_COORDINATES = {
    "ocean_s_coordinate",
    "ocean_s_coordinate_g1",
    "ocean_s_coordinate_g2",
    "atmosphere_hybrid_sigma_pressure_coordinate",
    "atmosphere_hybrid_height_coordinate",
    "ocean_double_sigma_coordinate",
    "ocean_sigma_z_coordinate",
    "ocean_sigma_coordinate",
    "atmosphere_sigma_coordinate",
    "atmosphere_ln_pressure_coordinate",
    "atmosphere_sleve_coordinate",
}


def attr_membership(attr_val, value_set, attr_type=str, modifier_fn=lambda x: x):
    """
    Helper function passed to netCDF4.Dataset.get_attributes_by_value
    Checks that `attr_val` exists, has the same type as `attr_type`,
    and is contained in `value_set`
    attr_val: The value of the attribute being checked
    attr_type: A type object that the `attr_val` is expected to have the same
               type as.  If the type is not the same, a warning is issued and
               the code attempts to cast `attr_val` to the expected type.
    value_set: The set against which membership for `attr_val` is tested
    modifier_fn: A function to apply to attr_val prior to applying the set
                 membership test

    """
    if attr_val is None:
        return False

    if not isinstance(attr_val, attr_type):
        warnings.warn(
            f"Attribute is of type {type(attr_val)!r}, {attr_type!r} expected. Attempting to cast to expected type.",
            stacklevel=2,
        )
        try:
            # if the expected type is str, try casting to unicode type
            # since str can't be instantiated
            if attr_type is str:
                new_attr_val = str(attr_val)
            else:
                new_attr_val = attr_type(attr_val)
        # catch casting errors
        except (ValueError, UnicodeEncodeError):
            warnings.warn(f"Could not cast to type {attr_type}", stacklevel=2)
            return False
    else:
        new_attr_val = attr_val

    try:
        is_in_set = modifier_fn(new_attr_val) in value_set
    except Exception as e:
        warnings.warn(
            f"Could not apply modifier function {modifier_fn} to value: {e.msg}",
            stacklevel=2,
        )
        return False

    return is_in_set


@lru_cache(128)
def is_dimensionless_standard_name(standard_name_table, standard_name):
    """
    Returns True if the units for the associated standard name are
    dimensionless.  Dimensionless standard names include those that have no
    units and units that are defined as constant units in the CF standard name
    table i.e. '1', or '1e-3'.
    """
    # standard_name must be string, so if it is not, it is *wrong* by default
    if not isinstance(standard_name, str):
        return False
    found_standard_name = standard_name_table.find(
        f".//entry[@id='{standard_name}']",
    )
    if found_standard_name is not None:
        canonical_units = Unit(found_standard_name.find("canonical_units").text)
        return canonical_units.is_dimensionless()
    # if the standard name is not found, assume we need units for the time being
    else:
        return False


def get_sea_names():
    """
    Returns a list of NODC sea names

    source of list: https://www.ncei.noaa.gov/resources/ocean-data-format-codes
    """
    global _SEA_NAMES
    if _SEA_NAMES is None:
        buf = {}
        with open(
            files("compliance_checker") / "data/seanames.csv",
        ) as f:
            reader = csv.reader(f)
            for code, sea_name in reader:
                buf[sea_name] = code
        _SEA_NAMES = buf
    return _SEA_NAMES


def is_unitless(nc, variable):
    """
    Returns true if the variable is unitless

    Note units of '1' are considered whole numbers or parts but still represent
    physical units and not the absence of units.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: Name of the variable
    """
    units = getattr(nc.variables[variable], "units", None)
    return units is None or units == ""


def is_geophysical(nc, variable):
    """
    Returns true if the dataset's variable is likely a geophysical variable

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: Name of the variable
    """
    ncvar = nc.variables[variable]

    if getattr(ncvar, "cf_role", None):
        return False

    # Check for axis
    if getattr(ncvar, "axis", None):
        return False

    standard_name_test = getattr(ncvar, "standard_name", "")
    unitless = is_unitless(nc, variable)

    if not isinstance(standard_name_test, str):
        warnings.warn(
            f"Variable {variable} has non string standard name, Attempting cast to string",
            stacklevel=2,
        )
        try:
            standard_name = str(standard_name_test)
        except ValueError:
            warnings.warn(
                "Unable to cast standard name to string, excluding from geophysical variables",
                stacklevel=2,
            )
    else:
        standard_name = standard_name_test

    # Is the standard name associated with coordinates
    if standard_name in {
        "time",
        "latitude",
        "longitude",
        "height",
        "depth",
        "altitude",
    }:
        return False

    if variable in get_coordinate_variables(nc):
        return False

    if variable in get_auxiliary_coordinate_variables(nc):
        return False

    if variable in get_forecast_metadata_variables(nc):
        return False

    # Is it dimensionless and unitless?
    if len(ncvar.shape) == 0 and unitless:
        return False

    # Is it a QC Flag?
    if "status_flag" in standard_name or hasattr(ncvar, "flag_meanings"):
        return False

    # Is it a §7.1 Cell Boundaries variable
    if variable in get_cell_boundary_variables(nc):
        return False

    if variable == get_climatology_variable(nc):
        return False

    # Is it a string but with no defined units?
    if hasattr(ncvar.dtype, "char") and ncvar.dtype.char == "S":
        return False
    elif ncvar.dtype is str:
        return False

    # Is it an instrument descriptor?
    if variable in get_instrument_variables(nc):
        return False

    # What about a platform descriptor?
    if variable in get_platform_variables(nc):
        return False

    # Skip count/index variables too
    if hasattr(ncvar, "sample_dimension") or hasattr(ncvar, "instance_dimension"):
        return False

    return True


@lru_cache
def get_coordinate_variables(nc):
    """
    Returns a list of variable names that identify as coordinate variables.

    A coordinate variable is a netCDF variable with exactly one dimension. The
    name of this dimension must be equivalent to the variable name.

    From CF §1.2 Terminology

    It is a one-dimensional variable with the same name as its dimension [e.g.,
    time(time) ], and it is defined as a numeric data type with values that are
    ordered monotonically. Missing values are not allowed in coordinate
    variables.

    :param netCDF4.Dataset nc: An open netCDF dataset
    """
    coord_vars = []
    for dimension in nc.dimensions:
        if dimension in nc.variables:
            # TODO: Handle string coordinate variables
            if nc.variables[dimension].dimensions == (dimension,):
                coord_vars.append(dimension)
    return coord_vars


def get_auxiliary_coordinate_variables(nc):
    """
    Returns a list of auxiliary coordinate variables

    An auxiliary coordinate variable is any netCDF variable that contains
    coordinate data, but is not a coordinate variable (in the sense of the term
    defined by CF).

    :param netCDf4.Dataset nc: An open netCDF dataset
    """
    aux_vars = []
    # get any variables referenced by the coordinates attribute
    for ncvar in nc.get_variables_by_attributes(
        coordinates=lambda x: isinstance(x, str),
    ):
        # split the coordinates into individual variable names
        referenced_variables = ncvar.coordinates.split(" ")
        # if the variable names exist, add them
        for referenced_variable in referenced_variables:
            if (
                referenced_variable in nc.variables
                and referenced_variable not in aux_vars
            ):
                aux_vars.append(referenced_variable)

    # axis variables are automatically in
    for variable in get_axis_variables(nc):
        if variable not in aux_vars:
            aux_vars.append(variable)

    # Last are any variables that define the common coordinate standard names
    coordinate_standard_names = [
        "time",
        "longitude",
        "latitude",
        "height",
        "depth",
        "altitude",
    ]
    coordinate_standard_names += DIMENSIONLESS_VERTICAL_COORDINATES

    # Some datasets like ROMS use multiple variables to define coordinates
    for ncvar in nc.get_variables_by_attributes(
        standard_name=lambda x: x in coordinate_standard_names,
    ):
        if ncvar.name not in aux_vars:
            aux_vars.append(ncvar.name)

    # Remove any that are purely coordinate variables
    ret_val = []
    for aux_var in aux_vars:
        if nc.variables[aux_var].dimensions == (aux_var,):
            continue
        ret_val.append(aux_var)

    return ret_val


def get_forecast_metadata_variables(nc):
    """
    Returns a list of variables that represent forecast reference time
    metadata.

    :param netCDF4.Dataset nc: An open netCDF4 Dataset.
    :rtype: list
    """
    forecast_metadata_standard_names = {
        "forecast_period",
        "forecast_reference_time",
    }
    forecast_metadata_variables = []
    for varname in nc.variables:
        standard_name = getattr(nc.variables[varname], "standard_name", None)
        if standard_name in forecast_metadata_standard_names:
            forecast_metadata_variables.append(varname)
    return forecast_metadata_variables


def get_cell_boundary_map(nc):
    """
    Returns a dictionary mapping a variable to its boundary variable. The
    returned dictionary maps a string variable name to the name of the boundary
    variable.

    :param netCDF4.Dataset nc: netCDF dataset
    """
    boundary_map = {}
    for variable in nc.get_variables_by_attributes(bounds=lambda x: x is not None):
        if variable.bounds in nc.variables:
            boundary_map[variable.name] = variable.bounds
    return boundary_map


def get_cell_boundary_variables(nc):
    """
    Returns a list of variable names for variables that represent cell
    boundaries through the `bounds` attribute

    :param netCDF4.Dataset nc: netCDF dataset
    """
    boundary_variables = []
    has_bounds = nc.get_variables_by_attributes(bounds=lambda x: x is not None)
    for var in has_bounds:
        if var.bounds in nc.variables:
            boundary_variables.append(var.bounds)
    return boundary_variables


def get_bounds_variables(nc):
    contains_bounds = nc.get_variables_by_attributes(bounds=lambda s: s in nc.variables)
    return {nc.variables[parent_var.bounds] for parent_var in contains_bounds}


def get_geophysical_variables(nc):
    """
    Returns a list of variable names for the variables detected as geophysical
    variables.

    :param netCDF4.Dataset nc: An open netCDF dataset
    """
    parameters = []
    for variable in nc.variables:
        if is_geophysical(nc, variable) and variable not in get_bounds_variables(nc):
            parameters.append(variable)
    return parameters


def get_z_variable(nc):
    """
    Returns the name of the variable that defines the Z axis or height/depth

    :param netCDF4.Dataset nc: netCDF dataset
    """
    z_variables = get_z_variables(nc)
    if not z_variables:
        return None

    # Priority is standard_name, units
    for var in z_variables:
        ncvar = nc.variables[var]
        if getattr(ncvar, "standard_name", None) in ("depth", "height", "altitude"):
            return var

    for var in z_variables:
        ncvar = nc.variables[var]
        units = getattr(ncvar, "units", None)
        if isinstance(units, str):
            if units_convertible(units, "bar"):
                return var
            if units_convertible(units, "m"):
                return var

    return z_variables[0]


def get_z_variables(nc):
    """
    Returns a list of all variables matching definitions for Z

    :param netcdf4.dataset nc: an open netcdf dataset object
    """
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
        units = getattr(coord_var, "units", None)
        positive = getattr(coord_var, "positive", None)
        standard_name = getattr(coord_var, "standard_name", None)
        axis = getattr(coord_var, "axis", None)
        # If there are no units, we can't identify it as a vertical coordinate
        # by checking pressure or positive
        if units is not None:
            if units_convertible(units, "bar"):
                z_variables.append(coord_name)
            elif isinstance(positive, str):
                if positive.lower() in ["up", "down"]:
                    z_variables.append(coord_name)
        # if axis='Z' we're good
        if coord_name not in z_variables and axis == "Z":
            z_variables.append(coord_name)
        if coord_name not in z_variables and standard_name in (
            "depth",
            "height",
            "altitude",
        ):
            z_variables.append(coord_name)
        if (
            coord_name not in z_variables
            and standard_name in DIMENSIONLESS_VERTICAL_COORDINATES
        ):
            z_variables.append(coord_name)

    return z_variables


def get_lat_variable(nc):
    """
    Returns the first variable matching latitude

    :param netcdf4.dataset nc: an open netcdf dataset object
    """
    latitudes = get_latitude_variables(nc)
    if latitudes:
        return latitudes[0]
    return None


def get_latitude_variables(nc):
    """
    Returns a list of all variables matching definitions for latitude

    :param netcdf4.dataset nc: an open netcdf dataset object
    """
    latitude_variables = []
    # standard_name takes precedence
    for variable in nc.get_variables_by_attributes(standard_name="latitude"):
        latitude_variables.append(variable.name)

    # Then axis
    for variable in nc.get_variables_by_attributes(axis="Y"):
        if not (
            variable.name in latitude_variables
            or getattr(variable, "standard_name", None)
            in {"projection_y_coordinate", "projection_y_angular_coordinate"}
        ):
            latitude_variables.append(variable.name)

    check_fn = partial(
        attr_membership,
        value_set=VALID_LAT_UNITS,
        modifier_fn=lambda s: s.lower(),
    )
    for variable in nc.get_variables_by_attributes(units=check_fn):
        if variable.name not in latitude_variables:
            latitude_variables.append(variable.name)

    return latitude_variables


def get_true_latitude_variables(nc):
    """
    Returns a list of variables defining true latitude.

    CF Chapter 4 refers to latitude as a coordinate variable that can also be
    used in non-standard coordinate systems like rotated pole and other
    projections. Chapter 5 refers to a concept of true latitude where the
    variable defines latitude in a standard projection.

    True latitude, for lack of a better definition, is simply latitude where
    the standard_name is latitude or the units are degrees_north.

    :param netCDF4.Dataset nc: An open netCDF dataset
    """
    lats = get_latitude_variables(nc)
    true_lats = []
    for lat in lats:
        standard_name = getattr(nc.variables[lat], "standard_name", None)
        units = getattr(nc.variables[lat], "units", None)
        if standard_name == "latitude":
            true_lats.append(lat)
        elif isinstance(units, str) and units.lower() in VALID_LAT_UNITS:
            true_lats.append(lat)
    return true_lats


def get_lon_variable(nc):
    """
    Returns the variable for longitude

    :param netCDF4.Dataset nc: netCDF dataset
    """
    longitudes = get_longitude_variables(nc)
    if longitudes:
        return longitudes[0]
    return None


def get_longitude_variables(nc):
    """
    Returns a list of all variables matching definitions for longitude

    :param netcdf4.dataset nc: an open netcdf dataset object
    """
    longitude_variables = []
    # standard_name takes precedence
    for variable in nc.get_variables_by_attributes(standard_name="longitude"):
        longitude_variables.append(variable.name)

    # Then axis
    for variable in nc.get_variables_by_attributes(axis="X"):
        if not (
            variable.name in longitude_variables
            or getattr(variable, "standard_name", None)
            in {"projection_x_coordinate", "projection_x_angular_coordinate"}
        ):
            longitude_variables.append(variable.name)

    check_fn = partial(
        attr_membership,
        value_set=VALID_LON_UNITS,
        modifier_fn=lambda s: s.lower(),
    )
    for variable in nc.get_variables_by_attributes(units=check_fn):
        if variable.name not in longitude_variables:
            longitude_variables.append(variable.name)

    return longitude_variables


def get_true_longitude_variables(nc):
    """
    Returns a list of variables defining true longitude.

    CF Chapter 4 refers to longitude as a coordinate variable that can also be
    used in non-standard coordinate systems like rotated pole and other
    projections. Chapter 5 refers to a concept of true longitude where the
    variable defines longitude in a standard projection.

    True longitude, for lack of a better definition, is simply longitude where
    the standard_name is longitude or the units are degrees_north.

    :param netCDF4.Dataset nc: An open netCDF dataset
    """
    lons = get_longitude_variables(nc)
    true_lons = []
    for lon in lons:
        standard_name = getattr(nc.variables[lon], "standard_name", None)
        units = getattr(nc.variables[lon], "units", None)
        if standard_name == "longitude":
            true_lons.append(lon)
        elif isinstance(units, str) and units.lower() in VALID_LON_UNITS:
            true_lons.append(lon)
    return true_lons


def get_platform_variables(nc):
    """
    Returns a list of platform variable NAMES

    :param netCDF4.Dataset nc: An open netCDF4 Dataset
    """
    candidates = []
    for variable in nc.variables:
        platform = getattr(nc.variables[variable], "platform", "")
        if platform and platform in nc.variables:
            if platform not in candidates:
                candidates.append(platform)

    platform = getattr(nc, "platform", "")
    if platform and platform in nc.variables:
        if platform not in candidates:
            candidates.append(platform)
    return candidates


def get_instrument_variables(nc):
    """
    Returns a list of instrument variables

    :param netCDF4.Dataset nc: An open netCDF4 Dataset
    """
    candidates = []
    for variable in nc.variables:
        instrument = getattr(nc.variables[variable], "instrument", "")
        if instrument and instrument in nc.variables:
            if instrument not in candidates:
                candidates.append(instrument)

    instrument = getattr(nc, "instrument", "")
    if instrument and instrument in nc.variables:
        if instrument not in candidates:
            candidates.append(instrument)
    return candidates


def get_time_variable(nc):
    """
    Returns the likeliest variable to be the time coordinate variable

    :param netCDF4.Dataset nc: An open netCDF4 Dataset
    """
    for var in nc.variables:
        if getattr(nc.variables[var], "axis", "") == "T":
            return var
    else:
        candidates = nc.get_variables_by_attributes(standard_name="time")
        if len(candidates) == 1:
            return candidates[0].name
        else:  # Look for a coordinate variable time
            for candidate in candidates:
                if candidate.dimensions == (candidate.name,):
                    return candidate.name

    # If we still haven't found the candidate
    time_variables = set(get_time_variables(nc))
    coordinate_variables = set(get_coordinate_variables(nc))
    if len(time_variables.intersection(coordinate_variables)) == 1:
        return list(time_variables.intersection(coordinate_variables))[0]

    auxiliary_coordinates = set(get_auxiliary_coordinate_variables(nc))
    if len(time_variables.intersection(auxiliary_coordinates)) == 1:
        return list(time_variables.intersection(auxiliary_coordinates))[0]
    return None


def get_time_variables(nc):
    """
    Returns a list of variables describing the time coordinate

    :param netCDF4.Dataset nc: An open netCDF4 Dataset
    """
    time_variables = set()
    for variable in nc.get_variables_by_attributes(standard_name="time"):
        time_variables.add(variable.name)

    for variable in nc.get_variables_by_attributes(axis="T"):
        if variable.name not in time_variables:
            time_variables.add(variable.name)

    regx = r"^(?:day|d|hour|hr|h|minute|min|second|s)s? since .*$"
    for variable in nc.get_variables_by_attributes(units=lambda x: isinstance(x, str)):
        if re.match(regx, variable.units) and variable.name not in time_variables:
            time_variables.add(variable.name)

    return time_variables


def get_axis_variables(nc):
    """
    Returns a list of variables that define an axis of the dataset

    :param netCDF4.Dataset nc: An open netCDF4 Dataset
    """
    axis_variables = []
    for ncvar in nc.get_variables_by_attributes(axis=lambda x: x is not None):
        axis_variables.append(ncvar.name)
    return axis_variables


def get_climatology_variable(nc):
    """
    Returns the variable describing climatology bounds if it exists.

    Climatology variables are similar to cell boundary variables that describe
    the climatology bounnc.

    See Example 7.8 in CF 1.6

    :param netCDF4.Dataset nc: An open netCDF4 Dataset
    :rtype: str or None
    """
    time = get_time_variable(nc)
    # If there's no time dimension there's no climatology bounds
    if not time:
        return None
    # Climatology variable is simply whatever time points to under the
    # `climatology` attribute.
    if hasattr(nc.variables[time], "climatology"):
        if nc.variables[time].climatology in nc.variables:
            return nc.variables[time].climatology
    return None


def _find_standard_name_modifier_variables(nc, return_deprecated=False):
    def match_modifier_variables(standard_name_str):
        if standard_name_str is None:
            return False
        if not return_deprecated:
            matches = re.search(r"^\w+ +\w+", standard_name_str)
        else:
            matches = re.search(
                r"^\w+ +(?:status_flag|number_of_observations)$",
                standard_name_str,
            )
        return bool(matches)

    return [
        var.name
        for var in nc.get_variables_by_attributes(
            standard_name=match_modifier_variables,
        )
    ]


def get_flag_variables(nc):
    """
    Returns a list of variables that are defined as flag variables

    :param netCDF4.Dataset nc: An open netCDF4 Dataset
    """
    flag_variables = []
    for name, ncvar in nc.variables.items():
        standard_name = getattr(ncvar, "standard_name", None)
        if isinstance(standard_name, str) and "status_flag" in standard_name:
            flag_variables.append(name)
        elif hasattr(ncvar, "flag_meanings"):
            flag_variables.append(name)
    return flag_variables


def get_grid_mapping_variables(nc):
    """
    Returns a list of grid mapping variables

    :param netCDF4.Dataset nc: An open netCDF4 Dataset
    """
    grid_mapping_variables = set()
    for ncvar in nc.get_variables_by_attributes(grid_mapping=lambda x: x is not None):
        if ncvar.grid_mapping in nc.variables:
            grid_mapping_variables.add(ncvar.grid_mapping)
    return grid_mapping_variables


def get_axis_map(nc, variable):
    """
    Returns an axis_map dictionary that contains an axis key and the coordinate
    names as values.

    For example::

        {'X': ['longitude'], 'Y': ['latitude'], 'T': ['time']}

    The axis C is for compressed coordinates like a reduced grid, and U is for
    unknown axis. This can sometimes be physical quantities representing a
    continuous discrete axis, like temperature or density.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: Variable name
    """
    all_coords = get_coordinate_variables(nc) + get_auxiliary_coordinate_variables(nc)

    latitudes = get_latitude_variables(nc)
    longitudes = get_longitude_variables(nc)
    times = get_time_variables(nc)
    heights = get_z_variables(nc)

    coordinates = getattr(nc.variables[variable], "coordinates", None)
    if not isinstance(coordinates, str):
        coordinates = []
    else:
        coordinates = coordinates.split(" ")

    # For example
    # {'x': ['longitude'], 'y': ['latitude'], 't': ['time']}
    axis_map = defaultdict(list)
    for coord_name in all_coords:
        axis = getattr(nc.variables[coord_name], "axis", None)
        if not axis or axis not in ("X", "Y", "Z", "T"):
            if is_compression_coordinate(nc, coord_name):
                axis = "C"
            elif coord_name in times:
                axis = "T"
            elif coord_name in longitudes:
                axis = "X"
            elif coord_name in latitudes:
                axis = "Y"
            elif coord_name in heights:
                axis = "Z"
            else:
                axis = "U"

        if coord_name in nc.variables[variable].dimensions:
            if coord_name not in axis_map[axis]:
                axis_map[axis].append(coord_name)

        elif coord_name in coordinates:
            if coord_name not in axis_map[axis]:
                axis_map[axis].append(coord_name)
    return axis_map


def is_coordinate_variable(nc, variable):
    """
    Returns True if the variable is a coordinate variable

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: Variable name
    """
    if variable not in nc.variables:
        return False
    return nc.variables[variable].dimensions == (variable,)


def is_compression_coordinate(nc, variable):
    """
    Returns True if the variable is a coordinate variable that defines a
    compression scheme.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: Variable name
    """
    # Must be a coordinate variable
    if not is_coordinate_variable(nc, variable):
        return False
    # must have a string attribute compress
    compress = getattr(nc.variables[variable], "compress", None)
    if not isinstance(compress, str):
        return False
    if not compress:
        return False
    # This should never happen or be allowed
    if variable in compress:
        return False
    # Must point to dimensions
    for dim in compress.split():
        if dim not in nc.dimensions:
            return False
    return True


def coordinate_dimension_matrix(nc):
    """
    Returns a dictionary of coordinates mapped to their dimensions

    :param netCDF4.Dataset nc: An open netCDF dataset
    """
    retval = {}
    x = get_lon_variable(nc)
    if x:
        retval["x"] = nc.variables[x].dimensions
    y = get_lat_variable(nc)
    if y:
        retval["y"] = nc.variables[y].dimensions

    z = get_z_variable(nc)
    if z:
        retval["z"] = nc.variables[z].dimensions

    t = get_time_variable(nc)
    if t:
        retval["t"] = nc.variables[t].dimensions
    return retval


def is_dataset_valid_ragged_array_repr_featureType(nc, feature_type: str):
    """
    Check if a data set is a valid representation of a ragged
    array structure. See inline comments.
    """

    # regardless of if compound type or not, must have a cf_role
    # variable; if compound, this will be the first part of the
    # feature_type as we'll have to search for one with profile_id
    # regardless; if single feature type, cf_role must match that
    # featureType
    cf_role_vars = nc.get_variables_by_attributes(cf_role=lambda x: x is not None)
    is_compound = False
    if feature_type.lower() in {"timeseriesprofile", "trajectoryprofile"}:
        is_compound = True
        ftype = feature_type.lower().split("profile")[0]
        if len(cf_role_vars) > 2:
            return False
    else:
        ftype = feature_type.lower()
        if len(cf_role_vars) > 1 or not ftype:
            return False

    cf_role_var = nc.get_variables_by_attributes(cf_role=f"{ftype}_id")[0]
    # if cf_role_var returns None, this should raise an error?
    if cf_role_var.cf_role.split("_id")[0].lower() != ftype:
        return False

    # now we'll check dimensions for singular feature types and/or
    # the first half of the compound featureType
    instance_dim = cf_role_var.dimensions
    if len(instance_dim) != 1:
        return False

    # Now we check for the presence of an index variable or count variable;
    # NOTE that if no index or count variables exist, we can't determine with
    # certainty that this is invalid, because single-instance data sets
    # are valid representations of the ragged array structures. Instead,
    # if the index/count variable is present, we check that only one of
    # each is present and that their dimensions are correct
    index_vars = nc.get_variables_by_attributes(
        instance_dimension=lambda x: x is not None,
    )
    count_vars = nc.get_variables_by_attributes(
        sample_dimension=lambda x: x is not None,
    )

    # if the featureType isn't compound, shouldn't have both count and index
    if index_vars and count_vars and not is_compound:
        return False

    # single featureType, checking for valid index variable
    elif index_vars and not is_compound:
        if len(index_vars) > 1:
            return False
        # the index variable's attr 'instance_dimension'
        # must be the same as the actual instance dimension,
        # which we get from the cf_role variable
        if index_vars[0].instance_dimension != instance_dim[0]:
            return False

    # single featureType, checking for valid count variable
    elif count_vars and not is_compound:
        if len(count_vars) > 1:
            return False
        # the count variable must have the same dimensions
        # as the instance variable, which has the instance
        # dimension as its dimension
        if count_vars[0].dimensions != instance_dim:
            return False

    # Now, if the featureType is compound, an index variable
    # must be present for the profile variable. To verify this, we will
    # check that the dimension of the index variable is the same dimension
    # that is present on the variable which has the attribute cf_role=profile_id.
    # The attribute of the index variable 'instance_dimension' should point to the
    # name of the dimension of the cf_role variable for either timeSeries or trajectory.
    # A count variable must also be present, and should have the same dimension,
    # but its attribute 'sample_dimension' must refer to the dimension, which is
    # DIFFERENT than the variable with the attribute cf_role=ftype, where ftype is the
    # first half of the compound featureType (so either timeseries or trajectory).
    # Thus, the dimension of the count variable must be the same dimension as the
    # dimension that all the other geophysical variables have.
    elif index_vars and count_vars and is_compound:
        if len(index_vars) > 1 or len(count_vars) > 1:
            return False

        profile_cf_role_vars = nc.get_variables_by_attributes(cf_role="profile_id")
        if len(profile_cf_role_vars) > 1:
            return False
        profile_cf_role_var = profile_cf_role_vars[0]

        # we first check the dimension of the index variable
        if index_vars[0].dimensions != profile_cf_role_var.dimensions:
            return False

        # the attribute 'instance_dimension' must point to the dimension
        # of the timeseries or trajectory cf_role var
        if index_vars[0].instance_dimension != cf_role_var.dimensions[0]:
            return False

        # get all geophysical dims
        geophysical_dims = [nc[v].dimensions for v in get_geophysical_variables(nc)]
        if len(set(geophysical_dims)) != 1:
            return False

        # check the dimension of the count var is the same and that the
        # sample_dimension attribute points to the same dimension that
        # the geophysical variables have
        if (
            count_vars[0].dimensions != profile_cf_role_var.dimensions
            or (count_vars[0].sample_dimension,) != geophysical_dims[0]
        ):
            return False

    else:
        return False

    return True


def resolve_ragged_array_dimension(ds: Dataset):
    # TODO: put in loop?
    ragged_variable = ds.get_variables_by_attributes(
        sample_dimension=lambda s: isinstance(s, str),
    )
    if ragged_variable:
        ragged_type = "sample_dimension"
    else:
        ragged_variable = ds.get_variables_by_attributes(
            instance_dimension=lambda s: isinstance(s, str),
        )
        ragged_type = "instance_dimension"
    if ragged_variable is None:
        raise ValueError("Could not find a ragged array related variable")
    return ragged_type


def is_variable_valid_ragged_array_repr_featureType(nc, variable: str) -> bool:
    """
    This method returns a boolean indicating whether the variable
    is a valid member of a contiguous ragged array representation
    or indexed ragged array representation of any of the three singular
    CF featureType types.

    For any ragged array representation, any DATA VARIABLE must have
    the sample dimension as its sole dimension. Additionally, for any
    featureType or compound featureType (e.g. timeSeriesProfile), the
    data variable must have the sample dimension as its dimension.
    """

    # Get all geophysical variables; should have only one
    # dimension in the set, and the dimension of the variable
    # should be equal.
    geo_vars = get_geophysical_variables(nc)
    dim_tuples = [nc.variables[v].dimensions for v in geo_vars]
    if len(set(dim_tuples)) < 1:
        return False

    # NOTE
    # Each dimension tuple - there should be only one - should only
    # have the sample dimension as its sole value. If there are more
    # than one, we assume the sample dimension is the first. Is this
    # an appropriate assumption?
    dim = dim_tuples[0]
    if len(dim) != 1:
        return False

    # this is the only thing we have to work with
    return nc.variables[variable].dimensions == dim


def is_point(nc, variable):
    """
    Returns true if the variable is a point feature type

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(o), y(o), z(o), t(o)
    # X(o)

    dims = nc.variables[variable].dimensions

    cmatrix = coordinate_dimension_matrix(nc)
    first_coord = None
    if "t" in cmatrix:
        first_coord = cmatrix["t"]
        if len(cmatrix["t"]) > 1:
            return False
    if "x" in cmatrix:
        if first_coord is None:
            first_coord = cmatrix["x"]
        if first_coord != cmatrix["x"]:
            return False
        if len(cmatrix["x"]) > 1:
            return False
    if "y" in cmatrix:
        if first_coord is None:
            first_coord = cmatrix["y"]
        if first_coord != cmatrix["y"]:
            return False
        if len(cmatrix["y"]) > 1:
            return False
    if "z" in cmatrix:
        if first_coord is None:
            first_coord = cmatrix["z"]
        if first_coord != cmatrix["z"]:
            return False
        if len(cmatrix["z"]) > 1:
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
    """
    Returns true if the variable is a time series feature type.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """

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
    if "x" in cmatrix:
        if len(cmatrix["x"]) != 0:
            return False
    if "y" in cmatrix:
        if len(cmatrix["y"]) != 0:
            return False
    if "z" in cmatrix:
        if len(cmatrix["z"]) != 0:
            return False

    return True


def is_multi_timeseries_orthogonal(nc, variable):
    """
    Returns true if the variable is a orthogonal multidimensional array
    representation of time series. For more information on what this means see
    CF 1.6 §H.2.1

    http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#_orthogonal_multidimensional_array_representation_of_time_series

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(i), y(i), z(i), t(o)
    # X(i, o)
    dims = nc.variables[variable].dimensions

    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "t"):
        if req not in cmatrix:
            return False
    if len(cmatrix["x"]) != 1 or cmatrix["x"] != cmatrix["y"]:
        return False
    if "z" in cmatrix and cmatrix["x"] != cmatrix["z"]:
        return False

    timevar = get_time_variable(nc)
    if cmatrix["t"] != (timevar,):
        return False

    i = cmatrix["x"][0]
    o = cmatrix["t"][0]
    if dims == (i, o):
        return True
    return False


def is_multi_timeseries_incomplete(nc, variable):
    """
    Returns true if the variable is an incomplete multidimensional array
    representation of time series. For more information on what this means see
    CF 1.6 §H.2.2

    http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#_incomplete_multidimensional_array_representation_of_time_series

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """

    # x(i), y(i), z(i), t(i, o)
    # X(i, o)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "t"):
        if req not in cmatrix:
            return False
    if len(cmatrix["x"]) != 1:
        return False
    if cmatrix["x"] != cmatrix["y"]:
        return False
    if len(cmatrix["t"]) != 2:
        return False
    if cmatrix["x"][0] != cmatrix["t"][0]:
        return False

    i = cmatrix["x"][0]
    o = cmatrix["t"][1]

    if dims == (i, o):
        return True
    return False


def isTimeSeries(nc, variable):
    """
    Attempt to consolidate all of the disparate timeseries checks.
    Is this being pragmatic or lazy? I'll argue pragmatic.
    Confidently pragmatic.
    """

    # first three check if variable is a valid multidimensional
    # representation, the last checks if it's a valid ragged array
    if (
        is_timeseries(nc, variable)
        or is_multi_timeseries_orthogonal(nc, variable)
        or is_multi_timeseries_incomplete(nc, variable)
    ):
        return True

    return False


def is_cf_trajectory(nc, variable):
    """
    Returns true if the variable is a CF trajectory feature type

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(i, o), y(i, o), z(i, o), t(i, o)
    # X(i, o)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "t"):
        if req not in cmatrix:
            return False
    if len(cmatrix["x"]) != 2:
        return False
    if cmatrix["x"] != cmatrix["y"]:
        return False
    if cmatrix["x"] != cmatrix["t"]:
        return False
    if "z" in cmatrix and cmatrix["x"] != cmatrix["z"]:
        return False
    if dims == cmatrix["x"]:
        return True
    return False


def is_single_trajectory(nc, variable):
    """
    Returns true if the variable is a single trajectory feature

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(o), y(o), z(o), t(o)
    # X(o)
    # cf_role must be trajectory
    dims = nc.variables[variable].dimensions

    cmatrix = coordinate_dimension_matrix(nc)
    # Time is required for trajectories
    if "t" not in cmatrix:
        return False
    # Coordinates changing with time are optional
    if "x" in cmatrix:
        if cmatrix["x"] != cmatrix["t"]:
            return False
    if "y" in cmatrix:
        if cmatrix["y"] != cmatrix["t"]:
            return False
    if "z" in cmatrix:
        if cmatrix["z"] != cmatrix["t"]:
            return False
    if dims != cmatrix["t"]:
        return False
    traj_ids = nc.get_variables_by_attributes(cf_role="trajectory_id")
    if len(traj_ids) != 1:
        return False
    return True


def isTrajectory(nc, variable):
    """
    Wrapper method for checking if a variable is detected as
    a trajectory featureType.
    """

    if is_cf_trajectory(nc, variable) or is_single_trajectory(nc, variable):
        return True

    return False


def is_profile_orthogonal(nc, variable):
    """
    Returns true if the variable is a orthogonal profile feature type

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # Every profile has the exact same depths, think thermister or ADCP
    # x(i), y(i), z(j), t(i)
    # X(i, j)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "z", "t"):
        if req not in cmatrix:
            return False
    if len(cmatrix["x"]) != 1:
        return False
    if cmatrix["x"] != cmatrix["y"]:
        return False
    if cmatrix["x"] != cmatrix["t"]:
        return False
    if len(cmatrix["z"]) != 1:
        return False

    i = cmatrix["x"][0]
    j = cmatrix["z"][0]

    if dims == (i, j):
        return True
    return False


def is_profile_incomplete(nc, variable):
    """
    Returns true if the variable is a incomplete profile feature type

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # Every profile may have different depths
    # x(i), y(i), z(i, j), t(i)
    # X(i, j)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "z", "t"):
        if req not in cmatrix:
            return False
    if len(cmatrix["x"]) != 1:
        return False
    if cmatrix["x"] != cmatrix["y"]:
        return False
    if cmatrix["x"] != cmatrix["t"]:
        return False
    if len(cmatrix["z"]) != 2:
        return False
    if cmatrix["z"][0] != cmatrix["x"][0]:
        return False

    i = cmatrix["x"][0]
    j = cmatrix["z"][1]

    if dims == (i, j):
        return True
    return False


def isProfile(nc, variable: str):
    """
    Per Ch 9 of the CF spec, profiles are a part of the CF Discrete
    Sampling Geometries. Profile data can be logically represented
    in a file one of four ways: orthogonal multidimensional array,
    incomplete multidimensional array, contiguous ragged array,
    and indexed ragged array. If the variable is found to be any
    valid for any one of these representations, return "profile"
    else return None.

    The very first part of this function attempts to use legacy code
    to test if the variable is a profile in the orthogonal multidimensional
    array representation or the incomplete multidimensional array
    representation. If neither of these are true, it moves on to testing
    for the contiguous ragged array and indexed ragged array representations.

    Parameters
    ----------
    nc      : netCDF4 Dataset
    variable: str name of variable

    Returns
    -------
    str or None
    """

    # NOTE
    # TODO
    # Does this take into account a single profile? This is a valid profile.

    # first check for orthogonal, incomplete
    if is_profile_orthogonal(nc, variable) or is_profile_incomplete(nc, variable):
        return True


def is_timeseries_profile_single_station(nc, variable):
    """
    Returns true if the variable is a time-series profile that represents a
    single station and each profile is the same length.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """

    # x, y, z(z), t(t)
    # X(t, z)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "z", "t"):
        if req not in cmatrix:
            return False
    if len(cmatrix["x"]) != 0:
        return False
    if cmatrix["x"] != cmatrix["y"]:
        return False

    z = get_z_variable(nc)
    if cmatrix["z"] != (z,):
        return False
    t = get_time_variable(nc)
    if cmatrix["t"] != (t,):
        return False

    if dims == (t, z):
        return True
    return False


def is_timeseries_profile_multi_station(nc, variable):
    """
    Returns true if the variable is a time-series profile that represents multiple stations with orthogonal time and depth

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(i), y(i), z(z), t(t)
    # X(i, t, z)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "z", "t"):
        if req not in cmatrix:
            return False
    if len(cmatrix["x"]) != 1:
        return False
    if cmatrix["x"] != cmatrix["y"]:
        return False
    i = cmatrix["x"][0]

    z = get_z_variable(nc)
    if cmatrix["z"] != (z,):
        return False
    t = get_time_variable(nc)
    if cmatrix["t"] != (t,):
        return False

    if dims == (i, t, z):
        return True
    return False


def is_timeseries_profile_single_ortho_time(nc, variable):
    """
    Returns true if the variable is a time-series profile that represents a
    single station with orthogonal time only.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x, y, z(t, j), t(t)
    # X(t, j)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "z", "t"):
        if req not in cmatrix:
            return False

    if len(cmatrix["x"]) != 0:
        return False
    if cmatrix["x"] != cmatrix["y"]:
        return False

    t = get_time_variable(nc)
    if cmatrix["t"] != (t,):
        return False

    if len(cmatrix["z"]) != 2:
        return False

    if cmatrix["z"][0] != t:
        return False

    j = cmatrix["z"][1]

    if dims == (t, j):
        return True
    return False


def is_timeseries_profile_multi_ortho_time(nc, variable):
    """
    Returns true if the variable is a time-series profile that represents a
    multi station with orthogonal time only.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(i), y(i), z(i, t, j), t(t)
    # X(i, t, j)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "z", "t"):
        if req not in cmatrix:
            return False

    if len(cmatrix["x"]) != 1:
        return False
    if cmatrix["x"] != cmatrix["y"]:
        return False

    t = get_time_variable(nc)
    if cmatrix["t"] != (t,):
        return False

    if len(cmatrix["z"]) != 3:
        return False

    if cmatrix["z"][1] != t:
        return False
    if cmatrix["z"][0] != cmatrix["x"][0]:
        return False

    i = cmatrix["x"][0]
    j = cmatrix["z"][2]

    if dims == (i, t, j):
        return True
    return False


def is_timeseries_profile_ortho_depth(nc, variable):
    """
    Returns true if the variable is a time-series profile with orthogonal depth
    only.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(i), y(i), z(z), t(i, j)
    # X(i, j, z)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "z", "t"):
        if req not in cmatrix:
            return False

    if len(cmatrix["x"]) != 1:
        return False
    if cmatrix["x"] != cmatrix["y"]:
        return False

    z = get_z_variable(nc)
    if cmatrix["z"] != (z,):
        return False

    i = cmatrix["x"][0]

    if len(cmatrix["t"]) != 2:
        return False
    if cmatrix["t"][0] != i:
        return False

    j = cmatrix["t"][1]

    if dims == (i, j, z):
        return True
    return False


def is_timeseries_profile_incomplete(nc, variable):
    """
    Returns true if the variable is a time-series profile incomplete depth and
    incomplete time.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(i), y(i), z(i, j, k), t(i, j)
    # X(i, j, k)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "z", "t"):
        if req not in cmatrix:
            return False

    if len(cmatrix["x"]) != 1:
        return False
    if cmatrix["x"] != cmatrix["y"]:
        return False
    i = cmatrix["x"][0]

    if len(cmatrix["t"]) != 2:
        return False
    if cmatrix["t"][0] != i:
        return False
    j = cmatrix["t"][1]

    if len(cmatrix["z"]) != 3:
        return False
    if cmatrix["z"][0] != i:
        return False
    if cmatrix["z"][1] != j:
        return False
    k = cmatrix["z"][2]

    if dims == (i, j, k):
        return True
    return False


def isTimeSeriesProfile(nc, variable):
    """
    Wrapper method.
    Verify if a variable matches with the timeSeriesProfile
    featureType. According to the CF specification, a data set
    with timeSeriesProfile features has two cf_role variables:
    one for the "station" (cf_role=timeseries_id) and one for the
    "profile" (cf_role=profile_id).
    """

    if (
        is_timeseries_profile_single_station(nc, variable)
        or is_timeseries_profile_multi_station(nc, variable)
        or is_timeseries_profile_single_ortho_time(nc, variable)
        or is_timeseries_profile_multi_ortho_time(nc, variable)
        or is_timeseries_profile_ortho_depth(nc, variable)
        or is_timeseries_profile_incomplete(nc, variable)
    ):
        return True

    return False


def is_trajectory_profile_orthogonal(nc, variable):
    """
    Returns true if the variable is a trajectory profile with orthogonal
    depths.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(i, o), y(i, o), z(z), t(i, o)
    # X(i, o, z)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "z", "t"):
        if req not in cmatrix:
            return False

    if len(cmatrix["x"]) != 2:
        return False
    if cmatrix["x"] != cmatrix["y"]:
        return False
    if cmatrix["x"] != cmatrix["t"]:
        return False

    i, o = cmatrix["x"]

    z = get_z_variable(nc)
    if cmatrix["z"] != (z,):
        return False

    if dims == (i, o, z):
        return True
    return False


def is_trajectory_profile_incomplete(nc, variable):
    """
    Returns true if the variable is a trajectory profile with incomplete
    depths.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(i, o), y(i, o), z(i, o, j), t(i, o)
    # X(i, o, j)
    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "z", "t"):
        if req not in cmatrix:
            return False

    if len(cmatrix["x"]) != 2:
        return False
    if cmatrix["x"] != cmatrix["y"]:
        return False
    if cmatrix["x"] != cmatrix["t"]:
        return False

    i, o = cmatrix["x"]

    if len(cmatrix["z"]) != 3:
        return False

    if cmatrix["z"][0] != i:
        return False
    if cmatrix["z"][1] != o:
        return False

    j = cmatrix["z"][2]

    if dims == (i, o, j):
        return True
    return False


def isTrajectoryProfile(nc, variable):
    """
    Wrapper method
    """

    # NOTE
    # does this take into account single trajectory profile?
    if is_trajectory_profile_orthogonal(
        nc,
        variable,
    ) or is_trajectory_profile_incomplete(nc, variable):
        return True

    return False


def is_2d_regular_grid(nc, variable):
    """
    Returns True if the variable is a 2D Regular grid.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(x), y(y), t(t)
    # X(t, y, x)

    if is_mapped_grid(nc, variable):
        return False

    dims = nc.variables[variable].dimensions

    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "t"):
        if req not in cmatrix:
            return False

    x = get_lon_variable(nc)
    y = get_lat_variable(nc)
    t = get_time_variable(nc)

    if cmatrix["x"] != (x,):
        return False
    if cmatrix["y"] != (y,):
        return False
    if cmatrix["t"] != (t,):
        return False

    # Relaxed dimension ordering
    if len(dims) == 3 and x in dims and y in dims and t in dims:
        return True
    return False


def is_2d_static_grid(nc, variable):
    """
    Returns True if the variable is a 2D Regular grid that does not vary with
    time.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(x), y(y)
    # X(y, x)

    if is_mapped_grid(nc, variable):
        return False

    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y"):
        if req not in cmatrix:
            return False

    x = get_lon_variable(nc)
    y = get_lat_variable(nc)

    if cmatrix["x"] != (x,):
        return False

    if cmatrix["y"] != (y,):
        return False

    if len(dims) != 2 or x not in dims or y not in dims:
        return False

    return True


def is_3d_regular_grid(nc, variable):
    """
    Returns True if the variable is a 3D Regular grid.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(x), y(y), z(z), t(t)
    # X(t, z, y, x)

    if is_mapped_grid(nc, variable):
        return False

    dims = nc.variables[variable].dimensions

    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "z", "t"):
        if req not in cmatrix:
            return False

    x = get_lon_variable(nc)
    y = get_lat_variable(nc)
    z = get_z_variable(nc)
    t = get_time_variable(nc)

    if cmatrix["x"] != (x,):
        return False
    if cmatrix["y"] != (y,):
        return False
    if cmatrix["z"] != (z,):
        return False
    if cmatrix["t"] != (t,):
        return False

    # Relaxed dimension ordering
    if len(dims) == 4 and x in dims and y in dims and t in dims and z in dims:
        return True
    return False


def is_3d_static_grid(nc, variable):
    """
    Returns True if the variable is a 2D Regular grid that does not vary with
    time.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(x), y(y), z(z)
    # X(z, y, x)

    if is_mapped_grid(nc, variable):
        return False

    dims = nc.variables[variable].dimensions
    cmatrix = coordinate_dimension_matrix(nc)

    for req in ("x", "y", "z"):
        if req not in cmatrix:
            return False

    x = get_lon_variable(nc)
    y = get_lat_variable(nc)
    z = get_z_variable(nc)

    if cmatrix["x"] != (x,):
        return False

    if cmatrix["y"] != (y,):
        return False

    if cmatrix["z"] != (z,):
        return False

    if len(dims) != 3 or x not in dims or y not in dims or z not in dims:
        return False

    return True


def is_mapped_grid(nc, variable):
    """
    Returns true if the feature-type of variable corresponds to a mapped grid
    type. Characterized by Appedix F of CF-1.6

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    # x(j, i), y(j, i), z?, t?
    # F(t?, z?, j, i)

    # The important and really defining characteristic of mapped grids is that
    # the true latitude and longitude coordinates are functions of (j,i) and
    # that the geophysical variables are also functions of (j,i) in their
    # dimensions.
    dims = nc.variables[variable].dimensions
    # For cases like ROMS, the coordinates are mapped using the coordinates attribute
    variable_coordinates = getattr(nc.variables[variable], "coordinates", "").split()

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

    comma_dimension = ",".join(dims)
    # Dimensions must be in the same order and the mapping coordinates i and j
    # must be in the same order
    if ",".join(x) not in comma_dimension:
        return False

    return True


def is_reduced_grid(nc, variable):
    """
    Returns True if the feature-type of the variable corresponds to a reduced
    horizontal grid.

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    axis_map = get_axis_map(nc, variable)

    if "X" not in axis_map:
        return False
    if "Y" not in axis_map:
        return False
    if "C" not in axis_map:
        return False

    compressed_coordinates = axis_map["C"]
    if len(compressed_coordinates) > 1:
        return False
    compressed_coordinate = axis_map["C"][0]
    for dim in nc.variables[compressed_coordinate].compress.split():
        if dim not in nc.dimensions:
            return False
    return True


def guess_feature_type(nc, variable):
    """
    Returns a string describing the feature type for this variable

    :param netCDF4.Dataset nc: An open netCDF dataset
    :param str variable: name of the variable to check
    """
    if is_point(nc, variable):
        return "point"
    if isProfile(nc, variable):
        return "profile"
    if isTimeSeries(nc, variable):
        return "timeseries"
    if isTrajectory(nc, variable):
        return "trajectory"
    if isTimeSeriesProfile(nc, variable):
        return "timeseriesprofile"
    if isTrajectoryProfile(nc, variable):
        return "trajectoryprofile"

    # TODO
    # consolidate below into "isGrid" ?
    if is_2d_regular_grid(nc, variable):
        return "2d-regular-grid"
    if is_2d_static_grid(nc, variable):
        return "2d-static-grid"
    if is_3d_regular_grid(nc, variable):
        return "3d-regular-grid"
    if is_3d_static_grid(nc, variable):
        return "3d-static-grid"
    if is_mapped_grid(nc, variable):
        return "mapped-grid"
    if is_reduced_grid(nc, variable):
        return "reduced-grid"


def units_convertible(units1, units2):
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
