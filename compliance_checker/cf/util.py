import itertools
import os
import posixpath
import sys
from functools import lru_cache
from importlib.resources import files
from pkgutil import get_data
from typing import Union

import requests
from cf_units import Unit
from lxml import etree
from netCDF4 import Dataset, Dimension, Group, Variable

from compliance_checker.cfutil import units_convertible

# copied from paegan
# paegan may depend on these later
_possiblet = {
    "time",
    "TIME",
    "Time",
    "t",
    "T",
    "ocean_time",
    "OCEAN_TIME",
    "jd",
    "JD",
    "dn",
    "DN",
    "times",
    "TIMES",
    "Times",
    "mt",
    "MT",
    "dt",
    "DT",
}
_possiblez = {
    "depth",
    "DEPTH",
    "depths",
    "DEPTHS",
    "height",
    "HEIGHT",
    "altitude",
    "ALTITUDE",
    "alt",
    "ALT",
    "Alt",
    "Altitude",
    "h",
    "H",
    "s_rho",
    "S_RHO",
    "s_w",
    "S_W",
    "z",
    "Z",
    "siglay",
    "SIGLAY",
    "siglev",
    "SIGLEV",
    "sigma",
    "SIGMA",
    "vertical",
    "VERTICAL",
    "lev",
    "LEV",
    "level",
    "LEVEL",
}
_possiblex = {
    "x",
    "X",
    "lon",
    "LON",
    "xlon",
    "XLON",
    "lonx",
    "lon_u",
    "LON_U",
    "lon_v",
    "LON_V",
    "lonc",
    "LONC",
    "Lon",
    "Longitude",
    "longitude",
    "LONGITUDE",
    "lon_rho",
    "LON_RHO",
    "lon_psi",
    "LON_PSI",
}
_possibley = {
    "y",
    "Y",
    "lat",
    "LAT",
    "ylat",
    "YLAT",
    "laty",
    "lat_u",
    "LAT_U",
    "lat_v",
    "LAT_V",
    "latc",
    "LATC",
    "Lat",
    "Latitude",
    "latitude",
    "LATITUDE",
    "lat_rho",
    "LAT_RHO",
    "lat_psi",
    "LAT_PSI",
}

_possibleaxis = _possiblet | _possiblez | _possiblex | _possibley


_possiblexunits = {
    "degrees_east",
    "degree_east",
    "degrees_E",
    "degree_E",
    "degreesE",
    "degreeE",
}

_possibleyunits = {
    "degrees_north",
    "degree_north",
    "degrees_N",
    "degree_N",
    "degreesN",
    "degreeN",
}

_possibletunits = {
    "day",
    "days",
    "d",
    "hour",
    "hours",
    "hr",
    "hrs",
    "h",
    "year",
    "years",
    "minute",
    "minutes",
    "m",
    "min",
    "mins",
    "second",
    "seconds",
    "s",
    "sec",
    "secs",
}

_possibleaxisunits = _possiblexunits | _possibleyunits | _possibletunits


def get_safe(dict_instance, keypath, default=None):
    """
    Returns a value with in a nested dict structure from a dot separated
    path expression such as "system.server.host" or a list of key entries
    @retval Value if found or None
    """
    try:
        obj = dict_instance
        keylist = keypath if isinstance(keypath, list) else keypath.split(".")
        for key in keylist:
            obj = obj[key]
        return obj
    except Exception:
        return default


class VariableReferenceError(Exception):
    """A variable to assign bad variable references to"""

    def __init__(self, name: str, dataset: Dataset = None):
        self.name = name
        self.dataset_path = dataset.filepath() if dataset is not None else None

    def __str__(self):
        return (
            f"Cannot find variable named {self.name} in dataset " f"{self}.dataset_path"
        )


class StandardNameTable:
    class NameEntry:
        def __init__(self, entrynode):
            self.canonical_units = self._get(entrynode, "canonical_units", True)
            self.grib = self._get(entrynode, "grib")
            self.amip = self._get(entrynode, "amip")
            self.description = self._get(entrynode, "description")

        def _get(self, entrynode, attrname, required=False):
            vals = entrynode.xpath(attrname)
            if len(vals) > 1:
                raise Exception(f"Multiple attrs ({attrname}) found")
            elif required and len(vals) == 0:
                raise Exception(f"Required attr ({attrname}) not found")

            return vals[0].text if len(vals) > 0 else None

    def __init__(self, cached_location=None):
        if cached_location:
            with open(cached_location, encoding="utf-8") as fp:
                resource_text = fp.read()
        elif os.environ.get("CF_STANDARD_NAME_TABLE") and os.path.exists(
            os.environ["CF_STANDARD_NAME_TABLE"],
        ):
            with open(
                os.environ["CF_STANDARD_NAME_TABLE"],
                encoding="utf-8",
            ) as fp:
                resource_text = fp.read()
        else:
            resource_text = get_data(
                "compliance_checker",
                "data/cf-standard-name-table.xml",
            )

        parser = etree.XMLParser(remove_blank_text=True)
        self._root = etree.fromstring(resource_text, parser)

        # generate and save a list of all standard names in file
        self._names = [node.get("id") for node in self._root.iter("entry")]
        self._aliases = [node.get("id") for node in self._root.iter("alias")]
        self._version = self._root.xpath("version_number")[0].text

    def __len__(self):
        return len(self._names) + len(self._aliases)

    def __getitem__(self, key):
        if not (key in self._names or key in self._aliases):
            raise KeyError(f"{key} not found in standard name table")

        if key in self._aliases:
            idx = self._aliases.index(key)
            entryids = self._root.xpath("alias")[idx].xpath("entry_id")

            if len(entryids) != 1:
                raise Exception(
                    f"Inconsistency in standard name table, could not lookup alias for {key}",
                )

            key = entryids[0].text

        if key not in self._names:
            raise KeyError(f"{key} not found in standard name table")

        idx = self._names.index(key)
        entry = self.NameEntry(self._root.xpath("entry")[idx])
        return entry

    def get(self, key, default=None):
        """
        Returns the item for the key or returns the default if it does not exist
        """
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        return key in self._names or key in self._aliases

    def __iter__(self):
        return iter(itertools.chain(self._names, self._aliases))


def download_cf_standard_name_table(version, location=None):
    """
    Downloads the specified CF standard name table version and saves it to file

    :param str version: CF standard name table version number (i.e 34)
    :param str location: Path/filename to write downloaded xml file to
    """

    if (
        location is None
    ):  # This case occurs when updating the packaged version from command line
        location = files("compliance_checker") / "data/cf-standard-name-table.xml"

    if version == "latest":
        url = "http://cfconventions.org/Data/cf-standard-names/current/src/cf-standard-name-table.xml"
    else:
        url = f"http://cfconventions.org/Data/cf-standard-names/{version}/src/cf-standard-name-table.xml"

    r = requests.get(url, allow_redirects=True)
    r.raise_for_status()

    print(
        f"Downloading cf-standard-names table version {version} from: {url}",
        file=sys.stderr,
    )
    with open(location, "wb") as f:
        f.write(r.content)


def create_cached_data_dir():
    """
    Returns the path to the data directory to download CF standard names.
    Use $XDG_DATA_HOME.
    """
    writable_directory = os.path.join(os.path.expanduser("~"), ".local", "share")
    data_directory = os.path.join(
        os.environ.get("XDG_DATA_HOME", writable_directory),
        "compliance-checker",
    )
    if not os.path.isdir(data_directory):
        os.makedirs(data_directory)

    return data_directory


def units_known(units):
    try:
        Unit(units)
    except ValueError:
        return False
    return True


def units_temporal(units):
    try:
        u = Unit(units)
    except ValueError:
        return False
    # IMPLEMENTATION CONFORMANCE REQUIRED 4.4 1/3
    # time units of a time_coordinate variable must contain a reference
    # time/date
    # IMPLEMENTATION CONFORMANCE REQUIRED 4.4 3/3
    # check that reference time seconds is not greater than or
    # equal to 60
    return u.is_time_reference()


def find_coord_vars(ncds):
    """
    Finds all coordinate variables in a dataset.

    A variable with the same name as a dimension is called a coordinate variable.
    """
    coord_vars = []

    for d in ncds.dimensions:
        if d in ncds.variables and ncds.variables[d].dimensions == (d,):
            coord_vars.append(ncds.variables[d])

    return coord_vars


def is_time_variable(varname, var):
    """
    Identifies if a variable is represents time
    """
    satisfied = varname.lower() == "time"
    satisfied |= getattr(var, "standard_name", "") == "time"
    satisfied |= getattr(var, "axis", "") == "T"
    satisfied |= units_convertible(
        "seconds since 1900-01-01",
        getattr(var, "units", ""),
    )
    return satisfied


def is_vertical_coordinate(var_name, var):
    """
    Determines if a variable is a vertical coordinate variable

    4.3
    A vertical coordinate will be identifiable by: units of pressure; or the presence of the positive attribute with a
    value of up or down (case insensitive).  Optionally, the vertical type may be indicated additionally by providing
    the standard_name attribute with an appropriate value, and/or the axis attribute with the value Z.
    """
    # Known name
    satisfied = var_name.lower() in _possiblez
    satisfied |= getattr(var, "standard_name", "") in _possiblez
    # Is the axis set to Z?
    satisfied |= getattr(var, "axis", "").lower() == "z"
    is_pressure = units_convertible(getattr(var, "units", "1"), "dbar")
    # Pressure defined or positive defined
    satisfied |= is_pressure
    if not is_pressure:
        satisfied |= getattr(var, "positive", "").lower() in ("up", "down")
    return satisfied


def string_from_var_type(variable):
    if isinstance(variable, str):
        return variable[:]
    elif variable.dtype.kind == "S":
        strip_char = variable.fill_value or b"\x00"
        return variable.tobytes().rstrip(strip_char).decode("utf-8")
    else:
        raise TypeError(
            f"Variable '{variable.name} has non-string/character' "
            f"dtype {variable.dtype}",
        )


def get_possible_label_variable_dimensions(variable: Variable) -> tuple[int, ...]:
    """
    Return dimensions if non-char variable, or return variable dimensions
    without trailing dimension if char variable, treating it as a label variable.
    """
    if variable.dtype.kind == "C" and len(variable.dimensions) > 0:
        return variable.dimensions[:-1]
    return variable.dimensions


@lru_cache
def maybe_lateral_reference_variable_or_dimension(
    group: Union[Group, Dataset],
    name: str,
    reference_type: Union[Variable, Dimension],
):

    def can_lateral_search(name):
        return not name.startswith(".") and posixpath.split(name)[0] == ""

    if reference_type == "variable":
        # first try to fetch any
        # can't set to None with .get
        try:
            maybe_var = group[name]
        except IndexError:
            maybe_var = None
        else:
            if isinstance(maybe_var, Variable):
                return maybe_var

        # alphanumeric string by itself, not a relative or absolute
        # search by proximity
        if posixpath.split(name)[0] == "" and not (
            name.startswith(".") or name.startswith("/")
        ):
            group_traverse = group
            while group_traverse.parent:
                group_traverse = group_traverse.parent
                # check_target = posixpath.join(group_traverse.path, name)
                try:
                    maybe_var = group_traverse[name]
                except IndexError:
                    maybe_var = None
                else:
                    if isinstance(maybe_var, Variable):
                        return maybe_var
        else:
            return VariableReferenceError(name)

            # can't find path relative to current group or absolute path
            # perform lateral search if we aren't in the root group


def reference_attr_variables(
    dataset: Dataset,
    attributes_string: str,
    split_by: str = None,
    reference_type: str = "variable",
    group: Union[Group, Dataset] = None,
):
    """
    Attempts to reference variables in the string, optionally splitting by
    a string
    """
    references, errors = [], []
    if attributes_string is None:
        return None, None
    elif reference_type == "variable":
        if split_by is None:
            return_val = dataset.variables.get(
                attributes_string,
                VariableReferenceError(attributes_string),
            )
            if not isinstance(return_val, VariableReferenceError):
                return return_val, None
            else:
                return None, return_val
        else:
            string_proc = attributes_string.split(split_by)
            for var_name in string_proc:
                if var_name in dataset.variables:
                    references.append(dataset.variables[var_name])
                else:
                    errors.append(VariableReferenceError(var_name))
            return references, errors
