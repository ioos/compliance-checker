import itertools
import requests
import os
from copy import deepcopy
from collections import defaultdict
from lxml import etree
from cf_units import Unit
from netCDF4 import Dimension, Variable
from pkgutil import get_data
from pkg_resources import resource_filename

# copied from paegan
# paegan may depend on these later
_possiblet = ["time", "TIME", "Time",
              "t", "T",
              "ocean_time", "OCEAN_TIME",
              "jd", "JD",
              "dn", "DN",
              "times", "TIMES", "Times",
              "mt", "MT",
              "dt", "DT",
              ]
_possiblez = ["depth", "DEPTH",
              "depths", "DEPTHS",
              "height", "HEIGHT",
              "altitude", "ALTITUDE",
              "alt", "ALT",
              "Alt", "Altitude",
              "h", "H",
              "s_rho", "S_RHO",
              "s_w", "S_W",
              "z", "Z",
              "siglay", "SIGLAY",
              "siglev", "SIGLEV",
              "sigma", "SIGMA",
              "vertical", "VERTICAL", "lev", "LEV", "level", "LEVEL"
              ]
_possiblex = ["x", "X",
              "lon", "LON",
              "xlon", "XLON",
              "lonx", "lonx",
              "lon_u", "LON_U",
              "lon_v", "LON_V",
              "lonc", "LONC",
              "Lon", "Longitude",
              "longitude", "LONGITUDE",
              "lon_rho", "LON_RHO",
              "lon_psi", "LON_PSI",

              ]
_possibley = ["y", "Y",
              "lat", "LAT",
              "ylat", "YLAT",
              "laty", "laty",
              "lat_u", "LAT_U",
              "lat_v", "LAT_V",
              "latc", "LATC",
              "Lat", "Latitude",
              "latitude", "LATITUDE",
              "lat_rho", "LAT_RHO",
              "lat_psi", "LAT_PSI",

              ]

_possibleaxis = _possiblet + _possiblez + _possiblex + _possibley


_possiblexunits = ['degrees_east',
                   'degree_east',
                   'degrees_E',
                   'degree_E',
                   'degreesE',
                   'degreeE'
                   ]

_possibleyunits = ['degrees_north',
                   'degree_north',
                   'degrees_N',
                   'degree_N',
                   'degreesN',
                   'degreeN'
                   ]

_possibletunits = ['day',
                   'days',
                   'd',
                   'hour',
                   'hours',
                   'hr',
                   'hrs',
                   'h',
                   'year',
                   'years',
                   'minute',
                   'minutes',
                   'm',
                   'min',
                   'mins',
                   'second',
                   'seconds',
                   's',
                   'sec',
                   'secs'
                   ]

_possibleaxisunits = _possiblexunits + _possibleyunits + _possibletunits


class DotDict(dict):
    """
    Subclass of dict that will recursively look up attributes with dot notation.
    This is primarily for working with JSON-style data in a cleaner way like javascript.
    Note that this will instantiate a number of child DotDicts when you first access attributes;
    do not use in performance-critical parts of your code.
    """

    def __dir__(self):
        return list(self.__dict__.keys()) + list(self.keys())

    def __getattr__(self, key):
        """ Make attempts to lookup by nonexistent attributes also attempt key lookups. """
        if key in self:
            return self[key]
        import sys
        import dis
        frame = sys._getframe(1)
        if '\x00%c' % dis.opmap['STORE_ATTR'] in frame.f_code.co_code:
            self[key] = DotDict()
            return self[key]

        raise AttributeError(key)

    def __setattr__(self, key, value):
        if key in dir(dict):
            raise AttributeError('%s conflicts with builtin.' % key)
        if isinstance(value, dict):
            self[key] = DotDict(value)
        else:
            self[key] = value

    def copy(self):
        return deepcopy(self)

    def get_safe(self, qual_key, default=None):
        """
        @brief Returns value of qualified key, such as "system.name" or None if not exists.
                If default is given, returns the default. No exception thrown.
        """
        value = get_safe(self, qual_key)
        if value is None:
            value = default
        return value

    @classmethod
    def fromkeys(cls, seq, value=None):
        return DotDict(dict.fromkeys(seq, value))


def get_safe(dict_instance, keypath, default=None):
    """
    Returns a value with in a nested dict structure from a dot separated
    path expression such as "system.server.host" or a list of key entries
    @retval Value if found or None
    """
    try:
        obj = dict_instance
        keylist = keypath if type(keypath) is list else keypath.split('.')
        for key in keylist:
            obj = obj[key]
        return obj
    except Exception:
        return default


class NCGraph(object):

    def __init__(self, ds, name, nc_object, self_reference_variables, reference_map=None):

        self.ds           = ds
        self.name         = name
        self.coords       = DotDict()
        self.dims         = DotDict()
        self.grid_mapping = DotDict()
        self.obj          = nc_object

        self.reference_variables = self_reference_variables
        self.reference_map = reference_map or {}

        self.reference_map[name] = self

        if isinstance(nc_object, Dimension):
            self._type = 'dim'

        elif isinstance(nc_object, Variable):
            self._type = 'var'

            self.get_references()

        else:
            raise TypeError("unknown type %s" % repr(type(nc_object)))

    def get_references(self):
        for dim in self.obj.dimensions:
            self.dims[dim] = self.get_dimension(dim)

        if hasattr(self.obj, 'coordinates'):
            coords = self.obj.coordinates.split(' ')
            for coord in coords:
                self.coords[coord] = self.get_coordinate(coord)

        if hasattr(self.obj, 'grid_mapping'):
            gm = self.obj.grid_mapping
            self.grid_mapping[gm] = self.get_grid_mapping(gm)

    def get_dimension(self, dim):
        if dim in self.reference_map:
            return self.reference_map[dim]
        return NCGraph(self.ds, dim, self.ds.dimensions[dim], self.reference_variables, self.reference_map)

    def get_coordinate(self, coord):
        if coord not in self.ds.variables:
            return
        if coord in self.reference_map:
            if self.name == coord:
                self.reference_variables.add(self.name)
            return self.reference_map[coord]
        return NCGraph(self.ds, coord, self.ds.variables[coord], self.reference_variables, self.reference_map)

    def get_grid_mapping(self, gm):
        if gm not in self.ds.variables:
            return
        if gm in self.reference_map:
            return self.reference_map[gm]
        return NCGraph(self.ds, gm, self.ds.variables[gm], self.reference_variables, self.reference_map)

    def __getattr__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        return getattr(self.obj, key)


class StandardNameTable(object):

    class NameEntry(object):

        def __init__(self, entrynode):
            self.canonical_units = self._get(entrynode, 'canonical_units', True)
            self.grib            = self._get(entrynode, 'grib')
            self.amip            = self._get(entrynode, 'amip')
            self.description     = self._get(entrynode, 'description')

        def _get(self, entrynode, attrname, required=False):
            vals = entrynode.xpath(attrname)
            if len(vals) > 1:
                raise Exception("Multiple attrs (%s) found" % attrname)
            elif required and len(vals) == 0:
                raise Exception("Required attr (%s) not found" % attrname)

            return vals[0].text

    def __init__(self, cached_location=None):
        if cached_location:
            try:
                with open(cached_location, 'r') as fp:
                    resource_text = fp.read()
            except Exception:
                resource_text = get_data("compliance_checker", "data/cf-standard-name-table.xml")
        else:
            resource_text = get_data("compliance_checker", "data/cf-standard-name-table.xml")

        parser = etree.XMLParser(remove_blank_text=True)
        self._root = etree.fromstring(resource_text, parser)

        # generate and save a list of all standard names in file
        self._names = [node.get('id') for node in self._root.iter('entry')]
        self._aliases = [node.get('id') for node in self._root.iter('alias')]
        self._version = self._root.xpath('version_number')[0].text

    def __len__(self):
        return len(self._names) + len(self._aliases)

    def __getitem__(self, key):
        if not (key in self._names or key in self._aliases):
            raise KeyError("%s not found in standard name table" % key)

        if key in self._aliases:
            idx = self._aliases.index(key)
            entryids = self._root.xpath('alias')[idx].xpath('entry_id')

            if len(entryids) != 1:
                raise Exception("Inconsistency in standard name table, could not lookup alias for %s" % key)

            key = entryids[0].text

        if key not in self._names:
            raise KeyError("%s not found in standard name table" % key)

        idx = self._names.index(key)
        entry = self.NameEntry(self._root.xpath('entry')[idx])
        return entry

    def get(self, key, default=None):
        '''
        Returns the item for the key or returns the default if it does not exist
        '''
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        return key in self._names or key in self._aliases

    def __iter__(self):
        return iter(itertools.chain(self._names, self._aliases))


def download_cf_standard_name_table(version, location=None):
    '''
    Downloads the specified CF standard name table version and saves it to file

    :param str version: CF standard name table version number (i.e 34)
    :param str location: Path/filename to write downloaded xml file to
    '''

    if location is None:  # This case occurs when updating the packaged version from command line
        location = resource_filename('compliance_checker', 'data/cf-standard-name-table.xml')

    url = "http://cfconventions.org/Data/cf-standard-names/{0}/src/cf-standard-name-table.xml".format(version)
    r = requests.get(url, allow_redirects=True)
    if r.status_code == 200:
        print("Downloading cf-standard-names table version {0} from: {1}".format(version, url))
        with open(location, 'wb') as f:
            f.write(r.content)
    else:
        r.raise_for_status()
    return

def create_cached_data_dir():
    '''
    Returns the path to the data directory to download CF standard names.
    Use $XDG_DATA_HOME.
    '''
    writable_directory = os.path.join(os.path.expanduser('~'), '.local', 'share')
    data_directory = os.path.join(os.environ.get("XDG_DATA_HOME", writable_directory),
                                  'compliance-checker')
    if not os.path.isdir(data_directory):
        os.makedirs(data_directory)

    return data_directory

def units_known(units):
    try:
        Unit(units)
    except ValueError:
        return False
    return True


def units_convertible(units1, units2, reftimeistime=True):
    """Return True if a Unit representing the string units1 can be converted
    to a Unit representing the string units2, else False."""
    try:
        u1 = Unit(units1)
        u2 = Unit(units2)
    except ValueError:
        return False
    return u1.is_convertible(u2)


def units_temporal(units):
    try:
        u = Unit(units)
    except ValueError:
        return False
    return u.is_time_reference()


def map_axes(dim_vars, reverse_map=False):
    """
    axis name       -> [dimension names]
    dimension name  -> [axis_name], length 0 if reverse_map
    """
    ret_val = defaultdict(list)
    axes = ['X', 'Y', 'Z', 'T']

    for k, v in dim_vars.items():
        axis = getattr(v, 'axis', '')
        if not axis:
            continue

        axis = axis.upper()
        if axis in axes:
            if reverse_map:
                ret_val[k].append(axis)
            else:
                ret_val[axis].append(k)

    return dict(ret_val)


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
    satisfied = varname.lower() == 'time'
    satisfied |= getattr(var, 'standard_name', '') == 'time'
    satisfied |= getattr(var, 'axis', '') == 'T'
    satisfied |= units_convertible('seconds since 1900-01-01', getattr(var, 'units', ''))
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
    satisfied |= getattr(var, 'standard_name', '') in _possiblez
    # Is the axis set to Z?
    satisfied |= getattr(var, 'axis', '').lower() == 'z'
    is_pressure = units_convertible(getattr(var, 'units', '1'), 'dbar')
    # Pressure defined or positive defined
    satisfied |= is_pressure
    if not is_pressure:
        satisfied |= getattr(var, 'positive', '').lower() in ('up', 'down')
    return satisfied
