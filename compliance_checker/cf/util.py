import os
import os.path
import itertools
from lxml import etree
try:
    from udunitspy import Unit, UdunitsError, Converter
except ImportError:
    pass #disabled as CF is not working and udunits won't install on centos/rhel yet
from netCDF4 import Dimension, Variable
from pkgutil import get_data

class DotDict(dict):
    """
    Subclass of dict that will recursively look up attributes with dot notation.
    This is primarily for working with JSON-style data in a cleaner way like javascript.
    Note that this will instantiate a number of child DotDicts when you first access attributes;
    do not use in performance-critical parts of your code.
    """

    def __dir__(self):
        return self.__dict__.keys() + self.keys()

    def __getattr__(self, key):
        """ Make attempts to lookup by nonexistent attributes also attempt key lookups. """
        if self.has_key(key):
            return self[key]
        import sys
        import dis
        frame = sys._getframe(1)
        if '\x00%c' % dis.opmap['STORE_ATTR'] in frame.f_code.co_code:
            self[key] = DotDict()
            return self[key]

        raise AttributeError(key)

    def __setattr__(self,key,value):
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

class NCGraph:
    def __init__(self, ds, name, nc_object):
        self.name         = name
        self.coords       = DotDict()
        self.dims         = DotDict()
        self.grid_mapping = DotDict()
        self.obj          = nc_object
        if isinstance(nc_object, Dimension):
            self._type = 'dim'
        elif isinstance(nc_object, Variable):
            self._type = 'var'
            for dim in nc_object.dimensions:
                self.dims[dim] = NCGraph(ds, dim, ds.dimensions[dim])
            if hasattr(nc_object, 'coordinates'):
                coords = nc_object.coordinates.split(' ')
                for coord in coords:
                    if coord in ds.variables:
                        self.coords[coord] = NCGraph(ds, coord, ds.variables[coord])
                    else:
                        self.coords[coord] = None
            if hasattr(nc_object, 'grid_mapping'):
                gm = nc_object.grid_mapping
                self.grid_mapping[gm] = None
                if gm in ds.variables:
                    self.grid_mapping[gm] = NCGraph(ds, gm, ds.variables[gm])

        else:
            raise TypeError("unknown type %s" % repr(type(nc_object)))

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
                raise StandardError("Multiple attrs (%s) found" % attrname)
            elif required and len(vals) == 0:
                raise StandardError("Required attr (%s) not found" % attrname)

            return vals[0].text

    def __init__(self, resource_name):
        resource_text = get_data("compliance_checker", "data/cf-standard-name-table.xml")
        parser = etree.XMLParser(remove_blank_text=True)
        self._root = etree.fromstring(resource_text, parser)

        # generate and save a list of all standard names in file
        self._names = [node.get('id') for node in self._root.iter('entry')]
        self._aliases = [node.get('id') for node in self._root.iter('alias')]

    def __len__(self):
        return len(self._names) + len(self._aliases)

    def __getitem__(self, key):
        if not (key in self._names or key in self._aliases):
            raise KeyError("%s not found in standard name table" % key)

        if key in self._aliases:
            idx = self._aliases.index(key)
            entryids = self._root.xpath('alias')[idx].xpath('entry_id')

            if len(entryids) != 1:
                raise StandardError("Inconsistency in standard name table, could not lookup alias for %s" % key)

            key = entryids[0].text

        if not key in self._names:
            raise KeyError("%s not found in standard name table" % key)

        idx = self._names.index(key)
        entry = self.NameEntry(self._root.xpath('entry')[idx])
        return entry

    def __contains__(self, key):
        return key in self._names or key in self._aliases

    def __iter__(self):
        return iter(itertools.chain(self._names, self._aliases))

def units_known(units):
    try:
        Unit(str(units))
    except UdunitsError:
        return False
    return True

def units_convertible(units1, units2, reftimeistime=True):
    try:
        Converter(str(units1), str(units2))
    except UdunitsError:
        return False

    return True

def units_temporal(units):
    r = False
    try:
        u = Unit('seconds since 1900-01-01')
        r = u.are_convertible(str(units))
    except UdunitsError:
        return False
    return r

