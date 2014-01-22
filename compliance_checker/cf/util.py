from netCDF4 import Dimension, Variable

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
