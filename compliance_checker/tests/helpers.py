import tempfile

from netCDF4 import Dataset


class MockNetCDF(Dataset):
    """
    Wrapper object around NetCDF Dataset to write data only to memory.
    """

    def __init__(self):
        # taken from test/tst_diskless.py NetCDF library
        # even though we aren't persisting data to disk, the constructor
        # requires a filename not currently in use by another Dataset object..
        tmp_filename = tempfile.NamedTemporaryFile(suffix=".nc", delete=True).name
        super(MockNetCDF, self).__init__(
            tmp_filename, "w", diskless=True, persist=False
        )


class MockTimeSeries(MockNetCDF):
    """
    Mock time series with time dimension and time, lon, lat, and depth
    variables defined
    """

    def __init__(self, default_fill_value=None):
        super(MockTimeSeries, self).__init__()
        self.createDimension("time", 500)
        for v in ("time", "lon", "lat", "depth"):
            self.createVariable(v, "d", ("time",), fill_value=default_fill_value)

        # give some applicable units
        self.variables["time"].units = "seconds since 2019-04-11T00:00:00"
        self.variables["time"].axis = "T"
        self.variables["lat"].units = "degree_north"
        self.variables["lat"].axis = "Y"
        self.variables["lon"].units = "degree_east"
        self.variables["lon"].axis = "X"
        self.variables["depth"].units = "meters"
        self.variables["depth"].axis = "Z"
        self.variables["depth"].positive = "down"


class MockVariable(object):
    """
    For mocking a dataset variable. Constructor optionally takes a NetCDF
    variable, the NetCDF attributes of which will be copied over to this
    object.
    """

    def __init__(self, copy_var=None):
        if copy_var is not None:
            self.name = copy_var.name
            self.dtype = copy_var.dtype
            self.dimensions = copy_var.dimensions
            self.ndim = copy_var.ndim
            for att in copy_var.ncattrs():
                setattr(self, att, getattr(copy_var, att))

    def ncattrs(self):
        return [
            att
            for att in vars(self)
            if att not in {"ndim", "name", "dtype", "dimensions"}
        ]
