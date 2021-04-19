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


class MockRaggedArrayRepr(MockNetCDF):
    """
    Class to construct a fake NetCDF dataset using the ragged
    array structure. User must specify whether the structure should
    be contiguous ragged array or incomplete ragged array.
    The user must also specify the featureType of the data set;
    compound featureTypes are acceptable, and all are case-insensitive.

    Users should note that if a compound featureType is given,
    only the valid-CF representation is available. This means
    profiles are organized by contiguous ragged array and an index
    variable is present to denote which "station" the profile belongs
    to.

    Data sets which follow the contiguous ragged array structure
    have a variable called the "count variable" denoting how many
    elements are in a particular instance. This variable must have the
    instance dimension as its dimension. The variable will have an
    attribute "sample_dimension" denoting which dimension is the
    sample dimension.

    Data sets which follow the indexed ragged array structure have
    a variable called the "index variable" denoting which instance
    a particular element belongs to. This variable will have an
    attribute "instance_dimension" denoting the sample dimension.
    The variable's dimension is the sample dimension.

    The netCDF4.Dataset interface available to the user will contain
    a very simple representation of the data model, yielding no
    extraneous dimensions/variables.

    Parameters
    ----------
    feature_type: str (
        profile, timeseries, trajectory,
        timeseriesprofile, trajectoryprofile)

    structure: str (contiguous|indexed)
    """

    def __init__(self, feature_type: str, structure="contiguous"):
        super(MockRaggedArrayRepr, self).__init__()

        if structure.lower() not in ("contiguous", "indexed"):
            raise ValueError("Must initialize MockRaggedArray as contiguous or indexed")

        if feature_type.lower() not in {
            "point",
            "profile",
            "timeseries",
            "trajectory",
            "timeseriesprofile",
            "trajectoryprofile",
        }:
            raise ValueError("Must initialize MockRaggedArray with valid featureType")

        is_compound = False
        if feature_type.lower() in {"timeseriesprofile", "trajectoryprofile"}:
            is_compound = True

        # the data will have 10 instances of whatever and 100
        # total elements
        self.createDimension("INSTANCE_DIMENSION", 10)
        self.createDimension("SAMPLE_DIMENSION", 100)

        # create a variable as the cf_role variable; if a compound
        # featureType is given, multiple variables are created
        if is_compound:
            # compound, need another dimension as well; this will
            # become the "station dimension"
            self.createDimension("STATION_DIMENSION", 5)

            # one variable for "timeseries" or "trajectory" which
            # has the station dimension and cf_role
            _var_name = feature_type.lower().split("profile")[0]
            self.createVariable(
                "{}_id_variable".format(_var_name),
                str,
                ("STATION_DIMENSION",),
                fill_value=None,
            )

            # set the cf_role
            self.variables["{}_id_variable".format(_var_name)].setncattr(
                "cf_role", "{}_id".format(_var_name)
            )

            # there will be one for the profile
            self.createVariable(
                "profile_id_variable", str, ("INSTANCE_DIMENSION",), fill_value=None
            )
            self.variables["profile_id_variable"].setncattr("cf_role", "profile_id")

            # will need a station index variable
            self.createVariable(
                "station_index_variable", int, ("INSTANCE_DIMENSION",), fill_value=None
            )

            self.variables["station_index_variable"].setncattr(
                "instance_dimension", "STATION_DIMENSION"
            )

            # also need counter variable, as compound featureTypes
            # are represented by having a contiguous repr for the
            # profiles and the indexed repr for the timeseries/trajectory
            # organization
            self.createVariable(
                "counter_var",
                "i",  # integer type
                ("INSTANCE_DIMENSION",),
                fill_value=None,
            )

            self.variables["counter_var"].setncattr(
                "sample_dimension", "SAMPLE_DIMENSION"
            )

        else:  # just a single featureType
            self.createVariable(
                "{}_id_variable".format(feature_type),
                str,
                ("INSTANCE_DIMENSION",),
                fill_value=None,
            )

            self.variables["{}_id_variable".format(feature_type)].setncattr(
                "cf_role", "{}_id".format(feature_type)
            )

            if structure == "contiguous":

                # create count variable
                self.createVariable(
                    "counter_var",
                    "i",  # integer type
                    ("INSTANCE_DIMENSION",),
                    fill_value=None,
                )

                self.variables["counter_var"].setncattr(
                    "sample_dimension", "SAMPLE_DIMENSION"
                )

            else:

                # create index variable
                self.createVariable(
                    "index_var",
                    "i",  # integer type
                    ("SAMPLE_DIMENSION",),
                    fill_value=None,
                )
                self.variables["index_var"].setncattr(
                    "instance_dimension", "INSTANCE_DIMENSION"
                )
