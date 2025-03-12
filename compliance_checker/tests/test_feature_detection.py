#!/usr/bin/env python
"""
compliance_checker/tests/test_feature_detection.py
"""

from unittest import TestCase

from netCDF4 import Dataset

from compliance_checker import cfutil as util
from compliance_checker.tests import resources
from compliance_checker.tests.helpers import MockRaggedArrayRepr


class TestFeatureDetection(TestCase):
    """
    Tests the feature type detection of cdftools
    """

    def test_point(self):
        """
        Ensures point detection works
        """
        with Dataset(resources.STATIC_FILES["point"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_point(nc, variable), f"{variable} is point"

    def test_timeseries(self):
        """
        Ensures timeseries detection works
        """
        with Dataset(resources.STATIC_FILES["timeseries"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_timeseries(nc, variable), f"{variable} is timeseries"

    def test_multi_timeseries_orthogonal(self):
        """
        Ensures multi-timeseries-orthogonal detection works
        """
        with Dataset(resources.STATIC_FILES["multi-timeseries-orthogonal"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_multi_timeseries_orthogonal(
                    nc,
                    variable,
                ), f"{variable} is multi-timeseries orthogonal"

    def test_multi_timeseries_incomplete(self):
        """
        Ensures multi-timeseries-incomplete detection works
        """
        with Dataset(resources.STATIC_FILES["multi-timeseries-incomplete"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_multi_timeseries_incomplete(
                    nc,
                    variable,
                ), f"{variable} is multi-timeseries incomplete"

    def test_trajectory(self):
        """
        Ensures trajectory detection works
        """
        with Dataset(resources.STATIC_FILES["trajectory"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_cf_trajectory(nc, variable), f"{variable} is trajectory"

    def test_trajectory_single(self):
        """
        Ensures trajectory-single detection works
        """
        with Dataset(resources.STATIC_FILES["trajectory-single"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_single_trajectory(
                    nc,
                    variable,
                ), f"{variable} is trajectory-single"

    def test_profile_orthogonal(self):
        """
        Ensures profile-orthogonal detection works
        """
        with Dataset(resources.STATIC_FILES["profile-orthogonal"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_profile_orthogonal(
                    nc,
                    variable,
                ), f"{variable} is profile-orthogonal"

    def test_profile_incomplete(self):
        """
        Ensures profile-incomplete detection works
        """
        with Dataset(resources.STATIC_FILES["profile-incomplete"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_profile_incomplete(
                    nc,
                    variable,
                ), f"{variable} is profile-incomplete"

    def test_timeseries_profile_single_station(self):
        """
        Ensures timeseries profile single station detection works
        """
        with Dataset(resources.STATIC_FILES["timeseries-profile-single-station"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_timeseries_profile_single_station(
                    nc,
                    variable,
                ), f"{variable} is timeseries-profile-single-station"

    def test_timeseries_profile_multi_station(self):
        """
        Ensures timeseries profile multi station detection works
        """
        with Dataset(resources.STATIC_FILES["timeseries-profile-multi-station"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_timeseries_profile_multi_station(
                    nc,
                    variable,
                ), f"{variable} is timeseries-profile-multi-station"

    def test_timeseries_profile_single_ortho_time(self):
        """
        Ensures timeseries profile single station ortho time detection works
        """
        with Dataset(
            resources.STATIC_FILES["timeseries-profile-single-ortho-time"],
        ) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_timeseries_profile_single_ortho_time(
                    nc,
                    variable,
                ), f"{variable} is timeseries-profile-single-ortho-time"

    def test_timeseries_profile_multi_ortho_time(self):
        """
        Ensures timeseries profile multi station ortho time detection works
        """
        with Dataset(
            resources.STATIC_FILES["timeseries-profile-multi-ortho-time"],
        ) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_timeseries_profile_multi_ortho_time(
                    nc,
                    variable,
                ), f"{variable} is timeseries-profile-multi-ortho-time"

    def test_timeseries_profile_ortho_depth(self):
        """
        Ensures timeseries profile ortho depth detection works
        """
        with Dataset(resources.STATIC_FILES["timeseries-profile-ortho-depth"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_timeseries_profile_ortho_depth(
                    nc,
                    variable,
                ), f"{variable} is timeseries-profile-ortho-depth"

    def test_timeseries_profile_incomplete(self):
        """
        Ensures timeseries profile station incomplete detection works
        """
        with Dataset(resources.STATIC_FILES["timeseries-profile-incomplete"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_timeseries_profile_incomplete(
                    nc,
                    variable,
                ), f"{variable} is timeseries-profile-incomplete"

    def test_trajectory_profile_orthogonal(self):
        """
        Ensures trajectory profile orthogonal detection works
        """
        with Dataset(resources.STATIC_FILES["trajectory-profile-orthogonal"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_trajectory_profile_orthogonal(
                    nc,
                    variable,
                ), f"{variable} is trajectory profile orthogonal"

    def test_trajectory_profile_incomplete(self):
        """
        Ensures trajectory profile incomplete detection works
        """
        with Dataset(resources.STATIC_FILES["trajectory-profile-incomplete"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_trajectory_profile_incomplete(
                    nc,
                    variable,
                ), f"{variable} is trajectory profile incomplete"

    def test_2d_regular_grid(self):
        """
        Ensures 2D Regular Grid detection works
        """
        with Dataset(resources.STATIC_FILES["2d-regular-grid"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_2d_regular_grid(
                    nc,
                    variable,
                ), f"{variable} is 2D regular grid"

    def test_2d_static_grid(self):
        """
        Ensures 2D Static Grid detection works
        """
        with Dataset(resources.STATIC_FILES["2d-static-grid"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_2d_static_grid(
                    nc,
                    variable,
                ), f"{variable} is a 2D static grid"

    def test_3d_regular_grid(self):
        """
        Ensures 2U Regular Grid detection works
        """
        with Dataset(resources.STATIC_FILES["3d-regular-grid"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_3d_regular_grid(
                    nc,
                    variable,
                ), f"{variable} is 3d regular grid"

    def test_3d_static_grid(self):
        """
        Ensures 3D Static Grid detection works
        """
        with Dataset(resources.STATIC_FILES["3d-static-grid"]) as nc:
            for variable in util.get_geophysical_variables(nc):
                assert util.is_3d_static_grid(
                    nc,
                    variable,
                ), f"{variable} is a 3D static grid"

    def test_boundaries(self):
        """
        Ensures that boundary variables are not listed as geophysical variables
        """
        with Dataset(resources.STATIC_FILES["grid-boundaries"]) as nc:
            assert "lat_bnds" not in util.get_geophysical_variables(nc)
            assert "lon_bnds" not in util.get_geophysical_variables(nc)
            assert "lat_bnds" in util.get_cell_boundary_variables(nc)
            assert "lon_bnds" in util.get_cell_boundary_variables(nc)

            boundary_map = util.get_cell_boundary_map(nc)
            assert boundary_map["lat"] == "lat_bnds"
            assert boundary_map["lon"] == "lon_bnds"

    def test_climatology(self):
        """
        Ensures that climatology variables are identified as climatology variables and not geophysical variables
        """
        with Dataset(resources.STATIC_FILES["climatology"]) as nc:
            geophysical_variables = util.get_geophysical_variables(nc)
            climatology_variable = util.get_climatology_variable(nc)
            assert "temperature" in geophysical_variables
            assert "climatology_bounds" not in geophysical_variables
            assert "climatology_bounds" == climatology_variable

    def test_grid_mapping(self):
        """
        Ensures that grid mapping variables are properly identified
        """
        with Dataset(resources.STATIC_FILES["rotated_pole_grid"]) as nc:
            grid_mapping = util.get_grid_mapping_variables(nc)
            coordinate_variables = util.get_coordinate_variables(nc)
            axis_variables = util.get_axis_variables(nc)

            assert "rotated_pole" in grid_mapping
            assert {"rlon", "rlat", "lev"} == set(coordinate_variables)
            assert {"rlon", "rlat", "lev"} == set(axis_variables)
            assert "lat" == util.get_lat_variable(nc)
            assert "lon" == util.get_lon_variable(nc)

    def test_auxiliary_coordinates(self):
        """
        Ensures variables are classified as auxiliary coordinate variables
        """
        with Dataset(resources.STATIC_FILES["bad_units"]) as nc:
            coordinate_variables = util.get_coordinate_variables(nc)
            assert {"time"} == set(coordinate_variables)

            aux_coord_vards = util.get_auxiliary_coordinate_variables(nc)
            assert {"lat", "lon"} == set(aux_coord_vards)

    def test_forecast_reference_metadata(self):
        """
        Tests variables used for forecast reference metadata to ensure they are
        not misclassified as geophysical variables.
        """
        with Dataset(resources.STATIC_FILES["forecast_reference"]) as nc:
            self.assertFalse(util.is_geophysical(nc, "forecast_reference_time"))
            self.assertFalse(util.is_geophysical(nc, "forecast_hour"))
            self.assertTrue(util.is_geophysical(nc, "air_temp"))
            self.assertFalse(util.is_geophysical(nc, "time"))

            assert len(util.get_coordinate_variables(nc)) == 3
            assert len(util.get_geophysical_variables(nc)) == 1

    def test_rotated_pole_grid(self):
        with Dataset(resources.STATIC_FILES["rotated_pole_grid"]) as nc:
            latitudes = util.get_latitude_variables(nc)
            assert sorted(latitudes) == ["lat", "rlat"]
            assert util.is_mapped_grid(nc, "temperature") is True

    def test_vertical_coords(self):
        with Dataset(resources.STATIC_FILES["vertical_coords"]) as nc:
            vertical = util.get_z_variables(nc)
            assert vertical == ["height"]

    def test_reduced_grid(self):
        with Dataset(resources.STATIC_FILES["reduced_horizontal_grid"]) as nc:
            assert util.guess_feature_type(nc, "PS") == "reduced-grid"

    def test_global_feature_detection(self):
        with Dataset(resources.STATIC_FILES["reduced_horizontal_grid"]) as nc:
            assert util.guess_feature_type(nc, "PS") == "reduced-grid"

        with Dataset(resources.STATIC_FILES["vertical_coords"]) as nc:
            assert util.guess_feature_type(nc, "temperature") == "point"

            axis_map = util.get_axis_map(nc, "temperature")
            assert axis_map["Z"] == ["height"]
            assert axis_map["T"] == ["time"]

        with Dataset(resources.STATIC_FILES["2d-regular-grid"]) as nc:
            assert util.guess_feature_type(nc, "temperature") == "2d-regular-grid"

            axis_map = util.get_axis_map(nc, "temperature")
            assert axis_map["T"] == ["time"]
            assert axis_map["Z"] == ["z"]
            assert axis_map["X"] == ["lon"]
            assert axis_map["Y"] == ["lat"]

        with Dataset(resources.STATIC_FILES["2dim"]) as nc:
            assert util.guess_feature_type(nc, "T") == "mapped-grid"

            axis_map = util.get_axis_map(nc, "T")
            assert axis_map["Z"] == ["lev"]
            assert axis_map["Y"] == ["yc", "lat"]
            assert axis_map["X"] == ["xc", "lon"]

        with Dataset(resources.STATIC_FILES["3d-regular-grid"]) as nc:
            assert util.guess_feature_type(nc, "temperature") == "3d-regular-grid"

            axis_map = util.get_axis_map(nc, "temperature")
            assert axis_map["T"] == ["time"]
            assert axis_map["Z"] == ["z"]
            assert axis_map["Y"] == ["lat"]
            assert axis_map["X"] == ["lon"]

        with Dataset(resources.STATIC_FILES["climatology"]) as nc:
            assert util.guess_feature_type(nc, "temperature") == "timeseries"

            axis_map = util.get_axis_map(nc, "temperature")
            assert axis_map["T"] == ["time"]
            assert axis_map["Z"] == []
            assert axis_map["Y"] == []
            assert axis_map["X"] == []

        with Dataset(resources.STATIC_FILES["index_ragged"]) as nc:
            assert util.guess_feature_type(nc, "temperature") == "trajectory"

            axis_map = util.get_axis_map(nc, "temperature")
            assert axis_map["T"] == ["time"]
            assert axis_map["Z"] == ["z"]
            assert axis_map["Y"] == ["lat"]
            assert axis_map["X"] == ["lon"]

        with Dataset(resources.STATIC_FILES["mapping"]) as nc:
            assert util.guess_feature_type(nc, "sea_surface_height") == "timeseries"

            axis_map = util.get_axis_map(nc, "sea_surface_height")
            assert axis_map["T"] == ["time"]
            assert axis_map["Z"] == []
            assert axis_map["Y"] == ["lat"]
            assert axis_map["X"] == ["lon"]

        with Dataset(resources.STATIC_FILES["rotated_pole_grid"]) as nc:
            assert util.guess_feature_type(nc, "temperature") == "mapped-grid"

            axis_map = util.get_axis_map(nc, "temperature")
            assert axis_map["T"] == []
            assert axis_map["Z"] == ["lev"]
            assert axis_map["Y"] == ["rlat", "lat"]
            assert axis_map["X"] == ["rlon", "lon"]

        with Dataset(resources.STATIC_FILES["rutgers"]) as nc:
            assert util.guess_feature_type(nc, "temperature") == "trajectory"

            axis_map = util.get_axis_map(nc, "temperature")
            assert axis_map["T"] == ["time"]
            assert axis_map["Z"] == ["depth"]
            assert axis_map["Y"] == ["lat"]
            assert axis_map["X"] == ["lon"]

        with Dataset(resources.STATIC_FILES["self-referencing-var"]) as nc:
            assert util.guess_feature_type(nc, "TEMP") == "point"

            axis_map = util.get_axis_map(nc, "TEMP")
            assert axis_map["T"] == ["TIME"]
            assert axis_map["Z"] == ["DEPTH"]
            assert axis_map["Y"] == []
            assert axis_map["X"] == []

        with Dataset(resources.STATIC_FILES["2d-static-grid"]) as nc:
            assert util.guess_feature_type(nc, "T") == "2d-static-grid"

            axis_map = util.get_axis_map(nc, "T")
            assert axis_map["X"] == ["lon"]
            assert axis_map["Y"] == ["lat"]
            assert axis_map["T"] == []
            assert axis_map["Z"] == []

        with Dataset(resources.STATIC_FILES["3d-static-grid"]) as nc:
            assert util.guess_feature_type(nc, "T") == "3d-static-grid"

            axis_map = util.get_axis_map(nc, "T")
            assert axis_map["X"] == ["lon"]
            assert axis_map["Y"] == ["lat"]
            assert axis_map["T"] == []
            assert axis_map["Z"] == ["depth"]

    def test_is_variable_valid_ragged_array_repr_featureType(self):
        nc = MockRaggedArrayRepr("timeseries", "indexed")

        # add a variable that isn't recognized as geophysical
        v = nc.createVariable("data1", "d", ("SAMPLE_DIMENSION",), fill_value=None)
        v.setncattr("cf_role", "blah")
        self.assertFalse(
            util.is_variable_valid_ragged_array_repr_featureType(nc, "data1"),
        )

        # add geophysical variable with correct dimension
        nc = MockRaggedArrayRepr("timeseries", "indexed")
        v = nc.createVariable("data1", "d", ("SAMPLE_DIMENSION",), fill_value=None)
        v.setncattr("standard_name", "sea_water_pressure")
        # test the variable
        self.assertTrue(
            util.is_variable_valid_ragged_array_repr_featureType(nc, "data1"),
        )

        # add good variable and another variable, this time with the improper dimension
        nc = MockRaggedArrayRepr("timeseries", "indexed")
        v = nc.createVariable("data1", "d", ("SAMPLE_DIMENSION",), fill_value=None)
        v.setncattr("standard_name", "sea_water_pressure")
        v2 = nc.createVariable("data2", "d", ("INSTANCE_DIMENSION",), fill_value=None)
        v2.setncattr("standard_name", "sea_water_salinity")

        # good variable should pass, second should fail
        self.assertTrue(
            util.is_variable_valid_ragged_array_repr_featureType(nc, "data1"),
        )
        self.assertFalse(
            util.is_variable_valid_ragged_array_repr_featureType(nc, "data2"),
        )

    def test_is_dataset_valid_ragged_array_repr_featureType(self):
        # first test single featureType

        # ----- timeseries, indexed ----- #

        nc = MockRaggedArrayRepr("timeseries", "indexed")
        self.assertTrue(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "timeseries"),
        )

        # we'll add another cf_role variable
        nc = MockRaggedArrayRepr("timeseries", "indexed")
        v = nc.createVariable("var2", "i", ("INSTANCE_DIMENSION",), fill_value=None)
        v.setncattr("cf_role", "yeetyeet_id")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "timeseries"),
        )

        # we'll add another index variable, also bad
        nc = MockRaggedArrayRepr("timeseries", "indexed")
        v = nc.createVariable("index_var2", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        v.setncattr("instance_dimension", "INSTANCE_DIMENSION")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "timeseries"),
        )

        # ----- timeseries, contiguous ----- #
        nc = MockRaggedArrayRepr("timeseries", "contiguous")
        self.assertTrue(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "timeseries"),
        )

        # add another cf_role var, bad
        nc = MockRaggedArrayRepr("timeseries", "contiguous")
        v = nc.createVariable("var2", "i", ("INSTANCE_DIMENSION",), fill_value=None)
        v.setncattr("cf_role", "yeetyeet_id")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "timeseries"),
        )

        # add another count variable, bad
        v = nc.createVariable(
            "count_var2",
            "i",
            ("INSTANCE_DIMENSION",),
            fill_value=None,
        )
        v.setncattr("sample_dimension", "SAMPLE_DIMENSION")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "timeseries"),
        )

        # ----- profile, indexed ----- #

        nc = MockRaggedArrayRepr("profile", "indexed")
        self.assertTrue(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "profile"),
        )

        # add another cf_role var
        nc = MockRaggedArrayRepr("profile", "indexed")
        v = nc.createVariable("var2", "i", ("INSTANCE_DIMENSION",), fill_value=None)
        v.setncattr("cf_role", "yeetyeet_id")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "profile"),
        )

        # we'll add another index variable, also bad
        nc = MockRaggedArrayRepr("profile", "indexed")
        v = nc.createVariable("index_var2", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        v.setncattr("instance_dimension", "INSTANCE_DIMENSION")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "profile"),
        )

        # ----- profile, contiguous ----- #
        nc = MockRaggedArrayRepr("profile", "contiguous")
        self.assertTrue(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "profile"),
        )

        # add another cf_role var
        nc = MockRaggedArrayRepr("profile", "contiguous")
        v = nc.createVariable("var2", "i", ("INSTANCE_DIMENSION",), fill_value=None)
        v.setncattr("cf_role", "yeetyeet_id")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "profile"),
        )

        # we'll add another count variable, also bad
        nc = MockRaggedArrayRepr("profile", "contiguous")
        v = nc.createVariable(
            "index_var2",
            "i",
            ("INSTANCE_DIMENSION",),
            fill_value=None,
        )
        v.setncattr("sample_dimension", "SAMPLE_DIMENSION")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "profile"),
        )

        # ----- trajectory, indexed ----- #
        nc = MockRaggedArrayRepr("trajectory", "indexed")
        self.assertTrue(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "trajectory"),
        )

        # add another cf_role var
        nc = MockRaggedArrayRepr("trajectory", "indexed")
        v = nc.createVariable("var2", "i", ("INSTANCE_DIMENSION",), fill_value=None)
        v.setncattr("cf_role", "yeetyeet_id")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "trajectory"),
        )

        # we'll add another index variable, also bad
        nc = MockRaggedArrayRepr("trajectory", "indexed")
        v = nc.createVariable("index_var2", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        v.setncattr("instance_dimension", "INSTANCE_DIMENSION")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "trajectory"),
        )

        # ----- trajectory, contiguous ----- #
        nc = MockRaggedArrayRepr("trajectory", "contiguous")
        self.assertTrue(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "trajectory"),
        )

        # add another cf_role var
        nc = MockRaggedArrayRepr("trajectory", "contiguous")
        v = nc.createVariable("var2", "i", ("INSTANCE_DIMENSION",), fill_value=None)
        v.setncattr("cf_role", "yeetyeet_id")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "trajectory"),
        )

        # we'll add another count variable, also bad
        nc = MockRaggedArrayRepr("trajectory", "contiguous")
        v = nc.createVariable(
            "index_var2",
            "i",
            ("INSTANCE_DIMENSION",),
            fill_value=None,
        )
        v.setncattr("sample_dimension", "SAMPLE_DIMENSION")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(nc, "trajectory"),
        )

        # ----- now test compound featureType ----- #

        # ----- timeSeriesProfile ----- #
        nc = MockRaggedArrayRepr("timeSeriesProfile")

        # NOTE
        # has no geophysical vars, so should (?) (will) fail
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(
                nc,
                "timeseriesprofile",
            ),
        )

        # add a geophysical variable and test again
        nc = MockRaggedArrayRepr("timeSeriesProfile")
        v1 = nc.createVariable("data1", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        v1.setncattr("standard_name", "pressure")
        self.assertTrue(
            util.is_dataset_valid_ragged_array_repr_featureType(
                nc,
                "timeseriesprofile",
            ),
        )

        nc = MockRaggedArrayRepr("timeSeriesProfile")
        v1 = nc.createVariable("data1", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        # add a third cf_role variable - this should fail
        v = nc.createVariable(
            "cf_role_var3",
            "i",
            ("INSTANCE_DIMENSION",),
            fill_value=None,
        )
        v.setncattr("cf_role", "yeetyeet_id")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(
                nc,
                "timeseriesprofile",
            ),
        )

        # set the index variable to have an incorrect attr
        nc = MockRaggedArrayRepr("timeSeriesProfile")
        v1 = nc.createVariable("data1", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        nc.variables["station_index_variable"].instance_dimension = "SIKE!"

        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(
                nc,
                "timeseriesprofile",
            ),
        )

        # change the sample_dimension attr on the count variable, bad
        nc = MockRaggedArrayRepr("timeSeriesProfile")
        v1 = nc.createVariable("data1", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        nc.variables["counter_var"].sample_dimension = "SIKE!"

        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(
                nc,
                "timeseriesprofile",
            ),
        )

        # give another geophysical data variable a different dimension
        nc = MockRaggedArrayRepr("timeSeriesProfile")
        v1 = nc.createVariable("data1", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        v1 = nc.createVariable(
            "data2",
            "i",
            ("STATION_DIMENSION",),
            fill_value=None,  # bad!
        )
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(
                nc,
                "timeseriesprofile",
            ),
        )

        # ----- trajectoryProfile ----- #
        nc = MockRaggedArrayRepr("trajectoryProfile")

        # NOTE
        # has no geophysical vars, so should (?) (will) fail
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(
                nc,
                "trajectoryprofile",
            ),
        )

        # add a geophysical variable and test again
        nc = MockRaggedArrayRepr("trajectoryProfile")
        v1 = nc.createVariable("data1", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        v1.setncattr("standard_name", "pressure")
        self.assertTrue(
            util.is_dataset_valid_ragged_array_repr_featureType(
                nc,
                "trajectoryprofile",
            ),
        )

        nc = MockRaggedArrayRepr("trajectoryProfile")
        v1 = nc.createVariable("data1", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        # add a third cf_role variable - this should fail
        v = nc.createVariable(
            "cf_role_var3",
            "i",
            ("INSTANCE_DIMENSION",),
            fill_value=None,
        )
        v.setncattr("cf_role", "yeetyeet_id")
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(
                nc,
                "trajectoryprofile",
            ),
        )

        # set the index variable to have an incorrect attr
        nc = MockRaggedArrayRepr("trajectoryProfile")
        v1 = nc.createVariable("data1", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        nc.variables["station_index_variable"].instance_dimension = "SIKE!"

        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(
                nc,
                "trajectoryprofile",
            ),
        )

        # change the sample_dimension attr on the count variable, bad
        nc = MockRaggedArrayRepr("trajectoryProfile")
        v1 = nc.createVariable("data1", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        nc.variables["counter_var"].sample_dimension = "SIKE!"

        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(
                nc,
                "trajectoryprofile",
            ),
        )

        # give another geophysical data variable a different dimension
        nc = MockRaggedArrayRepr("trajectoryProfile")
        v1 = nc.createVariable("data1", "i", ("SAMPLE_DIMENSION",), fill_value=None)
        v1 = nc.createVariable(
            "data2",
            "i",
            ("STATION_DIMENSION",),
            fill_value=None,  # bad!
        )
        self.assertFalse(
            util.is_dataset_valid_ragged_array_repr_featureType(
                nc,
                "trajectoryprofile",
            ),
        )
