#!/usr/bin/env python

import pytest

from compliance_checker.cf import util

# get current std names table version (it changes)
std_names = util.StandardNameTable()

# The meat and potatoes of most integration tests lies here
# To add a new integration test, simply add a new .CDL dataset somewhere in tests/data,
# and append this list with another tuple in the form:
# (file stem to test dataset in tests/data or subfolder : str, [expected messages]: list of strs )
# Note that this will also assert there was a missed compliance check with your dataset:
# assert scored < out_of
dataset_stem__expected_messages = [
    (
        "sldmb_43093_agg",
        [
            "attribute time:_CoordianteAxisType should begin with a letter and be composed of letters, digits, and underscores",
            "attribute lat:_CoordianteAxisType should begin with a letter and be composed of letters, digits, and underscores",
            "attribute lon:_CoordianteAxisType should begin with a letter and be composed of letters, digits, and underscores",
            "§2.6.2 global attribute history should exist and be a non-empty string",
            f"standard_name temperature is not defined in Standard Name Table v{std_names._version}. Possible close match(es): ['air_temperature', 'soil_temperature', 'snow_temperature']",
            "temperature's auxiliary coordinate specified by the coordinates attribute, precise_lat, is not a variable in this dataset",
            "temperature's auxiliary coordinate specified by the coordinates attribute, precise_lon, is not a variable in this dataset",
        ],
    ),
    (
        "usgs_dem_saipan",
        ['§2.6.1 Conventions global attribute does not contain "CF-1.9"'],
    ),
    (
        "l01-met",
        [
            "Attribute 'valid_range' (type: <class 'numpy.int16'>) and parent variable 'air_temperature_qc' (type: <class 'numpy.int8'>) must have equivalent datatypes",
            "Attribute 'valid_range' (type: <class 'numpy.int16'>) and parent variable 'barometric_pressure_qc' (type: <class 'numpy.int8'>) must have equivalent datatypes",
            "Attribute 'valid_range' (type: <class 'numpy.float64'>) and parent variable 'wind_gust' (type: <class 'numpy.float32'>) must have equivalent datatypes",
            "Attribute 'valid_range' (type: <class 'numpy.int16'>) and parent variable 'wind_gust_qc' (type: <class 'numpy.int8'>) must have equivalent datatypes",
            "Attribute 'valid_range' (type: <class 'numpy.float64'>) and parent variable 'wind_speed' (type: <class 'numpy.float32'>) must have equivalent datatypes",
            "Attribute 'valid_range' (type: <class 'numpy.int16'>) and parent variable 'wind_speed_qc' (type: <class 'numpy.int8'>) must have equivalent datatypes",
            "Attribute 'valid_range' (type: <class 'numpy.float64'>) and parent variable 'wind_direction' (type: <class 'numpy.float32'>) must have equivalent datatypes",
            "Attribute 'valid_range' (type: <class 'numpy.int16'>) and parent variable 'wind_direction_qc' (type: <class 'numpy.int8'>) must have equivalent datatypes",
            "Attribute 'valid_range' (type: <class 'numpy.int16'>) and parent variable 'visibility_qc' (type: <class 'numpy.int8'>) must have equivalent datatypes",
            '§2.6.1 Conventions global attribute does not contain "CF-1.8"',
            f"standard_name visibility is not defined in Standard Name Table v{std_names._version}. Possible close match(es): ['visibility_in_air']",
            'Standard name modifier "data_quality" for variable visibility_qc is not a valid modifier according to CF Appendix C',
            f"standard_name wind_direction is not defined in Standard Name Table v{std_names._version}. Possible close match(es): ['wind_to_direction', 'wind_from_direction', 'wind_gust_from_direction']",
            'Standard name modifier "data_quality" for variable wind_direction_qc is not a valid modifier according to CF Appendix C',
            f"standard_name wind_gust is not defined in Standard Name Table v{std_names._version}. Possible close match(es): ['y_wind_gust', 'x_wind_gust', 'wind_speed_of_gust']",
            'Standard name modifier "data_quality" for variable wind_gust_qc is not a valid modifier according to CF Appendix C',
            'Standard name modifier "data_quality" for variable air_temperature_qc is not a valid modifier according to CF Appendix C',
            f"standard_name use_wind is not defined in Standard Name Table v{std_names._version}. Possible close match(es): ['y_wind', 'x_wind']",
            f"standard_name barometric_pressure is not defined in Standard Name Table v{std_names._version}. Possible close match(es): ['air_pressure', 'reference_pressure', 'barometric_altitude']",
            'Standard name modifier "data_quality" for variable barometric_pressure_qc is not a valid modifier according to CF Appendix C',
            'Standard name modifier "data_quality" for variable wind_speed_qc is not a valid modifier according to CF Appendix C',
            f"standard_name barometric_pressure is not defined in Standard Name Table v{std_names._version}. Possible close match(es): ['air_pressure', 'reference_pressure', 'barometric_altitude']",
            "CF recommends latitude variable 'lat' to use units degrees_north",
            "CF recommends longitude variable 'lon' to use units degrees_east",
        ],
    ),
    ("sp041", ["lat_qc is not a variable in this dataset"]),
    (
        "3mf07",
        [
            "latitude:valid_min must be a numeric type not a string",
            "latitude:valid_max must be a numeric type not a string",
            "longitude:valid_min must be a numeric type not a string",
            "longitude:valid_max must be a numeric type not a string",
            "§2.6.2 references global attribute should be a non-empty string",
            "§2.6.2 comment global attribute should be a non-empty string",
            "dimensions for auxiliary coordinate variable z (z) are not a subset of dimensions for variable flag (profile)",
            "dimensions for auxiliary coordinate variable z (z) are not a subset of dimensions for variable haul (profile)",
        ],
    ),
    (
        "ooi_glider",
        [
            "§2.6.2 comment global attribute should be a non-empty string",
            "Attribute long_name or/and standard_name is highly recommended for variable deployment",
            "latitude variable 'latitude' should define standard_name='latitude' or axis='Y'",
            "longitude variable 'longitude' should define standard_name='longitude' or axis='X'",
        ],
    ),
    (
        "swan",
        [
            "global attribute _CoordSysBuilder should begin with a letter and be composed of letters, digits, and underscores",
            '§2.6.1 Conventions global attribute does not contain "CF-1.9"',
            'Units "hours since 2013-02-18T00:00:00Z" for variable time_offset must be convertible to canonical units "s"',
            "lon's axis attribute must be T, X, Y, or Z, currently x",
            "lat's axis attribute must be T, X, Y, or Z, currently y",
            "z's axis attribute must be T, X, Y, or Z, currently z",
            "z: vertical coordinates not defining pressure must include a positive attribute that is either 'up' or 'down'",
            "GRID is not a valid CF featureType. It must be one of point, timeseries, trajectory, profile, timeseriesprofile, trajectoryprofile",
        ],
    ),
    (
        "kibesillah",
        ["§2.6.2 global attribute title should exist and be a non-empty string"],
    ),
    (
        "pr_inundation",
        [
            'Units "degrees2" for variable area must be convertible to canonical units "m2"',
            "waterlevel's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are time (T), m (A), n (A) (with U: other/unknown; L: unlimited).",
            "velocity_x's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are time (T), Layer (Z), m (A), n (A) (with U: other/unknown; L: unlimited).",
            "velocity_y's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are time (T), Layer (Z), m (A), n (A) (with U: other/unknown; L: unlimited).",
            "tau_x's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are time (T), m (A), n (A) (with U: other/unknown; L: unlimited).",
            "tau_y's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are time (T), m (A), n (A) (with U: other/unknown; L: unlimited).",
            "§2.6.2 grid_depth:comment should be a non-empty string",
            "§2.6.2 depth:comment should be a non-empty string",
            "§2.6.2 institution global attribute should be a non-empty string",
            "§2.6.2 comment global attribute should be a non-empty string",
            'Units "degrees2" for variable area must be convertible to canonical units "m2"',
            "k: vertical coordinates not defining pressure must include a positive attribute that is either 'up' or 'down'",
            "grid_longitude has no coordinate associated with a variable identified as true latitude/longitude; its coordinate variable should also share a subset of grid_longitude's dimensions",
            "grid_latitude has no coordinate associated with a variable identified as true latitude/longitude; its coordinate variable should also share a subset of grid_latitude's dimensions",
            "time_bounds might be a cell boundary variable but there are no variables that define it as a boundary using the `bounds` attribute.",
        ],
    ),
    (
        "ww3",
        [
            "§2.6.2 global attribute title should exist and be a non-empty string",
            "§2.6.2 global attribute history should exist and be a non-empty string",
            "§2.6.1 Conventions field is not present",
            "Attribute long_name or/and standard_name is highly recommended for variable time",
            "Attribute long_name or/and standard_name is highly recommended for variable lon",
            "Attribute long_name or/and standard_name is highly recommended for variable lat",
            "latitude variable 'lat' should define standard_name='latitude' or axis='Y'",
            "longitude variable 'lon' should define standard_name='longitude' or axis='X'",
        ],
    ),
    (
        "glcfs",
        [
            # TODO: referenced/relative time is treated like time units
            'Units "hours since 2016-01-01T12:00:00Z" for variable time_offset must be convertible to canonical units "s"',
            f"standard_name cloud_cover is not defined in Standard Name Table v{std_names._version}. Possible close match(es): ['land_cover', 'land_cover_lccs', 'cloud_albedo']",
            f"standard_name dew_point is not defined in Standard Name Table v{std_names._version}. Possible close match(es): ['dew_point_depression', 'dew_point_temperature']",
            (
                "GRID is not a valid CF featureType. It must be one of point, timeseries, "
                "trajectory, profile, timeseriesprofile, trajectoryprofile"
            ),
            (
                "global attribute _CoordSysBuilder should begin with a letter and "
                "be composed of letters, digits, and underscores"
            ),
            'units for cl, "fraction" are not recognized by UDUNITS',
        ],
    ),
    (
        "bad_cf_role",
        [
            "§2.6.2 global attribute title should exist and be a non-empty string",
            "§2.6.2 global attribute history should exist and be a non-empty string",
            "§2.6.1 Conventions field is not present",
            "§9.5 The only acceptable values of cf_role for Discrete Geometry CF data sets are timeseries_id, profile_id, and trajectory_id",
        ],
    ),
    pytest.param(
        "ocos",
        [
            "AKs's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), s_w (Z), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            "AKt's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), s_w (Z), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            "AKv's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), s_w (Z), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            "latent's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            "lwrad's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            "salt's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), s_rho (Z), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            "sensible's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            "shflux's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            "swrad's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            "temp's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), s_rho (Z), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            "tke's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), s_w (Z), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            "u's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), s_rho (Z), eta_u (A), xi_u (A) (with U: other/unknown; L: unlimited).",
            "ubar's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), eta_u (A), xi_u (A) (with U: other/unknown; L: unlimited).",
            "v's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), s_rho (Z), eta_v (A), xi_v (A) (with U: other/unknown; L: unlimited).",
            "vbar's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), eta_v (A), xi_v (A) (with U: other/unknown; L: unlimited).",
            "w's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), s_w (Z), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            "zeta's spatio-temporal dimensions are not in the recommended order T, Z, Y, X and/or further dimensions are not located left of T, Z, Y, X. The dimensions (and their guessed types) are ocean_time (T), eta_rho (A), xi_rho (A) (with U: other/unknown; L: unlimited).",
            '§2.6.1 Conventions global attribute does not contain "CF-1.8"',
            "CF recommends latitude variable 'lat_rho' to use units degrees_north",
            "CF recommends latitude variable 'lat_u' to use units degrees_north",
            "CF recommends latitude variable 'lat_v' to use units degrees_north",
            "CF recommends latitude variable 'lat_psi' to use units degrees_north",
            "CF recommends longitude variable 'lon_rho' to use units degrees_east",
            "CF recommends longitude variable 'lon_u' to use units degrees_east",
            "CF recommends longitude variable 'lon_v' to use units degrees_east",
            "CF recommends longitude variable 'lon_psi' to use units degrees_east",
            "§4.3.3 The standard_name of `s_rho` must map to the correct computed_standard_name, `['altitude', 'height_above_geopotential_datum', 'height_above_mean_sea_level', 'height_above_reference_ellipsoid']`",
            "§4.3.3 The standard_name of `s_w` must map to the correct computed_standard_name, `['altitude', 'height_above_geopotential_datum', 'height_above_mean_sea_level', 'height_above_reference_ellipsoid']`",
        ],
        marks=pytest.mark.slowtest,
    ),
]

mult_msgs_diff = "Failed to find the following messages:\n{missing_msgs}\n\n\
        These were the messages captured:\n{found_msgs}\n\
            Please check wording and section names if messages have been altered since this test was written"


class TestCFIntegration:
    # --------------------------------------------------------------------------------
    # Helper Methods
    # --------------------------------------------------------------------------------

    def get_results(self, check_results, checksuite):
        """
        Returns a tuple of the value scored, possible, and a list of messages
        in the result set.\n
        cs: instance of CheckSuite object
        """
        aggregation = checksuite.build_structure(
            "cf",
            check_results["cf"][0],
            "test",
            1,
        )
        out_of = 0
        scored = 0
        results = aggregation["all_priorities"]
        for r in results:
            if isinstance(r.value, tuple):
                out_of += r.value[1]
                scored += r.value[0]
            else:
                out_of += 1
                scored += int(r.value)

        # Store the messages
        messages = []
        for r in results:
            messages.extend(r.msgs)

        return scored, out_of, messages

    # --------------------------------------------------------------------------------

    @pytest.mark.parametrize(
        "loaded_dataset,expected_messages",
        dataset_stem__expected_messages,
        indirect=[
            "loaded_dataset",
        ],  # must be specified to load this param at runtime, instead of at collection
    )
    def test_cf_integration(self, loaded_dataset, expected_messages, cs):
        # check_results = cs.run(loaded_dataset, [], "cf")
        check_results = cs.run_all(loaded_dataset, ["cf"], skip_checks=[])
        scored, out_of, messages = self.get_results(check_results, cs)

        assert scored < out_of

        assert all(m in messages for m in expected_messages), mult_msgs_diff.format(
            missing_msgs="\n".join([m for m in expected_messages if m not in messages]),
            found_msgs="\n".join(messages),
        )

    @pytest.mark.parametrize(
        "loaded_dataset,wrong_message",
        # From Github issue #845:\n
        # CF: incorrect errors for ragged array structure\n
        # \n
        # Structure needs to be correctly identified\n
        # http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html#_indexed_ragged_array_representation\n
        # https://github.com/ioos/compliance-checker/issues/845
        [
            ("index_ragged2", "are not a subset of dimensions for variable"),
            ("index_ragged2", "Unidentifiable feature for variable"),
        ],
        indirect=["loaded_dataset"],
    )
    def test_no_incorrect_errors(self, cs, loaded_dataset, wrong_message):
        # check_results = cs.run(loaded_dataset, [], True, "cf")
        check_results = cs.run_all(loaded_dataset, ["cf"], skip_checks=[])
        messages = self.get_results(check_results, cs)[-1]

        assert wrong_message not in "".join(messages)

    @pytest.mark.parametrize("loaded_dataset", ["fvcom"], indirect=True)
    def test_fvcom(self, cs, loaded_dataset):
        # check_results = cs.run(loaded_dataset, [], True, "cf")
        check_results = cs.run_all(loaded_dataset, ["cf"], skip_checks=[])
        scored, out_of, messages = self.get_results(check_results, cs)
        assert scored < out_of

        for msg in messages:
            if msg.startswith("dimensions for auxiliary coordinate variable siglay"):
                break
        # it's not clear to me what this is supposed to be doing -- this else clause is outside of the if
        else:
            raise AssertionError(
                '"dimensions for auxiliary coordinate variable siglay (node, siglay) '
                'are not a subset of dimensions for variable u (siglay, nele, time)"'
                " not in messages",
            )
        assert (
            '§2.6.1 Conventions global attribute does not contain "CF-1.8"'
        ) in messages

    @pytest.mark.parametrize(
        "loaded_dataset",
        ["NCEI_profile_template_v2.0_2016-09-22_181835.151325"],
        indirect=True,
    )
    def test_ncei_templates(self, cs, loaded_dataset):
        """
        Tests some of the NCEI NetCDF templates, which usually should get a
        perfect score.
        """
        # check_results = cs.run(loaded_dataset, [], "cf")
        check_results = cs.run_all(loaded_dataset, ["cf"], skip_checks=[])
        scored, out_of, messages = self.get_results(check_results, cs)
        assert scored < out_of
