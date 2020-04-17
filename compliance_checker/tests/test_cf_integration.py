#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from tempfile import gettempdir

import pytest

from netCDF4 import Dataset

from compliance_checker.cf import util
from compliance_checker.suite import CheckSuite
from compliance_checker.tests import BaseTestCase
from compliance_checker.tests.resources import STATIC_FILES


class TestCFIntegration(BaseTestCase):
    def setUp(self):
        """
        Initialize the dataset
        """
        self.cs = CheckSuite()
        self.cs.load_all_available_checkers()
        # get current std names table version (it changes)
        self._std_names = util.StandardNameTable()

    # --------------------------------------------------------------------------------
    # Helper Methods
    # --------------------------------------------------------------------------------

    def new_nc_file(self):
        """
        Make a new temporary netCDF file for the scope of the test
        """
        nc_file_path = os.path.join(gettempdir(), "example.nc")
        if os.path.exists(nc_file_path):
            raise IOError("File Exists: %s" % nc_file_path)
        nc = Dataset(nc_file_path, "w")
        self.addCleanup(os.remove, nc_file_path)
        self.addCleanup(nc.close)
        return nc

    def load_dataset(self, nc_dataset):
        """
        Return a loaded NC Dataset for the given path
        """
        if not isinstance(nc_dataset, str):
            raise ValueError("nc_dataset should be a string")

        nc_dataset = Dataset(nc_dataset, "r")
        self.addCleanup(nc_dataset.close)
        return nc_dataset

    def get_results(self, check_results):
        """
        Returns a tuple of the value scored, possible, and a list of messages
        in the result set.
        """
        aggregation = self.cs.build_structure("cf", check_results["cf"][0], "test", 1)
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

    def test_sldmb_43093_agg(self):
        dataset = self.load_dataset(STATIC_FILES["sldmb_43093_agg"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 7

        msgs = [
            u"attribute time:_CoordianteAxisType should begin with a letter and be composed of letters, digits, and underscores",
            u"attribute lat:_CoordianteAxisType should begin with a letter and be composed of letters, digits, and underscores",
            u"attribute lon:_CoordianteAxisType should begin with a letter and be composed of letters, digits, and underscores",
            u"§2.6.2 global attribute history should exist and be a non-empty string",
            u"standard_name temperature is not defined in Standard Name Table v49",
            u"temperature's auxiliary coordinate specified by the coordinates attribute, precise_lat, is not a variable in this dataset",
            u"temperature's auxiliary coordinate specified by the coordinates attribute, precise_lon, is not a variable in this dataset",
        ]
        assert all(m in messages for m in msgs)

    @pytest.mark.slowtest
    def test_ocos(self):
        dataset = self.load_dataset(STATIC_FILES["ocos"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert len(messages) == 63

        msgs = [
            u"zeta's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_rho, xi_rho",
            u"ubar's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_u, xi_u",
            u"vbar's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_v, xi_v",
            u"u's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_rho, eta_u, xi_u",
            u"v's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_rho, eta_v, xi_v",
            u"w's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_w, eta_rho, xi_rho",
            u"temp's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_rho, eta_rho, xi_rho",
            u"salt's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_rho, eta_rho, xi_rho",
            u"AKv's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_w, eta_rho, xi_rho",
            u"AKt's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_w, eta_rho, xi_rho",
            u"AKs's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_w, eta_rho, xi_rho",
            u"tke's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_w, eta_rho, xi_rho",
            u"shflux's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_rho, xi_rho",
            u"latent's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_rho, xi_rho",
            u"sensible's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_rho, xi_rho",
            u"lwrad's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_rho, xi_rho",
            u"swrad's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_rho, xi_rho",
            u'§2.6.1 Conventions global attribute does not contain "CF-1.6". The CF Checker only supports CF-1.6 at this time.',
            u"units (None) attribute of 's_w' must be a string compatible with UDUNITS",
            u"units (None) attribute of 's_rho' must be a string compatible with UDUNITS",
            u"units (None) attribute of 'Cs_w' must be a string compatible with UDUNITS",
            u"units (None) attribute of 'user' must be a string compatible with UDUNITS",
            u"units (None) attribute of 'Cs_r' must be a string compatible with UDUNITS",
            u"CF recommends latitude variable 'lat_rho' to use units degrees_north",
            u"CF recommends latitude variable 'lat_u' to use units degrees_north",
            u"CF recommends latitude variable 'lat_v' to use units degrees_north",
            u"CF recommends latitude variable 'lat_psi' to use units degrees_north",
            u"CF recommends longitude variable 'lon_rho' to use units degrees_east",
            u"CF recommends longitude variable 'lon_u' to use units degrees_east",
            u"CF recommends longitude variable 'lon_v' to use units degrees_east",
            u"CF recommends longitude variable 'lon_psi' to use units degrees_east",
            u"Unidentifiable feature for variable dt",
            u"Unidentifiable feature for variable dtfast",
            u"Unidentifiable feature for variable dstart",
            u"Unidentifiable feature for variable nl_tnu2",
            u"Unidentifiable feature for variable nl_visc2",
            u"Unidentifiable feature for variable Akt_bak",
            u"Unidentifiable feature for variable Akv_bak",
            u"Unidentifiable feature for variable Akk_bak",
            u"Unidentifiable feature for variable Akp_bak",
            u"Unidentifiable feature for variable rdrg",
            u"Unidentifiable feature for variable Zob",
            u"Unidentifiable feature for variable Zos",
            u"Unidentifiable feature for variable Znudg",
            u"Unidentifiable feature for variable M2nudg",
            u"Unidentifiable feature for variable M3nudg",
            u"Unidentifiable feature for variable Tnudg",
            u"Unidentifiable feature for variable FSobc_in",
            u"Unidentifiable feature for variable FSobc_out",
            u"Unidentifiable feature for variable M2obc_in",
            u"Unidentifiable feature for variable M2obc_out",
            u"Unidentifiable feature for variable Tobc_in",
            u"Unidentifiable feature for variable Tobc_out",
            u"Unidentifiable feature for variable M3obc_in",
            u"Unidentifiable feature for variable M3obc_out",
            u"Unidentifiable feature for variable rho0",
            u"Unidentifiable feature for variable xl",
            u"Unidentifiable feature for variable el",
            u"Unidentifiable feature for variable Tcline",
            u"Unidentifiable feature for variable hc",
            u"Unidentifiable feature for variable Cs_r",
            u"Unidentifiable feature for variable Cs_w",
            u"Unidentifiable feature for variable user",
        ]
        assert all([m in messages for m in msgs])

    def test_l01_met(self):
        dataset = self.load_dataset(STATIC_FILES["l01-met"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 16

        # The variable is supposed to be a status flag but it's mislabled
        msgs = [
            u"units for variable air_temperature_qc must be convertible to K currently they are 1",
            u"units for variable wind_speed_qc must be convertible to m s-1 currently they are 1",
            u"standard_name visibility is not defined in Standard Name Table v49",
            u"standard_name modifier data_quality for variable visibility_qc is not a valid modifier according to appendix C",
            u"standard_name wind_direction is not defined in Standard Name Table v49",
            u"standard_name modifier data_quality for variable wind_direction_qc is not a valid modifier according to appendix C",
            u"standard_name wind_gust is not defined in Standard Name Table v49",
            u"standard_name modifier data_quality for variable wind_gust_qc is not a valid modifier according to appendix C",
            u"standard_name modifier data_quality for variable air_temperature_qc is not a valid modifier according to appendix C",
            u"standard_name use_wind is not defined in Standard Name Table v49",
            u"standard_name barometric_pressure is not defined in Standard Name Table v49",
            u"standard_name modifier data_quality for variable barometric_pressure_qc is not a valid modifier according to appendix C",
            u"standard_name modifier data_quality for variable wind_speed_qc is not a valid modifier according to appendix C",
            u"standard_name barometric_pressure is not defined in Standard Name Table v49",
            u"CF recommends latitude variable 'lat' to use units degrees_north",
            u"CF recommends longitude variable 'lon' to use units degrees_east",
        ]

        assert all(m in messages for m in msgs)

    def test_usgs_dem_saipan(self):
        dataset = self.load_dataset(STATIC_FILES["usgs_dem_saipan"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 1

        msgs = [
            u'§2.6.1 Conventions global attribute does not contain "CF-1.6". The CF Checker only supports CF-1.6 at this time.'
        ]

        assert all(m in messages for m in msgs)

    def test_sp041(self):
        dataset = self.load_dataset(STATIC_FILES["sp041"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 3
        assert (u"lat_qc is not a variable in this dataset") in messages
        for i, msg in enumerate(messages):
            if msg.startswith("Different feature types"):
                break
        else:
            assert (
                False
            ), "'Different feature types discovered' was not found in the checker messages"

    def test_3mf07(self):
        """Load the 3mf07.nc file and run the CF check suite on it. There should be
        several variable/attribute combos which fail:
          - latitude:valid min
          - latitude:valid_max
          - longitude:valid_min
          - longitude:valid_max
          - references is an empty string
          - comment (global attr) is an empty string
          - z:dimensions are not a proper subset of dims for variable flag, haul
          - variable flag/haul has an unidentifiable feature"""

        dataset = self.load_dataset(STATIC_FILES["3mf07"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        msgs = [
            u"latitude:valid_min must be a numeric type not a string",
            u"latitude:valid_max must be a numeric type not a string",
            u"longitude:valid_min must be a numeric type not a string",
            u"longitude:valid_max must be a numeric type not a string",
            u"§2.6.2 references global attribute should be a non-empty string",
            u"§2.6.2 comment global attribute should be a non-empty string",
            u"dimensions for auxiliary coordinate variable z (z) are not a subset of dimensions for variable flag (profile)",
            u"dimensions for auxiliary coordinate variable z (z) are not a subset of dimensions for variable haul (profile)",
            u"Unidentifiable feature for variable flag",
            u"Unidentifiable feature for variable haul",
        ]

        assert scored < out_of
        assert all(m in messages for m in msgs)

    def test_ooi_glider(self):
        dataset = self.load_dataset(STATIC_FILES["ooi_glider"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 5

        msgs = [
            u"§2.6.2 comment global attribute should be a non-empty string",
            u"units (None) attribute of 'deployment' must be a string compatible with UDUNITS",
            u"Attribute long_name or/and standard_name is highly recommended for variable deployment",
            u"latitude variable 'latitude' should define standard_name='latitude' or axis='Y'",
            u"longitude variable 'longitude' should define standard_name='longitude' or axis='X'",
        ]

        assert all(m in messages for m in msgs)

    def test_swan(self):
        dataset = self.load_dataset(STATIC_FILES["swan"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 10

        msgs = [
            u"global attribute _CoordSysBuilder should begin with a letter and be composed of letters, digits, and underscores",
            u'§2.6.1 Conventions global attribute does not contain "CF-1.6". The CF Checker only supports CF-1.6 at this time.',
            u"units for variable time_offset must be convertible to s currently they are hours since 2013-02-18T00:00:00Z",
            u"units for variable time_run must be convertible to s currently they are hours since 2013-02-18 00:00:00.000 UTC",
            u"lon's axis attribute must be T, X, Y, or Z, currently x",
            "lat's axis attribute must be T, X, Y, or Z, currently y",
            u"z's axis attribute must be T, X, Y, or Z, currently z",
            u"z: vertical coordinates not defining pressure must include a positive attribute that is either 'up' or 'down'",
            u"GRID is not a valid CF featureType. It must be one of point, timeseries, trajectory, profile, timeseriesprofile, trajectoryprofile",
            u"Unidentifiable feature for variable time_offset",
        ]

        assert all(m in messages for m in msgs)

    def test_kibesillah(self):
        dataset = self.load_dataset(STATIC_FILES["kibesillah"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 1
        # test for global attributes (CF 2.6.2)
        assert (
            u"§2.6.2 global attribute title should exist and be a non-empty string"
        ) in messages

    def test_pr_inundation(self):
        dataset = self.load_dataset(STATIC_FILES["pr_inundation"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 21

        msgs = [
            u"waterlevel's dimensions are not in the recommended order T, Z, Y, X. They are time, m, n",
            u"velocity_x's dimensions are not in the recommended order T, Z, Y, X. They are time, Layer, m, n",
            u"velocity_y's dimensions are not in the recommended order T, Z, Y, X. They are time, Layer, m, n",
            u"tau_x's dimensions are not in the recommended order T, Z, Y, X. They are time, m, n",
            u"tau_y's dimensions are not in the recommended order T, Z, Y, X. They are time, m, n",
            u"§2.6.2 grid_depth:comment should be a non-empty string",
            u"§2.6.2 depth:comment should be a non-empty string",
            u"§2.6.2 institution global attribute should be a non-empty string",
            u"§2.6.2 comment global attribute should be a non-empty string",
            u"units (None) attribute of 'LayerInterf' must be a string compatible with UDUNITS",
            u"units (None) attribute of 'time_bounds' must be a string compatible with UDUNITS",
            u"units (None) attribute of 'Layer' must be a string compatible with UDUNITS",
            u"units for variable area must be convertible to m2 currently they are degrees2",
            u"k: vertical coordinates not defining pressure must include a positive attribute that is either 'up' or 'down'",
            u"grid_longitude is not associated with a coordinate defining true latitude and sharing a subset of dimensions",
            u"grid_longitude is not associated with a coordinate defining true longitude and sharing a subset of dimensions",
            u"grid_latitude is not associated with a coordinate defining true latitude and sharing a subset of dimensions",
            u"grid_latitude is not associated with a coordinate defining true longitude and sharing a subset of dimensions",
            u"time_bounds might be a cell boundary variable but there are no variables that define it as a boundary using the `bounds` attribute.",
            u"Unidentifiable feature for variable time_bounds",
            u"Unidentifiable feature for variable grid_depth",
        ]
        assert all(m in messages for m in msgs)

    def test_fvcom(self):
        dataset = self.load_dataset(STATIC_FILES["fvcom"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 40

        for msg in messages:
            if msg.startswith("dimensions for auxiliary coordinate variable siglay"):
                break
        # it's not clear to me what this is supposed to be doing -- this else clause is outside of the if
        else:
            raise AssertionError(
                u'"dimensions for auxiliary coordinate variable siglay (node, siglay) '
                'are not a subset of dimensions for variable u (siglay, nele, time)"'
                " not in messages"
            )
        assert (u"Unidentifiable feature for variable x") in messages
        assert (
            u"§2.6.1 Conventions global attribute does not contain "
            '"CF-1.6". The CF Checker only supports CF-1.6 '
            "at this time."
        ) in messages

    def test_ww3(self):
        dataset = self.load_dataset(STATIC_FILES["ww3"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 8

        msgs = [
            u"§2.6.2 global attribute title should exist and be a non-empty string",
            u"§2.6.2 global attribute history should exist and be a non-empty string",
            u"§2.6.1 Conventions field is not present",
            u"Attribute long_name or/and standard_name is highly recommended for variable time",
            u"Attribute long_name or/and standard_name is highly recommended for variable lon",
            u"Attribute long_name or/and standard_name is highly recommended for variable lat",
            u"latitude variable 'lat' should define standard_name='latitude' or axis='Y'",
            u"longitude variable 'lon' should define standard_name='longitude' or axis='X'",
        ]

        assert all(m in messages for m in msgs)

    def test_glcfs(self):
        dataset = self.load_dataset(STATIC_FILES["glcfs"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 14
        assert (
            u"units for variable time_offset must be convertible to s currently "
            "they are hours since 2016-01-01T12:00:00Z"
        ) in messages
        assert (
            u"standard_name cloud_cover is not defined in Standard Name Table v{}".format(
                self._std_names._version
            )
        ) in messages
        assert (
            u"standard_name dew_point is not defined in Standard Name Table v{}".format(
                self._std_names._version
            )
        ) in messages
        assert (
            u"GRID is not a valid CF featureType. It must be one of point, timeseries, "
            "trajectory, profile, timeseriesprofile, trajectoryprofile"
        ) in messages
        assert (
            u"global attribute _CoordSysBuilder should begin with a letter and "
            "be composed of letters, digits, and underscores"
        ) in messages
        assert u"source should be defined"
        assert (u'units for cl, "fraction" are not recognized by UDUNITS') in messages

    def test_ncei_templates(self):
        """
        Tests some of the NCEI NetCDF templates, which usually should get a
        perfect score.
        """
        dataset = self.load_dataset(STATIC_FILES["NCEI_profile_template_v2_0"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of

    def test_bad_cf_roles(self):
        """
        Tests the CF checker detects datasets with more than 2 defined cf_role
        variables.
        """
        dataset = self.load_dataset(STATIC_FILES["bad_cf_role"])
        check_results = self.cs.run(dataset, [], "cf")
        scored, out_of, messages = self.get_results(check_results)

        msgs = [
            u"§2.6.2 global attribute title should exist and be a non-empty string",
            u"§2.6.2 global attribute history should exist and be a non-empty string",
            u"§2.6.1 Conventions field is not present",
            u"Unidentifiable feature for variable T",
            u"§9.5 The only acceptable values of cf_role for Discrete Geometry CF data sets are timeseries_id, profile_id, and trajectory_id",
        ]

        assert scored < out_of
        assert all(m in messages for m in msgs)
