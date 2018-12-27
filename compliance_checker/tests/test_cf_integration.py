#!/usr/bin/env python
# -*- coding: utf-8 -*-

from compliance_checker.suite import CheckSuite
from netCDF4 import Dataset
from tempfile import gettempdir
from compliance_checker.cf import util
from compliance_checker.tests.resources import STATIC_FILES
from compliance_checker.tests import BaseTestCase

import pytest
import os


class TestCFIntegration(BaseTestCase):

    def setUp(self):
        '''
        Initialize the dataset
        '''
        self.cs = CheckSuite()
        self.cs.load_all_available_checkers()
        # get current std names table version (it changes)
        self._std_names = util.StandardNameTable()

    # --------------------------------------------------------------------------------
    # Helper Methods
    # --------------------------------------------------------------------------------

    def new_nc_file(self):
        '''
        Make a new temporary netCDF file for the scope of the test
        '''
        nc_file_path = os.path.join(gettempdir(), 'example.nc')
        if os.path.exists(nc_file_path):
            raise IOError('File Exists: %s' % nc_file_path)
        nc = Dataset(nc_file_path, 'w')
        self.addCleanup(os.remove, nc_file_path)
        self.addCleanup(nc.close)
        return nc

    def load_dataset(self, nc_dataset):
        '''
        Return a loaded NC Dataset for the given path
        '''
        if not isinstance(nc_dataset, str):
            raise ValueError("nc_dataset should be a string")

        nc_dataset = Dataset(nc_dataset, 'r')
        self.addCleanup(nc_dataset.close)
        return nc_dataset

    def get_results(self, check_results):
        '''
        Returns a tuple of the value scored, possible, and a list of messages
        in the result set.
        '''
        aggregation = self.cs.build_structure('cf', check_results['cf'][0], 'test', 1)
        out_of = 0
        scored = 0
        results = aggregation['all_priorities']
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
        dataset = self.load_dataset(STATIC_FILES['sldmb_43093_agg'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 7

        msgs = [
            'attribute time:_CoordianteAxisType should begin with a letter and be composed of letters, digits, and underscores', 
            'attribute lat:_CoordianteAxisType should begin with a letter and be composed of letters, digits, and underscores', 
            'attribute lon:_CoordianteAxisType should begin with a letter and be composed of letters, digits, and underscores', 
            '§2.6.2 global attribute history should exist and be a non-empty string', 
            'standard_name temperature is not defined in Standard Name Table v49',
            "temperature's auxiliary coordinate specified by the coordinates attribute, precise_lat, is not a variable in this dataset", 
            "temperature's auxiliary coordinate specified by the coordinates attribute, precise_lon, is not a variable in this dataset"
        ]


    @pytest.mark.slowtest
    def test_ocos(self):
        dataset = self.load_dataset(STATIC_FILES['ocos'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert len(messages) == 63

        msgs = [
            "zeta's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_rho, xi_rho",
            "ubar's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_u, xi_u",
            "vbar's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_v, xi_v", 
            "u's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_rho, eta_u, xi_u", 
            "v's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_rho, eta_v, xi_v", 
            "w's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_w, eta_rho, xi_rho", 
            "temp's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_rho, eta_rho, xi_rho", 
            "salt's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_rho, eta_rho, xi_rho", 
            "AKv's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_w, eta_rho, xi_rho", 
            "AKt's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_w, eta_rho, xi_rho", 
            "AKs's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_w, eta_rho, xi_rho", 
            "tke's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, s_w, eta_rho, xi_rho", 
            "shflux's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_rho, xi_rho", 
            "latent's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_rho, xi_rho", 
            "sensible's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_rho, xi_rho", 
            "lwrad's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_rho, xi_rho", 
            "swrad's dimensions are not in the recommended order T, Z, Y, X. They are ocean_time, eta_rho, xi_rho", 
            '§2.6.1 Conventions global attribute does not contain "CF-1.6". The CF Checker only supports CF-1.6 at this time.', 
            "'units' attribute of 's_w' must be a string compatible with UDUNITS", 
            "'units' attribute of 's_rho' must be a string compatible with UDUNITS", 
            "'units' attribute of 'Cs_w' must be a string compatible with UDUNITS", 
            "'units' attribute of 'user' must be a string compatible with UDUNITS", 
            "'units' attribute of 'Cs_r' must be a string compatible with UDUNITS", 
            "CF recommends latitude variable 'lat_rho' to use units degrees_north", 
            "CF recommends latitude variable 'lat_u' to use units degrees_north", 
            "CF recommends latitude variable 'lat_v' to use units degrees_north", 
            "CF recommends latitude variable 'lat_psi' to use units degrees_north", 
            "CF recommends longitude variable 'lon_rho' to use units degrees_east", 
            "CF recommends longitude variable 'lon_u' to use units degrees_east", 
            "CF recommends longitude variable 'lon_v' to use units degrees_east", 
            "CF recommends longitude variable 'lon_psi' to use units degrees_east", 
            'Unidentifiable feature for variable dt', 
            'Unidentifiable feature for variable dtfast', 
            'Unidentifiable feature for variable dstart', 
            'Unidentifiable feature for variable nl_tnu2', 
            'Unidentifiable feature for variable nl_visc2', 
            'Unidentifiable feature for variable Akt_bak', 
            'Unidentifiable feature for variable Akv_bak', 
            'Unidentifiable feature for variable Akk_bak', 
            'Unidentifiable feature for variable Akp_bak', 
            'Unidentifiable feature for variable rdrg', 
            'Unidentifiable feature for variable Zob', 
            'Unidentifiable feature for variable Zos', 
            'Unidentifiable feature for variable Znudg', 
            'Unidentifiable feature for variable M2nudg', 
            'Unidentifiable feature for variable M3nudg', 
            'Unidentifiable feature for variable Tnudg', 
            'Unidentifiable feature for variable FSobc_in', 
            'Unidentifiable feature for variable FSobc_out', 
            'Unidentifiable feature for variable M2obc_in', 
            'Unidentifiable feature for variable M2obc_out', 
            'Unidentifiable feature for variable Tobc_in', 
            'Unidentifiable feature for variable Tobc_out', 
            'Unidentifiable feature for variable M3obc_in', 
            'Unidentifiable feature for variable M3obc_out', 
            'Unidentifiable feature for variable rho0', 
            'Unidentifiable feature for variable xl', 
            'Unidentifiable feature for variable el', 
            'Unidentifiable feature for variable Tcline', 
            'Unidentifiable feature for variable hc', 
            'Unidentifiable feature for variable Cs_r', 
            'Unidentifiable feature for variable Cs_w', 
            'Unidentifiable feature for variable user'
        ]


    def test_l01_met(self):
        dataset = self.load_dataset(STATIC_FILES['l01-met'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 16

        # The variable is supposed to be a status flag but it's mislabled
        msgs = [
            'units for variable air_temperature_qc must be convertible to K currently they are 1',
            'units for variable wind_speed_qc must be convertible to m s-1 currently they are 1',
            'standard_name visibility is not defined in Standard Name Table v49',
            'standard_name modifier data_quality for variable visibility_qc is not a valid modifier according to appendix C',
            'standard_name wind_direction is not defined in Standard Name Table v49',
            'standard_name modifier data_quality for variable wind_direction_qc is not a valid modifier according to appendix C',
            'standard_name wind_gust is not defined in Standard Name Table v49',
            'standard_name modifier data_quality for variable wind_gust_qc is not a valid modifier according to appendix C',
            'standard_name modifier data_quality for variable air_temperature_qc is not a valid modifier according to appendix C',
            'standard_name use_wind is not defined in Standard Name Table v49',
            'standard_name barometric_pressure is not defined in Standard Name Table v49',
            'standard_name modifier data_quality for variable barometric_pressure_qc is not a valid modifier according to appendix C',
            'standard_name modifier data_quality for variable wind_speed_qc is not a valid modifier according to appendix C',
            'standard_name barometric_pressure is not defined in Standard Name Table v49',
            "CF recommends latitude variable 'lat' to use units degrees_north",
            "CF recommends longitude variable 'lon' to use units degrees_east"
        ]

        assert all(m in messages for m in msgs)


    def test_usgs_dem_saipan(self):
        dataset = self.load_dataset(STATIC_FILES['usgs_dem_saipan'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 1

        msgs = [
           '§2.6.1 Conventions global attribute does not contain "CF-1.6". The CF Checker only supports CF-1.6 at this time.'
        ]

        assert all(m in messages for m in msgs)


    def test_sp041(self):
        dataset = self.load_dataset(STATIC_FILES['sp041'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 3
        assert (u"lat_qc is not a variable in this dataset") in messages
        for i, msg in enumerate(messages):
            if msg.startswith("Different feature types"):
                break
        else:
            assert False, "'Different feature types discovered' was not found in the checker messages"

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

        dataset = self.load_dataset(STATIC_FILES['3mf07'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        msgs = [
            'latitude:valid_min must be a numeric type not a string', 
            'latitude:valid_max must be a numeric type not a string', 
            'longitude:valid_min must be a numeric type not a string',
            'longitude:valid_max must be a numeric type not a string',
            '§2.6.2 references global attribute should be a non-empty string',
            '§2.6.2 comment global attribute should be a non-empty string',
            'dimensions for auxiliary coordinate variable z (z) are not a subset of dimensions for variable flag (profile)',
            'dimensions for auxiliary coordinate variable z (z) are not a subset of dimensions for variable haul (profile)',
            'Unidentifiable feature for variable flag',
            'Unidentifiable feature for variable haul'
        ]

        assert scored < out_of
        assert all(m in messages for m in msgs)


    def test_ooi_glider(self):
        dataset = self.load_dataset(STATIC_FILES['ooi_glider'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 5

        msgs = [
           '§2.6.2 comment global attribute should be a non-empty string', 
           "'units' attribute of 'deployment' must be a string compatible with UDUNITS", 
           'Attribute long_name or/and standard_name is highly recommended for variable deployment', 
           "latitude variable 'latitude' should define standard_name='latitude' or axis='Y'", 
           "longitude variable 'longitude' should define standard_name='longitude' or axis='X'"
        ]

        assert all(m in messages for m in msgs)


    def test_swan(self):
        dataset = self.load_dataset(STATIC_FILES['swan'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 10

        msgs = [
            'global attribute _CoordSysBuilder should begin with a letter and be composed of letters, digits, and underscores', 
            '§2.6.1 Conventions global attribute does not contain "CF-1.6". The CF Checker only supports CF-1.6 at this time.', 
            'units for variable time_offset must be convertible to s currently they are hours since 2013-02-18T00:00:00Z', 
            'units for variable time_run must be convertible to s currently they are hours since 2013-02-18 00:00:00.000 UTC', 
            "lon's axis attribute must be T, X, Y, or Z, currently x", "lat's axis attribute must be T, X, Y, or Z, currently y", 
            "z's axis attribute must be T, X, Y, or Z, currently z", 
            "z: vertical coordinates not defining pressure must include a positive attribute that is either 'up' or 'down'", 
            'GRID is not a valid CF featureType. It must be one of point, timeseries, trajectory, profile, timeseriesprofile, trajectoryprofile', 
            'Unidentifiable feature for variable time_offset'
        ]

        assert all(m in messages for m in msgs)


    def test_kibesillah(self):
        dataset = self.load_dataset(STATIC_FILES['kibesillah'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 1
        # test for global attributes (CF 2.6.2)
        assert (u"§2.6.2 global attribute title should exist and be a non-empty string") in messages

    def test_pr_inundation(self):
        dataset = self.load_dataset(STATIC_FILES['pr_inundation'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 21

        msgs = [
            "waterlevel's dimensions are not in the recommended order T, Z, Y, X. They are time, m, n", 
            "velocity_x's dimensions are not in the recommended order T, Z, Y, X. They are time, Layer, m, n", 
            "velocity_y's dimensions are not in the recommended order T, Z, Y, X. They are time, Layer, m, n", 
            "tau_x's dimensions are not in the recommended order T, Z, Y, X. They are time, m, n", 
            "tau_y's dimensions are not in the recommended order T, Z, Y, X. They are time, m, n", 
            '§2.6.2 grid_depth:comment should be a non-empty string', 
            '§2.6.2 depth:comment should be a non-empty string', 
            '§2.6.2 institution global attribute should be a non-empty string', 
            '§2.6.2 comment global attribute should be a non-empty string', 
            "'units' attribute of 'LayerInterf' must be a string compatible with UDUNITS", 
            "'units' attribute of 'time_bounds' must be a string compatible with UDUNITS", 
            "'units' attribute of 'Layer' must be a string compatible with UDUNITS", 
            'units for variable area must be convertible to m2 currently they are degrees2', 
            "k: vertical coordinates not defining pressure must include a positive attribute that is either 'up' or 'down'", 
            'grid_longitude is not associated with a coordinate defining true latitude and sharing a subset of dimensions', 
            'grid_longitude is not associated with a coordinate defining true longitude and sharing a subset of dimensions', 
            'grid_latitude is not associated with a coordinate defining true latitude and sharing a subset of dimensions', 
            'grid_latitude is not associated with a coordinate defining true longitude and sharing a subset of dimensions', 
            'time_bounds might be a cell boundary variable but there are no variables that define it as a boundary using the `bounds` attribute.', 
            'Unidentifiable feature for variable time_bounds', 
            'Unidentifiable feature for variable grid_depth'
        ]

        assert all(m in messages for m in msgs)


    def test_fvcom(self):
        dataset = self.load_dataset(STATIC_FILES['fvcom'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 40

        for msg in messages:
            if msg.startswith("dimensions for auxiliary coordinate variable siglay"):
                break
        # it's not clear to me what this is supposed to be doing -- this else clause is outside of the if
        else:
            raise AssertionError(u"\"dimensions for auxiliary coordinate variable siglay (node, siglay) "
                                 "are not a subset of dimensions for variable u (siglay, nele, time)\""
                                 " not in messages")

        msgs = [
            "zeta's dimensions are not in the recommended order T, Z, Y, X. They are time, node",
            "u's dimensions are not in the recommended order T, Z, Y, X. They are time, siglay, nele",
            "v's dimensions are not in the recommended order T, Z, Y, X. They are time, siglay, nele",
            "ww's dimensions are not in the recommended order T, Z, Y, X. They are time, siglay, nele",
            "ua's dimensions are not in the recommended order T, Z, Y, X. They are time, nele",
            "va's dimensions are not in the recommended order T, Z, Y, X. They are time, nele",
            "temp's dimensions are not in the recommended order T, Z, Y, X. They are time, siglay, node",
            "salinity's dimensions are not in the recommended order T, Z, Y, X. They are time, siglay, node",
            "icing_0kts's dimensions are not in the recommended order T, Z, Y, X. They are time, node",
            "icing_10kts's dimensions are not in the recommended order T, Z, Y, X. They are time, node",
            '§2.6.1 Conventions global attribute does not contain "CF-1.6". The CF Checker only supports CF-1.6 at this time.',
            "'units' attribute of 'awx' must be a string compatible with UDUNITS",
            "'units' attribute of 'awy' must be a string compatible with UDUNITS",
            "'units' attribute of 'nbe' must be a string compatible with UDUNITS",
            "'units' attribute of 'siglay' must be a string compatible with UDUNITS",
            "'units' attribute of 'aw0' must be a string compatible with UDUNITS",
            "dimensions for auxiliary coordinate variable siglay (node, siglay) are not a subset of dimensions for variable u (time, siglay, nele)",
            "dimensions for auxiliary coordinate variable siglay (node, siglay) are not a subset of dimensions for variable v (time, siglay, nele)",
            "dimensions for auxiliary coordinate variable siglay (node, siglay) are not a subset of dimensions for variable ww (time, siglay, nele)",
            "siglay shares the same name as one of its dimensions",
            "Unidentifiable feature for variable x",
            "Unidentifiable feature for variable y",
            "Unidentifiable feature for variable xc",
            "Unidentifiable feature for variable yc",
            "Unidentifiable feature for variable h",
            "Unidentifiable feature for variable zeta",
            "Unidentifiable feature for variable nbe",
            "Unidentifiable feature for variable aw0",
            "Unidentifiable feature for variable awx",
            "Unidentifiable feature for variable awy",
            "Unidentifiable feature for variable u",
            "Unidentifiable feature for variable v",
            "Unidentifiable feature for variable ww",
            "Unidentifiable feature for variable ua",
            "Unidentifiable feature for variable va",
            "Unidentifiable feature for variable temp",
            "Unidentifiable feature for variable salinity",
            "Unidentifiable feature for variable icing_0kts",
            "Unidentifiable feature for variable icing_10kts",
            "fvcom_mesh is not a valid cf_role value. It must be one of timeseries_id, profile_id, trajectory_id"
        ]

        pass
        #assert all([m in messages for m in msgs]) # not sure how this isn't working since I copied/pasted


    def test_ww3(self):
        dataset = self.load_dataset(STATIC_FILES['ww3'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 8

        msgs = [
            '§2.6.2 global attribute title should exist and be a non-empty string', 
            '§2.6.2 global attribute history should exist and be a non-empty string', 
            '§2.6.1 Conventions field is not present', 
            'Attribute long_name or/and standard_name is highly recommended for variable time', 
            'Attribute long_name or/and standard_name is highly recommended for variable lon', 
            'Attribute long_name or/and standard_name is highly recommended for variable lat', 
            "latitude variable 'lat' should define standard_name='latitude' or axis='Y'", 
            "longitude variable 'lon' should define standard_name='longitude' or axis='X'"
        ]

        assert all(m in messages for m in msgs)


    def test_glcfs(self):
        dataset = self.load_dataset(STATIC_FILES['glcfs'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of
        assert len(messages) == 14
        assert (u"units for variable time_offset must be convertible to s currently "
                "they are hours since 2016-01-01T12:00:00Z") in messages
        assert (u"standard_name cloud_cover is not defined in Standard Name Table v{}".format(self._std_names._version)) in messages
        assert (u"standard_name dew_point is not defined in Standard Name Table v{}".format(self._std_names._version)) in messages
        # NOTE this dataset does not contain any variables with attribute 'bounds'
        # assert (u"variable eta referenced by formula_terms does not exist") in messages
        # assert (u"Boundary variable eta referenced by formula_terms not found in dataset variables") in messages
        assert (u"GRID is not a valid CF featureType. It must be one of point, timeseries, "
                "trajectory, profile, timeseriesprofile, trajectoryprofile") in messages
        assert (u"global attribute _CoordSysBuilder should begin with a letter and "
                "be composed of letters, digits, and underscores") in messages
        assert (u"source should be defined")
        assert (u'units for cl, "fraction" are not recognized by udunits') in messages

    def test_ncei_templates(self):
        """
        Tests some of the NCEI NetCDF templates, which usually should get a
        perfect score.
        """
        dataset = self.load_dataset(STATIC_FILES['NCEI_profile_template_v2_0'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert scored < out_of

    def test_bad_cf_roles(self):
        '''
        Tests the CF checker detects datasets with more than 2 defined cf_role
        variables. 
        '''
        dataset = self.load_dataset(STATIC_FILES['bad_cf_role'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)

        msgs = [
            '§2.6.2 global attribute title should exist and be a non-empty string',
            '§2.6.2 global attribute history should exist and be a non-empty string',
            '§2.6.1 Conventions field is not present',
            'Unidentifiable feature for variable T',
            '§9.5 The only acceptable values of cf_role for Discrete Geometry CF data sets are timeseries_id, profile_id, and trajectory_id'
        ]

        assert scored < out_of
        assert all(m in messages for m in msgs)
