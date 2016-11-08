#!/usr/bin/env python
# -*- coding: utf-8 -*-

from compliance_checker.suite import CheckSuite
from netCDF4 import Dataset
from tempfile import gettempdir
from compliance_checker.tests.resources import STATIC_FILES
from compliance_checker.tests import BaseTestCase

import pytest
import os
import re


class TestCFIntegration(BaseTestCase):

    def setUp(self):
        '''
        Initialize the dataset
        '''
        self.cs = CheckSuite()
        self.cs.load_all_available_checkers()

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
        assert (scored, out_of) == (141, 149)
        assert len(messages) == 8
        assert u'standard_name temperature is not defined in Standard Name Table v36' in messages
        assert (u'auxiliary coordinate specified by the coordinates attribute, precise_lat, '
                'is not a variable in this dataset') in messages
        assert (u'auxiliary coordinate specified by the coordinates attribute, precise_lon, '
                'is not a variable in this dataset') in messages

        assert (u'attribute time:_CoordianteAxisType should begin with a letter and be composed '
                'of letters, digits, and underscores') in messages

    @pytest.mark.slowtest
    def test_ocos(self):
        dataset = self.load_dataset(STATIC_FILES['ocos'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (1768, 1829)
        assert len(messages) == 61
        assert (u'units attribute is required for user') in messages
        assert (u'Unidentifiable feature for variable Akt_bak') in messages
        assert (u"zeta's dimensions are not in the recommended order T, Z, Y, X. They are "
                "ocean_time, eta_rho, xi_rho") in messages
        assert (u'Conventions global attribute does not contain "CF-1.6"') in messages
        assert (u"CF recommends latitude variable 'lat_psi' to use units degrees_north") in messages

    def test_l01_met(self):
        dataset = self.load_dataset(STATIC_FILES['l01-met'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (599, 615)
        assert len(messages) == 16

        # The variable is supposed to be a status flag but it's mislabled
        assert (u'units for variable air_temperature_qc must be convertible to K currently they are 1') in messages
        assert (u'standard_name barometric_pressure is not defined in Standard Name Table v36') in messages
        assert (u'standard_name use_wind is not defined in Standard Name Table v36') in messages
        assert (u'standard_name modifier data_quality is not a valid modifier according to appendix C') in messages
        assert (u"CF recommends latitude variable 'lat' to use units degrees_north") in messages

    def test_usgs_dem_saipan(self):
        dataset = self.load_dataset(STATIC_FILES['usgs_dem_saipan'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (111, 112)
        assert len(messages) == 1
        assert (u'Conventions global attribute does not contain "CF-1.6"') == messages[0]

    def test_sp041(self):
        dataset = self.load_dataset(STATIC_FILES['sp041'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (1311, 1315)
        assert len(messages) == 4

        assert (u"lat_qc is not a variable in this dataset") in messages
        assert (u"TrajectoryProfile is not a valid CF featureType. It must be one of point, "
                "timeSeries, trajectory, profile, timeSeriesProfile, trajectoryProfile") in messages
        for i, msg in enumerate(messages):
            if msg.startswith("Different feature types"):
                break
        else:
            assert False, "'Different feature types discovered' was not found in the checker messages"

    def test_3mf07(self):
        dataset = self.load_dataset(STATIC_FILES['3mf07'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (427, 433)
        assert len(messages) == 6
        assert (u"dimensions for auxiliary coordinate variable z (z) are not a subset of dimensions for "
                "variable flag (profile)") in messages
        assert (u"Unidentifiable feature for variable flag") in messages
        assert (u"references global attribute should be a non-empty string") in messages
        assert (u"comment global attribute should be a non-empty string") in messages

    def test_ooi_glider(self):
        dataset = self.load_dataset(STATIC_FILES['ooi_glider'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (610, 615)
        assert len(messages) == 5
        assert (u"units attribute is required for deployment") in messages
        assert (u"variable deployment's attribute standard_name must be a non-empty string or "
                "it should define a long_name attribute.") in messages
        assert (u"comment global attribute should be a non-empty string") in messages
        assert (u"latitude variable 'latitude' should define standard_name='latitude' or axis='Y'") in messages
        assert (u"longitude variable 'longitude' should define standard_name='longitude' or axis='X'") in messages

    def test_swan(self):
        dataset = self.load_dataset(STATIC_FILES['swan'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (370, 380)
        assert len(messages) == 10
        assert (u"units for variable time_offset must be convertible to s currently they are hours "
                "since 2013-02-18T00:00:00Z") in messages
        assert (u"axis attribute must be T, X, Y, or Z, currently y") in messages
        assert (u"vertical coordinates not defining pressure must include a positive attribute that "
                "is either 'up' or 'down'") in messages
        assert (u"GRID is not a valid CF featureType. It must be one of point, timeSeries, "
                "trajectory, profile, timeSeriesProfile, trajectoryProfile") in messages
        assert (u'Conventions global attribute does not contain "CF-1.6"') in messages

    def test_kibesillah(self):
        dataset = self.load_dataset(STATIC_FILES['kibesillah'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (201, 204)
        assert len(messages) == 3
        assert (u"source should be defined") in messages
        assert (u"references should be defined") in messages
        assert (u"global attribute title should exist and be a non-empty string") in messages
