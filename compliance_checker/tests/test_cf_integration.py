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
        assert (scored, out_of) == (137, 145)
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
        assert (scored, out_of) == (1840, 1841)
        assert len(messages) == 41
        assert (u'Unidentifiable feature for variable Akt_bak') in messages
        assert (u'Conventions global attribute does not contain "CF-1.6"') in messages
        assert (u"CF recommends latitude variable 'lat_psi' to use units degrees_north") in messages

    def test_l01_met(self):
        dataset = self.load_dataset(STATIC_FILES['l01-met'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (588, 602)
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
        assert (scored, out_of) == (109, 110)
        assert len(messages) == 1
        assert (u'Conventions global attribute does not contain "CF-1.6"') == messages[0]

    def test_sp041(self):
        dataset = self.load_dataset(STATIC_FILES['sp041'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (1188, 1192)
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
        assert (scored, out_of) == (418, 426)
        assert len(messages) == 10
        assert (u"dimensions for auxiliary coordinate variable z (z) are not a subset of dimensions for "
                "variable flag (profile)") in messages
        assert (u"Unidentifiable feature for variable flag") in messages
        assert (u"references global attribute should be a non-empty string") in messages
        assert (u"comment global attribute should be a non-empty string") in messages
        assert (u"latitude:valid_min must be a numeric type not a string") in messages
        assert (u"latitude:valid_max must be a numeric type not a string") in messages
        assert (u"longitude:valid_min must be a numeric type not a string") in messages
        assert (u"longitude:valid_max must be a numeric type not a string") in messages

    def test_ooi_glider(self):
        dataset = self.load_dataset(STATIC_FILES['ooi_glider'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (593, 597)
        assert len(messages) == 4
        assert (u"variable deployment's attribute standard_name must be a non-empty string or "
                "it should define a long_name attribute.") in messages
        assert (u"comment global attribute should be a non-empty string") in messages
        assert (u"latitude variable 'latitude' should define standard_name='latitude' or axis='Y'") in messages
        assert (u"longitude variable 'longitude' should define standard_name='longitude' or axis='X'") in messages

    def test_swan(self):
        dataset = self.load_dataset(STATIC_FILES['swan'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (363, 372)
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
        assert (scored, out_of) == (197, 200)
        assert len(messages) == 3
        assert (u"source should be defined") in messages
        assert (u"references should be defined") in messages
        assert (u"global attribute title should exist and be a non-empty string") in messages

    def test_pr_inundation(self):
        dataset = self.load_dataset(STATIC_FILES['pr_inundation'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (549, 557)
        assert len(messages) == 11
        assert (u"units for variable area must be convertible to m2 currently "
                "they are degrees2") in messages
        assert (u"vertical coordinates not defining pressure must include a positive "
                "attribute that is either 'up' or 'down'") in messages
        assert (u"Unidentifiable feature for variable time_bounds")in messages
        assert (u"grid_longitude's dimensions are not in the recommended order T, "
                "Z, Y, X. They are m, n, bounds4") in messages
        assert (u"depth:comment should be a non-empty string") in messages
        assert (u"institution global attribute should be a non-empty string") in messages
        assert (u"time_bounds might be a cell boundary variable but there are no variables that "
                "define it as a boundary using the `bounds` attribute.") in messages

    def test_fvcom(self):
        dataset = self.load_dataset(STATIC_FILES['fvcom'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (646, 651)
        assert len(messages) == 24

        for msg in messages:
            if msg.startswith("dimensions for auxiliary coordinate variable siglay"):
                break
        else:
            raise AssertionError(u"\"dimensions for auxiliary coordinate variable siglay (node, siglay) "
                                 "are not a subset of dimensions for variable u (siglay, nele, time)\""
                                 " not in messages")
        assert (u"Unidentifiable feature for variable x") in messages
        assert (u'Conventions global attribute does not contain "CF-1.6"') in messages
        assert (u"siglay shares the same name as one of its dimensions") in messages

    def test_ww3(self):
        dataset = self.load_dataset(STATIC_FILES['ww3'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (111, 121)
        assert len(messages) == 10

        assert (u"variable lat's attribute standard_name must be a non-empty string or it "
                "should define a long_name attribute.") in messages
        assert (u"Conventions field is not present") in messages
        assert (u"latitude variable 'lat' should define standard_name='latitude' or axis='Y'") in messages

    def test_glcfs(self):
        dataset = self.load_dataset(STATIC_FILES['glcfs'])
        check_results = self.cs.run(dataset, [], 'cf')
        scored, out_of, messages = self.get_results(check_results)
        assert (scored, out_of) == (330, 339)
        assert len(messages) == 10

        assert (u"units for variable time_offset must be convertible to s currently "
                "they are hours since 2016-01-01T12:00:00Z") in messages
        assert (u"standard_name cloud_cover is not defined in Standard Name Table v36") in messages
        assert (u"standard_name dew_point is not defined in Standard Name Table v36") in messages
        assert (u"variable eta referenced by formula_terms does not exist") in messages
        assert (u"GRID is not a valid CF featureType. It must be one of point, timeSeries, "
                "trajectory, profile, timeSeriesProfile, trajectoryProfile") in messages
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
        assert (scored, out_of) == (342, 346)
