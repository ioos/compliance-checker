#!/usr/bin/env python
'''
compliance_checker/tests/test_protocols.py

Unit tests that ensure the compliance checker can successfully identify protocol endpoints
'''
from unittest import TestCase
from compliance_checker.suite import CheckSuite

import pytest


@pytest.mark.integration
class TestProtocols(TestCase):

    def test_erddap(self):
        '''
        Tests that a connection can be made to ERDDAP's GridDAP
        '''
        url = 'http://coastwatch.pfeg.noaa.gov/erddap/griddap/osuChlaAnom'
        cs = CheckSuite()
        ds = cs.load_dataset(url)
        assert ds is not None

    def test_hyrax(self):
        '''
        Tests that a connection can be made to Hyrax
        '''
        url = 'http://ingria.coas.oregonstate.edu/opendap/hyrax/aggregated/ocean_time_aggregation.ncml'
        cs = CheckSuite()
        ds = cs.load_dataset(url)
        assert ds is not None

    def test_thredds(self):
        '''
        Tests that a connection can be made to a remote THREDDS endpoint
        '''
        url = 'http://data.ioos.us/thredds/dodsC/deployments/rutgers/ru24-20150105T1441/ru24-20150105T1441.nc3.nc'
        cs = CheckSuite()
        ds = cs.load_dataset(url)
        assert ds is not None

    def test_sos(self):
        '''
        Tests that a connection can be made to an SOS endpoint
        '''
        url = 'http://data.oceansmap.com/thredds/sos/caricoos_ag/VIA/VIA.ncml'
        cs = CheckSuite()
        ds = cs.load_dataset(url)
        assert ds is not None
