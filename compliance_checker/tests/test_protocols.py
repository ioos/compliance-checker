#!/usr/bin/env python
"""
compliance_checker/tests/test_protocols.py

Unit tests that ensure the compliance checker can successfully identify protocol endpoints
"""
from unittest import TestCase

import pytest

from compliance_checker.suite import CheckSuite


@pytest.mark.integration
class TestProtocols(TestCase):
    def test_netcdf_content_type(self):
        """
        Check that urls with Content-Type header of "application/x-netcdf" can
        successfully be read into memory for checks.
        """
        url = "https://gliders.ioos.us/erddap/tabledap/amelia-20180501T0000.ncCF?&time%3E=max(time)-1%20hour"
        cs = CheckSuite()
        ds = cs.load_dataset(url)
        assert ds is not None

    def test_erddap(self):
        """
        Tests that a connection can be made to ERDDAP's GridDAP
        """
        url = "http://coastwatch.pfeg.noaa.gov/erddap/griddap/osuChlaAnom"
        cs = CheckSuite()
        ds = cs.load_dataset(url)
        assert ds is not None

    def test_hyrax(self):
        """
        Tests that a connection can be made to Hyrax
        """
        url = "http://test.opendap.org:8080/opendap/ioos/mday_joinExist.ncml"
        cs = CheckSuite()
        ds = cs.load_dataset(url)
        assert ds is not None

    def test_thredds(self):
        """
        Tests that a connection can be made to a remote THREDDS endpoint
        """
        url = (
            "http://thredds.ucar.edu/thredds/dodsC/grib/NCEP/GFS/Global_0p25deg_ana/TP"
        )

        cs = CheckSuite()
        ds = cs.load_dataset(url)
        assert ds is not None

    def test_sos(self):
        """
        Tests that a connection can be made to an SOS endpoint
        """
        url = "https://thredds.aoos.org/thredds/sos/aoos/cruises/ecofoci/2dy12.nc"
        cs = CheckSuite()
        ds = cs.load_dataset(url)
        assert ds is not None
