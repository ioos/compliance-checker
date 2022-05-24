#!/usr/bin/env python
"""
compliance_checker/tests/test_protocols.py

Unit tests that ensure the compliance checker can successfully identify protocol endpoints
"""
import pytest

from compliance_checker.suite import CheckSuite


pytestmark = [pytest.mark.integration]


@pytest.mark.vcr()
def test_netcdf_content_type():
    """
    Check that urls with Content-Type header of "application/x-netcdf" can
    successfully be read into memory for checks.
    """
    url = "https://gliders.ioos.us/erddap/tabledap/amelia-20180501T0000.ncCF?&time%3E=max(time)-1%20hour"
    cs = CheckSuite()
    ds = cs.load_dataset(url)
    assert ds is not None


@pytest.mark.vcr()
def test_erddap():
    """
    Tests that a connection can be made to ERDDAP's GridDAP
    """
    url = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/osuChlaAnom"
    cs = CheckSuite()
    ds = cs.load_dataset(url)
    assert ds is not None


def test_hyrax():
    """
    Tests that a connection can be made to Hyrax
    """
    url = "http://ingria.coas.oregonstate.edu/opendap/hyrax/aggregated/ocean_time_aggregation.ncml"
    cs = CheckSuite()
    ds = cs.load_dataset(url)
    assert ds is not None


def test_thredds():
    """
    Tests that a connection can be made to a remote THREDDS endpoint
    """
    url = "http://thredds.ucar.edu/thredds/dodsC/grib/NCEP/GFS/Global_0p25deg_ana/TP"

    cs = CheckSuite()
    ds = cs.load_dataset(url)
    assert ds is not None


def test_sos():
    """
    Tests that a connection can be made to an SOS endpoint
    """
    url = "https://data.oceansmap.com/thredds/sos/caricoos_ag/VIA/VIA.ncml"
    cs = CheckSuite()
    ds = cs.load_dataset(url)
    assert ds is not None
