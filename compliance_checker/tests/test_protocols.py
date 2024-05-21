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
    url = "https://www.neracoos.org/erddap/griddap/WW3_EastCoast_latest"
    cs = CheckSuite()
    ds = cs.load_dataset(url)
    assert ds is not None


def test_hyrax():
    """
    Tests that a connection can be made to Hyrax
    """
    # Returns: error 405
    # url = "http://test.opendap.org:8080/opendap/ioos/mday_joinExist.ncml"
    # More direct file
    url = "http://test.opendap.org:8080/opendap/ioos/mday_joinExist.ncml.dap.nc4"
    cs = CheckSuite()
    ds = cs.load_dataset(url)
    assert ds is not None


def test_thredds():
    """
    Tests that a connection can be made to a remote THREDDS endpoint
    """
    # Returns: error 400
    # url = "http://thredds.ucar.edu/thredds/dodsC/grib/NCEP/GFS/Global_0p25deg_ana/TP"
    # Use a smaller dataset
    url = "https://thredds.ucar.edu/thredds/ncss/grid/grib/NCEP/GFS/Global_0p25deg_ana/TP?var=Temperature_altitude_above_msl&accept=netcdf3"

    cs = CheckSuite()
    ds = cs.load_dataset(url)
    assert ds is not None


@pytest.mark.skip(reason="The thredds endpoint is no longer serving SOS.")
def test_sos():
    """
    Tests that a connection can be made to an SOS endpoint
    """
    url = "https://thredds.aoos.org/thredds/dodsC/aoos/cruises/ecofoci/2dy12.nc"
    cs = CheckSuite()
    ds = cs.load_dataset(url)
    assert ds is not None
