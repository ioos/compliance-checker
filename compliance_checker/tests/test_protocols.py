#!/usr/bin/env python
"""
compliance_checker/tests/test_protocols.py

Unit tests that ensure the compliance checker can successfully identify protocol endpoints
"""
import platform

import pytest

from compliance_checker.protocols import zarr
from compliance_checker.suite import CheckSuite

from .conftest import datadir

pytestmark = [pytest.mark.integration]


@pytest.fixture
def cs():
    return CheckSuite()


@pytest.mark.vcr()
def test_netcdf_content_type(cs):
    """
    Check that urls with Content-Type header of "application/x-netcdf" can
    successfully be read into memory for checks.
    """
    url = "https://gliders.ioos.us/erddap/tabledap/amelia-20180501T0000.ncCF?&time%3E=max(time)-1%20hour"
    ds = cs.load_dataset(url)
    assert ds is not None


def test_erddap(cs):
    """
    Tests that a connection can be made to ERDDAP's GridDAP
    """
    url = "https://www.neracoos.org/erddap/griddap/WW3_EastCoast_latest"
    ds = cs.load_dataset(url)
    assert ds is not None


def test_hyrax(cs):
    """
    Tests that a connection can be made to Hyrax

    NOTE: This test creates a VCR file that is larger than GH's limit.
    We are not recording.
    """
    # Returns: error 405
    # url = "http://test.opendap.org:8080/opendap/ioos/mday_joinExist.ncml"
    # More direct file
    url = "http://test.opendap.org:8080/opendap/ioos/mday_joinExist.ncml.dap.nc4"
    ds = cs.load_dataset(url)
    assert ds is not None


def test_thredds(cs):
    """
    Tests that a connection can be made to a remote THREDDS endpoint

    NOTE: This test creates a VCR file that is larger than GH's limit.
    We are not recording.
    """
    # Returns: error 400
    # url = "http://thredds.ucar.edu/thredds/dodsC/grib/NCEP/GFS/Global_0p25deg_ana/TP"
    # Use a smaller dataset
    url = "https://thredds.ucar.edu/thredds/ncss/grid/grib/NCEP/GFS/Global_0p25deg_ana/TP?var=Temperature_altitude_above_msl&accept=netcdf3"

    cs = CheckSuite()
    ds = cs.load_dataset(url)
    assert ds is not None


str_dir = str(datadir.resolve()).replace("\\", "/")
file_url = "file://" + str_dir + "/trajectory.zarr#mode=nczarr,file"
s3_url = "s3://hrrrzarr/sfc/20210408/20210408_10z_anl.zarr#mode=nczarr,s3"
zip_url = "file://" + str_dir + "/zip.zarr#mode=nczarr,zip"
# replace slashes for windows compatibility
url_io = [
    ("s3://hrrrzarr/sfc/20210408/20210408_10z_anl.zarr", s3_url),
    (s3_url, s3_url),
    (datadir / "trajectory.zarr", file_url),
    ("file://" + str_dir + "/trajectory.zarr", file_url),
    (file_url, file_url),
    (datadir / "zip.zarr", zip_url),
    ("file://" + str_dir + "/zip.zarr", zip_url),
    (zip_url, zip_url),
]


@pytest.mark.skipif(
    platform.system() != "Linux",
    reason="NCZarr is not officially supported for your OS as of when this API was written",
)
@pytest.mark.parametrize("url_in, url_out", url_io)
def test_as_zarr(url_in, url_out):
    """
    Test is as_zurl can transform pointers to zarr datasets to valid nczarr urls,
    """
    assert zarr.as_zarr(url_in) == url_out
