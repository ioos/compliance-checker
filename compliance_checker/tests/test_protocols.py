#!/usr/bin/env python
"""
compliance_checker/tests/test_protocols.py

Unit tests that ensure the compliance checker can successfully identify protocol endpoints
"""
import pytest

from compliance_checker.protocols import zarr
from compliance_checker.suite import CheckSuite

from .conftest import datadir


id_url = {
    # Check that urls with Content-Type header of "application/x-netcdf" can
    # successfully be read into memory for checks.
    "netcdf_content_type": "https://gliders.ioos.us/erddap/tabledap/amelia-20180501T0000.ncCF?&time%3E=max(time)-1%20hour",
    # Tests that a connection can be made to ERDDAP's GridDAP
    "erddap": "http://coastwatch.pfeg.noaa.gov/erddap/griddap/osuChlaAnom",
    # Tests that a connection can be made to Hyrax
    "hyrax": "http://ingria.coas.oregonstate.edu/opendap/hyrax/aggregated/ocean_time_aggregation.ncml",
    # Tests that a connection can be made to a remote THREDDS endpoint
    "thredds": "http://thredds.ucar.edu/thredds/dodsC/grib/NCEP/GFS/Global_0p25deg_ana/TP",
    # Tests that a connection can be made to an SOS endpoint
    "sos": "https://data.oceansmap.com/thredds/sos/caricoos_ag/VIA/VIA.ncml",
}


class TestProtocols:
    @pytest.mark.integration
    @pytest.mark.slowtest
    @pytest.mark.parametrize("url", list(id_url.values()), ids=list(id_url.keys()))
    def test_connection(self, url):
        cs = CheckSuite()
        ds = cs.load_dataset(url)
        assert ds is not None

    # test that as_zurl can transform pointers to zarr datasets to valid nczarr urls
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

    @pytest.mark.parametrize("url_in,url_out", url_io)
    def test_as_zarr(self, url_in, url_out):
        assert zarr.as_zarr(url_in) == url_out
