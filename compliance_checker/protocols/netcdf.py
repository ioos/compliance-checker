#!/usr/bin/env python
"""
compliance_checker/protocols/netcdf.py

Functions to assist in determining if the URL points to a netCDF file
"""

import warnings

import requests


def is_netcdf(url):
    """
    Returns True if the URL points to a valid local netCDF file

    :param str url: Location of file on the file system
    """
    # Try an obvious exclusion of remote resources
    url = str(url)
    if url.startswith("http"):
        return False

    # If it's a known extension, give it a shot
    if url.endswith("nc"):
        return True

    try:
        # Brute force
        with open(url, "rb") as f:
            magic_number = f.read(4)
            if len(magic_number) < 4:
                return False
            if is_classic_netcdf(magic_number):
                return True
            elif is_hdf5(magic_number):
                return True
    except PermissionError:
        # open will fail for both a directory or a local url, either of which may be pointing to a Zarr dataset
        # directory
        return False
    except OSError:
        # local file url
        return False

    return False


def is_classic_netcdf(file_buffer):
    """
    Returns True if the contents of the byte array matches the magic number in
    netCDF files

    :param str file_buffer: Byte-array of the first 4 bytes of a file
    """
    # CDF.
    if file_buffer == b"\x43\x44\x46\x01":
        return True
    return False


def is_hdf5(file_buffer):
    """
    Returns True if the contents of the byte array matches the magic number in
    HDF5 files

    :param str file_buffer: Byte-array of the first 4 bytes of a file
    """
    # .HDF
    if file_buffer == b"\x89\x48\x44\x46":
        return True
    return False


def is_remote_netcdf(ds_str):
    """
    Check a remote path points to a NetCDF resource.

    Parameters
    ----------
    ds_str (str): remote path to a dataset

    Returns
    -------
    bool
    """

    # Some datasets do not support HEAD requests!  The vast majority will,
    # however, support GET requests
    try:
        head_req = requests.head(ds_str, allow_redirects=True, timeout=10)
        head_req.raise_for_status()
    except requests.exceptions.RequestException as e:
        warnings.warn(
            f"Received exception when making HEAD request to {ds_str}: {e}",
            stacklevel=2,
        )
        content_type = None
    else:
        content_type = head_req.headers.get("content-type")

    if content_type is None:
        return False

    # if the Content-Type header returned was "application/x-netcdf",
    # or a netCDF file (not OPeNDAP) we can open this into a Dataset
    # Add support for application/x-netcdf;ver=4
    return content_type.split(";")[0] == "application/x-netcdf"
