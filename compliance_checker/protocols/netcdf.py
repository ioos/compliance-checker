#!/usr/bin/env python
'''
compliance_checker/protocols/netcdf.py

Functions to assist in determining if the URL points to a netCDF file
'''


def is_netcdf(url):
    '''
    Returns True if the URL points to a valid local netCDF file

    :param str url: Location of file on the file system
    '''
    # Try an obvious exclusion of remote resources
    if url.startswith('http'):
        return False

    # If it's a known extension, give it a shot
    if url.endswith('nc'):
        return True

    # Brute force
    with open(url, 'rb') as f:
        magic_number = f.read(4)
        if len(magic_number) < 4:
            return False
        if is_classic_netcdf(magic_number):
            return True
        elif is_hdf5(magic_number):
            return True

        return False


def is_classic_netcdf(file_buffer):
    '''
    Returns True if the contents of the byte array matches the magic number in
    netCDF files

    :param str file_buffer: Byte-array of the first 4 bytes of a file
    '''
    # CDF.
    if file_buffer == b'\x43\x44\x46\x01':
        return True
    return False


def is_hdf5(file_buffer):
    '''
    Returns True if the contents of the byte array matches the magic number in
    HDF5 files

    :param str file_buffer: Byte-array of the first 4 bytes of a file
    '''
    # .HDF
    if file_buffer == b'\x89\x48\x44\x46':
        return True
    return False
