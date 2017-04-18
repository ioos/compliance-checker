#!/usr/bin/env python
'''
compliance_checker/protocols/cdl.py
'''
import os


def is_cdl(filename):
    '''
    Quick check for .cdl ascii file

    Example:
        netcdf sample_file {
        dimensions:
            name_strlen = 7 ;
            time = 96 ;
        variables:
            float lat ;
                lat:units = "degrees_north" ;
                lat:standard_name = "latitude" ;
                lat:long_name = "station latitude" ;
        etc...

    :param str filename: Absolute path of file to check
    :param str data: First chuck of data from file to check
    '''
    if os.path.splitext(filename)[-1] != '.cdl':
        return False

    with open(filename, 'rb') as f:
        data = f.read(32)
    if data.startswith(b'netcdf') or b'dimensions' in data:
        return True
    return False
