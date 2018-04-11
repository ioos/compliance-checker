#!/usr/bin/env python
'''
compliance_checker/protocols/opendap.py

Functions to assist in determining if the URL is an OPeNDAP endpoint
'''
import requests


def is_opendap(url):
    '''
    Returns True if the URL is a valid OPeNDAP URL

    :param str url: URL for a remote OPeNDAP endpoint
    '''
    # If the server replies to a Data Attribute Structure request
    das_url = url + '.das'
    response = requests.get(das_url, allow_redirects=True)
    if 'xdods-server' in response.headers:
        return True
    # Check if it is an access restricted ESGF thredds service
    if response.status_code == 401 and \
        'text/html' in response.headers['content-type'] and \
            'The following URL requires authentication:' in response.text:
        return True
    return False
