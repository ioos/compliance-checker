#!/usr/bin/env python
"""
compliance_checker/protocols/opendap.py

Functions to assist in determining if the URL is an OPeNDAP endpoint
"""
import urllib.parse
import urllib.request

import requests


def create_DAP_variable_str(url):
    """
    Create a URL-encoded string of variables for a given DAP dataset.
    Works on OPeNDAP datasets.

    Parameters
    ----------
    url (str): endpoint to *DAP dataset

    Returns
    -------
    str
    """

    # get dds
    with urllib.request.urlopen(f"{url}.dds") as resp:
        _str = resp.read().decode()[8:]

    # remove beginning and ending braces, split on newlines
    no_braces_newlines = list(
        filter(lambda x: "{" not in x and "}" not in x, _str.split("\n")),
    )

    # remove all the extra space used in the DDS string
    no_spaces = list(filter(None, (x.strip(" ") for x in no_braces_newlines)))

    # now need to split from type, grab only the variable and remove ;
    vars_only = [x.split(" ")[-1].strip(";") for x in no_spaces]

    # encode as proper URL characters
    varstr = urllib.parse.quote(",".join(vars_only))
    return varstr


def is_opendap(url):
    """
    Returns True if the URL is a valid OPeNDAP URL

    :param str url: URL for a remote OPeNDAP endpoint
    """
    # If the server replies to a Data Attribute Structure request
    if url.endswith("#fillmismatch"):
        das_url = url.replace("#fillmismatch", ".das")
    else:
        das_url = url + ".das"

    try:
        response = requests.get(das_url, allow_redirects=True)

        if "xdods-server" in response.headers:
            return True
        # Check if it is an access restricted ESGF thredds service
        if (
            response.status_code == 401
            and "text/html" in response.headers["content-type"]
            and "The following URL requires authentication:" in response.text
        ):
            return True
    except requests.exceptions.InvalidSchema:
        return False  # not opendap if url + ".das" isn't found
    return False
