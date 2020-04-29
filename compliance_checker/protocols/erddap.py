import io
import urllib.request
import compliance_checker.protocols.opendap as opendap

def is_tabledap(url):
    """
    Identify a dataset as an ERDDAP TableDAP dataset.

    Parameters
    ----------
    url (str) : URL to dataset

    Returns
    -------
    bool
    """

    if "tabledap" in url:
        return True
    return False

def get_tabledap_bytes(url, ftype):
    """
    ERDDAP TableDAP returns an OPeNDAP "sequence" response by default
    when no file extensions are provided. If a user wishes to get a dataset
    from an ERDDAP TableDAP URL, append the desired file extension and return
    a byte buffer object containing the data.

    Parameters
    ----------
    url (str)   : URL to TableDAP dataset
    ftype (str) : file format extension

    Return
    ------
    io.BytesIO buffer object
    """

    vstr = opendap.create_DAP_variable_str(url) # variable str from DDS
    _url = f'{".".join([url, ftype])}?{vstr}'
    with urllib.request.urlopen(_url, timeout=120) as resp:
        return io.BytesIO(resp.read())
