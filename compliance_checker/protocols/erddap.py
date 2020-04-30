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

    return "tabledap" in url
