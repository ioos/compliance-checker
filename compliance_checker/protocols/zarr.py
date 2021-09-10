import zipfile
from urllib.parse import urlparse
from pathlib import Path

# a valid Zarr dataset could be provided in any of the following forms:
"http://s3.amazonaws.com/bucket/dataset.zarr"

"/home/path/to/dataset.zarr"
"file:///home/path/to/dataset.zarr"
"file:///home/path/to/dataset.zarr#mode=nczarr,file"
"file:///home/path/to/dataset.zarr#mode=nczarr,zip"


def is_zarr(url):
    '''This check is only to be used once other protocols (is_netcdf) have come up empty\n
    Distinct from is_cdl etc in that it will return the appropriate URI '''
    if url.endswith("zarr"):
        return True

    if url.startswith('file:/'):
        return True

    if zipfile.is_zipfile(url):
        # if it's a folder or zip, assume it is a zarr
        return True
    
    if Path(url).is_dir():
        return True

    return False

def as_zarr(url):
    '''
    
    https://www.unidata.ucar.edu/blogs/developer/entry/overview-of-zarr-support-in
    '''
    pr = urlparse(str(url))
    zarr_url = Path(pr.path).resolve()
    mode = 'zip' if zipfile.is_zipfile(url) else 'file'

    zarr_url = f'{zarr_url.as_uri()}#mode=nczarr,{mode}'
    return zarr_url
