import zipfile
from urllib.parse import urlparse
from urllib.request import url2pathname
from pathlib import Path

# 


def is_zarr(url):
    '''This check is only to be used once other protocols (is_netcdf) have come up empty\n
    '''
    if url.endswith("zarr"):
        return True

    if url.startswith('file:/'):
        return True
    
    if url.lower().startswith('s3:/'):
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

    Distinct from is_cdl etc in that it will return the appropriate URI \n\n
    
    a valid Zarr dataset could be provided in any of the following forms:\n
    "http://s3.amazonaws.com/bucket/dataset.zarr"

    "/home/path/to/dataset.zarr"
    "file:///home/path/to/dataset.zarr"
    "file:///home/path/to/dataset.randomExt#mode=nczarr,file"
    "file:///home/path/to/dataset.zarr#mode=nczarr,zip"
    '''

    pr = urlparse(str(url))

    if '#mode=nczarr' in pr.fragment:
        if pr.netloc:
            return str(url) #already valid nczarr url
        elif pr.scheme == 'file':
            return str(url) #already valid nczarr url

    zarr_url = Path(url2pathname(pr.path)).resolve() #url2pathname necessary to avoid urlparse bug in windows

    if pr.netloc:
        mode = 's3'
    elif zipfile.is_zipfile(zarr_url):
        mode = 'zip'
    elif zarr_url.is_dir():
        mode = 'file'
    else:
        raise ValueError(f'Could not identify {url},\nif #mode=nczarr,zarr, please pass this explicitly\nValid url options are described here\nhttps://www.unidata.ucar.edu/blogs/developer/entry/overview-of-zarr-support-in')

    zarr_url = f'{zarr_url.as_uri()}#mode=nczarr,{mode}'
    return zarr_url
