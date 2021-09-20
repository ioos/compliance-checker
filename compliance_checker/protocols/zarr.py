from compliance_checker.protocols import netcdf
import zipfile
from zipfile import ZipFile
from urllib.parse import urlparse
from urllib.request import url2pathname
from pathlib import Path

# 


def is_zarr(url):
    '''
    '''

    if netcdf.is_netcdf(url):
        return False

    if '.zarr' in url:
        return True

    if urlparse(url).scheme in ('https','s3','file'):
        return True
    
    if zipfile.is_zipfile(url):
        if '.zmetadata' in ZipFile(url).namelist():
            return True
    
    if Path(url).is_dir():
        if (Path(url)/'.zmetadata').exists():
            return True

    return False

def as_zarr(url):
    '''
    Transform pointers to zarr datasets to valid nczarr urls, as described in
    https://www.unidata.ucar.edu/blogs/developer/entry/overview-of-zarr-support-in\n
    url: str or Path to valid zarr dataset\n
    Distinct from is_cdl etc in that it will return the appropriate URI \n\n
    
    A valid Zarr dataset could be provided in any of the following forms:\n
    "http://s3.amazonaws.com/bucket/dataset.zarr"\n
    "http://s3.amazonaws.com/bucket/dataset.zarr"#mode=nczarr,s3\n
    "/home/path/to/dataset.zarr"\n
    Path('/home/path/to/dataset.zarr')\n
    "file:///home/path/to/dataset.zarr"\n
    "file:///home/path/to/dataset.randomExt#mode=nczarr,file"
    "file:///home/path/to/dataset.zarr#mode=nczarr,zip"
    '''

    pr = urlparse(str(url))

    if 'mode=nczarr' in pr.fragment:
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

    url_base = url if mode=='s3' else zarr_url.as_uri()

    zarr_url = f'{url_base}#mode=nczarr,{mode}'
    return zarr_url

