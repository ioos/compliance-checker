import pytest
from pathlib import Path
from netCDF4 import Dataset

from compliance_checker.cf import util
from compliance_checker.suite import CheckSuite

import os
import subprocess

from pkg_resources import resource_filename

from functools import reduce
from operator import iconcat
def flatten_list(listOfLists):
    '''or list of tups'''
    return reduce(iconcat,listOfLists,[])
def glob_down(pth,suffix,lvls):
    '''globs "recursively" down (lvls: int) levels\n
    suffix in the form ".ipynb"\n
    pth: Path'''
    return flatten_list( [list(pth.glob(f'*{"/*"*lvl}{suffix}')) for lvl in range(lvls)] )

def generate_dataset(cdl_path, nc_path):
    subprocess.call([ "ncgen", "-o", str(nc_path), str(cdl_path) ])

def static_files(cdl_stem):
    """
    Returns the Path to a valid nc dataset\n
    replaces the old STATIC_FILES dict
    """
    datadir = Path()/'data'
    assert datadir.exists(), f'{datadir.resolve()} not found'

    cdl_path = glob_down(datadir,f'{cdl_stem}.cdl',3)[0].resolve()
    #PurePath object
    nc_path = cdl_path.parent/f'{cdl_path.stem}.nc'
    if not nc_path.exists():
        generate_dataset(cdl_path, nc_path)
        assert nc_path.exists(), f'ncgen CLI utility failed to produce {nc_path} from {cdl_path}'
    return nc_path

#---------Fixtures-----------

@pytest.fixture
def cs(scope='class'):
    """
    Initialize the dataset
    """
    cs = CheckSuite()
    cs.load_all_available_checkers()
    return cs

@pytest.fixture
def std_names(scope='class'):
    '''get current std names table version (it changes)'''
    _std_names = util.StandardNameTable()
    return _std_names

@pytest.fixture
def loaded_dataset(request):
    """
    Return a loaded NC Dataset for the given path\n
    nc_dataset_path parameterized for each test
    """
    nc_dataset_path = static_files(request.param)
    nc = Dataset(nc_dataset_path, "r")
    yield nc
    nc.close()

@pytest.fixture
def new_nc_file(tmpdir):
    """
    Make a new temporary netCDF file for the scope of the test
    """
    nc_file_path = os.path.join(tmpdir, "example.nc")
    if os.path.exists(nc_file_path):
        raise IOError("File Exists: %s" % nc_file_path)
    nc = Dataset(nc_file_path, "w")
    #no need for cleanup, built-in tmpdir fixture will handle it
    return nc

