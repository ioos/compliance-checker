import os
import subprocess
from importlib.resources import files
from itertools import chain

import pytest
from netCDF4 import Dataset

from compliance_checker.cf import util
from compliance_checker.suite import CheckSuite

datadir = files("compliance_checker").joinpath("tests/data").resolve()


def glob_down(pth, suffix, lvls):
    """globs down up to (lvls: int) levels of subfolders\n
    suffix in the form ".ipynb"\n
    pth: Path"""
    return list(chain(*[pth.glob(f'*{"/*"*lvl}{suffix}') for lvl in range(lvls)]))


def generate_dataset(cdl_path, nc_path):
    subprocess.call(["ncgen", "-4", "-o", str(nc_path), str(cdl_path)])


def static_files(cdl_stem):
    """
    Returns the Path to a valid nc dataset\n
    replaces the old STATIC_FILES dict
    """
    datadir = files("compliance_checker").joinpath("tests/data").resolve()
    assert datadir.exists(), f"{datadir} not found"

    cdl_paths = glob_down(datadir, f"{cdl_stem}.cdl", 3)
    assert (
        len(cdl_paths) > 0
    ), f"No file named {cdl_stem}.cdl found in {datadir} or its subfolders"
    assert (
        len(cdl_paths) == 1
    ), f"Multiple candidates found with the name {cdl_stem}.cdl:\n{cdl_paths}\nPlease reconcile naming conflict"
    cdl_path = cdl_paths[0]  # PurePath object

    nc_path = cdl_path.parent / f"{cdl_path.stem}.nc"
    if not nc_path.exists():
        generate_dataset(cdl_path, nc_path)
        assert (
            nc_path.exists()
        ), f"ncgen CLI utility failed to produce {nc_path} from {cdl_path}"
    return str(nc_path)


# ---------Fixtures-----------

# class scope:


@pytest.fixture
def cs(scope="class"):
    """
    Initialize the dataset
    """
    cs = CheckSuite()
    cs.load_all_available_checkers()
    return cs


@pytest.fixture
def std_names(scope="class"):
    """get current std names table version (it changes)"""
    _std_names = util.StandardNameTable()
    return _std_names


# func scope:


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
        raise OSError(f"File Exists: {nc_file_path}")
    nc = Dataset(nc_file_path, "w")
    # no need for cleanup, built-in tmpdir fixture will handle it
    return nc


@pytest.fixture
def tmp_txt_file(tmpdir):
    file_path = os.path.join(tmpdir, "output.txt")
    if os.path.exists(file_path):
        raise OSError(f"File Exists: {file_path}")

    return file_path


@pytest.fixture
def checksuite_setup():
    """For test_cli"""
    CheckSuite.checkers.clear()
    CheckSuite.load_all_available_checkers()
