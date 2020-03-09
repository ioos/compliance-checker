import pkg_resources
from functools import lru_cache

from netCDF4 import Dataset

try:
    __version__ = pkg_resources.get_distribution("compliance_checker").version
except Exception:
    __version__ = "unknown"

from contextlib import contextmanager
from tempfile import NamedTemporaryFile
from typing import BinaryIO, Generator


class MemoizedDataset(Dataset):
    """
    A NetCDF dataset which has its get_variables_by_attributes call memoized in
    order to speed up repeated calls to the function.  This should only really
    be used against netCDF Datasets opened in 'r' mode, as the attributes should
    not change upon reading the files.
    """

    @lru_cache(128)
    def get_variables_by_attributes(self, **kwargs):
        return super(MemoizedDataset, self).get_variables_by_attributes(**kwargs)

@contextmanager
def tempnc(data: BinaryIO) -> Generator[str, None, None]:
    """
    Create a temporary in-memory NetCDF file using a NamedTemporaryFile.
    Close the file automatically after scope is exited.

    Type aliasing and tempfile creation credit to @ocefpaf
    https://github.com/ioos/compliance-checker/pull/799#discussion_r411420587

    Parameters
    ----------
    data (bytes): raw bytes to store in the NamedTemporaryFile

    Returns
    -------
    context-managed generator
    """
    tmp = None
    try:
        tmp = NamedTemporaryFile(suffix=".nc", prefix="compliance-checker_")
        tmp.write(data)
        tmp.flush()
        yield tmp.name
    finally:
        if tmp is not None:
            tmp.close()
