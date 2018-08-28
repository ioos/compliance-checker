from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from netCDF4 import Dataset
try:
    from functools import lru_cache
# Fallback for Python < 3.2
except ImportError:
    from functools32 import lru_cache

class MemoizedDataset(Dataset):
    """
    A NetCDF dataset which has its get_variables_by_attributes call memoized in
    order to speed up repeated calls to the function.  This should only really
    be used against netCDF Datasets opened in 'r' mode, as the attributes should
    not change upon reading the files.
    """
    @lru_cache(128)
    def get_variables_by_attributes(self, **kwargs):
        return super(MemoizedDataset,
                     self).get_variables_by_attributes(**kwargs)
