from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from netCDF4 import Dataset
try:
    from functools import lru_cache
# Fallback for Python < 3.2
except ImportError:
    from repoze import lru_cache

class MemoizedDataset(Dataset):
    @lru_cache(128)
    def get_variables_by_attributes(self, **kwargs):
        return super(MemoizedDataset,
                     self).get_variables_by_attributes(**kwargs)
