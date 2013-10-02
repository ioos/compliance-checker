"""
Compliance Checker suite runner
"""

import inspect
import itertools
from netCDF4 import Dataset
from lxml import etree as ET
from wicken.netcdf_dogma import NetCDFDogma

class CheckSuite(object):

    def _get_checks(self, checkclass):
        meths = inspect.getmembers(checkclass, inspect.ismethod)
        return [x[1] for x in meths if x[0].startswith("check_")]

    def run(self, dataset, *args):

        all_checks = []
        map(all_checks.extend, (self._get_checks(a) for a in args))


        vals = [[v] if not isinstance(v, list) else v for v in [c(dataset) for c in all_checks]]

        # transform to scores

        return list(itertools.chain.from_iterable(vals))

    def load_dataset(self, ds_str, belief_map):
        ds = Dataset(ds_str)
        data_object = NetCDFDogma('ds', belief_map, ds)

        return data_object

