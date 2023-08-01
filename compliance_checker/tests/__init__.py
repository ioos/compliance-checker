import unittest
from pathlib import Path

from netCDF4 import Dataset


class BaseTestCase(unittest.TestCase):
    """
    Base test case compliance checker unit tests
    """

    def shortDescription(self):
        return None

    def __repr__(self):
        name = self.id()
        name = name.split(".")
        return "{} ( {} )".format(
            name[-1],
            ".".join(name[:-2]) + ":" + ".".join(name[-2:]),
        )

    __str__ = __repr__

    def load_dataset(self, nc_dataset):
        """
        Return a loaded NC Dataset for the given path
        """
        if not isinstance(nc_dataset, (str, Path)):
            raise ValueError("nc_dataset should be a valid path")

        nc_dataset = Dataset(nc_dataset, "r")
        self.addCleanup(nc_dataset.close)
        return nc_dataset

    def assert_result_is_good(self, result):
        if isinstance(result.value, bool):
            assert result.value is True
        else:
            assert result.value[0] == result.value[1]

    def assert_result_is_bad(self, result):
        if isinstance(result.value, bool):
            assert result.value is False
        else:
            assert result.value[0] != result.value[1]


class pytestBaseTest:
    """Base test case for migrating tests away from unittest and towards pytest only"""

    def shortDescription(self):
        return None

    def __repr__(self):
        name = self.id()
        name = name.split(".")
        return "{} ( {} )".format(
            name[-1],
            ".".join(name[:-2]) + ":" + ".".join(name[-2:]),
        )

    __str__ = __repr__

    def assert_result_is_good(self, result):
        if isinstance(result.value, bool):
            assert result.value is True
        else:
            assert result.value[0] == result.value[1]

    def assert_result_is_bad(self, result):
        if isinstance(result.value, bool):
            assert result.value is False
        else:
            assert result.value[0] != result.value[1]
