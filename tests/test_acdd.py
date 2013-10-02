import pytest

from compliance_checker.acdd import ACDDCheck
from netCDF4 import Dataset

"""
class TestACDDCheck:
    def __init__(self):
        self.ds = Dataset("/Users/asadeveloper/Downloads/hycomglobalnavy_2012120300.nc")
        self.acdd = ACDDCheck()

    def test_acdd(self):
        assert True == self.acdd.check_high(self.ds)
"""



def test_acdd():
    ds = Dataset("/Users/asadeveloper/Downloads/hycomglobalnavy_2012120300.nc")
    acdd = ACDDCheck()
    assert True == acdd.check_high(ds)
