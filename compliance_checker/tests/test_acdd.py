import pytest

from compliance_checker.acdd import ACDDBaseCheck
from netCDF4 import Dataset

# not updated


@pytest.mark.xfail
def test_acdd():

    ds = Dataset("/Users/asadeveloper/Downloads/hycomglobalnavy_2012120300.nc")
    acdd = ACDDCheck()
    assert True == acdd.check_high(ds)
