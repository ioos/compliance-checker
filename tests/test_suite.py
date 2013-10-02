from compliance_checker.suite import CheckSuite
from compliance_checker.acdd import ACDDCheck
from netCDF4 import Dataset
from pprint import pprint

def test_suite():
    cs = CheckSuite()
    ds = cs.load_dataset("/Users/asadeveloper/Downloads/hycomglobalnavy_2012120300.nc", ACDDCheck.beliefs())
    #ds = cs.load_dataset("/Users/asadeveloper/Downloads/hycom.ncml", ACDDCheck.beliefs)
    acdd = ACDDCheck()

    vals = cs.run(ds, acdd)

    pprint(vals)
    assert vals == []

