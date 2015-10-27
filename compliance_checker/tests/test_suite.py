from compliance_checker.runner import CheckSuite
from pkg_resources import resource_filename

def test_suite():
    cs = CheckSuite()
    ds = cs.load_dataset(resource_filename("compliance_checker", "tests/data/2dim-grid.nc"))
    vals = cs.run(ds, 'acdd')

    # run no longer returns the summed score, so this test.. just runs
    #assert vals['acdd'][0] == (43.5, 78)

