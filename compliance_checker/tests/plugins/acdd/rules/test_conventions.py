import xarray as xr
from xrlint.testing import RuleTest, RuleTester

from compliance_checker.plugins.acdd.rules.conventions import Conventions

valid_dataset_0 = xr.Dataset(attrs={"Conventions": "ACDD-1.3"})

invalid_dataset_0 = xr.Dataset()
invalid_dataset_1 = xr.Dataset(attrs={"Conventions": "ACDD-1.2"})
invalid_dataset_2 = xr.Dataset(attrs={"Conventions": 1.3})


ConventionsTest = RuleTester.define_test(
    "1.3_conventions",
    Conventions,
    valid=[
        RuleTest(dataset=valid_dataset_0),
    ],
    invalid=[
        RuleTest(
            dataset=invalid_dataset_0,
            expected=["Missing attribute 'Conventions'."],
        ),
        RuleTest(
            dataset=invalid_dataset_1,
            expected=[
                "Attribute 'Conventions' needs to contain 'ACDD-1.3' in addition to the current: 'ACDD-1.2'",
            ],
        ),
        RuleTest(
            dataset=invalid_dataset_2,
            expected=["Invalid attribute 'Conventions': 1.3"],
        ),
    ],
)
