import xarray as xr
from xrlint.testing import RuleTest, RuleTester

from compliance_checker.plugins.acdd.rules.no_id_blanks import NoBlanksInID

valid_dataset_0 = xr.Dataset(attrs={"id": "testind_dataset"})

invalid_dataset_0 = xr.Dataset()
invalid_dataset_1 = xr.Dataset(attrs={"id": "testing dataset"})


IdBlanksTest = RuleTester.define_test(
    "1.3_no_blanks_in_id",
    NoBlanksInID,
    valid=[RuleTest(dataset=valid_dataset_0)],
    invalid=[
        RuleTest(
            dataset=invalid_dataset_0,
            expected=["Missing attribute 'id'"],
        ),
        RuleTest(
            dataset=invalid_dataset_1,
            expected=["There should be no blanks in the id field"],
        ),
    ],
)
