import xarray as xr
from xrlint.testing import RuleTest, RuleTester

from compliance_checker.plugins.acdd.rules.iso_dates import IsoDates

valid_dataset_0 = xr.Dataset(attrs={"date_created": "2023-10-05T12:34:56Z"})
valid_dataset_1 = xr.Dataset(attrs={"date_modified": "2023-10-05"})
valid_dataset_2 = xr.Dataset(attrs={"date_issued": "2023-10-05T12:34:56+00:00"})
valid_dataset_3 = xr.Dataset(
    attrs={"date_metadata_modified": "2023-10-05T12:34:56-05:00"},
)
valid_dataset_4 = xr.Dataset()  # No date attributes


invalid_dataset_0 = xr.Dataset(attrs={"date_created": "10/05/2023"})
invalid_dataset_1 = xr.Dataset(attrs={"date_issued": "2023-13-05"})
invalid_dataset_2 = xr.Dataset(attrs={"date_modified": "2023-10-05 12:34:56"})
invalid_dataset_3 = xr.Dataset(attrs={"date_metadata_modified": "2023/10/05"})
invalid_dataset_4 = xr.Dataset(
    attrs={"date_created": "2023-10-05T25:34:56Z"},
)  # Invalid hour

IsoDatesTest = RuleTester.define_test(
    "1.3_dates_iso_format",
    IsoDates,
    valid=[
        RuleTest(dataset=valid_dataset_0),
        RuleTest(dataset=valid_dataset_1),
        RuleTest(dataset=valid_dataset_2),
        RuleTest(dataset=valid_dataset_3),
        RuleTest(dataset=valid_dataset_4),
    ],
    invalid=[
        RuleTest(
            dataset=invalid_dataset_0,
            expected=["Attribute 'date_created' is not in ISO format: '10/05/2023'"],
        ),
        RuleTest(
            dataset=invalid_dataset_1,
            expected=["Attribute 'date_issued' is not in ISO format: '2023-13-05'"],
        ),
        RuleTest(
            dataset=invalid_dataset_2,
            expected=[
                "Attribute 'date_modified' is not in ISO format: '2023-10-05 12:34:56'",
            ],
        ),
        RuleTest(
            dataset=invalid_dataset_3,
            expected=[
                "Attribute 'date_metadata_modified' is not in ISO format: '2023/10/05'",
            ],
        ),
        RuleTest(
            dataset=invalid_dataset_4,
            expected=[
                "Attribute 'date_created' is not in ISO format: '2023-10-05T25:34:56Z'",
            ],
        ),
    ],
)
