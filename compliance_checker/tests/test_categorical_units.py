"""
Test that categorical standard names with no canonical units
(e.g. soil_type, region, area_type) do not raise false units errors.
See: https://github.com/ioos/compliance-checker/issues/1219
"""
import numpy as np
import pytest
import netCDF4 as nc
import tempfile
import os
from compliance_checker.cf.cf_1_6 import CF1_6Check

def get_messages(ds):
    checker = CF1_6Check()
    results = checker.check_units(ds)
    messages = []
    for r in results:
        messages.extend(r.msgs)
    return messages

def test_soil_type_units_1_no_error():
    """units='1' for soil_type should not raise a canonical units error."""
    tmp = tempfile.NamedTemporaryFile(suffix=".nc", delete=False)
    tmp.close()
    try:
        ds = nc.Dataset(tmp.name, "w")
        ds.createDimension("x", 3)
        var = ds.createVariable("SOILTYP", "i4", ("x",))
        var.standard_name = "soil_type"
        var.units = "1"
        messages = get_messages(ds)
        bad = [m for m in messages if "canonical units" in m and "SOILTYP" in m]
        assert not bad, f"Unexpected canonical units errors: {bad}"
        ds.close()
    finally:
        os.unlink(tmp.name)

def test_soil_type_no_units_no_error():
    """Missing units for soil_type should not raise a units required error."""
    tmp = tempfile.NamedTemporaryFile(suffix=".nc", delete=False)
    tmp.close()
    try:
        ds = nc.Dataset(tmp.name, "w")
        ds.createDimension("x", 3)
        var = ds.createVariable("SOILTYP", "i4", ("x",))
        var.standard_name = "soil_type"
        # No units attribute intentionally
        messages = get_messages(ds)
        bad = [m for m in messages if "units attribute is required" in m and "SOILTYP" in m]
        assert not bad, f"Unexpected units required errors: {bad}"
        ds.close()
    finally:
        os.unlink(tmp.name)

def test_area_type_no_units_no_error():
    """area_type also has no canonical units — should not raise errors."""
    tmp = tempfile.NamedTemporaryFile(suffix=".nc", delete=False)
    tmp.close()
    try:
        ds = nc.Dataset(tmp.name, "w")
        ds.createDimension("x", 3)
        var = ds.createVariable("ATYPE", "i4", ("x",))
        var.standard_name = "area_type"
        var.units = "1"
        messages = get_messages(ds)
        bad = [m for m in messages if "canonical units" in m and "ATYPE" in m]
        assert not bad, f"Unexpected canonical units errors: {bad}"
        ds.close()
    finally:
        os.unlink(tmp.name)
