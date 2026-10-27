"""
compliance_checker/tests/test_util.py
"""

import pytest

from compliance_checker import util
from compliance_checker.tests.helpers import MockNetCDF


class TestUtils:
    """
    Test suite for utilities
    """

    @pytest.fixture
    def grouped_dataset(self):
        with MockNetCDF() as root:
            root.createVariable("temperature", "f4").units = "K"
            ocean = root.createGroup("ocean")
            ocean.createVariable("temperature", "f4").units = "K"
            ocean.createVariable("unlabelled", "f4")
            instruments = ocean.createGroup("instruments")
            instruments.createVariable("identifier", "i4").long_name = "Instrument identifier"
            atmosphere = root.createGroup("atmosphere")
            atmosphere.createVariable("temperature", "f4").units = "degC"
            root.createGroup("empty")
            yield root

    def test_group_order(self, grouped_dataset):
        assert [group.path for group in util.iter_groups(grouped_dataset)] == [
            "/",
            "/ocean",
            "/ocean/instruments",
            "/atmosphere",
            "/empty",
        ]

    def test_subtree(self, grouped_dataset):
        ocean = grouped_dataset.groups["ocean"]
        assert [group.path for group in util.iter_groups(ocean)] == [
            "/ocean",
            "/ocean/instruments",
        ]
        assert [variable.name for variable in util.iter_variables(ocean)] == [
            "temperature",
            "unlabelled",
            "identifier",
        ]

    def test_all_variables(self, grouped_dataset):
        assert [(variable.group().path, variable.name) for variable in util.iter_variables(grouped_dataset)] == [
            ("/", "temperature"),
            ("/ocean", "temperature"),
            ("/ocean", "unlabelled"),
            ("/ocean/instruments", "identifier"),
            ("/atmosphere", "temperature"),
        ]

    def test_attribute_filters(self, grouped_dataset):
        variables = list(util.iter_variables(grouped_dataset, units="K"))
        assert [variable.group().path for variable in variables] == ["/", "/ocean"]
        assert list(util.iter_variables(grouped_dataset, units="missing")) == []

    def test_callable_filters(self, grouped_dataset):
        variables = list(
            util.iter_variables(
                grouped_dataset,
                units=lambda value: value is None,
                long_name=lambda value: value is not None,
            ),
        )
        assert [(variable.group().path, variable.name) for variable in variables] == [
            ("/ocean/instruments", "identifier"),
        ]

    def test_empty_dataset(self):
        with MockNetCDF() as root:
            assert list(util.iter_groups(root)) == [root]
            assert list(util.iter_variables(root)) == []

    def test_datetime_is_iso(self):
        """
        Test that ISO 8601 dates are properly parsed, and non ISO 8601 dates
        are excluded
        """
        good_datetime = "2011-01-21T02:30:11Z"
        assert util.datetime_is_iso(good_datetime)[0]

        good_date = "2017-02-01"
        assert util.datetime_is_iso(good_date)[0]

        good_date = "20170201"
        assert util.datetime_is_iso(good_date)[0]

        good_datetime = "2017-09-19T23:06:17+00:00"
        assert util.datetime_is_iso(good_datetime)[0]

        good_datetime = "20170919T230617Z"
        assert util.datetime_is_iso(good_datetime)[0]

        good_date = "2017-W38"
        assert util.datetime_is_iso(good_date)[0]

        good_date = "2017-W38-2"
        assert util.datetime_is_iso(good_date)[0]

        good_date = "2017-262"
        assert util.datetime_is_iso(good_date)[0]

        # we need to fail on non-ISO 8601 compliant dates
        bad_datetime = "21 Dec 2015 10:02 PM"
        assert not util.datetime_is_iso(bad_datetime)[0]

        # Month first is not ISO-8601 compliant
        bad_datetime = "09192017T230617Z"
        assert not util.datetime_is_iso(bad_datetime)[0]

        bad_date = "09192017"
        assert not util.datetime_is_iso(bad_date)[0]
