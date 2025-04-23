#!/usr/bin/env python
"""
compliance_checker/tests/test_util.py
"""

from compliance_checker import util


class TestUtils:
    """
    Test suite for utilities
    """

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
