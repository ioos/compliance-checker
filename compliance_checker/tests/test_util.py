#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
compliance_checker/tests/test_util.py
'''
import unittest
from compliance_checker import util


class TestUtils(unittest.TestCase):
    '''
    Test suite for utilities
    '''

    def test_datetime_is_iso(self):
        """
        Test that ISO 8601 dates are properly parsed, and non ISO 8601 dates
        are excluded
        """
        good_datetime = '2011-01-21T02:30:11Z'
        self.assertTrue(util.datetime_is_iso(good_datetime)[0])
        # we need to fail on non-ISO 8601 compliant dates
        bad_datetime = '21 Dec 2015 10:02 PM'
        self.assertFalse(util.datetime_is_iso(bad_datetime)[0])

        good_date = '2017-02-01'
        self.assertTrue(util.datetime_is_iso(good_date))
