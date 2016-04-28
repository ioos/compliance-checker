#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Tests for base compliance checker class'''

from __future__ import unicode_literals
from unittest import TestCase
from netCDF4 import Dataset
from compliance_checker import base
import os

class TestBase(TestCase):
    '''
    Tests functionality of the base compliance checker class
    '''

    def setUp(self):
        self.acdd = base.BaseCheck()
        self.ds = Dataset(filename=os.devnull, mode='w', diskless=True)

    def tearDown(self):
        self.ds.close()

    def test_attr_presence(self):
        # attribute not present, should fail
        priority = base.BaseCheck.MEDIUM
        rv1, rv2, rv3, rv4 = [], [], [], []
        attr = 'test'
        base.attr_check(attr, self.ds, priority, rv1)
        assert rv1[0] == base.Result(priority, False, 'test',
                                  ['Attr test not present'])
        # test with empty string
        self.ds.test = ''
        base.attr_check(attr, self.ds, priority, rv2)
        assert rv2[0] == base.Result(priority, False, 'test',
                                ["Attr test is empty or completely whitespace"])
        # test with whitespace in the form of a space and a tab
        self.ds.test = ' 	'
        base.attr_check(attr, self.ds, priority, rv3)
        assert rv3[0] == base.Result(priority, False, 'test',
                                ["Attr test is empty or completely whitespace"])
        # test with actual string contents
        self.ds.test = 'abc 123'
        base.attr_check(attr, self.ds, priority, rv4)
        assert rv4[0] == base.Result(priority, True, 'test', [])

    def test_attr_in_valid_choices(self):
        """Tests attribute membership in a set"""
        rv1, rv2, rv3 = [], [], []
        priority = base.BaseCheck.MEDIUM
        valid_choices = ['a', 'b', 'c']
        attr = ('test', valid_choices)
        base.attr_check(attr, self.ds, priority, rv1)
        assert rv1[0] == base.Result(priority, (0, 2), 'test', ["Attr test not present"])
        self.ds.test = ''
        base.attr_check(attr, self.ds, priority, rv2)
        assert rv2[0] == base.Result(priority, (1, 2), 'test', ["Attr test present, but not in expected value list (%s)" % valid_choices])
        self.ds.test = 'a'
        base.attr_check(attr, self.ds, priority, rv3)
        assert rv3[0] == base.Result(priority, (2, 2), 'test', [])

    def test_attr_fn(self):
        """Test attribute against a checker function"""
        # simple test.  In an actual program, this use case would be covered
        rv1, rv2, rv3 = [], [], []
        priority = base.BaseCheck.MEDIUM

        def verify_dummy(ds):
            if ds.dummy + 'y' == 'dummyy':
                return base.ratable_result(True, 'dummy', [])
            else:
                return base.ratable_result(False, 'dummy', ['not "dummyy"'])

        attr = ('dummy', verify_dummy)
        base.attr_check(attr, self.ds, priority, rv1)
        assert rv1[0] == base.Result(priority, False, 'dummy',
                                  ['Attr dummy not present'])
        self.ds.dummy = 'doomy'
        base.attr_check(attr, self.ds, priority, rv2)
        assert rv2[0] == base.Result(priority, False, 'dummy', ['not "dummyy"'])
        self.ds.dummy = 'dummy'
        base.attr_check(attr, self.ds, priority, rv3)
        assert rv3[0] == base.Result(priority, True, 'dummy', [])
