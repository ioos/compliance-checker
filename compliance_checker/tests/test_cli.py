#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Tests for command line output and parsing

'''
import io
import sys
import json

from unittest import TestCase
from compliance_checker.runner import ComplianceChecker, CheckSuite
from compliance_checker.tests.resources import STATIC_FILES
import tempfile
import os

CheckSuite.load_all_available_checkers()


class TestCLI(TestCase):
    '''
    Tests various functions and aspects of the command line tool and runner
    '''

    def setUp(self):
        _, self.path = tempfile.mkstemp()

    def tearDown(self):
        if os.path.isfile(self.path):
            os.remove(self.path)

    def shortDescription(self):
        return None

    # override __str__ and __repr__ behavior to show a copy-pastable nosetest name for ion tests
    #  ion.module:TestClassName.test_function_name
    def __repr__(self):
        name = self.id()
        name = name.split('.')
        if name[0] not in ["ion", "pyon"]:
            return "%s (%s)" % (name[-1], '.'.join(name[:-1]))
        else:
            return "%s ( %s )" % (name[-1], '.'.join(name[:-2]) + ":" + '.'.join(name[-2:]))
    __str__ = __repr__

    def test_unicode_acdd_html(self):
        '''
        Tests that the checker is capable of producing HTML with unicode characters
        '''
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=STATIC_FILES['2dim'],
            verbose=0,
            criteria='strict',
            checker_names=['acdd'],
            output_filename=self.path,
            output_format='html'
        )

        assert os.stat(self.path).st_size > 0

    def test_unicode_cf_html(self):
        '''
        Tests that the CF checker can produce HTML output with unicode characters
        '''
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=STATIC_FILES['2dim'],
            verbose=0,
            criteria='strict',
            checker_names=['cf'],
            output_filename=self.path,
            output_format='html'
        )

        assert os.stat(self.path).st_size > 0

    def test_single_json_output(self):
        '''
        Tests that a suite can produce JSON output to a file
        '''
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=STATIC_FILES['conv_bad'],
            verbose=0,
            criteria='strict',
            checker_names=['cf'],
            output_filename=self.path,
            output_format='json'
        )

        assert os.stat(self.path).st_size > 0
        with open(self.path) as f:
            r = json.load(f)
            assert 'cf' in r

    def test_multiple_json_output(self):
        '''
        Tests that a suite can produce JSON output to a file
        '''
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=STATIC_FILES['conv_bad'],
            verbose=0,
            criteria='strict',
            checker_names=['acdd', 'cf'],
            output_filename=self.path,
            output_format='json'
        )

        assert os.stat(self.path).st_size > 0
        with open(self.path) as f:
            r = json.load(f)
            assert 'cf' in r
            assert 'acdd' in r

    def test_multiple_json_output_stdout(self):
        '''
        Tests that a suite can produce JSON output to stdout
        '''
        saved = sys.stdout
        try:
            fake_stdout = io.StringIO()
            sys.stdout = fake_stdout
            return_value, errors = ComplianceChecker.run_checker(
                ds_loc=STATIC_FILES['conv_bad'],
                verbose=0,
                criteria='strict',
                checker_names=['acdd', 'cf'],
                output_filename='-',
                output_format='json'
            )
            r = json.loads(fake_stdout.getvalue().strip())
            assert 'acdd' in r
            assert 'cf' in r
        finally:
            sys.stdout = saved

    def test_single_json_output_stdout(self):
        '''
        Tests that a suite can produce JSON output to stdout
        '''
        saved = sys.stdout
        try:
            fake_stdout = io.StringIO()
            sys.stdout = fake_stdout
            return_value, errors = ComplianceChecker.run_checker(
                ds_loc=STATIC_FILES['conv_bad'],
                verbose=0,
                criteria='strict',
                checker_names=['cf'],
                output_filename='-',
                output_format='json'
            )
            r = json.loads(fake_stdout.getvalue().strip())
            assert 'cf' in r
        finally:
            sys.stdout = saved

    def test_text_output(self):
        '''
        Tests that the 'text' output can be redirected to file with arguments
        to the command line
        '''
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=STATIC_FILES['conv_bad'],
            verbose=0,
            criteria='strict',
            checker_names=['acdd', 'cf'],
            output_filename=self.path,
            output_format='text'
        )

        assert os.stat(self.path).st_size > 0
