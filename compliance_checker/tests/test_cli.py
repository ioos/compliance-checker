#!/usr/bin/env python
"""
Tests for command line output and parsing

"""
import io
import json
import os
import sys
import tempfile
from argparse import Namespace
from unittest import TestCase

import pytest

from compliance_checker.runner import CheckSuite, ComplianceChecker
from compliance_checker.tests.resources import STATIC_FILES


class TestCLI(TestCase):
    """
    Tests various functions and aspects of the command line tool and runner
    """

    def setUp(self):
        self.fid, self.path = tempfile.mkstemp()
        # why is the class being written to
        CheckSuite.checkers.clear()
        CheckSuite.load_all_available_checkers()

    def tearDown(self):
        if os.path.isfile(self.path):
            os.close(self.fid)
            os.remove(self.path)

    def shortDescription(self):
        return None

    # override __str__ and __repr__ behavior to show a copy-pastable nosetest name for ion tests
    #  ion.module:TestClassName.test_function_name
    def __repr__(self):
        name = self.id()
        name = name.split(".")
        if name[0] not in ["ion", "pyon"]:
            return "{} ({})".format(name[-1], ".".join(name[:-1]))
        else:
            return "{} ( {} )".format(
                name[-1],
                ".".join(name[:-2]) + ":" + ".".join(name[-2:]),
            )

    __str__ = __repr__

    def test_unicode_acdd_html(self):
        """
        Tests that the checker is capable of producing HTML with unicode characters
        """
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=STATIC_FILES["2dim"],
            verbose=0,
            criteria="strict",
            checker_names=["acdd"],
            output_filename=self.path,
            output_format="html",
        )

        assert os.stat(self.path).st_size > 0

    def test_unicode_cf_html(self):
        """
        Tests that the CF checker can produce HTML output with unicode characters
        """
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=STATIC_FILES["2dim"],
            verbose=0,
            criteria="strict",
            checker_names=["cf"],
            output_filename=self.path,
            output_format="html",
        )

        assert os.stat(self.path).st_size > 0

    def test_single_json_output(self):
        """
        Tests that a suite can produce JSON output to a file
        """
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=STATIC_FILES["conv_bad"],
            verbose=0,
            criteria="strict",
            checker_names=["cf"],
            output_filename=self.path,
            output_format="json",
        )

        assert os.stat(self.path).st_size > 0
        with open(self.path) as f:
            r = json.load(f)
            assert "cf" in r

    def test_list_checks(self):
        """
        Tests listing of both old-style, deprecated checkers using .name
        attributes, and newer ones which use ._cc_spec and _cc_spec_version
        attributes
        """

        # hack: use argparse.Namespace to mock checker object with attributes
        # since SimpleNamespace is Python 3.3+ only
        CheckSuite.checkers.clear()
        # need to mock setuptools entrypoints here in order to load in checkers
        def checker_1():
            return Namespace(name="checker_1")

        def checker_2():
            return Namespace(_cc_spec="checker_2", _cc_spec_version="2.2")

        mock_checkers = [Namespace(resolve=checker_1), Namespace(resolve=checker_2)]
        with pytest.warns(DeprecationWarning):
            CheckSuite._load_checkers(mock_checkers)

        cs = CheckSuite()
        saved = sys.stdout
        try:
            # ugly!  consider refactoring to use py.test capsys
            fake_stdout = io.StringIO()
            sys.stdout = fake_stdout
            cs._print_suites()
            assert fake_stdout.getvalue() == " - checker_1:unknown\n - checker_2:2.2\n"
        finally:
            sys.stdout = saved
            fake_stdout.close()

    def test_multiple_json_output(self):
        """
        Tests that a suite can produce JSON output to a file
        """
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=STATIC_FILES["conv_bad"],
            verbose=0,
            criteria="strict",
            checker_names=["acdd", "cf"],
            output_filename=self.path,
            output_format="json",
        )

        assert os.stat(self.path).st_size > 0
        with open(self.path) as f:
            r = json.load(f)
            assert "cf" in r
            assert "acdd" in r

    def test_multiple_json_output_stdout(self):
        """
        Tests that a suite can produce JSON output to stdout
        """
        saved = sys.stdout
        try:
            fake_stdout = io.StringIO()
            sys.stdout = fake_stdout
            return_value, errors = ComplianceChecker.run_checker(
                ds_loc=STATIC_FILES["conv_bad"],
                verbose=0,
                criteria="strict",
                checker_names=["acdd", "cf"],
                output_filename="-",
                output_format="json",
            )
            r = json.loads(fake_stdout.getvalue().strip())
            assert "acdd" in r
            assert "cf" in r
        finally:
            sys.stdout = saved
            fake_stdout.close()

    def test_single_json_output_stdout(self):
        """
        Tests that a suite can produce JSON output to stdout
        """
        saved = sys.stdout
        try:
            fake_stdout = io.StringIO()
            sys.stdout = fake_stdout
            return_value, errors = ComplianceChecker.run_checker(
                ds_loc=STATIC_FILES["conv_bad"],
                verbose=0,
                criteria="strict",
                checker_names=["cf"],
                output_filename="-",
                output_format="json",
            )
            r = json.loads(fake_stdout.getvalue().strip())
            assert "cf" in r
        finally:
            sys.stdout = saved
            fake_stdout.close()

    def test_text_output(self):
        """
        Tests that the 'text' output can be redirected to file with arguments
        to the command line
        """
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=STATIC_FILES["conv_bad"],
            verbose=0,
            criteria="strict",
            checker_names=["acdd", "cf"],
            output_filename=self.path,
            output_format="text",
        )

        assert os.stat(self.path).st_size > 0

    def test_multi_checker_return_value(self):
        """
        Tests that any failure for multiple checkers results in a failure return
        status
        """
        # CF should pass here
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=STATIC_FILES["ncei_gold_point_1"],
            verbose=0,
            criteria="strict",
            checker_names=["cf:1.6"],
            output_filename=self.path,
            output_format="text",
        )
        self.assertTrue(return_value)

        # CF should pass, but ACDD will have some errors.  Overall return
        # status should be a failure
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=STATIC_FILES["ncei_gold_point_1"],
            verbose=0,
            criteria="strict",
            checker_names=["acdd", "cf"],
            output_filename=self.path,
            output_format="text",
        )
        self.assertFalse(return_value)
