#!/usr/bin/env python
"""
Tests for command line output and parsing

"""
import io
import json
import os
import sys
from argparse import Namespace

import pytest

from compliance_checker.runner import CheckSuite, ComplianceChecker

from .conftest import static_files


@pytest.mark.usefixtures("checksuite_setup")
class TestCLI:
    """
    Tests various functions and aspects of the command line tool and runner
    """

    def shortDescription(self):
        return None

    def test_unicode_acdd_html(self, tmp_txt_file):
        """
        Tests that the checker is capable of producing HTML with unicode characters
        """
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=static_files("2dim"),
            verbose=0,
            criteria="strict",
            checker_names=["acdd"],
            output_filename=tmp_txt_file,
            output_format="html",
        )

        assert os.stat(tmp_txt_file).st_size > 0

    def test_unicode_cf_html(self, tmp_txt_file):
        """
        Tests that the CF checker can produce HTML output with unicode characters
        """
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=static_files("2dim"),
            verbose=0,
            criteria="strict",
            checker_names=["cf"],
            output_filename=tmp_txt_file,
            output_format="html",
        )

        assert os.stat(tmp_txt_file).st_size > 0

    def test_single_json_output(self, tmp_txt_file):
        """
        Tests that a suite can produce JSON output to a file
        """
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=static_files("conv_bad"),
            verbose=0,
            criteria="strict",
            checker_names=["cf"],
            output_filename=tmp_txt_file,
            output_format="json",
        )

        assert os.stat(tmp_txt_file).st_size > 0
        with open(tmp_txt_file) as f:
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

        mock_checkers = [Namespace(load=checker_1), Namespace(load=checker_2)]
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

    def test_multiple_json_output(self, tmp_txt_file):
        """
        Tests that a suite can produce JSON output to a file
        """
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=static_files("conv_bad"),
            verbose=0,
            criteria="strict",
            checker_names=["acdd", "cf"],
            output_filename=tmp_txt_file,
            output_format="json",
        )

        assert os.stat(tmp_txt_file).st_size > 0
        with open(tmp_txt_file) as f:
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
                ds_loc=static_files("conv_bad"),
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
                ds_loc=static_files("conv_bad"),
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

    def test_text_output(self, tmp_txt_file):
        """
        Tests that the 'text' output can be redirected to file with arguments
        to the command line
        """
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=static_files("conv_bad"),
            verbose=0,
            criteria="strict",
            checker_names=["acdd", "cf"],
            output_filename=tmp_txt_file,
            output_format="text",
        )

        assert os.stat(tmp_txt_file).st_size > 0

    def test_multi_checker_return_value(self, tmp_txt_file):
        """
        Tests that any failure for multiple checkers results in a failure return
        status
        """
        # CF should pass here
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=static_files("ncei_gold_point_1"),
            verbose=0,
            criteria="strict",
            checker_names=["cf:1.6"],
            output_filename=tmp_txt_file,
            output_format="text",
        )
        assert return_value

        # CF should pass, but ACDD will have some errors.  Overall return
        # status should be a failure
        return_value, errors = ComplianceChecker.run_checker(
            ds_loc=static_files("ncei_gold_point_1"),
            verbose=0,
            criteria="strict",
            checker_names=["acdd", "cf"],
            output_filename=tmp_txt_file,
            output_format="text",
        )
        assert not return_value
