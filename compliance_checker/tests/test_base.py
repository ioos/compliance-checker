#!/usr/bin/env python
"""Tests for base compliance checker class"""

import os

from netCDF4 import Dataset

from compliance_checker import base


class TestBase:
    """
    Tests functionality of the base compliance checker class
    """

    def setup_method(self):
        self.acdd = base.BaseCheck()
        self.ds = Dataset(filename=os.devnull, mode="w", diskless=True)

    def teardown_method(self):
        self.ds.close()

    def test_attr_presence(self):
        # attribute not present, should fail
        priority = base.BaseCheck.MEDIUM
        rv1, rv2, rv3, rv4 = [], [], [], []
        attr = ("test", None)
        base.attr_check(attr, self.ds, priority, rv1)
        assert rv1[0] == base.Result(priority, False, "test", ["test not present"])
        # test with empty string
        self.ds.test = ""
        base.attr_check(attr, self.ds, priority, rv2)
        assert rv2[0] == base.Result(
            priority,
            False,
            "test",
            ["test is empty or completely whitespace"],
        )
        # test with whitespace in the form of a space and a tab
        self.ds.test = " 	"
        base.attr_check(attr, self.ds, priority, rv3)
        assert rv3[0] == base.Result(
            priority,
            False,
            "test",
            ["test is empty or completely whitespace"],
        )
        # test with actual string contents
        self.ds.test = "abc 123"
        base.attr_check(attr, self.ds, priority, rv4)
        assert rv4[0] == base.Result(priority, True, "test", [])

    def test_attr_in_valid_choices(self):
        """Tests attribute membership in a set"""
        rv1, rv2, rv3 = [], [], []
        priority = base.BaseCheck.MEDIUM
        valid_choices = ["a", "b", "c"]
        attr = ("test", valid_choices)
        base.attr_check(attr, self.ds, priority, rv1)
        assert rv1[0] == base.Result(priority, (0, 2), "test", ["test not present"])
        self.ds.test = ""
        base.attr_check(attr, self.ds, priority, rv2)
        assert rv2[0] == base.Result(
            priority,
            (1, 2),
            "test",
            [f"test present, but not in expected value list ({valid_choices})"],
        )
        self.ds.test = "a"
        base.attr_check(attr, self.ds, priority, rv3)
        assert rv3[0] == base.Result(priority, (2, 2), "test", [])

    def test_attr_fn(self):
        """Test attribute against a checker function"""

        # simple test.  In an actual program, this use case would be covered
        rv1, rv2, rv3 = [], [], []
        priority = base.BaseCheck.MEDIUM

        def verify_dummy(ds):
            """Sample function that will be called when passed into attr_check"""
            try:
                if ds.dummy + "y" == "dummyy":
                    return base.ratable_result(True, "dummy", [])
                else:
                    return base.ratable_result(False, "dummy", [ds.dummy + "y"])
            except AttributeError:
                return base.ratable_result(False, "dummy", [])

        attr = ("dummy", verify_dummy)
        base.attr_check(attr, self.ds, priority, rv1)
        assert rv1[0] == base.Result(priority, False, "dummy", [])
        self.ds.dummy = "doomy"
        base.attr_check(attr, self.ds, priority, rv2)
        assert rv2[0] == base.Result(priority, False, "dummy", ["doomyy"])
        self.ds.dummy = "dummy"
        base.attr_check(attr, self.ds, priority, rv3)
        assert rv3[0] == base.Result(priority, True, "dummy", [])

    def test_get_test_ctx(self):
        # acdd refers to a BaseCheck instance here -- perhaps the variable name
        # should reflect that?
        ctx = self.acdd.get_test_ctx(base.BaseCheck.HIGH, "Dummy Name")
        ctx.assert_true(1 + 1 == 2, "One plus one equals two")
        assert ctx.out_of == 1
        assert ctx.messages == []

        # ctx2 should be receive the same test context
        ctx2 = self.acdd.get_test_ctx(base.BaseCheck.HIGH, "Dummy Name")
        assert ctx is ctx2
        # will fail, obviously
        ctx2.assert_true(1 + 1 == 3, "One plus one equals three")
        assert ctx.out_of == 2
        assert ctx2.out_of == 2
        assert ctx2.messages == ["One plus one equals three"]

        ctx2 = self.acdd.get_test_ctx(base.BaseCheck.HIGH, "Test Name", "test_var_name")
        ctx3 = self.acdd.get_test_ctx(base.BaseCheck.HIGH, "Test Name", "test_var_name")
        # check that variable cache is working
        assert ctx3 is (
            self.acdd._defined_results["Test Name"]["test_var_name"][
                base.BaseCheck.HIGH
            ]
        )

    def test_email_validation(self):
        test_attr_name = "test"
        validator = base.EmailValidator()
        assert validator.validate(test_attr_name, "foo@bar.com")[0]
        bad_result = validator.validate(test_attr_name, "foo@@bar.com")
        assert not bad_result[0]
        assert bad_result[1] == ["test must be a valid email address"]

    def test_url_validation(self):
        """
        Test that URL validation works properly
        """
        test_attr_name = "test"
        # invalid URL
        test_url = "ssh://invalid_url"
        validator = base.UrlValidator()
        bad_result = validator.validate(test_attr_name, test_url)
        assert not bad_result[0]
        assert bad_result[1] == ["test must be a valid URL"]
        # valid URL
        test_url = "https://ioos.us"
        assert validator.validate(test_attr_name, test_url)[0]
        # test with CSV splitting rules, including checks with embedded commas,
        # which can appear in parts of URLs
        validator = base.UrlValidator(base.csv_splitter)
        url_multi_string = '"http://some-scientific-site.com/depth,temp",https://ioos.us/,http://google.com'
        assert validator.validate(test_attr_name, url_multi_string)[0]
        # add something that's invalid as a URL and check
        url_multi_string += ",noaa.ioos.webmaster@noaa.gov"
        bad_result = validator.validate(test_attr_name, url_multi_string)
        assert not bad_result[0]
        assert bad_result[1] == ["test must be a valid URL"]


class TestGenericFile:
    """
    Tests the GenericFile class.
    """

    def test_create_generic_file_success(self):
        path = "/tmp/test.txt"
        gf = base.GenericFile(path)
        assert gf.filepath() == path

    def test_create_generic_file_failure(self):
        gf = base.GenericFile("will not match")
        assert gf.filepath() != "do NOT MATCH"
