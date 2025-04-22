import os

import numpy as np
from netCDF4 import Dataset

from compliance_checker.acdd import ACDD1_1Check, ACDD1_3Check
from compliance_checker.tests import BaseTestCase
from compliance_checker.tests.helpers import MockTimeSeries
from compliance_checker.tests.resources import STATIC_FILES


def to_singleton_var(item):
    """
    Get the first value of a list if this implements iterator protocol and is
    not a string
    """
    return [
        x[0] if hasattr(x, "__iter__") and not isinstance(x, str) else x for x in item
    ]


def check_varset_nonintersect(group0, group1):
    """
    Returns true if both groups contain the same elements, regardless of
    order.

    :param list group0: A list of strings to compare
    :param list group1: A list of strings to compare
    """
    # Performs symmetric difference on two lists converted to sets
    return len(set(group0) ^ set(group1)) == 0


class TestACDD1_1(BaseTestCase):
    # TODO superclass this so ACDD1_3 can inherit?
    # Adapted using `pandas.read_html` from URL
    # http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery_1-1
    expected = {
        "Highly Recommended": ["title", "summary", "keywords"],
        "Highly Recommended Variable Attributes": [
            "long_name",
            "standard_name",
            "units",
            "coverage_content_type",
        ],
        "Recommended": [
            "id",
            "naming_authority",
            "keywords_vocabulary",
            "history",
            "comment",
            "date_created",
            "creator_name",
            "creator_url",
            "creator_email",
            "institution",
            "project",
            "processing_level",
            "geospatial_bounds",
            "geospatial_lat_min",
            "geospatial_lat_max",
            "geospatial_lon_min",
            "geospatial_lon_max",
            "geospatial_vertical_min",
            "geospatial_vertical_max",
            "time_coverage_start",
            "time_coverage_end",
            "time_coverage_duration",
            "time_coverage_resolution",
            "standard_name_vocabulary",
            "license",
        ],
        "Suggested": [
            "contributor_name",
            "contributor_role",
            "publisher_name",
            "publisher_url",
            "publisher_email",
            "date_modified",
            "date_issued",
            "geospatial_lat_units",
            "geospatial_lat_resolution",
            "geospatial_lon_units",
            "geospatial_lon_resolution",
            "geospatial_vertical_units",
            "geospatial_vertical_resolution",
            "geospatial_vertical_positive",
        ],
    }

    def setup_method(self):
        # Use the NCEI Gold Standard Point dataset for ACDD checks
        self.ds = self.load_dataset(STATIC_FILES["ncei_gold_point_1"])
        self.acdd = ACDD1_1Check()
        self.acdd_highly_recommended = to_singleton_var(self.acdd.high_rec_atts)
        self.acdd_recommended = to_singleton_var(self.acdd.rec_atts)
        self.acdd_suggested = to_singleton_var(self.acdd.sug_atts)

    def test_cc_meta(self):
        assert self.acdd._cc_spec == "acdd"
        assert self.acdd._cc_spec_version == "1.1"

    def test_highly_recommended(self):
        """
        Checks that all highly recommended attributes are present
        """
        assert check_varset_nonintersect(
            self.expected["Highly Recommended"],
            self.acdd_highly_recommended,
        )

        # Check the reference dataset, NCEI 1.1 Gold Standard Point
        missing = ["\"Conventions\" does not contain 'ACDD-1.3'"]
        results = self.acdd.check_high(self.ds)
        for result in results:
            if result.msgs and all(m in missing for m in result.msgs):
                # only the Conventions check should have failed
                self.assert_result_is_bad(result)
            self.assert_result_is_good(result)

        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        # TODO: DRY these dataset open/close, possibly with pytest features
        empty_ds = Dataset(os.devnull, "w", memory=True, diskless=True)
        results = self.acdd.check_high(empty_ds)
        for result in results:
            self.assert_result_is_bad(result)
        empty_ds.close()

    def test_recommended(self):
        """
        Checks that all recommended attributes are present
        """
        # 'geospatial_bounds' attribute currently has its own separate check
        # from the list of required atts
        assert check_varset_nonintersect(
            self.expected["Recommended"],
            self.acdd_recommended,
        )

        ncei_exceptions = [
            "geospatial_bounds not present",
            "time_coverage_duration not present",
            "time_coverage_resolution not present",
        ]
        results = self.acdd.check_recommended(self.ds)
        for result in results:
            if (result.msgs) and all(
                m in ncei_exceptions for m in result.msgs
            ):  # we're doing string comparisons, this is kind of hacky...
                self.assert_result_is_bad(result)
                continue

            self.assert_result_is_good(result)

        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        empty_ds = Dataset(os.devnull, "w", diskless=True, memory=True)

        results = self.acdd.check_recommended(empty_ds)
        for result in results:
            self.assert_result_is_bad(result)
        empty_ds.close()

    def test_suggested(self):
        """
        Checks that all suggested attributes are present
        """
        assert check_varset_nonintersect(
            self.expected["Suggested"],
            self.acdd_suggested,
        )

        # Attributes that are missing from NCEI but should be there
        missing = [
            "geospatial_lat_resolution not present",
            "geospatial_lon_resolution not present",
            "geospatial_vertical_resolution not present",
        ]

        results = self.acdd.check_suggested(self.ds)
        for result in results:
            if (result.msgs) and all(
                m in missing for m in result.msgs
            ):  # we're doing string comparisons, this is kind of hacky...
                self.assert_result_is_bad(result)
                continue

            self.assert_result_is_good(result)

        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        empty_ds = Dataset(os.devnull, "w", diskless=True, memory=True)

        results = self.acdd.check_recommended(empty_ds)
        for result in results:
            self.assert_result_is_bad(result)
        empty_ds.close()

    def test_acknowldegement_check(self):
        """Check both the British and American English spellings of 'acknowledgement'"""
        # Check British Spelling
        try:
            empty0 = Dataset(os.devnull, "w", diskless=True, memory=True)
            result = self.acdd.check_acknowledgment(empty0)
            self.assert_result_is_bad(result)

            empty0.acknowledgement = "Attribution goes here"
            result = self.acdd.check_acknowledgment(empty0)
            self.assert_result_is_good(result)
        finally:
            empty0.close()

        try:
            # Check American spelling
            empty1 = Dataset(os.devnull, "w", diskless=True, memory=True)
            result = self.acdd.check_acknowledgment(empty1)
            self.assert_result_is_bad(result)

            empty1.acknowledgment = "Attribution goes here"
            result = self.acdd.check_acknowledgment(empty1)
            self.assert_result_is_good(result)
        finally:
            empty1.close()


class TestACDD1_3(BaseTestCase):
    # Adapted using `pandas.read_html` from URL
    # http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery_1-3
    expected = {
        "Suggested": [
            "creator_type",
            "creator_institution",
            "publisher_type",
            "publisher_institution",
            "program",
            "contributor_name",
            "contributor_role",
            "geospatial_lat_units",
            "geospatial_lat_resolution",
            "geospatial_lon_units",
            "geospatial_lon_resolution",
            "geospatial_vertical_units",
            "geospatial_vertical_resolution",
            "date_modified",
            "date_issued",
            "date_metadata_modified",
            "product_version",
            "keywords_vocabulary",
            "platform",
            "platform_vocabulary",
            "instrument",
            "instrument_vocabulary",
            "metadata_link",
            "references",
        ],
        "Highly Recommended": ["title", "summary", "keywords", "Conventions"],
        "Recommended": [
            "id",
            "naming_authority",
            "history",
            "source",
            "processing_level",
            "comment",
            "license",
            "standard_name_vocabulary",
            "date_created",
            "creator_name",
            "creator_email",
            "creator_url",
            "institution",
            "project",
            "publisher_name",
            "publisher_email",
            "publisher_url",
            "geospatial_bounds",
            "geospatial_bounds_crs",
            "geospatial_bounds_vertical_crs",
            "geospatial_lat_min",
            "geospatial_lat_max",
            "geospatial_lon_min",
            "geospatial_lon_max",
            "geospatial_vertical_min",
            "geospatial_vertical_max",
            "geospatial_vertical_positive",
            "time_coverage_start",
            "time_coverage_end",
            "time_coverage_duration",
            "time_coverage_resolution",
        ],
        "Highly Recommended Variable Attributes": [
            "long_name",
            "standard_name",
            "units",
            "coverage_content_type",
        ],
    }

    def setup_method(self):
        # Use the NCEI Gold Standard Point dataset for ACDD checks
        self.ds = self.load_dataset(STATIC_FILES["ncei_gold_point_2"])
        self.acdd = ACDD1_3Check()
        self.acdd_highly_recommended = to_singleton_var(self.acdd.high_rec_atts)
        self.acdd_recommended = to_singleton_var(self.acdd.rec_atts)
        self.acdd_suggested = to_singleton_var(self.acdd.sug_atts)

    def test_cc_meta(self):
        assert self.acdd._cc_spec == "acdd"
        assert self.acdd._cc_spec_version == "1.3"

    def test_highly_recommended(self):
        """
        Checks that all highly recommended attributes are present
        """
        assert check_varset_nonintersect(
            self.expected["Highly Recommended"],
            self.acdd_highly_recommended,
        )

        results = self.acdd.check_high(self.ds)
        for result in results:
            # NODC 2.0 has a different value in the conventions field
            self.assert_result_is_good(result)

    def test_recommended(self):
        """
        Checks that all recommended attributes are present
        """
        assert check_varset_nonintersect(
            self.expected["Recommended"],
            self.acdd_recommended,
        )

        results = self.acdd.check_recommended(self.ds)
        ncei_exceptions = [
            "time_coverage_duration not present",
            "time_coverage_resolution not present",
        ]
        for result in results:
            if (result.msgs) and all(
                m in ncei_exceptions for m in result.msgs
            ):  # we're doing string comparisons, this is kind of hacky...
                self.assert_result_is_bad(result)
                continue

            self.assert_result_is_good(result)

    def test_suggested(self):
        """
        Checks that all suggested attributes are present
        """
        assert check_varset_nonintersect(
            self.expected["Suggested"],
            self.acdd_suggested,
        )

        results = self.acdd.check_suggested(self.ds)
        # NCEI does not require or suggest resolution attributes
        ncei_exceptions = [
            "geospatial_lat_resolution not present",
            "geospatial_lon_resolution not present",
            "geospatial_vertical_resolution not present",
        ]
        for result in results:
            if (result.msgs) and all(
                m in ncei_exceptions for m in result.msgs
            ):  # we're doing string comparisons, this is kind of hacky...
                self.assert_result_is_bad(result)
                continue
            self.assert_result_is_good(result)

    def test_variables(self):
        """
        Test that variables are checked for required attributes
        """
        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        empty_ds = Dataset(os.devnull, "w", diskless=True, memory=True)
        # self.addCleanup(empty_ds.close)

        # The dataset needs at least one variable to check that it's missing
        # all the required attributes.
        empty_ds.createDimension("time", 1)
        empty_ds.createVariable("fake", "float32", ("time",))

        # long_name
        results = self.acdd.check_var_long_name(self.ds)
        for result in results:
            self.assert_result_is_good(result)

        results = self.acdd.check_var_long_name(empty_ds)
        assert len(results) == 1
        for result in results:
            self.assert_result_is_bad(result)

        # standard_name
        results = self.acdd.check_var_standard_name(self.ds)
        for result in results:
            self.assert_result_is_good(result)

        results = self.acdd.check_var_standard_name(empty_ds)
        assert len(results) == 1
        for result in results:
            self.assert_result_is_bad(result)

        # units
        results = self.acdd.check_var_units(self.ds)
        for result in results:
            self.assert_result_is_good(result)

        results = self.acdd.check_var_units(empty_ds)
        assert len(results) == 1
        for result in results:
            self.assert_result_is_bad(result)

        empty_ds.close()

    def test_vertical_extents(self):
        """
        Test vertical extents are being checked
        """

        result = self.acdd.check_vertical_extents(self.ds)
        self.assert_result_is_good(result)

    def test_geospatial_bounds(self):
        """
        Test geospatial bounds are checked and provide a good error message
        """
        # Create an empty dataset that writes to /dev/null This acts as a
        # temporary netCDF file in-memory that never gets written to disk.
        empty_ds = Dataset(os.devnull, "w", diskless=True)

        # Misspelled WKT. Error message should include the attribute checked
        # and the value that was provided for easy troubleshooting
        empty_ds.geospatial_bounds = "POIT (-123.458000 38.048000)"
        results = self.acdd.check_recommended(empty_ds)
        for result in results:
            if result.variable_name == "geospatial_bounds":
                assert (
                    "Could not parse WKT from geospatial_bounds,"
                    f' possible bad value: "{empty_ds.geospatial_bounds}"'
                    in result.msgs
                )
        empty_ds.close()

    def test_time_extents(self):
        """
        Test that the time extents are being checked
        """
        result = self.acdd.check_time_extents(self.ds)
        self.assert_result_is_good(result)

        empty_ds = Dataset(os.devnull, "w", diskless=True)

        # The dataset needs at least one variable to check that it's missing
        # all the required attributes.
        empty_ds.createDimension("time", 1)
        time_var = empty_ds.createVariable("time", "float32", ("time",))
        time_var.units = "seconds since 1970-01-01 00:00:00 UTC"
        time_var[:] = [1451692800]  # 20160102T000000Z in seconds since epoch
        empty_ds.time_coverage_start = "20160102T000000Z"
        empty_ds.time_coverage_end = "20160102T000000Z"
        result = self.acdd.check_time_extents(empty_ds)
        self.assert_result_is_good(result)
        # try the same thing with time offsets
        time_var.units = "seconds since 1970-01-01 00:00:00-10:00"
        empty_ds.time_coverage_start = "20160102T000000-1000"
        empty_ds.time_coverage_end = "20160102T000000-1000"
        result = self.acdd.check_time_extents(empty_ds)
        self.assert_result_is_good(result)
        empty_ds.close()

    def test_check_lat_extents(self):
        """Test the check_lat_extents() method behaves expectedly"""

        # create dataset using MockDataset, give it lat/lon dimensions
        ds = MockTimeSeries()
        ds.variables["lat"][:] = np.linspace(
            -135.0,
            -130.0,
            num=500,
        )  # arbitrary, but matches time dim size

        # test no values, expect failure
        result = self.acdd.check_lat_extents(ds)
        self.assert_result_is_bad(result)

        # give integer geospatial_lat_max/min, test
        ds.setncattr("geospatial_lat_min", -135)
        ds.setncattr("geospatial_lat_max", -130)

        result = self.acdd.check_lat_extents(ds)
        self.assert_result_is_good(result)

        # give float geospatial_lat_min/max, test
        ds.setncattr("geospatial_lat_min", -135.0)
        ds.setncattr("geospatial_lat_max", -130.0)

        result = self.acdd.check_lat_extents(ds)
        self.assert_result_is_good(result)

        # give string (in number-form), test
        ds.setncattr("geospatial_lat_min", "-135.")
        ds.setncattr("geospatial_lat_max", "-130.")

        result = self.acdd.check_lat_extents(ds)
        self.assert_result_is_good(result)

        # give garbage string -- expect failure
        ds.setncattr("geospatial_lat_min", "bad")
        ds.setncattr("geospatial_lat_max", "val")

        result = self.acdd.check_lat_extents(ds)
        self.assert_result_is_bad(result)

    def test_check_lon_extents(self):
        """Test the check_lon_extents() method behaves expectedly"""

        # create dataset using MockDataset, give it lat/lon dimensions
        ds = MockTimeSeries()
        ds.variables["lon"][:] = np.linspace(65.0, 67.0, num=500)

        # test no values, expect failure
        # result = self.acdd.check_lon_extents(ds)
        # self.assert_result_is_bad(result)

        # give integer geospatial_lon_max/min, test
        ds.setncattr("geospatial_lon_min", 65)
        ds.setncattr("geospatial_lon_max", 67)

        result = self.acdd.check_lon_extents(ds)
        self.assert_result_is_good(result)

        # give float geospatial_lon_min/max, test
        ds.setncattr("geospatial_lon_min", 65.0)
        ds.setncattr("geospatial_lon_max", 67.0)

        result = self.acdd.check_lon_extents(ds)
        self.assert_result_is_good(result)

        # give string (in number-form), test
        ds.setncattr("geospatial_lon_min", "65.")
        ds.setncattr("geospatial_lon_max", "67.")

        result = self.acdd.check_lon_extents(ds)
        self.assert_result_is_good(result)

        # give garbage string -- expect failure
        ds.setncattr("geospatial_lon_min", "bad")
        ds.setncattr("geospatial_lon_max", "val")

        result = self.acdd.check_lon_extents(ds)
        self.assert_result_is_bad(result)

    def test_check_geospatial_vertical_max(self):
        ds = MockTimeSeries()
        ds.variables["depth"][:] = np.linspace(0.0, 30.0, num=500)

        # give integer geospatial_vertical_max/min, test
        ds.setncattr("geospatial_vertical_min", 0)
        ds.setncattr("geospatial_vertical_max", 30)

        result = self.acdd.check_vertical_extents(ds)
        self.assert_result_is_good(result)

        # give float geospatial_vertical_min/max, test
        ds.setncattr("geospatial_vertical_min", 0.0)
        ds.setncattr("geospatial_vertical_max", 30.0)

        result = self.acdd.check_vertical_extents(ds)
        self.assert_result_is_good(result)

        # give string (in number-form), test
        ds.setncattr("geospatial_vertical_min", "0.")
        ds.setncattr("geospatial_vertical_max", "30.")

        result = self.acdd.check_vertical_extents(ds)
        self.assert_result_is_good(result)

        # give garbage string -- expect failure
        ds.setncattr("geospatial_vertical_min", "bad")
        ds.setncattr("geospatial_vertical_max", "val")

        result = self.acdd.check_vertical_extents(ds)
        self.assert_result_is_bad(result)

        # all masked values mean that there are no valid array elements to get
        # the min/max of
        ds.setncattr("geospatial_vertical_min", 0.0)
        ds.setncattr("geospatial_vertical_max", 30.0)
        ds.variables["depth"][:] = np.ma.masked_all(ds.variables["depth"].shape)
        result = self.acdd.check_vertical_extents(ds)
        self.assert_result_is_bad(result)
