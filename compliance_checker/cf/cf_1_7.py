import logging
import os
import sqlite3
from warnings import warn

import numpy as np
import pyproj

import compliance_checker.cf.util as cfutil
from compliance_checker.base import BaseCheck, Result, TestCtx
from compliance_checker.cf.appendix_d import dimless_vertical_coordinates_1_7
from compliance_checker.cf.appendix_e import cell_methods17
from compliance_checker.cf.appendix_f import (
    ellipsoid_names17,
    grid_mapping_attr_types17,
    grid_mapping_dict17,
    horizontal_datum_names17,
    prime_meridian_names17,
)
from compliance_checker.cf.cf_1_6 import CF1_6Check
from compliance_checker.cf.cf_base import appendix_a_base

logger = logging.getLogger(__name__)


class CF1_7Check(CF1_6Check):
    """Implementation for CF v1.7. Inherits from CF1_6Check as most of the
    checks are the same."""

    # things that are specific to 1.7
    _cc_spec_version = "1.7"
    _cc_url = "http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html"

    appendix_a = appendix_a_base.copy()
    appendix_a.update(
        {
            "actual_range": {
                "Type": "N",
                "attr_loc": {"D", "C"},
                "cf_section": "2.5.1",
            },
            "comment": {
                "Type": "S",
                "attr_loc": {"G", "D", "C"},
                "cf_section": "2.6.2",
            },
            "external_variables": {
                "Type": "S",
                "attr_loc": {"G"},
                "cf_section": "2.6.3",
            },
            "scale_factor": {"Type": "N", "attr_loc": {"D", "C"}, "cf_section": "8.1"},
        },
    )

    def __init__(self, options=None):
        super().__init__(options)

        self.cell_methods = cell_methods17
        self.grid_mapping_dict = grid_mapping_dict17
        self.grid_mapping_attr_types = grid_mapping_attr_types17

    def check_external_variables(self, ds):
        """
        The global external_variables attribute is a blank-separated list of the
        names of variables which are named by attributes in the file but which
        are not present in the file. These variables are to be found in other
        files (called "external files") but CF does not provide conventions for
        identifying the files concerned. The only attribute for which CF
        standardises the use of external variables is cell_measures.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """
        external_vars_ctx = TestCtx(BaseCheck.MEDIUM, self.section_titles["2.6.3"])
        # IMPLEMENTATION CONFORMANCE 2.6.3 REQUIRED 2/2
        try:
            # Should external variables be stored in dataset state?
            external_var_names = set(ds.external_variables.strip().split())

            bad_external_var_names = external_var_names.intersection(ds.variables)

            if bad_external_var_names:
                external_vars_ctx.out_of += 1
                bad_msg = (
                    "Global attribute external_variables should not "
                    "have any variable names which are present in the dataset. "
                    "Currently, the following names appear in both external_variables "
                    f"and the dataset's variables: {bad_external_var_names}"
                )
                external_vars_ctx.messages.append(bad_msg)

        # string/global attributes are handled in Appendix A checks
        except (AttributeError, ValueError):
            pass

        return external_vars_ctx.to_result()

    def check_actual_range(self, ds):
        """
        Check the actual_range attribute of variables. As stated in
        section 2.5.1 of version 1.7, this convention defines a two-element
        vector attribute designed to describe the actual minimum and actual
        maximum values of variables containing numeric data. Conditions:
          - the fist value of the two-element vector must be equal to the
            minimum of the data, and the second element equal to the maximum
          - if the data is packed, the elements of actual_range should have
            the same data type as the *unpacked* data
          - if valid_range is specified, both elements of actual_range should
            be within valid_range

        If a variable does not have an actual_range attribute, let it pass;
        including this attribute is only suggested. However, if the user is
        specifying the actual_range, the Result will be considered
        high-priority."""

        ret_val = []

        for name, variable in ds.variables.items():
            msgs = []
            score = 0
            out_of = 0

            if not hasattr(variable, "actual_range"):
                continue  # having this attr is only suggested, no Result needed
            else:
                out_of += 1
                try:
                    if (
                        len(variable.actual_range) != 2
                    ):  # TODO is the attr also a numpy array? if so, .size
                        msgs.append(
                            f"actual_range of '{name}' must be 2 elements",
                        )
                        ret_val.append(
                            Result(  # putting result into list
                                BaseCheck.HIGH,
                                (score, out_of),
                                self.section_titles["2.5"],
                                msgs,
                            ),
                        )
                        continue  # no need to keep checking if already completely wrong
                    else:
                        score += 1
                except TypeError:  # in case it's just a single number
                    msgs.append(f"actual_range of '{name}' must be 2 elements")
                    ret_val.append(
                        Result(  # putting result into list
                            BaseCheck.HIGH,
                            (score, out_of),
                            self.section_titles["2.5"],
                            msgs,
                        ),
                    )
                    continue

                # check equality to existing min/max values
                # NOTE this is a data check
                # If every value is masked, a data check of actual_range isn't
                # appropriate, so skip.
                if not (hasattr(variable[:], "mask") and variable[:].mask.all()):
                    # if min/max values aren't close to actual_range bounds,
                    # fail.
                    out_of += 1
                    if not np.isclose(
                        variable.actual_range[0],
                        variable[:].min(),
                    ) or not np.isclose(variable.actual_range[1], variable[:].max()):
                        msgs.append(
                            f"actual_range elements of '{name}' inconsistent with its min/max values",
                        )
                    else:
                        score += 1

                # check that the actual range is within the valid range
                if hasattr(variable, "valid_range"):  # check within valid_range
                    out_of += 1
                    if (variable.actual_range[0] < variable.valid_range[0]) or (
                        variable.actual_range[1] > variable.valid_range[1]
                    ):
                        msgs.append(
                            f'"{name}"\'s actual_range must be within valid_range',
                        )
                    else:
                        score += 1

                # check the elements of the actual range have the appropriate
                # relationship to the valid_min and valid_max
                if hasattr(variable, "valid_min"):
                    out_of += 1
                    if variable.actual_range[0] < variable.valid_min:
                        msgs.append(
                            f'"{name}"\'s actual_range first element must be >= valid_min ({variable.valid_min})',
                        )
                    else:
                        score += 1
                if hasattr(variable, "valid_max"):
                    out_of += 1
                    if variable.actual_range[1] > variable.valid_max:
                        msgs.append(
                            f'"{name}"\'s actual_range second element must be <= valid_max ({variable.valid_max})',
                        )
                    else:
                        score += 1

            ret_val.append(
                Result(  # putting result into list
                    BaseCheck.HIGH,
                    (score, out_of),
                    self.section_titles["2.5"],
                    msgs,
                ),
            )
        return ret_val

    def check_cell_boundaries(self, ds):
        """
        Checks the dimensions of cell boundary variables to ensure they are CF compliant
        per section 7.1.

        This method extends the CF1_6Check method; please see the original method for the
        complete doc string.

        If any variable contains both a formula_terms attribute *and* a bounding variable,
        that bounds variable must also have a formula_terms attribute.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :returns list: List of results
        """

        # Note that test does not check monotonicity
        ret_val = []
        reasoning = []
        for variable_name, boundary_variable_name in cfutil.get_cell_boundary_map(
            ds,
        ).items():
            variable = ds.variables[variable_name]
            valid = True
            reasoning = []

            # 7.1 Required 1/5:
            # The type of the bounds attribute is a string whose value is a single variable name.
            # The specified variable must exist in the file.
            if boundary_variable_name not in ds.variables:
                valid = False
                reasoning.append(
                    f"Boundary variable {boundary_variable_name} referenced by {variable.name} not "
                    + "found in dataset variables",
                )
            else:
                boundary_variable = ds.variables[boundary_variable_name]

            # 7.1 Required 2/5:
            # The number of dimensions in the bounds variable should always be
            # the number of dimensions in the referring variable + 1
            if boundary_variable.ndim < 2:
                valid = False
                reasoning.append(
                    f"Boundary variable {boundary_variable.name} specified by {variable.name}"
                    + " should have at least two dimensions to enclose the base "
                    + "case of a one dimensionsal variable",
                )
            if boundary_variable.ndim != variable.ndim + 1:
                valid = False
                reasoning.append(
                    f"The number of dimensions of the variable {variable.name} is {variable.ndim}, but the "
                    f"number of dimensions of the boundary variable {boundary_variable.name} is {boundary_variable.ndim}. The boundary variable "
                    f"should have {variable.ndim + 1} dimensions",
                )
            if variable.dimensions[:] != boundary_variable.dimensions[: variable.ndim]:
                valid = False
                reasoning.append(
                    f"Boundary variable coordinates (for {variable.name}) are in improper order: {boundary_variable.dimensions}. Bounds-specific dimensions should be last"
                    "",
                )

            # 7.1 Required 2/5: continue
            # Ensure p vertices form a valid simplex given previous a...n
            # previous auxiliary coordinates
            if (
                ds.dimensions[boundary_variable.dimensions[-1]].size
                < len(boundary_variable.dimensions[:-1]) + 1
            ):
                valid = False
                reasoning.append(
                    f"Dimension {boundary_variable.name} of boundary variable (for {variable.name}) must have at least {len(variable.dimensions) + 1} elements to form a simplex/closed cell with previous dimensions {boundary_variable.dimensions[:-1]}.",
                )

            # 7.1 Required 3/5:
            # A boundary variable must be a numeric data type
            if boundary_variable.dtype.kind not in "biufc":
                valid = False
                reasoning.append(
                    f"Boundary variable {boundary_variable.name} specified by {variable.name}"
                    + "must be a numeric data type ",
                )

            # 7.1 Required 4/5:
            # If a boundary variable has units, standard_name, axis, positive, calendar, leap_month,
            # leap_year or month_lengths attributes, they must agree with those of its associated variable.
            if boundary_variable.__dict__.keys():
                for item in boundary_variable.__dict__.keys():
                    if hasattr(variable, item):
                        if getattr(variable, item) != getattr(boundary_variable, item):
                            valid = False
                            reasoning.append(
                                f"'{boundary_variable_name}' has attr '{item}' with value '{getattr(boundary_variable, item)}' that does not agree "
                                f"with its associated variable ('{variable.name}')'s attr value '{getattr(variable, item)}'"
                                "",
                            )

            # 7.1 Required 5/5:
            # check if formula_terms is present in the var; if so,
            # the bounds variable must also have a formula_terms attr
            if hasattr(variable, "formula_terms"):
                if not hasattr(boundary_variable, "formula_terms"):
                    valid = False
                    reasoning.append(
                        f"'{variable_name}' has 'formula_terms' attr, bounds variable '{boundary_variable_name}' must also have 'formula_terms'",
                    )

            # 7.1 Recommendations 2/2
            # Boundary variables should not have the _FillValue, missing_value, units, standard_name, axis,
            # positive, calendar, leap_month, leap_year or month_lengths attributes.
            attributes_to_check = {
                "_FillValue",
                "missing_value",
                "units",
                "standard_name",
                "axis",
                "positive",
                "calendar",
                "leap_month",
                "leap_year",
                "month_lengths",
            }
            if boundary_variable.__dict__.keys():
                lst1 = boundary_variable.__dict__.keys()
                lst2 = attributes_to_check
                unwanted_attributes = [value for value in lst1 if value in lst2]
                if unwanted_attributes:
                    valid = False
                    reasoning.append(
                        f"The Boundary variables '{boundary_variable_name}' should not have the attributes: '{unwanted_attributes}'",
                    )

            result = Result(
                BaseCheck.MEDIUM,
                valid,
                self.section_titles["7.1"],
                reasoning,
            )
            ret_val.append(result)
        return ret_val

    def check_cell_boundaries_interval(self, ds):
        """
        7.1 Cell Boundaries
        Recommendations: (1/2)
        The points specified by a coordinate or auxiliary coordinate variable
        should lie within, or on the boundary, of the cells specified by the
        associated boundary variable.
        """
        ret_val = []
        reasoning = []
        for variable_name, boundary_variable_name in cfutil.get_cell_boundary_map(
            ds,
        ).items():
            valid = True

            variable = ds.variables[variable_name]
            boundary_variable = ds.variables[boundary_variable_name]

            for ii in range(len(variable[:])):
                if abs(boundary_variable[ii][1]) >= abs(boundary_variable[ii][0]):
                    if not (
                        (abs(variable[ii]) >= abs(boundary_variable[ii][0]))
                        and (abs(variable[ii]) <= abs(boundary_variable[ii][1]))
                    ):
                        valid = False
                        reasoning.append(
                            f"The points specified by the coordinate variable {variable_name} ({variable[ii]})"
                            " lie outside the boundary of the cell specified by the "
                            f"associated boundary variable {boundary_variable_name} ({boundary_variable[ii]})",
                        )

                result = Result(
                    BaseCheck.MEDIUM,
                    valid,
                    self.section_titles["7.1"],
                    reasoning,
                )
                ret_val.append(result)
            return ret_val

    def check_cell_measures(self, ds):
        """
        A method to over-ride the CF1_6Check method. In CF 1.7, it is specified
        that variable referenced by cell_measures must be in the dataset OR
        referenced by the global attribute "external_variables", which represent
        all the variables used in the dataset but not found in the dataset.

        7.2 To indicate extra information about the spatial properties of a
        variable's grid cells, a cell_measures attribute may be defined for a
        variable. This is a string attribute comprising a list of
        blank-separated pairs of words of the form "measure: name". "area" and
        "volume" are the only defined measures.

        The "name" is the name of the variable containing the measure values,
        which we refer to as a "measure variable". The dimensions of the
        measure variable should be the same as or a subset of the dimensions of
        the variable to which they are related, but their order is not
        restricted.

        The variable must have a units attribute and may have other attributes
        such as a standard_name.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        variables = ds.get_variables_by_attributes(
            cell_measures=lambda c: c is not None,
        )
        try:
            external_variables_str = ds.getncattr("external_variables")
            if external_variables_str is not None:
                external_variables_names = set(external_variables_str.split(" "))
            else:
                external_variables_names = set()
        except (ValueError, AttributeError):
            external_variables_names = set()
        for var in variables:
            ret_val.append(
                self._cell_measures_core(
                    ds,
                    var,
                    external_variables_names,
                    "dataset or external variable",
                ),
            )

        return ret_val

    def _check_grid_mapping_attr_condition(self, attr, attr_name):
        """
        Evaluate a condition (or series of conditions) for a particular
        attribute. Implementation for CF-1.7.

        :param attr: attribute to teset condition for
        :param str attr_name: name of the attribute
        :rtype tuple
        :return two-tuple of (bool, str)
        """

        if attr_name == "geographic_crs_name":
            return self._evaluate_geographic_crs_name(attr)

        elif attr_name == "geoid_name":
            return self._evaluate_geoid_name(attr)

        elif attr_name == "geopotential_datum_name":
            return self._evaluate_geopotential_datum_name(attr)

        elif attr_name == "horizontal_datum_name":
            return self._evaluate_horizontal_datum_name(attr)

        elif attr_name == "prime_meridian_name":
            return self._evaluate_prime_meridian_name(attr)

        elif attr_name == "projected_crs_name":
            return self._evaluate_projected_crs_name(attr)

        elif attr_name == "reference_ellipsoid_name":
            return self._evaluate_reference_ellipsoid_name(attr)

        elif attr_name == "towgs84":
            return self._evaluate_towgs84(attr)

        else:  # invoke method from 1.6, as these names are all still valid
            return super()._check_grid_mapping_attr_condition(
                attr,
                attr_name,
            )

    def _check_gmattr_existence_condition_geoid_name_geoptl_datum_name(self, var):
        """
        Check to see if both geoid_name and geopotential_datum_name exist as attributes
        for `var`. They should not.

        :param netCDF4.Variable var
        :rtype tuple
        :return two-tuple (bool, str)
        """

        msg = "Both geoid_name and geopotential_datum_name cannot exist"

        if ("geoid_name" in var.ncattrs()) and (
            "geopotential_datum_name" in var.ncattrs()
        ):
            return (False, msg)

        else:
            return (True, msg)

    def _check_gmattr_existence_condition_ell_pmerid_hdatum(self, var):
        """
        If one of reference_ellipsoid_name, prime_meridian_name, or
        horizontal_datum_name are defined as grid_mapping attributes,
        they must all be defined.

        :param netCDF4.Variable var
        :rtype tuple
        :return two-tuple (bool, str)
        """

        msg = (
            "If any of reference_ellipsoid_name, prime_meridian_name, "
            "or horizontal_datum_name are defined, all must be defined."
        )

        _ncattrs = set(var.ncattrs())

        if any(
            x in _ncattrs
            for x in [
                "reference_ellipsoid_name",
                "prime_meridian_name",
                "horizontal_datum_name",
            ]
        ) and (
            not {
                "reference_ellipsoid_name",
                "prime_meridian_name",
                "horizontal_datum_name",
            }.issubset(_ncattrs)
        ):
            return (False, msg)

        else:
            return (True, msg)

    def _exec_query_str_with_params(self, qstr, argtuple):
        """
        Execute a query string in a database connection with the given argument
        tuple. Return a result set.

        :param str qstr: desired query to be executed
        :param tuple argtuple: tuple of arguments to be supplied to query
        :rtype set
        """

        proj_db_path = os.path.join(pyproj.datadir.get_data_dir(), "proj.db")
        with sqlite3.connect(proj_db_path) as conn:
            return conn.execute(qstr, argtuple)

    def _evaluate_geographic_crs_name(self, val):
        """
        Evaluate the condition for the geographic_crs_name attribute.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple of (bool, str)
        """

        query_str = (
            "SELECT 1 FROM geodetic_crs WHERE name = ? "
            "UNION ALL "  # need union in case contained in other tables
            "SELECT 1 FROM alias_name WHERE alt_name = ? "
            "AND table_name = 'geodetic_crs' LIMIT 1"
        )

        # try to find the value in the database
        res_set = self._exec_query_str_with_params(query_str, (val, val))

        # does it exist? if so, amt returned  be > 1
        return (
            len(res_set.fetchall()) > 0,
            "geographic_crs_name must correspond to a valid OGC WKT GEOGCS name",
        )

    def _evaluate_geoid_name(self, val):
        """
        Evaluate the condition for the geod_name attribute.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple of (bool, str)
        """

        query_str = (
            "SELECT 1 FROM vertical_datum WHERE name = ? "
            "UNION ALL "
            "SELECT 1 FROM alias_name WHERE alt_name = ? "
            "AND table_name = 'vertical_datum' LIMIT 1"
        )

        # try to find the value in the database
        res_set = self._exec_query_str_with_params(query_str, (val, val))

        return (
            len(res_set.fetchall()) > 0,
            "geoid_name must correspond to a valid OGC WKT VERT_DATUM name",
        )

    def _evaluate_geopotential_datum_name(self, val):
        """
        Evaluate the condition for the geogpotential_datum_name attribute.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple of (bool, str)
        """

        query_str = (
            "SELECT 1 FROM vertical_datum WHERE name = ? "
            "UNION ALL "
            "SELECT 1 FROM alias_name WHERE alt_name = ? "
            "AND table_name = 'vertical_datum' LIMIT 1"
        )

        # try to find the value in the database
        res_set = self._exec_query_str_with_params(query_str, (val, val))

        return (
            len(res_set.fetchall()) > 0,
            "geopotential_datum_name must correspond to a valid OGC WKT VERT_DATUM name",
        )

    def _evaluate_horizontal_datum_name(self, val):
        """
        Evaluate the condition for the horizontal_datum_name attribute.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple of (bool, str)
        """

        return (
            val in horizontal_datum_names17,
            (
                "{} must be a valid Horizontal Datum Name; "
                "see https://github.com/cf-convention/cf-conventions/wiki/Mapping-from-CF-Grid-Mapping-Attributes-to-CRS-WKT-Elements."
            ),
        )

    def _evaluate_prime_meridian_name(self, val):
        """
        Evaluate the condition for the prime_meridian_name.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple of (bool, str)
        """

        return (
            val in prime_meridian_names17,
            (
                "{} must be a valid Prime Meridian name; "
                "see https://github.com/cf-convention/cf-conventions/wiki/csv/prime_meridian.csv."
            ),
        )

    def _evaluate_projected_crs_name(self, val):
        """
        Evaluate the condition for the projected_crs attribute.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple of (bool, str)
        """

        query_str = (
            "SELECT 1 FROM projected_crs WHERE name = ? "
            "UNION ALL "
            "SELECT 1 FROM alias_name WHERE alt_name = ? "
            "AND table_name = 'projected_crs' LIMIT 1"
        )

        # try to find the value in the database
        res_set = self._exec_query_str_with_params(query_str, (val, val))

        return (
            len(res_set.fetchall()) > 0,
            "projected_crs_name must correspond to a valid OGC WKT PROJCS name",
        )

    def _evaluate_reference_ellipsoid_name(self, val):
        """
        Evaluate the condition for the reference_ellipsoid_name attribute.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple of (bool, str)
        """

        return (
            val in ellipsoid_names17,
            (
                "{} must be a valid Ellipsoid Name; "
                "see https://github.com/cf-convention/cf-conventions/wiki/csv/ellipsoid.csv."
            ),
        )

    def _evaluate_towgs84(self, val):
        """
        Evaluate the condition for the towgs84 attribute.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple of (bool, str)
        """

        msg = (
            "towgs84 must be an array of length 3, 6, or 7 of double-precision"
            " and correspond to anm OGC WKT TOWGS84 node"
        )

        # if not numpy type, return false
        if not getattr(val, "dtype", None):
            return (False, msg)

        # must be double-precision array
        elif val.dtype != np.float64:
            return (False, msg)

        # must be of length 3, 6, or 7
        elif not val.shape:  # single value
            return (False, msg)

        elif val.size not in (3, 6, 7):
            return (False, msg)

        else:
            return (True, msg)

    def check_grid_mapping(self, ds):
        prev_return = super().check_grid_mapping(ds)
        grid_mapping_variables = cfutil.get_grid_mapping_variables(ds)
        for var_name in sorted(grid_mapping_variables):
            var = ds.variables[var_name]
            test_ctx = self.get_test_ctx(
                BaseCheck.HIGH,
                self.section_titles["5.6"],
                var.name,
            )

            # TODO: check cases where crs_wkt provides part of a necessary
            #       grid_mapping attribute, or where a grid_mapping attribute
            #       overrides what has been provided in crs_wkt.
            # attempt to parse crs_wkt if it is present
            if "crs_wkt" in var.ncattrs():
                crs_wkt = var.crs_wkt
                if not isinstance(crs_wkt, str):
                    test_ctx.messages.append("crs_wkt attribute must be a string")
                    test_ctx.out_of += 1
                else:
                    try:
                        pyproj.CRS.from_wkt(crs_wkt)
                    except pyproj.exceptions.CRSError as crs_error:
                        test_ctx.messages.append(
                            f"Cannot parse crs_wkt attribute to CRS using Proj4. Proj4 error: {str(crs_error)}",
                        )
                    else:
                        test_ctx.score += 1
                    test_ctx.out_of += 1

            # existence_conditions
            exist_cond_1 = (
                self._check_gmattr_existence_condition_geoid_name_geoptl_datum_name(var)
            )
            test_ctx.assert_true(exist_cond_1[0], exist_cond_1[1])
            exist_cond_2 = self._check_gmattr_existence_condition_ell_pmerid_hdatum(var)
            test_ctx.assert_true(exist_cond_2[0], exist_cond_2[1])

            # handle vertical datum related grid_mapping attributes
            vert_datum_attrs = {}
            possible_vert_datum_attrs = {"geoid_name", "geopotential_datum_name"}
            vert_datum_attrs = possible_vert_datum_attrs.intersection(var.ncattrs())
            len_vdatum_name_attrs = len(vert_datum_attrs)
            # check that geoid_name and geopotential_datum_name are not both
            # present in the grid_mapping variable
            if len_vdatum_name_attrs == 2:
                test_ctx.out_of += 1
                test_ctx.messages.append(
                    "Cannot have both 'geoid_name' and "
                    "'geopotential_datum_name' attributes in "
                    f"grid mapping variable '{var.name}'",
                )
            elif len_vdatum_name_attrs == 1:
                # should be one or zero attrs
                proj_db_path = os.path.join(pyproj.datadir.get_data_dir(), "proj.db")
                try:
                    with sqlite3.connect(proj_db_path) as conn:
                        v_datum_attr = next(iter(vert_datum_attrs))
                        v_datum_value = getattr(var, v_datum_attr)
                        v_datum_str_valid = self._process_v_datum_str(
                            v_datum_value,
                            conn,
                        )

                        invalid_msg = (
                            f"Vertical datum value '{v_datum_value}' for "
                            f"attribute '{v_datum_attr}' in grid mapping "
                            f"variable '{var.name}' is not valid"
                        )
                        test_ctx.assert_true(v_datum_str_valid, invalid_msg)
                except sqlite3.Error as e:
                    # if we hit an error, skip the check
                    warn(
                        "Error occurred while trying to query "
                        f"Proj4 SQLite database at {proj_db_path}: {str(e)}",
                        stacklevel=2,
                    )
            prev_return[var.name] = test_ctx.to_result()

        return prev_return

    def check_standard_name_deprecated_modifiers(self, ds):
        """
        Not a standard check in that it won't raise pass/fail values,
        but instead warns upon finding deprecated CF standard name modifiers.
        :param netCDF4.Dataset ds: netCDF dataset
        """
        deprecated_var_names = cfutil._find_standard_name_modifier_variables(ds, True)
        if deprecated_var_names:
            warn(
                f"Deprecated standard_name modifiers found on variables {deprecated_var_names}",
                stacklevel=2,
            )

    def _process_v_datum_str(self, v_datum_str, conn):
        vdatum_query = """SELECT 1 FROM alias_name WHERE
                            table_name = 'vertical_datum' AND
                            alt_name = ?
                                UNION ALL
                            SELECT 1 FROM vertical_datum WHERE
                            name = ?
                            LIMIT 1"""
        res_set = conn.execute(vdatum_query, (v_datum_str, v_datum_str))
        return len(res_set.fetchall()) > 0

    def _check_dimensionless_vertical_coordinate_1_7(
        self,
        ds,
        vname,
        deprecated_units,
        ret_val,
        dim_vert_coords_dict,
    ):
        """
        Check that a dimensionless vertical coordinate variable is valid under
        CF-1.7.

        :param netCDF4.Dataset ds: open netCDF4 dataset
        :param str name: variable name
        :param list ret_val: array to append Results to
        :rtype None
        """
        variable = ds.variables[vname]
        standard_name = getattr(variable, "standard_name", None)
        formula_terms = getattr(variable, "formula_terms", None)
        # Skip the variable if it's dimensional
        correct_computed_std_name_ctx = TestCtx(
            BaseCheck.MEDIUM,
            self.section_titles["4.3"],
        )
        # IMPLEMENTATION CONFORMANCE 4.3.3 REQUIRED
        correct_computed_std_name_ctx.assert_true(
            not (formula_terms is None and hasattr(variable, "computed_standard_name")),
            f"Variable {vname} should have formula_terms attribute when "
            "computed_standard_name attribute is defined",
        )
        if formula_terms is None and standard_name not in dim_vert_coords_dict:
            return

        # assert that the computed_standard_name is maps to the standard_name correctly
        _comp_std_name = dim_vert_coords_dict[standard_name][1]
        correct_computed_std_name_ctx.assert_true(
            getattr(variable, "computed_standard_name", None) in _comp_std_name,
            f"§4.3.3 The standard_name of `{vname}` must map to the correct computed_standard_name, `{sorted(_comp_std_name)}`",
        )
        ret_val.append(correct_computed_std_name_ctx.to_result())

    def check_dimensionless_vertical_coordinates(self, ds):
        """
        Check the validity of dimensionless coordinates under CF

        CF §4.3.2 The units attribute is not required for dimensionless
        coordinates.

        The standard_name attribute associates a coordinate with its definition
        from Appendix D, Dimensionless Vertical Coordinates. The definition
        provides a mapping between the dimensionless coordinate values and
        dimensional values that can positively and uniquely indicate the
        location of the data.

        A new attribute, formula_terms, is used to associate terms in the
        definitions with variables in a netCDF file.  To maintain backwards
        compatibility with COARDS the use of these attributes is not required,
        but is strongly recommended.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []

        deprecated_units = ["level", "layer", "sigma_level"]

        # compose this function to use the results from the CF-1.6 check
        # and then extend it using a CF-1.7 addition
        ret_val.extend(
            self._check_dimensionless_vertical_coordinates(
                ds,
                deprecated_units,
                self._check_dimensionless_vertical_coordinate_1_6,
                dimless_vertical_coordinates_1_7,
            ),
        )

        ret_val.extend(
            self._check_dimensionless_vertical_coordinates(
                ds,
                deprecated_units,
                self._check_dimensionless_vertical_coordinate_1_7,
                dimless_vertical_coordinates_1_7,
            ),
        )

        return ret_val
