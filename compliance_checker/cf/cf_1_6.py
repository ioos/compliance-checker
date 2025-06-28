import difflib
import logging
from collections import defaultdict

import cftime
import numpy as np
import regex
from cf_units import Unit

import compliance_checker.cf.util as cfutil
from compliance_checker.base import BaseCheck, Result, TestCtx
from compliance_checker.cf import util
from compliance_checker.cf.appendix_c import valid_modifiers
from compliance_checker.cf.appendix_d import dimless_vertical_coordinates_1_6
from compliance_checker.cf.appendix_e import cell_methods16
from compliance_checker.cf.appendix_f import (
    grid_mapping_attr_types16,
    grid_mapping_dict16,
)
from compliance_checker.cf.cf_base import CFNCCheck, appendix_a_base

logger = logging.getLogger(__name__)


class CF1_6Check(CFNCCheck):
    """CF-1.6-specific implementation of CFBaseCheck; supports checking
    netCDF datasets.
    These checks are translated documents:
        https://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html"""

    register_checker = True
    _cc_spec = "cf"
    _cc_spec_version = "1.6"
    _cc_description = "Climate and Forecast Conventions (CF)"
    _cc_url = "http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html"
    _cc_display_headers = {3: "Errors", 2: "Warnings", 1: "Info"}
    appendix_a = appendix_a_base
    appendix_d_parametric_coords = dimless_vertical_coordinates_1_6
    _allowed_numeric_var_types = {
        np.character,
        np.bytes_,  # "|S1" dtype, byte array used as string
        np.int8,
        np.int16,
        np.int32,
        np.float32,
        np.float64,
    }

    def __init__(self, options=None):  # initialize with parent methods and data
        super().__init__(options)

        self.cell_methods = cell_methods16
        self.grid_mapping_dict = grid_mapping_dict16
        self.grid_mapping_attr_types = grid_mapping_attr_types16

    ###############################################################################
    # Chapter 2: NetCDF Files and Components
    ###############################################################################
    def check_filename(self, ds):
        """Checks that the filename ends with .nc"""
        # IMPLEMENTS CONFORMANCE 2.1
        filepath = ds.filepath()
        filename_suffix = TestCtx(BaseCheck.HIGH, self.section_titles["2.1"])
        filename_suffix.assert_true(
            filepath.endswith("nc"),
            f'Dataset path {filepath} must end with ".nc"',
        )
        return filename_suffix.to_result()

    def check_data_types(self, ds):
        """
        Checks the data type of all netCDF variables to ensure they are valid
        data types under CF.

        CF §2.2 The netCDF data types char, byte, short, int, float or real, and
        double are all acceptable

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """
        # IMPLEMENTS CONFORMANCE 2.2
        fails = []
        total = len(ds.variables)

        for k, v in ds.variables.items():
            if (
                v.dtype is not str
                and v.dtype.kind != "S"
                and v.dtype.type not in self._allowed_numeric_var_types
            ):
                fails.append(
                    f"The variable {k} failed because the datatype is {v.datatype}",
                )
        return Result(
            BaseCheck.HIGH,
            (total - len(fails), total),
            self.section_titles["2.2"],
            msgs=fails,
        )

    def check_child_attr_data_types(self, ds):
        """
        For any variables which contain any of the following attributes:
            - valid_min/valid_max
            - valid_range
            - scale_factor
            - add_offset
            - _FillValue
        the data type of the attribute must match the type of its parent variable as specified in the
        NetCDF User Guide (NUG) https://docs.unidata.ucar.edu/netcdf-c/current/attribute_conventions.html,
        referenced in the CF Conventions in Section 2.5.2
        (http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html#missing-data)

        :param netCDF4.Dataset ds: open netCDF dataset object
        :rtype: compliance_checker.base.Result
        """

        ctx = TestCtx(BaseCheck.MEDIUM, self.section_titles["2.5"])
        special_attrs = (
            "actual_range",
            "valid_min",
            "valid_max",
            "valid_range",
            "_FillValue",
        )

        for _var_name, var in ds.variables.items():
            for att_name in special_attrs:
                if att_name in var.ncattrs():
                    self._parent_var_attr_type_check(att_name, var, ctx)
        return ctx.to_result()

    # TODO: consider renaming to avoid confusion with non-underscore
    #       primary function version
    def _check_add_offset_scale_factor_type(self, variable, attr_name):
        """
        Reusable function for checking both add_offset and scale_factor.
        """

        msgs = []
        error_msg = (
            f"Variable {variable.name} and {attr_name} must be equivalent "
            f"data types or {variable.name} must be of type byte, short, or int "
            f"and {attr_name} must be float or double"
        )

        att = getattr(variable, attr_name, None)
        if not (isinstance(att, (np.number, float))):  # can't compare dtypes
            val = False

        else:
            val = (
                att.dtype == variable.dtype
            ) or (  # will short-circuit or if first condition is true
                isinstance(att, (np.float32, np.float64, float))
                and variable.dtype in (np.byte, np.short, np.int16, np.int32, int)
            )
        if not val:
            msgs.append(error_msg)

        return Result(BaseCheck.MEDIUM, val, self.section_titles["8.1"], msgs)

    def check_add_offset_scale_factor_type(self, ds):
        """
        If a variable has the attributes add_offset and scale_factor,
        check that the variables and attributes are of the same type
        OR that the variable is of type byte, short or int and the
        attributes are of type float or double.
        """

        results = []
        add_offset_vars = ds.get_variables_by_attributes(
            add_offset=lambda x: x is not None,
        )
        scale_factor_vars = ds.get_variables_by_attributes(
            scale_factor=lambda x: x is not None,
        )

        both = set(add_offset_vars).intersection(scale_factor_vars)
        both_msgs = []
        for both_var in sorted(both, key=lambda var: var.name):
            if both_var.scale_factor.dtype != both_var.add_offset.dtype:
                both_msgs.append(
                    "When both scale_factor and add_offset "
                    f"are supplied for variable {both_var.name}, "
                    "they must have the same type",
                )
        results.append(
            Result(
                BaseCheck.MEDIUM,
                not bool(both_msgs),
                self.section_titles["8.1"],
                both_msgs,
            ),
        )

        for _att_vars_tup in (
            ("add_offset", add_offset_vars),
            ("scale_factor", scale_factor_vars),
        ):
            results.extend(
                [
                    self._check_add_offset_scale_factor_type(
                        var,
                        _att_vars_tup[0],
                    )
                    for var in _att_vars_tup[1]
                ],
            )

        return results

    def check_naming_conventions(self, ds):
        """
        Checks the variable names to ensure they are valid CF variable names under CF.

        CF §2.3 Variable, dimension and attribute names should begin with a letter
        and be composed of letters, digits, and underscores.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """
        ret_val = []
        variable_naming = TestCtx(BaseCheck.MEDIUM, self.section_titles["2.3"])
        dimension_naming = TestCtx(BaseCheck.MEDIUM, self.section_titles["2.3"])
        attribute_naming = TestCtx(BaseCheck.MEDIUM, self.section_titles["2.3"])

        ignore_attributes = [
            "_FillValue",
            "DODS",
            "_ChunkSizes",
            "_Coordinate",
            "_Unsigned",
            "_Encoding",
        ]

        rname = regex.compile("^[A-Za-z][A-Za-z0-9_]*$")

        # IMPLEMENTATION CONFORMANCE 2.3 REQUIRED
        for name, variable in ds.variables.items():
            variable_naming.assert_true(
                rname.match(name) is not None,
                f"variable {name} should begin with a letter and be composed of "
                "letters, digits, and underscores",
            )

            # Keep track of all the attributes, we'll need to check them
            for attr in variable.ncattrs():
                if attr in ignore_attributes:
                    continue
                # Special attributes made by THREDDS
                if attr.startswith("DODS"):
                    continue
                # Ignore model produced attributes
                if attr.startswith("_Coordinate"):
                    continue
                attribute_naming.assert_true(
                    rname.match(attr) is not None,
                    f"attribute {name}:{attr} should begin with a letter and be composed of "
                    "letters, digits, and underscores",
                )

        ret_val.append(variable_naming.to_result())

        for dimension in ds.dimensions:
            dimension_naming.assert_true(
                rname.match(dimension) is not None,
                f"dimension {dimension} should begin with a latter and be composed of "
                "letters, digits, and underscores",
            )
        ret_val.append(dimension_naming.to_result())

        for global_attr in ds.ncattrs():
            # Special attributes made by THREDDS
            if global_attr.startswith("DODS"):
                continue
            if global_attr.startswith("EXTRA_DIMENSION"):
                continue
            attribute_naming.assert_true(
                rname.match(global_attr) is not None,
                f"global attribute {global_attr} should begin with a letter and be composed of "
                "letters, digits, and underscores",
            )
        ret_val.append(attribute_naming.to_result())

        return ret_val

    def check_names_unique(self, ds):
        """
        Checks the variable names for uniqueness regardless of case.

        CF §2.3 names should not be distinguished purely by case, i.e., if case
        is disregarded, no two names should be the same.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """
        fails = []
        total = len(ds.variables)
        names = defaultdict(int)

        # IMPLEMENTATION CONFORMANCE 2.3 RECOMMENDED
        for k in ds.variables:
            names[k.lower()] += 1

        fails = [
            f"Variables are not case sensitive. Duplicate variables named: {k}"
            for k, v in names.items()
            if v > 1
        ]
        return Result(
            BaseCheck.MEDIUM,
            (total - len(fails), total),
            self.section_titles["2.3"],
            msgs=fails,
        )

    def check_dimension_names(self, ds):
        """
        Checks variables contain no duplicate dimension names.

        CF §2.4 A variable may have any number of dimensions, including zero,
        and the dimensions must all have different names.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """
        fails = []
        total = len(ds.variables)

        for k, v in ds.variables.items():
            dims = defaultdict(int)
            for d in v.dimensions:
                dims[d] += 1

            # IMPLEMENTATION CONFORMANCE 2.4 REQUIRED
            for dimension, count in dims.items():
                if count > 1:
                    fails.append(
                        f"{k} has two or more dimensions named {dimension}",
                    )

        return Result(
            BaseCheck.HIGH,
            (total - len(fails), total),
            self.section_titles["2.4"],
            msgs=fails,
        )

    def check_dimension_order(self, ds):
        """
        Checks each variable's dimension order to ensure that the order is
        consistent and in order under CF §2.4

        CF §2.4 If any or all of the dimensions of a variable have the
        interpretations of "date or time" (T), "height or depth" (Z),
        "latitude" (Y), or "longitude" (X) then we recommend, those dimensions
        to appear in the relative order T, then Z, then Y, then X in the CDL
        definition corresponding to the file. All other dimensions should,
        whenever possible, be placed to the left of the spatiotemporal
        dimensions.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """
        valid_dimension_order = TestCtx(BaseCheck.MEDIUM, self.section_titles["2.4"])
        # Build a map from coordinate variable to axis
        coord_axis_map = self._get_coord_axis_map(ds)

        # Check each variable's dimension order, excluding climatology and
        # bounds variables
        any_clim = cfutil.get_climatology_variable(ds)
        any_bounds = cfutil.get_cell_boundary_variables(ds)
        for name, variable in ds.variables.items():
            # Skip bounds/climatology variables, as they should implicitly
            # have the same order except for the bounds specific dimension.
            # This is tested later in the respective checks
            if name in any_bounds or name == any_clim:
                continue

            # Skip strings/labels
            if hasattr(variable.dtype, "char") and variable.dtype.char == "S":
                continue
            elif variable.dtype == str:
                continue

            if variable.dimensions:
                dimension_order = self._get_dimension_order(ds, name, coord_axis_map)
                valid_dimension_order.assert_true(
                    self._dims_in_order(dimension_order),
                    "{}'s spatio-temporal dimensions are not in the "
                    "recommended order T, Z, Y, X and/or further dimensions "
                    "are not located left of T, Z, Y, X. The dimensions (and "
                    "their guessed types) are {} (with U: other/unknown; L: "
                    "unlimited).".format(
                        name,
                        self._get_pretty_dimension_order_with_type(
                            ds,
                            name,
                            dimension_order,
                        ),
                    ),
                )
        return valid_dimension_order.to_result()

    def check_fill_value_equal_missing_value(self, ds):
        """
        If both missing_value and _FillValue be used, they should have the same value.
        This according to CF §2.5.1 Recommendations:

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of Results
        """
        fails = []
        total = 0

        for variable in ds.variables.values():
            # If the variable have a defined _FillValue a defined missing_value check it.

            if hasattr(variable, "_FillValue") and hasattr(variable, "missing_value"):
                total = total + 1
                if not (
                    variable._FillValue == variable.missing_value
                    or (
                        np.isnan(variable._FillValue)
                        and np.isnan(variable.missing_value)
                    )
                ):
                    fails.append(
                        f"For the variable {variable.name} the missing_value must be equal to the _FillValue",
                    )

        return Result(
            BaseCheck.MEDIUM,
            (total - len(fails), total),
            self.section_titles["2.5"],
            msgs=fails,
        )

    def check_valid_range_and_valid_min_max_present(self, ds):
        """
        The valid_range attribute must not be present if the valid_min
        and/or valid_max attributes are present. This according to 2.5.1 Requirements.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of Results
        """
        fails = []
        total = 0

        for variable in ds.variables.values():
            if hasattr(variable, "valid_max") or hasattr(variable, "valid_min"):
                total += 1
                # if there's also valid_range in addition to
                # valid_min/valid_max, this is not compliant
                if hasattr(variable, "valid_range"):
                    fails.append(
                        f"For the variable {variable.name} the valid_range attribute must not be present "
                        "if the valid_min and/or valid_max attributes are present",
                    )
            # *Just* valid_range should be added to total as well
            elif hasattr(variable, "valid_range"):
                total += 1

        return Result(
            BaseCheck.MEDIUM,
            (total - len(fails), total),
            self.section_titles["2.5"],
            msgs=fails,
        )

    def check_fill_value_outside_valid_range(self, ds):
        """
        Checks each variable's _FillValue to ensure that it's in valid_range or
        between valid_min and valid_max according to CF §2.5.1

        CF §2.5.1 The _FillValue should be outside the range specified by
        valid_range (if used) for a variable.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of Results
        """
        valid_fill_range = TestCtx(BaseCheck.MEDIUM, self.section_titles["2.5"])

        for name, variable in ds.variables.items():
            # If the variable doesn't have a defined _FillValue don't check it.

            if not hasattr(variable, "_FillValue"):
                continue

            fill_value = variable._FillValue

            attrs = variable.ncattrs()

            if "valid_range" in attrs:
                if isinstance(variable.valid_range, str):
                    m = "§2.5.1 Fill Values should be outside the range specified by valid_range"  # subsection message
                    valid_fill_range.assert_true(
                        False,
                        f"{m};\n\t{name}:valid_range must be a numeric type not a string",
                    )
                    continue
                rmin, rmax = variable.valid_range
                spec_by = "valid_range"

            elif "valid_min" in attrs and "valid_max" in attrs:
                if isinstance(variable.valid_min, str):
                    valid_fill_range.assert_true(
                        False,
                        f"{name}:valid_min must be a numeric type not a string",
                    )
                if isinstance(variable.valid_max, str):
                    valid_fill_range.assert_true(
                        False,
                        f"{name}:valid_max must be a numeric type not a string",
                    )
                if isinstance(variable.valid_min, str) or isinstance(
                    variable.valid_max,
                    str,
                ):
                    continue
                rmin = variable.valid_min
                rmax = variable.valid_max
                spec_by = "valid_min/valid_max"
            else:
                continue

            if np.isnan(fill_value):
                valid = True
            else:
                valid = fill_value < rmin or fill_value > rmax

            valid_fill_range.assert_true(
                valid,
                f"{name}:_FillValue ({fill_value}) should be outside the range specified by {spec_by} ({rmin}, {rmax})"
                "",
            )

        return valid_fill_range.to_result()

    def check_convention_globals(self, ds):
        """
        Check the common global attributes are strings if they exist.

        CF §2.6.2 title/history global attributes, must be strings. Do not need
        to exist.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of Results
        """
        attrs = ["title", "history"]

        valid_globals = TestCtx(BaseCheck.MEDIUM, self.section_titles["2.6"])

        for attr in attrs:
            dataset_attr = getattr(ds, attr, None)
            is_string = isinstance(dataset_attr, str)
            valid_globals.assert_true(
                is_string and len(dataset_attr),
                f"§2.6.2 global attribute {attr} should exist and be a non-empty string"  # subsection message
                "",
            )
        return valid_globals.to_result()

    # IMPLEMENTATION CONFORMANCE 1.2
    def check_coordinate_variables_strict_monotonicity(self, ds):
        """
        Checks that data in coordinate variables is either monotonically
        increasing or decreasing
        """

        ret_val = []
        for coord_var_name in self._find_coord_vars(ds):
            coord_var = ds.variables[coord_var_name]
            arr_diff = np.diff(coord_var)
            monotonicity = TestCtx(BaseCheck.HIGH, self.section_titles["1.2"])
            monotonicity.assert_true(
                np.all(arr_diff > 0) or np.all(arr_diff < 0),
                f'Coordinate variable "{coord_var_name}" must be strictly monotonic',
            )
            ret_val.append(monotonicity.to_result())

        return ret_val

    def check_convention_possibly_var_attrs(self, ds):
        """
        Check variable and global attributes are strings for recommended attributes under CF §2.6.2

        CF §2.6.2 institution, source, references, and comment, either global
        or assigned to individual variables.  When an attribute appears both
        globally and as a variable attribute, the variable's version has
        precedence.  Must be strings.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of Results
        """
        # The attrs are optional and only needs to be a string and non-empty if it
        # exists.
        attrs = ["institution", "source", "references", "comment"]

        valid_attributes = TestCtx(BaseCheck.MEDIUM, self.section_titles["2.6"])

        attr_bin = set()
        # If the attribute is defined for any variable, check it and mark in
        # the set that we've seen it at least once.
        for name, variable in ds.variables.items():
            for attribute in variable.ncattrs():
                varattr = getattr(variable, attribute)
                if attribute in attrs:
                    is_string = isinstance(varattr, str)
                    valid_attributes.assert_true(
                        is_string and len(varattr) > 0,
                        f"§2.6.2 {name}:{attribute} should be a non-empty string" "",
                    )
                    attr_bin.add(attribute)

        # Check all the global attributes too and mark if we've seen them
        for attribute in ds.ncattrs():
            dsattr = getattr(ds, attribute)
            if attribute in attrs:
                is_string = isinstance(dsattr, str)
                valid_attributes.assert_true(
                    is_string and len(dsattr) > 0,
                    f"§2.6.2 {attribute} global attribute should be a non-empty string"
                    "",
                )
                attr_bin.add(attribute)
        return valid_attributes.to_result()

    ###############################################################################
    # Chapter 3: Description of the Data
    ###############################################################################

    def check_units(self, ds):
        """
        Check the units attribute for all variables to ensure they are CF
        compliant under CF §3.1

        CF §3.1 The units attribute is required for all variables that represent dimensional quantities
        (except for boundary variables defined in Section 7.1, "Cell Boundaries" and climatology variables
        defined in Section 7.4, "Climatological Statistics").

        Units are not required for dimensionless quantities. A variable with no units attribute is assumed
        to be dimensionless. However, a units attribute specifying a dimensionless unit may optionally be
        included.

        - units required
        - type must be recognized by udunits
        - if standard name specified, must be consistent with standard name table, must also be consistent with a
          specified cell_methods attribute if present

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """

        ret_val = []

        coordinate_variables = self._find_coord_vars(ds)
        auxiliary_coordinates = self._find_aux_coord_vars(ds)
        geophysical_variables = self._find_geophysical_vars(ds)
        modifier_variables = cfutil._find_standard_name_modifier_variables(ds)
        forecast_variables = cfutil.get_forecast_metadata_variables(ds)

        dimless_vert = {
            var.name
            for var in ds.get_variables_by_attributes(
                standard_name=lambda s: s in self.appendix_d_parametric_coords,
            )
            if not hasattr(var, "units")
        }
        # check anything remaining that has units
        # unit_containing =
        unit_required_variables = (
            set(
                coordinate_variables
                + auxiliary_coordinates
                + geophysical_variables
                + forecast_variables
                + modifier_variables,
            )  # standard names with modifiers require proper units, *except* for flags, where they should not be present
            - dimless_vert
        )

        for name in unit_required_variables:
            # For reduced horizontal grids, the compression index variable does
            # not require units.
            if cfutil.is_compression_coordinate(ds, name):
                continue

            variable = ds.variables[name]

            # Skip instance coordinate variables
            if getattr(variable, "cf_role", None) is not None:
                continue

            # Skip labels
            if (
                hasattr(variable.dtype, "char") and variable.dtype.char == "S"
            ) or variable.dtype == str:
                continue

            standard_name = getattr(variable, "standard_name", None)
            standard_name, standard_name_modifier = self._split_standard_name(
                standard_name,
            )

            units = getattr(variable, "units", None)

            valid_units = self._check_valid_cf_units(ds, name)
            ret_val.append(valid_units)

            units_attr_is_string = TestCtx(BaseCheck.MEDIUM, self.section_titles["3.1"])

            # side effects, but better than teasing out the individual result
            if units is not None and units_attr_is_string.assert_true(
                isinstance(units, str),
                f"units ({units}) attribute of '{variable.name}' must be a string compatible with UDUNITS",
            ):
                valid_udunits = self._check_valid_udunits(ds, name)
                ret_val.append(valid_udunits)
            ret_val.append(units_attr_is_string.to_result())

            if isinstance(standard_name, str):
                # CONFORMANCE 3.1 REQUIRED
                valid_standard_units = self._check_valid_standard_units(ds, name)
                ret_val.append(valid_standard_units)

        return ret_val

    def _check_valid_cf_units(self, ds, variable_name):
        """
        Checks that the variable contains units attribute, the attribute is a
        string and the value is not deprecated by CF

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str variable_name: Name of the variable to be checked
        :rtype:
        :return: List of results
        """

        # This list is straight from section 3
        deprecated = ["level", "layer", "sigma_level"]
        variable = ds.variables[variable_name]

        valid_units = TestCtx(BaseCheck.HIGH, self.section_titles["3.1"])
        units = getattr(variable, "units", None)
        standard_name_full = getattr(variable, "standard_name", None)
        standard_name, standard_name_modifier = self._split_standard_name(
            standard_name_full,
        )
        std_name_units_dimensionless = cfutil.is_dimensionless_standard_name(
            self._std_names._root,
            standard_name,
        )

        # 3) units are not deprecated
        valid_units.assert_true(
            units not in deprecated,
            f'units for {variable_name}, "{units}" are deprecated by CF 1.6',
        )
        # 4/5) Modifiers, if present, have the appropriate units, or none for
        #    status_flag
        if standard_name_modifier is not None:
            if standard_name_modifier not in valid_modifiers:
                # standard name modifier warning given elsewhere
                return valid_units.to_result()
            else:
                unit_type = valid_modifiers[standard_name_modifier]
        # no modifiers, just check against standard name canonical_units
        else:
            unit_type = "u"

        if unit_type == "u":
            try:
                reference = self._std_names[standard_name].canonical_units
            # if standard name isn't found, there won't be an associated units
            # but a standard name error will be raised elsewhere
            except KeyError:
                return valid_units.to_result()
        elif unit_type == "1":
            reference = "1"
        elif unit_type is None:
            valid_units.assert_true(
                units is None,
                f"units attribute for variable {variable_name} must be unset "
                "when status_flag standard name modifier is set",
            )
            return valid_units.to_result()

        # Is this even in the database? also, if there is no standard_name,
        # there's no way to know if it is dimensionless.
        should_be_dimensionless = (
            variable.dtype is str
            or (hasattr(variable.dtype, "char") and variable.dtype.char == "S")
            or std_name_units_dimensionless
            or standard_name is None
        )

        # 1) Units must exist
        valid_units.assert_true(
            should_be_dimensionless or units is not None,
            f"units attribute is required for {variable_name} when variable is not a dimensionless quantity",
        )

        # Don't bother checking the rest
        if units is None and not should_be_dimensionless:
            return valid_units.to_result()
        # 2) units attribute must be a string
        valid_units.assert_true(
            should_be_dimensionless or isinstance(units, str),
            f"units attribute for {variable_name} needs to be a string",
        )

        try:
            units_conv = Unit(units)
        except ValueError:
            valid_units.messages.append(
                f'Unit string "{units}" is not recognized by UDUnits',
            )
            valid_units.out_of += 1
            return valid_units
        else:
            valid_units.score += 1
            valid_units.out_of += 1

        # time and forecast_reference time have special unit handling rules
        # that use time relative to a reference point, despite canonical units
        # being expressed as "s"/seconds
        if standard_name not in {"time", "forecast_reference_time"}:
            valid_units.assert_true(
                units_conv.is_convertible(Unit(reference)),
                f'Units "{units}" for variable '
                f"{variable_name} must be convertible to "
                f'canonical units "{reference}"',
            )

        return valid_units.to_result()

    def _check_valid_udunits(self, ds, variable_name):
        """
        Checks that the variable's units are contained in UDUnits

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str variable_name: Name of the variable to be checked
        """
        variable = ds.variables[variable_name]

        units = getattr(variable, "units", None)
        standard_name = getattr(variable, "standard_name", None)
        standard_name, standard_name_modifier = self._split_standard_name(standard_name)
        std_name_units_dimensionless = cfutil.is_dimensionless_standard_name(
            self._std_names._root,
            standard_name,
        )

        # If the variable is supposed to be dimensionless, it automatically passes
        should_be_dimensionless = (
            variable.dtype is str
            or (hasattr(variable.dtype, "char") and variable.dtype.char == "S")
            or std_name_units_dimensionless
        )

        valid_udunits = TestCtx(BaseCheck.HIGH, self.section_titles["3.1"])
        are_udunits = units is not None and util.units_known(units)
        valid_udunits.assert_true(
            should_be_dimensionless or are_udunits or units is None,
            f'units for {variable_name}, "{units}" are not recognized by UDUNITS',
        )
        return valid_udunits.to_result()

    def _check_valid_standard_units(self, ds, variable_name):
        """
        Checks that the variable's units are appropriate for the standard name
        according to the CF standard name table and coordinate sections in CF
        1.6

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str variable_name: Name of the variable to be checked
        """
        variable = ds.variables[variable_name]
        units = getattr(variable, "units", None)
        standard_name = getattr(variable, "standard_name", None)

        valid_standard_units = TestCtx(BaseCheck.HIGH, self.section_titles["3.1"])

        # If the variable is supposed to be dimensionless, it automatically passes
        std_name_units_dimensionless = cfutil.is_dimensionless_standard_name(
            self._std_names._root,
            standard_name,
        )

        if std_name_units_dimensionless:
            return valid_standard_units.to_result()

        standard_name, standard_name_modifier = self._split_standard_name(standard_name)

        # Other standard_name modifiers have the same units as the
        # unmodified standard name or are not checked for units.

        # number_of_observations is a special case which always must be units
        # of "1"
        if standard_name_modifier == "number_of_observations":
            valid_standard_units.out_of += 1
            if units != "1":
                err_msg = (
                    f"When variable {variable_name} has a "
                    "standard name modifier of number_of_observations, "
                    "the specified units must be 1"
                )
                valid_standard_units.messages.append(err_msg)
            else:
                valid_standard_units.score += 1
            # number_of_observations should short circuit and not continue
            # on to further units checks
            return valid_standard_units.to_result()
        elif standard_name_modifier == "status_flag":
            # no units required - skip further checks
            return valid_standard_units.to_result()

        # This section represents the different cases where simple udunits
        # comparison isn't comprehensive enough to determine if the units are
        # appropriate under CF

        # UDUnits accepts "s" as a unit of time but it should be <unit> since <epoch>
        # TODO: forecast_reference_time.  Include upcoming merge.
        # IMPLEMENTATION CONFORMANCE 4.4 REQUIRED 1/2
        elif standard_name == "time":
            valid_standard_units.assert_true(
                cfutil.units_convertible(units, "seconds since 1970-01-01"),
                "time must be in a valid units format <unit> since <epoch> "
                f"not {units}",
            )

        # UDunits can't tell the difference between east and north facing coordinates
        elif standard_name == "latitude":
            # degrees is allowed if using a transformed grid
            allowed_units = cfutil.VALID_LAT_UNITS | {"degrees"}
            valid_standard_units.assert_true(
                (units.lower() if units is not None else None) in allowed_units,
                f'variables defining latitude ("{variable_name}") must use degrees_north '
                "or degrees if defining a transformed grid. Currently "
                f"{units}",
            )
        # UDunits can't tell the difference between east and north facing coordinates
        elif standard_name == "longitude":
            # degrees is allowed if using a transformed grid
            allowed_units = cfutil.VALID_LON_UNITS | {"degrees"}
            valid_standard_units.assert_true(
                (units.lower() if units is not None else None) in allowed_units,
                f'variables defining longitude ("{variable_name}") must use degrees_east '
                "or degrees if defining a transformed grid. Currently "
                f"{units}",
            )

        return valid_standard_units.to_result()

    def check_standard_name(self, ds):
        """
        Check a variables's standard_name attribute to ensure that it meets CF
        compliance.

        CF §3.3 A standard name is associated with a variable via the attribute
        standard_name which takes a string value comprised of a standard name
        optionally followed by one or more blanks and a standard name modifier

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []

        coord_vars = self._find_coord_vars(ds)
        aux_coord_vars = self._find_aux_coord_vars(ds)
        axis_vars = cfutil.get_axis_variables(ds)
        flag_vars = cfutil.get_flag_variables(ds)
        geophysical_vars = self._find_geophysical_vars(ds)

        variables_requiring_standard_names = (
            coord_vars + aux_coord_vars + axis_vars + flag_vars + geophysical_vars
        )
        for name in set(variables_requiring_standard_names):
            # Compression indices used in reduced horizontal grids or
            # compression schemes do not require attributes other than compress
            if cfutil.is_compression_coordinate(ds, name):
                continue

            ncvar = ds.variables[name]

            # §9 doesn't explicitly allow instance variables as coordinates but
            # it's loosely implied. Just in case, skip it.
            if hasattr(ncvar, "cf_role"):
                continue

            # Unfortunately, §6.1 allows for string types to be listed as
            # coordinates.
            if hasattr(ncvar.dtype, "char") and ncvar.dtype.char == "S":
                continue
            elif ncvar.dtype == str:
                continue

            standard_name = getattr(ncvar, "standard_name", None)
            standard_name, standard_name_modifier = self._split_standard_name(
                standard_name,
            )
            long_name = getattr(ncvar, "long_name", None)
            long_or_std_name = TestCtx(BaseCheck.HIGH, self.section_titles["3.3"])
            if long_name is not None:
                long_name_present = True
                long_or_std_name.assert_true(
                    isinstance(long_name, str),
                    f"Attribute long_name for variable {name} must be a string",
                )
            else:
                long_name_present = False
            # §1.3 The long_name and standard_name attributes are used to
            # describe the content of each variable. For backwards
            # compatibility with COARDS neither is required, but use of at
            # least one of them is strongly recommended.

            # If standard_name is not defined but long_name is, don't continue
            # the check for this variable
            # IMPLEMENTATION CONFORMANCE 3.3 REQUIRED 1, 2, 3 / 3
            if standard_name is not None:
                standard_name_present = True
                valid_std_name = TestCtx(BaseCheck.HIGH, self.section_titles["3.3"])
                valid_std_name.assert_true(
                    isinstance(standard_name, str),
                    f"Attribute standard_name for variable {name} must be a string",
                )
                valid_std_name.out_of += 1
                if standard_name not in self._std_names:
                    err_msg = "standard_name {} is not defined in Standard Name Table v{}.".format(
                        standard_name or "undefined",
                        self._std_names._version,
                    )
                    close_matches = difflib.get_close_matches(
                        standard_name,
                        self._std_names,
                    )
                    if close_matches:
                        err_msg += f" Possible close match(es): {close_matches}"
                    valid_std_name.messages.append(err_msg)
                else:
                    valid_std_name.score += 1

                ret_val.append(valid_std_name.to_result())

                # 2) optional - if modifiers, should be in table
                if standard_name_modifier is not None:
                    valid_modifier = TestCtx(BaseCheck.HIGH, self.section_titles["3.3"])
                    valid_modifier.assert_true(
                        standard_name_modifier in valid_modifiers,
                        f'Standard name modifier "{standard_name_modifier}" for variable {name} is not a valid modifier '
                        "according to CF Appendix C",
                    )

                    ret_val.append(valid_modifier.to_result())
            else:
                standard_name_present = False

            # IMPLEMENTATION CONFORMANCE 3 RECOMMENDED
            long_or_std_name.assert_true(
                long_name_present or standard_name_present,
                f"Attribute long_name or/and standard_name is highly recommended for variable {name}",
            )
            ret_val.append(long_or_std_name.to_result())
        return ret_val

    def check_ancillary_variables(self, ds):
        """
        Checks the ancillary_variable attribute for all variables to ensure
        they are CF compliant.

        CF §3.4 It is a string attribute whose value is a blank separated list
        of variable names.  The nature of the relationship between variables
        associated via ancillary_variables must be determined by other
        attributes. The variables listed by the ancillary_variables attribute
        will often have the standard name of the variable which points to them
        including a modifier (Appendix C, Standard Name Modifiers) to indicate
        the relationship.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []

        for ncvar in ds.get_variables_by_attributes(
            ancillary_variables=lambda x: x is not None,
        ):
            name = ncvar.name
            valid_ancillary = TestCtx(BaseCheck.HIGH, self.section_titles["3.4"])
            ancillary_variables = ncvar.ancillary_variables

            valid_ancillary.assert_true(
                isinstance(ancillary_variables, str),
                f"ancillary_variables attribute defined by {name} " "should be string",
            )

            # Can't perform the second check if it's not a string
            if not isinstance(ancillary_variables, str):
                ret_val.append(valid_ancillary.to_result())
                continue

            for ancillary_variable in ancillary_variables.split():
                valid_ancillary.assert_true(
                    ancillary_variable in ds.variables,
                    f"{ancillary_variable} is not a variable in this dataset",
                )

            ret_val.append(valid_ancillary.to_result())

        return ret_val

    def check_flags(self, ds):
        """
        Check the flag_values, flag_masks and flag_meanings attributes for
        variables to ensure they are CF compliant.

        CF §3.5 The attributes flag_values, flag_masks and flag_meanings are
        intended to make variables that contain flag values self describing.
        Status codes and Boolean (binary) condition flags may be expressed with
        different combinations of flag_values and flag_masks attribute
        definitions.

        The flag_values and flag_meanings attributes describe a status flag
        consisting of mutually exclusive coded values.

        The flag_meanings attribute is a string whose value is a blank
        separated list of descriptive words or phrases, one for each flag
        value. Each word or phrase should consist of characters from the
        alphanumeric set and the following five: '_', '-', '.', '+', '@'.

        The flag_masks and flag_meanings attributes describe a number of
        independent Boolean conditions using bit field notation by setting
        unique bits in each flag_masks value.

        The flag_masks, flag_values and flag_meanings attributes, used
        together, describe a blend of independent Boolean conditions and
        enumerated status codes. A flagged condition is identified by a bitwise
        AND of the variable value and each flag_masks value; a result that
        matches the flag_values value indicates a true condition.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []

        for name in cfutil.get_flag_variables(ds):
            variable = ds.variables[name]
            flag_values = getattr(variable, "flag_values", None)
            flag_masks = getattr(variable, "flag_masks", None)

            valid_flags_var = TestCtx(BaseCheck.HIGH, self.section_titles["3.5"])
            # Check that the variable defines mask or values
            valid_flags_var.assert_true(
                flag_values is not None or flag_masks is not None,
                f"{name} does not define either flag_masks or flag_values",
            )
            ret_val.append(valid_flags_var.to_result())

            valid_meanings = self._check_flag_meanings(ds, name)
            ret_val.append(valid_meanings)

            # check flag_values
            if flag_values is not None:
                valid_values = self._check_flag_values(ds, name)
                ret_val.append(valid_values)

            # check flag_masks
            if flag_masks is not None:
                valid_masks = self._check_flag_masks(ds, name)
                ret_val.append(valid_masks)

            if flag_values is not None and flag_masks is not None:
                vals_arr = np.array(flag_values, ndmin=1)
                masks_arr = np.array(flag_masks, ndmin=1)
                # IMPLEMENTATION CONFORMANCE 3.5 RECOMMENDED 1/1
                # If shapes aren't equal, we can't do proper elementwise
                # comparison
                if vals_arr.size != masks_arr.size or not (
                    np.issubdtype(vals_arr.dtype, np.integer)
                    and np.issubdtype(masks_arr.dtype, np.integer)
                ):
                    allv = False
                else:
                    allv = np.all(vals_arr & masks_arr == vals_arr)

                allvr = Result(BaseCheck.MEDIUM, allv, self.section_titles["3.5"])
                if not allvr.value:
                    allvr.msgs = [
                        f"flag masks and flag values for '{name}' combined don't equal flag values",
                    ]

                ret_val.append(allvr)

        return ret_val

    def _check_flag_values(self, ds, name):
        """
        Checks a variable's flag_values attribute for compliance under CF

        - flag_values exists as an array
        - unique elements in flag_values
        - flag_values si the same dtype as the variable
        - flag_values is the same length as flag_meanings

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Name of variable to check
        :rtype: compliance_checker.base.Result
        """
        variable = ds.variables[name]

        flag_values = getattr(variable, "flag_values", None)
        flag_meanings = getattr(variable, "flag_meanings", None)
        valid_values = TestCtx(BaseCheck.HIGH, self.section_titles["3.5"])

        # IMPLEMENTATION CONFORMANCE 3.5 REQUIRED 2/8
        valid_values.assert_true(
            hasattr(variable, "flag_meanings"),
            f"Variable {variable.name} must have attribute flag_meanings "
            "defined when flag_values attribute is present",
        )

        # the flag values must be independent, no repeating values
        flag_set = np.unique(flag_values)
        valid_values.assert_true(
            flag_set.size == np.array(flag_values).size,
            f"{name}'s flag_values must be independent and can not be repeated",
        )

        # IMPLEMENTATION CONFORMANCE 3.5 REQUIRED 1/8
        # the data type for flag_values should be the same as the variable
        flag_values_type = (
            flag_values.dtype.type
            if hasattr(flag_values, "dtype")
            else type(flag_values)
        )
        valid_values.assert_true(
            variable.dtype.type == flag_values_type,
            f"flag_values ({flag_values_type}) must be the same data type as {name} ({variable.dtype.type})"
            "",
        )

        # IMPLEMENTATION CONFORMANCE 3.5 REQUIRED 4/8
        if isinstance(flag_meanings, str):
            flag_meanings = flag_meanings.split()
            valid_values.assert_true(
                len(flag_meanings) == np.array(flag_values).size,
                f"{name}'s flag_meanings and flag_values should have the same "
                "number of elements.",
            )

        return valid_values.to_result()

    def _check_flag_masks(self, ds, name):
        """
        Check a variable's flag_masks attribute for compliance under CF

        - flag_masks exists as an array
        - flag_masks is the same dtype as the variable
        - variable's dtype can support bit-field
        - flag_masks is the same length as flag_meanings

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Variable name
        :rtype: compliance_checker.base.Result
        """
        variable = ds.variables[name]

        flag_masks = variable.flag_masks
        flag_meanings = getattr(variable, "flag_meanings", None)

        valid_masks = TestCtx(BaseCheck.HIGH, self.section_titles["3.5"])

        flag_masks_type = (
            flag_masks.dtype.type if hasattr(flag_masks, "dtype") else type(flag_masks)
        )
        valid_masks.assert_true(
            variable.dtype.type == flag_masks_type,
            f"flag_masks ({flag_masks_type}) must be the same data type as {name} ({variable.dtype.type})"
            "",
        )

        type_ok = (
            np.issubdtype(variable.dtype, np.integer)
            or np.issubdtype(variable.dtype, "S")
            or np.issubdtype(variable.dtype, "b")
        )

        valid_masks.assert_true(
            0 not in np.array(flag_masks),
            f"flag_masks for variable {variable.name} must "
            "not contain zero as an element",
        )

        valid_masks.assert_true(
            type_ok,
            f"{name}'s data type must be capable of bit-field expression",
        )

        if isinstance(flag_meanings, str):
            flag_meanings = flag_meanings.split()
            valid_masks.assert_true(
                # cast to array here as single element arrays are returned as
                # scalars from netCDF4 Python
                len(flag_meanings) == np.array(flag_masks).size,
                f"{name} flag_meanings and flag_masks should have the same "
                "number of elements.",
            )

        return valid_masks.to_result()

    def _check_flag_meanings(self, ds, name):
        """
        Check a variable's flag_meanings attribute for compliance under CF

        - flag_meanings exists
        - flag_meanings is a string
        - flag_meanings elements are valid strings

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Variable name
        :rtype: compliance_checker.base.Result
        """
        variable = ds.variables[name]
        flag_meanings = getattr(variable, "flag_meanings", None)
        valid_meanings = TestCtx(BaseCheck.HIGH, self.section_titles["3.5"])

        valid_meanings.assert_true(
            flag_meanings is not None,
            f"{name}'s flag_meanings attribute is required for flag variables",
        )

        valid_meanings.assert_true(
            isinstance(flag_meanings, str),
            f"{name}'s flag_meanings attribute must be a string",
        )

        # We can't perform any additional checks if it's not a string
        if not isinstance(flag_meanings, str):
            return valid_meanings.to_result()

        valid_meanings.assert_true(
            len(flag_meanings) > 0,
            f"{name}'s flag_meanings can't be empty",
        )

        # IMPLEMENTATION CONFORMANCE REQUIRED 3.5 3/8
        flag_regx = regex.compile(r"^[0-9A-Za-z_\-.+@]+$")
        meanings = flag_meanings.split()
        for meaning in meanings:
            if flag_regx.match(meaning) is None:
                valid_meanings.assert_true(
                    False,
                    f"{name}'s flag_meanings attribute defined an illegal flag meaning "
                    + f"{meaning}",
                )
        return valid_meanings.to_result()

    ###############################################################################
    # Chapter 4: Coordinate Types
    ###############################################################################

    def check_coordinate_types(self, ds):
        """
        Check the axis attribute of coordinate variables

        CF §4 The attribute axis may be attached to a coordinate variable and
        given one of the values X, Y, Z or T which stand for a longitude,
        latitude, vertical, or time axis respectively. Alternatively the
        standard_name attribute may be used for direct identification.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []

        for variable in ds.get_variables_by_attributes(axis=lambda x: x is not None):
            name = variable.name
            # Coordinate compressions should not be checked as a valid
            # coordinate, which they are not. They are a mechanism to project
            # an array of indices onto a 2-d grid containing valid coordinates.
            if cfutil.is_compression_coordinate(ds, name):
                continue

            variable = ds.variables[name]
            # Even though it's not allowed in CF 1.6, it is allowed in CF 1.7
            # and we see people do it, often.
            if hasattr(variable, "cf_role"):
                continue

            # §6.1 allows for labels to be referenced as auxiliary coordinate
            # variables, which should not be checked like the rest of the
            # coordinates.
            if hasattr(variable.dtype, "char") and variable.dtype.char == "S":
                continue
            elif variable.dtype == str:
                continue

            axis = getattr(variable, "axis", None)

            if axis is not None:
                valid_axis = self._check_axis(ds, name)
                ret_val.append(valid_axis)

        return ret_val

    def _check_axis(self, ds, name):
        """
        Checks that the axis attribute is a string and an allowed value, namely
        one of 'T', 'X', 'Y', or 'Z'.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :param str name: Name of the variable
        :rtype: compliance_checker.base.Result
        """
        allowed_axis = ["T", "X", "Y", "Z"]
        variable = ds.variables[name]
        axis = variable.axis

        valid_axis = TestCtx(BaseCheck.HIGH, self.section_titles["4"])
        axis_is_string = (isinstance(axis, str),)
        valid_axis.assert_true(
            axis_is_string and len(axis) > 0,
            f"{name}'s axis attribute must be a non-empty string",
        )

        # If axis isn't a string we can't continue any checks
        if not axis_is_string or len(axis) == 0:
            return valid_axis.to_result()

        valid_axis.assert_true(
            axis in allowed_axis,
            f"{name}'s axis attribute must be T, X, Y, or Z, " + f"currently {axis}",
        )

        return valid_axis.to_result()

    def check_latitude(self, ds):
        """
        Check variable(s) that define latitude and are defined correctly according to CF.

        CF §4.1 Variables representing latitude must always explicitly include
        the units attribute; there is no default value.  The recommended unit
        of latitude is degrees_north. Also acceptable are degree_north,
        degree_N, degrees_N, degreeN, and degreesN.

        Optionally, the latitude type may be indicated additionally by
        providing the standard_name attribute with the value latitude, and/or
        the axis attribute with the value Y.

        - Four checks per latitude variable
        - (H) latitude has units attribute
        - (M) latitude has an allowed units attribute
        - (L) latitude uses degrees_north (if not in rotated pole)
        - (M) latitude defines either standard_name or axis

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []

        allowed_lat_units = [
            "degrees_north",
            "degree_north",
            "degree_n",
            "degrees_n",
            "degreen",
            "degreesn",
        ]

        # Determine the grid mappings in this dataset
        grid_mapping = []
        grid_mapping_variables = cfutil.get_grid_mapping_variables(ds)
        for name in grid_mapping_variables:
            variable = ds.variables[name]
            grid_mapping_name = getattr(variable, "grid_mapping_name", None)
            if grid_mapping_name:
                grid_mapping.append(grid_mapping_name)

        latitude_variables = cfutil.get_latitude_variables(ds)
        for latitude in latitude_variables:
            variable = ds.variables[latitude]
            units = getattr(variable, "units", None)
            units_is_string = isinstance(units, str)
            standard_name = getattr(variable, "standard_name", None)
            axis = getattr(variable, "axis", None)

            # Check that latitude defines units
            valid_latitude = TestCtx(BaseCheck.HIGH, self.section_titles["4.1"])
            valid_latitude.assert_true(
                units is not None,
                f"latitude variable '{latitude}' must define units",
            )
            ret_val.append(valid_latitude.to_result())

            # Check that latitude uses allowed units
            allowed_units = TestCtx(BaseCheck.MEDIUM, self.section_titles["4.1"])
            if standard_name == "grid_latitude":
                e_n_units = cfutil.VALID_LAT_UNITS | cfutil.VALID_LON_UNITS
                # check that the units aren't in east and north degrees units,
                # but are convertible to angular units
                allowed_units.assert_true(
                    units not in e_n_units and Unit(units) == Unit("degree"),
                    f"Grid latitude variable '{latitude}' should use degree equivalent units without east or north components. "
                    f"Current units are {units}",
                )
            else:
                allowed_units.assert_true(
                    units_is_string and units.lower() in allowed_lat_units,
                    f"latitude variable '{latitude}' should define valid units for latitude"
                    "",
                )
            ret_val.append(allowed_units.to_result())

            # Check that latitude uses degrees_north
            if standard_name == "latitude" and units != "degrees_north":
                # This is only a recommendation and we won't penalize but we
                # will include a recommended action.
                msg = (
                    f"CF recommends latitude variable '{latitude}' to use units degrees_north"
                    ""
                )
                recommended_units = Result(
                    BaseCheck.LOW,
                    (1, 1),
                    self.section_titles["4.1"],
                    [msg],
                )
                ret_val.append(recommended_units)

            y_variables = ds.get_variables_by_attributes(axis="Y")
            # Check that latitude defines either standard_name or axis
            definition = TestCtx(BaseCheck.MEDIUM, self.section_titles["4.1"])
            definition.assert_true(
                standard_name == "latitude" or axis == "Y" or y_variables != [],
                f"latitude variable '{latitude}' should define standard_name='latitude' or axis='Y'"
                "",
            )
            ret_val.append(definition.to_result())

        return ret_val

    def check_longitude(self, ds):
        """
        Check variable(s) that define longitude and are defined correctly according to CF.

        CF §4.2 Variables representing longitude must always explicitly include
        the units attribute; there is no default value.  The recommended unit
        of longitude is degrees_east. Also acceptable are degree_east,
        degree_E, degrees_E, degreeE, and degreesE.

        Optionally, the longitude type may be indicated additionally by
        providing the standard_name attribute with the value longitude, and/or
        the axis attribute with the value X.

        - Four checks per longitude variable
        - (H) longitude has units attribute
        - (M) longitude has an allowed units attribute
        - (L) longitude uses degrees_east (if not in rotated pole)
        - (M) longitude defines either standard_name or axis

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """

        # TODO we already have a check_latitude... I'm sure we can make DRYer

        ret_val = []
        allowed_lon_units = [
            "degrees_east",
            "degree_east",
            "degree_e",
            "degrees_e",
            "degreee",
            "degreese",
        ]

        # Determine the grid mappings in this dataset
        grid_mapping = []
        grid_mapping_variables = cfutil.get_grid_mapping_variables(ds)
        for name in grid_mapping_variables:
            variable = ds.variables[name]
            grid_mapping_name = getattr(variable, "grid_mapping_name", None)
            if grid_mapping_name:
                grid_mapping.append(grid_mapping_name)

        longitude_variables = cfutil.get_longitude_variables(ds)
        for longitude in longitude_variables:
            variable = ds.variables[longitude]
            units = getattr(variable, "units", None)
            units_is_string = isinstance(units, str)
            standard_name = getattr(variable, "standard_name", None)
            axis = getattr(variable, "axis", None)

            # NOTE see docstring--should below be 4.1 or 4.2?
            # Check that longitude defines units
            valid_longitude = TestCtx(BaseCheck.HIGH, self.section_titles["4.2"])
            valid_longitude.assert_true(
                units is not None,
                f"longitude variable '{longitude}' must define units",
            )
            ret_val.append(valid_longitude.to_result())

            # Check that longitude uses allowed units
            allowed_units = TestCtx(BaseCheck.MEDIUM, self.section_titles["4.2"])
            if standard_name == "grid_longitude":
                e_n_units = cfutil.VALID_LAT_UNITS | cfutil.VALID_LON_UNITS
                # check that the units aren't in east and north degrees units,
                # but are convertible to angular units
                allowed_units.assert_true(
                    units not in e_n_units and Unit(units) == Unit("degree"),
                    f"Grid longitude variable '{longitude}' should use degree equivalent units without east or north components. "
                    f"Current units are {units}",
                )
            else:
                allowed_units.assert_true(
                    units_is_string and units.lower() in allowed_lon_units,
                    f"longitude variable '{longitude}' should define valid units for longitude"
                    "",
                )
            ret_val.append(allowed_units.to_result())

            # Check that longitude uses degrees_east
            if standard_name == "longitude" and units != "degrees_east":
                # This is only a recommendation and we won't penalize but we
                # will include a recommended action.
                msg = (
                    f"CF recommends longitude variable '{longitude}' to use units degrees_east"
                    ""
                )
                recommended_units = Result(
                    BaseCheck.LOW,
                    (1, 1),
                    self.section_titles["4.2"],
                    [msg],
                )
                ret_val.append(recommended_units)

            x_variables = ds.get_variables_by_attributes(axis="X")
            # Check that longitude defines either standard_name or axis
            definition = TestCtx(BaseCheck.MEDIUM, self.section_titles["4.2"])
            definition.assert_true(
                standard_name == "longitude" or axis == "X" or x_variables != [],
                f"longitude variable '{longitude}' should define standard_name='longitude' or axis='X'"
                "",
            )
            ret_val.append(definition.to_result())

        return ret_val

    def check_dimensional_vertical_coordinate(
        self,
        ds,
        dimless_vertical_coordinates=dimless_vertical_coordinates_1_6,
    ):
        """
        Check units for variables defining vertical position are valid under
        CF.

        CF §4.3.1 The units attribute for dimensional coordinates will be a string
        formatted as per the udunits.dat file.

        The acceptable units for vertical (depth or height) coordinate variables
        are:
        - units of pressure as listed in the file udunits.dat. For vertical axes
          the most commonly used of these include include bar, millibar,
          decibar, atmosphere (atm), pascal (Pa), and hPa.
        - units of length as listed in the file udunits.dat. For vertical axes
          the most commonly used of these include meter (metre, m), and
          kilometer (km).
        - other units listed in the file udunits.dat that may under certain
          circumstances reference vertical position such as units of density or
          temperature.

        Plural forms are also acceptable.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        z_variables = cfutil.get_z_variables(ds)
        # dimless_standard_names = [name for name, regx in dimless_vertical_coordinates]
        for name in z_variables:
            variable = ds.variables[name]
            standard_name = getattr(variable, "standard_name", None)
            units = getattr(variable, "units", None)
            positive = getattr(variable, "positive", None)
            # Skip the variable if it's dimensionless
            if (
                hasattr(variable, "formula_terms")
                or standard_name in dimless_vertical_coordinates
            ):
                continue

            valid_vertical_coord = TestCtx(BaseCheck.HIGH, self.section_titles["4.3"])
            valid_vertical_coord.assert_true(
                isinstance(units, str) and units,
                f"§4.3.1 {name}'s units must be defined for vertical coordinates, "
                "there is no default",
            )

            if not cfutil.units_convertible("bar", units):
                valid_vertical_coord.assert_true(
                    positive in ("up", "down"),
                    f"{name}: vertical coordinates not defining pressure must include "
                    "a positive attribute that is either 'up' or 'down'",
                )

            # _check_valid_standard_units, part of the Chapter 3 checks,
            # already verifies that this coordinate has valid units

            ret_val.append(valid_vertical_coord.to_result())

        return ret_val

    def _check_dimensionless_vertical_coordinate_1_6(
        self,
        ds,
        vname,
        deprecated_units,
        ret_val,
        dim_vert_coords_dict,
    ):
        """
        Check that a dimensionless vertical coordinate variable is valid under
        CF-1.6.

        :param netCDF4.Dataset ds: open netCDF4 dataset
        :param str name: variable name
        :param list ret_val: array to append Results to
        :rtype None
        """
        variable = ds.variables[vname]
        standard_name = getattr(variable, "standard_name", None)
        units = getattr(variable, "units", None)
        formula_terms = getattr(variable, "formula_terms", None)
        # Skip the variable if it's dimensional
        if formula_terms is None and standard_name not in dim_vert_coords_dict:
            return

        is_not_deprecated = TestCtx(BaseCheck.LOW, self.section_titles["4.3"])

        is_not_deprecated.assert_true(
            units not in deprecated_units,
            f"§4.3.2: units are deprecated by CF in variable {vname}: {units}" "",
        )

        # check the vertical coordinates
        ret_val.append(is_not_deprecated.to_result())
        ret_val.append(self._check_formula_terms(ds, vname, dim_vert_coords_dict))

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

        ret_val.extend(
            self._check_dimensionless_vertical_coordinates(
                ds,
                deprecated_units,
                self._check_dimensionless_vertical_coordinate_1_6,
                dimless_vertical_coordinates_1_6,
            ),
        )

        return ret_val

    def check_time_coordinate(self, ds):
        """
        Check variables defining time are valid under CF

        CF §4.4 Variables representing time must always explicitly include the
        units attribute; there is no default value.

        The units attribute takes a string value formatted as per the
        recommendations in the Udunits package.

        The acceptable units for time are listed in the udunits.dat file. The
        most commonly used of these strings (and their abbreviations) includes
        day (d), hour (hr, h), minute (min) and second (sec, s). Plural forms
        are also acceptable. The reference time string (appearing after the
        identifier since) may include date alone; date and time; or date, time,
        and time zone. The reference time is required. A reference time in year
        0 has a special meaning (see Section 7.4, "Climatological Statistics").

        Recommend that the unit year be used with caution. It is not a calendar
        year.  For similar reasons the unit month should also be used with
        caution.

        A time coordinate is identifiable from its units string alone.
        Optionally, the time coordinate may be indicated additionally by
        providing the standard_name attribute with an appropriate value, and/or
        the axis attribute with the value T.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """

        ret_val = []
        for name in cfutil.get_time_variables(ds):
            variable = ds.variables[name]
            # Has units
            has_units = hasattr(variable, "units")
            if not has_units:
                result = Result(
                    BaseCheck.HIGH,
                    False,
                    self.section_titles["4.4"],
                    [f"{name} does not have units"],
                )
                ret_val.append(result)
                continue
            # Correct and identifiable units
            # TODO: year zero climatological time warning
            result = Result(BaseCheck.HIGH, True, self.section_titles["4.4"])
            ret_val.append(result)
            correct_units = util.units_temporal(variable.units)
            reasoning = None
            if not correct_units:
                reasoning = [f"{name} does not have correct time units"]
                result = Result(
                    BaseCheck.HIGH,
                    correct_units,
                    self.section_titles["4.4"],
                    reasoning,
                )
                ret_val.append(result)
                continue
            # IMPLEMENTATION CONFORMANCE 4.4 RECOMMENDED 1/2
            if hasattr(variable, "climatology"):
                year_match = regex.match(r"\w+ since (?P<year>\d{1,4})", variable.units)
                # year should always exist at this point if it's been parsed as
                # valid date
                if int(year_match.group("year")) == 0:
                    message = (
                        f"Time coordinate variable {variable.name}'s "
                        "use of year 0 for climatological time is "
                        "deprecated"
                    )
                    result = Result(
                        BaseCheck.MEDIUM,
                        False,
                        self.section_titles["4.4"],
                        [message],
                    )
                    ret_val.append(result)
            # IMPLEMENTATION CONFORMANCE 4.4 RECOMMENDED 2/2
            # catch non-recommended months or years time interval
            if any(unit in variable.units for unit in ("months", "years")):
                message = f"Using relative time interval of months or years is not recommended for coordinate variable {variable.name}"
                result = Result(
                    BaseCheck.MEDIUM,
                    False,
                    self.section_titles["4.4"],
                    [message],
                )
                ret_val.append(result)
        return ret_val

    def check_calendar(self, ds):
        """
        Check the calendar attribute for variables defining time and ensure it
        is a valid calendar prescribed by CF.

        CF §4.4.1 In order to calculate a new date and time given a base date, base
        time and a time increment one must know what calendar to use.

        The values currently defined for calendar are:
        - gregorian or standard
        - proleptic_gregorian
        - noleap or 365_day
        - all_leap or 366_day
        - 360_day
        - julian
        - none

        The calendar attribute may be set to none in climate experiments that
        simulate a fixed time of year.
        The time of year is indicated by the date in the reference time of the
        units attribute.

        If none of the calendars defined above applies, a non-standard calendar
        can be defined. The lengths of each month are explicitly defined with
        the month_lengths attribute of the time axis.

        If leap years are included, then two other attributes of the time axis
        should also be defined:

        leap_year, leap_month

        The calendar attribute is not required when a non-standard calendar is
        being used. It is sufficient to define the calendar using the
        month_lengths attribute, along with leap_year, and leap_month as
        appropriate. However, the calendar attribute is allowed to take
        non-standard values and in that case defining the non-standard calendar
        using the appropriate attributes is required.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        standard_calendars = {
            "gregorian",
            "standard",
            "proleptic_gregorian",
            "noleap",
            "365_day",
            "all_leap",
            "366_day",
            "360_day",
            "julian",
            "none",
        }

        ret_val = []

        def check_standard_calendar_no_cross(time_var):
            """
            Check that the time variable does not cross the date
            1582-10-15 when standard or gregorian calendars are used
            """
            # Short-circuit if using months/years.
            if any(unit in time_var.units for unit in ("months", "years")):
                return Result(
                    BaseCheck.LOW,
                    False,
                    self.section_titles["4.4"],
                    [
                        "Miscellaneous failure when attempting to calculate crossover, possible malformed date",
                    ],
                )
            time_values = time_var[:]

            # IMPLEMENTATION CONFORMANCE 4.4.1 RECOMMENDED 2/2
            # Only get non-nan/FillValue times, as these are the only things
            # that make sense for conversion.  Furthermore, non-null checks
            # should be made for time coordinate variables anyways, so errors
            # should be caught where implemented there
            crossover_date = cftime.DatetimeGregorian(1582, 10, 15)
            crossover_date_value = cftime.date2num(
                crossover_date,
                time_var.units,
                calendar=time_var.calendar,
                has_year_zero=True,
            )
            # has_year_zero set to true in order to just check crossover,
            # actual year less than or equal to zero check handled elsewhere
            # when standard/Gregorian, or Julian calendars used.

            # Comparing cftime objects is awfully slow. Converting them toordinal makes this a bit faster.
            # https://github.com/ioos/compliance-checker/issues/1211
            crossover_1582 = np.any(time_values < crossover_date_value) and np.any(
                time_values >= crossover_date_value,
            )
            if not crossover_1582:
                reasoning = (
                    f"Variable {time_var.name} has standard or Gregorian "
                    "calendar and does not cross 1582-10-15T00:00Z"
                )
            else:
                reasoning = (
                    f"Variable {time_var.name} has time values "
                    "prior to 1582-10-15T00:00Z and utilizes "
                    "the standard or Gregorian calendar"
                )

            return Result(
                BaseCheck.LOW,
                not crossover_1582,
                self.section_titles["4.4"],
                [reasoning],
            )

        # if has a calendar, check that it is within the valid values
        # otherwise no calendar is valid

        # this will only fetch variables with time units defined
        for time_var_name in cfutil.get_time_variables(ds):
            if time_var_name not in {var.name for var in util.find_coord_vars(ds)}:
                continue
            time_var = ds.variables[time_var_name]
            if not hasattr(time_var, "calendar"):
                continue
            if time_var.calendar.lower() == "gregorian":
                reasoning = (
                    f"For time variable {time_var.name}, when using "
                    "the standard Gregorian calendar, the value "
                    '"standard" is preferred over "gregorian" for '
                    "the calendar attribute"
                )
                result = Result(
                    BaseCheck.LOW,
                    False,
                    self.section_titles["4.4.1"],
                    [reasoning],
                )
                ret_val.append(result)
                # check here and in the below case that time does not cross
                # thee date 1582-10-15 as requested by CF conformance
                ret_val.append(check_standard_calendar_no_cross(time_var))
            elif time_var.calendar == "standard":
                ret_val.append(check_standard_calendar_no_cross(time_var))
            # if a nonstandard calendar, then leap_years and leap_months must
            # must be present
            if time_var.calendar.lower() not in standard_calendars:
                result = self._check_leap_time(time_var)
            # passes if the calendar is valid, otherwise notify of invalid
            # calendar
            else:
                result = Result(BaseCheck.LOW, True, self.section_titles["4.4.1"])
            ret_val.append(result)

        return ret_val

    def _check_leap_time(self, time_variable):
        """
        Helper method to handle checking custom calendar leap time specifications
        """
        leap_time = TestCtx(BaseCheck.HIGH, self.section_titles["4.4"])
        leap_time.out_of = 1
        # IMPLEMENTATION CONFORMANCE 4.4.1 REQUIRED 2, 3 / 5
        if not hasattr(time_variable, "month_lengths") or not (
            hasattr(time_variable.month_lengths, "dtype")
            and np.issubdtype(time_variable.month_lengths.dtype, np.integer)
            and time_variable.month_lengths.size == 12
        ):
            leap_time.messages.append(
                f"For nonstandard calendar on variable {time_variable.name}, "
                "attribute month_lengths must be supplied as a 12-element "
                "integer array",
            )
            return leap_time.to_result()
        # If leap years are included, then attributes leap_month and
        # leap_year must be included.
        has_leap_year = hasattr(time_variable, "leap_year")
        # IMPLEMENTATION CONFORMANCE 4.4.1 REQUIRED 4,5/5
        if hasattr(time_variable, "leap_month"):
            leap_time.assert_true(
                (
                    np.isscalar(time_variable.leap_month)
                    and hasattr(time_variable.leap_month, "dtype")
                    and np.issubdtype(time_variable.leap_month.dtype, np.integer)
                    and 1 <= time_variable.leap_month <= 12
                ),
                "When attribute leap_month is supplied for variable "
                f"{time_variable.name}, the value must be a scalar integer "
                "between 1 and 12",
            )
            # IMPLEMENTATION CONFORMANCE 4.4.1 RECOMMENDED 1/2
            if not has_leap_year:
                leap_time.out_of += 1
                fail_message = (
                    f"For time variable {time_variable.name}, "
                    "attribute leap_year must be present if "
                    "leap_month attribute is defined"
                )
                leap_time.messages.append(fail_message)

        # IMPLEMENTATION CONFORMANCE 4.4.1 REQUIRED 5/5
        if has_leap_year:
            leap_time.assert_true(
                np.isscalar(time_variable.leap_year)
                and hasattr(time_variable.leap_year, "dtype"),
                "When attribute leap_year is supplied for variable "
                f"{time_variable.name}, the value must be a scalar "
                "integer",
            )
        return leap_time.to_result()

    ###############################################################################
    # Chapter 5: Coordinate Systems
    ###############################################################################

    def check_aux_coordinates(self, ds):
        """
        Chapter 5 paragraph 3

        The dimensions of an auxiliary coordinate variable must be a subset of
        the dimensions of the variable with which the coordinate is associated,
        with two exceptions. First, string-valued coordinates (Section 6.1,
        "Labels") have a dimension for maximum string length. Second, in the
        ragged array representations of data (Chapter 9, Discrete Sampling
        Geometries), special methods are needed to connect the data and
        coordinates.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """

        ret_val = []

        # for contiguous ragged array/indexed ragged array representations,
        # coordinates are not required to adhere to the same principles;
        # these representaitions can be identified by two attributes:

        # required for contiguous
        count_vars = ds.get_variables_by_attributes(
            sample_dimension=lambda x: x is not None,
        )

        # required for indexed
        index_vars = ds.get_variables_by_attributes(
            instance_dimension=lambda x: x is not None,
        )

        # if these attributes exist, we don't need to test
        # the coordinates
        if count_vars or index_vars:
            return ret_val

        geophysical_variables = self._find_geophysical_vars(ds)
        for name in geophysical_variables:
            variable = ds.variables[name]
            coordinates = getattr(variable, "coordinates", None)
            # We use a set so we can assert
            dim_set = set(variable.dimensions)
            # No auxiliary coordinates, no check
            if not isinstance(coordinates, str) or coordinates == "":
                continue

            valid_aux_coords = TestCtx(BaseCheck.HIGH, self.section_titles["5"])

            for aux_coord in coordinates.split():
                valid_aux_coords.assert_true(
                    aux_coord in ds.variables,
                    f"{name}'s auxiliary coordinate specified by the coordinates attribute, {aux_coord}, "
                    "is not a variable in this dataset"
                    "",
                )
                if aux_coord not in ds.variables:
                    continue

                # TODO CONFORMANCE: Partial implementation of labels
                # §6.1 Allows for "labels" to be referenced as coordinates
                if (
                    hasattr(ds.variables[aux_coord].dtype, "char")
                    and ds.variables[aux_coord].dtype.char == "S"
                ):
                    continue
                elif ds.variables[aux_coord].dtype == str:
                    continue

                aux_coord_dims = set(ds.variables[aux_coord].dimensions)
                valid_aux_coords.assert_true(
                    aux_coord_dims.issubset(dim_set),
                    "dimensions for auxiliary coordinate variable {} ({}) "
                    "are not a subset of dimensions for variable {} ({})"
                    "".format(
                        aux_coord,
                        ", ".join(aux_coord_dims),
                        name,
                        ", ".join(dim_set),
                    ),
                )
            ret_val.append(valid_aux_coords.to_result())
        return ret_val

    # IMPLEMENTATION Section 5 Coordinate Systems and Domain
    def check_coordinates_attribute_format(self, ds):
        """
        Checks that the `coordinates` attribute is a space-separated list of valid variable names,
        and all listed variables exist in the dataset.
        """
        ret_val = []

        geophysical_variables = self._find_geophysical_vars(ds)

        for var_name in geophysical_variables:
            var = ds.variables[var_name]
            coord_attr = getattr(var, "coordinates", None)

            if not coord_attr:
                continue  # No coordinates attribute, nothing to check

            check_coords_attrs_format = TestCtx(
                BaseCheck.HIGH,
                self.section_titles["5"],
            )

            # Check that it is a proper string
            check_coords_attrs_format.assert_true(
                isinstance(coord_attr, str),
                f"The 'coordinates' attribute of variable '{var_name}' must be a string.",
            )

            # Check for unexpected characters (e.g., commas or brackets)
            if isinstance(coord_attr, str):
                has_invalid_chars = any(
                    char in coord_attr for char in [",", "[", "]", "(", ")"]
                )
                check_coords_attrs_format.assert_true(
                    not has_invalid_chars,
                    f"The 'coordinates' attribute of variable '{var_name}' contains invalid characters: {coord_attr}",
                )

                # Check that each token is a valid variable name
                for token in coord_attr.split():
                    check_coords_attrs_format.assert_true(
                        token in ds.variables,
                        f"The 'coordinates' attribute of variable '{var_name}' references non-existent variable '{token}'.",
                    )

            ret_val.append(check_coords_attrs_format.to_result())

        return ret_val

    # IMPLEMENTATION Section 5.1 Independent Latitude, Longitude, Vertical, and Time Axes,
    def check_spatiotemporal_dims_have_coordinate_vars(self, ds):
        """
        Checks that spatial/temporal dimensions (time, lat, lon, height, ...)
        used in geophysical variables have proper coordinate variables.
        """

        # CF-recognized standard names for coordinate axes
        expected_standard_names = {
            "time": "time",
            "lat": "latitude",
            "latitude": "latitude",
            "lon": "longitude",
            "longitude": "longitude",
            "height": "height",
            "depth": "depth",
            "altitude": "altitude",
            "pressure": "air_pressure",
        }

        ret_val = []

        geophysical_variables = self._find_geophysical_vars(ds)
        for var_name in geophysical_variables:
            var = ds.variables[var_name]
            check_spatiotemporal_dims_coords = TestCtx(
                BaseCheck.HIGH,
                self.section_titles["5.1"],
            )

            for dim in var.dimensions:
                if dim in expected_standard_names:
                    if dim not in ds.variables:
                        check_spatiotemporal_dims_coords.assert_true(
                            False,
                            f"Dimension '{dim}' in variable '{var_name}' is expected to be a coordinate axis "
                            f"but no variable with that name exists.",
                        )
                        continue

                    coord_var = ds.variables[dim]
                    std_name = getattr(coord_var, "standard_name", None)

                    check_spatiotemporal_dims_coords.assert_true(
                        std_name == expected_standard_names[dim],
                        f"Coordinate variable '{dim}' should have standard_name='{expected_standard_names[dim]}', "
                        f"found: '{std_name}'",
                    )

            ret_val.append(check_spatiotemporal_dims_coords.to_result())
        return ret_val

    # IMPLEMENTATION Section 2.5.1 Coordinate Systems and Domain
    def check_invalid_coordinate_attr(self, ds):
        """
        Checks that a coordinate variable must not have the _FillValue or missing_value attributes.
        """
        ret_val = []
        for coord_var_name in cfutil.get_coordinate_variables(ds):
            coord_var = ds.variables[coord_var_name]
            valid_coords_attr = TestCtx(BaseCheck.HIGH, self.section_titles["2.5.1"])

            valid_coords_attr.assert_true(
                "_FillValue" not in coord_var.ncattrs(),
                f"The coordinate variable '{coord_var_name}' must not have the _FillValue attribute.",
            )

            valid_coords_attr.assert_true(
                "missing_value" not in coord_var.ncattrs(),
                f"The coordinate variable '{coord_var_name}' must not have the missing_value attribute.",
            )

            ret_val.append(valid_coords_attr.to_result())

        return ret_val

    def check_duplicate_axis(self, ds):
        """
        Checks that no variable contains two coordinates defining the same

        Chapter 5 paragraph 6

        If an axis attribute is attached to an auxiliary coordinate variable,
        it can be used by applications in the same way the `axis` attribute
        attached to a coordinate variable is used. However, it is not
        permissible for a [geophysical variable] to have both a coordinate
        variable and an auxiliary coordinate variable, or more than one of
        either type of variable, having an `axis` attribute with any given
        value e.g. there must be no more than one axis attribute for X for any
        [geophysical variable].

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        :return: List of results
        """

        ret_val = []
        geophysical_variables = self._find_geophysical_vars(ds)
        for name in geophysical_variables:
            no_duplicates = TestCtx(BaseCheck.HIGH, self.section_titles["5"])
            axis_map = cfutil.get_axis_map(ds, name)
            # For every coordinate associated with this variable, keep track of
            # which coordinates define an axis and assert that there are no
            # duplicate axis attributes defined in the set of associated
            # coordinates. axis_map includes coordinates that don't actually have
            # an axis attribute, so we need to ignore those here.
            for axis, coords in axis_map.items():
                coords = [c for c in coords if hasattr(ds.variables[c], "axis")]
                no_duplicates.assert_true(
                    len(coords) <= 1,
                    "'{}' has duplicate axis {} defined by [{}]".format(
                        name,
                        axis,
                        ", ".join(sorted(coords)),
                    ),
                )

            ret_val.append(no_duplicates.to_result())

        return ret_val

    def check_multi_dimensional_coords(self, ds):
        """
        Checks that no multidimensional coordinate shares a name with its
        dimensions.

        Chapter 5 paragraph 4

        We recommend that the name of a [multidimensional coordinate] should
        not match the name of any of its dimensions.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []

        # This can only apply to auxiliary coordinate variables
        for coord in self._find_aux_coord_vars(ds):
            variable = ds.variables[coord]
            if variable.ndim < 2:
                continue
            not_matching = TestCtx(BaseCheck.MEDIUM, self.section_titles["5"])

            not_matching.assert_true(
                coord not in variable.dimensions,
                f"{coord} shares the same name as one of its dimensions" "",
            )
            ret_val.append(not_matching.to_result())

        return ret_val

    # NOTE **********
    # IS THIS EVEN NEEDED ANYMORE?
    # ***************
    def check_grid_coordinates(self, ds):
        """
        5.6 When the coordinate variables for a horizontal grid are not
        longitude and latitude, it is required that the true latitude and
        longitude coordinates be supplied via the coordinates attribute.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        latitudes = cfutil.get_true_latitude_variables(ds)
        longitudes = cfutil.get_true_longitude_variables(ds)

        check_features = [
            "2d-regular-grid",
            "2d-static-grid",
            "3d-regular-grid",
            "3d-static-grid",
            "mapped-grid",
            "reduced-grid",
        ]

        # This one is tricky because there's a very subtle difference between
        # latitude as defined in Chapter 4 and "true" latitude as defined in
        # chapter 5.

        # For each geophysical variable that defines a grid, assert it is
        # associated with a true latitude or longitude coordinate.

        for variable in self._find_geophysical_vars(ds):
            # We use a set so we can do set-wise comparisons with coordinate
            # dimensions
            dimensions = set(ds.variables[variable].dimensions)
            # If it's not a grid, skip it
            if cfutil.guess_feature_type(ds, variable) not in check_features:
                continue
            has_coords = TestCtx(BaseCheck.HIGH, self.section_titles["5.6"])

            # axis_map is a defaultdict(list) mapping the axis to a list of
            # coordinate names. For example:
            # {'X': ['lon'], 'Y':['lat'], 'Z':['lev']}
            # The mapping comes from the dimensions of the variable and the
            # contents of the `coordinates` attribute only.
            axis_map = cfutil.get_axis_map(ds, variable)

            msg = (
                '{}\'s coordinate variable "{}" is not one of the variables identifying true '
                + "latitude/longitude and its dimensions are not a subset of {}'s dimensions"
            )

            alt = (
                "{} has no coordinate associated with a variable identified as true latitude/longitude; "
                "its coordinate variable should also share a subset of {}'s dimensions"
            )

            # Make sure we can find latitude and its dimensions are a subset
            _lat = None
            found_lat = False
            for lat in axis_map["Y"]:
                _lat = lat
                is_subset_dims = set(ds.variables[lat].dimensions).issubset(dimensions)

                if is_subset_dims and lat in latitudes:
                    found_lat = True
                    break
            if _lat:
                has_coords.assert_true(found_lat, msg.format(variable, _lat, variable))
            else:
                has_coords.assert_true(found_lat, alt.format(variable, variable))

            # Make sure we can find longitude and its dimensions are a subset
            _lon = None
            found_lon = False
            for lon in axis_map["X"]:
                _lon = lon
                is_subset_dims = set(ds.variables[lon].dimensions).issubset(dimensions)

                if is_subset_dims and lon in longitudes:
                    found_lon = True
                    break
            if _lon:
                has_coords.assert_true(found_lon, msg.format(variable, _lon, variable))
            else:
                has_coords.assert_true(found_lon, alt.format(variable, variable))

            ret_val.append(has_coords.to_result())
        return ret_val

    def check_reduced_horizontal_grid(self, ds):
        """
        5.3 A "reduced" longitude-latitude grid is one in which the points are
        arranged along constant latitude lines with the number of points on a
        latitude line decreasing toward the poles.

        Recommend that this type of gridded data be stored using the compression
        scheme described in Section 8.2, "Compression by Gathering". The
        compressed latitude and longitude auxiliary coordinate variables are
        identified by the coordinates attribute.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        # Create a set of coordinate variables defining `compress`
        lats = set(cfutil.get_latitude_variables(ds))
        lons = set(cfutil.get_longitude_variables(ds))

        for name in self._find_geophysical_vars(ds):
            coords = getattr(ds.variables[name], "coordinates", None)
            axis_map = cfutil.get_axis_map(ds, name)
            # If this variable has no coordinate that defines compression
            if "C" not in axis_map:
                continue

            valid_rgrid = TestCtx(BaseCheck.HIGH, self.section_titles["5.3"])
            # Make sure reduced grid features define coordinates
            valid_rgrid.assert_true(
                isinstance(coords, str) and coords,
                f"reduced grid feature {name} must define coordinates attribute" "",
            )
            # We can't check anything else if there are no defined coordinates
            if not isinstance(coords, str) and coords:
                continue

            coord_set = set(coords.split())

            # Make sure it's associated with valid lat and valid lon
            valid_rgrid.assert_true(
                len(coord_set.intersection(lons)) > 0,
                f"{name} must be associated with a valid longitude coordinate",
            )
            valid_rgrid.assert_true(
                len(coord_set.intersection(lats)) > 0,
                f"{name} must be associated with a valid latitude coordinate",
            )
            valid_rgrid.assert_true(
                len(axis_map["C"]) == 1,
                "{} can not be associated with more than one compressed coordinates: "
                "({})".format(name, ", ".join(axis_map["C"])),
            )

            for compressed_coord in axis_map["C"]:
                coord = ds.variables[compressed_coord]
                compress = getattr(coord, "compress", None)
                valid_rgrid.assert_true(
                    isinstance(compress, str) and compress,
                    f"compress attribute for compression coordinate {compressed_coord} must be a non-empty string"
                    "",
                )
                if not isinstance(compress, str):
                    continue
                for dim in compress.split():
                    valid_rgrid.assert_true(
                        dim in ds.dimensions,
                        f"dimension {dim} referenced by {compressed_coord}:compress must exist"
                        "",
                    )
            ret_val.append(valid_rgrid.to_result())

        return ret_val

    def _check_grid_mapping_attr_condition(self, attr, attr_name):
        """
        Evaluate a condition (or series of conditions) for a particular
        attribute. Implementation for CF-1.6.

        :param attr: attribute to teset condition for
        :param str attr_name: name of the attribute
        :rtype tuple
        :return two-tuple of (bool, str)
        """

        if attr_name == "latitude_of_projection_origin":
            return self._evaluate_latitude_of_projection_origin(attr)

        elif attr_name == "longitude_of_projection_origin":
            return self._evaluate_longitude_of_projection_origin(attr)

        elif attr_name == "longitude_of_central_meridian":
            return self._evaluate_longitude_of_central_meridian(attr)

        elif attr_name == "longitude_of_prime_meridian":
            return self._evaluate_longitude_of_prime_meridian(attr)

        elif attr_name == "scale_factor_at_central_meridian":
            return self._evaluate_scale_factor_at_central_meridian(attr)

        elif attr_name == "scale_factor_at_projection_origin":
            return self._evaluate_scale_factor_at_projection_origin(attr)

        elif attr_name == "standard_parallel":
            return self._evaluate_standard_parallel(attr)

        elif attr_name == "straight_vertical_longitude_from_pole":
            return self._evaluate_straight_vertical_longitude_from_pole(attr)

        else:
            raise NotImplementedError(
                f"Evaluation for {attr_name} not yet implemented",
            )

    def _evaluate_latitude_of_projection_origin(self, val):
        """
        Evaluate the condition for `latitude_of_projection_origin` attribute.
        Return result. Value must be -90 <= x <= 90.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple (bool, msg)
        """

        return (
            (val >= -90.0) and (val <= 90.0),
            "latitude_of_projection_origin must satisfy (-90 <= x <= 90)",
        )

    def _evaluate_longitude_of_projection_origin(self, val):
        """
        Evaluate the condition for `longitude_of_projection_origin` attribute.
        Return result.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple (bool, msg)
        """

        return (
            (val >= -180.0) and (val <= 180.0),
            "longitude_of_projection_origin must satisfy (-180 <= x <= 180)",
        )

    def _evaluate_longitude_of_central_meridian(self, val):
        """
        Evaluate the condition for `longitude_of_central_meridian` attribute.
        Return result.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple (bool, msg)
        """

        return (
            (val >= -180.0) and (val <= 180.0),
            "longitude_of_central_meridian must satisfy (-180 <= x <= 180)",
        )

    def _evaluate_longitude_of_prime_meridian(self, val):
        """
        Evaluate the condition for `longitude_of_prime_meridian` attribute.
        Return result.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple (bool, msg)
        """

        return (
            (val >= -180.0) and (val <= 180.0),
            "longitude_of_prime_meridian must satisfy (-180 <= x <= 180)",
        )

    def _evaluate_scale_factor_at_central_meridian(self, val):
        """
        Evaluate the condition for `scale_factor_at_central_meridian` attribute.
        Return result.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple (bool, msg)
        """

        return (val > 0.0, "scale_factor_at_central_meridian must be > 0.0")

    def _evaluate_scale_factor_at_projection_origin(self, val):
        """
        Evaluate the condition for `scale_factor_at_projection_origin` attribute.
        Return result.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple (bool, msg)
        """

        return (val > 0.0, "scale_factor_at_projection_origin must be > 0.0")

    def _evaluate_standard_parallel(self, val):
        """
        Evaluate the condition for `standard_parallel` attribute. Return result.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple (bool, msg)
        """

        return (
            (val >= -90.0) and (val <= 90),
            "standard_parallel must satisfy (-90 <= x <= 90)",
        )

    def _evaluate_straight_vertical_longitude_from_pole(self, val):
        """
        Evaluate the condition for `straight_vertical_longitude_from_pole`
        attribute. Return result.

        :param val: value to be tested
        :rtype tuple
        :return two-tuple (bool, msg)
        """

        return (
            (val >= -180.0) and (val <= 180),
            "straight_vertical_longitude_from_pole must satisfy (-180 <= x <= 180)",
        )

    ###############################################################################
    # Chapter 6: Labels and Alternative Coordinates
    ###############################################################################

    def check_geographic_region(self, ds):
        """
        6.1.1 When data is representative of geographic regions which can be identified by names but which have complex
        boundaries that cannot practically be specified using longitude and latitude boundary coordinates, a labeled
        axis should be used to identify the regions.

        Recommend that the names be chosen from the list of standardized region names whenever possible. To indicate
        that the label values are standardized the variable that contains the labels must be given the standard_name
        attribute with the value region.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        region_list = (
            [  # TODO maybe move this (and other info like it) into a config file?
                "africa",
                "antarctica",
                "arabian_sea",
                "aral_sea",
                "arctic_ocean",
                "asia",
                "atlantic_ocean",
                "australia",
                "baltic_sea",
                "barents_opening",
                "barents_sea",
                "beaufort_sea",
                "bellingshausen_sea",
                "bering_sea",
                "bering_strait",
                "black_sea",
                "canadian_archipelago",
                "caribbean_sea",
                "caspian_sea",
                "central_america",
                "chukchi_sea",
                "contiguous_united_states",
                "denmark_strait",
                "drake_passage",
                "east_china_sea",
                "english_channel",
                "eurasia",
                "europe",
                "faroe_scotland_channel",
                "florida_bahamas_strait",
                "fram_strait",
                "global",
                "global_land",
                "global_ocean",
                "great_lakes",
                "greenland",
                "gulf_of_alaska",
                "gulf_of_mexico",
                "hudson_bay",
                "iceland_faroe_channel",
                "indian_ocean",
                "indonesian_throughflow",
                "indo_pacific_ocean",
                "irish_sea",
                "lake_baykal",
                "lake_chad",
                "lake_malawi",
                "lake_tanganyika",
                "lake_victoria",
                "mediterranean_sea",
                "mozambique_channel",
                "north_america",
                "north_sea",
                "norwegian_sea",
                "pacific_equatorial_undercurrent",
                "pacific_ocean",
                "persian_gulf",
                "red_sea",
                "ross_sea",
                "sea_of_japan",
                "sea_of_okhotsk",
                "south_america",
                "south_china_sea",
                "southern_ocean",
                "taiwan_luzon_straits",
                "weddell_sea",
                "windward_passage",
                "yellow_sea",
            ]
        )

        for var in ds.get_variables_by_attributes(standard_name="region"):
            valid_region = TestCtx(BaseCheck.MEDIUM, self.section_titles["6.1"])
            region = var[:]
            if np.ma.isMA(region):
                region = region.data
            valid_region.assert_true(
                "".join(region.astype(str)).lower() in region_list,
                "6.1.1 '{}' specified by '{}' is not a valid region".format(
                    "".join(region.astype(str)),
                    var.name,
                ),
            )
            ret_val.append(valid_region.to_result())
        return ret_val

    ###############################################################################
    # Chapter 7: Data Representative of Cells
    ###############################################################################

    def check_cell_boundaries(self, ds):
        """
        Checks the dimensions of cell boundary variables to ensure they are CF compliant.

        7.1 To represent cells we add the attribute bounds to the appropriate coordinate variable(s). The value of bounds
        is the name of the variable that contains the vertices of the cell boundaries. We refer to this type of variable as
        a "boundary variable." A boundary variable will have one more dimension than its associated coordinate or auxiliary
        coordinate variable. The additional dimension should be the most rapidly varying one, and its size is the maximum
        number of cell vertices.

        Applications that process cell boundary data often times need to determine whether or not adjacent cells share an
        edge. In order to facilitate this type of processing the following restrictions are placed on the data in boundary
        variables:

        Bounds for 1-D coordinate variables

            For a coordinate variable such as lat(lat) with associated boundary variable latbnd(x,2), the interval endpoints
            must be ordered consistently with the associated coordinate, e.g., for an increasing coordinate, lat(1) > lat(0)
            implies latbnd(i,1) >= latbnd(i,0) for all i

            If adjacent intervals are contiguous, the shared endpoint must be represented identically in each instance where
            it occurs in the boundary variable. For example, if the intervals that contain grid points lat(i) and lat(i+1) are
            contiguous, then latbnd(i+1,0) = latbnd(i,1).

        Bounds for 2-D coordinate variables with 4-sided cells

            In the case where the horizontal grid is described by two-dimensional auxiliary coordinate variables in latitude
            lat(n,m) and longitude lon(n,m), and the associated cells are four-sided, then the boundary variables are given
            in the form latbnd(n,m,4) and lonbnd(n,m,4), where the trailing index runs over the four vertices of the cells.

        Bounds for multi-dimensional coordinate variables with p-sided cells

            In all other cases, the bounds should be dimensioned (...,n,p), where (...,n) are the dimensions of the auxiliary
            coordinate variables, and p the number of vertices of the cells. The vertices must be traversed anticlockwise in the
            lon-lat plane as viewed from above. The starting vertex is not specified.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
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
            if boundary_variable_name not in ds.variables:
                valid = False
                reasoning.append(
                    f"Boundary variable {boundary_variable_name} referenced by {variable.name} not "
                    + "found in dataset variables",
                )
            else:
                boundary_variable = ds.variables[boundary_variable_name]
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

            # ensure p vertices form a valid simplex given previous a...n
            # previous auxiliary coordinates
            if (
                ds.dimensions[boundary_variable.dimensions[-1]].size
                < len(boundary_variable.dimensions[:-1]) + 1
            ):
                valid = False
                reasoning.append(
                    f"Dimension {boundary_variable.name} of boundary variable (for {variable.name}) must have at least {len(variable.dimensions) + 1} elements to form a simplex/closed cell with previous dimensions {boundary_variable.dimensions[:-1]}.",
                )
            result = Result(
                BaseCheck.MEDIUM,
                valid,
                self.section_titles["7.1"],
                reasoning,
            )
            ret_val.append(result)
        return ret_val

    def _cell_measures_core(self, ds, var, external_set, variable_template):
        # IMPLEMENTATION CONFORMANCE REQUIRED 1/2
        reasoning = []
        search_str = (
            r"^(?P<measure_type>area|volume):\s+(?P<cell_measure_var_name>\w+)$"
        )
        search_res = regex.match(search_str, var.cell_measures)
        if not search_res:
            valid = False
            reasoning.append(
                f"The cell_measures attribute for variable {var.name} "
                "is formatted incorrectly. It should take the "
                "form of either 'area: cell_var' or "
                "'volume: cell_var' where cell_var is an existing name of "
                "a variable describing the cell measures.",
            )
        else:
            valid = True
            cell_measure_var_name = search_res.group("cell_measure_var_name")
            cell_measure_type = search_res.group("measure_type")
            # TODO: cache previous results
            if cell_measure_var_name not in set(ds.variables.keys()).union(
                external_set,
            ):
                valid = False
                reasoning.append(
                    f"Cell measure variable {cell_measure_var_name} referred to by "
                    f"{var.name} is not present in {variable_template}s".format(
                        cell_measure_var_name,
                        var.name,
                    ),
                )
            # CF 1.7+ assume external variables -- further checks can't be run here
            elif cell_measure_var_name in external_set:
                # can't test anything on an external var
                return Result(
                    BaseCheck.MEDIUM,
                    valid,
                    (self.section_titles["7.2"]),
                    reasoning,
                )

            else:
                cell_measure_var = ds.variables[cell_measure_var_name]
                if not hasattr(cell_measure_var, "units"):
                    valid = False
                    reasoning.append(
                        f"Cell measure variable {cell_measure_var_name} is required "
                        "to have units attribute defined",
                    )
                else:
                    # IMPLEMENTATION CONFORMANCE REQUIRED 2/2
                    # verify this combination {area: 'm2', volume: 'm3'}

                    # key is valid measure types, value is expected
                    # exponent
                    exponent_lookup = {"area": 2, "volume": 3}
                    exponent = exponent_lookup[search_res.group("measure_type")]
                    conversion_failure_msg = (
                        f'Variable "{cell_measure_var.name}" must have units which are convertible '
                        f'to UDUNITS "m{exponent}" when variable is referred to by a {variable_template} with '
                        f'cell_methods attribute with a measure type of "{cell_measure_type}".'
                    )
                    try:
                        cell_measure_units = Unit(cell_measure_var.units)
                    except ValueError:
                        valid = False
                        reasoning.append(conversion_failure_msg)
                    else:
                        if not cell_measure_units.is_convertible(Unit(f"m{exponent}")):
                            valid = False
                            reasoning.append(conversion_failure_msg)
                    if not set(cell_measure_var.dimensions).issubset(var.dimensions):
                        valid = False
                        reasoning.append(
                            f"Cell measure variable {cell_measure_var_name} must have "
                            "dimensions which are a subset of "
                            f"those defined in variable {var.name}.",
                        )

        return Result(BaseCheck.MEDIUM, valid, (self.section_titles["7.2"]), reasoning)

    def check_cell_measures(self, ds):
        """
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
        for var in variables:
            result = self._cell_measures_core(ds, var, set(), "dataset variable")
            ret_val.append(result)

        return ret_val

    def check_cell_methods(self, ds):
        """
        7.3 To describe the characteristic of a field that is represented by cell values, we define the cell_methods attribute
        of the variable. This is a string attribute comprising a list of blank-separated words of the form "name: method". Each
        "name: method" pair indicates that for an axis identified by name, the cell values representing the field have been
        determined or derived by the specified method.

        name can be a dimension of the variable, a scalar coordinate variable, a valid standard name, or the word "area"

        values of method should be selected from the list in Appendix E, Cell Methods, which includes point, sum, mean, maximum,
        minimum, mid_range, standard_deviation, variance, mode, and median. Case is not significant in the method name. Some
        methods (e.g., variance) imply a change of units of the variable, as is indicated in Appendix E, Cell Methods.

        Because the default interpretation for an intensive quantity differs from that of an extensive quantity and because this
        distinction may not be understood by some users of the data, it is recommended that every data variable include for each
        of its dimensions and each of its scalar coordinate variables the cell_methods information of interest (unless this
        information would not be meaningful). It is especially recommended that cell_methods be explicitly specified for each
        spatio-temporal dimension and each spatio-temporal scalar coordinate variable.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """

        ret_val = []
        # CONFORMANCE IMPLEMENTATION 7.3 1/3
        psep = regex.compile(
            r"(?P<vars>\w+: )+(?P<method>\w+) ?(?P<where>where (?P<wtypevar>\w+) "
            r"?(?P<over>over (?P<otypevar>\w+))?| ?)(?:\((?P<paren_contents>[^)]*)\))?",
        )

        for var in ds.get_variables_by_attributes(cell_methods=lambda x: x is not None):
            if not getattr(var, "cell_methods", ""):
                continue

            method = getattr(var, "cell_methods", "")

            valid_attribute = TestCtx(
                BaseCheck.HIGH,
                self.section_titles["7.3"],
            )  # changed from 7.1 to 7.3
            valid_attribute.assert_true(
                regex.match(psep, method) is not None,
                f'"{method}" is not a valid format for cell_methods attribute of "{var.name}"'
                "",
            )
            ret_val.append(valid_attribute.to_result())

            valid_cell_names = TestCtx(BaseCheck.MEDIUM, self.section_titles["7.3"])

            # check that the name is valid
            for match in regex.finditer(psep, method):
                # it is possible to have "var1: var2: ... varn: ...", so handle
                # that case
                for var_raw_str in match.captures("vars"):
                    # strip off the ' :' at the end of each match
                    var_str = var_raw_str[:-2]
                    if (
                        var_str in var.dimensions
                        or var_str == "area"
                        or var_str in getattr(var, "coordinates", "")
                    ):
                        valid = True
                    else:
                        valid = False

                    valid_cell_names.assert_true(
                        valid,
                        f"{var.name}'s cell_methods name component {var_str} does not match a dimension, "
                        "area or auxiliary coordinate",
                    )

            ret_val.append(valid_cell_names.to_result())

            # Checks if the method value of the 'name: method' pair is acceptable
            valid_cell_methods = TestCtx(BaseCheck.MEDIUM, self.section_titles["7.3"])

            for match in regex.finditer(psep, method):
                # CF section 7.3 - "Case is not significant in the method name."
                valid_cell_methods.assert_true(
                    match.group("method").lower() in self.cell_methods,
                    "{}:cell_methods contains an invalid method: {}"
                    "".format(var.name, match.group("method")),
                )

            ret_val.append(valid_cell_methods.to_result())

            for match in regex.finditer(psep, method):
                if match.group("paren_contents") is not None:
                    # split along spaces followed by words with a colon
                    # not sure what to do if a comment contains a colon!
                    ret_val.append(
                        self._check_cell_methods_paren_info(
                            match.group("paren_contents"),
                            var,
                        ).to_result(),
                    )

        return ret_val

    def _check_cell_methods_paren_info(self, paren_contents, var):
        """
        Checks that the spacing and/or comment info contained inside the
        parentheses in cell_methods is well-formed
        """
        # IMPLEMENTATION CONFORMANCE REQUIRED 3/3 - comment/paren contents
        valid_info = TestCtx(BaseCheck.MEDIUM, self.section_titles["7.3"])
        # if there are no colons, this is a simple comment
        # TODO: are empty comments considered valid?
        if ":" not in paren_contents:
            valid_info.out_of += 1
            valid_info.score += 1
            return valid_info
        # otherwise, split into k/v pairs
        kv_pair_pat = r"(\S+:)\s+(.*(?=\s+\w+:)|[^:]+$)\s*"
        # otherwise, we must split further with intervals coming
        # first, followed by non-standard comments
        # we need the count of the matches, and re.findall() only returns
        # groups if they are present and we wish to see if the entire match
        # object concatenated together is the same as the original string
        pmatches = list(regex.finditer(kv_pair_pat, paren_contents))
        for i, pmatch in enumerate(pmatches):
            keyword, val = pmatch.groups()
            if keyword == "interval:":
                valid_info.out_of += 2
                interval_matches = regex.match(
                    r"^\s*(?P<interval_number>\S+)\s+(?P<interval_units>\S+)\s*$",
                    val,
                )
                # attempt to get the number for the interval
                if not interval_matches:
                    valid_info.messages.append(
                        f'§7.3.3 {var.name}:cell_methods contains an interval specification that does not parse: "{val}". Should be in format "interval: <number> <units>"',
                    )
                else:
                    try:
                        float(interval_matches.group("interval_number"))
                    except ValueError:
                        valid_info.messages.append(
                            '§7.3.3 {}:cell_methods contains an interval value that does not parse as a numeric value: "{}".'.format(
                                var.name,
                                interval_matches.group("interval_number"),
                            ),
                        )
                    else:
                        valid_info.score += 1

                    # then the units
                    try:
                        Unit(interval_matches.group("interval_units"))
                    except ValueError:
                        valid_info.messages.append(
                            '§7.3.3 {}:cell_methods interval units "{}" is not parsable by UDUNITS.'.format(
                                var.name,
                                interval_matches.group("interval_units"),
                            ),
                        )
                    else:
                        valid_info.score += 1
            elif keyword == "comment:":
                # comments can't really be invalid, except
                # if they come first or aren't last, and
                # maybe if they contain colons embedded in the
                # comment string
                valid_info.out_of += 1
                if len(pmatches) == 1:
                    valid_info.messages.append(
                        f"§7.3.3 If there is no standardized information, the keyword comment: should be omitted for variable {var.name}",
                    )
                # otherwise check that the comment is the last
                # item in the parentheses
                elif i != len(pmatches) - 1:
                    valid_info.messages.append(
                        f'§7.3.3 The non-standard "comment:" element must come after any standard elements in cell_methods for variable {var.name}',
                    )
                #
                else:
                    valid_info.score += 1
            else:
                valid_info.out_of += 1
                valid_info.messages.append(
                    f'§7.3.3 Invalid cell_methods keyword "{keyword}" for variable {var.name}. Must be one of [interval, comment]',
                )

        # Ensure concatenated reconstructed matches are the same as the
        # original string.  If they're not, there's likely a formatting error
        valid_info.assert_true(
            "".join(m.group(0) for m in pmatches) == paren_contents,
            f"§7.3.3 Parenthetical content inside {var.name}:cell_methods is not well formed: {paren_contents}",
        )

        return valid_info

    def check_climatological_statistics(self, ds):
        """
        7.4 A climatological time coordinate variable does not have a bounds attribute. Instead, it has a climatology
        attribute, which names a variable with dimensions (n,2), n being the dimension of the climatological time axis.
        Using the units and calendar of the time coordinate variable, element (i,0) of the climatology variable specifies
        the beginning of the first subinterval and element (i,1) the end of the last subinterval used to evaluate the
        climatological statistics with index i in the time dimension. The time coordinates should be values that are
        representative of the climatological time intervals, such that an application which does not recognise climatological
        time will nonetheless be able to make a reasonable interpretation.

        A climatological axis may use different statistical methods to measure variation among years, within years, and within
        days. The methods which can be specified are those listed in Appendix E, Cell Methods and each entry in the cell_methods
        attribute may also contain non-standardised information in parentheses after the method. The value of the cell_method
        attribute must be in one of the following forms:
        - time: method1 within years   time: method2 over years
        - time: method1 within days    time: method2 over days
        - time: method1 within days    time: method2 over days   time: method3 over years

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """

        reasoning = []
        ret_val = []
        total_climate_count = 0
        valid_climate_count = 0
        all_clim_coord_var_names = []

        methods = [
            "point",  # TODO change to appendix import once cf1.7 merged
            "sum",
            "mean",
            "maximum",
            "minimum",
            "mid_range",
            "standard_deviation",
            "variance",
            "mode",
            "median",
        ]

        # find any climatology axis variables; any variables which contain climatological stats will use
        # these variables as coordinates
        clim_time_coord_vars = ds.get_variables_by_attributes(
            climatology=lambda s: s is not None,
        )

        # first, to determine whether or not we have a valid climatological time
        # coordinate variable, we need to make sure it has the attribute "climatology",
        # but not the attribute "bounds"
        time_vars = cfutil.get_time_variables(ds)
        for clim_coord_var in clim_time_coord_vars:
            climatology_ctx = TestCtx(BaseCheck.MEDIUM, self.section_titles["7.3"])
            # IMPLEMENTATION CONFORMANCE 7.4 REQUIRED 1/6
            if clim_coord_var.name not in time_vars:
                climatology_ctx.out_of += 1
                climatology_ctx.messages.append(
                    f"Variable {clim_coord_var.name} is not detected as a time "
                    "coordinate variable, but has climatology attribute",
                )
            # IMPLEMENTATION CONFORMANCE 7.4 REQUIRED
            if hasattr(clim_coord_var, "bounds"):
                climatology_ctx.out_of += 1
                climatology_ctx.messages.append(
                    f"Variable {clim_coord_var.name} has a climatology "
                    "attribute and cannot also have a bounds attribute.",
                )
                result = Result(
                    BaseCheck.MEDIUM,
                    False,
                    (self.section_titles["7.4"]),
                    reasoning,
                )

            # IMPLEMENTATION CONFORMANCE 7.4 REQUIRED 2/6
            # make sure the climatology variable referenced actually exists
            elif not isinstance(clim_coord_var.climatology, str):
                climatology_ctx.out_of += 1
                climatology_ctx.messages.append(
                    f"Variable {clim_coord_var.name} must have a climatology "
                    "attribute which is a string",
                )
                ret_val.append(climatology_ctx.to_result())
                continue
            elif clim_coord_var.climatology not in ds.variables:
                climatology_ctx.out_of += 1
                climatology_ctx.messages.append(
                    "Variable {} referenced in time's climatology attribute does not exist".format(
                        ds.variables["time"].climatology,
                    ),
                )
            else:
                clim_var = ds.variables[clim_coord_var.climatology]
                #
                # IMPLEMENTATION CONFORMANCE 7.4 REQUIRED 4/6
                if clim_var.dtype is str or not np.issubdtype(clim_var, np.number):
                    climatology_ctx.out_of += 1
                    climatology_ctx.messages.append(
                        f"Climatology variable {clim_var.name} is not a numeric type",
                    )
                # IMPLEMENTATION CONFORMANCE 7.4 REQUIRED 6/6
                if hasattr(clim_var, "_FillValue") or hasattr(
                    clim_var,
                    "missing_value",
                ):
                    climatology_ctx.out_of += 1
                    climatology_ctx.messages.append(
                        f"Climatology variable {clim_var.name} may not contain "
                        "attributes _FillValue or missing_value",
                    )

                # IMPLEMENTATION CONFORMANCE 7.4 REQUIRED 5/6
                for same_attr in ("units", "standard_name", "calendar"):
                    if hasattr(clim_var, same_attr):
                        climatology_ctx.assert_true(
                            getattr(clim_var, same_attr)
                            == getattr(clim_coord_var, same_attr, None),
                            f"Attribute {same_attr} must have the same value in both "
                            f"variables {clim_var.name} and {clim_coord_var.name}",
                        )
            ret_val.append(climatology_ctx.to_result())

            # check that coordinate bounds are in the proper order.
            # make sure last elements are boundary variable specific dimensions
            # IMPLEMENTATION CONFORMANCE 7.4 REQUIRED 3/6
            if (
                clim_coord_var.dimensions[:]
                != ds.variables[clim_coord_var.climatology].dimensions[
                    : clim_coord_var.ndim
                ]
            ):
                total_climate_count += 1
                reasoning.append(
                    f"Climatology variable coordinates are in improper order: {ds.variables[clim_coord_var.climatology].dimensions}. Bounds-specific dimensions should be last",
                )
                result = Result(
                    BaseCheck.MEDIUM,
                    (valid_climate_count, total_climate_count),
                    (self.section_titles["7.4"]),
                    reasoning,
                )
                ret_val.append(result)

            # IMPLEMENTATION CONFORMANCE 7.4 REQUIRED 3/6 - dim size of 2 for
            # climatology-specific dimension
            elif (
                ds.dimensions[
                    ds.variables[clim_coord_var.climatology].dimensions[-1]
                ].size
                != 2
            ):
                reasoning.append(
                    f'Climatology dimension "{ds.variables[clim_coord_var.climatology].name}" should only contain two elements',
                )
                total_climate_count += 1
                result = Result(
                    BaseCheck.MEDIUM,
                    (valid_climate_count, total_climate_count),
                    (self.section_titles["7.4"]),
                    reasoning,
                )
                ret_val.append(result)

            # passed all these checks, so we can add this clim_coord_var to our total list
            all_clim_coord_var_names.append(clim_coord_var.name)

        # for any variables which use a climatology time coordinate variable as a coordinate,
        # if they have a cell_methods attribute, it must comply with the form:
        #     time: method1 within years time: method2 over years
        #     time: method1 within days time: method2 over days
        #     time: method1 within days time: method2 over days time: method3 over years
        # optionally followed by parentheses for explaining additional
        # info, e.g.
        # "time: method1 within years time: method2 over years (sidereal years)"

        meth_regex = "(?:{})".format(
            "|".join(methods),
        )  # "or" comparison for the methods
        re_string = (
            rf"^time: {meth_regex} within (years|days)"  # regex string to test
            rf" time: {meth_regex} over \1(:?(?<=days) time: {meth_regex} over years)?"
            r"(?: \([^)]+\))?$"  # parenthesized comment section
        )

        # find any variables with a valid climatological cell_methods
        for cell_method_var in ds.get_variables_by_attributes(
            cell_methods=lambda s: s is not None,
        ):
            if any(
                dim in all_clim_coord_var_names for dim in cell_method_var.dimensions
            ):
                total_climate_count += 1
                if not regex.search(re_string, cell_method_var.cell_methods):
                    reasoning.append(
                        f'The "time: method within years/days over years/days" format is not correct in variable {cell_method_var.name}.',
                    )
                else:
                    valid_climate_count += 1

                result = Result(
                    BaseCheck.MEDIUM,
                    (valid_climate_count, total_climate_count),
                    (self.section_titles["7.4"]),
                    reasoning,
                )
                ret_val.append(result)

        return ret_val

    ###############################################################################
    # Chapter 8: Reduction of Dataset Size
    ###############################################################################

    def check_packed_data(self, ds):
        """
        8.1 Simple packing may be achieved through the use of the optional NUG defined attributes scale_factor and
        add_offset. After the data values of a variable have been read, they are to be multiplied by the scale_factor,
        and have add_offset added to them.

        The units of a variable should be representative of the unpacked data.

        If the scale_factor and add_offset attributes are of the same data type as the associated variable, the unpacked
        data is assumed to be of the same data type as the packed data. However, if the scale_factor and add_offset
        attributes are of a different data type from the variable (containing the packed data) then the unpacked data
        should match the type of these attributes, which must both be of type float or both be of type double. An additional
        restriction in this case is that the variable containing the packed data must be of type byte, short or int. It is
        not advised to unpack an int into a float as there is a potential precision loss.

        When data to be packed contains missing values the attributes that indicate missing values (_FillValue, valid_min,
        valid_max, valid_range) must be of the same data type as the packed data.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        for name, var in ds.variables.items():
            add_offset = getattr(var, "add_offset", None)
            scale_factor = getattr(var, "scale_factor", None)
            if not (add_offset or scale_factor):
                continue

            valid = True
            reasoning = []

            # if only one of these attributes is defined, assume they
            # are the same type (value doesn't matter here)
            if not add_offset:
                add_offset = scale_factor
            if not scale_factor:
                scale_factor = add_offset

            # IMPLEMENTATION CONFORMANCE 8.1 REQUIRED 1/3
            # scale_factor and add_offset same type
            if not isinstance(add_offset, type(scale_factor)):
                valid = False
                reasoning.append(
                    "Attributes add_offset and scale_factor have different data type.",
                )
            # IMPLEMENTATION CONFORMANCE 8.1 REQUIRED 2/3
            # scale_factor and add_offset must be floating point or double
            # if not the same type
            # FIXME: Check add_offset too.

            elif not isinstance(scale_factor, var.dtype.type):
                # Check both attributes are type float or double
                if not isinstance(scale_factor, (float, np.floating)):
                    valid = False
                    reasoning.append(
                        "Attributes add_offset and scale_factor are not of type float or double.",
                    )
                else:
                    # Check variable type is byte, short or int
                    if var.dtype.type not in [
                        int,
                        np.int8,
                        np.int16,
                        np.int32,
                        np.int64,
                    ]:
                        valid = False
                        # IMPLEMENTATION CONFORMANCE REQUIRED 3/3
                        # IMPLEMENTATION CONFORMANCE REQUIRED 3/3
                        reasoning.append(
                            "Variable is not of type byte, short, or int as required for different type add_offset/scale_factor.",
                        )

            result = Result(
                BaseCheck.MEDIUM,
                valid,
                self.section_titles["8.1"],
                reasoning,
            )
            ret_val.append(result)
            reasoning = []

            valid = True
            # test further with  _FillValue , valid_min , valid_max , valid_range
            if hasattr(var, "_FillValue"):
                if var._FillValue.dtype.type != var.dtype.type:
                    valid = False
                    reasoning.append(
                        f"Type of {name}:_FillValue attribute ({var._FillValue.dtype.name}) does not match variable type ({var.dtype.name})",
                    )
            if hasattr(var, "valid_min"):
                if var.valid_min.dtype.type != var.dtype.type:
                    valid = False
                    reasoning.append(
                        f"Type of {name}valid_min attribute ({var.valid_min.dtype.name}) does not match variable type ({var.dtype.name})",
                    )
            if hasattr(var, "valid_max"):
                if var.valid_max.dtype.type != var.dtype.type:
                    valid = False
                    reasoning.append(
                        f"Type of {name}:valid_max attribute ({var.valid_max.dtype.name}) does not match variable type ({var.dtype.name})",
                    )
            if hasattr(var, "valid_range"):
                if var.valid_range.dtype.type != var.dtype.type:
                    valid = False
                    reasoning.append(
                        f"Type of {name}:valid_range attribute ({var.valid_range.dtype.name}) does not match variable type ({var.dtype.name})",
                    )

            result = Result(
                BaseCheck.MEDIUM,
                valid,
                self.section_titles["8.1"],
                reasoning,
            )
            ret_val.append(result)

        return ret_val

    def check_compression_gathering(self, ds):
        """
        At the current time the netCDF interface does not provide for packing
        data. However a simple packing may be achieved through the use of the
        optional NUG defined attributes scale_factor and add_offset . After the
        data values of a variable have been read, they are to be multiplied by
        the scale_factor , and have add_offset added to them. If both
        attributes are present, the data are scaled before the offset is added.
        When scaled data are written, the application should first subtract the
        offset and then divide by the scale factor. The units of a variable
        should be representative of the unpacked data.

        This standard is more restrictive than the NUG with respect to the use
        of the scale_factor and add_offset attributes; ambiguities and
        precision problems related to data type conversions are resolved by
        these restrictions. If the scale_factor and add_offset attributes are
        of the same data type as the associated variable, the unpacked data is
        assumed to be of the same data type as the packed data. However, if the
        scale_factor and add_offset attributes are of a different data type
        from the variable (containing the packed data) then the unpacked data
        should match the type of these attributes, which must both be of type
        float or both be of type double . An additional restriction in this
        case is that the variable containing the packed data must be of type
        byte , short or int . It is not advised to unpack an int into a float
        as there is a potential precision loss.

        When data to be packed contains missing values the attributes that
        indicate missing values ( _FillValue , valid_min , valid_max ,
                                 valid_range ) must be of the same data type as
        the packed data. See Section 2.5.1, “Missing Data” for a discussion of
        how applications should treat variables that have attributes indicating
        both missing values and transformations defined by a scale and/or
        offset.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        for compress_var in ds.get_variables_by_attributes(
            compress=lambda s: s is not None,
        ):
            valid = True
            reasoning = []
            # puts the referenced variable being compressed into a set
            compress_set = set(compress_var.compress.split(" "))
            if compress_var.ndim != 1:
                valid = False
                reasoning.append(
                    f"Compression variable {compress_var.name} may only have one dimension",
                )
            # IMPLEMENTATION CONFORMANCE 8.2 REQUIRED 1/3
            # ensure compression variable is a proper index, and thus is an
            # signed or unsigned integer type of some sort
            if (compress_var.dtype is str) or (
                compress_var.dtype.kind not in {"i", "u"}
            ):
                valid = False
                reasoning.append(
                    f"Compression variable {compress_var.name} must be an integer type to form a proper array index",
                )
            # IMPLEMENTATION CONFORMANCE 8.2 REQUIRED 2/3
            # make sure all the variables referred to are contained by the
            # variables.
            if not compress_set.issubset(ds.dimensions):
                not_in_dims = sorted(compress_set.difference(ds.dimensions))
                valid = False
                reasoning.append(
                    f"The following dimensions referenced by the compress attribute of variable {compress_var.name} do not exist: {not_in_dims}",
                )
            # IMPLEMENTATION CONFORMANCE 8.2 REQUIRED 3/3
            # The values of the associated coordinate variable must be in the range
            # starting with 0 and going up to the product of the compressed dimension
            # sizes minus 1 (CDL index conventions).

            # Put the the values of the associated coordinate variable into a list
            coord_list_size = [
                item.size
                for item in ds.dimensions.values()
                if item.name in compress_set
            ]
            # get the upper limit of the dimenssion size
            upper_limit_size = np.prod(coord_list_size) - 1

            for coord_size in coord_list_size:
                if coord_size not in range(0, upper_limit_size):
                    valid = False
                    reasoning.append(
                        f"The dimenssion size {coord_size} referenced by the compress attribute is not "
                        "in the range (0, The product of the compressed dimension sizes minus 1)",
                    )
            result = Result(
                BaseCheck.MEDIUM,
                valid,
                self.section_titles["8.2"],
                reasoning,
            )
            ret_val.append(result)

        return ret_val

    ###############################################################################
    # Chapter 9: Discrete Sampling Geometries
    ###############################################################################

    def check_feature_type(self, ds):
        """
        Check the global attribute featureType for valid CF featureTypes

        9.4 A global attribute, featureType, is required for all Discrete Geometry representations except the orthogonal
        multidimensional array representation, for which it is highly recommended.

        The value assigned to the featureType attribute is case-insensitive.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """
        # Due to case insensitive requirement, we list the possible featuretypes
        # in lower case and check using the .lower() method
        feature_list = [
            "point",
            "timeseries",
            "trajectory",
            "profile",
            "timeseriesprofile",
            "trajectoryprofile",
        ]

        feature_type = getattr(ds, "featureType", None)
        valid_feature_type = TestCtx(
            BaseCheck.HIGH,
            "§9.1 Dataset contains a valid featureType",
        )
        valid_feature_type.assert_true(
            feature_type is None or feature_type.lower() in feature_list,
            "{} is not a valid CF featureType. It must be one of {}"
            "".format(feature_type, ", ".join(feature_list)),
        )
        return valid_feature_type.to_result()

    def check_cf_role(self, ds):
        """
        Check variables defining cf_role for legal cf_role values.

        §9.5 The only acceptable values of cf_role for Discrete Geometry CF
        data sets are timeseries_id, profile_id, and trajectory_id

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: compliance_checker.base.Result
        """
        valid_roles = ["timeseries_id", "profile_id", "trajectory_id"]
        variable_count = 0
        valid_cf_role = TestCtx(BaseCheck.HIGH, self.section_titles["9.5"])
        for variable in ds.get_variables_by_attributes(cf_role=lambda x: x is not None):
            variable_count += 1
            cf_role = variable.cf_role
            valid_cf_role.assert_true(
                cf_role in valid_roles,
                "{} is not a valid cf_role value. It must be one of {}"
                "".format(cf_role, ", ".join(valid_roles)),
            )
        if variable_count > 0:
            return valid_cf_role.to_result()

    def check_variable_features(self, ds):
        """
        Checks the variable feature types match the dataset featureType attribute.
        If more than one unique feature type is found, report this as an error.

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        feature_types_found = defaultdict(list)
        ret_val = []
        feature_list = {
            "point",
            "timeseries",
            "trajectory",
            "profile",
            "timeseriesprofile",
            "trajectoryprofile",
        }
        # Don't bother checking if it's not a legal featureType
        # if the featureType attribute doesn't exist
        feature_type = getattr(ds, "featureType", "")
        if feature_type is not None and feature_type.lower() not in feature_list:
            return []

        _feature = feature_type.lower()

        for name in self._find_geophysical_vars(ds):
            variable_feature = cfutil.guess_feature_type(ds, name)
            # If we can't figure it out, don't check it.
            if variable_feature is None:
                continue
            feature_types_found[variable_feature].append(name)
            matching_feature = TestCtx(BaseCheck.MEDIUM, self.section_titles["9.1"])
            matching_feature.assert_true(
                variable_feature.lower() == _feature,
                f"{name} is not a {_feature}, it is detected as a {variable_feature}"
                "",
            )
            ret_val.append(matching_feature.to_result())

        # create explanation of all of the different featureTypes
        # found in the dataset
        feature_description = ", ".join(
            [
                "{} ({})".format(ftr, ", ".join(vrs))
                for ftr, vrs in feature_types_found.items()
            ],
        )
        all_same_features = TestCtx(BaseCheck.HIGH, self.section_titles["9.1"])
        all_same_features.assert_true(
            len(feature_types_found) < 2,
            f"Different feature types discovered in this dataset: {feature_description}"
            "",
        )
        ret_val.append(all_same_features.to_result())

        return ret_val

    def check_hints(self, ds):
        """
        Checks for potentially mislabeled metadata and makes suggestions for how to correct

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []

        ret_val.extend(self._check_hint_bounds(ds))

        return ret_val

    def _check_hint_bounds(self, ds):
        """
        Checks for variables ending with _bounds, if they are not cell methods,
        make the recommendation

        :param netCDF4.Dataset ds: An open netCDF dataset
        :rtype: list
        :return: List of results
        """
        ret_val = []
        boundary_variables = cfutil.get_cell_boundary_variables(ds)
        for name in ds.variables:
            if name.endswith("_bounds") and name not in boundary_variables:
                msg = (
                    f"{name} might be a cell boundary variable but there are no variables that define it "
                    "as a boundary using the `bounds` attribute."
                )
                result = Result(BaseCheck.LOW, True, self.section_titles["7.1"], [msg])
                ret_val.append(result)

        return ret_val
