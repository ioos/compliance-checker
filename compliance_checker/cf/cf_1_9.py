import cftime
import numpy as np
import regex
from netCDF4 import Dataset

import compliance_checker.cf.util as cfutil
from compliance_checker.base import BaseCheck, Result, TestCtx
from compliance_checker.cf import util
from compliance_checker.cf.cf_1_8 import CF1_8Check
from compliance_checker.cf.util import (
    VariableReferenceError,
    get_coordinate_variables,
    reference_attr_variables,
)


class CF1_9Check(CF1_8Check):
    _cc_spec_version = "1.9"
    _cc_url = "http://cfconventions.org/Data/cf-conventions/cf-conventions-1.9/cf-conventions.html"
    _allowed_numeric_var_types = CF1_8Check._allowed_numeric_var_types.union(
        {np.ubyte, np.uint16, np.uint32, np.uint64},
    )

    def __init__(self, options=None):
        super().__init__(options)
        self.section_titles.update({"5.8": "§5.8 Domain Variables"})

    def check_time_coordinate_variable_has_calendar(self, ds):
        """
        Ensure that time coordinate variables have a calendar attribute
        """
        ret_val = []
        for name in cfutil.get_time_variables(ds):
            # DRY: get rid of time coordinate variable boilerplate
            if name not in {var.name for var in util.find_coord_vars(ds)}:
                continue
            time_var = ds.variables[name]
            if not hasattr(time_var, "calendar") or not isinstance(
                time_var.calendar,
                str,
            ):
                result = Result(
                    BaseCheck.MEDIUM,
                    True,
                    self.section_titles["4.4.1"],
                    [
                        f'Time coordinate variable "{name}" should have a '
                        'string valued attribute "calendar"',
                    ],
                )
                ret_val.append(result)
                continue
            if time_var.calendar.lower() in {"gregorian", "julian", "standard"}:
                try:
                    reference_year = cftime.num2date(
                        0,
                        time_var.units,
                        time_var.calendar,
                        has_year_zero=True,
                    ).year
                # will fail on months, certain other time specifications
                except ValueError:
                    continue
                if reference_year == 0:
                    reasoning = (
                        f'For time variable "{time_var.name}", when using '
                        "the Gregorian or Julian calendars, the use of year "
                        "zero is not recommended. Furthermore, the use of year "
                        "zero to signify a climatological variable as in COARDS "
                        "is deprecated in CF."
                    )
                    result = Result(
                        BaseCheck.MEDIUM,
                        False,
                        self.section_titles["4.4.1"],
                        [reasoning],
                    )
                elif reference_year < 0:
                    reasoning = (
                        f'For time variable "{time_var.name}", when using '
                        "the Gregorian or Julian calendars, negative "
                        "years are not allowed."
                    )
                    result = Result(
                        BaseCheck.MEDIUM,
                        False,
                        self.section_titles["4.4.1"],
                        [reasoning],
                    )
                else:
                    return ret_val

                ret_val.append(result)
        return ret_val

    def check_time_coordinate(self, ds):
        prev_return = super().check_time_coordinate(ds)
        # adds check for < 60 seconds
        seconds_regex = regex.compile(
            r"\w+ since \d{1,4}-\d{1,2}-\d{1,2}[ T]"
            r"\d{1,2}:\d{1,2}:(?P<seconds>\d{1,2})",
        )
        for name in cfutil.get_time_variables(ds):
            # DRY: get rid of time coordinate variable boilerplate
            if name not in {var.name for var in util.find_coord_vars(ds)}:
                continue
            time_var = ds.variables[name]
            test_ctx = self.get_test_ctx(
                BaseCheck.HIGH,
                self.section_titles["4.4"],
                name,
            )
            try:
                match = regex.match(seconds_regex, time_var.units)
            except AttributeError:
                # not much can be done if there are no units
                continue
            # only check if seconds are present and in a parseable format
            if match:
                test_ctx.assert_true(
                    match.group("seconds") is None or int(match.group("seconds")) < 60,
                    f'Time coordinate variable "{name}" must have '
                    "units with seconds less than 60",
                )
                prev_return.append(test_ctx.to_result())
        return prev_return

    def check_domain_variables(self, ds: Dataset):

        results = []
        domain_valid = TestCtx(BaseCheck.MEDIUM, self.section_titles["5.8"])
        is_ragged_array_repr = cfutil.is_dataset_valid_ragged_array_repr_featureType(
            ds,
            getattr(ds, "featureType", ""),
        )
        if is_ragged_array_repr:
            domain_valid.out_of += 1
            ragged_array_dim_variable, ragged_attr_name = (
                cfutil.resolve_ragged_array_dimension(ds)
            )
            dim_name = getattr(ragged_array_dim_variable, ragged_attr_name)
            referenced_dim = reference_attr_variables(
                ds,
                dim_name,
                reference_type="dimensions",
            )
            if isinstance(referenced_dim, VariableReferenceError):
                domain_valid.messages.append(
                    f"Found ragged array variable {ragged_array_dim_variable.name}, "
                    f"but dimension {dim_name} referenced from {ragged_attr_name} does not exist in file",
                )

            coord_var_reference_failures = []
            for coord_var in reference_attr_variables(ds, dim_name, " ", "dimensions"):
                if isinstance(coord_var, VariableReferenceError):
                    coord_var_reference_failures.append(coord_var)
                    domain_valid.messages.append(
                        f"Referenced coordinate variable {coord_var} does not exist in file",
                    )
                    continue
                # TODO: check for label variables
                if not set(
                    util.get_possible_label_variable_dimensions(coord_var),
                ).issubset({referenced_dim}):
                    domain_valid.messages.append(
                        f"Found ragged array variable {ragged_array_dim_variable.name}, "
                        f"but dimension {dim_name} referenced from {ragged_attr_name} does not exist in file",
                    )
                else:
                    domain_valid.score += 1
        else:
            # IMPLICIT CONFORMANCE REQUIRED 1/4
            # Domain variables must have a dimensions attribute
            for domain_var in (
                ds.variables[var]
                for var in ds.variables
                if "dimensions" in ds.variables[var].ncattrs()
            ):
                # store NetCDF variable attribute named "dimensions"
                dim_nc_attr = domain_var.getncattr("dimensions")
                # Every variable should have dimensionality (not to be confused with
                # dimensions NetCDF attribute here) even if scalar/empty, but non-empty
                # dimensions aren't allowed for domain variables

                domain_valid.out_of += 3
                if domain_var.dimensions:
                    domain_valid.messages.append(
                        "Domain variable "
                        f"{domain_var.name} should not have non-scalar/"
                        "non-empty variable dimensions",
                    )
                else:
                    domain_valid.score += 1

                # IMPLICIT CONFORMANCE REQUIRED 2/4
                # Aforementioned dimensions attribute is comprised of space separated
                # dimension names which must exist in the file

                domain_dims, dim_errors = reference_attr_variables(
                    ds,
                    dim_nc_attr,
                    " ",
                    "dimensions",
                )
                if dim_errors:
                    errors_str = ", ".join(dim_errors)
                    domain_valid.messages.append(
                        "Could not find the following "
                        "dimensions referenced in "
                        "dimensions attribute from "
                        "domain variable "
                        f"{domain_var.name}: {errors_str}",
                    )
                else:
                    domain_valid.score += 1
                if "coordinates" in domain_var.ncattrs():
                    domain_coord_vars, domain_coord_var_errors = (
                        reference_attr_variables(ds, domain_var.coordinates, " ")
                    )
                    if domain_coord_var_errors:
                        errors_str = ", ".join(
                            err.name for err in domain_coord_var_errors
                        )
                        domain_valid.messages.append(
                            "Could not find the following "
                            "variables referenced in "
                            "coordinates attribute from "
                            "domain variable "
                            f"{domain_var.name}: {errors_str}",
                        )
                    else:
                        domain_valid.score += 1
                        for coord_var in domain_coord_vars:
                            if coord_var not in get_coordinate_variables(ds):
                                continue
                            domain_valid.out_of += 1
                            domain_dims_names = {var.name for var in domain_dims}
                            variable_dim = util.get_possible_label_variable_dimensions(
                                coord_var,
                            )
                            if not (
                                set(
                                    util.get_possible_label_variable_dimensions(
                                        coord_var,
                                    ),
                                ).issubset(domain_dims_names)
                            ):
                                domain_valid.messages.append(
                                    "Could not find the following "
                                    "variables referenced in "
                                    "coordinates attribute from "
                                    "domain variable "
                                    f"{variable_dim}: {domain_dims_names}",
                                )
                            else:
                                domain_valid.score += 1

                # not in conformance docs, but mentioned as recommended anyways
                domain_valid.out_of += 1
                long_name = getattr(domain_var, "long_name", None)
                if long_name is None or not isinstance(long_name, str):
                    domain_valid.messages.append(
                        f"For domain variable {domain_var.name} "
                        f"it is recommended that attribute long_name be present and a string",
                    )
                else:
                    domain_valid.score += 1

                appendix_a_not_recommended_attrs = []
                for attr_name in domain_var.ncattrs():
                    if (
                        attr_name in self.appendix_a
                        and "D" not in self.appendix_a[attr_name]["attr_loc"]
                    ):
                        appendix_a_not_recommended_attrs.append(attr_name)

                domain_valid.out_of += 1
                if appendix_a_not_recommended_attrs:
                    domain_valid.messages.append(
                        f"The following attributes appear in variable {domain_var.name} "
                        "and CF Appendix A, but are not for use in domain variables: "
                        f"{appendix_a_not_recommended_attrs}",
                    )
                else:
                    # no errors occurred
                    domain_valid.score += 1

                # IMPLEMENTATION CONFORMANCE 5.8 REQUIRED 4/4
                # variables named by domain variable's cell_measures attributes must themselves be a subset
                # of dimensions named by domain variable's dimensions NetCDF attribute
                if hasattr(domain_var, "cell_measures"):
                    cell_measures_var_names = regex.findall(
                        r"\b(?:area|volume):\s+(\w+)",
                        domain_var.cell_measures,
                    )
                    # check exist
                    for var_name in cell_measures_var_names:
                        try:
                            cell_measures_variable = ds.variables[var_name]
                        except ValueError:
                            # TODO: what to do here?
                            continue
                        domain_valid.assert_true(
                            set(cell_measures_variable.dimensions).issubset(
                                {dd.name for dd in domain_dims},
                            ),
                            "Variables named in the cell_measures attributes must have dimensions with "
                            "values that are a subset of the referring domain variable's dimensions attribute",
                        )

            results.append(domain_valid.to_result())

        return results
