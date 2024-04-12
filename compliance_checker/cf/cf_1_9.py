import cftime
import numpy as np
import regex
from netCDF4 import Dataset

from compliance_checker import cfutil
from compliance_checker.base import BaseCheck, Result, TestCtx
from compliance_checker.cf import util
from compliance_checker.cf.cf_1_8 import CF1_8Check
from compliance_checker.cf.util import VariableReferenceError, reference_attr_variables


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

                    ret_val.append(result)
        return ret_val

    def check_time_coordinate(self, ds):
        prev_return = super().check_time_coordinate(ds)
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
            test_ctx.assert_true(
                match.group("seconds") is None or int(match.group("seconds")) < 60,
                f'Time coordinate variable "{name}" must have '
                "units with seconds less than 60",
            )
            prev_return.append(test_ctx.to_result())
        return prev_return

    def check_domain_variables(self, ds: Dataset):
        # Domain variables should have coordinates attribute, but should have
        # scalar dimensions
        results = []
        for domain_var in (
            var
            for var in ds.get_variables_by_attributes(
                coordinates=lambda c: c is not None,
            )
            # IMPLICIT
            if not var.dimensions
        ):
            domain_valid = TestCtx(BaseCheck.MEDIUM, self.section_titles["5.8"])
            domain_valid.out_of += 1
            domain_coord_vars = reference_attr_variables(
                ds,
                domain_var.coordinates,
                " ",
            )
            errors = [
                maybe_error.name
                for maybe_error in domain_coord_vars
                if isinstance(maybe_error, VariableReferenceError)
            ]
            if errors:
                errors_str = ", ".join(errors)
                domain_valid.messages.append(
                    "Could not find the following "
                    "variables referenced in "
                    "coordinates attribute from "
                    "domain variable "
                    f"{domain_var.name}: {errors_str}",
                )

            else:
                long_name = getattr(domain_var, "long_name", None)
                if long_name is None or not isinstance(long_name, str):
                    domain_valid.messages.append(
                        f"For domain variable {domain_var.name} "
                        f"it is recommended that attribute long_name be present and a string",
                    )
                    results.append(domain_valid.to_result())
                    continue
                appendix_a_not_recommended_attrs = []
                for attr_name in domain_var.ncattrs():
                    if "D" not in self.appendix_a[attr_name]["attr_loc"]:
                        appendix_a_not_recommended_attrs.append(attr_name)

                if appendix_a_not_recommended_attrs:
                    domain_valid.messages.append(
                        f"The following attributes appear in variable {domain_var.name} "
                        "and CF Appendix A, but are not for use in domain variables: "
                        f"{appendix_a_not_recommended_attrs}",
                    )

                # no errors occurred
                domain_valid.score += 1

            results.append(domain_valid.to_result())

        return results
