from functools import lru_cache

from compliance_checker.base import BaseCheck, TestCtx
from compliance_checker.cf.cf_1_10 import CF1_10Check


@lru_cache
def _temperature_standard_names(standard_name_table):
    re_ns = {"re": "http://exslt.org/regular-expressions"}
    temp_var_name_set = set(
        standard_name_table._root.xpath(
            "entry[re:test(canonical_units, " r"'(?:K|degree_C)(?:-?\d+)?')]/@id",
            namespaces=re_ns,
        ),
    )
    # need to include aliases for variable names that match temperature as well
    aliases = standard_name_table._root.findall("alias")
    all_names = temp_var_name_set.union(
        {
            alias.attrib["id"]
            for alias in aliases
            if alias.find("entry_id").text in temp_var_name_set
        },
    )
    return all_names


class CF1_11Check(CF1_10Check):
    _cc_spec_version = "1.11"
    _cc_url = "http://cfconventions.org/Data/cf-conventions/cf-conventions-1.11/cf-conventions.html"

    def __init__(self, options=None):
        super().__init__(options)
        self.section_titles.update(
            {
                "3.1.2": "§3.1.2 Temperature units",
            },
        )

    # IMPLEMENTATION CONFORMANCE 3.1 RECOMMENDED
    def check_temperature_units_metadata(self, ds):
        """Checks that units_metadata exists for variables with standard name of temperature"""
        temperature_variables = ds.get_variables_by_attributes(
            standard_name=lambda s: s in _temperature_standard_names(self._std_names),
        )
        if not temperature_variables:
            return []
        temperature_units_metadata_ctx = TestCtx(
            BaseCheck.MEDIUM,
            self.section_titles["3.1.2"],
        )
        for temperature_variable in temperature_variables:
            temperature_units_metadata_ctx.out_of += 1
            valid_temperature_units_metadata = [
                "temperature: difference",
                "temperature: on_scale",
                "temperature: unknown",
            ]
            if (
                getattr(temperature_variable, "units_metadata", None)
                not in valid_temperature_units_metadata
            ):
                temperature_units_metadata_ctx.messages.append(
                    f"Variable {temperature_variable.name} has a temperature related standard_name "
                    "and it is recommended that the units_metadata attribute is present and has one of the values "
                    f"{valid_temperature_units_metadata}",
                )
            else:
                temperature_units_metadata_ctx.score += 1

        return [temperature_units_metadata_ctx.to_result()]

    def check_time_units_metadata(self, ds):
        """Checks that units_metadata exists for time coordinates with specific calendar attributes"""
        time_variables = ds.get_variables_by_attributes(
            standard_name="time",
            calendar=lambda c: c in {"standard", "proleptic_gregorian", "julian"},
        )
        if not time_variables:
            return []

        time_units_metadata_ctx = self.get_test_ctx(
            BaseCheck.MEDIUM,
            self.section_titles["4.4"],
        )

        for time_variable in time_variables:
            time_units_metadata_ctx.out_of += 1
            valid_time_units_metadata = [
                "leap_seconds: none",
                "leap_seconds: utc",
                "leap_seconds: unknown",
            ]
            if (
                getattr(time_variable, "units_metadata", None)
                not in valid_time_units_metadata
            ):
                time_units_metadata_ctx.messages.append(
                    f"Variable {time_variable.name} has a calendar attribute of "
                    f"{time_variable.calendar} and it is recommended that the units_metadata attribute is present "
                    "and has one of the values "
                    f"{valid_time_units_metadata}",
                )
            else:
                time_units_metadata_ctx.score += 1

        return [time_units_metadata_ctx.to_result()]

    def check_single_cf_role(self, ds):
        test_ctx = self.get_test_ctx(
            BaseCheck.HIGH,
            self.section_titles["9.5"],
        )
        cf_role_var_names = [
            var.name
            for var in (ds.get_variables_by_attributes(cf_role=lambda x: x is not None))
        ]
        test_ctx.assert_true(
            len(cf_role_var_names) < 2,
            "There may only be one variable containing the cf_role attribute. "
            f"Currently the following variables have cf_role attributes: {cf_role_var_names}",
        )
        return test_ctx
