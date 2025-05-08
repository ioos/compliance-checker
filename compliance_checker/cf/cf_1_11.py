from functools import lru_cache

from compliance_checker.base import BaseCheck, TestCtx
from compliance_checker.cf.cf_1_10 import CF1_10Check


@lru_cache
def _temperature_standard_names(standard_name_table):
    re_ns = {"re": "http://exslt.org/regular-expressions"}
    return set(
        standard_name_table._root.xpath(
            "entry[re:test(canonical_units, " r"'K(-?\d+)?')]/@id",
            namespaces=re_ns,
        ),
    )


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
