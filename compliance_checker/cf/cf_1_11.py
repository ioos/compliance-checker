from functools import lru_cache

import numpy as np

from compliance_checker.base import BaseCheck, Result, TestCtx
from compliance_checker.cf.appendix_a import appendix_a
from compliance_checker.cf.cf_1_10 import CF1_10Check
from compliance_checker.cf.util import VariableReferenceError, reference_attr_variables


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

    def check_bounds_inherit_attributes(self, ds):
        """
        A boundary variable inherits the values of some attributes from its parent coordinate variable.
        If a coordinate variable has any of the attributes marked "BI" (for "inherit") in the "Use" column of <<attribute-appendix>>, they are assumed to apply to its bounds variable as well.
        It is recommended that BI attributes not be included on a boundary variable.
        If a BI attribute is included, it must also be present in the parent variable, and it must exactly match the parent attribute's data type and value.
        A boundary variable can only have inheritable attributes if they are also present on its parent coordinate variable.
        A bounds variable may have any of the attributes marked "BO" for ("own") in the "Use" column of <<attribute-appendix>>.
        These attributes take precedence over any corresponding attributes of the parent variable.
        In these cases, the parent variable's attribute does not apply to the bounds variable, regardless of whether the latter has its own attribute.
        """
        results = []

        appendix_a_bi_attrs = {
            attribute_name
            for attribute_name, data_dict in appendix_a.items()
            if "BI" in data_dict["Use"]
        }
        for parent_variable in ds.get_variables_by_attributes(
            bounds=lambda b: b is not None,
        ):
            parent_bi_attrs = set(parent_variable.ncattrs()) & appendix_a_bi_attrs
            bounds_variable = reference_attr_variables(ds, parent_variable.bounds)[0]
            # nonexistent bounds variable, skip
            if isinstance(bounds_variable, VariableReferenceError):
                continue

            bounds_bi_ctx = self.get_test_ctx(
                BaseCheck.MEDIUM,
                self.section_titles["7.1"],
            )
            bounds_ncattr_set = set(bounds_variable.ncattrs())
            # IMPLEMENTATION CONFORMANCE 7.3 REQUIRED 4
            # A boundary variable can only have inheritable attributes, i.e. any of those marked "BI" in the "Use" column of Appendix A, if they are also present on its parent coordinate variable.
            bounds_bi_attrs = bounds_ncattr_set & appendix_a_bi_attrs
            # failure case 1, BI attr is only in bounds variable
            bounds_bi_only_attrs = bounds_bi_attrs - parent_bi_attrs
            # If a boundary variable has an inheritable attribute then its data type and its value must be exactly the same as the parent variable’s attribute.
            bounds_bi_ctx.out_of += 1
            if bounds_bi_only_attrs:
                bounds_bi_ctx.messages.append(
                    f"Bounds variable {bounds_variable.name} has the following attributes which must appear on the parent variable {parent_variable.name}: "
                    f"{sorted(bounds_bi_only_attrs)}",
                )
            else:
                bounds_bi_ctx.score += 1
            # failure case 2, BI attrs are in both bounds and parent variable
            both_bi_attr = bounds_bi_attrs & parent_bi_attrs
            no_match_attrs, match_attrs = [], []
            for bi_attr in both_bi_attr:
                bounds_bi_ctx.out_of += 1
                parent_attr_val = getattr(parent_variable, bi_attr)
                bounds_attr_val = getattr(bounds_variable, bi_attr)
                # IMPLEMENTATION CONFORMANCE 7.3 REQUIRED 5
                # If a boundary variable has an inheritable attribute then its data type and its value must be exactly the same as the parent variable’s attribute.
                if (
                    type(parent_attr_val) is not type(bounds_attr_val)
                    or parent_attr_val != bounds_attr_val
                ):
                    no_match_attrs.append(bi_attr)
                else:
                    match_attrs.append(bi_attr)

                pass_bi_both = True
                if no_match_attrs:
                    pass_bi_both = False
                    bounds_bi_ctx.messages.append(
                        f"Bounds variable {bounds_variable.name} and parent variable {parent_variable.name} have the following non matching boundary related attributes: {sorted(no_match_attrs)}",
                    )

                if match_attrs:
                    pass_bi_both = False
                    bounds_bi_ctx.messages.append(
                        f"Bounds variable {bounds_variable.name} and parent variable {parent_variable.name} have the following matching attributes {sorted(match_attrs)}.  It is recommended that only the parent variable of the bounds variable contains these attributes",
                    )
                bounds_bi_ctx.score += pass_bi_both

            results.append(bounds_bi_ctx.to_result())
        return results

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
        return test_ctx.to_result()

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

    def _check_add_offset_scale_factor_type(self, variable, attr_name):
        """
        Reusable function for checking both add_offset and scale_factor.
        """

        msgs = []
        error_msg_template = (
            f"When attribute {attr_name} is of type {{}} variable "
            f"{variable.name} must be of type {{}}"
        )

        att = getattr(variable, attr_name, None)
        float_dtypes = (np.byte, np.ubyte, np.short, np.ushort)
        passing = True
        if not isinstance(att, (np.float32 | np.float64)):
            passing = False
            msgs.append(f"Attribute {attr_name} must be a float or double")
        elif isinstance(att, np.float32) and variable.dtype not in (
            np.byte,
            np.ubyte,
            np.short,
            np.ushort,
        ):
            passing = False
            msgs.append(
                error_msg_template.format(
                    "float",
                    "byte, unsigned byte, short, or unsigned short",
                ),
            )
        elif isinstance(att, np.float64) and variable.dtype not in float_dtypes + (
            np.int32,
            np.uint32,
        ):
            passing = False
            msgs.append(
                error_msg_template.format(
                    "double",
                    "byte, unsigned byte, short, unsigned short, int, or unsigned int",
                ),
            )
        return Result(
            BaseCheck.MEDIUM,
            (int(passing), 1),
            self.section_titles["8.1"],
            msgs,
        )
