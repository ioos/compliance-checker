from compliance_checker.base import BaseCheck
from compliance_checker.cf.cf_1_11 import CF1_11Check
from compliance_checker.util import iter_variables


class CF1_12Check(CF1_11Check):
    _cc_spec_version = "1.12"
    _cc_url = "https://cfconventions.org/Data/cf-conventions/cf-conventions-1.12/cf-conventions.html"

    def __init__(self, options=None):
        super().__init__(options)
        self.section_titles["1.3"] = "§1.3 Terminology"

    def check_string_variables_not_same_dimension_name(self, ds):
        """
        Checks that string variables don't have the same variable name as a dimension
        to avoid confusion with coordinate variables.
        """
        results = []
        for var in iter_variables(ds):
            if var.dtype is str:
                # IMPLEMENTATION CONFORMANCE 2.5 REQUIRED 1/1
                # However, in the CF specification itself, this requirement
                # resides in section 1.3
                var_str_test = self.get_test_ctx(BaseCheck.HIGH, self.section_titles["1.3"], var)
                var_str_test.out_of += 1
                # TODO: is this redundant with other coordinate variables we've already gathered?
                # Maybe use those instead?
                if var.ndim == 1 and var.dimensions[0] == var.name:
                    var_str_test.messages = [
                        f"One-dimensional variable length string variable {var.name} must not have the same name as its dimension's name to prevent confusion with coordinate variables"
                    ]
                else:
                    var_str_test.score += 1
                results.append(var_str_test.to_result())
        return results
