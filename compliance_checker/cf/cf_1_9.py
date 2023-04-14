from netCDF4 import Dataset

from compliance_checker.base import BaseCheck, TestCtx
from compliance_checker.cf.cf_1_8 import CF1_8Check
from compliance_checker.cf.util import VariableReferenceError, reference_attr_variables


class CF1_9Check(CF1_8Check):
    _cc_spec_version = "1.9"
    _cc_url = "http://cfconventions.org/Data/cf-conventions/cf-conventions-1.9/cf-conventions.html"

    def __init__(self, options=None):
        super(CF1_9Check, self).__init__(options)
        self.section_titles.update({"5.8": "§5.8 Domain Variables"})

    def check_domain_variables(self, ds: Dataset):
        # Domain variables should have coordinates attribute, but should have
        # scalar dimensions
        results = []
        for domain_var in (
            var
            for var in ds.get_variables_by_attributes(
                coordinates=lambda c: c is not None
            )
            if not var.dimensions
        ):
            domain_valid = TestCtx(BaseCheck.MEDIUM, self.section_titles["5.8"])
            domain_valid.out_of += 1
            domain_coord_vars = reference_attr_variables(
                ds, domain_var.coordinates, " "
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
                    f"{domain_var.name}: {errors_str}"
                )

            else:
                long_name = getattr(domain_var, "long_name", None)
                if long_name is None or not isinstance(long_name, str):
                    domain_valid.messages.append(
                        f"For domain variable {domain_var.name} "
                        f"it is recommended that attribute long_name be present and a string"
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
                        f"{appendix_a_not_recommended_attrs}"
                    )

                # no errors occurred
                domain_valid.score += 1

            results.append(domain_valid.to_result())

        return results
