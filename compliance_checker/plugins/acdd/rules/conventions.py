from xrlint.node import DatasetNode
from xrlint.rule import RuleContext, RuleOp

from compliance_checker.plugins.acdd.plugin import plugin


@plugin.define_rule(
    "1.3_conventions",
    version="1.3",
    docs_url="https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3",
)
class Conventions(RuleOp):
    def validate_dataset(self, ctx: RuleContext, node: DatasetNode):
        if "Conventions" not in node.dataset.attrs:
            ctx.report(
                "Missing attribute 'Conventions'.",
                suggestions=[
                    "Include 'Conventions' with 'ACDD-1.3' attribute in the dataset.",
                ],
            )
        else:
            value = node.dataset.attrs.get("Conventions")
            if not isinstance(value, str) and value:
                ctx.report(f"Invalid attribute 'Conventions': {value!r}")
            elif "ACDD-1.3" not in value:
                ctx.report(
                    f"Attribute 'Conventions' needs to contain 'ACDD-1.3' in addition to the current: {value!r}",
                    suggestions=["Include 'ACDD-1.3' in the 'Conventions' attribute."],
                )
