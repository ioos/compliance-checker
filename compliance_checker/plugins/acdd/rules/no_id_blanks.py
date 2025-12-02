from xrlint.node import DatasetNode
from xrlint.rule import RuleContext, RuleOp

from compliance_checker.plugins.acdd.plugin import plugin


@plugin.define_rule(
    "1.3_no_blanks_in_id",
    version="1.3",
    docs_url="https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3",
)
class NoBlanksInID(RuleOp):
    def validate_dataset(self, ctx: RuleContext, node: DatasetNode):
        try:
            value = node.dataset.attrs["id"]
        except KeyError:
            ctx.report(
                "Missing attribute 'id'",
                suggestions=["Include a non-blank 'id' attribute in the dataset."],
            )
            return
        if " " in value:
            ctx.report(
                "There should be no blanks in the id field",
                suggestions=["There should be no blanks in the id field"],
            )
