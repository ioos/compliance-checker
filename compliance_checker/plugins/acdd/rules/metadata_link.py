from xrlint.node import DatasetNode
from xrlint.rule import RuleContext, RuleOp

from compliance_checker.plugins.acdd.plugin import plugin


@plugin.define_rule(
    "1.3_metadata_link",
    version="1.3",
    docs_url="https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3",
)
class MetadataLink(RuleOp):
    def validate_dataset(self, ctx: RuleContext, node: DatasetNode):
        if "metadata_link" in node.dataset.attrs:
            value = node.dataset.attrs["metadata_link"]
            if "http" not in value:
                ctx.report(
                    f"Metadata URL should include http:// or https://: {value!r}",
                )
