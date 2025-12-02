from xrlint.node import DatasetNode
from xrlint.rule import RuleContext, RuleOp

from compliance_checker.plugins.acdd.plugin import plugin
from compliance_checker.util import datetime_is_iso


@plugin.define_rule(
    "1.3_dates_iso_format",
    version="1.3",
    docs_url="https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3",
)
class IsoDates(RuleOp):
    def validate_dataset(self, ctx: RuleContext, node: DatasetNode):
        for attr in (
            "date_created",
            "date_issued",
            "date_modified",
            "date_metadata_modified",
        ):
            if attr in node.dataset.attrs:
                value = node.dataset.attrs[attr]
                iso_check, _msg = datetime_is_iso(value)

                if not iso_check:
                    ctx.report(
                        f"Attribute '{attr}' is not in ISO format: {value!r}",
                        suggestions=[
                            "Change '{attr}' to be in ISO format (e.g. YYYY-MM-DDThh:mm:ssZ).",
                        ],
                    )
