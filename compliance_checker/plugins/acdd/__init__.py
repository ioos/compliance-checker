from xrlint.plugin import Plugin
from xrlint.util.importutil import import_submodules

rules_1_0 = {
    "1.0_attrs_highly_recommended": "error",
    "1.0_attrs_recommended": "warn",
    "1.0_attrs_suggested": "warn",
}

rules_1_1 = {
    "1.1_attrs_highly_recommended": "error",
    "1.1_attrs_recommended": "warn",
    "1.1_attrs_suggested": "warn",
}

rules_1_3 = {
    "1.3_conventions": "error",
    "1.3_attrs_highly_recommended": "error",
    "1.3_attrs_recommended": "warn",
    "1.3_attrs_suggested": "warn",
    "1.3_no_blanks_in_id": "warn",
    "1.3_metadata_link": "warn",
    "1.3_dates_iso_format": "warn",
}


def export_plugin() -> Plugin:
    from .plugin import plugin

    import_submodules("compliance_checker.plugins.acdd.rules")

    plugin.define_config("recommended", {"name": "recommended", "rules": rules_1_3})

    plugin.define_config("acdd_1.3", {"name": "ACDD 1.3", "rules": rules_1_3})

    plugin.define_config(
        "acdd_1.3_strict_reccomended",
        {
            "name": "ACDD 1.3 (strict recommended)",
            "rules": {
                **rules_1_3,
                "1.3_attrs_recommended": "error",
            },
        },
    )

    plugin.define_config(
        "acdd_1.3_strict",
        {"name": "ACDD 1.3 (strict)", "rules": dict.fromkeys(rules_1_3, "error")},
    )

    plugin.define_config(
        "acdd_1.3_warn",
        {"name": "ACDD 1.3 (as warnings)", "rules": dict.fromkeys(rules_1_3, "warn")},
    )

    plugin.define_config("acdd_1.1", {"name": "ACDD 1.1", "rules": rules_1_1})

    plugin.define_config("acdd_1.0", {"name": "ACDD 1.0", "rules": rules_1_0})

    return plugin
