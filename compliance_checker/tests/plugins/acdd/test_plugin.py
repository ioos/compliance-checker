from unittest import TestCase

from compliance_checker.plugins.acdd import export_plugin


class TestACDDPlugin(TestCase):
    def test_rules_complete(self):
        plugin = export_plugin()
        self.assertEqual(
            {
                "1.3_conventions",
                "1.0_attrs_suggested",
                "1.0_attrs_highly_recommended",
                "1.0_attrs_recommended",
                "1.3_attrs_recommended",
                "1.3_attrs_suggested",
                "1.3_attrs_highly_recommended",
            },
            set(plugin.rules.keys()),
        )

    def test_configs_complete(self):
        plugin = export_plugin()
        self.assertEqual(
            {
                "recommended",
                "acdd_1.3",
                "acdd_1.3_strict_reccomended",
                "acdd_1.3_strict",
                "acdd_1.3_warn",
                "acdd_1.1",
                "acdd_1.0",
            },
            set(plugin.configs.keys()),
        )
