import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_gold_standard_quality.py"
SPEC = importlib.util.spec_from_file_location("audit_gold_standard_quality", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)


class GoldStandardQualityAuditTests(unittest.TestCase):
    def test_quality_audit_passes_required_fields_and_trap_categories(self):
        gold = {
            "conclusions": [
                {
                    "gold_id": "pos_a",
                    "cancer_type": ["cancer"],
                    "mechanism_axis": "glycolysis lactate",
                    "expected_entities": ["glucose", "lactate"],
                    "expected_direction": "context_dependent",
                    "required_evidence_type": "reaction",
                    "source_refs": ["PMID:1"],
                    "polarity": "positive",
                },
                {
                    "gold_id": "neg_a",
                    "mechanism_axis": "ratio background",
                    "negative_trap_terms": ["ratio", "background"],
                    "source_refs": ["internal:test"],
                    "polarity": "negative",
                    "is_negative_control": True,
                },
            ]
        }
        requirements = {
            "minimum_positive_conclusions": 1,
            "minimum_negative_controls": 1,
            "required_positive_fields": [
                "gold_id",
                "cancer_type",
                "mechanism_axis",
                "expected_entities",
                "expected_direction",
                "required_evidence_type",
                "source_refs",
            ],
            "required_negative_fields": ["gold_id", "mechanism_axis", "negative_trap_terms", "source_refs"],
            "allowed_source_ref_prefixes": ["PMID:", "internal:"],
            "minimum_source_refs_per_positive": 1,
            "minimum_expected_entities_per_positive": 2,
            "negative_trap_categories": [{"category_id": "ratio", "terms": ["ratio"]}],
        }

        report = audit.audit_quality(gold, requirements)

        self.assertEqual(report["status"], "passed")
        self.assertEqual(report["positive_count"], 1)
        self.assertEqual(report["negative_control_count"], 1)

    def test_quality_audit_reports_count_field_source_and_category_failures(self):
        gold = {
            "conclusions": [
                {
                    "gold_id": "pos_bad",
                    "cancer_type": ["cancer"],
                    "mechanism_axis": "glycolysis",
                    "expected_entities": ["glucose"],
                    "expected_direction": "context_dependent",
                    "required_evidence_type": "reaction",
                    "source_refs": ["bad:ref"],
                    "polarity": "positive",
                }
            ]
        }
        requirements = {
            "minimum_positive_conclusions": 2,
            "minimum_negative_controls": 1,
            "required_positive_fields": ["gold_id", "source_refs", "expected_entities"],
            "required_negative_fields": ["gold_id", "source_refs"],
            "allowed_source_ref_prefixes": ["PMID:"],
            "minimum_source_refs_per_positive": 1,
            "minimum_expected_entities_per_positive": 2,
            "negative_trap_categories": [{"category_id": "ratio", "terms": ["ratio"]}],
        }

        report = audit.audit_quality(gold, requirements)
        codes = {row["code"] for row in report["failed_checks"]}

        self.assertEqual(report["status"], "failed")
        self.assertIn("positive_gold_count_low", codes)
        self.assertIn("negative_control_count_low", codes)
        self.assertIn("positive_gold_field_failures", codes)
        self.assertIn("negative_trap_category_gaps", codes)
        self.assertEqual(report["positive_field_failures"][0]["invalid_source_refs"], ["bad:ref"])

    def test_quality_audit_reports_positive_distribution_failures(self):
        gold = {
            "conclusions": [
                {
                    "gold_id": "pan_cancer_glycolysis",
                    "cancer_type": ["cancer"],
                    "mechanism_axis": "glycolysis lactate",
                    "expected_entities": ["glucose", "lactate"],
                    "expected_direction": "context_dependent",
                    "required_evidence_type": "reaction",
                    "source_refs": ["PMID:1"],
                    "polarity": "positive",
                },
                {
                    "gold_id": "pan_cancer_glutamine",
                    "cancer_type": ["cancer"],
                    "mechanism_axis": "glutamine metabolism",
                    "expected_entities": ["glutamine", "glutamate"],
                    "expected_direction": "context_dependent",
                    "required_evidence_type": "reaction",
                    "source_refs": ["PMID:1"],
                    "polarity": "positive",
                },
                {
                    "gold_id": "negative_ratio",
                    "mechanism_axis": "ratio background",
                    "negative_trap_terms": ["ratio", "background"],
                    "source_refs": ["internal:test"],
                    "polarity": "negative",
                    "is_negative_control": True,
                },
            ]
        }
        requirements = {
            "minimum_positive_conclusions": 2,
            "minimum_negative_controls": 1,
            "required_positive_fields": ["gold_id", "cancer_type", "mechanism_axis", "expected_entities", "source_refs"],
            "required_negative_fields": ["gold_id", "source_refs", "negative_trap_terms"],
            "allowed_source_ref_prefixes": ["PMID:", "internal:"],
            "minimum_source_refs_per_positive": 1,
            "minimum_expected_entities_per_positive": 2,
            "minimum_context_specific_positive_conclusions": 1,
            "maximum_pan_cancer_positive_fraction": 0.5,
            "minimum_unique_positive_cancer_contexts": 2,
            "minimum_unique_positive_mechanism_axes": 3,
            "minimum_unique_positive_source_refs": 2,
            "required_positive_theme_terms": [
                {"theme_id": "glycolysis", "terms": ["glycolysis"]},
                {"theme_id": "redox", "terms": ["redox", "glutathione"]},
            ],
            "negative_trap_categories": [{"category_id": "ratio", "terms": ["ratio"]}],
        }

        report = audit.audit_quality(gold, requirements)
        codes = {row["code"] for row in report["failed_checks"]}

        self.assertEqual(report["status"], "failed")
        self.assertIn("context_specific_positive_count_low", codes)
        self.assertIn("pan_cancer_positive_fraction_high", codes)
        self.assertIn("positive_cancer_context_diversity_low", codes)
        self.assertIn("positive_mechanism_axis_diversity_low", codes)
        self.assertIn("positive_source_ref_diversity_low", codes)
        self.assertIn("positive_required_theme_gaps", codes)
        self.assertEqual(report["positive_distribution"]["missing_required_themes"][0]["theme_id"], "redox")

    def test_quality_audit_reports_required_theme_depth_failures(self):
        gold = {
            "conclusions": [
                {
                    "gold_id": "pos_glycolysis",
                    "cancer_type": ["cancer"],
                    "mechanism_axis": "glycolysis lactate",
                    "expected_entities": ["glucose", "lactate"],
                    "expected_direction": "context_dependent",
                    "required_evidence_type": "reaction",
                    "source_refs": ["PMID:1"],
                    "polarity": "positive",
                },
                {
                    "gold_id": "pos_redox_a",
                    "cancer_type": ["breast cancer"],
                    "mechanism_axis": "redox glutathione",
                    "expected_entities": ["glutathione", "cysteine"],
                    "expected_direction": "context_dependent",
                    "required_evidence_type": "reaction",
                    "source_refs": ["PMID:2"],
                    "polarity": "positive",
                },
                {
                    "gold_id": "pos_redox_b",
                    "cancer_type": ["colon cancer"],
                    "mechanism_axis": "redox glutathione oxidative stress",
                    "expected_entities": ["glutathione", "nadph"],
                    "expected_direction": "context_dependent",
                    "required_evidence_type": "reaction",
                    "source_refs": ["PMID:3"],
                    "polarity": "positive",
                },
                {
                    "gold_id": "negative_ratio",
                    "mechanism_axis": "ratio background",
                    "negative_trap_terms": ["ratio", "background"],
                    "source_refs": ["internal:test"],
                    "polarity": "negative",
                    "is_negative_control": True,
                },
            ]
        }
        requirements = {
            "minimum_positive_conclusions": 3,
            "minimum_negative_controls": 1,
            "required_positive_fields": ["gold_id", "cancer_type", "mechanism_axis", "expected_entities", "source_refs"],
            "required_negative_fields": ["gold_id", "source_refs", "negative_trap_terms"],
            "allowed_source_ref_prefixes": ["PMID:", "internal:"],
            "minimum_source_refs_per_positive": 1,
            "minimum_expected_entities_per_positive": 2,
            "minimum_positive_controls_per_required_theme": 2,
            "required_positive_theme_terms": [
                {"theme_id": "glycolysis", "terms": ["glycolysis"]},
                {"theme_id": "redox", "terms": ["redox", "glutathione"]},
            ],
            "negative_trap_categories": [{"category_id": "ratio", "terms": ["ratio"]}],
        }

        report = audit.audit_quality(gold, requirements)
        by_code = {row["code"]: row for row in report["failed_checks"]}

        self.assertEqual(report["status"], "failed")
        self.assertIn("positive_required_theme_depth_low", by_code)
        self.assertEqual(by_code["positive_required_theme_depth_low"]["theme_ids"], ["glycolysis"])

    def test_quality_audit_passes_positive_distribution_requirements(self):
        gold = {
            "conclusions": [
                {
                    "gold_id": "pan_cancer_glycolysis",
                    "cancer_type": ["cancer"],
                    "mechanism_axis": "glycolysis lactate",
                    "expected_entities": ["glucose", "lactate"],
                    "expected_direction": "context_dependent",
                    "required_evidence_type": "reaction",
                    "source_refs": ["PMID:1"],
                    "polarity": "positive",
                },
                {
                    "gold_id": "brca_redox",
                    "cancer_type": ["breast cancer"],
                    "mechanism_axis": "redox glutathione",
                    "expected_entities": ["glutathione", "cysteine"],
                    "expected_direction": "context_dependent",
                    "required_evidence_type": "reaction",
                    "source_refs": ["PMID:2"],
                    "polarity": "positive",
                },
                {
                    "gold_id": "negative_ratio",
                    "mechanism_axis": "ratio background",
                    "negative_trap_terms": ["ratio", "background"],
                    "source_refs": ["internal:test"],
                    "polarity": "negative",
                    "is_negative_control": True,
                },
            ]
        }
        requirements = {
            "minimum_positive_conclusions": 2,
            "minimum_negative_controls": 1,
            "required_positive_fields": ["gold_id", "cancer_type", "mechanism_axis", "expected_entities", "source_refs"],
            "required_negative_fields": ["gold_id", "source_refs", "negative_trap_terms"],
            "allowed_source_ref_prefixes": ["PMID:", "internal:"],
            "minimum_source_refs_per_positive": 1,
            "minimum_expected_entities_per_positive": 2,
            "minimum_context_specific_positive_conclusions": 1,
            "maximum_pan_cancer_positive_fraction": 0.5,
            "minimum_unique_positive_cancer_contexts": 2,
            "minimum_unique_positive_mechanism_axes": 2,
            "minimum_unique_positive_source_refs": 2,
            "minimum_positive_controls_per_required_theme": 1,
            "required_positive_theme_terms": [
                {"theme_id": "glycolysis", "terms": ["glycolysis"]},
                {"theme_id": "redox", "terms": ["redox", "glutathione"]},
            ],
            "negative_trap_categories": [{"category_id": "ratio", "terms": ["ratio"]}],
        }

        report = audit.audit_quality(gold, requirements)

        self.assertEqual(report["status"], "passed")
        self.assertEqual(report["positive_distribution"]["context_specific_positive_count"], 1)
        self.assertEqual(report["positive_distribution"]["missing_required_themes"], [])


if __name__ == "__main__":
    unittest.main()
