import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_gold_standard_coverage.py"
SPEC = importlib.util.spec_from_file_location("audit_gold_standard_coverage", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)


class GoldStandardCoverageAuditTests(unittest.TestCase):
    def test_audit_reports_missing_cancer_theme_cells_and_negative_controls(self):
        gold = {
            "conclusions": [
                {
                    "gold_id": "pan_glycolysis",
                    "cancer_type": ["cancer"],
                    "mechanism_axis": "glycolysis lactate",
                    "expected_entities": ["glucose"],
                    "polarity": "positive",
                },
                {
                    "gold_id": "negative_ratio",
                    "mechanism_axis": "ratio",
                    "negative_trap_terms": ["ratio"],
                    "polarity": "negative",
                    "is_negative_control": True,
                },
            ]
        }
        requirements = {
            "dimensions": {
                "cancers": [
                    {"cancer_id": "pan_cancer", "terms": ["cancer"]},
                    {"cancer_id": "BRCA", "terms": ["breast cancer"]},
                ],
                "themes": [
                    {"theme_id": "glycolysis_lactate", "terms": ["glycolysis", "lactate"]},
                    {"theme_id": "glutamine_glutamate", "terms": ["glutamine", "glutamate"]},
                ],
            },
            "required_cells": [
                {"cancer_id": "pan_cancer", "theme_id": "glycolysis_lactate", "min_positive_controls": 1},
                {"cancer_id": "BRCA", "theme_id": "glutamine_glutamate", "min_positive_controls": 1},
            ],
            "minimum_negative_controls": 2,
        }

        report = audit.audit_coverage(gold, requirements)

        self.assertEqual(report["covered_cell_count"], 1)
        self.assertEqual(report["missing_cell_count"], 1)
        self.assertEqual(report["status"], "gaps_present")
        self.assertEqual(report["negative_control_requirement"]["status"], "missing")
        self.assertEqual(report["missing_cells"][0]["cancer_id"], "BRCA")

    def test_audit_marks_requirements_covered(self):
        gold = {
            "conclusions": [
                {
                    "gold_id": "brca_glutamine",
                    "cancer_type": ["breast cancer"],
                    "mechanism_axis": "glutamine glutamate metabolism",
                    "expected_entities": ["glutamine"],
                    "polarity": "positive",
                },
                {"gold_id": "neg1", "polarity": "negative", "is_negative_control": True, "negative_trap_terms": ["ratio"]},
            ]
        }
        requirements = {
            "dimensions": {
                "cancers": [{"cancer_id": "BRCA", "terms": ["breast cancer"]}],
                "themes": [{"theme_id": "glutamine_glutamate", "terms": ["glutamine", "glutamate"]}],
            },
            "required_cells": [{"cancer_id": "BRCA", "theme_id": "glutamine_glutamate", "min_positive_controls": 1}],
            "minimum_negative_controls": 1,
        }

        report = audit.audit_coverage(gold, requirements)

        self.assertEqual(report["status"], "covered")
        self.assertEqual(report["coverage_fraction"], 1.0)

    def test_pan_cancer_coverage_does_not_use_specific_cancer_rows(self):
        gold = {
            "conclusions": [
                {
                    "gold_id": "brca_nucleotide",
                    "cancer_type": ["breast cancer"],
                    "mechanism_axis": "nucleotide purine pyrimidine metabolism",
                    "expected_entities": ["purine"],
                    "polarity": "positive",
                },
                {"gold_id": "neg1", "polarity": "negative", "is_negative_control": True, "negative_trap_terms": ["ratio"]},
            ]
        }
        requirements = {
            "dimensions": {
                "cancers": [
                    {"cancer_id": "pan_cancer", "terms": ["cancer", "tumor", "carcinoma"]},
                    {"cancer_id": "BRCA", "terms": ["breast cancer"]},
                ],
                "themes": [{"theme_id": "nucleotide_metabolism", "terms": ["nucleotide", "purine", "pyrimidine"]}],
            },
            "required_cells": [
                {"cancer_id": "pan_cancer", "theme_id": "nucleotide_metabolism", "min_positive_controls": 1},
                {"cancer_id": "BRCA", "theme_id": "nucleotide_metabolism", "min_positive_controls": 1},
            ],
            "minimum_negative_controls": 1,
        }

        report = audit.audit_coverage(gold, requirements)

        by_cell = {(row["cancer_id"], row["theme_id"]): row for row in report["cells"]}
        self.assertEqual(by_cell[("pan_cancer", "nucleotide_metabolism")]["status"], "missing")
        self.assertEqual(by_cell[("BRCA", "nucleotide_metabolism")]["status"], "covered")

    def test_theme_coverage_uses_mechanism_fields_not_bridge_entities(self):
        gold = {
            "conclusions": [
                {
                    "gold_id": "pan_cancer_glutamine_only",
                    "cancer_type": ["cancer"],
                    "mechanism_axis": "glutamine glutamate metabolism",
                    "expected_entities": ["alpha-ketoglutarate"],
                    "polarity": "positive",
                },
                {"gold_id": "neg1", "polarity": "negative", "is_negative_control": True, "negative_trap_terms": ["ratio"]},
            ]
        }
        requirements = {
            "dimensions": {
                "cancers": [{"cancer_id": "pan_cancer", "terms": ["cancer", "tumor", "carcinoma"]}],
                "themes": [{"theme_id": "tca_cycle", "terms": ["tca", "tricarboxylic", "ketoglutarate"]}],
            },
            "required_cells": [{"cancer_id": "pan_cancer", "theme_id": "tca_cycle", "min_positive_controls": 1}],
            "minimum_negative_controls": 1,
        }

        report = audit.audit_coverage(gold, requirements)

        self.assertEqual(report["cells"][0]["status"], "missing")
        self.assertEqual(report["cells"][0]["matched_gold_ids"], [])

    def test_specificity_audit_flags_overbroad_gold_reuse_without_failing_coverage(self):
        gold = {
            "conclusions": [
                {
                    "gold_id": "coad_broad_core",
                    "cancer_type": ["colorectal cancer"],
                    "mechanism_axis": "glycolysis serine fatty acid glutathione metabolism",
                    "polarity": "positive",
                },
                {"gold_id": "neg1", "polarity": "negative", "is_negative_control": True, "negative_trap_terms": ["ratio"]},
            ]
        }
        requirements = {
            "recommended_max_cells_per_gold_id": 2,
            "dimensions": {
                "cancers": [{"cancer_id": "COAD_READ", "terms": ["colorectal cancer"]}],
                "themes": [
                    {"theme_id": "glycolysis_lactate", "terms": ["glycolysis"]},
                    {"theme_id": "serine_one_carbon", "terms": ["serine"]},
                    {"theme_id": "fatty_acid_lipid", "terms": ["fatty acid"]},
                ],
            },
            "required_cells": [
                {"cancer_id": "COAD_READ", "theme_id": "glycolysis_lactate", "min_positive_controls": 1},
                {"cancer_id": "COAD_READ", "theme_id": "serine_one_carbon", "min_positive_controls": 1},
                {"cancer_id": "COAD_READ", "theme_id": "fatty_acid_lipid", "min_positive_controls": 1},
            ],
            "minimum_negative_controls": 1,
        }

        report = audit.audit_coverage(gold, requirements)
        specificity = report["coverage_specificity"]

        self.assertEqual(report["status"], "covered")
        self.assertEqual(specificity["status"], "review_recommended")
        self.assertEqual(specificity["over_recommended_gold_ids"], ["coad_broad_core"])
        self.assertEqual(specificity["max_cells_per_gold_id_observed"], 3)
        self.assertEqual(specificity["multi_cell_gold_match_count"], 1)

    def test_specificity_audit_passes_when_gold_rows_are_cell_specific(self):
        gold = {
            "conclusions": [
                {
                    "gold_id": "brca_glycolysis",
                    "cancer_type": ["breast cancer"],
                    "mechanism_axis": "glycolysis lactate metabolism",
                    "polarity": "positive",
                },
                {
                    "gold_id": "brca_glutamine",
                    "cancer_type": ["breast cancer"],
                    "mechanism_axis": "glutamine glutamate metabolism",
                    "polarity": "positive",
                },
                {"gold_id": "neg1", "polarity": "negative", "is_negative_control": True, "negative_trap_terms": ["ratio"]},
            ]
        }
        requirements = {
            "recommended_max_cells_per_gold_id": 1,
            "dimensions": {
                "cancers": [{"cancer_id": "BRCA", "terms": ["breast cancer"]}],
                "themes": [
                    {"theme_id": "glycolysis_lactate", "terms": ["glycolysis"]},
                    {"theme_id": "glutamine_glutamate", "terms": ["glutamine"]},
                ],
            },
            "required_cells": [
                {"cancer_id": "BRCA", "theme_id": "glycolysis_lactate", "min_positive_controls": 1},
                {"cancer_id": "BRCA", "theme_id": "glutamine_glutamate", "min_positive_controls": 1},
            ],
            "minimum_negative_controls": 1,
        }

        report = audit.audit_coverage(gold, requirements)

        self.assertEqual(report["coverage_specificity"]["status"], "passed")
        self.assertEqual(report["coverage_specificity"]["single_cell_gold_fraction"], 1.0)
        self.assertEqual(report["coverage_specificity"]["over_recommended_count"], 0)

    def test_explicit_benchmark_coverage_cells_constrain_broad_gold_rows(self):
        gold = {
            "conclusions": [
                {
                    "gold_id": "coad_broad_core",
                    "cancer_type": ["colorectal cancer"],
                    "mechanism_axis": "glycolysis serine fatty acid metabolism",
                    "benchmark_coverage_cells": ["COAD_READ:glycolysis_lactate"],
                    "polarity": "positive",
                },
                {
                    "gold_id": "coad_serine_specific",
                    "cancer_type": ["colorectal cancer"],
                    "mechanism_axis": "serine one-carbon metabolism",
                    "benchmark_coverage_cells": [{"cancer_id": "COAD_READ", "theme_id": "serine_one_carbon"}],
                    "polarity": "positive",
                },
                {"gold_id": "neg1", "polarity": "negative", "is_negative_control": True, "negative_trap_terms": ["ratio"]},
            ]
        }
        requirements = {
            "recommended_max_cells_per_gold_id": 1,
            "dimensions": {
                "cancers": [{"cancer_id": "COAD_READ", "terms": ["colorectal cancer"]}],
                "themes": [
                    {"theme_id": "glycolysis_lactate", "terms": ["glycolysis"]},
                    {"theme_id": "serine_one_carbon", "terms": ["serine"]},
                    {"theme_id": "fatty_acid_lipid", "terms": ["fatty acid"]},
                ],
            },
            "required_cells": [
                {"cancer_id": "COAD_READ", "theme_id": "glycolysis_lactate", "min_positive_controls": 1},
                {"cancer_id": "COAD_READ", "theme_id": "serine_one_carbon", "min_positive_controls": 1},
                {"cancer_id": "COAD_READ", "theme_id": "fatty_acid_lipid", "min_positive_controls": 1},
            ],
            "minimum_negative_controls": 1,
        }

        report = audit.audit_coverage(gold, requirements)
        by_cell = {(row["cancer_id"], row["theme_id"]): row for row in report["cells"]}

        self.assertEqual(by_cell[("COAD_READ", "glycolysis_lactate")]["matched_gold_ids"], ["coad_broad_core"])
        self.assertEqual(by_cell[("COAD_READ", "serine_one_carbon")]["matched_gold_ids"], ["coad_serine_specific"])
        self.assertEqual(by_cell[("COAD_READ", "fatty_acid_lipid")]["matched_gold_ids"], [])
        self.assertEqual(report["coverage_specificity"]["over_recommended_count"], 0)


if __name__ == "__main__":
    unittest.main()
