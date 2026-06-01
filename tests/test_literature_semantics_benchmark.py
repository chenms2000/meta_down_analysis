import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_literature_semantics_benchmark.py"
SPEC = importlib.util.spec_from_file_location("run_literature_semantics_benchmark", SCRIPT)
benchmark = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = benchmark
SPEC.loader.exec_module(benchmark)


class LiteratureSemanticsBenchmarkTests(unittest.TestCase):
    def test_benchmark_passes_complex_semantics_cases(self):
        payload = {
            "minimum_status_accuracy": 1.0,
            "minimum_cue_recall": 1.0,
            "minimum_method_cue_recall": 1.0,
            "cases": [
                {
                    "case_id": "not_only_support",
                    "text": "Glucose was not only increased but also quantified by LC-MS in tumor tissue.",
                    "section": "results",
                    "expected_support_status": "support",
                    "expected_cues": ["non_negating_negation", "direct_assay"],
                    "expected_method_cues": ["patient_sample", "measurement_assay"],
                },
                {
                    "case_id": "correlation_not_causality",
                    "text": "Glutamine correlated with tumor burden, but correlation does not imply causality.",
                    "section": "discussion",
                    "expected_support_status": "uncertain",
                    "expected_cues": ["concession", "causal_weakening"],
                    "expected_method_cues": [],
                },
            ],
        }

        report = benchmark.benchmark_literature_semantics(payload)

        self.assertEqual(report["status"], "passed")
        self.assertEqual(report["metrics"]["status_accuracy"], 1.0)
        self.assertEqual(report["metrics"]["cue_recall"], 1.0)

    def test_benchmark_reports_failed_case_without_lowering_thresholds(self):
        payload = {
            "minimum_status_accuracy": 1.0,
            "minimum_cue_recall": 1.0,
            "minimum_method_cue_recall": 1.0,
            "cases": [
                {
                    "case_id": "wrong_expectation",
                    "text": "Lactate showed no significant difference between tumor and adjacent tissue.",
                    "section": "results",
                    "expected_support_status": "support",
                    "expected_cues": ["comparison_context"],
                    "expected_method_cues": [],
                }
            ],
        }

        report = benchmark.benchmark_literature_semantics(payload)

        self.assertEqual(report["status"], "failed")
        self.assertEqual(report["metrics"]["failed_case_count"], 1)
        self.assertEqual(report["failed_cases"][0]["case_id"], "wrong_expectation")
        self.assertEqual(report["failed_cases"][0]["observed_support_status"], "contradict")

    def test_benchmark_requires_minimum_case_count(self):
        payload = {
            "minimum_status_accuracy": 1.0,
            "minimum_cue_recall": 1.0,
            "minimum_method_cue_recall": 1.0,
            "minimum_case_count": 2,
            "cases": [
                {
                    "case_id": "single_case",
                    "text": "Patient tumor tissue lactate was quantified by LC-MS.",
                    "section": "results",
                    "expected_support_status": "support",
                    "expected_cues": ["direct_assay"],
                    "expected_method_cues": ["patient_sample", "measurement_assay"],
                }
            ],
        }

        report = benchmark.benchmark_literature_semantics(payload)

        self.assertEqual(report["status"], "failed")
        self.assertEqual(report["metrics"]["case_count"], 1)
        self.assertEqual(report["metrics"]["minimum_case_count"], 2)

    def test_benchmark_requires_declared_semantic_coverage(self):
        payload = {
            "minimum_status_accuracy": 1.0,
            "minimum_cue_recall": 1.0,
            "minimum_method_cue_recall": 1.0,
            "required_support_statuses": ["support", "background"],
            "required_semantic_cues": ["direct_assay", "context_boundary"],
            "required_method_cues": ["patient_sample", "animal_model"],
            "cases": [
                {
                    "case_id": "patient_direct_support",
                    "text": "Patient tumor tissue lactate was quantified by LC-MS.",
                    "section": "results",
                    "expected_support_status": "support",
                    "expected_cues": ["direct_assay"],
                    "expected_method_cues": ["patient_sample", "measurement_assay"],
                }
            ],
        }

        report = benchmark.benchmark_literature_semantics(payload)

        self.assertEqual(report["status"], "failed")
        self.assertEqual(report["metrics"]["required_coverage_missing_count"], 3)
        self.assertEqual(report["coverage"]["support_statuses"]["missing"], ["background"])
        self.assertEqual(report["coverage"]["semantic_cues"]["missing"], ["context_boundary"])
        self.assertEqual(report["coverage"]["method_cues"]["missing"], ["animal_model"])

    def test_benchmark_passes_when_declared_coverage_is_present(self):
        payload = {
            "minimum_status_accuracy": 1.0,
            "minimum_cue_recall": 1.0,
            "minimum_method_cue_recall": 1.0,
            "required_support_statuses": ["support", "uncertain"],
            "required_semantic_cues": ["direct_assay", "context_boundary"],
            "required_method_cues": ["patient_sample", "animal_model"],
            "cases": [
                {
                    "case_id": "patient_direct_support",
                    "text": "Patient tumor tissue lactate was quantified by LC-MS.",
                    "section": "results",
                    "expected_support_status": "support",
                    "expected_cues": ["direct_assay"],
                    "expected_method_cues": ["patient_sample", "measurement_assay"],
                },
                {
                    "case_id": "in_vitro_not_in_vivo_boundary",
                    "text": "Glutamine deprivation reduced proliferation in vitro but not in vivo.",
                    "section": "results",
                    "expected_support_status": "uncertain",
                    "expected_cues": ["context_boundary"],
                    "expected_method_cues": ["cell_line_model", "animal_model"],
                },
            ],
        }

        report = benchmark.benchmark_literature_semantics(payload)

        self.assertEqual(report["status"], "passed")
        self.assertEqual(report["metrics"]["required_coverage_missing_count"], 0)
        self.assertEqual(report["coverage"]["semantic_cues"]["status"], "covered")

    def test_benchmark_expands_generated_case_groups(self):
        payload = {
            "minimum_status_accuracy": 1.0,
            "minimum_cue_recall": 1.0,
            "minimum_method_cue_recall": 1.0,
            "minimum_case_count": 2,
            "minimum_generated_case_count": 2,
            "required_generated_groups": ["systematic_negation"],
            "generated_case_groups": [
                {
                    "group_id": "systematic_negation",
                    "section": "results",
                    "expected_support_status": "contradict",
                    "expected_cues": ["null_result", "negation"],
                    "expected_method_cues": ["patient_sample"],
                    "texts": [
                        "Lactate showed no significant difference between tumor and adjacent tissue.",
                        {
                            "text": "Glucose did not increase relative to matched normal samples.",
                            "expected_cues": ["null_result", "negation", "comparison_context"],
                        },
                    ],
                }
            ],
        }

        report = benchmark.benchmark_literature_semantics(payload)

        self.assertEqual(report["status"], "passed")
        self.assertEqual(report["metrics"]["manual_case_count"], 0)
        self.assertEqual(report["metrics"]["generated_case_count"], 2)
        self.assertEqual(report["metrics"]["generated_group_count"], 1)
        self.assertEqual(report["coverage"]["generated_groups"]["status"], "covered")
        self.assertTrue(all(row["generated_from_group"] == "systematic_negation" for row in report["cases"]))

    def test_benchmark_blocks_missing_generated_group_coverage(self):
        payload = {
            "minimum_status_accuracy": 1.0,
            "minimum_cue_recall": 1.0,
            "minimum_method_cue_recall": 1.0,
            "minimum_generated_case_count": 2,
            "required_generated_groups": ["systematic_negation", "context_boundary"],
            "generated_case_groups": [
                {
                    "group_id": "systematic_negation",
                    "section": "results",
                    "expected_support_status": "contradict",
                    "expected_cues": ["null_result", "negation"],
                    "expected_method_cues": ["patient_sample"],
                    "texts": ["Serine was not significantly elevated in clinical cohorts."],
                }
            ],
        }

        report = benchmark.benchmark_literature_semantics(payload)

        self.assertEqual(report["status"], "failed")
        self.assertEqual(report["metrics"]["generated_case_count"], 1)
        self.assertEqual(report["metrics"]["minimum_generated_case_count"], 2)
        self.assertEqual(report["coverage"]["generated_groups"]["missing"], ["context_boundary"])


if __name__ == "__main__":
    unittest.main()
