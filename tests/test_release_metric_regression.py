import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_release_metric_regression.py"
SPEC = importlib.util.spec_from_file_location("audit_release_metric_regression", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)


def batch_report():
    return {
        "status": {
            "passed": False,
            "failed_case_count": 0,
            "missing_case_count": 1,
            "benchmark_matrix_passed": False,
            "traceability_contract_passed": False,
        },
        "benchmark_matrix": {"registered_gap_count": 0, "executed_gap_count": 2},
        "scenario_traceability_contract": {"status": "failed", "failed_case_count": 1},
        "gold_quality": {
            "status": "passed",
            "positive_count": 20,
            "negative_control_count": 10,
            "positive_distribution": {
                "context_specific_positive_count": 12,
                "unique_positive_cancer_context_count": 7,
                "unique_positive_mechanism_axis_count": 18,
            },
        },
        "gold_coverage": {
            "required_cell_count": 28,
            "covered_cell_count": 28,
            "missing_cell_count": 0,
            "coverage_fraction": 1.0,
            "coverage_specificity": {
                "over_recommended_count": 1,
                "multi_cell_gold_match_count": 3,
                "single_cell_gold_fraction": 0.8,
                "max_cells_per_gold_id_observed": 3,
            },
            "status": "covered",
        },
        "literature_semantics": {
            "status": "passed",
            "metrics": {
                "case_count": 51,
                "manual_case_count": 36,
                "generated_case_count": 15,
                "generated_group_count": 7,
                "generated_failed_case_count": 0,
                "status_accuracy": 1.0,
                "cue_recall": 1.0,
                "method_cue_recall": 1.0,
                "required_coverage_missing_count": 0,
            },
        },
        "heldout_validation": {
            "status": "passed",
            "rerun_performed": True,
            "failed_checks": [],
            "warnings": [{"code": "source_heldout_random_auc"}, {"code": "other_warning"}],
        },
        "scenarios": [
            {
                "conclusion_evaluation": {
                    "metrics": {
                        "negative_trap_hits": 0,
                        "unsupported_top_claim_rate": 0.1,
                        "precision_at_1": 1.0,
                        "precision_at_3": 0.667,
                        "precision_at_5": 0.6,
                        "recall_at_1": 0.25,
                        "recall_at_3": 0.5,
                        "recall_at_5": 0.75,
                        "expected_topk_eligible": True,
                        "expected_precision_at_1": 1.0,
                        "expected_precision_at_3": 0.667,
                        "expected_precision_at_5": 0.6,
                        "expected_recall_at_1": 0.5,
                        "expected_recall_at_3": 1.0,
                        "expected_recall_at_5": 1.0,
                    }
                }
            },
            {
                "conclusion_evaluation": {
                    "metrics": {
                        "negative_trap_hits": 0,
                        "unsupported_top_claim_rate": 0.2,
                        "precision_at_1": 0.0,
                        "precision_at_3": 0.333,
                        "precision_at_5": 0.4,
                        "recall_at_1": 0.0,
                        "recall_at_3": 0.25,
                        "recall_at_5": 0.5,
                        "expected_topk_eligible": False,
                        "expected_precision_at_1": None,
                        "expected_precision_at_3": None,
                        "expected_precision_at_5": None,
                        "expected_recall_at_1": None,
                        "expected_recall_at_3": None,
                        "expected_recall_at_5": None,
                    }
                }
            },
        ],
    }


class ReleaseMetricRegressionTests(unittest.TestCase):
    def test_snapshot_extracts_stable_release_metrics(self):
        snap = audit.snapshot(
            batch_report(),
            {"status": "blocked", "blocker_count": 4},
            {"status": "incomplete", "satisfied_count": 3, "incomplete_count": 3},
        )
        metrics = snap["metrics"]

        self.assertEqual(metrics["gold_positive_count"], 20)
        self.assertEqual(metrics["gold_required_cell_count"], 28)
        self.assertEqual(metrics["gold_covered_cell_count"], 28)
        self.assertEqual(metrics["gold_missing_cell_count"], 0)
        self.assertEqual(metrics["gold_coverage_fraction"], 1.0)
        self.assertEqual(metrics["gold_specificity_over_recommended_count"], 1)
        self.assertEqual(metrics["gold_specificity_multi_cell_match_count"], 3)
        self.assertEqual(metrics["gold_specificity_single_cell_fraction"], 0.8)
        self.assertEqual(metrics["gold_specificity_max_cells_per_gold_id"], 3)
        self.assertEqual(metrics["literature_case_count"], 51)
        self.assertEqual(metrics["literature_generated_case_count"], 15)
        self.assertEqual(metrics["literature_generated_group_count"], 7)
        self.assertEqual(metrics["literature_generated_failed_case_count"], 0)
        self.assertEqual(metrics["literature_status_accuracy"], 1.0)
        self.assertEqual(metrics["benchmark_matrix_executed_gap_count"], 2)
        self.assertEqual(metrics["scenario_traceability_failed_case_count"], 1)
        self.assertEqual(metrics["scenario_traceability_passed"], 0)
        self.assertEqual(metrics["scenario_unsupported_top_claim_rate_max"], 0.2)
        self.assertEqual(metrics["scenario_precision_at_1_min"], 0.0)
        self.assertEqual(metrics["scenario_precision_at_3_min"], 0.333)
        self.assertEqual(metrics["scenario_precision_at_5_min"], 0.4)
        self.assertEqual(metrics["scenario_recall_at_1_min"], 0.0)
        self.assertEqual(metrics["scenario_recall_at_3_min"], 0.25)
        self.assertEqual(metrics["scenario_recall_at_5_min"], 0.5)
        self.assertEqual(metrics["scenario_expected_topk_eligible_count"], 1.0)
        self.assertEqual(metrics["scenario_expected_precision_at_1_min"], 1.0)
        self.assertEqual(metrics["scenario_expected_precision_at_3_min"], 0.667)
        self.assertEqual(metrics["scenario_expected_precision_at_5_min"], 0.6)
        self.assertEqual(metrics["scenario_expected_recall_at_1_min"], 0.5)
        self.assertEqual(metrics["scenario_expected_recall_at_3_min"], 1.0)
        self.assertEqual(metrics["scenario_expected_recall_at_5_min"], 1.0)
        self.assertEqual(metrics["heldout_warning_count"], 2)
        self.assertEqual(metrics["heldout_random_auc_warning_count"], 1)

    def test_compare_reports_snapshot_only_without_baseline(self):
        snap = {"metrics": {"literature_status_accuracy": 1.0}}

        comparison = audit.compare_snapshots(snap, {})

        self.assertEqual(comparison["status"], "snapshot_only")
        self.assertEqual(comparison["regression_count"], 0)

    def test_compare_detects_higher_is_better_and_lower_is_better_regressions(self):
        baseline = {
            "metrics": {
                "literature_status_accuracy": 1.0,
                "literature_generated_case_count": 15,
                "literature_generated_failed_case_count": 0,
                "gold_covered_cell_count": 28,
                "gold_missing_cell_count": 0,
                "gold_coverage_fraction": 1.0,
                "gold_specificity_over_recommended_count": 0,
                "gold_specificity_single_cell_fraction": 0.9,
                "batch_missing_case_count": 1,
                "scenario_traceability_failed_case_count": 0,
                "gold_positive_count": 20,
                "heldout_warning_count": 0,
                "scenario_precision_at_1_min": 1.0,
                "scenario_expected_recall_at_3_min": 1.0,
            }
        }
        current = {
            "metrics": {
                "literature_status_accuracy": 0.9,
                "literature_generated_case_count": 14,
                "literature_generated_failed_case_count": 1,
                "gold_covered_cell_count": 27,
                "gold_missing_cell_count": 1,
                "gold_coverage_fraction": 0.95,
                "gold_specificity_over_recommended_count": 1,
                "gold_specificity_single_cell_fraction": 0.8,
                "batch_missing_case_count": 2,
                "scenario_traceability_failed_case_count": 1,
                "gold_positive_count": 19,
                "heldout_warning_count": 1,
                "scenario_precision_at_1_min": 0.0,
                "scenario_expected_recall_at_3_min": 0.5,
            }
        }

        comparison = audit.compare_snapshots(current, baseline)
        by_metric = {row["metric"]: row for row in comparison["regressions"]}

        self.assertEqual(comparison["status"], "failed")
        self.assertIn("literature_status_accuracy", by_metric)
        self.assertIn("literature_generated_case_count", by_metric)
        self.assertIn("literature_generated_failed_case_count", by_metric)
        self.assertIn("gold_covered_cell_count", by_metric)
        self.assertIn("gold_missing_cell_count", by_metric)
        self.assertIn("gold_coverage_fraction", by_metric)
        self.assertIn("gold_specificity_over_recommended_count", by_metric)
        self.assertIn("gold_specificity_single_cell_fraction", by_metric)
        self.assertIn("batch_missing_case_count", by_metric)
        self.assertIn("scenario_traceability_failed_case_count", by_metric)
        self.assertIn("gold_positive_count", by_metric)
        self.assertIn("heldout_warning_count", by_metric)
        self.assertIn("scenario_precision_at_1_min", by_metric)
        self.assertIn("scenario_expected_recall_at_3_min", by_metric)

    def test_compare_passes_when_metrics_improve_or_hold(self):
        baseline = {"metrics": {"literature_status_accuracy": 0.9, "batch_missing_case_count": 2}}
        current = {"metrics": {"literature_status_accuracy": 1.0, "batch_missing_case_count": 1}}

        comparison = audit.compare_snapshots(current, baseline)

        self.assertEqual(comparison["status"], "passed")
        self.assertEqual(comparison["regressions"], [])


if __name__ == "__main__":
    unittest.main()
