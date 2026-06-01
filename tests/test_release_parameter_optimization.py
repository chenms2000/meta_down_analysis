import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_release_parameter_optimization.py"
SPEC = importlib.util.spec_from_file_location("audit_release_parameter_optimization", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)


class ReleaseParameterOptimizationTests(unittest.TestCase):
    def test_maps_release_gaps_to_parameter_actions_without_lowering_confidence(self):
        batch = {
            "fixture_audit": {
                "status": "failed",
                "missing_required": [
                    {
                        "case_id": "cscc_trait_score_csv",
                        "path": "missing.csv",
                        "expected_schema": {"expected_input_mode": "trait_score"},
                        "candidate_paths": [{"relative_path": "backup/cscc.csv"}],
                        "recovery_action": "restore fixture",
                        "recovery_rerun_plan": {
                            "case_id": "cscc_trait_score_csv",
                            "validation_sequence": [{"step_id": "fixture_schema_audit"}],
                        },
                    }
                ],
            },
            "benchmark_matrix": {
                "executed_gap_count": 2,
                "executed_gaps": ["cancers:cscc", "themes:arachidonate"],
            },
            "heldout_validation": {"status": "passed", "rerun_performed": False, "failed_checks": []},
            "scenario_traceability_contract": {"status": "passed", "failed_case_count": 0},
            "literature_semantics": {"status": "passed", "metrics": {"required_coverage_missing_count": 0}},
            "gold_coverage": {
                "status": "covered",
                "coverage_specificity": {
                    "over_recommended_count": 1,
                    "over_recommended_gold_ids": ["coad_broad_core"],
                },
            },
            "gold_quality": {"status": "passed"},
            "optimization_signals": {
                "by_code": [
                    {"code": "unsupported_claim_risk", "cases": ["case_a"]},
                    {"code": "gold_context_mismatch", "cases": ["case_b"]},
                ],
                "by_target": [
                    {
                        "bucket": "evidence_scoring",
                        "max_severity": "warning",
                        "signal_count": 1,
                        "cases": ["case_a"],
                        "codes": ["unsupported_claim_risk"],
                        "optimization_targets": ["evidence_scoring"],
                    }
                ],
                "by_dimension": [
                    {
                        "bucket": "input_modes:compound_name",
                        "max_severity": "warning",
                        "signal_count": 1,
                        "cases": ["case_a"],
                        "codes": ["unsupported_claim_risk"],
                        "optimization_targets": ["evidence_scoring"],
                    }
                ],
            },
        }
        readiness = {
            "status": "blocked",
            "blockers": [{"code": "fixture_audit_failed"}, {"code": "heldout_not_rerun"}],
        }
        objective = {
            "status": "incomplete",
            "remaining_blockers": [{"requirement_id": "cscc_real_csv_rerun_gate"}],
        }
        metrics = {
            "metrics": {
                "benchmark_matrix_executed_gap_count": 2,
                "heldout_rerun_performed": 0,
                "objective_incomplete_count": 3,
                "scenario_precision_at_1_min": 0.5,
                "scenario_recall_at_5_min": 0.25,
                "scenario_expected_topk_eligible_count": 4,
                "scenario_expected_recall_at_5_min": 0.75,
            }
        }

        report = audit.audit_parameter_optimization(batch, readiness, objective, metrics)
        by_id = {row["action_id"]: row for row in report["actions"]}

        self.assertEqual(report["status"], "action_required")
        self.assertIn("benchmark_matrix_execute_required_dimensions", by_id)
        self.assertIn("fixture_schema_restore_required_inputs", by_id)
        self.assertIn("heldout_rerun_enforce_current_validation", by_id)
        self.assertIn("evidence_support_requirement_tighten", by_id)
        self.assertIn("context_weight_recalibrate", by_id)
        self.assertIn("gold_standard_split_overbroad_controls", by_id)
        self.assertIn("objective_incomplete_release_block", by_id)
        self.assertIn("scenario_expected_topk_rank_calibration_review", by_id)
        self.assertTrue(by_id["benchmark_matrix_execute_required_dimensions"]["required"])
        self.assertIn("real_world_fixture_registry", by_id["fixture_schema_restore_required_inputs"]["parameter_targets"])
        self.assertEqual(
            by_id["fixture_schema_restore_required_inputs"]["fixture_recovery"][0]["expected_schema"]["expected_input_mode"],
            "trait_score",
        )
        self.assertEqual(
            by_id["fixture_schema_restore_required_inputs"]["fixture_recovery"][0]["candidate_paths"][0]["relative_path"],
            "backup/cscc.csv",
        )
        self.assertEqual(
            by_id["fixture_schema_restore_required_inputs"]["fixture_recovery"][0]["recovery_rerun_plan"]["case_id"],
            "cscc_trait_score_csv",
        )
        self.assertIn("heldout_transfer_quality_thresholds", by_id["heldout_rerun_enforce_current_validation"]["parameter_targets"])
        self.assertIn("source_family_holdout_strategy", by_id["heldout_rerun_enforce_current_validation"]["parameter_targets"])
        self.assertIn("literature_source_partitioning", by_id["heldout_rerun_enforce_current_validation"]["parameter_targets"])
        self.assertIn("context_weight", by_id["context_weight_recalibrate"]["parameter_targets"])
        self.assertFalse(by_id["gold_standard_split_overbroad_controls"]["required"])
        self.assertIn("coad_broad_core", by_id["gold_standard_split_overbroad_controls"]["evidence_refs"])
        self.assertFalse(by_id["scenario_expected_topk_rank_calibration_review"]["required"])
        self.assertIn("gold_entity_alias_map", by_id["scenario_expected_topk_rank_calibration_review"]["parameter_targets"])
        self.assertIn("Do not lower confidence thresholds", report["policy"]["confidence_policy"])
        self.assertEqual(report["release_metric_context"]["scenario_precision_at_1_min"], 0.5)
        self.assertEqual(report["release_metric_context"]["scenario_recall_at_5_min"], 0.25)
        self.assertEqual(report["release_metric_context"]["scenario_expected_topk_eligible_count"], 4)
        self.assertEqual(report["release_metric_context"]["scenario_expected_recall_at_5_min"], 0.75)
        matrix = report["diagnostics"]["optimization_signal_matrix"]
        self.assertEqual(matrix["top_targets"][0]["bucket"], "evidence_scoring")
        self.assertEqual(matrix["top_dimensions"][0]["bucket"], "input_modes:compound_name")

    def test_no_action_required_when_all_gates_are_clean(self):
        batch = {
            "fixture_audit": {"status": "passed"},
            "benchmark_matrix": {"executed_gap_count": 0, "executed_gaps": []},
            "heldout_validation": {"status": "passed", "rerun_performed": True, "warnings": []},
            "scenario_traceability_contract": {"status": "passed"},
            "literature_semantics": {"status": "passed", "metrics": {"required_coverage_missing_count": 0}},
            "gold_coverage": {"status": "covered"},
            "gold_quality": {"status": "passed"},
            "optimization_signals": {"by_code": []},
        }

        report = audit.audit_parameter_optimization(
            batch,
            {"status": "ready", "blockers": []},
            {"status": "complete", "remaining_blockers": []},
            {"metrics": {}},
        )

        self.assertEqual(report["status"], "no_action_required")
        self.assertEqual(report["actions"], [])

    def test_sparse_expected_topk_precision_does_not_trigger_rank_action_when_recall_is_complete(self):
        batch = {
            "fixture_audit": {"status": "passed"},
            "benchmark_matrix": {"executed_gap_count": 0, "executed_gaps": []},
            "heldout_validation": {"status": "passed", "rerun_performed": True, "warnings": []},
            "scenario_traceability_contract": {"status": "passed"},
            "literature_semantics": {"status": "passed", "metrics": {"required_coverage_missing_count": 0}},
            "gold_coverage": {"status": "covered"},
            "gold_quality": {"status": "passed"},
            "optimization_signals": {"by_code": []},
        }
        metrics = {
            "metrics": {
                "scenario_expected_topk_eligible_count": 6,
                "scenario_expected_recall_at_5_min": 1.0,
                "scenario_expected_precision_at_5_min": 0.2,
            }
        }

        report = audit.audit_parameter_optimization(
            batch,
            {"status": "ready", "blockers": []},
            {"status": "complete", "remaining_blockers": []},
            metrics,
        )

        by_id = {row["action_id"]: row for row in report["actions"]}
        self.assertNotIn("scenario_expected_topk_rank_calibration_review", by_id)
        self.assertEqual(report["status"], "no_action_required")

    def test_heldout_warning_becomes_recommended_calibration_action(self):
        batch = {
            "fixture_audit": {"status": "passed"},
            "benchmark_matrix": {"executed_gap_count": 0, "executed_gaps": []},
            "heldout_validation": {
                "status": "passed",
                "rerun_performed": True,
                "warnings": [{"code": "source_heldout_random_auc"}],
            },
            "scenario_traceability_contract": {"status": "passed"},
            "literature_semantics": {"status": "passed", "metrics": {"required_coverage_missing_count": 0}},
            "gold_coverage": {"status": "covered"},
            "gold_quality": {"status": "passed"},
            "optimization_signals": {"by_code": []},
        }

        report = audit.audit_parameter_optimization(
            batch,
            {"status": "ready", "blockers": []},
            {"status": "complete", "remaining_blockers": []},
            {"metrics": {"heldout_warning_count": 1}},
        )
        by_id = {row["action_id"]: row for row in report["actions"]}
        action = by_id["heldout_topk_global_rank_calibration_review"]

        self.assertEqual(report["status"], "action_required")
        self.assertFalse(action["required"])
        self.assertIn("global_auc_vs_topk_policy", action["parameter_targets"])
        self.assertIn("source_heldout_random_auc", action["evidence_refs"])

    def test_random_auc_warning_is_diagnostic_only_when_topk_transfer_passes(self):
        source_label = "literature:partition_0"
        batch = {
            "fixture_audit": {"status": "passed"},
            "benchmark_matrix": {"executed_gap_count": 0, "executed_gaps": []},
            "heldout": {
                "source_heldout": {
                    "summary": [
                        {
                            "source": source_label,
                            "validation_metrics": {
                                "heldout_lift_at_100": 4.5,
                                "heldout_unique_transfer_scores": 300,
                                "heldout_top_score_tie_fraction": 0.01,
                            },
                        }
                    ]
                }
            },
            "heldout_validation": {
                "status": "passed",
                "rerun_performed": True,
                "warnings": [
                    {
                        "code": "source_heldout_random_auc",
                        "metrics": {
                            "label": source_label,
                            "heldout_roc_auc": 0.41,
                            "diagnostics": {
                                "feature_transfer": {"shared_nonconstant_feature_count": 12},
                                "source_family": {"same_family_available": True},
                            },
                        },
                    }
                ],
            },
            "scenario_traceability_contract": {"status": "passed"},
            "literature_semantics": {"status": "passed", "metrics": {"required_coverage_missing_count": 0}},
            "gold_coverage": {"status": "covered"},
            "gold_quality": {"status": "passed"},
            "optimization_signals": {"by_code": []},
        }

        report = audit.audit_parameter_optimization(
            batch,
            {"status": "ready", "blockers": []},
            {"status": "complete", "remaining_blockers": []},
            {"metrics": {"heldout_warning_count": 1}},
        )

        by_id = {row["action_id"]: row for row in report["actions"]}
        self.assertNotIn("heldout_topk_global_rank_calibration_review", by_id)
        self.assertEqual(report["status"], "no_action_required")
        diagnostic = report["diagnostics"]["heldout_warning_diagnostics"][0]
        self.assertFalse(diagnostic["requires_action"])
        self.assertEqual(diagnostic["reason"], "global_auc_low_but_topk_transfer_and_feature_surface_passed")

    def test_negative_trap_becomes_required_ranking_recalibration_action(self):
        batch = {
            "fixture_audit": {"status": "passed"},
            "benchmark_matrix": {"executed_gap_count": 0, "executed_gaps": []},
            "heldout_validation": {"status": "passed", "rerun_performed": True},
            "scenario_traceability_contract": {"status": "passed"},
            "literature_semantics": {"status": "passed", "metrics": {"required_coverage_missing_count": 0}},
            "gold_coverage": {"status": "covered"},
            "gold_quality": {"status": "passed"},
            "optimization_signals": {
                "by_code": [{"code": "negative_trap_hit", "cases": ["trap_case"]}],
            },
        }

        report = audit.audit_parameter_optimization(
            batch,
            {"status": "ready", "blockers": []},
            {"status": "complete", "remaining_blockers": []},
            {"metrics": {}},
        )
        action = {row["action_id"]: row for row in report["actions"]}["negative_trap_blocking_recalibrate"]

        self.assertTrue(action["required"])
        self.assertIn("graph_propagation_penalty", action["parameter_targets"])
        self.assertIn("trap_case", action["evidence_refs"])

    def test_signal_matrix_summary_prioritizes_error_targets_and_dimensions(self):
        summary = audit.signal_matrix_summary(
            {
                "by_target": [
                    {
                        "bucket": "resolver",
                        "max_severity": "warning",
                        "signal_count": 4,
                        "cases": ["case_a"],
                        "codes": ["identity_ambiguity_review"],
                    },
                    {
                        "bucket": "database_schema",
                        "max_severity": "error",
                        "signal_count": 1,
                        "cases": ["fixture"],
                        "codes": ["required_fixture_missing"],
                    },
                ],
                "by_dimension": [
                    {
                        "bucket": "cancers:cscc",
                        "max_severity": "error",
                        "signal_count": 2,
                        "cases": ["fixture", "cancers:cscc"],
                        "codes": ["fixture_missing", "benchmark_matrix_gap"],
                        "optimization_targets": ["database_schema", "benchmark_runner"],
                    }
                ],
            }
        )

        self.assertEqual(summary["status"], "present")
        self.assertEqual(summary["top_targets"][0]["bucket"], "database_schema")
        self.assertEqual(summary["error_dimensions"][0]["bucket"], "cancers:cscc")
        self.assertIn("benchmark_runner", summary["error_dimensions"][0]["optimization_targets"])

    def test_markdown_summary_includes_fixture_recovery_details(self):
        report = {
            "status": "action_required",
            "action_count": 1,
            "policy": {"confidence_policy": "Do not lower confidence thresholds."},
            "actions": [
                {
                    "action_id": "fixture_schema_restore_required_inputs",
                    "category": "database_schema",
                    "required": True,
                    "rationale": "Restore fixture schema.",
                    "parameter_targets": ["fixture_identity_columns_any"],
                    "evidence_refs": ["cscc_trait_score_csv"],
                    "fixture_recovery": [
                        {
                            "case_id": "cscc_trait_score_csv",
                            "expected_schema": {
                                "expected_input_mode": "trait_score",
                                "identity_columns_any": ["trait"],
                                "effect_columns_any": ["cohen_d"],
                            },
                            "candidate_paths": [],
                            "recovery_rerun_plan": {"validation_sequence": [{"step_id": "fixture_schema_audit"}]},
                        }
                    ],
                }
            ],
        }

        markdown = audit.markdown_summary(report)

        self.assertIn("expected mode `trait_score`", markdown)
        self.assertIn("effect columns `cohen_d`", markdown)
        self.assertIn("Recovery rerun steps: `1`", markdown)

    def test_markdown_summary_includes_optimization_signal_matrix(self):
        report = {
            "status": "action_required",
            "action_count": 0,
            "policy": {"confidence_policy": "Do not lower confidence thresholds."},
            "actions": [],
            "diagnostics": {
                "optimization_signal_matrix": {
                    "top_targets": [
                        {
                            "bucket": "database_schema",
                            "max_severity": "error",
                            "codes": ["required_fixture_missing"],
                            "case_count": 1,
                        }
                    ],
                    "top_dimensions": [
                        {
                            "bucket": "cancers:cscc",
                            "max_severity": "error",
                            "optimization_targets": ["database_schema"],
                            "codes": ["benchmark_matrix_gap"],
                            "case_count": 2,
                        }
                    ],
                }
            },
        }

        markdown = audit.markdown_summary(report)

        self.assertIn("Optimization Signal Matrix", markdown)
        self.assertIn("| database_schema | error | required_fixture_missing | 1 |", markdown)
        self.assertIn("cancers:cscc", markdown)


if __name__ == "__main__":
    unittest.main()
