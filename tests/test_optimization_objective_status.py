import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_optimization_objective_status.py"
SPEC = importlib.util.spec_from_file_location("audit_optimization_objective_status", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)


def complete_batch_report():
    return {
        "gold_coverage": {"status": "covered"},
        "gold_quality": {
            "status": "passed",
            "positive_count": 20,
            "negative_control_count": 10,
            "positive_distribution": {
                "context_specific_positive_count": 12,
                "unique_positive_cancer_context_count": 7,
                "unique_positive_mechanism_axis_count": 18,
                "missing_required_themes": [],
            },
        },
        "literature_semantics": {
            "status": "passed",
            "metrics": {
                "case_count": 28,
                "minimum_case_count": 28,
                "status_accuracy": 1.0,
                "cue_recall": 1.0,
                "method_cue_recall": 1.0,
                "required_coverage_missing_count": 0,
            },
            "coverage": {"support_statuses": {"status": "covered"}},
        },
        "heldout_validation": {
            "status": "passed",
            "rerun_performed": True,
            "rerun_passed": {"temporal": True, "source": True},
            "failed_checks": [],
        },
        "benchmark_matrix": {
            "status": "passed",
            "registered_gap_count": 0,
            "executed_gap_count": 0,
            "registered_gaps": [],
            "executed_gaps": [],
            "dimensions": {
                "cancers": {"executed": ["cscc"]},
                "themes": {"executed": ["arachidonate"]},
            },
        },
        "scenario_traceability_contract": {"status": "passed", "failed_case_count": 0, "failure_codes": []},
        "status": {"passed": True, "traceability_contract_passed": True},
    }


class OptimizationObjectiveStatusTests(unittest.TestCase):
    def test_objective_status_complete_when_all_gates_are_proven(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cscc = root / "required_metabolism_fixture.csv"
            cscc.write_text("trait,cohen_d\nGCST1,1.0\n", encoding="utf-8")

            report = audit.audit_objective_status(
                root,
                complete_batch_report(),
                {"status": "ready", "blockers": []},
                {"status": "ready", "blockers": []},
                cscc,
                {"status": "no_action_required", "action_count": 0, "severity_counts": {"required": 0}},
            )

        self.assertEqual(report["status"], "complete")
        self.assertEqual(report["incomplete_count"], 0)

    def test_objective_status_incomplete_when_cscc_fixture_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            batch = complete_batch_report()
            batch["benchmark_matrix"] = {
                "status": "failed",
                "registered_gap_count": 0,
                "executed_gap_count": 2,
                "registered_gaps": [],
                "executed_gaps": ["cancers:cscc", "themes:arachidonate"],
                "dimensions": {
                    "cancers": {"executed": []},
                    "themes": {"executed": []},
                },
            }
            batch["status"] = {"passed": False}
            batch["fixture_audit"] = {
                "status": "failed",
                "missing_required": [
                    {
                        "case_id": "cscc_trait_score_csv",
                        "path": str(root / "missing.csv"),
                        "expected_schema": {"expected_input_mode": "trait_score"},
                        "candidate_paths": [{"relative_path": "backup/cscc.csv"}],
                        "recovery_action": "restore fixture",
                        "recovery_rerun_plan": {
                            "case_id": "cscc_trait_score_csv",
                            "validation_sequence": [{"step_id": "fixture_schema_audit"}],
                        },
                    }
                ],
            }

            report = audit.audit_objective_status(
                root,
                batch,
                {"status": "blocked", "blockers": [{"code": "benchmark_matrix_failed"}]},
                {"status": "blocked", "blockers": [{"code": "real_world_fixture_audit_failed"}]},
                root / "missing.csv",
                {
                    "status": "action_required",
                    "action_count": 1,
                    "severity_counts": {"required": 1},
                    "actions": [{"action_id": "benchmark_matrix_execute_required_dimensions"}],
                },
            )

        statuses = {row["requirement_id"]: row["status"] for row in report["requirements"]}
        blockers = {row["requirement_id"]: row["blockers"] for row in report["requirements"]}
        evidence = {row["requirement_id"]: row["evidence"] for row in report["requirements"]}

        self.assertEqual(report["status"], "incomplete")
        self.assertEqual(statuses["gold_standard_full_library_gate"], "satisfied")
        self.assertEqual(statuses["cscc_real_csv_rerun_gate"], "incomplete")
        self.assertEqual(statuses["release_batch_benchmark_matrix_gate"], "incomplete")
        self.assertIn("required cSCC CSV is missing or its scenario did not execute", blockers["cscc_real_csv_rerun_gate"])
        self.assertEqual(
            evidence["cscc_real_csv_rerun_gate"]["fixture_recovery"]["expected_schema"]["expected_input_mode"],
            "trait_score",
        )
        self.assertEqual(
            evidence["cscc_real_csv_rerun_gate"]["fixture_recovery"]["candidate_paths"][0]["relative_path"],
            "backup/cscc.csv",
        )
        self.assertEqual(
            evidence["cscc_real_csv_rerun_gate"]["fixture_recovery"]["recovery_rerun_plan"]["case_id"],
            "cscc_trait_score_csv",
        )

    def test_cscc_fixture_gate_is_not_coupled_to_arachidonate_theme(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cscc = root / "required_metabolism_fixture.csv"
            cscc.write_text("trait,cohen_d\nGCST1,1.0\n", encoding="utf-8")
            batch = complete_batch_report()
            batch["benchmark_matrix"]["dimensions"]["themes"]["executed"] = []

            report = audit.audit_objective_status(
                root,
                batch,
                {"status": "ready", "blockers": []},
                {"status": "ready", "blockers": []},
                cscc,
                {"status": "no_action_required", "action_count": 0, "severity_counts": {"required": 0}},
            )

        cscc_gate = next(row for row in report["requirements"] if row["requirement_id"] == "cscc_real_csv_rerun_gate")

        self.assertEqual(cscc_gate["status"], "satisfied")
        self.assertFalse(cscc_gate["evidence"]["arachidonate_executed_in_matrix"])
        self.assertEqual(
            cscc_gate["evidence"]["arachidonate_gate"],
            "covered_by_release_batch_benchmark_matrix_gate",
        )

    def test_objective_status_requires_scenario_traceability_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cscc = root / "required_metabolism_fixture.csv"
            cscc.write_text("trait,cohen_d\nGCST1,1.0\n", encoding="utf-8")
            batch = complete_batch_report()
            batch["scenario_traceability_contract"] = {
                "status": "failed",
                "failed_case_count": 1,
                "failure_codes": ["evidence_assertions_missing"],
            }
            batch["status"] = {"passed": False, "traceability_contract_passed": False}

            report = audit.audit_objective_status(
                root,
                batch,
                {"status": "blocked", "blockers": [{"code": "scenario_traceability_contract_failed"}]},
                {"status": "blocked", "blockers": [{"code": "objective_status_incomplete"}]},
                cscc,
                {
                    "status": "action_required",
                    "action_count": 1,
                    "severity_counts": {"required": 1},
                    "actions": [{"action_id": "scenario_traceability_contract_repair"}],
                },
            )

        matrix_gate = next(row for row in report["requirements"] if row["requirement_id"] == "release_batch_benchmark_matrix_gate")

        self.assertEqual(matrix_gate["status"], "incomplete")
        self.assertEqual(matrix_gate["evidence"]["scenario_traceability_status"], "failed")
        self.assertIn("evidence_assertions_missing", matrix_gate["evidence"]["scenario_traceability_failure_codes"])

    def test_objective_status_requires_parameter_optimization_audit_to_be_clean(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cscc = root / "required_metabolism_fixture.csv"
            cscc.write_text("trait,cohen_d\nGCST1,1.0\n", encoding="utf-8")

            report = audit.audit_objective_status(
                root,
                complete_batch_report(),
                {"status": "blocked", "blockers": [{"code": "parameter_optimization_action_required"}]},
                {"status": "blocked", "blockers": [{"code": "parameter_optimization_action_required"}]},
                cscc,
                {
                    "status": "action_required",
                    "action_count": 1,
                    "severity_counts": {"required": 1},
                    "actions": [{"action_id": "heldout_rerun_enforce_current_validation"}],
                },
            )

        parameter_gate = next(row for row in report["requirements"] if row["requirement_id"] == "parameter_optimization_gate")

        self.assertEqual(parameter_gate["status"], "incomplete")
        self.assertEqual(parameter_gate["evidence"]["required_action_count"], 1)
        self.assertEqual(parameter_gate["evidence"]["blocking_required_action_count"], 1)
        self.assertIn("heldout_rerun_enforce_current_validation", parameter_gate["evidence"]["action_ids"])
        self.assertIn("heldout_rerun_enforce_current_validation", parameter_gate["evidence"]["blocking_action_ids"])

    def test_objective_status_does_not_treat_self_release_guard_as_parameter_blocker(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cscc = root / "required_metabolism_fixture.csv"
            cscc.write_text("trait,cohen_d\nGCST1,1.0\n", encoding="utf-8")

            report = audit.audit_objective_status(
                root,
                complete_batch_report(),
                {"status": "ready", "blockers": []},
                {"status": "ready", "blockers": []},
                cscc,
                {
                    "status": "action_required",
                    "action_count": 1,
                    "severity_counts": {"required": 1},
                    "actions": [
                        {
                            "action_id": "objective_incomplete_release_block",
                            "required": True,
                        }
                    ],
                },
            )

        parameter_gate = next(row for row in report["requirements"] if row["requirement_id"] == "parameter_optimization_gate")

        self.assertEqual(parameter_gate["status"], "satisfied")
        self.assertEqual(parameter_gate["evidence"]["required_action_count"], 1)
        self.assertEqual(parameter_gate["evidence"]["blocking_required_action_count"], 0)
        self.assertEqual(parameter_gate["evidence"]["blocking_action_ids"], [])

    def test_objective_status_requires_executed_matrix_coverage_even_when_exploratory_matrix_passed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cscc = root / "required_metabolism_fixture.csv"
            cscc.write_text("trait,cohen_d\nGCST1,1.0\n", encoding="utf-8")
            batch = complete_batch_report()
            batch["benchmark_matrix"]["status"] = "passed"
            batch["benchmark_matrix"]["executed_gap_count"] = 1
            batch["benchmark_matrix"]["executed_gaps"] = ["cancers:cscc"]
            batch["status"] = {"passed": True, "traceability_contract_passed": True}

            report = audit.audit_objective_status(
                root,
                batch,
                {"status": "blocked", "blockers": [{"code": "benchmark_matrix_failed"}]},
                {"status": "blocked", "blockers": [{"code": "objective_status_incomplete"}]},
                cscc,
                {"status": "no_action_required", "action_count": 0, "severity_counts": {"required": 0}},
            )

        matrix_gate = next(row for row in report["requirements"] if row["requirement_id"] == "release_batch_benchmark_matrix_gate")

        self.assertEqual(matrix_gate["status"], "incomplete")
        self.assertEqual(matrix_gate["evidence"]["executed_gap_count"], 1)

    def test_markdown_summary_lists_blockers(self):
        report = {
            "status": "incomplete",
            "satisfied_count": 1,
            "incomplete_count": 1,
            "requirements": [
                {
                    "requirement_id": "cscc_real_csv_rerun_gate",
                    "status": "incomplete",
                    "blockers": ["missing csv"],
                    "evidence": {
                        "fixture_recovery": {
                            "case_id": "cscc_trait_score_csv",
                            "expected_schema": {
                                "expected_input_mode": "trait_score",
                                "identity_columns_any": ["trait"],
                                "effect_columns_any": ["cohen_d"],
                            },
                            "candidate_paths": [],
                            "recovery_rerun_plan": {"validation_sequence": [{"step_id": "fixture_schema_audit"}]},
                        }
                    },
                },
            ],
            "interpretation_boundaries": ["research only"],
        }

        markdown = audit.markdown_summary(report)

        self.assertIn("Status: `incomplete`", markdown)
        self.assertIn("cscc_real_csv_rerun_gate", markdown)
        self.assertIn("missing csv", markdown)
        self.assertIn("expected mode `trait_score`", markdown)
        self.assertIn("Recovery rerun steps: `1`", markdown)


if __name__ == "__main__":
    unittest.main()

