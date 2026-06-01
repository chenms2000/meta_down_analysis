import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_release_readiness.py"
SPEC = importlib.util.spec_from_file_location("audit_release_readiness", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)


def ready_report():
    return {
        "release_profile": {"profile": "strict"},
        "status": {"passed": True, "failed_case_count": 0, "missing_case_count": 0},
        "fixture_audit": {"status": "passed"},
        "benchmark_matrix": {"status": "passed", "registered_gap_count": 0, "executed_gap_count": 0},
        "scenario_traceability_contract": {"status": "passed", "failed_case_count": 0},
        "heldout_validation": {"status": "passed", "rerun_performed": True},
        "gold_coverage": {"status": "covered"},
        "gold_quality": {"status": "passed"},
        "literature_semantics": {"status": "passed"},
        "optimization_signals": {"by_code": []},
    }


class ReleaseReadinessAuditTests(unittest.TestCase):
    def test_ready_when_all_strict_gates_pass(self):
        report = audit.audit_readiness(ready_report())

        self.assertEqual(report["status"], "ready")
        self.assertEqual(report["blocker_count"], 0)
        self.assertEqual(report["gate_summary"]["release_profile"], "strict")
        self.assertEqual(report["gate_summary"]["benchmark_matrix_status"], "passed")
        self.assertEqual(report["gate_summary"]["scenario_traceability_status"], "passed")
        self.assertEqual(report["remediation_actions"], [])

    def test_blocks_on_missing_fixture_and_error_signal(self):
        payload = ready_report()
        payload["status"] = {"passed": False, "failed_case_count": 0, "missing_case_count": 1}
        payload["fixture_audit"] = {
            "status": "failed",
            "missing_required": [
                {
                    "case_id": "cscc_trait_score_csv",
                    "path": "missing.csv",
                    "expected_schema": {"expected_input_mode": "trait_score"},
                    "candidate_paths": [{"relative_path": "backup/cscc.csv"}],
                    "recovery_action": "restore fixture",
                    "recovery_rerun_plan": {"case_id": "cscc_trait_score_csv", "validation_sequence": [{"step_id": "fixture_schema_audit"}]},
                }
            ],
        }
        payload["optimization_signals"] = {
            "by_code": [{"code": "required_fixture_missing", "severity": "error", "cases": ["cscc_trait_score_csv"]}]
        }

        report = audit.audit_readiness(payload)
        codes = {row["code"] for row in report["blockers"]}

        self.assertEqual(report["status"], "blocked")
        self.assertIn("batch_benchmark_failed", codes)
        self.assertIn("fixture_audit_failed", codes)
        self.assertIn("error_optimization_signals", codes)
        actions = {row["code"]: row for row in report["remediation_actions"]}
        self.assertEqual(actions["fixture_audit_failed"]["priority"], "P0")
        self.assertIn("cscc_trait_score_csv", actions["fixture_audit_failed"]["evidence_refs"])
        self.assertEqual(
            actions["fixture_audit_failed"]["fixture_recovery"][0]["expected_schema"]["expected_input_mode"],
            "trait_score",
        )
        self.assertEqual(
            actions["fixture_audit_failed"]["fixture_recovery"][0]["candidate_paths"][0]["relative_path"],
            "backup/cscc.csv",
        )
        self.assertEqual(
            actions["fixture_audit_failed"]["fixture_recovery"][0]["recovery_rerun_plan"]["case_id"],
            "cscc_trait_score_csv",
        )

    def test_blocks_on_strict_benchmark_matrix_gap(self):
        payload = ready_report()
        payload["status"] = {"passed": False, "failed_case_count": 0, "missing_case_count": 1}
        payload["benchmark_matrix"] = {
            "status": "failed",
            "registered_gap_count": 0,
            "executed_gap_count": 2,
            "registered_gaps": [],
            "executed_gaps": ["cancers:cscc", "themes:arachidonate"],
            "cases": [{"case_id": "cscc_trait_score_csv", "status": "missing"}],
        }

        report = audit.audit_readiness(payload)
        codes = {row["code"] for row in report["blockers"]}

        self.assertEqual(report["status"], "blocked")
        self.assertIn("benchmark_matrix_failed", codes)
        self.assertEqual(report["gate_summary"]["benchmark_matrix_status"], "failed")
        self.assertEqual(report["gate_summary"]["benchmark_matrix_executed_gap_count"], 2)
        actions = {row["code"]: row for row in report["remediation_actions"]}
        self.assertEqual(actions["benchmark_matrix_failed"]["priority"], "P0")
        self.assertIn("cancers:cscc", actions["benchmark_matrix_failed"]["evidence_refs"])
        self.assertIn("themes:arachidonate", actions["benchmark_matrix_failed"]["evidence_refs"])
        self.assertNotIn("{'case_id': 'cscc_trait_score_csv', 'status': 'missing'}", actions["benchmark_matrix_failed"]["evidence_refs"])

    def test_strict_readiness_requires_zero_executed_matrix_gaps_even_if_matrix_status_passed(self):
        payload = ready_report()
        payload["benchmark_matrix"] = {
            "status": "passed",
            "registered_gap_count": 0,
            "executed_gap_count": 1,
            "registered_gaps": [],
            "executed_gaps": ["cancers:cscc"],
        }

        report = audit.audit_readiness(payload)

        self.assertEqual(report["status"], "blocked")
        self.assertIn("benchmark_matrix_failed", {row["code"] for row in report["blockers"]})

    def test_blocks_on_scenario_traceability_contract_gap(self):
        payload = ready_report()
        payload["status"] = {"passed": False, "failed_case_count": 0, "missing_case_count": 0}
        payload["scenario_traceability_contract"] = {
            "status": "failed",
            "failed_case_count": 1,
            "failed_cases": [
                {
                    "case_id": "unstructured_case",
                    "failures": [{"code": "evidence_assertions_missing"}],
                }
            ],
        }

        report = audit.audit_readiness(payload)
        codes = {row["code"] for row in report["blockers"]}
        actions = {row["code"]: row for row in report["remediation_actions"]}

        self.assertEqual(report["status"], "blocked")
        self.assertIn("scenario_traceability_contract_failed", codes)
        self.assertEqual(report["gate_summary"]["scenario_traceability_status"], "failed")
        self.assertIn("unstructured_case", actions["scenario_traceability_contract_failed"]["evidence_refs"])

    def test_blocks_when_strict_profile_not_used(self):
        payload = ready_report()
        payload["release_profile"] = {"profile": "exploratory"}

        report = audit.audit_readiness(payload)

        self.assertEqual(report["status"], "blocked")
        self.assertIn("strict_profile_not_used", {row["code"] for row in report["blockers"]})

    def test_allows_exploratory_when_explicitly_requested(self):
        payload = ready_report()
        payload["release_profile"] = {"profile": "exploratory"}

        report = audit.audit_readiness(payload, require_strict_profile=False)

        self.assertEqual(report["status"], "ready")

    def test_markdown_summary_includes_blockers_and_actions(self):
        payload = ready_report()
        payload["status"] = {"passed": False, "failed_case_count": 0, "missing_case_count": 1}
        payload["fixture_audit"] = {
            "status": "failed",
            "missing_required": [
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

        report = audit.audit_readiness(payload)
        markdown = audit.markdown_summary(report)

        self.assertIn("Status: `blocked`", markdown)
        self.assertIn("fixture_audit_failed", markdown)
        self.assertIn("Remediation Actions", markdown)
        self.assertIn("expected mode `trait_score`", markdown)
        self.assertIn("effect columns `cohen_d`", markdown)
        self.assertIn("Recovery rerun steps: `1`", markdown)


if __name__ == "__main__":
    unittest.main()
