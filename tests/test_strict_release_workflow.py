import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_strict_release_workflow.py"
SPEC = importlib.util.spec_from_file_location("run_strict_release_workflow", SCRIPT)
workflow = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = workflow
SPEC.loader.exec_module(workflow)


def write_registry(root: Path, fixture_path: str) -> Path:
    registry = {
        "schema_version": "test",
        "fixtures": [
            {
                "case_id": "cscc_trait_score_csv",
                "path": fixture_path,
                "fixture_type": "csv_trait_score",
                "required_for_strict_release": True,
                "min_rows": 1,
                "identity_columns_any": ["trait", "reported_trait"],
                "effect_columns_any": ["cohen_d", "log2FC"],
                "expected_input_mode": "trait_score",
            }
        ],
    }
    config_dir = root / "config"
    config_dir.mkdir(parents=True, exist_ok=True)
    path = config_dir / "real_world_fixture_registry.json"
    path.write_text(json.dumps(registry), encoding="utf-8")
    return path


class StrictReleaseWorkflowTests(unittest.TestCase):
    def test_missing_required_fixture_blocks_before_batch_benchmark(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_registry(root, "missing.csv")

            args = workflow.parse_args(
                [
                    "--workspace",
                    str(root),
                    "--release-id",
                    "test_release",
                    "--output-root",
                    "validation_reports",
                    "--fixture-registry",
                    "config/real_world_fixture_registry.json",
                ]
            )
            report = workflow.run_workflow(args)

            self.assertEqual(report["status"], "blocked")
            self.assertEqual(report["fixture_audit"]["status"], "failed")
            self.assertFalse(report["gate_summary"]["batch_benchmark_executed"])
            self.assertTrue(report["batch_benchmark"]["skipped"])
            self.assertFalse(report["post_release_audits"]["skipped"])
            self.assertEqual(report["post_release_audits"]["objective_status"]["status"], "incomplete")
            self.assertIn(report["post_release_audits"]["metric_regression"]["status"], {"snapshot_only", "passed"})
            self.assertEqual(report["post_release_audits"]["parameter_optimization"]["status"], "action_required")
            self.assertGreater(report["post_release_audits"]["parameter_optimization"]["blocking_required_action_count"], 0)
            blocker_codes = {row["code"] for row in report["blockers"]}
            self.assertIn("real_world_fixture_audit_failed", blocker_codes)
            self.assertIn("objective_status_incomplete", blocker_codes)
            self.assertIn("parameter_optimization_action_required", blocker_codes)
            self.assertTrue(Path(report["outputs"]["workflow_json"]).exists())
            self.assertTrue(Path(report["outputs"]["fixture_audit_json"]).exists())
            self.assertTrue(Path(report["outputs"]["objective_status_json"]).exists())
            self.assertTrue(Path(report["outputs"]["metric_regression_json"]).exists())
            self.assertTrue(Path(report["outputs"]["parameter_optimization_json"]).exists())
            objective_report = json.loads(Path(report["outputs"]["objective_status_json"]).read_text(encoding="utf-8"))
            readiness_gate = next(
                row
                for row in objective_report["requirements"]
                if row["requirement_id"] == "strict_release_readiness_gate"
            )
            self.assertIn("objective_status_incomplete", readiness_gate["evidence"]["workflow_blocker_codes"])
            self.assertIn("parameter_optimization_action_required", readiness_gate["evidence"]["workflow_blocker_codes"])

    def test_skip_post_audits_leaves_terminal_audits_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_registry(root, "missing.csv")

            args = workflow.parse_args(
                [
                    "--workspace",
                    str(root),
                    "--release-id",
                    "test_release",
                    "--output-root",
                    "validation_reports",
                    "--fixture-registry",
                    "config/real_world_fixture_registry.json",
                    "--skip-post-audits",
                ]
            )
            report = workflow.run_workflow(args)

            self.assertEqual(report["status"], "blocked")
            self.assertTrue(report["post_release_audits"]["skipped"])
            self.assertNotIn("objective_status_incomplete", {row["code"] for row in report["blockers"]})
            self.assertFalse(Path(report["outputs"]["objective_status_json"]).exists())
            self.assertFalse(Path(report["outputs"]["parameter_optimization_json"]).exists())

    def test_fixture_blockers_empty_when_fixture_audit_passes(self):
        report = {"status": "passed", "failed_required": []}

        self.assertEqual(workflow.fixture_blockers(report), [])

    def test_parameter_summary_excludes_self_release_guard_from_blocking_actions(self):
        summary = workflow.parameter_optimization_summary(
            {
                "status": "action_required",
                "action_count": 1,
                "severity_counts": {"required": 1},
                "actions": [{"action_id": "objective_incomplete_release_block", "required": True}],
            }
        )

        self.assertEqual(summary["required_action_count"], 1)
        self.assertEqual(summary["blocking_required_action_count"], 0)
        self.assertEqual(summary["blocking_action_ids"], [])

    def test_workflow_markdown_records_blockers_and_outputs(self):
        markdown = workflow.workflow_markdown(
            {
                "status": "blocked",
                "workflow_version": "test",
                "release_id": "release_a",
                "gate_summary": {"fixture_audit_status": "failed"},
                "blockers": [{"code": "real_world_fixture_audit_failed", "message": "fixture missing"}],
                "outputs": {"workflow_json": "report.json"},
            }
        )

        self.assertIn("Status: `blocked`", markdown)
        self.assertIn("real_world_fixture_audit_failed", markdown)
        self.assertIn("workflow_json", markdown)


if __name__ == "__main__":
    unittest.main()
