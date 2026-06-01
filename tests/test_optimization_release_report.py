import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_optimization_release_report.py"
SPEC = importlib.util.spec_from_file_location("build_optimization_release_report", SCRIPT)
reporter = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = reporter
SPEC.loader.exec_module(reporter)


class OptimizationReleaseReportTests(unittest.TestCase):
    def test_builds_markdown_and_json_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            release = "mvp_fake"
            run_id = "run_fake"
            phase_dir = root / "validation_reports" / release
            gate_dir = root / "learning_runs" / run_id / "reports"
            phase_dir.mkdir(parents=True)
            gate_dir.mkdir(parents=True)
            (phase_dir / "phase15_validation_report.json").write_text(
                json.dumps(
                    {
                        "report_hash": "phase_hash",
                        "summary": {"gate_status": "passed", "passed": 2, "total": 2},
                        "cases": [
                            {
                                "fixture": "common_metabolite_gold_set.json",
                                "case_id": "glucose_hmdb",
                                "status": "passed",
                                "validations": [
                                    {"endpoint": "/resolve", "status": "matched"},
                                    {"endpoint": "/precheck/metabolites", "status": "matched"},
                                ],
                            },
                            {
                                "fixture": "lcms_ambiguous_set.json",
                                "case_id": "hexose_mass_only",
                                "status": "passed",
                                "validations": [
                                    {"endpoint": "/precheck/metabolites", "status": "ambiguous"},
                                ],
                            },
                        ],
                        "optimization_summary": {
                            "entity_resolution": {"status": "passed", "resolver_case_count": 1, "precheck_case_count": 2, "abstention_case_count": 1},
                            "evidence_precision": {"status": "passed", "checks": {"p_literature_formula_contract": True}},
                            "ranking_calibration": {"status": "passed", "calibration_failure_count": 0, "traceability_failure_count": 0},
                            "release_readiness": {"failed_case_count": 0, "blocked_case_count": 0},
                        },
                        "evidence_validation": {"status": "passed", "support_row_count": 1, "sentence_count": 1},
                        "literature_evidence": {"status": "available", "evidence_candidate_count": 3, "supported_existing_edge_count": 2},
                    }
                ),
                encoding="utf-8",
            )
            llm_path = root / "validation_reports" / "llm_safe_adapter_regression_report.json"
            llm_path.parent.mkdir(parents=True, exist_ok=True)
            llm_path.write_text(
                json.dumps(
                    {
                        "report_hash": "llm_hash",
                        "summary": {"gate_status": "passed", "local_passed": 1, "adversarial_probe_passed": 3},
                    }
                ),
                encoding="utf-8",
            )
            (gate_dir / "release_gate_validation_report.json").write_text(
                json.dumps({"gate_status": "passed", "summary": {"passed": 4, "total": 4}}),
                encoding="utf-8",
            )
            name_dir = phase_dir / "name_audit"
            name_dir.mkdir()
            (name_dir / "compound_name_triage_summary.json").write_text(
                json.dumps(
                    {
                        "triage_version": "test",
                        "total_risk_rows": 3,
                        "policy_counts": {"auto_abstain_name_only": 1, "requires_identifier_or_orthogonal_feature": 2},
                        "release_policy": {
                            "auto_abstain_name_only": "abstain",
                            "requires_identifier_or_orthogonal_feature": "require evidence",
                        },
                    }
                ),
                encoding="utf-8",
            )
            (name_dir / "manual_adjudication_queue_summary.json").write_text(
                json.dumps(
                    {
                        "queue_version": "test",
                        "selected_rows": 2,
                        "sample_per_bucket": 1,
                        "bucket_counts": {"highest_risk_name_only": 1, "lipid_identifier_required": 1},
                        "output_csv": str(name_dir / "manual_adjudication_queue.csv"),
                        "interpretation": "This queue is a sampled manual follow-up target, not a completed adjudication result.",
                    }
                ),
                encoding="utf-8",
            )

            exit_code = reporter.main(["--workspace", str(root), "--release-id", release, "--run-id", run_id])
            self.assertEqual(exit_code, 0)
            output_json = phase_dir / "optimization_release_summary.json"
            output_md = phase_dir / "optimization_release_summary.md"
            self.assertTrue(output_json.exists())
            self.assertTrue(output_md.exists())
            payload = json.loads(output_json.read_text(encoding="utf-8"))
            self.assertEqual(payload["phase15_report_hash"], "phase_hash")
            self.assertEqual(payload["llm_report_hash"], "llm_hash")
            self.assertEqual(payload["compound_name_triage_summary"]["total_risk_rows"], 3)
            self.assertEqual(payload["fixture_benchmarks"][0][0], "common_metabolite_gold_set.json")
            self.assertEqual(payload["manual_adjudication_queue"]["selected_rows"], 2)
            markdown = output_md.read_text(encoding="utf-8")
            self.assertIn("entity_resolution", markdown)
            self.assertIn("Fixture Benchmarks", markdown)
            self.assertIn("Manual Adjudication Queue", markdown)
            self.assertIn("lcms_ambiguous_set.json", markdown)
            self.assertIn("compound_name_triage", markdown)
            self.assertIn("phase15_report_hash", markdown)


if __name__ == "__main__":
    unittest.main()
