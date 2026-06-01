import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_release_gate_validation.py"
SPEC = importlib.util.spec_from_file_location("run_release_gate_validation", SCRIPT)
release_gate = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = release_gate
SPEC.loader.exec_module(release_gate)


class ReleaseGateValidationTests(unittest.TestCase):
    def test_structured_prediction_gate_enforces_drug_and_generalized_bounds(self):
        structured = {
            "mode": "generalized",
            "high_confidence_themes": [
                {
                    "display_name": "glycolysis",
                    "confidence_tier": "high",
                    "matched_support_count": 1,
                    "related_input_count": 1,
                    "evidence_refs": [{"support_uid": "support1"}],
                    "claim_refs": {"traceability_passed": True},
                    "appendix": False,
                    "research_only": False,
                }
            ],
            "medium_confidence_themes": [],
            "exploratory_themes": [],
            "downgraded_but_supported_themes": [],
            "pathway_evidence": [
                {
                    "pathway_name": "Melanoma",
                    "confidence_tier": "low",
                    "evidence_refs": [{"support_uid": "support2"}],
                    "claim_refs": {"traceability_passed": True},
                    "appendix": True,
                    "research_only": True,
                }
            ],
            "disease_predictions": [
                {
                    "display_name": "Colon cancer",
                    "confidence_tier": "low",
                    "evidence_refs": [{"support_uid": "support3"}],
                    "claim_refs": {"traceability_passed": True},
                    "appendix": True,
                    "research_only": True,
                }
            ],
            "drug_hypotheses": [
                {
                    "drug_name": "Fixture drug",
                    "confidence_tier": "low",
                    "upstream_target_confidence": "low",
                    "evidence_refs": [{"support_uid": "support4"}],
                    "claim_refs": {"traceability_passed": True},
                    "research_only": True,
                    "appendix": True,
                    "excluded_from_primary_reason": "clinical_warning_or_research_only",
                    "clinical_warning": "Research-only hypothesis; not a treatment recommendation.",
                }
            ],
            "context_mismatch_results": [{"display_name": "wrong context", "calibration_status": "context_mismatch", "confidence_tier": "exploratory"}],
            "product_layers": {
                "primary_research_candidates": [{"display_name": "glycolysis", "result_type": "metabolic_theme", "confidence_tier": "high"}],
                "mechanistic_support": [],
                "appendix_overlay_only": [],
                "context_mismatch": [],
                "clinical_warning_or_research_only": [{"display_name": "Fixture drug", "result_type": "drug", "research_only": True, "appendix": True}],
                "excluded_from_primary_reason": [{"display_name": "Fixture drug", "result_type": "drug", "excluded_from_primary_reason": "clinical_warning_or_research_only"}],
            },
        }

        checks = release_gate.check_structured_prediction(structured, "fixture")
        self.assertTrue(all(check["passed"] for check in checks), checks)

        structured["drug_hypotheses"][0]["confidence_tier"] = "medium"
        failed = release_gate.check_structured_prediction(structured, "fixture")
        failed_names = {check["name"] for check in failed if not check["passed"]}
        self.assertIn("generalized_mode_no_strong_disease_or_drug_claims", failed_names)
        self.assertIn("drug_confidence_not_above_target_confidence", failed_names)

    def test_run_artifact_gate_checks_overlay_and_heldout_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_dir = root / "learning_runs" / "gate_run"
            views = run_dir / "views"
            reports = run_dir / "reports"
            views.mkdir(parents=True)
            pq.write_table(
                pa.Table.from_pandas(
                    pd.DataFrame(
                        [
                            {
                                "cell_type_id": "ct1",
                                "cell_type_name": "Epithelial",
                                "source_name": "CELLxGENE",
                                "source_record_id": "ds1",
                                "evidence_level": "cell_type_signature",
                                "license_id": "public_research:cellxgene",
                            }
                        ]
                    ),
                    preserve_index=False,
                ),
                views / "cell_context_overlay.parquet",
            )
            (reports / "temporal_holdout").mkdir(parents=True)
            (reports / "source_heldout").mkdir(parents=True)
            (reports / "temporal_holdout" / "temporal_holdout_validation_report.json").write_text(
                json.dumps({"metrics": {"heldout_rows": 5, "heldout_roc_auc": 0.61}}),
                encoding="utf-8",
            )
            (reports / "source_heldout" / "source_heldout_validation_report.json").write_text(
                json.dumps({"heldout_sources": [{"source": {"field": "source_name", "value": "literature"}, "metrics": {"heldout_rows": 4}}]}),
                encoding="utf-8",
            )
            validation_reports = root / "validation_reports"
            validation_reports.mkdir()
            (validation_reports / "llm_safe_adapter_regression_report.json").write_text(
                json.dumps({"summary": {"gate_status": "passed", "local_failed": 0, "adversarial_probe_failed": 0}}),
                encoding="utf-8",
            )

            report = release_gate.build_report(root, "gate_run", [])
            self.assertEqual(report["gate_status"], "passed")


if __name__ == "__main__":
    unittest.main()
