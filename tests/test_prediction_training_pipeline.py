import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover
    pq = None


ROOT = Path(__file__).resolve().parents[1]


def prediction_row(index: int, confidence: float, tier: str, label_support: int = 1):
    return {
        "rank": index,
        "pathway_uid": f"pathway_{index}",
        "display_name": f"Pathway {index}",
        "prediction_id": f"pathway_{index}:rank_{index}",
        "prediction": f"The model predicts Pathway {index}.",
        "prediction_task": "pathway_prediction",
        "result_type": "biochemical_process",
        "confidence_tier": tier,
        "calibrated_confidence": confidence,
        "calibration_status": "calibrated_in_scope" if tier in {"high", "medium"} else "exploratory_in_scope",
        "input_support_count": label_support,
        "matched_input_count": 3,
        "ambiguous_input_count": 0,
        "evidence_sources": ["input", "database_pathway"],
        "graph_distance": 1,
        "direction_consistency": 1.0,
        "is_directly_supported": True,
        "is_extrapolated": False,
        "needs_validation": tier != "high",
        "boundary": "Fixture boundary.",
        "ranker_score": confidence,
        "score": confidence,
        "literature_support_count": 1,
        "max_p_literature": 0.8,
        "claim_refs": {
            "traceability_passed": True,
            "edge_uids": [f"edge_{index}"],
            "evidence_ref_uids": [f"support_{index}"],
        },
        "calibration_components": {
            "raw_score": confidence,
            "node_degree": index,
        },
    }


@unittest.skipIf(pq is None, "pyarrow is required")
class PredictionTrainingPipelineTests(unittest.TestCase):
    def test_build_train_and_evaluate_from_analysis_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            analysis_json = tmp_path / "analysis_pack.json"
            rows = [
                prediction_row(1, 0.82, "high"),
                prediction_row(2, 0.66, "medium"),
                prediction_row(3, 0.52, "medium"),
                prediction_row(4, 0.32, "exploratory"),
                prediction_row(5, 0.2, "exploratory"),
                prediction_row(6, 0.08, "low", label_support=0),
            ]
            pack = {
                "analysis_pack": {
                    "contract_version": "analysis_pack.v1",
                    "input_summary": {"input_count": 3, "matched_count": 3, "ambiguous_count": 0, "unmatched_count": 0},
                    "pathway_rankings": rows,
                    "target_rankings": [],
                    "disease_rankings": [],
                    "prediction_model": {
                        "contract_version": "metabolic_prediction_model.v1",
                        "model_name": "fixture_rule_mvp",
                        "predictions": rows,
                    },
                    "release": {"release_id": "fixture_release"},
                    "determinism": {"analysis_pack_hash": "fixture_hash"},
                }
            }
            analysis_json.write_text(json.dumps(pack), encoding="utf-8")
            training_dir = tmp_path / "learning_runs" / "fixture_run" / "prediction_training"

            commands = [
                [
                    sys.executable,
                    str(ROOT / "scripts" / "build_prediction_training_set.py"),
                    "--workspace",
                    str(tmp_path),
                    "--run-id",
                    "fixture_run",
                    "--analysis-json",
                    str(analysis_json),
                    "--output-dir",
                    str(training_dir),
                    "--no-default-fixtures",
                ],
                [
                    sys.executable,
                    str(ROOT / "scripts" / "train_prediction_reranker.py"),
                    "--workspace",
                    str(tmp_path),
                    "--run-id",
                    "fixture_run",
                    "--training-dir",
                    str(training_dir),
                ],
                [
                    sys.executable,
                    str(ROOT / "scripts" / "evaluate_prediction_reranker.py"),
                    "--workspace",
                    str(tmp_path),
                    "--run-id",
                    "fixture_run",
                    "--training-dir",
                    str(training_dir),
                ],
            ]
            for command in commands:
                subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=True)

            examples = pq.read_table(training_dir / "training_examples.parquet").to_pandas()
            self.assertEqual(len(examples), 6)
            self.assertIn("weak_label", examples.columns)
            self.assertTrue((training_dir / "scored_training_examples.parquet").exists())
            report = json.loads((tmp_path / "learning_runs" / "fixture_run" / "reports" / "prediction_reranker_eval_report.json").read_text())
            self.assertEqual(report["task_count"], 1)
            self.assertEqual(report["missing_tasks"], [])


if __name__ == "__main__":
    unittest.main()
