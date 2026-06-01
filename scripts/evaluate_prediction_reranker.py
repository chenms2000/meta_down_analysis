"""Evaluate saved prediction rerankers on exported candidate examples."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import joblib
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from sklearn.metrics import average_precision_score, brier_score_loss, mean_absolute_error, roc_auc_score

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from train_prediction_reranker import (  # noqa: E402
    expected_calibration_error,
    ranking_metrics,
    resolved_label,
    score_frame,
    slugify,
)


DEFAULT_OUTPUT_ROOT = "learning_runs"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if hasattr(value, "item"):
        return value.item()
    return value


def write_parquet(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), path)


def evaluate_task(scored: pd.DataFrame) -> dict[str, Any]:
    labels = resolved_label(scored)
    label_score = (labels / 2.0).clip(0.0, 1.0)
    binary = (labels >= 2.0).astype(int)
    predicted = pd.to_numeric(scored["predicted_calibrated_confidence"], errors="coerce").fillna(0.0)
    metrics: dict[str, Any] = {
        "task": str(scored["prediction_task"].iloc[0]) if not scored.empty else "",
        "rows": int(len(scored)),
        "label_counts": json_safe(labels.value_counts().sort_index().to_dict()),
        "mae_label_score": float(mean_absolute_error(label_score, predicted)) if len(scored) else None,
        "brier_top_label": float(brier_score_loss(binary, predicted)) if len(set(binary.tolist())) > 1 else None,
        "ece_top_label": expected_calibration_error(binary.to_numpy(), predicted.to_numpy()),
        "overclaim_rate": float(((scored["predicted_confidence_tier"].isin(["high", "medium"])) & (labels <= 0)).mean()) if len(scored) else None,
    }
    if "context_mismatch" in scored.columns and len(scored):
        mismatch = pd.to_numeric(scored["context_mismatch"], errors="coerce").fillna(0).astype(int)
        metrics["context_mismatch_count"] = int(mismatch.sum())
        metrics["context_mismatch_overclaim_rate"] = float(((mismatch == 1) & scored["predicted_confidence_tier"].isin(["high", "medium"])).mean())
    if "appendix" in scored.columns and len(scored):
        appendix = pd.to_numeric(scored["appendix"], errors="coerce").fillna(0).astype(int)
        metrics["appendix_count"] = int(appendix.sum())
        metrics["appendix_overclaim_rate"] = float(((appendix == 1) & scored["predicted_confidence_tier"].isin(["high", "medium"])).mean())
    if len(set(binary.tolist())) > 1:
        metrics["roc_auc_top_label"] = float(roc_auc_score(binary, predicted))
        metrics["average_precision_top_label"] = float(average_precision_score(binary, predicted))
    metrics.update(ranking_metrics(scored))
    return metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate calibrated prediction rerankers.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--training-dir", default="")
    parser.add_argument("--model-dir", default="")
    parser.add_argument("--output-report", default="")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    run_dir = workspace / args.output_root / args.run_id
    training_dir = Path(args.training_dir).resolve() if args.training_dir else run_dir / "prediction_training"
    model_dir = Path(args.model_dir).resolve() if args.model_dir else run_dir / "models" / "prediction_reranker"
    report_path = Path(args.output_report).resolve() if args.output_report else run_dir / "reports" / "prediction_reranker_eval_report.json"

    examples = pq.read_table(training_dir / "training_examples.parquet").to_pandas()
    if examples.empty:
        raise RuntimeError(f"No training examples found in {training_dir / 'training_examples.parquet'}")
    scored_parts = []
    task_reports = []
    missing_tasks = []
    for task, task_frame in examples.groupby("prediction_task", dropna=False):
        task_name = str(task or "unknown_task")
        model_path = model_dir / slugify(task_name) / "model.joblib"
        if not model_path.exists():
            missing_tasks.append(task_name)
            continue
        bundle = joblib.load(model_path)
        scored = score_frame(
            task_frame.copy(),
            bundle["model"],
            bundle.get("calibrator"),
            bundle.get("numeric_feature_columns", []),
            bundle.get("categorical_feature_columns", []),
            bundle.get("feature_columns", []),
        )
        scored_parts.append(scored)
        task_report = evaluate_task(scored)
        task_report["model"] = str(model_path)
        task_report["calibration_status"] = bundle.get("calibration_status", "")
        task_reports.append(task_report)

    scored_all = pd.concat(scored_parts, ignore_index=True) if scored_parts else pd.DataFrame()
    scored_path = training_dir / "evaluated_training_examples.parquet"
    write_parquet(scored_path, scored_all)
    report = {
        "created_at_utc": utc_now(),
        "run_id": args.run_id,
        "training_examples": str(training_dir / "training_examples.parquet"),
        "scored_examples": str(scored_path),
        "model_dir": str(model_dir),
        "task_count": len(task_reports),
        "missing_tasks": missing_tasks,
        "tasks": task_reports,
        "notes": [
            "Metrics are against the current weak or human utility labels in the examples table.",
            "Use a separate human-reviewed table before interpreting these as biological performance.",
        ],
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(json_safe(report), ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"run_id": args.run_id, "task_count": len(task_reports), "missing_tasks": missing_tasks, "report": str(report_path)}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
