"""Train task-specific prediction rerankers and calibration layers."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import average_precision_score, brier_score_loss, mean_absolute_error, roc_auc_score


DEFAULT_OUTPUT_ROOT = "learning_runs"
TARGET_COLUMNS = {"label", "weak_label", "human_label"}
TIER_THRESHOLDS = (("high", 0.75), ("medium", 0.45), ("exploratory", 0.18), ("low", 0.0))
SLUG_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def slugify(value: str) -> str:
    slug = SLUG_RE.sub("_", value.strip())[:80].strip("_")
    return slug or "unknown_task"


def stable_bucket(value: str, buckets: int = 10) -> int:
    return int(hashlib.sha256(value.encode("utf-8")).hexdigest()[:8], 16) % buckets


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(number) or math.isinf(number):
        return default
    return number


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [json_safe(item) for item in value]
    if hasattr(value, "item"):
        return value.item()
    return value


def write_parquet(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), path)


def confidence_tier(score: float) -> str:
    for tier, threshold in TIER_THRESHOLDS:
        if score >= threshold:
            return tier
    return "low"


def build_feature_matrix(
    frame: pd.DataFrame,
    numeric_columns: list[str],
    categorical_columns: list[str],
    feature_columns: list[str] | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    numeric_pieces = []
    for column in numeric_columns:
        if column in TARGET_COLUMNS or column not in frame.columns:
            continue
        numeric_pieces.append(pd.to_numeric(frame[column], errors="coerce").fillna(0.0).rename(column))
    if numeric_pieces:
        matrix = pd.concat(numeric_pieces, axis=1)
    else:
        matrix = pd.DataFrame(index=frame.index)
    present_categoricals = [column for column in categorical_columns if column in frame.columns]
    if present_categoricals:
        dummies = pd.get_dummies(frame[present_categoricals].fillna("").astype(str), prefix=present_categoricals)
        matrix = pd.concat([matrix, dummies.astype(float)], axis=1)
    if feature_columns is not None:
        for column in feature_columns:
            if column not in matrix.columns:
                matrix[column] = 0.0
        return matrix[feature_columns], feature_columns
    columns = matrix.columns.tolist()
    return matrix, columns


def resolved_label(frame: pd.DataFrame) -> pd.Series:
    labels = pd.to_numeric(frame.get("label", 0), errors="coerce")
    if "human_label" in frame.columns:
        human = pd.to_numeric(frame["human_label"], errors="coerce")
        if human.notna().any():
            labels = human.combine_first(labels)
    return labels.fillna(0.0).clip(0.0, 2.0)


def ensure_split(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    if "split" not in result.columns or result["split"].fillna("").eq("").all():
        result["split"] = [
            "validation" if stable_bucket(f"{row.analysis_id}:{row.prediction_task}:{row.candidate_id}") == 1 else "train"
            for row in result.itertuples(index=False)
        ]
    if (result["split"] == "train").sum() == 0:
        result.loc[result.index, "split"] = "train"
    if (result["split"] == "validation").sum() == 0 and len(result) > 4:
        validation_index = result.index[::5]
        result.loc[validation_index, "split"] = "validation"
        result.loc[result.index.difference(validation_index), "split"] = "train"
    return result


def calibrate(raw_valid: np.ndarray, y_binary: np.ndarray) -> tuple[Any | None, str]:
    if len(raw_valid) >= 6 and len(set(y_binary.astype(int).tolist())) > 1:
        calibrator = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
        calibrator.fit(raw_valid, y_binary)
        return calibrator, "isotonic_validation"
    return None, "insufficient_validation_classes"


def apply_calibrator(raw: np.ndarray, calibrator: Any | None) -> np.ndarray:
    raw = np.clip(raw, 0.0, 1.0)
    if calibrator is None:
        return raw
    return np.clip(calibrator.predict(raw), 0.0, 1.0)


def expected_calibration_error(y_true: np.ndarray, y_score: np.ndarray, bins: int = 10) -> float | None:
    if len(y_true) == 0:
        return None
    ece = 0.0
    for lower in np.linspace(0.0, 1.0, bins, endpoint=False):
        upper = lower + 1.0 / bins
        mask = (y_score >= lower) & (y_score < upper if upper < 1.0 else y_score <= upper)
        if not np.any(mask):
            continue
        ece += float(np.mean(mask)) * abs(float(np.mean(y_score[mask])) - float(np.mean(y_true[mask])))
    return float(ece)


def ranking_metrics(frame: pd.DataFrame, score_column: str = "predicted_calibrated_confidence", k: int = 5) -> dict[str, Any]:
    if frame.empty:
        return {"precision_at_5": None, "ndcg_at_10": None, "mrr": None, "group_count": 0}
    precision_values = []
    ndcg_values = []
    mrr_values = []
    for _, group in frame.groupby("analysis_id", dropna=False):
        ranked = group.sort_values(score_column, ascending=False).reset_index(drop=True)
        labels = resolved_label(ranked).to_numpy()
        binary = (labels >= 2.0).astype(float)
        top = binary[: min(k, len(binary))]
        if len(top):
            precision_values.append(float(np.mean(top)))
        gains = (np.power(2.0, labels[:10]) - 1.0) / np.log2(np.arange(2, min(10, len(labels)) + 2))
        ideal_labels = np.sort(labels)[::-1]
        ideal = (np.power(2.0, ideal_labels[:10]) - 1.0) / np.log2(np.arange(2, min(10, len(ideal_labels)) + 2))
        ideal_sum = float(np.sum(ideal))
        ndcg_values.append(float(np.sum(gains) / ideal_sum) if ideal_sum > 0 else 0.0)
        positive_positions = np.where(binary > 0)[0]
        mrr_values.append(float(1.0 / (positive_positions[0] + 1)) if len(positive_positions) else 0.0)
    return {
        "precision_at_5": float(np.mean(precision_values)) if precision_values else None,
        "ndcg_at_10": float(np.mean(ndcg_values)) if ndcg_values else None,
        "mrr": float(np.mean(mrr_values)) if mrr_values else None,
        "group_count": int(frame["analysis_id"].nunique(dropna=False)),
    }


def score_frame(
    frame: pd.DataFrame,
    model: Any,
    calibrator: Any | None,
    numeric_columns: list[str],
    categorical_columns: list[str],
    feature_columns: list[str],
) -> pd.DataFrame:
    matrix, _ = build_feature_matrix(frame, numeric_columns, categorical_columns, feature_columns)
    result = frame.copy()
    raw = np.clip(model.predict(matrix), 0.0, 1.0)
    calibrated = apply_calibrator(raw, calibrator)
    result["predicted_raw_score"] = np.round(raw, 6)
    result["predicted_calibrated_confidence"] = np.round(calibrated, 6)
    result["predicted_confidence_tier"] = [confidence_tier(float(score)) for score in calibrated]
    return result


def train_task(task_frame: pd.DataFrame, schema: dict[str, Any]) -> tuple[dict[str, Any], pd.DataFrame]:
    task_frame = ensure_split(task_frame)
    labels = resolved_label(task_frame)
    task_frame = task_frame.copy()
    task_frame["label_score"] = (labels / 2.0).clip(0.0, 1.0)
    task_frame["label_binary_top"] = (labels >= 2.0).astype(int)
    numeric_columns = [column for column in schema.get("numeric_feature_columns", []) if column not in TARGET_COLUMNS]
    categorical_columns = schema.get("categorical_feature_columns", [])
    train_frame = task_frame[task_frame["split"] == "train"].copy()
    valid_frame = task_frame[task_frame["split"] == "validation"].copy()
    if train_frame.empty:
        train_frame = task_frame.copy()
    if valid_frame.empty:
        valid_frame = task_frame.copy()

    x_train, feature_columns = build_feature_matrix(train_frame, numeric_columns, categorical_columns)
    y_train = train_frame["label_score"].astype(float).to_numpy()
    if len(train_frame) < 6 or len(set(np.round(y_train, 6).tolist())) <= 1 or x_train.shape[1] == 0:
        model: Any = DummyRegressor(strategy="mean")
        model_family = "dummy_mean_regressor"
    else:
        model = HistGradientBoostingRegressor(
            max_iter=180,
            learning_rate=0.05,
            max_leaf_nodes=31,
            l2_regularization=0.02,
            random_state=13,
        )
        model_family = "hist_gradient_boosting_regressor"
    model.fit(x_train, y_train)

    x_valid, _ = build_feature_matrix(valid_frame, numeric_columns, categorical_columns, feature_columns)
    valid_raw = np.clip(model.predict(x_valid), 0.0, 1.0)
    y_valid_binary = valid_frame["label_binary_top"].astype(int).to_numpy()
    calibrator, calibration_status = calibrate(valid_raw, y_valid_binary)
    scored = score_frame(task_frame, model, calibrator, numeric_columns, categorical_columns, feature_columns)
    valid_scored = scored[scored["split"] == "validation"].copy()
    if valid_scored.empty:
        valid_scored = scored
    y_valid_score = valid_scored["label_score"].astype(float).to_numpy()
    y_binary = valid_scored["label_binary_top"].astype(int).to_numpy()
    predicted = valid_scored["predicted_calibrated_confidence"].astype(float).to_numpy()

    metrics: dict[str, Any] = {
        "task": str(task_frame["prediction_task"].iloc[0]),
        "model_family": model_family,
        "calibration_status": calibration_status,
        "rows": int(len(task_frame)),
        "train_rows": int(len(train_frame)),
        "validation_rows": int(len(valid_scored)),
        "label_counts": json_safe(task_frame["label"].value_counts().sort_index().to_dict()),
        "mae_label_score": float(mean_absolute_error(y_valid_score, predicted)) if len(predicted) else None,
        "brier_top_label": float(brier_score_loss(y_binary, predicted)) if len(set(y_binary.tolist())) > 1 else None,
        "ece_top_label": expected_calibration_error(y_binary, predicted),
        "overclaim_rate": float(((valid_scored["predicted_confidence_tier"].isin(["high", "medium"])) & (valid_scored["label"].astype(float) <= 0)).mean())
        if len(valid_scored)
        else None,
    }
    if len(set(y_binary.tolist())) > 1:
        metrics["roc_auc_top_label"] = float(roc_auc_score(y_binary, predicted))
        metrics["average_precision_top_label"] = float(average_precision_score(y_binary, predicted))
    metrics.update(ranking_metrics(valid_scored))
    bundle = {
        "model": model,
        "calibrator": calibrator,
        "task": metrics["task"],
        "model_family": model_family,
        "calibration_status": calibration_status,
        "feature_columns": feature_columns,
        "numeric_feature_columns": numeric_columns,
        "categorical_feature_columns": categorical_columns,
        "tier_thresholds": TIER_THRESHOLDS,
        "metrics": metrics,
    }
    return bundle, scored


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train calibrated prediction rerankers from exported training examples.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--training-dir", default="")
    parser.add_argument("--min-task-rows", type=int, default=1)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    run_dir = workspace / args.output_root / args.run_id
    training_dir = Path(args.training_dir).resolve() if args.training_dir else run_dir / "prediction_training"
    model_dir = run_dir / "models" / "prediction_reranker"
    report_dir = run_dir / "reports"

    examples = pq.read_table(training_dir / "training_examples.parquet").to_pandas()
    schema = read_json(training_dir / "feature_schema.json")
    if examples.empty:
        raise RuntimeError(f"No training examples found in {training_dir / 'training_examples.parquet'}")

    scored_parts = []
    tasks = []
    for task, task_frame in examples.groupby("prediction_task", dropna=False):
        task_name = str(task or "unknown_task")
        if len(task_frame) < args.min_task_rows:
            continue
        bundle, scored = train_task(task_frame.copy(), schema)
        task_slug = slugify(task_name)
        task_model_dir = model_dir / task_slug
        task_model_dir.mkdir(parents=True, exist_ok=True)
        joblib.dump(bundle, task_model_dir / "model.joblib")
        feature_schema = {
            "task": task_name,
            "numeric_feature_columns": bundle["numeric_feature_columns"],
            "categorical_feature_columns": bundle["categorical_feature_columns"],
            "expanded_feature_columns": bundle["feature_columns"],
            "tier_thresholds": TIER_THRESHOLDS,
        }
        (task_model_dir / "feature_schema.json").write_text(json.dumps(feature_schema, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        tasks.append({"task": task_name, "slug": task_slug, "metrics": bundle["metrics"], "model": str(task_model_dir / "model.joblib")})
        scored_parts.append(scored)

    scored_all = pd.concat(scored_parts, ignore_index=True) if scored_parts else pd.DataFrame()
    write_parquet(training_dir / "scored_training_examples.parquet", scored_all)
    report = {
        "created_at_utc": utc_now(),
        "run_id": args.run_id,
        "training_examples": str(training_dir / "training_examples.parquet"),
        "scored_training_examples": str(training_dir / "scored_training_examples.parquet"),
        "task_count": len(tasks),
        "tasks": tasks,
        "notes": [
            "Reranker confidence is calibrated against weak utility labels, not biological truth.",
            "Tasks with too little label diversity use a constant baseline until more reviewed examples exist.",
        ],
    }
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "prediction_reranker_training_report.json").write_text(
        json.dumps(json_safe(report), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"run_id": args.run_id, "task_count": len(tasks), "scored_rows": len(scored_all)}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
