"""Run temporal hold-out validation for generalized learning runs.

The validation retrains a temporary priority model after removing candidates
backed by recent publications, then scores those held-out rows to estimate
time-forward transfer. It does not overwrite the baseline ranker.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import joblib
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import average_precision_score, mean_absolute_error, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from train_priority_ranker import stable_bucket, training_sample_weights  # noqa: E402


DEFAULT_OUTPUT_ROOT = "learning_runs"
POSITIVE_LABELS = {"strong_positive", "medium_positive", "weak_positive"}
TEMPORAL_NUMERIC_COLUMNS = [
    "subject_graph_log",
    "object_graph_log",
    "subject_mention_log",
    "object_mention_log",
    "subject_support_log",
    "object_support_log",
    "subject_overlay_log",
    "object_overlay_log",
    "graph_degree_sum_log",
    "graph_degree_min_log",
    "graph_degree_max_log",
    "graph_degree_balance",
    "canonical_edge_sum_log",
    "direct_graph_support_indicator",
    "entity_semantic_overlap",
    "metabolic_theme_overlap",
    "cell_context_overlap",
    "signature_entity_overlap",
    "subject_pathway_log",
    "object_pathway_log",
    "shared_pathway_count",
    "shared_pathway_jaccard",
    "metabolite_target_pathway_bridge",
    "endpoint_overlay_any",
    "endpoint_support_min_log",
    "source_independent_prior_score",
    "overlay_confidence",
    "overlay_score",
]
TEMPORAL_CATEGORICAL_COLUMNS = ["subject_type", "object_type", "predicate", "relation_family", "subject_object_family"]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_parquet(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pq.read_table(path).to_pandas()


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(number) or math.isinf(number):
        return default
    return number


def build_temporal_feature_matrix(frame: pd.DataFrame, feature_columns: list[str] | None = None) -> tuple[pd.DataFrame, list[str]]:
    numeric = frame.copy()
    for column in TEMPORAL_NUMERIC_COLUMNS:
        values = numeric[column] if column in numeric.columns else pd.Series(0.0, index=numeric.index)
        numeric[column] = pd.to_numeric(values, errors="coerce").fillna(0.0)
    pieces = [numeric[TEMPORAL_NUMERIC_COLUMNS]]
    categorical_columns = [column for column in TEMPORAL_CATEGORICAL_COLUMNS if column in numeric.columns]
    categoricals = pd.get_dummies(numeric[categorical_columns].fillna(""), prefix=categorical_columns)
    pieces.append(categoricals.astype(float))
    matrix = pd.concat(pieces, axis=1)
    if feature_columns is not None:
        for column in feature_columns:
            if column not in matrix.columns:
                matrix[column] = 0.0
        matrix = matrix[feature_columns]
        return matrix, feature_columns
    columns = matrix.columns.tolist()
    return matrix, columns


def train_temporal_model(labels: pd.DataFrame) -> tuple[HistGradientBoostingRegressor, list[str], dict[str, Any]]:
    trainable = labels[labels["weak_label"] != "abstain"].copy()
    trainable["label_score"] = pd.to_numeric(trainable["label_score"], errors="coerce").fillna(0.0)
    trainable["validation_bucket"] = trainable["candidate_uid"].astype(str).map(stable_bucket)
    train_frame = trainable[trainable["validation_bucket"] != 0].copy()
    valid_frame = trainable[trainable["validation_bucket"] == 0].copy()
    if train_frame.empty or valid_frame.empty:
        train_frame = trainable.sample(frac=0.8, random_state=13) if len(trainable) > 5 else trainable
        valid_frame = trainable.drop(train_frame.index)
        if valid_frame.empty:
            valid_frame = train_frame

    x_train, feature_columns = build_temporal_feature_matrix(train_frame)
    x_valid, _ = build_temporal_feature_matrix(valid_frame, feature_columns)
    y_train = train_frame["label_score"].astype(float).to_numpy()
    y_valid = valid_frame["label_score"].astype(float).to_numpy()
    weights, weight_diagnostics = training_sample_weights(train_frame, source_balanced=True)

    model = HistGradientBoostingRegressor(
        max_iter=180,
        learning_rate=0.05,
        max_leaf_nodes=31,
        l2_regularization=0.02,
        random_state=13,
    )
    model.fit(x_train, y_train, sample_weight=weights)
    predicted = np.clip(model.predict(x_valid), 0.0, 1.0)
    metrics: dict[str, Any] = {
        "train_rows": int(len(train_frame)),
        "validation_rows": int(len(valid_frame)),
        "mae": float(mean_absolute_error(y_valid, predicted)),
        "mean_prediction_positive": float(np.mean(predicted[y_valid >= 0.5])) if np.any(y_valid >= 0.5) else None,
        "mean_prediction_negative": float(np.mean(predicted[y_valid < 0.5])) if np.any(y_valid < 0.5) else None,
        "feature_mode": "temporal_no_direct_literature_or_label_reason",
        "source_balanced": True,
        "feature_count": int(len(feature_columns)),
        "weight_diagnostics": weight_diagnostics,
    }
    binary = (y_valid >= 0.5).astype(int)
    if len(set(binary.tolist())) > 1:
        metrics["roc_auc"] = float(roc_auc_score(binary, predicted))
        metrics["average_precision"] = float(average_precision_score(binary, predicted))
    return model, feature_columns, metrics


def score_with_model(labels: pd.DataFrame, model: Any, feature_columns: list[str]) -> pd.DataFrame:
    matrix, _ = build_temporal_feature_matrix(labels, feature_columns)
    scored = labels.copy()
    scored["temporal_transfer_score"] = np.clip(model.predict(matrix), 0.0, 1.0).round(6)
    return scored.sort_values(["temporal_transfer_score", "candidate_uid"], ascending=[False, True]).reset_index(drop=True)


def precision_at_k(frame: pd.DataFrame, k: int) -> float | None:
    if frame.empty:
        return None
    top = frame.head(min(k, len(frame)))
    if top.empty:
        return None
    return float(top["heldout_positive"].mean())


def lift_at_k(frame: pd.DataFrame, k: int, prevalence: float | None) -> float | None:
    precision = precision_at_k(frame, k)
    if precision is None or prevalence is None or prevalence <= 0:
        return None
    return float(precision / prevalence)


def rank_fraction(ranks: pd.Series, denominator: int, cutoff_fraction: float) -> float | None:
    if denominator <= 0 or ranks.empty:
        return None
    return float((ranks <= max(1, int(denominator * cutoff_fraction))).mean())


def validation_metrics(scored_all: pd.DataFrame, heldout_scored: pd.DataFrame) -> dict[str, Any]:
    positives = heldout_scored[heldout_scored["heldout_positive"]].copy()
    labels = heldout_scored["heldout_positive"].astype(int).to_numpy()
    scores = pd.to_numeric(heldout_scored["temporal_transfer_score"], errors="coerce").fillna(0.0).to_numpy()
    prevalence = float(labels.mean()) if len(labels) else None
    metrics: dict[str, Any] = {
        "heldout_rows": int(len(heldout_scored)),
        "heldout_positive_rows": int(labels.sum()),
        "heldout_nonpositive_rows": int(len(labels) - labels.sum()),
        "heldout_positive_prevalence": prevalence,
        "heldout_unique_transfer_scores": int(pd.Series(scores).nunique()) if len(scores) else 0,
        "heldout_mean_transfer_score": float(np.mean(scores)) if len(scores) else None,
        "heldout_positive_mean_transfer_score": float(positives["temporal_transfer_score"].mean()) if not positives.empty else None,
        "heldout_positive_median_rank": float(positives["global_rank"].median()) if not positives.empty else None,
        "heldout_positive_top_1pct_fraction": rank_fraction(positives["global_rank"], len(scored_all), 0.01),
        "heldout_positive_top_5pct_fraction": rank_fraction(positives["global_rank"], len(scored_all), 0.05),
        "heldout_positive_top_10pct_fraction": rank_fraction(positives["global_rank"], len(scored_all), 0.10),
        "heldout_precision_at_50": precision_at_k(heldout_scored, 50),
        "heldout_precision_at_100": precision_at_k(heldout_scored, 100),
        "heldout_precision_at_500": precision_at_k(heldout_scored, 500),
        "heldout_lift_at_50": lift_at_k(heldout_scored, 50, prevalence),
        "heldout_lift_at_100": lift_at_k(heldout_scored, 100, prevalence),
        "heldout_lift_at_500": lift_at_k(heldout_scored, 500, prevalence),
    }
    if len(set(labels.tolist())) > 1:
        metrics["heldout_roc_auc"] = float(roc_auc_score(labels, scores))
        metrics["heldout_average_precision"] = float(average_precision_score(labels, scores))
    return metrics


def top_candidate_table(scored_all: pd.DataFrame, heldout_mask: pd.Series, top_n: int) -> pd.DataFrame:
    frame = scored_all.loc[heldout_mask.reindex(scored_all.index, fill_value=False)].copy()
    if frame.empty:
        return frame
    columns = [
        "global_rank",
        "candidate_uid",
        "publication_year",
        "subject_name",
        "predicate",
        "object_name",
        "source_kind",
        "weak_label",
        "label_score",
        "temporal_transfer_score",
        "label_reason",
        "paper_id",
        "source_ref",
    ]
    return frame[[column for column in columns if column in frame.columns]].head(top_n)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run temporal hold-out validation for a learning run.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--min-train-rows", type=int, default=100)
    parser.add_argument("--top-n", type=int, default=500)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    run_dir = workspace / args.output_root / args.run_id
    label_dir = run_dir / "weak_labels"
    report_dir = run_dir / "reports" / "temporal_holdout"

    labels = read_parquet(label_dir / "weak_labels.parquet")
    if labels.empty:
        raise FileNotFoundError(f"No weak labels found for run: {args.run_id}")
    if "temporal_holdout" not in labels.columns:
        raise ValueError("Weak labels do not include temporal_holdout.")

    heldout = labels["temporal_holdout"].fillna(False).astype(bool)
    heldout_positive = heldout & labels["weak_label"].isin(POSITIVE_LABELS)
    train_labels = labels.loc[~heldout & (labels["weak_label"] != "abstain")].copy()
    heldout_labels = labels.loc[heldout].copy()
    if int(heldout_positive.sum()) == 0:
        raise ValueError("Temporal hold-out has no positive labels.")
    if len(train_labels) < args.min_train_rows:
        raise ValueError(f"Training set too small after temporal hold-out: {len(train_labels)} rows")

    model, feature_columns, train_metrics = train_temporal_model(train_labels)
    scored_all = score_with_model(labels, model, feature_columns)
    scored_all["global_rank"] = np.arange(1, len(scored_all) + 1)
    scored_all["heldout_temporal"] = scored_all["temporal_holdout"].fillna(False).astype(bool)
    scored_all["heldout_positive"] = scored_all["heldout_temporal"] & scored_all["weak_label"].isin(POSITIVE_LABELS)
    heldout_scored = scored_all[scored_all["heldout_temporal"]].copy()
    heldout_scored = heldout_scored.sort_values(["temporal_transfer_score", "candidate_uid"], ascending=[False, True])
    metrics = validation_metrics(scored_all, heldout_scored)

    report_dir.mkdir(parents=True, exist_ok=True)
    top_candidates = top_candidate_table(scored_all, scored_all["heldout_temporal"], args.top_n)
    top_candidates.to_csv(report_dir / "temporal_top_candidates.csv", index=False, encoding="utf-8-sig")
    heldout_scored.to_csv(report_dir / "temporal_all_heldout_scored.csv", index=False, encoding="utf-8-sig")
    joblib.dump(
        {
            "model": model,
            "feature_columns": feature_columns,
            "heldout": {"field": "temporal_holdout", "value": True},
            "feature_mode": "temporal_no_direct_literature_or_label_reason",
            "source_balanced": True,
        },
        report_dir / "temporal_temporary_model.joblib",
    )

    year_counts = {}
    if "publication_year" in labels.columns:
        year_counts = labels.loc[heldout, "publication_year"].fillna("").astype(str).value_counts().sort_index().to_dict()

    summary = {
        "created_at_utc": utc_now(),
        "run_id": args.run_id,
        "scope": "temporal_holdout_validation",
        "heldout": {"field": "temporal_holdout", "value": True},
        "input_rows": {
            "weak_labels": int(len(labels)),
            "train_rows": int(len(train_labels)),
            "heldout_rows": int(len(heldout_labels)),
            "heldout_positive_rows": int(heldout_positive.sum()),
        },
        "heldout_year_counts": year_counts,
        "training_metrics": train_metrics,
        "validation_metrics": metrics,
        "outputs": {
            "top_candidates": str(report_dir / "temporal_top_candidates.csv"),
            "all_heldout_scored": str(report_dir / "temporal_all_heldout_scored.csv"),
            "temporary_model": str(report_dir / "temporal_temporary_model.joblib"),
            "summary": str(report_dir / "temporal_holdout_validation_report.json"),
        },
        "notes": [
            "Temporary temporal models are validation artifacts and do not replace the frozen baseline ranker.",
            "Held-out labels are weak supervision rows backed by publications mapped to 2024 or later.",
            "The temporal feature mode excludes direct literature evidence counts, p_literature, source_kind, and label_reason to reduce future-evidence leakage, with source-balanced weights during temporary training.",
            "No canonical graph, normalized store, baseline model, or baseline rankings are modified.",
        ],
    }
    (report_dir / "temporal_holdout_validation_report.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
