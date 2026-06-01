"""Run source-held-out validation for generalized learning runs.

The validation retrains a temporary priority ranker after removing one evidence
source from weak supervision, then scores the held-out candidates to check
whether their relations still rank well from transferable graph/literature
features. It does not overwrite the baseline model or ranking outputs.
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
from train_priority_ranker import build_feature_matrix, stable_bucket, training_sample_weights  # noqa: E402


DEFAULT_OUTPUT_ROOT = "learning_runs"
LITERATURE_SOURCE_PARTITION_BUCKETS = 5
POSITIVE_LABELS = {"strong_positive", "medium_positive", "weak_positive"}
NEGATIVE_LABELS = {"negative_shuffled", "negative_refute"}


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


def normalize_source_name(value: Any) -> str:
    text = str(value or "").strip()
    return text if text else "unknown"


def source_kind_family_rows(labels: pd.DataFrame, heldout_mask: pd.Series, source_kinds: set[str]) -> int:
    if "source_kind" not in labels.columns:
        return 0
    trainable = labels.loc[~heldout_mask].copy()
    if "weak_label" in trainable.columns:
        trainable = trainable[trainable["weak_label"] != "abstain"]
    return int(trainable["source_kind"].astype(str).isin(source_kinds).sum())


def add_source_names(labels: pd.DataFrame, drug_overlay: pd.DataFrame, literature: pd.DataFrame) -> pd.DataFrame:
    result = labels.copy()
    result["source_name"] = result["source_kind"].astype(str)
    result["source_version"] = ""
    result["source_license_id"] = ""

    if not drug_overlay.empty and "source_record_id" in drug_overlay.columns:
        drug_lookup = drug_overlay[
            [
                column
                for column in ["source_record_id", "source_name", "source_version", "license_id"]
                if column in drug_overlay.columns
            ]
        ].copy()
        drug_lookup = drug_lookup.rename(
            columns={
                "source_record_id": "source_ref",
                "source_name": "overlay_source_name",
                "source_version": "overlay_source_version",
                "license_id": "overlay_license_id",
            }
        ).drop_duplicates("source_ref")
        result = result.merge(drug_lookup, on="source_ref", how="left")
        result.loc[result["overlay_source_name"].notna(), "source_name"] = result.loc[
            result["overlay_source_name"].notna(), "overlay_source_name"
        ]
        result.loc[result["overlay_source_version"].notna(), "source_version"] = result.loc[
            result["overlay_source_version"].notna(), "overlay_source_version"
        ]
        result.loc[result["overlay_license_id"].notna(), "source_license_id"] = result.loc[
            result["overlay_license_id"].notna(), "overlay_license_id"
        ]
        result = result.drop(columns=["overlay_source_name", "overlay_source_version", "overlay_license_id"], errors="ignore")

    if not literature.empty and "support_uid" in literature.columns:
        lit_lookup = literature[
            [
                column
                for column in ["support_uid", "source_release", "license_id", "config_hash"]
                if column in literature.columns
            ]
        ].copy()
        lit_lookup = lit_lookup.rename(
            columns={
                "support_uid": "source_ref",
                "source_release": "lit_source_release",
                "license_id": "lit_license_id",
                "config_hash": "lit_config_hash",
            }
        ).drop_duplicates("source_ref")
        result = result.merge(lit_lookup, on="source_ref", how="left")
        lit_mask = result["lit_license_id"].notna()
        result.loc[lit_mask, "source_name"] = "literature:" + result.loc[lit_mask, "lit_license_id"].astype(str)
        result.loc[result["lit_source_release"].notna(), "source_version"] = result.loc[
            result["lit_source_release"].notna(), "lit_source_release"
        ]
        result.loc[result["lit_license_id"].notna(), "source_license_id"] = result.loc[
            result["lit_license_id"].notna(), "lit_license_id"
        ]
        result = result.drop(columns=["lit_source_release", "lit_license_id", "lit_config_hash"], errors="ignore")

    result["source_name"] = result["source_name"].map(normalize_source_name)
    result["source_partition"] = result["source_name"]
    result["source_partition_strategy"] = "source_name"
    if {"source_kind", "source_name"}.issubset(result.columns):
        lit_mask = result["source_kind"].astype(str).eq("literature_relation")
        if lit_mask.any():
            seed = pd.Series("", index=result.index, dtype=object)
            for column in ("paper_id", "source_ref", "candidate_uid"):
                if column not in result.columns:
                    continue
                values = result[column].fillna("").astype(str)
                seed = seed.where(seed.astype(str).str.len() > 0, values)
            bucket = seed.loc[lit_mask].map(lambda value: stable_bucket(str(value), LITERATURE_SOURCE_PARTITION_BUCKETS))
            result.loc[lit_mask, "source_partition"] = (
                result.loc[lit_mask, "source_name"].astype(str) + ":partition_" + bucket.astype(str)
            )
            result.loc[lit_mask, "source_partition_strategy"] = "literature_stable_bucket"
    return result


def infer_source_specs(labels: pd.DataFrame, min_positive_rows: int) -> list[dict[str, str]]:
    specs: list[dict[str, str]] = []
    candidates = labels[labels["source_kind"] != "shuffled_negative"].copy()
    for source_name, group in candidates.groupby("source_name", dropna=False):
        positive_rows = int(group["weak_label"].isin(POSITIVE_LABELS).sum())
        if positive_rows >= min_positive_rows:
            source_value = normalize_source_name(source_name)
            heldout_mask = source_mask(labels, "source_name", source_value)
            source_kinds = {str(kind) for kind in group.get("source_kind", pd.Series(dtype=object)).dropna().unique()}
            if source_kind_family_rows(labels, heldout_mask, source_kinds) > 0:
                specs.append({"field": "source_name", "value": source_value})
                continue

            partition_specs: list[dict[str, str]] = []
            if "source_partition" in group.columns:
                for partition, partition_group in group.groupby("source_partition", dropna=False):
                    partition_value = normalize_source_name(partition)
                    partition_positive = int(partition_group["weak_label"].isin(POSITIVE_LABELS).sum())
                    if partition_positive < min_positive_rows:
                        continue
                    partition_mask = source_mask(labels, "source_partition", partition_value)
                    if source_kind_family_rows(labels, partition_mask, source_kinds) <= 0:
                        continue
                    partition_specs.append({"field": "source_partition", "value": partition_value})
            specs.extend(partition_specs or [{"field": "source_name", "value": source_value}])
    return specs


def source_mask(labels: pd.DataFrame, field: str, value: str) -> pd.Series:
    if field not in labels.columns:
        raise ValueError(f"Unknown held-out source field: {field}")
    return labels[field].astype(str).str.casefold() == value.casefold()


def train_transfer_model(labels: pd.DataFrame) -> tuple[HistGradientBoostingRegressor, list[str], dict[str, Any]]:
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

    x_train, feature_columns = build_feature_matrix(train_frame, feature_mode="no_leakage")
    x_valid, _ = build_feature_matrix(valid_frame, feature_columns, feature_mode="no_leakage")
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
        "feature_mode": "no_leakage",
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
    matrix, _ = build_feature_matrix(labels, feature_columns, feature_mode="no_leakage")
    scored = labels.copy()
    scored["transfer_score"] = np.clip(model.predict(matrix), 0.0, 1.0).round(6)
    return scored.sort_values(["transfer_score", "candidate_uid"], ascending=[False, True]).reset_index(drop=True)


def source_family_diagnostics(train_labels: pd.DataFrame, heldout_labels: pd.DataFrame) -> dict[str, Any]:
    train_counts = {
        str(key): int(value)
        for key, value in train_labels.get("source_kind", pd.Series(dtype=object)).value_counts(dropna=False).sort_index().items()
    }
    heldout_counts = {
        str(key): int(value)
        for key, value in heldout_labels.get("source_kind", pd.Series(dtype=object)).value_counts(dropna=False).sort_index().items()
    }
    heldout_kinds = sorted(heldout_counts)
    missing_kinds = [kind for kind in heldout_kinds if kind not in train_counts]
    same_family_train_rows = int(sum(train_counts.get(kind, 0) for kind in heldout_kinds))
    return {
        "train_source_kind_counts": train_counts,
        "heldout_source_kind_counts": heldout_counts,
        "heldout_source_kinds": heldout_kinds,
        "missing_heldout_source_kinds_in_train": missing_kinds,
        "same_family_train_rows": same_family_train_rows,
        "same_family_available": same_family_train_rows > 0 and not missing_kinds,
    }


def feature_transfer_diagnostics(
    train_labels: pd.DataFrame,
    heldout_labels: pd.DataFrame,
    feature_columns: list[str],
) -> dict[str, Any]:
    train_matrix, _ = build_feature_matrix(train_labels, feature_columns, feature_mode="no_leakage")
    heldout_matrix, _ = build_feature_matrix(heldout_labels, feature_columns, feature_mode="no_leakage")
    train_variance = train_matrix.var(numeric_only=True)
    heldout_variance = heldout_matrix.var(numeric_only=True)
    train_nonconstant = {str(column) for column, value in train_variance.items() if safe_float(value) > 0.0}
    heldout_nonconstant = {str(column) for column, value in heldout_variance.items() if safe_float(value) > 0.0}
    shared_nonconstant = sorted(train_nonconstant & heldout_nonconstant)
    heldout_only = sorted(heldout_nonconstant - train_nonconstant)
    train_only = sorted(train_nonconstant - heldout_nonconstant)
    return {
        "feature_count": int(len(feature_columns)),
        "train_nonconstant_feature_count": int(len(train_nonconstant)),
        "heldout_nonconstant_feature_count": int(len(heldout_nonconstant)),
        "shared_nonconstant_feature_count": int(len(shared_nonconstant)),
        "heldout_only_nonconstant_feature_count": int(len(heldout_only)),
        "train_only_nonconstant_feature_count": int(len(train_only)),
        "shared_nonconstant_features_sample": shared_nonconstant[:20],
        "heldout_only_nonconstant_features_sample": heldout_only[:20],
        "train_only_nonconstant_features_sample": train_only[:20],
    }


def rank_fraction(ranks: pd.Series, denominator: int, cutoff_fraction: float) -> float | None:
    if denominator <= 0 or ranks.empty:
        return None
    return float((ranks <= max(1, int(denominator * cutoff_fraction))).mean())


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


def validation_metrics(scored_all: pd.DataFrame, heldout_scored: pd.DataFrame) -> dict[str, Any]:
    positives = heldout_scored[heldout_scored["heldout_positive"]].copy()
    labels = heldout_scored["heldout_positive"].astype(int).to_numpy()
    scores = pd.to_numeric(heldout_scored["transfer_score"], errors="coerce").fillna(0.0).to_numpy()
    prevalence = float(labels.mean()) if len(labels) else None
    unique_scores = int(pd.Series(scores).nunique()) if len(scores) else 0
    top_score_fraction = float((scores == np.max(scores)).mean()) if len(scores) else None
    metrics: dict[str, Any] = {
        "heldout_rows": int(len(heldout_scored)),
        "heldout_positive_rows": int(labels.sum()),
        "heldout_nonpositive_rows": int(len(labels) - labels.sum()),
        "heldout_positive_prevalence": prevalence,
        "heldout_unique_transfer_scores": unique_scores,
        "heldout_top_score_tie_fraction": top_score_fraction,
        "heldout_mean_transfer_score": float(np.mean(scores)) if len(scores) else None,
        "heldout_positive_mean_transfer_score": float(positives["transfer_score"].mean()) if not positives.empty else None,
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
        "subject_name",
        "predicate",
        "object_name",
        "source_name",
        "source_partition",
        "source_kind",
        "weak_label",
        "label_score",
        "transfer_score",
        "label_reason",
        "source_ref",
    ]
    return frame[[column for column in columns if column in frame.columns]].head(top_n)


def run_one_holdout(
    labels: pd.DataFrame,
    field: str,
    value: str,
    run_dir: Path,
    min_train_rows: int,
    top_n: int,
) -> dict[str, Any]:
    heldout = source_mask(labels, field, value)
    heldout_positive = heldout & labels["weak_label"].isin(POSITIVE_LABELS)
    train_labels = labels.loc[~heldout & (labels["weak_label"] != "abstain")].copy()
    heldout_labels = labels.loc[heldout].copy()
    if int(heldout_positive.sum()) == 0:
        raise ValueError(f"Held-out source has no positive labels: {field}={value}")
    if len(train_labels) < min_train_rows:
        raise ValueError(f"Training set too small after holding out {field}={value}: {len(train_labels)} rows")

    model, feature_columns, train_metrics = train_transfer_model(train_labels)
    scored_all = score_with_model(labels, model, feature_columns)
    scored_all["global_rank"] = np.arange(1, len(scored_all) + 1)
    scored_all["heldout_source"] = source_mask(scored_all, field, value)
    scored_all["heldout_positive"] = scored_all["heldout_source"] & scored_all["weak_label"].isin(POSITIVE_LABELS)
    heldout_scored = scored_all[scored_all["heldout_source"]].copy()
    heldout_scored = heldout_scored.sort_values(["transfer_score", "candidate_uid"], ascending=[False, True])
    metrics = validation_metrics(scored_all, heldout_scored)
    source_diagnostics = source_family_diagnostics(train_labels, heldout_labels)
    feature_diagnostics = feature_transfer_diagnostics(train_labels, heldout_labels, feature_columns)

    slug = f"{field}_{value}".replace(":", "_").replace("/", "_").replace("\\", "_").replace(" ", "_")
    output_dir = run_dir / "reports" / "source_heldout"
    output_dir.mkdir(parents=True, exist_ok=True)
    top_candidates = top_candidate_table(scored_all, scored_all["heldout_source"], top_n)
    top_candidates.to_csv(output_dir / f"{slug}_top_candidates.csv", index=False, encoding="utf-8-sig")
    heldout_scored.to_csv(output_dir / f"{slug}_all_heldout_scored.csv", index=False, encoding="utf-8-sig")
    joblib.dump(
        {
            "model": model,
            "feature_columns": feature_columns,
            "heldout": {"field": field, "value": value},
            "feature_mode": "no_leakage",
            "source_balanced": True,
        },
        output_dir / f"{slug}_temporary_model.joblib",
    )

    return {
        "heldout": {"field": field, "value": value},
        "train_rows": int(len(train_labels)),
        "heldout_rows": int(len(heldout_labels)),
        "heldout_positive_rows": int(heldout_positive.sum()),
        "training_metrics": train_metrics,
        "validation_metrics": metrics,
        "source_family_diagnostics": source_diagnostics,
        "feature_transfer_diagnostics": feature_diagnostics,
        "outputs": {
            "top_candidates": str(output_dir / f"{slug}_top_candidates.csv"),
            "all_heldout_scored": str(output_dir / f"{slug}_all_heldout_scored.csv"),
            "temporary_model": str(output_dir / f"{slug}_temporary_model.joblib"),
        },
    }


def parse_heldout_sources(values: list[str]) -> list[dict[str, str]]:
    specs: list[dict[str, str]] = []
    for value in values:
        if "=" not in value:
            specs.append({"field": "source_name", "value": value})
            continue
        field, source_value = value.split("=", 1)
        specs.append({"field": field.strip(), "value": source_value.strip()})
    return specs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run source-held-out validation for a learning run.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--heldout-source",
        action="append",
        default=[],
        help="Source to hold out. Use 'DrugCentral' or 'source_kind=drug_target_overlay'. Can be repeated.",
    )
    parser.add_argument("--min-positive-rows", type=int, default=50)
    parser.add_argument("--min-train-rows", type=int, default=100)
    parser.add_argument("--top-n", type=int, default=500)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    run_dir = workspace / args.output_root / args.run_id
    views_dir = run_dir / "views"
    label_dir = run_dir / "weak_labels"
    report_dir = run_dir / "reports" / "source_heldout"

    labels = read_parquet(label_dir / "weak_labels.parquet")
    if labels.empty:
        raise FileNotFoundError(f"No weak labels found for run: {args.run_id}")
    drug_overlay = read_parquet(views_dir / "drug_target_overlay.parquet")
    literature = read_parquet(views_dir / "literature_relations.parquet")
    labels = add_source_names(labels, drug_overlay, literature)

    specs = parse_heldout_sources(args.heldout_source)
    if not specs:
        specs = infer_source_specs(labels, args.min_positive_rows)
    if not specs:
        raise ValueError("No held-out sources met the minimum positive-row threshold.")

    source_counts = (
        labels[labels["source_kind"] != "shuffled_negative"]
        .groupby(["source_name", "source_partition", "source_partition_strategy", "source_kind", "weak_label"], dropna=False)
        .size()
        .reset_index(name="rows")
        .sort_values(["source_name", "source_partition", "source_kind", "rows"], ascending=[True, True, True, False])
    )
    report_dir.mkdir(parents=True, exist_ok=True)
    source_counts.to_csv(report_dir / "source_label_distribution.csv", index=False, encoding="utf-8-sig")

    validations = [
        run_one_holdout(
            labels=labels,
            field=spec["field"],
            value=spec["value"],
            run_dir=run_dir,
            min_train_rows=args.min_train_rows,
            top_n=args.top_n,
        )
        for spec in specs
    ]

    summary = {
        "created_at_utc": utc_now(),
        "run_id": args.run_id,
        "scope": "source_heldout_validation",
        "input_rows": {
            "weak_labels": int(len(labels)),
            "drug_target_overlay": int(len(drug_overlay)),
            "literature_relations": int(len(literature)),
        },
        "heldout_sources": validations,
        "outputs": {
            "source_label_distribution": str(report_dir / "source_label_distribution.csv"),
            "summary": str(report_dir / "source_heldout_validation_report.json"),
        },
        "notes": [
            "Temporary held-out models are validation artifacts and do not replace the frozen baseline ranker.",
            "Held-out labels are still weak supervision, but grouped by evidence source to test cross-source transfer.",
            "Validation scores use an uncapped transfer model in no_leakage feature mode with source-balanced weights; precision@k is reported with prevalence-normalized lift and tie diagnostics.",
            "No canonical graph, normalized store, baseline model, or baseline rankings are modified.",
        ],
    }
    (report_dir / "source_heldout_validation_report.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
