"""Run proxy validation and debiased ranking for generalized learning outputs.

The script evaluates weak-supervision separation, pseudo-negative ordering, and
popularity/hub bias. It also emits debiased pathway, target, and drug rankings without
modifying the source learning run outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score


DEFAULT_OUTPUT_ROOT = "learning_runs"
GENERIC_TERMS = {
    "neoplasm",
    "carcinoma",
    "cancer",
    "tumor",
    "tumour",
    "disease",
    "malignant neoplasm",
    "cell",
    "cells",
    "protein",
}
HUB_LABELS = {
    "tp53",
    "tumor protein p53",
    "egfr",
    "epidermal growth factor receptor",
    "myc",
    "stat3",
    "signal transducer and activator of transcription 3",
    "tnf",
    "tumor necrosis factor",
    "mtor",
    "mechanistic target of rapamycin kinase",
    "akt1",
    "mapk1",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def stable_bucket(value: str, buckets: int = 5) -> int:
    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % buckets


def read_parquet(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pq.read_table(path).to_pandas()


def parse_exclude_terms(value: str) -> list[str]:
    return [term.strip().casefold() for term in str(value or "").split(",") if term.strip()]


def filter_excluded(frame: pd.DataFrame, exclude_terms: list[str]) -> pd.DataFrame:
    if frame.empty or not exclude_terms:
        return frame
    string_columns = [column for column in frame.columns if frame[column].dtype == object]
    if not string_columns:
        return frame
    mask = pd.Series(False, index=frame.index)
    for column in string_columns:
        text = frame[column].fillna("").astype(str).str.casefold()
        for term in exclude_terms:
            mask = mask | text.str.contains(term, regex=False)
    return frame.loc[~mask].copy()


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(number) or math.isinf(number):
        return default
    return number


def bounded(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return min(upper, max(lower, value))


def label_text(value: Any) -> str:
    return str(value or "").casefold().strip()


def contains_generic(value: Any) -> bool:
    text = label_text(value)
    if not text:
        return False
    if text in GENERIC_TERMS:
        return True
    return any(text == term or text.endswith(" " + term) for term in GENERIC_TERMS)


def contains_hub(value: Any) -> bool:
    text = label_text(value)
    if not text:
        return False
    return any(text == label or label in text for label in HUB_LABELS)


def normalize_series(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").fillna(0.0)
    max_value = float(numeric.max()) if len(numeric) else 0.0
    if max_value <= 0:
        return numeric
    return numeric / max_value


def precision_at_k(labels: np.ndarray, scores: np.ndarray, k: int) -> float | None:
    if len(labels) == 0:
        return None
    k = min(k, len(labels))
    if k <= 0:
        return None
    order = np.argsort(-scores)[:k]
    return float(np.mean(labels[order]))


def weak_label_holdout_metrics(relations: pd.DataFrame) -> dict[str, Any]:
    frame = relations[relations["weak_label"] != "abstain"].copy()
    frame["bucket"] = frame["candidate_uid"].astype(str).map(stable_bucket)
    valid = frame[frame["bucket"] == 0].copy()
    if valid.empty:
        valid = frame.copy()
    labels = (pd.to_numeric(valid["label_score"], errors="coerce").fillna(0.0) >= 0.5).astype(int).to_numpy()
    scores = pd.to_numeric(valid["priority_score"], errors="coerce").fillna(0.0).to_numpy()
    metrics: dict[str, Any] = {
        "validation_rows": int(len(valid)),
        "positive_rows": int(labels.sum()),
        "negative_rows": int(len(labels) - labels.sum()),
        "precision_at_50": precision_at_k(labels, scores, 50),
        "precision_at_100": precision_at_k(labels, scores, 100),
        "precision_at_500": precision_at_k(labels, scores, 500),
    }
    if len(set(labels.tolist())) > 1:
        metrics["roc_auc"] = float(roc_auc_score(labels, scores))
        metrics["average_precision"] = float(average_precision_score(labels, scores))
        precision, recall, thresholds = precision_recall_curve(labels, scores)
        if len(thresholds):
            f1 = 2 * precision[:-1] * recall[:-1] / np.maximum(precision[:-1] + recall[:-1], 1e-12)
            idx = int(np.argmax(f1))
            metrics["best_f1"] = float(f1[idx])
            metrics["best_threshold"] = float(thresholds[idx])
    return metrics


def pseudo_negative_pairwise_metrics(relations: pd.DataFrame, sample_limit: int = 200_000) -> dict[str, Any]:
    positives = relations[relations["label_score"] >= 0.5].copy()
    negatives = relations[relations["weak_label"] == "negative_shuffled"].copy()
    if positives.empty or negatives.empty:
        return {"pair_count": 0}
    positives["pair_key"] = positives["subject_uid"].astype(str) + "|" + positives["predicate"].astype(str) + "|" + positives["object_type"].astype(str)
    negatives["pair_key"] = negatives["subject_uid"].astype(str) + "|" + negatives["predicate"].astype(str) + "|" + negatives["object_type"].astype(str)
    merged = positives[["pair_key", "priority_score"]].merge(
        negatives[["pair_key", "priority_score"]],
        on="pair_key",
        suffixes=("_positive", "_negative"),
    )
    if len(merged) > sample_limit:
        merged = merged.sample(sample_limit, random_state=13)
    if merged.empty:
        return {"pair_count": 0}
    margin = merged["priority_score_positive"] - merged["priority_score_negative"]
    return {
        "pair_count": int(len(merged)),
        "pairwise_accuracy": float((margin > 0).mean()),
        "mean_margin": float(margin.mean()),
        "min_margin": float(margin.min()),
    }


def enrich_entity_features(frame: pd.DataFrame, entities: pd.DataFrame, uid_column: str) -> pd.DataFrame:
    feature_columns = [
        "entity_uid",
        "graph_degree",
        "canonical_edge_count",
        "literature_mention_count",
        "literature_article_count",
        "literature_support_count",
    ]
    lookup = entities[[column for column in feature_columns if column in entities.columns]].copy()
    return frame.merge(lookup, left_on=uid_column, right_on="entity_uid", how="left").drop(columns=["entity_uid"], errors="ignore")


def bias_correlations(frame: pd.DataFrame, uid_column: str, entities: pd.DataFrame) -> dict[str, Any]:
    if frame.empty:
        return {}
    enriched = enrich_entity_features(frame, entities, uid_column)
    columns = [
        "candidate_count",
        "positive_signal_count",
        "graph_degree",
        "canonical_edge_count",
        "literature_mention_count",
        "literature_article_count",
        "literature_support_count",
    ]
    result: dict[str, Any] = {}
    for column in columns:
        if column not in enriched.columns:
            continue
        pair = enriched[["priority_score", column]].copy()
        pair[column] = pd.to_numeric(pair[column], errors="coerce")
        pair = pair.dropna()
        if len(pair) < 3 or pair[column].nunique() < 2:
            continue
        result[f"spearman_priority_vs_{column}"] = float(pair["priority_score"].corr(pair[column], method="spearman"))
    top = enriched.head(100)
    result["top100_generic_counterpart_fraction"] = float(top.get("top_counterpart_name", pd.Series(dtype=str)).map(contains_generic).mean()) if not top.empty else None
    result["top100_hub_label_fraction"] = float(
        (
            top.get("display_name", pd.Series(dtype=str)).map(contains_hub)
            | top.get("top_counterpart_name", pd.Series(dtype=str)).map(contains_hub)
        ).mean()
    ) if not top.empty else None
    return result


def debias_rankings(frame: pd.DataFrame, entities: pd.DataFrame, uid_column: str, kind: str) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()
    result = enrich_entity_features(frame, entities, uid_column)
    for column in ("priority_score", "mean_priority_score", "candidate_count", "positive_signal_count", "graph_degree", "literature_mention_count", "literature_support_count"):
        if column in result.columns:
            result[column] = pd.to_numeric(result[column], errors="coerce").fillna(0.0)

    raw = result["priority_score"]
    candidate_count = result.get("candidate_count", pd.Series(0.0, index=result.index))
    positive_signal_count = result.get("positive_signal_count", pd.Series(0.0, index=result.index))
    graph_degree = result.get("graph_degree", pd.Series(0.0, index=result.index))
    mention_count = result.get("literature_mention_count", pd.Series(0.0, index=result.index))
    support_count = result.get("literature_support_count", pd.Series(0.0, index=result.index))

    evidence_density = positive_signal_count / np.sqrt(np.maximum(candidate_count, 1.0))
    density_norm = normalize_series(evidence_density)
    popularity_penalty = 1.0 + 0.10 * np.log1p(candidate_count) + 0.04 * np.log1p(graph_degree) + 0.05 * np.log1p(mention_count)
    generic_penalty = result.get("top_counterpart_name", pd.Series("", index=result.index)).map(lambda value: 0.65 if contains_generic(value) else 1.0)
    hub_penalty = (
        result.get("display_name", pd.Series("", index=result.index)).map(lambda value: 0.72 if contains_hub(value) else 1.0)
        * result.get("top_counterpart_name", pd.Series("", index=result.index)).map(lambda value: 0.82 if contains_hub(value) else 1.0)
    )
    support_bonus = 0.85 + 0.15 * normalize_series(support_count)
    density_bonus = 0.55 + 0.45 * density_norm
    debiased = raw * density_bonus * support_bonus * generic_penalty * hub_penalty / popularity_penalty

    result[f"{kind}_uid"] = result[uid_column]
    result["raw_priority_score"] = raw.round(6)
    result["debiased_score"] = debiased.round(6)
    result["evidence_density"] = evidence_density.round(6)
    result["popularity_penalty"] = popularity_penalty.round(6)
    result["generic_penalty"] = generic_penalty.round(6)
    result["hub_penalty"] = hub_penalty.round(6)
    result["support_bonus"] = support_bonus.round(6)
    result["debias_notes"] = [
        ";".join(note for note in notes if note)
        for notes in zip(
            result.get("top_counterpart_name", pd.Series("", index=result.index)).map(lambda value: "generic_counterpart" if contains_generic(value) else ""),
            result.get("display_name", pd.Series("", index=result.index)).map(lambda value: "hub_label" if contains_hub(value) else ""),
            result.get("top_counterpart_name", pd.Series("", index=result.index)).map(lambda value: "hub_counterpart" if contains_hub(value) else ""),
        )
    ]
    order_cols = [
        f"{kind}_uid",
        "display_name",
        "debiased_score",
        "raw_priority_score",
        "mean_priority_score",
        "candidate_count",
        "positive_signal_count",
        "evidence_density",
        "popularity_penalty",
        "generic_penalty",
        "hub_penalty",
        "support_bonus",
        "top_counterpart_name",
        "top_label_reason",
        "debias_notes",
    ]
    result = result.sort_values(["debiased_score", "positive_signal_count", "evidence_density"], ascending=False)
    return result[[column for column in order_cols if column in result.columns]]


def validation_failures(
    pathways: pd.DataFrame,
    targets: pd.DataFrame,
    drugs: pd.DataFrame,
    debiased_pathways: pd.DataFrame,
    debiased_targets: pd.DataFrame,
    debiased_drugs: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    def add(kind: str, frame: pd.DataFrame, uid_column: str) -> None:
        if frame.empty:
            return
        for idx, row in frame.head(200).reset_index(drop=True).iterrows():
            issues: list[str] = []
            if contains_generic(row.get("top_counterpart_name", "")):
                issues.append("generic_counterpart")
            if contains_hub(row.get("display_name", "")) or contains_hub(row.get("top_counterpart_name", "")):
                issues.append("hub_label_or_counterpart")
            candidate_count = safe_float(row.get("candidate_count"))
            positives = safe_float(row.get("positive_signal_count"))
            if candidate_count >= 100 and positives / max(candidate_count, 1.0) < 0.08:
                issues.append("high_volume_low_density")
            if safe_float(row.get("priority_score")) >= 0.9 and positives <= 1 and candidate_count >= 5:
                issues.append("high_score_sparse_positive_signal")
            if issues:
                rows.append(
                    {
                        "kind": kind,
                        "rank": int(idx + 1),
                        "entity_uid": row.get(uid_column, ""),
                        "display_name": row.get("display_name", ""),
                        "priority_score": safe_float(row.get("priority_score")),
                        "candidate_count": candidate_count,
                        "positive_signal_count": positives,
                        "top_counterpart_name": row.get("top_counterpart_name", ""),
                        "issues": ";".join(issues),
                    }
                )

    add("pathway", pathways, "pathway_uid")
    add("target", targets, "target_uid")
    add("drug", drugs, "drug_uid")

    for kind, raw, debiased, uid_column in (
        ("pathway", pathways, debiased_pathways, "pathway_uid"),
        ("target", targets, debiased_targets, "target_uid"),
        ("drug", drugs, debiased_drugs, "drug_uid"),
    ):
        if raw.empty or debiased.empty:
            continue
        raw_top = set(raw.head(50)[uid_column].astype(str))
        debiased_top = set(debiased.head(50)[f"{kind}_uid"].astype(str))
        dropped = raw_top - debiased_top
        for uid in sorted(dropped)[:50]:
            raw_row = raw[raw[uid_column].astype(str) == uid].iloc[0]
            rows.append(
                {
                    "kind": kind,
                    "rank": int(raw.index[raw[uid_column].astype(str) == uid][0] + 1),
                    "entity_uid": uid,
                    "display_name": raw_row.get("display_name", ""),
                    "priority_score": safe_float(raw_row.get("priority_score")),
                    "candidate_count": safe_float(raw_row.get("candidate_count")),
                    "positive_signal_count": safe_float(raw_row.get("positive_signal_count")),
                    "top_counterpart_name": raw_row.get("top_counterpart_name", ""),
                    "issues": "dropped_from_debiased_top50",
                }
            )
    return pd.DataFrame(rows)


def top_overlap(raw: pd.DataFrame, debiased: pd.DataFrame, raw_uid: str, debiased_uid: str, k: int) -> float | None:
    if raw.empty or debiased.empty:
        return None
    a = set(raw.head(k)[raw_uid].astype(str))
    b = set(debiased.head(k)[debiased_uid].astype(str))
    if not a:
        return None
    return float(len(a & b) / len(a))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run proxy validation and debiasing for a generalized learning run.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--exclude-terms",
        default="",
        help="Comma-separated terms to omit from validation/debiasing outputs and metrics.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    run_dir = workspace / args.output_root / args.run_id
    ranking_dir = run_dir / "rankings"
    views_dir = run_dir / "views"
    report_dir = run_dir / "reports"

    entities = read_parquet(views_dir / "entities.parquet")
    relations = read_parquet(ranking_dir / "relation_priorities.parquet")
    pathways = read_parquet(ranking_dir / "pathway_priorities.parquet")
    targets = read_parquet(ranking_dir / "target_priorities.parquet")
    drugs = read_parquet(ranking_dir / "drug_priorities.parquet")
    exclude_terms = parse_exclude_terms(args.exclude_terms)
    relations = filter_excluded(relations, exclude_terms)
    pathways = filter_excluded(pathways, exclude_terms)
    targets = filter_excluded(targets, exclude_terms)
    drugs = filter_excluded(drugs, exclude_terms)

    holdout = weak_label_holdout_metrics(relations)
    pairwise = pseudo_negative_pairwise_metrics(relations)
    pathway_bias = bias_correlations(pathways, "pathway_uid", entities)
    target_bias = bias_correlations(targets, "target_uid", entities)
    drug_bias = bias_correlations(drugs, "drug_uid", entities)

    debiased_pathways = debias_rankings(pathways, entities, "pathway_uid", "pathway")
    debiased_targets = debias_rankings(targets, entities, "target_uid", "target")
    debiased_drugs = debias_rankings(drugs, entities, "drug_uid", "drug")
    failures = validation_failures(pathways, targets, drugs, debiased_pathways, debiased_targets, debiased_drugs)

    report_dir.mkdir(parents=True, exist_ok=True)
    debiased_pathways.to_csv(report_dir / "debiased_top_pathways.csv", index=False, encoding="utf-8-sig")
    debiased_targets.to_csv(report_dir / "debiased_top_targets.csv", index=False, encoding="utf-8-sig")
    debiased_drugs.to_csv(report_dir / "debiased_top_drugs.csv", index=False, encoding="utf-8-sig")
    failures.to_csv(report_dir / "validation_failures.csv", index=False, encoding="utf-8-sig")

    report = {
        "created_at_utc": utc_now(),
        "run_id": args.run_id,
        "scope": "generalized_proxy_validation_and_debiasing",
        "exclude_term_count": len(exclude_terms),
        "input_rows": {
            "entities": int(len(entities)),
            "relations": int(len(relations)),
            "pathways": int(len(pathways)),
            "targets": int(len(targets)),
            "drugs": int(len(drugs)),
        },
        "proxy_validation": {
            "weak_label_holdout": holdout,
            "pseudo_negative_pairwise": pairwise,
            "top50_overlap_after_debias": {
                "pathways": top_overlap(pathways, debiased_pathways, "pathway_uid", "pathway_uid", 50),
                "targets": top_overlap(targets, debiased_targets, "target_uid", "target_uid", 50),
                "drugs": top_overlap(drugs, debiased_drugs, "drug_uid", "drug_uid", 50),
            },
        },
        "bias_diagnostics": {
            "pathways": pathway_bias,
            "targets": target_bias,
            "drugs": drug_bias,
        },
        "outputs": {
            "debiased_top_pathways": str(report_dir / "debiased_top_pathways.csv"),
            "debiased_top_targets": str(report_dir / "debiased_top_targets.csv"),
            "debiased_top_drugs": str(report_dir / "debiased_top_drugs.csv"),
            "validation_failures": str(report_dir / "validation_failures.csv"),
        },
        "notes": [
            "Proxy metrics use weak labels and pseudo-negatives, not manual gold labels.",
            "Debiased scores reduce popularity, hub, and generic-label effects; they are still research-prioritization scores.",
            "No canonical graph, normalized store, or learning source table is modified.",
        ],
    }
    (report_dir / "proxy_validation_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"proxy_validation": report["proxy_validation"], "outputs": report["outputs"]}, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
