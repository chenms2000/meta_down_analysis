"""Train a weak-supervision priority ranker from candidate labels."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import average_precision_score, mean_absolute_error, roc_auc_score


DEFAULT_OUTPUT_ROOT = "learning_runs"
NUMERIC_COLUMNS = [
    "p_literature",
    "literature_sentence_count",
    "literature_distinct_article_count",
    "supported_existing_edge_count",
    "overlay_confidence",
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
    "literature_support_score",
    "independent_article_score",
    "curated_edge_score",
    "overlay_score",
    "biomedbert_relevance_mean",
    "biomedbert_relevance_max",
    "biomedbert_relevance_min",
    "biomedbert_relevance_scored_sentence_count",
    "biomedbert_relevance_evidence_candidate_fraction",
    "biomedbert_weighted_literature_score",
    "biomedbert_relevance_coverage",
    "scispacy_entity_candidate_count",
    "scispacy_entity_candidate_log",
    "scispacy_unique_surface_count",
    "scispacy_label_count",
    "scispacy_candidate_sentence_fraction",
    "scispacy_candidate_density",
    "scispacy_evidence_complexity_score",
]
CATEGORICAL_COLUMNS = ["subject_type", "object_type", "predicate", "relation_family", "subject_object_family", "source_kind", "label_reason"]
FEATURE_MODES = {"full", "no_leakage"}
LEAKAGE_CATEGORICAL_COLUMNS = {"source_kind", "label_reason"}

DISEASE_MODEL_PATHWAY_TERMS = (
    "cancer",
    "carcinoma",
    "melanoma",
    "glioma",
    "glioblastoma",
    "leukemia",
    "leukaemia",
    "lymphoma",
    "myeloma",
    "sarcoma",
    "adenoma",
    "tumor",
    "tumour",
    "metastasis",
    "metastatic",
    "disease",
    "disorder",
    "syndrome",
    "deficiency",
    "defective",
    "cell line",
    "cancer cells",
    "carcinoma cells",
    "tumor cells",
    "tumour cells",
)

PRIMARY_METABOLIC_ALLOW_TERMS = (
    "glycolysis",
    "gluconeogenesis",
    "tca",
    "citric acid",
    "tricarboxylic",
    "oxidative phosphorylation",
    "fatty acid",
    "beta oxidation",
    "carnitine",
    "amino acid",
    "glutamine",
    "glutamate",
    "arginine",
    "nitric oxide",
    "purine",
    "pyrimidine",
    "nucleotide",
    "one carbon",
    "methionine",
    "folate",
    "tryptophan",
    "kynurenine",
    "glutathione",
    "redox",
    "ferroptosis",
    "lipid",
    "phospholipid",
    "sphingolipid",
    "glycosylation",
    "glycan",
    "uremic",
    "organic anion",
    "osmolyte",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def stable_bucket(value: str, buckets: int = 5) -> int:
    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % buckets


def write_parquet(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), path)


def feature_columns_for_mode(feature_mode: str) -> tuple[list[str], list[str]]:
    if feature_mode not in FEATURE_MODES:
        raise ValueError(f"Unknown feature mode: {feature_mode}")
    categorical_columns = CATEGORICAL_COLUMNS
    if feature_mode == "no_leakage":
        categorical_columns = [column for column in CATEGORICAL_COLUMNS if column not in LEAKAGE_CATEGORICAL_COLUMNS]
    return NUMERIC_COLUMNS, categorical_columns


def build_feature_matrix(
    frame: pd.DataFrame,
    feature_columns: list[str] | None = None,
    feature_mode: str = "full",
) -> tuple[pd.DataFrame, list[str]]:
    numeric_columns, categorical_columns = feature_columns_for_mode(feature_mode)
    numeric = frame.copy()
    for column in numeric_columns:
        values = numeric[column] if column in numeric.columns else pd.Series(0.0, index=numeric.index)
        numeric[column] = pd.to_numeric(values, errors="coerce").fillna(0.0)
    pieces = [numeric[numeric_columns]]
    available_categoricals = [column for column in categorical_columns if column in numeric.columns]
    categoricals = pd.get_dummies(numeric[available_categoricals].fillna(""), prefix=available_categoricals)
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


def training_sample_weights(train_frame: pd.DataFrame, source_balanced: bool) -> tuple[np.ndarray, dict[str, Any]]:
    weights = 1.0 + 0.5 * (train_frame["weak_label"] == "strong_positive").astype(float).to_numpy()
    diagnostics: dict[str, Any] = {
        "source_balanced": bool(source_balanced),
        "strong_positive_multiplier": 1.5,
    }
    if not source_balanced:
        diagnostics["min_weight"] = float(np.min(weights)) if len(weights) else None
        diagnostics["max_weight"] = float(np.max(weights)) if len(weights) else None
        return weights, diagnostics

    group_columns = [column for column in ("source_kind", "weak_label") if column in train_frame.columns]
    if not group_columns:
        diagnostics["balance_group_columns"] = []
        diagnostics["min_weight"] = float(np.min(weights)) if len(weights) else None
        diagnostics["max_weight"] = float(np.max(weights)) if len(weights) else None
        return weights, diagnostics

    group_sizes = train_frame.groupby(group_columns, dropna=False)["candidate_uid"].transform("size").astype(float)
    group_count = int(train_frame[group_columns].drop_duplicates().shape[0])
    target_group_size = max(1.0, len(train_frame) / max(group_count, 1))
    balance = np.sqrt(target_group_size / np.maximum(group_sizes.to_numpy(), 1.0))
    balance = np.clip(balance, 0.35, 6.0)
    weights = weights * balance
    diagnostics.update(
        {
            "balance_group_columns": group_columns,
            "balance_group_count": group_count,
            "balance_target_group_size": float(target_group_size),
            "min_weight": float(np.min(weights)) if len(weights) else None,
            "max_weight": float(np.max(weights)) if len(weights) else None,
            "mean_weight": float(np.mean(weights)) if len(weights) else None,
        }
    )
    return weights, diagnostics


def pathway_appendix_classification(display_name: Any) -> tuple[bool, str]:
    text = str(display_name or "").casefold()
    if not text:
        return False, ""
    has_disease_model = any(term in text for term in DISEASE_MODEL_PATHWAY_TERMS)
    if not has_disease_model:
        return False, ""
    has_primary_metabolic = any(term in text for term in PRIMARY_METABOLIC_ALLOW_TERMS)
    if has_primary_metabolic and not any(term in text for term in (" in cancer", " cancer cell", "tumor cell", "cell line")):
        return False, ""
    return True, "disease_or_model_specific_pathway_appendix_in_generalized_report"


def train_model(
    labels: pd.DataFrame,
    feature_mode: str = "full",
    source_balanced: bool = False,
) -> tuple[HistGradientBoostingRegressor, list[str], dict[str, Any]]:
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

    x_train, feature_columns = build_feature_matrix(train_frame, feature_mode=feature_mode)
    x_valid, _ = build_feature_matrix(valid_frame, feature_columns, feature_mode=feature_mode)
    y_train = train_frame["label_score"].astype(float).to_numpy()
    y_valid = valid_frame["label_score"].astype(float).to_numpy()
    weights, weight_diagnostics = training_sample_weights(train_frame, source_balanced=source_balanced)

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
        "feature_mode": feature_mode,
        "source_balanced": bool(source_balanced),
        "feature_count": int(len(feature_columns)),
        "weight_diagnostics": weight_diagnostics,
        "mae": float(mean_absolute_error(y_valid, predicted)),
        "mean_prediction_positive": float(np.mean(predicted[y_valid >= 0.5])) if np.any(y_valid >= 0.5) else None,
        "mean_prediction_negative": float(np.mean(predicted[y_valid < 0.5])) if np.any(y_valid < 0.5) else None,
    }
    binary = (y_valid >= 0.5).astype(int)
    if len(set(binary.tolist())) > 1:
        metrics["roc_auc"] = float(roc_auc_score(binary, predicted))
        metrics["average_precision"] = float(average_precision_score(binary, predicted))
    return model, feature_columns, metrics


def predict_all(
    labels: pd.DataFrame,
    model: HistGradientBoostingRegressor,
    feature_columns: list[str],
    feature_mode: str = "full",
) -> pd.DataFrame:
    matrix, _ = build_feature_matrix(labels, feature_columns, feature_mode=feature_mode)
    result = labels.copy()
    result["priority_score_raw"] = np.clip(model.predict(matrix), 0.0, 1.0)
    caps = {
        "strong_positive": 1.0,
        "medium_positive": 0.82,
        "weak_positive": 0.62,
        "abstain": 0.25,
        "negative_shuffled": 0.05,
        "negative_refute": 0.05,
    }
    result["evidence_cap"] = result["weak_label"].map(caps).fillna(0.25).astype(float)
    result["priority_score"] = np.minimum(result["priority_score_raw"], result["evidence_cap"])
    result["priority_score"] = result["priority_score"].round(6)
    result["priority_score_raw"] = result["priority_score_raw"].round(6)
    return result.sort_values(["priority_score", "label_score"], ascending=False)


def aggregate_priorities(predictions: pd.DataFrame, entity_kind: str) -> pd.DataFrame:
    if entity_kind == "drug":
        frame = predictions[predictions["subject_type"] == "drug"].copy()
        frame["entity_uid"] = frame["subject_uid"]
        frame["entity_name"] = frame["subject_name"]
        frame["counterpart_uid"] = frame["object_uid"]
        frame["counterpart_name"] = frame["object_name"]
    else:
        object_frame = predictions[predictions["object_type"] == entity_kind].copy()
        subject_frame = predictions[predictions["subject_type"] == entity_kind].copy()
        object_frame["entity_uid"] = object_frame["object_uid"]
        object_frame["entity_name"] = object_frame["object_name"]
        object_frame["counterpart_uid"] = object_frame["subject_uid"]
        object_frame["counterpart_name"] = object_frame["subject_name"]
        subject_frame["entity_uid"] = subject_frame["subject_uid"]
        subject_frame["entity_name"] = subject_frame["subject_name"]
        subject_frame["counterpart_uid"] = subject_frame["object_uid"]
        subject_frame["counterpart_name"] = subject_frame["object_name"]
        frame = pd.concat([object_frame, subject_frame], ignore_index=True)
    if frame.empty:
        return pd.DataFrame(columns=[f"{entity_kind}_uid", "display_name", "priority_score", "candidate_count"])
    grouped = []
    for uid, group in frame.groupby("entity_uid"):
        best = group.sort_values("priority_score", ascending=False).iloc[0]
        grouped.append(
            {
                f"{entity_kind}_uid": uid,
                "display_name": best.get("entity_name", ""),
                "priority_score": float(group["priority_score"].max()),
                "priority_score_raw": float(pd.to_numeric(group.get("priority_score_raw", group["priority_score"]), errors="coerce").fillna(0.0).max()),
                "mean_priority_score": float(group["priority_score"].mean()),
                "candidate_count": int(len(group)),
                "positive_signal_count": int((group["label_score"] >= 0.5).sum()),
                "top_source_kind": best.get("source_kind", ""),
                "top_label_reason": best.get("label_reason", ""),
                "top_counterpart_uid": best.get("counterpart_uid", ""),
                "top_counterpart_name": best.get("counterpart_name", ""),
            }
        )
    result = pd.DataFrame(grouped)
    if entity_kind == "pathway" and not result.empty:
        appendix_flags = [pathway_appendix_classification(value) for value in result["display_name"]]
        result["appendix"] = [flag for flag, _reason in appendix_flags]
        result["downgrade_reason"] = [reason for _flag, reason in appendix_flags]
        result.loc[result["appendix"], "priority_score"] = result.loc[result["appendix"], "priority_score"].clip(upper=0.25)
        result.loc[result["appendix"], "mean_priority_score"] = result.loc[result["appendix"], "mean_priority_score"].clip(upper=0.25)
        result["display_tier"] = np.where(result["appendix"], "appendix", "primary")
    return result.sort_values(["priority_score", "positive_signal_count"], ascending=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train a priority ranker from weak labels.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--feature-mode",
        choices=sorted(FEATURE_MODES),
        default="full",
        help="Feature set to use. no_leakage removes source_kind and label_reason categorical shortcuts.",
    )
    parser.add_argument(
        "--source-balanced",
        action="store_true",
        help="Apply inverse group-size weights over source_kind and weak_label groups.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    run_dir = workspace / args.output_root / args.run_id
    label_dir = run_dir / "weak_labels"
    ranking_dir = run_dir / "rankings"
    model_dir = run_dir / "models"
    report_dir = run_dir / "reports"

    labels = pq.read_table(label_dir / "weak_labels.parquet").to_pandas()
    model, feature_columns, metrics = train_model(labels, feature_mode=args.feature_mode, source_balanced=args.source_balanced)
    predictions = predict_all(labels, model, feature_columns, feature_mode=args.feature_mode)

    relation_priority_columns = [
        "candidate_uid",
        "subject_uid",
        "subject_type",
        "subject_name",
        "predicate",
        "object_uid",
        "object_type",
        "object_name",
        "source_kind",
        "weak_label",
        "weak_label_class",
        "label_score",
        "priority_score",
        "priority_score_raw",
        "evidence_cap",
        "label_reason",
        "group_id",
        "paper_id",
        "publication_year",
        "split",
        "temporal_holdout",
        "context_holdout",
        "source_ref",
        "biomedbert_relevance_mean",
        "biomedbert_relevance_max",
        "biomedbert_weighted_literature_score",
        "biomedbert_relevance_coverage",
        "scispacy_entity_candidate_log",
        "scispacy_candidate_density",
        "scispacy_evidence_complexity_score",
    ]
    relation_priorities = predictions[[column for column in relation_priority_columns if column in predictions.columns]].copy()
    pathway_priorities = aggregate_priorities(predictions, "pathway")
    target_priorities = aggregate_priorities(predictions, "target")
    drug_priorities = aggregate_priorities(predictions, "drug")

    write_parquet(ranking_dir / "relation_priorities.parquet", relation_priorities)
    write_parquet(ranking_dir / "pathway_priorities.parquet", pathway_priorities)
    write_parquet(ranking_dir / "target_priorities.parquet", target_priorities)
    write_parquet(ranking_dir / "drug_priorities.parquet", drug_priorities)

    model_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "model": model,
            "feature_columns": feature_columns,
            "feature_mode": args.feature_mode,
            "source_balanced": bool(args.source_balanced),
        },
        model_dir / "priority_ranker.joblib",
    )

    report = {
        "created_at_utc": utc_now(),
        "run_id": args.run_id,
        "feature_mode": args.feature_mode,
        "source_balanced": bool(args.source_balanced),
        "metrics": metrics,
        "rows": {
            "weak_labels": int(len(labels)),
            "relation_priorities": int(len(relation_priorities)),
            "pathway_priorities": int(len(pathway_priorities)),
            "target_priorities": int(len(target_priorities)),
            "drug_priorities": int(len(drug_priorities)),
        },
        "outputs": {
            "relation_priorities": str(ranking_dir / "relation_priorities.parquet"),
            "pathway_priorities": str(ranking_dir / "pathway_priorities.parquet"),
            "target_priorities": str(ranking_dir / "target_priorities.parquet"),
            "drug_priorities": str(ranking_dir / "drug_priorities.parquet"),
            "model": str(model_dir / "priority_ranker.joblib"),
        },
        "notes": [
            "Priority scores are research-ranking outputs, not truth probabilities.",
            "Metrics are proxy validation over weak labels and pseudo-negatives.",
            "scispaCy and BiomedBERT/PubMedBERT values are used only as preprocessing, relevance, and ranking features; they are not truth labels.",
            "no_leakage feature mode removes source_kind and label_reason categorical shortcuts.",
            "source_balanced training reweights source_kind x weak_label groups to reduce source-volume dominance.",
        ],
    }
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "priority_ranker_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"metrics": metrics, "rows": report["rows"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
