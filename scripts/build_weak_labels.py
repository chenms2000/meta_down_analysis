"""Build weak supervision labels for research-prioritization candidates."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


DEFAULT_OUTPUT_ROOT = "learning_runs"
POSITIVE_CLASSES = {"confirm", "support_direction", "support", "supported_existing_edge"}
WEAK_CLASSES = {"novel_candidate", "candidate", "uncertain"}
NEGATIVE_CLASSES = {"conflict", "refute", "refuted", "opposite_direction"}
CONTEXT_MISMATCH_CLASSES = {"context_mismatch", "wrong_context", "out_of_context"}
TOKEN_RE = re.compile(r"[a-z0-9][a-z0-9_\-]{1,}")
METABOLIC_THEME_TERMS = {
    "glycolysis",
    "lactate",
    "pyruvate",
    "tca",
    "citrate",
    "glutamine",
    "glutamate",
    "fatty",
    "lipid",
    "carnitine",
    "serine",
    "folate",
    "one",
    "carbon",
    "purine",
    "pyrimidine",
    "redox",
    "glutathione",
    "oxidative",
    "phosphorylation",
}
CELL_CONTEXT_TERMS = {
    "cell",
    "cells",
    "cellline",
    "cell_line",
    "epithelial",
    "fibroblast",
    "immune",
    "macrophage",
    "tcell",
    "bcell",
    "stromal",
    "endothelial",
    "tumor",
    "tumour",
    "cancer",
    "lineage",
    "tissue",
}
DISEASE_TERMS = {"disease", "cancer", "tumor", "tumour", "carcinoma", "leukemia", "lymphoma"}
PATHWAY_TERMS = {"pathway", "process", "reaction", "metabolic", "metabolism"}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def stable_hash(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(number) or math.isinf(number):
        return default
    return number


def safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def nullable_int(value: Any) -> int | None:
    try:
        number = int(float(value))
    except (TypeError, ValueError):
        return None
    if number <= 0:
        return None
    return number


def log1p(value: Any) -> float:
    return math.log1p(max(0.0, safe_float(value)))


def token_set(*values: Any) -> set[str]:
    tokens: set[str] = set()
    for value in values:
        if isinstance(value, (list, tuple, set)):
            tokens.update(token_set(*value))
            continue
        tokens.update(TOKEN_RE.findall(str(value or "").casefold()))
    return {token for token in tokens if len(token) > 1}


def jaccard(left: set[str], right: set[str]) -> float:
    if not left or not right:
        return 0.0
    return len(left & right) / len(left | right)


def overlap_count(left: set[str], right: set[str], vocabulary: set[str] | None = None) -> int:
    if vocabulary is not None:
        left = left & vocabulary
        right = right & vocabulary
    return len(left & right)


def uid_set(value: Any) -> set[str]:
    if isinstance(value, (list, tuple, set)):
        return {str(item) for item in value if str(item or "")}
    text = str(value or "")
    return {item.strip() for item in re.split(r"[;|,\s]+", text) if item.strip()}


def relation_family(subject_type: str, object_type: str, predicate: str) -> str:
    subject = str(subject_type or "").casefold()
    obj = str(object_type or "").casefold()
    pred_tokens = token_set(predicate)
    if obj == "disease" or pred_tokens & DISEASE_TERMS:
        return f"{subject}_disease"
    if obj == "pathway" or pred_tokens & PATHWAY_TERMS:
        return f"{subject}_pathway"
    if subject == "metabolite" and obj in {"gene", "target", "protein"}:
        return f"metabolite_{obj}"
    return f"{subject}_{obj}".strip("_") or "unknown"


def as_list(value: Any) -> list[str]:
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        value = value.tolist()
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value if str(item or "")]
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    text = str(value)
    return [text] if text else []


def write_parquet(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), path)


def parse_exclude_terms(value: str) -> list[str]:
    return [term.strip().casefold() for term in str(value or "").split(",") if term.strip()]


def stable_bucket(value: str, buckets: int = 100) -> int:
    return int(stable_hash(value)[:12], 16) % buckets


def split_for_group(group_id: str, temporal: bool = False, context_holdout: bool = False) -> str:
    if temporal:
        return "temporal_holdout"
    if context_holdout:
        return "context_holdout"
    bucket = stable_bucket(group_id or "ungrouped", buckets=100)
    if bucket < 15:
        return "test"
    if bucket < 30:
        return "validation"
    return "train"


def five_class_weak_label(label: str, reason: str = "") -> str:
    normalized = str(label or "").casefold()
    reason_text = str(reason or "").casefold()
    if normalized in CONTEXT_MISMATCH_CLASSES or "context_mismatch" in reason_text:
        return "context_mismatch"
    if normalized in NEGATIVE_CLASSES or normalized.startswith("negative") or "refute" in reason_text:
        return "contradicted"
    if normalized in {"strong_positive", "medium_positive", "weak_positive"}:
        return "supported"
    if normalized in {"negative_shuffled"}:
        return "unrelated"
    return "insufficient_evidence"


def first_group_id(literature: dict[str, Any], overlay: dict[str, Any]) -> tuple[str, str]:
    pmids = as_list(literature.get("pmids")) or as_list(literature.get("pmid"))
    if pmids:
        return "paper_id", sorted(pmids)[0]
    for key in ("study_id", "dataset_id", "source_record_id", "support_uid"):
        value = literature.get(key, overlay.get(key, ""))
        if value:
            return key, str(value)
    return "candidate_uid", ""


def publication_year_fields(literature: dict[str, Any], overlay: dict[str, Any]) -> tuple[int | None, int | None, int | None]:
    min_year = nullable_int(literature.get("min_publication_year", overlay.get("min_publication_year", "")))
    max_year = nullable_int(literature.get("max_publication_year", overlay.get("max_publication_year", "")))
    year = nullable_int(literature.get("publication_year", literature.get("year", overlay.get("publication_year", ""))))
    if year is None:
        year = max_year or min_year
    if min_year is None:
        min_year = year
    if max_year is None:
        max_year = year
    return year, min_year, max_year


def filter_excluded_labels(labels: pd.DataFrame, exclude_terms: list[str]) -> tuple[pd.DataFrame, int]:
    if labels.empty or not exclude_terms:
        return labels, 0
    string_columns = [column for column in labels.columns if labels[column].dtype == object]
    if not string_columns:
        return labels, 0
    mask = pd.Series(False, index=labels.index)
    for column in string_columns:
        text = labels[column].fillna("").astype(str).str.casefold()
        for term in exclude_terms:
            mask = mask | text.str.contains(term, regex=False)
    return labels.loc[~mask].copy(), int(mask.sum())


def entity_feature_lookup(entities: pd.DataFrame) -> dict[str, dict[str, Any]]:
    columns = [
        "entity_uid",
        "node_type",
        "canonical_name",
        "display_name",
        "document",
        "graph_degree",
        "canonical_edge_count",
        "literature_mention_count",
        "literature_article_count",
        "literature_support_count",
        "prediction_overlay_count",
        "pathway_context_count",
        "pathway_context_uids",
    ]
    frame = entities[[column for column in columns if column in entities.columns]].copy()
    return {str(row["entity_uid"]): row for row in frame.to_dict(orient="records")}


def entities_by_type(entities: pd.DataFrame) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for node_type, group in entities.groupby("node_type"):
        values = group["entity_uid"].astype(str).tolist()
        if values:
            result[str(node_type)] = values
    return result


def relation_label(row: pd.Series) -> tuple[str, float, str]:
    support_class = str(row.get("support_class", "") or "").casefold()
    p_literature = safe_float(row.get("p_literature"))
    articles = safe_int(row.get("distinct_article_count"))
    existing_edges = safe_int(row.get("supported_existing_edge_count"))
    sentences = safe_int(row.get("evidence_sentence_count"))

    if support_class in NEGATIVE_CLASSES:
        return "negative_refute", 0.0, "literature_refute"
    if support_class in POSITIVE_CLASSES and p_literature >= 0.85 and articles >= 2 and existing_edges > 0:
        return "strong_positive", 1.0, "curated_and_literature"
    if support_class in POSITIVE_CLASSES and p_literature >= 0.70 and sentences >= 2:
        return "medium_positive", 0.78, "literature_support"
    if support_class in WEAK_CLASSES and p_literature >= 0.75 and articles >= 2:
        return "weak_positive", 0.58, "novel_literature_candidate"
    if p_literature >= 0.90 and articles >= 3:
        return "weak_positive", 0.55, "high_literature_probability"
    return "abstain", 0.0, "insufficient_or_context_only"


def feature_row(
    candidate_uid: str,
    subject_uid: str,
    subject_type: str,
    object_uid: str,
    object_type: str,
    predicate: str,
    source_kind: str,
    entity_features: dict[str, dict[str, Any]],
    literature: dict[str, Any] | None = None,
    overlay: dict[str, Any] | None = None,
    weak_label: str = "abstain",
    label_score: float = 0.0,
    label_reason: str = "",
) -> dict[str, Any]:
    literature = literature or {}
    overlay = overlay or {}
    subject = entity_features.get(subject_uid, {})
    obj = entity_features.get(object_uid, {})
    subject_tokens = token_set(subject.get("display_name"), subject.get("canonical_name"), subject.get("document"))
    object_tokens = token_set(obj.get("display_name"), obj.get("canonical_name"), obj.get("document"))
    p_literature = safe_float(literature.get("p_literature"))
    distinct_articles = safe_int(literature.get("distinct_article_count"))
    sentence_count = safe_int(literature.get("evidence_sentence_count"))
    existing_edge_count = safe_int(literature.get("supported_existing_edge_count"))
    overlay_confidence = safe_float(overlay.get("confidence"))
    group_key, group_value = first_group_id(literature, overlay)
    context_value = " ".join(str(literature.get(key, overlay.get(key, ""))) for key in ("context", "context_terms", "disease", "organ", "tissue", "cell_type"))
    context_tokens = token_set(context_value, overlay.get("cell_state", ""), overlay.get("cell_line", ""))
    overlay_signature_tokens = token_set(
        overlay.get("target_symbols", ""),
        overlay.get("target_uids", ""),
        overlay.get("pathway_terms", ""),
        overlay.get("metabolite_terms", ""),
        overlay.get("context_terms", ""),
    )
    context_holdout = bool(context_value and stable_bucket(context_value, buckets=20) == 0)
    publication_year, min_publication_year, max_publication_year = publication_year_fields(literature, overlay)
    temporal_holdout = bool(publication_year and publication_year >= 2024)
    weak_label_class = five_class_weak_label(weak_label, label_reason)
    group_id = f"{group_key}:{group_value or candidate_uid}"
    existing_edge_count = safe_int(literature.get("supported_existing_edge_count", overlay.get("supported_existing_edge_count", 0)))
    subject_pathways = uid_set(subject.get("pathway_context_uids", ""))
    object_pathways = uid_set(obj.get("pathway_context_uids", ""))
    shared_pathways = subject_pathways & object_pathways
    pathway_union = subject_pathways | object_pathways
    subject_pathway_count = safe_int(subject.get("pathway_context_count", len(subject_pathways)))
    object_pathway_count = safe_int(obj.get("pathway_context_count", len(object_pathways)))
    graph_degree_min = min(safe_float(subject.get("graph_degree", 0)), safe_float(obj.get("graph_degree", 0)))
    graph_degree_max = max(safe_float(subject.get("graph_degree", 0)), safe_float(obj.get("graph_degree", 0)))
    family = relation_family(subject_type, object_type, predicate)
    metabolite_target_bridge = (
        family in {"metabolite_gene", "metabolite_target", "metabolite_protein"}
        and len(shared_pathways) > 0
    )
    endpoint_support_min_log = log1p(
        min(safe_float(subject.get("literature_support_count", 0)), safe_float(obj.get("literature_support_count", 0)))
    )
    graph_degree_min_log = log1p(graph_degree_min)
    shared_pathway_jaccard = float(len(shared_pathways) / len(pathway_union)) if pathway_union else 0.0
    metabolic_theme_signal = float(overlap_count(subject_tokens | overlay_signature_tokens, object_tokens | overlay_signature_tokens, METABOLIC_THEME_TERMS) > 0)
    source_independent_prior = min(
        1.0,
        0.18 * min(1.0, endpoint_support_min_log / 5.0)
        + 0.16 * min(1.0, graph_degree_min_log / 5.0)
        + 0.16 * min(1.0, jaccard(subject_tokens, object_tokens) * 2.5)
        + 0.18 * min(1.0, shared_pathway_jaccard * 5.0)
        + 0.12 * (1.0 if metabolite_target_bridge else 0.0)
        + 0.10 * (1.0 if safe_float(subject.get("prediction_overlay_count", 0)) > 0 or safe_float(obj.get("prediction_overlay_count", 0)) > 0 else 0.0)
        + 0.10 * metabolic_theme_signal,
    )
    return {
        "candidate_uid": candidate_uid,
        "subject_uid": subject_uid,
        "subject_type": subject_type,
        "subject_name": subject.get("display_name", overlay.get("drug_name", "")),
        "predicate": predicate,
        "relation_family": family,
        "subject_object_family": f"{subject_type}:{object_type}",
        "object_uid": object_uid,
        "object_type": object_type,
        "object_name": obj.get("display_name", overlay.get("target_name", "")),
        "source_kind": source_kind,
        "weak_label": weak_label,
        "weak_label_class": weak_label_class,
        "label_score": float(label_score),
        "label_reason": label_reason,
        "group_key": group_key,
        "group_id": group_id,
        "paper_id": group_value if group_key == "paper_id" else "",
        "publication_year": publication_year,
        "min_publication_year": min_publication_year,
        "max_publication_year": max_publication_year,
        "study_id": literature.get("study_id", overlay.get("study_id", "")),
        "dataset_id": literature.get("dataset_id", overlay.get("dataset_id", "")),
        "split": split_for_group(group_id, temporal=temporal_holdout, context_holdout=context_holdout),
        "temporal_holdout": temporal_holdout,
        "context_holdout": context_holdout,
        "p_literature": p_literature,
        "literature_sentence_count": sentence_count,
        "literature_distinct_article_count": distinct_articles,
        "supported_existing_edge_count": existing_edge_count,
        "overlay_confidence": overlay_confidence,
        "subject_graph_log": log1p(subject.get("graph_degree", 0)),
        "object_graph_log": log1p(obj.get("graph_degree", 0)),
        "subject_mention_log": log1p(subject.get("literature_mention_count", 0)),
        "object_mention_log": log1p(obj.get("literature_mention_count", 0)),
        "subject_support_log": log1p(subject.get("literature_support_count", 0)),
        "object_support_log": log1p(obj.get("literature_support_count", 0)),
        "subject_overlay_log": log1p(subject.get("prediction_overlay_count", 0)),
        "object_overlay_log": log1p(obj.get("prediction_overlay_count", 0)),
        "graph_degree_sum_log": log1p(safe_float(subject.get("graph_degree", 0)) + safe_float(obj.get("graph_degree", 0))),
        "graph_degree_min_log": graph_degree_min_log,
        "graph_degree_max_log": log1p(graph_degree_max),
        "graph_degree_balance": float(graph_degree_min / graph_degree_max) if graph_degree_max > 0 else 0.0,
        "canonical_edge_sum_log": log1p(safe_float(subject.get("canonical_edge_count", 0)) + safe_float(obj.get("canonical_edge_count", 0))),
        "direct_graph_support_indicator": 1.0 if existing_edge_count > 0 else 0.0,
        "entity_semantic_overlap": jaccard(subject_tokens, object_tokens),
        "metabolic_theme_overlap": float(overlap_count(subject_tokens | overlay_signature_tokens, object_tokens | overlay_signature_tokens, METABOLIC_THEME_TERMS)),
        "cell_context_overlap": float(overlap_count(context_tokens | overlay_signature_tokens, subject_tokens | object_tokens, CELL_CONTEXT_TERMS)),
        "signature_entity_overlap": float(len(overlay_signature_tokens & (subject_tokens | object_tokens))),
        "subject_pathway_log": log1p(subject_pathway_count),
        "object_pathway_log": log1p(object_pathway_count),
        "shared_pathway_count": float(len(shared_pathways)),
        "shared_pathway_jaccard": shared_pathway_jaccard,
        "metabolite_target_pathway_bridge": 1.0 if metabolite_target_bridge else 0.0,
        "endpoint_overlay_any": 1.0 if safe_float(subject.get("prediction_overlay_count", 0)) > 0 or safe_float(obj.get("prediction_overlay_count", 0)) > 0 else 0.0,
        "endpoint_support_min_log": endpoint_support_min_log,
        "source_independent_prior_score": source_independent_prior,
        "literature_support_score": log1p(sentence_count) * p_literature,
        "independent_article_score": log1p(distinct_articles),
        "curated_edge_score": log1p(existing_edge_count),
        "overlay_score": overlay_confidence,
        "source_ref": literature.get("support_uid", overlay.get("source_record_id", "")),
    }


def build_literature_candidates(
    literature: pd.DataFrame,
    entity_features: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for item in literature.to_dict(orient="records"):
        subject_uid = str(item.get("subject_uid", ""))
        object_uid = str(item.get("object_uid", ""))
        if not subject_uid or not object_uid:
            continue
        weak_label, label_score, reason = relation_label(pd.Series(item))
        candidate_uid = "lit_" + stable_hash("|".join([subject_uid, str(item.get("predicate", "")), object_uid, str(item.get("support_uid", ""))]))[:24]
        rows.append(
            feature_row(
                candidate_uid=candidate_uid,
                subject_uid=subject_uid,
                subject_type=str(item.get("subject_type", "")),
                object_uid=object_uid,
                object_type=str(item.get("object_type", "")),
                predicate=str(item.get("predicate", "")),
                source_kind="literature_relation",
                entity_features=entity_features,
                literature=item,
                weak_label=weak_label,
                label_score=label_score,
                label_reason=reason,
            )
        )
    return rows


def build_canonical_graph_candidates(
    graph_edges: pd.DataFrame,
    entity_features: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if graph_edges.empty:
        return rows
    for item in graph_edges.to_dict(orient="records"):
        subject_uid = str(item.get("subject_uid", "") or "")
        object_uid = str(item.get("object_uid", "") or "")
        edge_uid = str(item.get("edge_uid", "") or "")
        if not subject_uid or not object_uid or not edge_uid:
            continue
        evidence_level = str(item.get("evidence_level", "") or "").casefold()
        if "curated" in evidence_level:
            weak_label, label_score = "medium_positive", 0.70
        else:
            weak_label, label_score = "weak_positive", 0.55
        predicate = str(item.get("edge_type", item.get("predicate", "")) or "")
        row = feature_row(
            candidate_uid="graph_" + stable_hash(edge_uid)[:24],
            subject_uid=subject_uid,
            subject_type=str(item.get("subject_type", "")),
            object_uid=object_uid,
            object_type=str(item.get("object_type", "")),
            predicate=predicate,
            source_kind="canonical_graph_edge",
            entity_features=entity_features,
            overlay={
                "source_record_id": edge_uid,
                "source_name": item.get("source_name", ""),
                "source_version": item.get("source_release", ""),
                "license_id": item.get("license_id", ""),
                "supported_existing_edge_count": 1,
            },
            weak_label=weak_label,
            label_score=label_score,
            label_reason="canonical_graph_edge",
        )
        graph_prior = safe_float(row.get("source_independent_prior_score"))
        row["label_score"] = round(min(0.86, max(0.50, 0.50 + 0.36 * graph_prior + (0.06 if "curated" in evidence_level else 0.0))), 6)
        row["weak_label"] = "medium_positive" if row["label_score"] >= 0.65 else "weak_positive"
        rows.append(row)
    return rows


def build_drug_overlay_candidates(
    drug: pd.DataFrame,
    entity_features: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if drug.empty:
        return rows
    for item in drug.to_dict(orient="records"):
        target_uid = str(item.get("target_uid", "") or "")
        drug_id = str(item.get("drug_id", "") or "")
        if not target_uid or not drug_id:
            continue
        confidence = safe_float(item.get("confidence"))
        if confidence >= 0.75:
            weak_label, label_score, reason = "medium_positive", 0.72, "curated_drug_target_overlay"
        elif confidence >= 0.45:
            weak_label, label_score, reason = "weak_positive", 0.52, "drug_target_overlay"
        else:
            weak_label, label_score, reason = "abstain", 0.0, "low_confidence_overlay"
        candidate_uid = "drug_" + stable_hash("|".join([drug_id, target_uid, str(item.get("source_record_id", ""))]))[:24]
        rows.append(
            feature_row(
                candidate_uid=candidate_uid,
                subject_uid=drug_id,
                subject_type="drug",
                object_uid=target_uid,
                object_type="target",
                predicate="targets",
                source_kind="drug_target_overlay",
                entity_features=entity_features,
                overlay=item,
                weak_label=weak_label,
                label_score=label_score,
                label_reason=reason,
            )
        )
    return rows


def build_negative_candidates(
    positives: list[dict[str, Any]],
    by_type: dict[str, list[str]],
    entity_features: dict[str, dict[str, Any]],
    negatives_per_positive: float,
    seed: int,
) -> list[dict[str, Any]]:
    rng = random.Random(seed)
    rows: list[dict[str, Any]] = []
    carry = 0.0
    for positive in positives:
        if positive["weak_label"] not in {"strong_positive", "medium_positive", "weak_positive"}:
            continue
        carry += negatives_per_positive
        while carry >= 1.0:
            carry -= 1.0
            object_type = str(positive.get("object_type", ""))
            pool = by_type.get(object_type, [])
            if len(pool) < 2:
                break
            object_uid = rng.choice(pool)
            tries = 0
            while object_uid == positive.get("object_uid") and tries < 8:
                object_uid = rng.choice(pool)
                tries += 1
            if object_uid == positive.get("object_uid"):
                break
            candidate_uid = "neg_" + stable_hash("|".join([positive["subject_uid"], positive["predicate"], object_uid, positive["candidate_uid"]]))[:24]
            split_context = {
                "pmids": [positive.get("paper_id", "")] if positive.get("paper_id") else [],
                "publication_year": positive.get("publication_year", ""),
                "min_publication_year": positive.get("min_publication_year", ""),
                "max_publication_year": positive.get("max_publication_year", ""),
                "study_id": positive.get("study_id", ""),
                "dataset_id": positive.get("dataset_id", ""),
            }
            rows.append(
                feature_row(
                    candidate_uid=candidate_uid,
                    subject_uid=positive["subject_uid"],
                    subject_type=positive["subject_type"],
                    object_uid=object_uid,
                    object_type=object_type,
                    predicate=positive["predicate"],
                    source_kind="shuffled_negative",
                    entity_features=entity_features,
                    literature=split_context,
                    weak_label="negative_shuffled",
                    label_score=0.0,
                    label_reason="same_type_random_mismatch",
                )
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build weak labels for prioritization learning.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--negatives-per-positive", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument(
        "--exclude-terms",
        default="",
        help="Comma-separated terms to remove from weak-label candidates before training.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    run_dir = workspace / args.output_root / args.run_id
    views_dir = run_dir / "views"
    output_dir = run_dir / "weak_labels"

    entities = pq.read_table(views_dir / "entities.parquet").to_pandas()
    literature = pq.read_table(views_dir / "literature_relations.parquet").to_pandas()
    drug = pq.read_table(views_dir / "drug_target_overlay.parquet").to_pandas()
    graph_edges_path = views_dir / "canonical_graph_edges.parquet"
    graph_edges = pq.read_table(graph_edges_path).to_pandas() if graph_edges_path.exists() else pd.DataFrame()

    feature_lookup = entity_feature_lookup(entities)
    by_type = entities_by_type(entities)
    literature_rows = build_literature_candidates(literature, feature_lookup)
    graph_rows = build_canonical_graph_candidates(graph_edges, feature_lookup)
    drug_rows = build_drug_overlay_candidates(drug, feature_lookup)
    positive_for_negative = literature_rows + graph_rows + drug_rows
    negative_rows = build_negative_candidates(positive_for_negative, by_type, feature_lookup, args.negatives_per_positive, args.seed)
    all_rows = literature_rows + graph_rows + drug_rows + negative_rows
    labels = pd.DataFrame(all_rows)
    exclude_terms = parse_exclude_terms(args.exclude_terms)
    labels, excluded_rows = filter_excluded_labels(labels, exclude_terms)

    write_parquet(output_dir / "weak_labels.parquet", labels)
    summary = {
        "created_at_utc": utc_now(),
        "run_id": args.run_id,
        "rows": int(len(labels)),
        "excluded_rows": excluded_rows,
        "exclude_term_count": len(exclude_terms),
        "label_counts": labels["weak_label"].value_counts().to_dict() if not labels.empty else {},
        "weak_label_class_counts": labels["weak_label_class"].value_counts().to_dict() if not labels.empty and "weak_label_class" in labels else {},
        "split_counts": labels["split"].value_counts().to_dict() if not labels.empty and "split" in labels else {},
        "grouping": {
            "primary": "PMID/paper_id when available, then study_id, dataset_id, source_record_id, support_uid",
            "default_split": "70% train / 15% validation / 15% test by stable group hash",
            "ood_splits": ["temporal_holdout", "context_holdout"],
        },
        "source_counts": labels["source_kind"].value_counts().to_dict() if not labels.empty else {},
        "outputs": {
            "weak_labels": {"path": str(output_dir / "weak_labels.parquet"), "rows": int(len(labels))},
        },
        "notes": [
            "Weak labels are training signals for prioritization, not truth labels.",
            "Abstain rows should not be used as positive training examples.",
            "Shuffled negatives are pseudo-negatives and can contain rare false negatives.",
        ],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "weak_label_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary["label_counts"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
