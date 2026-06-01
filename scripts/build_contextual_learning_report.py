"""Build a context-aware learning report from MVP learning outputs.

The report combines:

- the current metabolite analysis for a records/context pair
- unsupervised entity embeddings and neighbors
- weak-supervision priority rankings
- prediction overlays and literature evidence references

Outputs are research-prioritization artifacts only; this script does not mutate
canonical release tables.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from metabo_service import MetaboService, load_context_arg, load_records_arg


DEFAULT_RELEASE_ID = "mvp_20260513T002254"
DEFAULT_OUTPUT_ROOT = "learning_runs"
TOKEN_RE = re.compile(r"[A-Za-z][A-Za-z0-9_\-]{1,}")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


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


def normalize_text(value: Any) -> str:
    return " ".join(TOKEN_RE.findall(str(value or "").casefold()))


def context_terms_from_context(context: Any) -> list[str]:
    values: list[str] = []
    if isinstance(context, dict):
        for key in ("cancer_type", "tissue", "cell_type", "cell_state", "cell_line", "context_terms"):
            value = context.get(key)
            if isinstance(value, list):
                values.extend(str(item) for item in value)
            elif value:
                values.append(str(value))
    elif context:
        values.append(str(context))
    terms: list[str] = []
    seen: set[str] = set()
    for value in values:
        normalized = normalize_text(value)
        if normalized and normalized not in seen:
            terms.append(normalized)
            seen.add(normalized)
        for token in normalized.split():
            if len(token) >= 4 and token not in seen:
                terms.append(token)
                seen.add(token)
    return terms


def context_match_score(*values: Any, context_terms: list[str]) -> tuple[float, str]:
    haystack = normalize_text(" ".join(str(value or "") for value in values))
    if not haystack or not context_terms:
        return 0.0, ""
    matched = [term for term in context_terms if term and term in haystack]
    if not matched:
        return 0.0, ""
    phrase_hits = [term for term in matched if " " in term]
    score = min(1.0, 0.2 * len(matched) + 0.35 * len(phrase_hits))
    return score, ";".join(matched[:12])


def evidence_strength(evidence_refs: Any) -> tuple[float, list[str], list[str]]:
    if not isinstance(evidence_refs, list):
        return 0.0, [], []
    max_p = 0.0
    pmids: list[str] = []
    support_classes: list[str] = []
    for ref in evidence_refs:
        if not isinstance(ref, dict):
            continue
        max_p = max(max_p, safe_float(ref.get("p_literature")))
        for pmid in ref.get("pmids") or []:
            text = str(pmid)
            if text and text not in pmids:
                pmids.append(text)
        support_class = str(ref.get("support_class") or "")
        if support_class and support_class not in support_classes:
            support_classes.append(support_class)
    score = min(1.0, max_p * (0.65 + 0.08 * math.log1p(len(pmids))))
    return score, pmids[:20], support_classes[:8]


def service_rank_score(rank: int) -> float:
    return 1.0 / math.log2(rank + 2.0)


def normalize_series(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").fillna(0.0)
    max_value = float(values.max()) if len(values) else 0.0
    if max_value <= 0:
        return values
    return values / max_value


def read_parquet(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pq.read_table(path).to_pandas()


def semicolon(values: list[Any]) -> str:
    return ";".join(str(value) for value in values if str(value))


def build_pathway_context(
    service_output: dict[str, Any],
    learned_pathways: pd.DataFrame,
    context_terms: list[str],
    top_n: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    learned = {
        str(row["pathway_uid"]): row
        for row in learned_pathways.to_dict(orient="records")
        if str(row.get("pathway_uid", ""))
    }
    service_rows = service_output.get("rankings", {}).get("pathways", []) or []
    seen: set[str] = set()
    for idx, row in enumerate(service_rows, start=1):
        uid = str(row.get("pathway_uid", ""))
        if not uid:
            continue
        seen.add(uid)
        learned_row = learned.get(uid, {})
        evidence_score, pmids, support_classes = evidence_strength(row.get("evidence_refs"))
        overlap_score = min(1.0, safe_float(row.get("coverage")) + 0.08 * safe_int(row.get("overlap_count")))
        context_score, matched_context = context_match_score(
            row.get("name"),
            learned_row.get("top_counterpart_name", ""),
            semicolon(pmids),
            context_terms=context_terms,
        )
        learned_score = safe_float(learned_row.get("priority_score"))
        p_value = safe_float(row.get("p_value"), default=1.0)
        p_value_score = min(1.0, -math.log10(max(p_value, 1e-300)) / 8.0) if p_value > 0 else 1.0
        score = (
            0.30 * service_rank_score(idx)
            + 0.18 * overlap_score
            + 0.16 * evidence_score
            + 0.16 * learned_score
            + 0.10 * p_value_score
            + 0.10 * context_score
        )
        rows.append(
            {
                "pathway_uid": uid,
                "pathway_name": row.get("name", ""),
                "contextual_score": round(float(score), 6),
                "analysis_rank": idx,
                "analysis_overlap_count": safe_int(row.get("overlap_count")),
                "analysis_coverage": round(safe_float(row.get("coverage")), 6),
                "analysis_p_value": p_value,
                "learned_priority_score": learned_score,
                "literature_evidence_score": round(evidence_score, 6),
                "matched_context_terms": matched_context,
                "support_classes": semicolon(support_classes),
                "pmids": semicolon(pmids),
                "top_learned_counterpart": learned_row.get("top_counterpart_name", ""),
                "evidence_tier": evidence_tier(score, evidence_score, learned_score),
            }
        )
    for uid, learned_row in learned.items():
        if uid in seen:
            continue
        context_score, matched_context = context_match_score(
            learned_row.get("display_name", ""),
            learned_row.get("top_counterpart_name", ""),
            context_terms=context_terms,
        )
        if context_score <= 0:
            continue
        learned_score = safe_float(learned_row.get("priority_score"))
        score = min(0.34, 0.22 * learned_score + 0.14 * context_score)
        rows.append(
            {
                "pathway_uid": uid,
                "pathway_name": learned_row.get("display_name", ""),
                "contextual_score": round(float(score), 6),
                "analysis_rank": 0,
                "analysis_overlap_count": 0,
                "analysis_coverage": 0.0,
                "analysis_p_value": "",
                "learned_priority_score": learned_score,
                "literature_evidence_score": 0.0,
                "matched_context_terms": matched_context,
                "support_classes": "",
                "pmids": "",
                "top_learned_counterpart": learned_row.get("top_counterpart_name", ""),
                "evidence_tier": "hypothesis",
            }
        )
    frame = pd.DataFrame(rows)
    return frame.sort_values("contextual_score", ascending=False).head(top_n).reset_index(drop=True)


def build_target_context(
    service_output: dict[str, Any],
    learned_targets: pd.DataFrame,
    drug_predictions: list[dict[str, Any]],
    context_terms: list[str],
    top_n: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    learned = {
        str(row["target_uid"]): row
        for row in learned_targets.to_dict(orient="records")
        if str(row.get("target_uid", ""))
    }
    drug_target_hits: dict[str, float] = {}
    drug_names: dict[str, list[str]] = {}
    for drug in drug_predictions:
        for target in drug.get("matched_targets") or []:
            uid = str(target.get("target_uid", ""))
            if not uid:
                continue
            drug_target_hits[uid] = max(drug_target_hits.get(uid, 0.0), safe_float(drug.get("score")))
            drug_names.setdefault(uid, []).append(str(drug.get("drug_name", "")))

    service_rows = service_output.get("rankings", {}).get("targets", []) or []
    max_service_score = max([safe_float(row.get("score")) for row in service_rows] or [0.0])
    seen: set[str] = set()
    for idx, row in enumerate(service_rows, start=1):
        uid = str(row.get("target_uid", ""))
        if not uid:
            continue
        seen.add(uid)
        learned_row = learned.get(uid, {})
        evidence_score, pmids, support_classes = evidence_strength(row.get("evidence_refs"))
        service_score = safe_float(row.get("score")) / max_service_score if max_service_score > 0 else service_rank_score(idx)
        learned_score = safe_float(learned_row.get("priority_score"))
        drug_hit_score = min(1.0, drug_target_hits.get(uid, 0.0))
        context_score, matched_context = context_match_score(
            row.get("display_name"),
            learned_row.get("top_counterpart_name", ""),
            semicolon(support_classes),
            context_terms=context_terms,
        )
        score = (
            0.32 * service_score
            + 0.20 * evidence_score
            + 0.18 * learned_score
            + 0.15 * drug_hit_score
            + 0.10 * service_rank_score(idx)
            + 0.05 * context_score
        )
        rows.append(
            {
                "target_uid": uid,
                "target_name": row.get("display_name", ""),
                "contextual_score": round(float(score), 6),
                "analysis_rank": idx,
                "analysis_score": safe_float(row.get("score")),
                "learned_priority_score": learned_score,
                "literature_support_count": safe_int(row.get("literature_support_count")),
                "max_p_literature": safe_float(row.get("max_p_literature")),
                "literature_evidence_score": round(evidence_score, 6),
                "drug_overlay_score": round(drug_hit_score, 6),
                "matched_drugs": semicolon(sorted(set(drug_names.get(uid, [])))[:10]),
                "matched_context_terms": matched_context,
                "support_classes": semicolon(support_classes),
                "pmids": semicolon(pmids),
                "top_learned_counterpart": learned_row.get("top_counterpart_name", ""),
                "evidence_tier": evidence_tier(score, evidence_score, learned_score),
            }
        )
    for uid, learned_row in learned.items():
        if uid in seen:
            continue
        context_score, matched_context = context_match_score(
            learned_row.get("display_name", ""),
            learned_row.get("top_counterpart_name", ""),
            context_terms=context_terms,
        )
        drug_hit_score = min(1.0, drug_target_hits.get(uid, 0.0))
        if context_score <= 0 and drug_hit_score <= 0:
            continue
        learned_score = safe_float(learned_row.get("priority_score"))
        score = min(0.34, 0.22 * learned_score + 0.18 * drug_hit_score + 0.12 * context_score)
        rows.append(
            {
                "target_uid": uid,
                "target_name": learned_row.get("display_name", ""),
                "contextual_score": round(float(score), 6),
                "analysis_rank": 0,
                "analysis_score": 0.0,
                "learned_priority_score": learned_score,
                "literature_support_count": 0,
                "max_p_literature": 0.0,
                "literature_evidence_score": 0.0,
                "drug_overlay_score": round(drug_hit_score, 6),
                "matched_drugs": semicolon(sorted(set(drug_names.get(uid, [])))[:10]),
                "matched_context_terms": matched_context,
                "support_classes": "",
                "pmids": "",
                "top_learned_counterpart": learned_row.get("top_counterpart_name", ""),
                "evidence_tier": "hypothesis",
            }
        )
    frame = pd.DataFrame(rows)
    return frame.sort_values("contextual_score", ascending=False).head(top_n).reset_index(drop=True)


def build_drug_context(
    service_drugs: list[dict[str, Any]],
    learned_drugs: pd.DataFrame,
    contextual_targets: pd.DataFrame,
    top_n: int,
) -> pd.DataFrame:
    target_scores = {
        str(row["target_uid"]): safe_float(row["contextual_score"])
        for row in contextual_targets.to_dict(orient="records")
    }
    learned = {
        str(row["drug_uid"]): row
        for row in learned_drugs.to_dict(orient="records")
        if str(row.get("drug_uid", ""))
    }
    seen: set[str] = set()
    rows: list[dict[str, Any]] = []
    for idx, row in enumerate(service_drugs, start=1):
        uid = str(row.get("drug_id", ""))
        if not uid:
            continue
        seen.add(uid)
        learned_row = learned.get(uid, {})
        matched_targets = row.get("matched_targets") or []
        target_uids = [str(target.get("target_uid", "")) for target in matched_targets if target.get("target_uid")]
        target_names = [str(target.get("display_name", "")) for target in matched_targets if target.get("display_name")]
        target_context_score = max([target_scores.get(uid, 0.0) for uid in target_uids] or [0.0])
        learned_score = safe_float(learned_row.get("priority_score"))
        overlay_refs = row.get("evidence_refs") or []
        overlay_score = 0.75 if overlay_refs else 0.0
        service_score = safe_float(row.get("score"))
        score = (
            0.40 * min(1.0, service_score)
            + 0.25 * learned_score
            + 0.25 * target_context_score
            + 0.10 * overlay_score
        )
        rows.append(
            {
                "drug_uid": uid,
                "drug_name": row.get("drug_name", ""),
                "contextual_score": round(float(score), 6),
                "prediction_rank": idx,
                "prediction_score": service_score,
                "learned_priority_score": learned_score,
                "target_context_score": round(float(target_context_score), 6),
                "overlay_evidence_score": overlay_score,
                "matched_targets": semicolon(target_names),
                "mechanisms": semicolon(row.get("mechanisms") or []),
                "source_refs": semicolon([ref.get("source_record_id", "") for ref in overlay_refs if isinstance(ref, dict)]),
                "evidence_tier": "medium" if score >= 0.45 and overlay_refs else "hypothesis",
            }
        )
    for uid, learned_row in learned.items():
        if uid in seen:
            continue
        target_uid = str(learned_row.get("top_counterpart_uid", ""))
        target_context_score = target_scores.get(target_uid, 0.0)
        if target_context_score <= 0:
            continue
        learned_score = safe_float(learned_row.get("priority_score"))
        score = 0.35 * learned_score + 0.35 * target_context_score
        rows.append(
            {
                "drug_uid": uid,
                "drug_name": learned_row.get("display_name", ""),
                "contextual_score": round(float(score), 6),
                "prediction_rank": 0,
                "prediction_score": 0.0,
                "learned_priority_score": learned_score,
                "target_context_score": round(float(target_context_score), 6),
                "overlay_evidence_score": 0.0,
                "matched_targets": learned_row.get("top_counterpart_name", ""),
                "mechanisms": "",
                "source_refs": "",
                "evidence_tier": "hypothesis",
            }
        )
    frame = pd.DataFrame(rows)
    return frame.sort_values("contextual_score", ascending=False).head(top_n).reset_index(drop=True)


def evidence_tier(score: float, evidence_score: float, learned_score: float) -> str:
    if score >= 0.75 and evidence_score >= 0.65:
        return "strong"
    if score >= 0.45 and (evidence_score >= 0.35 or learned_score >= 0.5):
        return "medium"
    return "hypothesis"


def matched_metabolite_seeds(service_output: dict[str, Any]) -> list[dict[str, str]]:
    seeds: list[dict[str, str]] = []
    display_by_uid = {
        str(row.get("metabolite_uid", "")): str(row.get("display_name") or row.get("input_name") or row.get("metabolite_uid", ""))
        for row in service_output.get("analysis_features", {}).get("matched", []) or []
        if row.get("metabolite_uid")
    }
    for row in service_output.get("analysis_pack", {}).get("matched", []) or []:
        uid = str(row.get("metabolite_uid", ""))
        if uid:
            candidate = (row.get("candidates") or [{}])[0]
            display_name = display_by_uid.get(uid) or str(candidate.get("display_name") or uid)
            seeds.append({"entity_uid": uid, "display_name": display_name, "seed_type": "matched_metabolite"})
    return seeds


def embedding_neighbors(embedding_frame: pd.DataFrame, seeds: list[dict[str, str]], top_n: int) -> pd.DataFrame:
    if embedding_frame.empty or not seeds:
        return pd.DataFrame()
    dim_cols = [column for column in embedding_frame.columns if column.startswith("dim_")]
    if not dim_cols:
        return pd.DataFrame()
    matrix = embedding_frame[dim_cols].to_numpy(dtype=np.float32)
    uid_to_idx = {str(uid): idx for idx, uid in enumerate(embedding_frame["entity_uid"].astype(str).tolist())}
    rows: list[dict[str, Any]] = []
    for seed in seeds:
        seed_uid = seed["entity_uid"]
        seed_idx = uid_to_idx.get(seed_uid)
        if seed_idx is None:
            continue
        vector = matrix[seed_idx]
        scores = matrix @ vector
        order = np.argpartition(-scores, min(top_n + 1, len(scores) - 1))[: top_n + 1]
        order = order[np.argsort(-scores[order])]
        rank = 0
        for idx in order:
            if int(idx) == seed_idx:
                continue
            rank += 1
            item = embedding_frame.iloc[int(idx)]
            rows.append(
                {
                    "seed_uid": seed_uid,
                    "seed_name": seed["display_name"],
                    "seed_type": seed["seed_type"],
                    "neighbor_uid": item["entity_uid"],
                    "neighbor_type": item["node_type"],
                    "neighbor_name": item["display_name"],
                    "rank": rank,
                    "cosine_similarity": round(float(scores[int(idx)]), 6),
                }
            )
            if rank >= top_n:
                break
    return pd.DataFrame(rows)


def md_table(frame: pd.DataFrame, columns: list[str], max_rows: int = 10) -> str:
    if frame.empty:
        return "_No rows._\n"
    subset = frame[columns].head(max_rows).copy()
    headers = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |"
    rows = []
    for item in subset.to_dict(orient="records"):
        values = [str(item.get(column, "")).replace("\n", " ")[:160] for column in columns]
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join([headers, sep, *rows]) + "\n"


def write_report(
    path: Path,
    run_id: str,
    release_id: str,
    records: list[Any],
    context: Any,
    pathways: pd.DataFrame,
    targets: pd.DataFrame,
    drugs: pd.DataFrame,
    neighbors: pd.DataFrame,
    analysis: dict[str, Any],
) -> None:
    matched_count = len(analysis.get("analysis_pack", {}).get("matched", []) or [])
    ambiguous_count = len(analysis.get("analysis_pack", {}).get("ambiguous", []) or [])
    unmatched_count = len(analysis.get("analysis_pack", {}).get("unmatched", []) or [])
    lines = [
        f"# Contextual Learning Report: {run_id}",
        "",
        f"- Release: `{release_id}`",
        f"- Created: `{utc_now()}`",
        f"- Input records: `{len(records)}`",
        f"- Matched / ambiguous / unmatched: `{matched_count}` / `{ambiguous_count}` / `{unmatched_count}`",
        f"- Context: `{json.dumps(context, ensure_ascii=False)}`",
        "",
        "This report combines current metabolite analysis, unsupervised entity embeddings, weak-supervision priorities, literature evidence, and prediction overlays. Scores are research-prioritization signals, not truth probabilities.",
        "",
        "## Top Pathways",
        md_table(pathways, ["pathway_name", "contextual_score", "analysis_rank", "analysis_overlap_count", "learned_priority_score", "evidence_tier", "support_classes", "top_learned_counterpart"]),
        "## Top Targets",
        md_table(targets, ["target_name", "contextual_score", "analysis_rank", "analysis_score", "learned_priority_score", "drug_overlay_score", "evidence_tier", "matched_drugs"]),
        "## Top Drugs",
        md_table(drugs, ["drug_name", "contextual_score", "prediction_rank", "prediction_score", "learned_priority_score", "target_context_score", "mechanisms", "matched_targets"]),
        "## Similarity Neighbors",
        "",
        "These are exploratory nearest neighbors from the lightweight unsupervised MVP embedding. They are useful for browsing, but noisier than the pathway/target/drug rankings.",
        md_table(neighbors, ["seed_name", "neighbor_name", "neighbor_type", "rank", "cosine_similarity"], max_rows=20),
        "## Notes",
        "- `strong` and `medium` tiers indicate stronger local evidence, not validated biological truth.",
        "- `hypothesis` rows are useful for exploration but should not be promoted to graph facts.",
        "- Context-only Cellosaurus rows are intentionally weak; they help orient the selected context but do not prove cell-line-specific metabolism.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build context-aware learning report from a learning run.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", default=DEFAULT_RELEASE_ID)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--records-json", required=True)
    parser.add_argument("--context-json", default="")
    parser.add_argument("--context-mode", default="soft", choices=["soft", "hard"])
    parser.add_argument("--output-prefix", default="context")
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--top-n", type=int, default=40)
    parser.add_argument("--neighbor-top-n", type=int, default=12)
    parser.add_argument("--max-paths", type=int, default=25)
    parser.add_argument("--max-hops", type=int, default=4)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    run_dir = workspace / args.output_root / args.run_id
    ranking_dir = run_dir / "rankings"
    embedding_dir = run_dir / "embeddings"
    report_dir = run_dir / "reports"

    records = load_records_arg(args.records_json)
    context = load_context_arg(args.context_json)
    context_terms = context_terms_from_context(context)
    service = MetaboService(workspace, release_id=args.release_id)
    analysis = service.analyze_metabolites(
        records,
        max_paths=args.max_paths,
        max_hops=args.max_hops,
        context=context,
        context_mode=args.context_mode,
    )
    predictions = analysis.get("predictions", {})
    drug_predictions = predictions.get("drug_rankings", []) or []

    learned_pathways = read_parquet(ranking_dir / "pathway_priorities.parquet")
    learned_targets = read_parquet(ranking_dir / "target_priorities.parquet")
    learned_drugs = read_parquet(ranking_dir / "drug_priorities.parquet")
    embeddings = read_parquet(embedding_dir / "entity_embeddings.parquet")

    contextual_pathways = build_pathway_context(analysis, learned_pathways, context_terms, args.top_n)
    contextual_targets = build_target_context(analysis, learned_targets, drug_predictions, context_terms, args.top_n)
    contextual_drugs = build_drug_context(drug_predictions, learned_drugs, contextual_targets, args.top_n)
    neighbors = embedding_neighbors(embeddings, matched_metabolite_seeds(analysis), args.neighbor_top_n)

    report_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_prefix
    contextual_pathways.to_csv(report_dir / f"{prefix}_top_pathways.csv", index=False, encoding="utf-8-sig")
    contextual_targets.to_csv(report_dir / f"{prefix}_top_targets.csv", index=False, encoding="utf-8-sig")
    contextual_drugs.to_csv(report_dir / f"{prefix}_top_drugs.csv", index=False, encoding="utf-8-sig")
    neighbors.to_csv(report_dir / f"{prefix}_entity_neighbors.csv", index=False, encoding="utf-8-sig")
    write_report(
        report_dir / f"{prefix}_context_report.md",
        args.run_id,
        args.release_id,
        records,
        context,
        contextual_pathways,
        contextual_targets,
        contextual_drugs,
        neighbors,
        analysis,
    )
    manifest = {
        "created_at_utc": utc_now(),
        "run_id": args.run_id,
        "release_id": args.release_id,
        "records_json": str(Path(args.records_json)),
        "context_json": str(Path(args.context_json)) if args.context_json else "",
        "context_terms": context_terms,
        "outputs": {
            "context_report": str(report_dir / f"{prefix}_context_report.md"),
            "top_pathways": str(report_dir / f"{prefix}_top_pathways.csv"),
            "top_targets": str(report_dir / f"{prefix}_top_targets.csv"),
            "top_drugs": str(report_dir / f"{prefix}_top_drugs.csv"),
            "entity_neighbors": str(report_dir / f"{prefix}_entity_neighbors.csv"),
        },
        "analysis_hash": analysis.get("determinism", {}).get("response_hash", ""),
        "notes": [
            "Contextual scores are research-prioritization signals only.",
            "The script reads learning outputs and frozen release artifacts; it does not write canonical graph facts.",
        ],
    }
    (report_dir / f"{prefix}_context_report_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest["outputs"], ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
