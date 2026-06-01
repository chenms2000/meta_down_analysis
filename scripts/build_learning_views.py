"""Build release-local learning views for discovery and prioritization.

This script creates read-only derived views under learning_runs/<run_id>/views.
It does not update canonical graph or normalized release tables.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


DEFAULT_RELEASE_ID = "mvp_20260513T002254"
DEFAULT_OUTPUT_ROOT = "learning_runs"
TOKEN_RE = re.compile(r"[A-Za-z][A-Za-z0-9_\-]{1,}")
YEAR_RE = re.compile(r"((?:19|20)[0-9]{2})")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def normalize_text(value: Any) -> str:
    return " ".join(TOKEN_RE.findall(str(value or "").lower()))


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(number) or math.isinf(number):
        return default
    return number


def list_len(value: Any) -> int:
    if isinstance(value, (list, tuple, set)):
        return len(value)
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return 0
    return 1


def as_list(value: Any) -> list[str]:
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        value = value.tolist()
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value if str(item or "")]
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    text = str(value)
    return [text] if text else []


def publication_year(value: Any) -> int | None:
    match = YEAR_RE.search(str(value or ""))
    if not match:
        return None
    year = int(match.group(1))
    if 1800 <= year <= 2100:
        return year
    return None


def counter_terms(counter: Counter[str], limit: int = 24) -> str:
    terms: list[str] = []
    for key, count in counter.most_common(limit):
        term = normalize_text(key.replace("_", " "))
        if not term:
            continue
        repeats = max(1, min(3, int(math.log1p(count)) + 1))
        terms.extend([term] * repeats)
    return " ".join(terms)


def read_csv_if_exists(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def write_parquet(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pandas(frame, preserve_index=False)
    pq.write_table(table, path)


def load_nodes(graph_dir: Path) -> pd.DataFrame:
    path = graph_dir / "nodes.parquet"
    columns = ["node_uid", "node_type", "canonical_name", "display_name", "primary_external_id", "source_release", "license_id"]
    return pq.read_table(path, columns=columns).to_pandas()


def stable_bucket(value: str, buckets: int = 100) -> int:
    import hashlib

    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()
    return int(digest[:12], 16) % buckets


def add_pathway_context(pathway_context: dict[str, set[str]], uid: Any, pathway_uid: Any) -> None:
    uid_text = str(uid or "")
    pathway_text = str(pathway_uid or "")
    if uid_text and pathway_text:
        pathway_context[uid_text].add(pathway_text)


def aggregate_graph_edges(
    graph_dir: Path,
    batch_size: int,
) -> tuple[Counter[str], dict[str, Counter[str]], Counter[str], dict[str, set[str]]]:
    path = graph_dir / "edges.parquet"
    degree: Counter[str] = Counter()
    edge_terms: dict[str, Counter[str]] = defaultdict(Counter)
    canonical_edge_count: Counter[str] = Counter()
    pathway_context: dict[str, set[str]] = defaultdict(set)
    target_gene_map: dict[str, set[str]] = defaultdict(set)
    gene_pathway_map: dict[str, set[str]] = defaultdict(set)
    parquet = pq.ParquetFile(path)
    columns = ["subject_uid", "subject_type", "object_uid", "object_type", "edge_type", "predicate", "evidence_level", "source_name"]
    for batch in parquet.iter_batches(columns=columns, batch_size=batch_size):
        frame = batch.to_pandas()
        for uid, count in frame["subject_uid"].value_counts(dropna=True).items():
            degree[str(uid)] += int(count)
            canonical_edge_count[str(uid)] += int(count)
        for uid, count in frame["object_uid"].value_counts(dropna=True).items():
            degree[str(uid)] += int(count)
            canonical_edge_count[str(uid)] += int(count)

        grouped = frame.groupby(["subject_uid", "edge_type"], dropna=True).size()
        for (uid, edge_type), count in grouped.items():
            edge_terms[str(uid)][f"out_{edge_type}"] += int(count)
        grouped = frame.groupby(["object_uid", "edge_type"], dropna=True).size()
        for (uid, edge_type), count in grouped.items():
            edge_terms[str(uid)][f"in_{edge_type}"] += int(count)
        grouped = frame.groupby(["subject_uid", "predicate"], dropna=True).size()
        for (uid, predicate), count in grouped.items():
            edge_terms[str(uid)][f"predicate_{predicate}"] += int(count)

        pathway_edges = frame[frame["object_type"].astype(str) == "pathway"]
        for row in pathway_edges.itertuples(index=False):
            subject_uid = str(getattr(row, "subject_uid", "") or "")
            subject_type = str(getattr(row, "subject_type", "") or "")
            pathway_uid = str(getattr(row, "object_uid", "") or "")
            add_pathway_context(pathway_context, subject_uid, pathway_uid)
            if subject_type == "gene":
                gene_pathway_map[subject_uid].add(pathway_uid)

        target_gene_edges = frame[
            (frame["subject_type"].astype(str) == "target")
            & (frame["object_type"].astype(str) == "gene")
        ]
        for row in target_gene_edges.itertuples(index=False):
            target_gene_map[str(getattr(row, "subject_uid", "") or "")].add(str(getattr(row, "object_uid", "") or ""))

    for target_uid, gene_uids in target_gene_map.items():
        for gene_uid in gene_uids:
            pathway_context[target_uid].update(gene_pathway_map.get(gene_uid, set()))
    return degree, edge_terms, canonical_edge_count, pathway_context


def load_canonical_graph_edges(graph_dir: Path, max_rows_per_edge_type: int, batch_size: int) -> pd.DataFrame:
    path = graph_dir / "edges.parquet"
    if not path.exists() or max_rows_per_edge_type <= 0:
        return pd.DataFrame()
    parquet = pq.ParquetFile(path)
    columns = [
        "edge_uid",
        "edge_type",
        "subject_uid",
        "subject_type",
        "predicate",
        "object_uid",
        "object_type",
        "source_name",
        "source_record_id",
        "evidence_level",
        "source_release",
        "license_id",
    ]
    keep_edge_types = {
        "target_associated_with_disease",
        "gene_involved_in_pathway",
        "metabolite_participates_in_pathway",
        "protein_participates_in_reaction",
        "reaction_in_pathway",
        "pathway_parent_of_pathway",
    }
    parts: list[pd.DataFrame] = []
    for batch in parquet.iter_batches(columns=columns, batch_size=batch_size):
        frame = batch.to_pandas()
        frame = frame[frame["edge_type"].astype(str).isin(keep_edge_types)].copy()
        if frame.empty:
            continue
        frame["bucket"] = frame["edge_uid"].astype(str).map(stable_bucket)
        sampled = (
            frame.sort_values(["edge_type", "bucket", "edge_uid"])
            .groupby("edge_type", group_keys=False)
            .head(max_rows_per_edge_type)
        )
        parts.append(sampled.drop(columns=["bucket"]))
    if not parts:
        return pd.DataFrame(columns=columns)
    combined = pd.concat(parts, ignore_index=True)
    combined["bucket"] = combined["edge_uid"].astype(str).map(stable_bucket)
    combined = (
        combined.sort_values(["edge_type", "bucket", "edge_uid"])
        .groupby("edge_type", group_keys=False)
        .head(max_rows_per_edge_type)
        .drop_duplicates("edge_uid")
        .drop(columns=["bucket"])
        .reset_index(drop=True)
    )
    return combined


def aggregate_sentence_mentions(literature_dir: Path) -> tuple[Counter[str], Counter[str], dict[str, Counter[str]]]:
    path = literature_dir / "sentence_mentions.parquet"
    if not path.exists():
        return Counter(), Counter(), defaultdict(Counter)
    columns = ["entity_uid", "article_uid", "surface", "display_name", "entity_type", "resolution_confidence"]
    mentions = pq.read_table(path, columns=columns).to_pandas()
    mention_count: Counter[str] = Counter(mentions["entity_uid"].astype(str))
    article_count: Counter[str] = Counter(
        mentions[["entity_uid", "article_uid"]].drop_duplicates()["entity_uid"].astype(str)
    )
    mention_terms: dict[str, Counter[str]] = defaultdict(Counter)
    term_frame = mentions[["entity_uid", "surface", "display_name", "entity_type"]].fillna("")
    for row in term_frame.itertuples(index=False):
        uid = str(row.entity_uid)
        for value in (row.surface, row.display_name, row.entity_type):
            text = normalize_text(value)
            if text:
                mention_terms[uid][text] += 1
    return mention_count, article_count, mention_terms


def load_literature_relations(literature_dir: Path) -> tuple[pd.DataFrame, Counter[str], dict[str, Counter[str]]]:
    path = literature_dir / "literature_edge_support.parquet"
    if not path.exists():
        return pd.DataFrame(), Counter(), defaultdict(Counter)
    frame = pq.read_table(path).to_pandas()
    for column in ("p_literature", "raw_score_max", "calibrated_prob_max", "evidence_sentence_count", "distinct_article_count"):
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0.0)
    frame["supported_existing_edge_count"] = frame.get("supported_existing_edge_uids", pd.Series(dtype=object)).apply(list_len)
    frame["evidence_relation_count"] = frame.get("evidence_relation_uids", pd.Series(dtype=object)).apply(list_len)
    frame["sentence_count"] = frame.get("sentence_uids", pd.Series(dtype=object)).apply(list_len)
    frame["pmid_count"] = frame.get("pmids", pd.Series(dtype=object)).apply(list_len)
    frame["pmcid_count"] = frame.get("pmcids", pd.Series(dtype=object)).apply(list_len)

    support_count: Counter[str] = Counter()
    support_terms: dict[str, Counter[str]] = defaultdict(Counter)
    for row in frame.itertuples(index=False):
        subject_uid = str(getattr(row, "subject_uid"))
        object_uid = str(getattr(row, "object_uid"))
        support_count[subject_uid] += 1
        support_count[object_uid] += 1
        support_class = str(getattr(row, "support_class", ""))
        predicate = str(getattr(row, "predicate", ""))
        subject_type = str(getattr(row, "subject_type", ""))
        object_type = str(getattr(row, "object_type", ""))
        term = " ".join(filter(None, [support_class, predicate, subject_type, object_type]))
        if term:
            support_terms[subject_uid][term] += 1
            support_terms[object_uid][term] += 1
    keep_columns = [
        "support_uid",
        "subject_uid",
        "subject_type",
        "predicate",
        "object_uid",
        "object_type",
        "pmids",
        "pmcids",
        "sentence_uids",
        "support_class",
        "evidence_sentence_count",
        "distinct_article_count",
        "p_literature",
        "raw_score_max",
        "calibrated_prob_max",
        "supported_existing_edge_count",
        "evidence_relation_count",
        "sentence_count",
        "pmid_count",
        "pmcid_count",
        "source_release",
        "license_id",
        "config_hash",
    ]
    return frame[[column for column in keep_columns if column in frame.columns]].copy(), support_count, support_terms


def load_article_year_lookup(normalized_dir: Path) -> dict[str, int]:
    path = normalized_dir / "articles.parquet"
    if not path.exists():
        return {}
    columns = ["article_uid", "pmid", "pmcid", "pub_date"]
    available = [column for column in columns if column in pq.ParquetFile(path).schema.names]
    frame = pq.read_table(path, columns=available).to_pandas()
    lookup: dict[str, int] = {}
    if "pub_date" not in frame.columns:
        return lookup
    for row in frame.to_dict(orient="records"):
        year = publication_year(row.get("pub_date"))
        if year is None:
            continue
        for key in ("pmid", "pmcid", "article_uid"):
            value = str(row.get(key) or "")
            if value:
                lookup[value] = year
    return lookup


def add_publication_years(literature: pd.DataFrame, article_year_lookup: dict[str, int]) -> pd.DataFrame:
    if literature.empty:
        return literature
    result = literature.copy()
    years_by_row: list[list[int]] = []
    for row in result.to_dict(orient="records"):
        years = sorted(
            {
                article_year_lookup[token]
                for key in ("pmids", "pmcids", "article_uids")
                for token in as_list(row.get(key))
                if token in article_year_lookup
            }
        )
        years_by_row.append(years)
    result["publication_years"] = years_by_row
    result["min_publication_year"] = [min(years) if years else None for years in years_by_row]
    result["max_publication_year"] = [max(years) if years else None for years in years_by_row]
    result["publication_year"] = result["max_publication_year"]
    result["publication_year_source"] = ["normalized_store_articles" if years else "" for years in years_by_row]
    return result


def overlay_terms_and_counts(overlay_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, Counter[str], dict[str, Counter[str]]]:
    drug = read_csv_if_exists(overlay_dir / "drug_targets.csv")
    cell = read_csv_if_exists(overlay_dir / "cell_type_signatures.csv")
    overlay_count: Counter[str] = Counter()
    overlay_terms: dict[str, Counter[str]] = defaultdict(Counter)

    if not drug.empty:
        drug["confidence"] = pd.to_numeric(drug.get("confidence", 0.0), errors="coerce").fillna(0.0)
        for row in drug.itertuples(index=False):
            uid = str(getattr(row, "target_uid", "") or "")
            if not uid:
                continue
            overlay_count[uid] += 1
            values = [
                getattr(row, "drug_name", ""),
                getattr(row, "mechanism", ""),
                getattr(row, "target_symbol", ""),
                getattr(row, "target_name", ""),
                getattr(row, "context_terms", ""),
                getattr(row, "source_name", ""),
            ]
            for value in values:
                text = normalize_text(value)
                if text:
                    overlay_terms[uid][f"drug_overlay {text}"] += 1

    if not cell.empty:
        cell["confidence"] = pd.to_numeric(cell.get("confidence", 0.0), errors="coerce").fillna(0.0)
        for row in cell.itertuples(index=False):
            for column in ("target_uids", "target_uid"):
                value = str(getattr(row, column, "") or "")
                for uid in re.split(r"[;|,\s]+", value):
                    uid = uid.strip()
                    if not uid:
                        continue
                    overlay_count[uid] += 1
                    overlay_terms[uid]["cell_context_overlay"] += 1
    return drug, cell, overlay_count, overlay_terms


def build_entities(
    nodes: pd.DataFrame,
    degree: Counter[str],
    graph_terms: dict[str, Counter[str]],
    canonical_edge_count: Counter[str],
    mention_count: Counter[str],
    article_count: Counter[str],
    mention_terms: dict[str, Counter[str]],
    support_count: Counter[str],
    support_terms: dict[str, Counter[str]],
    overlay_count: Counter[str],
    overlay_terms: dict[str, Counter[str]],
    pathway_context: dict[str, set[str]],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for node in nodes.itertuples(index=False):
        uid = str(node.node_uid)
        name_bits = " ".join(
            normalize_text(value)
            for value in (node.canonical_name, node.display_name, node.primary_external_id, node.node_type)
            if normalize_text(value)
        )
        doc = " ".join(
            bit
            for bit in (
                name_bits,
                counter_terms(graph_terms.get(uid, Counter())),
                counter_terms(mention_terms.get(uid, Counter())),
                counter_terms(support_terms.get(uid, Counter())),
                counter_terms(overlay_terms.get(uid, Counter())),
            )
            if bit
        )
        rows.append(
            {
                "entity_uid": uid,
                "node_type": node.node_type,
                "canonical_name": node.canonical_name,
                "display_name": node.display_name,
                "primary_external_id": node.primary_external_id,
                "graph_degree": int(degree.get(uid, 0)),
                "canonical_edge_count": int(canonical_edge_count.get(uid, 0)),
                "literature_mention_count": int(mention_count.get(uid, 0)),
                "literature_article_count": int(article_count.get(uid, 0)),
                "literature_support_count": int(support_count.get(uid, 0)),
                "prediction_overlay_count": int(overlay_count.get(uid, 0)),
                "pathway_context_count": int(len(pathway_context.get(uid, set()))),
                "pathway_context_uids": ";".join(sorted(pathway_context.get(uid, set()))[:256]),
                "source_release": node.source_release,
                "license_id": node.license_id,
                "document": doc,
            }
        )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build learning views from frozen metabolism release artifacts.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", default=DEFAULT_RELEASE_ID)
    parser.add_argument("--run-id", default="")
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--edge-batch-size", type=int, default=250_000)
    parser.add_argument("--max-graph-positive-rows-per-edge-type", type=int, default=25_000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    release_id = args.release_id
    run_id = args.run_id or f"learn_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
    run_dir = workspace / args.output_root / run_id
    views_dir = run_dir / "views"

    graph_dir = workspace / "graph_projection" / release_id
    normalized_dir = workspace / "normalized_store" / release_id
    literature_dir = workspace / "literature_evidence" / release_id
    overlay_dir = workspace / "manual_sources" / "prediction_overlays" / release_id

    nodes = load_nodes(graph_dir)
    degree, graph_terms, canonical_edge_count, pathway_context = aggregate_graph_edges(graph_dir, args.edge_batch_size)
    canonical_graph_edges = load_canonical_graph_edges(
        graph_dir,
        max_rows_per_edge_type=args.max_graph_positive_rows_per_edge_type,
        batch_size=args.edge_batch_size,
    )
    mention_count, article_count, mention_terms = aggregate_sentence_mentions(literature_dir)
    literature_relations, support_count, support_terms = load_literature_relations(literature_dir)
    literature_relations = add_publication_years(literature_relations, load_article_year_lookup(normalized_dir))
    drug_overlay, cell_overlay, overlay_count, overlay_terms = overlay_terms_and_counts(overlay_dir)
    entities = build_entities(
        nodes,
        degree,
        graph_terms,
        canonical_edge_count,
        mention_count,
        article_count,
        mention_terms,
        support_count,
        support_terms,
        overlay_count,
        overlay_terms,
        pathway_context,
    )

    write_parquet(views_dir / "entities.parquet", entities)
    write_parquet(views_dir / "literature_relations.parquet", literature_relations)
    write_parquet(views_dir / "drug_target_overlay.parquet", drug_overlay)
    write_parquet(views_dir / "cell_context_overlay.parquet", cell_overlay)
    write_parquet(views_dir / "canonical_graph_edges.parquet", canonical_graph_edges)

    manifest = {
        "created_at_utc": utc_now(),
        "release_id": release_id,
        "run_id": run_id,
        "workspace": str(workspace),
        "views_dir": str(views_dir),
        "inputs": {
            "nodes": str(graph_dir / "nodes.parquet"),
            "edges": str(graph_dir / "edges.parquet"),
            "articles": str(normalized_dir / "articles.parquet"),
            "sentence_mentions": str(literature_dir / "sentence_mentions.parquet"),
            "literature_edge_support": str(literature_dir / "literature_edge_support.parquet"),
            "drug_targets": str(overlay_dir / "drug_targets.csv"),
            "cell_type_signatures": str(overlay_dir / "cell_type_signatures.csv"),
        },
        "outputs": {
            "entities": {"path": str(views_dir / "entities.parquet"), "rows": int(len(entities))},
            "literature_relations": {"path": str(views_dir / "literature_relations.parquet"), "rows": int(len(literature_relations))},
            "drug_target_overlay": {"path": str(views_dir / "drug_target_overlay.parquet"), "rows": int(len(drug_overlay))},
            "cell_context_overlay": {"path": str(views_dir / "cell_context_overlay.parquet"), "rows": int(len(cell_overlay))},
            "canonical_graph_edges": {"path": str(views_dir / "canonical_graph_edges.parquet"), "rows": int(len(canonical_graph_edges))},
        },
        "notes": [
            "Derived learning views only; canonical release tables are not modified.",
            "Prediction overlays remain research-prioritization evidence, not graph facts.",
        ],
    }
    views_dir.mkdir(parents=True, exist_ok=True)
    (views_dir / "learning_view_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest["outputs"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
