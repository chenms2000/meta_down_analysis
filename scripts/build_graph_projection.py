"""Build derived graph projection and resolver indexes from a normalized release.

The normalized store remains the source of truth. This script emits a replayable
derived layer for graph analysis, sparse-matrix style algorithms, and exact
entity resolution.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import sys
import time
import unicodedata
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

try:
    import pyarrow as pa
    import pyarrow.dataset as ds
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover - handled in main for friendlier errors.
    pa = None
    pq = None
    ds = None

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from build_normalized_store import json_dumps, stable_uid  # noqa: E402


DEFAULT_NORMALIZED_ROOT = "normalized_store"
DEFAULT_OUTPUT_ROOT = "graph_projection"
DEFAULT_RELEASE_RE = re.compile(r"mvp_[0-9T]+$")


def require_arrow() -> None:
    if pa is None or pq is None:
        raise RuntimeError("pyarrow is required to build graph_projection parquet outputs.")


def utc_now() -> str:
    return dt.datetime.now(dt.UTC).isoformat()


def safe_float(value: Any, default: float = 1.0) -> float:
    if value is None:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def normalize_lookup_key(value: str) -> str:
    """Normalize exact text lookup keys without doing fuzzy matching."""
    text = unicodedata.normalize("NFKC", str(value or "")).casefold().strip()
    text = re.sub(r"\s+", " ", text)
    return text


def normalize_id_token(value: str) -> str:
    text = normalize_lookup_key(value)
    text = text.replace("_", ":")
    text = re.sub(r"\s+", "", text)
    return text


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def table_path(root: Path, name: str) -> Path:
    path = root / f"{name}.parquet"
    if not path.exists():
        raise FileNotFoundError(f"Missing normalized table: {path}")
    return path


def has_table(root: Path, name: str) -> bool:
    return (root / f"{name}.parquet").exists()


def latest_normalized_release_id(normalized_root: Path) -> str:
    candidates: list[tuple[float, str]] = []
    if not normalized_root.exists():
        raise FileNotFoundError(f"Normalized root does not exist: {normalized_root}")
    for path in normalized_root.iterdir():
        if path.is_dir() and DEFAULT_RELEASE_RE.fullmatch(path.name):
            candidates.append((path.stat().st_mtime, path.name))
    if not candidates:
        raise FileNotFoundError(f"No mvp_* normalized releases found in {normalized_root}")
    return sorted(candidates)[-1][1]


def iter_table_rows(path: Path, columns: list[str], batch_size: int = 50_000) -> Iterator[dict[str, Any]]:
    parquet_file = pq.ParquetFile(path)
    for batch in parquet_file.iter_batches(columns=columns, batch_size=batch_size):
        for row in batch.to_pylist():
            yield row


def schema_for(name: str) -> pa.Schema:
    list_str = pa.list_(pa.string())
    schemas = {
        "nodes": pa.schema(
            [
                ("node_uid", pa.string()),
                ("node_idx", pa.int64()),
                ("node_type", pa.string()),
                ("canonical_name", pa.string()),
                ("display_name", pa.string()),
                ("primary_external_id", pa.string()),
                ("external_xrefs", list_str),
                ("source_priority", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
            ]
        ),
        "node_index": pa.schema(
            [
                ("node_uid", pa.string()),
                ("node_idx", pa.int64()),
                ("node_type", pa.string()),
            ]
        ),
        "edges": pa.schema(
            [
                ("edge_uid", pa.string()),
                ("edge_type", pa.string()),
                ("subject_uid", pa.string()),
                ("subject_type", pa.string()),
                ("predicate", pa.string()),
                ("object_uid", pa.string()),
                ("object_type", pa.string()),
                ("weight", pa.float64()),
                ("source_table", pa.string()),
                ("source_name", pa.string()),
                ("source_record_id", pa.string()),
                ("evidence_level", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
                ("metadata_json", pa.string()),
            ]
        ),
        "sparse_edges": pa.schema(
            [
                ("edge_uid", pa.string()),
                ("source_idx", pa.int64()),
                ("target_idx", pa.int64()),
                ("edge_type_id", pa.int64()),
                ("weight", pa.float64()),
            ]
        ),
        "edge_type_index": pa.schema(
            [
                ("edge_type_id", pa.int64()),
                ("edge_type", pa.string()),
                ("predicate", pa.string()),
                ("subject_type", pa.string()),
                ("object_type", pa.string()),
            ]
        ),
        "resolver_index": pa.schema(
            [
                ("lookup_uid", pa.string()),
                ("entity_uid", pa.string()),
                ("entity_type", pa.string()),
                ("namespace", pa.string()),
                ("lookup_key", pa.string()),
                ("raw_value", pa.string()),
                ("match_field", pa.string()),
                ("rank", pa.float64()),
                ("source_table", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
            ]
        ),
    }
    return schemas[name]


class ParquetWriter:
    def __init__(self, path: Path, schema_name: str):
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            path.unlink()
        self.path = path
        self.schema_name = schema_name
        self.schema = schema_for(schema_name)
        self.writer = pq.ParquetWriter(path, self.schema, compression="snappy")
        self.rows = 0

    def write(self, rows: Iterable[dict[str, Any]]) -> None:
        batch = list(rows)
        if not batch:
            return
        self.writer.write_table(pa.Table.from_pylist(batch, schema=self.schema))
        self.rows += len(batch)

    def close(self) -> dict[str, Any]:
        self.writer.close()
        return {
            "table": self.schema_name,
            "path": str(self.path),
            "rows": self.rows,
            "bytes": self.path.stat().st_size if self.path.exists() else 0,
            "sha256": sha256_file(self.path) if self.path.exists() else "",
        }


@dataclass
class ProjectionContext:
    normalized_dir: Path
    output_dir: Path
    release_id: str
    created_at_utc: str = field(default_factory=utc_now)
    node_index: dict[str, int] = field(default_factory=dict)
    node_type: dict[str, str] = field(default_factory=dict)
    edge_type_ids: dict[str, int] = field(default_factory=dict)
    edge_type_meta: dict[str, tuple[str, str, str]] = field(default_factory=dict)
    skipped_edges: dict[str, int] = field(default_factory=dict)
    table_audit: list[dict[str, Any]] = field(default_factory=list)

    def writer(self, relative_path: str, schema_name: str) -> ParquetWriter:
        return ParquetWriter(self.output_dir / relative_path, schema_name)

    def register_edge_type(self, edge_type: str, predicate: str, subject_type: str, object_type: str) -> int:
        if edge_type not in self.edge_type_ids:
            self.edge_type_ids[edge_type] = len(self.edge_type_ids)
            self.edge_type_meta[edge_type] = (predicate, subject_type, object_type)
        return self.edge_type_ids[edge_type]

    def note_skipped(self, edge_type: str) -> None:
        self.skipped_edges[edge_type] = self.skipped_edges.get(edge_type, 0) + 1


NODE_SPECS = [
    {
        "table": "metabolites",
        "uid": "metabolite_uid",
        "type": "metabolite",
        "name": "canonical_name",
        "display": "canonical_name",
        "primary": None,
        "xrefs": "external_xrefs",
        "priority": "source_priority",
    },
    {
        "table": "genes",
        "uid": "gene_uid",
        "type": "gene",
        "name": "symbol",
        "display": "symbol",
        "primary": "ensembl_gene_id",
        "xrefs": "external_xrefs",
        "priority": "biotype",
    },
    {
        "table": "pathways",
        "uid": "pathway_uid",
        "type": "pathway",
        "name": "name",
        "display": "name",
        "primary": "primary_external_id",
        "xrefs": "external_xrefs",
        "priority": "source_name",
    },
    {
        "table": "reactions",
        "uid": "reaction_uid",
        "type": "reaction",
        "name": "name",
        "display": "name",
        "primary": "primary_external_id",
        "xrefs": "external_xrefs",
        "priority": "source_name",
    },
    {
        "table": "diseases",
        "uid": "disease_uid",
        "type": "disease",
        "name": "name",
        "display": "name",
        "primary": "primary_external_id",
        "xrefs": "external_xrefs",
        "priority": "source_priority",
    },
    {
        "table": "cell_types",
        "uid": "cell_type_uid",
        "type": "cell_type",
        "name": "name",
        "display": "name",
        "primary": "primary_external_id",
        "xrefs": "external_xrefs",
        "priority": "source_priority",
        "optional": True,
    },
    {
        "table": "cell_states",
        "uid": "cell_state_uid",
        "type": "cell_state",
        "name": "name",
        "display": "name",
        "primary": "primary_external_id",
        "xrefs": "external_xrefs",
        "priority": "state_category",
        "optional": True,
    },
    {
        "table": "tissues",
        "uid": "tissue_uid",
        "type": "tissue",
        "name": "name",
        "display": "name",
        "primary": "primary_external_id",
        "xrefs": "external_xrefs",
        "priority": "source_priority",
        "optional": True,
    },
    {
        "table": "targets",
        "uid": "target_uid",
        "type": "target",
        "name": "approved_symbol",
        "display": "preferred_name",
        "primary": "target_external_id",
        "xrefs": "external_xrefs",
        "priority": "target_type",
    },
]


def primary_id(row: dict[str, Any], spec: dict[str, Any]) -> str:
    column = spec.get("primary")
    value = row.get(column) if column else ""
    if value and spec["table"] == "genes" and str(value).startswith("ENSG"):
        return f"ENSEMBL:{value}"
    if value and spec["table"] == "targets" and str(value).startswith("ENSG"):
        return f"OPENTARGETS:{value}"
    if value:
        return str(value)
    xrefs = row.get(spec["xrefs"]) or []
    return str(xrefs[0]) if xrefs else ""


def build_nodes(ctx: ProjectionContext) -> None:
    started = time.time()
    nodes_writer = ctx.writer("nodes.parquet", "nodes")
    index_writer = ctx.writer("node_index.parquet", "node_index")
    rows: list[dict[str, Any]] = []
    index_rows: list[dict[str, Any]] = []
    node_idx = 0

    def emit_node(node_row: dict[str, Any]) -> None:
        nonlocal node_idx, rows, index_rows
        node_uid = node_row["node_uid"]
        if not node_uid or node_uid in ctx.node_index:
            return
        node_row["node_idx"] = node_idx
        ctx.node_index[node_uid] = node_idx
        ctx.node_type[node_uid] = node_row["node_type"]
        rows.append(node_row)
        index_rows.append({"node_uid": node_uid, "node_idx": node_idx, "node_type": node_row["node_type"]})
        node_idx += 1
        if len(rows) >= 100_000:
            nodes_writer.write(rows)
            index_writer.write(index_rows)
            rows = []
            index_rows = []

    for spec in NODE_SPECS:
        if spec.get("optional") and not has_table(ctx.normalized_dir, spec["table"]):
            continue
        columns = [
            spec["uid"],
            spec["name"],
            spec["display"],
            spec["xrefs"],
            spec["priority"],
            "source_release",
            "license_id",
        ]
        if spec.get("primary"):
            columns.append(spec["primary"])
        for row in iter_table_rows(table_path(ctx.normalized_dir, spec["table"]), columns):
            node_uid = row.get(spec["uid"], "")
            if not node_uid:
                continue
            name = row.get(spec["name"], "") or primary_id(row, spec)
            display = row.get(spec["display"], "") or name
            emit_node(
                {
                "node_uid": node_uid,
                "node_type": spec["type"],
                "canonical_name": str(name or ""),
                "display_name": str(display or ""),
                "primary_external_id": primary_id(row, spec),
                "external_xrefs": sorted(set(str(x) for x in (row.get(spec["xrefs"]) or []) if x)),
                "source_priority": str(row.get(spec["priority"], "") or ""),
                "source_release": str(row.get("source_release", "") or ""),
                "license_id": str(row.get("license_id", "") or ""),
                "node_idx": -1,
                }
            )

    for row in iter_table_rows(
        table_path(ctx.normalized_dir, "reaction_participants"),
        ["participant_external_id", "participant_type", "source_release", "license_id"],
    ):
        participant_type = str(row.get("participant_type", "") or "").lower()
        external_id = str(row.get("participant_external_id", "") or "").strip()
        if participant_type != "protein" or not external_id:
            continue
        node_uid = stable_uid("protein", external_id)
        emit_node(
            {
                "node_uid": node_uid,
                "node_idx": -1,
                "node_type": "protein",
                "canonical_name": external_id,
                "display_name": external_id,
                "primary_external_id": external_id,
                "external_xrefs": [external_id],
                "source_priority": "derived_reactome_participant",
                "source_release": str(row.get("source_release", "") or ""),
                "license_id": str(row.get("license_id", "") or ""),
            }
        )

    for row in iter_table_rows(
        table_path(ctx.normalized_dir, "target_disease_edges"),
        ["disease_uid", "disease_external_id", "source_release", "license_id"],
    ):
        disease_uid = row.get("disease_uid", "")
        disease_external_id = str(row.get("disease_external_id", "") or "").strip()
        if not disease_uid or disease_uid in ctx.node_index or not disease_external_id:
            continue
        emit_node(
            {
                "node_uid": disease_uid,
                "node_idx": -1,
                "node_type": "disease",
                "canonical_name": disease_external_id,
                "display_name": disease_external_id,
                "primary_external_id": disease_external_id,
                "external_xrefs": [disease_external_id],
                "source_priority": "derived_opentargets_association_endpoint",
                "source_release": str(row.get("source_release", "") or ""),
                "license_id": str(row.get("license_id", "") or ""),
            }
        )
    nodes_writer.write(rows)
    index_writer.write(index_rows)
    audit = nodes_writer.close()
    audit["seconds"] = round(time.time() - started, 3)
    ctx.table_audit.append(audit)
    audit = index_writer.close()
    audit["seconds"] = round(time.time() - started, 3)
    ctx.table_audit.append(audit)


@dataclass(frozen=True)
class EdgeSpec:
    table: str
    subject_col: str
    object_col: str
    subject_type: str
    object_type: str
    edge_type: str
    predicate: str
    predicate_col: str = ""
    weight_col: str = ""
    source_record_col: str = "source_record_id"
    evidence_col: str = "evidence_level"
    metadata_cols: tuple[str, ...] = ()
    optional: bool = False

    def columns(self) -> list[str]:
        cols = {
            "edge_uid",
            self.subject_col,
            self.object_col,
            "source_name",
            "source_release",
            "license_id",
            "parser_hash",
        }
        if self.predicate_col:
            cols.add(self.predicate_col)
        if self.weight_col:
            cols.add(self.weight_col)
        if self.source_record_col:
            cols.add(self.source_record_col)
        if self.evidence_col:
            cols.add(self.evidence_col)
        cols.update(self.metadata_cols)
        return sorted(cols)


EDGE_SPECS = [
    EdgeSpec(
        table="metabolite_pathway_edges",
        subject_col="metabolite_uid",
        object_col="pathway_uid",
        subject_type="metabolite",
        object_type="pathway",
        edge_type="metabolite_participates_in_pathway",
        predicate="participates_in",
        predicate_col="predicate",
        metadata_cols=("metabolite_external_id", "pathway_external_id", "evidence_code", "species"),
    ),
    EdgeSpec(
        table="metabolite_identity_edges",
        subject_col="subject_metabolite_uid",
        object_col="object_metabolite_uid",
        subject_type="metabolite",
        object_type="metabolite",
        edge_type="metabolite_same_as_metabolite",
        predicate="same_as",
        predicate_col="predicate",
        source_record_col="source_record_id",
        evidence_col="evidence_level",
        metadata_cols=("subject_external_id", "object_external_id", "match_basis", "linked_name"),
        optional=True,
    ),
    EdgeSpec(
        table="gene_pathway_edges",
        subject_col="gene_uid",
        object_col="pathway_uid",
        subject_type="gene",
        object_type="pathway",
        edge_type="gene_involved_in_pathway",
        predicate="involved_in",
        predicate_col="predicate",
        metadata_cols=("gene_external_id", "pathway_external_id", "evidence_code", "species"),
    ),
    EdgeSpec(
        table="target_gene_edges",
        subject_col="target_uid",
        object_col="gene_uid",
        subject_type="target",
        object_type="gene",
        edge_type="target_maps_to_gene",
        predicate="targets",
        predicate_col="predicate",
        source_record_col="",
        evidence_col="",
        metadata_cols=("target_external_id", "gene_external_id"),
    ),
    EdgeSpec(
        table="target_disease_edges",
        subject_col="target_uid",
        object_col="disease_uid",
        subject_type="target",
        object_type="disease",
        edge_type="target_associated_with_disease",
        predicate="associated_with",
        predicate_col="predicate",
        weight_col="association_score",
        metadata_cols=(
            "target_external_id",
            "disease_external_id",
            "evidence_count",
            "aggregation_type",
            "aggregation_value",
            "current_novelty",
            "score_components_json",
        ),
    ),
    EdgeSpec(
        table="disease_parent_edges",
        subject_col="child_disease_uid",
        object_col="parent_disease_uid",
        subject_type="disease",
        object_type="disease",
        edge_type="disease_is_a_disease",
        predicate="is_a",
        predicate_col="predicate",
        source_record_col="",
        evidence_col="",
        metadata_cols=("child_external_id", "parent_external_id"),
    ),
    EdgeSpec(
        table="cell_type_parent_edges",
        subject_col="child_cell_type_uid",
        object_col="parent_cell_type_uid",
        subject_type="cell_type",
        object_type="cell_type",
        edge_type="cell_type_is_a_cell_type",
        predicate="is_a",
        predicate_col="predicate",
        source_record_col="",
        evidence_col="",
        metadata_cols=("child_external_id", "parent_external_id"),
        optional=True,
    ),
    EdgeSpec(
        table="cell_state_parent_edges",
        subject_col="child_cell_state_uid",
        object_col="parent_cell_state_uid",
        subject_type="cell_state",
        object_type="cell_state",
        edge_type="cell_state_is_a_cell_state",
        predicate="is_a",
        predicate_col="predicate",
        source_record_col="",
        evidence_col="",
        metadata_cols=("child_external_id", "parent_external_id"),
        optional=True,
    ),
    EdgeSpec(
        table="tissue_parent_edges",
        subject_col="child_tissue_uid",
        object_col="parent_tissue_uid",
        subject_type="tissue",
        object_type="tissue",
        edge_type="tissue_is_a_tissue",
        predicate="is_a",
        predicate_col="predicate",
        source_record_col="",
        evidence_col="",
        metadata_cols=("child_external_id", "parent_external_id"),
        optional=True,
    ),
    EdgeSpec(
        table="pathway_hierarchy_edges",
        subject_col="parent_pathway_uid",
        object_col="child_pathway_uid",
        subject_type="pathway",
        object_type="pathway",
        edge_type="pathway_parent_of_pathway",
        predicate="parent_of",
        predicate_col="predicate",
        source_record_col="",
        evidence_col="",
        metadata_cols=("parent_external_id", "child_external_id"),
    ),
    EdgeSpec(
        table="reactions",
        subject_col="reaction_uid",
        object_col="pathway_uid",
        subject_type="reaction",
        object_type="pathway",
        edge_type="reaction_in_pathway",
        predicate="in_pathway",
        source_record_col="primary_external_id",
        evidence_col="",
        metadata_cols=("primary_external_id", "pathway_external_id", "species"),
    ),
]


def build_edge_row(ctx: ProjectionContext, spec: EdgeSpec, row: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]] | None:
    subject_uid = row.get(spec.subject_col, "")
    object_uid = row.get(spec.object_col, "")
    if not subject_uid or not object_uid or subject_uid not in ctx.node_index or object_uid not in ctx.node_index:
        ctx.note_skipped(spec.edge_type)
        return None
    predicate = row.get(spec.predicate_col, "") if spec.predicate_col else ""
    predicate = predicate or spec.predicate
    weight = safe_float(row.get(spec.weight_col), 1.0) if spec.weight_col else 1.0
    edge_type_id = ctx.register_edge_type(spec.edge_type, predicate, spec.subject_type, spec.object_type)
    metadata = {key: row.get(key) for key in spec.metadata_cols if key in row}
    edge_uid = row.get("edge_uid", "") or stable_uid("edge", spec.edge_type, subject_uid, object_uid, json_dumps(metadata))
    edge = {
        "edge_uid": edge_uid,
        "edge_type": spec.edge_type,
        "subject_uid": subject_uid,
        "subject_type": spec.subject_type,
        "predicate": predicate,
        "object_uid": object_uid,
        "object_type": spec.object_type,
        "weight": weight,
        "source_table": spec.table,
        "source_name": str(row.get("source_name", "") or ""),
        "source_record_id": str(row.get(spec.source_record_col, "") or "") if spec.source_record_col else "",
        "evidence_level": str(row.get(spec.evidence_col, "") or "") if spec.evidence_col else "",
        "source_release": str(row.get("source_release", "") or ""),
        "license_id": str(row.get("license_id", "") or ""),
        "parser_hash": str(row.get("parser_hash", "") or ""),
        "metadata_json": json_dumps(metadata),
    }
    sparse = {
        "edge_uid": edge_uid,
        "source_idx": ctx.node_index[subject_uid],
        "target_idx": ctx.node_index[object_uid],
        "edge_type_id": edge_type_id,
        "weight": weight,
    }
    return edge, sparse


def reaction_participant_edges(ctx: ProjectionContext) -> Iterator[tuple[dict[str, Any], dict[str, Any]]]:
    path = table_path(ctx.normalized_dir, "reaction_participants")
    columns = [
        "edge_uid",
        "reaction_uid",
        "participant_uid",
        "participant_external_id",
        "participant_type",
        "role",
        "source_name",
        "source_record_id",
        "source_release",
        "license_id",
        "parser_hash",
    ]
    for row in iter_table_rows(path, columns):
        participant_type = str(row.get("participant_type", "") or "entity").lower()
        if not row.get("participant_uid") and participant_type == "protein" and row.get("participant_external_id"):
            row = {**row, "participant_uid": stable_uid("protein", row["participant_external_id"])}
        subject_type = participant_type if participant_type in {"metabolite", "gene", "target", "disease", "pathway", "protein"} else "entity"
        spec = EdgeSpec(
            table="reaction_participants",
            subject_col="participant_uid",
            object_col="reaction_uid",
            subject_type=subject_type,
            object_type="reaction",
            edge_type=f"{subject_type}_participates_in_reaction",
            predicate="participates_in_reaction",
            metadata_cols=("participant_external_id", "participant_type", "role"),
        )
        edge = build_edge_row(ctx, spec, row)
        if edge:
            yield edge


def build_edges(ctx: ProjectionContext) -> None:
    started = time.time()
    edge_writer = ctx.writer("edges.parquet", "edges")
    sparse_writer = ctx.writer("sparse_edges.parquet", "sparse_edges")
    edge_rows: list[dict[str, Any]] = []
    sparse_rows: list[dict[str, Any]] = []

    def write_edge_pair(pair: tuple[dict[str, Any], dict[str, Any]]) -> None:
        nonlocal edge_rows, sparse_rows
        edge, sparse = pair
        edge_rows.append(edge)
        sparse_rows.append(sparse)
        if len(edge_rows) >= 100_000:
            edge_writer.write(edge_rows)
            sparse_writer.write(sparse_rows)
            edge_rows = []
            sparse_rows = []

    for spec in EDGE_SPECS:
        if spec.optional and not has_table(ctx.normalized_dir, spec.table):
            continue
        for row in iter_table_rows(table_path(ctx.normalized_dir, spec.table), spec.columns()):
            pair = build_edge_row(ctx, spec, row)
            if pair:
                write_edge_pair(pair)
    for pair in reaction_participant_edges(ctx):
        write_edge_pair(pair)

    edge_writer.write(edge_rows)
    sparse_writer.write(sparse_rows)
    audit = edge_writer.close()
    audit["seconds"] = round(time.time() - started, 3)
    ctx.table_audit.append(audit)
    audit = sparse_writer.close()
    audit["seconds"] = round(time.time() - started, 3)
    ctx.table_audit.append(audit)

    type_writer = ctx.writer("edge_type_index.parquet", "edge_type_index")
    type_rows = []
    for edge_type, edge_type_id in sorted(ctx.edge_type_ids.items(), key=lambda item: item[1]):
        predicate, subject_type, object_type = ctx.edge_type_meta[edge_type]
        type_rows.append(
            {
                "edge_type_id": edge_type_id,
                "edge_type": edge_type,
                "predicate": predicate,
                "subject_type": subject_type,
                "object_type": object_type,
            }
        )
    type_writer.write(type_rows)
    audit = type_writer.close()
    audit["seconds"] = round(time.time() - started, 3)
    ctx.table_audit.append(audit)


def lookup_row(
    entity_uid: str,
    entity_type: str,
    namespace: str,
    raw_value: str,
    match_field: str,
    rank: float,
    source_table: str,
    source_release: str,
    license_id: str,
    key: str | None = None,
) -> dict[str, Any] | None:
    raw = str(raw_value or "").strip()
    lookup_key = key if key is not None else normalize_lookup_key(raw)
    if not raw or not lookup_key:
        return None
    namespace = namespace.upper()
    return {
        "lookup_uid": stable_uid("lookup", entity_uid, entity_type, namespace, lookup_key, match_field),
        "entity_uid": entity_uid,
        "entity_type": entity_type,
        "namespace": namespace,
        "lookup_key": lookup_key,
        "raw_value": raw,
        "match_field": match_field,
        "rank": float(rank),
        "source_table": source_table,
        "source_release": str(source_release or ""),
        "license_id": str(license_id or ""),
    }


def xref_lookup_rows(
    entity_uid: str,
    entity_type: str,
    xref_source: str,
    xref_id: str,
    xref_key: str,
    source_table: str,
    source_release: str,
    license_id: str,
    rank: float = 95.0,
) -> Iterator[dict[str, Any]]:
    source = str(xref_source or "").upper()
    xid = str(xref_id or "").strip()
    xkey = str(xref_key or "").strip()
    if not source or not xid:
        return
    variants: set[tuple[str, str, str]] = {
        (source, xid, normalize_id_token(xid)),
        ("ID", xkey or f"{source}:{xid}", normalize_id_token(xkey or f"{source}:{xid}")),
    }
    if source == "PUBCHEM.COMPOUND":
        variants.add(("PUBCHEM", xid, normalize_id_token(xid)))
        variants.add(("CID", xid, normalize_id_token(xid)))
    if source in {"NCBI.GENE", "ENTREZ"}:
        variants.add(("NCBI", xid, normalize_id_token(xid)))
        variants.add(("ENTREZ", xid, normalize_id_token(xid)))
    if source == "ENSEMBL":
        variants.add(("ENSEMBL_GENE", xid, normalize_id_token(xid)))
    if source in {"MONDO", "DOID", "MESH", "EFO", "HP", "ORPHANET", "OTAR", "NCIT", "ONCOTREE"}:
        variants.add(("OPENTARGETS_DISEASE", xid, normalize_id_token(xid)))
    for namespace, raw, key in sorted(variants):
        row = lookup_row(entity_uid, entity_type, namespace, raw, "xref", rank, source_table, source_release, license_id, key=key)
        if row:
            yield row


def flush_lookup(writer: ParquetWriter, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if len(rows) >= 100_000:
        writer.write(rows)
        return []
    return rows


def write_node_text_resolver_rows(ctx: ProjectionContext, writer: ParquetWriter) -> None:
    rows: list[dict[str, Any]] = []
    node_tables = {
        "metabolites": ("metabolite_uid", "metabolite", "canonical_name", "synonyms"),
        "genes": ("gene_uid", "gene", "symbol", "aliases"),
        "pathways": ("pathway_uid", "pathway", "name", ""),
        "reactions": ("reaction_uid", "reaction", "name", ""),
        "diseases": ("disease_uid", "disease", "name", "aliases"),
        "cell_types": ("cell_type_uid", "cell_type", "name", "aliases"),
        "cell_states": ("cell_state_uid", "cell_state", "name", "aliases"),
        "tissues": ("tissue_uid", "tissue", "name", "aliases"),
        "targets": ("target_uid", "target", "approved_symbol", ""),
    }
    for table, (uid_col, entity_type, name_col, synonyms_col) in node_tables.items():
        if not has_table(ctx.normalized_dir, table):
            continue
        columns = [uid_col, name_col, "source_release", "license_id"]
        if synonyms_col:
            columns.append(synonyms_col)
        if table == "metabolites":
            columns.extend(["inchikey", "external_xrefs"])
        if table == "genes":
            columns.extend(["ensembl_gene_id", "entrez_gene_id", "uniprot_ids", "external_xrefs"])
        if table in {"diseases", "cell_types", "cell_states", "tissues"}:
            columns.extend(["primary_external_id", "external_xrefs"])
        if table == "targets":
            columns.extend(["target_external_id", "preferred_name", "external_xrefs"])
        if table in {"pathways", "reactions"}:
            columns.extend(["primary_external_id", "external_xrefs"])
        for row in iter_table_rows(table_path(ctx.normalized_dir, table), columns):
            entity_uid = row.get(uid_col, "")
            source_release = row.get("source_release", "")
            license_id = row.get("license_id", "")
            seen: set[tuple[str, str, str]] = set()

            def add(namespace: str, raw: str, match_field: str, rank: float, key: str | None = None) -> None:
                out = lookup_row(entity_uid, entity_type, namespace, raw, match_field, rank, table, source_release, license_id, key)
                if out and (out["namespace"], out["lookup_key"], out["match_field"]) not in seen:
                    seen.add((out["namespace"], out["lookup_key"], out["match_field"]))
                    rows.append(out)

            add("TEXT", row.get(name_col, ""), name_col, 100.0)
            if table == "targets":
                add("TEXT", row.get("preferred_name", ""), "preferred_name", 92.0)
            for synonym in row.get(synonyms_col) or []:
                add("TEXT", synonym, "synonym", 80.0)
            if table == "metabolites":
                add("INCHIKEY", row.get("inchikey", ""), "inchikey", 100.0, key=normalize_id_token(row.get("inchikey", "")))
            for xref in row.get("external_xrefs") or []:
                if ":" not in str(xref):
                    continue
                source, xid = str(xref).split(":", 1)
                for out in xref_lookup_rows(entity_uid, entity_type, source, xid, str(xref), table, source_release, license_id):
                    if (out["namespace"], out["lookup_key"], out["match_field"]) not in seen:
                        seen.add((out["namespace"], out["lookup_key"], out["match_field"]))
                        rows.append(out)
            if table == "genes":
                add("ENSEMBL", row.get("ensembl_gene_id", ""), "ensembl_gene_id", 100.0, key=normalize_id_token(row.get("ensembl_gene_id", "")))
                add("NCBI", row.get("entrez_gene_id", ""), "entrez_gene_id", 98.0, key=normalize_id_token(row.get("entrez_gene_id", "")))
                for uniprot in row.get("uniprot_ids") or []:
                    add("UNIPROT", uniprot, "uniprot_id", 92.0, key=normalize_id_token(uniprot))
            if table == "targets":
                add("OPENTARGETS", row.get("target_external_id", ""), "target_external_id", 100.0, key=normalize_id_token(row.get("target_external_id", "")))
            if table in {"pathways", "reactions", "diseases", "cell_types", "cell_states", "tissues"}:
                primary = row.get("primary_external_id", "")
                if primary and ":" in str(primary):
                    source, xid = str(primary).split(":", 1)
                    for out in xref_lookup_rows(entity_uid, entity_type, source, xid, primary, table, source_release, license_id, rank=100.0):
                        if (out["namespace"], out["lookup_key"], out["match_field"]) not in seen:
                            seen.add((out["namespace"], out["lookup_key"], out["match_field"]))
                            rows.append(out)
            rows = flush_lookup(writer, rows)
    writer.write(rows)


def write_xref_table_resolver_rows(ctx: ProjectionContext, writer: ParquetWriter) -> None:
    rows: list[dict[str, Any]] = []
    xref_tables = {
        "metabolite_xrefs": ("metabolite_uid", "metabolite"),
        "gene_xrefs": ("gene_uid", "gene"),
        "disease_xrefs": ("disease_uid", "disease"),
        "cell_type_xrefs": ("cell_type_uid", "cell_type"),
        "cell_state_xrefs": ("cell_state_uid", "cell_state"),
        "tissue_xrefs": ("tissue_uid", "tissue"),
        "target_xrefs": ("target_uid", "target"),
    }
    for table, (uid_col, entity_type) in xref_tables.items():
        if not has_table(ctx.normalized_dir, table):
            continue
        columns = [uid_col, "xref_source", "xref_id", "xref_key", "source_release", "license_id"]
        for row in iter_table_rows(table_path(ctx.normalized_dir, table), columns):
            for out in xref_lookup_rows(
                row.get(uid_col, ""),
                entity_type,
                row.get("xref_source", ""),
                row.get("xref_id", ""),
                row.get("xref_key", ""),
                table,
                row.get("source_release", ""),
                row.get("license_id", ""),
            ):
                rows.append(out)
            rows = flush_lookup(writer, rows)
    writer.write(rows)


def write_edge_endpoint_resolver_rows(ctx: ProjectionContext, writer: ParquetWriter) -> None:
    rows: list[dict[str, Any]] = []
    for row in iter_table_rows(
        table_path(ctx.normalized_dir, "target_disease_edges"),
        ["disease_uid", "disease_external_id", "source_release", "license_id"],
    ):
        disease_uid = row.get("disease_uid", "")
        external_id = str(row.get("disease_external_id", "") or "").strip()
        if not disease_uid or not external_id:
            continue
        if ":" in external_id:
            source, xid = external_id.split(":", 1)
            for out in xref_lookup_rows(
                disease_uid,
                "disease",
                source,
                xid,
                external_id,
                "target_disease_edges",
                row.get("source_release", ""),
                row.get("license_id", ""),
                rank=90.0,
            ):
                rows.append(out)
        else:
            out = lookup_row(
                disease_uid,
                "disease",
                "OPENTARGETS_DISEASE",
                external_id,
                "disease_external_id",
                90.0,
                "target_disease_edges",
                row.get("source_release", ""),
                row.get("license_id", ""),
                key=normalize_id_token(external_id),
            )
            if out:
                rows.append(out)
        rows = flush_lookup(writer, rows)

    for row in iter_table_rows(
        table_path(ctx.normalized_dir, "reaction_participants"),
        ["participant_external_id", "participant_type", "source_release", "license_id"],
    ):
        if str(row.get("participant_type", "") or "").lower() != "protein":
            continue
        external_id = str(row.get("participant_external_id", "") or "").strip()
        if not external_id:
            continue
        protein_uid = stable_uid("protein", external_id)
        if ":" in external_id:
            source, xid = external_id.split(":", 1)
        else:
            source, xid = "UNIPROT", external_id
        for out in xref_lookup_rows(
            protein_uid,
            "protein",
            source,
            xid,
            external_id,
            "reaction_participants",
            row.get("source_release", ""),
            row.get("license_id", ""),
            rank=85.0,
        ):
            rows.append(out)
        rows = flush_lookup(writer, rows)
    writer.write(rows)


def build_resolver_index(ctx: ProjectionContext) -> None:
    started = time.time()
    writer = ctx.writer("resolver_index.parquet", "resolver_index")
    write_node_text_resolver_rows(ctx, writer)
    write_xref_table_resolver_rows(ctx, writer)
    write_edge_endpoint_resolver_rows(ctx, writer)
    audit = writer.close()
    audit["seconds"] = round(time.time() - started, 3)
    ctx.table_audit.append(audit)


def write_manifest(ctx: ProjectionContext) -> Path:
    manifest = {
        "release_id": ctx.release_id,
        "created_at_utc": ctx.created_at_utc,
        "normalized_dir": str(ctx.normalized_dir),
        "output_dir": str(ctx.output_dir),
        "tables": ctx.table_audit,
        "node_count": len(ctx.node_index),
        "edge_type_count": len(ctx.edge_type_ids),
        "skipped_edges": ctx.skipped_edges,
        "notes": [
            "Graph projection is derived from normalized parquet tables; normalized_store remains the source of truth.",
            "Article and sentence evidence remain in normalized_store until mention/relation edges are materialized.",
        ],
    }
    path = ctx.output_dir / "graph_manifest.json"
    path.write_text(json_dumps(manifest) + "\n", encoding="utf-8")
    return path


def build_graph_projection(normalized_root: Path, output_root: Path, release_id: str | None = None) -> Path:
    require_arrow()
    resolved_release = release_id or latest_normalized_release_id(normalized_root)
    normalized_dir = normalized_root / resolved_release
    if not normalized_dir.exists():
        raise FileNotFoundError(f"Normalized release does not exist: {normalized_dir}")
    output_dir = output_root / resolved_release
    output_dir.mkdir(parents=True, exist_ok=True)
    ctx = ProjectionContext(normalized_dir=normalized_dir, output_dir=output_dir, release_id=resolved_release)
    build_nodes(ctx)
    build_edges(ctx)
    build_resolver_index(ctx)
    return write_manifest(ctx)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build graph projection and resolver indexes from normalized_store.")
    parser.add_argument("--normalized-root", default=DEFAULT_NORMALIZED_ROOT)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--release-id", default="")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    workspace = Path.cwd()
    normalized_root = (workspace / args.normalized_root).resolve()
    output_root = (workspace / args.output_root).resolve()
    manifest = build_graph_projection(normalized_root, output_root, args.release_id or None)
    print(f"[graph_projection] done: {manifest}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
