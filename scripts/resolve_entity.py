"""Deterministic local resolver over graph_projection resolver_index.parquet."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

try:
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.dataset as ds
except Exception:  # pragma: no cover - handled in main.
    pa = None
    pc = None
    ds = None

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from build_graph_projection import (  # noqa: E402
    latest_normalized_release_id,
    normalize_id_token,
    normalize_lookup_key,
)


ID_PREFIX_ALIASES = {
    "CHEBI": ["CHEBI"],
    "HMDB": ["HMDB"],
    "KEGG": ["KEGG"],
    "CID": ["CID", "PUBCHEM", "PUBCHEM.COMPOUND"],
    "PUBCHEM": ["PUBCHEM", "PUBCHEM.COMPOUND", "CID"],
    "PUBCHEM.COMPOUND": ["PUBCHEM.COMPOUND", "PUBCHEM", "CID"],
    "ENSEMBL": ["ENSEMBL", "ENSEMBL_GENE"],
    "ENSEMBL_GENE": ["ENSEMBL_GENE", "ENSEMBL"],
    "NCBI": ["NCBI", "NCBI.GENE", "ENTREZ"],
    "NCBI.GENE": ["NCBI.GENE", "NCBI", "ENTREZ"],
    "ENTREZ": ["ENTREZ", "NCBI", "NCBI.GENE"],
    "MONDO": ["MONDO", "OPENTARGETS_DISEASE"],
    "DOID": ["DOID", "OPENTARGETS_DISEASE"],
    "MESH": ["MESH", "OPENTARGETS_DISEASE"],
    "EFO": ["EFO", "OPENTARGETS_DISEASE"],
    "HP": ["HP", "OPENTARGETS_DISEASE"],
    "ORPHANET": ["ORPHANET", "OPENTARGETS_DISEASE"],
    "OTAR": ["OTAR", "OPENTARGETS_DISEASE"],
    "NCIT": ["NCIT", "OPENTARGETS_DISEASE"],
    "ONCOTREE": ["ONCOTREE", "OPENTARGETS_DISEASE"],
    "CL": ["CL"],
    "UBERON": ["UBERON"],
    "METABO_STATE": ["METABO_STATE"],
    "UNIPROT": ["UNIPROT"],
}


def require_arrow() -> None:
    if pa is None or pc is None or ds is None:
        raise RuntimeError("pyarrow is required for resolver queries.")


def resolver_index_path(graph_root: Path, release_id: str | None) -> Path:
    resolved = release_id or latest_normalized_release_id(graph_root)
    path = graph_root / resolved / "resolver_index.parquet"
    if not path.exists():
        raise FileNotFoundError(f"Resolver index does not exist: {path}")
    return path


def candidate_queries(query: str) -> list[tuple[str, str, int]]:
    raw = str(query or "").strip()
    text_key = normalize_lookup_key(raw)
    id_key = normalize_id_token(raw)
    candidates: list[tuple[str, str, int]] = []

    def add(namespace: str, key: str, priority: int) -> None:
        if key:
            item = (namespace.upper(), key, priority)
            if item not in candidates:
                candidates.append(item)

    prefix_match = re.match(r"^([A-Za-z][A-Za-z0-9_.-]*):(.+)$", id_key)
    if raw.casefold().startswith("metabo_state:"):
        add("METABO_STATE", normalize_id_token(raw.split(":", 1)[1]), 120)
    if prefix_match:
        prefix = prefix_match.group(1).upper()
        value = prefix_match.group(2)
        for namespace in ID_PREFIX_ALIASES.get(prefix, [prefix]):
            add(namespace, value, 120)
        add("ID", id_key, 90)
    if re.match(r"^hmdb[0-9]{7}$", id_key):
        add("HMDB", id_key, 120)
    if re.match(r"^chebi:[0-9]+$", id_key):
        add("CHEBI", id_key.split(":", 1)[1], 120)
    if re.match(r"^cl:[0-9]+$", id_key):
        add("CL", id_key.split(":", 1)[1], 120)
    if re.match(r"^uberon:[0-9]+$", id_key):
        add("UBERON", id_key.split(":", 1)[1], 120)
    if re.match(r"^c[0-9]{5}$", id_key):
        add("KEGG", id_key, 110)
    if re.match(r"^ensg[0-9]+", id_key):
        add("ENSEMBL", id_key, 120)
        add("ENSEMBL_GENE", id_key, 115)
        add("OPENTARGETS", id_key, 90)
    if re.match(r"^[a-z][0-9][a-z0-9]{3}[0-9]$", id_key) or re.match(r"^[a-z0-9]{10}$", id_key):
        add("UNIPROT", id_key, 100)
    if re.match(r"^[a-z0-9]{14}-[a-z0-9]{10}-[a-z0-9]$", id_key):
        add("INCHIKEY", id_key, 120)
    if re.match(r"^[0-9]+$", id_key):
        add("CID", id_key, 85)
        add("PUBCHEM", id_key, 80)
        add("NCBI", id_key, 70)
        add("ENTREZ", id_key, 70)
    add("TEXT", text_key, 10)
    return candidates


def resolve_entities(
    index_path: Path,
    query: str,
    entity_type: str | None = None,
    namespace: str | None = None,
    limit: int = 10,
) -> list[dict[str, Any]]:
    require_arrow()
    if namespace:
        key = normalize_lookup_key(query) if namespace.upper() == "TEXT" else normalize_id_token(query)
        candidates = [(namespace.upper(), key, 100)]
    else:
        candidates = candidate_queries(query)

    dataset = ds.dataset(index_path, format="parquet")
    rows: list[dict[str, Any]] = []
    for ns, key, priority in candidates:
        filt = (ds.field("namespace") == ns) & (ds.field("lookup_key") == key)
        if entity_type:
            filt = filt & (ds.field("entity_type") == entity_type.lower())
        table = dataset.to_table(filter=filt)
        for row in table.to_pylist():
            row["_query_priority"] = priority
            rows.append(row)

    deduped: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    for row in rows:
        dedupe_key = (row["entity_uid"], row["namespace"], row["lookup_key"], row["match_field"])
        prev = deduped.get(dedupe_key)
        if prev is None or (row["_query_priority"], row["rank"]) > (prev["_query_priority"], prev["rank"]):
            deduped[dedupe_key] = row

    ordered = sorted(
        deduped.values(),
        key=lambda row: (-row["_query_priority"], -row["rank"], row["entity_type"], row["entity_uid"], row["raw_value"]),
    )
    for row in ordered:
        row.pop("_query_priority", None)
    return ordered[:limit]


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Resolve entity names or IDs against a graph_projection resolver index.")
    parser.add_argument("query")
    parser.add_argument("--graph-root", default="graph_projection")
    parser.add_argument("--release-id", default="")
    parser.add_argument(
        "--entity-type",
        choices=["metabolite", "gene", "disease", "pathway", "reaction", "target", "protein", "cell_type", "cell_state", "tissue"],
    )
    parser.add_argument("--namespace", help="Force a resolver namespace, e.g. HMDB, CHEBI, TEXT, ENSEMBL, MONDO.")
    parser.add_argument("--limit", type=int, default=10)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    graph_root = (Path.cwd() / args.graph_root).resolve()
    path = resolver_index_path(graph_root, args.release_id or None)
    rows = resolve_entities(path, args.query, args.entity_type, args.namespace, args.limit)
    print(json.dumps({"query": args.query, "resolver_index": str(path), "matches": rows}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
