"""Promote reviewed compound identity overlays into normalized graph tables.

This script turns confirmed PubChem/InChIKey overlay compounds into normalized
metabolite nodes, xrefs, identity bridge edges, and copied pathway/reaction
links from exact local same-identity metabolites. It is intentionally
conservative: pathway/reaction links are created only when an overlay compound
has an exact InChIKey or stable xref match to an existing local metabolite.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any

try:
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover - handled by require_arrow.
    pa = None
    pc = None
    pq = None

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from build_normalized_store import json_dumps, stable_uid  # noqa: E402


DEFAULT_NORMALIZED_ROOT = "normalized_store"
DEFAULT_OVERLAY_ROOT = "manual_sources/compound_identity_overlays"
OVERLAY_SOURCE_NAME = "compound_identity_overlay"
OVERLAY_LICENSE_ID = "local:european_trait_identity_review"


def require_arrow() -> None:
    if pa is None or pc is None or pq is None:
        raise RuntimeError("pyarrow is required to promote compound identity overlays.")


def utc_now() -> str:
    return dt.datetime.now(dt.UTC).isoformat()


def parser_hash() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_overlay(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return [
            {str(key or "").strip().casefold().replace(" ", "_"): str(value or "").strip() for key, value in row.items()}
            for row in reader
        ]


def split_values(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    out: list[str] = []
    seen: set[str] = set()
    for part in re.split(r"[|;]", text):
        cleaned = str(part or "").strip()
        key = cleaned.casefold()
        if cleaned and key not in seen:
            seen.add(key)
            out.append(cleaned)
    return out


def dedupe(values: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        cleaned = str(value or "").strip()
        if cleaned and cleaned not in seen:
            seen.add(cleaned)
            out.append(cleaned)
    return out


def maybe_float(value: Any) -> float | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def infer_charge(formula: str) -> int:
    text = str(formula or "").strip()
    if text.endswith("-"):
        return -1
    if text.endswith("+"):
        return 1
    match = re.search(r"([+-])(\d+)$", text)
    if not match:
        return 0
    value = int(match.group(2))
    return value if match.group(1) == "+" else -value


def normalize_xref_source(source: str) -> str:
    key = str(source or "").strip().upper()
    aliases = {
        "CID": "PUBCHEM.COMPOUND",
        "PUBCHEM": "PUBCHEM.COMPOUND",
        "PUBCHEM_CID": "PUBCHEM.COMPOUND",
        "KEGG": "KEGG.COMPOUND",
    }
    return aliases.get(key, key)


def overlay_xrefs(row: dict[str, str]) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    cid = str(row.get("pubchem_cid") or "").strip()
    if cid:
        pairs.append(("PUBCHEM.COMPOUND", cid))
    inchikey = str(row.get("inchikey") or "").strip()
    if inchikey:
        pairs.append(("INCHIKEY", inchikey))
    for source, field_name in (
        ("HMDB", "hmdb_id"),
        ("CHEBI", "chebi_id"),
        ("KEGG.COMPOUND", "kegg_id"),
    ):
        for value in split_values(row.get(field_name)):
            pairs.append((source, value.split(":", 1)[1] if ":" in value and value.upper().startswith(source.split(".")[0]) else value))
    seen: set[tuple[str, str]] = set()
    out: list[tuple[str, str]] = []
    for source, xid in pairs:
        normalized = normalize_xref_source(source)
        cleaned = str(xid or "").strip()
        key = (normalized, cleaned.upper() if normalized in {"INCHIKEY", "HMDB"} else cleaned)
        if cleaned and key not in seen:
            seen.add(key)
            out.append((normalized, cleaned))
    return out


def xref_key(source: str, xid: str) -> str:
    return f"{normalize_xref_source(source)}:{str(xid or '').strip()}"


def table_path(normalized_dir: Path, table: str) -> Path:
    return normalized_dir / f"{table}.parquet"


def read_table_if_exists(path: Path, schema: pa.Schema | None = None) -> pa.Table:
    if path.exists():
        return pq.read_table(path)
    if schema is None:
        raise FileNotFoundError(path)
    return pa.Table.from_pylist([], schema=schema)


def filter_out_values(table: pa.Table, column: str, values: set[str]) -> pa.Table:
    if not values or column not in table.column_names or table.num_rows == 0:
        return table
    field_type = table.schema.field(column).type
    mask = pc.invert(pc.is_in(table[column], value_set=pa.array(sorted(values), type=field_type)))
    return table.filter(mask)


def write_replaced_table(path: Path, base_table: pa.Table, new_rows: list[dict[str, Any]], schema: pa.Schema | None = None) -> dict[str, Any]:
    schema = schema or base_table.schema
    new_table = pa.Table.from_pylist(new_rows, schema=schema) if new_rows else pa.Table.from_pylist([], schema=schema)
    output = pa.concat_tables([base_table.cast(schema), new_table], promote_options="default")
    pq.write_table(output, path, compression="snappy")
    return {
        "table": path.stem,
        "path": str(path),
        "rows": output.num_rows,
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def identity_edge_schema() -> pa.Schema:
    return pa.schema(
        [
            ("edge_uid", pa.string()),
            ("subject_metabolite_uid", pa.string()),
            ("object_metabolite_uid", pa.string()),
            ("subject_external_id", pa.string()),
            ("object_external_id", pa.string()),
            ("predicate", pa.string()),
            ("source_name", pa.string()),
            ("source_record_id", pa.string()),
            ("evidence_level", pa.string()),
            ("match_basis", pa.string()),
            ("linked_name", pa.string()),
            ("source_release", pa.string()),
            ("license_id", pa.string()),
            ("parser_hash", pa.string()),
        ]
    )


def build_existing_indexes(
    metabolites: list[dict[str, Any]],
    xrefs: list[dict[str, Any]],
    overlay_uids: set[str],
) -> tuple[dict[str, list[dict[str, Any]]], dict[tuple[str, str], list[dict[str, Any]]]]:
    by_inchikey: dict[str, list[dict[str, Any]]] = {}
    by_xref: dict[tuple[str, str], list[dict[str, Any]]] = {}

    def add_xref(source: str, xid: str, metabolite: dict[str, Any]) -> None:
        source = normalize_xref_source(source)
        xid = str(xid or "").strip()
        if source and xid:
            by_xref.setdefault((source, xid.upper() if source in {"INCHIKEY", "HMDB"} else xid), []).append(metabolite)

    metabolites_by_uid = {row.get("metabolite_uid"): row for row in metabolites if row.get("metabolite_uid") not in overlay_uids}
    for row in metabolites_by_uid.values():
        inchikey = str(row.get("inchikey") or "").strip().upper()
        if inchikey:
            by_inchikey.setdefault(inchikey, []).append(row)
        for raw_xref in row.get("external_xrefs") or []:
            if ":" not in str(raw_xref):
                continue
            source, xid = str(raw_xref).split(":", 1)
            add_xref(source, xid, row)
    for row in xrefs:
        metabolite_uid = row.get("metabolite_uid")
        metabolite = metabolites_by_uid.get(metabolite_uid)
        if not metabolite:
            continue
        add_xref(row.get("xref_source", ""), row.get("xref_id", ""), metabolite)
    return by_inchikey, by_xref


def promote_overlay(
    normalized_dir: Path,
    overlay_path: Path,
    release_id: str,
) -> dict[str, Any]:
    require_arrow()
    overlay_rows = read_overlay(overlay_path)
    phash = parser_hash()
    started = utc_now()
    overlay_uids = {
        str(row.get("metabolite_uid") or f"manual_pubchem_cid_{row.get('pubchem_cid')}").strip()
        for row in overlay_rows
        if str(row.get("pubchem_cid") or row.get("metabolite_uid") or "").strip()
    }
    stats: dict[str, Any] = {
        "overlay_path": str(overlay_path),
        "input_rows": len(overlay_rows),
        "metabolite_rows_added": 0,
        "xref_rows_added": 0,
        "identity_edges_added": 0,
        "pathway_edges_added": 0,
        "reaction_participant_rows_added": 0,
        "linked_overlay_compounds": 0,
        "unlinked_overlay_compounds": 0,
        "created_at_utc": started,
        "tables": [],
    }
    if not overlay_rows:
        return stats

    metabolites_path = table_path(normalized_dir, "metabolites")
    xrefs_path = table_path(normalized_dir, "metabolite_xrefs")
    pathway_path = table_path(normalized_dir, "metabolite_pathway_edges")
    participants_path = table_path(normalized_dir, "reaction_participants")
    identity_path = table_path(normalized_dir, "metabolite_identity_edges")

    metabolites_table = read_table_if_exists(metabolites_path)
    xrefs_table = read_table_if_exists(xrefs_path)
    pathway_table = read_table_if_exists(pathway_path)
    participants_table = read_table_if_exists(participants_path)
    identity_table = read_table_if_exists(identity_path, identity_edge_schema())

    base_metabolites = filter_out_values(metabolites_table, "metabolite_uid", overlay_uids)
    base_xrefs = filter_out_values(xrefs_table, "metabolite_uid", overlay_uids)
    base_pathways = filter_out_values(pathway_table, "metabolite_uid", overlay_uids)
    base_participants = filter_out_values(participants_table, "participant_uid", overlay_uids)
    base_identity = filter_out_values(identity_table, "subject_metabolite_uid", overlay_uids)

    metabolites = base_metabolites.to_pylist()
    xrefs = base_xrefs.to_pylist()
    by_inchikey, by_xref = build_existing_indexes(metabolites, xrefs, overlay_uids)
    pathway_by_metabolite: dict[str, list[dict[str, Any]]] = {}
    for row in base_pathways.to_pylist():
        pathway_by_metabolite.setdefault(str(row.get("metabolite_uid") or ""), []).append(row)
    participants_by_metabolite: dict[str, list[dict[str, Any]]] = {}
    for row in base_participants.to_pylist():
        if str(row.get("participant_type") or "").lower() == "metabolite":
            participants_by_metabolite.setdefault(str(row.get("participant_uid") or ""), []).append(row)

    new_metabolites: list[dict[str, Any]] = []
    new_xrefs: list[dict[str, Any]] = []
    new_identity_edges: list[dict[str, Any]] = []
    new_pathway_edges: list[dict[str, Any]] = []
    new_participants: list[dict[str, Any]] = []

    for row in overlay_rows:
        cid = str(row.get("pubchem_cid") or "").strip()
        metabolite_uid = str(row.get("metabolite_uid") or f"manual_pubchem_cid_{cid}").strip()
        if not metabolite_uid:
            continue
        canonical_name = row.get("canonical_name") or row.get("name") or f"PubChem CID {cid}"
        synonyms = dedupe(split_values(row.get("synonyms")) + split_values(row.get("iupac_name")))
        xrefs_for_row = overlay_xrefs(row)
        external_xrefs = sorted({xref_key(source, xid) for source, xid in xrefs_for_row})
        primary_external_id = f"PUBCHEM.COMPOUND:{cid}" if cid else (external_xrefs[0] if external_xrefs else metabolite_uid)
        formula = str(row.get("formula") or "").strip()
        checksum = hashlib.sha256(json_dumps({"source": OVERLAY_SOURCE_NAME, **row}).encode("utf-8")).hexdigest()
        new_metabolites.append(
            {
                "metabolite_uid": metabolite_uid,
                "canonical_name": canonical_name,
                "synonyms": synonyms,
                "formula": formula,
                "exact_mass": maybe_float(row.get("exact_mass")),
                "charge": infer_charge(formula),
                "inchikey": str(row.get("inchikey") or "").strip(),
                "smiles": str(row.get("smiles") or "").strip(),
                "external_xrefs": external_xrefs,
                "source_priority": "reviewed_pubchem_overlay",
                "source_release": release_id,
                "license_id": row.get("license_id") or OVERLAY_LICENSE_ID,
                "checksum": checksum,
                "parser_hash": phash,
            }
        )
        for source, xid in xrefs_for_row:
            key = xref_key(source, xid)
            new_xrefs.append(
                {
                    "xref_uid": stable_uid("xref", metabolite_uid, source, xid, OVERLAY_SOURCE_NAME),
                    "metabolite_uid": metabolite_uid,
                    "xref_source": source,
                    "xref_id": xid,
                    "xref_key": key,
                    "source_name": OVERLAY_SOURCE_NAME,
                    "source_release": release_id,
                    "license_id": row.get("license_id") or OVERLAY_LICENSE_ID,
                    "parser_hash": phash,
                }
            )

        linked: dict[str, tuple[dict[str, Any], str]] = {}
        inchikey = str(row.get("inchikey") or "").strip().upper()
        for linked_row in by_inchikey.get(inchikey, []):
            linked[linked_row["metabolite_uid"]] = (linked_row, "inchikey_full")
        for source, xid in xrefs_for_row:
            key = (source, xid.upper() if source in {"INCHIKEY", "HMDB"} else xid)
            for linked_row in by_xref.get(key, []):
                linked.setdefault(linked_row["metabolite_uid"], (linked_row, f"xref:{source}"))

        if linked:
            stats["linked_overlay_compounds"] += 1
        else:
            stats["unlinked_overlay_compounds"] += 1
        for linked_uid, (linked_row, basis) in sorted(linked.items()):
            linked_primary = (linked_row.get("external_xrefs") or [linked_uid])[0]
            new_identity_edges.append(
                {
                    "edge_uid": stable_uid("edge", "metabolite_same_as", metabolite_uid, linked_uid, basis),
                    "subject_metabolite_uid": metabolite_uid,
                    "object_metabolite_uid": linked_uid,
                    "subject_external_id": primary_external_id,
                    "object_external_id": linked_primary,
                    "predicate": "same_as",
                    "source_name": OVERLAY_SOURCE_NAME,
                    "source_record_id": f"{primary_external_id}|same_as|{linked_primary}",
                    "evidence_level": "reviewed_identity_bridge",
                    "match_basis": basis,
                    "linked_name": str(linked_row.get("canonical_name") or ""),
                    "source_release": release_id,
                    "license_id": row.get("license_id") or OVERLAY_LICENSE_ID,
                    "parser_hash": phash,
                }
            )
            for pathway_edge in pathway_by_metabolite.get(linked_uid, []):
                pathway_uid = pathway_edge.get("pathway_uid")
                if not pathway_uid:
                    continue
                new_pathway_edges.append(
                    {
                        "edge_uid": stable_uid("edge", "overlay_metabolite_pathway", metabolite_uid, pathway_uid, linked_uid),
                        "metabolite_uid": metabolite_uid,
                        "pathway_uid": pathway_uid,
                        "metabolite_external_id": primary_external_id,
                        "pathway_external_id": pathway_edge.get("pathway_external_id", ""),
                        "predicate": "participates_in",
                        "source_name": OVERLAY_SOURCE_NAME,
                        "source_record_id": f"{primary_external_id}|same_as:{linked_uid}|{pathway_edge.get('pathway_external_id', '')}",
                        "evidence_level": "inferred_from_same_identity",
                        "evidence_code": f"identity_bridge:{basis}",
                        "species": pathway_edge.get("species", ""),
                        "source_release": release_id,
                        "license_id": f"{row.get('license_id') or OVERLAY_LICENSE_ID}+{pathway_edge.get('license_id', '')}",
                        "parser_hash": phash,
                    }
                )
            for participant in participants_by_metabolite.get(linked_uid, []):
                reaction_uid = participant.get("reaction_uid")
                if not reaction_uid:
                    continue
                new_participants.append(
                    {
                        "edge_uid": stable_uid("edge", "overlay_reaction_participant", metabolite_uid, reaction_uid, linked_uid, participant.get("role", "")),
                        "reaction_uid": reaction_uid,
                        "participant_uid": metabolite_uid,
                        "participant_external_id": primary_external_id,
                        "participant_type": "metabolite",
                        "role": participant.get("role", ""),
                        "source_name": OVERLAY_SOURCE_NAME,
                        "source_record_id": f"{primary_external_id}|same_as:{linked_uid}|{participant.get('source_record_id', '')}",
                        "source_release": release_id,
                        "license_id": f"{row.get('license_id') or OVERLAY_LICENSE_ID}+{participant.get('license_id', '')}",
                        "parser_hash": phash,
                    }
                )

    stats["metabolite_rows_added"] = len(new_metabolites)
    stats["xref_rows_added"] = len(new_xrefs)
    stats["identity_edges_added"] = len(new_identity_edges)
    stats["pathway_edges_added"] = len(new_pathway_edges)
    stats["reaction_participant_rows_added"] = len(new_participants)

    stats["tables"].append(write_replaced_table(metabolites_path, base_metabolites, new_metabolites))
    stats["tables"].append(write_replaced_table(xrefs_path, base_xrefs, new_xrefs))
    stats["tables"].append(write_replaced_table(identity_path, base_identity, new_identity_edges, identity_edge_schema()))
    stats["tables"].append(write_replaced_table(pathway_path, base_pathways, new_pathway_edges))
    stats["tables"].append(write_replaced_table(participants_path, base_participants, new_participants))
    update_manifest(normalized_dir, stats)
    return stats


def update_manifest(normalized_dir: Path, stats: dict[str, Any]) -> None:
    manifest_path = normalized_dir / "normalized_manifest.json"
    if not manifest_path.exists():
        return
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    table_updates = {row["table"]: row for row in stats.get("tables", [])}
    existing_tables = {row.get("table"): row for row in manifest.get("tables", [])}
    for table_name, update in table_updates.items():
        if table_name in existing_tables:
            existing_tables[table_name].update(update)
        else:
            existing_tables[table_name] = {**update, "seconds": 0.0}
    manifest["tables"] = list(existing_tables.values())
    notes = list(manifest.get("notes") or [])
    notes = [note for note in notes if note.get("component") != OVERLAY_SOURCE_NAME]
    notes.append(
        {
            "note_uid": stable_uid("note", OVERLAY_SOURCE_NAME, stats.get("created_at_utc", "")),
            "component": OVERLAY_SOURCE_NAME,
            "severity": "info",
            "created_at_utc": stats.get("created_at_utc", utc_now()),
            "message": (
                f"Promoted {stats.get('metabolite_rows_added', 0)} reviewed PubChem overlay compounds; "
                f"added {stats.get('identity_edges_added', 0)} same-identity bridges, "
                f"{stats.get('pathway_edges_added', 0)} inferred pathway links, and "
                f"{stats.get('reaction_participant_rows_added', 0)} inferred reaction participant links."
            ),
        }
    )
    manifest["notes"] = notes
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Promote compound identity overlays into normalized graph tables.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", default="mvp_20260513T002254")
    parser.add_argument("--normalized-root", default=DEFAULT_NORMALIZED_ROOT)
    parser.add_argument("--overlay-root", default=DEFAULT_OVERLAY_ROOT)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    workspace = Path(args.workspace).resolve()
    normalized_dir = (workspace / args.normalized_root / args.release_id).resolve()
    overlay_path = (workspace / args.overlay_root / args.release_id / "compound_identity_overlay.csv").resolve()
    stats = promote_overlay(normalized_dir, overlay_path, args.release_id)
    print(json.dumps(stats, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
