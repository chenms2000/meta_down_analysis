"""Build a PubChem CID property cache for CIDs already referenced by normalized_store.

This is a derived open-core cache. It never changes canonical metabolite identity;
it only materializes PubChem Compound Extras fields for the CID xrefs already
present in a frozen normalized release.
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import re
import sys
import time
import unicodedata
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

try:
    import pyarrow as pa
    import pyarrow.dataset as ds
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover - handled by require_arrow.
    pa = None
    ds = None
    pq = None

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from build_graph_projection import latest_normalized_release_id, normalize_id_token, normalize_lookup_key  # noqa: E402
from build_normalized_store import json_dumps, stable_uid  # noqa: E402


DEFAULT_NORMALIZED_ROOT = "normalized_store"
DEFAULT_RAW_ROOT = "raw_lake"
DEFAULT_OUTPUT_ROOT = "pubchem_cid_cache"
PUBCHEM_SOURCE_ID = "pubchem_compound_extras"
PROPERTY_FILES = {
    "title": "CID-Title.gz",
    "inchi_key": "CID-InChI-Key.gz",
    "smiles": "CID-SMILES.gz",
    "mass": "CID-Mass.gz",
    "synonyms": "CID-Synonym-filtered.gz",
}


def require_arrow() -> None:
    if pa is None or ds is None or pq is None:
        raise RuntimeError("pyarrow is required to build PubChem CID cache parquet outputs.")


def utc_now() -> str:
    return dt.datetime.now(dt.UTC).isoformat()


def maybe_float(value: str) -> float | None:
    if value in {None, ""}:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parser_hash() -> str:
    text = Path(__file__).read_text(encoding="utf-8")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def schema_for(name: str) -> pa.Schema:
    list_str = pa.list_(pa.string())
    schemas = {
        "cid_metabolite_links": pa.schema(
            [
                ("pubchem_cid", pa.string()),
                ("metabolite_uid", pa.string()),
                ("xref_uid", pa.string()),
                ("xref_key", pa.string()),
                ("normalized_source_name", pa.string()),
                ("normalized_source_release", pa.string()),
                ("normalized_license_id", pa.string()),
                ("normalized_parser_hash", pa.string()),
            ]
        ),
        "cid_properties": pa.schema(
            [
                ("pubchem_cid", pa.string()),
                ("metabolite_uids", list_str),
                ("title", pa.string()),
                ("formula", pa.string()),
                ("monoisotopic_mass", pa.float64()),
                ("exact_mass", pa.float64()),
                ("inchi", pa.string()),
                ("inchikey", pa.string()),
                ("isomeric_smiles", pa.string()),
                ("synonyms", list_str),
                ("synonym_count", pa.int64()),
                ("fields_present", list_str),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "cid_missing": pa.schema(
            [
                ("pubchem_cid", pa.string()),
                ("metabolite_uids", list_str),
                ("reason", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "cid_lookup_index": pa.schema(
            [
                ("lookup_uid", pa.string()),
                ("pubchem_cid", pa.string()),
                ("metabolite_uids", list_str),
                ("namespace", pa.string()),
                ("lookup_key", pa.string()),
                ("raw_value", pa.string()),
                ("match_field", pa.string()),
                ("rank", pa.float64()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
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
        self.writer = pq.ParquetWriter(path, schema_for(schema_name), compression="snappy")
        self.rows = 0

    def write(self, rows: Iterable[dict[str, Any]]) -> None:
        batch = list(rows)
        if not batch:
            return
        self.writer.write_table(pa.Table.from_pylist(batch, schema=schema_for(self.schema_name)))
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass
class PubChemProperty:
    pubchem_cid: str
    metabolite_uids: set[str] = field(default_factory=set)
    title: str = ""
    formula: str = ""
    monoisotopic_mass: float | None = None
    exact_mass: float | None = None
    inchi: str = ""
    inchikey: str = ""
    isomeric_smiles: str = ""
    synonyms: list[str] = field(default_factory=list)
    synonym_count: int = 0

    def fields_present(self) -> list[str]:
        fields = []
        if self.title:
            fields.append("title")
        if self.formula:
            fields.append("formula")
        if self.monoisotopic_mass is not None:
            fields.append("monoisotopic_mass")
        if self.exact_mass is not None:
            fields.append("exact_mass")
        if self.inchi:
            fields.append("inchi")
        if self.inchikey:
            fields.append("inchikey")
        if self.isomeric_smiles:
            fields.append("isomeric_smiles")
        if self.synonyms:
            fields.append("synonyms")
        return fields

    def has_any_property(self) -> bool:
        return bool(self.fields_present())


@dataclass
class CacheContext:
    normalized_dir: Path
    raw_dir: Path
    output_dir: Path
    release_id: str
    max_synonyms_per_cid: int
    include_synonyms: bool
    created_at_utc: str = field(default_factory=utc_now)
    parser_hash: str = field(default_factory=parser_hash)
    table_audit: list[dict[str, Any]] = field(default_factory=list)
    file_audit: list[dict[str, Any]] = field(default_factory=list)

    def writer(self, filename: str, schema_name: str) -> ParquetWriter:
        return ParquetWriter(self.output_dir / filename, schema_name)


def table_path(root: Path, name: str) -> Path:
    path = root / f"{name}.parquet"
    if not path.exists():
        raise FileNotFoundError(f"Missing normalized table: {path}")
    return path


def raw_file_metadata(normalized_dir: Path) -> dict[str, dict[str, Any]]:
    dataset = ds.dataset(table_path(normalized_dir, "raw_file_manifest"), format="parquet")
    table = dataset.to_table(filter=ds.field("source_id") == PUBCHEM_SOURCE_ID)
    metadata: dict[str, dict[str, Any]] = {}
    for row in table.to_pylist():
        path = Path(row.get("path", ""))
        metadata[path.name] = row
    return metadata


def read_pubchem_links(ctx: CacheContext) -> tuple[set[str], dict[str, set[str]]]:
    dataset = ds.dataset(table_path(ctx.normalized_dir, "metabolite_xrefs"), format="parquet")
    table = dataset.to_table(
        columns=[
            "xref_uid",
            "metabolite_uid",
            "xref_source",
            "xref_id",
            "xref_key",
            "source_name",
            "source_release",
            "license_id",
            "parser_hash",
        ],
        filter=ds.field("xref_source") == "PUBCHEM.COMPOUND",
    )
    writer = ctx.writer("cid_metabolite_links.parquet", "cid_metabolite_links")
    links: list[dict[str, Any]] = []
    target_cids: set[str] = set()
    cid_to_metabolites: dict[str, set[str]] = defaultdict(set)
    for row in table.to_pylist():
        cid = str(row.get("xref_id", "") or "").strip()
        metabolite_uid = str(row.get("metabolite_uid", "") or "")
        if cid and cid.isdigit() and int(cid) > 0:
            target_cids.add(cid)
            cid_to_metabolites[cid].add(metabolite_uid)
        links.append(
            {
                "pubchem_cid": cid,
                "metabolite_uid": metabolite_uid,
                "xref_uid": str(row.get("xref_uid", "") or ""),
                "xref_key": str(row.get("xref_key", "") or ""),
                "normalized_source_name": str(row.get("source_name", "") or ""),
                "normalized_source_release": str(row.get("source_release", "") or ""),
                "normalized_license_id": str(row.get("license_id", "") or ""),
                "normalized_parser_hash": str(row.get("parser_hash", "") or ""),
            }
        )
        if len(links) >= 100_000:
            writer.write(links)
            links = []
    writer.write(links)
    audit = writer.close()
    ctx.table_audit.append(audit)
    return target_cids, cid_to_metabolites


def ensure_properties(target_cids: set[str], cid_to_metabolites: dict[str, set[str]]) -> dict[str, PubChemProperty]:
    return {
        cid: PubChemProperty(pubchem_cid=cid, metabolite_uids=set(cid_to_metabolites.get(cid, set())))
        for cid in sorted(target_cids, key=lambda x: int(x))
    }


def iter_gzip_tsv(path: Path) -> Iterator[list[str]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            yield line.rstrip("\n").split("\t")


def audit_file(ctx: CacheContext, filename: str, rows_scanned: int, matched_rows: int, seconds: float) -> None:
    path = ctx.raw_dir / filename
    manifest_meta = raw_file_metadata(ctx.normalized_dir).get(filename, {})
    row = {
        "file": filename,
        "path": str(path),
        "bytes": manifest_meta.get("bytes", path.stat().st_size if path.exists() else 0),
        "sha256": manifest_meta.get("sha256", ""),
        "rows_scanned": rows_scanned,
        "matched_rows": matched_rows,
        "seconds": round(seconds, 3),
    }
    ctx.file_audit.append(row)
    print(
        f"[pubchem_cid_cache] scanned {filename}: rows={rows_scanned} matched={matched_rows} seconds={row['seconds']}",
        flush=True,
    )


def scan_title(ctx: CacheContext, properties: dict[str, PubChemProperty]) -> None:
    started = time.time()
    path = ctx.raw_dir / PROPERTY_FILES["title"]
    scanned = matched = 0
    for parts in iter_gzip_tsv(path):
        scanned += 1
        if len(parts) < 2:
            continue
        prop = properties.get(parts[0])
        if prop is None:
            continue
        prop.title = parts[1]
        matched += 1
    audit_file(ctx, path.name, scanned, matched, time.time() - started)


def scan_inchi_key(ctx: CacheContext, properties: dict[str, PubChemProperty]) -> None:
    started = time.time()
    path = ctx.raw_dir / PROPERTY_FILES["inchi_key"]
    scanned = matched = 0
    for parts in iter_gzip_tsv(path):
        scanned += 1
        if len(parts) < 3:
            continue
        prop = properties.get(parts[0])
        if prop is None:
            continue
        prop.inchi = parts[1]
        prop.inchikey = parts[2]
        matched += 1
    audit_file(ctx, path.name, scanned, matched, time.time() - started)


def scan_smiles(ctx: CacheContext, properties: dict[str, PubChemProperty]) -> None:
    started = time.time()
    path = ctx.raw_dir / PROPERTY_FILES["smiles"]
    scanned = matched = 0
    for parts in iter_gzip_tsv(path):
        scanned += 1
        if len(parts) < 2:
            continue
        prop = properties.get(parts[0])
        if prop is None:
            continue
        prop.isomeric_smiles = parts[1]
        matched += 1
    audit_file(ctx, path.name, scanned, matched, time.time() - started)


def scan_mass(ctx: CacheContext, properties: dict[str, PubChemProperty]) -> None:
    started = time.time()
    path = ctx.raw_dir / PROPERTY_FILES["mass"]
    scanned = matched = 0
    for parts in iter_gzip_tsv(path):
        scanned += 1
        if len(parts) < 4:
            continue
        prop = properties.get(parts[0])
        if prop is None:
            continue
        prop.formula = parts[1]
        prop.monoisotopic_mass = maybe_float(parts[2])
        prop.exact_mass = maybe_float(parts[3])
        matched += 1
    audit_file(ctx, path.name, scanned, matched, time.time() - started)


def scan_synonyms(ctx: CacheContext, properties: dict[str, PubChemProperty]) -> None:
    started = time.time()
    path = ctx.raw_dir / PROPERTY_FILES["synonyms"]
    scanned = matched = 0
    for parts in iter_gzip_tsv(path):
        scanned += 1
        if len(parts) < 2:
            continue
        prop = properties.get(parts[0])
        if prop is None:
            continue
        matched += 1
        prop.synonym_count += 1
        synonym = parts[1].strip()
        if synonym and len(prop.synonyms) < ctx.max_synonyms_per_cid and synonym not in prop.synonyms:
            prop.synonyms.append(synonym)
    audit_file(ctx, path.name, scanned, matched, time.time() - started)


def write_property_tables(ctx: CacheContext, properties: dict[str, PubChemProperty]) -> None:
    source_release = ctx.release_id
    license_id = f"open_core:{PUBCHEM_SOURCE_ID}"
    prop_writer = ctx.writer("cid_properties.parquet", "cid_properties")
    missing_writer = ctx.writer("cid_missing.parquet", "cid_missing")
    prop_rows: list[dict[str, Any]] = []
    missing_rows: list[dict[str, Any]] = []
    for cid in sorted(properties, key=lambda x: int(x)):
        prop = properties[cid]
        base = {
            "pubchem_cid": cid,
            "metabolite_uids": sorted(prop.metabolite_uids),
            "source_release": source_release,
            "license_id": license_id,
            "parser_hash": ctx.parser_hash,
        }
        if prop.has_any_property():
            prop_rows.append(
                {
                    **base,
                    "title": prop.title,
                    "formula": prop.formula,
                    "monoisotopic_mass": prop.monoisotopic_mass,
                    "exact_mass": prop.exact_mass,
                    "inchi": prop.inchi,
                    "inchikey": prop.inchikey,
                    "isomeric_smiles": prop.isomeric_smiles,
                    "synonyms": prop.synonyms,
                    "synonym_count": prop.synonym_count,
                    "fields_present": prop.fields_present(),
                }
            )
        else:
            missing_rows.append({**base, "reason": "no_pubchem_extra_property_rows_matched"})
        if len(prop_rows) >= 100_000:
            prop_writer.write(prop_rows)
            prop_rows = []
        if len(missing_rows) >= 100_000:
            missing_writer.write(missing_rows)
            missing_rows = []
    prop_writer.write(prop_rows)
    missing_writer.write(missing_rows)
    ctx.table_audit.append(prop_writer.close())
    ctx.table_audit.append(missing_writer.close())


def lookup_rows_for_property(ctx: CacheContext, prop: PubChemProperty) -> Iterator[dict[str, Any]]:
    source_release = ctx.release_id
    license_id = f"open_core:{PUBCHEM_SOURCE_ID}"
    metabolite_uids = sorted(prop.metabolite_uids)

    def row(namespace: str, raw_value: str, match_field: str, rank: float, key: str | None = None) -> dict[str, Any] | None:
        raw = str(raw_value or "").strip()
        if not raw:
            return None
        lookup_key = key if key is not None else normalize_lookup_key(raw)
        if not lookup_key:
            return None
        return {
            "lookup_uid": stable_uid("pubchem_lookup", prop.pubchem_cid, namespace, lookup_key, match_field),
            "pubchem_cid": prop.pubchem_cid,
            "metabolite_uids": metabolite_uids,
            "namespace": namespace,
            "lookup_key": lookup_key,
            "raw_value": raw,
            "match_field": match_field,
            "rank": float(rank),
            "source_release": source_release,
            "license_id": license_id,
            "parser_hash": ctx.parser_hash,
        }

    seed_rows = [
        row("CID", prop.pubchem_cid, "pubchem_cid", 100.0, key=normalize_id_token(prop.pubchem_cid)),
        row("PUBCHEM", prop.pubchem_cid, "pubchem_cid", 100.0, key=normalize_id_token(prop.pubchem_cid)),
        row("PUBCHEM.COMPOUND", prop.pubchem_cid, "pubchem_cid", 100.0, key=normalize_id_token(prop.pubchem_cid)),
        row("ID", f"PUBCHEM.COMPOUND:{prop.pubchem_cid}", "pubchem_cid", 95.0, key=f"pubchem.compound:{prop.pubchem_cid}"),
        row("INCHIKEY", prop.inchikey, "inchikey", 95.0, key=normalize_id_token(prop.inchikey)),
        row("TEXT", prop.title, "title", 70.0),
    ]
    for candidate in seed_rows:
        if candidate:
            yield candidate
    seen_synonyms: set[str] = set()
    for synonym in prop.synonyms:
        key = normalize_lookup_key(synonym)
        if not key or key in seen_synonyms:
            continue
        seen_synonyms.add(key)
        candidate = row("TEXT", synonym, "synonym", 50.0, key=key)
        if candidate:
            yield candidate


def write_lookup_index(ctx: CacheContext, properties: dict[str, PubChemProperty]) -> None:
    writer = ctx.writer("cid_lookup_index.parquet", "cid_lookup_index")
    rows: list[dict[str, Any]] = []
    for cid in sorted(properties, key=lambda x: int(x)):
        prop = properties[cid]
        if not prop.has_any_property():
            continue
        rows.extend(lookup_rows_for_property(ctx, prop))
        if len(rows) >= 100_000:
            writer.write(rows)
            rows = []
    writer.write(rows)
    ctx.table_audit.append(writer.close())


def write_manifest(ctx: CacheContext, properties: dict[str, PubChemProperty], target_cids: set[str]) -> Path:
    matched = sum(1 for prop in properties.values() if prop.has_any_property())
    manifest = {
        "release_id": ctx.release_id,
        "created_at_utc": ctx.created_at_utc,
        "normalized_dir": str(ctx.normalized_dir),
        "raw_dir": str(ctx.raw_dir),
        "output_dir": str(ctx.output_dir),
        "target_cid_count": len(target_cids),
        "matched_cid_count": matched,
        "missing_cid_count": len(target_cids) - matched,
        "include_synonyms": ctx.include_synonyms,
        "max_synonyms_per_cid": ctx.max_synonyms_per_cid,
        "source_files": ctx.file_audit,
        "tables": ctx.table_audit,
        "notes": [
            "Cache is restricted to PubChem CIDs already present in normalized_store metabolite_xrefs.",
            "PubChem properties are an open-core lookup enhancement and do not define canonical metabolite identity.",
        ],
    }
    path = ctx.output_dir / "pubchem_cid_cache_manifest.json"
    path.write_text(json_dumps(manifest) + "\n", encoding="utf-8")
    return path


def build_pubchem_cid_cache(
    normalized_root: Path,
    raw_root: Path,
    output_root: Path,
    release_id: str | None,
    max_synonyms_per_cid: int = 20,
    include_synonyms: bool = True,
) -> Path:
    require_arrow()
    resolved_release = release_id or latest_normalized_release_id(normalized_root)
    normalized_dir = normalized_root / resolved_release
    raw_dir = raw_root / PUBCHEM_SOURCE_ID / resolved_release
    if not normalized_dir.exists():
        raise FileNotFoundError(f"Normalized release does not exist: {normalized_dir}")
    if not raw_dir.exists():
        raise FileNotFoundError(f"PubChem raw directory does not exist: {raw_dir}")
    output_dir = output_root / resolved_release
    output_dir.mkdir(parents=True, exist_ok=True)
    ctx = CacheContext(
        normalized_dir=normalized_dir,
        raw_dir=raw_dir,
        output_dir=output_dir,
        release_id=resolved_release,
        max_synonyms_per_cid=max_synonyms_per_cid,
        include_synonyms=include_synonyms,
    )
    target_cids, cid_to_metabolites = read_pubchem_links(ctx)
    properties = ensure_properties(target_cids, cid_to_metabolites)
    scan_title(ctx, properties)
    scan_inchi_key(ctx, properties)
    scan_smiles(ctx, properties)
    scan_mass(ctx, properties)
    if include_synonyms:
        scan_synonyms(ctx, properties)
    write_property_tables(ctx, properties)
    write_lookup_index(ctx, properties)
    return write_manifest(ctx, properties, target_cids)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build PubChem CID property cache for a normalized release.")
    parser.add_argument("--normalized-root", default=DEFAULT_NORMALIZED_ROOT)
    parser.add_argument("--raw-root", default=DEFAULT_RAW_ROOT)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--release-id", default="")
    parser.add_argument("--max-synonyms-per-cid", type=int, default=20)
    parser.add_argument("--skip-synonyms", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    workspace = Path.cwd()
    manifest = build_pubchem_cid_cache(
        normalized_root=(workspace / args.normalized_root).resolve(),
        raw_root=(workspace / args.raw_root).resolve(),
        output_root=(workspace / args.output_root).resolve(),
        release_id=args.release_id or None,
        max_synonyms_per_cid=max(0, args.max_synonyms_per_cid),
        include_synonyms=not args.skip_synonyms,
    )
    print(f"[pubchem_cid_cache] done: {manifest}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
