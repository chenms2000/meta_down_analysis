"""Build compound-aware matching indexes for uploaded metabolite tables.

This derived layer is optimized for deterministic compound resolution. It keeps
identity evidence separate by mode: external identifiers, normalized names,
InChIKey layers, and formula/exact-mass candidates.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
import sys
import time
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover - handled by require_arrow.
    pa = None
    pq = None

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from build_graph_projection import latest_normalized_release_id, normalize_id_token, normalize_lookup_key  # noqa: E402
from build_normalized_store import stable_uid  # noqa: E402


DEFAULT_NORMALIZED_ROOT = "normalized_store"
DEFAULT_PUBCHEM_ROOT = "pubchem_cid_cache"
DEFAULT_OUTPUT_ROOT = "compound_match_index"
DEFAULT_MANUAL_FEATURE_ROOT = "manual_sources/compound_features"
DEFAULT_MANUAL_COMPOUND_ROOT = "manual_sources/compound_identity_overlays"


def require_arrow() -> None:
    if pa is None or pq is None:
        raise RuntimeError("pyarrow is required to build compound match indexes.")


def utc_now() -> str:
    return dt.datetime.now(dt.UTC).isoformat()


def parser_hash() -> str:
    text = Path(__file__).read_text(encoding="utf-8")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_table_rows(path: Path, columns: list[str], batch_size: int = 50_000) -> Iterator[dict[str, Any]]:
    parquet_file = pq.ParquetFile(path)
    available = set(parquet_file.schema_arrow.names)
    selected = [column for column in columns if column in available]
    for batch in parquet_file.iter_batches(columns=selected, batch_size=batch_size):
        for row in batch.to_pylist():
            yield row


def maybe_float(value: Any) -> float | None:
    if value in {None, ""}:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if number == number else None


def normalize_formula(value: str) -> str:
    raw = re.sub(r"\s+", "", str(value or ""))
    if not raw:
        return ""
    tokens = re.findall(r"([A-Z][a-z]?)([0-9]*)", raw)
    if not tokens or "".join(f"{element}{count}" for element, count in tokens) != raw:
        return raw
    counts: dict[str, int] = {}
    for element, count_text in tokens:
        counts[element] = counts.get(element, 0) + int(count_text or "1")
    order = []
    if "C" in counts:
        order.append("C")
    if "H" in counts:
        order.append("H")
    order.extend(sorted(element for element in counts if element not in {"C", "H"}))
    return "".join(f"{element}{counts[element] if counts[element] != 1 else ''}" for element in order)


def inchikey_connectivity(value: str) -> str:
    key = normalize_id_token(value).upper()
    return key.split("-", 1)[0] if "-" in key else key[:14]


def schema_for(name: str) -> pa.Schema:
    schemas = {
        "compound_identifier_index": pa.schema(
            [
                ("lookup_uid", pa.string()),
                ("metabolite_uid", pa.string()),
                ("namespace", pa.string()),
                ("lookup_key", pa.string()),
                ("raw_value", pa.string()),
                ("match_field", pa.string()),
                ("rank", pa.float64()),
                ("source_table", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "compound_name_index": pa.schema(
            [
                ("lookup_uid", pa.string()),
                ("metabolite_uid", pa.string()),
                ("name_key", pa.string()),
                ("raw_value", pa.string()),
                ("match_field", pa.string()),
                ("rank", pa.float64()),
                ("source_table", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "compound_inchikey_index": pa.schema(
            [
                ("lookup_uid", pa.string()),
                ("metabolite_uid", pa.string()),
                ("inchikey", pa.string()),
                ("connectivity_key", pa.string()),
                ("raw_value", pa.string()),
                ("match_field", pa.string()),
                ("rank", pa.float64()),
                ("source_table", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "compound_formula_mass_index": pa.schema(
            [
                ("lookup_uid", pa.string()),
                ("metabolite_uid", pa.string()),
                ("formula", pa.string()),
                ("formula_key", pa.string()),
                ("exact_mass", pa.float64()),
                ("mass_mda_bucket", pa.int64()),
                ("match_field", pa.string()),
                ("rank", pa.float64()),
                ("source_table", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "compound_rt_index": pa.schema(
            [
                ("lookup_uid", pa.string()),
                ("metabolite_uid", pa.string()),
                ("rt", pa.float64()),
                ("rt_unit", pa.string()),
                ("method_key", pa.string()),
                ("source_name", pa.string()),
                ("source_record_id", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "compound_ms2_index": pa.schema(
            [
                ("lookup_uid", pa.string()),
                ("metabolite_uid", pa.string()),
                ("precursor_mz", pa.float64()),
                ("adduct", pa.string()),
                ("ion_mode", pa.string()),
                ("fragment_mz", pa.float64()),
                ("fragment_intensity", pa.float64()),
                ("source_name", pa.string()),
                ("source_record_id", pa.string()),
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
class CompoundIndexContext:
    normalized_dir: Path
    pubchem_dir: Path
    output_dir: Path
    manual_feature_dir: Path
    manual_compound_dir: Path
    release_id: str
    parser_hash: str = field(default_factory=parser_hash)
    created_at_utc: str = field(default_factory=utc_now)
    table_audit: list[dict[str, Any]] = field(default_factory=list)
    feature_import_stats: dict[str, Any] = field(default_factory=dict)
    identity_overlay_import_stats: dict[str, Any] = field(default_factory=dict)

    def writer(self, name: str) -> ParquetWriter:
        return ParquetWriter(self.output_dir / f"{name}.parquet", name)


def parse_xref(raw_xref: str) -> tuple[str, str] | None:
    text = str(raw_xref or "").strip()
    if ":" not in text:
        return None
    namespace, value = text.split(":", 1)
    namespace = namespace.strip().upper()
    value = value.strip()
    if not namespace or not value:
        return None
    return namespace, value


def identifier_rows(
    metabolite_uid: str,
    namespace: str,
    value: str,
    match_field: str,
    rank: float,
    source_table: str,
    source_release: str,
    license_id: str,
    parser_hash_value: str,
) -> Iterator[dict[str, Any]]:
    namespace = str(namespace or "").upper()
    value = str(value or "").strip()
    if not metabolite_uid or not namespace or not value:
        return
    variants = {(namespace, value, normalize_id_token(value))}
    xref_key = f"{namespace}:{value}"
    variants.add(("ID", xref_key, normalize_id_token(xref_key)))
    if namespace == "PUBCHEM.COMPOUND":
        variants.add(("PUBCHEM", value, normalize_id_token(value)))
        variants.add(("CID", value, normalize_id_token(value)))
    if namespace in {"HMDB", "CHEBI", "KEGG", "INCHIKEY"}:
        variants.add((namespace, value, normalize_id_token(value)))
    for ns, raw_value, lookup_key in variants:
        if not lookup_key:
            continue
        yield {
            "lookup_uid": stable_uid("compound_identifier_lookup", metabolite_uid, ns, lookup_key, match_field),
            "metabolite_uid": metabolite_uid,
            "namespace": ns,
            "lookup_key": lookup_key,
            "raw_value": raw_value,
            "match_field": match_field,
            "rank": float(rank),
            "source_table": source_table,
            "source_release": source_release,
            "license_id": license_id,
            "parser_hash": parser_hash_value,
        }


def name_row(
    metabolite_uid: str,
    value: str,
    match_field: str,
    rank: float,
    source_table: str,
    source_release: str,
    license_id: str,
    parser_hash_value: str,
) -> dict[str, Any] | None:
    raw = str(value or "").strip()
    key = normalize_lookup_key(raw)
    if not metabolite_uid or not raw or not key:
        return None
    return {
        "lookup_uid": stable_uid("compound_name_lookup", metabolite_uid, key, match_field),
        "metabolite_uid": metabolite_uid,
        "name_key": key,
        "raw_value": raw,
        "match_field": match_field,
        "rank": float(rank),
        "source_table": source_table,
        "source_release": source_release,
        "license_id": license_id,
        "parser_hash": parser_hash_value,
    }


def inchikey_row(
    metabolite_uid: str,
    value: str,
    match_field: str,
    rank: float,
    source_table: str,
    source_release: str,
    license_id: str,
    parser_hash_value: str,
) -> dict[str, Any] | None:
    raw = str(value or "").strip()
    key = normalize_id_token(raw).upper()
    if not metabolite_uid or not re.match(r"^[A-Z0-9]{14}-[A-Z0-9]{10}-[A-Z0-9]$", key):
        return None
    return {
        "lookup_uid": stable_uid("compound_inchikey_lookup", metabolite_uid, key, match_field),
        "metabolite_uid": metabolite_uid,
        "inchikey": key,
        "connectivity_key": inchikey_connectivity(key),
        "raw_value": raw,
        "match_field": match_field,
        "rank": float(rank),
        "source_table": source_table,
        "source_release": source_release,
        "license_id": license_id,
        "parser_hash": parser_hash_value,
    }


def formula_mass_row(
    metabolite_uid: str,
    formula: str,
    exact_mass: Any,
    match_field: str,
    rank: float,
    source_table: str,
    source_release: str,
    license_id: str,
    parser_hash_value: str,
) -> dict[str, Any] | None:
    mass = maybe_float(exact_mass)
    formula_key = normalize_formula(formula)
    if not metabolite_uid or mass is None or mass <= 0:
        return None
    return {
        "lookup_uid": stable_uid("compound_formula_mass_lookup", metabolite_uid, formula_key, f"{mass:.8f}", match_field),
        "metabolite_uid": metabolite_uid,
        "formula": str(formula or ""),
        "formula_key": formula_key,
        "exact_mass": mass,
        "mass_mda_bucket": int(round(mass * 1000.0)),
        "match_field": match_field,
        "rank": float(rank),
        "source_table": source_table,
        "source_release": source_release,
        "license_id": license_id,
        "parser_hash": parser_hash_value,
    }


def split_overlay_values(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    parts = re.split(r"[|;]", text)
    out: list[str] = []
    seen: set[str] = set()
    for part in parts:
        cleaned = str(part or "").strip()
        key = cleaned.casefold()
        if cleaned and key not in seen:
            seen.add(key)
            out.append(cleaned)
    return out


def iter_identity_overlay_rows(compound_dir: Path) -> Iterator[tuple[Path, dict[str, Any]]]:
    for suffix, reader in (
        (".csv", read_delimited_rows),
        (".tsv", read_delimited_rows),
        (".jsonl", read_jsonl_rows),
    ):
        path = compound_dir / f"compound_identity_overlay{suffix}"
        for row in reader(path):
            yield path, row


def build_from_metabolites(ctx: CompoundIndexContext) -> None:
    started = time.time()
    id_writer = ctx.writer("compound_identifier_index")
    name_writer = ctx.writer("compound_name_index")
    inchikey_writer = ctx.writer("compound_inchikey_index")
    mass_writer = ctx.writer("compound_formula_mass_index")
    id_rows: list[dict[str, Any]] = []
    name_rows: list[dict[str, Any]] = []
    inchikey_rows: list[dict[str, Any]] = []
    mass_rows: list[dict[str, Any]] = []
    indexed_identifier_keys: set[tuple[str, str]] = set()
    overlay_stats: dict[str, Any] = {
        "manual_compound_dir_exists": ctx.manual_compound_dir.exists(),
        "input_rows": 0,
        "imported_rows": 0,
        "imported_identifier_rows": 0,
        "imported_name_rows": 0,
        "imported_inchikey_rows": 0,
        "imported_formula_mass_rows": 0,
        "skipped_existing_cid": 0,
        "skipped_missing_identity": 0,
        "formats": [],
    }

    def flush(force: bool = False) -> None:
        nonlocal id_rows, name_rows, inchikey_rows, mass_rows
        if force or len(id_rows) >= 100_000:
            id_writer.write(id_rows)
            id_rows = []
        if force or len(name_rows) >= 100_000:
            name_writer.write(name_rows)
            name_rows = []
        if force or len(inchikey_rows) >= 100_000:
            inchikey_writer.write(inchikey_rows)
            inchikey_rows = []
        if force or len(mass_rows) >= 100_000:
            mass_writer.write(mass_rows)
            mass_rows = []

    def append_identifier_rows(
        metabolite_uid: str,
        namespace: str,
        value: str,
        match_field: str,
        rank: float,
        source_table: str,
        source_release: str,
        license_id: str,
    ) -> int:
        rows = list(
            identifier_rows(
                metabolite_uid,
                namespace,
                value,
                match_field,
                rank,
                source_table,
                source_release,
                license_id,
                ctx.parser_hash,
            )
        )
        for item in rows:
            indexed_identifier_keys.add((str(item.get("namespace") or "").upper(), str(item.get("lookup_key") or "")))
        id_rows.extend(rows)
        return len(rows)

    metabolites_path = ctx.normalized_dir / "metabolites.parquet"
    if metabolites_path.exists():
        for row in iter_table_rows(
            metabolites_path,
            [
                "metabolite_uid",
                "canonical_name",
                "synonyms",
                "formula",
                "exact_mass",
                "inchikey",
                "external_xrefs",
                "source_release",
                "license_id",
            ],
        ):
            metabolite_uid = row.get("metabolite_uid", "")
            source_release = str(row.get("source_release", "") or ctx.release_id)
            license_id = str(row.get("license_id", "") or "")
            canonical = name_row(
                metabolite_uid,
                row.get("canonical_name", ""),
                "canonical_name",
                100.0,
                "metabolites",
                source_release,
                license_id,
                ctx.parser_hash,
            )
            if canonical:
                name_rows.append(canonical)
            for synonym in row.get("synonyms") or []:
                item = name_row(
                    metabolite_uid,
                    synonym,
                    "synonym",
                    80.0,
                    "metabolites",
                    source_release,
                    license_id,
                    ctx.parser_hash,
                )
                if item:
                    name_rows.append(item)
            inchikey = inchikey_row(
                metabolite_uid,
                row.get("inchikey", ""),
                "inchikey",
                100.0,
                "metabolites",
                source_release,
                license_id,
                ctx.parser_hash,
            )
            if inchikey:
                inchikey_rows.append(inchikey)
            mass = formula_mass_row(
                metabolite_uid,
                row.get("formula", ""),
                row.get("exact_mass"),
                "formula_exact_mass",
                100.0,
                "metabolites",
                source_release,
                license_id,
                ctx.parser_hash,
            )
            if mass:
                mass_rows.append(mass)
            for xref in row.get("external_xrefs") or []:
                parsed = parse_xref(xref)
                if parsed:
                    append_identifier_rows(
                        metabolite_uid,
                        parsed[0],
                        parsed[1],
                        "external_xref",
                        90.0,
                        "metabolites",
                        source_release,
                        license_id,
                    )
            flush()

    xrefs_path = ctx.normalized_dir / "metabolite_xrefs.parquet"
    if xrefs_path.exists():
        for row in iter_table_rows(
            xrefs_path,
            ["metabolite_uid", "xref_source", "xref_id", "source_release", "license_id"],
        ):
            append_identifier_rows(
                row.get("metabolite_uid", ""),
                row.get("xref_source", ""),
                row.get("xref_id", ""),
                "metabolite_xref",
                95.0,
                "metabolite_xrefs",
                str(row.get("source_release", "") or ctx.release_id),
                str(row.get("license_id", "") or ""),
            )
            flush()

    pubchem_path = ctx.pubchem_dir / "cid_properties.parquet"
    if pubchem_path.exists():
        for row in iter_table_rows(
            pubchem_path,
            [
                "pubchem_cid",
                "metabolite_uids",
                "title",
                "formula",
                "exact_mass",
                "inchikey",
                "synonyms",
                "source_release",
                "license_id",
            ],
        ):
            source_release = str(row.get("source_release", "") or ctx.release_id)
            license_id = str(row.get("license_id", "") or "")
            for metabolite_uid in row.get("metabolite_uids") or []:
                append_identifier_rows(
                    metabolite_uid,
                    "PUBCHEM.COMPOUND",
                    row.get("pubchem_cid", ""),
                    "pubchem_cid",
                    95.0,
                    "cid_properties",
                    source_release,
                    license_id,
                )
                title = name_row(
                    metabolite_uid,
                    row.get("title", ""),
                    "pubchem_title",
                    85.0,
                    "cid_properties",
                    source_release,
                    license_id,
                    ctx.parser_hash,
                )
                if title:
                    name_rows.append(title)
                for synonym in row.get("synonyms") or []:
                    item = name_row(
                        metabolite_uid,
                        synonym,
                        "pubchem_synonym",
                        70.0,
                        "cid_properties",
                        source_release,
                        license_id,
                        ctx.parser_hash,
                    )
                    if item:
                        name_rows.append(item)
                inchikey = inchikey_row(
                    metabolite_uid,
                    row.get("inchikey", ""),
                    "pubchem_inchikey",
                    90.0,
                    "cid_properties",
                    source_release,
                    license_id,
                    ctx.parser_hash,
                )
                if inchikey:
                    inchikey_rows.append(inchikey)
                mass = formula_mass_row(
                    metabolite_uid,
                    row.get("formula", ""),
                    row.get("exact_mass"),
                    "pubchem_formula_exact_mass",
                    90.0,
                    "cid_properties",
                    source_release,
                    license_id,
                    ctx.parser_hash,
                )
                if mass:
                    mass_rows.append(mass)
                flush()

    for overlay_path, row in iter_identity_overlay_rows(ctx.manual_compound_dir):
        if overlay_path.name not in overlay_stats["formats"]:
            overlay_stats["formats"].append(overlay_path.name)
        overlay_stats["input_rows"] += 1
        pubchem_cid = clean_identifier_value("CID", row.get("pubchem_cid") or row.get("cid") or row.get("pubchem"))
        cid_key = normalize_id_token(pubchem_cid)
        if cid_key and ("CID", cid_key) in indexed_identifier_keys:
            overlay_stats["skipped_existing_cid"] += 1
            continue
        inchikey_value = row.get("inchikey") or row.get("inchi_key") or ""
        canonical_name = row.get("canonical_name") or row.get("reviewed_name") or row.get("name") or row.get("title") or ""
        metabolite_uid = (
            str(row.get("metabolite_uid") or "").strip()
            or stable_uid("manual_pubchem_metabolite", pubchem_cid or inchikey_value or canonical_name)
        )
        if not (pubchem_cid or inchikey_value or canonical_name):
            overlay_stats["skipped_missing_identity"] += 1
            continue
        source_release = str(row.get("source_release") or ctx.release_id)
        license_id = str(row.get("license_id") or "local:compound_identity_overlay")
        source_table = "manual_compound_identity_overlay"

        if pubchem_cid:
            overlay_stats["imported_identifier_rows"] += append_identifier_rows(
                metabolite_uid,
                "PUBCHEM.COMPOUND",
                pubchem_cid,
                "manual_pubchem_cid",
                99.0,
                source_table,
                source_release,
                license_id,
            )
        for namespace, *fields in (
            ("HMDB", "hmdb_id", "hmdb"),
            ("CHEBI", "chebi_id", "chebi"),
            ("KEGG", "kegg_id", "kegg"),
        ):
            for field_name in fields:
                for value in split_overlay_values(row.get(field_name)):
                    overlay_stats["imported_identifier_rows"] += append_identifier_rows(
                        metabolite_uid,
                        namespace,
                        clean_identifier_value(namespace, value),
                        f"manual_{field_name}",
                        99.0,
                        source_table,
                        source_release,
                        license_id,
                    )
        if inchikey_value:
            overlay_stats["imported_identifier_rows"] += append_identifier_rows(
                metabolite_uid,
                "INCHIKEY",
                str(inchikey_value),
                "manual_inchikey",
                99.0,
                source_table,
                source_release,
                license_id,
            )
            inchikey = inchikey_row(
                metabolite_uid,
                str(inchikey_value),
                "manual_inchikey",
                99.0,
                source_table,
                source_release,
                license_id,
                ctx.parser_hash,
            )
            if inchikey:
                inchikey_rows.append(inchikey)
                overlay_stats["imported_inchikey_rows"] += 1

        canonical = name_row(
            metabolite_uid,
            str(canonical_name),
            "manual_canonical_name",
            99.0,
            source_table,
            source_release,
            license_id,
            ctx.parser_hash,
        )
        if canonical:
            name_rows.append(canonical)
            overlay_stats["imported_name_rows"] += 1
        synonym_values: list[str] = []
        for field_name in ("synonyms", "aliases", "mapped_names", "candidate_names", "iupac_name"):
            synonym_values.extend(split_overlay_values(row.get(field_name)))
        seen_name_keys = {normalize_lookup_key(str(canonical_name))}
        for synonym in synonym_values:
            name_key = normalize_lookup_key(synonym)
            if not name_key or name_key in seen_name_keys:
                continue
            seen_name_keys.add(name_key)
            item = name_row(
                metabolite_uid,
                synonym,
                "manual_synonym",
                85.0,
                source_table,
                source_release,
                license_id,
                ctx.parser_hash,
            )
            if item:
                name_rows.append(item)
                overlay_stats["imported_name_rows"] += 1
        mass = formula_mass_row(
            metabolite_uid,
            row.get("formula") or row.get("molecular_formula") or "",
            row.get("exact_mass") or row.get("monoisotopic_mass"),
            "manual_formula_exact_mass",
            90.0,
            source_table,
            source_release,
            license_id,
            ctx.parser_hash,
        )
        if mass:
            mass_rows.append(mass)
            overlay_stats["imported_formula_mass_rows"] += 1
        overlay_stats["imported_rows"] += 1
        flush()

    flush(force=True)
    for writer in (id_writer, name_writer, inchikey_writer, mass_writer):
        audit = writer.close()
        audit["seconds"] = round(time.time() - started, 3)
        ctx.table_audit.append(audit)
    overlay_stats["formats"] = sorted(overlay_stats["formats"])
    ctx.identity_overlay_import_stats = overlay_stats


def normalize_method_key(value: str) -> str:
    return normalize_lookup_key(value or "default")


def read_csv_rows(path: Path) -> Iterator[dict[str, str]]:
    if not path.exists():
        return
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            yield {str(key or "").strip().casefold().replace(" ", "_"): str(value or "").strip() for key, value in row.items()}


def read_delimited_rows(path: Path) -> Iterator[dict[str, str]]:
    if not path.exists():
        return
    delimiter = "\t" if path.suffix.casefold() == ".tsv" else ","
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            yield {str(key or "").strip().casefold().replace(" ", "_"): str(value or "").strip() for key, value in row.items()}


def read_jsonl_rows(path: Path) -> Iterator[dict[str, Any]]:
    if not path.exists():
        return
    with path.open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            payload = json.loads(line)
            if isinstance(payload, dict):
                yield {str(key or "").strip().casefold().replace(" ", "_"): value for key, value in payload.items()}


def parse_peak_string(value: str) -> list[tuple[float, float]]:
    peaks: list[tuple[float, float]] = []
    text = str(value or "").replace(";", " ").replace(",", " ")
    for token in text.split():
        if not token:
            continue
        if ":" in token:
            mz_text, intensity_text = token.split(":", 1)
        elif "|" in token:
            mz_text, intensity_text = token.split("|", 1)
        else:
            mz_text, intensity_text = token, "1"
        mz = maybe_float(mz_text)
        intensity = maybe_float(intensity_text)
        if mz and mz > 0:
            peaks.append((mz, intensity if intensity is not None and intensity > 0 else 1.0))
    return peaks


def normalize_msp_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value or "").strip().casefold()).strip("_")


def parse_msp_rows(path: Path) -> Iterator[dict[str, Any]]:
    if not path.exists():
        return

    def emit(block: dict[str, Any], peaks: list[tuple[float, float]]) -> dict[str, Any] | None:
        if not block and not peaks:
            return None
        row = {normalize_msp_key(key): value for key, value in block.items()}
        if peaks:
            row["peaks"] = " ".join(f"{mz}:{intensity}" for mz, intensity in peaks)
        aliases = {
            "precursormz": "precursor_mz",
            "precursor_type": "adduct",
            "precursortype": "adduct",
            "ionmode": "ion_mode",
            "db": "source_record_id",
            "database_id": "source_record_id",
            "hmdb": "hmdb_id",
            "chebi": "chebi_id",
            "pubchem": "pubchem_cid",
            "pubchem_cid": "pubchem_cid",
            "kegg": "kegg_id",
        }
        for source, target in aliases.items():
            if source in row and target not in row:
                row[target] = row[source]
        row.setdefault("source_name", path.stem)
        return row

    block: dict[str, Any] = {}
    peaks: list[tuple[float, float]] = []
    in_peaks = False
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                parsed = emit(block, peaks)
                if parsed:
                    yield parsed
                block = {}
                peaks = []
                in_peaks = False
                continue
            if ":" in line and not in_peaks:
                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip()
                block[key] = value
                if normalize_msp_key(key) in {"num_peaks", "numpeaks"}:
                    in_peaks = True
                continue
            parts = line.replace(";", " ").replace(",", " ").split()
            if len(parts) >= 2:
                mz = maybe_float(parts[0])
                intensity = maybe_float(parts[1])
                if mz and mz > 0:
                    peaks.append((mz, intensity if intensity is not None and intensity > 0 else 1.0))
    parsed = emit(block, peaks)
    if parsed:
        yield parsed


def iter_feature_rows(feature_dir: Path, stem: str) -> Iterator[tuple[Path, dict[str, Any]]]:
    for suffix, reader in (
        (".csv", read_delimited_rows),
        (".tsv", read_delimited_rows),
        (".jsonl", read_jsonl_rows),
    ):
        path = feature_dir / f"{stem}{suffix}"
        for row in reader(path):
            yield path, row


def build_feature_uid_resolver(ctx: CompoundIndexContext) -> dict[tuple[str, str], set[str]]:
    resolver: dict[tuple[str, str], set[str]] = {}

    def add(namespace: str, value: str, metabolite_uid: str) -> None:
        namespace = str(namespace or "").upper()
        key = normalize_id_token(value)
        if namespace and key and metabolite_uid:
            resolver.setdefault((namespace, key), set()).add(str(metabolite_uid))

    id_path = ctx.output_dir / "compound_identifier_index.parquet"
    if id_path.exists():
        for row in iter_table_rows(id_path, ["namespace", "lookup_key", "metabolite_uid"]):
            add(row.get("namespace", ""), row.get("lookup_key", ""), row.get("metabolite_uid", ""))

    inchikey_path = ctx.output_dir / "compound_inchikey_index.parquet"
    if inchikey_path.exists():
        for row in iter_table_rows(inchikey_path, ["inchikey", "metabolite_uid"]):
            add("INCHIKEY", row.get("inchikey", ""), row.get("metabolite_uid", ""))
    return resolver


def clean_identifier_value(namespace: str, value: Any) -> str:
    text = str(value or "").strip()
    if ":" not in text:
        return text
    prefix, suffix = text.split(":", 1)
    prefix_key = normalize_id_token(prefix).upper()
    namespace_key = normalize_id_token(namespace).upper()
    aliases = {
        "PUBCHEM": {"PUBCHEM", "PUBCHEMCOMPOUND", "CID"},
        "CID": {"PUBCHEM", "PUBCHEMCOMPOUND", "CID"},
        "CHEBI": {"CHEBI"},
        "HMDB": {"HMDB"},
        "KEGG": {"KEGG", "KEGGCOMPOUND"},
        "INCHIKEY": {"INCHIKEY", "INCHI"},
    }
    if prefix_key in aliases.get(namespace_key, {namespace_key}):
        return suffix.strip()
    return text


def feature_metabolite_uids(row: dict[str, Any], resolver: dict[tuple[str, str], set[str]]) -> list[str]:
    direct = str(row.get("metabolite_uid") or "").strip()
    if direct:
        return [direct]
    candidates: set[str] = set()
    fields = [
        ("hmdb_id", "HMDB"),
        ("hmdb", "HMDB"),
        ("chebi_id", "CHEBI"),
        ("chebi", "CHEBI"),
        ("pubchem_cid", "CID"),
        ("cid", "CID"),
        ("pubchem", "CID"),
        ("kegg_id", "KEGG"),
        ("kegg", "KEGG"),
        ("inchikey", "INCHIKEY"),
        ("inchi_key", "INCHIKEY"),
    ]
    for field_name, namespace in fields:
        value = row.get(field_name)
        if value in {None, ""}:
            continue
        key = normalize_id_token(clean_identifier_value(namespace, value))
        candidates.update(resolver.get((namespace, key), set()))
    return sorted(candidates)


def build_optional_feature_indexes(ctx: CompoundIndexContext) -> None:
    started = time.time()
    rt_writer = ctx.writer("compound_rt_index")
    ms2_writer = ctx.writer("compound_ms2_index")
    rt_rows: list[dict[str, Any]] = []
    ms2_rows: list[dict[str, Any]] = []
    feature_dir = ctx.manual_feature_dir
    resolver = build_feature_uid_resolver(ctx)
    stats = {
        "feature_dir_exists": feature_dir.exists(),
        "rt_input_rows": 0,
        "rt_imported_rows": 0,
        "ms2_input_records": 0,
        "ms2_imported_peaks": 0,
        "skipped_unresolved": 0,
        "skipped_ambiguous": 0,
        "formats": [],
    }

    for rt_path, row in iter_feature_rows(feature_dir, "compound_rt"):
        if rt_path.name not in stats["formats"]:
            stats["formats"].append(rt_path.name)
        stats["rt_input_rows"] += 1
        metabolite_uids = feature_metabolite_uids(row, resolver)
        if not metabolite_uids:
            stats["skipped_unresolved"] += 1
            continue
        if len(metabolite_uids) > 1:
            stats["skipped_ambiguous"] += 1
            continue
        for metabolite_uid in metabolite_uids:
            rt = maybe_float(row.get("rt") or row.get("retention_time"))
            if not metabolite_uid or rt is None:
                continue
            method_key = normalize_method_key(row.get("method_key") or row.get("method") or row.get("lc_method") or "default")
            source_record_id = row.get("source_record_id") or stable_uid("rt_source", metabolite_uid, method_key, f"{rt:.6f}")
            rt_rows.append(
                {
                    "lookup_uid": stable_uid("compound_rt_lookup", metabolite_uid, method_key, f"{rt:.6f}", source_record_id),
                    "metabolite_uid": metabolite_uid,
                    "rt": rt,
                    "rt_unit": row.get("rt_unit") or "min",
                    "method_key": method_key,
                    "source_name": row.get("source_name") or "manual_compound_features",
                    "source_record_id": source_record_id,
                    "source_release": row.get("source_release") or ctx.release_id,
                    "license_id": row.get("license_id") or "local:compound_features",
                    "parser_hash": ctx.parser_hash,
                }
            )
            stats["rt_imported_rows"] += 1

    ms2_rows_by_source: list[tuple[Path, dict[str, Any]]] = list(iter_feature_rows(feature_dir, "compound_ms2"))
    for path in sorted(set([feature_dir / "compound_ms2.msp", *feature_dir.glob("*.msp")])):
        if path.exists():
            ms2_rows_by_source.extend((path, row) for row in parse_msp_rows(path))

    for ms2_path, row in ms2_rows_by_source:
        if ms2_path.name not in stats["formats"]:
            stats["formats"].append(ms2_path.name)
        stats["ms2_input_records"] += 1
        metabolite_uids = feature_metabolite_uids(row, resolver)
        if not metabolite_uids:
            stats["skipped_unresolved"] += 1
            continue
        if len(metabolite_uids) > 1:
            stats["skipped_ambiguous"] += 1
            continue
        metabolite_uid = metabolite_uids[0]
        precursor_mz = maybe_float(row.get("precursor_mz"))
        adduct = row.get("adduct", "")
        ion_mode = row.get("ion_mode") or row.get("polarity") or ""
        source_record_id = row.get("source_record_id") or stable_uid("ms2_source", metabolite_uid, precursor_mz, adduct, row.get("peaks", ""))
        peaks: list[tuple[float, float]] = []
        fragment_mz = maybe_float(row.get("fragment_mz"))
        if fragment_mz:
            peaks.append((fragment_mz, maybe_float(row.get("fragment_intensity")) or 1.0))
        peaks.extend(parse_peak_string(row.get("peaks", "")))
        for mz, intensity in peaks:
            if not metabolite_uid or mz <= 0:
                continue
            ms2_rows.append(
                {
                    "lookup_uid": stable_uid("compound_ms2_lookup", metabolite_uid, source_record_id, f"{mz:.6f}"),
                    "metabolite_uid": metabolite_uid,
                    "precursor_mz": precursor_mz,
                    "adduct": adduct,
                    "ion_mode": ion_mode,
                    "fragment_mz": mz,
                    "fragment_intensity": intensity,
                    "source_name": row.get("source_name") or "manual_compound_features",
                    "source_record_id": source_record_id,
                    "source_release": row.get("source_release") or ctx.release_id,
                    "license_id": row.get("license_id") or "local:compound_features",
                    "parser_hash": ctx.parser_hash,
                }
            )
            stats["ms2_imported_peaks"] += 1

    rt_writer.write(rt_rows)
    ms2_writer.write(ms2_rows)
    for writer in (rt_writer, ms2_writer):
        audit = writer.close()
        audit["seconds"] = round(time.time() - started, 3)
        ctx.table_audit.append(audit)
    stats["formats"] = sorted(stats["formats"])
    ctx.feature_import_stats = stats


def build_compound_match_index(
    normalized_root: Path,
    pubchem_root: Path,
    output_root: Path,
    release_id: str | None = None,
    manual_feature_root: Path | None = None,
    manual_compound_root: Path | None = None,
) -> Path:
    require_arrow()
    resolved_release = release_id or latest_normalized_release_id(normalized_root)
    normalized_dir = normalized_root / resolved_release
    pubchem_dir = pubchem_root / resolved_release
    output_dir = output_root / resolved_release
    manual_feature_dir = (manual_feature_root / resolved_release) if manual_feature_root else (output_dir / "__no_manual_features__")
    manual_compound_dir = (manual_compound_root / resolved_release) if manual_compound_root else (output_dir / "__no_manual_compounds__")
    output_dir.mkdir(parents=True, exist_ok=True)
    ctx = CompoundIndexContext(
        normalized_dir=normalized_dir,
        pubchem_dir=pubchem_dir,
        output_dir=output_dir,
        manual_feature_dir=manual_feature_dir,
        manual_compound_dir=manual_compound_dir,
        release_id=resolved_release,
    )
    build_from_metabolites(ctx)
    build_optional_feature_indexes(ctx)
    manifest = {
        "release_id": resolved_release,
        "created_at_utc": ctx.created_at_utc,
        "normalized_dir": str(normalized_dir),
        "pubchem_dir": str(pubchem_dir),
        "manual_feature_dir": str(manual_feature_dir),
        "manual_compound_dir": str(manual_compound_dir),
        "output_dir": str(output_dir),
        "parser_hash": ctx.parser_hash,
        "identity_overlay_import": ctx.identity_overlay_import_stats,
        "feature_import": ctx.feature_import_stats,
        "tables": ctx.table_audit,
        "notes": [
            "Compound match index is derived from normalized_store and PubChem CID cache.",
            "Identity resolution remains deterministic and abstains on ambiguous candidates.",
            "Reviewed compound identity overlays can add confirmed PubChem/HMDB/ChEBI/InChIKey identities that are not yet in the normalized graph.",
            "Optional RT/MS2 features can be imported from compound_rt.csv/.tsv/.jsonl and compound_ms2.csv/.tsv/.jsonl/.msp under manual_feature_dir.",
        ],
    }
    manifest_path = output_dir / "compound_match_index_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, sort_keys=True, indent=2), encoding="utf-8")
    return manifest_path


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build compound-aware matching indexes.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--normalized-root", default=DEFAULT_NORMALIZED_ROOT)
    parser.add_argument("--pubchem-root", default=DEFAULT_PUBCHEM_ROOT)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--manual-feature-root", default=DEFAULT_MANUAL_FEATURE_ROOT)
    parser.add_argument("--manual-compound-root", default=DEFAULT_MANUAL_COMPOUND_ROOT)
    parser.add_argument("--release-id", default="")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    workspace = Path(args.workspace).resolve()
    manifest_path = build_compound_match_index(
        normalized_root=(workspace / args.normalized_root).resolve(),
        pubchem_root=(workspace / args.pubchem_root).resolve(),
        output_root=(workspace / args.output_root).resolve(),
        release_id=args.release_id or None,
        manual_feature_root=(workspace / args.manual_feature_root).resolve(),
        manual_compound_root=(workspace / args.manual_compound_root).resolve(),
    )
    print(json.dumps({"compound_match_index_manifest": str(manifest_path)}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
