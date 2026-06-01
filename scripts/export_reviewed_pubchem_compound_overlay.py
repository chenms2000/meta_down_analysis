"""Export reviewed European PubChem identities into a local compound overlay.

The overlay is intentionally narrow: it promotes only strict, reviewed PubChem
identities that were absent from the local compound index at curation time.
The compound match index builder can import this CSV without pretending these
records already have full graph/pathway connectivity.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any


DEFAULT_REVIEW_PATH = "manual_sources/european_trait_identity/european_trait_identity_review.csv"
DEFAULT_OUTPUT_ROOT = "manual_sources/compound_identity_overlays"
DEFAULT_CACHE_PATH = "raw_lake/European/pubchem_identity_autofill_cache.json"
ACCEPTED_REVIEW_STATUSES = {"reviewed", "auto_reviewed", "accepted", "confirmed"}


def normalize_key(value: Any) -> str:
    return str(value or "").strip().casefold().replace(" ", "_")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return [
            {normalize_key(key): str(value or "").strip() for key, value in row.items()}
            for row in reader
        ]


def split_values(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    values: list[str] = []
    seen: set[str] = set()
    for part in re.split(r"[|;]", text):
        cleaned = part.strip()
        key = cleaned.casefold()
        if cleaned and key not in seen:
            seen.add(key)
            values.append(cleaned)
    return values


def normalize_cid(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    match = re.search(r"(\d+)", text)
    return match.group(1) if match else ""


def load_cache(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def cache_bundle(cache: dict[str, Any], cid: str) -> dict[str, Any]:
    bundle = cache.get(f"cid_bundle::{cid}") or {}
    return bundle if isinstance(bundle, dict) else {}


def bundle_properties(bundle: dict[str, Any]) -> dict[str, Any]:
    properties = bundle.get("properties") or {}
    return properties if isinstance(properties, dict) else {}


def bundle_synonyms(bundle: dict[str, Any]) -> list[str]:
    synonyms = bundle.get("synonyms") or []
    return [str(item or "").strip() for item in synonyms if str(item or "").strip()]


def extract_external_ids(values: list[str]) -> tuple[list[str], list[str], list[str]]:
    hmdb: list[str] = []
    chebi: list[str] = []
    kegg: list[str] = []
    for value in values:
        upper = value.upper()
        if re.fullmatch(r"HMDB\d{7}", upper):
            hmdb.append(upper)
        elif upper.startswith("HMDB:"):
            hmdb.append(upper.split(":", 1)[1])
        elif upper.startswith("CHEBI:"):
            chebi.append(upper.split(":", 1)[1])
        elif upper.startswith("KEGG:"):
            kegg.append(upper.split(":", 1)[1])
    return sorted(set(hmdb)), sorted(set(chebi)), sorted(set(kegg))


def is_name_like(value: str) -> bool:
    text = str(value or "").strip()
    if not text:
        return False
    upper = value.upper()
    if ":" in value:
        return False
    if re.fullmatch(r"HMDB\d{7}", upper):
        return False
    if upper.startswith(("HMDB:", "CHEBI:", "KEGG:", "PUBCHEM:", "CID:")):
        return False
    if re.fullmatch(r"\d+", value):
        return False
    if re.fullmatch(r"\d{2,7}-\d{2}-\d", value):
        return False
    if re.fullmatch(r"[A-Z0-9]{14}-[A-Z0-9]{10}-[A-Z0-9]", upper):
        return False
    registry_prefixes = (
        "AKOS",
        "AS-",
        "BDBM",
        "CAS-",
        "CCG-",
        "CHEMBL",
        "DB-",
        "DTXCID",
        "DTXSID",
        "EBC-",
        "EN",
        "GTPL",
        "HMS",
        "HY-",
        "MFCD",
        "NCGC",
        "NSC",
        "NS0",
        "Q",
        "SBB",
        "SCHEMBL",
        "T3D",
        "Tox21".upper(),
        "UNII-",
    )
    if upper.startswith(registry_prefixes):
        return False
    if re.fullmatch(r"[A-Z0-9_.-]{4,}", value) and not re.search(r"[a-z]", value):
        return False
    return True


def append_unique(target: list[str], values: list[str], limit: int = 80) -> None:
    seen = {value.casefold() for value in target}
    for value in values:
        cleaned = str(value or "").strip()
        key = cleaned.casefold()
        if cleaned and key not in seen:
            target.append(cleaned)
            seen.add(key)
        if len(target) >= limit:
            return


def export_overlay(
    review_path: Path,
    output_path: Path,
    pubchem_cache_path: Path,
    release_id: str,
    include_all_strict: bool = False,
) -> dict[str, Any]:
    rows = read_csv_rows(review_path)
    cache = load_cache(pubchem_cache_path)
    grouped: dict[str, dict[str, Any]] = {}
    stats = {
        "review_rows": len(rows),
        "eligible_rows": 0,
        "written_rows": 0,
        "skipped_non_strict": 0,
        "skipped_not_absent": 0,
        "skipped_no_cid": 0,
        "output_path": str(output_path),
    }

    for row in rows:
        if row.get("review_status") not in ACCEPTED_REVIEW_STATUSES or row.get("identity_scope") != "strict_identity":
            stats["skipped_non_strict"] += 1
            continue
        if not include_all_strict and row.get("local_cid_status") != "local_cid_absent":
            stats["skipped_not_absent"] += 1
            continue
        cids = [normalize_cid(value) for value in split_values(row.get("reviewed_pubchem_cid") or row.get("pubchem_cids"))]
        cids = [cid for cid in cids if cid]
        if not cids:
            stats["skipped_no_cid"] += 1
            continue
        stats["eligible_rows"] += 1
        for cid in cids:
            bundle = cache_bundle(cache, cid)
            properties = bundle_properties(bundle)
            pubchem_synonyms = bundle_synonyms(bundle)
            record = grouped.setdefault(
                cid,
                {
                    "metabolite_uid": f"manual_pubchem_cid_{cid}",
                    "canonical_name": "",
                    "synonyms": [],
                    "pubchem_cid": cid,
                    "inchikey": "",
                    "formula": "",
                    "exact_mass": "",
                    "iupac_name": "",
                    "hmdb_id": [],
                    "chebi_id": [],
                    "kegg_id": [],
                    "source_name": "European_identity_review",
                    "source_record_id": f"European_identity_review:{cid}",
                    "source_release": release_id,
                    "license_id": "local:european_trait_identity_review",
                    "evidence_source": [],
                    "evidence_url": [],
                    "accession_ids": [],
                    "notes": [],
                },
            )
            canonical = (
                row.get("reviewed_name")
                or properties.get("Title")
                or row.get("component_name")
                or row.get("reported_trait")
                or ""
            )
            if canonical and not record["canonical_name"]:
                record["canonical_name"] = canonical
            record["inchikey"] = record["inchikey"] or row.get("reviewed_inchikey") or properties.get("InChIKey") or ""
            record["formula"] = record["formula"] or properties.get("MolecularFormula") or ""
            record["iupac_name"] = record["iupac_name"] or properties.get("IUPACName") or ""
            synonym_candidates: list[str] = []
            for field_name in ("reviewed_name", "component_name", "mapped_names", "candidate_names", "pubchem_titles"):
                synonym_candidates.extend(split_values(row.get(field_name)))
            synonym_candidates.extend([str(properties.get("Title") or ""), str(properties.get("IUPACName") or "")])
            synonym_candidates.extend(pubchem_synonyms)
            hmdb_ids, chebi_ids, kegg_ids = extract_external_ids(synonym_candidates)
            append_unique(record["hmdb_id"], hmdb_ids)
            append_unique(record["chebi_id"], chebi_ids)
            append_unique(record["kegg_id"], kegg_ids)
            append_unique(record["synonyms"], [value for value in synonym_candidates if is_name_like(value)], limit=40)
            append_unique(record["evidence_source"], split_values(row.get("evidence_source")))
            append_unique(record["evidence_url"], split_values(row.get("evidence_url")))
            append_unique(record["accession_ids"], split_values(row.get("accession_id")))
            append_unique(record["notes"], split_values(row.get("notes")), limit=20)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "metabolite_uid",
        "canonical_name",
        "synonyms",
        "pubchem_cid",
        "inchikey",
        "formula",
        "exact_mass",
        "iupac_name",
        "hmdb_id",
        "chebi_id",
        "kegg_id",
        "source_name",
        "source_record_id",
        "source_release",
        "license_id",
        "evidence_source",
        "evidence_url",
        "accession_ids",
        "notes",
    ]
    with output_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for cid in sorted(grouped, key=lambda value: int(value)):
            record = grouped[cid]
            writer.writerow(
                {
                    key: "|".join(value) if isinstance(value, list) else value
                    for key, value in record.items()
                }
            )
            stats["written_rows"] += 1
    return stats


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export reviewed PubChem identities into compound_identity_overlay.csv.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", default="")
    parser.add_argument("--review-path", default=DEFAULT_REVIEW_PATH)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--pubchem-cache-path", default=DEFAULT_CACHE_PATH)
    parser.add_argument("--include-all-strict", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    workspace = Path(args.workspace).resolve()
    release_id = args.release_id or "mvp_20260513T002254"
    stats = export_overlay(
        review_path=(workspace / args.review_path).resolve(),
        output_path=(workspace / args.output_root / release_id / "compound_identity_overlay.csv").resolve(),
        pubchem_cache_path=(workspace / args.pubchem_cache_path).resolve(),
        release_id=release_id,
        include_all_strict=bool(args.include_all_strict),
    )
    print(json.dumps(stats, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
