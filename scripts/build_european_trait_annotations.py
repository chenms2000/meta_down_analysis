"""Build European GWAS trait annotation hints for metabolite resolution.

The input European.csv uses GWAS accessions (GCST...) and reportedTrait strings
that often need lightweight normalization before the local compound resolver can
use them. This script creates a small reusable annotation CSV, optionally using
PubChem PUG REST for names that are not already covered by the local compound
name index.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import pyarrow.dataset as ds
except Exception:  # pragma: no cover
    ds = None

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from build_graph_projection import normalize_id_token, normalize_lookup_key  # noqa: E402
from metabo_service import (  # noqa: E402
    clean_reported_trait_name,
    normalize_gwas_accession,
    reported_trait_candidate_names,
)


DEFAULT_RELEASE_ID = "mvp_20260513T002254"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
REVIEW_PRESERVE_FIELDS = [
    "review_status",
    "identity_scope",
    "reviewed_name",
    "reviewed_pubchem_cid",
    "reviewed_hmdb_id",
    "reviewed_chebi_id",
    "reviewed_inchikey",
    "evidence_source",
    "evidence_url",
]
IDENTITY_REVIEW_STATUSES = {"reviewed", "auto_reviewed", "accepted", "confirmed"}
EUROPEAN_SOURCE_METADATA_FIELDS = [
    "summary_statistics_url",
    "pubmed_id",
    "paper_title",
    "journal",
    "publication_date",
    "efo_traits",
    "bg_traits",
    "initial_sample_description",
    "discovery_sample_ancestry",
]
BIOCHEMICAL_ABBREVIATIONS = {
    "7-hoca": ["7-alpha-hydroxy-3-oxo-4-cholestenoic acid", "7-alpha-hydroxy-3-oxo-4-cholestenoate"],
    "13-hode": ["13-hydroxyoctadecadienoic acid", "13-hydroxyoctadecadienoate"],
    "9-hode": ["9-hydroxyoctadecadienoic acid", "9-hydroxyoctadecadienoate"],
    "hode": ["hydroxyoctadecadienoic acid", "hydroxyoctadecadienoate"],
    "hete": ["hydroxyeicosatetraenoic acid", "hydroxyeicosatetraenoate"],
    "dpa": ["docosapentaenoic acid", "docosapentaenoate"],
    "dha": ["docosahexaenoic acid", "docosahexaenoate"],
    "epa": ["eicosapentaenoic acid", "eicosapentaenoate"],
    "ida": ["iminodiacetic acid", "iminodiacetate"],
}


def dedupe(values: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        text = re.sub(r"\s+", " ", str(value or "").strip())
        if not text:
            continue
        key = text.casefold()
        if key in seen:
            continue
        seen.add(key)
        out.append(text)
    return out


def acid_variants(name: str) -> list[str]:
    variants: list[str] = []
    replacements = {
        "docosahexaenoate": "docosahexaenoic acid",
        "eicosapentaenoate": "eicosapentaenoic acid",
        "arachidonate": "arachidonic acid",
        "linoleate": "linoleic acid",
        "linolenate": "linolenic acid",
        "oleate": "oleic acid",
        "palmitate": "palmitic acid",
        "stearate": "stearic acid",
        "myristate": "myristic acid",
        "laurate": "lauric acid",
        "hydrocinnamate": "hydrocinnamic acid",
        "tartarate": "tartaric acid",
        "tartronate": "tartronic acid",
    }
    lowered = name.casefold()
    for src, dst in replacements.items():
        if src in lowered:
            variants.append(re.sub(src, dst, name, flags=re.IGNORECASE))
    if lowered.endswith("ate") and len(name) > 5:
        variants.append(f"{name[:-3]}ic acid")
    return variants


def lipid_abbreviation_variants(name: str) -> list[str]:
    variants: list[str] = []
    text = name.strip()
    replacements = {
        "gpc": "glycerophosphocholine",
        "gps": "glycerophosphoserine",
        "gpi": "glycerophosphoinositol",
        "gpe": "glycerophosphoethanolamine",
        "gpg": "glycerophosphoglycerol",
    }
    for src, dst in replacements.items():
        if re.search(rf"\b{src}\b", text, flags=re.IGNORECASE):
            variants.append(re.sub(rf"\b{src}\b", dst, text, flags=re.IGNORECASE))
    return variants


def split_outside_parentheses_any(text: str, separators: list[str]) -> list[str]:
    parts: list[str] = []
    current: list[str] = []
    depth = 0
    index = 0
    lowered = text.casefold()
    while index < len(text):
        char = text[index]
        if char == "(":
            depth += 1
            current.append(char)
            index += 1
            continue
        if char == ")" and depth:
            depth -= 1
            current.append(char)
            index += 1
            continue
        matched = ""
        if depth == 0:
            for separator in separators:
                if lowered.startswith(separator.casefold(), index):
                    matched = separator
                    break
        if matched:
            part = "".join(current).strip()
            if part:
                parts.append(part)
            current = []
            index += len(matched)
            continue
        current.append(char)
        index += 1
    part = "".join(current).strip()
    if part:
        parts.append(part)
    return parts or [text]


def composite_trait_variants(name: str) -> list[str]:
    text = name.strip()
    variants: list[str] = []
    if re.search(r"\s+\+\s+|/", text):
        for part in split_outside_parentheses_any(text, [" + ", "+", "/"]):
            cleaned = clean_reported_trait_name(part)
            if cleaned and cleaned != text:
                variants.append(cleaned)
    return variants


def biochemical_abbreviation_variants(name: str) -> list[str]:
    variants: list[str] = []
    text = name.strip()
    lowered = text.casefold()
    for token, expansions in BIOCHEMICAL_ABBREVIATIONS.items():
        if re.search(rf"(?<![a-z0-9]){re.escape(token)}(?![a-z0-9])", lowered):
            variants.extend(expansions)
            for expansion in expansions:
                variants.append(re.sub(re.escape(token), expansion, text, flags=re.IGNORECASE))
    return variants


def stereochemistry_relaxed_variants(name: str) -> list[str]:
    variants: list[str] = []
    text = name.strip()
    relaxed = re.sub(r"(?<![A-Za-z0-9])\(?[RSZE]\)?\s*-\s*", "", text, flags=re.IGNORECASE)
    relaxed = re.sub(r"\b\d+[RSZE](?:,\d+[RSZE])*[-,]?", "", relaxed, flags=re.IGNORECASE)
    relaxed = re.sub(r"\s+", " ", relaxed).strip(" ,-")
    if relaxed and relaxed != text:
        variants.append(relaxed)
    bracket_choice = re.sub(r"\s*\[[^\]]*\bor\b[^\]]*\]", "", text, flags=re.IGNORECASE).strip()
    if bracket_choice and bracket_choice != text:
        variants.append(bracket_choice)
    return variants


def expanded_trait_candidates(reported_trait: Any) -> list[str]:
    seeds = reported_trait_candidate_names(reported_trait)
    variants: list[str] = []
    for seed in seeds:
        variants.append(seed)
        no_semicolon = seed.split(";", 1)[0].strip()
        if no_semicolon != seed:
            variants.append(no_semicolon)
        no_number_suffix = re.sub(r"\s+\(\d+\)$", "", seed).strip()
        if no_number_suffix != seed:
            variants.append(no_number_suffix)
        no_parentheses = clean_reported_trait_name(re.sub(r"\([^)]*\)", "", seed))
        if no_parentheses != seed:
            variants.append(no_parentheses)
        if re.search(r"\bDHA\b", seed):
            variants.extend(["DHA", "docosahexaenoic acid"])
        if re.search(r"\bEPA\b", seed):
            variants.extend(["EPA", "eicosapentaenoic acid"])
        composite_variants = composite_trait_variants(seed)
        variants.extend(composite_variants)
        variants.extend(acid_variants(seed))
        variants.extend(lipid_abbreviation_variants(seed))
        variants.extend(biochemical_abbreviation_variants(seed))
        variants.extend(stereochemistry_relaxed_variants(seed))
        for composite in composite_variants:
            variants.extend(acid_variants(composite))
            variants.extend(lipid_abbreviation_variants(composite))
            variants.extend(biochemical_abbreviation_variants(composite))
            variants.extend(stereochemistry_relaxed_variants(composite))
    return dedupe(variants)


def unresolved_review_bucket(reported_trait: Any, candidate_names: Any = "") -> tuple[str, str, str]:
    text = clean_reported_trait_name(reported_trait)
    lowered = text.casefold()
    candidates = str(candidate_names or "")
    if re.search(r"\bX[- ]?\d+\b", text, flags=re.IGNORECASE):
        return (
            "platform_x_code",
            "manual_source_mapping_required",
            "Platform-internal metabolite code; do not infer chemical identity without the original study supplement or vendor mapping.",
        )
    if re.search(r"bilirubin degradation product.*C\d+H\d+|(?:conjugate|glucuronide) of C\d+H\d+", text, flags=re.IGNORECASE):
        return (
            "formula_or_unknown_conjugate",
            "review_queue_only",
            "Formula/class description is not enough to assign a unique structure.",
        )
    if re.search(r"\s+\+\s+|/|\bto\b|\bratio\b", text, flags=re.IGNORECASE):
        return (
            "composite_or_ratio",
            "component_review",
            "Trait may contain multiple components or a ratio; resolve components separately and preserve GCST-level provenance.",
        )
    if re.search(
        r"\b(?:gpc|gpe|gpi|gps|gpg)\b|\d{1,2}:\d|sphingo|ceramide|stearoyl|oleoyl|linoleoyl|palmitoyl|arachidon",
        lowered,
    ):
        return (
            "lipid_shorthand",
            "lipid_identifier_review",
            "Lipid shorthand may be chain-specific or isomeric; prefer LIPID MAPS/HMDB/InChIKey before promotion.",
        )
    if re.search(r"carnitine|\bC\d+-DC\b|\(C\d", text, flags=re.IGNORECASE):
        return (
            "acylcarnitine_shorthand",
            "identifier_or_alias_review",
            "Acylcarnitine shorthand needs identifier or curated synonym confirmation.",
        )
    if re.search(r"\(\d+\)|\bor\b|\[[^\]]*\bor[^\]]*\]", text, flags=re.IGNORECASE):
        return (
            "positional_or_isomer_ambiguous",
            "review_queue_only",
            "Positional/isomer label is ambiguous; keep as review item unless an orthogonal identifier is available.",
        )
    if re.search(r"\b(?:hode|hete|hoca|dpa|dha|epa|ida)\b", candidates, flags=re.IGNORECASE):
        return (
            "biochemical_abbrev",
            "alias_rule_candidate",
            "Biochemical abbreviation expanded but still did not match locally.",
        )
    if re.search(r"ate\b|acid\b|sulfate|glucuronide|glycine|lysine|proline", text, flags=re.IGNORECASE):
        return (
            "alias_variant_missing",
            "alias_rule_or_manual_identifier",
            "Likely synonym/ionization/stereochemistry gap; review before adding deterministic alias.",
        )
    return (
        "manual_review",
        "review_queue_only",
        "No safe automatic rule identified.",
    )


def write_unresolved_review(path: Path, output_rows: list[dict[str, Any]], existing_review: dict[str, dict[str, Any]]) -> dict[str, int]:
    review_rows: list[dict[str, Any]] = []
    counts: dict[str, int] = {}
    for row in output_rows:
        if row.get("resolution_status") != "unresolved":
            continue
        accession = str(row.get("accession_id", "") or "")
        preserved = existing_review.get(normalize_gwas_accession(accession), {})
        bucket, action, note = unresolved_review_bucket(row.get("reported_trait"), row.get("candidate_names"))
        counts[bucket] = counts.get(bucket, 0) + 1
        review_rows.append(
            {
                "accession_id": accession,
                "reported_trait": row.get("reported_trait", ""),
                "summary_statistics_url": row.get("summary_statistics_url", ""),
                "pubmed_id": row.get("pubmed_id", ""),
                "paper_title": row.get("paper_title", ""),
                "review_bucket": bucket,
                "recommended_action": action,
                "review_status": preserved.get("review_status") or "needs_review",
                "identity_scope": preserved.get("identity_scope", ""),
                "candidate_names": row.get("candidate_names", ""),
                "mapped_names": row.get("mapped_names", ""),
                "pubchem_query": row.get("pubchem_query", ""),
                "notes": note,
                "reviewed_name": preserved.get("reviewed_name", ""),
                "reviewed_pubchem_cid": preserved.get("reviewed_pubchem_cid", ""),
                "reviewed_hmdb_id": preserved.get("reviewed_hmdb_id", ""),
                "reviewed_chebi_id": preserved.get("reviewed_chebi_id", ""),
                "reviewed_inchikey": preserved.get("reviewed_inchikey", ""),
                "evidence_source": preserved.get("evidence_source", ""),
                "evidence_url": preserved.get("evidence_url", ""),
            }
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "accession_id",
            "reported_trait",
            "summary_statistics_url",
            "pubmed_id",
            "paper_title",
            "review_bucket",
            "recommended_action",
            "review_status",
            "identity_scope",
            "candidate_names",
            "mapped_names",
            "pubchem_query",
            "notes",
            "reviewed_name",
            "reviewed_pubchem_cid",
            "reviewed_hmdb_id",
            "reviewed_chebi_id",
            "reviewed_inchikey",
            "evidence_source",
            "evidence_url",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(review_rows)
    return counts


def read_european_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        counts: dict[str, int] = {}
        unique_header: list[str] = []
        for raw in header:
            name = str(raw or "").strip() or "unnamed"
            count = counts.get(name, 0)
            counts[name] = count + 1
            unique_header.append(name if count == 0 else f"{name}.{count}")
        rows = []
        for row in reader:
            if not any(str(cell).strip() for cell in row):
                continue
            item = {unique_header[index]: str(cell or "").strip() for index, cell in enumerate(row) if index < len(unique_header)}
            rows.append(item)
        return rows


def load_parquet_column(path: Path, column: str) -> list[Any]:
    if ds is None or not path.exists():
        return []
    return ds.dataset(str(path), format="parquet").to_table(columns=[column]).column(column).to_pylist()


def load_local_indexes(compound_dir: Path) -> tuple[set[str], set[str]]:
    name_keys = set(str(value or "") for value in load_parquet_column(compound_dir / "compound_name_index.parquet", "name_key"))
    identifier_path = compound_dir / "compound_identifier_index.parquet"
    local_cids: set[str] = set()
    if ds is not None and identifier_path.exists():
        table = ds.dataset(str(identifier_path), format="parquet").to_table(columns=["namespace", "lookup_key"])
        for namespace, lookup_key in zip(table.column("namespace").to_pylist(), table.column("lookup_key").to_pylist()):
            if str(namespace or "").upper() == "CID":
                local_cids.add(normalize_id_token(lookup_key))
    return name_keys, local_cids


def load_cache(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def save_cache(path: Path, cache: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cache, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def split_review_values(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    delimiter = "|" if "|" in text else ";"
    return [part.strip() for part in text.split(delimiter) if part.strip()]


def load_existing_review(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    review: dict[str, dict[str, Any]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            accession = normalize_gwas_accession(row.get("accession_id"))
            if not accession:
                continue
            preserved = {field: str(row.get(field) or "").strip() for field in REVIEW_PRESERVE_FIELDS}
            if any(preserved.values()):
                review[accession] = preserved
    return review


def load_identity_review(path: Path) -> dict[str, list[dict[str, Any]]]:
    if not path.exists():
        return {}
    review: dict[str, list[dict[str, Any]]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            accession = normalize_gwas_accession(row.get("accession_id"))
            if not accession:
                continue
            item = {field: str(row.get(field) or "").strip() for field in REVIEW_PRESERVE_FIELDS}
            item["component_name"] = str(row.get("component_name") or "").strip()
            item["relationship_to_trait"] = str(row.get("relationship_to_trait") or "").strip()
            if any(item.values()):
                review.setdefault(accession, []).append(item)
    return review


def reviewed_identity_rows(rows: list[dict[str, Any]], identity_scope: str | None = None) -> list[dict[str, Any]]:
    reviewed = [
        row
        for row in rows
        if str(row.get("review_status") or "").strip().casefold() in IDENTITY_REVIEW_STATUSES
    ]
    if identity_scope is not None:
        reviewed = [
            row
            for row in reviewed
            if str(row.get("identity_scope") or "").strip().casefold() == identity_scope
        ]
    return reviewed


def joined_review_field(rows: list[dict[str, Any]], field: str) -> str:
    return "|".join(dedupe([str(row.get(field) or "").strip() for row in rows]))


def european_source_metadata(row: dict[str, Any]) -> dict[str, str]:
    return {
        "summary_statistics_url": str(row.get("summaryStatistics") or "").strip(),
        "pubmed_id": str(row.get("pubmedId") or "").strip(),
        "paper_title": str(row.get("title") or "").strip(),
        "journal": str(row.get("journal") or "").strip(),
        "publication_date": str(row.get("publicationDate") or "").strip(),
        "efo_traits": str(row.get("efoTraits") or "").strip(),
        "bg_traits": str(row.get("bgTraits") or "").strip(),
        "initial_sample_description": str(row.get("initialSampleDescription") or "").strip(),
        "discovery_sample_ancestry": str(row.get("discoverySampleAncestry") or "").strip(),
    }


def fetch_json(url: str, timeout: int) -> tuple[dict[str, Any] | None, str]:
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "metabo-european-trait-annotation/1.0"})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8")), ""
    except urllib.error.HTTPError as exc:
        return None, f"http_{exc.code}"
    except Exception as exc:
        return None, type(exc).__name__


def pubchem_name_lookup(name: str, cache: dict[str, Any], timeout: int, sleep_seconds: float) -> dict[str, Any]:
    key = f"name::{name.casefold()}"
    if key in cache:
        return cache[key]
    quoted = urllib.parse.quote(name)
    url = f"{PUBCHEM_BASE}/name/{quoted}/cids/JSON"
    payload, error = fetch_json(url, timeout)
    cids = []
    if payload:
        cids = [str(cid) for cid in payload.get("IdentifierList", {}).get("CID", [])[:5]]
    result = {"query": name, "url": url, "cids": cids, "error": error}
    cache[key] = result
    if sleep_seconds:
        time.sleep(sleep_seconds)
    return result


def pubchem_properties(cids: list[str], cache: dict[str, Any], timeout: int, sleep_seconds: float) -> list[dict[str, Any]]:
    missing = [cid for cid in cids if f"cid::{cid}" not in cache]
    for cid in missing:
        url = f"{PUBCHEM_BASE}/cid/{urllib.parse.quote(cid)}/property/Title,MolecularFormula,InChIKey,IUPACName/JSON"
        payload, error = fetch_json(url, timeout)
        props = {}
        if payload:
            rows = payload.get("PropertyTable", {}).get("Properties", [])
            props = rows[0] if rows else {}
        cache[f"cid::{cid}"] = {"cid": cid, "url": url, "properties": props, "error": error}
        if sleep_seconds:
            time.sleep(sleep_seconds)
    return [cache.get(f"cid::{cid}", {}).get("properties", {}) for cid in cids]


def build_annotations(args: argparse.Namespace) -> dict[str, Any]:
    workspace = Path(args.workspace).resolve()
    release_id = args.release_id
    european_path = workspace / args.european_csv
    output_path = workspace / args.output
    review_output_path = workspace / args.review_output
    identity_review_path = workspace / args.identity_review
    cache_path = workspace / args.cache
    compound_dir = workspace / args.compound_root / release_id

    rows = read_european_rows(european_path)
    name_keys, local_cids = load_local_indexes(compound_dir)
    cache = load_cache(cache_path)
    existing_review = load_existing_review(review_output_path)
    identity_review = load_identity_review(identity_review_path)
    output_rows: list[dict[str, Any]] = []

    for index, row in enumerate(rows, start=1):
        accession = normalize_gwas_accession(row.get("accessionId") or row.get("accessionId.1"))
        reported_trait = row.get("reportedTrait", "")
        source_metadata = european_source_metadata(row)
        preserved_review = existing_review.get(accession, {})
        identity_rows = reviewed_identity_rows(identity_review.get(accession, []))
        strict_identity_rows = reviewed_identity_rows(identity_review.get(accession, []), "strict_identity")
        reviewed_names = split_review_values(preserved_review.get("reviewed_name"))
        reviewed_names.extend(row.get("reviewed_name") or row.get("component_name") or "" for row in strict_identity_rows)
        reviewed_cids = split_review_values(preserved_review.get("reviewed_pubchem_cid"))
        for review_row in strict_identity_rows:
            reviewed_cids.extend(split_review_values(review_row.get("reviewed_pubchem_cid")))
        reviewed_hmdb_ids = dedupe(
            [
                *split_review_values(preserved_review.get("reviewed_hmdb_id")),
                *[value for item in strict_identity_rows for value in split_review_values(item.get("reviewed_hmdb_id"))],
            ]
        )
        reviewed_chebi_ids = dedupe(
            [
                *split_review_values(preserved_review.get("reviewed_chebi_id")),
                *[value for item in strict_identity_rows for value in split_review_values(item.get("reviewed_chebi_id"))],
            ]
        )
        reviewed_inchikeys = dedupe(
            [
                *split_review_values(preserved_review.get("reviewed_inchikey")),
                *[value for item in strict_identity_rows for value in split_review_values(item.get("reviewed_inchikey"))],
            ]
        )
        candidates = expanded_trait_candidates(reported_trait)
        candidates = dedupe([*reviewed_names, *candidates])
        local_name_hits = [name for name in candidates if normalize_lookup_key(name) in name_keys]
        pubchem_hits: list[dict[str, Any]] = []
        pubchem_cids: list[str] = dedupe(reviewed_cids)
        pubchem_titles: list[str] = []
        pubchem_formulas: list[str] = []
        pubchem_inchikeys: list[str] = []
        if args.online and not local_name_hits:
            for candidate in candidates[: args.max_queries_per_trait]:
                hit = pubchem_name_lookup(candidate, cache, args.timeout, args.sleep_seconds)
                if hit.get("cids"):
                    pubchem_hits.append(hit)
                    for cid in hit.get("cids", []):
                        if cid not in pubchem_cids:
                            pubchem_cids.append(cid)
                    break
            props = pubchem_properties(pubchem_cids[: args.max_cids_per_trait], cache, args.timeout, args.sleep_seconds) if pubchem_cids else []
            for prop in props:
                title = str(prop.get("Title") or "").strip()
                formula = str(prop.get("MolecularFormula") or "").strip()
                inchikey = str(prop.get("InChIKey") or "").strip()
                if title:
                    pubchem_titles.append(title)
                if formula:
                    pubchem_formulas.append(formula)
                if inchikey:
                    pubchem_inchikeys.append(inchikey)
        mapped_names = dedupe([*local_name_hits, *pubchem_titles, *candidates])
        local_cid_hits = [cid for cid in pubchem_cids if normalize_id_token(cid) in local_cids]
        if local_name_hits:
            status = "local_name"
        elif local_cid_hits:
            status = "pubchem_cid_in_local_index"
        elif pubchem_cids:
            status = "pubchem_only"
        else:
            status = "unresolved"
        output_rows.append(
            {
                "accession_id": accession,
                "reported_trait": reported_trait,
                **source_metadata,
                "resolution_status": status,
                "mapped_names": "|".join(mapped_names),
                "local_name_hits": "|".join(local_name_hits),
                "pubchem_cids": "|".join(pubchem_cids),
                "local_pubchem_cid_hits": "|".join(local_cid_hits),
                "pubchem_titles": "|".join(dedupe(pubchem_titles)),
                "pubchem_formulas": "|".join(dedupe(pubchem_formulas)),
                "pubchem_inchikeys": "|".join(dedupe(pubchem_inchikeys)),
                "candidate_names": "|".join(candidates),
                "pubchem_query": pubchem_hits[0].get("query", "") if pubchem_hits else "",
                "source_name": "GWAS Catalog European metabolite summary + PubChem PUG REST",
                "manual_review_status": joined_review_field(identity_rows, "review_status"),
                "manual_identity_scope": joined_review_field(identity_rows, "identity_scope"),
                "reviewed_name": joined_review_field(strict_identity_rows, "reviewed_name"),
                "reviewed_pubchem_cid": "|".join(dedupe(reviewed_cids)),
                "reviewed_hmdb_id": "|".join(reviewed_hmdb_ids),
                "reviewed_chebi_id": "|".join(reviewed_chebi_ids),
                "reviewed_inchikey": "|".join(reviewed_inchikeys),
                "evidence_source": joined_review_field(identity_rows, "evidence_source"),
                "evidence_url": joined_review_field(identity_rows, "evidence_url"),
            }
        )
        if args.cache_every and index % args.cache_every == 0:
            save_cache(cache_path, cache)
            print(f"[european_trait_annotations] processed {index}/{len(rows)}", flush=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "accession_id",
            "reported_trait",
            *EUROPEAN_SOURCE_METADATA_FIELDS,
            "resolution_status",
            "mapped_names",
            "local_name_hits",
            "pubchem_cids",
            "local_pubchem_cid_hits",
            "pubchem_titles",
            "pubchem_formulas",
            "pubchem_inchikeys",
            "candidate_names",
            "pubchem_query",
            "source_name",
            "manual_review_status",
            "manual_identity_scope",
            "reviewed_name",
            "reviewed_pubchem_cid",
            "reviewed_hmdb_id",
            "reviewed_chebi_id",
            "reviewed_inchikey",
            "evidence_source",
            "evidence_url",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)
    save_cache(cache_path, cache)
    review_bucket_counts = write_unresolved_review(review_output_path, output_rows, existing_review)

    summary = {
        "annotation_version": "european_trait_annotations.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "release_id": release_id,
        "input_rows": len(rows),
        "output": str(output_path),
        "review_output": str(review_output_path),
        "identity_review": str(identity_review_path),
        "cache": str(cache_path),
        "online": bool(args.online),
        "sources": [
            str(european_path),
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON",
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Title,MolecularFormula,InChIKey,IUPACName/JSON",
        ],
        "source_metadata_fields": EUROPEAN_SOURCE_METADATA_FIELDS,
        "manual_identity_reviewed_count": sum(1 for rows in identity_review.values() for row in reviewed_identity_rows(rows)),
        "manual_strict_identity_count": sum(
            1 for rows in identity_review.values() for row in reviewed_identity_rows(rows, "strict_identity")
        ),
        "status_counts": {
            status: sum(1 for row in output_rows if row["resolution_status"] == status)
            for status in sorted({row["resolution_status"] for row in output_rows})
        },
        "unresolved_review_bucket_counts": review_bucket_counts,
    }
    manifest_path = output_path.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build European trait annotation hints for metabolite matching.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", default=DEFAULT_RELEASE_ID)
    parser.add_argument("--compound-root", type=Path, default=Path("compound_match_index"))
    parser.add_argument("--european-csv", default="raw_lake/European/European.csv")
    parser.add_argument("--output", default="raw_lake/European/European_trait_annotations.csv")
    parser.add_argument("--review-output", default="raw_lake/European/European_unresolved_review.csv")
    parser.add_argument("--identity-review", default="manual_sources/european_trait_identity/european_trait_identity_review.csv")
    parser.add_argument("--cache", default="raw_lake/European/pubchem_name_cache.json")
    parser.add_argument("--online", action="store_true")
    parser.add_argument("--timeout", type=int, default=20)
    parser.add_argument("--sleep-seconds", type=float, default=0.12)
    parser.add_argument("--max-queries-per-trait", type=int, default=8)
    parser.add_argument("--max-cids-per-trait", type=int, default=3)
    parser.add_argument("--cache-every", type=int, default=25)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    summary = build_annotations(parse_args(argv))
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
