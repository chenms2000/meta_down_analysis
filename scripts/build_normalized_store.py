#!/usr/bin/env python3
"""Build the MVP normalized store from the versioned raw lake.

The builder keeps external accessions as xrefs and uses deterministic internal
UIDs for canonical rows. It is intentionally batch-oriented: every run writes a
new normalized release folder from immutable raw inputs and records an audit
manifest beside the generated tables.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time
import zipfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Iterator
from xml.etree import ElementTree as ET

try:
    import pyarrow as pa
    import pyarrow.dataset as ds
    import pyarrow.parquet as pq
except ImportError:  # pragma: no cover - exercised only on minimal machines.
    pa = None
    ds = None
    pq = None


csv.field_size_limit(min(sys.maxsize, 2_147_483_647))

DEFAULT_RELEASE_RE = re.compile(r"(?P<release>mvp_[0-9T]+)\.manifest\.json$")

BRIDGEDB_CODE_TO_SOURCE = {
    "Ca": "CAS",
    "Ce": "CHEBI",
    "Ch": "HMDB",
    "Ck": "KEGG.COMPOUND",
    "Cks": "KNAPSACK",
    "Cl": "CHEMBL.COMPOUND",
    "Cpc": "PUBCHEM.COMPOUND",
    "Cs": "CHEMSPIDER",
    "Dr": "DRUGBANK",
    "Ect": "EPA.COMPTOX",
    "Gpl": "GUIDE_TO_PHARMACOLOGY",
    "Ik": "INCHIKEY",
    "Kd": "KEGG.DRUG",
    "Lm": "LIPIDMAPS",
    "Wd": "WIKIDATA",
}
BRIDGEDB_DIRECT_ANCHOR_CODES = ("Ce", "Ch", "Lm", "Cpc", "Ck", "Kd", "Ik")
GTF_ATTR_RE = re.compile(r'(\S+) "([^"]*)"')
SENTENCE_RE = re.compile(r"[^.!?\n]+(?:[.!?]+|$)", re.M)
INTRO_RE = re.compile(r"\n\s*(?:1[.]?\s+)?Introduction\s*\n", re.I)
DEFAULT_SCISPACY_MODEL = "en_core_sci_sm"
ONCOLOGY_DISEASE_RE = re.compile(
    r"\b("
    r"cancer|tumou?r|neoplasm|neoplastic|malignan|carcinoma|sarcoma|"
    r"leukemia|leukaemia|lymphoma|melanoma|glioma|glioblastoma|"
    r"myeloma|adenoma|blastoma|metasta|mesothelioma"
    r")\b",
    re.I,
)
NCIT_NON_DISEASE_CONCEPT_NAME_RE = re.compile(
    r"\b("
    r"agent|therapeutic procedure|therapy agent|chemotherapy|radiation therapy|"
    r"radiotherapy|immunotherapy|screening|measurement|assay|test|criteria|"
    r"risk group|gene mutation|wt allele|allele|oncogene|gene|protein|receptor|"
    r"ligand|factor|kinase|antigen|enzyme|pathway|process|route of administration|"
    r"microenvironment|immunity|angiogenesis|invasion|carcinogenesis|tumou?rigenesis|"
    r"progression|regression|suppression|fibroblast|question|consortium|atlas"
    r")\b",
    re.I,
)
NCIT_NON_DISEASE_CONCEPT_DESCRIPTION_RE = re.compile(
    r"\b(encoded by|encodes|protein|therapeutic agent|therapeutic procedure|"
    r"determination of the amount|route of administration)\b",
    re.I,
)


def utc_now() -> str:
    return dt.datetime.now(dt.UTC).isoformat()


def short_hash(value: str, length: int = 20) -> str:
    return hashlib.sha1(value.encode("utf-8")).hexdigest()[:length]


def stable_uid(prefix: str, *parts: object) -> str:
    payload = "|".join(str(part).strip() for part in parts if part is not None and str(part).strip())
    return f"{prefix}_{short_hash(payload)}"


def parser_hash() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]


def json_dumps(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def dedupe(values: Iterable[str | None]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if value is None:
            continue
        cleaned = str(value).strip()
        if not cleaned or cleaned == "-":
            continue
        if cleaned not in seen:
            seen.add(cleaned)
            out.append(cleaned)
    return out


def maybe_float(value: Any) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text == "-":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def maybe_int(value: Any) -> int | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text == "-":
        return None
    try:
        return int(text)
    except ValueError:
        return None


def text_open(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def iter_tsv(path: Path) -> Iterator[dict[str, str]]:
    with text_open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield {key: (value if value is not None else "") for key, value in row.items()}


def parse_gtf_attrs(value: str) -> dict[str, str]:
    return {match.group(1): match.group(2) for match in GTF_ATTR_RE.finditer(value)}


def normalize_source_name(source: str) -> str:
    text = source.strip().strip('"').strip()
    low = text.lower().replace("_", ".")
    mapping = {
        "chebi": "CHEBI",
        "pubchem.compound": "PUBCHEM.COMPOUND",
        "pubchem": "PUBCHEM.COMPOUND",
        "hmdb": "HMDB",
        "kegg.compound": "KEGG.COMPOUND",
        "kegg.drug": "KEGG.DRUG",
        "lipidmaps": "LIPIDMAPS",
        "inchikey": "INCHIKEY",
        "inchi.key": "INCHIKEY",
        "cas": "CAS",
        "chemspider": "CHEMSPIDER",
        "chembl.compound": "CHEMBL.COMPOUND",
        "drugbank": "DRUGBANK",
        "epa.comptox": "EPA.COMPTOX",
        "comptox": "EPA.COMPTOX",
        "wikidata": "WIKIDATA",
        "guide.to.pharmacology": "GUIDE_TO_PHARMACOLOGY",
        "knapsack": "KNAPSACK",
        "mesh": "MESH",
        "doid": "DOID",
        "mondo": "MONDO",
        "efo": "EFO",
        "ncit": "NCIT",
        "nci": "NCIT",
        "umls": "UMLS",
        "umls.cui": "UMLS",
        "umls_cui": "UMLS",
        "orphanet": "ORPHANET",
        "ensembl": "ENSEMBL",
        "ncbigene": "NCBI.GENE",
        "ncbi.gene": "NCBI.GENE",
        "hgnc": "HGNC",
        "mim": "MIM",
        "uniprot": "UNIPROT",
        "refseq": "REFSEQ",
        "reactome": "REACTOME",
        "wikipathways": "WIKIPATHWAYS",
        "opentargets": "OPENTARGETS",
        "oncotree": "ONCOTREE",
        "cellosaurus": "CELLOSAURUS",
        "cl": "CL",
        "cell.ontology": "CL",
        "cellontology": "CL",
        "uberon": "UBERON",
        "fma": "FMA",
        "bto": "BTO",
        "metabo.state": "METABO_STATE",
        "metabo_state": "METABO_STATE",
    }
    return mapping.get(low, text.upper())


def normalize_xref(raw: str) -> tuple[str, str, str] | None:
    token = raw.strip()
    if not token:
        return None
    token = token.split(" {", 1)[0].split(" !", 1)[0].strip()
    if ":" not in token:
        return None
    source, external_id = token.split(":", 1)
    source = normalize_source_name(source)
    external_id = external_id.strip()
    if not external_id:
        return None
    return source, external_id, f"{source}:{external_id}"


def normalize_prefixed_id(value: str) -> tuple[str, str, str] | None:
    text = value.strip()
    if not text:
        return None
    if ":" in text:
        return normalize_xref(text)
    for prefix in ("MONDO", "DOID", "EFO", "HP", "GO", "CHEBI", "NCIT", "ORPHANET", "CL", "UBERON"):
        marker = f"{prefix}_"
        if text.upper().startswith(marker):
            return normalize_xref(f"{prefix}:{text[len(marker):]}")
    if text.startswith("ENSG"):
        return "ENSEMBL", text.split(".", 1)[0], f"ENSEMBL:{text.split('.', 1)[0]}"
    if text.startswith("WP"):
        return "WIKIPATHWAYS", text, f"WIKIPATHWAYS:{text}"
    if text.startswith("R-HSA-"):
        return "REACTOME", text, f"REACTOME:{text}"
    return None


def xref_key(source: str, external_id: str) -> str:
    return f"{normalize_source_name(source)}:{external_id.strip()}"


def bridgedb_xref(code: str, external_id: str) -> str | None:
    value = str(external_id or "").strip()
    if not value:
        return None
    source = BRIDGEDB_CODE_TO_SOURCE.get(str(code or "").strip())
    if not source:
        return None
    if source == "CHEBI":
        return chebi_key(value.split(":", 1)[1] if value.upper().startswith("CHEBI:") else value)
    if source == "PUBCHEM.COMPOUND" and not re.fullmatch(r"\d+", value):
        return None
    return f"{source}:{value}"


def chebi_key(value: str | int) -> str:
    text = str(value).strip()
    if text.upper().startswith("CHEBI:"):
        return f"CHEBI:{text.split(':', 1)[1]}"
    return f"CHEBI:{text}"


def safe_sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def latest_release_id(manifest_dir: Path) -> str:
    candidates: list[tuple[float, str]] = []
    for path in manifest_dir.glob("*.manifest.json"):
        match = DEFAULT_RELEASE_RE.search(path.name)
        if match:
            candidates.append((path.stat().st_mtime, match.group("release")))
    if not candidates:
        raise FileNotFoundError(f"No mvp_*.manifest.json files found in {manifest_dir}")
    return sorted(candidates)[-1][1]


def iter_obo_terms(path: Path) -> Iterator[dict[str, list[str]]]:
    current: dict[str, list[str]] | None = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "[Term]":
                if current:
                    yield current
                current = {}
                continue
            if line.startswith("["):
                if current:
                    yield current
                current = None
                continue
            if current is None or not line.strip() or ": " not in line:
                continue
            key, value = line.split(": ", 1)
            current.setdefault(key, []).append(value)
    if current:
        yield current


def quoted_value(line: str) -> str:
    match = re.search(r'"([^"]+)"', line)
    return match.group(1) if match else line


def arrow_schema(name: str):
    if pa is None:
        return None
    list_str = pa.list_(pa.string())
    schemas: dict[str, pa.Schema] = {
        "raw_source_registry": pa.schema(
            [
                ("source_id", pa.string()),
                ("source_name", pa.string()),
                ("layer", pa.string()),
                ("official_url", pa.string()),
                ("license_policy", pa.string()),
                ("release_id", pa.string()),
                ("source_release", pa.string()),
                ("local_root", pa.string()),
                ("raw_file_count", pa.int64()),
                ("raw_bytes", pa.int64()),
                ("checksum_digest", pa.string()),
                ("created_at_utc", pa.string()),
            ]
        ),
        "raw_file_manifest": pa.schema(
            [
                ("file_uid", pa.string()),
                ("release_id", pa.string()),
                ("source_id", pa.string()),
                ("source_name", pa.string()),
                ("path", pa.string()),
                ("url", pa.string()),
                ("bytes", pa.int64()),
                ("sha256", pa.string()),
                ("status", pa.string()),
                ("local_exists", pa.bool_()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
            ]
        ),
        "metabolites": pa.schema(
            [
                ("metabolite_uid", pa.string()),
                ("canonical_name", pa.string()),
                ("synonyms", list_str),
                ("formula", pa.string()),
                ("exact_mass", pa.float64()),
                ("charge", pa.int64()),
                ("inchikey", pa.string()),
                ("smiles", pa.string()),
                ("external_xrefs", list_str),
                ("source_priority", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("checksum", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "metabolite_xrefs": pa.schema(
            [
                ("xref_uid", pa.string()),
                ("metabolite_uid", pa.string()),
                ("xref_source", pa.string()),
                ("xref_id", pa.string()),
                ("xref_key", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "hmdb_metabolite_details": pa.schema(
            [
                ("hmdb_detail_uid", pa.string()),
                ("metabolite_uid", pa.string()),
                ("hmdb_accession", pa.string()),
                ("name", pa.string()),
                ("status", pa.string()),
                ("description", pa.string()),
                ("chemical_formula", pa.string()),
                ("average_molecular_weight", pa.float64()),
                ("monoisotopic_molecular_weight", pa.float64()),
                ("inchikey", pa.string()),
                ("smiles", pa.string()),
                ("external_xrefs", list_str),
                ("synonyms", list_str),
                ("biospecimen_locations", list_str),
                ("tissue_locations", list_str),
                ("cellular_locations", list_str),
                ("pathway_refs", list_str),
                ("protein_refs", list_str),
                ("disease_refs", list_str),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("checksum", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "genes": pa.schema(
            [
                ("gene_uid", pa.string()),
                ("ensembl_gene_id", pa.string()),
                ("entrez_gene_id", pa.string()),
                ("symbol", pa.string()),
                ("aliases", list_str),
                ("description", pa.string()),
                ("taxon", pa.string()),
                ("biotype", pa.string()),
                ("chromosome", pa.string()),
                ("start", pa.int64()),
                ("end", pa.int64()),
                ("strand", pa.string()),
                ("uniprot_ids", list_str),
                ("external_xrefs", list_str),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("checksum", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "gene_xrefs": pa.schema(
            [
                ("xref_uid", pa.string()),
                ("gene_uid", pa.string()),
                ("xref_source", pa.string()),
                ("xref_id", pa.string()),
                ("xref_key", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "pathways": pa.schema(
            [
                ("pathway_uid", pa.string()),
                ("name", pa.string()),
                ("species", pa.string()),
                ("source_name", pa.string()),
                ("primary_external_id", pa.string()),
                ("hierarchy_path", pa.string()),
                ("external_xrefs", list_str),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("checksum", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "reactions": pa.schema(
            [
                ("reaction_uid", pa.string()),
                ("name", pa.string()),
                ("species", pa.string()),
                ("primary_external_id", pa.string()),
                ("pathway_uid", pa.string()),
                ("pathway_external_id", pa.string()),
                ("source_name", pa.string()),
                ("external_xrefs", list_str),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("checksum", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "pathway_hierarchy_edges": pa.schema(
            [
                ("edge_uid", pa.string()),
                ("parent_pathway_uid", pa.string()),
                ("child_pathway_uid", pa.string()),
                ("parent_external_id", pa.string()),
                ("child_external_id", pa.string()),
                ("predicate", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "reaction_participants": pa.schema(
            [
                ("edge_uid", pa.string()),
                ("reaction_uid", pa.string()),
                ("participant_uid", pa.string()),
                ("participant_external_id", pa.string()),
                ("participant_type", pa.string()),
                ("role", pa.string()),
                ("source_name", pa.string()),
                ("source_record_id", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "gene_pathway_edges": pa.schema(
            [
                ("edge_uid", pa.string()),
                ("gene_uid", pa.string()),
                ("pathway_uid", pa.string()),
                ("gene_external_id", pa.string()),
                ("pathway_external_id", pa.string()),
                ("predicate", pa.string()),
                ("source_name", pa.string()),
                ("source_record_id", pa.string()),
                ("evidence_level", pa.string()),
                ("evidence_code", pa.string()),
                ("species", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "metabolite_pathway_edges": pa.schema(
            [
                ("edge_uid", pa.string()),
                ("metabolite_uid", pa.string()),
                ("pathway_uid", pa.string()),
                ("metabolite_external_id", pa.string()),
                ("pathway_external_id", pa.string()),
                ("predicate", pa.string()),
                ("source_name", pa.string()),
                ("source_record_id", pa.string()),
                ("evidence_level", pa.string()),
                ("evidence_code", pa.string()),
                ("species", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "diseases": pa.schema(
            [
                ("disease_uid", pa.string()),
                ("primary_external_id", pa.string()),
                ("name", pa.string()),
                ("aliases", list_str),
                ("description", pa.string()),
                ("parents", list_str),
                ("external_xrefs", list_str),
                ("source_priority", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("checksum", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "disease_xrefs": pa.schema(
            [
                ("xref_uid", pa.string()),
                ("disease_uid", pa.string()),
                ("xref_source", pa.string()),
                ("xref_id", pa.string()),
                ("xref_key", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "disease_parent_edges": pa.schema(
            [
                ("edge_uid", pa.string()),
                ("child_disease_uid", pa.string()),
                ("parent_disease_uid", pa.string()),
                ("child_external_id", pa.string()),
                ("parent_external_id", pa.string()),
                ("predicate", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "cell_types": pa.schema(
            [
                ("cell_type_uid", pa.string()),
                ("primary_external_id", pa.string()),
                ("name", pa.string()),
                ("aliases", list_str),
                ("description", pa.string()),
                ("parents", list_str),
                ("external_xrefs", list_str),
                ("source_priority", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("checksum", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "cell_type_xrefs": pa.schema(
            [
                ("xref_uid", pa.string()),
                ("cell_type_uid", pa.string()),
                ("xref_source", pa.string()),
                ("xref_id", pa.string()),
                ("xref_key", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "cell_type_parent_edges": pa.schema(
            [
                ("edge_uid", pa.string()),
                ("child_cell_type_uid", pa.string()),
                ("parent_cell_type_uid", pa.string()),
                ("child_external_id", pa.string()),
                ("parent_external_id", pa.string()),
                ("predicate", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "cell_states": pa.schema(
            [
                ("cell_state_uid", pa.string()),
                ("primary_external_id", pa.string()),
                ("name", pa.string()),
                ("aliases", list_str),
                ("description", pa.string()),
                ("parents", list_str),
                ("external_xrefs", list_str),
                ("state_category", pa.string()),
                ("source_priority", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("checksum", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "cell_state_xrefs": pa.schema(
            [
                ("xref_uid", pa.string()),
                ("cell_state_uid", pa.string()),
                ("xref_source", pa.string()),
                ("xref_id", pa.string()),
                ("xref_key", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "cell_state_parent_edges": pa.schema(
            [
                ("edge_uid", pa.string()),
                ("child_cell_state_uid", pa.string()),
                ("parent_cell_state_uid", pa.string()),
                ("child_external_id", pa.string()),
                ("parent_external_id", pa.string()),
                ("predicate", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "tissues": pa.schema(
            [
                ("tissue_uid", pa.string()),
                ("primary_external_id", pa.string()),
                ("name", pa.string()),
                ("aliases", list_str),
                ("description", pa.string()),
                ("parents", list_str),
                ("external_xrefs", list_str),
                ("source_priority", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("checksum", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "tissue_xrefs": pa.schema(
            [
                ("xref_uid", pa.string()),
                ("tissue_uid", pa.string()),
                ("xref_source", pa.string()),
                ("xref_id", pa.string()),
                ("xref_key", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "tissue_parent_edges": pa.schema(
            [
                ("edge_uid", pa.string()),
                ("child_tissue_uid", pa.string()),
                ("parent_tissue_uid", pa.string()),
                ("child_external_id", pa.string()),
                ("parent_external_id", pa.string()),
                ("predicate", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "targets": pa.schema(
            [
                ("target_uid", pa.string()),
                ("target_external_id", pa.string()),
                ("preferred_name", pa.string()),
                ("approved_symbol", pa.string()),
                ("target_type", pa.string()),
                ("gene_uid", pa.string()),
                ("tractability_flags_json", pa.string()),
                ("external_xrefs", list_str),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("checksum", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "target_xrefs": pa.schema(
            [
                ("xref_uid", pa.string()),
                ("target_uid", pa.string()),
                ("xref_source", pa.string()),
                ("xref_id", pa.string()),
                ("xref_key", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "target_gene_edges": pa.schema(
            [
                ("edge_uid", pa.string()),
                ("target_uid", pa.string()),
                ("gene_uid", pa.string()),
                ("target_external_id", pa.string()),
                ("gene_external_id", pa.string()),
                ("predicate", pa.string()),
                ("source_name", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "target_disease_edges": pa.schema(
            [
                ("edge_uid", pa.string()),
                ("target_uid", pa.string()),
                ("disease_uid", pa.string()),
                ("target_external_id", pa.string()),
                ("disease_external_id", pa.string()),
                ("predicate", pa.string()),
                ("association_score", pa.float64()),
                ("evidence_count", pa.int64()),
                ("aggregation_type", pa.string()),
                ("aggregation_value", pa.string()),
                ("current_novelty", pa.float64()),
                ("score_components_json", pa.string()),
                ("source_name", pa.string()),
                ("source_record_id", pa.string()),
                ("evidence_level", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "articles": pa.schema(
            [
                ("article_uid", pa.string()),
                ("pmid", pa.string()),
                ("pmcid", pa.string()),
                ("doi", pa.string()),
                ("title", pa.string()),
                ("abstract", pa.string()),
                ("journal", pa.string()),
                ("pub_date", pa.string()),
                ("article_type", pa.string()),
                ("source_path", pa.string()),
                ("source_bytes", pa.int64()),
                ("checksum", pa.string()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "sentences": pa.schema(
            [
                ("sentence_uid", pa.string()),
                ("article_uid", pa.string()),
                ("pmid", pa.string()),
                ("pmcid", pa.string()),
                ("section", pa.string()),
                ("sentence_text", pa.string()),
                ("text_hash", pa.string()),
                ("start_offset", pa.int64()),
                ("end_offset", pa.int64()),
                ("language", pa.string()),
                ("segmenter", pa.string()),
                ("parser_model", pa.string()),
                ("source_release", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "sentence_entity_candidates": pa.schema(
            [
                ("candidate_uid", pa.string()),
                ("sentence_uid", pa.string()),
                ("article_uid", pa.string()),
                ("pmid", pa.string()),
                ("pmcid", pa.string()),
                ("section", pa.string()),
                ("mention_text", pa.string()),
                ("normalized_surface", pa.string()),
                ("start_offset", pa.int64()),
                ("end_offset", pa.int64()),
                ("scispacy_label", pa.string()),
                ("candidate_source", pa.string()),
                ("parser_model", pa.string()),
                ("confidence_proxy", pa.float64()),
                ("source_release", pa.string()),
                ("license_id", pa.string()),
                ("parser_hash", pa.string()),
            ]
        ),
        "build_notes": pa.schema(
            [
                ("note_uid", pa.string()),
                ("severity", pa.string()),
                ("component", pa.string()),
                ("message", pa.string()),
                ("created_at_utc", pa.string()),
            ]
        ),
    }
    return schemas[name]


class TableWriter:
    def __init__(self, output_dir: Path, name: str):
        self.name = name
        self.schema = arrow_schema(name)
        self.count = 0
        self._writer = None
        self._json_handle = None
        if pa is not None:
            self.path = output_dir / f"{name}.parquet"
            if self.path.exists():
                self.path.unlink()
            self._writer = pq.ParquetWriter(self.path, self.schema, compression="snappy")
        else:
            self.path = output_dir / f"{name}.jsonl.gz"
            if self.path.exists():
                self.path.unlink()
            self._json_handle = gzip.open(self.path, "wt", encoding="utf-8")

    def write(self, rows: Iterable[dict[str, Any]]) -> None:
        batch = list(rows)
        if not batch:
            return
        if self._writer is not None:
            table = pa.Table.from_pylist(batch, schema=self.schema)
            self._writer.write_table(table)
        else:
            assert self._json_handle is not None
            for row in batch:
                self._json_handle.write(json_dumps(row) + "\n")
        self.count += len(batch)

    def close(self) -> dict[str, Any]:
        if self._writer is not None:
            self._writer.close()
        if self._json_handle is not None:
            self._json_handle.close()
        return {"table": self.name, "path": str(self.path), "rows": self.count}


@dataclass
class BuildContext:
    workspace: Path
    raw_root: Path
    manifest_dir: Path
    output_root: Path
    release_id: str
    manifest: dict[str, Any]
    catalog: dict[str, Any]
    parser_hash: str
    created_at_utc: str = field(default_factory=utc_now)
    table_audit: list[dict[str, Any]] = field(default_factory=list)
    notes: list[dict[str, Any]] = field(default_factory=list)
    source_by_id: dict[str, dict[str, Any]] = field(default_factory=dict)
    source_file_rows: dict[str, list[dict[str, Any]]] = field(default_factory=lambda: defaultdict(list))
    source_checksum: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        catalog_sources = {source["id"]: source for source in self.catalog.get("sources", [])}
        manifest_sources = {source["id"]: source for source in self.manifest.get("sources", [])}
        self.source_by_id = {**catalog_sources, **manifest_sources}
        for row in self.manifest.get("files", []):
            self.source_file_rows[row.get("source_id", "")].append(row)
        for source_id, rows in self.source_file_rows.items():
            checksums = sorted(row.get("sha256", "") for row in rows if row.get("sha256"))
            self.source_checksum[source_id] = hashlib.sha256("".join(checksums).encode("utf-8")).hexdigest()

    @property
    def output_dir(self) -> Path:
        return self.output_root / self.release_id

    def source_dir(self, source_id: str) -> Path:
        return self.raw_root / source_id / self.release_id

    def source_release(self, source_id: str) -> str:
        if source_id == "hmdb":
            return "HMDB 5.0"
        if source_id == "rhea":
            return "rhea:current"
        if source_id == "oncotree":
            return "oncotree:latest_stable"
        if source_id == "ncit":
            obo_path = self.source_dir("ncit") / "ncit" / "ncit.obo"
            if obo_path.exists():
                try:
                    with obo_path.open("r", encoding="utf-8", errors="replace") as handle:
                        for _ in range(50):
                            line = handle.readline()
                            if not line:
                                break
                            if line.startswith("data-version:"):
                                return "ncit:" + line.split(":", 1)[1].strip()
                except OSError:
                    pass
            return "ncit:current"
        if source_id in {"cell_ontology", "uberon", "efo"}:
            for obo_path in sorted(self.source_dir(source_id).glob("*.obo")):
                try:
                    with obo_path.open("r", encoding="utf-8", errors="replace") as handle:
                        for _ in range(80):
                            line = handle.readline()
                            if not line:
                                break
                            if line.startswith("data-version:"):
                                return f"{source_id}:{line.split(':', 1)[1].strip()}"
                except OSError:
                    continue
            return f"{source_id}:current"
        if source_id == "opentargets_core":
            for row in self.source_file_rows.get(source_id, []):
                url = row.get("url", "")
                match = re.search(r"/platform/([^/]+)/", url)
                if match:
                    return f"opentargets:{match.group(1)}"
        if source_id == "ensembl_human":
            for row in self.source_file_rows.get(source_id, []):
                match = re.search(r"GRCh38\.([0-9]+)\.", row.get("path", ""))
                if match:
                    return f"ensembl:{match.group(1)}"
        return self.release_id

    def license_id(self, source_id: str) -> str:
        source = self.source_by_id.get(source_id, {})
        return f"{source.get('layer', 'open_core')}:{source_id}"

    def checksum(self, source_id: str) -> str:
        return self.source_checksum.get(source_id, "")

    def note(self, severity: str, component: str, message: str) -> None:
        self.notes.append(
            {
                "note_uid": stable_uid("note", severity, component, message),
                "severity": severity,
                "component": component,
                "message": message,
                "created_at_utc": utc_now(),
            }
        )

    def writer(self, name: str) -> TableWriter:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        return TableWriter(self.output_dir, name)

    def close_writer(self, writer: TableWriter, started: float) -> None:
        audit = writer.close()
        audit["seconds"] = round(time.time() - started, 3)
        audit["parser_hash"] = self.parser_hash
        self.table_audit.append(audit)

    def discover_unmanifested_local_sources(self) -> None:
        for source_id in sorted(self.source_by_id):
            if source_id in self.source_file_rows:
                continue
            root = self.source_dir(source_id)
            if not root.exists():
                continue
            rows = []
            for path in sorted(p for p in root.rglob("*") if p.is_file()):
                rel = path.relative_to(self.workspace)
                rows.append(
                    {
                        "bytes": path.stat().st_size,
                        "path": str(rel),
                        "sha256": safe_sha256_file(path),
                        "source_id": source_id,
                        "status": "manual_existing",
                        "url": self.source_by_id[source_id].get("official_url", ""),
                    }
                )
            if rows:
                self.source_file_rows[source_id] = rows
                checksums = sorted(row["sha256"] for row in rows)
                self.source_checksum[source_id] = hashlib.sha256("".join(checksums).encode("utf-8")).hexdigest()


def load_manifest(workspace: Path, manifest_dir: Path, release_id: str | None) -> tuple[str, dict[str, Any]]:
    resolved_release = release_id or latest_release_id(manifest_dir)
    path = manifest_dir / f"{resolved_release}.manifest.json"
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")
    return resolved_release, json.loads(path.read_text(encoding="utf-8"))


def load_catalog(workspace: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    catalog_ref = str(manifest.get("catalog", "config/source_catalog.toml")).replace("\\", "/")
    catalog_path = workspace / catalog_ref
    with catalog_path.open("rb") as handle:
        import tomllib

        return tomllib.load(handle)


def build_raw_registry(ctx: BuildContext) -> None:
    started = time.time()
    source_writer = ctx.writer("raw_source_registry")
    file_writer = ctx.writer("raw_file_manifest")
    source_rows = []
    file_rows = []
    source_ids = sorted({source["id"] for source in ctx.manifest.get("sources", [])} | set(ctx.source_file_rows))
    for source_id in source_ids:
        source = ctx.source_by_id.get(source_id, {"id": source_id})
        files = ctx.source_file_rows.get(source_id, [])
        source_rows.append(
            {
                "source_id": source_id,
                "source_name": source.get("name", ""),
                "layer": source.get("layer", ""),
                "official_url": source.get("official_url", ""),
                "license_policy": source.get("license_policy", ""),
                "release_id": ctx.release_id,
                "source_release": ctx.source_release(source_id),
                "local_root": str(ctx.source_dir(source_id)),
                "raw_file_count": len(files),
                "raw_bytes": sum(int(row.get("bytes") or 0) for row in files),
                "checksum_digest": ctx.checksum(source_id),
                "created_at_utc": ctx.created_at_utc,
            }
        )
        for row in files:
            local_path = ctx.workspace / row.get("path", "")
            file_rows.append(
                {
                    "file_uid": stable_uid("rawfile", ctx.release_id, row.get("path", "")),
                    "release_id": ctx.release_id,
                    "source_id": source_id,
                    "source_name": source.get("name", ""),
                    "path": row.get("path", ""),
                    "url": row.get("url", ""),
                    "bytes": int(row.get("bytes") or 0),
                    "sha256": row.get("sha256", ""),
                    "status": row.get("status", ""),
                    "local_exists": local_path.exists(),
                    "source_release": ctx.source_release(source_id),
                    "license_id": ctx.license_id(source_id),
                }
            )
    source_writer.write(source_rows)
    file_writer.write(file_rows)
    ctx.close_writer(source_writer, started)
    ctx.close_writer(file_writer, started)


def parse_chebi_obo_xrefs(path: Path) -> dict[str, set[str]]:
    xrefs_by_chebi: dict[str, set[str]] = defaultdict(set)
    for term in iter_obo_terms(path):
        ids = term.get("id", [])
        if not ids or not ids[0].startswith("CHEBI:"):
            continue
        chebi = ids[0]
        xrefs_by_chebi[chebi].add(chebi)
        for alt_id in term.get("alt_id", []):
            if alt_id.startswith("CHEBI:"):
                xrefs_by_chebi[chebi].add(alt_id)
        for raw_xref in term.get("xref", []):
            parsed = normalize_xref(raw_xref)
            if parsed:
                xrefs_by_chebi[chebi].add(parsed[2])
    return xrefs_by_chebi


def write_xref_rows(
    writer: TableWriter,
    entity_field: str,
    entity_uid: str,
    xrefs: Iterable[str],
    source_name: str,
    source_release: str,
    license_id: str,
    parser_hash_value: str,
) -> None:
    rows = []
    for xref in sorted(set(xrefs)):
        parsed = normalize_xref(xref)
        if not parsed:
            continue
        source, external_id, key = parsed
        rows.append(
            {
                "xref_uid": stable_uid("xref", entity_uid, key),
                entity_field: entity_uid,
                "xref_source": source,
                "xref_id": external_id,
                "xref_key": key,
                "source_name": source_name,
                "source_release": source_release,
                "license_id": license_id,
                "parser_hash": parser_hash_value,
            }
        )
    writer.write(rows)


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def child_by_name(node: ET.Element, name: str) -> ET.Element | None:
    for child in list(node):
        if local_name(child.tag) == name:
            return child
    return None


def children_by_name(node: ET.Element, name: str) -> list[ET.Element]:
    return [child for child in list(node) if local_name(child.tag) == name]


def children_path(node: ET.Element, *names: str) -> list[ET.Element]:
    nodes = [node]
    for name in names:
        next_nodes: list[ET.Element] = []
        for current in nodes:
            next_nodes.extend(children_by_name(current, name))
        nodes = next_nodes
        if not nodes:
            break
    return nodes


def child_text_local(node: ET.Element, name: str) -> str:
    child = child_by_name(node, name)
    return (child.text or "").strip() if child is not None and child.text else ""


def child_text_path(node: ET.Element, *names: str) -> str:
    current = node
    for name in names:
        child = child_by_name(current, name)
        if child is None:
            return ""
        current = child
    return (current.text or "").strip() if current.text else ""


def child_texts_path(node: ET.Element, *names: str) -> list[str]:
    return dedupe((item.text or "").strip() for item in children_path(node, *names) if item.text)


def iter_hmdb_metabolite_elements(path: Path) -> Iterator[ET.Element]:
    def emit_from(handle) -> Iterator[ET.Element]:
        for _, elem in ET.iterparse(handle, events=("end",)):
            if local_name(elem.tag) == "metabolite":
                yield elem
                elem.clear()

    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            members = [name for name in archive.namelist() if name.lower().endswith(".xml")]
            preferred = [name for name in members if "metabolite" in Path(name).name.lower()]
            selected = (preferred or members)
            if not selected:
                return
            with archive.open(selected[0]) as handle:
                yield from emit_from(handle)
    elif path.suffix.lower() == ".gz":
        with gzip.open(path, "rb") as handle:
            yield from emit_from(handle)
    else:
        with path.open("rb") as handle:
            yield from emit_from(handle)


def find_hmdb_metabolite_files(root: Path) -> list[Path]:
    if not root.exists():
        return []
    allowed_suffixes = {".zip", ".xml", ".gz"}
    candidates = [
        path
        for path in root.rglob("*")
        if path.is_file()
        and path.suffix.lower() in allowed_suffixes
        and "metabolite" in path.name.lower()
        and "spectra" not in path.name.lower()
    ]
    full_markers = ("hmdb_metabolites", "all_metabolites", "metabolites.xml")
    full = [path for path in candidates if any(marker in path.name.lower() for marker in full_markers)]
    if full:
        return sorted(full)
    return sorted(candidates)


def hmdb_metabolite_payload(elem: ET.Element) -> dict[str, Any]:
    accession = child_text_local(elem, "accession")
    secondary = child_texts_path(elem, "secondary_accessions", "accession")
    chebi_id = child_text_local(elem, "chebi_id")
    pubchem_id = child_text_local(elem, "pubchem_compound_id")
    kegg_id = child_text_local(elem, "kegg_id")
    cas_id = child_text_local(elem, "cas_registry_number")
    xrefs = [f"HMDB:{accession}"]
    xrefs.extend(f"HMDB:{value}" for value in secondary)
    if chebi_id:
        xrefs.append(f"CHEBI:{chebi_id.split(':')[-1]}")
    if pubchem_id:
        xrefs.append(f"PUBCHEM.COMPOUND:{pubchem_id}")
    if kegg_id:
        xrefs.append(f"KEGG.COMPOUND:{kegg_id}")
    if cas_id:
        xrefs.append(f"CAS:{cas_id}")

    pathway_refs = []
    for pathway in children_path(elem, "pathways", "pathway"):
        parts = dedupe(
            [
                child_text_local(pathway, "name"),
                f"SMPDB:{child_text_local(pathway, 'smpdb_id')}" if child_text_local(pathway, "smpdb_id") else "",
                f"KEGG.PATHWAY:{child_text_local(pathway, 'kegg_map_id')}" if child_text_local(pathway, "kegg_map_id") else "",
            ]
        )
        if parts:
            pathway_refs.append("|".join(parts))

    protein_refs = []
    for protein in children_path(elem, "protein_associations", "protein"):
        parts = dedupe(
            [
                child_text_local(protein, "protein_accession"),
                child_text_local(protein, "name"),
                child_text_local(protein, "gene_name"),
                f"UNIPROT:{child_text_local(protein, 'uniprot_id')}" if child_text_local(protein, "uniprot_id") else "",
            ]
        )
        if parts:
            protein_refs.append("|".join(parts))

    disease_refs = []
    for disease in children_path(elem, "diseases", "disease"):
        parts = dedupe([child_text_local(disease, "name"), child_text_local(disease, "omim_id")])
        if parts:
            disease_refs.append("|".join(parts))

    return {
        "accession": accession,
        "name": child_text_local(elem, "name"),
        "status": child_text_local(elem, "status"),
        "description": child_text_local(elem, "description"),
        "chemical_formula": child_text_local(elem, "chemical_formula"),
        "average_molecular_weight": maybe_float(child_text_local(elem, "average_molecular_weight")),
        "monoisotopic_molecular_weight": maybe_float(
            child_text_local(elem, "monisotopic_molecular_weight")
            or child_text_local(elem, "monoisotopic_molecular_weight")
        ),
        "inchikey": child_text_local(elem, "inchikey"),
        "smiles": child_text_local(elem, "smiles"),
        "external_xrefs": dedupe(xrefs),
        "synonyms": child_texts_path(elem, "synonyms", "synonym"),
        "biospecimen_locations": child_texts_path(elem, "biospecimen_locations", "biospecimen"),
        "tissue_locations": child_texts_path(elem, "tissue_locations", "tissue"),
        "cellular_locations": child_texts_path(elem, "cellular_locations", "cellular"),
        "pathway_refs": dedupe(pathway_refs),
        "protein_refs": dedupe(protein_refs),
        "disease_refs": dedupe(disease_refs),
    }


def ingest_hmdb_metabolites(
    ctx: BuildContext,
    metabolite_writer: TableWriter,
    xref_writer: TableWriter,
    hmdb_writer: TableWriter,
    uid_by_xref: dict[str, str],
    seen_metabolites: set[str],
    seen_xrefs: set[tuple[str, str]],
) -> None:
    hmdb_dir = ctx.source_dir("hmdb")
    files = find_hmdb_metabolite_files(hmdb_dir)
    if not files:
        ctx.note(
            "info",
            "hmdb",
            "HMDB parser is enabled, but no HMDB metabolite XML/ZIP file was found under raw_lake/hmdb/<release_id>/.",
        )
        return

    source_release = ctx.source_release("hmdb")
    license_id = ctx.license_id("hmdb")
    checksum = ctx.checksum("hmdb")
    new_metabolite_rows = []
    xref_rows = []
    detail_rows = []
    parsed_count = 0

    for file_path in files:
        for elem in iter_hmdb_metabolite_elements(file_path):
            payload = hmdb_metabolite_payload(elem)
            accession = payload["accession"]
            if not accession:
                continue
            parsed_count += 1
            metabolite_uid = ""
            for xref in payload["external_xrefs"]:
                metabolite_uid = uid_by_xref.get(xref, "")
                if metabolite_uid:
                    break
            if not metabolite_uid:
                metabolite_uid = stable_uid("metabolite", "HMDB", accession)
            for xref in payload["external_xrefs"]:
                uid_by_xref[xref] = metabolite_uid

            if metabolite_uid not in seen_metabolites:
                new_metabolite_rows.append(
                    {
                        "metabolite_uid": metabolite_uid,
                        "canonical_name": payload["name"] or accession,
                        "synonyms": payload["synonyms"],
                        "formula": payload["chemical_formula"],
                        "exact_mass": payload["monoisotopic_molecular_weight"],
                        "charge": None,
                        "inchikey": payload["inchikey"],
                        "smiles": payload["smiles"],
                        "external_xrefs": payload["external_xrefs"],
                        "source_priority": "HMDB",
                        "source_release": source_release,
                        "license_id": license_id,
                        "checksum": checksum,
                        "parser_hash": ctx.parser_hash,
                    }
                )
                seen_metabolites.add(metabolite_uid)

            for xref in payload["external_xrefs"]:
                parsed = normalize_xref(xref)
                if not parsed:
                    continue
                key = (metabolite_uid, parsed[2])
                if key in seen_xrefs:
                    continue
                seen_xrefs.add(key)
                xref_rows.append(
                    {
                        "xref_uid": stable_uid("xref", metabolite_uid, parsed[2]),
                        "metabolite_uid": metabolite_uid,
                        "xref_source": parsed[0],
                        "xref_id": parsed[1],
                        "xref_key": parsed[2],
                        "source_name": "HMDB",
                        "source_release": source_release,
                        "license_id": license_id,
                        "parser_hash": ctx.parser_hash,
                    }
                )

            detail_rows.append(
                {
                    "hmdb_detail_uid": stable_uid("hmdbdetail", accession),
                    "metabolite_uid": metabolite_uid,
                    "hmdb_accession": accession,
                    "name": payload["name"],
                    "status": payload["status"],
                    "description": payload["description"],
                    "chemical_formula": payload["chemical_formula"],
                    "average_molecular_weight": payload["average_molecular_weight"],
                    "monoisotopic_molecular_weight": payload["monoisotopic_molecular_weight"],
                    "inchikey": payload["inchikey"],
                    "smiles": payload["smiles"],
                    "external_xrefs": payload["external_xrefs"],
                    "synonyms": payload["synonyms"],
                    "biospecimen_locations": payload["biospecimen_locations"],
                    "tissue_locations": payload["tissue_locations"],
                    "cellular_locations": payload["cellular_locations"],
                    "pathway_refs": payload["pathway_refs"],
                    "protein_refs": payload["protein_refs"],
                    "disease_refs": payload["disease_refs"],
                    "source_release": source_release,
                    "license_id": license_id,
                    "checksum": checksum,
                    "parser_hash": ctx.parser_hash,
                }
            )

            if len(new_metabolite_rows) >= 10_000:
                metabolite_writer.write(new_metabolite_rows)
                new_metabolite_rows = []
            if len(xref_rows) >= 25_000:
                xref_writer.write(xref_rows)
                xref_rows = []
            if len(detail_rows) >= 10_000:
                hmdb_writer.write(detail_rows)
                detail_rows = []

    metabolite_writer.write(new_metabolite_rows)
    xref_writer.write(xref_rows)
    hmdb_writer.write(detail_rows)
    ctx.note(
        "info",
        "hmdb",
        f"Parsed {parsed_count} HMDB metabolite records from {len(files)} local file(s) as a non-commercial research enhancement layer.",
    )


def find_first_existing(paths: Iterable[Path | None]) -> Path | None:
    for path in paths:
        if path and path.exists():
            return path
    return None


def bridgedb_java_paths() -> tuple[Path | None, Path | None, Path | None]:
    java_home = os.environ.get("JAVA_HOME", "")
    derby_home = os.environ.get("DERBY_HOME", "")
    java_candidates = [
        Path(os.environ["BRIDGEDB_JAVA"]) if os.environ.get("BRIDGEDB_JAVA") else None,
        Path(java_home) / "bin" / ("java.exe" if os.name == "nt" else "java") if java_home else None,
        Path("D:/java/bin/java.exe"),
        Path(shutil.which("java")) if shutil.which("java") else None,
    ]
    javac_candidates = [
        Path(os.environ["BRIDGEDB_JAVAC"]) if os.environ.get("BRIDGEDB_JAVAC") else None,
        Path(java_home) / "bin" / ("javac.exe" if os.name == "nt" else "javac") if java_home else None,
        Path("D:/java/bin/javac.exe"),
        Path(shutil.which("javac")) if shutil.which("javac") else None,
    ]
    derby_candidates = [
        Path(os.environ["BRIDGEDB_DERBY_JAR"]) if os.environ.get("BRIDGEDB_DERBY_JAR") else None,
        Path(derby_home) / "lib" / "derby.jar" if derby_home else None,
        Path("D:/java/Derby/db-derby-10.14.2.0-bin/lib/derby.jar"),
    ]
    java_root = Path("D:/java")
    if java_root.exists():
        derby_candidates.extend(sorted(java_root.glob("Derby/**/lib/derby.jar")))
    return find_first_existing(java_candidates), find_first_existing(javac_candidates), find_first_existing(derby_candidates)


def compile_bridgedb_exporter(ctx: BuildContext, javac_path: Path, derby_jar: Path) -> Path:
    source = ctx.workspace / "scripts" / "BridgeDbDerbyExport.java"
    classes = ctx.workspace / "tmp" / "bridgedb_java_classes"
    classes.mkdir(parents=True, exist_ok=True)
    class_file = classes / "BridgeDbDerbyExport.class"
    if not class_file.exists() or class_file.stat().st_mtime < source.stat().st_mtime:
        result = subprocess.run(
            [
                str(javac_path),
                "-cp",
                str(derby_jar),
                "-d",
                str(classes),
                str(source),
            ],
            cwd=ctx.workspace,
            text=True,
            capture_output=True,
            timeout=120,
        )
        if result.returncode != 0:
            raise RuntimeError(f"javac failed: {result.stderr.strip() or result.stdout.strip()}")
    return classes


def extract_bridgedb_database(ctx: BuildContext, bridge_file: Path) -> Path:
    target = ctx.workspace / "tmp" / "bridgedb_derby" / bridge_file.stem
    database = target / "database"
    if (database / "service.properties").exists():
        return database
    target.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(bridge_file) as archive:
        archive.extractall(target)
    if not (database / "service.properties").exists():
        raise FileNotFoundError(f"BridgeDb Derby service.properties not found after extracting {bridge_file}")
    return database


def run_bridgedb_exporter(
    ctx: BuildContext,
    bridge_file: Path,
    mode: str,
    code_csv: str = "",
) -> Iterator[list[str]]:
    java_path, javac_path, derby_jar = bridgedb_java_paths()
    if java_path is None or javac_path is None or derby_jar is None:
        raise FileNotFoundError("Java, javac, or derby.jar was not found. Set BRIDGEDB_JAVA, BRIDGEDB_JAVAC, or BRIDGEDB_DERBY_JAR if needed.")
    classes = compile_bridgedb_exporter(ctx, javac_path, derby_jar)
    database = extract_bridgedb_database(ctx, bridge_file)
    classpath = os.pathsep.join([str(classes), str(derby_jar)])
    command = [str(java_path), "-cp", classpath, "BridgeDbDerbyExport", str(database), mode]
    if code_csv:
        command.append(code_csv)
    proc = subprocess.Popen(
        command,
        cwd=ctx.workspace,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    assert proc.stdout is not None
    header_seen = False
    for raw_line in proc.stdout:
        line = raw_line.rstrip("\n")
        if not header_seen:
            header_seen = True
            continue
        yield line.split("\t")
    assert proc.stderr is not None
    stderr = proc.stderr.read()
    return_code = proc.wait()
    if return_code != 0:
        raise RuntimeError(f"BridgeDb Derby exporter failed with exit {return_code}: {stderr.strip()}")


def read_bridgedb_info(ctx: BuildContext, bridge_file: Path) -> dict[str, str]:
    try:
        rows = list(run_bridgedb_exporter(ctx, bridge_file, "info"))
    except Exception:
        return {}
    if not rows:
        return {}
    values = rows[0] + [""] * 6
    return {
        "schema_version": values[0],
        "build_date": values[1],
        "data_source_name": values[2],
        "data_source_version": values[3],
        "data_type": values[4],
        "series": values[5],
    }


def ingest_bridgedb_metabolites(
    ctx: BuildContext,
    xref_writer: TableWriter,
    uid_by_xref: dict[str, str],
    seen_xrefs: set[tuple[str, str]],
) -> None:
    bridge_dir = ctx.source_dir("bridgedb_metabolites")
    bridge_files = sorted(bridge_dir.glob("metabolites_*.bridge"))
    if not bridge_files:
        return
    latest = bridge_files[-1]
    try:
        with zipfile.ZipFile(latest) as archive:
            service = archive.read("database/service.properties").decode("latin1", errors="replace")
    except (KeyError, zipfile.BadZipFile, OSError) as exc:
        ctx.note("warning", "bridgedb_metabolites", f"Could not inspect BridgeDb file {latest.name}: {exc}")
        return
    if "derby.serviceProtocol" not in service:
        ctx.note("warning", "bridgedb_metabolites", f"BridgeDb file {latest.name} did not look like a Derby database.")
        return

    source_info = read_bridgedb_info(ctx, latest)
    source_release = f"bridgedb:{source_info.get('build_date') or latest.stem.rsplit('_', 1)[-1]}"
    license_id = ctx.license_id("bridgedb_metabolites")
    anchor_uid_by_xref = dict(uid_by_xref)
    scanned = 0
    added = 0
    conflicts = 0
    skipped = 0
    xref_rows: list[dict[str, Any]] = []
    try:
        for parts in run_bridgedb_exporter(ctx, latest, "links", ",".join(BRIDGEDB_DIRECT_ANCHOR_CODES)):
            scanned += 1
            values = parts + [""] * 5
            left_key = bridgedb_xref(values[1], values[0])
            right_key = bridgedb_xref(values[3], values[2])
            if not left_key or not right_key:
                skipped += 1
                continue
            left_uid = anchor_uid_by_xref.get(left_key, "")
            right_uid = anchor_uid_by_xref.get(right_key, "")
            if left_uid and right_uid and left_uid != right_uid:
                conflicts += 1
                continue
            candidate_pairs: list[tuple[str, str]] = []
            if left_uid:
                candidate_pairs.append((left_uid, right_key))
            if right_uid:
                candidate_pairs.append((right_uid, left_key))
            if left_uid and right_uid and left_uid == right_uid:
                candidate_pairs.append((left_uid, left_key))
                candidate_pairs.append((left_uid, right_key))
            if not candidate_pairs:
                skipped += 1
                continue
            for uid, xref in candidate_pairs:
                parsed = normalize_xref(xref)
                if not parsed:
                    continue
                key = (uid, parsed[2])
                if key in seen_xrefs:
                    continue
                seen_xrefs.add(key)
                uid_by_xref[parsed[2]] = uid
                xref_rows.append(
                    {
                        "xref_uid": stable_uid("xref", uid, parsed[2]),
                        "metabolite_uid": uid,
                        "xref_source": parsed[0],
                        "xref_id": parsed[1],
                        "xref_key": parsed[2],
                        "source_name": "BridgeDb",
                        "source_release": source_release,
                        "license_id": license_id,
                        "parser_hash": ctx.parser_hash,
                    }
                )
                added += 1
                if len(xref_rows) >= 50_000:
                    xref_writer.write(xref_rows)
                    xref_rows = []
        xref_writer.write(xref_rows)
        ctx.note(
            "info",
            "bridgedb_metabolites",
            "Parsed "
            f"{latest.name} with Java/Derby; scanned {scanned} direct anchor links, added {added} xrefs, "
            f"skipped {skipped}, conflicts {conflicts}. Source {source_info.get('data_source_version', source_release)}.",
        )
    except Exception as exc:
        ctx.note("warning", "bridgedb_metabolites", f"BridgeDb Derby extraction failed for {latest.name}: {exc}")


def build_metabolites(ctx: BuildContext, include_pubchem_enrichment: bool) -> dict[str, str]:
    started = time.time()
    metabolite_writer = ctx.writer("metabolites")
    xref_writer = ctx.writer("metabolite_xrefs")
    hmdb_writer = ctx.writer("hmdb_metabolite_details")
    chebi_dir = ctx.source_dir("chebi")
    lipid_dir = ctx.source_dir("lipidmaps") / "rest"
    chebi_release = ctx.source_release("chebi")
    chebi_license = ctx.license_id("chebi")
    lipid_release = ctx.source_release("lipidmaps")
    lipid_license = ctx.license_id("lipidmaps")
    uid_by_xref: dict[str, str] = {}
    seen_metabolites: set[str] = set()
    seen_xrefs: set[tuple[str, str]] = set()

    compounds: dict[str, dict[str, str]] = {}
    for row in iter_tsv(chebi_dir / "compounds.tsv.gz"):
        compounds[row["id"]] = row

    synonyms: dict[str, list[str]] = defaultdict(list)
    for row in iter_tsv(chebi_dir / "names.tsv.gz"):
        if row.get("language_code") and row.get("language_code") != "en":
            continue
        synonyms[row["compound_id"]].append(row.get("ascii_name") or row.get("name") or "")

    chemical: dict[str, dict[str, str]] = {}
    for row in iter_tsv(chebi_dir / "chemical_data.tsv.gz"):
        current = chemical.get(row["compound_id"])
        preferred = row.get("is_autogenerated", "").lower() == "false"
        if current is None or preferred:
            chemical[row["compound_id"]] = row

    structures: dict[str, dict[str, str]] = {}
    for row in iter_tsv(chebi_dir / "structures.tsv.gz"):
        current = structures.get(row["compound_id"])
        preferred = row.get("default_structure", "").lower() == "true"
        if current is None or preferred:
            structures[row["compound_id"]] = row

    obo_xrefs = parse_chebi_obo_xrefs(chebi_dir / "ontology" / "chebi.obo")

    database_xrefs: dict[str, set[str]] = defaultdict(set)
    for row in iter_tsv(chebi_dir / "database_accession.tsv.gz"):
        key = chebi_key(row["compound_id"])
        accession = row.get("accession_number", "").strip()
        if not accession:
            continue
        source_id = row.get("source_id", "").strip()
        xref = f"CHEBI.DB_SOURCE_{source_id}:{accession}" if source_id else f"CHEBI.DB_SOURCE:{accession}"
        database_xrefs[key].add(xref)

    batch = []
    for compound_id, compound in compounds.items():
        chebi = compound.get("chebi_accession") or chebi_key(compound_id)
        uid = stable_uid("metabolite", chebi)
        uid_by_xref[chebi] = uid
        xrefs = set(obo_xrefs.get(chebi, set()))
        xrefs.add(chebi)
        xrefs.update(database_xrefs.get(chebi, set()))
        for x in xrefs:
            uid_by_xref[x] = uid
        chem = chemical.get(compound_id, {})
        struct = structures.get(compound_id, {})
        canonical = compound.get("ascii_name") or compound.get("name") or chebi
        syns = dedupe([s for s in synonyms.get(compound_id, []) if s and s != canonical])
        external_xrefs = sorted(xrefs)
        row = {
            "metabolite_uid": uid,
            "canonical_name": canonical,
            "synonyms": syns,
            "formula": chem.get("formula") or "",
            "exact_mass": maybe_float(chem.get("monoisotopic_mass") or chem.get("mass")),
            "charge": maybe_int(chem.get("charge")),
            "inchikey": struct.get("standard_inchi_key") or "",
            "smiles": struct.get("smiles") or "",
            "external_xrefs": external_xrefs,
            "source_priority": "ChEBI",
            "source_release": chebi_release,
            "license_id": chebi_license,
            "checksum": ctx.checksum("chebi"),
            "parser_hash": ctx.parser_hash,
        }
        batch.append(row)
        seen_metabolites.add(uid)
        for x in external_xrefs:
            seen_xrefs.add((uid, x))
        if len(batch) >= 25_000:
            metabolite_writer.write(batch)
            batch = []
    metabolite_writer.write(batch)

    xref_rows = []
    for uid, x in sorted(seen_xrefs):
        parsed = normalize_xref(x)
        if not parsed:
            continue
        source, external_id, key = parsed
        xref_rows.append(
            {
                "xref_uid": stable_uid("xref", uid, key),
                "metabolite_uid": uid,
                "xref_source": source,
                "xref_id": external_id,
                "xref_key": key,
                "source_name": "ChEBI",
                "source_release": chebi_release,
                "license_id": chebi_license,
                "parser_hash": ctx.parser_hash,
            }
        )
        if len(xref_rows) >= 50_000:
            xref_writer.write(xref_rows)
            xref_rows = []
    xref_writer.write(xref_rows)

    lipid_new_rows = []
    lipid_xref_rows = []
    for lipid_file in sorted(lipid_dir.glob("LM*.tsv")):
        with lipid_file.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            first_line = handle.readline()
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                lm_id = row.get("lm_id", "").strip()
                if not lm_id:
                    continue
                candidate_xrefs = [
                    f"LIPIDMAPS:{lm_id}",
                    f"CHEBI:{row.get('chebi_id', '').strip()}" if row.get("chebi_id", "").strip() else "",
                    f"KEGG.COMPOUND:{row.get('kegg_id', '').strip()}" if row.get("kegg_id", "").strip() else "",
                    f"HMDB:{row.get('hmdb_id', '').strip()}" if row.get("hmdb_id", "").strip() else "",
                    f"PUBCHEM.COMPOUND:{row.get('pubchem_cid', '').strip()}" if row.get("pubchem_cid", "").strip() else "",
                ]
                resolved_uid = ""
                for candidate in candidate_xrefs:
                    if candidate and candidate in uid_by_xref:
                        resolved_uid = uid_by_xref[candidate]
                        break
                if not resolved_uid:
                    resolved_uid = stable_uid("metabolite", "LIPIDMAPS", lm_id)
                    uid_by_xref[f"LIPIDMAPS:{lm_id}"] = resolved_uid
                    if resolved_uid not in seen_metabolites:
                        lipid_new_rows.append(
                            {
                                "metabolite_uid": resolved_uid,
                                "canonical_name": row.get("name") or lm_id,
                                "synonyms": dedupe((row.get("synonyms") or "").split("|")),
                                "formula": row.get("formula") or "",
                                "exact_mass": maybe_float(row.get("exactmass")),
                                "charge": None,
                                "inchikey": row.get("inchi_key") or "",
                                "smiles": row.get("smiles") or "",
                                "external_xrefs": dedupe(candidate_xrefs),
                                "source_priority": "LIPID MAPS",
                                "source_release": first_line.strip() or lipid_release,
                                "license_id": lipid_license,
                                "checksum": ctx.checksum("lipidmaps"),
                                "parser_hash": ctx.parser_hash,
                            }
                        )
                        seen_metabolites.add(resolved_uid)
                for candidate in candidate_xrefs:
                    parsed = normalize_xref(candidate) if candidate else None
                    if not parsed:
                        continue
                    uid_by_xref[parsed[2]] = resolved_uid
                    key = (resolved_uid, parsed[2])
                    if key in seen_xrefs:
                        continue
                    seen_xrefs.add(key)
                    lipid_xref_rows.append(
                        {
                            "xref_uid": stable_uid("xref", resolved_uid, parsed[2]),
                            "metabolite_uid": resolved_uid,
                            "xref_source": parsed[0],
                            "xref_id": parsed[1],
                            "xref_key": parsed[2],
                            "source_name": "LIPID MAPS",
                            "source_release": first_line.strip() or lipid_release,
                            "license_id": lipid_license,
                            "parser_hash": ctx.parser_hash,
                        }
                    )
                if len(lipid_new_rows) >= 10_000:
                    metabolite_writer.write(lipid_new_rows)
                    lipid_new_rows = []
                if len(lipid_xref_rows) >= 25_000:
                    xref_writer.write(lipid_xref_rows)
                    lipid_xref_rows = []
    metabolite_writer.write(lipid_new_rows)
    xref_writer.write(lipid_xref_rows)

    ingest_hmdb_metabolites(
        ctx,
        metabolite_writer,
        xref_writer,
        hmdb_writer,
        uid_by_xref,
        seen_metabolites,
        seen_xrefs,
    )
    ingest_bridgedb_metabolites(ctx, xref_writer, uid_by_xref, seen_xrefs)
    if include_pubchem_enrichment:
        ctx.note(
            "info",
            "pubchem_compound_extras",
            "PubChem enrichment is implemented as a separate derived cache. Run scripts/build_pubchem_cid_cache.py after the normalized build; PubChem CIDs are retained here as metabolite xrefs.",
        )
    else:
        ctx.note(
            "info",
            "pubchem_compound_extras",
            "PubChem Compound Extras are represented here through PubChem CID xrefs. Use scripts/build_pubchem_cid_cache.py to materialize CID title, formula, mass, InChIKey, SMILES, and synonym properties.",
        )
    ctx.close_writer(metabolite_writer, started)
    ctx.close_writer(xref_writer, started)
    ctx.close_writer(hmdb_writer, started)
    return uid_by_xref


def build_genes(ctx: BuildContext) -> tuple[dict[str, str], dict[str, str]]:
    started = time.time()
    gene_writer = ctx.writer("genes")
    xref_writer = ctx.writer("gene_xrefs")
    ensembl_dir = ctx.source_dir("ensembl_human")
    ncbi_dir = ctx.source_dir("ncbi_gene_human")
    genes: dict[str, dict[str, Any]] = {}
    xrefs: dict[str, set[str]] = defaultdict(set)
    uid_by_ensembl: dict[str, str] = {}
    uid_by_entrez: dict[str, str] = {}

    gtf_files = sorted(ensembl_dir.glob("*.gtf.gz"))
    if gtf_files:
        with gzip.open(gtf_files[0], "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9 or parts[2] != "gene":
                    continue
                attrs = parse_gtf_attrs(parts[8])
                ensembl_id = attrs.get("gene_id", "").split(".", 1)[0]
                if not ensembl_id:
                    continue
                uid = stable_uid("gene", "ENSEMBL", ensembl_id)
                uid_by_ensembl[ensembl_id] = uid
                xrefs[uid].add(f"ENSEMBL:{ensembl_id}")
                genes[uid] = {
                    "gene_uid": uid,
                    "ensembl_gene_id": ensembl_id,
                    "entrez_gene_id": "",
                    "symbol": attrs.get("gene_name", ""),
                    "aliases": [],
                    "description": "",
                    "taxon": "9606",
                    "biotype": attrs.get("gene_biotype", ""),
                    "chromosome": parts[0],
                    "start": maybe_int(parts[3]),
                    "end": maybe_int(parts[4]),
                    "strand": parts[6],
                    "uniprot_ids": [],
                    "external_xrefs": [],
                    "source_release": ctx.source_release("ensembl_human"),
                    "license_id": ctx.license_id("ensembl_human"),
                    "checksum": ctx.checksum("ensembl_human"),
                    "parser_hash": ctx.parser_hash,
                }

    gene_info = ncbi_dir / "Homo_sapiens.gene_info.gz"
    if gene_info.exists():
        for row in iter_tsv(gene_info):
            if row.get("#tax_id") != "9606":
                continue
            entrez = row.get("GeneID", "").strip()
            dbx = dedupe((row.get("dbXrefs") or "").split("|"))
            ensembl_ids = [x.split(":", 1)[1].split(".", 1)[0] for x in dbx if x.startswith("Ensembl:")]
            uid = ""
            for ensembl_id in ensembl_ids:
                uid = uid_by_ensembl.get(ensembl_id, "")
                if uid:
                    break
            if not uid:
                uid = stable_uid("gene", "NCBI.GENE", entrez)
            if entrez:
                uid_by_entrez[entrez] = uid
                xrefs[uid].add(f"NCBI.GENE:{entrez}")
            for ensembl_id in ensembl_ids:
                uid_by_ensembl[ensembl_id] = uid
                xrefs[uid].add(f"ENSEMBL:{ensembl_id}")
            for raw in dbx:
                parsed = normalize_xref(raw)
                if parsed:
                    xrefs[uid].add(parsed[2])
            aliases = dedupe((row.get("Synonyms") or "").split("|") + (row.get("Other_designations") or "").split("|"))
            gene_row = genes.get(
                uid,
                {
                    "gene_uid": uid,
                    "ensembl_gene_id": ensembl_ids[0] if ensembl_ids else "",
                    "entrez_gene_id": entrez,
                    "symbol": "",
                    "aliases": [],
                    "description": "",
                    "taxon": "9606",
                    "biotype": row.get("type_of_gene", ""),
                    "chromosome": row.get("chromosome", ""),
                    "start": None,
                    "end": None,
                    "strand": "",
                    "uniprot_ids": [],
                    "external_xrefs": [],
                    "source_release": ctx.source_release("ncbi_gene_human"),
                    "license_id": ctx.license_id("ncbi_gene_human"),
                    "checksum": ctx.checksum("ncbi_gene_human"),
                    "parser_hash": ctx.parser_hash,
                },
            )
            gene_row["entrez_gene_id"] = gene_row.get("entrez_gene_id") or entrez
            gene_row["symbol"] = gene_row.get("symbol") or row.get("Symbol") or row.get("Symbol_from_nomenclature_authority")
            gene_row["description"] = gene_row.get("description") or row.get("description") or row.get("Full_name_from_nomenclature_authority")
            gene_row["aliases"] = dedupe(list(gene_row.get("aliases", [])) + aliases)
            genes[uid] = gene_row

    gene2ensembl = ncbi_dir / "gene2ensembl.gz"
    if gene2ensembl.exists():
        for row in iter_tsv(gene2ensembl):
            if row.get("#tax_id") != "9606":
                continue
            entrez = row.get("GeneID", "").strip()
            ensembl_gene = row.get("Ensembl_gene_identifier", "").split(".", 1)[0]
            uid = uid_by_entrez.get(entrez) or uid_by_ensembl.get(ensembl_gene)
            if not uid:
                continue
            if entrez:
                uid_by_entrez[entrez] = uid
                xrefs[uid].add(f"NCBI.GENE:{entrez}")
            if ensembl_gene and ensembl_gene != "-":
                uid_by_ensembl[ensembl_gene] = uid
                xrefs[uid].add(f"ENSEMBL:{ensembl_gene}")
            for field_name, source_name in [
                ("RNA_nucleotide_accession.version", "REFSEQ"),
                ("protein_accession.version", "REFSEQ"),
                ("Ensembl_rna_identifier", "ENSEMBL.TRANSCRIPT"),
                ("Ensembl_protein_identifier", "ENSEMBL.PROTEIN"),
            ]:
                value = row.get(field_name, "").strip()
                if value and value != "-":
                    xrefs[uid].add(f"{source_name}:{value}")

    gene_rows = []
    xref_rows = []
    for uid, row in sorted(genes.items()):
        row["external_xrefs"] = sorted(xrefs.get(uid, set()))
        row["uniprot_ids"] = sorted({x.split(":", 1)[1] for x in row["external_xrefs"] if x.startswith("UNIPROT:")})
        gene_rows.append(row)
        for x in row["external_xrefs"]:
            parsed = normalize_xref(x)
            if not parsed:
                continue
            xref_rows.append(
                {
                    "xref_uid": stable_uid("xref", uid, parsed[2]),
                    "gene_uid": uid,
                    "xref_source": parsed[0],
                    "xref_id": parsed[1],
                    "xref_key": parsed[2],
                    "source_name": "Ensembl/NCBI Gene",
                    "source_release": f"{ctx.source_release('ensembl_human')}+{ctx.source_release('ncbi_gene_human')}",
                    "license_id": f"{ctx.license_id('ensembl_human')}|{ctx.license_id('ncbi_gene_human')}",
                    "parser_hash": ctx.parser_hash,
                }
            )
        if len(gene_rows) >= 25_000:
            gene_writer.write(gene_rows)
            gene_rows = []
        if len(xref_rows) >= 50_000:
            xref_writer.write(xref_rows)
            xref_rows = []
    gene_writer.write(gene_rows)
    xref_writer.write(xref_rows)
    ctx.close_writer(gene_writer, started)
    ctx.close_writer(xref_writer, started)
    return uid_by_ensembl, uid_by_entrez


def strip_reactome_version(value: str) -> str:
    text = str(value or "").strip()
    if text.startswith("REACTOME:"):
        text = text.split(":", 1)[1]
    return text.split(".", 1)[0]


def iter_rhea_tsv_rows(ctx: BuildContext, filename: str) -> Iterator[dict[str, str]]:
    rhea_dir = ctx.source_dir("rhea")
    for candidate in [
        rhea_dir / filename,
        rhea_dir / "tsv" / filename,
    ]:
        if candidate.exists():
            yield from iter_tsv(candidate)
            return

    for archive_path in [
        rhea_dir / "rhea-tsv.tar.gz",
        rhea_dir / "tsv" / "rhea-tsv.tar.gz",
    ]:
        if not archive_path.exists():
            continue
        try:
            with tarfile.open(archive_path, "r:gz") as archive:
                member = next((item for item in archive.getmembers() if Path(item.name).name == filename), None)
                if member is None:
                    continue
                extracted = archive.extractfile(member)
                if extracted is None:
                    continue
                with io.TextIOWrapper(extracted, encoding="utf-8", errors="replace", newline="") as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    for row in reader:
                        yield {key: (value if value is not None else "") for key, value in row.items()}
                return
        except (OSError, tarfile.TarError):
            ctx.note("warning", "rhea", f"Could not read {filename} from {archive_path}.")
            return


def build_pathways(
    ctx: BuildContext,
    metabolite_uid_by_xref: dict[str, str],
    gene_uid_by_ensembl: dict[str, str],
    gene_uid_by_entrez: dict[str, str],
) -> tuple[dict[str, str], dict[str, str]]:
    started = time.time()
    pathway_writer = ctx.writer("pathways")
    reaction_writer = ctx.writer("reactions")
    hierarchy_writer = ctx.writer("pathway_hierarchy_edges")
    participant_writer = ctx.writer("reaction_participants")
    gene_edge_writer = ctx.writer("gene_pathway_edges")
    metabolite_edge_writer = ctx.writer("metabolite_pathway_edges")
    reactome_dir = ctx.source_dir("reactome")
    wiki_dir = ctx.source_dir("wikipathways")
    reactome_release = ctx.source_release("reactome")
    reactome_license = ctx.license_id("reactome")
    rhea_release = ctx.source_release("rhea")
    rhea_license = ctx.license_id("rhea")
    wiki_release = ctx.source_release("wikipathways")
    wiki_license = ctx.license_id("wikipathways")
    pathway_uid_by_ext: dict[str, str] = {}
    reaction_uid_by_ext: dict[str, str] = {}
    pathways: dict[str, dict[str, Any]] = {}
    reactions: dict[str, dict[str, Any]] = {}

    path_file = reactome_dir / "ReactomePathways.txt"
    if path_file.exists():
        with path_file.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for parts in reader:
                if len(parts) < 3 or parts[2] != "Homo sapiens":
                    continue
                rid, name, species = parts[0], parts[1], parts[2]
                uid = stable_uid("pathway", "REACTOME", rid)
                pathway_uid_by_ext[rid] = uid
                pathways[uid] = {
                    "pathway_uid": uid,
                    "name": name,
                    "species": species,
                    "source_name": "Reactome",
                    "primary_external_id": rid,
                    "hierarchy_path": "",
                    "external_xrefs": [f"REACTOME:{rid}"],
                    "source_release": reactome_release,
                    "license_id": reactome_license,
                    "checksum": ctx.checksum("reactome"),
                    "parser_hash": ctx.parser_hash,
                }

    hierarchy_rows = []
    rel_file = reactome_dir / "ReactomePathwaysRelation.txt"
    if rel_file.exists():
        with rel_file.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for parent, child, *_ in reader:
                if parent in pathway_uid_by_ext and child in pathway_uid_by_ext:
                    hierarchy_rows.append(
                        {
                            "edge_uid": stable_uid("edge", "reactome_parent", parent, child),
                            "parent_pathway_uid": pathway_uid_by_ext[parent],
                            "child_pathway_uid": pathway_uid_by_ext[child],
                            "parent_external_id": parent,
                            "child_external_id": child,
                            "predicate": "parent_of",
                            "source_name": "Reactome",
                            "source_release": reactome_release,
                            "license_id": reactome_license,
                            "parser_hash": ctx.parser_hash,
                        }
                    )
                if len(hierarchy_rows) >= 50_000:
                    hierarchy_writer.write(hierarchy_rows)
                    hierarchy_rows = []
    hierarchy_writer.write(hierarchy_rows)

    reaction_participant_rows = []
    reactome_uniprot_ids: set[str] = set()
    rxn_file = reactome_dir / "reactome_reaction_exporter.txt"
    if rxn_file.exists():
        with rxn_file.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                pathway_id = row.get("pathway_id", "")
                reaction_id = row.get("reaction_id", "")
                if not pathway_id.startswith("R-HSA-") or not reaction_id.startswith("R-HSA-"):
                    continue
                pathway_uid = pathway_uid_by_ext.get(pathway_id, "")
                if not pathway_uid:
                    continue
                reaction_uid = reaction_uid_by_ext.get(reaction_id) or stable_uid("reaction", "REACTOME", reaction_id)
                reaction_uid_by_ext[reaction_id] = reaction_uid
                reactions.setdefault(
                    reaction_uid,
                    {
                        "reaction_uid": reaction_uid,
                        "name": row.get("reaction_name", ""),
                        "species": "Homo sapiens",
                        "primary_external_id": reaction_id,
                        "pathway_uid": pathway_uid,
                        "pathway_external_id": pathway_id,
                        "source_name": "Reactome",
                        "external_xrefs": [f"REACTOME:{reaction_id}"],
                        "source_release": reactome_release,
                        "license_id": reactome_license,
                        "checksum": ctx.checksum("reactome"),
                        "parser_hash": ctx.parser_hash,
                    },
                )
                uniprot = row.get("uniprot_acc", "").strip()
                if uniprot:
                    reactome_uniprot_ids.add(uniprot)
                    reaction_participant_rows.append(
                        {
                            "edge_uid": stable_uid("edge", "reaction_uniprot", reaction_id, uniprot, row.get("role_in_reaction", "")),
                            "reaction_uid": reaction_uid,
                            "participant_uid": "",
                            "participant_external_id": f"UNIPROT:{uniprot}",
                            "participant_type": "protein",
                            "role": row.get("role_in_reaction", ""),
                            "source_name": "Reactome",
                            "source_record_id": f"{pathway_id}|{reaction_id}|{uniprot}",
                            "source_release": reactome_release,
                            "license_id": reactome_license,
                            "parser_hash": ctx.parser_hash,
                        }
                    )
                if len(reaction_participant_rows) >= 50_000:
                    participant_writer.write(reaction_participant_rows)
                    reaction_participant_rows = []
    participant_writer.write(reaction_participant_rows)

    rhea_reaction_uids_by_master: dict[str, set[str]] = defaultdict(set)
    rhea_xref_count = 0
    for row in iter_rhea_tsv_rows(ctx, "rhea2reactome.tsv"):
        reactome_id = strip_reactome_version(row.get("ID", ""))
        reaction_uid = reaction_uid_by_ext.get(reactome_id)
        if not reaction_uid or reaction_uid not in reactions:
            continue
        master_id = str(row.get("MASTER_ID") or row.get("RHEA_ID") or "").strip()
        rhea_id = str(row.get("RHEA_ID") or "").strip()
        if master_id:
            rhea_reaction_uids_by_master[master_id].add(reaction_uid)
        xrefs = set(reactions[reaction_uid].get("external_xrefs") or [])
        before = len(xrefs)
        if master_id:
            xrefs.add(f"RHEA:{master_id}")
        if rhea_id:
            xrefs.add(f"RHEA:{rhea_id}")
        if len(xrefs) > before:
            reactions[reaction_uid]["external_xrefs"] = sorted(xrefs)
            rhea_xref_count += len(xrefs) - before

    rhea_participant_rows = []
    rhea_enzyme_edges = 0
    if rhea_reaction_uids_by_master and reactome_uniprot_ids:
        for row in iter_rhea_tsv_rows(ctx, "rhea2uniprot_sprot.tsv"):
            uniprot = str(row.get("ID") or "").strip()
            master_id = str(row.get("MASTER_ID") or row.get("RHEA_ID") or "").strip()
            if not uniprot or uniprot not in reactome_uniprot_ids or not master_id:
                continue
            for reaction_uid in sorted(rhea_reaction_uids_by_master.get(master_id, set())):
                rhea_participant_rows.append(
                    {
                        "edge_uid": stable_uid("edge", "rhea_uniprot", reaction_uid, uniprot, master_id),
                        "reaction_uid": reaction_uid,
                        "participant_uid": "",
                        "participant_external_id": f"UNIPROT:{uniprot}",
                        "participant_type": "protein",
                        "role": "enzyme_or_transporter",
                        "source_name": "Rhea",
                        "source_record_id": f"RHEA:{master_id}|UNIPROT:{uniprot}",
                        "source_release": rhea_release,
                        "license_id": rhea_license,
                        "parser_hash": ctx.parser_hash,
                    }
                )
                rhea_enzyme_edges += 1
                if len(rhea_participant_rows) >= 50_000:
                    participant_writer.write(rhea_participant_rows)
                    rhea_participant_rows = []
    participant_writer.write(rhea_participant_rows)
    if rhea_xref_count or rhea_enzyme_edges:
        ctx.note(
            "info",
            "rhea",
            f"Integrated Rhea cross-references into {len(rhea_reaction_uids_by_master)} Reactome-linked master reactions; added {rhea_xref_count} reaction xrefs and {rhea_enzyme_edges} human protein-reaction support edges.",
        )

    gene_edge_rows = []
    for reactome_name, gene_map, prefix in [
        ("Ensembl2Reactome.txt", gene_uid_by_ensembl, "ENSEMBL"),
        ("NCBI2Reactome.txt", gene_uid_by_entrez, "NCBI.GENE"),
    ]:
        file_path = reactome_dir / reactome_name
        if not file_path.exists():
            continue
        with file_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for parts in reader:
                if len(parts) < 6 or parts[5] != "Homo sapiens":
                    continue
                gene_id, pathway_id, _, _, evidence_code, species = parts[:6]
                pathway_uid = pathway_uid_by_ext.get(pathway_id)
                gene_uid = gene_map.get(gene_id.split(".", 1)[0])
                if not pathway_uid or not gene_uid:
                    continue
                gene_edge_rows.append(
                    {
                        "edge_uid": stable_uid("edge", "gene_pathway", gene_uid, pathway_uid, reactome_name),
                        "gene_uid": gene_uid,
                        "pathway_uid": pathway_uid,
                        "gene_external_id": f"{prefix}:{gene_id}",
                        "pathway_external_id": pathway_id,
                        "predicate": "involved_in",
                        "source_name": "Reactome",
                        "source_record_id": f"{gene_id}|{pathway_id}",
                        "evidence_level": "curated" if evidence_code != "IEA" else "inferred",
                        "evidence_code": evidence_code,
                        "species": species,
                        "source_release": reactome_release,
                        "license_id": reactome_license,
                        "parser_hash": ctx.parser_hash,
                    }
                )
                if len(gene_edge_rows) >= 100_000:
                    gene_edge_writer.write(gene_edge_rows)
                    gene_edge_rows = []
    gene_edge_writer.write(gene_edge_rows)

    metabolite_edge_rows = []
    chebi_reactome = reactome_dir / "ChEBI2Reactome.txt"
    if chebi_reactome.exists():
        with chebi_reactome.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for parts in reader:
                if len(parts) < 6 or parts[5] != "Homo sapiens":
                    continue
                chebi_id, pathway_id, _, _, evidence_code, species = parts[:6]
                metabolite_uid = metabolite_uid_by_xref.get(chebi_key(chebi_id))
                pathway_uid = pathway_uid_by_ext.get(pathway_id)
                if not metabolite_uid or not pathway_uid:
                    continue
                metabolite_edge_rows.append(
                    {
                        "edge_uid": stable_uid("edge", "metabolite_pathway", metabolite_uid, pathway_uid, "reactome"),
                        "metabolite_uid": metabolite_uid,
                        "pathway_uid": pathway_uid,
                        "metabolite_external_id": chebi_key(chebi_id),
                        "pathway_external_id": pathway_id,
                        "predicate": "participates_in",
                        "source_name": "Reactome",
                        "source_record_id": f"{chebi_id}|{pathway_id}",
                        "evidence_level": "curated" if evidence_code != "IEA" else "inferred",
                        "evidence_code": evidence_code,
                        "species": species,
                        "source_release": reactome_release,
                        "license_id": reactome_license,
                        "parser_hash": ctx.parser_hash,
                    }
                )
                if len(metabolite_edge_rows) >= 100_000:
                    metabolite_edge_writer.write(metabolite_edge_rows)
                    metabolite_edge_rows = []
    metabolite_edge_writer.write(metabolite_edge_rows)

    wiki_files = sorted(wiki_dir.glob("*gmt-Homo_sapiens.gmt"))
    wiki_gene_edges = []
    for wiki_file in wiki_files:
        with wiki_file.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for parts in reader:
                if len(parts) < 3:
                    continue
                desc = parts[0].split("%")
                if len(desc) < 4 or desc[3] != "Homo sapiens":
                    continue
                name, release, wp_id, species = desc[:4]
                uid = stable_uid("pathway", "WIKIPATHWAYS", wp_id)
                pathway_uid_by_ext[wp_id] = uid
                pathways.setdefault(
                    uid,
                    {
                        "pathway_uid": uid,
                        "name": name,
                        "species": species,
                        "source_name": "WikiPathways",
                        "primary_external_id": wp_id,
                        "hierarchy_path": "",
                        "external_xrefs": [f"WIKIPATHWAYS:{wp_id}", parts[1]],
                        "source_release": release or wiki_release,
                        "license_id": wiki_license,
                        "checksum": ctx.checksum("wikipathways"),
                        "parser_hash": ctx.parser_hash,
                    },
                )
                for entrez in parts[2:]:
                    gene_uid = gene_uid_by_entrez.get(entrez)
                    if not gene_uid:
                        continue
                    wiki_gene_edges.append(
                        {
                            "edge_uid": stable_uid("edge", "gene_pathway", gene_uid, uid, "wikipathways"),
                            "gene_uid": gene_uid,
                            "pathway_uid": uid,
                            "gene_external_id": f"NCBI.GENE:{entrez}",
                            "pathway_external_id": wp_id,
                            "predicate": "involved_in",
                            "source_name": "WikiPathways",
                            "source_record_id": f"{wp_id}|{entrez}",
                            "evidence_level": "curated",
                            "evidence_code": "pathway_member",
                            "species": species,
                            "source_release": release or wiki_release,
                            "license_id": wiki_license,
                            "parser_hash": ctx.parser_hash,
                        }
                    )
                if len(wiki_gene_edges) >= 50_000:
                    gene_edge_writer.write(wiki_gene_edges)
                    wiki_gene_edges = []
    gene_edge_writer.write(wiki_gene_edges)

    pathway_writer.write(pathways.values())
    reaction_writer.write(reactions.values())
    ctx.close_writer(pathway_writer, started)
    ctx.close_writer(reaction_writer, started)
    ctx.close_writer(hierarchy_writer, started)
    ctx.close_writer(participant_writer, started)
    ctx.close_writer(gene_edge_writer, started)
    ctx.close_writer(metabolite_edge_writer, started)
    return pathway_uid_by_ext, reaction_uid_by_ext


def disease_uid_for_external(external_key: str) -> str:
    return stable_uid("disease", external_key)


def parse_obo_disease_terms(path: Path, source_priority: str, ctx: BuildContext) -> list[dict[str, Any]]:
    rows = []
    for term in iter_obo_terms(path):
        ids = term.get("id", [])
        names = term.get("name", [])
        if not ids or not names or term.get("is_obsolete", ["false"])[0] == "true":
            continue
        primary = ids[0]
        parsed_primary = normalize_xref(primary)
        if not parsed_primary:
            continue
        xrefs = {parsed_primary[2]}
        for alt_id in term.get("alt_id", []):
            parsed = normalize_xref(alt_id)
            if parsed:
                xrefs.add(parsed[2])
        for raw_xref in term.get("xref", []):
            parsed = normalize_xref(raw_xref)
            if parsed:
                xrefs.add(parsed[2])
        aliases = [quoted_value(raw) for raw in term.get("synonym", [])]
        parents = []
        for parent in term.get("is_a", []):
            parsed = normalize_xref(parent)
            if parsed:
                parents.append(parsed[2])
        rows.append(
            {
                "primary_external_id": parsed_primary[2],
                "name": names[0],
                "aliases": dedupe(aliases),
                "description": quoted_value(term.get("def", [""])[0]) if term.get("def") else "",
                "parents": dedupe(parents),
                "external_xrefs": sorted(xrefs),
                "source_priority": source_priority,
            }
        )
    return rows


def is_oncology_disease_term(name: str, aliases: Iterable[str], subsets: Iterable[str] = ()) -> bool:
    if any("oncotree" in str(subset).casefold() for subset in subsets):
        return True
    labels = [name, *list(aliases)]
    return any(ONCOLOGY_DISEASE_RE.search(str(label or "")) for label in labels)


def is_ncit_non_disease_concept(name: str, description: str = "") -> bool:
    return bool(
        NCIT_NON_DISEASE_CONCEPT_NAME_RE.search(name)
        or NCIT_NON_DISEASE_CONCEPT_DESCRIPTION_RE.search(description)
    )


def parse_ncit_cancer_terms(path: Path, ctx: BuildContext) -> list[dict[str, Any]]:
    rows = []
    for term in iter_obo_terms(path):
        ids = term.get("id", [])
        names = term.get("name", [])
        if not ids or not names or term.get("is_obsolete", ["false"])[0] == "true":
            continue
        parsed_primary = normalize_xref(ids[0])
        if not parsed_primary or parsed_primary[0] != "NCIT":
            continue
        aliases = [quoted_value(raw) for raw in term.get("synonym", [])]
        if not is_oncology_disease_term(names[0], aliases, term.get("subset", [])):
            continue
        description = quoted_value(term.get("def", [""])[0]) if term.get("def") else ""
        if is_ncit_non_disease_concept(names[0], description):
            continue
        xrefs = {parsed_primary[2]}
        for alt_id in term.get("alt_id", []):
            parsed = normalize_xref(alt_id)
            if parsed:
                xrefs.add(parsed[2])
        for raw_xref in term.get("xref", []):
            parsed = normalize_xref(raw_xref)
            if parsed:
                xrefs.add(parsed[2])
        parents = []
        for parent in term.get("is_a", []):
            parsed = normalize_xref(parent)
            if parsed:
                parents.append(parsed[2])
        rows.append(
            {
                "primary_external_id": parsed_primary[2],
                "name": names[0],
                "aliases": dedupe(aliases),
                "description": description,
                "parents": dedupe(parents),
                "external_xrefs": sorted(xrefs),
                "source_priority": "NCIt",
            }
        )
    return rows


def iter_oncotree_nodes(payload: Any) -> Iterator[dict[str, Any]]:
    if isinstance(payload, list):
        for item in payload:
            yield from iter_oncotree_nodes(item)
        return
    if not isinstance(payload, dict):
        return
    if "code" in payload:
        yield payload
    children = payload.get("children")
    if isinstance(children, dict):
        for child in children.values():
            yield from iter_oncotree_nodes(child)
    elif isinstance(children, list):
        for child in children:
            yield from iter_oncotree_nodes(child)


def parse_oncotree_tumor_types(path: Path) -> list[dict[str, Any]]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return []
    rows = []
    for node in iter_oncotree_nodes(payload):
        code = str(node.get("code") or "").strip()
        if not code or code == "TISSUE":
            continue
        name = str(node.get("name") or "").strip()
        main_type = str(node.get("mainType") or "").strip()
        tissue = str(node.get("tissue") or "").strip()
        canonical = main_type if main_type and not is_oncology_disease_term(name, [], []) else name
        aliases = dedupe([name, main_type, tissue])
        if not is_oncology_disease_term(canonical, aliases, []):
            continue
        xrefs = [f"ONCOTREE:{code}"]
        refs = node.get("externalReferences") or {}
        if isinstance(refs, dict):
            for raw_source, values in refs.items():
                source = "NCIT" if str(raw_source).upper() == "NCI" else normalize_source_name(str(raw_source))
                for value in values or []:
                    parsed = normalize_xref(f"{source}:{value}")
                    if parsed:
                        xrefs.append(parsed[2])
        parent = str(node.get("parent") or "").strip()
        parents = [f"ONCOTREE:{parent}"] if parent and parent != "TISSUE" else []
        rows.append(
            {
                "primary_external_id": f"ONCOTREE:{code}",
                "name": canonical or name or code,
                "aliases": aliases,
                "description": f"OncoTree code {code}; tissue={tissue}" if tissue else f"OncoTree code {code}",
                "parents": parents,
                "external_xrefs": dedupe(xrefs),
                "source_priority": "OncoTree",
            }
        )
    return rows


def collect_texts(node: ET.Element, path: str) -> list[str]:
    found = node.findall(path)
    values = []
    for item in found:
        if item.text:
            values.append(item.text.strip())
    return values


def build_diseases(ctx: BuildContext) -> dict[str, str]:
    started = time.time()
    disease_writer = ctx.writer("diseases")
    xref_writer = ctx.writer("disease_xrefs")
    parent_writer = ctx.writer("disease_parent_edges")
    disease_dir = ctx.source_dir("disease_ontologies")
    disease_rows_by_uid: dict[str, dict[str, Any]] = {}
    external_to_uid: dict[str, str] = {}
    xrefs_by_uid: dict[str, set[str]] = defaultdict(set)

    def add_disease(row: dict[str, Any], source_id: str) -> str:
        uid = ""
        for x in row["external_xrefs"]:
            if x in external_to_uid:
                uid = external_to_uid[x]
                break
        if not uid:
            uid = disease_uid_for_external(row["primary_external_id"])
        existing = disease_rows_by_uid.get(uid)
        if existing is None or existing.get("source_priority") != "MONDO":
            disease_rows_by_uid[uid] = {
                "disease_uid": uid,
                "primary_external_id": existing.get("primary_external_id") if existing and existing.get("source_priority") == "MONDO" else row["primary_external_id"],
                "name": row.get("name") or (existing or {}).get("name", ""),
                "aliases": dedupe((existing or {}).get("aliases", []) + row.get("aliases", [])),
                "description": row.get("description") or (existing or {}).get("description", ""),
                "parents": dedupe((existing or {}).get("parents", []) + row.get("parents", [])),
                "external_xrefs": [],
                "source_priority": row.get("source_priority", source_id),
                "source_release": ctx.source_release(source_id),
                "license_id": ctx.license_id(source_id),
                "checksum": ctx.checksum(source_id),
                "parser_hash": ctx.parser_hash,
            }
        else:
            existing["aliases"] = dedupe(existing.get("aliases", []) + row.get("aliases", []))
            existing["parents"] = dedupe(existing.get("parents", []) + row.get("parents", []))
        for x in row["external_xrefs"]:
            external_to_uid[x] = uid
            xrefs_by_uid[uid].add(x)
        return uid

    mondo_path = disease_dir / "mondo" / "mondo.obo"
    if mondo_path.exists():
        for row in parse_obo_disease_terms(mondo_path, "MONDO", ctx):
            add_disease(row, "disease_ontologies")

    doid_path = disease_dir / "doid" / "doid.obo"
    if doid_path.exists():
        for row in parse_obo_disease_terms(doid_path, "DOID", ctx):
            add_disease(row, "disease_ontologies")

    mesh_needed = {x.split(":", 1)[1] for x in external_to_uid if x.startswith("MESH:")}
    desc_path = disease_dir / "desc2026.gz"
    if desc_path.exists():
        with gzip.open(desc_path, "rb") as handle:
            for _, elem in ET.iterparse(handle, events=("end",)):
                if elem.tag != "DescriptorRecord":
                    continue
                mesh_id = (elem.findtext("DescriptorUI") or "").strip()
                name = (elem.findtext("DescriptorName/String") or "").strip()
                tree_numbers = collect_texts(elem, "TreeNumberList/TreeNumber")
                is_disease_tree = any(t.startswith("C") for t in tree_numbers)
                if mesh_id and name and (mesh_id in mesh_needed or is_disease_tree):
                    key = f"MESH:{mesh_id}"
                    uid = external_to_uid.get(key) or disease_uid_for_external(key)
                    aliases = collect_texts(elem, "ConceptList/Concept/TermList/Term/String")
                    row = {
                        "primary_external_id": key,
                        "name": name,
                        "aliases": dedupe(aliases),
                        "description": "",
                        "parents": [],
                        "external_xrefs": [key],
                        "source_priority": "MeSH",
                    }
                    add_disease(row, "disease_ontologies")
                    external_to_uid[key] = uid
                elem.clear()

    ncit_path = ctx.source_dir("ncit") / "ncit" / "ncit.obo"
    if ncit_path.exists():
        ncit_rows = parse_ncit_cancer_terms(ncit_path, ctx)
        for row in ncit_rows:
            add_disease(row, "ncit")
        ctx.note("info", "ncit", f"Parsed {len(ncit_rows)} oncology-focused NCIt terms into the disease vocabulary.")

    oncotree_path = ctx.source_dir("oncotree") / "tumorTypes.latest_stable.json"
    if oncotree_path.exists():
        oncotree_rows = parse_oncotree_tumor_types(oncotree_path)
        for row in oncotree_rows:
            add_disease(row, "oncotree")
        ctx.note("info", "oncotree", f"Parsed {len(oncotree_rows)} OncoTree tumor type terms into the disease vocabulary.")

    ot_disease_path = ctx.source_dir("opentargets_core") / "disease" / "disease.parquet"
    if pa is not None and ot_disease_path.exists():
        table = pq.read_table(ot_disease_path, columns=["id", "name", "description", "dbXRefs", "parents", "synonyms"])
        for row in table.to_pylist():
            ot_id = row.get("id", "")
            parsed = normalize_prefixed_id(ot_id)
            if not parsed:
                continue
            xrefs = {parsed[2]}
            for raw in row.get("dbXRefs") or []:
                parsed_x = normalize_xref(raw)
                if parsed_x:
                    xrefs.add(parsed_x[2])
            synonyms_obj = row.get("synonyms") or {}
            aliases: list[str] = []
            if isinstance(synonyms_obj, dict):
                for values in synonyms_obj.values():
                    aliases.extend(values or [])
            parents = []
            for parent in row.get("parents") or []:
                parsed_parent = normalize_prefixed_id(parent)
                if parsed_parent:
                    parents.append(parsed_parent[2])
            add_disease(
                {
                    "primary_external_id": parsed[2],
                    "name": row.get("name", ""),
                    "aliases": dedupe(aliases),
                    "description": row.get("description", ""),
                    "parents": dedupe(parents),
                    "external_xrefs": sorted(xrefs),
                    "source_priority": "Open Targets",
                },
                "opentargets_core",
            )

    disease_rows = []
    xref_rows = []
    parent_rows = []
    for uid, row in sorted(disease_rows_by_uid.items()):
        row["external_xrefs"] = sorted(xrefs_by_uid.get(uid, set()))
        disease_rows.append(row)
        for x in row["external_xrefs"]:
            parsed = normalize_xref(x)
            if not parsed:
                continue
            xref_rows.append(
                {
                    "xref_uid": stable_uid("xref", uid, parsed[2]),
                    "disease_uid": uid,
                    "xref_source": parsed[0],
                    "xref_id": parsed[1],
                    "xref_key": parsed[2],
                    "source_name": row["source_priority"],
                    "source_release": row["source_release"],
                    "license_id": row["license_id"],
                    "parser_hash": ctx.parser_hash,
                }
            )
        for parent_key in row.get("parents", []):
            parent_uid = external_to_uid.get(parent_key)
            if not parent_uid:
                continue
            parent_rows.append(
                {
                    "edge_uid": stable_uid("edge", "disease_parent", uid, parent_uid),
                    "child_disease_uid": uid,
                    "parent_disease_uid": parent_uid,
                    "child_external_id": row["primary_external_id"],
                    "parent_external_id": parent_key,
                    "predicate": "is_a",
                    "source_name": row["source_priority"],
                    "source_release": row["source_release"],
                    "license_id": row["license_id"],
                    "parser_hash": ctx.parser_hash,
                }
            )
        if len(disease_rows) >= 20_000:
            disease_writer.write(disease_rows)
            disease_rows = []
        if len(xref_rows) >= 50_000:
            xref_writer.write(xref_rows)
            xref_rows = []
        if len(parent_rows) >= 50_000:
            parent_writer.write(parent_rows)
            parent_rows = []
    disease_writer.write(disease_rows)
    xref_writer.write(xref_rows)
    parent_writer.write(parent_rows)
    ctx.close_writer(disease_writer, started)
    ctx.close_writer(xref_writer, started)
    ctx.close_writer(parent_writer, started)
    return external_to_uid


CELL_STATE_SEEDS = [
    {
        "primary_external_id": "METABO_STATE:malignant",
        "name": "malignant",
        "aliases": ["malignant state", "tumor cell state", "cancer cell state"],
        "description": "Cell state indicating malignant or neoplastic behavior in a cancer context.",
        "parents": [],
        "state_category": "tumor_status",
    },
    {
        "primary_external_id": "METABO_STATE:proliferating",
        "name": "proliferating",
        "aliases": ["cycling", "cell cycle active", "mitotic"],
        "description": "Cell state indicating active proliferation or cell-cycle activity.",
        "parents": [],
        "state_category": "proliferation",
    },
    {
        "primary_external_id": "METABO_STATE:hypoxic",
        "name": "hypoxic",
        "aliases": ["hypoxia-associated", "low oxygen"],
        "description": "Cell state associated with hypoxia or low-oxygen response.",
        "parents": [],
        "state_category": "microenvironment",
    },
    {
        "primary_external_id": "METABO_STATE:emt_like",
        "name": "EMT-like",
        "aliases": ["epithelial mesenchymal transition-like", "mesenchymal-like"],
        "description": "Cell state resembling epithelial-to-mesenchymal transition programs.",
        "parents": [],
        "state_category": "transition",
    },
    {
        "primary_external_id": "METABO_STATE:stem_like",
        "name": "stem-like",
        "aliases": ["stemness", "cancer stem-like"],
        "description": "Cell state associated with stemness or progenitor-like programs.",
        "parents": [],
        "state_category": "differentiation",
    },
    {
        "primary_external_id": "METABO_STATE:glycolytic",
        "name": "glycolytic",
        "aliases": ["glycolysis-high", "glycolytic state"],
        "description": "Cell state associated with elevated glycolysis or lactate-linked metabolic programs.",
        "parents": [],
        "state_category": "metabolic",
    },
]


def ontology_uid_for_external(uid_prefix: str, external_key: str) -> str:
    return stable_uid(uid_prefix, external_key)


def parse_obo_axis_terms(
    path: Path,
    allowed_primary_sources: set[str],
    source_priority: str,
    required_xref_sources: set[str] | None = None,
) -> list[dict[str, Any]]:
    rows = []
    for term in iter_obo_terms(path):
        ids = term.get("id", [])
        names = term.get("name", [])
        if not ids or not names or term.get("is_obsolete", ["false"])[0] == "true":
            continue
        parsed_primary = normalize_xref(ids[0])
        if not parsed_primary or parsed_primary[0] not in allowed_primary_sources:
            continue
        xrefs = {parsed_primary[2]}
        for alt_id in term.get("alt_id", []):
            parsed = normalize_xref(alt_id)
            if parsed:
                xrefs.add(parsed[2])
        for raw_xref in term.get("xref", []):
            parsed = normalize_xref(raw_xref)
            if parsed:
                xrefs.add(parsed[2])
        if required_xref_sources and not any(x.split(":", 1)[0] in required_xref_sources for x in xrefs if ":" in x):
            continue
        parents = []
        for parent in term.get("is_a", []):
            parsed = normalize_xref(parent)
            if parsed:
                parents.append(parsed[2])
        aliases = [quoted_value(raw) for raw in term.get("synonym", [])]
        rows.append(
            {
                "primary_external_id": parsed_primary[2],
                "name": names[0],
                "aliases": dedupe(aliases),
                "description": quoted_value(term.get("def", [""])[0]) if term.get("def") else "",
                "parents": dedupe(parents),
                "external_xrefs": sorted(xrefs),
                "source_priority": source_priority,
            }
        )
    return rows


def write_context_axis_tables(
    ctx: BuildContext,
    table_name: str,
    xref_table_name: str,
    parent_table_name: str,
    uid_prefix: str,
    uid_col: str,
    child_col: str,
    parent_col: str,
    rows_by_uid: dict[str, dict[str, Any]],
    external_to_uid: dict[str, str],
    xrefs_by_uid: dict[str, set[str]],
    started: float,
) -> dict[str, str]:
    entity_writer = ctx.writer(table_name)
    xref_writer = ctx.writer(xref_table_name)
    parent_writer = ctx.writer(parent_table_name)
    entity_rows = []
    xref_rows = []
    parent_rows = []
    for uid, row in sorted(rows_by_uid.items()):
        row["external_xrefs"] = sorted(xrefs_by_uid.get(uid, set()))
        entity_rows.append(row)
        for x in row["external_xrefs"]:
            parsed = normalize_xref(x)
            if not parsed:
                continue
            xref_rows.append(
                {
                    "xref_uid": stable_uid("xref", uid, parsed[2]),
                    uid_col: uid,
                    "xref_source": parsed[0],
                    "xref_id": parsed[1],
                    "xref_key": parsed[2],
                    "source_name": row["source_priority"],
                    "source_release": row["source_release"],
                    "license_id": row["license_id"],
                    "parser_hash": ctx.parser_hash,
                }
            )
        for parent_key in row.get("parents", []):
            parent_uid = external_to_uid.get(parent_key)
            if not parent_uid:
                continue
            parent_rows.append(
                {
                    "edge_uid": stable_uid("edge", parent_table_name, uid, parent_uid),
                    child_col: uid,
                    parent_col: parent_uid,
                    "child_external_id": row["primary_external_id"],
                    "parent_external_id": parent_key,
                    "predicate": "is_a",
                    "source_name": row["source_priority"],
                    "source_release": row["source_release"],
                    "license_id": row["license_id"],
                    "parser_hash": ctx.parser_hash,
                }
            )
        if len(entity_rows) >= 20_000:
            entity_writer.write(entity_rows)
            entity_rows = []
        if len(xref_rows) >= 50_000:
            xref_writer.write(xref_rows)
            xref_rows = []
        if len(parent_rows) >= 50_000:
            parent_writer.write(parent_rows)
            parent_rows = []
    entity_writer.write(entity_rows)
    xref_writer.write(xref_rows)
    parent_writer.write(parent_rows)
    ctx.close_writer(entity_writer, started)
    ctx.close_writer(xref_writer, started)
    ctx.close_writer(parent_writer, started)
    return external_to_uid


def build_context_axis(
    ctx: BuildContext,
    table_name: str,
    xref_table_name: str,
    parent_table_name: str,
    uid_prefix: str,
    uid_col: str,
    child_col: str,
    parent_col: str,
    source_rows: Iterable[tuple[dict[str, Any], str]],
) -> dict[str, str]:
    started = time.time()
    rows_by_uid: dict[str, dict[str, Any]] = {}
    external_to_uid: dict[str, str] = {}
    xrefs_by_uid: dict[str, set[str]] = defaultdict(set)

    def add_row(row: dict[str, Any], source_id: str) -> str:
        uid = ""
        for xref in row["external_xrefs"]:
            if xref in external_to_uid:
                uid = external_to_uid[xref]
                break
        if not uid:
            uid = ontology_uid_for_external(uid_prefix, row["primary_external_id"])
        existing = rows_by_uid.get(uid)
        if existing is None:
            payload = {
                uid_col: uid,
                "primary_external_id": row["primary_external_id"],
                "name": row.get("name", ""),
                "aliases": dedupe(row.get("aliases", [])),
                "description": row.get("description", ""),
                "parents": dedupe(row.get("parents", [])),
                "external_xrefs": [],
                "source_priority": row.get("source_priority", source_id),
                "source_release": ctx.source_release(source_id),
                "license_id": ctx.license_id(source_id),
                "checksum": ctx.checksum(source_id),
                "parser_hash": ctx.parser_hash,
            }
            if table_name == "cell_states":
                payload["state_category"] = row.get("state_category", "")
            rows_by_uid[uid] = payload
        else:
            existing["aliases"] = dedupe(existing.get("aliases", []) + row.get("aliases", []))
            existing["parents"] = dedupe(existing.get("parents", []) + row.get("parents", []))
            if not existing.get("description"):
                existing["description"] = row.get("description", "")
            if table_name == "cell_states" and not existing.get("state_category"):
                existing["state_category"] = row.get("state_category", "")
        for xref in row["external_xrefs"]:
            external_to_uid[xref] = uid
            xrefs_by_uid[uid].add(xref)
        return uid

    added = 0
    for row, source_id in source_rows:
        add_row(row, source_id)
        added += 1
    ctx.note("info", table_name, f"Parsed {added} source rows into {table_name}.")
    return write_context_axis_tables(
        ctx,
        table_name,
        xref_table_name,
        parent_table_name,
        uid_prefix,
        uid_col,
        child_col,
        parent_col,
        rows_by_uid,
        external_to_uid,
        xrefs_by_uid,
        started,
    )


def cell_type_source_rows(ctx: BuildContext) -> Iterator[tuple[dict[str, Any], str]]:
    source_dir = ctx.source_dir("cell_ontology")
    for filename in ("cl.obo", "cl-basic.obo"):
        path = source_dir / filename
        if path.exists():
            for row in parse_obo_axis_terms(path, {"CL"}, "Cell Ontology"):
                yield row, "cell_ontology"
            return


def cell_state_source_rows() -> Iterator[tuple[dict[str, Any], str]]:
    for row in CELL_STATE_SEEDS:
        payload = {
            **row,
            "external_xrefs": [row["primary_external_id"]],
            "source_priority": "local controlled vocabulary",
        }
        yield payload, "manual_context_axes"


def tissue_source_rows(ctx: BuildContext) -> Iterator[tuple[dict[str, Any], str]]:
    uberon_dir = ctx.source_dir("uberon")
    for filename in ("uberon.obo", "basic.obo", "uberon-basic.obo"):
        path = uberon_dir / filename
        if path.exists():
            for row in parse_obo_axis_terms(path, {"UBERON"}, "UBERON"):
                yield row, "uberon"
            break
    efo_path = ctx.source_dir("efo") / "efo.obo"
    if efo_path.exists():
        for row in parse_obo_axis_terms(efo_path, {"EFO"}, "EFO", required_xref_sources={"UBERON"}):
            yield row, "efo"


def build_context_axes(ctx: BuildContext) -> None:
    build_context_axis(
        ctx,
        "cell_types",
        "cell_type_xrefs",
        "cell_type_parent_edges",
        "cell_type",
        "cell_type_uid",
        "child_cell_type_uid",
        "parent_cell_type_uid",
        cell_type_source_rows(ctx),
    )
    build_context_axis(
        ctx,
        "cell_states",
        "cell_state_xrefs",
        "cell_state_parent_edges",
        "cell_state",
        "cell_state_uid",
        "child_cell_state_uid",
        "parent_cell_state_uid",
        cell_state_source_rows(),
    )
    build_context_axis(
        ctx,
        "tissues",
        "tissue_xrefs",
        "tissue_parent_edges",
        "tissue",
        "tissue_uid",
        "child_tissue_uid",
        "parent_tissue_uid",
        tissue_source_rows(ctx),
    )


def parquet_files(path: Path) -> list[Path]:
    return sorted(p for p in path.glob("*.parquet") if p.is_file())


def iter_parquet_rows(
    files: list[Path],
    columns: list[str],
    batch_size: int = 10_000,
    ctx: BuildContext | None = None,
    component: str = "parquet",
) -> Iterator[dict[str, Any]]:
    if pa is None:
        return
    for path in files:
        try:
            parquet_file = pq.ParquetFile(path)
            for batch in parquet_file.iter_batches(columns=columns, batch_size=batch_size):
                for row in batch.to_pylist():
                    yield row
        except Exception as exc:
            if ctx is not None:
                ctx.note("warning", component, f"Skipped parquet partition {path.name}: {exc}")
            continue


def build_opentargets(
    ctx: BuildContext,
    gene_uid_by_ensembl: dict[str, str],
    disease_uid_by_xref: dict[str, str],
) -> dict[str, str]:
    started = time.time()
    target_writer = ctx.writer("targets")
    xref_writer = ctx.writer("target_xrefs")
    target_gene_writer = ctx.writer("target_gene_edges")
    target_disease_writer = ctx.writer("target_disease_edges")
    uid_by_target: dict[str, str] = {}
    ot_dir = ctx.source_dir("opentargets_core")
    if pa is None:
        ctx.note("warning", "opentargets_core", "pyarrow is not installed; Open Targets parquet parsing was skipped.")
        for writer in [target_writer, xref_writer, target_gene_writer, target_disease_writer]:
            ctx.close_writer(writer, started)
        return uid_by_target

    target_rows = []
    xref_rows = []
    target_gene_rows = []
    target_files = parquet_files(ot_dir / "target")
    for row in iter_parquet_rows(
        target_files,
        ["id", "approvedSymbol", "approvedName", "biotype", "dbXrefs", "proteinIds", "tractability"],
        ctx=ctx,
        component="opentargets_core.target",
    ):
        target_id = row.get("id", "")
        if not target_id:
            continue
        target_uid = stable_uid("target", "OPENTARGETS", target_id)
        uid_by_target[target_id] = target_uid
        gene_uid = gene_uid_by_ensembl.get(target_id, "")
        xrefs = {f"OPENTARGETS:{target_id}"}
        if target_id.startswith("ENSG"):
            xrefs.add(f"ENSEMBL:{target_id}")
        for raw in row.get("dbXrefs") or []:
            if isinstance(raw, dict):
                parsed = normalize_xref(f"{raw.get('source', '')}:{raw.get('id', '')}")
            else:
                parsed = normalize_xref(str(raw))
            if parsed:
                xrefs.add(parsed[2])
        for protein in row.get("proteinIds") or []:
            if isinstance(protein, dict) and protein.get("id"):
                source = normalize_source_name(protein.get("source", "UNIPROT"))
                xrefs.add(f"{source}:{protein['id']}")
        tractability = {}
        for item in row.get("tractability") or []:
            if isinstance(item, dict):
                key = f"{item.get('modality')}:{item.get('id')}"
                tractability[key] = item.get("value")
        target_rows.append(
            {
                "target_uid": target_uid,
                "target_external_id": target_id,
                "preferred_name": row.get("approvedName", ""),
                "approved_symbol": row.get("approvedSymbol", ""),
                "target_type": row.get("biotype", ""),
                "gene_uid": gene_uid,
                "tractability_flags_json": json_dumps(tractability),
                "external_xrefs": sorted(xrefs),
                "source_release": ctx.source_release("opentargets_core"),
                "license_id": ctx.license_id("opentargets_core"),
                "checksum": ctx.checksum("opentargets_core"),
                "parser_hash": ctx.parser_hash,
            }
        )
        for x in sorted(xrefs):
            parsed = normalize_xref(x)
            if not parsed:
                continue
            xref_rows.append(
                {
                    "xref_uid": stable_uid("xref", target_uid, parsed[2]),
                    "target_uid": target_uid,
                    "xref_source": parsed[0],
                    "xref_id": parsed[1],
                    "xref_key": parsed[2],
                    "source_name": "Open Targets",
                    "source_release": ctx.source_release("opentargets_core"),
                    "license_id": ctx.license_id("opentargets_core"),
                    "parser_hash": ctx.parser_hash,
                }
            )
        if gene_uid:
            target_gene_rows.append(
                {
                    "edge_uid": stable_uid("edge", "target_gene", target_uid, gene_uid),
                    "target_uid": target_uid,
                    "gene_uid": gene_uid,
                    "target_external_id": target_id,
                    "gene_external_id": f"ENSEMBL:{target_id}",
                    "predicate": "targets",
                    "source_name": "Open Targets",
                    "source_release": ctx.source_release("opentargets_core"),
                    "license_id": ctx.license_id("opentargets_core"),
                    "parser_hash": ctx.parser_hash,
                }
            )
        if len(target_rows) >= 10_000:
            target_writer.write(target_rows)
            target_rows = []
        if len(xref_rows) >= 50_000:
            xref_writer.write(xref_rows)
            xref_rows = []
        if len(target_gene_rows) >= 10_000:
            target_gene_writer.write(target_gene_rows)
            target_gene_rows = []
    target_writer.write(target_rows)
    xref_writer.write(xref_rows)
    target_gene_writer.write(target_gene_rows)

    edge_rows = []
    assoc_files = parquet_files(ot_dir / "association_overall_direct")
    for row in iter_parquet_rows(
        assoc_files,
        ["diseaseId", "targetId", "aggregationType", "aggregationValue", "associationScore", "evidenceCount", "currentNovelty"],
        batch_size=50_000,
        ctx=ctx,
        component="opentargets_core.association_overall_direct",
    ):
        target_id = row.get("targetId", "")
        disease_id = row.get("diseaseId", "")
        target_uid = uid_by_target.get(target_id) or stable_uid("target", "OPENTARGETS", target_id)
        parsed_disease = normalize_prefixed_id(disease_id)
        disease_key = parsed_disease[2] if parsed_disease else disease_id
        disease_uid = disease_uid_by_xref.get(disease_key) or disease_uid_for_external(disease_key)
        score = maybe_float(row.get("associationScore"))
        edge_rows.append(
            {
                "edge_uid": stable_uid("edge", "target_disease", target_id, disease_key, row.get("aggregationType", ""), row.get("aggregationValue", "")),
                "target_uid": target_uid,
                "disease_uid": disease_uid,
                "target_external_id": target_id,
                "disease_external_id": disease_key,
                "predicate": "associated_with",
                "association_score": score,
                "evidence_count": maybe_int(row.get("evidenceCount")),
                "aggregation_type": row.get("aggregationType", ""),
                "aggregation_value": row.get("aggregationValue", ""),
                "current_novelty": maybe_float(row.get("currentNovelty")),
                "score_components_json": json_dumps({"associationScore": score, "evidenceCount": row.get("evidenceCount")}),
                "source_name": "Open Targets",
                "source_record_id": f"{target_id}|{disease_key}|{row.get('aggregationType', '')}|{row.get('aggregationValue', '')}",
                "evidence_level": "curated/inferred",
                "source_release": ctx.source_release("opentargets_core"),
                "license_id": ctx.license_id("opentargets_core"),
                "parser_hash": ctx.parser_hash,
            }
        )
        if len(edge_rows) >= 100_000:
            target_disease_writer.write(edge_rows)
            edge_rows = []
    target_disease_writer.write(edge_rows)
    ctx.close_writer(target_writer, started)
    ctx.close_writer(xref_writer, started)
    ctx.close_writer(target_gene_writer, started)
    ctx.close_writer(target_disease_writer, started)
    return uid_by_target


def parse_article_metadata(text: str, source_path: Path) -> dict[str, str]:
    def line_value(label: str) -> str:
        match = re.search(rf"^{re.escape(label)}:\s*(.+)$", text, re.M)
        return match.group(1).strip() if match else ""

    pmcid = line_value("PMCID") or source_path.stem.split(".", 1)[0]
    pmid = line_value("PMID")
    doi = line_value("DOI")
    journal = line_value("NLM Title Abbreviation") or line_value("Journal ID")
    article_type = line_value("Subjects")
    pub_date = line_value("Publication date") or line_value("Electronic publication date")

    lines = text.splitlines()
    title = ""
    for idx, line in enumerate(lines):
        if line.startswith("Subjects:"):
            for candidate in lines[idx + 1 : idx + 12]:
                candidate = candidate.strip()
                if candidate:
                    title = candidate
                    break
            break
    if not title:
        title = source_path.stem

    intro = INTRO_RE.search(text)
    abstract = ""
    if intro:
        before_intro = text[: intro.start()]
        marker = max(before_intro.rfind("=============================="), before_intro.rfind("================"))
        if marker >= 0:
            abstract = before_intro[marker:].split("\n", 1)[-1].strip()
        else:
            abstract = before_intro[-4000:].strip()
    abstract = re.sub(r"\s+", " ", abstract).strip()
    if abstract and title and abstract.startswith(title):
        abstract = abstract[len(title) :].strip()
    return {
        "pmcid": pmcid,
        "pmid": pmid,
        "doi": doi,
        "journal": journal,
        "article_type": article_type,
        "pub_date": pub_date,
        "title": title,
        "abstract": abstract,
    }


def load_scispacy_pipeline(model_name: str) -> Any:
    try:
        import spacy  # type: ignore
    except ImportError as exc:  # pragma: no cover - depends on optional env.
        raise RuntimeError("scispaCy parsing requires optional dependency 'spacy' and the scispaCy model package.") from exc
    try:
        nlp = spacy.load(model_name)
    except Exception as exc:  # pragma: no cover - depends on optional env.
        raise RuntimeError(f"Unable to load scispaCy model '{model_name}'. Install scispaCy and the model before using --sentence-parser scispacy.") from exc
    if "parser" not in nlp.pipe_names and "senter" not in nlp.pipe_names and "sentencizer" not in nlp.pipe_names:
        nlp.add_pipe("sentencizer")
    return nlp


def scispacy_sentence_spans(content: str, nlp: Any) -> list[tuple[int, int, str]]:
    doc = nlp(content)
    spans: list[tuple[int, int, str]] = []
    for sent in doc.sents:
        sentence = re.sub(r"\s+", " ", sent.text).strip()
        if len(sentence) < 8:
            continue
        spans.append((int(sent.start_char), int(sent.end_char), sentence))
    return spans


def rule_sentence_spans(content: str) -> list[tuple[int, int, str]]:
    spans: list[tuple[int, int, str]] = []
    for match in SENTENCE_RE.finditer(content):
        sentence = re.sub(r"\s+", " ", match.group(0)).strip()
        if len(sentence) < 8:
            continue
        spans.append((match.start(), match.end(), sentence))
    return spans


def section_spans(meta: dict[str, str], text: str, scope: str) -> list[tuple[str, str]]:
    sections: list[tuple[str, str]] = []
    if meta.get("title"):
        sections.append(("title", meta["title"]))
    if scope in {"abstract", "full"} and meta.get("abstract"):
        sections.append(("abstract", meta["abstract"]))
    if scope == "full":
        intro = INTRO_RE.search(text)
        body = text[intro.start() :] if intro else text
        sections.append(("full_text", body))
    return sections


def sentence_rows_from_spans(
    meta: dict[str, str],
    text: str,
    article_uid: str,
    release_id: str,
    phash: str,
    section: str,
    content: str,
    spans: list[tuple[int, int, str]],
    segmenter: str,
    parser_model: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    cursor = 0
    for start, _end, sentence in spans:
        offset = text.find(sentence[:80], cursor)
        if offset < 0:
            content_offset = text.find(content[:80])
            offset = (content_offset if content_offset >= 0 else 0) + start
        end_offset = offset + len(sentence)
        cursor = max(cursor, end_offset)
        text_hash = hashlib.sha256(sentence.encode("utf-8")).hexdigest()
        rows.append(
            {
                "sentence_uid": stable_uid("sent", article_uid, section, offset, text_hash[:16]),
                "article_uid": article_uid,
                "pmid": meta.get("pmid", ""),
                "pmcid": meta.get("pmcid", ""),
                "section": section,
                "sentence_text": sentence,
                "text_hash": text_hash,
                "start_offset": offset,
                "end_offset": end_offset,
                "language": "en",
                "segmenter": segmenter,
                "parser_model": parser_model,
                "source_release": release_id,
                "parser_hash": phash,
            }
        )
    return rows


def scispacy_candidate_rows(
    meta: dict[str, str],
    article_uid: str,
    release_id: str,
    phash: str,
    section: str,
    content: str,
    sentence_rows: list[dict[str, Any]],
    nlp: Any,
    parser_model: str,
) -> list[dict[str, Any]]:
    if not sentence_rows:
        return []
    doc = nlp(content)
    candidates: list[dict[str, Any]] = []
    used_offsets: set[tuple[str, int, int, str]] = set()
    for ent in doc.ents:
        mention = re.sub(r"\s+", " ", ent.text).strip()
        if len(mention) < 2:
            continue
        sent_text = re.sub(r"\s+", " ", getattr(ent, "sent", ent).text).strip()
        matched = sentence_entity_offset(sentence_rows, mention, sent_text)
        if matched is None:
            continue
        matched_row, start_offset, end_offset = matched
        sentence_uid = str(matched_row["sentence_uid"])
        dedupe_key = (sentence_uid, start_offset, end_offset, str(ent.label_ or ""))
        if dedupe_key in used_offsets:
            continue
        used_offsets.add(dedupe_key)
        candidates.append(
            {
                "candidate_uid": stable_uid("entcand", sentence_uid, mention, start_offset, end_offset, ent.label_),
                "sentence_uid": sentence_uid,
                "article_uid": article_uid,
                "pmid": meta.get("pmid", ""),
                "pmcid": meta.get("pmcid", ""),
                "section": section,
                "mention_text": mention,
                "normalized_surface": re.sub(r"\s+", " ", mention).strip().casefold(),
                "start_offset": start_offset,
                "end_offset": end_offset,
                "scispacy_label": str(ent.label_ or ""),
                "candidate_source": "scispacy",
                "parser_model": parser_model,
                "confidence_proxy": 0.5,
                "source_release": release_id,
                "license_id": "local_articles:articles_collect",
                "parser_hash": phash,
            }
        )
    return candidates


def sentence_entity_offset(sentence_rows: list[dict[str, Any]], mention: str, preferred_sentence: str) -> tuple[dict[str, Any], int, int] | None:
    preferred_norm = re.sub(r"\s+", " ", preferred_sentence).strip()
    for row in sentence_rows:
        sent_start = int(row["start_offset"])
        sent_end = int(row["end_offset"])
        sentence_text = str(row.get("sentence_text") or "")
        if preferred_norm and sentence_text != preferred_norm:
            continue
        match = re.search(re.escape(mention), sentence_text)
        if match:
            return row, sent_start + match.start(), min(sent_start + match.end(), sent_end)
    for row in sentence_rows:
        sent_start = int(row["start_offset"])
        sent_end = int(row["end_offset"])
        sentence_text = str(row.get("sentence_text") or "")
        match = re.search(re.escape(mention), sentence_text)
        if match:
            return row, sent_start + match.start(), min(sent_start + match.end(), sent_end)
    return None


def sentence_and_candidate_rows_for_article(
    meta: dict[str, str],
    text: str,
    article_uid: str,
    scope: str,
    release_id: str,
    phash: str,
    sentence_parser: str,
    scispacy_model: str,
    nlp: Any | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], int]:
    sentence_rows: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    fallback_count = 0
    for section, content in section_spans(meta, text, scope):
        segmenter = "rules"
        parser_model = "regex"
        spans = rule_sentence_spans(content)
        use_scispacy = sentence_parser in {"scispacy", "hybrid"} and nlp is not None
        if use_scispacy:
            try:
                scispacy_spans = scispacy_sentence_spans(content, nlp)
            except Exception:
                if sentence_parser == "scispacy":
                    raise
                fallback_count += 1
                scispacy_spans = []
            if scispacy_spans:
                spans = scispacy_spans
                segmenter = "scispacy"
                parser_model = scispacy_model
            elif sentence_parser == "scispacy":
                spans = []
                segmenter = "scispacy"
                parser_model = scispacy_model
            elif sentence_parser == "hybrid":
                fallback_count += 1
        section_sentence_rows = sentence_rows_from_spans(meta, text, article_uid, release_id, phash, section, content, spans, segmenter, parser_model)
        sentence_rows.extend(section_sentence_rows)
        if segmenter == "scispacy" and nlp is not None:
            candidate_rows.extend(scispacy_candidate_rows(meta, article_uid, release_id, phash, section, content, section_sentence_rows, nlp, parser_model))
    return sentence_rows, candidate_rows, fallback_count


def sentence_rows_for_article(meta: dict[str, str], text: str, article_uid: str, scope: str, release_id: str, phash: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for section, content in section_spans(meta, text, scope):
        rows.extend(sentence_rows_from_spans(meta, text, article_uid, release_id, phash, section, content, rule_sentence_spans(content), "rules", "regex"))
    return rows


def build_articles(
    ctx: BuildContext,
    articles_dir: Path,
    max_articles: int,
    sentence_scope: str,
    sentence_parser: str = "rules",
    scispacy_model: str = DEFAULT_SCISPACY_MODEL,
) -> None:
    started = time.time()
    article_writer = ctx.writer("articles")
    sentence_writer = ctx.writer("sentences")
    candidate_writer = ctx.writer("sentence_entity_candidates")
    files = sorted(articles_dir.glob("*.txt"))
    if max_articles > 0:
        files = files[:max_articles]
    article_rows = []
    sentence_rows = []
    candidate_rows = []
    scispacy_fallback_count = 0
    nlp = None
    if sentence_parser in {"scispacy", "hybrid"}:
        try:
            nlp = load_scispacy_pipeline(scispacy_model)
        except RuntimeError as exc:
            if sentence_parser == "scispacy":
                raise
            scispacy_fallback_count += 1
            ctx.note("warning", "articles_collect", f"scispaCy unavailable for hybrid sentence parsing; falling back to rules. {exc}")
    for index, path in enumerate(files, start=1):
        raw = path.read_bytes()
        checksum = hashlib.sha256(raw).hexdigest()
        text = raw.decode("utf-8", errors="replace")
        meta = parse_article_metadata(text, path)
        anchor = meta.get("pmid") or meta.get("pmcid") or path.stem
        article_uid = stable_uid("article", anchor)
        article_rows.append(
            {
                "article_uid": article_uid,
                "pmid": meta.get("pmid", ""),
                "pmcid": meta.get("pmcid", ""),
                "doi": meta.get("doi", ""),
                "title": meta.get("title", ""),
                "abstract": meta.get("abstract", ""),
                "journal": meta.get("journal", ""),
                "pub_date": meta.get("pub_date", ""),
                "article_type": meta.get("article_type", ""),
                "source_path": str(path),
                "source_bytes": len(raw),
                "checksum": checksum,
                "source_release": ctx.release_id,
                "license_id": "local_articles:articles_collect",
                "parser_hash": ctx.parser_hash,
            }
        )
        if sentence_parser == "rules":
            article_sentence_rows = sentence_rows_for_article(meta, text, article_uid, sentence_scope, ctx.release_id, ctx.parser_hash)
            article_candidate_rows: list[dict[str, Any]] = []
            article_fallbacks = 0
        else:
            article_sentence_rows, article_candidate_rows, article_fallbacks = sentence_and_candidate_rows_for_article(
                meta,
                text,
                article_uid,
                sentence_scope,
                ctx.release_id,
                ctx.parser_hash,
                sentence_parser,
                scispacy_model,
                nlp,
            )
        scispacy_fallback_count += article_fallbacks
        sentence_rows.extend(article_sentence_rows)
        candidate_rows.extend(article_candidate_rows)
        if len(article_rows) >= 1000:
            article_writer.write(article_rows)
            article_rows = []
        if len(sentence_rows) >= 25_000:
            sentence_writer.write(sentence_rows)
            sentence_rows = []
        if len(candidate_rows) >= 25_000:
            candidate_writer.write(candidate_rows)
            candidate_rows = []
        if index % 5000 == 0:
            print(f"[articles] parsed {index}/{len(files)}", flush=True)
    article_writer.write(article_rows)
    sentence_writer.write(sentence_rows)
    candidate_writer.write(candidate_rows)
    ctx.close_writer(article_writer, started)
    ctx.close_writer(sentence_writer, started)
    ctx.close_writer(candidate_writer, started)
    if sentence_parser in {"scispacy", "hybrid"}:
        ctx.note(
            "info",
            "articles_collect",
            f"Sentence parser mode={sentence_parser}; scispacy_model={scispacy_model}; scispacy_fallback_count={scispacy_fallback_count}. Entity candidates are candidate spans only, not canonical facts.",
        )
    if sentence_scope != "full":
        ctx.note(
            "info",
            "articles_collect",
            "Article rows include full-text path, byte size, and checksum. Sentence table was built for title/abstract scope; rerun with --article-sentence-scope full for full-text sentence materialization.",
        )


def write_notes_and_manifest(ctx: BuildContext, args: argparse.Namespace) -> None:
    started = time.time()
    notes_writer = ctx.writer("build_notes")
    notes_writer.write(ctx.notes)
    ctx.close_writer(notes_writer, started)
    manifest = {
        "release_id": ctx.release_id,
        "created_at_utc": ctx.created_at_utc,
        "parser_hash": ctx.parser_hash,
        "raw_manifest": str(ctx.manifest_dir / f"{ctx.release_id}.manifest.json"),
        "output_dir": str(ctx.output_dir),
        "args": vars(args),
        "tables": ctx.table_audit,
        "notes": ctx.notes,
    }
    (ctx.output_dir / "normalized_manifest.json").write_text(json_dumps(manifest) + "\n", encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build MVP normalized store from raw_lake.")
    parser.add_argument("--workspace", default=".", help="Project workspace root.")
    parser.add_argument("--raw-root", default="raw_lake")
    parser.add_argument("--manifest-dir", default="manifests")
    parser.add_argument("--output-root", default="normalized_store")
    parser.add_argument("--release-id", default="", help="Raw release id. Defaults to latest mvp_*.manifest.json.")
    parser.add_argument("--articles-dir", default="articles_collect")
    parser.add_argument("--max-articles", type=int, default=0, help="Limit local article parsing for smoke tests.")
    parser.add_argument(
        "--article-sentence-scope",
        choices=["none", "abstract", "full"],
        default="abstract",
        help="Materialize no sentences, title/abstract sentences, or full-text sentences.",
    )
    parser.add_argument("--skip-articles", action="store_true")
    parser.add_argument("--pubchem-enrich", action="store_true", help="Reserved for a cached PubChem CID property index.")
    parser.add_argument(
        "--sentence-parser",
        choices=["rules", "scispacy", "hybrid"],
        default="rules",
        help="Sentence parser for local article materialization. hybrid tries scispaCy and falls back to rules.",
    )
    parser.add_argument("--scispacy-model", default=DEFAULT_SCISPACY_MODEL, help="spaCy/scispaCy model name for --sentence-parser scispacy|hybrid.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    workspace = Path(args.workspace).resolve()
    raw_root = (workspace / args.raw_root).resolve()
    manifest_dir = (workspace / args.manifest_dir).resolve()
    output_root = (workspace / args.output_root).resolve()
    release_id, manifest = load_manifest(workspace, manifest_dir, args.release_id or None)
    catalog = load_catalog(workspace, manifest)
    ctx = BuildContext(
        workspace=workspace,
        raw_root=raw_root,
        manifest_dir=manifest_dir,
        output_root=output_root,
        release_id=release_id,
        manifest=manifest,
        catalog=catalog,
        parser_hash=parser_hash(),
    )
    ctx.output_dir.mkdir(parents=True, exist_ok=True)
    ctx.discover_unmanifested_local_sources()

    print(f"[mvpbuild] release={release_id} output={ctx.output_dir}", flush=True)
    build_raw_registry(ctx)
    metabolite_uid_by_xref = build_metabolites(ctx, include_pubchem_enrichment=args.pubchem_enrich)
    gene_uid_by_ensembl, gene_uid_by_entrez = build_genes(ctx)
    build_pathways(ctx, metabolite_uid_by_xref, gene_uid_by_ensembl, gene_uid_by_entrez)
    disease_uid_by_xref = build_diseases(ctx)
    build_context_axes(ctx)
    build_opentargets(ctx, gene_uid_by_ensembl, disease_uid_by_xref)
    if not args.skip_articles and args.article_sentence_scope != "none":
        build_articles(
            ctx,
            (workspace / args.articles_dir).resolve(),
            args.max_articles,
            args.article_sentence_scope,
            args.sentence_parser,
            args.scispacy_model,
        )
    elif not args.skip_articles:
        build_articles(
            ctx,
            (workspace / args.articles_dir).resolve(),
            args.max_articles,
            "none",
            args.sentence_parser,
            args.scispacy_model,
        )
    else:
        ctx.note("info", "articles_collect", "Article parsing skipped by --skip-articles.")
    write_notes_and_manifest(ctx, args)
    print(f"[mvpbuild] done: {ctx.output_dir / 'normalized_manifest.json'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
