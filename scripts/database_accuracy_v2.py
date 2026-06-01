#!/usr/bin/env python
"""Build accuracy-oriented v2 fact stores from the current frozen release.

The v2 store is intentionally parallel to normalized_store/graph_projection.
It separates exact chemicals, classes, traits, identity decisions, relations,
and evidence assertions so downstream analysis can stop treating every matched
surface as the same kind of metabolite fact.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import tarfile
from pathlib import Path
from typing import Any, Iterable

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except Exception as exc:  # pragma: no cover
    raise SystemExit("pyarrow is required for database_accuracy_v2") from exc

try:  # XML semantic sources are optional, but lxml is available in the app runtime.
    from lxml import etree
except Exception:  # pragma: no cover
    etree = None


DEFAULT_NORMALIZED_ROOT = "normalized_store"
DEFAULT_LITERATURE_ROOT = "literature_evidence"
DEFAULT_OUTPUT_ROOT = "database_accuracy_store"

POOL_OR_CLASS_RE = re.compile(
    r"\b(?:lipid|sphingo|sphingosine|sphingomyelin|ceramide|acylcarnitine|"
    r"phosphatidyl|lysophosphatidyl|plasmalogen|sterol|steroid|bile acid|"
    r"fatty acid|gpc|gpe|gpi|gps|gpg|pool|class|family)\b",
    flags=re.IGNORECASE,
)
LIPID_SPECIES_RE = re.compile(r"\b[A-Z]{1,5}\(?\d{1,2}:\d|d\d{1,2}:\d|\d{1,2}:\d", flags=re.IGNORECASE)


def short_hash(value: Any, length: int = 20) -> str:
    return hashlib.sha256(str(value).encode("utf-8")).hexdigest()[:length]


def stable_uid(prefix: str, *parts: object) -> str:
    payload = "|".join(str(part).strip() for part in parts if part is not None and str(part).strip())
    return f"{prefix}_{short_hash(payload)}"


def json_dumps(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def normalize_accession(value: Any) -> str:
    token = re.sub(r"[^0-9A-Za-z]+", "", str(value or "")).upper()
    if token.startswith("CGST"):
        token = "GCST" + token[4:]
    return token


def normalize_name(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip())


def name_key(value: Any) -> str:
    return re.sub(r"[^0-9a-z]+", " ", normalize_name(value).casefold()).strip()


def chebi_key(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    match = re.search(r"(\d+)$", text)
    return match.group(1) if match else ""


def strip_reactome_version(value: Any) -> str:
    text = str(value or "").strip()
    return re.sub(r"\.\d+$", "", text)


def split_values(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    delimiter = "|" if "|" in text else ";"
    return [part.strip() for part in text.split(delimiter) if part.strip()]


def european_source_path(workspace: Path, filename: str) -> Path:
    for directory in ("European", "European_point"):
        path = workspace / "raw_lake" / directory / filename
        if path.exists():
            return path
    return workspace / "raw_lake" / "European" / filename


def is_ratio_text(value: Any) -> bool:
    text = normalize_name(value).casefold()
    if not text:
        return False
    return bool(re.search(r"\bratio\b", text) or " to " in text or "/" in text)


def is_class_text(value: Any) -> bool:
    return bool(POOL_OR_CLASS_RE.search(str(value or "")))


def classify_trait_type(reported_trait: str, annotation: dict[str, Any]) -> str:
    manual_scope = " ".join(split_values(annotation.get("manual_identity_scope")))
    text = f"{reported_trait} {annotation.get('mapped_names', '')} {manual_scope}"
    if is_ratio_text(text):
        return "metabolite_ratio"
    if is_class_text(text):
        return "class_trait"
    if "x-" in text.casefold() or "unknown" in text.casefold():
        return "unknown_trait"
    if reported_trait or annotation.get("mapped_names"):
        return "metabolite_level"
    return "unknown_trait"


def entity_granularity(row: dict[str, Any]) -> str:
    label = " ".join([str(row.get("canonical_name") or ""), " ".join(row.get("synonyms") or [])])
    if not (row.get("inchikey") or row.get("external_xrefs")):
        return "unknown"
    if is_class_text(label) and not LIPID_SPECIES_RE.search(label):
        return "unknown"
    if LIPID_SPECIES_RE.search(label):
        return "lipid_species"
    return "exact_compound"


def name_risk(value: Any) -> str:
    text = normalize_name(value)
    if not text:
        return "high"
    if is_ratio_text(text):
        return "high"
    if is_class_text(text) and not LIPID_SPECIES_RE.search(text):
        return "high"
    if len(text) <= 3:
        return "medium"
    return "low"


def table_schema(name: str) -> pa.Schema:
    list_str = pa.list_(pa.string())
    schemas = {
        "source_records": pa.schema(
            [
                ("source_record_uid", pa.string()),
                ("source_name", pa.string()),
                ("source_type", pa.string()),
                ("source_version", pa.string()),
                ("download_date", pa.string()),
                ("license_id", pa.string()),
                ("file_path", pa.string()),
                ("checksum", pa.string()),
                ("parser_name", pa.string()),
                ("parser_hash", pa.string()),
                ("raw_record_id", pa.string()),
                ("raw_payload_json", pa.string()),
            ]
        ),
        "chemical_entities": pa.schema(
            [
                ("chemical_uid", pa.string()),
                ("entity_granularity", pa.string()),
                ("canonical_name", pa.string()),
                ("formula", pa.string()),
                ("monoisotopic_mass", pa.float64()),
                ("inchi", pa.string()),
                ("inchikey", pa.string()),
                ("inchikey14", pa.string()),
                ("smiles", pa.string()),
                ("charge", pa.int64()),
                ("source_priority", pa.string()),
                ("identity_status", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "chemical_xrefs": pa.schema(
            [
                ("chemical_uid", pa.string()),
                ("xref_source", pa.string()),
                ("xref_id", pa.string()),
                ("xref_type", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "chemical_names": pa.schema(
            [
                ("chemical_uid", pa.string()),
                ("name", pa.string()),
                ("name_type", pa.string()),
                ("language", pa.string()),
                ("source_record_uid", pa.string()),
                ("name_risk", pa.string()),
            ]
        ),
        "metabolite_classes": pa.schema(
            [
                ("class_uid", pa.string()),
                ("class_name", pa.string()),
                ("class_type", pa.string()),
                ("parent_class_uid", pa.string()),
                ("definition", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "chemical_class_members": pa.schema(
            [
                ("class_uid", pa.string()),
                ("chemical_uid", pa.string()),
                ("membership_type", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "trait_entities": pa.schema(
            [
                ("trait_uid", pa.string()),
                ("accession_id", pa.string()),
                ("trait_name", pa.string()),
                ("trait_type", pa.string()),
                ("trait_source", pa.string()),
                ("reported_trait", pa.string()),
                ("summary_statistics_url", pa.string()),
                ("source_record_uid", pa.string()),
                ("identity_scope", pa.string()),
            ]
        ),
        "trait_components": pa.schema(
            [
                ("trait_uid", pa.string()),
                ("component_role", pa.string()),
                ("component_uid", pa.string()),
                ("component_entity_type", pa.string()),
                ("component_name", pa.string()),
                ("direction_semantics", pa.string()),
                ("curation_status", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "input_features": pa.schema(
            [
                ("input_feature_uid", pa.string()),
                ("input_row_id", pa.string()),
                ("feature_label", pa.string()),
                ("feature_type", pa.string()),
                ("effect_value", pa.float64()),
                ("effect_label", pa.string()),
                ("direction", pa.string()),
                ("pvalue", pa.float64()),
                ("padj", pa.float64()),
                ("sample_context_uid", pa.string()),
                ("raw_input_json", pa.string()),
            ]
        ),
        "identity_candidates": pa.schema(
            [
                ("input_feature_uid", pa.string()),
                ("candidate_uid", pa.string()),
                ("candidate_entity_type", pa.string()),
                ("candidate_name", pa.string()),
                ("match_basis", pa.string()),
                ("score", pa.float64()),
                ("margin", pa.float64()),
                ("rank", pa.int64()),
                ("false_match_risk", pa.string()),
                ("blocking_reason", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "identity_decisions": pa.schema(
            [
                ("decision_uid", pa.string()),
                ("input_feature_uid", pa.string()),
                ("decision_status", pa.string()),
                ("accepted_entity_uid", pa.string()),
                ("accepted_entity_type", pa.string()),
                ("decision_rule", pa.string()),
                ("decision_confidence", pa.float64()),
                ("review_status", pa.string()),
                ("reviewer", pa.string()),
                ("reviewed_at", pa.string()),
                ("decision_notes", pa.string()),
            ]
        ),
        "reaction_entities": pa.schema(
            [
                ("reaction_uid", pa.string()),
                ("reaction_name", pa.string()),
                ("reaction_source", pa.string()),
                ("equation", pa.string()),
                ("compartment_uid", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "reaction_participants_v2": pa.schema(
            [
                ("reaction_uid", pa.string()),
                ("chemical_uid", pa.string()),
                ("chemical_source_id", pa.string()),
                ("physical_entity_id", pa.string()),
                ("physical_entity_name", pa.string()),
                ("participant_role", pa.string()),
                ("stoichiometry", pa.float64()),
                ("compartment_uid", pa.string()),
                ("directionality", pa.string()),
                ("relation_source", pa.string()),
                ("semantic_source_uri", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "reaction_catalysts": pa.schema(
            [
                ("reaction_uid", pa.string()),
                ("gene_uid", pa.string()),
                ("protein_uid", pa.string()),
                ("catalyst_role", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "reaction_equations": pa.schema(
            [
                ("reaction_uid", pa.string()),
                ("source_reaction_id", pa.string()),
                ("equation_source", pa.string()),
                ("directionality", pa.string()),
                ("left_to_right_uid", pa.string()),
                ("right_to_left_uid", pa.string()),
                ("bidirectional_uid", pa.string()),
                ("equation_text", pa.string()),
                ("reaction_smiles", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "reaction_side_participants": pa.schema(
            [
                ("reaction_uid", pa.string()),
                ("chemical_uid", pa.string()),
                ("chemical_source_id", pa.string()),
                ("side", pa.string()),
                ("participant_role", pa.string()),
                ("stoichiometry", pa.float64()),
                ("directionality", pa.string()),
                ("relation_source", pa.string()),
                ("physical_entity_id", pa.string()),
                ("physical_entity_name", pa.string()),
                ("semantic_source_uri", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "reaction_xrefs": pa.schema(
            [
                ("reaction_uid", pa.string()),
                ("xref_source", pa.string()),
                ("xref_id", pa.string()),
                ("direction", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "enzyme_reaction_links": pa.schema(
            [
                ("reaction_uid", pa.string()),
                ("protein_uid", pa.string()),
                ("protein_source_id", pa.string()),
                ("gene_uid", pa.string()),
                ("enzyme_role", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "reaction_publication_links": pa.schema(
            [
                ("reaction_uid", pa.string()),
                ("pmid", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "pathway_modules": pa.schema(
            [
                ("module_uid", pa.string()),
                ("module_name", pa.string()),
                ("parent_pathway_uid", pa.string()),
                ("module_type", pa.string()),
                ("definition", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "module_members": pa.schema(
            [
                ("module_uid", pa.string()),
                ("member_uid", pa.string()),
                ("member_type", pa.string()),
                ("member_role", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "evidence_contexts": pa.schema(
            [
                ("context_uid", pa.string()),
                ("cancer_type_uid", pa.string()),
                ("tissue_uid", pa.string()),
                ("cell_type_uid", pa.string()),
                ("cell_state_uid", pa.string()),
                ("species", pa.string()),
                ("model_system", pa.string()),
                ("comparison", pa.string()),
                ("context_text", pa.string()),
            ]
        ),
        "evidence_sentences": pa.schema(
            [
                ("sentence_uid", pa.string()),
                ("pmid", pa.string()),
                ("pmcid", pa.string()),
                ("section", pa.string()),
                ("sentence_text", pa.string()),
                ("sentence_hash", pa.string()),
                ("source_record_uid", pa.string()),
            ]
        ),
        "evidence_assertions": pa.schema(
            [
                ("assertion_uid", pa.string()),
                ("subject_uid", pa.string()),
                ("subject_type", pa.string()),
                ("predicate", pa.string()),
                ("object_uid", pa.string()),
                ("object_type", pa.string()),
                ("polarity", pa.string()),
                ("direction", pa.string()),
                ("evidence_class", pa.string()),
                ("support_status", pa.string()),
                ("context_uid", pa.string()),
                ("method", pa.string()),
                ("source_record_uid", pa.string()),
                ("sentence_uid", pa.string()),
                ("confidence_components_json", pa.string()),
            ]
        ),
        "fact_candidates": pa.schema(
            [
                ("candidate_fact_uid", pa.string()),
                ("candidate_fact_type", pa.string()),
                ("subject_uid", pa.string()),
                ("subject_type", pa.string()),
                ("predicate", pa.string()),
                ("object_uid", pa.string()),
                ("object_type", pa.string()),
                ("relation_source", pa.string()),
                ("source_record_uid", pa.string()),
                ("context_uid", pa.string()),
                ("direction", pa.string()),
                ("role", pa.string()),
                ("evidence_assertion_uids", list_str),
                ("blocking_reasons", list_str),
                ("readiness_score", pa.float64()),
                ("metadata_json", pa.string()),
            ]
        ),
        "mechanism_ready_facts": pa.schema(
            [
                ("fact_uid", pa.string()),
                ("fact_type", pa.string()),
                ("subject_uid", pa.string()),
                ("subject_type", pa.string()),
                ("predicate", pa.string()),
                ("object_uid", pa.string()),
                ("object_type", pa.string()),
                ("context_uid", pa.string()),
                ("direction", pa.string()),
                ("role", pa.string()),
                ("readiness_tier", pa.string()),
                ("allowed_claim_scope", pa.string()),
                ("blocking_reasons", list_str),
                ("supporting_assertion_uids", list_str),
                ("source_record_uids", list_str),
                ("identity_decision_uids", list_str),
                ("confidence_components_json", pa.string()),
                ("boundary_text", pa.string()),
                ("metadata_json", pa.string()),
            ]
        ),
    }
    return schemas[name]


def read_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    return pq.read_table(path).to_pylist()


def write_table(path: Path, rows: list[dict[str, Any]], schema_name: str) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    schema = table_schema(schema_name)
    table = pa.Table.from_pylist(rows, schema=schema) if rows else pa.Table.from_pylist([], schema=schema)
    pq.write_table(table, path, compression="snappy")
    return {
        "table": schema_name,
        "path": str(path),
        "rows": len(rows),
        "schema_hash": short_hash(str(schema), 16),
    }


def read_tsv_from_tar(path: Path, member_name: str) -> list[list[str]]:
    if not path.exists():
        return []
    with tarfile.open(path, "r:gz") as archive:
        try:
            member = archive.getmember(member_name)
        except KeyError:
            return []
        handle = archive.extractfile(member)
        if handle is None:
            return []
        text = handle.read().decode("utf-8", errors="replace").splitlines()
    return [row for row in csv.reader(text, delimiter="\t") if row]


def tsv_records(rows: list[list[str]]) -> list[dict[str, str]]:
    if not rows:
        return []
    header = [str(cell or "").strip() for cell in rows[0]]
    return [
        {header[index]: str(cell or "").strip() for index, cell in enumerate(row) if index < len(header)}
        for row in rows[1:]
    ]


def split_reaction_smiles(value: str) -> tuple[list[str], list[str]]:
    if ">>" not in value:
        return [], []
    left, right = value.split(">>", 1)
    return (
        [part.strip() for part in left.split(".") if part.strip()],
        [part.strip() for part in right.split(".") if part.strip()],
    )


def open_maybe_gzip(path: Path):
    if path.suffix.casefold() == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def xml_local_name(tag: Any) -> str:
    text = str(tag or "")
    if "}" in text:
        return text.rsplit("}", 1)[-1]
    return text.rsplit("#", 1)[-1]


def xml_attr(element: Any, local_name: str) -> str:
    target = local_name.casefold()
    for key, value in getattr(element, "attrib", {}).items():
        if xml_local_name(key).casefold() == target:
            return str(value or "")
    return ""


def xml_resource(element: Any) -> str:
    return xml_attr(element, "resource") or xml_attr(element, "about") or xml_attr(element, "ID")


def resource_key(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if text.startswith("#"):
        return text[1:]
    if "/" in text or "#" in text:
        return re.split(r"[/#]", text)[-1]
    return text


def resource_refs(element: Any, child_names: set[str]) -> list[str]:
    refs = []
    for child in element:
        if xml_local_name(child.tag).casefold() in child_names:
            ref = xml_resource(child)
            if ref:
                refs.append(ref)
    return refs


def child_text(element: Any, child_names: set[str]) -> str:
    for child in element:
        if xml_local_name(child.tag).casefold() in child_names and child.text:
            return str(child.text).strip()
    return ""


def extract_chebi_keys_from_text(value: Any) -> list[str]:
    text = str(value or "")
    return [match.group(1) for match in re.finditer(r"CHEBI[:_/ ]+(\d+)", text, flags=re.IGNORECASE)]


def extract_reactome_ids_from_text(value: Any) -> list[str]:
    return [strip_reactome_version(match.group(0)) for match in re.finditer(r"R-[A-Z]+-\d+(?:\.\d+)?", str(value or ""))]


def extract_rhea_ids_from_text(value: Any) -> list[str]:
    ids = []
    for match in re.finditer(r"(?:RHEA[:_/ -]*)?(\d{3,})", str(value or ""), flags=re.IGNORECASE):
        token = match.group(1)
        if token:
            ids.append(token)
    return ids


def semantic_source_paths(workspace: Path, release_id: str) -> dict[str, list[Path]]:
    rhea_root = workspace / "raw_lake" / "rhea" / release_id
    reactome_root = workspace / "raw_lake" / "reactome" / release_id
    return {
        "rhea_biopax": [
            rhea_root / "rhea.rdf.gz",
            rhea_root / "rhea.rdf",
            rhea_root / "rhea-biopax.owl.gz",
            rhea_root / "rhea-biopax.owl",
        ],
        "reactome_sbml": [
            reactome_root / "Homo_sapiens.sbml",
            reactome_root / "Homo_sapiens.sbml.gz",
            reactome_root / "homo_sapiens.3.1.sbml",
            reactome_root / "homo_sapiens.3.1.sbml.gz",
        ],
        "reactome_biopax": [
            reactome_root / "biopax-level3.owl.gz",
            reactome_root / "biopax-level3.owl",
            reactome_root / "Homo_sapiens.owl.gz",
            reactome_root / "Homo_sapiens.owl",
        ],
    }


def source_uid(source_name: str, release_id: str, raw_record_id: str) -> str:
    return stable_uid("src", source_name, release_id, raw_record_id)


def source_record(source_name: str, source_type: str, release_id: str, path: Path, raw_record_id: str, row: dict[str, Any] | None = None) -> dict[str, Any]:
    return {
        "source_record_uid": source_uid(source_name, release_id, raw_record_id),
        "source_name": source_name,
        "source_type": source_type,
        "source_version": release_id,
        "download_date": "",
        "license_id": str((row or {}).get("license_id") or ""),
        "file_path": str(path),
        "checksum": str((row or {}).get("checksum") or ""),
        "parser_name": "database_accuracy_v2",
        "parser_hash": short_hash(Path(__file__).read_text(encoding="utf-8"), 16),
        "raw_record_id": raw_record_id,
        "raw_payload_json": json_dumps(row or {}),
    }


def add_unique(target: list[dict[str, Any]], seen: set[tuple[Any, ...]], key: tuple[Any, ...], row: dict[str, Any]) -> None:
    if key in seen:
        return
    seen.add(key)
    target.append(row)


def add_reaction_side_fact_rows(
    *,
    reaction_uid: str,
    reaction_external_id: str,
    chemical_uid: str,
    chemical_source_id: str,
    side: str,
    participant_role: str,
    stoichiometry: float | None,
    directionality: str,
    relation_source: str,
    physical_entity_id: str,
    physical_entity_name: str,
    semantic_source_uri: str,
    source_record_uid: str,
    reaction_side_participants: list[dict[str, Any]],
    participants_v2: list[dict[str, Any]],
    side_seen: set[tuple[str, str, str, str]],
    participant_seen: set[tuple[str, str, str, str]],
) -> None:
    add_unique(
        reaction_side_participants,
        side_seen,
        (reaction_uid, chemical_uid, side, source_record_uid),
        {
            "reaction_uid": reaction_uid,
            "chemical_uid": chemical_uid,
            "chemical_source_id": chemical_source_id,
            "side": side,
            "participant_role": participant_role,
            "stoichiometry": stoichiometry,
            "directionality": directionality or "unknown",
            "relation_source": relation_source,
            "physical_entity_id": physical_entity_id,
            "physical_entity_name": physical_entity_name,
            "semantic_source_uri": semantic_source_uri,
            "source_record_uid": source_record_uid,
        },
    )
    add_unique(
        participants_v2,
        participant_seen,
        (reaction_uid, chemical_uid, participant_role, source_record_uid),
        {
            "reaction_uid": reaction_uid,
            "chemical_uid": chemical_uid,
            "chemical_source_id": chemical_source_id,
            "physical_entity_id": physical_entity_id or f"{reaction_external_id}:{side}",
            "physical_entity_name": physical_entity_name,
            "participant_role": participant_role,
            "stoichiometry": stoichiometry,
            "compartment_uid": "",
            "directionality": directionality or "unknown",
            "relation_source": relation_source,
            "semantic_source_uri": semantic_source_uri,
            "source_record_uid": source_record_uid,
        },
    )


def add_unresolved_semantic_participant_row(
    *,
    reaction_uid: str,
    reaction_external_id: str,
    side: str,
    participant_role: str,
    stoichiometry: float | None,
    directionality: str,
    relation_source: str,
    physical_entity_id: str,
    physical_entity_name: str,
    semantic_source_uri: str,
    source_record_uid: str,
    participants_v2: list[dict[str, Any]],
    participant_seen: set[tuple[str, str, str, str]],
) -> None:
    add_unique(
        participants_v2,
        participant_seen,
        (reaction_uid, physical_entity_id or physical_entity_name, participant_role, source_record_uid),
        {
            "reaction_uid": reaction_uid,
            "chemical_uid": "",
            "chemical_source_id": "",
            "physical_entity_id": physical_entity_id or f"{reaction_external_id}:{side}:unresolved",
            "physical_entity_name": physical_entity_name,
            "participant_role": participant_role,
            "stoichiometry": stoichiometry,
            "compartment_uid": "",
            "directionality": directionality or "unknown",
            "relation_source": relation_source,
            "semantic_source_uri": semantic_source_uri,
            "source_record_uid": source_record_uid,
        },
    )


def parse_biopax_semantic_sides(
    path: Path,
    *,
    release_id: str,
    chebi_to_chemical: dict[str, str],
    reaction_by_external_base: dict[str, dict[str, Any]],
    rhea_to_reaction_uids: dict[str, list[str]],
    direction_by_rhea: dict[str, dict[str, str]],
    reaction_equations: list[dict[str, Any]],
    reaction_side_participants: list[dict[str, Any]],
    participants_v2: list[dict[str, Any]],
    reaction_xrefs: list[dict[str, Any]],
    equation_seen: set[tuple[str, str]],
    side_seen: set[tuple[str, str, str, str]],
    participant_seen: set[tuple[str, str, str, str]],
    xref_seen: set[tuple[str, str, str]],
) -> int:
    if etree is None or not path.exists():
        return 0
    try:
        with open_maybe_gzip(path) as handle:
            tree = etree.parse(handle)
    except Exception:
        return 0

    root = tree.getroot()
    elements_by_key: dict[str, Any] = {}
    xrefs_by_key: dict[str, dict[str, str]] = {}
    for element in root.iter():
        about = xml_attr(element, "about") or xml_attr(element, "ID")
        key = resource_key(about)
        if key:
            elements_by_key.setdefault(key, element)
            elements_by_key.setdefault(about, element)
        local = xml_local_name(element.tag).casefold()
        if local.endswith("xref"):
            db = child_text(element, {"db"})
            xid = child_text(element, {"id"})
            if not db:
                db = xml_attr(element, "db")
            if not xid:
                xid = xml_attr(element, "id")
            if key and (db or xid):
                xrefs_by_key[key] = {"db": db, "id": xid}

    def element_display_name(element: Any) -> str:
        return child_text(element, {"display-name", "displayname", "standard-name", "name"}) or xml_attr(element, "displayName")

    def xrefs_for_element(element: Any) -> list[dict[str, str]]:
        xrefs = []
        for ref in resource_refs(element, {"xref"}):
            info = xrefs_by_key.get(resource_key(ref)) or xrefs_by_key.get(ref)
            if info:
                xrefs.append(info)
        text = " ".join([str(element.get(key) or "") for key in element.attrib])
        if getattr(element, "text", None):
            text += f" {element.text}"
        for child in element.iter():
            if getattr(child, "text", None):
                text += f" {child.text}"
        for chebi in extract_chebi_keys_from_text(text):
            xrefs.append({"db": "ChEBI", "id": f"CHEBI:{chebi}"})
        return xrefs

    def chebi_keys_for_entity_ref(ref: str) -> tuple[list[str], str, str]:
        key = resource_key(ref)
        element = elements_by_key.get(key)
        if element is None:
            element = elements_by_key.get(ref)
        if element is None:
            return [], key, ""
        name = element_display_name(element)
        keys: list[str] = []
        for info in xrefs_for_element(element):
            db = str(info.get("db") or "").casefold()
            xid = str(info.get("id") or "")
            if "chebi" in db or "chebi" in xid.casefold():
                keys.extend(extract_chebi_keys_from_text(xid) or [chebi_key(xid)])
        for child_ref in resource_refs(element, {"entityreference", "entity-reference", "component"}):
            child_keys, _child_id, child_name = chebi_keys_for_entity_ref(child_ref)
            keys.extend(child_keys)
            if not name:
                name = child_name
        return sorted({key for key in keys if key}), key, name

    relation_source = "rhea_biopax_explicit_side" if "rhea" in path.name.casefold() else "reactome_biopax_conversion"
    source_name = f"raw_lake.{relation_source}"
    promoted_count = 0
    reaction_locals = {"biochemicalreaction", "conversion", "transport", "transportwithbiochemicalreaction"}
    for element in root.iter():
        if xml_local_name(element.tag).casefold() not in reaction_locals:
            continue
        reaction_uri = xml_attr(element, "about") or xml_attr(element, "ID")
        reaction_key = resource_key(reaction_uri)
        reaction_xref_infos = xrefs_for_element(element)
        rhea_ids: list[str] = []
        reactome_ids: list[str] = []
        for info in reaction_xref_infos:
            db = str(info.get("db") or "").casefold()
            xid = str(info.get("id") or "")
            if "rhea" in db or "rhea" in xid.casefold():
                rhea_ids.extend(extract_rhea_ids_from_text(xid))
            if "reactome" in db or "reactome" in xid.casefold():
                reactome_ids.extend(extract_reactome_ids_from_text(xid))
        rhea_ids.extend(extract_rhea_ids_from_text(reaction_uri))
        reactome_ids.extend(extract_reactome_ids_from_text(reaction_uri))
        reaction_uids: set[str] = set()
        direction = child_text(element, {"conversiondirection", "conversion-direction"}).upper().replace("_", "-")
        for rhea_id in sorted(set(rhea_ids)):
            reaction_uids.update(rhea_to_reaction_uids.get(rhea_id, []))
            if not direction:
                direction = direction_by_rhea.get(rhea_id, {}).get("direction", "")
        for reactome_id in sorted(set(reactome_ids)):
            reaction = reaction_by_external_base.get(strip_reactome_version(reactome_id), {})
            if reaction.get("reaction_uid"):
                reaction_uids.add(str(reaction["reaction_uid"]))
        if not reaction_uids:
            continue
        if direction in {"LEFT-TO-RIGHT", "LEFT_TO_RIGHT"}:
            direction = "LR"
        elif direction in {"RIGHT-TO-LEFT", "RIGHT_TO_LEFT"}:
            direction = "RL"
        elif direction in {"REVERSIBLE", "BIDIRECTIONAL", "BIDIRECTIONAL-REACTION"}:
            direction = "BI"
        direction = direction or "unknown"
        left_refs = resource_refs(element, {"left"})
        right_refs = resource_refs(element, {"right"})
        source_id = f"{reaction_key}|{'|'.join(sorted(set(rhea_ids + reactome_ids)))}"
        for reaction_uid in sorted(reaction_uids):
            for rhea_id in sorted(set(rhea_ids)):
                src_uid = source_uid(source_name + ".xref", release_id, f"{reaction_uid}|RHEA|{rhea_id}|{reaction_key}")
                add_unique(
                    reaction_xrefs,
                    xref_seen,
                    (reaction_uid, "RHEA", rhea_id),
                    {"reaction_uid": reaction_uid, "xref_source": "RHEA", "xref_id": rhea_id, "direction": direction, "source_record_uid": src_uid},
                )
            add_unique(
                reaction_equations,
                equation_seen,
                (reaction_uid, f"{relation_source}:{source_id}"),
                {
                    "reaction_uid": reaction_uid,
                    "source_reaction_id": sorted(set([f"RHEA:{x}" for x in rhea_ids] + reactome_ids))[0] if (rhea_ids or reactome_ids) else reaction_key,
                    "equation_source": relation_source,
                    "directionality": direction,
                    "left_to_right_uid": "",
                    "right_to_left_uid": "",
                    "bidirectional_uid": "",
                    "equation_text": "",
                    "reaction_smiles": "",
                    "source_record_uid": source_uid(source_name + ".reaction", release_id, source_id),
                },
            )
            for side, role, refs in (("left", "substrate", left_refs), ("right", "product", right_refs)):
                for ref in refs:
                    chebi_keys, entity_key, entity_name = chebi_keys_for_entity_ref(ref)
                    if not chebi_keys:
                        add_unresolved_semantic_participant_row(
                            reaction_uid=reaction_uid,
                            reaction_external_id=reaction_key,
                            side=side,
                            participant_role=role,
                            stoichiometry=None,
                            directionality=direction,
                            relation_source=relation_source,
                            physical_entity_id=entity_key,
                            physical_entity_name=entity_name,
                            semantic_source_uri=ref,
                            source_record_uid=source_uid(source_name + ".participant", release_id, f"{reaction_uid}|{side}|{entity_key}|unresolved"),
                            participants_v2=participants_v2,
                            participant_seen=participant_seen,
                        )
                        continue
                    for key in chebi_keys:
                        chemical_uid = chebi_to_chemical.get(key, "")
                        src_uid = source_uid(source_name + ".participant", release_id, f"{reaction_uid}|{side}|CHEBI:{key}|{entity_key}")
                        if not chemical_uid:
                            add_unresolved_semantic_participant_row(
                                reaction_uid=reaction_uid,
                                reaction_external_id=reaction_key,
                                side=side,
                                participant_role=role,
                                stoichiometry=None,
                                directionality=direction,
                                relation_source=relation_source,
                                physical_entity_id=entity_key,
                                physical_entity_name=entity_name,
                                semantic_source_uri=ref,
                                source_record_uid=src_uid,
                                participants_v2=participants_v2,
                                participant_seen=participant_seen,
                            )
                            continue
                        add_reaction_side_fact_rows(
                            reaction_uid=reaction_uid,
                            reaction_external_id=reaction_key,
                            chemical_uid=chemical_uid,
                            chemical_source_id=f"CHEBI:{key}",
                            side=side,
                            participant_role=role,
                            stoichiometry=None,
                            directionality=direction,
                            relation_source=relation_source,
                            physical_entity_id=entity_key,
                            physical_entity_name=entity_name,
                            semantic_source_uri=ref,
                            source_record_uid=src_uid,
                            reaction_side_participants=reaction_side_participants,
                            participants_v2=participants_v2,
                            side_seen=side_seen,
                            participant_seen=participant_seen,
                        )
                        promoted_count += 1
    return promoted_count


def parse_rhea_rdf_semantic_sides(
    path: Path,
    *,
    release_id: str,
    chebi_to_chemical: dict[str, str],
    rhea_to_reaction_uids: dict[str, list[str]],
    direction_by_rhea: dict[str, dict[str, str]],
    reaction_equations: list[dict[str, Any]],
    reaction_side_participants: list[dict[str, Any]],
    participants_v2: list[dict[str, Any]],
    reaction_xrefs: list[dict[str, Any]],
    equation_seen: set[tuple[str, str]],
    side_seen: set[tuple[str, str, str, str]],
    participant_seen: set[tuple[str, str, str, str]],
    xref_seen: set[tuple[str, str, str]],
) -> int:
    if etree is None or not path.exists():
        return 0
    target_rhea_ids = {str(rhea_id) for rhea_id, reaction_uids in rhea_to_reaction_uids.items() if reaction_uids}
    if not target_rhea_ids:
        return 0
    reactions: dict[str, dict[str, Any]] = {}
    try:
        with open_maybe_gzip(path) as handle:
            context = etree.iterparse(handle, events=("end",), tag="{*}Description", recover=True)
            for _event, element in context:
                accession = child_text(element, {"accession"})
                if accession.startswith("RHEA:"):
                    rhea_id = chebi_key(accession)
                    if rhea_id in target_rhea_ids:
                        substrate_refs = resource_refs(element, {"substrates"})
                        product_refs = resource_refs(element, {"products"})
                        reversible_refs = resource_refs(element, {"substratesorproducts"})
                        reactions[rhea_id] = {
                            "rhea_id": rhea_id,
                            "equation": child_text(element, {"equation"}),
                            "substrates": [resource_key(ref) for ref in substrate_refs],
                            "products": [resource_key(ref) for ref in product_refs],
                            "substrates_or_products": [resource_key(ref) for ref in reversible_refs],
                        }
                element.clear()
                parent = element.getparent()
                while parent is not None and element.getprevious() is not None:
                    del parent[0]
    except Exception:
        return 0

    target_side_keys = {
        side_key
        for reaction in reactions.values()
        for side_key in (
            list(reaction.get("substrates") or [])
            + list(reaction.get("products") or [])
            + list(reaction.get("substrates_or_products") or [])
        )
    }
    side_participants: dict[str, list[str]] = {}
    try:
        with open_maybe_gzip(path) as handle:
            context = etree.iterparse(handle, events=("end",), tag="{*}Description", recover=True)
            for _event, element in context:
                key = resource_key(xml_attr(element, "about"))
                if key in target_side_keys:
                    contains_refs = resource_refs(element, {"contains", "contains1"})
                    if contains_refs:
                        side_participants.setdefault(key, []).extend(resource_key(ref) for ref in contains_refs)
                element.clear()
                parent = element.getparent()
                while parent is not None and element.getprevious() is not None:
                    del parent[0]
    except Exception:
        return 0

    target_participant_keys = {participant for participants in side_participants.values() for participant in participants}
    participant_compound: dict[str, str] = {}
    try:
        with open_maybe_gzip(path) as handle:
            context = etree.iterparse(handle, events=("end",), tag="{*}Description", recover=True)
            for _event, element in context:
                key = resource_key(xml_attr(element, "about"))
                if key in target_participant_keys:
                    compound_refs = resource_refs(element, {"compound"})
                    if compound_refs:
                        participant_compound[key] = resource_key(compound_refs[0])
                element.clear()
                parent = element.getparent()
                while parent is not None and element.getprevious() is not None:
                    del parent[0]
    except Exception:
        return 0

    target_compound_keys = set(participant_compound.values())
    compound_chebi: dict[str, str] = {}
    try:
        with open_maybe_gzip(path) as handle:
            context = etree.iterparse(handle, events=("end",), tag="{*}Description", recover=True)
            for _event, element in context:
                key = resource_key(xml_attr(element, "about"))
                if key in target_compound_keys:
                    accession = child_text(element, {"accession"})
                    if accession.startswith("CHEBI:"):
                        compound_chebi[key] = chebi_key(accession)
                element.clear()
                parent = element.getparent()
                while parent is not None and element.getprevious() is not None:
                    del parent[0]
    except Exception:
        return 0

    relation_source = "rhea_rdf_explicit_side"
    source_name = "raw_lake.rhea_rdf_explicit_side"
    promoted_count = 0
    for rhea_id, reaction_info in reactions.items():
        reaction_uids = rhea_to_reaction_uids.get(rhea_id, [])
        if not reaction_uids:
            continue
        direction = direction_by_rhea.get(rhea_id, {}).get("direction") or ("BI" if reaction_info.get("substrates_or_products") else "unknown")
        side_roles: list[tuple[str, str, str]] = []
        for side_key in reaction_info.get("substrates") or []:
            side_roles.append((side_key, "left" if side_key.endswith("_L") else "substrate", "substrate"))
        for side_key in reaction_info.get("products") or []:
            side_roles.append((side_key, "right" if side_key.endswith("_R") else "product", "product"))
        for side_key in reaction_info.get("substrates_or_products") or []:
            side = "left" if side_key.endswith("_L") else "right" if side_key.endswith("_R") else "unknown"
            role = "substrate" if side == "left" else "product" if side == "right" else "participant"
            side_roles.append((side_key, side, role))
        source_uri = f"http://rdf.rhea-db.org/{rhea_id}"
        for reaction_uid in reaction_uids:
            add_unique(
                reaction_xrefs,
                xref_seen,
                (reaction_uid, "RHEA", rhea_id),
                {
                    "reaction_uid": reaction_uid,
                    "xref_source": "RHEA",
                    "xref_id": rhea_id,
                    "direction": direction,
                    "source_record_uid": source_uid(source_name + ".xref", release_id, f"{reaction_uid}|RHEA|{rhea_id}"),
                },
            )
            add_unique(
                reaction_equations,
                equation_seen,
                (reaction_uid, f"{relation_source}:{rhea_id}"),
                {
                    "reaction_uid": reaction_uid,
                    "source_reaction_id": f"RHEA:{rhea_id}",
                    "equation_source": relation_source,
                    "directionality": direction,
                    "left_to_right_uid": direction_by_rhea.get(rhea_id, {}).get("lr", ""),
                    "right_to_left_uid": direction_by_rhea.get(rhea_id, {}).get("rl", ""),
                    "bidirectional_uid": direction_by_rhea.get(rhea_id, {}).get("bi", ""),
                    "equation_text": str(reaction_info.get("equation") or ""),
                    "reaction_smiles": "",
                    "source_record_uid": source_uid(source_name + ".reaction", release_id, rhea_id),
                },
            )
            for side_key, side, role in side_roles:
                for participant_key in side_participants.get(side_key, []):
                    compound_key = participant_compound.get(participant_key, "")
                    chebi_id = compound_chebi.get(compound_key, "")
                    entity_uri = f"http://rdf.rhea-db.org/{participant_key}"
                    if not chebi_id:
                        add_unresolved_semantic_participant_row(
                            reaction_uid=reaction_uid,
                            reaction_external_id=f"RHEA:{rhea_id}",
                            side=side,
                            participant_role=role,
                            stoichiometry=None,
                            directionality=direction,
                            relation_source=relation_source,
                            physical_entity_id=participant_key,
                            physical_entity_name="",
                            semantic_source_uri=entity_uri,
                            source_record_uid=source_uid(source_name + ".participant", release_id, f"{reaction_uid}|{side}|{participant_key}|unresolved"),
                            participants_v2=participants_v2,
                            participant_seen=participant_seen,
                        )
                        continue
                    chemical_uid = chebi_to_chemical.get(chebi_id, "")
                    src_uid = source_uid(source_name + ".participant", release_id, f"{reaction_uid}|{side}|CHEBI:{chebi_id}|{participant_key}")
                    if not chemical_uid:
                        add_unresolved_semantic_participant_row(
                            reaction_uid=reaction_uid,
                            reaction_external_id=f"RHEA:{rhea_id}",
                            side=side,
                            participant_role=role,
                            stoichiometry=None,
                            directionality=direction,
                            relation_source=relation_source,
                            physical_entity_id=participant_key,
                            physical_entity_name="",
                            semantic_source_uri=entity_uri,
                            source_record_uid=src_uid,
                            participants_v2=participants_v2,
                            participant_seen=participant_seen,
                        )
                        continue
                    add_reaction_side_fact_rows(
                        reaction_uid=reaction_uid,
                        reaction_external_id=f"RHEA:{rhea_id}",
                        chemical_uid=chemical_uid,
                        chemical_source_id=f"CHEBI:{chebi_id}",
                        side=side,
                        participant_role=role,
                        stoichiometry=None,
                        directionality=direction,
                        relation_source=relation_source,
                        physical_entity_id=participant_key,
                        physical_entity_name="",
                        semantic_source_uri=entity_uri,
                        source_record_uid=src_uid,
                        reaction_side_participants=reaction_side_participants,
                        participants_v2=participants_v2,
                        side_seen=side_seen,
                        participant_seen=participant_seen,
                    )
                    promoted_count += 1
    return promoted_count


def parse_sbml_semantic_sides(
    path: Path,
    *,
    release_id: str,
    chebi_to_chemical: dict[str, str],
    reaction_by_external_base: dict[str, dict[str, Any]],
    reaction_equations: list[dict[str, Any]],
    reaction_side_participants: list[dict[str, Any]],
    participants_v2: list[dict[str, Any]],
    equation_seen: set[tuple[str, str]],
    side_seen: set[tuple[str, str, str, str]],
    participant_seen: set[tuple[str, str, str, str]],
) -> int:
    if etree is None or not path.exists():
        return 0
    try:
        with open_maybe_gzip(path) as handle:
            tree = etree.parse(handle)
    except Exception:
        return 0

    species: dict[str, dict[str, Any]] = {}
    for element in tree.getroot().iter():
        if xml_local_name(element.tag).casefold() != "species":
            continue
        species_id = xml_attr(element, "id")
        if not species_id:
            continue
        name = xml_attr(element, "name")
        payload = " ".join([species_id, name])
        for child in element.iter():
            payload += f" {xml_attr(child, 'resource')} {child.text or ''}"
        species[species_id] = {
            "id": species_id,
            "name": name,
            "chebi_keys": sorted(set(extract_chebi_keys_from_text(payload))),
            "source_uri": f"{path.name}#{species_id}",
        }

    promoted_count = 0
    relation_source = "reactome_sbml_species_reference"
    source_name = "raw_lake.reactome_sbml_species_reference"
    for element in tree.getroot().iter():
        if xml_local_name(element.tag).casefold() != "reaction":
            continue
        reaction_id = xml_attr(element, "id")
        reaction_name = xml_attr(element, "name")
        reactome_ids = extract_reactome_ids_from_text(f"{reaction_id} {reaction_name}")
        reaction_uids = {
            str(reaction_by_external_base.get(strip_reactome_version(reactome_id), {}).get("reaction_uid") or "")
            for reactome_id in reactome_ids
        }
        reaction_uids = {uid for uid in reaction_uids if uid}
        if not reaction_uids:
            continue
        direction = "BI" if str(xml_attr(element, "reversible")).casefold() == "true" else "LR"
        for reaction_uid in sorted(reaction_uids):
            add_unique(
                reaction_equations,
                equation_seen,
                (reaction_uid, f"{relation_source}:{reaction_id}"),
                {
                    "reaction_uid": reaction_uid,
                    "source_reaction_id": sorted(set(reactome_ids))[0] if reactome_ids else reaction_id,
                    "equation_source": relation_source,
                    "directionality": direction,
                    "left_to_right_uid": "",
                    "right_to_left_uid": "",
                    "bidirectional_uid": "",
                    "equation_text": reaction_name,
                    "reaction_smiles": "",
                    "source_record_uid": source_uid(source_name + ".reaction", release_id, reaction_id),
                },
            )
            for side_container, side, role in (("listofreactants", "left", "substrate"), ("listofproducts", "right", "product")):
                for container in element:
                    if xml_local_name(container.tag).casefold() != side_container:
                        continue
                    for species_ref in container:
                        if xml_local_name(species_ref.tag).casefold() not in {"speciesreference", "modifierspeciesreference"}:
                            continue
                        species_id = xml_attr(species_ref, "species")
                        info = species.get(species_id, {})
                        stoich_text = xml_attr(species_ref, "stoichiometry")
                        try:
                            stoich = float(stoich_text) if stoich_text else None
                        except ValueError:
                            stoich = None
                        chebi_keys = info.get("chebi_keys") or []
                        if not chebi_keys:
                            add_unresolved_semantic_participant_row(
                                reaction_uid=reaction_uid,
                                reaction_external_id=reaction_id,
                                side=side,
                                participant_role=role,
                                stoichiometry=stoich,
                                directionality=direction,
                                relation_source=relation_source,
                                physical_entity_id=species_id,
                                physical_entity_name=str(info.get("name") or ""),
                                semantic_source_uri=str(info.get("source_uri") or ""),
                                source_record_uid=source_uid(source_name + ".participant", release_id, f"{reaction_id}|{side}|{species_id}|unresolved"),
                                participants_v2=participants_v2,
                                participant_seen=participant_seen,
                            )
                            continue
                        for key in chebi_keys:
                            chemical_uid = chebi_to_chemical.get(key, "")
                            src_uid = source_uid(source_name + ".participant", release_id, f"{reaction_id}|{side}|{species_id}|CHEBI:{key}")
                            if not chemical_uid:
                                add_unresolved_semantic_participant_row(
                                    reaction_uid=reaction_uid,
                                    reaction_external_id=reaction_id,
                                    side=side,
                                    participant_role=role,
                                    stoichiometry=stoich,
                                    directionality=direction,
                                    relation_source=relation_source,
                                    physical_entity_id=species_id,
                                    physical_entity_name=str(info.get("name") or ""),
                                    semantic_source_uri=str(info.get("source_uri") or ""),
                                    source_record_uid=src_uid,
                                    participants_v2=participants_v2,
                                    participant_seen=participant_seen,
                                )
                                continue
                            add_reaction_side_fact_rows(
                                reaction_uid=reaction_uid,
                                reaction_external_id=reaction_id,
                                chemical_uid=chemical_uid,
                                chemical_source_id=f"CHEBI:{key}",
                                side=side,
                                participant_role=role,
                                stoichiometry=stoich,
                                directionality=direction,
                                relation_source=relation_source,
                                physical_entity_id=species_id,
                                physical_entity_name=str(info.get("name") or ""),
                                semantic_source_uri=str(info.get("source_uri") or ""),
                                source_record_uid=src_uid,
                                reaction_side_participants=reaction_side_participants,
                                participants_v2=participants_v2,
                                side_seen=side_seen,
                                participant_seen=participant_seen,
                            )
                            promoted_count += 1
    return promoted_count


def read_european_annotations(workspace: Path) -> dict[str, dict[str, Any]]:
    path = european_source_path(workspace, "European_trait_annotations.csv")
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return {normalize_accession(row.get("accession_id")): dict(row) for row in csv.DictReader(handle) if normalize_accession(row.get("accession_id"))}


def read_european_traits(workspace: Path) -> dict[str, dict[str, Any]]:
    path = european_source_path(workspace, "European.csv")
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = {}
        for row in reader:
            accession = normalize_accession(row.get("accession_id") or row.get("accessionId") or row.get("trait"))
            if accession:
                rows[accession] = dict(row)
        return rows


def trait_identity_scope(trait_type: str) -> str:
    return {
        "metabolite_ratio": "ratio",
        "class_trait": "class_level",
        "metabolite_level": "exact_chemical",
    }.get(trait_type, "unknown")


def component_role(index: int, trait_type: str) -> tuple[str, str]:
    if trait_type == "metabolite_ratio":
        if index == 0:
            return "numerator", "same_direction"
        if index == 1:
            return "denominator", "inverse_direction"
    return "component", "undefined"


def class_type_for_name(value: str) -> str:
    text = value.casefold()
    if "sphingo" in text or "ceramide" in text:
        return "lipid_class"
    if "pool" in text:
        return "pathway_pool"
    if "redox" in text:
        return "redox_pair"
    return "metabolite_family"


def build_entity_store(workspace: Path, release_id: str, normalized_dir: Path, output_dir: Path) -> list[dict[str, Any]]:
    metabolites_path = normalized_dir / "metabolites.parquet"
    xrefs_path = normalized_dir / "metabolite_xrefs.parquet"
    metabolites = read_rows(metabolites_path)
    xrefs = read_rows(xrefs_path)
    annotations = read_european_annotations(workspace)
    traits = read_european_traits(workspace)

    source_records: list[dict[str, Any]] = []
    source_seen: set[tuple[Any, ...]] = set()
    chemical_entities: list[dict[str, Any]] = []
    chemical_xrefs: list[dict[str, Any]] = []
    chemical_names: list[dict[str, Any]] = []
    metabolite_classes: list[dict[str, Any]] = []
    class_members: list[dict[str, Any]] = []
    trait_entities: list[dict[str, Any]] = []
    trait_components: list[dict[str, Any]] = []
    class_seen: set[tuple[Any, ...]] = set()
    member_seen: set[tuple[Any, ...]] = set()
    name_to_chemical: dict[str, str] = {}

    for row in metabolites:
        chemical_uid = str(row.get("metabolite_uid") or "")
        if not chemical_uid:
            continue
        src_uid = source_uid("normalized_store.metabolites", release_id, chemical_uid)
        add_unique(source_records, source_seen, (src_uid,), source_record("normalized_store.metabolites", "curated_db", release_id, metabolites_path, chemical_uid, row))
        synonyms = [normalize_name(value) for value in (row.get("synonyms") or []) if normalize_name(value)]
        chemical_entities.append(
            {
                "chemical_uid": chemical_uid,
                "entity_granularity": entity_granularity(row),
                "canonical_name": str(row.get("canonical_name") or ""),
                "formula": str(row.get("formula") or ""),
                "monoisotopic_mass": row.get("exact_mass"),
                "inchi": "",
                "inchikey": str(row.get("inchikey") or ""),
                "inchikey14": str(row.get("inchikey") or "")[:14],
                "smiles": str(row.get("smiles") or ""),
                "charge": row.get("charge"),
                "source_priority": str(row.get("source_priority") or ""),
                "identity_status": "canonical" if entity_granularity(row) != "unknown" else "review",
                "source_record_uid": src_uid,
            }
        )
        for name, name_type in [(row.get("canonical_name"), "canonical"), *[(syn, "synonym") for syn in synonyms]]:
            cleaned = normalize_name(name)
            if not cleaned:
                continue
            chemical_names.append(
                {
                    "chemical_uid": chemical_uid,
                    "name": cleaned,
                    "name_type": "class_name" if name_risk(cleaned) == "high" and is_class_text(cleaned) else name_type,
                    "language": "",
                    "source_record_uid": src_uid,
                    "name_risk": name_risk(cleaned),
                }
            )
            if name_risk(cleaned) == "low":
                name_to_chemical.setdefault(name_key(cleaned), chemical_uid)
            if is_class_text(cleaned) and not LIPID_SPECIES_RE.search(cleaned):
                class_uid = stable_uid("class", name_key(cleaned))
                add_unique(
                    metabolite_classes,
                    class_seen,
                    (class_uid,),
                    {
                        "class_uid": class_uid,
                        "class_name": cleaned,
                        "class_type": class_type_for_name(cleaned),
                        "parent_class_uid": "",
                        "definition": "Class-level surface separated from exact chemical identity.",
                        "source_record_uid": src_uid,
                    },
                )
                add_unique(
                    class_members,
                    member_seen,
                    (class_uid, chemical_uid),
                    {
                        "class_uid": class_uid,
                        "chemical_uid": chemical_uid,
                        "membership_type": "possible_member",
                        "source_record_uid": src_uid,
                    },
                )

    for row in xrefs:
        chemical_uid = str(row.get("metabolite_uid") or "")
        if not chemical_uid:
            continue
        raw_id = str(row.get("xref_key") or row.get("xref_id") or "")
        src_uid = source_uid("normalized_store.metabolite_xrefs", release_id, raw_id or chemical_uid)
        add_unique(source_records, source_seen, (src_uid,), source_record("normalized_store.metabolite_xrefs", "curated_db", release_id, xrefs_path, raw_id or chemical_uid, row))
        chemical_xrefs.append(
            {
                "chemical_uid": chemical_uid,
                "xref_source": str(row.get("xref_source") or ""),
                "xref_id": str(row.get("xref_id") or ""),
                "xref_type": "primary" if raw_id and raw_id in (row.get("xref_key") or "") else "secondary",
                "source_record_uid": src_uid,
            }
        )

    accession_rows = {**traits, **annotations}
    for accession in sorted(accession_rows):
        annotation = annotations.get(accession, {})
        trait_row = traits.get(accession, {})
        reported_trait = str(annotation.get("reported_trait") or trait_row.get("reported_trait") or trait_row.get("reportedTrait") or "")
        trait_type = classify_trait_type(reported_trait, annotation)
        trait_uid = stable_uid("trait", accession)
        src_uid = source_uid("raw_lake.European", release_id, accession)
        add_unique(
            source_records,
            source_seen,
            (src_uid,),
            source_record(
                "raw_lake.European",
                "curated_db",
                release_id,
                european_source_path(workspace, "European_trait_annotations.csv").parent,
                accession,
                {**trait_row, **annotation},
            ),
        )
        trait_entities.append(
            {
                "trait_uid": trait_uid,
                "accession_id": accession,
                "trait_name": reported_trait or accession,
                "trait_type": trait_type,
                "trait_source": "GWAS",
                "reported_trait": reported_trait,
                "summary_statistics_url": str(annotation.get("summary_statistics_url") or trait_row.get("summary_statistics_url") or trait_row.get("summaryStatistics") or ""),
                "source_record_uid": src_uid,
                "identity_scope": trait_identity_scope(trait_type),
            }
        )
        names = split_values(annotation.get("mapped_names")) or split_values(annotation.get("candidate_names"))
        if not names and reported_trait:
            names = re.split(r"\s+(?:to|/)\s+|\s*/\s*", reported_trait, flags=re.IGNORECASE) if is_ratio_text(reported_trait) else [reported_trait]
        for index, component_name in enumerate([normalize_name(name) for name in names if normalize_name(name)]):
            role, direction = component_role(index, trait_type)
            component_uid = name_to_chemical.get(name_key(component_name), "")
            component_type = "chemical" if component_uid else ("class" if is_class_text(component_name) else "unknown")
            if component_type == "class":
                component_uid = stable_uid("class", name_key(component_name))
                add_unique(
                    metabolite_classes,
                    class_seen,
                    (component_uid,),
                    {
                        "class_uid": component_uid,
                        "class_name": component_name,
                        "class_type": class_type_for_name(component_name),
                        "parent_class_uid": "",
                        "definition": "Trait component is class-level and requires review before exact identity use.",
                        "source_record_uid": src_uid,
                    },
                )
            trait_components.append(
                {
                    "trait_uid": trait_uid,
                    "component_role": role,
                    "component_uid": component_uid,
                    "component_entity_type": component_type,
                    "component_name": component_name,
                    "direction_semantics": direction,
                    "curation_status": "rule_inferred" if component_uid else "review_required",
                    "source_record_uid": src_uid,
                }
            )

    audits = []
    root = output_dir / "entity_store"
    audits.append(write_table(root / "source_records.parquet", source_records, "source_records"))
    audits.append(write_table(root / "chemical_entities.parquet", chemical_entities, "chemical_entities"))
    audits.append(write_table(root / "chemical_xrefs.parquet", chemical_xrefs, "chemical_xrefs"))
    audits.append(write_table(root / "chemical_names.parquet", chemical_names, "chemical_names"))
    audits.append(write_table(root / "metabolite_classes.parquet", metabolite_classes, "metabolite_classes"))
    audits.append(write_table(root / "chemical_class_members.parquet", class_members, "chemical_class_members"))
    audits.append(write_table(root / "trait_entities.parquet", trait_entities, "trait_entities"))
    audits.append(write_table(root / "trait_components.parquet", trait_components, "trait_components"))
    return audits


def feature_type_for_trait(trait_type: str) -> str:
    if trait_type == "metabolite_ratio":
        return "ratio"
    if trait_type == "class_trait":
        return "lipid_class"
    if trait_type in {"metabolite_level", "unknown_trait"}:
        return "trait_score"
    return "unknown_peak"


def decision_for_feature(feature_type: str, trait_uid: str, class_uid: str, chemical_uid: str) -> tuple[str, str, str, str, float, str]:
    if feature_type == "ratio":
        if trait_uid:
            return "accepted_trait", trait_uid, "trait", "ratio_not_exact_abundance", 0.8, "needs_review"
        return "ambiguous", "", "", "ratio_not_exact_abundance;missing_trait_entity", 0.0, "needs_review"
    if feature_type == "lipid_class":
        if class_uid:
            return "accepted_class", class_uid, "class", "class_name_not_exact_identity", 0.65, "needs_review"
        return "ambiguous", "", "", "class_name_not_exact_identity", 0.0, "needs_review"
    if feature_type == "trait_score":
        if trait_uid:
            return "accepted_trait", trait_uid, "trait", "trait_score_not_direct_abundance", 0.7, "needs_review"
        return "ambiguous", "", "", "trait_score_without_registered_trait", 0.0, "needs_review"
    if chemical_uid:
        return "accepted_exact", chemical_uid, "chemical", "stable_source_entity", 0.9, "auto"
    return "unmatched", "", "", "no_candidate_entity", 0.0, "needs_review"


def build_identity_resolution_store(workspace: Path, release_id: str, output_dir: Path) -> list[dict[str, Any]]:
    entity_root = output_dir / "entity_store"
    traits = read_rows(entity_root / "trait_entities.parquet")
    components = read_rows(entity_root / "trait_components.parquet")
    classes = read_rows(entity_root / "metabolite_classes.parquet")
    chemicals = read_rows(entity_root / "chemical_entities.parquet")
    names = read_rows(entity_root / "chemical_names.parquet")
    trait_by_uid = {row["trait_uid"]: row for row in traits}
    class_by_name = {name_key(row.get("class_name")): row for row in classes}
    chemical_by_name = {
        name_key(row.get("name")): row
        for row in names
        if row.get("name_risk") == "low" and row.get("chemical_uid")
    }
    chemical_by_uid = {row["chemical_uid"]: row for row in chemicals}

    input_features: list[dict[str, Any]] = []
    identity_candidates: list[dict[str, Any]] = []
    identity_decisions: list[dict[str, Any]] = []

    for trait in traits:
        trait_uid = str(trait.get("trait_uid") or "")
        accession = str(trait.get("accession_id") or trait_uid)
        feature_uid = stable_uid("input_feature", "trait", accession)
        feature_type = feature_type_for_trait(str(trait.get("trait_type") or ""))
        input_features.append(
            {
                "input_feature_uid": feature_uid,
                "input_row_id": accession,
                "feature_label": str(trait.get("trait_name") or accession),
                "feature_type": feature_type,
                "effect_value": None,
                "effect_label": "",
                "direction": "",
                "pvalue": None,
                "padj": None,
                "sample_context_uid": "",
                "raw_input_json": json_dumps(trait),
            }
        )
        identity_candidates.append(
            {
                "input_feature_uid": feature_uid,
                "candidate_uid": trait_uid,
                "candidate_entity_type": "trait",
                "candidate_name": str(trait.get("trait_name") or accession),
                "match_basis": "trait_annotation",
                "score": 100.0,
                "margin": 100.0,
                "rank": 1,
                "false_match_risk": "low" if feature_type == "trait_score" else "medium",
                "blocking_reason": "ratio_not_exact_abundance" if feature_type == "ratio" else ("class_name_not_exact_identity" if feature_type == "lipid_class" else ""),
                "source_record_uid": str(trait.get("source_record_uid") or ""),
            }
        )
        component_rows = [row for row in components if row.get("trait_uid") == trait_uid]
        for rank, component in enumerate(component_rows, start=2):
            component_uid = str(component.get("component_uid") or "")
            component_type = str(component.get("component_entity_type") or "unknown")
            if not component_uid and component_type == "unknown":
                match_basis = "trait_component_unresolved"
                risk = "high"
            elif component_type == "chemical":
                match_basis = "trait_component_chemical"
                risk = "medium" if feature_type == "ratio" else "low"
            else:
                match_basis = "trait_component_class"
                risk = "high"
            identity_candidates.append(
                {
                    "input_feature_uid": feature_uid,
                    "candidate_uid": component_uid,
                    "candidate_entity_type": component_type,
                    "candidate_name": str(component.get("component_name") or ""),
                    "match_basis": match_basis,
                    "score": 60.0 if component_uid else 0.0,
                    "margin": 0.0,
                    "rank": rank,
                    "false_match_risk": risk,
                    "blocking_reason": "ratio_component_not_abundance_seed" if feature_type == "ratio" else ("" if component_uid else "missing_component_identity"),
                    "source_record_uid": str(component.get("source_record_uid") or ""),
                }
            )
        status, accepted_uid, accepted_type, rule, confidence, review = decision_for_feature(feature_type, trait_uid, "", "")
        identity_decisions.append(
            {
                "decision_uid": stable_uid("decision", feature_uid, status, accepted_uid, rule),
                "input_feature_uid": feature_uid,
                "decision_status": status,
                "accepted_entity_uid": accepted_uid,
                "accepted_entity_type": accepted_type,
                "decision_rule": rule,
                "decision_confidence": confidence,
                "review_status": review,
                "reviewer": "",
                "reviewed_at": "",
                "decision_notes": "Trait inputs are retained as traits and never promoted to exact abundance metabolites.",
            }
        )

    for class_row in classes:
        class_name = str(class_row.get("class_name") or "")
        if not class_name:
            continue
        class_uid = str(class_row.get("class_uid") or "")
        feature_uid = stable_uid("input_feature", "class", class_uid)
        input_features.append(
            {
                "input_feature_uid": feature_uid,
                "input_row_id": class_uid,
                "feature_label": class_name,
                "feature_type": "lipid_class",
                "effect_value": None,
                "effect_label": "",
                "direction": "",
                "pvalue": None,
                "padj": None,
                "sample_context_uid": "",
                "raw_input_json": json_dumps(class_row),
            }
        )
        identity_candidates.append(
            {
                "input_feature_uid": feature_uid,
                "candidate_uid": class_uid,
                "candidate_entity_type": "class",
                "candidate_name": class_name,
                "match_basis": "class_name",
                "score": 80.0,
                "margin": 80.0,
                "rank": 1,
                "false_match_risk": "high",
                "blocking_reason": "class_name_not_exact_identity",
                "source_record_uid": str(class_row.get("source_record_uid") or ""),
            }
        )
        status, accepted_uid, accepted_type, rule, confidence, review = decision_for_feature("lipid_class", "", class_uid, "")
        identity_decisions.append(
            {
                "decision_uid": stable_uid("decision", feature_uid, status, accepted_uid, rule),
                "input_feature_uid": feature_uid,
                "decision_status": status,
                "accepted_entity_uid": accepted_uid,
                "accepted_entity_type": accepted_type,
                "decision_rule": rule,
                "decision_confidence": confidence,
                "review_status": review,
                "reviewer": "",
                "reviewed_at": "",
                "decision_notes": "Class-level surface is retained as class evidence and cannot become exact chemical identity.",
            }
        )

    for chemical_uid, chemical in chemical_by_uid.items():
        label = str(chemical.get("canonical_name") or chemical_uid)
        feature_uid = stable_uid("input_feature", "chemical", chemical_uid)
        input_features.append(
            {
                "input_feature_uid": feature_uid,
                "input_row_id": chemical_uid,
                "feature_label": label,
                "feature_type": "direct_metabolite",
                "effect_value": None,
                "effect_label": "",
                "direction": "",
                "pvalue": None,
                "padj": None,
                "sample_context_uid": "",
                "raw_input_json": json_dumps(chemical),
            }
        )
        identity_candidates.append(
            {
                "input_feature_uid": feature_uid,
                "candidate_uid": chemical_uid,
                "candidate_entity_type": "chemical",
                "candidate_name": label,
                "match_basis": "stable_source_entity",
                "score": 100.0 if chemical.get("identity_status") == "canonical" else 0.0,
                "margin": 100.0 if chemical.get("identity_status") == "canonical" else 0.0,
                "rank": 1,
                "false_match_risk": "low" if chemical.get("identity_status") == "canonical" else "high",
                "blocking_reason": "" if chemical.get("identity_status") == "canonical" else "high_risk_name",
                "source_record_uid": str(chemical.get("source_record_uid") or ""),
            }
        )
        status = "accepted_exact" if chemical.get("identity_status") == "canonical" else "rejected"
        identity_decisions.append(
            {
                "decision_uid": stable_uid("decision", feature_uid, status, chemical_uid if status == "accepted_exact" else "", "stable_source_entity" if status == "accepted_exact" else "high_risk_name"),
                "input_feature_uid": feature_uid,
                "decision_status": status,
                "accepted_entity_uid": chemical_uid if status == "accepted_exact" else "",
                "accepted_entity_type": "chemical" if status == "accepted_exact" else "",
                "decision_rule": "stable_source_entity" if status == "accepted_exact" else "high_risk_name",
                "decision_confidence": 0.9 if status == "accepted_exact" else 0.0,
                "review_status": "auto" if status == "accepted_exact" else "needs_review",
                "reviewer": "",
                "reviewed_at": "",
                "decision_notes": "",
            }
        )

    root = output_dir / "identity_resolution_store"
    return [
        write_table(root / "input_features.parquet", input_features, "input_features"),
        write_table(root / "identity_candidates.parquet", identity_candidates, "identity_candidates"),
        write_table(root / "identity_decisions.parquet", identity_decisions, "identity_decisions"),
    ]


def build_relation_store(workspace: Path, release_id: str, normalized_dir: Path, output_dir: Path) -> list[dict[str, Any]]:
    reactions_path = normalized_dir / "reactions.parquet"
    participants_path = normalized_dir / "reaction_participants.parquet"
    pathways_path = normalized_dir / "pathways.parquet"
    xrefs_path = normalized_dir / "metabolite_xrefs.parquet"
    reactome_pe_path = workspace / "raw_lake" / "reactome" / release_id / "ChEBI2Reactome_PE_Reactions.txt"
    reactome_exporter_path = workspace / "raw_lake" / "reactome" / release_id / "reactome_reaction_exporter.txt"
    reaction_pmids_path = workspace / "raw_lake" / "reactome" / release_id / "ReactionPMIDS.txt"
    rhea_tsv_path = workspace / "raw_lake" / "rhea" / release_id / "rhea-tsv.tar.gz"
    reactions = read_rows(reactions_path)
    participants = read_rows(participants_path)
    pathways = read_rows(pathways_path)
    xrefs = read_rows(xrefs_path)
    modules_by_pathway: dict[str, str] = {}
    reaction_by_external = {str(row.get("primary_external_id") or ""): row for row in reactions}
    reaction_by_external_base = {strip_reactome_version(row.get("primary_external_id")): row for row in reactions}
    reaction_to_pathway = {str(row.get("reaction_uid") or ""): str(row.get("pathway_uid") or "") for row in reactions}
    chebi_to_chemical: dict[str, str] = {}
    for row in xrefs:
        chemical_uid = str(row.get("metabolite_uid") or "")
        if not chemical_uid:
            continue
        if str(row.get("xref_source") or "").casefold() == "chebi":
            key = chebi_key(row.get("xref_id") or row.get("xref_key"))
            if key:
                chebi_to_chemical.setdefault(key, chemical_uid)
    reaction_entities = []
    participants_v2 = []
    reaction_equations = []
    reaction_side_participants = []
    reaction_xrefs = []
    enzyme_reaction_links = []
    reaction_publication_links = []
    pathway_modules = []
    module_members = []
    participant_seen: set[tuple[str, str, str, str]] = set()
    side_seen: set[tuple[str, str, str, str]] = set()
    equation_seen: set[tuple[str, str]] = set()
    xref_seen: set[tuple[str, str, str]] = set()
    enzyme_seen: set[tuple[str, str, str]] = set()
    publication_seen: set[tuple[str, str]] = set()
    module_member_seen: set[tuple[str, str, str, str]] = set()
    rhea_to_reaction_uids: dict[str, list[str]] = {}

    for row in reactions:
        reaction_uid = str(row.get("reaction_uid") or "")
        if not reaction_uid:
            continue
        src_uid = source_uid("normalized_store.reactions", release_id, reaction_uid)
        reaction_entities.append(
            {
                "reaction_uid": reaction_uid,
                "reaction_name": str(row.get("name") or ""),
                "reaction_source": str(row.get("source_name") or ""),
                "equation": "",
                "compartment_uid": "",
                "source_record_uid": src_uid,
            }
        )
        pathway_uid = str(row.get("pathway_uid") or "")
        if pathway_uid:
            module_uid = modules_by_pathway.setdefault(pathway_uid, stable_uid("module", pathway_uid))
            module_members.append(
                {
                    "module_uid": module_uid,
                    "member_uid": reaction_uid,
                    "member_type": "reaction",
                    "member_role": "core",
                    "source_record_uid": src_uid,
                }
            )

    for row in pathways:
        pathway_uid = str(row.get("pathway_uid") or "")
        if not pathway_uid:
            continue
        module_uid = modules_by_pathway.setdefault(pathway_uid, stable_uid("module", pathway_uid))
        pathway_modules.append(
            {
                "module_uid": module_uid,
                "module_name": str(row.get("name") or ""),
                "parent_pathway_uid": pathway_uid,
                "module_type": "pathway_module",
                "definition": "Initial v2 module projected from curated pathway; refine into mechanism modules during curation.",
                "source_record_uid": source_uid("normalized_store.pathways", release_id, pathway_uid),
            }
        )

    for row in participants:
        participant_uid = str(row.get("participant_uid") or "")
        reaction_uid = str(row.get("reaction_uid") or "")
        if not participant_uid or not reaction_uid:
            continue
        role = str(row.get("role") or "participant").strip().casefold()
        participant_role = "substrate" if role in {"input", "left", "substrate"} else ("product" if role in {"output", "right", "product"} else role)
        src_uid = source_uid("normalized_store.reaction_participants", release_id, str(row.get("edge_uid") or f"{reaction_uid}|{participant_uid}"))
        if str(row.get("participant_type") or "") == "metabolite":
            add_unique(
                participants_v2,
                participant_seen,
                (reaction_uid, participant_uid, participant_role, src_uid),
                {
                    "reaction_uid": reaction_uid,
                    "chemical_uid": participant_uid,
                    "chemical_source_id": str(row.get("participant_external_id") or ""),
                    "physical_entity_id": "",
                    "physical_entity_name": "",
                    "participant_role": participant_role,
                    "stoichiometry": None,
                    "compartment_uid": "",
                    "directionality": "unknown",
                    "relation_source": "normalized_store.reaction_participants",
                    "semantic_source_uri": "",
                    "source_record_uid": src_uid,
                },
            )
            pathway_uid = reaction_to_pathway.get(reaction_uid, "")
            module_uid = modules_by_pathway.get(pathway_uid, "")
            if module_uid:
                add_unique(
                    module_members,
                    module_member_seen,
                    (module_uid, participant_uid, "chemical", src_uid),
                    {
                        "module_uid": module_uid,
                        "member_uid": participant_uid,
                        "member_type": "chemical",
                        "member_role": "supporting",
                        "source_record_uid": src_uid,
                    }
                )

    if reactome_pe_path.exists():
        with reactome_pe_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for raw in reader:
                if len(raw) < 8:
                    continue
                chebi_id, physical_entity_id, physical_entity_name, reaction_external_id, _url, _reaction_name, evidence_code, species = raw[:8]
                if species != "Homo sapiens":
                    continue
                chemical_uid = chebi_to_chemical.get(chebi_key(chebi_id), "")
                reaction = reaction_by_external.get(reaction_external_id, {})
                reaction_uid = str(reaction.get("reaction_uid") or "")
                if not chemical_uid or not reaction_uid:
                    continue
                src_uid = source_uid(
                    "raw_lake.reactome.ChEBI2Reactome_PE_Reactions",
                    release_id,
                    f"{chebi_id}|{physical_entity_id}|{reaction_external_id}",
                )
                add_unique(
                    participants_v2,
                    participant_seen,
                    (reaction_uid, chemical_uid, "participant", src_uid),
                    {
                        "reaction_uid": reaction_uid,
                        "chemical_uid": chemical_uid,
                        "chemical_source_id": f"CHEBI:{chebi_key(chebi_id)}",
                        "physical_entity_id": physical_entity_id,
                        "physical_entity_name": physical_entity_name,
                        "participant_role": "participant",
                        "stoichiometry": None,
                        "compartment_uid": "",
                        "directionality": "unknown",
                        "relation_source": "reactome_pe_participation",
                        "semantic_source_uri": "",
                        "source_record_uid": src_uid,
                    },
                )
                pathway_uid = reaction_to_pathway.get(reaction_uid, "")
                module_uid = modules_by_pathway.get(pathway_uid, "")
                if module_uid:
                    add_unique(
                        module_members,
                        module_member_seen,
                        (module_uid, chemical_uid, "chemical", src_uid),
                        {
                            "module_uid": module_uid,
                            "member_uid": chemical_uid,
                            "member_type": "chemical",
                            "member_role": "supporting",
                            "source_record_uid": src_uid,
                        }
                    )

    direction_rows = tsv_records(read_tsv_from_tar(rhea_tsv_path, "tsv/rhea-directions.tsv"))
    direction_by_master = {
        str(row.get("RHEA_ID_MASTER") or ""): {
            "lr": str(row.get("RHEA_ID_LR") or ""),
            "rl": str(row.get("RHEA_ID_RL") or ""),
            "bi": str(row.get("RHEA_ID_BI") or ""),
        }
        for row in direction_rows
        if row.get("RHEA_ID_MASTER")
    }
    direction_by_rhea: dict[str, dict[str, str]] = {}
    for master, mapping in direction_by_master.items():
        if master:
            direction_by_rhea[master] = {**mapping, "direction": "UN", "master": master}
        for direction, rhea_id in (("LR", mapping.get("lr", "")), ("RL", mapping.get("rl", "")), ("BI", mapping.get("bi", ""))):
            if rhea_id:
                direction_by_rhea[rhea_id] = {**mapping, "direction": direction, "master": master}

    chebi_smiles_rows = read_tsv_from_tar(rhea_tsv_path, "tsv/rhea-chebi-smiles.tsv")
    smiles_to_chebi: dict[str, list[str]] = {}
    for row in chebi_smiles_rows:
        if len(row) < 2:
            continue
        chebi_id, smiles = row[0].strip(), row[1].strip()
        if smiles:
            smiles_to_chebi.setdefault(smiles, []).append(chebi_id)

    rhea_smiles: dict[str, str] = {}
    for row in read_tsv_from_tar(rhea_tsv_path, "tsv/rhea-reaction-smiles.tsv"):
        if len(row) >= 2 and row[0].strip().isdigit():
            rhea_smiles[row[0].strip()] = row[1].strip()

    for row in tsv_records(read_tsv_from_tar(rhea_tsv_path, "tsv/rhea2reactome.tsv")):
        rhea_id = str(row.get("RHEA_ID") or "").strip()
        reactome_id = strip_reactome_version(row.get("ID"))
        reaction = reaction_by_external_base.get(reactome_id, {})
        reaction_uid = str(reaction.get("reaction_uid") or "")
        if not rhea_id or not reaction_uid:
            continue
        direction = str(row.get("DIRECTION") or direction_by_rhea.get(rhea_id, {}).get("direction") or "")
        master_id = str(row.get("MASTER_ID") or direction_by_rhea.get(rhea_id, {}).get("master") or "")
        mapping = direction_by_master.get(master_id, {})
        src_uid = source_uid("raw_lake.rhea.rhea2reactome", release_id, f"{rhea_id}|{reactome_id}")
        add_unique(
            reaction_xrefs,
            xref_seen,
            (reaction_uid, "RHEA", rhea_id),
            {
                "reaction_uid": reaction_uid,
                "xref_source": "RHEA",
                "xref_id": rhea_id,
                "direction": direction,
                "source_record_uid": src_uid,
            },
        )
        rhea_to_reaction_uids.setdefault(rhea_id, []).append(reaction_uid)
        smiles = rhea_smiles.get(rhea_id, "")
        if not smiles:
            continue
        add_unique(
            reaction_equations,
            equation_seen,
            (reaction_uid, rhea_id),
            {
                "reaction_uid": reaction_uid,
                "source_reaction_id": f"RHEA:{rhea_id}",
                "equation_source": "Rhea reaction SMILES",
                "directionality": direction or "unknown",
                "left_to_right_uid": mapping.get("lr", ""),
                "right_to_left_uid": mapping.get("rl", ""),
                "bidirectional_uid": mapping.get("bi", ""),
                "equation_text": smiles,
                "reaction_smiles": smiles,
                "source_record_uid": source_uid("raw_lake.rhea.rhea-reaction-smiles", release_id, rhea_id),
            },
        )
        left_smiles, right_smiles = split_reaction_smiles(smiles)
        for side, role, parts in (("left", "substrate", left_smiles), ("right", "product", right_smiles)):
            for part in parts:
                chebi_ids = smiles_to_chebi.get(part, [])
                for chebi_id in chebi_ids:
                    chemical_uid = chebi_to_chemical.get(chebi_key(chebi_id), "")
                    if not chemical_uid:
                        continue
                    side_src_uid = source_uid("raw_lake.rhea.rhea-reaction-smiles.participant", release_id, f"{rhea_id}|{side}|{chebi_id}|{part}")
                    add_unique(
                        reaction_side_participants,
                        side_seen,
                        (reaction_uid, chemical_uid, side, side_src_uid),
                        {
                            "reaction_uid": reaction_uid,
                            "chemical_uid": chemical_uid,
                            "chemical_source_id": f"CHEBI:{chebi_key(chebi_id)}",
                            "side": side,
                            "participant_role": role,
                            "stoichiometry": None,
                            "directionality": direction or "unknown",
                            "relation_source": "rhea_smiles_exact_component",
                            "physical_entity_id": f"RHEA:{rhea_id}:{side}",
                            "physical_entity_name": part,
                            "semantic_source_uri": "",
                            "source_record_uid": side_src_uid,
                        },
                    )
                    add_unique(
                        participants_v2,
                        participant_seen,
                        (reaction_uid, chemical_uid, role, side_src_uid),
                        {
                            "reaction_uid": reaction_uid,
                            "chemical_uid": chemical_uid,
                            "chemical_source_id": f"CHEBI:{chebi_key(chebi_id)}",
                            "physical_entity_id": f"RHEA:{rhea_id}:{side}",
                            "physical_entity_name": part,
                            "participant_role": role,
                            "stoichiometry": None,
                            "compartment_uid": "",
                            "directionality": direction or "unknown",
                            "relation_source": "rhea_smiles_exact_component",
                            "semantic_source_uri": "",
                            "source_record_uid": side_src_uid,
                        },
                    )

    semantic_paths = semantic_source_paths(workspace, release_id)
    for semantic_path in semantic_paths["rhea_biopax"]:
        if semantic_path.exists():
            if "rhea.rdf" in semantic_path.name.casefold():
                parse_rhea_rdf_semantic_sides(
                    semantic_path,
                    release_id=release_id,
                    chebi_to_chemical=chebi_to_chemical,
                    rhea_to_reaction_uids=rhea_to_reaction_uids,
                    direction_by_rhea=direction_by_rhea,
                    reaction_equations=reaction_equations,
                    reaction_side_participants=reaction_side_participants,
                    participants_v2=participants_v2,
                    reaction_xrefs=reaction_xrefs,
                    equation_seen=equation_seen,
                    side_seen=side_seen,
                    participant_seen=participant_seen,
                    xref_seen=xref_seen,
                )
            else:
                parse_biopax_semantic_sides(
                    semantic_path,
                    release_id=release_id,
                    chebi_to_chemical=chebi_to_chemical,
                    reaction_by_external_base=reaction_by_external_base,
                    rhea_to_reaction_uids=rhea_to_reaction_uids,
                    direction_by_rhea=direction_by_rhea,
                    reaction_equations=reaction_equations,
                    reaction_side_participants=reaction_side_participants,
                    participants_v2=participants_v2,
                    reaction_xrefs=reaction_xrefs,
                    equation_seen=equation_seen,
                    side_seen=side_seen,
                    participant_seen=participant_seen,
                    xref_seen=xref_seen,
                )
            break

    for semantic_path in semantic_paths["reactome_sbml"]:
        if semantic_path.exists():
            parse_sbml_semantic_sides(
                semantic_path,
                release_id=release_id,
                chebi_to_chemical=chebi_to_chemical,
                reaction_by_external_base=reaction_by_external_base,
                reaction_equations=reaction_equations,
                reaction_side_participants=reaction_side_participants,
                participants_v2=participants_v2,
                equation_seen=equation_seen,
                side_seen=side_seen,
                participant_seen=participant_seen,
            )
            break

    for semantic_path in semantic_paths["reactome_biopax"]:
        if semantic_path.exists():
            parse_biopax_semantic_sides(
                semantic_path,
                release_id=release_id,
                chebi_to_chemical=chebi_to_chemical,
                reaction_by_external_base=reaction_by_external_base,
                rhea_to_reaction_uids=rhea_to_reaction_uids,
                direction_by_rhea=direction_by_rhea,
                reaction_equations=reaction_equations,
                reaction_side_participants=reaction_side_participants,
                participants_v2=participants_v2,
                reaction_xrefs=reaction_xrefs,
                equation_seen=equation_seen,
                side_seen=side_seen,
                participant_seen=participant_seen,
                xref_seen=xref_seen,
            )
            break

    for row in read_tsv_from_tar(rhea_tsv_path, "tsv/rhea2uniprot_sprot.tsv"):
        if not row or str(row[0]).startswith("RHEA_ID"):
            continue
        if len(row) < 4:
            continue
        rhea_id, direction, _master_id, protein_id = [cell.strip() for cell in row[:4]]
        for reaction_uid in rhea_to_reaction_uids.get(rhea_id, []):
            if not reaction_uid or not protein_id:
                continue
            src_uid = source_uid("raw_lake.rhea.rhea2uniprot_sprot", release_id, f"{rhea_id}|{protein_id}")
            add_unique(
                enzyme_reaction_links,
                enzyme_seen,
                (reaction_uid, protein_id, src_uid),
                {
                    "reaction_uid": reaction_uid,
                    "protein_uid": "",
                    "protein_source_id": protein_id,
                    "gene_uid": "",
                    "enzyme_role": "enzyme",
                    "source_record_uid": src_uid,
                },
            )

    if reactome_exporter_path.exists():
        with reactome_exporter_path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                reaction_id = strip_reactome_version(row.get("reaction_id"))
                reaction = reaction_by_external_base.get(reaction_id, {})
                reaction_uid = str(reaction.get("reaction_uid") or "")
                uniprot = str(row.get("uniprot_acc") or "").strip()
                if not reaction_uid or not uniprot:
                    continue
                role = str(row.get("role_in_reaction") or "enzyme").strip() or "enzyme"
                src_uid = source_uid("raw_lake.reactome.reactome_reaction_exporter", release_id, f"{reaction_id}|{uniprot}|{role}")
                add_unique(
                    enzyme_reaction_links,
                    enzyme_seen,
                    (reaction_uid, uniprot, src_uid),
                    {
                        "reaction_uid": reaction_uid,
                        "protein_uid": "",
                        "protein_source_id": uniprot,
                        "gene_uid": "",
                        "enzyme_role": role,
                        "source_record_uid": src_uid,
                    },
                )

    if reaction_pmids_path.exists():
        with reaction_pmids_path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) < 2:
                    continue
                reaction_id = strip_reactome_version(row[0])
                reaction = reaction_by_external_base.get(reaction_id, {})
                reaction_uid = str(reaction.get("reaction_uid") or "")
                pmid = str(row[1] or "").strip()
                if not reaction_uid or not pmid:
                    continue
                src_uid = source_uid("raw_lake.reactome.ReactionPMIDS", release_id, f"{reaction_id}|{pmid}")
                add_unique(
                    reaction_publication_links,
                    publication_seen,
                    (reaction_uid, pmid),
                    {"reaction_uid": reaction_uid, "pmid": pmid, "source_record_uid": src_uid},
                )

    root = output_dir / "relation_store"
    return [
        write_table(root / "reaction_entities.parquet", reaction_entities, "reaction_entities"),
        write_table(root / "reaction_participants_v2.parquet", participants_v2, "reaction_participants_v2"),
        write_table(root / "reaction_catalysts.parquet", [], "reaction_catalysts"),
        write_table(root / "reaction_equations.parquet", reaction_equations, "reaction_equations"),
        write_table(root / "reaction_side_participants.parquet", reaction_side_participants, "reaction_side_participants"),
        write_table(root / "reaction_xrefs.parquet", reaction_xrefs, "reaction_xrefs"),
        write_table(root / "enzyme_reaction_links.parquet", enzyme_reaction_links, "enzyme_reaction_links"),
        write_table(root / "reaction_publication_links.parquet", reaction_publication_links, "reaction_publication_links"),
        write_table(root / "pathway_modules.parquet", pathway_modules, "pathway_modules"),
        write_table(root / "module_members.parquet", module_members, "module_members"),
    ]


def support_status_for_class(value: str) -> str:
    if value == "conflict_candidate":
        return "contradict"
    if value == "unresolved":
        return "uncertain"
    if value in {"confirm", "support_direction"}:
        return "support"
    return "background"


def evidence_class_for_support(value: str) -> str:
    if value in {"confirm", "support_direction"}:
        return "omics_association"
    if value == "conflict_candidate":
        return "omics_association"
    return "literature_background"


def build_evidence_store(release_id: str, normalized_dir: Path, literature_dir: Path, output_dir: Path) -> list[dict[str, Any]]:
    sentences_path = normalized_dir / "sentences.parquet"
    support_path = literature_dir / "literature_edge_support.parquet"
    sentences = read_rows(sentences_path)
    supports = read_rows(support_path)
    evidence_sentences = []
    evidence_contexts = [
        {
            "context_uid": "context_unspecified",
            "cancer_type_uid": "",
            "tissue_uid": "",
            "cell_type_uid": "",
            "cell_state_uid": "",
            "species": "",
            "model_system": "",
            "comparison": "",
            "context_text": "Unspecified context; cannot support context-specific mechanism without further curation.",
        }
    ]
    evidence_assertions = []
    for row in sentences:
        sentence_uid = str(row.get("sentence_uid") or "")
        if not sentence_uid:
            continue
        evidence_sentences.append(
            {
                "sentence_uid": sentence_uid,
                "pmid": str(row.get("pmid") or ""),
                "pmcid": str(row.get("pmcid") or ""),
                "section": str(row.get("section") or ""),
                "sentence_text": str(row.get("sentence_text") or ""),
                "sentence_hash": str(row.get("text_hash") or ""),
                "source_record_uid": source_uid("normalized_store.sentences", release_id, sentence_uid),
            }
        )
    for row in supports:
        support_uid = str(row.get("support_uid") or "")
        sentence_uids = row.get("sentence_uids") or [""]
        support_class = str(row.get("support_class") or "")
        confidence = {
            "p_literature": row.get("p_literature"),
            "raw_score_max": row.get("raw_score_max"),
            "calibrated_prob_max": row.get("calibrated_prob_max"),
            "support_class": support_class,
            "support_status_set": row.get("support_status_set") or [],
        }
        for sentence_uid in sentence_uids or [""]:
            evidence_assertions.append(
                {
                    "assertion_uid": stable_uid("assertion", support_uid, sentence_uid),
                    "subject_uid": str(row.get("subject_uid") or ""),
                    "subject_type": str(row.get("subject_type") or ""),
                    "predicate": str(row.get("predicate") or ""),
                    "object_uid": str(row.get("object_uid") or ""),
                    "object_type": str(row.get("object_type") or ""),
                    "polarity": ",".join(row.get("polarity_set") or []),
                    "direction": "mixed",
                    "evidence_class": evidence_class_for_support(support_class),
                    "support_status": support_status_for_class(support_class),
                    "context_uid": "context_unspecified",
                    "method": "text_mined",
                    "source_record_uid": source_uid("literature_evidence.literature_edge_support", release_id, support_uid),
                    "sentence_uid": str(sentence_uid or ""),
                    "confidence_components_json": json_dumps(confidence),
                }
            )

    root = output_dir / "evidence_store"
    return [
        write_table(root / "evidence_contexts.parquet", evidence_contexts, "evidence_contexts"),
        write_table(root / "evidence_sentences.parquet", evidence_sentences, "evidence_sentences"),
        write_table(root / "evidence_assertions.parquet", evidence_assertions, "evidence_assertions"),
    ]


def build_analysis_view(output_dir: Path) -> list[dict[str, Any]]:
    entity_root = output_dir / "entity_store"
    identity_root = output_dir / "identity_resolution_store"
    relation_root = output_dir / "relation_store"
    evidence_root = output_dir / "evidence_store"
    chemicals = {row["chemical_uid"]: row for row in read_rows(entity_root / "chemical_entities.parquet")}
    reactions = {row["reaction_uid"]: row for row in read_rows(relation_root / "reaction_entities.parquet")}
    participants = read_rows(relation_root / "reaction_participants_v2.parquet")
    module_members = read_rows(relation_root / "module_members.parquet")
    enzyme_links = read_rows(relation_root / "enzyme_reaction_links.parquet")
    publication_links = read_rows(relation_root / "reaction_publication_links.parquet")
    assertions = read_rows(evidence_root / "evidence_assertions.parquet")
    decisions = read_rows(identity_root / "identity_decisions.parquet")
    accepted_exact_decisions: dict[str, list[str]] = {}
    for row in decisions:
        if row.get("decision_status") == "accepted_exact" and row.get("accepted_entity_type") == "chemical":
            accepted_exact_decisions.setdefault(str(row.get("accepted_entity_uid") or ""), []).append(
                str(row.get("decision_uid") or row.get("input_feature_uid") or "")
            )
    assertions_by_subject: dict[str, list[str]] = {}
    for row in assertions:
        if row.get("support_status") == "support" and row.get("evidence_class") in {"curated_database", "direct_assay", "omics_association"}:
            assertions_by_subject.setdefault(str(row.get("subject_uid") or ""), []).append(str(row.get("assertion_uid") or ""))
    enzymes_by_reaction: dict[str, list[dict[str, Any]]] = {}
    for row in enzyme_links:
        enzymes_by_reaction.setdefault(str(row.get("reaction_uid") or ""), []).append(row)
    pmids_by_reaction: dict[str, list[str]] = {}
    for row in publication_links:
        pmid = str(row.get("pmid") or "")
        if pmid:
            pmids_by_reaction.setdefault(str(row.get("reaction_uid") or ""), []).append(pmid)
    fact_candidates = []
    facts = []
    for row in participants:
        chemical_uid = str(row.get("chemical_uid") or "")
        reaction_uid = str(row.get("reaction_uid") or "")
        chemical = chemicals.get(chemical_uid, {})
        reaction = reactions.get(reaction_uid, {})
        relation_source = str(row.get("relation_source") or "curated_reaction_participant")
        semantic_source_uri = str(row.get("semantic_source_uri") or "")
        semantic_sources = {
            "rhea_rdf_explicit_side",
            "rhea_biopax_explicit_side",
            "reactome_sbml_species_reference",
            "reactome_biopax_conversion",
        }
        blocking_reasons = []
        if relation_source in semantic_sources and not chemical_uid:
            blocking_reasons.append("unresolved_semantic_participant")
        if not chemical:
            blocking_reasons.append("chemical_not_in_entity_store")
        if not accepted_exact_decisions.get(chemical_uid):
            blocking_reasons.append("identity_decision_not_accepted_exact")
        if not row.get("source_record_uid"):
            blocking_reasons.append("missing_relation_source_record")
        if not reaction:
            blocking_reasons.append("reaction_not_in_relation_store")
        role = str(row.get("participant_role") or "participant")
        direction = str(row.get("directionality") or "unknown")
        if role in {"substrate", "product"}:
            candidate_type = "reaction_side_participant"
            claim_scope = "bidirectional_reaction_fact" if direction == "BI" else "directional_reaction_fact"
            if relation_source in semantic_sources:
                direction_gate = "explicit_semantic_side_assignment"
                boundary = (
                    "Exact chemical side is derived from explicit RDF/BioPAX/SBML reaction semantics. This supports substrate/product wording for the reaction, "
                    "but still does not prove sample-level flux or pathway activation without input direction and module consistency."
                )
                readiness_score = 0.96
            else:
                direction_gate = "side_assignment_from_rhea_smiles"
                boundary = (
                    "Exact chemical side is derived from a Rhea directed equation. This supports substrate/product wording for the reaction, "
                    "but still does not prove sample-level flux or pathway activation without input direction and module consistency."
                )
                readiness_score = 0.92
        elif role == "participant":
            candidate_type = "reaction_participation"
            claim_scope = "role_unknown_reaction_fact"
            direction_gate = "role_unknown_allowed_with_boundary"
            boundary = (
                "Exact chemical participates in a curated reaction, but substrate/product role is unknown. "
                "Use only participation wording."
            )
            readiness_score = 0.78
        else:
            candidate_type = "reaction_participation"
            claim_scope = "exact_reaction_fact"
            direction_gate = "source_role_preserved"
            boundary = (
                "Exact chemical has a curated reaction role from source data. This supports reaction-level wording, "
                "not pathway activation or flux direction by itself."
            )
            readiness_score = 0.85
        candidate_uid = stable_uid("fact_candidate", "reaction_participation", reaction_uid, chemical_uid, row.get("source_record_uid"))
        candidate = {
            "candidate_fact_uid": candidate_uid,
            "candidate_fact_type": candidate_type,
            "subject_uid": chemical_uid,
            "subject_type": "chemical",
            "predicate": "participates_in_reaction",
            "object_uid": reaction_uid,
            "object_type": "reaction",
            "relation_source": relation_source,
            "source_record_uid": str(row.get("source_record_uid") or ""),
            "context_uid": "context_unspecified",
            "direction": direction,
            "role": role,
            "evidence_assertion_uids": assertions_by_subject.get(chemical_uid, [])[:20],
            "blocking_reasons": blocking_reasons,
            "readiness_score": readiness_score if not blocking_reasons else 0.0,
            "metadata_json": json_dumps(
                {
                    "chemical_source_id": row.get("chemical_source_id"),
                    "physical_entity_id": row.get("physical_entity_id"),
                    "physical_entity_name": row.get("physical_entity_name"),
                    "relation_source": relation_source,
                    "semantic_source_uri": semantic_source_uri,
                    "reaction_name": reaction.get("reaction_name", ""),
                    "enzyme_protein_source_ids": [link.get("protein_source_id", "") for link in enzymes_by_reaction.get(reaction_uid, [])[:10]],
                    "pmids": pmids_by_reaction.get(reaction_uid, [])[:10],
                    "entity_granularity": chemical.get("entity_granularity", ""),
                    "fact_gate": "identity+relation_source+reaction_store",
                }
            ),
        }
        fact_candidates.append(candidate)
        if blocking_reasons:
            continue
        facts.append(
            {
                "fact_uid": stable_uid("fact", "reaction_participation", reaction_uid, chemical_uid, row.get("participant_role"), row.get("source_record_uid")),
                "fact_type": candidate_type,
                "subject_uid": chemical_uid,
                "subject_type": "chemical",
                "predicate": "participates_in_reaction",
                "object_uid": reaction_uid,
                "object_type": "reaction",
                "context_uid": "context_unspecified",
                "direction": direction,
                "role": role,
                "readiness_tier": "high",
                "allowed_claim_scope": claim_scope,
                "blocking_reasons": [],
                "supporting_assertion_uids": assertions_by_subject.get(chemical_uid, [])[:20],
                "source_record_uids": [str(row.get("source_record_uid") or "")],
                "identity_decision_uids": accepted_exact_decisions.get(chemical_uid, [])[:5],
                "confidence_components_json": json_dumps(
                    {
                        "identity_gate": "accepted_exact",
                        "relation_gate": relation_source,
                        "context_gate": "unspecified_context",
                        "direction_gate": direction_gate,
                        "readiness_score": readiness_score,
                    }
                ),
                "boundary_text": boundary,
                "metadata_json": json_dumps(
                    {
                        "participant_role": row.get("participant_role"),
                        "chemical_source_id": row.get("chemical_source_id"),
                        "physical_entity_id": row.get("physical_entity_id"),
                        "physical_entity_name": row.get("physical_entity_name"),
                        "relation_source": relation_source,
                        "semantic_source_uri": semantic_source_uri,
                        "reaction_name": reaction.get("reaction_name", ""),
                        "enzyme_protein_source_ids": [link.get("protein_source_id", "") for link in enzymes_by_reaction.get(reaction_uid, [])[:10]],
                        "pmids": pmids_by_reaction.get(reaction_uid, [])[:10],
                        "entity_granularity": chemical.get("entity_granularity"),
                    }
                ),
            }
        )
    for row in enzyme_links:
        reaction_uid = str(row.get("reaction_uid") or "")
        reaction = reactions.get(reaction_uid, {})
        fact_candidates.append(
            {
                "candidate_fact_uid": stable_uid("fact_candidate", "enzyme_reaction_link", reaction_uid, row.get("protein_source_id"), row.get("source_record_uid")),
                "candidate_fact_type": "enzyme_reaction_link",
                "subject_uid": reaction_uid,
                "subject_type": "reaction",
                "predicate": "has_enzyme_or_protein_role",
                "object_uid": str(row.get("protein_uid") or row.get("protein_source_id") or ""),
                "object_type": "protein",
                "relation_source": "Rhea/Reactome enzyme reaction link",
                "source_record_uid": str(row.get("source_record_uid") or ""),
                "context_uid": "context_unspecified",
                "direction": "unknown",
                "role": str(row.get("enzyme_role") or "enzyme"),
                "evidence_assertion_uids": [],
                "blocking_reasons": [] if reaction else ["reaction_not_in_relation_store"],
                "readiness_score": 0.7 if reaction else 0.0,
                "metadata_json": json_dumps({"reaction_name": reaction.get("reaction_name", ""), "claim_scope": "enzyme_reaction_fact"}),
            }
        )
    for row in module_members:
        if row.get("member_type") != "chemical":
            continue
        chemical_uid = str(row.get("member_uid") or "")
        fact_candidates.append(
            {
                "candidate_fact_uid": stable_uid("fact_candidate", "module_membership", row.get("module_uid"), chemical_uid, row.get("source_record_uid")),
                "candidate_fact_type": "module_membership",
                "subject_uid": chemical_uid,
                "subject_type": "chemical",
                "predicate": "member_of_module",
                "object_uid": str(row.get("module_uid") or ""),
                "object_type": "module",
                "relation_source": "pathway_projected_module",
                "source_record_uid": str(row.get("source_record_uid") or ""),
                "context_uid": "context_unspecified",
                "direction": "unknown",
                "role": str(row.get("member_role") or "supporting"),
                "evidence_assertion_uids": assertions_by_subject.get(chemical_uid, [])[:20],
                "blocking_reasons": ["pathway_projected_module_candidate_only"],
                "readiness_score": 0.25,
                "metadata_json": json_dumps({"fact_gate": "candidate_only_until_curated_mechanism_module"}),
            }
        )
    for row in assertions:
        support_status = str(row.get("support_status") or "")
        evidence_class = str(row.get("evidence_class") or "")
        blocking = []
        if support_status != "support":
            blocking.append(f"evidence_status_{support_status or 'unknown'}")
        if evidence_class not in {"direct_assay", "curated_database", "omics_association"}:
            blocking.append(f"evidence_class_{evidence_class or 'unknown'}")
        fact_candidates.append(
            {
                "candidate_fact_uid": stable_uid("fact_candidate", "evidence_assertion", row.get("assertion_uid")),
                "candidate_fact_type": "evidence_assertion",
                "subject_uid": str(row.get("subject_uid") or ""),
                "subject_type": str(row.get("subject_type") or ""),
                "predicate": str(row.get("predicate") or ""),
                "object_uid": str(row.get("object_uid") or ""),
                "object_type": str(row.get("object_type") or ""),
                "relation_source": str(row.get("method") or "evidence_assertion"),
                "source_record_uid": str(row.get("source_record_uid") or ""),
                "context_uid": str(row.get("context_uid") or "context_unspecified"),
                "direction": str(row.get("direction") or "unknown"),
                "role": "",
                "evidence_assertion_uids": [str(row.get("assertion_uid") or "")],
                "blocking_reasons": blocking or ["evidence_assertion_candidate_requires_relation_anchor"],
                "readiness_score": 0.4 if not blocking else 0.0,
                "metadata_json": json_dumps({"evidence_class": evidence_class, "support_status": support_status}),
            }
        )
    return [
        write_table(output_dir / "analysis_view" / "fact_candidates.parquet", fact_candidates, "fact_candidates"),
        write_table(output_dir / "analysis_view" / "mechanism_ready_facts.parquet", facts, "mechanism_ready_facts"),
    ]


def build_database_accuracy_store(
    workspace: Path,
    release_id: str,
    normalized_root: str = DEFAULT_NORMALIZED_ROOT,
    literature_root: str = DEFAULT_LITERATURE_ROOT,
    output_root: str = DEFAULT_OUTPUT_ROOT,
) -> Path:
    normalized_dir = workspace / normalized_root / release_id
    literature_dir = workspace / literature_root / release_id
    output_dir = workspace / output_root / release_id
    if not normalized_dir.exists():
        raise FileNotFoundError(f"normalized release not found: {normalized_dir}")
    audits = []
    audits.extend(build_entity_store(workspace, release_id, normalized_dir, output_dir))
    audits.extend(build_identity_resolution_store(workspace, release_id, output_dir))
    audits.extend(build_relation_store(workspace, release_id, normalized_dir, output_dir))
    audits.extend(build_evidence_store(release_id, normalized_dir, literature_dir, output_dir))
    audits.extend(build_analysis_view(output_dir))
    manifest = {
        "contract_version": "database_accuracy_store.v2",
        "release_id": release_id,
        "parser_name": "database_accuracy_v2",
        "parser_hash": short_hash(Path(__file__).read_text(encoding="utf-8"), 16),
        "table_count": len(audits),
        "total_rows": sum(int(row.get("rows") or 0) for row in audits),
        "source": {
            "normalized_dir": str(normalized_dir),
            "literature_dir": str(literature_dir),
        },
        "tables": audits,
        "notes": [
            "Trait, ratio, class, and exact chemical entities are split before graph projection.",
            "Graph projection remains a derived view; entity_store/relation_store/evidence_store are the accuracy-oriented fact layers.",
        ],
    }
    manifest_path = output_dir / "database_accuracy_store_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest["manifest_hash"] = short_hash(manifest, 16)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Build database accuracy v2 fact stores from a frozen normalized release.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--normalized-root", default=DEFAULT_NORMALIZED_ROOT)
    parser.add_argument("--literature-root", default=DEFAULT_LITERATURE_ROOT)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    manifest = build_database_accuracy_store(
        Path(args.workspace).resolve(),
        args.release_id,
        normalized_root=args.normalized_root,
        literature_root=args.literature_root,
        output_root=args.output_root,
    )
    print(json.dumps({"database_accuracy_store_manifest": str(manifest)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
