#!/usr/bin/env python3
"""Build a deterministic Phase 2-lite sentence evidence overlay.

The builder intentionally keeps literature-derived facts out of the curated
graph. It materializes normalized sentence mentions, high-precision rule-based
relation candidates, and an aggregated p_literature support table that can be
used by later scoring layers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import time
import unicodedata
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator

try:
    import pyarrow as pa
    import pyarrow.dataset as ds
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover - checked at runtime.
    pa = None
    ds = None
    pq = None


DEFAULT_NORMALIZED_ROOT = "normalized_store"
DEFAULT_OUTPUT_ROOT = "literature_evidence"
DEFAULT_ARTICLES_DIR = "articles_collect"
DEFAULT_PRECISION_FILTER_CONFIG = "config/evidence_precision_filters.json"
DEFAULT_MIN_RELATION_PROB = 0.20
DEFAULT_MIN_SUPPORT_PROB = 0.25
DEFAULT_MAX_MENTIONS_PER_SENTENCE = 14
DEFAULT_BATCH_SIZE = 10_000
DEFAULT_RELEVANCE_BATCH_SIZE = 16
DEFAULT_RELEVANCE_MODEL = "microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext"
BUILDER_VERSION = "phase2_lite.literature_evidence.20260513"
PRECISION_FILTER_SCHEMA_VERSION = "evidence_precision_filters.v1"
SUPPORTED_ASSERTION_STATUSES = {"support", "contradict", "uncertain", "background"}

TOKEN_RE = re.compile(r"[A-Za-z0-9]+")
RELEASE_RE = re.compile(r"^mvp_[0-9T]+$")
NEGATION_RE = re.compile(r"\b(no|not|neither|without|failed to|did not|lack(?:ed|s|ing)?|absence of|little effect|no effect)\b", re.I)
NON_NEGATING_NEGATION_RE = re.compile(r"\b(not only|not merely|not just|not simply)\b", re.I)
HEDGED_NEGATION_RE = re.compile(r"\b(not necessarily|not always|not uniformly|not exclusively|not consistently)\b", re.I)
NULL_RESULT_RE = re.compile(
    r"\b(no significant (?:difference|change|association|correlation)|"
    r"no (?:detectable |measurable |appreciable )?(?:effect|change|difference|association|correlation)|"
    r"failed to (?:show|demonstrate|detect|observe|confirm|rescue|suppress)|"
    r"did not (?:show|demonstrate|detect|observe|confirm|increase|decrease|alter|rescue|suppress)|"
    r"not (?:significantly )?(?:associated|correlated|linked|elevated|increased|decreased|changed|altered|higher|lower))\b",
    re.I,
)
UNCERTAINTY_RE = re.compile(
    r"\b(may|might|could|suggests?|potential|putative|hypothes(?:is|ized)|possible|possibly|"
    r"unclear|inconclusive|cannot exclude|remains unknown|will be required|future work)\b",
    re.I,
)
BACKGROUND_CONTEXT_RE = re.compile(
    r"\b(review|summarizes?|background|overview|known to|has been implicated|previous studies|prior studies|"
    r"has been reported|have been reported|is reported to|are reported to|canonical|hallmark)\b",
    re.I,
)
DIRECT_ASSAY_RE = re.compile(
    r"\b(measured|quantified|detected|validated|confirmed|knockdown|overexpression|inhibitor|"
    r"isotope|tracing|13c|lc[- ]?ms|ms/ms|western blot|qpcr|single[- ]cell|spatial metabolomics)\b",
    re.I,
)
CONCESSION_RE = re.compile(r"\b(although|though|whereas|while|despite|however|nevertheless|nonetheless|but)\b", re.I)
COMPARISON_CONTEXT_RE = re.compile(r"\b(compared (?:with|to)|relative to|versus|vs\.?|tumou?r[- ]?vs[- ]?normal|adjacent|matched normal|control(?:s)?)\b", re.I)
CAUSAL_WEAKENING_RE = re.compile(
    r"\b(correlation does not (?:imply|prove) caus(?:e|ality)|does not prove caus(?:e|ality)|"
    r"not sufficient|insufficient to|not required|independent of|dispensable|not necessary|necessary but not sufficient|"
    r"cannot establish caus(?:e|ality)|remains to be determined|"
    r"(?:no longer|not) significant after (?:adjustment|adjusting|controlling)|"
    r"adjust(?:ed|ing) for .{0,80}(?:abolished|attenuated|weakened))\b",
    re.I,
)
ALTERNATIVE_EXPLANATION_RE = re.compile(
    r"\b(rather than|instead of|secondary to|(?:not )?attributable to|(?:not )?explained by|confounded by|"
    r"driven by .{0,80} rather than)\b",
    re.I,
)
WEAK_OBSERVATION_RE = re.compile(r"\b(trend(?:ed)? toward|appears? to|tended to|nominal(?:ly)?|borderline|not statistically significant)\b", re.I)
CONFLICT_RE = re.compile(r"\b(conflicting|inconsistent|mixed results|discordant|controversial|contradictory)\b", re.I)
CONTEXT_BOUNDARY_RE = re.compile(
    r"\b(only in (?:cell lines|cells|mice|mouse|xenografts?|animal models?)|"
    r"in vitro but not in vivo|cell lines? but not (?:patient|tumou?r|clinical)|"
    r"not observed in (?:patient|tumou?r|clinical|in vivo)|absent in (?:patient|tumou?r|clinical|in vivo)|"
    r"only when|unless|except in|restricted to|"
    r"context[- ]dependent|subtype[- ]specific)\b",
    re.I,
)
METHOD_CUE_PATTERNS = {
    "patient_sample": re.compile(
        r"\b(patient(?:s)?|tumou?r tissue|tumou?r regions?|samples?|clinical sample|clinical cohorts?|"
        r"cohorts?|biopsy|adjacent tissue|matched normal|plasma|serum)\b",
        re.I,
    ),
    "cell_line_model": re.compile(r"\b(cell line|cells?|in vitro|culture)\b", re.I),
    "animal_model": re.compile(r"\b(mouse|mice|murine|xenograft|in vivo)\b", re.I),
    "measurement_assay": re.compile(
        r"\b(measured|quantified|detected|lc[- ]?ms|ms/ms|mass spectrometry|metabolomics|"
        r"western blot|qpcr|single[- ]cell|spatial metabolomics|isotope|tracing|13c)\b",
        re.I,
    ),
    "intervention_assay": re.compile(r"\b(knockdown|knockout|overexpression|inhibitor|inhibition|silencing|crispr|siRNA|shRNA|treated with|depletion|deprivation|withdrawal)\b", re.I),
}

ONCOLOGY_KEYWORDS = {
    "cancer",
    "cancers",
    "tumor",
    "tumors",
    "tumour",
    "tumours",
    "carcinoma",
    "adenocarcinoma",
    "sarcoma",
    "melanoma",
    "leukemia",
    "leukaemia",
    "lymphoma",
    "glioma",
    "glioblastoma",
    "neoplasm",
    "neoplasms",
    "metastasis",
    "metastatic",
    "oncogenic",
    "malignancy",
    "malignant",
}

GENERIC_SURFACES = {
    "acid",
    "activation",
    "activity",
    "analysis",
    "assay",
    "background",
    "binding",
    "body",
    "cancer",
    "cell",
    "cells",
    "control",
    "controls",
    "data",
    "disease",
    "effect",
    "effects",
    "expression",
    "factor",
    "family",
    "growth",
    "group",
    "human",
    "level",
    "levels",
    "line",
    "lines",
    "metabolism",
    "method",
    "model",
    "models",
    "pathway",
    "patient",
    "patients",
    "protein",
    "response",
    "results",
    "role",
    "sample",
    "samples",
    "study",
    "target",
    "therapy",
    "tissue",
    "tissues",
    "tumor",
    "tumors",
}

LOW_PRECISION_DISEASE_SURFACES = {
    "abnormal",
    "angiogenesis",
    "braf",
    "breast cancer cell",
    "c myc",
    "cafs",
    "cancer cell",
    "cancer immunotherapy",
    "cancer progression",
    "cancer risk",
    "cancer therapy",
    "cancer treatment",
    "carcinogenesis",
    "ccnd1",
    "central",
    "chemotherapy",
    "chemo",
    "chip",
    "cyclin d1",
    "decreasing",
    "erbb2",
    "expansion",
    "favorable",
    "gastric",
    "gli1",
    "grade",
    "her2",
    "high risk",
    "igf2bp2",
    "igf2bp3",
    "intestinal",
    "intratumoral",
    "invasion",
    "irradiation",
    "kras",
    "ldha",
    "lkb1",
    "many",
    "measured",
    "mixed",
    "moderate",
    "negative",
    "normal",
    "oncogene",
    "oncogenesis",
    "oncology",
    "ovary",
    "pcos",
    "pfkfb3",
    "positive",
    "poor",
    "primary",
    "radiation",
    "radiotherapy",
    "recurrent",
    "runx1",
    "screening",
    "site",
    "stage",
    "targeting",
    "tcga",
    "tet1",
    "the cancer genome atlas",
    "tumor cell",
    "tumor expansion",
    "tumor immunity",
    "tumor necrosis",
    "tumor microenvironment",
    "tumor necrosis factor",
    "tumor progression",
    "tumor regression",
    "tumor suppression",
    "tumor suppressor",
    "tumor tissue",
    "tumorigenesis",
    "with tumor",
    "ythdf2",
    "brca1",
    "tp53",
    "cell cycle arrest",
    "cancer metastasis",
    "tumor metastasis",
    "metastases",
    "cancer related",
    "early stage",
    "interact",
    "deep",
}

LOW_PRECISION_DISEASE_LABEL_RE = re.compile(
    r"\b("
    r"agent|therapeutic procedure|therapy agent|chemotherapy|radiation therapy|"
    r"radiotherapy|immunotherapy|screening|measurement|assay|test|criteria|"
    r"risk group|gene mutation|wt allele|allele|oncogene|gene|protein|receptor|"
    r"ligand|factor|kinase|antigen|enzyme|pathway|process|route of administration|"
    r"microenvironment|immunity|angiogenesis|invasion|carcinogenesis|tumou?rigenesis|"
    r"progression|regression|suppression|regulation|fibroblast|question|consortium|atlas"
    r")\b",
    re.I,
)

GENE_SYMBOL_BLOCKLIST = {
    "ACE",
    "ALL",
    "AND",
    "ARE",
    "AS",
    "AT",
    "CAT",
    "COX",
    "DNA",
    "FOR",
    "GAP",
    "HAD",
    "HAS",
    "IN",
    "INH",
    "ITS",
    "MAP",
    "MET",
    "NO",
    "NOT",
    "OR",
    "PR",
    "RAN",
    "SET",
    "WAS",
}

DEFAULT_PRECISION_FILTERS = {
    "schema_version": PRECISION_FILTER_SCHEMA_VERSION,
    "surface_blocklist": [
        {"surface": "can", "entity_types": ["gene", "target"], "matched_fields": ["alias"], "reason": "English modal verb."},
        {"surface": "via", "entity_types": ["gene", "target"], "matched_fields": ["alias"], "reason": "English preposition."},
        {"surface": "hcc", "entity_types": ["gene", "target"], "matched_fields": ["alias"], "reason": "Disease abbreviation, not gene alias."},
        {"surface": "protein kinase", "entity_types": ["gene", "target", "pathway"], "reason": "Generic protein family phrase."},
        {"surface": "kinase", "entity_types": ["gene", "target", "pathway"], "reason": "Generic protein family word."},
        {"surface": "proteins", "entity_types": ["metabolite"], "matched_fields": ["synonym"], "reason": "Generic biomolecule class."},
        {"surface": "mrna", "entity_types": ["metabolite"], "matched_fields": ["synonym"], "reason": "Generic transcript class."},
    ],
    "surface_downweight": [],
    "matched_field_downweight": [],
}

SECTION_WEIGHTS = {
    "title": 1.05,
    "abstract": 1.00,
    "results": 1.00,
    "conclusion": 1.00,
    "conclusions": 1.00,
    "discussion": 0.75,
    "introduction": 0.65,
    "background": 0.60,
    "full_text": 0.70,
    "methods": 0.20,
    "materials and methods": 0.20,
    "references": 0.05,
}

POLARITY_OPPOSITES = {
    "increase": "decrease",
    "decrease": "increase",
    "activation": "inhibition",
    "inhibition": "activation",
}


def require_arrow() -> None:
    if pa is None or ds is None or pq is None:
        raise RuntimeError("pyarrow is required for literature evidence building.")


def stable_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def content_hash(value: Any) -> str:
    return hashlib.sha256(stable_json(value).encode("utf-8")).hexdigest()


def short_hash(value: str, length: int = 20) -> str:
    return hashlib.sha1(value.encode("utf-8")).hexdigest()[:length]


def stable_uid(prefix: str, *parts: object) -> str:
    payload = "|".join(str(part).strip() for part in parts if part is not None and str(part).strip())
    return f"{prefix}_{short_hash(payload)}"


def parser_hash() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]


def normalize_text(value: Any) -> str:
    text = unicodedata.normalize("NFKC", str(value or ""))
    text = text.replace("鈥", "'").replace("怣", "m").replace("恈", "c").replace("恡", "t")
    text = re.sub(r"[\u2010-\u2015\-_/]+", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip().casefold()


def iter_tokens_with_offsets(text: str) -> list[tuple[str, int, int]]:
    return [(match.group(0).casefold(), match.start(), match.end()) for match in TOKEN_RE.finditer(text)]


def surface_tokens(surface: str) -> tuple[str, ...]:
    return tuple(token for token, _start, _end in iter_tokens_with_offsets(normalize_text(surface)))


def clamp_unit(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        number = default
    if math.isnan(number) or math.isinf(number):
        number = default
    return min(1.0, max(0.0, number))


def independent_probability_union(values: Iterable[float]) -> float:
    remaining = 1.0
    for value in values:
        remaining *= 1.0 - clamp_unit(value)
    return clamp_unit(1.0 - remaining)


def parquet_row_count(path: Path) -> int:
    if not path.exists():
        return 0
    return int(pq.ParquetFile(path).metadata.num_rows)


def latest_release_id(root: Path) -> str:
    candidates = [path.name for path in root.iterdir() if path.is_dir() and RELEASE_RE.match(path.name)] if root.exists() else []
    if not candidates:
        raise FileNotFoundError(f"No mvp_* release directories found in {root}")
    return sorted(candidates)[-1]


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def deep_copy_json(value: Any) -> Any:
    return json.loads(json.dumps(value, ensure_ascii=False))


def merge_precision_filters(config: dict[str, Any]) -> dict[str, Any]:
    merged = deep_copy_json(DEFAULT_PRECISION_FILTERS)
    for key in ("surface_blocklist", "surface_downweight", "matched_field_downweight"):
        merged[key].extend(config.get(key, []) or [])
    merged["schema_version"] = str(config.get("schema_version") or merged["schema_version"])
    if config.get("purpose"):
        merged["purpose"] = config["purpose"]
    return merged


def load_precision_filters(workspace: Path, config_path: str | None) -> tuple[dict[str, Any], str]:
    if not config_path:
        return merge_precision_filters({}), ""
    path = Path(config_path)
    if not path.is_absolute():
        path = workspace / path
    if not path.exists():
        return merge_precision_filters({}), str(path)
    return read_json(path), str(path)


def precision_filter_stats(filters: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": filters.get("schema_version", ""),
        "surface_blocklist_count": len(filters.get("surface_blocklist", []) or []),
        "surface_downweight_count": len(filters.get("surface_downweight", []) or []),
        "matched_field_downweight_count": len(filters.get("matched_field_downweight", []) or []),
        "generic_blocked_surface_count": 0,
        "precision_blocked_surface_count": 0,
        "precision_downweighted_surface_count": 0,
        "field_downweighted_surface_count": 0,
        "blocklist_hits_by_surface": Counter(),
        "downweight_hits_by_surface": Counter(),
        "field_downweight_hits_by_field": Counter(),
    }


def normalized_surface_key(surface: Any) -> str:
    return " ".join(surface_tokens(str(surface or "")))


def rule_matches(rule: dict[str, Any], normalized: str, entity_type: str, matched_field: str) -> bool:
    rule_surface = normalized_surface_key(rule.get("surface", ""))
    if rule_surface and rule_surface != normalized:
        return False
    entity_types = set(str(item) for item in (rule.get("entity_types") or []))
    if entity_types and "*" not in entity_types and entity_type not in entity_types:
        return False
    matched_fields = set(str(item) for item in (rule.get("matched_fields") or []))
    if matched_fields and "*" not in matched_fields and matched_field not in matched_fields:
        return False
    return True


def is_precision_blocked(filters: dict[str, Any], normalized: str, entity_type: str, matched_field: str) -> bool:
    return any(rule_matches(rule, normalized, entity_type, matched_field) for rule in filters.get("surface_blocklist", []) or [])


def precision_confidence_multiplier(filters: dict[str, Any], normalized: str, entity_type: str, matched_field: str) -> tuple[float, bool, bool]:
    multiplier = 1.0
    surface_hit = False
    field_hit = False
    for rule in filters.get("surface_downweight", []) or []:
        if rule_matches(rule, normalized, entity_type, matched_field):
            multiplier *= float(rule.get("confidence_multiplier", 1.0) or 1.0)
            surface_hit = True
    for rule in filters.get("matched_field_downweight", []) or []:
        if rule_matches(rule, "", entity_type, matched_field):
            multiplier *= float(rule.get("confidence_multiplier", 1.0) or 1.0)
            field_hit = True
    return max(0.0, min(1.0, multiplier)), surface_hit, field_hit


def finalize_precision_stats(stats: dict[str, Any], filters: dict[str, Any], path: str, config_hash: str) -> dict[str, Any]:
    return {
        **{key: value for key, value in stats.items() if not isinstance(value, Counter)},
        "config_path": path,
        "config_hash": config_hash,
        "blocklist_hits_by_surface": dict(sorted(stats["blocklist_hits_by_surface"].items())),
        "downweight_hits_by_surface": dict(sorted(stats["downweight_hits_by_surface"].items())),
        "field_downweight_hits_by_field": dict(sorted(stats["field_downweight_hits_by_field"].items())),
        "notes": "Precision filters only affect the literature evidence overlay lexicon; curated graph rows are unchanged.",
        "active_filter_hash": content_hash(filters)[:16],
    }


def table_dataset(path: Path) -> ds.Dataset:
    if not path.exists():
        raise FileNotFoundError(f"Missing parquet table: {path}")
    return ds.dataset(path, format="parquet")


def iter_table_rows(path: Path, columns: list[str] | None = None, batch_size: int = DEFAULT_BATCH_SIZE) -> Iterator[dict[str, Any]]:
    parquet_file = pq.ParquetFile(path)
    for batch in parquet_file.iter_batches(columns=columns, batch_size=batch_size):
        for row in batch.to_pylist():
            yield row


def safe_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value if item is not None and str(item).strip()]
    if isinstance(value, tuple):
        return [str(item) for item in value if item is not None and str(item).strip()]
    text = str(value).strip()
    return [text] if text else []


def is_generic_surface(surface: str, entity_type: str) -> bool:
    normalized = normalize_text(surface)
    tokens = surface_tokens(surface)
    if not normalized or not tokens:
        return True
    if len(normalized) > 96 or len(tokens) > 10:
        return True
    if not any(any(ch.isalpha() for ch in token) for token in tokens):
        return True
    if normalized in GENERIC_SURFACES:
        return True
    if len(tokens) == 1 and len(tokens[0]) < 4 and entity_type not in {"gene", "target"}:
        return True
    if entity_type in {"gene", "target"}:
        upper = str(surface or "").strip().upper()
        if upper in GENE_SYMBOL_BLOCKLIST:
            return True
        if len(tokens) == 1 and len(tokens[0]) < 3:
            return True
    return False


def is_low_precision_disease_surface(surface: Any, display_name: Any, matched_field: str) -> bool:
    normalized = normalized_surface_key(surface)
    if normalized in LOW_PRECISION_DISEASE_SURFACES:
        return True
    label = normalize_text(display_name)
    if LOW_PRECISION_DISEASE_LABEL_RE.search(label):
        return True
    if matched_field == "alias" and normalized in GENERIC_SURFACES:
        return True
    return False


def is_oncology_label(name: Any, aliases: Iterable[str] = ()) -> bool:
    text = normalize_text(" ".join([str(name or ""), *[str(alias or "") for alias in aliases]]))
    return any(keyword in text.split() or keyword in text for keyword in ONCOLOGY_KEYWORDS)


@dataclass(frozen=True)
class LexiconEntry:
    tokens: tuple[str, ...]
    surface: str
    normalized_surface: str
    entity_uid: str
    entity_type: str
    display_name: str
    matched_field: str
    resolution_confidence: float
    source_release: str
    license_id: str


@dataclass
class Mention:
    mention_uid: str
    sentence_uid: str
    article_uid: str
    pmid: str
    pmcid: str
    section: str
    entity_uid: str
    entity_type: str
    display_name: str
    surface: str
    normalized_surface: str
    start_offset: int
    end_offset: int
    matched_field: str
    resolution_confidence: float
    source_release: str
    license_id: str
    parser_hash: str
    config_hash: str

    def as_row(self) -> dict[str, Any]:
        return {
            "mention_uid": self.mention_uid,
            "sentence_uid": self.sentence_uid,
            "article_uid": self.article_uid,
            "pmid": self.pmid,
            "pmcid": self.pmcid,
            "section": self.section,
            "entity_uid": self.entity_uid,
            "entity_type": self.entity_type,
            "display_name": self.display_name,
            "surface": self.surface,
            "normalized_surface": self.normalized_surface,
            "start_offset": self.start_offset,
            "end_offset": self.end_offset,
            "matched_field": self.matched_field,
            "resolution_confidence": round(self.resolution_confidence, 6),
            "source_release": self.source_release,
            "license_id": self.license_id,
            "parser_hash": self.parser_hash,
            "config_hash": self.config_hash,
        }


@dataclass(frozen=True)
class Trigger:
    phrase: str
    polarity: str
    raw_score: float
    rule_id: str


MENTION_SCHEMA = pa.schema(
    [
        ("mention_uid", pa.string()),
        ("sentence_uid", pa.string()),
        ("article_uid", pa.string()),
        ("pmid", pa.string()),
        ("pmcid", pa.string()),
        ("section", pa.string()),
        ("entity_uid", pa.string()),
        ("entity_type", pa.string()),
        ("display_name", pa.string()),
        ("surface", pa.string()),
        ("normalized_surface", pa.string()),
        ("start_offset", pa.int64()),
        ("end_offset", pa.int64()),
        ("matched_field", pa.string()),
        ("resolution_confidence", pa.float64()),
        ("source_release", pa.string()),
        ("license_id", pa.string()),
        ("parser_hash", pa.string()),
        ("config_hash", pa.string()),
    ]
) if pa is not None else None

RELATION_SCHEMA = pa.schema(
    [
        ("relation_uid", pa.string()),
        ("subject_uid", pa.string()),
        ("subject_type", pa.string()),
        ("predicate", pa.string()),
        ("object_uid", pa.string()),
        ("object_type", pa.string()),
        ("polarity", pa.string()),
        ("cancer_context", pa.string()),
        ("sentence_uid", pa.string()),
        ("pmid", pa.string()),
        ("pmcid", pa.string()),
        ("section", pa.string()),
        ("trigger_phrase", pa.string()),
        ("extraction_rule_id", pa.string()),
        ("support_status", pa.string()),
        ("raw_score", pa.float64()),
        ("calibrated_prob", pa.float64()),
        ("score_components_json", pa.string()),
        ("license_id", pa.string()),
        ("source_release", pa.string()),
        ("parser_hash", pa.string()),
        ("config_hash", pa.string()),
    ]
) if pa is not None else None

SUPPORT_SCHEMA = pa.schema(
    [
        ("support_uid", pa.string()),
        ("subject_uid", pa.string()),
        ("subject_type", pa.string()),
        ("predicate", pa.string()),
        ("object_uid", pa.string()),
        ("object_type", pa.string()),
        ("polarity_set", pa.list_(pa.string())),
        ("support_status_set", pa.list_(pa.string())),
        ("support_class", pa.string()),
        ("supported_existing_edge_uids", pa.list_(pa.string())),
        ("evidence_relation_uids", pa.list_(pa.string())),
        ("sentence_uids", pa.list_(pa.string())),
        ("pmids", pa.list_(pa.string())),
        ("pmcids", pa.list_(pa.string())),
        ("evidence_sentence_count", pa.int64()),
        ("distinct_article_count", pa.int64()),
        ("p_literature", pa.float64()),
        ("raw_score_max", pa.float64()),
        ("calibrated_prob_max", pa.float64()),
        ("score_components_json", pa.string()),
        ("license_id", pa.string()),
        ("source_release", pa.string()),
        ("parser_hash", pa.string()),
        ("config_hash", pa.string()),
    ]
) if pa is not None else None

RELEVANCE_SCHEMA = pa.schema(
    [
        ("sentence_uid", pa.string()),
        ("article_uid", pa.string()),
        ("pmid", pa.string()),
        ("pmcid", pa.string()),
        ("section", pa.string()),
        ("model_name", pa.string()),
        ("relevance_mode", pa.string()),
        ("relevance_score", pa.float64()),
        ("decision", pa.string()),
        ("score_components_json", pa.string()),
        ("source_release", pa.string()),
        ("parser_hash", pa.string()),
        ("config_hash", pa.string()),
    ]
) if pa is not None else None


class ParquetBatchWriter:
    def __init__(self, path: Path, schema: pa.Schema):
        self.path = path
        self.schema = schema
        self.writer: pq.ParquetWriter | None = None
        self.rows = 0

    def write(self, rows: list[dict[str, Any]]) -> None:
        if not rows:
            return
        table = pa.Table.from_pylist(rows, schema=self.schema)
        if self.writer is None:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.writer = pq.ParquetWriter(self.path, self.schema)
        self.writer.write_table(table)
        self.rows += len(rows)

    def close(self) -> None:
        if self.writer is None:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            pq.write_table(pa.Table.from_pylist([], schema=self.schema), self.path)
            return
        self.writer.close()


class PubMedBertRelevanceScorer:
    """Embedding-similarity relevance scorer; it filters/reranks but does not extract facts."""

    anchors = [
        "tumor metabolism mechanism metabolic reprogramming cancer metabolite pathway gene regulation",
        "cancer cells glycolysis lactate glutamine lipid metabolism metabolic flux assay",
        "metabolite abundance altered in tumor tissue pathway enzyme transporter mechanism",
    ]

    def __init__(self, model_name: str):
        self.model_name = model_name
        try:
            import torch  # type: ignore
            from transformers import AutoModel, AutoTokenizer  # type: ignore
        except ImportError as exc:  # pragma: no cover - depends on optional env.
            raise RuntimeError("PubMedBERT relevance mode requires optional dependencies 'torch' and 'transformers'.") from exc
        self.torch = torch
        try:
            self.tokenizer = AutoTokenizer.from_pretrained(model_name)
            self.model = AutoModel.from_pretrained(model_name)
        except Exception as exc:  # pragma: no cover - depends on optional env/model cache.
            raise RuntimeError(f"Unable to load relevance model '{model_name}'. Install/cache the model before using --relevance-mode pubmedbert.") from exc
        self.model.eval()
        with torch.no_grad():
            self.anchor_embedding = torch.nn.functional.normalize(self._embed(self.anchors).mean(dim=0, keepdim=True), p=2, dim=1)

    def _embed(self, texts: list[str]):
        torch = self.torch
        encoded = self.tokenizer(texts, padding=True, truncation=True, max_length=256, return_tensors="pt")
        with torch.no_grad():
            output = self.model(**encoded)
            hidden = output.last_hidden_state
            mask = encoded["attention_mask"].unsqueeze(-1).to(hidden.dtype)
            pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1.0)
            return torch.nn.functional.normalize(pooled, p=2, dim=1)

    def score_many(self, texts: list[str]) -> list[tuple[float, dict[str, Any]]]:
        torch = self.torch
        cleaned = [re.sub(r"\s+", " ", str(text or "")).strip() for text in texts]
        results: list[tuple[float, dict[str, Any]]] = []
        nonempty_indices = [idx for idx, text in enumerate(cleaned) if text]
        scores_by_index: dict[int, tuple[float, dict[str, Any]]] = {}
        if nonempty_indices:
            embeddings = self._embed([cleaned[idx] for idx in nonempty_indices])
            cosines = torch.matmul(embeddings, self.anchor_embedding.T).squeeze(dim=1).tolist()
            for idx, cosine_value in zip(nonempty_indices, cosines):
                cosine = float(cosine_value)
                relevance_score = clamp_unit((cosine + 1.0) / 2.0)
                scores_by_index[idx] = (
                    relevance_score,
                    {
                        "scorer": "pubmedbert_embedding_similarity",
                        "cosine_to_tumor_metabolism_anchor": round(cosine, 6),
                        "formula": "(cosine(sentence_embedding,tumor_metabolism_anchor)+1)/2",
                        "fact_boundary": "relevance score is used only for filtering/reranking and probability weighting",
                    },
                )
        for idx, text in enumerate(cleaned):
            results.append(scores_by_index.get(idx, (0.0, {"scorer": "pubmedbert_embedding_similarity", "empty_text": not bool(text)})))
        return results

    def score(self, text: str) -> tuple[float, dict[str, Any]]:
        return self.score_many([text])[0]


def relevance_decision(score: float) -> str:
    if score >= 0.65:
        return "evidence_candidate"
    if score >= 0.45:
        return "background"
    return "skip"


def relevance_row(row: dict[str, Any], scorer: PubMedBertRelevanceScorer, cfg_hash: str) -> tuple[dict[str, Any], float]:
    score, components = scorer.score(str(row.get("sentence_text") or ""))
    result = {
        "sentence_uid": str(row.get("sentence_uid") or ""),
        "article_uid": str(row.get("article_uid") or ""),
        "pmid": str(row.get("pmid") or ""),
        "pmcid": str(row.get("pmcid") or ""),
        "section": str(row.get("section") or ""),
        "model_name": scorer.model_name,
        "relevance_mode": "pubmedbert",
        "relevance_score": round(score, 6),
        "decision": relevance_decision(score),
        "score_components_json": stable_json(components),
        "source_release": str(row.get("source_release") or ""),
        "parser_hash": str(row.get("parser_hash") or ""),
        "config_hash": cfg_hash,
    }
    return result, score


def relevance_rows(rows: list[dict[str, Any]], scorer: PubMedBertRelevanceScorer, cfg_hash: str) -> list[tuple[dict[str, Any], float]]:
    scored = scorer.score_many([str(row.get("sentence_text") or "") for row in rows])
    output: list[tuple[dict[str, Any], float]] = []
    for row, (score, components) in zip(rows, scored):
        output.append(
            (
                {
                    "sentence_uid": str(row.get("sentence_uid") or ""),
                    "article_uid": str(row.get("article_uid") or ""),
                    "pmid": str(row.get("pmid") or ""),
                    "pmcid": str(row.get("pmcid") or ""),
                    "section": str(row.get("section") or ""),
                    "model_name": scorer.model_name,
                    "relevance_mode": "pubmedbert",
                    "relevance_score": round(score, 6),
                    "decision": relevance_decision(score),
                    "score_components_json": stable_json(components),
                    "source_release": str(row.get("source_release") or ""),
                    "parser_hash": str(row.get("parser_hash") or ""),
                    "config_hash": cfg_hash,
                },
                score,
            )
        )
    return output


def apply_relevance_weight(row: dict[str, Any], relevance_score: float | None, min_prob: float) -> dict[str, Any] | None:
    if relevance_score is None:
        return row
    original_prob = clamp_unit(row.get("calibrated_prob", 0.0))
    multiplier = 0.5 + 0.5 * clamp_unit(relevance_score)
    final_prob = clamp_unit(original_prob * multiplier)
    if final_prob < min_prob:
        return None
    components = read_score_components(row)
    components["relevance_weight"] = {
        "relevance_score": round(clamp_unit(relevance_score), 6),
        "multiplier": round(multiplier, 6),
        "original_calibrated_prob": round(original_prob, 6),
        "formula": "rule_relation_prob*(0.5+0.5*relevance_score)",
        "fact_boundary": "model relevance only downweights/reranks rule-extracted relation candidates",
    }
    row = dict(row)
    row["calibrated_prob"] = round(final_prob, 6)
    row["score_components_json"] = stable_json(components)
    return row


def read_score_components(row: dict[str, Any]) -> dict[str, Any]:
    try:
        value = json.loads(str(row.get("score_components_json") or "{}"))
    except json.JSONDecodeError:
        return {}
    return value if isinstance(value, dict) else {}


def add_surface(
    surfaces: dict[tuple[str, str], list[tuple[float, LexiconEntry]]],
    precision_filters: dict[str, Any],
    precision_stats: dict[str, Any],
    entity_type: str,
    entity_uid: str,
    display_name: str,
    surface: Any,
    matched_field: str,
    base_confidence: float,
    source_release: str,
    license_id: str,
) -> None:
    text = str(surface or "").strip()
    if is_generic_surface(text, entity_type):
        precision_stats["generic_blocked_surface_count"] += 1
        return
    if entity_type == "disease" and is_low_precision_disease_surface(text, display_name, matched_field):
        precision_stats["generic_blocked_surface_count"] += 1
        return
    tokens = surface_tokens(text)
    if not tokens:
        return
    normalized = " ".join(tokens)
    if is_precision_blocked(precision_filters, normalized, entity_type, matched_field):
        precision_stats["precision_blocked_surface_count"] += 1
        precision_stats["blocklist_hits_by_surface"][f"{entity_type}:{matched_field}:{normalized}"] += 1
        return
    multiplier, surface_downweighted, field_downweighted = precision_confidence_multiplier(
        precision_filters,
        normalized,
        entity_type,
        matched_field,
    )
    if surface_downweighted:
        precision_stats["precision_downweighted_surface_count"] += 1
        precision_stats["downweight_hits_by_surface"][f"{entity_type}:{matched_field}:{normalized}"] += 1
    if field_downweighted:
        precision_stats["field_downweighted_surface_count"] += 1
        precision_stats["field_downweight_hits_by_field"][f"{entity_type}:{matched_field}"] += 1
    adjusted_confidence = round(max(0.0, min(1.0, base_confidence * multiplier)), 6)
    entry = LexiconEntry(
        tokens=tokens,
        surface=text,
        normalized_surface=normalized,
        entity_uid=entity_uid,
        entity_type=entity_type,
        display_name=display_name,
        matched_field=matched_field,
        resolution_confidence=adjusted_confidence,
        source_release=source_release,
        license_id=license_id,
    )
    field_bonus = {
        "canonical_name": 0.08,
        "symbol": 0.08,
        "approved_symbol": 0.08,
        "name": 0.05,
        "preferred_name": 0.04,
        "synonym": 0.0,
        "alias": 0.0,
    }.get(matched_field, 0.0)
    token_bonus = min(0.06, 0.015 * len(tokens))
    priority = min(1.0, adjusted_confidence + field_bonus + token_bonus)
    surfaces[(entity_type, normalized)].append((priority, entry))


def collapse_surfaces(surfaces: dict[tuple[str, str], list[tuple[float, LexiconEntry]]]) -> tuple[list[LexiconEntry], dict[str, int]]:
    entries: list[LexiconEntry] = []
    stats = {"surface_keys": len(surfaces), "ambiguous_surface_keys": 0, "retained_surface_keys": 0}
    for (_entity_type, _surface), candidates in sorted(surfaces.items()):
        candidates.sort(key=lambda item: (-item[0], item[1].entity_uid, item[1].matched_field))
        top_score, top = candidates[0]
        runner_up = candidates[1][0] if len(candidates) > 1 else -1.0
        unique_uids = {entry.entity_uid for _score, entry in candidates}
        if len(unique_uids) > 1 and top_score - runner_up < 0.12:
            stats["ambiguous_surface_keys"] += 1
            continue
        entries.append(top)
        stats["retained_surface_keys"] += 1
    entries.sort(key=lambda entry: (entry.tokens, entry.entity_type, entry.entity_uid))
    return entries, stats


def load_graph_limited_uids(normalized_dir: Path) -> tuple[set[str], set[str]]:
    metabolite_uids: set[str] = set()
    gene_uids: set[str] = set()
    met_edges = normalized_dir / "metabolite_pathway_edges.parquet"
    if met_edges.exists():
        for row in iter_table_rows(met_edges, ["metabolite_uid"]):
            if row.get("metabolite_uid"):
                metabolite_uids.add(row["metabolite_uid"])
    gene_edges = normalized_dir / "gene_pathway_edges.parquet"
    if gene_edges.exists():
        for row in iter_table_rows(gene_edges, ["gene_uid"]):
            if row.get("gene_uid"):
                gene_uids.add(row["gene_uid"])
    return metabolite_uids, gene_uids


def load_article_license_map(normalized_dir: Path) -> dict[str, str]:
    path = normalized_dir / "articles.parquet"
    if not path.exists():
        return {}
    article_license: dict[str, str] = {}
    for row in iter_table_rows(path, ["article_uid", "license_id"]):
        article_uid = str(row.get("article_uid") or "")
        if article_uid:
            article_license[article_uid] = str(row.get("license_id") or "local_articles:articles_collect")
    return article_license


def build_lexicon(
    normalized_dir: Path,
    limit_to_graph_entities: bool = True,
    precision_filters: dict[str, Any] | None = None,
    precision_stats: dict[str, Any] | None = None,
) -> tuple[dict[str, list[LexiconEntry]], dict[str, Any]]:
    precision_filters = precision_filters or merge_precision_filters({})
    precision_stats = precision_stats or precision_filter_stats(precision_filters)
    surfaces: dict[tuple[str, str], list[tuple[float, LexiconEntry]]] = defaultdict(list)
    graph_metabolites, graph_genes = load_graph_limited_uids(normalized_dir) if limit_to_graph_entities else (set(), set())
    target_gene_uids: set[str] = set()

    targets_path = normalized_dir / "targets.parquet"
    if targets_path.exists():
        for row in iter_table_rows(targets_path, ["target_uid", "preferred_name", "approved_symbol", "gene_uid", "source_release", "license_id"]):
            target_uid = str(row.get("target_uid") or "")
            if not target_uid:
                continue
            gene_uid = str(row.get("gene_uid") or "")
            if gene_uid:
                target_gene_uids.add(gene_uid)
            label = str(row.get("approved_symbol") or row.get("preferred_name") or target_uid)
            add_surface(surfaces, precision_filters, precision_stats, "target", target_uid, label, row.get("approved_symbol"), "approved_symbol", 0.94, row.get("source_release", ""), row.get("license_id", ""))
            add_surface(surfaces, precision_filters, precision_stats, "target", target_uid, label, row.get("preferred_name"), "preferred_name", 0.86, row.get("source_release", ""), row.get("license_id", ""))

    metabolites_path = normalized_dir / "metabolites.parquet"
    if metabolites_path.exists():
        columns = ["metabolite_uid", "canonical_name", "synonyms", "source_release", "license_id"]
        for row in iter_table_rows(metabolites_path, columns):
            metabolite_uid = str(row.get("metabolite_uid") or "")
            if limit_to_graph_entities and metabolite_uid not in graph_metabolites:
                continue
            label = str(row.get("canonical_name") or metabolite_uid)
            add_surface(surfaces, precision_filters, precision_stats, "metabolite", metabolite_uid, label, label, "canonical_name", 0.93, row.get("source_release", ""), row.get("license_id", ""))
            for synonym in safe_list(row.get("synonyms"))[:24]:
                add_surface(surfaces, precision_filters, precision_stats, "metabolite", metabolite_uid, label, synonym, "synonym", 0.82, row.get("source_release", ""), row.get("license_id", ""))

    genes_path = normalized_dir / "genes.parquet"
    if genes_path.exists():
        columns = ["gene_uid", "symbol", "aliases", "biotype", "taxon", "source_release", "license_id"]
        for row in iter_table_rows(genes_path, columns):
            gene_uid = str(row.get("gene_uid") or "")
            if limit_to_graph_entities and gene_uid not in graph_genes and gene_uid not in target_gene_uids:
                continue
            if str(row.get("taxon") or "") not in {"", "9606"}:
                continue
            label = str(row.get("symbol") or gene_uid)
            add_surface(surfaces, precision_filters, precision_stats, "gene", gene_uid, label, label, "symbol", 0.94, row.get("source_release", ""), row.get("license_id", ""))
            for alias in safe_list(row.get("aliases"))[:12]:
                add_surface(surfaces, precision_filters, precision_stats, "gene", gene_uid, label, alias, "alias", 0.80, row.get("source_release", ""), row.get("license_id", ""))

    pathways_path = normalized_dir / "pathways.parquet"
    if pathways_path.exists():
        for row in iter_table_rows(pathways_path, ["pathway_uid", "name", "species", "source_release", "license_id"]):
            if str(row.get("species") or "") not in {"", "Homo sapiens"}:
                continue
            pathway_uid = str(row.get("pathway_uid") or "")
            label = str(row.get("name") or pathway_uid)
            add_surface(surfaces, precision_filters, precision_stats, "pathway", pathway_uid, label, label, "name", 0.90, row.get("source_release", ""), row.get("license_id", ""))

    diseases_path = normalized_dir / "diseases.parquet"
    if diseases_path.exists():
        for row in iter_table_rows(diseases_path, ["disease_uid", "name", "aliases", "source_release", "license_id"]):
            aliases = safe_list(row.get("aliases"))
            if not is_oncology_label(row.get("name"), aliases):
                continue
            disease_uid = str(row.get("disease_uid") or "")
            label = str(row.get("name") or disease_uid)
            add_surface(surfaces, precision_filters, precision_stats, "disease", disease_uid, label, label, "name", 0.91, row.get("source_release", ""), row.get("license_id", ""))
            for alias in aliases[:16]:
                add_surface(surfaces, precision_filters, precision_stats, "disease", disease_uid, label, alias, "alias", 0.80, row.get("source_release", ""), row.get("license_id", ""))

    entries, stats = collapse_surfaces(surfaces)
    index: dict[str, list[LexiconEntry]] = defaultdict(list)
    for entry in entries:
        index[entry.tokens[0]].append(entry)
    for bucket in index.values():
        bucket.sort(key=lambda entry: (-len(entry.tokens), entry.entity_type, entry.entity_uid))
    stats.update(
        {
            "entry_count": len(entries),
            "first_token_count": len(index),
            "limited_to_graph_entities": limit_to_graph_entities,
            "graph_metabolite_uid_count": len(graph_metabolites),
            "graph_gene_uid_count": len(graph_genes),
        }
    )
    return dict(index), stats


def match_mentions(row: dict[str, Any], lexicon: dict[str, list[LexiconEntry]], phash: str, cfg_hash: str, max_mentions: int) -> list[Mention]:
    text = str(row.get("sentence_text") or "")
    tokens = iter_tokens_with_offsets(text)
    mentions: list[Mention] = []
    occupied: set[tuple[int, int, str, str]] = set()
    occupied_type_spans: list[tuple[int, int, str]] = []
    seen_entity_spans: set[tuple[str, int, int]] = set()
    for idx, (token, start, _end) in enumerate(tokens):
        for entry in lexicon.get(token, []):
            width = len(entry.tokens)
            if idx + width > len(tokens):
                continue
            if tuple(tok for tok, _s, _e in tokens[idx : idx + width]) != entry.tokens:
                continue
            span_start = start
            span_end = tokens[idx + width - 1][2]
            span_key = (span_start, span_end, entry.entity_type, entry.entity_uid)
            if any(
                existing_type == entry.entity_type and not (span_end <= existing_start or span_start >= existing_end)
                for existing_start, existing_end, existing_type in occupied_type_spans
            ):
                continue
            if span_key in occupied or (entry.entity_uid, span_start, span_end) in seen_entity_spans:
                continue
            surface = text[span_start:span_end]
            mention_uid = stable_uid("mention", row.get("sentence_uid", ""), entry.entity_uid, span_start, span_end)
            mentions.append(
                Mention(
                    mention_uid=mention_uid,
                    sentence_uid=str(row.get("sentence_uid") or ""),
                    article_uid=str(row.get("article_uid") or ""),
                    pmid=str(row.get("pmid") or ""),
                    pmcid=str(row.get("pmcid") or ""),
                    section=str(row.get("section") or ""),
                    entity_uid=entry.entity_uid,
                    entity_type=entry.entity_type,
                    display_name=entry.display_name,
                    surface=surface,
                    normalized_surface=entry.normalized_surface,
                    start_offset=span_start,
                    end_offset=span_end,
                    matched_field=entry.matched_field,
                    resolution_confidence=entry.resolution_confidence,
                    source_release=str(row.get("source_release") or entry.source_release),
                    license_id=str(row.get("_article_license_id") or "local_articles:articles_collect"),
                    parser_hash=phash,
                    config_hash=cfg_hash,
                )
            )
            occupied.add(span_key)
            occupied_type_spans.append((span_start, span_end, entry.entity_type))
            seen_entity_spans.add((entry.entity_uid, span_start, span_end))
    mentions.sort(key=lambda item: (item.start_offset, -(item.end_offset - item.start_offset), item.entity_type, item.entity_uid))
    if len(mentions) > max_mentions:
        mentions = sorted(mentions, key=lambda item: (-item.resolution_confidence, item.start_offset, item.entity_uid))[:max_mentions]
        mentions.sort(key=lambda item: (item.start_offset, item.entity_type, item.entity_uid))
    return mentions


def find_trigger(text: str, patterns: list[tuple[str, str, float, str]]) -> Trigger | None:
    normalized = normalize_text(text)
    for pattern, polarity, score, rule_id in patterns:
        match = re.search(pattern, normalized)
        if match:
            return Trigger(match.group(0), polarity, score, rule_id)
    return None


CHANGE_PATTERNS = [
    (r"\b(elevated|increased|upregulated|overexpressed|higher|accumulated|enriched)\b", "increase", 0.78, "changed_in_cancer.increase"),
    (r"\b(decreased|reduced|downregulated|lower|depleted|suppressed)\b", "decrease", 0.78, "changed_in_cancer.decrease"),
]
ASSOCIATION_PATTERNS = [
    (r"\b(associated with|correlated with|linked to|predictive of|marker of|biomarker for|poor prognosis|overall survival|progression free survival)\b", "association", 0.70, "association.explicit"),
]
REGULATION_PATTERNS = [
    (r"\b(activates|activated|promotes|promoted|enhances|enhanced|induces|induced|drives|stabilizes|upregulates|upregulated)\b", "activation", 0.74, "regulates.activation"),
    (r"\b(inhibits|inhibited|suppresses|suppressed|blocks|blocked|reduces|reduced|downregulates|downregulated)\b", "inhibition", 0.74, "regulates.inhibition"),
    (r"\b(regulates|regulated|modulates|modulated|controls|controlled)\b", "association", 0.66, "regulates.generic"),
]


def mention_context(mentions: list[Mention]) -> str:
    diseases = [mention.display_name for mention in mentions if mention.entity_type == "disease"]
    return "; ".join(sorted(set(diseases))) if diseases else "oncology_keyword"


def section_weight(section: str) -> float:
    return SECTION_WEIGHTS.get(normalize_text(section), 0.70)


def oncology_context_weight(text: str, mentions: list[Mention]) -> float:
    if any(mention.entity_type == "disease" for mention in mentions):
        return 1.0
    normalized = normalize_text(text)
    return 0.85 if any(keyword in normalized for keyword in ONCOLOGY_KEYWORDS) else 0.0


def direction_weight(polarity: str) -> float:
    return 1.0 if polarity in {"increase", "decrease", "activation", "inhibition"} else 0.82


def pair_distance_weight(subject: Mention, obj: Mention) -> float:
    distance = pair_distance(subject, obj)
    if distance <= 80:
        return 1.0
    if distance <= 160:
        return 0.88
    return 0.72


def mention_density_weight(mentions: list[Mention]) -> float:
    if len(mentions) <= 6:
        return 1.0
    return round(max(0.72, 1.0 - 0.035 * (len(mentions) - 6)), 6)


def certainty_context_weight(text: str) -> float:
    weight = 1.0
    semantics = assertion_semantics(text)
    if UNCERTAINTY_RE.search(text) or WEAK_OBSERVATION_RE.search(text):
        weight *= 0.72
    if semantics["support_status"] == "contradict":
        weight *= 0.70
    if CAUSAL_WEAKENING_RE.search(text):
        weight *= 0.68
    if CONFLICT_RE.search(text):
        weight *= 0.70
    if CONTEXT_BOUNDARY_RE.search(text):
        weight *= 0.74
    if CONCESSION_RE.search(text) and not DIRECT_ASSAY_RE.search(text):
        weight *= 0.84
    if BACKGROUND_CONTEXT_RE.search(text):
        weight *= 0.86
    if COMPARISON_CONTEXT_RE.search(text):
        weight *= 1.03
    if DIRECT_ASSAY_RE.search(text):
        weight *= 1.06
    return round(max(0.55, min(1.08, weight)), 6)


def has_effective_negation(text: str) -> bool:
    normalized = normalize_text(text)
    if NON_NEGATING_NEGATION_RE.search(normalized) and not (NULL_RESULT_RE.search(normalized) or HEDGED_NEGATION_RE.search(normalized)):
        without_idiom = NON_NEGATING_NEGATION_RE.sub("", normalized)
        return bool(NEGATION_RE.search(without_idiom))
    return bool(NEGATION_RE.search(normalized))


def method_cues(text: str) -> list[str]:
    return [cue for cue, pattern in METHOD_CUE_PATTERNS.items() if pattern.search(text)]


def semantic_cues(text: str, section: str = "") -> list[str]:
    normalized_section = normalize_text(section)
    cues: list[str] = []
    if NON_NEGATING_NEGATION_RE.search(text):
        cues.append("non_negating_negation")
    if HEDGED_NEGATION_RE.search(text):
        cues.append("hedged_negation")
    if NULL_RESULT_RE.search(text):
        cues.append("null_result")
    if has_effective_negation(text):
        cues.append("negation")
    if CAUSAL_WEAKENING_RE.search(text):
        cues.append("causal_weakening")
    if ALTERNATIVE_EXPLANATION_RE.search(text):
        cues.extend(["alternative_explanation", "causal_weakening"])
    if UNCERTAINTY_RE.search(text):
        cues.append("uncertainty")
    if WEAK_OBSERVATION_RE.search(text):
        cues.append("weak_observation")
    if CONFLICT_RE.search(text):
        cues.append("conflict")
    if CONTEXT_BOUNDARY_RE.search(text):
        cues.append("context_boundary")
    if CONCESSION_RE.search(text):
        cues.append("concession")
    if BACKGROUND_CONTEXT_RE.search(text) or normalized_section in {"background", "review"}:
        cues.append("background")
    if COMPARISON_CONTEXT_RE.search(text):
        cues.append("comparison_context")
    if DIRECT_ASSAY_RE.search(text):
        cues.append("direct_assay")
    return cues


def assertion_semantics(text: str, section: str = "") -> dict[str, Any]:
    cues = set(semantic_cues(text, section))
    if "null_result" in cues:
        status = "contradict"
        reason = "null_result"
    elif {"conflict", "context_boundary"} & cues:
        status = "uncertain"
        reason = "conflict_or_context_boundary"
    elif "causal_weakening" in cues:
        status = "uncertain"
        reason = "causal_weakening"
    elif "negation" in cues and "hedged_negation" not in cues and "weak_observation" not in cues:
        status = "contradict"
        reason = "effective_negation"
    elif {"uncertainty", "weak_observation", "hedged_negation"} & cues:
        status = "uncertain"
        reason = "uncertain_or_weak_language"
    elif "concession" in cues and "direct_assay" not in cues and "comparison_context" not in cues:
        status = "uncertain"
        reason = "unsupported_concession"
    elif "background" in cues:
        status = "background"
        reason = "background_or_review_context"
    else:
        status = "support"
        reason = "direct_or_contextual_support"
    return {
        "support_status": status,
        "semantic_cues": sorted(cues),
        "method_cues": method_cues(text),
        "decision_reason": reason,
    }


def sentence_assertion_status(text: str, section: str = "") -> str:
    return str(assertion_semantics(text, section).get("support_status") or "uncertain")


def relation_probability(trigger: Trigger, subject: Mention, obj: Mention, section: str, text: str, mentions: list[Mention]) -> tuple[float, dict[str, Any]]:
    semantics = assertion_semantics(text, section)
    support_status = str(semantics["support_status"])
    components = {
        "base_rule_confidence": round(trigger.raw_score, 6),
        "section_weight": round(section_weight(section), 6),
        "entity_resolution_confidence": round(min(subject.resolution_confidence, obj.resolution_confidence), 6),
        "oncology_context_weight": round(oncology_context_weight(text, mentions), 6),
        "direction_clarity_weight": round(direction_weight(trigger.polarity), 6),
        "pair_distance_weight": round(pair_distance_weight(subject, obj), 6),
        "mention_density_weight": round(mention_density_weight(mentions), 6),
        "certainty_context_weight": round(certainty_context_weight(text), 6),
        "support_status": support_status,
        "support_status_reason": semantics["decision_reason"],
        "semantic_cues": semantics["semantic_cues"],
        "method_cues": semantics["method_cues"],
        "formula": "base_rule_confidence*section_weight*entity_resolution_confidence*oncology_context_weight*direction_clarity_weight*pair_distance_weight*mention_density_weight*certainty_context_weight",
    }
    prob = trigger.raw_score
    prob *= components["section_weight"]
    prob *= components["entity_resolution_confidence"]
    prob *= components["oncology_context_weight"]
    prob *= components["direction_clarity_weight"]
    prob *= components["pair_distance_weight"]
    prob *= components["mention_density_weight"]
    prob *= components["certainty_context_weight"]
    return clamp_unit(prob), components


def pair_distance(left: Mention, right: Mention) -> int:
    if left.end_offset <= right.start_offset:
        return right.start_offset - left.end_offset
    if right.end_offset <= left.start_offset:
        return left.start_offset - right.end_offset
    return 0


def relation_row(
    subject: Mention,
    predicate: str,
    obj: Mention,
    trigger: Trigger,
    text: str,
    mentions: list[Mention],
    min_prob: float,
) -> dict[str, Any] | None:
    if subject.entity_uid == obj.entity_uid:
        return None
    if pair_distance(subject, obj) > 240:
        return None
    prob, components = relation_probability(trigger, subject, obj, subject.section, text, mentions)
    if prob < min_prob:
        return None
    relation_uid = stable_uid(
        "litrel",
        subject.sentence_uid,
        subject.entity_uid,
        predicate,
        obj.entity_uid,
        trigger.rule_id,
        trigger.polarity,
    )
    license_id = subject.license_id if subject.license_id == obj.license_id else "mixed:" + content_hash(sorted([subject.license_id, obj.license_id]))[:12]
    source_release = subject.source_release if subject.source_release == obj.source_release else "mixed:" + content_hash(sorted([subject.source_release, obj.source_release]))[:12]
    return {
        "relation_uid": relation_uid,
        "subject_uid": subject.entity_uid,
        "subject_type": subject.entity_type,
        "predicate": predicate,
        "object_uid": obj.entity_uid,
        "object_type": obj.entity_type,
        "polarity": trigger.polarity,
        "cancer_context": mention_context(mentions),
        "sentence_uid": subject.sentence_uid,
        "pmid": subject.pmid,
        "pmcid": subject.pmcid,
        "section": subject.section,
        "trigger_phrase": trigger.phrase,
        "extraction_rule_id": trigger.rule_id,
        "support_status": sentence_assertion_status(text, subject.section),
        "raw_score": round(trigger.raw_score, 6),
        "calibrated_prob": round(prob, 6),
        "score_components_json": stable_json(components),
        "license_id": license_id,
        "source_release": source_release,
        "parser_hash": subject.parser_hash,
        "config_hash": subject.config_hash,
    }


def extract_relations(text: str, mentions: list[Mention], min_prob: float) -> list[dict[str, Any]]:
    if not mentions or oncology_context_weight(text, mentions) <= 0.0:
        return []
    by_type: dict[str, list[Mention]] = defaultdict(list)
    for mention in mentions:
        by_type[mention.entity_type].append(mention)
    rows: list[dict[str, Any]] = []

    change_trigger = find_trigger(text, CHANGE_PATTERNS)
    association_trigger = find_trigger(text, ASSOCIATION_PATTERNS)
    regulation_trigger = find_trigger(text, REGULATION_PATTERNS)

    if change_trigger:
        for metabolite in by_type.get("metabolite", []):
            for disease in by_type.get("disease", []):
                row = relation_row(metabolite, "metabolite_changed_in_cancer", disease, change_trigger, text, mentions, min_prob)
                if row:
                    rows.append(row)

    if association_trigger:
        for metabolite in by_type.get("metabolite", []):
            for disease in by_type.get("disease", []):
                row = relation_row(metabolite, "metabolite_associated_with_disease", disease, association_trigger, text, mentions, min_prob)
                if row:
                    rows.append(row)
        for target in by_type.get("target", []):
            for disease in by_type.get("disease", []):
                row = relation_row(target, "target_associated_with_disease", disease, association_trigger, text, mentions, min_prob)
                if row:
                    rows.append(row)
        for pathway in by_type.get("pathway", []):
            for disease in by_type.get("disease", []):
                row = relation_row(pathway, "pathway_associated_with_disease", disease, association_trigger, text, mentions, min_prob)
                if row:
                    rows.append(row)

    if regulation_trigger:
        for metabolite in by_type.get("metabolite", []):
            for gene in by_type.get("gene", []) + by_type.get("target", []):
                row = relation_row(metabolite, "metabolite_regulates_gene", gene, regulation_trigger, text, mentions, min_prob)
                if row:
                    rows.append(row)
        for gene in by_type.get("gene", []):
            for pathway in by_type.get("pathway", []):
                row = relation_row(gene, "gene_regulates_metabolic_process", pathway, regulation_trigger, text, mentions, min_prob)
                if row:
                    rows.append(row)
        for target in by_type.get("target", []):
            for pathway in by_type.get("pathway", []):
                row = relation_row(target, "gene_regulates_metabolic_process", pathway, regulation_trigger, text, mentions, min_prob)
                if row:
                    row["subject_type"] = "target"
                    rows.append(row)
        for target in by_type.get("target", []):
            for disease in by_type.get("disease", []):
                row = relation_row(target, "target_associated_with_disease", disease, regulation_trigger, text, mentions, min_prob)
                if row:
                    rows.append(row)

    deduped = {row["relation_uid"]: row for row in rows}
    return [deduped[key] for key in sorted(deduped)]


def read_existing_edges(normalized_dir: Path, relation_rows: list[dict[str, Any]]) -> dict[tuple[str, str, str], list[str]]:
    wanted: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for row in relation_rows:
        predicate = row["predicate"]
        if predicate == "target_associated_with_disease" and row["subject_type"] == "target":
            wanted["target_disease_edges"].add((row["subject_uid"], row["object_uid"]))
        elif predicate == "gene_regulates_metabolic_process" and row["subject_type"] == "gene":
            wanted["gene_pathway_edges"].add((row["subject_uid"], row["object_uid"]))

    existing: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    specs = {
        "target_disease_edges": ("target_uid", "disease_uid", "target_associated_with_disease"),
        "gene_pathway_edges": ("gene_uid", "pathway_uid", "gene_regulates_metabolic_process"),
    }
    for table_name, pairs in wanted.items():
        if not pairs:
            continue
        path = normalized_dir / f"{table_name}.parquet"
        if not path.exists():
            continue
        subject_col, object_col, predicate = specs[table_name]
        subject_values = sorted({pair[0] for pair in pairs})
        object_values = sorted({pair[1] for pair in pairs})
        dataset = table_dataset(path)
        filt = ds.field(subject_col).isin(subject_values) & ds.field(object_col).isin(object_values)
        table = dataset.to_table(columns=["edge_uid", subject_col, object_col], filter=filt)
        wanted_set = set(pairs)
        for row in table.to_pylist():
            pair = (row.get(subject_col, ""), row.get(object_col, ""))
            if pair in wanted_set:
                existing[(pair[0], predicate, pair[1])].append(row.get("edge_uid", ""))
    return {key: sorted(set(value)) for key, value in existing.items()}


def has_polarity_conflict(polarities: set[str]) -> bool:
    return any(POLARITY_OPPOSITES.get(polarity) in polarities for polarity in polarities)


def normalized_support_statuses(rows: list[dict[str, Any]]) -> set[str]:
    statuses = {str(row.get("support_status") or "support") for row in rows}
    return {status if status in SUPPORTED_ASSERTION_STATUSES else "uncertain" for status in statuses}


def evidence_article_key(row: dict[str, Any]) -> str:
    return str(row.get("pmid") or row.get("pmcid") or row.get("sentence_uid") or row.get("relation_uid") or "")


def article_clustered_probability(rows: list[dict[str, Any]]) -> tuple[float, list[float], dict[str, float]]:
    by_article: dict[str, float] = {}
    for row in rows:
        key = evidence_article_key(row)
        by_article[key] = max(by_article.get(key, 0.0), clamp_unit(row.get("calibrated_prob", 0.0)))
    article_probabilities = [round(value, 6) for _key, value in sorted(by_article.items())]
    return independent_probability_union(article_probabilities), article_probabilities, by_article


def aggregate_support(
    relation_rows: list[dict[str, Any]],
    normalized_dir: Path,
    min_support_prob: float,
    phash: str,
    cfg_hash: str,
) -> list[dict[str, Any]]:
    existing_edges = read_existing_edges(normalized_dir, relation_rows)
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in relation_rows:
        grouped[(row["subject_uid"], row["predicate"], row["object_uid"])].append(row)
    support_rows: list[dict[str, Any]] = []
    for key in sorted(grouped):
        rows = sorted(grouped[key], key=lambda row: row["relation_uid"])
        subject_uid, predicate, object_uid = key
        polarities = {row.get("polarity", "") for row in rows if row.get("polarity")}
        support_statuses = normalized_support_statuses(rows)
        sentence_probabilities = [round(float(row.get("calibrated_prob") or 0.0), 6) for row in rows]
        p_literature, article_probabilities, by_article = article_clustered_probability(rows)
        existing = existing_edges.get(key, [])
        has_direct_support = "support" in support_statuses
        if "contradict" in support_statuses or has_polarity_conflict(polarities):
            support_class = "conflict_candidate"
        elif not has_direct_support and support_statuses <= {"background"}:
            support_class = "background_candidate"
        elif not has_direct_support:
            support_class = "unresolved"
        elif existing:
            support_class = "support_direction" if any(polarity not in {"association"} for polarity in polarities) else "confirm"
        elif p_literature >= min_support_prob:
            support_class = "novel_candidate"
        else:
            support_class = "unresolved"
        sentence_uids = sorted({row.get("sentence_uid", "") for row in rows if row.get("sentence_uid")})
        pmids = sorted({row.get("pmid", "") for row in rows if row.get("pmid")})
        pmcids = sorted({row.get("pmcid", "") for row in rows if row.get("pmcid")})
        score_components = {
            "p_literature_formula": "1-prod(1-max_sentence_probability_per_article)",
            "sentence_probabilities": sentence_probabilities,
            "article_probabilities": article_probabilities,
            "article_count_for_probability": len(by_article),
            "support_statuses": sorted(support_statuses),
            "support_class_rule": "direct support can confirm existing edges or form novel candidates; contradiction conflicts; background/uncertain evidence cannot core-support claims",
            "min_support_prob": min_support_prob,
        }
        first = rows[0]
        support_rows.append(
            {
                "support_uid": stable_uid("litsup", subject_uid, predicate, object_uid),
                "subject_uid": subject_uid,
                "subject_type": first.get("subject_type", ""),
                "predicate": predicate,
                "object_uid": object_uid,
                "object_type": first.get("object_type", ""),
                "polarity_set": sorted(polarities),
                "support_status_set": sorted(support_statuses),
                "support_class": support_class,
                "supported_existing_edge_uids": existing,
                "evidence_relation_uids": [row["relation_uid"] for row in rows],
                "sentence_uids": sentence_uids,
                "pmids": pmids,
                "pmcids": pmcids,
                "evidence_sentence_count": len(sentence_uids),
                "distinct_article_count": len(set(pmids) | set(pmcids)),
                "p_literature": round(p_literature, 6),
                "raw_score_max": round(max(float(row.get("raw_score") or 0.0) for row in rows), 6),
                "calibrated_prob_max": round(max(float(row.get("calibrated_prob") or 0.0) for row in rows), 6),
                "score_components_json": stable_json(score_components),
                "license_id": first.get("license_id", ""),
                "source_release": first.get("source_release", ""),
                "parser_hash": phash,
                "config_hash": cfg_hash,
            }
        )
    return support_rows


def write_parquet(path: Path, rows: list[dict[str, Any]], schema: pa.Schema) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=schema), path)


def build_evidence(args: argparse.Namespace) -> dict[str, Any]:
    require_arrow()
    workspace = Path(args.workspace).resolve()
    normalized_root = (workspace / args.normalized_root).resolve()
    release_id = args.release_id or latest_release_id(normalized_root)
    normalized_dir = normalized_root / release_id
    output_dir = (workspace / args.output_root / release_id).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    phash = parser_hash()
    precision_filters, precision_config_path = load_precision_filters(workspace, args.precision_filter_config)
    precision_config_hash = content_hash(precision_filters)[:16]
    precision_stats = precision_filter_stats(precision_filters)
    config = {
        "builder_version": BUILDER_VERSION,
        "min_relation_prob": args.min_relation_prob,
        "min_support_prob": args.min_support_prob,
        "max_mentions_per_sentence": args.max_mentions_per_sentence,
        "max_sentences": args.max_sentences,
        "relevance_mode": args.relevance_mode,
        "relevance_model": args.relevance_model if args.relevance_mode == "pubmedbert" else "",
        "relevance_batch_size": args.relevance_batch_size if args.relevance_mode == "pubmedbert" else 0,
        "relevance_formula": "rule_relation_prob*(0.5+0.5*relevance_score) when relevance_mode=pubmedbert",
        "limit_to_graph_entities": not args.no_graph_entity_limit,
        "precision_filter_schema_version": precision_filters.get("schema_version", ""),
        "precision_filter_config_path": precision_config_path,
        "precision_filter_config_hash": precision_config_hash,
    }
    cfg_hash = content_hash(config)[:16]
    started = time.time()
    relevance_scorer = PubMedBertRelevanceScorer(args.relevance_model) if args.relevance_mode == "pubmedbert" else None

    lexicon, lexicon_stats = build_lexicon(
        normalized_dir,
        limit_to_graph_entities=not args.no_graph_entity_limit,
        precision_filters=precision_filters,
        precision_stats=precision_stats,
    )
    article_license_by_uid = load_article_license_map(normalized_dir)
    mention_writer = ParquetBatchWriter(output_dir / "sentence_mentions.parquet", MENTION_SCHEMA)
    relevance_writer = ParquetBatchWriter(output_dir / "sentence_relevance.parquet", RELEVANCE_SCHEMA)
    mention_batch: list[dict[str, Any]] = []
    relevance_batch: list[dict[str, Any]] = []
    relation_rows: list[dict[str, Any]] = []
    sentence_count = 0
    sentence_with_mentions = 0
    relevance_scored_count = 0
    relevance_decision_counter: Counter[str] = Counter()
    relation_counter: Counter[str] = Counter()
    mention_counter: Counter[str] = Counter()

    sentence_path = normalized_dir / "sentences.parquet"
    columns = ["sentence_uid", "article_uid", "pmid", "pmcid", "section", "sentence_text", "source_release", "parser_hash"]

    def process_sentence(row: dict[str, Any], relevance_score: float | None) -> None:
        nonlocal mention_batch, sentence_with_mentions
        mentions = match_mentions(row, lexicon, phash, cfg_hash, args.max_mentions_per_sentence)
        if mentions:
            sentence_with_mentions += 1
        for mention in mentions:
            mention_counter[mention.entity_type] += 1
            mention_batch.append(mention.as_row())
        if len(mention_batch) >= args.batch_size:
            mention_writer.write(mention_batch)
            mention_batch = []
        extracted = extract_relations(str(row.get("sentence_text") or ""), mentions, args.min_relation_prob)
        if relevance_score is not None:
            extracted = [weighted for relation in extracted if (weighted := apply_relevance_weight(relation, relevance_score, args.min_relation_prob)) is not None]
        for relation in extracted:
            relation_counter[relation["predicate"]] += 1
        relation_rows.extend(extracted)

    def process_relevance_rows(rows: list[dict[str, Any]]) -> None:
        nonlocal relevance_batch, relevance_scored_count
        if not rows:
            return
        assert relevance_scorer is not None
        for row, (rel_row, relevance_score) in zip(rows, relevance_rows(rows, relevance_scorer, cfg_hash)):
            relevance_batch.append(rel_row)
            relevance_scored_count += 1
            relevance_decision_counter[rel_row["decision"]] += 1
            if len(relevance_batch) >= args.batch_size:
                relevance_writer.write(relevance_batch)
                relevance_batch = []
            process_sentence(row, relevance_score)

    pending_relevance_rows: list[dict[str, Any]] = []
    for row in iter_table_rows(sentence_path, columns=columns, batch_size=args.batch_size):
        if args.max_sentences > 0 and sentence_count >= args.max_sentences:
            break
        sentence_count += 1
        row["_article_license_id"] = article_license_by_uid.get(str(row.get("article_uid") or ""), "local_articles:articles_collect")
        if relevance_scorer is None:
            process_sentence(row, None)
            continue
        pending_relevance_rows.append(row)
        if len(pending_relevance_rows) >= args.relevance_batch_size:
            process_relevance_rows(pending_relevance_rows)
            pending_relevance_rows = []
    if relevance_scorer is not None:
        process_relevance_rows(pending_relevance_rows)
    mention_writer.write(mention_batch)
    mention_writer.close()
    relevance_writer.write(relevance_batch)
    relevance_writer.close()

    relation_rows.sort(key=lambda row: row["relation_uid"])
    write_parquet(output_dir / "relation_candidates.parquet", relation_rows, RELATION_SCHEMA)
    support_rows = aggregate_support(relation_rows, normalized_dir, args.min_support_prob, phash, cfg_hash)
    support_rows.sort(key=lambda row: row["support_uid"])
    write_parquet(output_dir / "literature_edge_support.parquet", support_rows, SUPPORT_SCHEMA)

    support_counter = Counter(row["support_class"] for row in support_rows)
    normalized_manifest = read_json(normalized_dir / "normalized_manifest.json")
    articles_dir = (workspace / args.articles_dir).resolve()
    metrics = {
        "sentence_count_scanned": sentence_count,
        "sentence_with_mentions": sentence_with_mentions,
        "sentence_mention_count": mention_writer.rows,
        "normalized_mention_rate": round(1.0 if mention_writer.rows else 0.0, 6),
        "relation_candidate_count": len(relation_rows),
        "evidence_candidate_count": len(relation_rows),
        "edge_support_count": len(support_rows),
        "relevance_scored_sentence_count": relevance_scored_count,
        "relevance_count_by_decision": dict(sorted(relevance_decision_counter.items())),
        "supported_existing_edge_count": support_counter["confirm"] + support_counter["support_direction"],
        "novel_candidate_count": support_counter["novel_candidate"],
        "conflict_candidate_count": support_counter["conflict_candidate"],
        "unresolved_candidate_count": support_counter["unresolved"],
        "mention_count_by_type": dict(sorted(mention_counter.items())),
        "relation_count_by_predicate": dict(sorted(relation_counter.items())),
        "support_count_by_class": dict(sorted(support_counter.items())),
    }
    manifest = {
        "builder_version": BUILDER_VERSION,
        "release_id": release_id,
        "source_release": release_id,
        "normalized_manifest": str(normalized_dir / "normalized_manifest.json"),
        "normalized_manifest_hash": content_hash(normalized_manifest)[:16] if normalized_manifest else "",
        "articles_dir": str(articles_dir),
        "articles_collect_file_count": len(list(articles_dir.glob("*.txt"))) if articles_dir.exists() else 0,
        "article_table_row_count": len(article_license_by_uid),
        "output_dir": str(output_dir),
        "config": config,
        "parser_hash": phash,
        "config_hash": cfg_hash,
        "metrics": metrics,
        "lexicon": lexicon_stats,
        "precision_filters": finalize_precision_stats(precision_stats, precision_filters, precision_config_path, precision_config_hash),
        "tables": [
            {"table": "sentence_mentions", "path": str(output_dir / "sentence_mentions.parquet"), "rows": mention_writer.rows},
            {"table": "sentence_relevance", "path": str(output_dir / "sentence_relevance.parquet"), "rows": relevance_writer.rows},
            {"table": "relation_candidates", "path": str(output_dir / "relation_candidates.parquet"), "rows": len(relation_rows)},
            {"table": "literature_edge_support", "path": str(output_dir / "literature_edge_support.parquet"), "rows": len(support_rows)},
        ],
        "seconds": round(time.time() - started, 3),
    }
    manifest["manifest_hash"] = content_hash({key: value for key, value in manifest.items() if key != "seconds"})
    (output_dir / "literature_evidence_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, sort_keys=True, indent=2), encoding="utf-8")
    return manifest


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Phase 2-lite literature evidence overlay from normalized sentences.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--normalized-root", default=DEFAULT_NORMALIZED_ROOT)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--articles-dir", default=DEFAULT_ARTICLES_DIR)
    parser.add_argument("--release-id", default="")
    parser.add_argument("--max-sentences", type=int, default=0, help="Limit sentence rows for smoke tests. 0 means all.")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--max-mentions-per-sentence", type=int, default=DEFAULT_MAX_MENTIONS_PER_SENTENCE)
    parser.add_argument("--min-relation-prob", type=float, default=DEFAULT_MIN_RELATION_PROB)
    parser.add_argument("--min-support-prob", type=float, default=DEFAULT_MIN_SUPPORT_PROB)
    parser.add_argument("--precision-filter-config", default=DEFAULT_PRECISION_FILTER_CONFIG)
    parser.add_argument("--no-graph-entity-limit", action="store_true", help="Allow all canonical entities into the mention lexicon.")
    parser.add_argument("--relevance-mode", choices=["none", "pubmedbert"], default="none", help="Optional sentence relevance scorer for filtering/reranking.")
    parser.add_argument("--relevance-model", default=DEFAULT_RELEVANCE_MODEL, help="Transformer model name used when --relevance-mode pubmedbert.")
    parser.add_argument("--relevance-batch-size", type=int, default=DEFAULT_RELEVANCE_BATCH_SIZE, help="Sentence batch size for PubMedBERT CPU/GPU inference.")
    args = parser.parse_args(argv)
    if args.relevance_batch_size < 1:
        parser.error("--relevance-batch-size must be >= 1")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    manifest = build_evidence(args)
    print(
        json.dumps(
            {
                "release_id": manifest["release_id"],
                "output_dir": manifest["output_dir"],
                "metrics": manifest["metrics"],
                "manifest_hash": manifest["manifest_hash"],
            },
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
