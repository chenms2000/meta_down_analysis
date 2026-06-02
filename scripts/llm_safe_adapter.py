"""Local deterministic safety layer for LLM-facing explanations.

This module intentionally does not call an external LLM. It defines the
adapter contract, produces a small deterministic explanation bundle, and
guards any adapter-like output against unsupported or mutating claims.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import re
import socket
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib import request as urllib_request
from urllib.error import HTTPError, URLError


INPUT_CONTRACT_VERSION = "llm_safe_adapter.input.v1"
OUTPUT_CONTRACT_VERSION = "llm_safe_adapter.output.v1"
ADAPTER_VERSION = "llm_safe_adapter.local.20260513"
EXTERNAL_ADAPTER_VERSION = "llm_safe_adapter.external.20260513"
LLM_NARRATOR_PROMPT_VERSION = "llm_safe_adapter.prompt.20260513"
DEFAULT_OPENAI_COMPATIBLE_CHAT_ENDPOINT = "https://api.openai.com/v1/chat/completions"

ALLOWED_READ_ENDPOINTS = {
    "/analyze/metabolites",
    "/analyze/differential-table",
    "/evidence",
    "/subgraph",
    "/releases",
}

SOURCE_TYPES = {
    "analysis_pack",
    "evidence",
    "subgraph",
    "release",
}

EXPLANATION_FIELDS = {
    "narrative_summary",
    "evidence_digest",
    "question_to_spec",
    "candidate_suggestion",
    "insufficient_evidence",
}

ALLOWED_TOP_LEVEL_OUTPUT_FIELDS = {
    "contract_version",
    "adapter_version",
    "status",
    "guard",
    "determinism",
    *EXPLANATION_FIELDS,
}

ALLOWED_STATUSES = {"ok", "insufficient_evidence", "blocked_by_guard"}

REFERENCE_SCALAR_KEYS = {
    "analysis_pack_hash",
    "edge_uid",
    "entity_uid",
    "input_id",
    "mention_uid",
    "metabolite_uid",
    "node_uid",
    "path_id",
    "pathway_uid",
    "pmcid",
    "pmid",
    "primary_external_id",
    "reaction_uid",
    "relation_uid",
    "release_id",
    "sentence_uid",
    "subject_uid",
    "support_uid",
    "target_uid",
    "terminal_node_uid",
    "object_uid",
    "disease_uid",
    "gene_uid",
    "protein_uid",
    "manifest_hash",
}

REFERENCE_LIST_KEYS = {
    "duplicate_matched_metabolite_uids",
    "edge_uids",
    "evidence_ref_uids",
    "evidence_relation_uids",
    "mention_uids",
    "node_uids",
    "pmcids",
    "pmids",
    "relation_uids",
    "sentence_uids",
    "source_records",
    "supported_existing_edge_uids",
}

FORBIDDEN_OUTPUT_KEYS = {
    "calibrated_prob",
    "canonical_entity_uid",
    "curated_graph_write",
    "entity_resolution_decision",
    "graph_mutation",
    "new_edge",
    "new_edges",
    "new_entity",
    "new_entities",
    "new_fact",
    "new_facts",
    "new_p_final",
    "new_score",
    "new_scores",
    "p_final",
    "raw_score",
    "resolved_entity_uid",
    "resolver_decision",
    "score",
    "score_components",
    "write_actions",
}

FORBIDDEN_KEY_SUBSTRINGS = (
    "curated_graph",
    "entity_resolution_decision",
    "graph_write",
    "new_score",
    "p_final",
)

FORBIDDEN_TEXT_PATTERNS = [
    (
        "mutates_scores_or_graph",
        re.compile(
            r"\b(update|updated|write|overwrite|change|changed|set|raise|raised|lower|lowered|modify|modified|mutate|mutated|save|saved)\b"
            r".{0,80}\b(p_final|score|curated graph|graph|release|edge)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "entity_resolution_decision",
        re.compile(
            r"\b(i|we|adapter|llm)\s+(resolved|merged|mapped|canonicalized|selected)\b"
            r"|\bshould\s+(resolve|merge|map|canonicalize|select)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "overstrong_claim",
        re.compile(
            r"\b(proves|proof|proved|definitively|guarantees|causes|cures|is causal|confirmed new fact|"
            r"activated|activation of|key driver|effective treatment|patients should receive)\b"
            r"|证明了|证明|表明该疾病|说明该药物有效|关键驱动|该通路被激活|预测患者适合|适合某治疗",
            re.IGNORECASE,
        ),
    ),
    (
        "curation_action",
        re.compile(r"\b(add|insert|write|curate|commit)\b.{0,80}\b(edge|entity|fact|graph|source of truth)\b", re.IGNORECASE),
    ),
]

ID_TOKEN_RE = re.compile(
    r"\b(?:met|metabolite|gene|protein|target|disease|pathway|reaction|edge|litsup|litrel|sent|mention)_[A-Za-z0-9][A-Za-z0-9_:-]*\b"
    r"|\bPMID[: ]?\d+\b"
    r"|\bPMC\d+\b"
)

NUMERIC_CLAIM_RE = re.compile(r"(?<![A-Za-z_])\d+(?:\.\d+)?%?")
LANGUAGE_KEYS = {"text", "reason", "suggestion", "digest", "narrative", "claim"}
IGNORED_OUTPUT_PATH_PREFIXES = ("$.guard", "$.determinism")
STRUCTURED_PREDICTION_GROUPS = (
    "high_confidence_themes",
    "medium_confidence_themes",
    "exploratory_themes",
    "downgraded_but_supported_themes",
    "pathway_evidence",
    "target_predictions",
    "disease_predictions",
    "drug_hypotheses",
    "context_mismatch_results",
)
PREDICTION_REF_REQUIRED_FIELDS = {
    "prediction_id",
    "confidence_tier",
    "evidence_refs",
    "claim_refs",
    "appendix",
    "research_only",
}


class LLMAdapterError(RuntimeError):
    def __init__(self, code: str, message: str, detail: dict[str, Any] | None = None):
        super().__init__(message)
        self.code = code
        self.detail = detail or {}

    def as_issue(self) -> dict[str, Any]:
        return {"code": self.code, "path": "$.adapter_backend", "detail": {"message": str(self), **self.detail}}


@dataclass(frozen=True)
class ExternalLLMConfig:
    enabled: bool = False
    provider: str = "openai_compatible"
    endpoint: str = DEFAULT_OPENAI_COMPATIBLE_CHAT_ENDPOINT
    model: str = ""
    api_key: str = ""
    timeout_seconds: float = 30.0
    max_output_tokens: int = 1600
    proxy_url: str = ""

    @classmethod
    def from_env(
        cls,
        enabled: bool | None = None,
        provider: str = "",
        endpoint: str = "",
        model: str = "",
        api_key_env: str = "",
        timeout_seconds: float | None = None,
        max_output_tokens: int | None = None,
        proxy_url: str = "",
    ) -> "ExternalLLMConfig":
        env_enabled = str(os.environ.get("LLM_SAFE_ADAPTER_ENABLE_EXTERNAL", "")).casefold() in {"1", "true", "yes", "on"}
        key_env = api_key_env or os.environ.get("LLM_SAFE_ADAPTER_API_KEY_ENV", "LLM_SAFE_ADAPTER_API_KEY")
        api_key = os.environ.get(key_env) or os.environ.get("OPENAI_API_KEY", "")
        timeout = timeout_seconds
        if timeout is None:
            try:
                timeout = float(os.environ.get("LLM_SAFE_ADAPTER_TIMEOUT_SECONDS", "30"))
            except ValueError:
                timeout = 30.0
        tokens = max_output_tokens
        if tokens is None:
            try:
                tokens = int(os.environ.get("LLM_SAFE_ADAPTER_MAX_OUTPUT_TOKENS", "1600"))
            except ValueError:
                tokens = 1600
        return cls(
            enabled=env_enabled if enabled is None else bool(enabled),
            provider=provider or os.environ.get("LLM_SAFE_ADAPTER_PROVIDER", "openai_compatible"),
            endpoint=endpoint or os.environ.get("LLM_SAFE_ADAPTER_ENDPOINT", DEFAULT_OPENAI_COMPATIBLE_CHAT_ENDPOINT),
            model=model or os.environ.get("LLM_SAFE_ADAPTER_MODEL", ""),
            api_key=api_key,
            timeout_seconds=float(timeout),
            max_output_tokens=int(tokens),
            proxy_url=proxy_url or os.environ.get("LLM_SAFE_ADAPTER_PROXY", ""),
        )

    def validate_ready(self) -> None:
        if not self.enabled:
            raise LLMAdapterError("external_llm_disabled", "External LLM backend is disabled.")
        if self.provider != "openai_compatible":
            raise LLMAdapterError("unsupported_llm_provider", "Only openai_compatible provider is supported.", {"provider": self.provider})
        if not self.endpoint:
            raise LLMAdapterError("missing_llm_endpoint", "External LLM endpoint is not configured.")
        if not self.model:
            raise LLMAdapterError("missing_llm_model", "External LLM model is not configured.")
        if not self.api_key:
            raise LLMAdapterError("missing_llm_api_key", "External LLM API key is not configured.")
        try:
            self.api_key.encode("latin-1")
        except UnicodeEncodeError as exc:
            raise LLMAdapterError(
                "invalid_llm_api_key_characters",
                "External LLM API key contains non-HTTP-header characters.",
            ) from exc

    def sanitized(self) -> dict[str, Any]:
        return {
            "enabled": self.enabled,
            "provider": self.provider,
            "endpoint": self.endpoint,
            "model": self.model,
            "timeout_seconds": self.timeout_seconds,
            "max_output_tokens": self.max_output_tokens,
            "proxy_configured": bool(self.proxy_url),
            "proxy_url": self.proxy_url,
            "api_key_configured": bool(self.api_key),
        }


@dataclass(frozen=True)
class BackendResult:
    output: dict[str, Any]
    audit: dict[str, Any]


def stable_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def content_hash(value: Any) -> str:
    return hashlib.sha256(stable_json(value).encode("utf-8")).hexdigest()


def guard_policy_hash() -> str:
    return content_hash(
        {
            "allowed_read_endpoints": sorted(ALLOWED_READ_ENDPOINTS),
            "allowed_output_fields": sorted(ALLOWED_TOP_LEVEL_OUTPUT_FIELDS),
            "forbidden_output_keys": sorted(FORBIDDEN_OUTPUT_KEYS),
            "forbidden_key_substrings": sorted(FORBIDDEN_KEY_SUBSTRINGS),
            "forbidden_text_patterns": [name for name, _pattern in FORBIDDEN_TEXT_PATTERNS],
        }
    )


def is_scalar(value: Any) -> bool:
    return value is None or isinstance(value, (str, int, float, bool))


def walk_json(value: Any, path: str = "$"):
    yield path, value
    if isinstance(value, dict):
        for key, child in value.items():
            yield from walk_json(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from walk_json(child, f"{path}[{index}]")


def source_payloads(input_pack: dict[str, Any]) -> dict[str, Any]:
    analysis_pack = input_pack.get("analysis_pack")
    evidence = input_pack.get("evidence")
    subgraph = input_pack.get("subgraph")
    release = input_pack.get("release")
    if not release and isinstance(analysis_pack, dict):
        release = analysis_pack.get("release")
    return {
        "analysis_pack": analysis_pack if isinstance(analysis_pack, dict) else {},
        "evidence": evidence if isinstance(evidence, dict) else {},
        "subgraph": subgraph if isinstance(subgraph, dict) else {},
        "release": release if isinstance(release, dict) else {},
    }


def add_catalog_ref(catalog: dict[str, Any], source_type: str, ref_id: Any, path: str) -> None:
    if ref_id is None:
        return
    text = str(ref_id).strip()
    if not text:
        return
    catalog[source_type]["ref_ids"].add(text)
    catalog[source_type]["path_by_ref"][text].add(path)
    if source_type == "release" and not text.startswith("release:"):
        release_ref = f"release:{text}"
        catalog[source_type]["ref_ids"].add(release_ref)
        catalog[source_type]["path_by_ref"][release_ref].add(path)


def build_source_catalog(input_pack: dict[str, Any]) -> dict[str, Any]:
    catalog: dict[str, Any] = {
        source_type: {"paths": set(), "ref_ids": set(), "path_by_ref": defaultdict(set)}
        for source_type in SOURCE_TYPES
    }
    for source_type, payload in source_payloads(input_pack).items():
        root_path = f"$.{source_type}"
        catalog[source_type]["paths"].add(root_path)
        add_catalog_ref(catalog, source_type, f"{source_type}:root", root_path)
        for path, value in walk_json(payload, root_path):
            catalog[source_type]["paths"].add(path)
            key = path.rsplit(".", 1)[-1]
            if "[" in key:
                key = key.split("[", 1)[0]
            if key in REFERENCE_SCALAR_KEYS and is_scalar(value):
                add_catalog_ref(catalog, source_type, value, path)
                if key == "pmid" and value:
                    add_catalog_ref(catalog, source_type, f"PMID:{value}", path)
                if key == "pmcid" and value:
                    add_catalog_ref(catalog, source_type, f"PMCID:{value}", path)
                if key == "release_id" and value:
                    add_catalog_ref(catalog, source_type, f"release:{value}", path)
            elif key in REFERENCE_LIST_KEYS and isinstance(value, list):
                for item in value:
                    if is_scalar(item):
                        add_catalog_ref(catalog, source_type, item, path)
                        if key == "pmids" and item:
                            add_catalog_ref(catalog, source_type, f"PMID:{item}", path)
                        if key == "pmcids" and item:
                            add_catalog_ref(catalog, source_type, f"PMCID:{item}", path)
    return catalog


def all_catalog_refs(catalog: dict[str, Any]) -> set[str]:
    refs: set[str] = set()
    for source in catalog.values():
        refs.update(source["ref_ids"])
    return refs


def ref_path(catalog: dict[str, Any], source_type: str, ref_id: str, fallback_path: str) -> str:
    paths = sorted(catalog.get(source_type, {}).get("path_by_ref", {}).get(ref_id, []))
    return paths[0] if paths else fallback_path


def source_ref(catalog: dict[str, Any], source_type: str, ref_id: str, fallback_path: str) -> dict[str, str]:
    return {
        "source_type": source_type,
        "ref_id": ref_id,
        "path": ref_path(catalog, source_type, ref_id, fallback_path),
    }


def release_source_ref(input_pack: dict[str, Any], catalog: dict[str, Any]) -> dict[str, str]:
    release = source_payloads(input_pack).get("release", {})
    release_id = str(release.get("release_id") or "root")
    ref_id = f"release:{release_id}" if release_id != "root" else "release:root"
    if ref_id not in catalog["release"]["ref_ids"]:
        ref_id = "release:root"
    return source_ref(catalog, "release", ref_id, "$.release")


def validate_adapter_input(input_pack: dict[str, Any]) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    if input_pack.get("contract_version") != INPUT_CONTRACT_VERSION:
        issues.append(
            {
                "code": "input_contract_version_mismatch",
                "path": "$.contract_version",
                "detail": {
                    "expected": INPUT_CONTRACT_VERSION,
                    "actual": input_pack.get("contract_version", ""),
                },
            }
        )
    for index, row in enumerate(input_pack.get("source_endpoints", []) or []):
        endpoint = row.get("endpoint") if isinstance(row, dict) else row
        if endpoint not in ALLOWED_READ_ENDPOINTS:
            issues.append(
                {
                    "code": "unauthorized_input_endpoint",
                    "path": f"$.source_endpoints[{index}]",
                    "detail": {"endpoint": endpoint, "allowed": sorted(ALLOWED_READ_ENDPOINTS)},
                }
            )
    for key in input_pack:
        if key.startswith("resolve") or key in {"resolver", "entity_resolution", "curated_graph"}:
            issues.append({"code": "unauthorized_input_payload", "path": f"$.{key}", "detail": {"key": key}})
    return issues


def path_is_ignored(path: str) -> bool:
    return any(path == prefix or path.startswith(f"{prefix}.") or path.startswith(f"{prefix}[") for prefix in IGNORED_OUTPUT_PATH_PREFIXES)


def iter_language_blocks(output: dict[str, Any]):
    for path, value in walk_json(output):
        if path_is_ignored(path):
            continue
        if not isinstance(value, dict):
            continue
        texts = []
        for key in LANGUAGE_KEYS:
            text = value.get(key)
            if isinstance(text, str) and text.strip():
                texts.append((key, text))
        for key, text in texts:
            yield f"{path}.{key}", text, value.get("source_refs")


def normalize_id_token(token: str) -> list[str]:
    token = token.strip()
    if token.upper().startswith("PMID "):
        return [token, f"PMID:{token.split(None, 1)[1]}"]
    if token.upper().startswith("PMID:"):
        return [token, token.split(":", 1)[1]]
    if token.upper().startswith("PMC"):
        return [token, f"PMCID:{token}"]
    return [token]


def validate_source_refs(source_refs: Any, path: str, catalog: dict[str, Any]) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    if not isinstance(source_refs, list) or not source_refs:
        return [{"code": "language_without_source_refs", "path": path, "detail": {}}]
    for index, ref in enumerate(source_refs):
        ref_path_value = f"{path}.source_refs[{index}]"
        if not isinstance(ref, dict):
            issues.append({"code": "invalid_source_ref", "path": ref_path_value, "detail": {"reason": "source ref must be an object"}})
            continue
        source_type = ref.get("source_type")
        ref_id = str(ref.get("ref_id") or "")
        cited_path = ref.get("path")
        if source_type not in SOURCE_TYPES:
            issues.append({"code": "invalid_source_type", "path": ref_path_value, "detail": {"source_type": source_type}})
            continue
        if not ref_id:
            issues.append({"code": "missing_source_ref_id", "path": ref_path_value, "detail": {"source_type": source_type}})
            continue
        if ref_id not in catalog[source_type]["ref_ids"]:
            issues.append(
                {
                    "code": "source_ref_not_in_input",
                    "path": ref_path_value,
                    "detail": {"source_type": source_type, "ref_id": ref_id},
                }
            )
        if cited_path and cited_path not in catalog[source_type]["paths"]:
            issues.append(
                {
                    "code": "source_ref_path_not_in_input",
                    "path": ref_path_value,
                    "detail": {"source_type": source_type, "path": cited_path},
                }
            )
    return issues


def validate_language_text(text: str, path: str, catalog: dict[str, Any], input_text: str) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    all_refs = all_catalog_refs(catalog)
    for issue_code, pattern in FORBIDDEN_TEXT_PATTERNS:
        if pattern.search(text):
            issues.append({"code": issue_code, "path": path, "detail": {"text": text}})
    for token in ID_TOKEN_RE.findall(text):
        alternatives = normalize_id_token(token)
        if not any(alt in all_refs for alt in alternatives):
            issues.append({"code": "text_ref_not_in_input", "path": path, "detail": {"ref_id": token}})
    for token in NUMERIC_CLAIM_RE.findall(text):
        normalized = token.rstrip("%")
        if len(normalized.replace(".", "")) < 2 and "." not in normalized and "%" not in token:
            continue
        if token not in input_text and normalized not in input_text:
            issues.append({"code": "unsupported_numeric_claim", "path": path, "detail": {"number": token}})
    return issues


def prediction_row_id(row: dict[str, Any], group: str, index: int) -> str:
    for key in (
        "prediction_id",
        "theme_id",
        "state_id",
        "pathway_id",
        "pathway_uid",
        "target_uid",
        "disease_uid",
        "drug_id",
        "result_id",
        "result_uid",
    ):
        value = str(row.get(key) or "").strip()
        if value:
            return value
    name = str(row.get("display_name") or row.get("drug_name") or "").strip()
    if name:
        return f"{group}:{name}"
    return f"{group}:{index}"


def iter_structured_prediction_rows(input_pack: dict[str, Any]) -> list[dict[str, Any]]:
    structured = source_payloads(input_pack)["analysis_pack"].get("structured_prediction") or {}
    if not isinstance(structured, dict):
        return []
    rows: list[dict[str, Any]] = []
    for group in STRUCTURED_PREDICTION_GROUPS:
        values = structured.get(group) or []
        if not isinstance(values, list):
            continue
        for index, value in enumerate(values):
            if not isinstance(value, dict):
                continue
            row = copy.deepcopy(value)
            row["_structured_group"] = group
            row["_prediction_id"] = prediction_row_id(row, group, index)
            rows.append(row)
    return rows


def structured_prediction_catalog(input_pack: dict[str, Any]) -> dict[str, Any]:
    rows = iter_structured_prediction_rows(input_pack)
    by_id: dict[str, dict[str, Any]] = {}
    display_names: set[str] = set()
    for row in rows:
        row_id = str(row.get("_prediction_id") or "")
        if row_id:
            by_id[row_id] = row
        for key in ("display_name", "drug_name", "pathway_name"):
            name = str(row.get(key) or "").strip()
            if name:
                display_names.add(name.casefold())
    return {"rows": rows, "by_id": by_id, "display_names": display_names}


def prediction_ref_from_row(row: dict[str, Any]) -> dict[str, Any]:
    allowed_evidence_keys = {
        "ref_type",
        "edge_uid",
        "support_uid",
        "support_class",
        "pmids",
        "pmcids",
        "sentence_uids",
        "relation_uids",
        "source_name",
        "source_record_id",
        "source_release",
        "license_id",
    }
    compact_evidence_refs = [
        {key: value for key, value in ref.items() if key in allowed_evidence_keys}
        for ref in (row.get("evidence_refs", []) or [])
        if isinstance(ref, dict)
    ]
    return {
        "prediction_id": row.get("_prediction_id") or prediction_row_id(row, str(row.get("_structured_group") or "prediction"), 0),
        "display_name": row.get("display_name") or row.get("drug_name") or row.get("pathway_name") or "",
        "confidence_tier": row.get("confidence_tier") or row.get("confidence") or row.get("display_confidence") or "low",
        "evidence_refs": compact_evidence_refs,
        "claim_refs": row.get("claim_refs", {}) or {},
        "appendix": bool(row.get("appendix")),
        "appendix_reason": row.get("appendix_reason", ""),
        "research_only": bool(row.get("research_only")),
    }


def validate_prediction_refs(block: dict[str, Any], path: str, prediction_catalog: dict[str, Any]) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    rows_by_id = prediction_catalog.get("by_id", {})
    if not rows_by_id:
        return issues
    reported_names = [str(item).strip().casefold() for item in block.get("reported_display_names", []) or [] if str(item).strip()]
    reported_ids = [str(item).strip() for item in block.get("reported_prediction_ids", []) or [] if str(item).strip()]
    prediction_refs = block.get("prediction_refs")
    has_prediction_claim = bool(reported_names or reported_ids or block.get("claim"))
    if not has_prediction_claim:
        return issues
    if not isinstance(prediction_refs, list) or not prediction_refs:
        issues.append({"code": "prediction_claim_without_structured_prediction_refs", "path": path, "detail": {}})
        return issues

    allowed_names = prediction_catalog.get("display_names", set())
    for name in reported_names:
        if name not in allowed_names:
            issues.append({"code": "unsupported_claim_not_in_structured_prediction", "path": path, "detail": {"display_name": name}})
    for prediction_id in reported_ids:
        if prediction_id not in rows_by_id:
            issues.append({"code": "unsupported_claim_not_in_structured_prediction", "path": path, "detail": {"prediction_id": prediction_id}})

    for index, ref in enumerate(prediction_refs):
        ref_path = f"{path}.prediction_refs[{index}]"
        if not isinstance(ref, dict):
            issues.append({"code": "invalid_prediction_ref", "path": ref_path, "detail": {"reason": "prediction ref must be an object"}})
            continue
        missing = sorted(field for field in PREDICTION_REF_REQUIRED_FIELDS if field not in ref)
        if missing:
            issues.append({"code": "prediction_ref_missing_required_fields", "path": ref_path, "detail": {"missing": missing}})
        prediction_id = str(ref.get("prediction_id") or "")
        source_row = rows_by_id.get(prediction_id)
        if not source_row:
            issues.append({"code": "prediction_ref_not_in_structured_prediction", "path": ref_path, "detail": {"prediction_id": prediction_id}})
            continue
        expected_tier = source_row.get("confidence_tier") or source_row.get("confidence") or source_row.get("display_confidence")
        if expected_tier and ref.get("confidence_tier") != expected_tier:
            issues.append(
                {
                    "code": "prediction_ref_confidence_tier_mismatch",
                    "path": ref_path,
                    "detail": {"expected": expected_tier, "actual": ref.get("confidence_tier")},
                }
            )
        if not isinstance(ref.get("evidence_refs"), list) or not ref.get("evidence_refs"):
            issues.append({"code": "prediction_ref_missing_evidence_refs", "path": ref_path, "detail": {"prediction_id": prediction_id}})
        if not isinstance(ref.get("claim_refs"), dict) or not ref.get("claim_refs"):
            issues.append({"code": "prediction_ref_missing_claim_refs", "path": ref_path, "detail": {"prediction_id": prediction_id}})
    return issues


def guard_adapter_output(input_pack: dict[str, Any], output: dict[str, Any]) -> dict[str, Any]:
    catalog = build_source_catalog(input_pack)
    prediction_catalog = structured_prediction_catalog(input_pack)
    input_text = stable_json(input_pack)
    issues = validate_adapter_input(input_pack)

    if output.get("contract_version") != OUTPUT_CONTRACT_VERSION:
        issues.append(
            {
                "code": "output_contract_version_mismatch",
                "path": "$.contract_version",
                "detail": {
                    "expected": OUTPUT_CONTRACT_VERSION,
                    "actual": output.get("contract_version", ""),
                },
            }
        )
    if output.get("status") not in ALLOWED_STATUSES:
        issues.append({"code": "invalid_output_status", "path": "$.status", "detail": {"status": output.get("status")}})

    unknown_top = sorted(set(output) - ALLOWED_TOP_LEVEL_OUTPUT_FIELDS)
    for key in unknown_top:
        issues.append({"code": "unauthorized_output_field", "path": f"$.{key}", "detail": {"field": key}})

    for path, value in walk_json(output):
        if path_is_ignored(path):
            continue
        if not isinstance(value, dict):
            continue
        for key in value:
            normalized_key = str(key).casefold()
            if normalized_key in FORBIDDEN_OUTPUT_KEYS or any(fragment in normalized_key for fragment in FORBIDDEN_KEY_SUBSTRINGS):
                issues.append({"code": "forbidden_output_field", "path": f"{path}.{key}", "detail": {"field": key}})
        issues.extend(validate_prediction_refs(value, path, prediction_catalog))

    language_block_count = 0
    for path, text, source_refs in iter_language_blocks(output):
        language_block_count += 1
        issues.extend(validate_source_refs(source_refs, path.rsplit(".", 1)[0], catalog))
        issues.extend(validate_language_text(text, path, catalog, input_text))

    if output.get("status") == "insufficient_evidence" and "insufficient_evidence" not in output:
        issues.append({"code": "missing_insufficient_evidence_payload", "path": "$.insufficient_evidence", "detail": {}})
    if output.get("status") == "ok" and not any(field in output for field in ("narrative_summary", "evidence_digest", "question_to_spec", "candidate_suggestion")):
        issues.append({"code": "missing_explanation_payload", "path": "$", "detail": {}})

    return {
        "passed": not issues,
        "guard_version": ADAPTER_VERSION,
        "policy_hash": guard_policy_hash(),
        "language_block_count": language_block_count,
        "issue_count": len(issues),
        "issues": issues,
    }


def has_explainable_material(input_pack: dict[str, Any]) -> bool:
    payloads = source_payloads(input_pack)
    analysis = payloads["analysis_pack"]
    evidence = payloads["evidence"]
    subgraph = payloads["subgraph"]
    if any(analysis.get(key) for key in ("input_summary", "matched", "ambiguous", "unmatched", "pathway_rankings", "target_rankings", "disease_rankings")):
        return True
    if any(evidence.get(key) for key in ("support", "relation_candidates", "sentences", "mentions")):
        return True
    if any(subgraph.get(key) for key in ("nodes", "edges")):
        return True
    return False


def first_ranking(analysis_pack: dict[str, Any]) -> tuple[str, dict[str, Any], str] | None:
    ranking_specs = [
        ("pathway_rankings", "pathway_uid", "pathway"),
        ("target_rankings", "target_uid", "target"),
        ("disease_rankings", "disease_uid", "disease"),
    ]
    for field, uid_key, label in ranking_specs:
        rows = analysis_pack.get(field) or []
        if rows:
            return uid_key, rows[0], label
    return None


def build_narrative_summary(input_pack: dict[str, Any], catalog: dict[str, Any]) -> dict[str, Any]:
    analysis = source_payloads(input_pack)["analysis_pack"]
    release_ref = release_source_ref(input_pack, catalog)
    sections: list[dict[str, Any]] = []
    summary = analysis.get("input_summary") or {}
    prediction_model = analysis.get("prediction_model") or {}
    assessment = prediction_model.get("assessment") or {}
    analysis_ref = source_ref(catalog, "analysis_pack", "analysis_pack:root", "$.analysis_pack")
    structured_rows = iter_structured_prediction_rows(input_pack)

    def prediction_names(rows: Any, limit: int = 3) -> list[str]:
        names = []
        for row in (rows or [])[:limit]:
            name = str(row.get("display_name") or row.get("result_uid") or row.get("prediction_id") or "").strip()
            if name:
                names.append(name)
        return names

    def structured_rows_by_tier(tiers: set[str], limit: int = 3) -> list[dict[str, Any]]:
        selected = []
        for row in structured_rows:
            tier = str(row.get("confidence_tier") or row.get("confidence") or row.get("display_confidence") or "")
            if tier in tiers and row.get("evidence_refs") and row.get("claim_refs"):
                selected.append(row)
            if len(selected) >= limit:
                break
        return selected

    def structured_names(rows: list[dict[str, Any]]) -> list[str]:
        return [
            str(row.get("display_name") or row.get("drug_name") or row.get("pathway_name") or row.get("_prediction_id") or "")
            for row in rows
            if str(row.get("display_name") or row.get("drug_name") or row.get("pathway_name") or row.get("_prediction_id") or "").strip()
        ]

    def structured_ids(rows: list[dict[str, Any]]) -> list[str]:
        return [str(row.get("_prediction_id") or "") for row in rows if str(row.get("_prediction_id") or "")]

    def structured_refs(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
        return [prediction_ref_from_row(row) for row in rows]

    if summary:
        sections.append(
            {
                "section_id": "precision_input_parsing",
                "text": (
                    "The precision parsing layer reports "
                    f"{summary.get('matched_count', 0)} matched, "
                    f"{summary.get('ambiguous_count', 0)} ambiguous, and "
                    f"{summary.get('unmatched_count', 0)} unmatched input rows; prediction calibration level is "
                    f"{assessment.get('level', 'not_assessed')}."
                ),
                "source_refs": [analysis_ref, release_ref],
            }
        )
    high_structured = structured_rows_by_tier({"high"}, limit=3)
    high_names = structured_names(high_structured) or ([] if structured_rows else prediction_names(prediction_model.get("high_confidence") or [], limit=3))
    if high_names:
        sections.append(
            {
                "section_id": "high_confidence_predictions",
                "reported_display_names": high_names,
                "reported_prediction_ids": structured_ids(high_structured),
                "prediction_refs": structured_refs(high_structured),
                "text": "High-confidence rows are calibrated model predictions with stronger input support and shorter evidence paths.",
                "source_refs": [analysis_ref],
            }
        )
    else:
        sections.append(
            {
                "section_id": "high_confidence_predictions",
                "text": "No high-confidence prediction is present in the current prediction model output.",
                "source_refs": [analysis_ref],
            }
        )
    medium_structured = structured_rows_by_tier({"medium"}, limit=3)
    medium_names = structured_names(medium_structured) or ([] if structured_rows else prediction_names(prediction_model.get("medium_confidence") or [], limit=3))
    sections.append(
        {
            "section_id": "medium_confidence_predictions",
            "reported_display_names": medium_names,
            "reported_prediction_ids": structured_ids(medium_structured),
            "prediction_refs": structured_refs(medium_structured),
            "text": "Medium-confidence predictions can guide candidate mechanism review, but should stay tied to evidence and context.",
            "source_refs": [analysis_ref],
        }
    )
    exploratory_structured = structured_rows_by_tier({"exploratory", "low"}, limit=3)
    exploratory_names = structured_names(exploratory_structured) or ([] if structured_rows else prediction_names(prediction_model.get("exploratory") or [], limit=3))
    sections.append(
        {
            "section_id": "exploratory_predictions",
            "reported_display_names": exploratory_names,
            "reported_prediction_ids": structured_ids(exploratory_structured),
            "prediction_refs": structured_refs(exploratory_structured),
            "text": "Exploratory predictions are retained for discovery, but they require follow-up validation before strong interpretation.",
            "source_refs": [analysis_ref],
        }
    )
    boundary_text = (
        "Review ambiguous or unmatched inputs before trusting calibrated predictions."
        if (summary.get("ambiguous_count") or summary.get("unmatched_count") or summary.get("invalid_count"))
        else "Narration may describe model predictions and confidence tiers, but must not add treatment, diagnosis, or uncited causal claims."
    )
    sections.append(
        {
            "section_id": "calibration_boundaries",
            "text": boundary_text,
            "source_refs": [analysis_ref, release_ref],
        }
    )
    if analysis.get("blocked_reasons"):
        sections.append(
            {
                "section_id": "blocked_reasons",
                "text": "The analysis pack includes blocked reasons, so the adapter should explain the block instead of inventing results.",
                "source_refs": [analysis_ref],
            }
        )
    return {"sections": sections}


def build_evidence_digest(input_pack: dict[str, Any], catalog: dict[str, Any]) -> dict[str, Any]:
    evidence = source_payloads(input_pack)["evidence"]
    items: list[dict[str, Any]] = []
    for index, support in enumerate(evidence.get("support", [])[:3]):
        support_uid = str(support.get("support_uid") or f"support_{index}")
        subject = support.get("subject_uid", "")
        obj = support.get("object_uid", "")
        pmids = [str(item) for item in support.get("pmids", []) if item]
        pmcids = [str(item) for item in support.get("pmcids", []) if item]
        sentence_uids = [str(item) for item in support.get("sentence_uids", []) if item]
        article_bits = []
        if pmids:
            article_bits.append(f"PMID:{pmids[0]}")
        if pmcids:
            article_bits.append(pmcids[0])
        if sentence_uids:
            article_bits.append(f"sentence {sentence_uids[0]}")
        detail = ", ".join(article_bits) if article_bits else "the cited evidence rows"
        refs = [source_ref(catalog, "evidence", support_uid, "$.evidence.support")]
        if sentence_uids:
            refs.append(source_ref(catalog, "evidence", sentence_uids[0], "$.evidence.sentences"))
        items.append(
            {
                "digest_id": f"evidence_{index + 1}",
                "text": f"Evidence support {support_uid} links {subject} to {obj} through {detail}.",
                "source_refs": refs,
            }
        )
    if not items and evidence:
        items.append(
            {
                "digest_id": "evidence_not_found",
                "text": "No sentence-level evidence support is present in the adapter evidence input.",
                "source_refs": [source_ref(catalog, "evidence", "evidence:root", "$.evidence")],
            }
        )
    return {"items": items}


def infer_question_endpoint(question: str) -> str:
    q = question.casefold()
    if any(token in q for token in ("evidence", "pmid", "pmcid", "sentence")):
        return "/evidence"
    if any(token in q for token in ("subgraph", "edge", "neighbor", "path")):
        return "/subgraph"
    if "release" in q or "manifest" in q:
        return "/releases"
    return "/analyze/metabolites"


def build_question_to_spec(input_pack: dict[str, Any], catalog: dict[str, Any]) -> dict[str, Any]:
    question = str((input_pack.get("request") or {}).get("question") or "")
    endpoint = infer_question_endpoint(question)
    return {
        "items": [
            {
                "spec_id": "question_spec_1",
                "endpoint": endpoint,
                "method": "read",
                "params": {},
                "text": f"The question maps to the read-only {endpoint} contract; scoring and entity decisions stay out of adapter scope.",
                "source_refs": [
                    source_ref(catalog, "analysis_pack", "analysis_pack:root", "$.analysis_pack"),
                    release_source_ref(input_pack, catalog),
                ],
            }
        ]
    }


def build_candidate_suggestion(input_pack: dict[str, Any], catalog: dict[str, Any]) -> dict[str, Any]:
    analysis = source_payloads(input_pack)["analysis_pack"]
    items: list[dict[str, Any]] = []
    for bucket in ("ambiguous", "unmatched"):
        for row in (analysis.get(bucket) or [])[:2]:
            input_id = str(row.get("input_id") or f"{bucket}_row")
            items.append(
                {
                    "suggestion_id": f"review_{input_id}",
                    "review_state": "needs_human_review",
                    "text": f"Review input row {input_id}; it remains {bucket} in the analysis pack and the adapter leaves entity resolution unchanged.",
                    "source_refs": [source_ref(catalog, "analysis_pack", input_id, "$.analysis_pack")],
                }
            )
    return {"items": items}


def build_insufficient_evidence(input_pack: dict[str, Any], catalog: dict[str, Any]) -> dict[str, Any]:
    return {
        "reason": "The adapter input contains no analysis rows, evidence support, or subgraph slice to narrate; returning insufficient_evidence.",
        "source_refs": [release_source_ref(input_pack, catalog)],
    }


def build_local_adapter_output(input_pack: dict[str, Any]) -> dict[str, Any]:
    catalog = build_source_catalog(input_pack)
    if not has_explainable_material(input_pack):
        output = {
            "contract_version": OUTPUT_CONTRACT_VERSION,
            "adapter_version": ADAPTER_VERSION,
            "status": "insufficient_evidence",
            "insufficient_evidence": build_insufficient_evidence(input_pack, catalog),
        }
        return finalize_adapter_output(input_pack, output)

    output = {
        "contract_version": OUTPUT_CONTRACT_VERSION,
        "adapter_version": ADAPTER_VERSION,
        "status": "ok",
        "narrative_summary": build_narrative_summary(input_pack, catalog),
        "evidence_digest": build_evidence_digest(input_pack, catalog),
        "question_to_spec": build_question_to_spec(input_pack, catalog),
        "candidate_suggestion": build_candidate_suggestion(input_pack, catalog),
    }
    return finalize_adapter_output(input_pack, output)


def finalize_adapter_output(input_pack: dict[str, Any], output: dict[str, Any]) -> dict[str, Any]:
    bare_output = copy.deepcopy(output)
    guard = guard_adapter_output(input_pack, bare_output)
    final = copy.deepcopy(bare_output)
    final["guard"] = guard
    final["determinism"] = {
        "adapter_input_hash": content_hash(input_pack),
        "adapter_output_hash": content_hash(bare_output),
        "guard_policy_hash": guard_policy_hash(),
    }
    final["determinism"]["response_hash"] = content_hash({key: value for key, value in final.items() if key != "determinism"})
    return final


def llm_narrator_system_prompt() -> str:
    return (
        "You are a narrator over a frozen metabolism-oncology analysis contract. "
        "You are not a resolver, curator, graph writer, or source of truth. "
        "Read only the supplied llm_safe_adapter.input.v1 JSON. "
        "Return only one valid llm_safe_adapter.output.v1 JSON object. "
        "Do not write markdown, code fences, commentary, apologies, or prose outside the JSON object. "
        "All keys must remain in English exactly as requested. "
        "Allowed top-level explanation fields are narrative_summary, evidence_digest, "
        "question_to_spec, candidate_suggestion, and insufficient_evidence. "
        "Every natural-language text, reason, suggestion, digest, narrative, or claim must include source_refs "
        "whose ref_id and path exist in the input. "
        "Use analysis_pack.structured_prediction as the primary prediction JSON when present, and "
        "analysis_pack.prediction_model only for fields already present in the frozen input. "
        "When any output section names or claims a prediction, include reported_prediction_ids and prediction_refs; "
        "each prediction_ref must copy prediction_id, confidence_tier, evidence_refs, claim_refs, "
        "appendix, and research_only from analysis_pack.structured_prediction. "
        "Prediction narration must preserve prediction_task, confidence_tier, calibrated_confidence, "
        "calibration_status, evidence_sources, evidence_refs, claim_refs, appendix/research_only status, and boundary. "
        "The model is allowed to make pathway, target, phenotype, state, and drug-hypothesis predictions, "
        "but the narration must preserve each prediction's confidence tier and validation boundary. "
        "Write bounded narrative sections when possible: precision_input_parsing, high_confidence_predictions, "
        "medium_confidence_predictions, exploratory_predictions, and calibration_boundaries. "
        "Do not resolve entities, choose canonical IDs, add facts, add graph edges, update scores, "
        "change p_final, or turn novel/conflict candidates into curated truth. "
        "Do not say proves, causes, cures, key driver, pathway activated, effective treatment, "
        "or patient treatment suitability. "
        "If support is absent, return status insufficient_evidence. "
        "Treat display names as labels, not as biomedical claims."
    )


def build_external_llm_messages(input_pack: dict[str, Any]) -> list[dict[str, str]]:
    return [
        {"role": "system", "content": llm_narrator_system_prompt()},
        {
            "role": "user",
            "content": (
                "Produce a source-bound explanation JSON for this adapter input. "
                "Return JSON only, with no markdown fences. "
                "Use this shape when evidence is sufficient: "
                '{"contract_version":"llm_safe_adapter.output.v1","adapter_version":"external_llm",'
                '"status":"ok","narrative_summary":{"sections":[{"section_id":"high_confidence_predictions",'
                '"text":"...","reported_prediction_ids":[],"prediction_refs":[],"source_refs":[{"source_type":"analysis_pack","ref_id":"analysis_pack:root",'
                '"path":"$.analysis_pack"}]}]},"question_to_spec":{"items":[{"spec_id":"mapping",'
                '"endpoint":"/analyze/metabolites","method":"read","params":{},'
                '"text":"...","source_refs":[{"source_type":"release","ref_id":"release:'
                f'{input_pack.get("release", {}).get("release_id", "")}","path":"$.release.release_id"}}'
                "]}}]}}. "
                "If evidence is insufficient, return the same contract with status insufficient_evidence and "
                "an insufficient_evidence.reason plus source_refs. "
                "Do not include score, p_final, resolved_entity_uid, new facts, or any field outside the allowed contract. "
                "Do not invent predictions outside analysis_pack.structured_prediction when it exists.\n\n"
                f"{stable_json(input_pack)}"
            ),
        },
    ]


def build_external_llm_request(input_pack: dict[str, Any], config: ExternalLLMConfig) -> dict[str, Any]:
    return {
        "model": config.model,
        "messages": build_external_llm_messages(input_pack),
        "temperature": 0,
        "max_tokens": config.max_output_tokens,
        "response_format": {"type": "json_object"},
    }


def llm_text_narrator_system_prompt() -> str:
    return (
        "你是一个证据校准型代谢机制预测系统的报告层。"
        "你只能根据 frozen analysis/evidence/release JSON 写中文说明。"
        "analysis_pack.prediction_model 是预测依据：必须保留 prediction_task、confidence_tier、"
        "calibrated_confidence、calibration_status、evidence_sources 和 boundary。"
        "允许描述模型预测的代谢主题、通路、靶点、表型状态和候选干预方向，"
        "但必须按高置信、中置信、探索性、低置信分层表达。"
        "不要做实体匹配决定，不要修改分数，不要新增数据库事实，不要给诊断、疗效或治疗建议。"
        "不得使用“证明了、说明药物有效、关键驱动、通路被激活、预测患者适合治疗”等强结论表达。"
        "输出普通中文文本即可，不要输出 JSON，不要 markdown 代码块。"
    )


def compact_text_narrator_input(input_pack: dict[str, Any]) -> dict[str, Any]:
    analysis = input_pack.get("analysis_pack", {}) if isinstance(input_pack, dict) else {}
    prediction = input_pack.get("prediction_pack", {}) if isinstance(input_pack, dict) else {}
    prediction_model = analysis.get("prediction_model", {})

    def compact_resolution_rows(rows: Any, limit: int = 5) -> list[dict[str, Any]]:
        compact_rows: list[dict[str, Any]] = []
        for row in (rows or [])[:limit]:
            candidate = (row.get("candidates") or [{}])[0] if isinstance(row, dict) else {}
            record = row.get("record", {}) if isinstance(row, dict) else {}
            compact_rows.append(
                {
                    "input_id": row.get("input_id"),
                    "input": record.get("name") or record.get("metabolite") or record.get("HMDB") or record.get("ChEBI"),
                    "status": row.get("resolution_status"),
                    "matched_name": candidate.get("display_name") or row.get("metabolite_uid"),
                    "external_id": candidate.get("primary_external_id"),
                    "top_score": row.get("top_score"),
                    "top_margin": row.get("top_margin"),
                }
            )
        return compact_rows

    def compact_ranking_rows(rows: Any, id_key: str, limit: int = 5) -> list[dict[str, Any]]:
        return [
            {
                "rank": row.get("rank"),
                "name": row.get("display_name") or row.get("name"),
                "id": row.get("primary_external_id") or row.get(id_key),
                "score": row.get("score"),
                "prediction_task": row.get("prediction_task"),
                "result_type": row.get("result_type"),
                "confidence_tier": row.get("confidence_tier"),
                "calibrated_confidence": row.get("calibrated_confidence"),
                "calibration_status": row.get("calibration_status"),
                "input_support_count": row.get("input_support_count"),
                "graph_distance": row.get("graph_distance"),
                "boundary": row.get("boundary"),
                "traceable": row.get("claim_refs", {}).get("traceability_passed"),
                "support_classes": row.get("support_classes", [])[:3],
                "pmids": row.get("supported_pmids", [])[:3],
            }
            for row in (rows or [])[:limit]
        ]

    return {
        "question": input_pack.get("request", {}).get("question", ""),
        "release_id": input_pack.get("release", {}).get("release_id", ""),
        "input_summary": analysis.get("input_summary", {}),
        "prediction_assessment": prediction_model.get("assessment", {}),
        "high_confidence": prediction_model.get("high_confidence", [])[:5],
        "medium_confidence": prediction_model.get("medium_confidence", [])[:5],
        "exploratory": prediction_model.get("exploratory", [])[:5],
        "low_confidence": prediction_model.get("low_confidence", [])[:5],
        "quality_warnings": analysis.get("quality_warnings", []),
        "blocked_reasons": analysis.get("blocked_reasons", []),
        "matched": compact_resolution_rows(analysis.get("matched", []), limit=5),
        "ambiguous": compact_resolution_rows(analysis.get("ambiguous", []), limit=5),
        "unmatched": compact_resolution_rows(analysis.get("unmatched", []), limit=5),
        "pathway_rankings": compact_ranking_rows(analysis.get("pathway_rankings", []), "pathway_uid", limit=5),
        "target_rankings": compact_ranking_rows(analysis.get("target_rankings", []), "target_uid", limit=5),
        "disease_rankings": compact_ranking_rows(analysis.get("disease_rankings", []), "disease_uid", limit=5),
        "literature_evidence_pack": analysis.get("literature_evidence_pack", {}),
        "prediction_status": {
            "drug_count": len(prediction.get("drug_rankings", []) or []),
            "cell_type_count": len(prediction.get("cell_type_rankings", []) or []),
            "blocked_reasons": prediction.get("blocked_reasons", [])[:3],
        },
    }


def build_external_text_llm_messages(input_pack: dict[str, Any]) -> list[dict[str, str]]:
    compact_input = compact_text_narrator_input(input_pack)
    return [
        {"role": "system", "content": llm_text_narrator_system_prompt()},
        {
            "role": "user",
            "content": (
                "请根据下面的结构化预测摘要写中文解释。必须按五块组织："
                "1. 精准输入解析；2. 高置信预测；3. 中置信预测；"
                "4. 探索性预测；5. 置信度边界和下一步验证。"
                "每个预测都要保留 confidence_tier 和 boundary 的含义。不要输出 JSON。\n\n"
                f"{stable_json(compact_input)}"
            ),
        },
    ]


def build_external_text_llm_request(input_pack: dict[str, Any], config: ExternalLLMConfig) -> dict[str, Any]:
    return {
        "model": config.model,
        "messages": build_external_text_llm_messages(input_pack),
        "temperature": 0,
        "max_tokens": config.max_output_tokens,
    }


def build_external_text_output(input_pack: dict[str, Any], text: str) -> dict[str, Any]:
    catalog = build_source_catalog(input_pack)
    analysis_ref = source_ref(catalog, "analysis_pack", "analysis_pack:root", "$.analysis_pack")
    release_ref = release_source_ref(input_pack, catalog)
    status = "ok" if has_explainable_material(input_pack) else "insufficient_evidence"
    field = "narrative_summary" if status == "ok" else "insufficient_evidence"
    payload = {
        "reason" if status == "insufficient_evidence" else "sections": [
            {
                "section_id": "external_text_summary",
                "text": text.strip(),
                "source_refs": [analysis_ref, release_ref],
            }
        ]
    }
    if status == "insufficient_evidence":
        payload = {"reason": text.strip(), "source_refs": [analysis_ref, release_ref]}
    return {
        "contract_version": OUTPUT_CONTRACT_VERSION,
        "adapter_version": f"{EXTERNAL_ADAPTER_VERSION}.text",
        "status": status,
        field: payload,
    }


def external_text_section_output(input_pack: dict[str, Any], section_id: str, text: str) -> dict[str, Any]:
    catalog = build_source_catalog(input_pack)
    analysis_ref = source_ref(catalog, "analysis_pack", "analysis_pack:root", "$.analysis_pack")
    release_ref = release_source_ref(input_pack, catalog)
    return {
        "contract_version": OUTPUT_CONTRACT_VERSION,
        "adapter_version": f"{EXTERNAL_ADAPTER_VERSION}.chunked_text",
        "status": "ok",
        "narrative_summary": {
            "sections": [
                {
                    "section_id": section_id,
                    "text": text.strip(),
                    "source_refs": [analysis_ref, release_ref],
                }
            ]
        },
    }


def chunk_fallback_section(input_pack: dict[str, Any], section_id: str, title: str, issue_codes: list[str]) -> dict[str, Any]:
    reason = ", ".join(issue_codes) if issue_codes else "chunk_guard_failed"
    text = (
        f"{title} 分块未通过外部文本安全校验，已保留为结构化结果审阅项；"
        f"原因：{reason}。该分块不用于强化结论。"
    )
    return external_text_section_output(input_pack, section_id, text)["narrative_summary"]["sections"][0]


def compact_feature_rows(features: Any, predicate: Any, limit: int = 12) -> list[dict[str, Any]]:
    rows = []
    for row in features or []:
        if not isinstance(row, dict) or not predicate(row):
            continue
        ratio = row.get("ratio_trait") or {}
        component = row.get("ratio_component") or {}
        class_seed = row.get("class_seed") or {}
        cluster = row.get("identity_cluster") or {}
        rows.append(
            {
                "input_name": row.get("input_name"),
                "display_name": row.get("display_name"),
                "seed_class": row.get("seed_class"),
                "direction": row.get("direction"),
                "ratio_label": ratio.get("label"),
                "ratio_role": component.get("role"),
                "ratio_signed_direction": component.get("signed_ratio_direction"),
                "ratio_boundary": ratio.get("interpretation_boundary") or row.get("ratio_effect_semantics"),
                "class_ids": class_seed.get("class_ids"),
                "chain_constraints": class_seed.get("chain_constraints"),
                "allowed_claim_level": class_seed.get("allowed_claim_level"),
                "class_boundary": class_seed.get("interpretation_boundary"),
                "cluster_basis": cluster.get("cluster_basis"),
                "candidate_count": cluster.get("candidate_count"),
                "candidate_names": cluster.get("candidate_names"),
                "identity_boundary": cluster.get("interpretation_boundary") or row.get("identity_effect_semantics"),
            }
        )
        if len(rows) >= limit:
            break
    return rows


def compact_evidence_rows(evidence: dict[str, Any], limit: int = 6) -> list[dict[str, Any]]:
    rows = []
    for row in (evidence.get("support") or [])[:limit]:
        rows.append(
            {
                "support_uid": row.get("support_uid"),
                "subject_uid": row.get("subject_uid"),
                "object_uid": row.get("object_uid"),
                "pmids": row.get("pmids", [])[:3],
                "pmcids": row.get("pmcids", [])[:3],
                "sentence_uids": row.get("sentence_uids", [])[:3],
                "confidence": row.get("confidence"),
                "evidence_type": row.get("evidence_type"),
            }
        )
    return rows


def build_external_text_chunks(input_pack: dict[str, Any]) -> list[dict[str, Any]]:
    compact = compact_text_narrator_input(input_pack)
    analysis = input_pack.get("analysis_pack", {}) if isinstance(input_pack, dict) else {}
    evidence = input_pack.get("evidence", {}) if isinstance(input_pack, dict) else {}
    features = (analysis.get("analysis_features") or {}).get("matched", [])
    chunks: list[dict[str, Any]] = []

    def add(chunk_id: str, title: str, payload: dict[str, Any], instruction: str) -> None:
        if any(value not in ({}, [], None, "") for value in payload.values()):
            chunks.append({"chunk_id": chunk_id, "title": title, "payload": payload, "instruction": instruction})

    add(
        "precision_input_parsing",
        "精准输入解析",
        {
            "input_summary": compact.get("input_summary"),
            "quality_warnings": compact.get("quality_warnings"),
            "blocked_reasons": compact.get("blocked_reasons"),
            "matched": compact.get("matched"),
            "ambiguous": compact.get("ambiguous"),
            "unmatched": compact.get("unmatched"),
            "differential_table_selection": analysis.get("differential_table_selection"),
            "trait_score_selection": analysis.get("trait_score_selection"),
        },
        "说明输入解析质量、选择阈值、未匹配和低权重信号；不要改写实体身份或显著性。",
    )
    add(
        "ratio_component",
        "ratio component",
        {"ratio_features": compact_feature_rows(features, lambda row: row.get("ratio_component") or row.get("ratio_trait"))},
        "解释 numerator/denominator 或复合比值的相对平衡边界；不要声称组件丰度真实升降。",
    )
    add(
        "class_seed",
        "class seed",
        {"class_features": compact_feature_rows(features, lambda row: row.get("class_seed"))},
        "解释类别、池子、链长或异构体约束；默认降级为类别层或候选脂质物种线索。",
    )
    add(
        "identity_cluster",
        "identity cluster",
        {"identity_features": compact_feature_rows(features, lambda row: row.get("identity_cluster"))},
        "解释软身份候选集合和聚类依据；不要强行唯一化候选。",
    )
    add(
        "pathway_rankings",
        "通路排序",
        {"pathway_rankings": compact.get("pathway_rankings"), "prediction_assessment": compact.get("prediction_assessment")},
        "按置信度和边界解释通路排序；只能作为研究机制线索。",
    )
    add(
        "target_disease_rankings",
        "靶点和疾病关联",
        {"target_rankings": compact.get("target_rankings"), "disease_rankings": compact.get("disease_rankings")},
        "解释靶点/疾病关联的证据边界；不得给临床诊断或治疗建议。",
    )
    add(
        "evidence_digest",
        "证据摘要",
        {"literature_evidence_pack": compact.get("literature_evidence_pack"), "evidence_support": compact_evidence_rows(evidence)},
        "概括证据来源、PMID/句子级支持和缺口；不要新增文献事实。",
    )
    add(
        "low_confidence_appendix",
        "低置信附录",
        {
            "exploratory": compact.get("exploratory"),
            "low_confidence": compact.get("low_confidence"),
            "prediction_status": compact.get("prediction_status"),
        },
        "把探索性和低置信结果作为验证优先级，不得升级为强结论。",
    )
    return chunks


def build_external_chunk_text_llm_request(chunk: dict[str, Any], config: ExternalLLMConfig) -> dict[str, Any]:
    return {
        "model": config.model,
        "messages": [
            {"role": "system", "content": llm_text_narrator_system_prompt()},
            {
                "role": "user",
                "content": (
                    f"请只写“{chunk['title']}”这一块的中文解释，控制在 2-4 段。"
                    f"{chunk['instruction']} 必须保留 confidence_tier、boundary、候选集合或解释边界的含义。"
                    "不要输出 JSON，不要 markdown 代码块，不要新增事实。\n\n"
                    f"{stable_json(chunk['payload'])}"
                ),
            },
        ],
        "temperature": 0,
        "max_tokens": min(max(int(config.max_output_tokens or 1600), 256), 900),
    }


def build_external_final_text_llm_request(sections: list[dict[str, Any]], config: ExternalLLMConfig) -> dict[str, Any]:
    compact_sections = [
        {"section_id": row.get("section_id"), "text": row.get("text", "")[:1800]}
        for row in sections
        if row.get("text")
    ]
    return {
        "model": config.model,
        "messages": [
            {"role": "system", "content": llm_text_narrator_system_prompt()},
            {
                "role": "user",
                "content": (
                    "下面是已经通过安全校验的分块解释。请写一个总述，按："
                    "核心发现、主要不确定性、验证优先级 三段组织。"
                    "只能总结这些分块，不要新增实体、数字、因果或临床建议。不要输出 JSON。\n\n"
                    f"{stable_json(compact_sections)}"
                ),
            },
        ],
        "temperature": 0,
        "max_tokens": min(max(int(config.max_output_tokens or 1600), 256), 900),
    }


def build_external_chunked_text_output(input_pack: dict[str, Any], sections: list[dict[str, Any]]) -> dict[str, Any]:
    if not sections:
        return {
            "contract_version": OUTPUT_CONTRACT_VERSION,
            "adapter_version": f"{EXTERNAL_ADAPTER_VERSION}.chunked_text",
            "status": "insufficient_evidence",
            "insufficient_evidence": build_insufficient_evidence(input_pack, build_source_catalog(input_pack)),
        }
    return {
        "contract_version": OUTPUT_CONTRACT_VERSION,
        "adapter_version": f"{EXTERNAL_ADAPTER_VERSION}.chunked_text",
        "status": "ok",
        "narrative_summary": {"sections": sections},
    }


def default_openai_compatible_transport(config: ExternalLLMConfig, request_payload: dict[str, Any]) -> dict[str, Any]:
    data = json.dumps(request_payload, ensure_ascii=False).encode("utf-8")
    req = urllib_request.Request(
        config.endpoint,
        data=data,
        headers={
            "Authorization": f"Bearer {config.api_key}",
            "Content-Type": "application/json",
        },
        method="POST",
    )
    try:
        if config.proxy_url:
            opener = urllib_request.build_opener(
                urllib_request.ProxyHandler({"http": config.proxy_url, "https": config.proxy_url})
            )
            response_context = opener.open(req, timeout=config.timeout_seconds)  # noqa: S310 - explicit user-enabled endpoint.
        else:
            response_context = urllib_request.urlopen(req, timeout=config.timeout_seconds)  # noqa: S310 - explicit user-enabled endpoint.
        with response_context as response:
            raw = response.read().decode("utf-8")
    except HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")[:1000]
        raise LLMAdapterError("external_llm_http_error", "External LLM request failed.", {"status": exc.code, "body": body}) from exc
    except (TimeoutError, socket.timeout) as exc:
        raise LLMAdapterError("external_llm_transport_error", "External LLM transport failed.", {"reason": str(exc)}) from exc
    except URLError as exc:
        raise LLMAdapterError("external_llm_transport_error", "External LLM transport failed.", {"reason": str(exc.reason)}) from exc
    except OSError as exc:
        raise LLMAdapterError("external_llm_transport_error", "External LLM transport failed.", {"reason": str(exc)}) from exc
    try:
        parsed = json.loads(raw)
    except ValueError as exc:
        raise LLMAdapterError("external_llm_invalid_response", "External LLM response was not JSON.") from exc
    if not isinstance(parsed, dict):
        raise LLMAdapterError("external_llm_invalid_response", "External LLM response JSON must be an object.")
    return parsed


def strip_json_fence(text: str) -> str:
    stripped = text.strip()
    if stripped.startswith("```"):
        stripped = re.sub(r"^```(?:json)?\s*", "", stripped, flags=re.IGNORECASE)
        stripped = re.sub(r"\s*```$", "", stripped)
    return stripped.strip()


def extract_json_object_text(text: str) -> str:
    stripped = strip_json_fence(text)
    if stripped.startswith("{"):
        return stripped
    fence = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", stripped, flags=re.IGNORECASE | re.DOTALL)
    if fence:
        return fence.group(1).strip()
    decoder = json.JSONDecoder()
    for match in re.finditer(r"\{", stripped):
        candidate = stripped[match.start() :]
        try:
            _parsed, end = decoder.raw_decode(candidate)
        except ValueError:
            continue
        return candidate[:end].strip()
    return stripped


def extract_chat_completion_content(response: dict[str, Any]) -> str:
    choices = response.get("choices")
    if isinstance(choices, list) and choices:
        first = choices[0] if isinstance(choices[0], dict) else {}
        message = first.get("message") if isinstance(first.get("message"), dict) else {}
        content = message.get("content")
        if isinstance(content, str) and content.strip():
            return content
        if isinstance(content, list):
            parts = []
            for item in content:
                if isinstance(item, dict) and isinstance(item.get("text"), str):
                    parts.append(item["text"])
            if parts:
                return "".join(parts)
        reasoning_content = message.get("reasoning_content")
        if isinstance(reasoning_content, str) and reasoning_content.strip():
            return reasoning_content
        choice_text = first.get("text")
        if isinstance(choice_text, str) and choice_text.strip():
            return choice_text
        raise LLMAdapterError(
            "external_llm_empty_content",
            "External LLM response contained no assistant content.",
            {
                "choice_keys": sorted(first.keys()),
                "finish_reason": first.get("finish_reason", ""),
                "message_keys": sorted(message.keys()),
            },
        )
    output_text = response.get("output_text")
    if isinstance(output_text, str):
        return output_text
    raise LLMAdapterError("external_llm_missing_content", "External LLM response did not contain message content.")


def parse_external_llm_output(response: dict[str, Any]) -> dict[str, Any]:
    content = extract_json_object_text(extract_chat_completion_content(response))
    try:
        parsed = json.loads(content)
    except ValueError as exc:
        preview = content[:500].replace("\n", " ")
        raise LLMAdapterError(
            "external_llm_output_not_json",
            "External LLM output was not valid JSON.",
            {"output_preview": preview},
        ) from exc
    if not isinstance(parsed, dict):
        raise LLMAdapterError("external_llm_output_not_object", "External LLM output must be a JSON object.")
    return parsed


class LocalNarratorBackend:
    backend_name = "local"

    def generate(self, input_pack: dict[str, Any]) -> BackendResult:
        output = build_local_adapter_output(input_pack)
        return BackendResult(
            output=output,
            audit={
                "backend": self.backend_name,
                "adapter_version": ADAPTER_VERSION,
                "adapter_input_hash": content_hash(input_pack),
                "adapter_output_hash": content_hash(output),
            },
        )


class ExternalLLMNarratorBackend:
    backend_name = "external_llm"

    def __init__(self, config: ExternalLLMConfig | None = None, transport: Any | None = None):
        self.config = config or ExternalLLMConfig.from_env()
        self.transport = transport or default_openai_compatible_transport

    def generate(self, input_pack: dict[str, Any]) -> BackendResult:
        self.config.validate_ready()
        request_payload = build_external_text_llm_request(input_pack, self.config)
        response = self.transport(self.config, request_payload)
        if not isinstance(response, dict):
            raise LLMAdapterError("external_llm_invalid_response", "External LLM transport returned a non-object response.")
        text = extract_chat_completion_content(response).strip()
        if not text:
            raise LLMAdapterError("external_llm_empty_content", "External LLM response contained no assistant content.")
        output = build_external_text_output(input_pack, text)
        output["adapter_version"] = f"{EXTERNAL_ADAPTER_VERSION}.text_fixed_contract"
        audit = {
            "backend": self.backend_name,
            "adapter_version": output.get("adapter_version", f"{EXTERNAL_ADAPTER_VERSION}.text_fixed_contract"),
            "provider": self.config.provider,
            "endpoint": self.config.endpoint,
            "model": self.config.model,
            "prompt_version": f"{LLM_NARRATOR_PROMPT_VERSION}.text_fixed_contract",
            "prompt_hash": content_hash(
                {
                    "prompt_version": f"{LLM_NARRATOR_PROMPT_VERSION}.text_fixed_contract",
                    "messages": request_payload["messages"],
                }
            ),
            "request_hash": content_hash({key: value for key, value in request_payload.items() if key != "messages"}),
            "adapter_input_hash": content_hash(input_pack),
            "raw_response_hash": content_hash(response),
            "raw_output_hash": content_hash(output),
        }
        return BackendResult(output=output, audit=audit)


class ExternalLLMTextNarratorBackend:
    backend_name = "external_text"

    def __init__(self, config: ExternalLLMConfig | None = None, transport: Any | None = None):
        self.config = config or ExternalLLMConfig.from_env()
        self.transport = transport or default_openai_compatible_transport

    def generate(self, input_pack: dict[str, Any]) -> BackendResult:
        self.config.validate_ready()
        request_payload = build_external_text_llm_request(input_pack, self.config)
        response = self.transport(self.config, request_payload)
        if not isinstance(response, dict):
            raise LLMAdapterError("external_llm_invalid_response", "External LLM transport returned a non-object response.")
        text = extract_chat_completion_content(response).strip()
        if not text:
            raise LLMAdapterError("external_llm_empty_content", "External LLM response contained no assistant content.")
        output = build_external_text_output(input_pack, text)
        audit = {
            "backend": self.backend_name,
            "adapter_version": output.get("adapter_version", f"{EXTERNAL_ADAPTER_VERSION}.text"),
            "provider": self.config.provider,
            "endpoint": self.config.endpoint,
            "model": self.config.model,
            "prompt_version": f"{LLM_NARRATOR_PROMPT_VERSION}.text",
            "prompt_hash": content_hash(
                {
                    "prompt_version": f"{LLM_NARRATOR_PROMPT_VERSION}.text",
                    "messages": request_payload["messages"],
                }
            ),
            "request_hash": content_hash({key: value for key, value in request_payload.items() if key != "messages"}),
            "adapter_input_hash": content_hash(input_pack),
            "raw_response_hash": content_hash(response),
            "raw_output_hash": content_hash(output),
        }
        return BackendResult(output=output, audit=audit)


class ExternalLLMChunkedTextNarratorBackend:
    backend_name = "external_text_chunked"

    def __init__(self, config: ExternalLLMConfig | None = None, transport: Any | None = None):
        self.config = config or ExternalLLMConfig.from_env()
        self.transport = transport or default_openai_compatible_transport

    def generate(self, input_pack: dict[str, Any]) -> BackendResult:
        self.config.validate_ready()
        chunks = build_external_text_chunks(input_pack)
        if not chunks:
            output = build_external_chunked_text_output(input_pack, [])
            return BackendResult(
                output=output,
                audit={
                    "backend": self.backend_name,
                    "adapter_version": output.get("adapter_version", f"{EXTERNAL_ADAPTER_VERSION}.chunked_text"),
                    "provider": self.config.provider,
                    "endpoint": self.config.endpoint,
                    "model": self.config.model,
                    "prompt_version": f"{LLM_NARRATOR_PROMPT_VERSION}.chunked_text",
                    "chunk_count": 0,
                    "adapter_input_hash": content_hash(input_pack),
                },
            )

        sections: list[dict[str, Any]] = []
        chunk_audits: list[dict[str, Any]] = []
        prompt_material: list[dict[str, Any]] = []
        raw_response_hashes: list[str] = []
        for chunk in chunks:
            request_payload = build_external_chunk_text_llm_request(chunk, self.config)
            prompt_material.append({"chunk_id": chunk["chunk_id"], "messages": request_payload["messages"]})
            response = self.transport(self.config, request_payload)
            if not isinstance(response, dict):
                raise LLMAdapterError("external_llm_invalid_response", "External LLM transport returned a non-object response.")
            raw_response_hashes.append(content_hash(response))
            text = extract_chat_completion_content(response).strip()
            if not text:
                raise LLMAdapterError("external_llm_empty_content", "External LLM response contained no assistant content.")
            section_id = f"chunk_{chunk['chunk_id']}"
            section_output = external_text_section_output(input_pack, section_id, text)
            guard = guard_adapter_output(input_pack, section_output)
            issue_codes = [str(issue.get("code") or "") for issue in guard.get("issues", []) if issue.get("code")]
            if guard.get("passed"):
                sections.extend(section_output["narrative_summary"]["sections"])
            else:
                sections.append(chunk_fallback_section(input_pack, section_id, str(chunk.get("title") or chunk["chunk_id"]), issue_codes))
            chunk_audits.append(
                {
                    "chunk_id": chunk["chunk_id"],
                    "guard_passed": bool(guard.get("passed")),
                    "guard_issue_codes": issue_codes,
                    "request_hash": content_hash({key: value for key, value in request_payload.items() if key != "messages"}),
                    "raw_response_hash": raw_response_hashes[-1],
                }
            )

        if sections:
            final_request = build_external_final_text_llm_request(sections, self.config)
            prompt_material.append({"chunk_id": "final_summary", "messages": final_request["messages"]})
            response = self.transport(self.config, final_request)
            if not isinstance(response, dict):
                raise LLMAdapterError("external_llm_invalid_response", "External LLM transport returned a non-object response.")
            raw_response_hashes.append(content_hash(response))
            final_text = extract_chat_completion_content(response).strip()
            final_output = external_text_section_output(input_pack, "chunk_final_summary", final_text)
            final_guard = guard_adapter_output(input_pack, final_output)
            final_issue_codes = [str(issue.get("code") or "") for issue in final_guard.get("issues", []) if issue.get("code")]
            if final_guard.get("passed"):
                sections = [*final_output["narrative_summary"]["sections"], *sections]
            chunk_audits.append(
                {
                    "chunk_id": "final_summary",
                    "guard_passed": bool(final_guard.get("passed")),
                    "guard_issue_codes": final_issue_codes,
                    "request_hash": content_hash({key: value for key, value in final_request.items() if key != "messages"}),
                    "raw_response_hash": raw_response_hashes[-1],
                }
            )

        output = build_external_chunked_text_output(input_pack, sections)
        audit = {
            "backend": self.backend_name,
            "adapter_version": output.get("adapter_version", f"{EXTERNAL_ADAPTER_VERSION}.chunked_text"),
            "provider": self.config.provider,
            "endpoint": self.config.endpoint,
            "model": self.config.model,
            "prompt_version": f"{LLM_NARRATOR_PROMPT_VERSION}.chunked_text",
            "prompt_hash": content_hash(
                {
                    "prompt_version": f"{LLM_NARRATOR_PROMPT_VERSION}.chunked_text",
                    "chunks": prompt_material,
                }
            ),
            "chunk_count": len(chunks),
            "passed_chunk_count": sum(1 for row in chunk_audits if row["chunk_id"] != "final_summary" and row["guard_passed"]),
            "chunk_audits": chunk_audits,
            "adapter_input_hash": content_hash(input_pack),
            "raw_response_hashes": raw_response_hashes,
            "raw_output_hash": content_hash(output),
        }
        return BackendResult(output=output, audit=audit)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build or guard a local LLM safe adapter explanation bundle.")
    parser.add_argument("--input-pack", required=True, help="Path to an llm_safe_adapter.input.v1 JSON file.")
    parser.add_argument("--adapter-output", default="", help="Optional adapter output JSON to guard instead of generating local output.")
    parser.add_argument("--output", default="", help="Optional path to write the generated or guarded output JSON.")
    args = parser.parse_args(argv)

    input_pack = read_json(Path(args.input_pack))
    if args.adapter_output:
        adapter_output = read_json(Path(args.adapter_output))
        result = copy.deepcopy(adapter_output)
        result["guard"] = guard_adapter_output(input_pack, adapter_output)
    else:
        result = build_local_adapter_output(input_pack)

    text = json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2)
    if args.output:
        Path(args.output).write_text(text + "\n", encoding="utf-8")
    else:
        sys.stdout.write(text + "\n")
    return 0 if result.get("guard", {}).get("passed") else 2


if __name__ == "__main__":
    raise SystemExit(main())
