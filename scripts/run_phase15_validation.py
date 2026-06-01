"""Phase 1.5 validation gate for the metabolism-oncology MVP.

The runner is intentionally deterministic: it reads fixed fixtures, queries a
frozen release through the read-only service layer, and writes a stable JSON
report without wall-clock timestamps.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

try:
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover - optional for fake-service unit tests.
    pq = None

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from metabo_service import (  # noqa: E402
    DEFAULT_COMPOUND_MATCH_ROOT,
    DEFAULT_GRAPH_ROOT,
    DEFAULT_NORMALIZED_ROOT,
    DEFAULT_PUBCHEM_ROOT,
    MetaboService,
    ServiceConfig,
    content_hash,
    stable_json,
)


REPORT_NAME = "phase15_validation_report.json"
RUNNER_VERSION = "phase1.5.validation.optimization.20260520"
DEFAULT_FIXTURE_DIR = "tests/fixtures/validation"
DEFAULT_OUTPUT_ROOT = "validation_reports"
DEFAULT_LITERATURE_ROOT = "literature_evidence"
RT_MS2_TABLES = {
    "compound_rt.csv": "compound_rt_index.parquet",
    "compound_ms2.csv": "compound_ms2_index.parquet",
}
ANALYSIS_PACK_REQUIRED_FIELDS = {
    "input_summary",
    "matched",
    "ambiguous",
    "unmatched",
    "biological_entity_pools",
    "quality_warnings",
    "pathway_rankings",
    "target_rankings",
    "disease_rankings",
    "prediction_model",
    "top_explanation_paths",
    "literature_evidence_pack",
    "evidence_refs",
    "determinism",
    "release",
    "blocked_reasons",
}
CONFIDENCE_TIERS = {"high", "medium", "exploratory", "low"}
RANKING_CALIBRATION_FIELDS = {
    "prediction_task",
    "confidence_tier",
    "calibrated_confidence",
    "calibration_status",
    "boundary",
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_fixtures(fixture_dir: Path) -> list[dict[str, Any]]:
    fixtures = []
    for path in sorted(fixture_dir.glob("*.json")):
        payload = read_json(path)
        payload["_fixture_file"] = path.name
        fixtures.append(payload)
    return fixtures


def parquet_row_count(path: Path) -> int:
    if not path.exists() or pq is None:
        return 0
    return int(pq.ParquetFile(path).metadata.num_rows)


def literature_evidence_summary(literature_root: Path | None, release_id: str) -> dict[str, Any]:
    if literature_root is None:
        return {
            "status": "not_configured",
            "evidence_candidate_count": 0,
            "normalized_mention_rate": 0.0,
            "supported_existing_edge_count": 0,
            "novel_candidate_count": 0,
            "conflict_candidate_count": 0,
        }
    evidence_dir = literature_root / release_id
    manifest_path = evidence_dir / "literature_evidence_manifest.json"
    if not evidence_dir.exists():
        return {
            "status": "missing",
            "path": str(evidence_dir),
            "evidence_candidate_count": 0,
            "normalized_mention_rate": 0.0,
            "supported_existing_edge_count": 0,
            "novel_candidate_count": 0,
            "conflict_candidate_count": 0,
        }
    manifest = read_json(manifest_path) if manifest_path.exists() else {}
    metrics = manifest.get("metrics", {}) if isinstance(manifest.get("metrics"), dict) else {}
    if not metrics:
        support_counts = {"confirm": 0, "support_direction": 0, "novel_candidate": 0, "conflict_candidate": 0}
        support_path = evidence_dir / "literature_edge_support.parquet"
        if support_path.exists() and pq is not None:
            table = pq.read_table(support_path, columns=["support_class"])
            for row in table.to_pylist():
                support_counts[row.get("support_class", "")] = support_counts.get(row.get("support_class", ""), 0) + 1
        metrics = {
            "evidence_candidate_count": parquet_row_count(evidence_dir / "relation_candidates.parquet"),
            "normalized_mention_rate": 1.0 if parquet_row_count(evidence_dir / "sentence_mentions.parquet") else 0.0,
            "supported_existing_edge_count": support_counts.get("confirm", 0) + support_counts.get("support_direction", 0),
            "novel_candidate_count": support_counts.get("novel_candidate", 0),
            "conflict_candidate_count": support_counts.get("conflict_candidate", 0),
        }
    return {
        "status": "available",
        "path": str(evidence_dir),
        "manifest_hash": manifest.get("manifest_hash", ""),
        "evidence_candidate_count": int(metrics.get("evidence_candidate_count", metrics.get("relation_candidate_count", 0)) or 0),
        "normalized_mention_rate": float(metrics.get("normalized_mention_rate", 0.0) or 0.0),
        "supported_existing_edge_count": int(metrics.get("supported_existing_edge_count", 0) or 0),
        "novel_candidate_count": int(metrics.get("novel_candidate_count", 0) or 0),
        "conflict_candidate_count": int(metrics.get("conflict_candidate_count", 0) or 0),
        "details": {
            "sentence_mention_count": int(metrics.get("sentence_mention_count", 0) or 0),
            "edge_support_count": int(metrics.get("edge_support_count", 0) or 0),
            "unresolved_candidate_count": int(metrics.get("unresolved_candidate_count", 0) or 0),
            "relation_count_by_predicate": metrics.get("relation_count_by_predicate", {}),
            "support_count_by_class": metrics.get("support_count_by_class", {}),
        },
    }


def missing_manual_feature_tables(service: Any, fixture: dict[str, Any]) -> list[str]:
    missing = []
    compound_dir = Path(getattr(service, "compound_dir", ""))
    for table_name in fixture.get("manual_feature_tables", []):
        parquet_name = RT_MS2_TABLES.get(table_name, table_name)
        if parquet_row_count(compound_dir / parquet_name) <= 0:
            missing.append(table_name)
    return missing


def status_matches(actual: str, expected: str) -> bool:
    expected_set = {
        "matched": {"matched"},
        "ambiguous": {"ambiguous"},
        "unmatched": {"unmatched"},
        "invalid": {"invalid"},
        "matched_or_ambiguous": {"matched", "ambiguous"},
        "ambiguous_or_unmatched": {"ambiguous", "unmatched"},
    }.get(str(expected or ""), {str(expected or "")})
    return str(actual or "") in expected_set


def text_contains_label(value: Any, label: str) -> bool:
    needle = str(label or "").casefold()
    if not needle:
        return True
    if isinstance(value, str):
        return needle in value.casefold()
    if isinstance(value, list):
        return any(text_contains_label(item, label) for item in value)
    if isinstance(value, dict):
        return any(text_contains_label(item, label) for item in value.values())
    return False


def endpoint_hashes_equal(first: dict[str, Any], second: dict[str, Any]) -> bool:
    return first.get("determinism", {}).get("response_hash") == second.get("determinism", {}).get("response_hash")


def endpoint_hash_snapshot(response: dict[str, Any]) -> dict[str, Any]:
    determinism = response.get("determinism", {})
    return {
        "request_hash": determinism.get("request_hash", ""),
        "response_hash": determinism.get("response_hash", ""),
    }


def independent_probability_union(values: list[float]) -> float:
    remaining = 1.0
    for value in values:
        probability = min(0.999999, max(0.0, float(value or 0.0)))
        remaining *= 1.0 - probability
    return 1.0 - remaining


P_LITERATURE_SCORE_NOTES = {
    "p_literature(edge)=1-prod(1-p_sentence_i)",
    "p_literature(edge)=1-prod(1-max_sentence_probability_per_article)",
}

P_LITERATURE_COMPONENT_FORMULAS = {
    "1-prod(1-p_sentence_i)",
    "1-prod(1-max_sentence_probability_per_article)",
}


def summarize_candidate(candidate: dict[str, Any] | None) -> dict[str, Any] | None:
    if not candidate:
        return None
    return {
        "entity_uid": candidate.get("entity_uid", ""),
        "entity_type": candidate.get("entity_type", ""),
        "display_name": candidate.get("display_name", ""),
        "primary_external_id": candidate.get("primary_external_id", ""),
        "score": candidate.get("score", 0.0),
        "score_components": candidate.get("score_components", {}),
    }


def summarize_ranking(rows: list[dict[str, Any]], id_key: str, limit: int = 5) -> list[dict[str, Any]]:
    summary = []
    for row in rows[:limit]:
        summary.append(
            {
                id_key: row.get(id_key, ""),
                "name": row.get("name") or row.get("display_name", ""),
                "primary_external_id": row.get("primary_external_id", ""),
                "score": row.get("score", 0.0),
                "score_components": row.get("score_components", {}),
            }
        )
    return summary


def summarize_evidence_support(rows: list[dict[str, Any]], limit: int = 5) -> list[dict[str, Any]]:
    return [
        {
            "support_uid": row.get("support_uid", ""),
            "support_class": row.get("support_class", ""),
            "predicate": row.get("predicate", ""),
            "subject_uid": row.get("subject_uid", ""),
            "object_uid": row.get("object_uid", ""),
            "p_literature": row.get("p_literature", 0.0),
            "pmids": row.get("pmids", [])[:5],
            "pmcids": row.get("pmcids", [])[:5],
            "sentence_uids": row.get("sentence_uids", [])[:5],
        }
        for row in rows[:limit]
    ]


def summarize_relation_candidates(rows: list[dict[str, Any]], limit: int = 5) -> list[dict[str, Any]]:
    return [
        {
            "relation_uid": row.get("relation_uid", ""),
            "predicate": row.get("predicate", ""),
            "subject_uid": row.get("subject_uid", ""),
            "object_uid": row.get("object_uid", ""),
            "sentence_uid": row.get("sentence_uid", ""),
            "pmid": row.get("pmid", ""),
            "pmcid": row.get("pmcid", ""),
            "calibrated_prob": row.get("calibrated_prob", 0.0),
        }
        for row in rows[:limit]
    ]


def query_from_record(record: Any) -> tuple[str, str | None] | None:
    if isinstance(record, str):
        return record, None
    if not isinstance(record, dict):
        return None
    normalized = {str(key).strip().casefold().replace(" ", "_"): value for key, value in record.items()}
    for key, prefix in [
        ("hmdb", "HMDB"),
        ("hmdb_id", "HMDB"),
        ("chebi", "CHEBI"),
        ("chebi_id", "CHEBI"),
        ("pubchem_cid", "PubChem"),
        ("cid", "PubChem"),
        ("pubchem", "PubChem"),
        ("kegg", "KEGG"),
        ("kegg_id", "KEGG"),
    ]:
        value = normalized.get(key)
        if value is None or str(value).strip() == "":
            continue
        text = str(value).strip()
        if ":" not in text:
            text = f"{prefix}:{text}"
        return text, "metabolite"
    for key in ("name", "metabolite", "compound"):
        value = normalized.get(key)
        if value is not None and str(value).strip():
            return str(value).strip(), "metabolite"
    return None


def build_check(name: str, passed: bool, details: dict[str, Any] | None = None) -> dict[str, Any]:
    return {"name": name, "passed": bool(passed), "details": details or {}}


def evidence_request_from_input(value: Any) -> dict[str, Any]:
    allowed = ("edge_uid", "subject_uid", "predicate", "object_uid", "entity_uid", "sentence_uid", "relation_uid", "limit")
    if not isinstance(value, dict):
        return {"edge_uid": str(value or "")}
    request = {}
    for key in allowed:
        if key in value and value[key] not in {None, ""}:
            request[key] = value[key]
    return request


def validate_evidence(service: Any, case: dict[str, Any]) -> dict[str, Any]:
    expected = case.get("expected", {})
    request = evidence_request_from_input(case.get("input", {}))
    first = service.evidence(**request)
    second = service.evidence(**request)
    support_rows = first.get("support", [])
    relation_rows = first.get("relation_candidates", [])
    sentences = first.get("sentences", [])
    mentions = first.get("mentions", [])
    checks = [
        build_check("/evidence_idempotency", endpoint_hashes_equal(first, second), endpoint_hash_snapshot(first)),
    ]
    if expected.get("status"):
        checks.append(
            build_check(
                "evidence_status",
                status_matches(first.get("status", ""), expected["status"]),
                {"actual": first.get("status", ""), "expected": expected["status"]},
            )
        )

    required_support_classes = set(expected.get("support_classes") or [])
    if required_support_classes:
        actual = {row.get("support_class", "") for row in support_rows}
        checks.append(
            build_check(
                "evidence_support_classes",
                required_support_classes.issubset(actual),
                {"actual": sorted(actual), "expected": sorted(required_support_classes)},
            )
        )

    if expected.get("requires_sentence_trace"):
        support_sentence_uids = {uid for row in support_rows for uid in row.get("sentence_uids", [])}
        returned_sentence_uids = {row.get("sentence_uid", "") for row in sentences}
        relation_sentence_uids = {row.get("sentence_uid", "") for row in relation_rows if row.get("sentence_uid")}
        mention_sentence_uids = {row.get("sentence_uid", "") for row in mentions}
        checks.append(
            build_check(
                "evidence_claim_traceability",
                bool(support_rows)
                and bool(support_sentence_uids)
                and support_sentence_uids.issubset(returned_sentence_uids)
                and relation_sentence_uids.issubset(returned_sentence_uids)
                and (not expected.get("requires_mentions") or bool(mention_sentence_uids & returned_sentence_uids)),
                {
                    "support_sentence_uids": sorted(support_sentence_uids)[:10],
                    "returned_sentence_uids": sorted(returned_sentence_uids)[:10],
                    "relation_sentence_uids": sorted(relation_sentence_uids)[:10],
                    "mention_count": len(mentions),
                },
            )
        )

    if expected.get("requires_relation_candidate_provenance"):
        required_fields = ("license_id", "parser_hash", "config_hash")
        checks.append(
            build_check(
                "relation_candidate_provenance",
                bool(relation_rows)
                and all((row.get("pmid") or row.get("pmcid")) and all(row.get(field) for field in required_fields) for row in relation_rows),
                {
                    "relation_candidate_count": len(relation_rows),
                    "required_fields": ["pmid_or_pmcid", *required_fields],
                    "sample": summarize_relation_candidates(relation_rows, limit=3),
                },
            )
        )

    if expected.get("requires_p_literature_contract"):
        formula = first.get("score_notes", {}).get("p_literature", "")
        contract_rows = []
        for row in support_rows:
            components = row.get("score_components", {})
            component_formula = components.get("p_literature_formula", "")
            if component_formula == "1-prod(1-max_sentence_probability_per_article)":
                probabilities = [float(value or 0.0) for value in components.get("article_probabilities", [])]
            else:
                probabilities = [float(value or 0.0) for value in components.get("sentence_probabilities", [])]
            expected_probability = independent_probability_union(probabilities) if probabilities else float(row.get("p_literature") or 0.0)
            actual_probability = float(row.get("p_literature") or 0.0)
            contract_rows.append(
                {
                    "support_uid": row.get("support_uid", ""),
                    "has_formula": component_formula in P_LITERATURE_COMPONENT_FORMULAS,
                    "actual": round(actual_probability, 6),
                    "expected": round(expected_probability, 6),
                    "delta": abs(actual_probability - expected_probability),
                }
            )
        checks.append(
            build_check(
                "p_literature_formula_contract",
                formula in P_LITERATURE_SCORE_NOTES
                and bool(contract_rows)
                and all(row["has_formula"] and row["delta"] <= 1e-6 for row in contract_rows),
                {"score_note": formula, "support_rows": contract_rows[:5]},
            )
        )

    if expected.get("must_not_enter_graph_scoring"):
        disallowed_classes = set(expected.get("support_classes") or ["novel_candidate", "conflict_candidate"])
        target_rows = [row for row in support_rows if row.get("support_class") in disallowed_classes]
        attached_edge_uids = [uid for row in target_rows for uid in row.get("supported_existing_edge_uids", [])]
        target_support_uids = {row.get("support_uid") for row in target_rows}
        scoring_refs = []
        if hasattr(service, "subgraph"):
            seeds = [request.get("subject_uid"), request.get("object_uid")]
            for seed in sorted({str(seed) for seed in seeds if seed}):
                subgraph = service.subgraph([seed], max_hops=1)
                for edge in subgraph.get("edges", []):
                    for support in edge.get("literature_support", []):
                        if support.get("support_class") in disallowed_classes or support.get("support_uid") in target_support_uids:
                            scoring_refs.append({"edge_uid": edge.get("edge_uid", ""), "support_uid": support.get("support_uid", "")})
        checks.append(
            build_check(
                "novel_conflict_not_in_curated_scoring",
                bool(target_rows) and not attached_edge_uids and not scoring_refs,
                {
                    "support_classes": sorted(disallowed_classes),
                    "support_uids": [row.get("support_uid", "") for row in target_rows[:10]],
                    "supported_existing_edge_uids": attached_edge_uids[:10],
                    "subgraph_scoring_refs": scoring_refs[:10],
                },
            )
        )

    return {
        "endpoint": "/evidence",
        "query": request,
        "status": first.get("status", ""),
        "hashes": endpoint_hash_snapshot(first),
        "support_count": len(support_rows),
        "relation_candidate_count": len(relation_rows),
        "sentence_count": len(sentences),
        "mention_count": len(mentions),
        "support_snapshot": summarize_evidence_support(support_rows),
        "relation_candidates_snapshot": summarize_relation_candidates(relation_rows),
        "checks": checks,
    }


def validate_resolver(service: Any, record: Any, expected: dict[str, Any]) -> dict[str, Any] | None:
    query = query_from_record(record)
    if query is None:
        return None
    query_text, entity_type = query
    first = service.resolve(query_text, entity_type=entity_type)
    second = service.resolve(query_text, entity_type=entity_type)
    candidates = first.get("candidates", [])
    top_candidate = candidates[0] if candidates else None
    checks = [
        build_check("resolver_idempotency", endpoint_hashes_equal(first, second), endpoint_hash_snapshot(first)),
    ]
    if expected.get("status") and expected["status"] not in {"ambiguous", "unmatched", "invalid", "ambiguous_or_unmatched"}:
        checks.append(build_check("resolver_status", status_matches(first.get("status", ""), expected["status"]), {"actual": first.get("status", ""), "expected": expected["status"]}))
    if expected.get("entity_type") and top_candidate:
        checks.append(
            build_check(
                "resolver_entity_type",
                top_candidate.get("entity_type") == expected["entity_type"],
                {"actual": top_candidate.get("entity_type"), "expected": expected["entity_type"]},
            )
        )
    if expected.get("label") and top_candidate:
        checks.append(
            build_check(
                "resolver_label",
                text_contains_label(top_candidate, expected["label"]),
                {"expected_label": expected["label"], "top_candidate": summarize_candidate(top_candidate)},
            )
        )
    return {
        "endpoint": "/resolve",
        "query": query_text,
        "status": first.get("status", ""),
        "hashes": endpoint_hash_snapshot(first),
        "top_candidate": summarize_candidate(top_candidate),
        "checks": checks,
    }


def validate_precheck(service: Any, record: Any, expected: dict[str, Any]) -> dict[str, Any]:
    first = service.precheck_metabolites([record])
    second = service.precheck_metabolites([record])
    summary = first.get("summary", {})
    row = None
    actual_status = "unmatched"
    for bucket in ("matched", "ambiguous", "unmatched", "invalid"):
        rows = first.get(bucket, [])
        if rows:
            row = rows[0]
            actual_status = bucket
            break
    candidates = (row or {}).get("resolution", {}).get("candidates", [])
    top_candidate = candidates[0] if candidates else None
    checks = [
        build_check("precheck_idempotency", endpoint_hashes_equal(first, second), endpoint_hash_snapshot(first)),
    ]
    if expected.get("status"):
        checks.append(
            build_check(
                "precheck_status",
                status_matches(actual_status, expected["status"]),
                {"actual": actual_status, "expected": expected["status"], "summary": summary},
            )
        )
    if expected.get("entity_type") and top_candidate:
        checks.append(
            build_check(
                "precheck_entity_type",
                top_candidate.get("entity_type") == expected["entity_type"],
                {"actual": top_candidate.get("entity_type"), "expected": expected["entity_type"]},
            )
        )
    if expected.get("label") and top_candidate:
        checks.append(
            build_check(
                "precheck_label",
                text_contains_label(top_candidate, expected["label"]),
                {"expected_label": expected["label"], "top_candidate": summarize_candidate(top_candidate)},
            )
        )
    required_components = expected.get("required_components") or []
    if required_components and top_candidate:
        components = top_candidate.get("score_components", {})
        checks.append(
            build_check(
                "required_score_components",
                all(float(components.get(component, 0.0) or 0.0) > 0.0 for component in required_components),
                {"required": required_components, "actual": components},
            )
        )
    return {
        "endpoint": "/precheck/metabolites",
        "status": actual_status,
        "summary": summary,
        "hashes": endpoint_hash_snapshot(first),
        "top_candidate": summarize_candidate(top_candidate),
        "checks": checks,
    }


def analysis_pack_hash(pack: dict[str, Any]) -> str:
    return content_hash({key: value for key, value in pack.items() if key != "determinism"})


def analysis_pack_ranking_rows(pack: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for key in ("pathway_rankings", "target_rankings", "disease_rankings"):
        rows.extend(pack.get(key, []) or [])
    return rows


def analysis_pack_traceable_claims(pack: dict[str, Any]) -> tuple[bool, list[dict[str, Any]]]:
    failures = []
    for row in analysis_pack_ranking_rows(pack):
        refs = row.get("claim_refs", {})
        if not refs.get("traceability_passed"):
            failures.append(
                {
                    "rank": row.get("rank"),
                    "id": row.get("pathway_uid") or row.get("target_uid") or row.get("disease_uid", ""),
                    "display_name": row.get("display_name", ""),
                    "claim_refs": refs,
                }
            )
    return not failures, failures


def analysis_pack_ranking_calibration_failures(pack: dict[str, Any]) -> list[dict[str, Any]]:
    failures = []
    for row in analysis_pack_ranking_rows(pack):
        missing = sorted(field for field in RANKING_CALIBRATION_FIELDS if field not in row)
        tier = str(row.get("confidence_tier") or "")
        boundary = str(row.get("boundary") or "").strip()
        calibrated_confidence = row.get("calibrated_confidence")
        confidence_is_numeric = isinstance(calibrated_confidence, int | float)
        confidence_in_range = confidence_is_numeric and 0.0 <= float(calibrated_confidence) <= 1.0
        if missing or tier not in CONFIDENCE_TIERS or not boundary or not confidence_in_range:
            failures.append(
                {
                    "rank": row.get("rank"),
                    "id": row.get("pathway_uid") or row.get("target_uid") or row.get("disease_uid", ""),
                    "display_name": row.get("display_name", ""),
                    "missing_fields": missing,
                    "confidence_tier": tier,
                    "calibrated_confidence": calibrated_confidence,
                    "has_boundary": bool(boundary),
                }
            )
    return failures


def validate_analysis_pack(first: dict[str, Any], second: dict[str, Any], expected: dict[str, Any]) -> list[dict[str, Any]]:
    pack = first.get("analysis_pack", {})
    repeated_pack = second.get("analysis_pack", {})
    missing_fields = sorted(ANALYSIS_PACK_REQUIRED_FIELDS - set(pack))
    checks = [
        build_check(
            "analysis_pack_contract_fields",
            bool(pack) and not missing_fields,
            {"missing_fields": missing_fields, "contract_version": pack.get("contract_version", "")},
        )
    ]
    if not pack:
        return checks

    determinism = pack.get("determinism", {})
    repeated_determinism = repeated_pack.get("determinism", {})
    expected_pack_hash = analysis_pack_hash(pack)
    checks.append(
        build_check(
            "analysis_pack_idempotency",
            determinism.get("analysis_pack_hash") == repeated_determinism.get("analysis_pack_hash") == expected_pack_hash,
            {
                "actual": determinism.get("analysis_pack_hash", ""),
                "repeated": repeated_determinism.get("analysis_pack_hash", ""),
                "expected": expected_pack_hash,
            },
        )
    )
    checks.append(
        build_check(
            "analysis_pack_release_trace",
            pack.get("release", {}).get("release_id") == first.get("release", {}).get("release_id"),
            {
                "pack_release": pack.get("release", {}).get("release_id", ""),
                "response_release": first.get("release", {}).get("release_id", ""),
            },
        )
    )
    traceable, failures = analysis_pack_traceable_claims(pack)
    checks.append(
        build_check(
            "analysis_pack_ranking_claim_traceability",
            bool(analysis_pack_ranking_rows(pack)) and traceable,
            {"failure_count": len(failures), "failures": failures[:5]},
        )
    )
    calibration_failures = analysis_pack_ranking_calibration_failures(pack)
    checks.append(
        build_check(
            "analysis_pack_ranking_calibration_fields",
            bool(analysis_pack_ranking_rows(pack)) and not calibration_failures,
            {"failure_count": len(calibration_failures), "failures": calibration_failures[:5]},
        )
    )
    prediction_model = pack.get("prediction_model", {})
    checks.append(
        build_check(
            "analysis_pack_prediction_model_present",
            isinstance(prediction_model, dict)
            and bool(prediction_model)
            and ("assessment" in prediction_model or "predictions" in prediction_model or "high_confidence" in stable_json(prediction_model)),
            {"keys": sorted(prediction_model)[:20] if isinstance(prediction_model, dict) else []},
        )
    )
    input_summary = pack.get("input_summary", {})
    expected_min_counts = {
        "matched": expected.get("expected_matched_min"),
        "ambiguous": expected.get("expected_ambiguous_min"),
        "unmatched": expected.get("expected_unmatched_min"),
    }
    for bucket, minimum in expected_min_counts.items():
        if minimum is None:
            continue
        actual = len(pack.get(bucket, [])) if bucket != "matched" else int(input_summary.get("matched_count", 0) or 0)
        checks.append(
            build_check(
                f"analysis_pack_{bucket}_min",
                actual >= int(minimum),
                {"actual": actual, "expected_min": int(minimum)},
            )
        )
    expected_direction_counts = expected.get("expected_direction_counts") or {}
    if expected_direction_counts:
        actual_directions = input_summary.get("feature_summary", {}).get("directions", {})
        checks.append(
            build_check(
                "analysis_pack_direction_counts",
                all(int(actual_directions.get(direction, 0) or 0) >= int(minimum) for direction, minimum in expected_direction_counts.items()),
                {"actual": actual_directions, "expected_min": expected_direction_counts},
            )
        )
    if expected.get("expected_with_fold_change_min") is not None:
        actual = int(input_summary.get("feature_summary", {}).get("with_fold_change", 0) or 0)
        checks.append(
            build_check(
                "analysis_pack_fold_change_columns",
                actual >= int(expected["expected_with_fold_change_min"]),
                {"actual": actual, "expected_min": int(expected["expected_with_fold_change_min"])},
            )
        )
    if expected.get("expected_with_p_value_min") is not None:
        actual = int(input_summary.get("feature_summary", {}).get("with_p_value", 0) or 0)
        checks.append(
            build_check(
                "analysis_pack_p_value_columns",
                actual >= int(expected["expected_with_p_value_min"]),
                {"actual": actual, "expected_min": int(expected["expected_with_p_value_min"])},
            )
        )
    required_sections = expected.get("required_ranking_sections") or []
    if expected.get("requires_rankings_non_empty") and not required_sections:
        required_sections = ["pathway_rankings", "target_rankings", "disease_rankings"]
    if required_sections:
        counts = {section: len(pack.get(section, []) or []) for section in required_sections}
        checks.append(
            build_check(
                "analysis_pack_rankings_non_empty",
                all(count > 0 for count in counts.values()),
                {"ranking_counts": counts},
            )
        )
    if expected.get("requires_evidence_refs"):
        checks.append(
            build_check(
                "analysis_pack_evidence_refs_present",
                bool(pack.get("evidence_refs")),
                {"evidence_ref_count": len(pack.get("evidence_refs", []) or [])},
            )
        )
    expected_warning_codes = set(expected.get("expected_quality_warning_codes") or [])
    if expected_warning_codes:
        actual_warning_codes = {row.get("code", "") for row in pack.get("quality_warnings", [])}
        checks.append(
            build_check(
                "analysis_pack_expected_quality_warning_codes",
                expected_warning_codes.issubset(actual_warning_codes),
                {"actual": sorted(actual_warning_codes), "expected": sorted(expected_warning_codes)},
            )
        )
    if expected.get("requires_review_queue"):
        checks.append(
            build_check(
                "analysis_pack_review_queue",
                bool(pack.get("ambiguous") or pack.get("unmatched")),
                {"ambiguous": len(pack.get("ambiguous", [])), "unmatched": len(pack.get("unmatched", []))},
            )
        )
    if expected.get("requires_quality_warnings"):
        checks.append(
            build_check(
                "analysis_pack_quality_warnings",
                bool(pack.get("quality_warnings")),
                {"quality_warnings": pack.get("quality_warnings", [])[:5]},
            )
        )
    if expected.get("requires_duplicate_warning"):
        warning_codes = {row.get("code", "") for row in pack.get("quality_warnings", [])}
        checks.append(
            build_check(
                "analysis_pack_duplicate_warning",
                "duplicate_matched_metabolites" in warning_codes,
                {"warning_codes": sorted(warning_codes)},
            )
        )
    if expected.get("requires_literature_evidence_pack"):
        literature = pack.get("literature_evidence_pack", {})
        checks.append(
            build_check(
                "analysis_pack_literature_evidence_pack",
                int(literature.get("support_count", 0) or 0) > 0 and bool(literature.get("evidence_refs")),
                {
                    "support_count": literature.get("support_count", 0),
                    "support_classes": literature.get("support_classes", []),
                },
            )
        )
    return checks


def validate_analyze(
    service: Any,
    records: list[Any],
    expected: dict[str, Any],
    max_paths: int,
    max_hops: int,
) -> dict[str, Any]:
    first = service.analyze_metabolites(records, max_paths=max_paths, max_hops=max_hops)
    second = service.analyze_metabolites(records, max_paths=max_paths, max_hops=max_hops)
    rankings = first.get("rankings", {})
    pathways = rankings.get("pathways", [])
    targets = rankings.get("targets", [])
    diseases = rankings.get("diseases", [])
    directional = first.get("directional_enrichment", {})
    propagation = first.get("propagation", {})
    explanation_paths = first.get("explanation_paths", [])
    compressed_paths = first.get("compressed_explanation_paths", {})
    checks = [
        build_check("analyze_idempotency", endpoint_hashes_equal(first, second), endpoint_hash_snapshot(first)),
        build_check("propagation_contract", propagation.get("formula") == "pi=(1-alpha)y+alpha*W*pi", {"formula": propagation.get("formula", "")}),
    ]
    checks.extend(validate_analysis_pack(first, second, expected))
    if expected.get("top_pathway_label"):
        label = expected["top_pathway_label"]
        checks.append(
            build_check(
                "pathway_label_in_top_rankings",
                any(text_contains_label(row.get("name", ""), label) for row in pathways[:50]),
                {"expected_label": label, "top_pathway": pathways[0].get("name", "") if pathways else ""},
            )
        )
    if expected.get("directional_enrichment"):
        direction = expected["directional_enrichment"]
        summary = directional.get("summary", {})
        direction_rows = directional.get(direction, [])
        checks.append(
            build_check(
                "directional_enrichment",
                int(summary.get(direction, 0) or 0) > 0 and bool(direction_rows),
                {"direction": direction, "summary": summary},
            )
        )
    if expected.get("requires_explanation_paths"):
        checks.append(build_check("explanation_paths_present", bool(explanation_paths), {"count": len(explanation_paths)}))
        checks.append(
            build_check(
                "compressed_path_signatures_present",
                int(compressed_paths.get("signature_count", 0) or 0) > 0,
                {"signature_count": compressed_paths.get("signature_count", 0)},
            )
        )
    checks.append(
        build_check(
            "propagation_rankings_present",
            bool(pathways) and (bool(targets) or bool(diseases) or propagation.get("graph", {}).get("node_count", 0) > 0),
            {
                "pathways": len(pathways),
                "targets": len(targets),
                "diseases": len(diseases),
                "propagation_graph": propagation.get("graph", {}),
            },
        )
    )
    return {
        "endpoint": "/analyze/metabolites",
        "hashes": endpoint_hash_snapshot(first),
        "precheck_summary": first.get("precheck", {}).get("summary", {}),
        "analysis_feature_summary": first.get("analysis_features", {}).get("summary", {}),
        "analysis_pack": {
            "contract_version": first.get("analysis_pack", {}).get("contract_version", ""),
            "input_summary": first.get("analysis_pack", {}).get("input_summary", {}),
            "quality_warnings": first.get("analysis_pack", {}).get("quality_warnings", [])[:10],
            "blocked_reasons": first.get("analysis_pack", {}).get("blocked_reasons", []),
            "hash": first.get("analysis_pack", {}).get("determinism", {}).get("analysis_pack_hash", ""),
        },
        "rankings_snapshot": {
            "pathways": summarize_ranking(pathways, "pathway_uid"),
            "targets": summarize_ranking(targets, "target_uid"),
            "diseases": summarize_ranking(diseases, "disease_uid"),
        },
        "propagation": {
            "formula": propagation.get("formula", ""),
            "alpha": propagation.get("alpha"),
            "iterations": propagation.get("iterations"),
            "converged": propagation.get("converged"),
            "graph": propagation.get("graph", {}),
        },
        "explanation_path_count": len(explanation_paths),
        "compressed_explanation_paths": {
            "mode": compressed_paths.get("mode", ""),
            "input_path_count": compressed_paths.get("input_path_count", 0),
            "signature_count": compressed_paths.get("signature_count", 0),
            "signatures": (compressed_paths.get("signatures") or [])[:5],
        },
        "checks": checks,
    }


def validate_case(service: Any, fixture: dict[str, Any], case: dict[str, Any], max_paths: int, max_hops: int) -> dict[str, Any]:
    expected = case.get("expected", {})
    endpoint = case.get("endpoint") or fixture.get("endpoint") or case.get("case_type") or fixture.get("case_type")
    if endpoint in {"/evidence", "evidence"}:
        validation = validate_evidence(service, case)
        checks = validation.get("checks", [])
        passed = all(check["passed"] for check in checks)
        return {
            "fixture": fixture.get("_fixture_file", ""),
            "case_id": case.get("case_id", ""),
            "status": "passed" if passed else "failed",
            "validations": [validation],
            "checks": checks,
        }

    missing_features = missing_manual_feature_tables(service, fixture)
    if missing_features:
        return {
            "fixture": fixture.get("_fixture_file", ""),
            "case_id": case.get("case_id", ""),
            "status": "blocked",
            "block_reason": "missing_manual_feature_tables",
            "missing_manual_feature_tables": missing_features,
            "checks": [],
        }

    record_or_records = case.get("input")
    validations = []
    if isinstance(record_or_records, list):
        case_max_paths = int(case.get("max_paths") or fixture.get("max_paths") or max_paths)
        case_max_hops = int(case.get("max_hops") or fixture.get("max_hops") or max_hops)
        validations.append(validate_analyze(service, record_or_records, expected, max_paths=case_max_paths, max_hops=case_max_hops))
    else:
        resolver_validation = validate_resolver(service, record_or_records, expected)
        if resolver_validation:
            validations.append(resolver_validation)
        validations.append(validate_precheck(service, record_or_records, expected))

    checks = [check for validation in validations for check in validation.get("checks", [])]
    passed = all(check["passed"] for check in checks)
    return {
        "fixture": fixture.get("_fixture_file", ""),
        "case_id": case.get("case_id", ""),
        "status": "passed" if passed else "failed",
        "validations": validations,
        "checks": checks,
    }


def report_summary(case_results: list[dict[str, Any]]) -> dict[str, Any]:
    counts = {
        "passed": sum(1 for row in case_results if row.get("status") == "passed"),
        "failed": sum(1 for row in case_results if row.get("status") == "failed"),
        "blocked": sum(1 for row in case_results if row.get("status") == "blocked"),
    }
    counts["total"] = len(case_results)
    return {
        **counts,
        "gate_passed": counts["failed"] == 0,
        "gate_status": "failed" if counts["failed"] else "pass_with_known_blocks" if counts["blocked"] else "passed",
    }


def validation_checks(case_results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [check for case in case_results for check in case.get("checks", [])]


def validations_by_endpoint(case_results: list[dict[str, Any]], endpoint: str) -> list[dict[str, Any]]:
    return [
        validation
        for case in case_results
        for validation in case.get("validations", [])
        if validation.get("endpoint") == endpoint
    ]


def check_passed(checks: list[dict[str, Any]], name: str) -> bool | None:
    matched = [check for check in checks if check.get("name") == name]
    if not matched:
        return None
    return all(check.get("passed") for check in matched)


def optimization_validation_summary(case_results: list[dict[str, Any]]) -> dict[str, Any]:
    checks = validation_checks(case_results)
    prechecks = validations_by_endpoint(case_results, "/precheck/metabolites")
    resolver_checks = [check for check in checks if str(check.get("name", "")).startswith("resolver_")]
    precheck_checks = [check for check in checks if str(check.get("name", "")).startswith("precheck_")]
    precheck_status_counts: dict[str, int] = {}
    for validation in prechecks:
        status = str(validation.get("status") or "unknown")
        precheck_status_counts[status] = precheck_status_counts.get(status, 0) + 1
    calibration_failures = [
        failure
        for check in checks
        if check.get("name") == "analysis_pack_ranking_calibration_fields"
        for failure in check.get("details", {}).get("failures", [])
    ]
    traceability_failures = [
        failure
        for check in checks
        if check.get("name") == "analysis_pack_ranking_claim_traceability"
        for failure in check.get("details", {}).get("failures", [])
    ]
    evidence_precision_checks = [
        "evidence_claim_traceability",
        "relation_candidate_provenance",
        "p_literature_formula_contract",
        "novel_conflict_not_in_curated_scoring",
    ]
    return {
        "entity_resolution": {
            "status": "passed" if resolver_checks + precheck_checks and all(check.get("passed") for check in resolver_checks + precheck_checks) else "failed",
            "resolver_case_count": len(validations_by_endpoint(case_results, "/resolve")),
            "precheck_case_count": len(prechecks),
            "precheck_status_counts": precheck_status_counts,
            "abstention_case_count": precheck_status_counts.get("ambiguous", 0) + precheck_status_counts.get("unmatched", 0),
            "label_checks_passed": check_passed(checks, "precheck_label"),
        },
        "evidence_precision": {
            "status": "passed"
            if all(check_passed(checks, name) is not False for name in evidence_precision_checks)
            else "failed",
            "checks": {name: check_passed(checks, name) for name in evidence_precision_checks},
        },
        "ranking_calibration": {
            "status": "passed"
            if check_passed(checks, "analysis_pack_ranking_calibration_fields") is not False
            and check_passed(checks, "analysis_pack_prediction_model_present") is not False
            else "failed",
            "calibration_failure_count": len(calibration_failures),
            "calibration_failures": calibration_failures[:5],
            "traceability_failure_count": len(traceability_failures),
            "traceability_failures": traceability_failures[:5],
        },
        "release_readiness": {
            "failed_case_count": sum(1 for row in case_results if row.get("status") == "failed"),
            "blocked_case_count": sum(1 for row in case_results if row.get("status") == "blocked"),
            "known_block_reasons": sorted({row.get("block_reason", "") for row in case_results if row.get("status") == "blocked"}),
        },
    }


def evidence_validation_summary(case_results: list[dict[str, Any]]) -> dict[str, Any]:
    validations = [
        validation
        for case in case_results
        for validation in case.get("validations", [])
        if validation.get("endpoint") == "/evidence"
    ]
    if not validations:
        return {
            "status": "not_run",
            "contract_passed": None,
            "evidence_case_count": 0,
            "evidence_idempotency_hashes": [],
            "evidence_claim_traceability": None,
        }
    checks = [check for validation in validations for check in validation.get("checks", [])]
    traceability_checks = [check for check in checks if check.get("name") == "evidence_claim_traceability"]
    contract_passed = all(check.get("passed") for check in checks)
    return {
        "status": "passed" if contract_passed else "failed",
        "contract_passed": contract_passed,
        "evidence_case_count": len(validations),
        "evidence_idempotency_hashes": [validation.get("hashes", {}) for validation in validations],
        "evidence_claim_traceability": all(check.get("passed") for check in traceability_checks) if traceability_checks else None,
        "support_row_count": sum(int(validation.get("support_count", 0) or 0) for validation in validations),
        "relation_candidate_count": sum(int(validation.get("relation_candidate_count", 0) or 0) for validation in validations),
        "sentence_count": sum(int(validation.get("sentence_count", 0) or 0) for validation in validations),
        "mention_count": sum(int(validation.get("mention_count", 0) or 0) for validation in validations),
    }


def build_report(
    service: Any,
    fixtures: list[dict[str, Any]],
    fixture_dir: Path,
    max_paths: int,
    max_hops: int,
    literature_root: Path | None = None,
) -> dict[str, Any]:
    case_results = []
    for fixture in fixtures:
        for case in fixture.get("cases", []):
            case_results.append(validate_case(service, fixture, case, max_paths=max_paths, max_hops=max_hops))
    report = {
        "runner_version": RUNNER_VERSION,
        "release_id": getattr(service, "release_id", ""),
        "fixture_dir": str(fixture_dir.as_posix()),
        "fixture_files": [fixture.get("_fixture_file", "") for fixture in fixtures],
        "service_release": service.release_meta() if hasattr(service, "release_meta") else {"release_id": getattr(service, "release_id", "")},
        "config": service.config.as_dict() if hasattr(getattr(service, "config", None), "as_dict") else {},
        "literature_evidence": literature_evidence_summary(literature_root, getattr(service, "release_id", "")),
        "evidence_validation": evidence_validation_summary(case_results),
        "optimization_summary": optimization_validation_summary(case_results),
        "summary": report_summary(case_results),
        "cases": case_results,
    }
    report["report_hash"] = content_hash(report)
    return report


def write_report(report: dict[str, Any], output_root: Path, release_id: str) -> Path:
    output_dir = output_root / release_id
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / REPORT_NAME
    path.write_text(json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2), encoding="utf-8")
    return path


def build_service(args: argparse.Namespace) -> MetaboService:
    workspace = Path(args.workspace).resolve()
    config = ServiceConfig(
        max_paths=args.max_paths,
        max_hops=args.max_hops,
        max_propagation_nodes=args.max_propagation_nodes,
        max_propagation_edges=args.max_propagation_edges,
        propagation_beam_per_type=args.propagation_beam_per_type,
        propagation_retained_mass=args.propagation_retained_mass,
        propagation_degree_penalty=not args.disable_propagation_degree_penalty,
    )
    return MetaboService(
        workspace=workspace,
        normalized_root=(workspace / args.normalized_root).resolve(),
        graph_root=(workspace / args.graph_root).resolve(),
        pubchem_root=(workspace / args.pubchem_root).resolve(),
        compound_root=(workspace / args.compound_root).resolve(),
        literature_root=(workspace / args.literature_root).resolve(),
        release_id=args.release_id or None,
        config=config,
    )


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 1.5 validation fixtures against a frozen release.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", default="")
    parser.add_argument("--fixture-dir", default=DEFAULT_FIXTURE_DIR)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--literature-root", default=DEFAULT_LITERATURE_ROOT)
    parser.add_argument("--normalized-root", default=DEFAULT_NORMALIZED_ROOT)
    parser.add_argument("--graph-root", default=DEFAULT_GRAPH_ROOT)
    parser.add_argument("--pubchem-root", default=DEFAULT_PUBCHEM_ROOT)
    parser.add_argument("--compound-root", default=DEFAULT_COMPOUND_MATCH_ROOT)
    parser.add_argument("--max-paths", type=int, default=25)
    parser.add_argument("--max-hops", type=int, default=4)
    parser.add_argument("--max-propagation-nodes", type=int, default=1000)
    parser.add_argument("--max-propagation-edges", type=int, default=5000)
    parser.add_argument("--propagation-beam-per-type", type=int, default=40)
    parser.add_argument("--propagation-retained-mass", type=float, default=0.95)
    parser.add_argument("--disable-propagation-degree-penalty", action="store_true")
    parser.add_argument("--no-write-report", action="store_true")
    parser.add_argument("--fail-on-blocked", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    workspace = Path(args.workspace).resolve()
    fixture_dir = (workspace / args.fixture_dir).resolve()
    output_root = (workspace / args.output_root).resolve()
    literature_root = (workspace / args.literature_root).resolve()
    service = build_service(args)
    fixtures = load_fixtures(fixture_dir)
    report = build_report(
        service,
        fixtures,
        fixture_dir.relative_to(workspace) if fixture_dir.is_relative_to(workspace) else fixture_dir,
        args.max_paths,
        args.max_hops,
        literature_root=literature_root,
    )
    report_path = None
    if not args.no_write_report:
        report_path = write_report(report, output_root, service.release_id)
    print(
        json.dumps(
            {
                "release_id": service.release_id,
                "gate_status": report["summary"]["gate_status"],
                "summary": report["summary"],
                "report_hash": report["report_hash"],
                "report_path": str(report_path) if report_path else "",
            },
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
        )
    )
    if report["summary"]["failed"] > 0:
        return 1
    if args.fail_on_blocked and report["summary"]["blocked"] > 0:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
