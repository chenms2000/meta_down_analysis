"""Build candidate-level training rows from metabolic prediction packs.

The output is a weak-supervision training table for the prediction reranker.
It is a derived research artifact and does not update canonical graph tables.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


DEFAULT_OUTPUT_ROOT = "learning_runs"
DEFAULT_RUN_PREFIX = "prediction_reranker"
RANKING_FIELDS = {
    "pathway_rankings": ("pathway_uid", "pathway_prediction"),
    "target_rankings": ("target_uid", "target_prediction"),
    "disease_rankings": ("disease_uid", "phenotype_context_prediction"),
}
MODEL_ONLY_TASKS = {"metabolic_theme_prediction", "metabolic_state_prediction", "drug_hypothesis_prediction"}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def stable_hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, ensure_ascii=False, sort_keys=True, default=str).encode("utf-8")).hexdigest()


def stable_bucket(value: str, buckets: int = 10) -> int:
    return int(hashlib.sha256(value.encode("utf-8")).hexdigest()[:8], 16) % buckets


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(number) or math.isinf(number):
        return default
    return number


def safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def list_len(value: Any) -> int:
    if isinstance(value, (list, tuple, set)):
        return len(value)
    return 0 if value in (None, "") else 1


def write_parquet(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), path)


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if hasattr(value, "item"):
        return value.item()
    return value


def json_load(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def is_analysis_pack(value: Any) -> bool:
    return isinstance(value, dict) and (
        value.get("contract_version") == "analysis_pack.v1"
        or "prediction_model" in value
        or any(field in value for field in RANKING_FIELDS)
    )


def iter_analysis_packs(value: Any, source_path: Path, context: dict[str, Any] | None = None) -> Iterable[tuple[dict[str, Any], dict[str, Any]]]:
    context = dict(context or {})
    if isinstance(value, dict):
        next_context = dict(context)
        for key in ("case_id", "fixture_version", "release_id"):
            if key in value and key not in next_context:
                next_context[key] = value[key]
        if is_analysis_pack(value):
            yield value, {"source_path": str(source_path), **next_context}
            return
        analysis_pack = value.get("analysis_pack")
        if is_analysis_pack(analysis_pack):
            yield analysis_pack, {"source_path": str(source_path), **next_context}
        for key, child in value.items():
            if key == "analysis_pack":
                continue
            if isinstance(child, (dict, list)):
                yield from iter_analysis_packs(child, source_path, next_context)
    elif isinstance(value, list):
        for child in value:
            yield from iter_analysis_packs(child, source_path, context)


def discover_json_files(paths: list[str]) -> list[Path]:
    files: list[Path] = []
    for raw in paths:
        path = Path(raw).resolve()
        if path.is_dir():
            files.extend(sorted(path.rglob("*.json")))
        elif path.exists():
            files.append(path)
    return files


def prediction_id(row: dict[str, Any], uid_key: str = "") -> str:
    return str(
        row.get("prediction_id")
        or row.get(uid_key)
        or row.get("result_uid")
        or row.get("theme_id")
        or row.get("state_id")
        or row.get("display_name")
        or stable_hash(row)[:16]
    )


def weak_label(row: dict[str, Any]) -> tuple[int, str]:
    tier = str(row.get("confidence_tier") or "").lower()
    status = str(row.get("calibration_status") or "")
    result_type = str(row.get("result_type") or "")
    evidence_sources = set(row.get("evidence_sources") or [])
    input_support = safe_int(row.get("input_support_count"))
    graph_distance = safe_int(row.get("graph_distance"), default=9)
    literature_support = safe_int(row.get("literature_support_count"))
    traceable = bool((row.get("claim_refs") or {}).get("traceability_passed"))
    direct = bool(row.get("is_directly_supported"))
    confidence = safe_float(row.get("calibrated_confidence"))

    high_evidence = input_support > 0 and graph_distance <= 1 and (
        direct or traceable or "input_theme" in evidence_sources or "input_state_signature" in evidence_sources
    )
    if tier == "high" and high_evidence:
        return 2, "high_tier_direct_or_traceable"
    if tier == "medium" and high_evidence and confidence >= 0.45:
        return 2, "medium_tier_strong_support"
    if tier in {"high", "medium"} and input_support > 0:
        return 1, "supported_candidate"
    if tier == "exploratory" and (input_support > 0 or literature_support > 0 or status in {"overlay_informed", "graph_extrapolation"}):
        return 1, "exploratory_supported_candidate"
    if result_type in {"drug", "disease", "cell_context_overlay"} and input_support <= 0:
        return 0, "context_or_overlay_without_input_support"
    if graph_distance >= 3 and input_support <= 0:
        return 0, "distant_without_input_support"
    if tier == "low":
        return 0, "low_confidence_prediction"
    return 0, "weak_or_missing_prediction_signal"


def split_for(analysis_id: str, task: str, candidate_id: str, group_id: str = "") -> str:
    group = group_id or analysis_id or f"{task}:{candidate_id}"
    bucket = stable_bucket(group, buckets=100)
    if bucket < 15:
        return "test"
    if bucket < 30:
        return "validation"
    return "train"


def row_from_prediction(
    pack: dict[str, Any],
    row: dict[str, Any],
    source_context: dict[str, Any],
    uid_key: str = "",
    analysis_id: str = "",
    ranking_field: str = "",
) -> dict[str, Any]:
    input_summary = pack.get("input_summary") or {}
    release = pack.get("release") or {}
    components = row.get("calibration_components") or {}
    evidence_sources = list(row.get("evidence_sources") or [])
    label, label_reason = weak_label(row)
    task = str(row.get("prediction_task") or "")
    candidate_id = prediction_id(row, uid_key)
    group_key = "analysis_id"
    group_id = str(source_context.get("paper_id") or source_context.get("pmid") or "")
    if group_id:
        group_key = "paper_id"
    elif source_context.get("study_id"):
        group_key = "study_id"
        group_id = str(source_context.get("study_id"))
    elif source_context.get("dataset_id"):
        group_key = "dataset_id"
        group_id = str(source_context.get("dataset_id"))
    else:
        group_id = analysis_id
    split = split_for(analysis_id, task, candidate_id, group_id)
    raw_score = safe_float(components.get("raw_score", row.get("score", row.get("ranker_score", row.get("calibrated_confidence")))))
    support_classes = row.get("support_classes") or []
    return {
        "analysis_id": analysis_id,
        "source_path": source_context.get("source_path", ""),
        "case_id": source_context.get("case_id", ""),
        "release_id": release.get("release_id") or pack.get("release_id") or source_context.get("release_id", ""),
        "model_name": (pack.get("prediction_model") or {}).get("model_name", ""),
        "prediction_task": task,
        "ranking_field": ranking_field,
        "candidate_id": candidate_id,
        "candidate_name": row.get("display_name") or row.get("prediction") or candidate_id,
        "candidate_external_id": row.get("primary_external_id", ""),
        "result_type": row.get("result_type", ""),
        "rank": safe_int(row.get("rank")),
        "raw_score": raw_score,
        "ranker_score": safe_float(row.get("ranker_score", row.get("calibrated_confidence"))),
        "current_calibrated_confidence": safe_float(row.get("calibrated_confidence")),
        "current_confidence_tier": row.get("confidence_tier", ""),
        "current_calibration_status": row.get("calibration_status", ""),
        "input_support_count": safe_int(row.get("input_support_count")),
        "matched_input_count": safe_int(row.get("matched_input_count", input_summary.get("matched_count"))),
        "ambiguous_input_count": safe_int(row.get("ambiguous_input_count", input_summary.get("ambiguous_count"))),
        "unmatched_input_count": safe_int(input_summary.get("unmatched_count")),
        "direction_consistency": safe_float(row.get("direction_consistency")),
        "graph_distance": safe_int(row.get("graph_distance"), default=9),
        "node_degree": safe_int(components.get("node_degree")),
        "literature_support_count": safe_int(row.get("literature_support_count")),
        "max_p_literature": safe_float(row.get("max_p_literature")),
        "evidence_source_count": len(set(evidence_sources)),
        "has_input_evidence": int(any(src in evidence_sources for src in ("input", "input_theme", "input_state_signature"))),
        "has_database_evidence": int(any(src in evidence_sources for src in ("database", "database_pathway"))),
        "has_graph_evidence": int("graph_propagation" in evidence_sources),
        "has_literature_evidence": int(any("literature" in src for src in evidence_sources)),
        "has_overlay_evidence": int(any("overlay" in src for src in evidence_sources)),
        "context_match_score": safe_float(row.get("context_match_score", row.get("path_confidence", 0.0))),
        "traceability_passed": int(bool((row.get("claim_refs") or {}).get("traceability_passed"))),
        "edge_ref_count": list_len((row.get("claim_refs") or {}).get("edge_uids")),
        "evidence_ref_count": list_len((row.get("claim_refs") or {}).get("evidence_ref_uids")),
        "support_class_count": list_len(support_classes),
        "is_directly_supported": int(bool(row.get("is_directly_supported"))),
        "is_extrapolated": int(bool(row.get("is_extrapolated"))),
        "context_mismatch": int(bool(row.get("context_mismatch"))),
        "appendix": int(bool(row.get("appendix"))),
        "needs_validation": int(bool(row.get("needs_validation"))),
        "boundary": row.get("boundary", ""),
        "weak_label": label,
        "weak_label_reason": label_reason,
        "group_key": group_key,
        "group_id": group_id,
        "human_label": None,
        "label": label,
        "label_source": "weak_rule_v1",
        "split": split,
    }


def rows_from_pack(pack: dict[str, Any], source_context: dict[str, Any]) -> list[dict[str, Any]]:
    analysis_id = (
        (pack.get("determinism") or {}).get("analysis_pack_hash")
        or stable_hash({"source": source_context, "pack": pack})[:24]
    )
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for ranking_field, (uid_key, default_task) in RANKING_FIELDS.items():
        for row in pack.get(ranking_field) or []:
            if not row.get("prediction_task"):
                row = {**row, "prediction_task": default_task}
            exported = row_from_prediction(pack, row, source_context, uid_key, analysis_id, ranking_field)
            seen.add((exported["prediction_task"], exported["candidate_id"]))
            rows.append(exported)

    prediction_model = pack.get("prediction_model") or {}
    for row in prediction_model.get("predictions") or []:
        task = str(row.get("prediction_task") or "")
        candidate_id = prediction_id(row)
        if (task, candidate_id) in seen and task not in MODEL_ONLY_TASKS:
            continue
        exported = row_from_prediction(pack, row, source_context, "", analysis_id, "prediction_model.predictions")
        seen.add((exported["prediction_task"], exported["candidate_id"]))
        rows.append(exported)
    return rows


def build_service(workspace: Path, release_id: str) -> Any:
    script = workspace / "scripts" / "metabo_service.py"
    spec = importlib.util.spec_from_file_location("metabo_service_for_training_export", script)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import metabo_service.py from {script}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    args = argparse.Namespace(
        workspace=str(workspace),
        normalized_root=module.DEFAULT_NORMALIZED_ROOT,
        graph_root=module.DEFAULT_GRAPH_ROOT,
        pubchem_root=module.DEFAULT_PUBCHEM_ROOT,
        compound_root=module.DEFAULT_COMPOUND_MATCH_ROOT,
        literature_root=module.DEFAULT_LITERATURE_ROOT,
        prediction_root=module.DEFAULT_PREDICTION_OVERLAY_ROOT,
        release_id=release_id,
        enable_external_llm=False,
        llm_provider="",
        llm_endpoint="",
        llm_model="",
        llm_api_key_env="",
        llm_timeout_seconds=None,
        llm_max_output_tokens=None,
        llm_proxy="",
    )
    return module.build_service_from_args(args)


def packs_from_fixture_files(args: argparse.Namespace, workspace: Path) -> Iterable[tuple[dict[str, Any], dict[str, Any]]]:
    fixture_paths = discover_json_files(args.fixture_json)
    if not fixture_paths:
        return
    service = build_service(workspace, args.release_id)
    for path in fixture_paths:
        fixture = json_load(path)
        cases = fixture.get("cases") if isinstance(fixture, dict) else []
        for index, case in enumerate(cases or []):
            if args.max_cases and index >= args.max_cases:
                break
            records = case.get("input")
            if not isinstance(records, list):
                continue
            analyzed = service.analyze_metabolites(
                records,
                max_paths=args.max_paths,
                max_hops=args.max_hops,
                context=case.get("context"),
                context_mode=case.get("context_mode", "soft"),
            )
            pack = analyzed.get("analysis_pack") or {}
            if is_analysis_pack(pack):
                yield pack, {
                    "source_path": str(path.resolve()),
                    "case_id": case.get("case_id", f"case_{index + 1}"),
                    "fixture_version": fixture.get("fixture_version", ""),
                    "release_id": analyzed.get("release", {}).get("release_id", fixture.get("release_id", "")),
                }


def default_fixture_paths(workspace: Path) -> list[str]:
    return [
        str(workspace / "tests" / "fixtures" / "validation" / "analysis_pack_contract_cases.json"),
        str(workspace / "tests" / "fixtures" / "validation" / "real_analysis_fixture_cases.json"),
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export prediction-model training examples from analysis packs.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--run-id", default="")
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--output-dir", default="")
    parser.add_argument("--analysis-json", action="append", default=[], help="Analysis-pack JSON file or directory. Can be repeated.")
    parser.add_argument("--fixture-json", action="append", default=[], help="Fixture JSON with case inputs to analyze before export. Can be repeated.")
    parser.add_argument("--release-id", default="mvp_20260513T002254")
    parser.add_argument("--max-paths", type=int, default=8)
    parser.add_argument("--max-hops", type=int, default=4)
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument("--no-default-fixtures", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    run_id = args.run_id or f"{DEFAULT_RUN_PREFIX}_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
    output_dir = Path(args.output_dir).resolve() if args.output_dir else workspace / args.output_root / run_id / "prediction_training"
    analysis_paths = list(args.analysis_json)
    fixture_paths = list(args.fixture_json)
    if not analysis_paths and not fixture_paths and not args.no_default_fixtures:
        fixture_paths = default_fixture_paths(workspace)
        args.fixture_json = fixture_paths

    rows: list[dict[str, Any]] = []
    pack_count = 0
    for path in discover_json_files(analysis_paths):
        for pack, context in iter_analysis_packs(json_load(path), path):
            pack_rows = rows_from_pack(pack, context)
            rows.extend(pack_rows)
            pack_count += 1
    for pack, context in packs_from_fixture_files(args, workspace) or []:
        pack_rows = rows_from_pack(pack, context)
        rows.extend(pack_rows)
        pack_count += 1

    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame = frame.sort_values(["analysis_id", "prediction_task", "rank", "candidate_id"]).reset_index(drop=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_parquet(output_dir / "training_examples.parquet", frame)

    numeric_columns = [
        column
        for column in frame.columns
        if column
        not in {
            "label",
            "weak_label",
            "human_label",
            "analysis_id",
            "source_path",
            "case_id",
            "release_id",
            "model_name",
            "prediction_task",
            "ranking_field",
            "candidate_id",
            "candidate_name",
            "candidate_external_id",
            "result_type",
            "current_confidence_tier",
            "current_calibration_status",
            "boundary",
            "weak_label_reason",
            "group_key",
            "group_id",
            "label_source",
            "split",
        }
        and pd.api.types.is_numeric_dtype(frame[column])
    ]
    categorical_columns = ["prediction_task", "ranking_field", "result_type", "current_confidence_tier", "current_calibration_status"]
    schema = {
        "schema_version": "prediction_training_examples.v1",
        "label_definition": {
            "2": "should enter high-confidence or top-priority predictions",
            "1": "medium-confidence or candidate mechanism",
            "0": "low-value, noisy, distant, or under-supported prediction",
        },
        "numeric_feature_columns": numeric_columns,
        "categorical_feature_columns": categorical_columns,
        "target_columns": ["label", "weak_label", "human_label"],
    }
    manifest = {
        "created_at_utc": utc_now(),
        "run_id": run_id,
        "analysis_pack_count": pack_count,
        "row_count": int(len(frame)),
        "task_counts": json_safe(frame["prediction_task"].value_counts().to_dict()) if not frame.empty else {},
        "label_counts": json_safe(frame["label"].value_counts().sort_index().to_dict()) if not frame.empty else {},
        "outputs": {
            "training_examples": str(output_dir / "training_examples.parquet"),
            "feature_schema": str(output_dir / "feature_schema.json"),
            "export_manifest": str(output_dir / "export_manifest.json"),
        },
        "notes": [
            "Labels are weak-supervision utility labels, not biological truth labels.",
            "Human-reviewed labels can be added in human_label and used by downstream trainers.",
        ],
    }
    (output_dir / "feature_schema.json").write_text(json.dumps(schema, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output_dir / "export_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"run_id": run_id, "rows": len(frame), "analysis_packs": pack_count, "output_dir": str(output_dir)}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
