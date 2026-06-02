"""Run release-level benchmark checks across input modes and gold conclusions.

This runner is intentionally read-only for release artifacts. It exercises the
service with representative input scenarios, records conclusion-evaluation
metrics, checks for missing real-world fixtures such as the cSCC CSV, and can
optionally rerun temporal/source held-out validation.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from metabo_service import MetaboService, latest_release_id  # noqa: E402
from audit_real_world_fixtures import expected_schema, find_candidate_paths, recovery_rerun_plan  # noqa: E402
from audit_gold_standard_coverage import audit_coverage  # noqa: E402
from audit_gold_standard_quality import audit_quality  # noqa: E402
from run_literature_semantics_benchmark import benchmark_literature_semantics  # noqa: E402


RUNNER_VERSION = "release_batch_benchmark.v5"
DEFAULT_OUTPUT_ROOT = "validation_reports"
DEFAULT_CSSC_CSV = "required_metabolism_fixture.csv"
DEFAULT_RUN_ID = "full_context_engine_20260519"
DEFAULT_GOLD_PATH = "config/gold_standard_conclusions.json"
DEFAULT_GOLD_COVERAGE_REQUIREMENTS = "config/gold_standard_coverage_requirements.json"
DEFAULT_GOLD_QUALITY_REQUIREMENTS = "config/gold_standard_quality_requirements.json"
DEFAULT_LITERATURE_SEMANTICS_BENCHMARK = "config/literature_semantics_benchmark.json"
DEFAULT_REAL_WORLD_FIXTURE_REGISTRY = "config/real_world_fixture_registry.json"
DEFAULT_MAX_HELDOUT_AGE_HOURS = 168
DEFAULT_HELDOUT_RERUN_TIMEOUT_SECONDS = 900
DEFAULT_SCENARIO_TIMEOUT_SECONDS = 300
DEFAULT_MIN_TEMPORAL_HELDOUT_ROWS = 100
DEFAULT_MIN_SOURCE_HELDOUT_SOURCES = 2
DEFAULT_MIN_SOURCE_HELDOUT_ROWS = 100
DEFAULT_MIN_HELDOUT_LIFT_AT_100 = 1.0
DEFAULT_MIN_HELDOUT_ROC_AUC = 0.5
DEFAULT_MIN_HELDOUT_UNIQUE_TRANSFER_SCORES = 2
DEFAULT_MAX_HELDOUT_TOP_SCORE_TIE_FRACTION = 0.99
DEFAULT_REQUIRED_REAL_FIXTURE_CASES = ["cscc_trait_score_csv"]
DEFAULT_REQUIRED_BENCHMARK_CANCERS = ["pan_cancer", "brca", "coad_read", "luad", "hcc", "melanoma", "cscc"]
DEFAULT_REQUIRED_BENCHMARK_THEMES = [
    "glycolysis_lactate",
    "tca_cycle",
    "glutamine",
    "serine_one_carbon",
    "fatty_acid_lipid",
    "bile_acid",
    "redox_glutathione",
    "amino_acid",
    "arachidonate",
]
DEFAULT_REQUIRED_BENCHMARK_INPUT_MODES = [
    "stable_identifier",
    "compound_name",
    "trait_score",
    "lcms_feature",
    "lipid_class",
    "ratio_trait",
]
TRACEABILITY_CONTRACT_VERSION = "scenario_traceability_contract.v1"
OPTIMIZATION_SIGNAL_TARGETS = {
    "execution_error": ["benchmark_runner"],
    "expected_topk_recall_gap": ["ranking_calibration", "mechanism_fact_priority"],
    "false_negative_proxy": ["benchmark_runner", "entity_coverage", "relation_parser", "gold_standard"],
    "false_positive_trap": ["negative_trap_block_rule", "context_model", "graph_propagation", "ranking_calibration"],
    "fixture_missing": ["real_world_fixture_registry", "database_schema"],
    "fixture_schema_invalid": ["database_schema", "input_parser"],
    "gold_context_mismatch": ["context_model", "gold_standard"],
    "identity_ambiguity_review": ["resolver", "input_parser"],
    "identity_unmatched_review": ["resolver", "synonym_index"],
    "evidence_not_effective": ["evidence_assertion_linker", "literature_weight"],
    "negative_trap_hit": ["negative_trap_block_rule", "ranking_calibration"],
    "unsupported_claim_risk": ["evidence_scoring", "output_contract"],
    "context_drift_review": ["context_model", "graph_propagation"],
    "trait_score_confidence_overclaim": ["ranking_calibration", "output_contract", "scenario_traceability_contract"],
    "benchmark_matrix_gap": ["benchmark_runner", "real_world_fixture_registry"],
    "gold_coverage_gap": ["gold_standard", "coverage_requirements"],
    "gold_quality_gap": ["gold_standard", "negative_traps"],
    "heldout_validation_gap": ["heldout_validation", "ranker_calibration"],
    "literature_semantics_regression": ["literature_parser", "support_status_classifier"],
    "required_fixture_missing": ["real_world_fixture_registry", "database_schema"],
    "scenario_child_failed": ["benchmark_runner", "analysis_runtime"],
    "scenario_timeout": ["benchmark_runner", "analysis_runtime", "batching_cache"],
    "scenario_traceability_contract_failed": ["output_contract", "analysis_pack_contract"],
}
DEFAULT_CSSC_BENCHMARK_TAGS = {
    "cancer_ids": ["cscc"],
    "theme_ids": ["amino_acid"],
    "input_modes": ["trait_score"],
}
CSSC_IDENTITY_COLUMNS = [
    "trait",
    "reported_trait",
    "trait_name",
    "name",
    "metabolite",
    "compound",
    "hmdb",
    "chebi",
    "mz",
]
CSSC_EFFECT_COLUMNS = [
    "cohen_d",
    "log2fc",
    "log2FC",
    "effect_size",
    "trait_score_diff",
    "score_diff",
    "delta",
    "diff",
    "statistic",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def parse_utc_datetime(value: Any) -> datetime | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def slugify(value: str, fallback: str = "scenario") -> str:
    slug = "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in str(value or "").strip())
    slug = "_".join(part for part in slug.split("_") if part)
    return slug[:80] or fallback


def progress_event(
    events: list[dict[str, Any]],
    event: str,
    *,
    case_id: str = "",
    progress_log: Path | None = None,
    emit: bool = False,
    **fields: Any,
) -> dict[str, Any]:
    row = {"created_at_utc": utc_now(), "event": event}
    if case_id:
        row["case_id"] = case_id
    row.update(fields)
    events.append(row)
    if progress_log is not None:
        progress_log.parent.mkdir(parents=True, exist_ok=True)
        with progress_log.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
    if emit:
        print(json.dumps(row, ensure_ascii=False, sort_keys=True), flush=True)
    return row


def apply_release_profile(args: argparse.Namespace) -> argparse.Namespace:
    if bool(getattr(args, "strict_release", False)):
        args.require_real_fixtures = True
        args.require_heldout_rerun = True
        args.rerun_heldout = True
        args.release_profile = "strict"
    else:
        args.release_profile = "exploratory"
    return args


def release_profile_summary(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "profile": getattr(args, "release_profile", "exploratory"),
        "strict_release": bool(getattr(args, "strict_release", False)),
        "require_real_fixtures": bool(getattr(args, "require_real_fixtures", False)),
        "rerun_heldout": bool(getattr(args, "rerun_heldout", False)),
        "require_heldout_rerun": bool(getattr(args, "require_heldout_rerun", False)),
        "heldout_rerun_timeout_seconds": int(
            getattr(args, "heldout_rerun_timeout_seconds", DEFAULT_HELDOUT_RERUN_TIMEOUT_SECONDS)
        ),
        "scenario_timeout_seconds": int(getattr(args, "scenario_timeout_seconds", DEFAULT_SCENARIO_TIMEOUT_SECONDS)),
        "max_heldout_age_hours": int(getattr(args, "max_heldout_age_hours", DEFAULT_MAX_HELDOUT_AGE_HOURS)),
        "required_fixture_cases": list(getattr(args, "required_fixture_case", []) or []),
    }


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value) if str(value).strip() else default
    except (TypeError, ValueError):
        return default


def sum_numeric_values(payload: dict[str, Any]) -> int:
    total = 0
    for value in payload.values():
        try:
            total += int(value)
        except (TypeError, ValueError):
            continue
    return total


def parse_optional_float(value: Any) -> float | None:
    try:
        text = str(value).strip()
        return float(text) if text else None
    except (TypeError, ValueError):
        return None


def column_lookup(columns: list[str]) -> dict[str, str]:
    return {str(column).strip().casefold(): str(column) for column in columns if str(column).strip()}


def first_present_column(columns: list[str], candidates: list[str]) -> str:
    lookup = column_lookup(columns)
    for candidate in candidates:
        column = lookup.get(candidate.casefold())
        if column:
            return column
    return ""


def infer_input_mode(columns: list[str]) -> str:
    lookup = column_lookup(columns)
    if "trait" in lookup or "reported_trait" in lookup or "trait_name" in lookup:
        return "trait_score"
    if "hmdb" in lookup or "chebi" in lookup:
        return "stable_identifier"
    if "mz" in lookup:
        return "lcms_feature"
    if any(column in lookup for column in ("name", "metabolite", "compound")):
        return "compound_name"
    return "unknown"


def validate_fixture_csv_rows(
    rows: list[dict[str, Any]],
    columns: list[str],
    identity_columns: list[str],
    effect_columns: list[str],
    expected_input_mode: str = "",
    min_rows: int = 1,
) -> dict[str, Any]:
    identity_column = first_present_column(columns, identity_columns)
    effect_column = first_present_column(columns, effect_columns)
    numeric_effect_rows = 0
    if effect_column:
        numeric_effect_rows = sum(1 for row in rows if parse_optional_float(row.get(effect_column)) is not None)
    input_mode = infer_input_mode(columns)
    errors: list[dict[str, Any]] = []
    if len(rows) < min_rows:
        errors.append({"code": "row_count_below_minimum", "message": f"Observed {len(rows)} rows, required {min_rows}."})
    if not identity_column:
        errors.append(
            {
                "code": "missing_identity_column",
                "message": "CSV must include a trait, compound, stable identifier, or mz column for input mode detection.",
            }
        )
    if not effect_column:
        errors.append(
            {
                "code": "missing_effect_column",
                "message": "CSV must include an effect/ranking column such as cohen_d, log2FC, trait_score_diff, delta, or diff.",
            }
        )
    elif numeric_effect_rows <= 0:
        errors.append({"code": "effect_column_not_numeric", "message": f"Effect column {effect_column} has no numeric values."})
    if expected_input_mode and input_mode != expected_input_mode:
        errors.append({"code": "input_mode_mismatch", "message": f"Expected {expected_input_mode}, observed {input_mode}."})
    return {
        "status": "passed" if not errors else "failed",
        "columns": columns,
        "row_count": len(rows),
        "identity_column": identity_column,
        "effect_column": effect_column,
        "numeric_effect_rows": numeric_effect_rows,
        "input_mode": input_mode,
        "expected_input_mode": expected_input_mode,
        "errors": errors,
    }


def validate_cssc_csv_rows(rows: list[dict[str, Any]], columns: list[str]) -> dict[str, Any]:
    return validate_fixture_csv_rows(rows, columns, CSSC_IDENTITY_COLUMNS, CSSC_EFFECT_COLUMNS, "trait_score", 1)


def fixture_benchmark_tags(case_id: str, fixture: dict[str, Any]) -> dict[str, Any]:
    tags = fixture.get("benchmark_tags") or {}
    if isinstance(tags, dict) and tags:
        return tags
    if case_id == "cscc_trait_score_csv":
        return dict(DEFAULT_CSSC_BENCHMARK_TAGS)
    return {}


def scenario_payload(
    case_id: str,
    records: list[dict[str, Any]],
    context: dict[str, Any] | None = None,
    expected_gold_ids: list[str] | None = None,
    expected_identity_behavior: dict[str, Any] | None = None,
    benchmark_tags: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "case_id": case_id,
        "records": records,
        "context": context or {},
        "expected_gold_ids": expected_gold_ids or [],
        "expected_identity_behavior": expected_identity_behavior or {},
        "benchmark_tags": benchmark_tags or {},
        "status": "ready",
    }


def default_scenarios() -> list[dict[str, Any]]:
    return [
        scenario_payload(
            "direct_hmdb_glucose_epithelial_tumor",
            [{"HMDB": "HMDB0000122", "log2FC": 1.0, "padj": 0.01, "direction": "up"}],
            {"cell_type": "epithelial", "context_terms": ["tumor", "epithelial"]},
            ["pan_cancer_glycolysis_lactate_axis"],
            benchmark_tags={
                "cancer_ids": ["pan_cancer"],
                "theme_ids": ["glycolysis_lactate"],
                "input_modes": ["stable_identifier"],
            },
        ),
        scenario_payload(
            "lipid_class_review_only",
            [{"name": "Sphingomyelin", "log2FC": 1.0, "padj": 0.02, "direction": "up"}],
            {"context_terms": ["tumor", "lipid metabolism"]},
            [],
            {"allow_ambiguity": True, "reason": "lipid class or pool must not be promoted to exact compound identity"},
            {
                "cancer_ids": ["pan_cancer"],
                "theme_ids": ["fatty_acid_lipid"],
                "input_modes": ["lipid_class", "compound_name"],
            },
        ),
        scenario_payload(
            "ratio_trait_review_only",
            [{"trait": "GCST90201000", "reported_trait": "Tryptophan to Pyruvate ratio", "log2FC": 1.0, "padj": 0.01}],
            {"context_terms": ["tumor", "amino acid metabolism"]},
            [],
            {"allow_ambiguity": True, "reason": "ratio trait must remain observation-only and not become abundance seed"},
            {
                "cancer_ids": ["pan_cancer"],
                "theme_ids": ["amino_acid"],
                "input_modes": ["ratio_trait", "trait_score"],
            },
        ),
        scenario_payload(
            "mass_only_ambiguous_review",
            [{"mz": 180.063, "adduct": "[M+H]+", "ppm_tolerance": 10, "log2FC": 1.0}],
            {"context_terms": ["tumor"]},
            [],
            {"allow_ambiguity": True, "reason": "mass-only input should abstain unless orthogonal evidence is present"},
            {
                "cancer_ids": ["pan_cancer"],
                "theme_ids": ["lcms_ambiguity"],
                "input_modes": ["lcms_feature"],
            },
        ),
        scenario_payload(
            "pan_cancer_arachidonate_eicosanoid",
            [{"name": "Arachidonic acid", "log2FC": 0.8, "padj": 0.03, "direction": "up"}],
            {"context_terms": ["tumor", "arachidonate", "eicosanoid"]},
            ["pan_cancer_arachidonate_eicosanoid_axis"],
            benchmark_tags={
                "cancer_ids": ["pan_cancer"],
                "theme_ids": ["arachidonate"],
                "input_modes": ["compound_name"],
            },
        ),
        scenario_payload(
            "brca_glucose_lactate_glutamine",
            [
                {"name": "Glucose", "log2FC": 1.0, "padj": 0.01, "direction": "up"},
                {"name": "Lactate", "log2FC": 1.1, "padj": 0.02, "direction": "up"},
                {"name": "Glutamine", "log2FC": -0.8, "padj": 0.03, "direction": "down"},
            ],
            {"cancer_type": "breast cancer", "tissue": "breast", "context_terms": ["BRCA", "breast cancer", "tumor"]},
            ["brca_glycolysis_glutamine_lactate_axis", "pan_cancer_glycolysis_lactate_axis", "pan_cancer_glutamine_glutamate_axis"],
            benchmark_tags={
                "cancer_ids": ["brca"],
                "theme_ids": ["glycolysis_lactate", "glutamine"],
                "input_modes": ["compound_name"],
            },
        ),
        scenario_payload(
            "coad_read_serine_glutathione_lipid",
            [
                {"name": "Serine", "log2FC": 0.9, "padj": 0.02, "direction": "up"},
                {"name": "Glycine", "log2FC": 0.6, "padj": 0.04, "direction": "up"},
                {"name": "Glutathione", "log2FC": -0.7, "padj": 0.04, "direction": "down"},
                {"name": "Fatty acid", "log2FC": 0.5, "padj": 0.05, "direction": "up"},
            ],
            {"cancer_type": "colorectal cancer", "tissue": "colon", "context_terms": ["COAD", "READ", "CRC", "colorectal cancer"]},
            [
                "coad_read_colorectal_core_metabolism_axis",
                "coad_read_serine_one_carbon_axis",
                "coad_read_redox_glutathione_axis",
                "coad_read_fatty_acid_lipid_axis",
            ],
            benchmark_tags={
                "cancer_ids": ["coad_read"],
                "theme_ids": ["serine_one_carbon", "redox_glutathione", "fatty_acid_lipid"],
                "input_modes": ["compound_name", "lipid_class"],
            },
        ),
        scenario_payload(
            "luad_lactate_glutamine",
            [
                {"name": "Lactate", "log2FC": 1.0, "padj": 0.01, "direction": "up"},
                {"name": "Glutamine", "log2FC": -0.7, "padj": 0.04, "direction": "down"},
                {"name": "Glutamate", "log2FC": 0.8, "padj": 0.02, "direction": "up"},
            ],
            {"cancer_type": "lung adenocarcinoma", "tissue": "lung", "context_terms": ["LUAD", "lung adenocarcinoma", "tumor"]},
            ["luad_lung_glucose_glutamine_axis", "pan_cancer_glycolysis_lactate_axis", "pan_cancer_glutamine_glutamate_axis"],
            benchmark_tags={
                "cancer_ids": ["luad"],
                "theme_ids": ["glycolysis_lactate", "glutamine"],
                "input_modes": ["compound_name"],
            },
        ),
        scenario_payload(
            "hcc_bile_acid_fatty_acid_glutamine",
            [
                {"name": "Bile acid", "log2FC": 0.9, "padj": 0.03, "direction": "up"},
                {"name": "Fatty acid", "log2FC": 0.7, "padj": 0.02, "direction": "up"},
                {"name": "Glutamine", "log2FC": -0.6, "padj": 0.04, "direction": "down"},
            ],
            {"cancer_type": "hepatocellular carcinoma", "tissue": "liver", "context_terms": ["HCC", "hepatocellular carcinoma", "liver"]},
            ["hcc_bile_acid_cholesterol_axis", "hcc_bile_acid_lipid_glutamine_axis"],
            {"allow_ambiguity": True, "reason": "bile acid and fatty acid are class-level inputs, not exact compound identities"},
            {
                "cancer_ids": ["hcc"],
                "theme_ids": ["bile_acid", "fatty_acid_lipid", "glutamine"],
                "input_modes": ["compound_name", "lipid_class"],
            },
        ),
        scenario_payload(
            "melanoma_glutamine_tca_lipid",
            [
                {"name": "Glutamine", "log2FC": -0.8, "padj": 0.03, "direction": "down"},
                {"name": "Citrate", "log2FC": 0.6, "padj": 0.04, "direction": "up"},
                {"name": "Fatty acid", "log2FC": 0.5, "padj": 0.05, "direction": "up"},
            ],
            {"cancer_type": "melanoma", "tissue": "skin", "context_terms": ["melanoma", "skin", "tumor"]},
            [
                "melanoma_glutamine_glutamate_axis",
                "melanoma_glutamine_tca_fatty_acid_axis",
                "melanoma_fatty_acid_lipid_axis",
            ],
            benchmark_tags={
                "cancer_ids": ["melanoma"],
                "theme_ids": ["glutamine", "tca_cycle", "fatty_acid_lipid"],
                "input_modes": ["compound_name", "lipid_class"],
            },
        ),
    ]


def load_cssc_scenario(path: Path, top_n: int) -> dict[str, Any]:
    if not path.exists():
        return {
            "case_id": "cscc_trait_score_csv",
            "status": "missing",
            "path": str(path),
            "message": "cSCC CSV was not found in the current workspace; this scenario is reported but not executed.",
            "benchmark_tags": dict(DEFAULT_CSSC_BENCHMARK_TAGS),
        }
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        columns = list(reader.fieldnames or [])
    schema_audit = validate_cssc_csv_rows(rows, columns)
    if schema_audit["status"] != "passed":
        return {
            "case_id": "cscc_trait_score_csv",
            "status": "failed",
            "path": str(path),
            "error_type": "invalid_fixture_schema",
            "error_message": "cSCC CSV failed schema/input-mode validation.",
            "fixture_schema": schema_audit,
        }
    effect_column = str(schema_audit.get("effect_column") or "")
    selected = sorted(rows, key=lambda row: abs(safe_float(row.get(effect_column))), reverse=True)[:top_n]
    return {
        "case_id": "cscc_trait_score_csv",
        "status": "ready",
        "path": str(path),
        "records": selected,
        "input_rows": len(rows),
        "selected_rows": len(selected),
        "fixture_schema": schema_audit,
        "context": {
            "cancer_type": "cutaneous squamous cell carcinoma",
            "tissue": "skin",
            "cell_type": "epithelial",
            "comparison": "Tumor vs Adjacent",
            "context_terms": ["cSCC", "cutaneous squamous cell carcinoma", "epithelial", "tumor", "adjacent"],
        },
        "expected_gold_ids": ["cscc_amino_acid_metabolism_axis"],
        "benchmark_tags": dict(DEFAULT_CSSC_BENCHMARK_TAGS),
    }


def registered_fixture_scenario(workspace: Path, fixture: dict[str, Any]) -> dict[str, Any]:
    case_id = str(fixture.get("case_id") or "")
    path = workspace / str(fixture.get("path") or "")
    if not path.exists():
        return {
            "case_id": case_id,
            "status": "missing",
            "path": str(path),
            "message": "Registered real-world fixture was not found in the current workspace; this scenario is reported but not executed.",
            "expected_schema": expected_schema(fixture),
            "candidate_paths": find_candidate_paths(workspace, fixture),
            "recovery_action": "Restore the registered fixture path or update the registry to one of the candidate paths after schema validation.",
            "recovery_rerun_plan": recovery_rerun_plan(fixture),
            "benchmark_tags": fixture_benchmark_tags(case_id, fixture),
            "fixture_registry": {
                "fixture_type": fixture.get("fixture_type", ""),
                "required_for_strict_release": bool(fixture.get("required_for_strict_release")),
            },
        }
    fixture_type = str(fixture.get("fixture_type") or "")
    if not fixture_type.startswith("csv"):
        return {
            "case_id": case_id,
            "status": "failed",
            "path": str(path),
            "error_type": "unsupported_fixture_type",
            "error_message": f"Unsupported registered fixture type: {fixture_type}",
        }
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        columns = list(reader.fieldnames or [])
    schema_audit = validate_fixture_csv_rows(
        rows,
        columns,
        [str(column) for column in fixture.get("identity_columns_any", []) or []],
        [str(column) for column in fixture.get("effect_columns_any", []) or []],
        str(fixture.get("expected_input_mode") or ""),
        int(fixture.get("min_rows") or 1),
    )
    if schema_audit["status"] != "passed":
        return {
            "case_id": case_id,
            "status": "failed",
            "path": str(path),
            "error_type": "invalid_fixture_schema",
            "error_message": "Registered real-world CSV failed schema/input-mode validation.",
            "fixture_schema": schema_audit,
        }
    effect_column = str(schema_audit.get("effect_column") or "")
    top_n = int(fixture.get("top_n") or 80)
    selected = sorted(rows, key=lambda row: abs(safe_float(row.get(effect_column))), reverse=True)[:top_n]
    return {
        "case_id": case_id,
        "status": "ready",
        "path": str(path),
        "records": selected,
        "input_rows": len(rows),
        "selected_rows": len(selected),
        "fixture_schema": schema_audit,
        "fixture_registry": {
            "fixture_type": fixture_type,
            "required_for_strict_release": bool(fixture.get("required_for_strict_release")),
        },
        "context": fixture.get("context") or {},
        "expected_gold_ids": [str(gold_id) for gold_id in fixture.get("expected_gold_ids", []) or []],
        "expected_identity_behavior": fixture.get("expected_identity_behavior") or {},
        "benchmark_tags": fixture_benchmark_tags(case_id, fixture),
    }


def registry_required_cases(registry: dict[str, Any]) -> list[str]:
    return [
        str(fixture.get("case_id") or "")
        for fixture in registry.get("fixtures", []) or []
        if isinstance(fixture, dict) and fixture.get("required_for_strict_release") and str(fixture.get("case_id") or "")
    ]


def load_registry_scenarios(workspace: Path, registry: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        registered_fixture_scenario(workspace, fixture)
        for fixture in registry.get("fixtures", []) or []
        if isinstance(fixture, dict) and str(fixture.get("case_id") or "")
    ]


def dedupe_preserve(values: list[str]) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        if value and value not in seen:
            seen.add(value)
            result.append(value)
    return result


def scenario_expected_topk_metrics(evaluation: dict[str, Any], expected_gold_ids: list[str]) -> dict[str, Any]:
    expected = [gold_id for gold_id in expected_gold_ids if gold_id]
    metrics: dict[str, Any] = {
        "expected_topk_eligible": bool(expected),
        "expected_gold_count": len(expected),
    }
    detail: dict[str, Any] = {
        "eligible": bool(expected),
        "expected_gold_ids": expected,
        "by_k": {},
    }
    per_gold = {str(row.get("gold_id") or ""): row for row in evaluation.get("per_gold", []) or [] if isinstance(row, dict)}
    ranked_gold = [
        str(row.get("gold_id") or "")
        for row in evaluation.get("top_ranked_gold_conclusions", []) or []
        if isinstance(row, dict) and str(row.get("gold_id") or "")
    ]
    ranked_candidates = [
        str(row.get("candidate_id") or "")
        for row in evaluation.get("top_ranked_candidates", []) or []
        if isinstance(row, dict) and str(row.get("candidate_id") or "")
    ]
    for k in (1, 3, 5):
        gold_window = ranked_gold[:k]
        candidate_window = ranked_candidates[:k]
        candidate_hits: set[str] = set()
        gold_hits: set[str] = set()
        if ranked_gold:
            gold_hits = {gold_id for gold_id in expected if gold_id in gold_window}
            candidate_hits = set(gold_hits)
            denominator = min(k, len(gold_window))
        else:
            for gold_id in expected:
                for match in (per_gold.get(gold_id) or {}).get("top_matches", []) or []:
                    candidate_id = str(match.get("candidate_id") or "")
                    if candidate_id and candidate_id in candidate_window:
                        candidate_hits.add(candidate_id)
                        gold_hits.add(gold_id)
            denominator = min(k, len(candidate_window))
        metrics[f"expected_precision_at_{k}"] = round(len(candidate_hits) / denominator, 6) if denominator else None
        metrics[f"expected_recall_at_{k}"] = round(len(gold_hits) / len(expected), 6) if expected else None
        detail["by_k"][str(k)] = {
            "candidate_hit_count": len(candidate_hits),
            "matched_gold_ids": sorted(gold_hits),
            "missing_gold_ids": [gold_id for gold_id in expected if gold_id not in gold_hits],
            "ranked_gold_ids": gold_window,
            "ranked_candidate_ids": candidate_window,
            "ranking_basis": "top_ranked_gold_conclusions" if ranked_gold else "top_ranked_candidates",
        }
    return {"metrics": metrics, "detail": detail}


def scenario_prediction_quality(pack: dict[str, Any], input_summary: dict[str, Any]) -> dict[str, Any]:
    analysis_mode = str(input_summary.get("analysis_mode") or "")
    weak_trait_score_high: list[dict[str, Any]] = []
    context_mismatch_predictions: list[dict[str, Any]] = []
    high_or_medium_context_mismatch = 0
    for ranking_name in ("pathway_rankings", "target_rankings", "disease_rankings"):
        for row in pack.get(ranking_name, []) or []:
            confidence_tier = str(row.get("confidence_tier") or "")
            display_name = str(row.get("display_name") or "")
            if row.get("context_mismatch"):
                if confidence_tier in {"high", "medium"}:
                    high_or_medium_context_mismatch += 1
                if len(context_mismatch_predictions) < 20:
                    context_mismatch_predictions.append(
                        {
                            "ranking": ranking_name,
                            "rank": row.get("rank"),
                            "display_name": display_name,
                            "confidence_tier": confidence_tier,
                            "calibration_status": row.get("calibration_status", ""),
                            "downgrade_reason": row.get("downgrade_reason", ""),
                        }
                    )
            if (
                ranking_name == "pathway_rankings"
                and analysis_mode == "two_group_trait_comparison_table"
                and confidence_tier == "high"
            ):
                overlap_count = int(row.get("overlap_count") or 0)
                literature_support_count = int(row.get("literature_support_count") or 0)
                significant_support_count = int(row.get("significant_support_count") or 0)
                if overlap_count < 2 or significant_support_count < 2 or literature_support_count <= 0:
                    weak_trait_score_high.append(
                        {
                            "rank": row.get("rank"),
                            "display_name": display_name,
                            "overlap_count": overlap_count,
                            "literature_support_count": literature_support_count,
                            "significant_support_count": significant_support_count,
                            "downgrade_reason": row.get("downgrade_reason", ""),
                        }
                    )
    return {
        "analysis_mode": analysis_mode,
        "high_confidence_trait_score_weak_support_count": len(weak_trait_score_high),
        "high_confidence_trait_score_weak_support": weak_trait_score_high[:20],
        "context_mismatch_prediction_count": len(context_mismatch_predictions),
        "high_or_medium_context_mismatch_prediction_count": high_or_medium_context_mismatch,
        "context_mismatch_predictions": context_mismatch_predictions,
    }


def summarize_analysis(case: dict[str, Any], analyzed: dict[str, Any]) -> dict[str, Any]:
    pack = analyzed.get("analysis_pack") or {}
    evaluation = pack.get("conclusion_evaluation") or {}
    input_summary = pack.get("input_summary") or {}
    database_accuracy = pack.get("database_accuracy") or {}
    identity_summary = database_accuracy.get("identity_decision_summary") or {}
    evidence_summary = database_accuracy.get("evidence_assertion_summary") or {}
    mechanism_facts = database_accuracy.get("mechanism_ready_facts") or []
    determinism = pack.get("determinism") or {}
    warnings = pack.get("quality_warnings") or []
    per_gold = evaluation.get("per_gold", []) or []
    per_gold_status_counts: dict[str, int] = {}
    for row in per_gold:
        status = str(row.get("gold_match_status") or "unknown")
        per_gold_status_counts[status] = per_gold_status_counts.get(status, 0) + 1
    per_gold_by_id = {str(row.get("gold_id") or ""): row for row in per_gold}
    expected_gold_ids = [str(gold_id) for gold_id in case.get("expected_gold_ids", []) or [] if str(gold_id)]
    expected_statuses = {
        gold_id: str((per_gold_by_id.get(gold_id) or {}).get("gold_match_status") or "missing_from_gold_standard")
        for gold_id in expected_gold_ids
    }
    expected_match_statuses = {"exact_match", "partial_match"}
    expected_matched = sum(1 for status in expected_statuses.values() if status in expected_match_statuses)
    evaluation_metrics = dict(evaluation.get("metrics", {}) or {})
    expected_topk = scenario_expected_topk_metrics(evaluation, expected_gold_ids)
    evaluation_metrics.update(expected_topk["metrics"])
    scenario_expected_metrics = {
        "expected_gold_ids": expected_gold_ids,
        "expected_count": len(expected_gold_ids),
        "expected_matched": expected_matched,
        "expected_recall": round(expected_matched / len(expected_gold_ids), 6) if expected_gold_ids else None,
        "expected_statuses": expected_statuses,
        "missing_expected_gold_ids": [
            gold_id for gold_id, status in expected_statuses.items() if status not in expected_match_statuses
        ],
    }
    return {
        "case_id": case.get("case_id", ""),
        "status": "passed" if (evaluation.get("release_gate") or {}).get("passed", True) else "failed",
        "benchmark_tags": case.get("benchmark_tags") or {},
        "input_summary": {
            key: input_summary.get(key)
            for key in (
                "analysis_mode",
                "input_count",
                "matched_count",
                "ambiguous_count",
                "unmatched_count",
                "invalid_count",
                "expanded_candidate_count",
                "genetic_exposure_count",
            )
        },
        "quality_warning_codes": [warning.get("code", "") for warning in warnings],
        "prediction_quality": scenario_prediction_quality(pack, input_summary),
        "expected_identity_behavior": case.get("expected_identity_behavior") or {},
        "scenario_expected_metrics": scenario_expected_metrics,
        "traceability": {
            "contract_version": TRACEABILITY_CONTRACT_VERSION,
            "input_summary_present": bool(input_summary),
            "identity_decision_summary": identity_summary,
            "identity_decision_count": sum_numeric_values(identity_summary),
            "mechanism_ready_fact_count": len(mechanism_facts) if isinstance(mechanism_facts, list) else 0,
            "evidence_assertion_summary": evidence_summary,
            "evidence_assertion_count": int(evidence_summary.get("assertion_count") or 0),
            "evidence_support_statuses": sorted((evidence_summary.get("by_support_status") or {}).keys()),
            "conclusion_evaluation_present": bool(evaluation),
            "evaluation_metrics_present": bool(evaluation.get("metrics")),
            "release_gate_present": bool(evaluation.get("release_gate")),
            "analysis_pack_hash_present": bool(determinism.get("analysis_pack_hash")),
            "response_hash_present": bool((analyzed.get("determinism") or {}).get("response_hash")),
        },
        "conclusion_evaluation": {
            "gold_standard": evaluation.get("gold_standard", {}),
            "metrics": evaluation_metrics,
            "expected_topk": expected_topk["detail"],
            "top_ranked_gold_conclusions": evaluation.get("top_ranked_gold_conclusions", [])[:10],
            "release_gate": evaluation.get("release_gate", {}),
            "bug_queue": evaluation.get("bug_queue", []),
            "per_gold_status_counts": dict(sorted(per_gold_status_counts.items())),
            "matched_gold": [
                {
                    "gold_id": row.get("gold_id", ""),
                    "status": row.get("gold_match_status", ""),
                }
                for row in evaluation.get("per_gold", []) or []
                if row.get("gold_match_status") not in {
                    "unsupported",
                    "not_applicable_context",
                    "appendix_only",
                    "negative_control_not_triggered",
                }
            ][:20],
        },
        "determinism": {
            "analysis_pack_hash": determinism.get("analysis_pack_hash", ""),
            "response_hash": (analyzed.get("determinism") or {}).get("response_hash", ""),
        },
    }


def audit_scenario_traceability(scenarios: list[dict[str, Any]]) -> dict[str, Any]:
    failed_cases: list[dict[str, Any]] = []
    checked_count = 0
    skipped_count = 0
    for case in scenarios:
        if case.get("status") == "missing":
            skipped_count += 1
            continue
        checked_count += 1
        trace = case.get("traceability") or {}
        summary = case.get("input_summary") or {}
        failures: list[dict[str, Any]] = []
        input_count = int(summary.get("input_count") or 0)
        if not trace.get("input_summary_present"):
            failures.append({"code": "input_summary_missing", "message": "Scenario did not expose an input summary."})
        if input_count <= 0:
            failures.append({"code": "input_summary_empty", "message": "Scenario input summary has no input rows."})
        if int(trace.get("identity_decision_count") or 0) <= 0 and input_count > 0:
            failures.append({"code": "identity_decisions_missing", "message": "Scenario has inputs but no identity-decision summary."})
        if "mechanism_ready_fact_count" not in trace:
            failures.append({"code": "mechanism_facts_unreported", "message": "Scenario did not report mechanism-ready fact count."})
        if int(trace.get("evidence_assertion_count") or 0) <= 0:
            failures.append({"code": "evidence_assertions_missing", "message": "Scenario did not expose structured evidence assertions."})
        if not trace.get("conclusion_evaluation_present"):
            failures.append({"code": "conclusion_evaluation_missing", "message": "Scenario did not expose conclusion evaluation."})
        if not trace.get("evaluation_metrics_present"):
            failures.append({"code": "evaluation_metrics_missing", "message": "Scenario did not expose conclusion evaluation metrics."})
        if not trace.get("release_gate_present"):
            failures.append({"code": "release_gate_missing", "message": "Scenario did not expose conclusion release-gate status."})
        if failures:
            failed_cases.append({"case_id": case.get("case_id", ""), "failures": failures})
    failure_codes = sorted({failure["code"] for case in failed_cases for failure in case.get("failures", [])})
    return {
        "contract_version": TRACEABILITY_CONTRACT_VERSION,
        "status": "passed" if not failed_cases else "failed",
        "checked_count": checked_count,
        "skipped_missing_count": skipped_count,
        "failed_case_count": len(failed_cases),
        "failure_codes": failure_codes,
        "failed_cases": failed_cases,
        "interpretation": (
            "Every executed scenario must expose auditable input summaries, identity decisions, mechanism facts, "
            "structured evidence assertions, conclusion metrics, and release-gate state."
        ),
    }


def attach_optimization_targets(signals: list[dict[str, Any]]) -> list[dict[str, Any]]:
    for row in signals:
        row["optimization_targets"] = sorted(set(OPTIMIZATION_SIGNAL_TARGETS.get(str(row.get("code") or ""), ["manual_review"])))
    return signals


def case_optimization_signals(case: dict[str, Any]) -> list[dict[str, Any]]:
    if case.get("status") == "missing":
        return attach_optimization_targets([
            {
                "code": "fixture_missing",
                "severity": "warning",
                "recommended_action": "Restore the real-world fixture or register an explicit unavailable-fixture waiver.",
            }
        ])
    if case.get("status") == "failed" and case.get("error_type") == "invalid_fixture_schema":
        return attach_optimization_targets([
            {
                "code": "fixture_schema_invalid",
                "severity": "error",
                "recommended_action": "Inspect CSV columns and input mode; require identity and numeric effect/ranking fields before analysis.",
            }
        ])
    if case.get("status") == "failed" and not case.get("conclusion_evaluation"):
        return attach_optimization_targets([
            {
                "code": "execution_error",
                "severity": "error",
                "recommended_action": "Inspect the exception path before interpreting benchmark metrics.",
            }
        ])

    signals: list[dict[str, Any]] = []
    summary = case.get("input_summary") or {}
    metrics = ((case.get("conclusion_evaluation") or {}).get("metrics") or {})
    status_counts = ((case.get("conclusion_evaluation") or {}).get("per_gold_status_counts") or {})
    input_count = int(summary.get("input_count") or 0)
    ambiguous_count = int(summary.get("ambiguous_count") or 0)
    unmatched_count = int(summary.get("unmatched_count") or 0)
    expected_identity = case.get("expected_identity_behavior") or {}
    expected_ambiguity = bool(expected_identity.get("allow_ambiguity"))
    if input_count and ambiguous_count / input_count >= 0.5 and not expected_ambiguity:
        signals.append(
            {
                "code": "identity_ambiguity_review",
                "severity": "warning",
                "recommended_action": "Review resolver schema, exact identifiers, lipid/class handling, and ambiguity abstention rules.",
            }
        )
    if input_count and unmatched_count / input_count >= 0.5:
        signals.append(
            {
                "code": "identity_unmatched_review",
                "severity": "warning",
                "recommended_action": "Inspect synonym/index coverage and input parser mode detection.",
            }
        )
    if int(metrics.get("negative_trap_hits") or 0) > 0:
        signals.append(
            {
                "code": "false_positive_trap",
                "severity": "error",
                "recommended_action": "Block release and inspect schema, resolver, parser, and graph propagation weights.",
            }
        )
    expected_metrics = case.get("scenario_expected_metrics") or {}
    expected_total = expected_metrics.get("expected_count")
    expected_matched = expected_metrics.get("expected_matched")
    if expected_total is not None and int(expected_total or 0) > 0 and int(expected_matched or 0) < int(expected_total or 0):
        signals.append(
            {
                "code": "false_negative_proxy",
                "severity": "info",
                "recommended_action": "Classify missing scenario-expected gold hits by entity coverage, context mismatch, relation parser gap, or gold-set overbreadth.",
            }
        )
    expected_recall_at_5 = metrics.get("expected_recall_at_5")
    if metrics.get("expected_topk_eligible") and expected_recall_at_5 is not None and float(expected_recall_at_5 or 0.0) < 1.0:
        signals.append(
            {
                "code": "expected_topk_recall_gap",
                "severity": "warning",
                "recommended_action": "Inspect expected-top-k missing gold IDs, candidate source priority, alias coverage, and mechanism-fact ranking.",
            }
        )
    unsupported_rate = metrics.get("unsupported_top_claim_rate")
    if unsupported_rate is not None and float(unsupported_rate or 0.0) > 0:
        signals.append(
            {
                "code": "unsupported_claim_risk",
                "severity": "warning",
                "recommended_action": "Ensure top claims have evidence refs or mechanism traces; keep observation-only rows appendix-only.",
            }
        )
    context_accuracy = metrics.get("context_accuracy")
    if context_accuracy is not None and float(context_accuracy or 0.0) < 1.0:
        signals.append(
            {
                "code": "context_drift_review",
                "severity": "warning",
                "recommended_action": "Inspect disease/tissue/cell-type terms and downgrade out-of-context graph propagation.",
            }
        )
    prediction_quality = case.get("prediction_quality") or {}
    if int(prediction_quality.get("context_mismatch_prediction_count") or 0) > 0:
        signals.append(
            {
                "code": "context_drift_review",
                "severity": "warning",
                "recommended_action": "Review context-mismatch prediction rows and keep incompatible disease/tissue/cell results appendix-only.",
            }
        )
    if int(prediction_quality.get("high_confidence_trait_score_weak_support_count") or 0) > 0:
        signals.append(
            {
                "code": "trait_score_confidence_overclaim",
                "severity": "error",
                "recommended_action": "Cap trait-score pathway confidence unless high-confidence rows have multiple stable seeds or literature support.",
            }
        )
    if int(status_counts.get("context_mismatch") or 0) > 0:
        signals.append(
            {
                "code": "gold_context_mismatch",
                "severity": "warning",
                "recommended_action": "Separate broad pan-cancer support from requested cancer/tissue/cell context.",
            }
        )
    evidence_usefulness = metrics.get("evidence_usefulness")
    matched_gold = (case.get("conclusion_evaluation") or {}).get("matched_gold") or []
    if evidence_usefulness is not None and float(evidence_usefulness or 0.0) <= 0.0 and matched_gold:
        signals.append(
            {
                "code": "evidence_not_effective",
                "severity": "warning",
                "recommended_action": "Verify that literature/evidence assertions are structurally linked and not merely present as text.",
            }
        )
    bug_codes = {
        str(item.get("code") or "")
        for item in ((case.get("conclusion_evaluation") or {}).get("bug_queue") or [])
        if item.get("code")
    }
    for code in sorted(bug_codes):
        signals.append(
            {
                "code": code,
                "severity": "error",
                "recommended_action": "Inspect the release gate bug queue entry for the matched candidate and provenance trace.",
            }
        )
    attach_optimization_targets(signals)
    deduped = {row["code"]: row for row in signals}
    return [deduped[key] for key in sorted(deduped)]


def case_signal_dimensions(case: dict[str, Any]) -> dict[str, list[str]]:
    tags = case.get("benchmark_tags") or {}
    return {
        "cancers": [str(value) for value in tags.get("cancer_ids", []) or [] if str(value)],
        "themes": [str(value) for value in tags.get("theme_ids", []) or [] if str(value)],
        "input_modes": [str(value) for value in tags.get("input_modes", []) or [] if str(value)],
    }


def add_signal_bucket(
    buckets: dict[str, dict[str, Any]],
    bucket_id: str,
    signal: dict[str, Any],
    case_id: str,
) -> None:
    entry = buckets.setdefault(
        bucket_id,
        {
            "bucket": bucket_id,
            "signal_count": 0,
            "cases": [],
            "codes": [],
            "optimization_targets": [],
            "max_severity": "",
        },
    )
    entry["signal_count"] += 1
    if case_id and case_id not in entry["cases"]:
        entry["cases"].append(case_id)
    code = str(signal.get("code") or "")
    if code and code not in entry["codes"]:
        entry["codes"].append(code)
    for target in signal.get("optimization_targets") or []:
        target_text = str(target)
        if target_text and target_text not in entry["optimization_targets"]:
            entry["optimization_targets"].append(target_text)
    if signal.get("severity") == "error" or not entry["max_severity"]:
        entry["max_severity"] = str(signal.get("severity") or "")


def aggregate_optimization_signals(scenarios: list[dict[str, Any]]) -> dict[str, Any]:
    by_code: dict[str, dict[str, Any]] = {}
    by_case: dict[str, list[dict[str, Any]]] = {}
    by_dimension: dict[str, dict[str, Any]] = {}
    by_target: dict[str, dict[str, Any]] = {}
    for case in scenarios:
        signals = case_optimization_signals(case)
        case_id = str(case.get("case_id") or "")
        dimensions = case_signal_dimensions(case)
        if signals:
            by_case[case_id] = signals
        for signal in signals:
            code = str(signal.get("code") or "unknown")
            entry = by_code.setdefault(
                code,
                {
                    "code": code,
                    "severity": signal.get("severity", ""),
                    "case_count": 0,
                    "cases": [],
                    "recommended_action": signal.get("recommended_action", ""),
                    "optimization_targets": signal.get("optimization_targets", []),
                },
            )
            entry["case_count"] += 1
            entry["cases"].append(case_id)
            for dimension_name, values in dimensions.items():
                for value in values:
                    add_signal_bucket(by_dimension, f"{dimension_name}:{value}", signal, case_id)
            for target in signal.get("optimization_targets") or []:
                add_signal_bucket(by_target, str(target), signal, case_id)
    return {
        "by_code": sorted(by_code.values(), key=lambda row: (row.get("severity") != "error", row["code"])),
        "by_case": by_case,
        "by_dimension": sorted(by_dimension.values(), key=lambda row: (row.get("max_severity") != "error", row["bucket"])),
        "by_target": sorted(by_target.values(), key=lambda row: (row.get("max_severity") != "error", row["bucket"])),
    }


def gold_coverage_signal(gold_coverage: dict[str, Any]) -> dict[str, Any] | None:
    if not gold_coverage or gold_coverage.get("status") == "covered":
        return None
    missing_count = int(gold_coverage.get("missing_cell_count") or 0)
    negative_status = (gold_coverage.get("negative_control_requirement") or {}).get("status", "")
    cases = [f"{row.get('cancer_id')}:{row.get('theme_id')}" for row in gold_coverage.get("missing_cells", [])[:10]]
    if negative_status != "covered":
        cases.append("negative_controls")
    return {
        "code": "gold_coverage_gap",
        "severity": "warning",
        "case_count": missing_count + (1 if negative_status != "covered" else 0),
        "cases": cases,
        "recommended_action": "Add reviewed positive controls for missing cancer/theme cells or update coverage requirements with an explicit rationale.",
    }


def gold_quality_signal(gold_quality: dict[str, Any]) -> dict[str, Any] | None:
    if not gold_quality or gold_quality.get("status") == "passed":
        return None
    failed = gold_quality.get("failed_checks") or []
    cases = [str(row.get("code") or "gold_quality_gap") for row in failed]
    return {
        "code": "gold_quality_gap",
        "severity": "error",
        "case_count": len(cases) or 1,
        "cases": cases[:10] or ["gold_quality"],
        "recommended_action": "Repair gold-standard quantity, required fields, source refs, or negative-trap category coverage before release.",
    }


def literature_semantics_signal(literature_semantics: dict[str, Any]) -> dict[str, Any] | None:
    if not literature_semantics or literature_semantics.get("status") == "passed":
        return None
    failed_cases = [str(row.get("case_id") or "") for row in literature_semantics.get("failed_cases", [])]
    metrics = literature_semantics.get("metrics") or {}
    return {
        "code": "literature_semantics_regression",
        "severity": "error",
        "case_count": int(metrics.get("failed_case_count") or len(failed_cases) or 1),
        "cases": failed_cases[:10] or ["literature_semantics_benchmark"],
        "recommended_action": "Fix assertion semantics before release; do not lower thresholds or promote weak/negated/background literature into core support.",
    }


def audit_fixture_availability(
    scenarios: list[dict[str, Any]],
    required_fixture_cases: list[str],
    require_real_fixtures: bool,
) -> dict[str, Any]:
    by_case = {str(case.get("case_id") or ""): case for case in scenarios}
    required = [case_id for case_id in required_fixture_cases if case_id]
    missing_required = []
    available_required = []
    for case_id in required:
        case = by_case.get(case_id)
        if not case or case.get("status") == "missing":
            missing_required.append(
                {
                    "case_id": case_id,
                    "path": str((case or {}).get("path") or ""),
                    "message": str((case or {}).get("message") or "Required real-world fixture case is not registered."),
                    "expected_schema": (case or {}).get("expected_schema") or {},
                    "candidate_paths": (case or {}).get("candidate_paths") or [],
                    "recovery_action": str((case or {}).get("recovery_action") or ""),
                    "recovery_rerun_plan": (case or {}).get("recovery_rerun_plan") or {},
                }
            )
        else:
            available_required.append({"case_id": case_id, "status": case.get("status", ""), "path": str(case.get("path") or "")})
    failed = require_real_fixtures and bool(missing_required)
    return {
        "status": "failed" if failed else "passed",
        "require_real_fixtures": require_real_fixtures,
        "required_fixture_cases": required,
        "available_required_count": len(available_required),
        "missing_required_count": len(missing_required),
        "available_required": available_required,
        "missing_required": missing_required,
        "interpretation": (
            "Missing real-world fixtures are warnings in exploratory runs; with --require-real-fixtures they block release."
        ),
    }


def fixture_audit_signal(fixture_audit: dict[str, Any]) -> dict[str, Any] | None:
    if not fixture_audit or fixture_audit.get("status") == "passed":
        return None
    missing = fixture_audit.get("missing_required") or []
    return {
        "code": "required_fixture_missing",
        "severity": "error",
        "case_count": len(missing) or 1,
        "cases": [str(row.get("case_id") or "") for row in missing[:10]] or ["real_world_fixture"],
        "recommended_action": "Restore the required real-world fixture file and rerun the release benchmark with --require-real-fixtures.",
    }


def as_optional_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def heldout_quality_findings(
    scope: str,
    metrics: dict[str, Any],
    label: str = "",
    diagnostics: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    findings: list[dict[str, Any]] = []
    if not isinstance(metrics, dict) or not metrics:
        return findings

    diagnostics = diagnostics or {}
    heldout_rows = int(metrics.get("heldout_rows") or 0)
    positive_rows = int(metrics.get("heldout_positive_rows") or 0)
    nonpositive_rows = int(metrics.get("heldout_nonpositive_rows") or 0)
    context = {
        "scope": scope,
        "label": label,
        "heldout_rows": heldout_rows,
        "heldout_positive_rows": positive_rows,
        "diagnostics": diagnostics,
    }
    if heldout_rows <= 0 or positive_rows <= 0:
        return findings

    if nonpositive_rows > 0:
        lift_at_100 = as_optional_float(metrics.get("heldout_lift_at_100"))
        if lift_at_100 is not None and lift_at_100 <= DEFAULT_MIN_HELDOUT_LIFT_AT_100:
            findings.append(
                {
                    "code": f"{scope}_heldout_no_topk_lift",
                    "severity": "error",
                    "message": f"{scope} held-out validation does not beat top-k prevalence baseline.",
                    "metrics": {**context, "heldout_lift_at_100": lift_at_100, "minimum_exclusive": DEFAULT_MIN_HELDOUT_LIFT_AT_100},
                }
            )
        roc_auc = as_optional_float(metrics.get("heldout_roc_auc"))
        if roc_auc is not None and roc_auc <= DEFAULT_MIN_HELDOUT_ROC_AUC:
            severity = "error" if lift_at_100 is not None and lift_at_100 <= DEFAULT_MIN_HELDOUT_LIFT_AT_100 else "warning"
            findings.append(
                {
                    "code": f"{scope}_heldout_random_auc",
                    "severity": severity,
                    "message": f"{scope} held-out validation has random-or-worse ROC-AUC.",
                    "metrics": {**context, "heldout_roc_auc": roc_auc, "minimum_exclusive": DEFAULT_MIN_HELDOUT_ROC_AUC},
                }
            )

    unique_scores = as_optional_float(metrics.get("heldout_unique_transfer_scores"))
    if unique_scores is not None and unique_scores < DEFAULT_MIN_HELDOUT_UNIQUE_TRANSFER_SCORES:
        findings.append(
            {
                "code": f"{scope}_heldout_score_collapse",
                "severity": "error",
                "message": f"{scope} held-out validation collapsed to too few transfer scores.",
                "metrics": {
                    **context,
                    "heldout_unique_transfer_scores": unique_scores,
                    "minimum": DEFAULT_MIN_HELDOUT_UNIQUE_TRANSFER_SCORES,
                },
            }
        )
    tie_fraction = as_optional_float(metrics.get("heldout_top_score_tie_fraction"))
    if tie_fraction is not None and tie_fraction >= DEFAULT_MAX_HELDOUT_TOP_SCORE_TIE_FRACTION:
        findings.append(
            {
                "code": f"{scope}_heldout_score_tie_collapse",
                "severity": "error",
                "message": f"{scope} held-out validation top scores are dominated by ties.",
                "metrics": {
                    **context,
                    "heldout_top_score_tie_fraction": tie_fraction,
                    "maximum_exclusive": DEFAULT_MAX_HELDOUT_TOP_SCORE_TIE_FRACTION,
                },
            }
        )
    return findings


def audit_heldout_validation(
    heldout: dict[str, Any],
    rerun: dict[str, Any] | None,
    max_age_hours: int,
    require_rerun: bool,
    min_temporal_rows: int,
    min_source_count: int,
    min_source_rows: int,
) -> dict[str, Any]:
    now = datetime.now(timezone.utc)
    failed_checks: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []
    temporal = heldout.get("temporal_holdout") or {}
    source = heldout.get("source_heldout") or {}
    expected_rerun_names = ("temporal", "source")
    invalid_rerun_names = [name for name, payload in (rerun or {}).items() if not isinstance(payload, dict)]
    rerun_payloads = {name: payload for name, payload in (rerun or {}).items() if isinstance(payload, dict)}
    rerun_passed = {name: bool(payload.get("passed")) for name, payload in rerun_payloads.items()}
    rerun_performed = all(name in rerun_payloads for name in expected_rerun_names)
    rerun_details = {
        name: {
            "command": list(payload.get("command") or []),
            "started_at_utc": payload.get("started_at_utc", ""),
            "completed_at_utc": payload.get("completed_at_utc", ""),
            "duration_seconds": payload.get("duration_seconds"),
            "returncode": payload.get("returncode"),
            "passed": bool(payload.get("passed")),
            "timed_out": bool(payload.get("timed_out")),
            "timeout_seconds": payload.get("timeout_seconds"),
            "report_path": payload.get("report_path", ""),
        }
        for name, payload in rerun_payloads.items()
    }
    report_freshness: dict[str, dict[str, Any]] = {}

    if require_rerun and not rerun_performed:
        failed_checks.append(
            {
                "code": "heldout_rerun_required",
                "message": "Release benchmark was run without complete temporal and source --rerun-heldout results while --require-heldout-rerun was set.",
                "missing": [name for name in expected_rerun_names if name not in rerun_payloads],
            }
        )
    for name in invalid_rerun_names:
        failed_checks.append(
            {
                "code": f"heldout_rerun_invalid:{name}",
                "message": f"{name} held-out rerun payload is not a structured object.",
            }
        )
    for name, passed in rerun_passed.items():
        if not passed:
            timed_out = bool((rerun_payloads.get(name) or {}).get("timed_out"))
            code = f"heldout_rerun_timeout:{name}" if timed_out else f"heldout_rerun_failed:{name}"
            message = f"{name} held-out rerun timed out." if timed_out else f"{name} held-out rerun failed."
            failed_checks.append({"code": code, "message": message})

    for name, payload in (("temporal", temporal), ("source", source)):
        if not payload.get("available"):
            failed_checks.append({"code": f"{name}_heldout_missing", "message": f"{name} held-out report is missing."})
            continue
        created = parse_utc_datetime(payload.get("created_at_utc"))
        if created is None:
            failed_checks.append({"code": f"{name}_heldout_unparseable_timestamp", "message": f"{name} held-out timestamp is missing or invalid."})
            report_freshness[name] = {
                "available": True,
                "created_at_utc": payload.get("created_at_utc", ""),
                "age_hours": None,
                "status": "unparseable_timestamp",
            }
            continue
        age_hours = (now - created).total_seconds() / 3600.0
        payload["age_hours"] = round(age_hours, 3)
        freshness_status = "fresh"
        if max_age_hours > 0 and age_hours > max_age_hours:
            freshness_status = "stale"
            failed_checks.append(
                {
                    "code": f"{name}_heldout_stale",
                    "message": f"{name} held-out report is {age_hours:.1f} hours old, above {max_age_hours} hours.",
                }
            )
        report_freshness[name] = {
            "available": True,
            "created_at_utc": payload.get("created_at_utc", ""),
            "age_hours": round(age_hours, 3),
            "max_age_hours": max_age_hours,
            "status": freshness_status,
        }

    temporal_metrics = temporal.get("validation_metrics") or {}
    temporal_rows = int(temporal_metrics.get("heldout_rows") or 0)
    if temporal.get("available") and temporal_rows < min_temporal_rows:
        failed_checks.append(
            {
                "code": "temporal_heldout_too_small",
                "message": f"Temporal held-out rows {temporal_rows} are below required {min_temporal_rows}.",
            }
        )
    if temporal.get("available") and int(temporal_metrics.get("heldout_positive_rows") or 0) <= 0:
        failed_checks.append({"code": "temporal_heldout_no_positives", "message": "Temporal held-out report has no positive rows."})
    for finding in heldout_quality_findings("temporal", temporal_metrics):
        if finding.get("severity") == "warning":
            warnings.append(finding)
        else:
            failed_checks.append(finding)

    source_summary = source.get("summary") or []
    source_count = int(source.get("source_count") or 0)
    if source.get("available") and source_count < min_source_count:
        failed_checks.append(
            {
                "code": "source_heldout_source_count_low",
                "message": f"Source held-out source count {source_count} is below required {min_source_count}.",
            }
        )
    small_sources = [
        str(row.get("source") or "unknown")
        for row in source_summary
        if int(row.get("heldout_rows") or 0) < min_source_rows
    ]
    if source.get("available") and small_sources:
        warnings.append(
            {
                "code": "source_heldout_small_sources",
                "message": f"Some held-out sources have fewer than {min_source_rows} rows.",
                "sources": small_sources[:10],
            }
        )
    for row in source_summary:
        source_label = str(row.get("source") or "unknown")
        diagnostics = {
            "source_family": row.get("source_family_diagnostics") or {},
            "feature_transfer": row.get("feature_transfer_diagnostics") or {},
        }
        for finding in heldout_quality_findings("source", row.get("validation_metrics") or {}, source_label, diagnostics):
            if finding.get("severity") == "warning":
                warnings.append(finding)
            else:
                failed_checks.append(finding)

    return {
        "status": "passed" if not failed_checks else "failed",
        "rerun_performed": rerun_performed,
        "rerun_passed": rerun_passed,
        "rerun_details": rerun_details,
        "report_freshness": report_freshness,
        "require_rerun": require_rerun,
        "max_age_hours": max_age_hours,
        "min_temporal_rows": min_temporal_rows,
        "min_source_count": min_source_count,
        "min_source_rows": min_source_rows,
        "failed_checks": failed_checks,
        "warnings": warnings,
    }


def heldout_validation_signal(heldout_validation: dict[str, Any]) -> dict[str, Any] | None:
    if not heldout_validation or heldout_validation.get("status") == "passed":
        return None
    failed = heldout_validation.get("failed_checks") or []
    cases = [str(row.get("code") or "heldout_validation_gap") for row in failed]
    return {
        "code": "heldout_validation_gap",
        "severity": "error",
        "case_count": len(cases) or 1,
        "cases": cases[:10] or ["heldout_validation"],
        "recommended_action": "Rerun temporal/source held-out validation or repair missing/stale/undersized held-out reports before release.",
    }


def tag_values(case: dict[str, Any], key: str) -> list[str]:
    tags = case.get("benchmark_tags") or {}
    value = tags.get(key) if isinstance(tags, dict) else []
    if isinstance(value, list):
        return [str(item) for item in value if str(item)]
    return [str(value)] if str(value) else []


def coverage_dimension(required: list[str], registered: set[str], executed: set[str]) -> dict[str, Any]:
    required_set = set(required)
    return {
        "required": sorted(required_set),
        "registered": sorted(required_set & registered),
        "executed": sorted(required_set & executed),
        "missing_registered": sorted(required_set - registered),
        "missing_executed": sorted(required_set - executed),
        "registered_status": "covered" if not (required_set - registered) else "gaps_present",
        "executed_status": "covered" if not (required_set - executed) else "gaps_present",
    }


def audit_benchmark_matrix(
    scenarios: list[dict[str, Any]],
    required_cancers: list[str],
    required_themes: list[str],
    required_input_modes: list[str],
    require_executed_coverage: bool,
) -> dict[str, Any]:
    registered_cancers: set[str] = set()
    registered_themes: set[str] = set()
    registered_input_modes: set[str] = set()
    executed_cancers: set[str] = set()
    executed_themes: set[str] = set()
    executed_input_modes: set[str] = set()
    matrix_cases = []
    for case in scenarios:
        cancers = set(tag_values(case, "cancer_ids"))
        themes = set(tag_values(case, "theme_ids"))
        input_modes = set(tag_values(case, "input_modes"))
        registered_cancers.update(cancers)
        registered_themes.update(themes)
        registered_input_modes.update(input_modes)
        executed = case.get("status") != "missing"
        if executed:
            executed_cancers.update(cancers)
            executed_themes.update(themes)
            executed_input_modes.update(input_modes)
        matrix_cases.append(
            {
                "case_id": case.get("case_id", ""),
                "status": case.get("status", ""),
                "cancer_ids": sorted(cancers),
                "theme_ids": sorted(themes),
                "input_modes": sorted(input_modes),
                "executed": executed,
            }
        )
    dimensions = {
        "cancers": coverage_dimension(required_cancers, registered_cancers, executed_cancers),
        "themes": coverage_dimension(required_themes, registered_themes, executed_themes),
        "input_modes": coverage_dimension(required_input_modes, registered_input_modes, executed_input_modes),
    }
    registered_gaps = [
        f"{name}:{value}"
        for name, dimension in dimensions.items()
        for value in dimension["missing_registered"]
    ]
    executed_gaps = [
        f"{name}:{value}"
        for name, dimension in dimensions.items()
        for value in dimension["missing_executed"]
    ]
    failed = bool(registered_gaps or (require_executed_coverage and executed_gaps))
    return {
        "status": "passed" if not failed else "failed",
        "require_executed_coverage": require_executed_coverage,
        "registered_gap_count": len(registered_gaps),
        "executed_gap_count": len(executed_gaps),
        "registered_gaps": registered_gaps,
        "executed_gaps": executed_gaps,
        "dimensions": dimensions,
        "cases": matrix_cases,
        "interpretation": (
            "Registered coverage includes unavailable real-world fixtures; executed coverage counts only scenarios that actually ran. "
            "Strict release requires executed matrix coverage, while exploratory runs report execution gaps for fixture follow-up."
        ),
    }


def benchmark_matrix_signal(benchmark_matrix: dict[str, Any]) -> dict[str, Any] | None:
    if not benchmark_matrix or benchmark_matrix.get("status") == "passed":
        return None
    cases = list(benchmark_matrix.get("registered_gaps") or []) + list(benchmark_matrix.get("executed_gaps") or [])
    return {
        "code": "benchmark_matrix_gap",
        "severity": "error",
        "case_count": len(cases) or 1,
        "cases": cases[:10] or ["benchmark_matrix"],
        "recommended_action": "Add or restore benchmark scenarios so required cancer, theme, and input-mode dimensions are registered and executed for strict release.",
    }


def traceability_contract_signal(traceability_contract: dict[str, Any]) -> dict[str, Any] | None:
    if not traceability_contract or traceability_contract.get("status") == "passed":
        return None
    failed_cases = traceability_contract.get("failed_cases") or []
    return {
        "code": "scenario_traceability_contract_failed",
        "severity": "error",
        "case_count": len(failed_cases) or 1,
        "cases": [str(row.get("case_id") or "") for row in failed_cases[:10]] or ["scenario_traceability"],
        "recommended_action": (
            "Repair scenario output contracts so benchmark conclusions link to input summaries, identity decisions, "
            "mechanism facts, evidence assertions, and conclusion-evaluation metrics."
        ),
    }


def merge_report_signals(scenario_signals: dict[str, Any], extra_signals: list[dict[str, Any]]) -> dict[str, Any]:
    by_case = dict(scenario_signals.get("by_case") or {})
    by_code = {str(row.get("code") or ""): dict(row) for row in scenario_signals.get("by_code", []) or []}
    by_dimension = {str(row.get("bucket") or ""): dict(row) for row in scenario_signals.get("by_dimension", []) or []}
    by_target = {str(row.get("bucket") or ""): dict(row) for row in scenario_signals.get("by_target", []) or []}
    for signal in extra_signals:
        if not signal:
            continue
        signal = dict(signal)
        attach_optimization_targets([signal])
        code = str(signal.get("code") or "")
        if not code:
            continue
        if code in by_code:
            existing = by_code[code]
            existing["case_count"] = int(existing.get("case_count") or 0) + int(signal.get("case_count") or 0)
            existing["cases"] = list(existing.get("cases") or []) + list(signal.get("cases") or [])
            existing_targets = set(existing.get("optimization_targets") or [])
            existing["optimization_targets"] = sorted(existing_targets | set(signal.get("optimization_targets") or []))
        else:
            by_code[code] = dict(signal)
        target_case_refs = [str(case_ref) for case_ref in signal.get("cases") or [] if str(case_ref)] or [code]
        for target in signal.get("optimization_targets") or []:
            for case_ref in target_case_refs:
                add_signal_bucket(by_target, str(target), signal, case_ref)
        for case_ref in signal.get("cases") or []:
            case_text = str(case_ref)
            if case_text.startswith(("cancers:", "themes:", "input_modes:")):
                add_signal_bucket(by_dimension, case_text, signal, case_text)
    return {
        "by_case": by_case,
        "by_code": sorted(by_code.values(), key=lambda row: (row.get("severity") != "error", row.get("code", ""))),
        "by_dimension": sorted(
            [row for key, row in by_dimension.items() if key],
            key=lambda row: (row.get("max_severity") != "error", row.get("bucket", "")),
        ),
        "by_target": sorted(
            [row for key, row in by_target.items() if key],
            key=lambda row: (row.get("max_severity") != "error", row.get("bucket", "")),
        ),
    }


def summarize_existing_heldout(workspace: Path, run_id: str) -> dict[str, Any]:
    report_root = workspace / "learning_runs" / run_id / "reports"
    temporal = read_json(report_root / "temporal_holdout" / "temporal_holdout_validation_report.json")
    source = read_json(report_root / "source_heldout" / "source_heldout_validation_report.json")
    release_gate = read_json(report_root / "release_gate_validation_report.json")
    return {
        "run_id": run_id,
        "temporal_holdout": {
            "available": bool(temporal),
            "created_at_utc": temporal.get("created_at_utc", ""),
            "validation_metrics": temporal.get("validation_metrics", {}),
        },
        "source_heldout": {
            "available": bool(source),
            "created_at_utc": source.get("created_at_utc", ""),
            "source_count": len(source.get("heldout_sources", []) or []),
            "summary": [
                {
                    "source": ((row.get("heldout") or {}).get("value") or ""),
                    "heldout_rows": row.get("heldout_rows"),
                    "validation_metrics": row.get("validation_metrics", {}),
                    "source_family_diagnostics": row.get("source_family_diagnostics", {}),
                    "feature_transfer_diagnostics": row.get("feature_transfer_diagnostics", {}),
                }
                for row in (source.get("heldout_sources", []) or [])[:10]
            ],
        },
        "release_gate": {
            "available": bool(release_gate),
            "gate_status": release_gate.get("gate_status", ""),
            "summary": release_gate.get("summary", {}),
        },
    }


def run_subprocess(command: list[str], workspace: Path, timeout_seconds: int) -> dict[str, Any]:
    started_dt = datetime.now(timezone.utc).replace(microsecond=0)
    started = started_dt.isoformat().replace("+00:00", "Z")
    env = os.environ.copy()
    env.setdefault("LOKY_MAX_CPU_COUNT", "1")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("NUMEXPR_NUM_THREADS", "1")
    try:
        completed = subprocess.run(
            command,
            cwd=workspace,
            capture_output=True,
            text=True,
            env=env,
            timeout=timeout_seconds if timeout_seconds > 0 else None,
        )
    except subprocess.TimeoutExpired as exc:
        completed_dt = datetime.now(timezone.utc).replace(microsecond=0)
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        return {
            "command": command,
            "started_at_utc": started,
            "completed_at_utc": completed_dt.isoformat().replace("+00:00", "Z"),
            "duration_seconds": round((completed_dt - started_dt).total_seconds(), 3),
            "returncode": None,
            "stdout": str(stdout)[-4000:],
            "stderr": str(stderr)[-4000:],
            "passed": False,
            "timed_out": True,
            "timeout_seconds": timeout_seconds,
        }
    completed_dt = datetime.now(timezone.utc).replace(microsecond=0)
    return {
        "command": command,
        "started_at_utc": started,
        "completed_at_utc": completed_dt.isoformat().replace("+00:00", "Z"),
        "duration_seconds": round((completed_dt - started_dt).total_seconds(), 3),
        "returncode": completed.returncode,
        "stdout": completed.stdout[-4000:],
        "stderr": completed.stderr[-4000:],
        "passed": completed.returncode == 0,
        "timed_out": False,
        "timeout_seconds": timeout_seconds,
    }


def rerun_heldout(workspace: Path, run_id: str, top_n: int, timeout_seconds: int) -> dict[str, Any]:
    python = sys.executable
    report_root = workspace / "learning_runs" / run_id / "reports"
    temporal_report = report_root / "temporal_holdout" / "temporal_holdout_validation_report.json"
    source_report = report_root / "source_heldout" / "source_heldout_validation_report.json"
    temporal = run_subprocess(
        [python, str(SCRIPT_DIR / "run_temporal_holdout_validation.py"), "--workspace", ".", "--run-id", run_id, "--top-n", str(top_n)],
        workspace,
        timeout_seconds,
    )
    temporal["report_path"] = str(temporal_report)
    source = run_subprocess(
        [python, str(SCRIPT_DIR / "run_source_heldout_validation.py"), "--workspace", ".", "--run-id", run_id, "--top-n", str(top_n)],
        workspace,
        timeout_seconds,
    )
    source["report_path"] = str(source_report)
    return {
        "temporal": temporal,
        "source": source,
    }


def scenario_timeout_result(case: dict[str, Any], timeout_seconds: int, duration_seconds: float, stdout: str = "", stderr: str = "") -> dict[str, Any]:
    return {
        "case_id": case.get("case_id", ""),
        "status": "failed",
        "error_type": "scenario_timeout",
        "error_message": f"Scenario exceeded timeout of {timeout_seconds} seconds.",
        "duration_seconds": round(duration_seconds, 3),
        "timeout_seconds": timeout_seconds,
        "timed_out": True,
        "benchmark_tags": case.get("benchmark_tags") or {},
        "input_summary": {"input_count": len(case.get("records", []) or [])},
        "expected_identity_behavior": case.get("expected_identity_behavior") or {},
        "scenario_expected_metrics": {
            "expected_gold_ids": list(case.get("expected_gold_ids", []) or []),
            "expected_count": len(case.get("expected_gold_ids", []) or []),
            "expected_matched": 0,
            "expected_recall": 0.0 if case.get("expected_gold_ids") else None,
            "expected_statuses": {},
            "missing_expected_gold_ids": list(case.get("expected_gold_ids", []) or []),
        },
        "traceability": {
            "contract_version": TRACEABILITY_CONTRACT_VERSION,
            "input_summary_present": True,
            "identity_decision_count": 0,
            "mechanism_ready_fact_count": 0,
            "evidence_assertion_count": 0,
            "evidence_support_statuses": [],
            "conclusion_evaluation_present": False,
            "evaluation_metrics_present": False,
            "release_gate_present": False,
            "analysis_pack_hash_present": False,
            "response_hash_present": False,
        },
        "conclusion_evaluation": {
            "metrics": {
                "gold_positive_count": len(case.get("expected_gold_ids", []) or []),
                "gold_positive_matched": 0,
                "negative_trap_hits": 0,
                "unsupported_top_claim_rate": 1.0,
                "evidence_usefulness": 0.0,
            },
            "bug_queue": [
                {
                    "code": "scenario_timeout",
                    "message": "Scenario analysis timed out and must be profiled or batched before strict release.",
                }
            ],
            "per_gold_status_counts": {},
            "matched_gold": [],
        },
        "stdout_tail": stdout[-4000:],
        "stderr_tail": stderr[-4000:],
    }


def scenario_child_failure_result(case: dict[str, Any], completed: subprocess.CompletedProcess[str], duration_seconds: float) -> dict[str, Any]:
    return {
        "case_id": case.get("case_id", ""),
        "status": "failed",
        "error_type": "scenario_child_failed",
        "error_message": f"Scenario subprocess exited with code {completed.returncode}.",
        "duration_seconds": round(duration_seconds, 3),
        "returncode": completed.returncode,
        "benchmark_tags": case.get("benchmark_tags") or {},
        "input_summary": {"input_count": len(case.get("records", []) or [])},
        "expected_identity_behavior": case.get("expected_identity_behavior") or {},
        "conclusion_evaluation": {
            "metrics": {"negative_trap_hits": 0, "unsupported_top_claim_rate": 1.0, "evidence_usefulness": 0.0},
            "bug_queue": [{"code": "scenario_child_failed", "message": "Scenario subprocess failed before producing a valid summary."}],
            "per_gold_status_counts": {},
            "matched_gold": [],
        },
        "stdout_tail": (completed.stdout or "")[-4000:],
        "stderr_tail": (completed.stderr or "")[-4000:],
    }


def run_scenario_child_main(args: argparse.Namespace) -> int:
    payload = read_json(Path(args.single_scenario_input))
    workspace = Path(payload.get("workspace") or args.workspace).resolve()
    release_id = str(payload.get("release_id") or args.release_id or latest_release_id(workspace / "graph_projection"))
    case = payload.get("case") or {}
    service = MetaboService(workspace, release_id=release_id)
    try:
        analyzed = service.analyze_metabolites(
            case.get("records", []),
            max_paths=int(payload.get("max_paths") or args.max_paths),
            max_hops=int(payload.get("max_hops") or args.max_hops),
            context=case.get("context") or {},
        )
        result = summarize_analysis(case, analyzed)
    except Exception as exc:  # pragma: no cover - child diagnostics are integration exercised.
        result = {
            "case_id": case.get("case_id", ""),
            "status": "failed",
            "error_type": type(exc).__name__,
            "error_message": str(exc),
            "benchmark_tags": case.get("benchmark_tags") or {},
            "input_summary": {"input_count": len(case.get("records", []) or [])},
        }
    write_json(Path(args.single_scenario_output), result)
    return 0


def run_ready_scenario_isolated(
    case: dict[str, Any],
    *,
    index: int,
    workspace: Path,
    release_id: str,
    max_paths: int,
    max_hops: int,
    timeout_seconds: int,
    run_dir: Path,
) -> dict[str, Any]:
    run_dir.mkdir(parents=True, exist_ok=True)
    slug = slugify(str(case.get("case_id") or f"scenario_{index}"))
    input_path = run_dir / f"{index:03d}_{slug}_input.json"
    output_path = run_dir / f"{index:03d}_{slug}_output.json"
    write_json(
        input_path,
        {
            "workspace": str(workspace),
            "release_id": release_id,
            "max_paths": max_paths,
            "max_hops": max_hops,
            "case": case,
        },
    )
    command = [
        sys.executable,
        str(SCRIPT_DIR / "run_release_batch_benchmark.py"),
        "--workspace",
        str(workspace),
        "--release-id",
        release_id,
        "--max-paths",
        str(max_paths),
        "--max-hops",
        str(max_hops),
        "--single-scenario-input",
        str(input_path),
        "--single-scenario-output",
        str(output_path),
    ]
    started = datetime.now(timezone.utc)
    try:
        completed = subprocess.run(
            command,
            cwd=workspace,
            capture_output=True,
            text=True,
            timeout=timeout_seconds if timeout_seconds > 0 else None,
        )
    except subprocess.TimeoutExpired as exc:
        duration = (datetime.now(timezone.utc) - started).total_seconds()
        return scenario_timeout_result(case, timeout_seconds, duration, str(exc.stdout or ""), str(exc.stderr or ""))
    duration = (datetime.now(timezone.utc) - started).total_seconds()
    if completed.returncode != 0:
        return scenario_child_failure_result(case, completed, duration)
    result = read_json(output_path)
    if not result:
        return scenario_child_failure_result(case, completed, duration)
    result["duration_seconds"] = round(duration, 3)
    result["timeout_seconds"] = timeout_seconds
    result["timed_out"] = False
    result["scenario_run_artifacts"] = {"input": str(input_path), "output": str(output_path)}
    return result


def benchmark_status(
    scenarios: list[dict[str, Any]],
    heldout: dict[str, Any],
    rerun: dict[str, Any] | None,
    literature_semantics: dict[str, Any] | None = None,
    heldout_validation: dict[str, Any] | None = None,
    fixture_audit: dict[str, Any] | None = None,
    benchmark_matrix: dict[str, Any] | None = None,
    traceability_contract: dict[str, Any] | None = None,
) -> dict[str, Any]:
    failed_cases = [case for case in scenarios if case.get("status") == "failed"]
    missing_cases = [case for case in scenarios if case.get("status") == "missing"]
    heldout_available = bool((heldout.get("temporal_holdout") or {}).get("available")) and bool((heldout.get("source_heldout") or {}).get("available"))
    rerun_failed = []
    if rerun:
        rerun_failed = [name for name, payload in rerun.items() if not payload.get("passed")]
    literature_semantics_passed = not literature_semantics or literature_semantics.get("status") == "passed"
    heldout_validation_passed = not heldout_validation or heldout_validation.get("status") == "passed"
    fixture_audit_passed = not fixture_audit or fixture_audit.get("status") == "passed"
    benchmark_matrix_passed = not benchmark_matrix or benchmark_matrix.get("status") == "passed"
    traceability_contract_passed = not traceability_contract or traceability_contract.get("status") == "passed"
    return {
        "passed": (
            not failed_cases
            and heldout_available
            and not rerun_failed
            and literature_semantics_passed
            and heldout_validation_passed
            and fixture_audit_passed
            and benchmark_matrix_passed
            and traceability_contract_passed
        ),
        "failed_case_count": len(failed_cases),
        "missing_case_count": len(missing_cases),
        "heldout_reports_available": heldout_available,
        "rerun_failed": rerun_failed,
        "literature_semantics_passed": literature_semantics_passed,
        "heldout_validation_passed": heldout_validation_passed,
        "fixture_audit_passed": fixture_audit_passed,
        "benchmark_matrix_passed": benchmark_matrix_passed,
        "traceability_contract_passed": traceability_contract_passed,
        "interpretation": (
            "Missing optional real-world fixture cases are reported for follow-up but do not fail the automated gate; "
            "failed executed cases, missing held-out reports, held-out rerun failures, literature semantics regressions, "
            "required benchmark matrix gaps, traceability contract gaps, or required real-world fixture gaps do fail it."
        ),
    }


def write_markdown(path: Path, report: dict[str, Any]) -> None:
    lines = [
        "# Release Batch Benchmark",
        "",
        f"Release: `{report.get('release_id', '')}`",
        f"Created: `{report.get('created_at_utc', '')}`",
        f"Status: `{'passed' if report.get('status', {}).get('passed') else 'failed'}`",
        f"Profile: `{(report.get('release_profile') or {}).get('profile', '')}`",
        f"Progress log: `{report.get('progress_log', '')}`",
        "",
        "## Scenario Summary",
        "",
        "| case | status | matched | ambiguous | identity expectation | negative traps | gold recall | expected recall |",
        "| --- | --- | ---: | ---: | --- | ---: | ---: | ---: |",
    ]
    for case in report.get("scenarios", []):
        summary = case.get("input_summary") or {}
        metrics = ((case.get("conclusion_evaluation") or {}).get("metrics") or {})
        identity_expectation = "allow_ambiguity" if (case.get("expected_identity_behavior") or {}).get("allow_ambiguity") else ""
        lines.append(
            "| {case} | {status} | {matched} | {ambiguous} | {identity_expectation} | {traps} | {recall} | {expected_recall} |".format(
                case=case.get("case_id", ""),
                status=case.get("status", ""),
                matched=summary.get("matched_count", ""),
                ambiguous=summary.get("ambiguous_count", ""),
                identity_expectation=identity_expectation,
                traps=metrics.get("negative_trap_hits", ""),
                recall=metrics.get("gold_recall", ""),
                expected_recall=(case.get("scenario_expected_metrics") or {}).get("expected_recall", ""),
            )
        )
    prediction_quality_cases = [
        case
        for case in report.get("scenarios", [])
        if case.get("prediction_quality")
        and (
            int((case.get("prediction_quality") or {}).get("high_confidence_trait_score_weak_support_count") or 0) > 0
            or int((case.get("prediction_quality") or {}).get("context_mismatch_prediction_count") or 0) > 0
        )
    ]
    if prediction_quality_cases:
        lines.extend(
            [
                "",
                "## Prediction Quality Review",
                "",
                "| case | mode | weak trait-score high confidence | context mismatch predictions | high/medium context mismatch |",
                "| --- | --- | ---: | ---: | ---: |",
            ]
        )
        for case in prediction_quality_cases:
            quality = case.get("prediction_quality") or {}
            lines.append(
                "| {case} | {mode} | {weak} | {mismatch} | {high_medium} |".format(
                    case=case.get("case_id", ""),
                    mode=quality.get("analysis_mode", ""),
                    weak=quality.get("high_confidence_trait_score_weak_support_count", 0),
                    mismatch=quality.get("context_mismatch_prediction_count", 0),
                    high_medium=quality.get("high_or_medium_context_mismatch_prediction_count", 0),
                )
            )
        weak_rows = [
            (case.get("case_id", ""), row)
            for case in prediction_quality_cases
            for row in (case.get("prediction_quality") or {}).get("high_confidence_trait_score_weak_support", [])[:10]
        ]
        if weak_rows:
            lines.extend(["", "Weak high-confidence trait-score pathway rows:"])
            for case_id, row in weak_rows[:12]:
                lines.append(
                    f"- `{case_id}` rank `{row.get('rank', '')}` {row.get('display_name', '')}: "
                    f"overlap={row.get('overlap_count', 0)}, significant={row.get('significant_support_count', 0)}, "
                    f"literature={row.get('literature_support_count', 0)}"
                )
        mismatch_rows = [
            (case.get("case_id", ""), row)
            for case in prediction_quality_cases
            for row in (case.get("prediction_quality") or {}).get("context_mismatch_predictions", [])[:10]
        ]
        if mismatch_rows:
            lines.extend(["", "Context-mismatch prediction rows:"])
            for case_id, row in mismatch_rows[:12]:
                lines.append(
                    f"- `{case_id}` {row.get('ranking', '')} rank `{row.get('rank', '')}` "
                    f"{row.get('display_name', '')}: tier=`{row.get('confidence_tier', '')}`, "
                    f"status=`{row.get('calibration_status', '')}`"
                )
    events = report.get("progress_events") or []
    if events:
        lines.extend(["", "## Progress Events", ""])
        for row in events[-20:]:
            event_case = f" `{row.get('case_id')}`" if row.get("case_id") else ""
            status = f" status=`{row.get('status')}`" if row.get("status") else ""
            duration = f" duration=`{row.get('duration_seconds')}`" if row.get("duration_seconds") is not None else ""
            error = f" error=`{row.get('error_type')}`" if row.get("error_type") else ""
            lines.append(f"- `{row.get('event')}`{event_case}{status}{duration}{error}")
    missing = [case for case in report.get("scenarios", []) if case.get("status") == "missing"]
    if missing:
        lines.extend(["", "## Missing Fixtures", ""])
        for case in missing:
            lines.append(f"- `{case.get('case_id')}`: {case.get('message', '')} Path: `{case.get('path', '')}`")
    fixture_audit = report.get("fixture_audit") or {}
    if fixture_audit:
        lines.extend(
            [
                "",
                "## Fixture Audit",
                "",
                f"- Status: `{fixture_audit.get('status', '')}`",
                f"- Require real fixtures: `{fixture_audit.get('require_real_fixtures')}`",
                f"- Available required: `{fixture_audit.get('available_required_count', 0)}/{len(fixture_audit.get('required_fixture_cases') or [])}`",
            ]
        )
        missing_required = fixture_audit.get("missing_required") or []
        if missing_required:
            lines.extend(["", "Missing required real-world fixtures:"])
            for row in missing_required[:12]:
                lines.append(f"- `{row.get('case_id')}`: {row.get('message', '')} Path: `{row.get('path', '')}`")
    benchmark_matrix = report.get("benchmark_matrix") or {}
    if benchmark_matrix:
        lines.extend(
            [
                "",
                "## Benchmark Matrix",
                "",
                f"- Status: `{benchmark_matrix.get('status', '')}`",
                f"- Require executed coverage: `{benchmark_matrix.get('require_executed_coverage')}`",
                f"- Registered gaps: `{benchmark_matrix.get('registered_gap_count', 0)}`",
                f"- Executed gaps: `{benchmark_matrix.get('executed_gap_count', 0)}`",
            ]
        )
        for name, dimension in (benchmark_matrix.get("dimensions") or {}).items():
            lines.append(
                f"- `{name}` registered: `{len(dimension.get('registered') or [])}/{len(dimension.get('required') or [])}`, "
                f"executed: `{len(dimension.get('executed') or [])}/{len(dimension.get('required') or [])}`"
            )
        gaps = benchmark_matrix.get("executed_gaps") or []
        if gaps:
            lines.extend(["", "Executed benchmark matrix gaps:"])
            for gap in gaps[:12]:
                lines.append(f"- `{gap}`")
    traceability_contract = report.get("scenario_traceability_contract") or {}
    if traceability_contract:
        lines.extend(
            [
                "",
                "## Scenario Traceability Contract",
                "",
                f"- Status: `{traceability_contract.get('status', '')}`",
                f"- Checked cases: `{traceability_contract.get('checked_count', 0)}`",
                f"- Failed cases: `{traceability_contract.get('failed_case_count', 0)}`",
            ]
        )
        for case in (traceability_contract.get("failed_cases") or [])[:12]:
            codes = ", ".join(str(row.get("code") or "") for row in case.get("failures", []) or [])
            lines.append(f"- `{case.get('case_id', '')}`: {codes}")
    signals = ((report.get("optimization_signals") or {}).get("by_code") or [])
    if signals:
        lines.extend(
            [
                "",
                "## Optimization Queue",
                "",
                "| signal | severity | targets | cases | recommended action |",
                "| --- | --- | --- | ---: | --- |",
            ]
        )
        for signal in signals:
            cases = ", ".join(str(case) for case in signal.get("cases", [])[:6])
            if len(signal.get("cases", []) or []) > 6:
                cases += ", ..."
            lines.append(
                "| {code} | {severity} | {targets} | {count} | {action} |".format(
                    code=signal.get("code", ""),
                    severity=signal.get("severity", ""),
                    targets=", ".join(str(target) for target in signal.get("optimization_targets", [])[:6]),
                    count=signal.get("case_count", 0),
                    action=f"{signal.get('recommended_action', '')} Cases: {cases}",
                )
            )
        target_rows = ((report.get("optimization_signals") or {}).get("by_target") or [])
        if target_rows:
            lines.extend(
                [
                    "",
                    "Optimization targets:",
                    "",
                    "| target | severity | signals | cases |",
                    "| --- | --- | --- | ---: |",
                ]
            )
            for row in target_rows[:20]:
                lines.append(
                    "| {target} | {severity} | {codes} | {cases} |".format(
                        target=row.get("bucket", ""),
                        severity=row.get("max_severity", ""),
                        codes=", ".join(str(code) for code in row.get("codes", [])[:6]),
                        cases=len(row.get("cases") or []),
                    )
                )
        dimension_rows = ((report.get("optimization_signals") or {}).get("by_dimension") or [])
        if dimension_rows:
            lines.extend(
                [
                    "",
                    "Optimization dimensions:",
                    "",
                    "| dimension | severity | targets | signals | cases |",
                    "| --- | --- | --- | --- | ---: |",
                ]
            )
            for row in dimension_rows[:20]:
                lines.append(
                    "| {dimension} | {severity} | {targets} | {codes} | {cases} |".format(
                        dimension=row.get("bucket", ""),
                        severity=row.get("max_severity", ""),
                        targets=", ".join(str(target) for target in row.get("optimization_targets", [])[:6]),
                        codes=", ".join(str(code) for code in row.get("codes", [])[:6]),
                        cases=len(row.get("cases") or []),
                    )
                )
    gold_coverage = report.get("gold_coverage") or {}
    if gold_coverage:
        lines.extend(
            [
                "",
                "## Gold Coverage",
                "",
                f"- Status: `{gold_coverage.get('status', '')}`",
                f"- Covered cells: `{gold_coverage.get('covered_cell_count', 0)}/{gold_coverage.get('required_cell_count', 0)}`",
                f"- Coverage fraction: `{gold_coverage.get('coverage_fraction')}`",
                f"- Negative controls: `{(gold_coverage.get('negative_control_requirement') or {}).get('observed')}/{(gold_coverage.get('negative_control_requirement') or {}).get('minimum')}`",
            ]
        )
        missing = gold_coverage.get("missing_cells") or []
        if missing:
            lines.extend(["", "Missing gold coverage cells:"])
            for row in missing[:12]:
                lines.append(
                    f"- `{row.get('cancer_id')}` / `{row.get('theme_id')}`: "
                    f"{row.get('matched_positive_count', 0)}/{row.get('min_positive_controls', 1)} positive controls"
                )
    gold_quality = report.get("gold_quality") or {}
    if gold_quality:
        lines.extend(
            [
                "",
                "## Gold Quality",
                "",
                f"- Status: `{gold_quality.get('status', '')}`",
                f"- Positive conclusions: `{gold_quality.get('positive_count', 0)}/{gold_quality.get('minimum_positive_conclusions', 0)}`",
                f"- Negative controls: `{gold_quality.get('negative_control_count', 0)}/{gold_quality.get('minimum_negative_controls', 0)}`",
                f"- Negative trap categories covered: `{len(gold_quality.get('negative_trap_categories') or []) - len(gold_quality.get('missing_negative_trap_categories') or [])}/{len(gold_quality.get('negative_trap_categories') or [])}`",
            ]
        )
        failed = gold_quality.get("failed_checks") or []
        if failed:
            lines.extend(["", "Failed gold quality checks:"])
            for row in failed[:12]:
                lines.append(f"- `{row.get('code')}`")
    literature_semantics = report.get("literature_semantics") or {}
    if literature_semantics:
        metrics = literature_semantics.get("metrics") or {}
        lines.extend(
            [
                "",
                "## Literature Semantics",
                "",
                f"- Status: `{literature_semantics.get('status', '')}`",
                f"- Cases: `{metrics.get('passed_case_count', 0)}/{metrics.get('case_count', 0)}`",
                f"- Status accuracy: `{metrics.get('status_accuracy')}`",
                f"- Cue recall: `{metrics.get('cue_recall')}`",
                f"- Method cue recall: `{metrics.get('method_cue_recall')}`",
            ]
        )
        failed = literature_semantics.get("failed_cases") or []
        if failed:
            lines.extend(["", "Failed literature semantics cases:"])
            for row in failed[:12]:
                lines.append(
                    f"- `{row.get('case_id')}`: expected `{row.get('expected_support_status')}`, "
                    f"observed `{row.get('observed_support_status')}`"
                )
    lines.extend(
        [
            "",
            "## Held-Out Validation",
            "",
            f"- Temporal available: `{((report.get('heldout') or {}).get('temporal_holdout') or {}).get('available')}`",
            f"- Source available: `{((report.get('heldout') or {}).get('source_heldout') or {}).get('available')}`",
            f"- Release gate: `{((report.get('heldout') or {}).get('release_gate') or {}).get('gate_status', '')}`",
        ]
    )
    heldout_validation = report.get("heldout_validation") or {}
    if heldout_validation:
        lines.extend(
            [
                "",
                "## Held-Out Audit",
                "",
                f"- Status: `{heldout_validation.get('status', '')}`",
                f"- Rerun performed: `{heldout_validation.get('rerun_performed')}`",
                f"- Require rerun: `{heldout_validation.get('require_rerun')}`",
                f"- Max age hours: `{heldout_validation.get('max_age_hours')}`",
            ]
        )
        failed = heldout_validation.get("failed_checks") or []
        if failed:
            lines.extend(["", "Failed held-out checks:"])
            for row in failed[:12]:
                lines.append(f"- `{row.get('code')}`: {row.get('message', '')}")
        warnings = heldout_validation.get("warnings") or []
        if warnings:
            lines.extend(["", "Held-out warnings:"])
            for row in warnings[:12]:
                lines.append(f"- `{row.get('code')}`: {row.get('message', '')}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_benchmark(args: argparse.Namespace) -> dict[str, Any]:
    args = apply_release_profile(args)
    workspace = Path(args.workspace).resolve()
    release_id = args.release_id or latest_release_id(workspace / "graph_projection")
    output_dir = workspace / args.output_root / release_id / "batch_benchmark"
    progress_log = Path(args.progress_log).resolve() if args.progress_log else output_dir / "release_batch_benchmark_progress.jsonl"
    if progress_log.exists():
        progress_log.write_text("", encoding="utf-8")
    progress_events: list[dict[str, Any]] = []
    progress_event(
        progress_events,
        "batch_started",
        progress_log=progress_log,
        emit=not args.quiet_progress,
        release_id=release_id,
        scenario_timeout_seconds=args.scenario_timeout_seconds,
    )
    fixture_registry = read_json(workspace / args.fixture_registry)
    registry_scenarios = load_registry_scenarios(workspace, fixture_registry)
    args.required_fixture_case = dedupe_preserve(list(args.required_fixture_case or []) + registry_required_cases(fixture_registry))
    scenarios = default_scenarios()
    if registry_scenarios:
        scenarios.extend(registry_scenarios)
    else:
        scenarios.append(load_cssc_scenario((workspace / args.cscc_csv).resolve(), args.cscc_top_n))

    results: list[dict[str, Any]] = []
    scenario_run_dir = output_dir / "scenario_runs"
    for index, case in enumerate(scenarios, start=1):
        case_id = str(case.get("case_id") or f"scenario_{index}")
        if case.get("status") != "ready":
            progress_event(
                progress_events,
                "scenario_skipped",
                case_id=case_id,
                progress_log=progress_log,
                emit=not args.quiet_progress,
                status=case.get("status", ""),
            )
            results.append(case)
            continue
        progress_event(
            progress_events,
            "scenario_started",
            case_id=case_id,
            progress_log=progress_log,
            emit=not args.quiet_progress,
            index=index,
            total=len(scenarios),
            input_count=len(case.get("records", []) or []),
        )
        result = run_ready_scenario_isolated(
            case,
            index=index,
            workspace=workspace,
            release_id=release_id,
            max_paths=args.max_paths,
            max_hops=args.max_hops,
            timeout_seconds=args.scenario_timeout_seconds,
            run_dir=scenario_run_dir,
        )
        results.append(result)
        progress_event(
            progress_events,
            "scenario_completed",
            case_id=case_id,
            progress_log=progress_log,
            emit=not args.quiet_progress,
            status=result.get("status", ""),
            duration_seconds=result.get("duration_seconds"),
            timed_out=bool(result.get("timed_out")),
            error_type=result.get("error_type", ""),
        )

    fixture_audit = audit_fixture_availability(results, args.required_fixture_case, args.require_real_fixtures)
    progress_event(
        progress_events,
        "heldout_rerun_started" if args.rerun_heldout else "heldout_rerun_skipped",
        progress_log=progress_log,
        emit=not args.quiet_progress,
        timeout_seconds=args.heldout_rerun_timeout_seconds,
    )
    rerun = (
        rerun_heldout(workspace, args.run_id, args.heldout_top_n, args.heldout_rerun_timeout_seconds)
        if args.rerun_heldout
        else None
    )
    if args.rerun_heldout:
        progress_event(
            progress_events,
            "heldout_rerun_completed",
            progress_log=progress_log,
            emit=not args.quiet_progress,
            temporal_passed=bool(((rerun or {}).get("temporal") or {}).get("passed")),
            source_passed=bool(((rerun or {}).get("source") or {}).get("passed")),
        )
    heldout = summarize_existing_heldout(workspace, args.run_id)
    heldout_validation = audit_heldout_validation(
        heldout,
        rerun,
        args.max_heldout_age_hours,
        args.require_heldout_rerun,
        args.min_temporal_heldout_rows,
        args.min_source_heldout_sources,
        args.min_source_heldout_rows,
    )
    gold_coverage = audit_coverage(read_json(workspace / args.gold_path), read_json(workspace / args.gold_coverage_requirements))
    gold_quality = audit_quality(read_json(workspace / args.gold_path), read_json(workspace / args.gold_quality_requirements))
    literature_semantics = benchmark_literature_semantics(read_json(workspace / args.literature_semantics_benchmark))
    benchmark_matrix = audit_benchmark_matrix(
        results,
        list(args.required_benchmark_cancer or []),
        list(args.required_benchmark_theme or []),
        list(args.required_benchmark_input_mode or []),
        bool(args.require_real_fixtures),
    )
    scenario_traceability_contract = audit_scenario_traceability(results)
    scenario_signals = aggregate_optimization_signals(results)
    report = {
        "report_version": RUNNER_VERSION,
        "created_at_utc": utc_now(),
        "workspace": str(workspace),
        "release_id": release_id,
        "run_id": args.run_id,
        "release_profile": release_profile_summary(args),
        "progress_log": str(progress_log),
        "progress_events": progress_events,
        "fixture_registry": {
            "path": str((workspace / args.fixture_registry).resolve()),
            "schema_version": fixture_registry.get("schema_version", ""),
            "fixture_count": len(fixture_registry.get("fixtures", []) or []),
        },
        "scenarios": results,
        "fixture_audit": fixture_audit,
        "heldout": heldout,
        "heldout_rerun": rerun,
        "heldout_validation": heldout_validation,
        "gold_coverage": gold_coverage,
        "gold_quality": gold_quality,
        "literature_semantics": literature_semantics,
        "benchmark_matrix": benchmark_matrix,
        "scenario_traceability_contract": scenario_traceability_contract,
        "optimization_signals": merge_report_signals(
            scenario_signals,
            [
                gold_coverage_signal(gold_coverage),
                gold_quality_signal(gold_quality),
                literature_semantics_signal(literature_semantics),
                benchmark_matrix_signal(benchmark_matrix),
                traceability_contract_signal(scenario_traceability_contract),
                heldout_validation_signal(heldout_validation),
                fixture_audit_signal(fixture_audit),
            ],
        ),
        "status": benchmark_status(
            results,
            heldout,
            rerun,
            literature_semantics,
            heldout_validation,
            fixture_audit,
            benchmark_matrix,
            scenario_traceability_contract,
        ),
    }
    json_path = output_dir / "release_batch_benchmark_report.json"
    md_path = output_dir / "release_batch_benchmark_report.md"
    write_json(json_path, report)
    write_markdown(md_path, report)
    report["outputs"] = {"json": str(json_path), "markdown": str(md_path)}
    write_json(json_path, report)
    write_markdown(md_path, report)
    return report


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run release-level benchmark scenarios and held-out report checks.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", default="")
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--cscc-csv", default=DEFAULT_CSSC_CSV)
    parser.add_argument("--cscc-top-n", type=int, default=80)
    parser.add_argument("--fixture-registry", default=DEFAULT_REAL_WORLD_FIXTURE_REGISTRY)
    parser.add_argument("--strict-release", action="store_true")
    parser.add_argument("--require-real-fixtures", action="store_true")
    parser.add_argument("--required-fixture-case", action="append", default=list(DEFAULT_REQUIRED_REAL_FIXTURE_CASES))
    parser.add_argument("--max-paths", type=int, default=10)
    parser.add_argument("--max-hops", type=int, default=3)
    parser.add_argument("--rerun-heldout", action="store_true")
    parser.add_argument("--heldout-top-n", type=int, default=200)
    parser.add_argument("--heldout-rerun-timeout-seconds", type=int, default=DEFAULT_HELDOUT_RERUN_TIMEOUT_SECONDS)
    parser.add_argument("--scenario-timeout-seconds", type=int, default=DEFAULT_SCENARIO_TIMEOUT_SECONDS)
    parser.add_argument("--progress-log", default="")
    parser.add_argument("--quiet-progress", action="store_true")
    parser.add_argument("--single-scenario-input", default=argparse.SUPPRESS)
    parser.add_argument("--single-scenario-output", default=argparse.SUPPRESS)
    parser.add_argument("--require-heldout-rerun", action="store_true")
    parser.add_argument("--max-heldout-age-hours", type=int, default=DEFAULT_MAX_HELDOUT_AGE_HOURS)
    parser.add_argument("--min-temporal-heldout-rows", type=int, default=DEFAULT_MIN_TEMPORAL_HELDOUT_ROWS)
    parser.add_argument("--min-source-heldout-sources", type=int, default=DEFAULT_MIN_SOURCE_HELDOUT_SOURCES)
    parser.add_argument("--min-source-heldout-rows", type=int, default=DEFAULT_MIN_SOURCE_HELDOUT_ROWS)
    parser.add_argument("--gold-path", default=DEFAULT_GOLD_PATH)
    parser.add_argument("--gold-coverage-requirements", default=DEFAULT_GOLD_COVERAGE_REQUIREMENTS)
    parser.add_argument("--gold-quality-requirements", default=DEFAULT_GOLD_QUALITY_REQUIREMENTS)
    parser.add_argument("--literature-semantics-benchmark", default=DEFAULT_LITERATURE_SEMANTICS_BENCHMARK)
    parser.add_argument("--required-benchmark-cancer", action="append", default=list(DEFAULT_REQUIRED_BENCHMARK_CANCERS))
    parser.add_argument("--required-benchmark-theme", action="append", default=list(DEFAULT_REQUIRED_BENCHMARK_THEMES))
    parser.add_argument("--required-benchmark-input-mode", action="append", default=list(DEFAULT_REQUIRED_BENCHMARK_INPUT_MODES))
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if hasattr(args, "single_scenario_input"):
        return run_scenario_child_main(args)
    report = run_benchmark(args)
    print(
        json.dumps(
            {
                "passed": report["status"]["passed"],
                "report": report["outputs"]["json"],
                "markdown": report["outputs"]["markdown"],
                "failed_case_count": report["status"]["failed_case_count"],
                "missing_case_count": report["status"]["missing_case_count"],
                "release_profile": (report.get("release_profile") or {}).get("profile", ""),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0 if report["status"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

