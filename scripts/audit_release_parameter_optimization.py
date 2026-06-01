"""Build a release-level parameter and algorithm optimization queue.

This audit turns benchmark/readiness/objective evidence into concrete schema,
parameter, or algorithm follow-up actions. It is deliberately conservative:
actions must improve evidence structure or calibration and must not lower
confidence thresholds to hide errors.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AUDIT_VERSION = "release_parameter_optimization.v1"
DEFAULT_BATCH_REPORT = "validation_reports/mvp_20260513T002254/batch_benchmark/release_batch_benchmark_report.json"
DEFAULT_READINESS_REPORT = "validation_reports/mvp_20260513T002254/batch_benchmark/release_readiness_audit.json"
DEFAULT_OBJECTIVE_REPORT = "validation_reports/mvp_20260513T002254/optimization_objective_status.json"
DEFAULT_METRIC_SNAPSHOT = "validation_reports/mvp_20260513T002254/release_metric_snapshot.json"

MANDATORY_VALIDATION_SUITES = [
    "gold_set",
    "negative_traps",
    "source_heldout",
    "temporal_heldout",
    "release_batch_benchmark",
    "scenario_traceability_contract",
    "cscc_real_fixture_regression",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def finding_codes(rows: list[dict[str, Any]]) -> set[str]:
    return {str(row.get("code") or "") for row in rows if str(row.get("code") or "")}


def add_action(
    actions: dict[str, dict[str, Any]],
    action_id: str,
    category: str,
    trigger: str,
    parameter_targets: list[str],
    required: bool,
    rationale: str,
    evidence_refs: list[str] | None = None,
    validation_suites: list[str] | None = None,
    details: dict[str, Any] | None = None,
) -> None:
    row = actions.setdefault(
        action_id,
        {
            "action_id": action_id,
            "category": category,
            "trigger": trigger,
            "parameter_targets": [],
            "required": required,
            "rationale": rationale,
            "evidence_refs": [],
            "validation_suites": [],
        },
    )
    row["required"] = bool(row.get("required") or required)
    row["parameter_targets"] = sorted(set(list(row.get("parameter_targets") or []) + parameter_targets))
    row["evidence_refs"] = sorted(set(list(row.get("evidence_refs") or []) + list(evidence_refs or [])))
    suites = list(validation_suites or MANDATORY_VALIDATION_SUITES)
    row["validation_suites"] = sorted(set(list(row.get("validation_suites") or []) + suites))
    if details:
        for key, value in details.items():
            if isinstance(value, list):
                existing = row.setdefault(key, [])
                if isinstance(existing, list):
                    existing.extend(value)
            elif isinstance(value, dict):
                existing = row.setdefault(key, {})
                if isinstance(existing, dict):
                    existing.update(value)
            else:
                row[key] = value


def signal_refs(signal: dict[str, Any]) -> list[str]:
    return [str(case) for case in signal.get("cases", []) or [] if str(case)]


def signal_matrix_summary(optimization_signals: dict[str, Any]) -> dict[str, Any]:
    by_target = [row for row in optimization_signals.get("by_target", []) or [] if isinstance(row, dict)]
    by_dimension = [row for row in optimization_signals.get("by_dimension", []) or [] if isinstance(row, dict)]

    def normalize_row(row: dict[str, Any]) -> dict[str, Any]:
        return {
            "bucket": str(row.get("bucket") or ""),
            "max_severity": str(row.get("max_severity") or ""),
            "signal_count": int(row.get("signal_count") or 0),
            "case_count": len(row.get("cases") or []),
            "codes": [str(code) for code in row.get("codes", []) or [] if str(code)],
            "cases": [str(case) for case in row.get("cases", []) or [] if str(case)][:12],
            "optimization_targets": [str(target) for target in row.get("optimization_targets", []) or [] if str(target)],
        }

    def sort_key(row: dict[str, Any]) -> tuple[int, int, str]:
        severity_rank = 0 if row.get("max_severity") == "error" else 1
        return (severity_rank, -int(row.get("signal_count") or 0), str(row.get("bucket") or ""))

    targets = sorted((normalize_row(row) for row in by_target), key=sort_key)
    dimensions = sorted((normalize_row(row) for row in by_dimension), key=sort_key)
    return {
        "status": "present" if targets or dimensions else "absent",
        "target_count": len(targets),
        "dimension_count": len(dimensions),
        "error_targets": [row for row in targets if row.get("max_severity") == "error"],
        "warning_targets": [row for row in targets if row.get("max_severity") == "warning"],
        "error_dimensions": [row for row in dimensions if row.get("max_severity") == "error"],
        "warning_dimensions": [row for row in dimensions if row.get("max_severity") == "warning"],
        "top_targets": targets[:10],
        "top_dimensions": dimensions[:10],
        "interpretation": (
            "Optimization targets and dimensions are derived from release benchmark signals; they indicate which schema, "
            "resolver, parser, context, graph, ranking, or fixture surfaces should be changed before rerunning gates."
        ),
    }


def fixture_recovery_details(missing_required: list[Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in missing_required:
        if not isinstance(row, dict):
            continue
        rows.append(
            {
                "case_id": str(row.get("case_id") or ""),
                "path": str(row.get("path") or ""),
                "expected_schema": row.get("expected_schema") or {},
                "candidate_paths": row.get("candidate_paths") or [],
                "recovery_action": str(row.get("recovery_action") or ""),
                "recovery_rerun_plan": row.get("recovery_rerun_plan") or {},
            }
        )
    return rows


def source_heldout_summary_by_label(batch_report: dict[str, Any]) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    summary = ((batch_report.get("heldout") or {}).get("source_heldout") or {}).get("summary") or []
    for row in summary:
        if not isinstance(row, dict):
            continue
        label = str(row.get("source") or "")
        if label:
            rows[label] = row
    return rows


def heldout_warning_requires_action(
    warning: dict[str, Any],
    source_summary_by_label: dict[str, dict[str, Any]],
) -> tuple[bool, dict[str, Any]]:
    code = str(warning.get("code") or "")
    label = str((warning.get("metrics") or {}).get("label") or "")
    diagnostic = {
        "code": code,
        "label": label,
        "requires_action": True,
        "reason": "unclassified_heldout_warning",
    }
    if code != "source_heldout_random_auc":
        return True, diagnostic

    summary_row = source_summary_by_label.get(label, {})
    validation_metrics = summary_row.get("validation_metrics") or {}
    warning_metrics = warning.get("metrics") or {}
    diagnostics = warning_metrics.get("diagnostics") or {}
    feature_transfer = diagnostics.get("feature_transfer") or summary_row.get("feature_transfer_diagnostics") or {}
    source_family = diagnostics.get("source_family") or summary_row.get("source_family_diagnostics") or {}
    lift_at_100 = safe_float(validation_metrics.get("heldout_lift_at_100"), 0.0)
    unique_scores = safe_float(validation_metrics.get("heldout_unique_transfer_scores"), 0.0)
    tie_fraction = safe_float(validation_metrics.get("heldout_top_score_tie_fraction"), 1.0)
    shared_features = safe_float(feature_transfer.get("shared_nonconstant_feature_count"), 0.0)
    same_family_available = bool(source_family.get("same_family_available"))
    topk_transfer_passed = lift_at_100 > 1.0 and unique_scores >= 2.0 and tie_fraction < 0.99
    transfer_surface_available = same_family_available and shared_features > 0.0
    diagnostic.update(
        {
            "heldout_lift_at_100": lift_at_100,
            "heldout_unique_transfer_scores": unique_scores,
            "heldout_top_score_tie_fraction": tie_fraction,
            "same_family_available": same_family_available,
            "shared_nonconstant_feature_count": shared_features,
            "topk_transfer_passed": topk_transfer_passed,
            "transfer_surface_available": transfer_surface_available,
        }
    )
    if topk_transfer_passed and transfer_surface_available:
        diagnostic["requires_action"] = False
        diagnostic["reason"] = "global_auc_low_but_topk_transfer_and_feature_surface_passed"
        return False, diagnostic
    diagnostic["reason"] = "global_auc_low_without_sufficient_topk_or_transfer_surface"
    return True, diagnostic


def audit_parameter_optimization(
    batch_report: dict[str, Any],
    readiness_report: dict[str, Any],
    objective_report: dict[str, Any],
    metric_snapshot: dict[str, Any],
) -> dict[str, Any]:
    actions: dict[str, dict[str, Any]] = {}
    diagnostics: dict[str, Any] = {}
    batch_status = batch_report.get("status") or {}
    matrix = batch_report.get("benchmark_matrix") or {}
    fixture = batch_report.get("fixture_audit") or {}
    heldout = batch_report.get("heldout_validation") or {}
    traceability = batch_report.get("scenario_traceability_contract") or {}
    literature = batch_report.get("literature_semantics") or {}
    gold_quality = batch_report.get("gold_quality") or {}
    gold_coverage = batch_report.get("gold_coverage") or {}
    gold_specificity = gold_coverage.get("coverage_specificity") or {}
    metrics = metric_snapshot.get("metrics") or {}
    readiness_blocker_codes = finding_codes(readiness_report.get("blockers", []) or [])
    objective_status = str(objective_report.get("status") or "")
    signal_matrix = signal_matrix_summary(batch_report.get("optimization_signals") or {})
    if signal_matrix.get("status") == "present":
        diagnostics["optimization_signal_matrix"] = signal_matrix

    executed_gaps = [str(gap) for gap in matrix.get("executed_gaps", []) or [] if str(gap)]
    if int(matrix.get("executed_gap_count") or 0) > 0:
        add_action(
            actions,
            "benchmark_matrix_execute_required_dimensions",
            "benchmark_runner",
            "required benchmark dimensions are registered but not executed",
            ["required_benchmark_cancers", "required_benchmark_themes", "real_world_fixture_registry"],
            True,
            "Release-level calibration cannot rely on registered-only coverage; missing executed dimensions must be restored or added.",
            executed_gaps,
            ["release_batch_benchmark", "scenario_traceability_contract", "cscc_real_fixture_regression"],
        )

    missing_required = fixture.get("missing_required") or []
    if fixture.get("status") == "failed" or "fixture_audit_failed" in readiness_blocker_codes:
        add_action(
            actions,
            "fixture_schema_restore_required_inputs",
            "database_schema",
            "required real-world fixture missing or invalid",
            ["real_world_fixture_registry", "fixture_identity_columns_any", "fixture_effect_columns_any"],
            True,
            "Required fixtures are part of the calibration surface; repair registry/schema rather than lowering evidence gates.",
            [str(row.get("case_id") or "") for row in missing_required if isinstance(row, dict)],
            ["release_batch_benchmark", "cscc_real_fixture_regression"],
            {"fixture_recovery": fixture_recovery_details(missing_required)},
        )

    if heldout.get("status") != "passed" or not heldout.get("rerun_performed") or "heldout_not_rerun" in readiness_blocker_codes:
        add_action(
            actions,
            "heldout_rerun_enforce_current_validation",
            "heldout_validation",
            "held-out validation missing, stale, failed, or not rerun in strict path",
            [
                "require_heldout_rerun",
                "max_heldout_age_hours",
                "min_temporal_heldout_rows",
                "min_source_heldout_sources",
                "heldout_transfer_quality_thresholds",
                "source_specific_feature_coverage",
                "source_family_holdout_strategy",
                "literature_source_partitioning",
                "tie_collapse_guard",
            ],
            True,
            "Parameter changes require current temporal/source held-out checks with measurable transfer signal before becoming defaults.",
            [str(row.get("code") or "") for row in heldout.get("failed_checks", []) or [] if isinstance(row, dict)],
            ["source_heldout", "temporal_heldout"],
        )
    source_summary = source_heldout_summary_by_label(batch_report)
    actionable_warnings: list[str] = []
    heldout_warning_diagnostics: list[dict[str, Any]] = []
    for row in heldout.get("warnings", []) or []:
        if not isinstance(row, dict) or not str(row.get("code") or ""):
            continue
        requires_action, diagnostic = heldout_warning_requires_action(row, source_summary)
        heldout_warning_diagnostics.append(diagnostic)
        if requires_action:
            actionable_warnings.append(str(row.get("code") or ""))
    if heldout_warning_diagnostics:
        diagnostics["heldout_warning_diagnostics"] = heldout_warning_diagnostics
    if actionable_warnings or (not heldout_warning_diagnostics and safe_float(metrics.get("heldout_warning_count")) > 0):
        add_action(
            actions,
            "heldout_topk_global_rank_calibration_review",
            "heldout_validation",
            "held-out warnings remain after strict pass",
            [
                "global_auc_vs_topk_policy",
                "source_partition_bucket_count",
                "topk_ranker_calibration",
                "literature_partition_score_monotonicity",
            ],
            False,
            "Passing top-k transfer with low global AUC should stay non-blocking but remain tracked for ranker calibration.",
            actionable_warnings,
            ["source_heldout", "temporal_heldout", "release_batch_benchmark"],
        )

    if traceability.get("status") not in {"", "passed"}:
        add_action(
            actions,
            "scenario_traceability_contract_repair",
            "output_contract",
            "executed scenarios failed to expose auditable structured evidence",
            ["analysis_pack_contract", "identity_decision_summary", "evidence_assertion_summary", "mechanism_ready_fact_count"],
            True,
            "Benchmark results must be traceable to inputs, identity decisions, mechanism facts, and evidence assertions.",
            [str(row.get("case_id") or "") for row in traceability.get("failed_cases", []) or [] if isinstance(row, dict)],
            ["scenario_traceability_contract", "release_batch_benchmark"],
        )

    if literature.get("status") not in {"", "passed"} or safe_float(metrics.get("literature_required_coverage_missing_count")) > 0:
        add_action(
            actions,
            "literature_parser_semantics_expand",
            "literature_parser",
            "literature semantics benchmark failed or cue coverage is incomplete",
            ["support_status_classifier", "semantic_cue_coverage", "method_cue_coverage", "background_penalty"],
            True,
            "Complex sentences must enter scoring only as structured assertions with support/contradict/uncertain/background status.",
            [str(row.get("case_id") or "") for row in literature.get("failed_cases", []) or [] if isinstance(row, dict)],
            ["release_batch_benchmark", "negative_traps"],
        )

    if gold_coverage.get("status") not in {"", "covered"} or gold_quality.get("status") not in {"", "passed"}:
        add_action(
            actions,
            "gold_standard_expand_or_repair",
            "gold_standard",
            "gold coverage or quality gate failed",
            ["gold_standard_conclusions", "negative_trap_categories", "coverage_requirements"],
            True,
            "Calibration needs broad, context-specific positives and negative traps rather than single-example tuning.",
            [str(row.get("code") or "") for row in gold_quality.get("failed_checks", []) or [] if isinstance(row, dict)],
            ["gold_set", "negative_traps", "release_batch_benchmark"],
        )

    if int(gold_specificity.get("over_recommended_count") or 0) > 0:
        add_action(
            actions,
            "gold_standard_split_overbroad_controls",
            "gold_standard",
            "gold standard controls are reused across too many cancer-theme cells",
            ["gold_standard_conclusions", "coverage_requirements", "mechanism_axis_specificity"],
            False,
            "Overbroad gold controls should be split into narrower disease-theme conclusions so benchmark coverage reflects precise calibration rather than broad wording.",
            [str(gold_id) for gold_id in gold_specificity.get("over_recommended_gold_ids", []) or [] if str(gold_id)],
            ["gold_set", "release_batch_benchmark", "scenario_traceability_contract"],
        )

    for signal in ((batch_report.get("optimization_signals") or {}).get("by_code") or []):
        code = str(signal.get("code") or "")
        refs = signal_refs(signal)
        if code in {"identity_ambiguity_review", "identity_unmatched_review"}:
            add_action(
                actions,
                f"resolver_{code}",
                "resolver",
                code,
                ["stable_identifier_priority", "name_ambiguity_abstention", "lipid_class_exact_match_guard", "mz_rt_ms2_orthogonal_support"],
                False,
                "Resolver calibration should improve identity evidence and abstention behavior, not promote weak matches.",
                refs,
            )
        elif code in {"false_positive_trap", "negative_trap_hit"}:
            add_action(
                actions,
                "negative_trap_blocking_recalibrate",
                "ranking",
                code,
                ["negative_trap_block_rule", "context_weight", "graph_propagation_penalty", "generic_pathway_penalty"],
                True,
                "False-positive traps must block release and trigger schema/parser/ranking inspection.",
                refs,
            )
        elif code == "unsupported_claim_risk":
            add_action(
                actions,
                "evidence_support_requirement_tighten",
                "evidence_scoring",
                code,
                ["top_claim_evidence_ref_requirement", "observation_appendix_policy", "evidence_assertion_linking"],
                False,
                "Top claims should require evidence refs or mechanism traces; observation-only claims stay appendix-only.",
                refs,
            )
        elif code in {"context_drift_review", "gold_context_mismatch"}:
            add_action(
                actions,
                "context_weight_recalibrate",
                "context_model",
                code,
                ["context_weight", "context_mismatch_cap", "tissue_cell_type_penalty", "pan_cancer_to_context_penalty"],
                False,
                "Context drift should be fixed by disease/tissue/cell-type calibration and propagation penalties.",
                refs,
            )
        elif code == "evidence_not_effective":
            add_action(
                actions,
                "evidence_assertion_linkage_repair",
                "evidence_scoring",
                code,
                ["evidence_usefulness_minimum", "assertion_to_fact_linker", "literature_weight"],
                False,
                "Literature evidence must have measurable effects through structured assertion links.",
                refs,
            )
        elif code == "false_negative_proxy":
            add_action(
                actions,
                "false_negative_error_classification",
                "benchmark_runner",
                code,
                ["entity_coverage_queue", "relation_parser_gap_queue", "gold_context_scope_review"],
                False,
                "Missing expected gold hits should be classified before changing ranking parameters.",
                refs,
            )
        elif code == "expected_topk_recall_gap":
            add_action(
                actions,
                "scenario_expected_topk_rank_calibration_review",
                "ranking_calibration",
                code,
                ["mechanism_fact_priority", "gold_entity_alias_map", "candidate_source_priority", "pathway_auxiliary_weight"],
                False,
                "Scenario-level expected gold should be recovered near the top-k ranking without promoting appendix or weak-evidence candidates.",
                refs,
                ["release_batch_benchmark", "scenario_traceability_contract", "gold_set"],
            )

    expected_eligible = float(metrics.get("scenario_expected_topk_eligible_count") or 0.0)
    expected_recall_at_5 = metrics.get("scenario_expected_recall_at_5_min")
    expected_precision_at_5 = metrics.get("scenario_expected_precision_at_5_min")
    if expected_eligible > 0 and expected_recall_at_5 is not None and float(expected_recall_at_5 or 0.0) < 1.0:
        add_action(
            actions,
            "scenario_expected_topk_rank_calibration_review",
            "ranking_calibration",
            "scenario expected gold top-k recall is below target",
            ["mechanism_fact_priority", "gold_entity_alias_map", "candidate_source_priority", "pathway_auxiliary_weight"],
            False,
            "Scenario-level expected gold should be recovered near the top-k ranking without promoting appendix or weak-evidence candidates.",
            [
                f"scenario_expected_recall_at_5_min={expected_recall_at_5}",
                f"scenario_expected_precision_at_5_min={expected_precision_at_5}",
            ],
            ["release_batch_benchmark", "scenario_traceability_contract", "gold_set"],
        )

    if objective_status and objective_status != "complete":
        add_action(
            actions,
            "objective_incomplete_release_block",
            "release_gate",
            "optimization objective audit is incomplete",
            ["strict_release_workflow", "objective_completion_gate"],
            True,
            "The default path must remain blocked until every objective requirement is proven by current reports.",
            [str(row.get("requirement_id") or "") for row in objective_report.get("remaining_blockers", []) or [] if isinstance(row, dict)],
            MANDATORY_VALIDATION_SUITES,
        )

    severity_counts = {
        "required": sum(1 for row in actions.values() if row.get("required")),
        "recommended": sum(1 for row in actions.values() if not row.get("required")),
    }
    ordered = sorted(actions.values(), key=lambda row: (not row.get("required"), row.get("category", ""), row.get("action_id", "")))
    return {
        "audit_version": AUDIT_VERSION,
        "created_at_utc": utc_now(),
        "status": "action_required" if ordered else "no_action_required",
        "action_count": len(ordered),
        "severity_counts": severity_counts,
        "actions": ordered,
        "diagnostics": diagnostics,
        "release_metric_context": {
            "batch_passed": metrics.get("batch_passed"),
            "benchmark_matrix_executed_gap_count": metrics.get("benchmark_matrix_executed_gap_count"),
            "heldout_rerun_performed": metrics.get("heldout_rerun_performed"),
            "objective_incomplete_count": metrics.get("objective_incomplete_count"),
            "scenario_traceability_failed_case_count": metrics.get("scenario_traceability_failed_case_count"),
            "scenario_negative_trap_hits": metrics.get("scenario_negative_trap_hits"),
            "scenario_unsupported_top_claim_rate_max": metrics.get("scenario_unsupported_top_claim_rate_max"),
            "scenario_precision_at_1_min": metrics.get("scenario_precision_at_1_min"),
            "scenario_precision_at_3_min": metrics.get("scenario_precision_at_3_min"),
            "scenario_precision_at_5_min": metrics.get("scenario_precision_at_5_min"),
            "scenario_recall_at_1_min": metrics.get("scenario_recall_at_1_min"),
            "scenario_recall_at_3_min": metrics.get("scenario_recall_at_3_min"),
            "scenario_recall_at_5_min": metrics.get("scenario_recall_at_5_min"),
            "scenario_expected_topk_eligible_count": metrics.get("scenario_expected_topk_eligible_count"),
            "scenario_expected_precision_at_1_min": metrics.get("scenario_expected_precision_at_1_min"),
            "scenario_expected_precision_at_3_min": metrics.get("scenario_expected_precision_at_3_min"),
            "scenario_expected_precision_at_5_min": metrics.get("scenario_expected_precision_at_5_min"),
            "scenario_expected_recall_at_1_min": metrics.get("scenario_expected_recall_at_1_min"),
            "scenario_expected_recall_at_3_min": metrics.get("scenario_expected_recall_at_3_min"),
            "scenario_expected_recall_at_5_min": metrics.get("scenario_expected_recall_at_5_min"),
            "gold_specificity_over_recommended_count": metrics.get("gold_specificity_over_recommended_count"),
            "gold_specificity_single_cell_fraction": metrics.get("gold_specificity_single_cell_fraction"),
        },
        "policy": {
            "confidence_policy": "Do not lower confidence thresholds to pass gates; optimize schema, parser, resolver, context weights, graph propagation penalties, or benchmark coverage.",
            "validation_required_before_default": MANDATORY_VALIDATION_SUITES,
        },
    }


def markdown_summary(report: dict[str, Any]) -> str:
    lines = [
        "# Release Parameter Optimization Audit",
        "",
        f"Status: `{report.get('status', '')}`",
        f"Actions: `{report.get('action_count', 0)}`",
        "",
        "## Policy",
        "",
        f"- {((report.get('policy') or {}).get('confidence_policy') or '')}",
        "",
        "## Actions",
        "",
    ]
    for row in report.get("actions", []) or []:
        required = "required" if row.get("required") else "recommended"
        targets = ", ".join(row.get("parameter_targets") or [])
        refs = ", ".join(row.get("evidence_refs") or [])
        lines.append(f"- `{required}` `{row.get('action_id')}` ({row.get('category')}): {row.get('rationale')}")
        if targets:
            lines.append(f"  - Targets: `{targets}`")
        if refs:
            lines.append(f"  - Evidence: `{refs}`")
        for recovery in row.get("fixture_recovery") or []:
            schema = recovery.get("expected_schema") or {}
            identity_cols = ", ".join(schema.get("identity_columns_any") or [])
            effect_cols = ", ".join(schema.get("effect_columns_any") or [])
            lines.append(
                f"  - Fixture `{recovery.get('case_id', '')}` expected mode `{schema.get('expected_input_mode', '')}`, "
                f"identity columns `{identity_cols}`, effect columns `{effect_cols}`, "
                f"candidate paths `{len(recovery.get('candidate_paths') or [])}`."
            )
            rerun_steps = (recovery.get("recovery_rerun_plan") or {}).get("validation_sequence") or []
            if rerun_steps:
                lines.append(f"  - Recovery rerun steps: `{len(rerun_steps)}`.")
    matrix = ((report.get("diagnostics") or {}).get("optimization_signal_matrix") or {})
    if matrix:
        lines.extend(["", "## Optimization Signal Matrix", ""])
        top_targets = matrix.get("top_targets") or []
        if top_targets:
            lines.extend(["| target | severity | signals | cases |", "| --- | --- | --- | ---: |"])
            for row in top_targets[:12]:
                lines.append(
                    "| {target} | {severity} | {codes} | {cases} |".format(
                        target=row.get("bucket", ""),
                        severity=row.get("max_severity", ""),
                        codes=", ".join(row.get("codes") or []),
                        cases=row.get("case_count", 0),
                    )
                )
        top_dimensions = matrix.get("top_dimensions") or []
        if top_dimensions:
            lines.extend(["", "| dimension | severity | targets | signals | cases |", "| --- | --- | --- | --- | ---: |"])
            for row in top_dimensions[:12]:
                lines.append(
                    "| {dimension} | {severity} | {targets} | {codes} | {cases} |".format(
                        dimension=row.get("bucket", ""),
                        severity=row.get("max_severity", ""),
                        targets=", ".join(row.get("optimization_targets") or []),
                        codes=", ".join(row.get("codes") or []),
                        cases=row.get("case_count", 0),
                    )
                )
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit release metrics into parameter/schema/algorithm optimization actions.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--batch-report", default=DEFAULT_BATCH_REPORT)
    parser.add_argument("--readiness-report", default=DEFAULT_READINESS_REPORT)
    parser.add_argument("--objective-report", default=DEFAULT_OBJECTIVE_REPORT)
    parser.add_argument("--metric-snapshot", default=DEFAULT_METRIC_SNAPSHOT)
    parser.add_argument("--output", default="")
    parser.add_argument("--markdown-output", default="")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    workspace = Path(args.workspace).resolve()
    report = audit_parameter_optimization(
        read_json(workspace / args.batch_report),
        read_json(workspace / args.readiness_report),
        read_json(workspace / args.objective_report),
        read_json(workspace / args.metric_snapshot),
    )
    if args.output:
        write_json(workspace / args.output, report)
    if args.markdown_output:
        path = workspace / args.markdown_output
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(markdown_summary(report), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if report.get("status") == "no_action_required" else 1


if __name__ == "__main__":
    raise SystemExit(main())
