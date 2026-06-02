"""Audit completion status for the tumor metabolism optimization objective.

This is a meta-audit over the release evidence. It does not rerun analysis; it
checks whether the current reports prove each explicit optimization requirement
or whether evidence is still missing/incomplete.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AUDIT_VERSION = "optimization_objective_status.v1"
DEFAULT_RELEASE_ID = "mvp_20260513T002254"
DEFAULT_BATCH_REPORT = "validation_reports/mvp_20260513T002254/batch_benchmark/release_batch_benchmark_report.json"
DEFAULT_READINESS_REPORT = "validation_reports/mvp_20260513T002254/batch_benchmark/release_readiness_audit.json"
DEFAULT_WORKFLOW_REPORT = "validation_reports/mvp_20260513T002254/strict_release_workflow/strict_release_workflow_report.json"
DEFAULT_PARAMETER_REPORT = "validation_reports/mvp_20260513T002254/release_parameter_optimization_audit.json"
DEFAULT_CSSC_CSV = "required_metabolism_fixture.csv"
NON_BLOCKING_PARAMETER_ACTION_IDS = {"objective_incomplete_release_block"}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def gate(status: bool, satisfied: str = "satisfied", incomplete: str = "incomplete") -> str:
    return satisfied if status else incomplete


def requirement_row(
    requirement_id: str,
    status: str,
    evidence: dict[str, Any],
    blockers: list[str] | None = None,
    interpretation: str = "",
) -> dict[str, Any]:
    return {
        "requirement_id": requirement_id,
        "status": status,
        "evidence": evidence,
        "blockers": blockers or [],
        "interpretation": interpretation,
    }


def has_executed_matrix_value(batch_report: dict[str, Any], dimension: str, value: str) -> bool:
    matrix = batch_report.get("benchmark_matrix") or {}
    dimensions = matrix.get("dimensions") or {}
    return value in set((dimensions.get(dimension) or {}).get("executed") or [])


def finding_codes(rows: list[dict[str, Any]]) -> list[str]:
    return [str(row.get("code") or "") for row in rows if str(row.get("code") or "")]


def remediation_refs(rows: list[dict[str, Any]]) -> list[str]:
    refs: list[str] = []
    for row in rows:
        refs.extend(str(ref) for ref in row.get("evidence_refs", []) or [] if str(ref))
    return sorted(set(refs))


def blocking_parameter_actions(parameter_report: dict[str, Any]) -> list[dict[str, Any]]:
    actions = parameter_report.get("actions", []) or []
    required_count = int((parameter_report.get("severity_counts") or {}).get("required") or 0)
    has_explicit_required_flags = any(isinstance(row, dict) and "required" in row for row in actions)
    return [
        row
        for row in actions
        if isinstance(row, dict)
        and (row.get("required") or (required_count > 0 and not has_explicit_required_flags))
        and str(row.get("action_id") or "") not in NON_BLOCKING_PARAMETER_ACTION_IDS
    ]


def cscc_fixture_recovery(batch_report: dict[str, Any]) -> dict[str, Any]:
    missing_required = ((batch_report.get("fixture_audit") or {}).get("missing_required") or [])
    for row in missing_required:
        if isinstance(row, dict) and str(row.get("case_id") or "") == "cscc_trait_score_csv":
            return {
                "case_id": str(row.get("case_id") or ""),
                "path": str(row.get("path") or ""),
                "expected_schema": row.get("expected_schema") or {},
                "candidate_paths": row.get("candidate_paths") or [],
                "recovery_action": str(row.get("recovery_action") or ""),
                "recovery_rerun_plan": row.get("recovery_rerun_plan") or {},
            }
    return {}


def audit_objective_status(
    workspace: Path,
    batch_report: dict[str, Any],
    readiness_report: dict[str, Any],
    workflow_report: dict[str, Any],
    cscc_csv: Path,
    parameter_report: dict[str, Any] | None = None,
) -> dict[str, Any]:
    gold_coverage = batch_report.get("gold_coverage") or {}
    gold_quality = batch_report.get("gold_quality") or {}
    literature = batch_report.get("literature_semantics") or {}
    heldout = batch_report.get("heldout_validation") or {}
    matrix = batch_report.get("benchmark_matrix") or {}
    traceability = batch_report.get("scenario_traceability_contract") or {}
    batch_status = batch_report.get("status") or {}
    readiness_status = readiness_report.get("status", "")
    workflow_status = workflow_report.get("status", "")
    parameter_report = parameter_report or {}
    parameter_status = str(parameter_report.get("status") or "")
    blocking_actions = blocking_parameter_actions(parameter_report)
    blocking_action_ids = [str(row.get("action_id") or "") for row in blocking_actions]

    gold_ok = gold_coverage.get("status") == "covered" and gold_quality.get("status") == "passed"
    gold_distribution = gold_quality.get("positive_distribution") or {}
    literature_metrics = literature.get("metrics") or {}
    literature_ok = (
        literature.get("status") == "passed"
        and int(literature_metrics.get("required_coverage_missing_count") or 0) == 0
        and int(literature_metrics.get("case_count") or 0) >= int(literature_metrics.get("minimum_case_count") or 1)
    )
    heldout_ok = heldout.get("status") == "passed" and bool(heldout.get("rerun_performed"))
    registered_matrix_ok = matrix.get("registered_gap_count") == 0
    executed_matrix_ok = matrix.get("status") == "passed" and int(matrix.get("executed_gap_count") or 0) == 0
    traceability_ok = traceability.get("status") == "passed"
    cscc_exists = cscc_csv.exists()
    cscc_executed = has_executed_matrix_value(batch_report, "cancers", "cscc")
    arachidonate_executed = has_executed_matrix_value(batch_report, "themes", "arachidonate")
    cscc_recovery = cscc_fixture_recovery(batch_report)
    cscc_ok = cscc_exists and cscc_executed
    batch_ok = bool(batch_status.get("passed")) and registered_matrix_ok and executed_matrix_ok and traceability_ok
    parameter_ok = not blocking_actions
    readiness_ok = readiness_status == "ready" and workflow_status == "ready"

    requirements = [
        requirement_row(
            "gold_standard_full_library_gate",
            gate(gold_ok),
            {
                "gold_coverage_status": gold_coverage.get("status", ""),
                "gold_quality_status": gold_quality.get("status", ""),
                "positive_count": gold_quality.get("positive_count"),
                "negative_control_count": gold_quality.get("negative_control_count"),
                "context_specific_positive_count": gold_distribution.get("context_specific_positive_count"),
                "unique_positive_cancer_context_count": gold_distribution.get("unique_positive_cancer_context_count"),
                "unique_positive_mechanism_axis_count": gold_distribution.get("unique_positive_mechanism_axis_count"),
                "missing_required_themes": gold_distribution.get("missing_required_themes", []),
            },
            [] if gold_ok else ["gold coverage or quality gate is not satisfied"],
            "Gold set must be broad enough for release-level benchmark calibration, not just a small seed set.",
        ),
        requirement_row(
            "cscc_real_csv_rerun_gate",
            gate(cscc_ok),
            {
                "csv_path": str(cscc_csv),
                "csv_exists": cscc_exists,
                "cscc_executed_in_matrix": cscc_executed,
                "arachidonate_executed_in_matrix": arachidonate_executed,
                "arachidonate_gate": "covered_by_release_batch_benchmark_matrix_gate",
                "fixture_recovery": cscc_recovery,
            },
            [] if cscc_ok else ["required cSCC CSV is missing or its scenario did not execute"],
            "The user-specified cSCC trait-score input must be present and executed before final completion.",
        ),
        requirement_row(
            "literature_semantics_systematic_gate",
            gate(literature_ok),
            {
                "status": literature.get("status", ""),
                "case_count": literature_metrics.get("case_count"),
                "minimum_case_count": literature_metrics.get("minimum_case_count"),
                "status_accuracy": literature_metrics.get("status_accuracy"),
                "cue_recall": literature_metrics.get("cue_recall"),
                "method_cue_recall": literature_metrics.get("method_cue_recall"),
                "required_coverage_missing_count": literature_metrics.get("required_coverage_missing_count"),
                "coverage": literature.get("coverage", {}),
            },
            [] if literature_ok else ["literature semantics benchmark or required cue coverage failed"],
            "Literature evidence must be structured and stress-tested for support, contradiction, uncertainty, and background semantics.",
        ),
        requirement_row(
            "temporal_source_heldout_rerun_gate",
            gate(heldout_ok),
            {
                "status": heldout.get("status", ""),
                "rerun_performed": heldout.get("rerun_performed"),
                "rerun_passed": heldout.get("rerun_passed", {}),
                "failed_checks": heldout.get("failed_checks", []),
            },
            [] if heldout_ok else ["temporal/source held-out validation was not rerun or did not pass"],
            "Strict validation requires current temporal/source held-out reruns, not only stale reports.",
        ),
        requirement_row(
            "release_batch_benchmark_matrix_gate",
            gate(batch_ok),
            {
                "batch_passed": batch_status.get("passed"),
                "registered_gap_count": matrix.get("registered_gap_count"),
                "executed_gap_count": matrix.get("executed_gap_count"),
                "registered_gaps": matrix.get("registered_gaps", []),
                "executed_gaps": matrix.get("executed_gaps", []),
                "matrix_status": matrix.get("status", ""),
                "scenario_traceability_status": traceability.get("status", ""),
                "scenario_traceability_failed_case_count": traceability.get("failed_case_count"),
                "scenario_traceability_failure_codes": traceability.get("failure_codes", []),
            },
            [] if batch_ok else ["release batch benchmark did not pass with executed matrix and traceability coverage"],
            "Batch benchmark must cover required cancers, metabolism themes, input modes, and structured scenario traceability at release level.",
        ),
        requirement_row(
            "parameter_optimization_gate",
            gate(parameter_ok),
            {
                "status": parameter_status,
                "action_count": parameter_report.get("action_count"),
                "required_action_count": int((parameter_report.get("severity_counts") or {}).get("required") or 0),
                "blocking_required_action_count": len(blocking_actions),
                "action_ids": [str(row.get("action_id") or "") for row in parameter_report.get("actions", []) or []],
                "blocking_action_ids": blocking_action_ids,
                "non_blocking_action_ids": sorted(NON_BLOCKING_PARAMETER_ACTION_IDS),
            },
            [] if parameter_ok else ["release parameter optimization audit still has blocking required actions"],
            "Parameter changes must be driven by release metrics and classified error types, not manual impressions or lower confidence thresholds.",
        ),
        requirement_row(
            "strict_release_readiness_gate",
            gate(readiness_ok),
            {
                "readiness_status": readiness_status,
                "workflow_status": workflow_status,
                "readiness_blocker_codes": finding_codes(readiness_report.get("blockers", []) or []),
                "workflow_blocker_codes": finding_codes(workflow_report.get("blockers", []) or []),
                "readiness_remediation_refs": remediation_refs(readiness_report.get("remediation_actions", []) or []),
            },
            [] if readiness_ok else ["strict release readiness or workflow is blocked"],
            "No default release path should be considered complete while strict readiness is blocked.",
        ),
    ]
    incomplete = [row for row in requirements if row["status"] != "satisfied"]
    return {
        "audit_version": AUDIT_VERSION,
        "created_at_utc": utc_now(),
        "workspace": str(workspace),
        "status": "complete" if not incomplete else "incomplete",
        "satisfied_count": len(requirements) - len(incomplete),
        "incomplete_count": len(incomplete),
        "requirements": requirements,
        "remaining_blockers": [
            {"requirement_id": row["requirement_id"], "blockers": row["blockers"], "evidence": row["evidence"]}
            for row in incomplete
        ],
        "interpretation_boundaries": [
            "This is a research interpretation and release-readiness audit, not clinical decision support.",
            "Manual adjudication was not performed at this stage; current validation is based on deterministic fixtures, provenance checks, held-out proxy tasks, and conservative abstention rules.",
        ],
    }


def markdown_summary(report: dict[str, Any]) -> str:
    lines = [
        "# Optimization Objective Status",
        "",
        f"Status: `{report.get('status', '')}`",
        f"Satisfied: `{report.get('satisfied_count', 0)}`",
        f"Incomplete: `{report.get('incomplete_count', 0)}`",
        "",
        "## Requirements",
        "",
    ]
    for row in report.get("requirements", []) or []:
        lines.append(f"- `{row.get('requirement_id')}`: `{row.get('status')}`")
        for blocker in row.get("blockers") or []:
            lines.append(f"  - Blocker: {blocker}")
        recovery = (row.get("evidence") or {}).get("fixture_recovery") or {}
        if recovery:
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
    boundaries = report.get("interpretation_boundaries") or []
    if boundaries:
        lines.extend(["", "## Boundaries", ""])
        for boundary in boundaries:
            lines.append(f"- {boundary}")
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit objective-level optimization completion status.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--batch-report", default=DEFAULT_BATCH_REPORT)
    parser.add_argument("--readiness-report", default=DEFAULT_READINESS_REPORT)
    parser.add_argument("--workflow-report", default=DEFAULT_WORKFLOW_REPORT)
    parser.add_argument("--parameter-report", default=DEFAULT_PARAMETER_REPORT)
    parser.add_argument("--cscc-csv", default=DEFAULT_CSSC_CSV)
    parser.add_argument("--output", default="")
    parser.add_argument("--markdown-output", default="")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    workspace = Path(args.workspace).resolve()
    report = audit_objective_status(
        workspace,
        read_json(workspace / args.batch_report),
        read_json(workspace / args.readiness_report),
        read_json(workspace / args.workflow_report),
        workspace / args.cscc_csv,
        read_json(workspace / args.parameter_report),
    )
    if args.output:
        write_json(workspace / args.output, report)
    if args.markdown_output:
        output = workspace / args.markdown_output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(markdown_summary(report), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())

