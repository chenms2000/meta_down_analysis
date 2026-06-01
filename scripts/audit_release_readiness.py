"""Audit a release batch benchmark report for publish/readiness status."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AUDIT_VERSION = "release_readiness_audit.v1"
DEFAULT_REPORT = "validation_reports/mvp_20260513T002254/batch_benchmark/release_batch_benchmark_report.json"

ACTION_CATALOG = {
    "strict_profile_not_used": {
        "priority": "P0",
        "action": "Rerun release batch benchmark with --strict-release so real fixtures and held-out reruns are enforced.",
    },
    "batch_benchmark_failed": {
        "priority": "P0",
        "action": "Inspect failing gate summaries in the batch benchmark report and rerun after repairing upstream blockers.",
    },
    "fixture_audit_failed": {
        "priority": "P0",
        "action": "Restore required real-world fixture files, verify fixture schema, then rerun strict release benchmark.",
    },
    "benchmark_matrix_failed": {
        "priority": "P0",
        "action": "Restore or add benchmark scenarios until required cancer, theme, and input-mode dimensions are executed in strict release.",
    },
    "scenario_traceability_contract_failed": {
        "priority": "P0",
        "action": "Repair benchmark scenario outputs so each executed case exposes input summaries, identity decisions, mechanism facts, evidence assertions, and conclusion metrics.",
    },
    "error_optimization_signals": {
        "priority": "P0",
        "action": "Resolve all error-severity optimization signals before publication.",
    },
    "heldout_audit_failed": {
        "priority": "P0",
        "action": "Rerun or repair temporal/source held-out validation until held-out audit passes.",
    },
    "heldout_not_rerun": {
        "priority": "P0",
        "action": "Rerun strict release benchmark with temporal/source held-out validation enabled.",
    },
    "gold_coverage_not_covered": {
        "priority": "P0",
        "action": "Add or repair gold-standard conclusions for missing cancer/theme coverage cells.",
    },
    "gold_quality_failed": {
        "priority": "P0",
        "action": "Repair gold-standard quantity, fields, source refs, or negative-trap category coverage.",
    },
    "literature_semantics_failed": {
        "priority": "P0",
        "action": "Fix literature assertion semantics; do not lower semantic benchmark thresholds.",
    },
    "scenario_execution_failed": {
        "priority": "P0",
        "action": "Inspect scenario exception details and repair parser/resolver/service execution before interpreting benchmark metrics.",
    },
    "fixture_missing": {
        "priority": "P1",
        "action": "Restore optional real-world fixture or document an explicit unavailable-fixture waiver for exploratory runs.",
    },
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def add_blocker(blockers: list[dict[str, Any]], code: str, message: str, evidence: Any = None) -> None:
    row = {"code": code, "message": message}
    if evidence is not None:
        row["evidence"] = evidence
    blockers.append(row)


def add_warning(warnings: list[dict[str, Any]], code: str, message: str, evidence: Any = None) -> None:
    row = {"code": code, "message": message}
    if evidence is not None:
        row["evidence"] = evidence
    warnings.append(row)


def signal_rows(report: dict[str, Any], severity: str) -> list[dict[str, Any]]:
    return [
        row
        for row in ((report.get("optimization_signals") or {}).get("by_code") or [])
        if str(row.get("severity") or "") == severity
    ]


def action_for_code(code: str) -> dict[str, str]:
    default = {"priority": "P1", "action": "Review the cited gate evidence and add a concrete remediation before release."}
    return dict(ACTION_CATALOG.get(code, default))


def actions_from_findings(blockers: list[dict[str, Any]], warnings: list[dict[str, Any]]) -> list[dict[str, Any]]:
    actions: dict[str, dict[str, Any]] = {}
    for finding in blockers + warnings:
        code = str(finding.get("code") or "")
        catalog = action_for_code(code)
        action = actions.setdefault(
            code,
            {
                "code": code,
                "priority": catalog["priority"],
                "action": catalog["action"],
                "finding_count": 0,
                "evidence_refs": [],
            },
        )
        action["finding_count"] += 1
        evidence = finding.get("evidence")
        if isinstance(evidence, dict):
            if evidence.get("missing_required"):
                fixture_recovery = action.setdefault("fixture_recovery", [])
                for row in evidence.get("missing_required") or []:
                    if not isinstance(row, dict):
                        continue
                    action["evidence_refs"].append(str(row.get("case_id") or ""))
                    fixture_recovery.append(
                        {
                            "case_id": str(row.get("case_id") or ""),
                            "path": str(row.get("path") or ""),
                            "expected_schema": row.get("expected_schema") or {},
                            "candidate_paths": row.get("candidate_paths") or [],
                            "recovery_action": str(row.get("recovery_action") or ""),
                            "recovery_rerun_plan": row.get("recovery_rerun_plan") or {},
                        }
                    )
            if evidence.get("cases"):
                action["evidence_refs"].extend(str(case) for case in evidence.get("cases") or [] if not isinstance(case, dict))
            if evidence.get("registered_gaps"):
                action["evidence_refs"].extend(str(gap) for gap in evidence.get("registered_gaps") or [])
            if evidence.get("executed_gaps"):
                action["evidence_refs"].extend(str(gap) for gap in evidence.get("executed_gaps") or [])
            if evidence.get("failed_cases"):
                action["evidence_refs"].extend(
                    str(row.get("case_id") or "")
                    for row in evidence.get("failed_cases") or []
                    if isinstance(row, dict)
                )
        elif isinstance(evidence, list):
            for row in evidence:
                if isinstance(row, dict):
                    action["evidence_refs"].extend(str(case) for case in row.get("cases", []) or [])
        action["evidence_refs"] = sorted(set(ref for ref in action["evidence_refs"] if ref))
    return sorted(actions.values(), key=lambda row: (row["priority"], row["code"]))


def markdown_summary(report: dict[str, Any]) -> str:
    lines = [
        "# Release Readiness Audit",
        "",
        f"Status: `{report.get('status', '')}`",
        f"Blockers: `{report.get('blocker_count', 0)}`",
        f"Warnings: `{report.get('warning_count', 0)}`",
        "",
        "## Gate Summary",
        "",
    ]
    for key, value in (report.get("gate_summary") or {}).items():
        lines.append(f"- `{key}`: `{value}`")
    blockers = report.get("blockers") or []
    if blockers:
        lines.extend(["", "## Blockers", ""])
        for row in blockers:
            lines.append(f"- `{row.get('code')}`: {row.get('message', '')}")
    actions = report.get("remediation_actions") or []
    if actions:
        lines.extend(["", "## Remediation Actions", ""])
        for row in actions:
            refs = ", ".join(row.get("evidence_refs") or [])
            suffix = f" Evidence: `{refs}`" if refs else ""
            lines.append(f"- `{row.get('priority')}` `{row.get('code')}`: {row.get('action')}{suffix}")
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
    warnings = report.get("warnings") or []
    if warnings:
        lines.extend(["", "## Warnings", ""])
        for row in warnings:
            lines.append(f"- `{row.get('code')}`: {row.get('message', '')}")
    return "\n".join(lines) + "\n"


def audit_readiness(report: dict[str, Any], require_strict_profile: bool = True) -> dict[str, Any]:
    blockers: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []
    if not report:
        add_blocker(blockers, "report_missing", "Release batch benchmark report is missing or empty.")
        result = {
            "audit_version": AUDIT_VERSION,
            "created_at_utc": utc_now(),
            "status": "blocked",
            "blockers": blockers,
            "warnings": warnings,
            "gate_summary": {},
        }
        result["remediation_actions"] = actions_from_findings(blockers, warnings)
        return result

    profile = report.get("release_profile") or {}
    status = report.get("status") or {}
    fixture = report.get("fixture_audit") or {}
    benchmark_matrix = report.get("benchmark_matrix") or {}
    traceability_contract = report.get("scenario_traceability_contract") or {}
    heldout = report.get("heldout_validation") or {}
    gold_coverage = report.get("gold_coverage") or {}
    gold_quality = report.get("gold_quality") or {}
    literature = report.get("literature_semantics") or {}
    error_signals = signal_rows(report, "error")
    warning_signals = signal_rows(report, "warning")

    if require_strict_profile and profile.get("profile") != "strict":
        add_blocker(blockers, "strict_profile_not_used", "Readiness requires a strict release profile.", profile)
    if not status.get("passed"):
        add_blocker(blockers, "batch_benchmark_failed", "Release batch benchmark status is not passed.", status)
    if int(status.get("failed_case_count") or 0) > 0:
        add_blocker(blockers, "scenario_execution_failed", "One or more benchmark scenarios failed.", status)
    if fixture.get("status") != "passed":
        add_blocker(blockers, "fixture_audit_failed", "Required real-world fixture audit failed.", fixture)
    strict_matrix_gap = int(benchmark_matrix.get("executed_gap_count") or 0) > 0
    if require_strict_profile and (benchmark_matrix.get("status") != "passed" or strict_matrix_gap):
        add_blocker(blockers, "benchmark_matrix_failed", "Strict readiness requires executed benchmark matrix coverage.", benchmark_matrix)
    if traceability_contract.get("status") not in {"", "passed"}:
        add_blocker(
            blockers,
            "scenario_traceability_contract_failed",
            "Executed benchmark scenarios did not satisfy the traceability contract.",
            traceability_contract,
        )
    if heldout.get("status") != "passed":
        add_blocker(blockers, "heldout_audit_failed", "Temporal/source held-out audit failed.", heldout)
    if require_strict_profile and not heldout.get("rerun_performed"):
        add_blocker(blockers, "heldout_not_rerun", "Strict readiness requires temporal/source held-out rerun.", heldout)
    if gold_coverage.get("status") != "covered":
        add_blocker(blockers, "gold_coverage_not_covered", "Gold standard coverage audit is not fully covered.", gold_coverage)
    if gold_quality.get("status") != "passed":
        add_blocker(blockers, "gold_quality_failed", "Gold standard quality audit failed.", gold_quality)
    if literature.get("status") != "passed":
        add_blocker(blockers, "literature_semantics_failed", "Literature semantics benchmark failed.", literature)
    if error_signals:
        add_blocker(blockers, "error_optimization_signals", "Error-severity optimization signals remain.", error_signals)
    for signal in warning_signals:
        add_warning(warnings, str(signal.get("code") or "warning_signal"), "Warning-severity optimization signal remains.", signal)

    result = {
        "audit_version": AUDIT_VERSION,
        "created_at_utc": utc_now(),
        "status": "ready" if not blockers else "blocked",
        "blocker_count": len(blockers),
        "warning_count": len(warnings),
        "blockers": blockers,
        "warnings": warnings,
        "gate_summary": {
            "release_profile": profile.get("profile", ""),
            "batch_passed": bool(status.get("passed")),
            "failed_case_count": int(status.get("failed_case_count") or 0),
            "missing_case_count": int(status.get("missing_case_count") or 0),
            "fixture_status": fixture.get("status", ""),
            "benchmark_matrix_status": benchmark_matrix.get("status", ""),
            "benchmark_matrix_registered_gap_count": int(benchmark_matrix.get("registered_gap_count") or 0),
            "benchmark_matrix_executed_gap_count": int(benchmark_matrix.get("executed_gap_count") or 0),
            "scenario_traceability_status": traceability_contract.get("status", ""),
            "scenario_traceability_failed_case_count": int(traceability_contract.get("failed_case_count") or 0),
            "heldout_status": heldout.get("status", ""),
            "heldout_rerun_performed": bool(heldout.get("rerun_performed")),
            "gold_coverage_status": gold_coverage.get("status", ""),
            "gold_quality_status": gold_quality.get("status", ""),
            "literature_semantics_status": literature.get("status", ""),
            "error_signal_count": len(error_signals),
            "warning_signal_count": len(warning_signals),
        },
    }
    result["remediation_actions"] = actions_from_findings(blockers, warnings)
    return result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit release readiness from a batch benchmark report.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--report", default=DEFAULT_REPORT)
    parser.add_argument("--output", default="")
    parser.add_argument("--markdown-output", default="")
    parser.add_argument("--allow-exploratory", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    workspace = Path(args.workspace).resolve()
    report_path = workspace / args.report
    report = audit_readiness(read_json(report_path), require_strict_profile=not args.allow_exploratory)
    report["input_report"] = str(report_path)
    if args.output:
        write_json(workspace / args.output, report)
    if args.markdown_output:
        output = workspace / args.markdown_output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(markdown_summary(report), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if report["status"] == "ready" else 1


if __name__ == "__main__":
    raise SystemExit(main())
