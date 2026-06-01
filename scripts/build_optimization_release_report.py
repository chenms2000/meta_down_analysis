#!/usr/bin/env python3
"""Build a paper-ready optimization summary for a frozen release."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


REPORT_VERSION = "optimization_release_report.20260520"
DEFAULT_OUTPUT_ROOT = "validation_reports"
DEFAULT_LLM_REPORT = "validation_reports/llm_safe_adapter_regression_report.json"
DEFAULT_RELEASE_GATE_ROOT = "learning_runs"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def maybe_read_json(path: Path) -> dict[str, Any]:
    return read_json(path) if path.is_file() else {}


def status_label(value: Any) -> str:
    text = str(value or "").strip()
    return text or "not_available"


def markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(item).replace("\n", " ") for item in row) + " |")
    return "\n".join(lines)


def phase15_path(workspace: Path, release_id: str) -> Path:
    return workspace / DEFAULT_OUTPUT_ROOT / release_id / "phase15_validation_report.json"


def release_gate_path(workspace: Path, run_id: str) -> Path:
    return workspace / DEFAULT_RELEASE_GATE_ROOT / run_id / "reports" / "release_gate_validation_report.json"


def name_audit_path(workspace: Path, release_id: str) -> Path:
    return workspace / DEFAULT_OUTPUT_ROOT / release_id / "name_audit" / "compound_name_audit_report.json"


def name_audit_triage_path(workspace: Path, release_id: str) -> Path:
    return workspace / DEFAULT_OUTPUT_ROOT / release_id / "name_audit" / "compound_name_triage_summary.json"


def manual_adjudication_queue_path(workspace: Path, release_id: str) -> Path:
    return workspace / DEFAULT_OUTPUT_ROOT / release_id / "name_audit" / "manual_adjudication_queue_summary.json"


def temporal_holdout_path(workspace: Path, run_id: str) -> Path:
    return workspace / DEFAULT_RELEASE_GATE_ROOT / run_id / "reports" / "temporal_holdout" / "temporal_holdout_validation_report.json"


def source_heldout_path(workspace: Path, run_id: str) -> Path:
    return workspace / DEFAULT_RELEASE_GATE_ROOT / run_id / "reports" / "source_heldout" / "source_heldout_validation_report.json"


def resolve_path(value: str, workspace: Path) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (workspace / path).resolve()


def optimization_rows(phase15: dict[str, Any]) -> list[list[Any]]:
    optimization = phase15.get("optimization_summary") or {}
    entity = optimization.get("entity_resolution") or {}
    evidence = optimization.get("evidence_precision") or {}
    ranking = optimization.get("ranking_calibration") or {}
    readiness = optimization.get("release_readiness") or {}
    return [
        [
            "entity_resolution",
            status_label(entity.get("status")),
            f"resolver={entity.get('resolver_case_count', 0)}, precheck={entity.get('precheck_case_count', 0)}, abstention={entity.get('abstention_case_count', 0)}",
        ],
        [
            "literature_evidence_precision",
            status_label(evidence.get("status")),
            ", ".join(f"{key}={value}" for key, value in sorted((evidence.get("checks") or {}).items())),
        ],
        [
            "ranking_calibration",
            status_label(ranking.get("status")),
            f"calibration_failures={ranking.get('calibration_failure_count', 0)}, traceability_failures={ranking.get('traceability_failure_count', 0)}",
        ],
        [
            "release_readiness",
            "passed" if not readiness.get("failed_case_count") and not readiness.get("blocked_case_count") else "review_required",
            f"failed={readiness.get('failed_case_count', 0)}, blocked={readiness.get('blocked_case_count', 0)}",
        ],
    ]


def validation_rows(phase15: dict[str, Any], llm: dict[str, Any], release_gate: dict[str, Any]) -> list[list[Any]]:
    summary = phase15.get("summary") or {}
    evidence_validation = phase15.get("evidence_validation") or {}
    literature = phase15.get("literature_evidence") or {}
    llm_summary = llm.get("summary") or {}
    gate_summary = release_gate.get("summary") or {}
    return [
        [
            "phase15_fixture_gate",
            status_label(summary.get("gate_status")),
            f"{summary.get('passed', 0)}/{summary.get('total', 0)} passed",
        ],
        [
            "evidence_contract",
            status_label(evidence_validation.get("status")),
            f"support={evidence_validation.get('support_row_count', 0)}, sentences={evidence_validation.get('sentence_count', 0)}",
        ],
        [
            "literature_manifest",
            status_label(literature.get("status")),
            f"candidates={literature.get('evidence_candidate_count', 0)}, supported_edges={literature.get('supported_existing_edge_count', 0)}",
        ],
        [
            "llm_safe_adapter",
            status_label(llm_summary.get("gate_status")),
            f"local_passed={llm_summary.get('local_passed', 0)}, adversarial_passed={llm_summary.get('adversarial_probe_passed', 0)}",
        ],
        [
            "learning_release_gate",
            status_label(release_gate.get("gate_status")),
            f"{gate_summary.get('passed', 0)}/{gate_summary.get('total', 0)} passed",
        ],
    ]


def fixture_benchmark_rows(phase15: dict[str, Any]) -> list[list[Any]]:
    """Summarize fixed validation fixtures as proxy benchmark rows."""
    grouped: dict[str, dict[str, Any]] = {}
    for case in phase15.get("cases") or []:
        fixture = str(case.get("fixture") or "unknown")
        row = grouped.setdefault(
            fixture,
            {
                "case_count": 0,
                "passed": 0,
                "failed": 0,
                "blocked": 0,
                "resolver": 0,
                "precheck": 0,
                "abstain": 0,
            },
        )
        row["case_count"] += 1
        status = str(case.get("status") or "")
        if status in {"passed", "failed", "blocked"}:
            row[status] += 1
        for validation in case.get("validations") or []:
            endpoint = validation.get("endpoint")
            if endpoint == "/resolve":
                row["resolver"] += 1
            if endpoint == "/precheck/metabolites":
                row["precheck"] += 1
                if validation.get("status") in {"ambiguous", "unmatched", "invalid"}:
                    row["abstain"] += 1
    return [
        [
            fixture,
            values["case_count"],
            f"{values['passed']}/{values['case_count']}",
            values["failed"],
            values["blocked"],
            values["resolver"],
            values["precheck"],
            values["abstain"],
        ]
        for fixture, values in sorted(grouped.items())
    ]


def global_validation_rows(
    name_audit: dict[str, Any],
    name_triage: dict[str, Any],
    temporal: dict[str, Any],
    source: dict[str, Any],
) -> list[list[Any]]:
    temporal_metrics = temporal.get("validation_metrics") or {}
    source_rows = source.get("heldout_sources") or []
    total_source_heldout_rows = sum(int(row.get("heldout_rows") or 0) for row in source_rows)
    source_names = [
        str((row.get("heldout") or {}).get("value") or "")
        for row in source_rows
        if (row.get("heldout") or {}).get("value")
    ]
    audit_summary = name_audit.get("summary") or {}
    action_counts = audit_summary.get("by_recommended_action") or {}
    triage_counts = name_triage.get("policy_counts") or {}
    return [
        [
            "compound_name_audit",
            "available" if name_audit else "missing",
            f"risk_rows={name_audit.get('risk_row_count', 0)}, online_checks={name_audit.get('online_check_count', 0)}, actions={len(action_counts)}",
        ],
        [
            "compound_name_triage",
            "available" if name_triage else "missing",
            f"risk_rows={name_triage.get('total_risk_rows', 0)}, policies={len(triage_counts)}, auto_abstain={triage_counts.get('auto_abstain_name_only', 0)}, require_id_or_feature={triage_counts.get('requires_identifier_or_orthogonal_feature', 0)}",
        ],
        [
            "temporal_holdout",
            "available" if temporal else "missing",
            f"rows={temporal_metrics.get('heldout_rows', 0)}, roc_auc={temporal_metrics.get('heldout_roc_auc', 'n/a')}, lift@50={temporal_metrics.get('heldout_lift_at_50', 'n/a')}",
        ],
        [
            "source_heldout",
            "available" if source else "missing",
            f"sources={len(source_rows)} ({', '.join(source_names)}), rows={total_source_heldout_rows}",
        ],
    ]


def source_heldout_rows(source: dict[str, Any]) -> list[list[Any]]:
    rows = []
    for item in source.get("heldout_sources") or []:
        metrics = item.get("validation_metrics") or {}
        heldout = item.get("heldout") or {}
        rows.append(
            [
                heldout.get("value", ""),
                item.get("heldout_rows", 0),
                item.get("heldout_positive_rows", 0),
                metrics.get("heldout_roc_auc", "n/a"),
                metrics.get("heldout_average_precision", "n/a"),
                metrics.get("heldout_positive_top_5pct_fraction", "n/a"),
                metrics.get("heldout_top_score_tie_fraction", "n/a"),
            ]
        )
    return rows


def manual_queue_rows(queue: dict[str, Any]) -> list[list[Any]]:
    return [[key, value] for key, value in sorted((queue.get("bucket_counts") or {}).items())]


def build_report(
    workspace: Path,
    release_id: str,
    run_id: str,
    phase15_report: Path,
    llm_report: Path,
    release_gate_report: Path,
    name_audit_report: Path,
    name_triage_report: Path,
    manual_queue_report: Path,
    temporal_report: Path,
    source_report: Path,
) -> dict[str, Any]:
    phase15 = maybe_read_json(phase15_report)
    llm = maybe_read_json(llm_report)
    release_gate = maybe_read_json(release_gate_report)
    name_audit = maybe_read_json(name_audit_report)
    name_triage = maybe_read_json(name_triage_report)
    manual_queue = maybe_read_json(manual_queue_report)
    temporal = maybe_read_json(temporal_report)
    source = maybe_read_json(source_report)
    return {
        "report_version": REPORT_VERSION,
        "release_id": release_id,
        "run_id": run_id,
        "inputs": {
            "phase15_report": str(phase15_report),
            "llm_regression_report": str(llm_report),
            "release_gate_report": str(release_gate_report),
            "compound_name_audit_report": str(name_audit_report),
            "compound_name_triage_report": str(name_triage_report),
            "manual_adjudication_queue_report": str(manual_queue_report),
            "temporal_holdout_report": str(temporal_report),
            "source_heldout_report": str(source_report),
        },
        "phase15_summary": phase15.get("summary", {}),
        "optimization_summary": phase15.get("optimization_summary", {}),
        "fixture_benchmarks": fixture_benchmark_rows(phase15),
        "evidence_validation": phase15.get("evidence_validation", {}),
        "literature_evidence": phase15.get("literature_evidence", {}),
        "llm_safe_adapter_summary": llm.get("summary", {}),
        "release_gate_summary": release_gate.get("summary", {}),
        "release_gate_status": release_gate.get("gate_status", "not_available"),
        "compound_name_audit_summary": name_audit.get("summary", {}),
        "compound_name_audit_counts": {
            "risk_row_count": name_audit.get("risk_row_count", 0),
            "online_check_count": name_audit.get("online_check_count", 0),
        },
        "compound_name_triage_summary": {
            "triage_version": name_triage.get("triage_version", ""),
            "total_risk_rows": name_triage.get("total_risk_rows", 0),
            "policy_counts": name_triage.get("policy_counts", {}),
            "release_policy": name_triage.get("release_policy", {}),
        },
        "manual_adjudication_queue": {
            "queue_version": manual_queue.get("queue_version", ""),
            "selected_rows": manual_queue.get("selected_rows", 0),
            "sample_per_bucket": manual_queue.get("sample_per_bucket", 0),
            "bucket_counts": manual_queue.get("bucket_counts", {}),
            "output_csv": manual_queue.get("output_csv", ""),
            "interpretation": manual_queue.get("interpretation", ""),
        },
        "temporal_holdout_validation": temporal.get("validation_metrics", {}),
        "temporal_holdout_input_rows": temporal.get("input_rows", {}),
        "source_heldout_validation": source.get("heldout_sources", []),
        "source_heldout_input_rows": source.get("input_rows", {}),
        "phase15_report_hash": phase15.get("report_hash", ""),
        "llm_report_hash": llm.get("report_hash", ""),
    }


def render_markdown(report: dict[str, Any]) -> str:
    release_id = report.get("release_id", "")
    phase15 = {
        "summary": report.get("phase15_summary", {}),
        "optimization_summary": report.get("optimization_summary", {}),
        "evidence_validation": report.get("evidence_validation", {}),
        "literature_evidence": report.get("literature_evidence", {}),
    }
    llm = {"summary": report.get("llm_safe_adapter_summary", {})}
    release_gate = {
        "gate_status": report.get("release_gate_status", "not_available"),
        "summary": report.get("release_gate_summary", {}),
    }
    name_audit = {
        "summary": report.get("compound_name_audit_summary", {}),
        **(report.get("compound_name_audit_counts") or {}),
    }
    name_triage = report.get("compound_name_triage_summary") or {}
    temporal = {
        "validation_metrics": report.get("temporal_holdout_validation", {}),
        "input_rows": report.get("temporal_holdout_input_rows", {}),
    }
    source = {"heldout_sources": report.get("source_heldout_validation", [])}
    manual_queue = report.get("manual_adjudication_queue") or {}
    lines = [
        "# Optimization Release Summary",
        "",
        f"Release: `{release_id}`",
        f"Learning run: `{report.get('run_id', '')}`",
        "",
        "## Optimization Gates",
        "",
        markdown_table(["area", "status", "details"], optimization_rows(phase15)),
        "",
        "## Validation Evidence",
        "",
        markdown_table(["gate", "status", "summary"], validation_rows(phase15, llm, release_gate)),
        "",
        "## Fixture Benchmarks",
        "",
        markdown_table(
            ["fixture", "cases", "passed", "failed", "blocked", "resolver_checks", "precheck_checks", "abstain_or_ambiguous"],
            report.get("fixture_benchmarks") or [],
        ),
        "",
        "## Global Validation",
        "",
        markdown_table(["area", "status", "summary"], global_validation_rows(name_audit, name_triage, temporal, source)),
        "",
        "## Compound Name Risk Handling",
        "",
        markdown_table(
            ["policy", "count", "release handling"],
            [
                [policy, count, (name_triage.get("release_policy") or {}).get(policy, "")]
                for policy, count in sorted((name_triage.get("policy_counts") or {}).items())
            ],
        ),
        "",
        "## Manual Adjudication Queue",
        "",
        markdown_table(
            ["field", "value"],
            [
                ["selected_rows", manual_queue.get("selected_rows", 0)],
                ["sample_per_bucket", manual_queue.get("sample_per_bucket", 0)],
                ["output_csv", manual_queue.get("output_csv", "")],
                ["interpretation", manual_queue.get("interpretation", "")],
            ],
        ),
        "",
        markdown_table(["review_bucket", "selected"], manual_queue_rows(manual_queue)),
        "",
        "## Source-Held-Out Details",
        "",
        markdown_table(
            ["source", "heldout_rows", "positive_rows", "roc_auc", "average_precision", "top_5pct_positive_fraction", "top_score_tie_fraction"],
            source_heldout_rows(source),
        ),
        "",
        "## Reproducibility",
        "",
        markdown_table(
            ["artifact", "value"],
            [
                ["phase15_report_hash", report.get("phase15_report_hash", "")],
                ["llm_report_hash", report.get("llm_report_hash", "")],
                ["phase15_report", report.get("inputs", {}).get("phase15_report", "")],
                ["llm_regression_report", report.get("inputs", {}).get("llm_regression_report", "")],
                ["release_gate_report", report.get("inputs", {}).get("release_gate_report", "")],
                ["compound_name_audit_report", report.get("inputs", {}).get("compound_name_audit_report", "")],
                ["compound_name_triage_report", report.get("inputs", {}).get("compound_name_triage_report", "")],
                ["manual_adjudication_queue_report", report.get("inputs", {}).get("manual_adjudication_queue_report", "")],
                ["temporal_holdout_report", report.get("inputs", {}).get("temporal_holdout_report", "")],
                ["source_heldout_report", report.get("inputs", {}).get("source_heldout_report", "")],
            ],
        ),
        "",
        "## Interpretation Boundary",
        "",
        "- This summary describes an optional research interpretation plugin, not the core algorithm.",
        "- Confidence tiers are research-prioritization strata, not biological truth probabilities.",
        "- Fixture benchmarks are proxy resolver/ambiguity checks; they are not a full biological gold standard.",
        "- Compound-name risk rows were not manually adjudicated one by one in this release; they are handled through deterministic triage policies that abstain, require identifiers/orthogonal features, warn on formatting, or monitor/downweight ambiguous surfaces.",
        "- Literature evidence remains an overlay at this stage; novel/conflict evidence is not promoted to canonical facts before a future human-adjudicated optimization cycle.",
        "- Manual adjudication is temporarily unavailable for this release summary and is listed as a follow-up optimization step; current validation relies on deterministic fixtures, provenance contracts, held-out proxy tasks, and conservative abstention.",
        "- LLM output is accepted only when source-bound and guard-validated.",
        "",
    ]
    return "\n".join(lines)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build an optimization release summary.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--run-id", default="")
    parser.add_argument("--phase15-report", default="")
    parser.add_argument("--llm-report", default=DEFAULT_LLM_REPORT)
    parser.add_argument("--release-gate-report", default="")
    parser.add_argument("--compound-name-audit-report", default="")
    parser.add_argument("--compound-name-triage-report", default="")
    parser.add_argument("--manual-adjudication-queue-report", default="")
    parser.add_argument("--temporal-holdout-report", default="")
    parser.add_argument("--source-heldout-report", default="")
    parser.add_argument("--output-json", default="")
    parser.add_argument("--output-md", default="")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    workspace = Path(args.workspace).resolve()
    phase15_report = resolve_path(args.phase15_report, workspace) if args.phase15_report else phase15_path(workspace, args.release_id)
    llm_report = resolve_path(args.llm_report, workspace)
    release_gate_report = (
        resolve_path(args.release_gate_report, workspace)
        if args.release_gate_report
        else release_gate_path(workspace, args.run_id)
        if args.run_id
        else Path("")
    )
    name_audit_report = (
        resolve_path(args.compound_name_audit_report, workspace)
        if args.compound_name_audit_report
        else name_audit_path(workspace, args.release_id)
    )
    name_triage_report = (
        resolve_path(args.compound_name_triage_report, workspace)
        if args.compound_name_triage_report
        else name_audit_triage_path(workspace, args.release_id)
    )
    manual_queue_report = (
        resolve_path(args.manual_adjudication_queue_report, workspace)
        if args.manual_adjudication_queue_report
        else manual_adjudication_queue_path(workspace, args.release_id)
    )
    temporal_report = (
        resolve_path(args.temporal_holdout_report, workspace)
        if args.temporal_holdout_report
        else temporal_holdout_path(workspace, args.run_id)
        if args.run_id
        else Path("")
    )
    source_report = (
        resolve_path(args.source_heldout_report, workspace)
        if args.source_heldout_report
        else source_heldout_path(workspace, args.run_id)
        if args.run_id
        else Path("")
    )
    report = build_report(
        workspace,
        args.release_id,
        args.run_id,
        phase15_report,
        llm_report,
        release_gate_report,
        name_audit_report,
        name_triage_report,
        manual_queue_report,
        temporal_report,
        source_report,
    )
    default_dir = workspace / DEFAULT_OUTPUT_ROOT / args.release_id
    output_json = resolve_path(args.output_json, workspace) if args.output_json else default_dir / "optimization_release_summary.json"
    output_md = resolve_path(args.output_md, workspace) if args.output_md else default_dir / "optimization_release_summary.md"
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    output_md.write_text(render_markdown(report), encoding="utf-8")
    print(json.dumps({"release_id": args.release_id, "output_json": str(output_json), "output_md": str(output_md)}, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
