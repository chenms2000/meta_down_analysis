"""Run the strict release workflow with fixture-first failure gates.

The batch benchmark can rerun expensive held-out validations. This wrapper
audits required real-world fixtures first, then runs the strict benchmark and
readiness audit only when the required input surface is present and valid.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from audit_real_world_fixtures import audit_registry, read_json as read_fixture_json, write_json  # noqa: E402
from audit_release_readiness import audit_readiness, markdown_summary  # noqa: E402
from metabo_service import latest_release_id  # noqa: E402


WORKFLOW_VERSION = "strict_release_workflow.v1"
DEFAULT_OUTPUT_ROOT = "validation_reports"
DEFAULT_RUN_ID = "full_context_engine_20260519"
DEFAULT_FIXTURE_REGISTRY = "config/real_world_fixture_registry.json"
DEFAULT_CSSC_CSV = "trait_score_diff_cSCC_Epi_cSCC_Tumor_vs_Adjacent.csv"
NON_BLOCKING_PARAMETER_ACTION_IDS = {"objective_incomplete_release_block"}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def run_command(command: list[str], workspace: Path) -> dict[str, Any]:
    completed = subprocess.run(command, cwd=workspace, capture_output=True, text=True, encoding="utf-8")
    return {
        "command": command,
        "returncode": completed.returncode,
        "stdout_tail": completed.stdout[-4000:],
        "stderr_tail": completed.stderr[-4000:],
        "passed": completed.returncode == 0,
    }


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


def parameter_optimization_summary(parameter_report: dict[str, Any]) -> dict[str, Any]:
    blocking_actions = blocking_parameter_actions(parameter_report)
    return {
        "status": parameter_report.get("status", ""),
        "action_count": parameter_report.get("action_count", 0),
        "required_action_count": (parameter_report.get("severity_counts") or {}).get("required", 0),
        "blocking_required_action_count": len(blocking_actions),
        "blocking_action_ids": [str(row.get("action_id") or "") for row in blocking_actions],
        "non_blocking_action_ids": sorted(NON_BLOCKING_PARAMETER_ACTION_IDS),
    }


def workflow_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Strict Release Workflow",
        "",
        f"Status: `{report.get('status', '')}`",
        f"Workflow version: `{report.get('workflow_version', '')}`",
        f"Release ID: `{report.get('release_id', '')}`",
        "",
        "## Gates",
        "",
    ]
    for key, value in (report.get("gate_summary") or {}).items():
        lines.append(f"- `{key}`: `{value}`")
    blockers = report.get("blockers") or []
    if blockers:
        lines.extend(["", "## Blockers", ""])
        for row in blockers:
            lines.append(f"- `{row.get('code')}`: {row.get('message', '')}")
    outputs = report.get("outputs") or {}
    if outputs:
        lines.extend(["", "## Outputs", ""])
        for key, value in outputs.items():
            lines.append(f"- `{key}`: `{value}`")
    return "\n".join(lines) + "\n"


def fixture_blockers(fixture_report: dict[str, Any]) -> list[dict[str, Any]]:
    if fixture_report.get("status") == "passed":
        return []
    failed_required = fixture_report.get("failed_required") or []
    return [
        {
            "code": "real_world_fixture_audit_failed",
            "message": "Required real-world fixture audit failed before strict benchmark execution.",
            "evidence": {
                "failed_required_count": len(failed_required),
                "case_ids": [str(row.get("case_id") or "") for row in failed_required],
            },
        }
    ]


def append_blocker_once(blockers: list[dict[str, Any]], blocker: dict[str, Any]) -> None:
    code = str(blocker.get("code") or "")
    if code and code in {str(row.get("code") or "") for row in blockers}:
        return
    blockers.append(blocker)


def build_batch_command(args: argparse.Namespace, release_id: str) -> list[str]:
    command = [
        sys.executable,
        str(SCRIPT_DIR / "run_release_batch_benchmark.py"),
        "--workspace",
        ".",
        "--release-id",
        release_id,
        "--run-id",
        args.run_id,
        "--output-root",
        args.output_root,
        "--fixture-registry",
        args.fixture_registry,
        "--strict-release",
        "--heldout-top-n",
        str(args.heldout_top_n),
        "--max-paths",
        str(args.max_paths),
        "--max-hops",
        str(args.max_hops),
        "--scenario-timeout-seconds",
        str(args.scenario_timeout_seconds),
        "--heldout-rerun-timeout-seconds",
        str(args.heldout_rerun_timeout_seconds),
    ]
    if args.cscc_csv:
        command.extend(["--cscc-csv", args.cscc_csv])
    return command


def build_objective_command(args: argparse.Namespace, release_id: str) -> list[str]:
    output_root = Path(args.output_root) / release_id
    return [
        sys.executable,
        str(SCRIPT_DIR / "audit_optimization_objective_status.py"),
        "--workspace",
        ".",
        "--batch-report",
        str(output_root / "batch_benchmark" / "release_batch_benchmark_report.json"),
        "--readiness-report",
        str(output_root / "batch_benchmark" / "release_readiness_audit.json"),
        "--workflow-report",
        str(output_root / "strict_release_workflow" / "strict_release_workflow_report.json"),
        "--parameter-report",
        str(output_root / "release_parameter_optimization_audit.json"),
        "--cscc-csv",
        args.cscc_csv or DEFAULT_CSSC_CSV,
        "--output",
        str(output_root / "optimization_objective_status.json"),
        "--markdown-output",
        str(output_root / "optimization_objective_status.md"),
    ]


def build_metric_regression_command(args: argparse.Namespace, release_id: str) -> list[str]:
    output_root = Path(args.output_root) / release_id
    command = [
        sys.executable,
        str(SCRIPT_DIR / "audit_release_metric_regression.py"),
        "--workspace",
        ".",
        "--batch-report",
        str(output_root / "batch_benchmark" / "release_batch_benchmark_report.json"),
        "--readiness-report",
        str(output_root / "batch_benchmark" / "release_readiness_audit.json"),
        "--objective-report",
        str(output_root / "optimization_objective_status.json"),
        "--snapshot-output",
        str(output_root / "release_metric_snapshot.json"),
        "--output",
        str(output_root / "release_metric_regression_report.json"),
        "--markdown-output",
        str(output_root / "release_metric_regression_report.md"),
    ]
    if args.metric_baseline:
        command.extend(["--baseline", args.metric_baseline])
    return command


def build_parameter_optimization_command(args: argparse.Namespace, release_id: str) -> list[str]:
    output_root = Path(args.output_root) / release_id
    return [
        sys.executable,
        str(SCRIPT_DIR / "audit_release_parameter_optimization.py"),
        "--workspace",
        ".",
        "--batch-report",
        str(output_root / "batch_benchmark" / "release_batch_benchmark_report.json"),
        "--readiness-report",
        str(output_root / "batch_benchmark" / "release_readiness_audit.json"),
        "--objective-report",
        str(output_root / "optimization_objective_status.json"),
        "--metric-snapshot",
        str(output_root / "release_metric_snapshot.json"),
        "--output",
        str(output_root / "release_parameter_optimization_audit.json"),
        "--markdown-output",
        str(output_root / "release_parameter_optimization_audit.md"),
    ]


def post_audit_status(payload: dict[str, Any], status_key: str, success_statuses: set[str]) -> str:
    if not payload:
        return ""
    status = str(payload.get(status_key) or "")
    return status if status in success_statuses or status else ""


def run_workflow(args: argparse.Namespace) -> dict[str, Any]:
    workspace = Path(args.workspace).resolve()
    release_id = args.release_id or latest_release_id(workspace / "graph_projection")
    workflow_dir = workspace / args.output_root / release_id / "strict_release_workflow"
    batch_dir = workspace / args.output_root / release_id / "batch_benchmark"
    workflow_json = workflow_dir / "strict_release_workflow_report.json"
    workflow_md = workflow_dir / "strict_release_workflow_report.md"
    fixture_output = workflow_dir / "real_world_fixture_audit.json"
    readiness_output = batch_dir / "release_readiness_audit.json"
    readiness_md = batch_dir / "release_readiness_audit.md"
    batch_report_path = batch_dir / "release_batch_benchmark_report.json"
    objective_output = workspace / args.output_root / release_id / "optimization_objective_status.json"
    objective_md = workspace / args.output_root / release_id / "optimization_objective_status.md"
    metric_snapshot_output = workspace / args.output_root / release_id / "release_metric_snapshot.json"
    metric_regression_output = workspace / args.output_root / release_id / "release_metric_regression_report.json"
    metric_regression_md = workspace / args.output_root / release_id / "release_metric_regression_report.md"
    parameter_optimization_output = workspace / args.output_root / release_id / "release_parameter_optimization_audit.json"
    parameter_optimization_md = workspace / args.output_root / release_id / "release_parameter_optimization_audit.md"

    fixture_registry = read_fixture_json(workspace / args.fixture_registry)
    fixture_report = audit_registry(workspace, fixture_registry)
    write_json(fixture_output, fixture_report)

    blockers = fixture_blockers(fixture_report)
    batch_result: dict[str, Any] | None = None
    readiness_report: dict[str, Any] | None = None
    readiness_written = False
    batch_report_updated = False

    if not blockers or args.continue_on_fixture_failure:
        batch_command = build_batch_command(args, release_id)
        batch_started_at = time.time()
        batch_result = run_command(batch_command, workspace)
        batch_report_updated = batch_report_path.exists() and batch_report_path.stat().st_mtime >= batch_started_at - 0.001
        batch_report = read_json(batch_report_path) if batch_report_updated else {}
        if batch_report:
            readiness_report = audit_readiness(batch_report, require_strict_profile=True)
            readiness_report["input_report"] = str(batch_report_path)
            write_json(readiness_output, readiness_report)
            readiness_md.parent.mkdir(parents=True, exist_ok=True)
            readiness_md.write_text(markdown_summary(readiness_report), encoding="utf-8")
            readiness_written = True
            if readiness_report.get("status") != "ready":
                blockers.append(
                    {
                        "code": "release_readiness_blocked",
                        "message": "Release readiness audit blocked publication after strict benchmark.",
                        "evidence": {"blocker_count": readiness_report.get("blocker_count", 0)},
                    }
                )
        else:
            blockers.append(
                {
                    "code": "batch_report_missing",
                    "message": "Strict batch benchmark did not produce a readable report.",
                    "evidence": {"path": str(batch_report_path)},
                }
            )
        if batch_result and not batch_result.get("passed"):
            blockers.append(
                {
                    "code": "strict_batch_benchmark_failed",
                    "message": "Strict batch benchmark command returned a non-zero exit code.",
                    "evidence": {"returncode": batch_result.get("returncode")},
                }
            )

    status = "ready" if not blockers else "blocked"
    report = {
        "workflow_version": WORKFLOW_VERSION,
        "created_at_utc": utc_now(),
        "status": status,
        "workspace": str(workspace),
        "release_id": release_id,
        "run_id": args.run_id,
        "fixture_registry": str((workspace / args.fixture_registry).resolve()),
        "continue_on_fixture_failure": bool(args.continue_on_fixture_failure),
        "gate_summary": {
            "fixture_audit_status": fixture_report.get("status", ""),
            "fixture_required_failed_count": fixture_report.get("required_failed_count", 0),
            "batch_benchmark_executed": batch_result is not None,
            "batch_benchmark_returncode": None if batch_result is None else batch_result.get("returncode"),
            "batch_report_updated": batch_report_updated,
            "readiness_audit_status": "" if readiness_report is None else readiness_report.get("status", ""),
            "readiness_audit_written": readiness_written,
        },
        "blocker_count": len(blockers),
        "blockers": blockers,
        "fixture_audit": fixture_report,
        "batch_benchmark": batch_result or {"skipped": True, "reason": "fixture_audit_failed"},
        "readiness_audit": readiness_report or {"skipped": True},
        "post_release_audits": {"skipped": True},
        "outputs": {
            "workflow_json": str(workflow_json),
            "workflow_markdown": str(workflow_md),
            "fixture_audit_json": str(fixture_output),
            "batch_report_json": str(batch_report_path),
            "readiness_json": str(readiness_output),
            "readiness_markdown": str(readiness_md),
            "objective_status_json": str(objective_output),
            "objective_status_markdown": str(objective_md),
            "metric_snapshot_json": str(metric_snapshot_output),
            "metric_regression_json": str(metric_regression_output),
            "metric_regression_markdown": str(metric_regression_md),
            "parameter_optimization_json": str(parameter_optimization_output),
            "parameter_optimization_markdown": str(parameter_optimization_md),
        },
    }
    write_json(workflow_json, report)
    workflow_md.parent.mkdir(parents=True, exist_ok=True)
    workflow_md.write_text(workflow_markdown(report), encoding="utf-8")

    post_audits: dict[str, Any] = {"skipped": bool(args.skip_post_audits)}
    if not args.skip_post_audits:
        objective_result = run_command(build_objective_command(args, release_id), workspace)
        objective_report = read_json(objective_output)
        metric_result = run_command(build_metric_regression_command(args, release_id), workspace)
        metric_report = read_json(metric_regression_output)
        parameter_result = run_command(build_parameter_optimization_command(args, release_id), workspace)
        parameter_report = read_json(parameter_optimization_output)
        post_audits = {
            "skipped": False,
            "objective_status": {
                "command": objective_result,
                "status": objective_report.get("status", ""),
                "satisfied_count": objective_report.get("satisfied_count", 0),
                "incomplete_count": objective_report.get("incomplete_count", 0),
            },
            "metric_regression": {
                "command": metric_result,
                "status": (metric_report.get("comparison") or {}).get("status", ""),
                "regression_count": (metric_report.get("comparison") or {}).get("regression_count", 0),
            },
            "parameter_optimization": {
                "command": parameter_result,
                **parameter_optimization_summary(parameter_report),
            },
        }
        report["post_release_audits"] = post_audits
        report["gate_summary"]["objective_status"] = post_audits["objective_status"]["status"]
        report["gate_summary"]["metric_regression_status"] = post_audits["metric_regression"]["status"]
        report["gate_summary"]["metric_regression_count"] = post_audits["metric_regression"]["regression_count"]
        report["gate_summary"]["parameter_optimization_status"] = post_audits["parameter_optimization"]["status"]
        report["gate_summary"]["parameter_optimization_required_action_count"] = post_audits["parameter_optimization"]["required_action_count"]
        report["gate_summary"]["parameter_optimization_blocking_required_action_count"] = post_audits["parameter_optimization"][
            "blocking_required_action_count"
        ]
        # Post audits make the terminal release contract explicit: the objective
        # audit must be complete, and metric regressions block only when a
        # baseline comparison detects deterioration.
        if post_audits["objective_status"]["status"] != "complete":
            append_blocker_once(
                blockers,
                {
                    "code": "objective_status_incomplete",
                    "message": "Optimization objective status audit is not complete.",
                    "evidence": {
                        "status": post_audits["objective_status"]["status"],
                        "incomplete_count": post_audits["objective_status"]["incomplete_count"],
                    },
                }
            )
            report["status"] = "blocked"
        if post_audits["metric_regression"]["status"] == "failed":
            append_blocker_once(
                blockers,
                {
                    "code": "metric_regression_failed",
                    "message": "Release metric regression audit detected a deterioration against the baseline.",
                    "evidence": {"regression_count": post_audits["metric_regression"]["regression_count"]},
                }
            )
            report["status"] = "blocked"
        if int(post_audits["parameter_optimization"]["blocking_required_action_count"] or 0) > 0:
            append_blocker_once(
                blockers,
                {
                    "code": "parameter_optimization_action_required",
                    "message": "Release parameter optimization audit has blocking required actions.",
                    "evidence": {
                        "status": post_audits["parameter_optimization"]["status"],
                        "required_action_count": post_audits["parameter_optimization"]["required_action_count"],
                        "blocking_required_action_count": post_audits["parameter_optimization"][
                            "blocking_required_action_count"
                        ],
                        "blocking_action_ids": post_audits["parameter_optimization"]["blocking_action_ids"],
                    },
                },
            )
            report["status"] = "blocked"
        report["blockers"] = blockers
        report["blocker_count"] = len(blockers)
        write_json(workflow_json, report)
        workflow_md.write_text(workflow_markdown(report), encoding="utf-8")

        # Refresh terminal audits against the final workflow report so their
        # evidence reflects blockers introduced by post-audit gates themselves.
        objective_result = run_command(build_objective_command(args, release_id), workspace)
        objective_report = read_json(objective_output)
        post_audits["objective_status"] = {
            "command": objective_result,
            "status": objective_report.get("status", ""),
            "satisfied_count": objective_report.get("satisfied_count", 0),
            "incomplete_count": objective_report.get("incomplete_count", 0),
        }
        metric_result = run_command(build_metric_regression_command(args, release_id), workspace)
        metric_report = read_json(metric_regression_output)
        parameter_result = run_command(build_parameter_optimization_command(args, release_id), workspace)
        parameter_report = read_json(parameter_optimization_output)
        post_audits["metric_regression"] = {
            "command": metric_result,
            "status": (metric_report.get("comparison") or {}).get("status", ""),
            "regression_count": (metric_report.get("comparison") or {}).get("regression_count", 0),
        }
        post_audits["parameter_optimization"] = {
            "command": parameter_result,
            **parameter_optimization_summary(parameter_report),
        }
        report["post_release_audits"] = post_audits
        report["gate_summary"]["objective_status"] = post_audits["objective_status"]["status"]
        report["gate_summary"]["metric_regression_status"] = post_audits["metric_regression"]["status"]
        report["gate_summary"]["metric_regression_count"] = post_audits["metric_regression"]["regression_count"]
        report["gate_summary"]["parameter_optimization_status"] = post_audits["parameter_optimization"]["status"]
        report["gate_summary"]["parameter_optimization_required_action_count"] = post_audits["parameter_optimization"]["required_action_count"]
        report["gate_summary"]["parameter_optimization_blocking_required_action_count"] = post_audits["parameter_optimization"][
            "blocking_required_action_count"
        ]
        if post_audits["objective_status"]["status"] != "complete":
            append_blocker_once(
                blockers,
                {
                    "code": "objective_status_incomplete",
                    "message": "Optimization objective status audit is not complete.",
                    "evidence": {
                        "status": post_audits["objective_status"]["status"],
                        "incomplete_count": post_audits["objective_status"]["incomplete_count"],
                    },
                },
            )
            report["status"] = "blocked"
        if post_audits["metric_regression"]["status"] == "failed":
            append_blocker_once(
                blockers,
                {
                    "code": "metric_regression_failed",
                    "message": "Release metric regression audit detected a deterioration against the baseline.",
                    "evidence": {"regression_count": post_audits["metric_regression"]["regression_count"]},
                },
            )
            report["status"] = "blocked"
        if int(post_audits["parameter_optimization"]["blocking_required_action_count"] or 0) > 0:
            append_blocker_once(
                blockers,
                {
                    "code": "parameter_optimization_action_required",
                    "message": "Release parameter optimization audit has blocking required actions.",
                    "evidence": {
                        "status": post_audits["parameter_optimization"]["status"],
                        "required_action_count": post_audits["parameter_optimization"]["required_action_count"],
                        "blocking_required_action_count": post_audits["parameter_optimization"][
                            "blocking_required_action_count"
                        ],
                        "blocking_action_ids": post_audits["parameter_optimization"]["blocking_action_ids"],
                    },
                },
            )
            report["status"] = "blocked"
        report["blockers"] = blockers
        report["blocker_count"] = len(blockers)
        write_json(workflow_json, report)
        workflow_md.write_text(workflow_markdown(report), encoding="utf-8")
    return report


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run fixture-first strict release workflow.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", default="")
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--fixture-registry", default=DEFAULT_FIXTURE_REGISTRY)
    parser.add_argument("--cscc-csv", default="")
    parser.add_argument("--heldout-top-n", type=int, default=200)
    parser.add_argument("--max-paths", type=int, default=10)
    parser.add_argument("--max-hops", type=int, default=3)
    parser.add_argument("--scenario-timeout-seconds", type=int, default=300)
    parser.add_argument("--heldout-rerun-timeout-seconds", type=int, default=900)
    parser.add_argument("--continue-on-fixture-failure", action="store_true")
    parser.add_argument("--metric-baseline", default="")
    parser.add_argument("--skip-post-audits", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    report = run_workflow(parse_args(argv))
    print(
        json.dumps(
            {
                "status": report["status"],
                "workflow_report": report["outputs"]["workflow_json"],
                "fixture_audit_status": report["gate_summary"]["fixture_audit_status"],
                "batch_benchmark_executed": report["gate_summary"]["batch_benchmark_executed"],
                "readiness_audit_status": report["gate_summary"]["readiness_audit_status"],
                "objective_status": report["gate_summary"].get("objective_status", ""),
                "metric_regression_status": report["gate_summary"].get("metric_regression_status", ""),
                "blocker_count": report["blocker_count"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0 if report["status"] == "ready" else 1


if __name__ == "__main__":
    raise SystemExit(main())
