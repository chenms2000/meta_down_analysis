"""D5 external LLM pilot gate for the safe adapter.

This script wraps the D4 regression runner with explicit external-LLM readiness
and promotion criteria. It never calls an external model unless --run-external
is supplied and the external backend is explicitly enabled/configured.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from llm_safe_adapter import ExternalLLMConfig, content_hash, guard_policy_hash  # noqa: E402
from run_llm_safe_adapter_regression import (  # noqa: E402
    DEFAULT_FIXTURE_PATH,
    DEFAULT_OUTPUT_ROOT,
    build_report as build_regression_report,
    build_service,
    read_json,
)


PILOT_VERSION = "llm_safe_adapter.external_pilot.20260514"
REPORT_NAME = "llm_safe_adapter_external_pilot_report.json"


def readiness(llm_config: ExternalLLMConfig, run_external: bool) -> dict[str, Any]:
    missing = []
    if not llm_config.enabled:
        missing.append("external_llm_not_enabled")
    if not llm_config.model:
        missing.append("missing_model")
    if not llm_config.api_key:
        missing.append("missing_api_key")
    if not llm_config.endpoint:
        missing.append("missing_endpoint")
    return {
        "run_external_requested": bool(run_external),
        "ready": bool(run_external) and not missing,
        "missing": missing,
        "config": llm_config.sanitized(),
    }


def external_case_statuses(regression_report: dict[str, Any]) -> list[str]:
    return [case.get("external", {}).get("status", "unknown") for case in regression_report.get("cases", [])]


def external_guard_issue_counts(regression_report: dict[str, Any]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for case in regression_report.get("cases", []):
        counts.update(case.get("external", {}).get("guard_issue_codes", []) or [])
    return dict(sorted(counts.items()))


def score_pilot(regression_report: dict[str, Any], ready: dict[str, Any], min_external_pass_rate: float) -> dict[str, Any]:
    statuses = external_case_statuses(regression_report)
    total = len(statuses)
    passed = sum(1 for status in statuses if status == "passed")
    pass_rate = 0.0 if total == 0 else passed / total
    adversarial_failed = int(regression_report.get("summary", {}).get("adversarial_probe_failed", 0) or 0)
    if not ready.get("run_external_requested"):
        gate_status = "not_run"
        gate_passed = True
    elif not ready.get("ready"):
        gate_status = "not_ready"
        gate_passed = False
    elif adversarial_failed:
        gate_status = "failed_adversarial_regression"
        gate_passed = False
    elif pass_rate >= min_external_pass_rate:
        gate_status = "passed"
        gate_passed = True
    else:
        gate_status = "needs_prompt_iteration"
        gate_passed = False
    return {
        "external_case_count": total,
        "external_passed": passed,
        "external_pass_rate": round(pass_rate, 6),
        "min_external_pass_rate": float(min_external_pass_rate),
        "external_status_counts": dict(sorted(Counter(statuses).items())),
        "external_guard_issue_counts": external_guard_issue_counts(regression_report),
        "adversarial_probe_failed": adversarial_failed,
        "gate_status": gate_status,
        "gate_passed": gate_passed,
    }


def build_pilot_report(args: argparse.Namespace) -> dict[str, Any]:
    fixture_path = Path(args.fixture)
    fixture = read_json(fixture_path)
    service = build_service(args)
    ready = readiness(service.llm_config, run_external=bool(args.run_external))
    regression_report = build_regression_report(
        service,
        fixture,
        fixture_path,
        service.llm_config,
        run_external=bool(args.run_external),
        max_cases=args.max_cases,
    )
    score = score_pilot(regression_report, ready, min_external_pass_rate=args.min_external_pass_rate)
    report = {
        "pilot_version": PILOT_VERSION,
        "release_id": getattr(service, "release_id", fixture.get("release_id", "")),
        "fixture_path": str(fixture_path.as_posix()),
        "readiness": ready,
        "guard_policy_hash": guard_policy_hash(),
        "regression_report_hash": regression_report.get("report_hash", ""),
        "regression_summary": regression_report.get("summary", {}),
        "score": score,
        "cases": [
            {
                "case_id": case.get("case_id", ""),
                "local_status": case.get("local", {}).get("status", ""),
                "external_status": case.get("external", {}).get("status", ""),
                "external_guard_issue_codes": case.get("external", {}).get("guard_issue_codes", []),
                "adapter_input_hash": case.get("input_summary", {}).get("adapter_input_hash", ""),
            }
            for case in regression_report.get("cases", [])
        ],
    }
    report["report_hash"] = content_hash(report)
    return report


def write_report(report: dict[str, Any], output_root: Path) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    path = output_root / REPORT_NAME
    path.write_text(json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2), encoding="utf-8")
    return path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run D5 external LLM pilot/readiness checks for the safe adapter.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--normalized-root", default="normalized_store")
    parser.add_argument("--graph-root", default="graph_projection")
    parser.add_argument("--pubchem-root", default="pubchem_cid_cache")
    parser.add_argument("--compound-root", default="compound_match_index")
    parser.add_argument("--literature-root", default="literature_evidence")
    parser.add_argument("--release-id", default="")
    parser.add_argument("--fixture", default=DEFAULT_FIXTURE_PATH)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--write-report", action="store_true")
    parser.add_argument("--max-cases", type=int, default=None)
    parser.add_argument("--run-external", action="store_true")
    parser.add_argument("--enable-external-llm", action="store_true")
    parser.add_argument("--llm-provider", default="")
    parser.add_argument("--llm-endpoint", default="")
    parser.add_argument("--llm-model", default="")
    parser.add_argument("--llm-api-key-env", default="")
    parser.add_argument("--llm-timeout-seconds", type=float, default=None)
    parser.add_argument("--llm-max-output-tokens", type=int, default=None)
    parser.add_argument("--min-external-pass-rate", type=float, default=1.0)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    report = build_pilot_report(args)
    if args.write_report:
        write_report(report, Path(args.output_root))
    sys.stdout.write(json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2) + "\n")
    return 0 if report.get("score", {}).get("gate_passed") else 2


if __name__ == "__main__":
    raise SystemExit(main())
