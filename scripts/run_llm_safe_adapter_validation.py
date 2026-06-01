"""Validation runner for the local LLM safe adapter layer."""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from llm_safe_adapter import (  # noqa: E402
    ADAPTER_VERSION,
    build_local_adapter_output,
    content_hash,
    guard_adapter_output,
)


RUNNER_VERSION = "llm_safe_adapter.validation.20260513"
DEFAULT_FIXTURE_DIR = "tests/fixtures/llm_safe_adapter"
DEFAULT_OUTPUT_ROOT = "validation_reports"
REPORT_NAME = "llm_safe_adapter_validation_report.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_fixtures(fixture_dir: Path) -> list[dict[str, Any]]:
    fixtures = []
    for path in sorted(fixture_dir.glob("*.json")):
        payload = read_json(path)
        if payload.get("fixture_version") != "llm_safe_adapter.20260513":
            continue
        payload["_fixture_file"] = path.name
        fixtures.append(payload)
    return fixtures


def build_check(name: str, passed: bool, detail: dict[str, Any] | None = None) -> dict[str, Any]:
    return {"name": name, "passed": bool(passed), "detail": detail or {}}


def issue_codes(guard: dict[str, Any]) -> set[str]:
    return {str(issue.get("code") or "") for issue in guard.get("issues", [])}


def validate_case(fixture: dict[str, Any], case: dict[str, Any]) -> dict[str, Any]:
    input_pack = case.get("input_pack", {})
    expected = case.get("expected", {})
    if "adapter_output" in case:
        output = copy.deepcopy(case["adapter_output"])
        guard = guard_adapter_output(input_pack, output)
    else:
        output = build_local_adapter_output(input_pack)
        guard = output.get("guard", {})

    codes = issue_codes(guard)
    required_issue_codes = set(expected.get("required_issue_codes") or [])
    required_outputs = set(expected.get("required_outputs") or [])
    checks = [
        build_check(
            "guard_expectation",
            bool(guard.get("passed")) == bool(expected.get("guard_passed", True)),
            {"actual": bool(guard.get("passed")), "expected": bool(expected.get("guard_passed", True)), "issue_codes": sorted(codes)},
        ),
        build_check(
            "status_expectation",
            output.get("status") == expected.get("status", output.get("status")),
            {"actual": output.get("status"), "expected": expected.get("status", output.get("status"))},
        ),
        build_check(
            "required_issue_codes",
            required_issue_codes.issubset(codes),
            {"actual": sorted(codes), "required": sorted(required_issue_codes)},
        ),
        build_check(
            "required_output_fields",
            required_outputs.issubset(set(output)),
            {"actual": sorted(set(output)), "required": sorted(required_outputs)},
        ),
    ]
    if expected.get("requires_source_bound_language"):
        checks.append(
            build_check(
                "source_bound_language",
                "language_without_source_refs" not in codes and int(guard.get("language_block_count", 0) or 0) > 0,
                {"language_block_count": guard.get("language_block_count", 0), "issue_codes": sorted(codes)},
            )
        )

    passed = all(check["passed"] for check in checks)
    return {
        "fixture": fixture.get("_fixture_file", ""),
        "case_id": case.get("case_id", ""),
        "status": "passed" if passed else "failed",
        "adapter_output_status": output.get("status", ""),
        "guard_passed": bool(guard.get("passed")),
        "guard_issue_codes": sorted(codes),
        "output_hash": content_hash(output),
        "checks": checks,
    }


def report_summary(case_results: list[dict[str, Any]]) -> dict[str, Any]:
    passed = sum(1 for row in case_results if row.get("status") == "passed")
    failed = sum(1 for row in case_results if row.get("status") == "failed")
    return {
        "passed": passed,
        "failed": failed,
        "total": len(case_results),
        "gate_passed": failed == 0,
        "gate_status": "passed" if failed == 0 else "failed",
    }


def build_report(fixtures: list[dict[str, Any]], fixture_dir: Path) -> dict[str, Any]:
    cases = []
    for fixture in fixtures:
        for case in fixture.get("cases", []):
            cases.append(validate_case(fixture, case))
    report = {
        "runner_version": RUNNER_VERSION,
        "adapter_version": ADAPTER_VERSION,
        "fixture_dir": str(fixture_dir.as_posix()),
        "fixture_files": [fixture.get("_fixture_file", "") for fixture in fixtures],
        "summary": report_summary(cases),
        "cases": cases,
    }
    report["report_hash"] = content_hash(report)
    return report


def write_report(report: dict[str, Any], output_root: Path) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    path = output_root / REPORT_NAME
    path.write_text(json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2), encoding="utf-8")
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run LLM safe adapter validation fixtures.")
    parser.add_argument("--fixture-dir", default=DEFAULT_FIXTURE_DIR)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--write-report", action="store_true")
    args = parser.parse_args(argv)

    fixture_dir = Path(args.fixture_dir)
    fixtures = load_fixtures(fixture_dir)
    report = build_report(fixtures, fixture_dir)
    if args.write_report:
        write_report(report, Path(args.output_root))
    sys.stdout.write(json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2) + "\n")
    return 0 if report.get("summary", {}).get("gate_passed") else 2


if __name__ == "__main__":
    raise SystemExit(main())
