"""Run release-level checks for structured literature assertion semantics."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from build_literature_evidence import assertion_semantics  # noqa: E402


BENCHMARK_VERSION = "literature_semantics_benchmark.v1"
DEFAULT_BENCHMARK_PATH = "config/literature_semantics_benchmark.json"
DEFAULT_OUTPUT_ROOT = "validation_reports"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def expected_list(row: dict[str, Any], key: str) -> list[str]:
    value = row.get(key) or []
    if isinstance(value, list):
        return [str(item) for item in value if str(item)]
    return [str(value)] if str(value) else []


def required_list(benchmark: dict[str, Any], key: str) -> list[str]:
    return sorted(set(expected_list(benchmark, key)))


def generated_cases(benchmark: dict[str, Any]) -> list[dict[str, Any]]:
    """Expand systematic semantic stress groups into benchmark cases."""
    expanded: list[dict[str, Any]] = []
    groups = [row for row in benchmark.get("generated_case_groups", []) or [] if isinstance(row, dict)]
    for group in groups:
        group_id = str(group.get("group_id") or "").strip()
        if not group_id:
            continue
        texts = group.get("texts") or []
        if not isinstance(texts, list):
            continue
        for index, item in enumerate(texts, start=1):
            if isinstance(item, dict):
                text = str(item.get("text") or "")
                section = str(item.get("section") or group.get("section") or "")
                expected_status = str(item.get("expected_support_status") or group.get("expected_support_status") or "")
                expected_cues = expected_list(item, "expected_cues") or expected_list(group, "expected_cues")
                expected_method_cues = expected_list(item, "expected_method_cues") or expected_list(group, "expected_method_cues")
                forbidden_cues = expected_list(item, "forbidden_cues") or expected_list(group, "forbidden_cues")
            else:
                text = str(item)
                section = str(group.get("section") or "")
                expected_status = str(group.get("expected_support_status") or "")
                expected_cues = expected_list(group, "expected_cues")
                expected_method_cues = expected_list(group, "expected_method_cues")
                forbidden_cues = expected_list(group, "forbidden_cues")
            if not text.strip():
                continue
            case = {
                "case_id": f"{group_id}__{index:02d}",
                "text": text,
                "section": section,
                "expected_support_status": expected_status,
                "expected_cues": expected_cues,
                "expected_method_cues": expected_method_cues,
                "forbidden_cues": forbidden_cues,
                "generated_from_group": group_id,
            }
            expanded.append(case)
    return expanded


def evaluate_case(row: dict[str, Any]) -> dict[str, Any]:
    observed = assertion_semantics(str(row.get("text") or ""), str(row.get("section") or ""))
    expected_status = str(row.get("expected_support_status") or "")
    expected_cues = set(expected_list(row, "expected_cues"))
    expected_method_cues = set(expected_list(row, "expected_method_cues"))
    forbidden_cues = set(expected_list(row, "forbidden_cues"))
    observed_cues = set(observed.get("semantic_cues") or [])
    observed_method_cues = set(observed.get("method_cues") or [])
    missing_cues = sorted(expected_cues - observed_cues)
    missing_method_cues = sorted(expected_method_cues - observed_method_cues)
    forbidden_hits = sorted(forbidden_cues & observed_cues)
    status_ok = observed.get("support_status") == expected_status
    cues_ok = not missing_cues and not forbidden_hits
    method_cues_ok = not missing_method_cues
    return {
        "case_id": str(row.get("case_id") or ""),
        "section": str(row.get("section") or ""),
        "expected_support_status": expected_status,
        "observed_support_status": observed.get("support_status"),
        "expected_cues": sorted(expected_cues),
        "observed_cues": sorted(observed_cues),
        "missing_cues": missing_cues,
        "forbidden_cue_hits": forbidden_hits,
        "expected_method_cues": sorted(expected_method_cues),
        "observed_method_cues": sorted(observed_method_cues),
        "missing_method_cues": missing_method_cues,
        "decision_reason": observed.get("decision_reason", ""),
        "generated_from_group": str(row.get("generated_from_group") or ""),
        "status_ok": status_ok,
        "cues_ok": cues_ok,
        "method_cues_ok": method_cues_ok,
        "passed": status_ok and cues_ok and method_cues_ok,
    }


def safe_fraction(numerator: int, denominator: int) -> float | None:
    if denominator <= 0:
        return None
    return round(numerator / denominator, 6)


def coverage_report(required: list[str], covered: set[str]) -> dict[str, Any]:
    required_set = set(required)
    return {
        "required": sorted(required_set),
        "covered": sorted(required_set & covered),
        "missing": sorted(required_set - covered),
        "status": "covered" if not (required_set - covered) else "gaps_present",
    }


def benchmark_literature_semantics(benchmark: dict[str, Any]) -> dict[str, Any]:
    manual_cases = [row for row in benchmark.get("cases", []) or [] if isinstance(row, dict)]
    generated = generated_cases(benchmark)
    cases = manual_cases + generated
    results = [evaluate_case(row) for row in cases]
    expected_cue_total = sum(len(row["expected_cues"]) for row in results)
    missing_cue_total = sum(len(row["missing_cues"]) + len(row["forbidden_cue_hits"]) for row in results)
    expected_method_total = sum(len(row["expected_method_cues"]) for row in results)
    missing_method_total = sum(len(row["missing_method_cues"]) for row in results)
    status_accuracy = safe_fraction(sum(1 for row in results if row["status_ok"]), len(results))
    cue_recall = safe_fraction(expected_cue_total - missing_cue_total, expected_cue_total)
    method_cue_recall = safe_fraction(expected_method_total - missing_method_total, expected_method_total)
    min_status = float(benchmark.get("minimum_status_accuracy", 1.0))
    min_cue = float(benchmark.get("minimum_cue_recall", 0.95))
    min_method = float(benchmark.get("minimum_method_cue_recall", 0.90))
    min_cases = int(benchmark.get("minimum_case_count") or 1)
    min_generated_cases = int(benchmark.get("minimum_generated_case_count") or 0)
    required_generated_groups = required_list(benchmark, "required_generated_groups")
    generated_groups = {
        str(row.get("generated_from_group") or "")
        for row in results
        if str(row.get("generated_from_group") or "")
    }
    support_status_coverage = coverage_report(
        required_list(benchmark, "required_support_statuses"),
        {str(row["expected_support_status"]) for row in results if str(row["expected_support_status"])},
    )
    semantic_cue_coverage = coverage_report(
        required_list(benchmark, "required_semantic_cues"),
        {cue for row in results for cue in row["expected_cues"]},
    )
    method_cue_coverage = coverage_report(
        required_list(benchmark, "required_method_cues"),
        {cue for row in results for cue in row["expected_method_cues"]},
    )
    coverage = {
        "support_statuses": support_status_coverage,
        "semantic_cues": semantic_cue_coverage,
        "method_cues": method_cue_coverage,
        "generated_groups": coverage_report(required_generated_groups, generated_groups),
    }
    missing_required_coverage = (
        support_status_coverage["missing"]
        + semantic_cue_coverage["missing"]
        + method_cue_coverage["missing"]
        + coverage["generated_groups"]["missing"]
    )
    generated_failed_count = sum(1 for row in results if row["generated_from_group"] and not row["passed"])
    metrics = {
        "case_count": len(results),
        "manual_case_count": len(manual_cases),
        "generated_case_count": len(generated),
        "generated_group_count": len(generated_groups),
        "generated_failed_case_count": generated_failed_count,
        "passed_case_count": sum(1 for row in results if row["passed"]),
        "failed_case_count": sum(1 for row in results if not row["passed"]),
        "status_accuracy": status_accuracy,
        "cue_recall": cue_recall,
        "method_cue_recall": method_cue_recall,
        "minimum_status_accuracy": min_status,
        "minimum_cue_recall": min_cue,
        "minimum_method_cue_recall": min_method,
        "minimum_case_count": min_cases,
        "minimum_generated_case_count": min_generated_cases,
        "required_coverage_missing_count": len(missing_required_coverage),
    }
    passed = (
        bool(results)
        and len(results) >= min_cases
        and len(generated) >= min_generated_cases
        and status_accuracy is not None
        and cue_recall is not None
        and method_cue_recall is not None
        and status_accuracy >= min_status
        and cue_recall >= min_cue
        and method_cue_recall >= min_method
        and not [row for row in results if not row["passed"]]
        and not missing_required_coverage
    )
    return {
        "benchmark_version": BENCHMARK_VERSION,
        "created_at_utc": utc_now(),
        "schema_version": benchmark.get("schema_version", ""),
        "metrics": metrics,
        "coverage": coverage,
        "failed_cases": [row for row in results if not row["passed"]],
        "cases": results,
        "status": "passed" if passed else "failed",
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run literature assertion semantics benchmark.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--benchmark-path", default=DEFAULT_BENCHMARK_PATH)
    parser.add_argument("--output", default="")
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--release-id", default="manual")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    workspace = Path(args.workspace).resolve()
    report = benchmark_literature_semantics(read_json(workspace / args.benchmark_path))
    output = args.output
    if not output:
        output = str(Path(args.output_root) / args.release_id / "literature_semantics_benchmark_report.json")
    output_path = workspace / output
    write_json(output_path, report)
    report["output"] = str(output_path)
    print(json.dumps({key: report[key] for key in ("status", "metrics", "output")}, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
