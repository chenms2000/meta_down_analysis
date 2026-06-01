"""Create and compare release-level metric snapshots.

The release reports are rich but nested. This script extracts stable metrics
that can be compared across optimization cycles to catch regressions when a
change improves one topic while damaging another.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AUDIT_VERSION = "release_metric_regression.v1"
DEFAULT_BATCH_REPORT = "validation_reports/mvp_20260513T002254/batch_benchmark/release_batch_benchmark_report.json"
DEFAULT_READINESS_REPORT = "validation_reports/mvp_20260513T002254/batch_benchmark/release_readiness_audit.json"
DEFAULT_OBJECTIVE_REPORT = "validation_reports/mvp_20260513T002254/optimization_objective_status.json"


METRIC_DIRECTIONS = {
    "batch_failed_case_count": "lower",
    "batch_missing_case_count": "lower",
    "benchmark_matrix_registered_gap_count": "lower",
    "benchmark_matrix_executed_gap_count": "lower",
    "scenario_traceability_failed_case_count": "lower",
    "scenario_traceability_passed": "higher",
    "gold_positive_count": "higher",
    "gold_negative_control_count": "higher",
    "gold_context_specific_positive_count": "higher",
    "gold_unique_cancer_context_count": "higher",
    "gold_unique_mechanism_axis_count": "higher",
    "gold_required_cell_count": "higher",
    "gold_covered_cell_count": "higher",
    "gold_missing_cell_count": "lower",
    "gold_coverage_fraction": "higher",
    "gold_specificity_over_recommended_count": "lower",
    "gold_specificity_multi_cell_match_count": "lower",
    "gold_specificity_single_cell_fraction": "higher",
    "gold_specificity_max_cells_per_gold_id": "lower",
    "literature_case_count": "higher",
    "literature_manual_case_count": "higher",
    "literature_generated_case_count": "higher",
    "literature_generated_group_count": "higher",
    "literature_generated_failed_case_count": "lower",
    "literature_status_accuracy": "higher",
    "literature_cue_recall": "higher",
    "literature_method_cue_recall": "higher",
    "literature_required_coverage_missing_count": "lower",
    "heldout_failed_check_count": "lower",
    "heldout_warning_count": "lower",
    "heldout_random_auc_warning_count": "lower",
    "readiness_blocker_count": "lower",
    "objective_satisfied_count": "higher",
    "objective_incomplete_count": "lower",
    "scenario_negative_trap_hits": "lower",
    "scenario_unsupported_top_claim_rate_max": "lower",
    "scenario_precision_at_1_min": "higher",
    "scenario_precision_at_3_min": "higher",
    "scenario_precision_at_5_min": "higher",
    "scenario_recall_at_1_min": "higher",
    "scenario_recall_at_3_min": "higher",
    "scenario_recall_at_5_min": "higher",
    "scenario_expected_precision_at_1_min": "higher",
    "scenario_expected_precision_at_3_min": "higher",
    "scenario_expected_precision_at_5_min": "higher",
    "scenario_expected_recall_at_1_min": "higher",
    "scenario_expected_recall_at_3_min": "higher",
    "scenario_expected_recall_at_5_min": "higher",
    "scenario_expected_topk_eligible_count": "higher",
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


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def bool_metric(value: Any) -> int:
    return 1 if bool(value) else 0


def max_scenario_metric(batch: dict[str, Any], key: str) -> float:
    values = []
    for case in batch.get("scenarios") or []:
        metrics = ((case.get("conclusion_evaluation") or {}).get("metrics") or {})
        if key in metrics:
            values.append(safe_float(metrics.get(key)))
    return max(values) if values else 0.0


def sum_scenario_metric(batch: dict[str, Any], key: str) -> float:
    total = 0.0
    for case in batch.get("scenarios") or []:
        metrics = ((case.get("conclusion_evaluation") or {}).get("metrics") or {})
        total += safe_float(metrics.get(key))
    return total


def min_scenario_metric(batch: dict[str, Any], key: str) -> float:
    values = []
    for case in batch.get("scenarios") or []:
        metrics = ((case.get("conclusion_evaluation") or {}).get("metrics") or {})
        if key in metrics and metrics.get(key) is not None:
            values.append(safe_float(metrics.get(key)))
    return min(values) if values else 0.0


def min_expected_scenario_metric(batch: dict[str, Any], key: str) -> float:
    values = []
    for case in batch.get("scenarios") or []:
        metrics = ((case.get("conclusion_evaluation") or {}).get("metrics") or {})
        if not metrics.get("expected_topk_eligible"):
            continue
        if key in metrics and metrics.get(key) is not None:
            values.append(safe_float(metrics.get(key)))
    return min(values) if values else 0.0


def expected_topk_eligible_count(batch: dict[str, Any]) -> float:
    return float(
        sum(
            1
            for case in batch.get("scenarios") or []
            if (((case.get("conclusion_evaluation") or {}).get("metrics") or {}).get("expected_topk_eligible"))
        )
    )


def extract_metrics(batch: dict[str, Any], readiness: dict[str, Any], objective: dict[str, Any]) -> dict[str, float]:
    status = batch.get("status") or {}
    matrix = batch.get("benchmark_matrix") or {}
    traceability = batch.get("scenario_traceability_contract") or {}
    gold_quality = batch.get("gold_quality") or {}
    gold_coverage = batch.get("gold_coverage") or {}
    gold_specificity = gold_coverage.get("coverage_specificity") or {}
    gold_distribution = gold_quality.get("positive_distribution") or {}
    literature_metrics = (batch.get("literature_semantics") or {}).get("metrics") or {}
    heldout = batch.get("heldout_validation") or {}
    return {
        "batch_passed": bool_metric(status.get("passed")),
        "batch_failed_case_count": safe_float(status.get("failed_case_count")),
        "batch_missing_case_count": safe_float(status.get("missing_case_count")),
        "benchmark_matrix_passed": bool_metric(status.get("benchmark_matrix_passed")),
        "benchmark_matrix_registered_gap_count": safe_float(matrix.get("registered_gap_count")),
        "benchmark_matrix_executed_gap_count": safe_float(matrix.get("executed_gap_count")),
        "scenario_traceability_passed": bool_metric(status.get("traceability_contract_passed")),
        "scenario_traceability_failed_case_count": safe_float(traceability.get("failed_case_count")),
        "gold_quality_passed": bool_metric(gold_quality.get("status") == "passed"),
        "gold_positive_count": safe_float(gold_quality.get("positive_count")),
        "gold_negative_control_count": safe_float(gold_quality.get("negative_control_count")),
        "gold_context_specific_positive_count": safe_float(gold_distribution.get("context_specific_positive_count")),
        "gold_unique_cancer_context_count": safe_float(gold_distribution.get("unique_positive_cancer_context_count")),
        "gold_unique_mechanism_axis_count": safe_float(gold_distribution.get("unique_positive_mechanism_axis_count")),
        "gold_required_cell_count": safe_float(gold_coverage.get("required_cell_count")),
        "gold_covered_cell_count": safe_float(gold_coverage.get("covered_cell_count")),
        "gold_missing_cell_count": safe_float(gold_coverage.get("missing_cell_count")),
        "gold_coverage_fraction": safe_float(gold_coverage.get("coverage_fraction")),
        "gold_specificity_over_recommended_count": safe_float(gold_specificity.get("over_recommended_count")),
        "gold_specificity_multi_cell_match_count": safe_float(gold_specificity.get("multi_cell_gold_match_count")),
        "gold_specificity_single_cell_fraction": safe_float(gold_specificity.get("single_cell_gold_fraction")),
        "gold_specificity_max_cells_per_gold_id": safe_float(gold_specificity.get("max_cells_per_gold_id_observed")),
        "literature_passed": bool_metric((batch.get("literature_semantics") or {}).get("status") == "passed"),
        "literature_case_count": safe_float(literature_metrics.get("case_count")),
        "literature_manual_case_count": safe_float(literature_metrics.get("manual_case_count")),
        "literature_generated_case_count": safe_float(literature_metrics.get("generated_case_count")),
        "literature_generated_group_count": safe_float(literature_metrics.get("generated_group_count")),
        "literature_generated_failed_case_count": safe_float(literature_metrics.get("generated_failed_case_count")),
        "literature_status_accuracy": safe_float(literature_metrics.get("status_accuracy")),
        "literature_cue_recall": safe_float(literature_metrics.get("cue_recall")),
        "literature_method_cue_recall": safe_float(literature_metrics.get("method_cue_recall")),
        "literature_required_coverage_missing_count": safe_float(literature_metrics.get("required_coverage_missing_count")),
        "heldout_passed": bool_metric(heldout.get("status") == "passed"),
        "heldout_rerun_performed": bool_metric(heldout.get("rerun_performed")),
        "heldout_failed_check_count": float(len(heldout.get("failed_checks") or [])),
        "heldout_warning_count": float(len(heldout.get("warnings") or [])),
        "heldout_random_auc_warning_count": float(
            sum(1 for row in heldout.get("warnings") or [] if isinstance(row, dict) and row.get("code") == "source_heldout_random_auc")
        ),
        "readiness_ready": bool_metric(readiness.get("status") == "ready"),
        "readiness_blocker_count": safe_float(readiness.get("blocker_count")),
        "objective_complete": bool_metric(objective.get("status") == "complete"),
        "objective_satisfied_count": safe_float(objective.get("satisfied_count")),
        "objective_incomplete_count": safe_float(objective.get("incomplete_count")),
        "scenario_negative_trap_hits": sum_scenario_metric(batch, "negative_trap_hits"),
        "scenario_unsupported_top_claim_rate_max": max_scenario_metric(batch, "unsupported_top_claim_rate"),
        "scenario_precision_at_1_min": min_scenario_metric(batch, "precision_at_1"),
        "scenario_precision_at_3_min": min_scenario_metric(batch, "precision_at_3"),
        "scenario_precision_at_5_min": min_scenario_metric(batch, "precision_at_5"),
        "scenario_recall_at_1_min": min_scenario_metric(batch, "recall_at_1"),
        "scenario_recall_at_3_min": min_scenario_metric(batch, "recall_at_3"),
        "scenario_recall_at_5_min": min_scenario_metric(batch, "recall_at_5"),
        "scenario_expected_topk_eligible_count": expected_topk_eligible_count(batch),
        "scenario_expected_precision_at_1_min": min_expected_scenario_metric(batch, "expected_precision_at_1"),
        "scenario_expected_precision_at_3_min": min_expected_scenario_metric(batch, "expected_precision_at_3"),
        "scenario_expected_precision_at_5_min": min_expected_scenario_metric(batch, "expected_precision_at_5"),
        "scenario_expected_recall_at_1_min": min_expected_scenario_metric(batch, "expected_recall_at_1"),
        "scenario_expected_recall_at_3_min": min_expected_scenario_metric(batch, "expected_recall_at_3"),
        "scenario_expected_recall_at_5_min": min_expected_scenario_metric(batch, "expected_recall_at_5"),
    }


def snapshot(batch: dict[str, Any], readiness: dict[str, Any], objective: dict[str, Any]) -> dict[str, Any]:
    metrics = extract_metrics(batch, readiness, objective)
    return {
        "snapshot_version": AUDIT_VERSION,
        "created_at_utc": utc_now(),
        "metrics": metrics,
        "metric_directions": {key: METRIC_DIRECTIONS.get(key, "equal_or_higher") for key in metrics},
        "source_status": {
            "batch_status": (batch.get("status") or {}).get("passed"),
            "readiness_status": readiness.get("status", ""),
            "objective_status": objective.get("status", ""),
        },
    }


def compare_snapshots(current: dict[str, Any], baseline: dict[str, Any], tolerance: float = 0.0) -> dict[str, Any]:
    if not baseline:
        return {
            "status": "snapshot_only",
            "regression_count": 0,
            "regressions": [],
            "interpretation": "No baseline snapshot was supplied; current metrics were recorded for future comparison.",
        }
    current_metrics = current.get("metrics") or {}
    baseline_metrics = baseline.get("metrics") or {}
    regressions = []
    for name, current_value in sorted(current_metrics.items()):
        if name not in baseline_metrics:
            continue
        baseline_value = safe_float(baseline_metrics.get(name))
        observed = safe_float(current_value)
        direction = METRIC_DIRECTIONS.get(name, "higher")
        regressed = False
        if direction == "lower":
            regressed = observed > baseline_value + tolerance
        else:
            regressed = observed < baseline_value - tolerance
        if regressed:
            regressions.append(
                {
                    "metric": name,
                    "direction": direction,
                    "baseline": baseline_value,
                    "current": observed,
                    "delta": round(observed - baseline_value, 6),
                }
            )
    return {
        "status": "passed" if not regressions else "failed",
        "regression_count": len(regressions),
        "regressions": regressions,
        "interpretation": "Metric regressions are directional: lower-is-better counts must not increase; higher-is-better coverage/accuracy metrics must not decrease.",
    }


def markdown_summary(report: dict[str, Any]) -> str:
    comparison = report.get("comparison") or {}
    lines = [
        "# Release Metric Regression",
        "",
        f"Status: `{comparison.get('status', '')}`",
        f"Regressions: `{comparison.get('regression_count', 0)}`",
        "",
        "## Metrics",
        "",
    ]
    for key, value in sorted((report.get("current_snapshot") or {}).get("metrics", {}).items()):
        lines.append(f"- `{key}`: `{value}`")
    regressions = comparison.get("regressions") or []
    if regressions:
        lines.extend(["", "## Regressions", ""])
        for row in regressions:
            lines.append(
                f"- `{row.get('metric')}`: baseline `{row.get('baseline')}`, current `{row.get('current')}`, direction `{row.get('direction')}`"
            )
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create and compare release metric snapshots.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--batch-report", default=DEFAULT_BATCH_REPORT)
    parser.add_argument("--readiness-report", default=DEFAULT_READINESS_REPORT)
    parser.add_argument("--objective-report", default=DEFAULT_OBJECTIVE_REPORT)
    parser.add_argument("--baseline", default="")
    parser.add_argument("--output", default="")
    parser.add_argument("--snapshot-output", default="")
    parser.add_argument("--markdown-output", default="")
    parser.add_argument("--tolerance", type=float, default=0.0)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    workspace = Path(args.workspace).resolve()
    current = snapshot(
        read_json(workspace / args.batch_report),
        read_json(workspace / args.readiness_report),
        read_json(workspace / args.objective_report),
    )
    baseline = read_json(workspace / args.baseline) if args.baseline else {}
    comparison = compare_snapshots(current, baseline, args.tolerance)
    report = {
        "audit_version": AUDIT_VERSION,
        "created_at_utc": utc_now(),
        "current_snapshot": current,
        "baseline_path": str(workspace / args.baseline) if args.baseline else "",
        "comparison": comparison,
    }
    if args.snapshot_output:
        write_json(workspace / args.snapshot_output, current)
    if args.output:
        write_json(workspace / args.output, report)
    if args.markdown_output:
        path = workspace / args.markdown_output
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(markdown_summary(report), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if comparison["status"] in {"passed", "snapshot_only"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
