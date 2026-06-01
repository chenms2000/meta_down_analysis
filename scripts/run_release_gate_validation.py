"""Final release-gate checks for research-prioritization outputs.

This runner validates artifacts that are easy to regress late in the pipeline:
held-out reports, cell overlay provenance, generalized-mode claim strength, drug
confidence caps, context mismatch downgrades, and structured prediction shape.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq


RUNNER_VERSION = "release_gate.validation.20260519"
DEFAULT_OUTPUT_ROOT = "learning_runs"
TIER_RANK = {"low": 0, "exploratory": 1, "medium": 2, "high": 3}
DISEASE_MODEL_TERMS = (
    "cancer",
    "carcinoma",
    "melanoma",
    "glioma",
    "leukemia",
    "lymphoma",
    "tumor",
    "tumour",
    "disease",
    "cell line",
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_check(name: str, passed: bool, details: dict[str, Any] | None = None) -> dict[str, Any]:
    return {"name": name, "passed": bool(passed), "details": details or {}}


def tier_rank(value: Any) -> int:
    return TIER_RANK.get(str(value or "low"), 0)


def structured_from_analysis(path: Path) -> dict[str, Any]:
    payload = read_json(path)
    if "structured_prediction" in payload:
        return payload["structured_prediction"]
    if "analysis_pack" in payload and "structured_prediction" in payload["analysis_pack"]:
        return payload["analysis_pack"]["structured_prediction"]
    return {}


def check_structured_prediction(structured: dict[str, Any], source: str) -> list[dict[str, Any]]:
    if not structured:
        return [build_check("structured_prediction_present", False, {"source": source})]
    checks = [build_check("structured_prediction_present", True, {"source": source})]
    mode = structured.get("mode", "")
    disease_rows = structured.get("disease_predictions") or []
    drug_rows = structured.get("drug_hypotheses") or []
    product_layers = structured.get("product_layers") or {}
    required_layers = [
        "primary_research_candidates",
        "mechanistic_support",
        "appendix_overlay_only",
        "context_mismatch",
        "clinical_warning_or_research_only",
        "excluded_from_primary_reason",
    ]
    missing_layers = [key for key in required_layers if key not in product_layers]
    checks.append(
        build_check(
            "structured_prediction_product_layers_present",
            not missing_layers,
            {"missing_layers": missing_layers},
        )
    )
    primary_rows = product_layers.get("primary_research_candidates") or []
    primary_bad_rows = [
        row.get("display_name", row.get("prediction_id", ""))
        for row in primary_rows
        if row.get("result_type") in {"drug", "disease"}
        or row.get("appendix")
        or row.get("research_only")
        or row.get("excluded_from_primary_reason")
    ]
    checks.append(
        build_check(
            "primary_research_candidates_exclude_drug_disease_appendix",
            not primary_bad_rows,
            {"failure_count": len(primary_bad_rows), "sample": primary_bad_rows[:5]},
        )
    )
    claim_required_groups = [
        "high_confidence_themes",
        "medium_confidence_themes",
        "exploratory_themes",
        "downgraded_but_supported_themes",
        "pathway_evidence",
        "target_predictions",
        "disease_predictions",
        "drug_hypotheses",
    ]
    required_fields = ["confidence_tier", "evidence_refs", "claim_refs", "appendix", "research_only"]
    claim_field_failures = []
    for group in claim_required_groups:
        for index, row in enumerate(structured.get(group, []) or []):
            if not isinstance(row, dict):
                claim_field_failures.append(f"{group}[{index}]:not_object")
                continue
            missing = []
            for field in required_fields:
                if field == "confidence_tier" and ("confidence_tier" in row or "confidence" in row or "display_confidence" in row):
                    continue
                if field not in row:
                    missing.append(field)
            if missing:
                claim_field_failures.append(f"{group}[{index}] missing {','.join(missing)}")
    checks.append(
        build_check(
            "structured_prediction_claim_rows_have_required_bindings",
            not claim_field_failures,
            {"failure_count": len(claim_field_failures), "sample": claim_field_failures[:5]},
        )
    )
    if mode == "generalized":
        strong_rows = [
            row
            for row in [*disease_rows, *drug_rows]
            if tier_rank(row.get("confidence_tier") or row.get("display_confidence")) >= tier_rank("medium")
        ]
        checks.append(
            build_check(
                "generalized_mode_no_strong_disease_or_drug_claims",
                not strong_rows,
                {"strong_row_count": len(strong_rows)},
            )
        )

    theme_rows = [
        row
        for key in (
            "high_confidence_themes",
            "medium_confidence_themes",
            "exploratory_themes",
            "downgraded_but_supported_themes",
        )
        for row in structured.get(key, []) or []
    ]
    support_failures = [
        row.get("display_name", "")
        for row in theme_rows
        if int(row.get("matched_support_count") or 0) > int(row.get("related_input_count") or 0)
    ]
    checks.append(
        build_check(
            "matched_support_count_not_above_related_input_count",
            not support_failures,
            {"failure_count": len(support_failures), "sample": support_failures[:5]},
        )
    )

    context_rows = structured.get("context_mismatch_results") or []
    context_failures = [
        row.get("display_name", row.get("prediction_id", ""))
        for row in context_rows
        if row.get("calibration_status") != "context_mismatch"
        and tier_rank(row.get("confidence_tier")) > tier_rank("exploratory")
    ]
    checks.append(
        build_check(
            "context_mismatch_results_are_downgraded",
            not context_failures,
            {"failure_count": len(context_failures), "sample": context_failures[:5]},
        )
    )

    pathway_failures = []
    for row in structured.get("pathway_evidence", []) or []:
        name = str(row.get("pathway_name") or "").casefold()
        if any(term in name for term in DISEASE_MODEL_TERMS) and not row.get("appendix"):
            pathway_failures.append(row.get("pathway_name", ""))
    checks.append(
        build_check(
            "disease_model_pathways_are_appendix",
            not pathway_failures,
            {"failure_count": len(pathway_failures), "sample": pathway_failures[:5]},
        )
    )

    drug_failures = []
    clinical_failures = []
    exclusion_failures = []
    for row in drug_rows:
        if tier_rank(row.get("confidence_tier") or row.get("display_confidence")) > tier_rank(row.get("upstream_target_confidence")):
            drug_failures.append(row.get("drug_name", ""))
        warning = str(row.get("clinical_warning") or "").casefold()
        if not row.get("research_only") or not row.get("appendix") or "not a treatment recommendation" not in warning:
            clinical_failures.append(row.get("drug_name", ""))
        if not row.get("excluded_from_primary_reason"):
            exclusion_failures.append(row.get("drug_name", ""))
    checks.append(
        build_check(
            "drug_confidence_not_above_target_confidence",
            not drug_failures,
            {"failure_count": len(drug_failures), "sample": drug_failures[:5]},
        )
    )
    checks.append(
        build_check(
            "drug_outputs_remain_research_only",
            not clinical_failures,
            {"failure_count": len(clinical_failures), "sample": clinical_failures[:5]},
        )
    )
    checks.append(
        build_check(
            "drug_outputs_have_primary_exclusion_reason",
            not exclusion_failures,
            {"failure_count": len(exclusion_failures), "sample": exclusion_failures[:5]},
        )
    )
    return checks


def check_cell_overlay(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "views" / "cell_context_overlay.parquet"
    if not path.exists():
        return build_check("cell_overlay_nonempty_with_provenance", False, {"path": str(path), "reason": "missing"})
    frame = pq.read_table(path).to_pandas()
    required = ["source_name", "source_record_id", "evidence_level", "license_id"]
    missing_columns = [column for column in required if column not in frame.columns]
    missing_values = []
    if not frame.empty and not missing_columns:
        for column in required:
            if frame[column].fillna("").astype(str).str.strip().eq("").any():
                missing_values.append(column)
    return build_check(
        "cell_overlay_nonempty_with_provenance",
        not frame.empty and not missing_columns and not missing_values,
        {
            "path": str(path),
            "rows": int(len(frame)),
            "missing_columns": missing_columns,
            "columns_with_blank_values": missing_values,
        },
    )


def check_heldout_reports(run_dir: Path) -> list[dict[str, Any]]:
    temporal = run_dir / "reports" / "temporal_holdout" / "temporal_holdout_validation_report.json"
    source = run_dir / "reports" / "source_heldout" / "source_heldout_validation_report.json"
    checks = []
    if temporal.exists():
        report = read_json(temporal)
        metrics = report.get("validation_metrics") or report.get("metrics") or {}
        checks.append(
            build_check(
                "temporal_heldout_report_metrics_readable",
                "heldout_rows" in metrics and any(key in metrics for key in ("heldout_roc_auc", "heldout_average_precision")),
                {"path": str(temporal), "metric_keys": sorted(metrics)[:20]},
            )
        )
    else:
        checks.append(build_check("temporal_heldout_report_metrics_readable", False, {"path": str(temporal), "reason": "missing"}))
    if source.exists():
        report = read_json(source)
        heldouts = report.get("heldout_sources") or []
        readable = bool(heldouts) and all(
            isinstance(row.get("validation_metrics") or row.get("metrics"), dict) for row in heldouts
        )
        checks.append(
            build_check(
                "source_heldout_report_metrics_readable",
                readable,
                {"path": str(source), "heldout_count": len(heldouts)},
            )
        )
    else:
        checks.append(build_check("source_heldout_report_metrics_readable", False, {"path": str(source), "reason": "missing"}))
    return checks


def check_llm_regression_report(workspace: Path) -> dict[str, Any]:
    path = workspace / "validation_reports" / "llm_safe_adapter_regression_report.json"
    if not path.exists():
        return build_check("llm_report_no_unsupported_claim_regression", False, {"path": str(path), "reason": "missing"})
    report = read_json(path)
    summary = report.get("summary") or {}
    passed = summary.get("gate_status") == "passed" or (
        int(summary.get("local_failed") or summary.get("failed") or 0) == 0
        and int(summary.get("adversarial_probe_failed") or 0) == 0
    )
    return build_check(
        "llm_report_no_unsupported_claim_regression",
        passed,
        {"path": str(path), "summary": summary},
    )


def build_report(workspace: Path, run_id: str, analysis_jsons: list[Path]) -> dict[str, Any]:
    run_dir = workspace / DEFAULT_OUTPUT_ROOT / run_id
    checks = [check_cell_overlay(run_dir), *check_heldout_reports(run_dir), check_llm_regression_report(workspace)]
    for path in analysis_jsons:
        checks.extend(check_structured_prediction(structured_from_analysis(path), str(path)))
    failed = [check for check in checks if not check["passed"]]
    return {
        "runner_version": RUNNER_VERSION,
        "workspace": str(workspace),
        "run_id": run_id,
        "gate_status": "passed" if not failed else "failed",
        "summary": {"total": len(checks), "passed": len(checks) - len(failed), "failed": len(failed)},
        "checks": checks,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run final release-gate checks for a learning run.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--analysis-json", action="append", default=[], help="Analysis JSON containing structured_prediction. Repeatable.")
    parser.add_argument("--output", default="")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    workspace = Path(args.workspace).resolve()
    report = build_report(workspace, args.run_id, [Path(value).resolve() for value in args.analysis_json])
    output = Path(args.output).resolve() if args.output else workspace / DEFAULT_OUTPUT_ROOT / args.run_id / "reports" / "release_gate_validation_report.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"gate_status": report["gate_status"], "summary": report["summary"], "report_path": str(output)}, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if report["gate_status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
