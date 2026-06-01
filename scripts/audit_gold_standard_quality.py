"""Audit gold-standard conclusion quality and maturity gates."""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DEFAULT_GOLD_PATH = "config/gold_standard_conclusions.json"
DEFAULT_REQUIREMENTS_PATH = "config/gold_standard_quality_requirements.json"
AUDIT_VERSION = "gold_standard_quality_audit.v1"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def normalize_text(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").casefold()).strip()


def as_list(value: Any) -> list[Any]:
    if isinstance(value, list):
        return value
    if value in {None, ""}:
        return []
    return [value]


def rows(gold: dict[str, Any]) -> list[dict[str, Any]]:
    return [row for row in (gold.get("conclusions", []) or gold.get("entries", []) or []) if isinstance(row, dict)]


def is_negative(row: dict[str, Any]) -> bool:
    return bool(row.get("is_negative_control") is True or normalize_text(row.get("polarity")) == "negative")


def positive_rows(gold: dict[str, Any]) -> list[dict[str, Any]]:
    return [row for row in rows(gold) if not is_negative(row)]


def negative_rows(gold: dict[str, Any]) -> list[dict[str, Any]]:
    return [row for row in rows(gold) if is_negative(row)]


def field_present(row: dict[str, Any], field: str) -> bool:
    value = row.get(field)
    if isinstance(value, list):
        return bool([item for item in value if str(item).strip()])
    return bool(str(value or "").strip())


def row_text(row: dict[str, Any]) -> str:
    values: list[str] = []
    for value in row.values():
        if isinstance(value, list):
            values.extend(str(item) for item in value)
        elif isinstance(value, dict):
            values.extend(str(item) for item in value.values())
        else:
            values.append(str(value))
    return normalize_text(" ".join(values))


def first_nonempty(value: Any) -> str:
    for item in as_list(value):
        text = normalize_text(item)
        if text:
            return text
    return ""


def source_ref_valid(ref: Any, allowed_prefixes: list[str]) -> bool:
    text = str(ref or "").strip()
    return bool(text and any(text.startswith(prefix) for prefix in allowed_prefixes))


def field_audit(
    selected_rows: list[dict[str, Any]],
    required_fields: list[str],
    allowed_prefixes: list[str],
    min_source_refs: int,
    min_expected_entities: int,
    require_expected_entities: bool,
) -> list[dict[str, Any]]:
    failures: list[dict[str, Any]] = []
    for row in selected_rows:
        gold_id = str(row.get("gold_id") or row.get("id") or "")
        missing = [field for field in required_fields if not field_present(row, field)]
        source_refs = as_list(row.get("source_refs"))
        invalid_refs = [str(ref) for ref in source_refs if not source_ref_valid(ref, allowed_prefixes)]
        if len(source_refs) < min_source_refs:
            missing.append("source_refs:min_count")
        if require_expected_entities and len(as_list(row.get("expected_entities"))) < min_expected_entities:
            missing.append("expected_entities:min_count")
        if missing or invalid_refs:
            failures.append(
                {
                    "gold_id": gold_id,
                    "missing_or_short_fields": sorted(set(missing)),
                    "invalid_source_refs": invalid_refs,
                }
            )
    return failures


def category_coverage(negative_controls: list[dict[str, Any]], categories: list[dict[str, Any]]) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for category in categories:
        category_id = str(category.get("category_id") or "")
        terms = [normalize_text(term) for term in as_list(category.get("terms")) if normalize_text(term)]
        matched = []
        for row in negative_controls:
            text = row_text(row)
            if any(term in text for term in terms):
                matched.append(str(row.get("gold_id") or row.get("id") or ""))
        results.append(
            {
                "category_id": category_id,
                "terms": terms,
                "matched_negative_count": len(set(matched)),
                "matched_gold_ids": sorted(set(matched)),
                "status": "covered" if matched else "missing",
            }
        )
    return results


def positive_distribution_audit(positives: list[dict[str, Any]], requirements: dict[str, Any]) -> dict[str, Any]:
    positive_count = len(positives)
    pan_cancer_rows = [
        row
        for row in positives
        if str(row.get("gold_id") or "").startswith("pan_cancer_")
        or first_nonempty(row.get("cancer_type")) in {"cancer", "tumor", "carcinoma", "pan-cancer", "pan cancer"}
    ]
    context_specific_rows = [row for row in positives if row not in pan_cancer_rows]
    cancer_contexts = {
        first_nonempty(row.get("cancer_type"))
        for row in positives
        if first_nonempty(row.get("cancer_type"))
    }
    mechanism_axes = {
        normalize_text(row.get("mechanism_axis") or row.get("mechanism") or row.get("theme"))
        for row in positives
        if normalize_text(row.get("mechanism_axis") or row.get("mechanism") or row.get("theme"))
    }
    source_refs = {
        str(ref).strip()
        for row in positives
        for ref in as_list(row.get("source_refs"))
        if str(ref).strip()
    }
    theme_requirements = [
        row
        for row in requirements.get("required_positive_theme_terms", []) or []
        if isinstance(row, dict)
    ]
    theme_coverage = []
    for requirement in theme_requirements:
        theme_id = str(requirement.get("theme_id") or "")
        terms = [normalize_text(term) for term in as_list(requirement.get("terms")) if normalize_text(term)]
        matched = [
            str(row.get("gold_id") or row.get("id") or "")
            for row in positives
            if any(term in normalize_text(row.get("mechanism_axis") or row.get("mechanism") or row.get("theme")) for term in terms)
        ]
        theme_coverage.append(
            {
                "theme_id": theme_id,
                "terms": terms,
                "matched_positive_count": len(set(matched)),
                "matched_gold_ids": sorted(set(matched)),
                "status": "covered" if matched else "missing",
            }
        )
    pan_fraction = round(len(pan_cancer_rows) / positive_count, 6) if positive_count else None
    return {
        "positive_count": positive_count,
        "pan_cancer_positive_count": len(pan_cancer_rows),
        "context_specific_positive_count": len(context_specific_rows),
        "pan_cancer_positive_fraction": pan_fraction,
        "unique_positive_cancer_context_count": len(cancer_contexts),
        "unique_positive_cancer_contexts": sorted(cancer_contexts),
        "unique_positive_mechanism_axis_count": len(mechanism_axes),
        "unique_positive_mechanism_axes": sorted(mechanism_axes),
        "unique_positive_source_ref_count": len(source_refs),
        "unique_positive_source_refs": sorted(source_refs),
        "required_theme_coverage": theme_coverage,
        "missing_required_themes": [row for row in theme_coverage if row["status"] != "covered"],
    }


def audit_quality(gold: dict[str, Any], requirements: dict[str, Any]) -> dict[str, Any]:
    positives = positive_rows(gold)
    negatives = negative_rows(gold)
    allowed_prefixes = [str(prefix) for prefix in requirements.get("allowed_source_ref_prefixes", []) or []]
    min_positive = int(requirements.get("minimum_positive_conclusions") or 0)
    min_negative = int(requirements.get("minimum_negative_controls") or 0)
    positive_field_failures = field_audit(
        positives,
        [str(field) for field in requirements.get("required_positive_fields", []) or []],
        allowed_prefixes,
        int(requirements.get("minimum_source_refs_per_positive") or 0),
        int(requirements.get("minimum_expected_entities_per_positive") or 0),
        True,
    )
    negative_field_failures = field_audit(
        negatives,
        [str(field) for field in requirements.get("required_negative_fields", []) or []],
        allowed_prefixes,
        1,
        0,
        False,
    )
    category_rows = category_coverage(negatives, requirements.get("negative_trap_categories", []) or [])
    missing_categories = [row for row in category_rows if row["status"] != "covered"]
    distribution = positive_distribution_audit(positives, requirements)
    failed_checks: list[dict[str, Any]] = []
    if len(positives) < min_positive:
        failed_checks.append({"code": "positive_gold_count_low", "observed": len(positives), "minimum": min_positive})
    if len(negatives) < min_negative:
        failed_checks.append({"code": "negative_control_count_low", "observed": len(negatives), "minimum": min_negative})
    if positive_field_failures:
        failed_checks.append({"code": "positive_gold_field_failures", "count": len(positive_field_failures)})
    if negative_field_failures:
        failed_checks.append({"code": "negative_gold_field_failures", "count": len(negative_field_failures)})
    if missing_categories:
        failed_checks.append({"code": "negative_trap_category_gaps", "count": len(missing_categories)})
    min_context_specific = int(requirements.get("minimum_context_specific_positive_conclusions") or 0)
    if distribution["context_specific_positive_count"] < min_context_specific:
        failed_checks.append(
            {
                "code": "context_specific_positive_count_low",
                "observed": distribution["context_specific_positive_count"],
                "minimum": min_context_specific,
            }
        )
    max_pan_fraction = requirements.get("maximum_pan_cancer_positive_fraction")
    if max_pan_fraction is not None and distribution["pan_cancer_positive_fraction"] is not None:
        max_pan_fraction_float = float(max_pan_fraction)
        if distribution["pan_cancer_positive_fraction"] > max_pan_fraction_float:
            failed_checks.append(
                {
                    "code": "pan_cancer_positive_fraction_high",
                    "observed": distribution["pan_cancer_positive_fraction"],
                    "maximum": max_pan_fraction_float,
                }
            )
    min_cancer_contexts = int(requirements.get("minimum_unique_positive_cancer_contexts") or 0)
    if distribution["unique_positive_cancer_context_count"] < min_cancer_contexts:
        failed_checks.append(
            {
                "code": "positive_cancer_context_diversity_low",
                "observed": distribution["unique_positive_cancer_context_count"],
                "minimum": min_cancer_contexts,
            }
        )
    min_mechanism_axes = int(requirements.get("minimum_unique_positive_mechanism_axes") or 0)
    if distribution["unique_positive_mechanism_axis_count"] < min_mechanism_axes:
        failed_checks.append(
            {
                "code": "positive_mechanism_axis_diversity_low",
                "observed": distribution["unique_positive_mechanism_axis_count"],
                "minimum": min_mechanism_axes,
            }
        )
    min_source_refs = int(requirements.get("minimum_unique_positive_source_refs") or 0)
    if distribution["unique_positive_source_ref_count"] < min_source_refs:
        failed_checks.append(
            {
                "code": "positive_source_ref_diversity_low",
                "observed": distribution["unique_positive_source_ref_count"],
                "minimum": min_source_refs,
            }
        )
    if distribution["missing_required_themes"]:
        failed_checks.append(
            {
                "code": "positive_required_theme_gaps",
                "count": len(distribution["missing_required_themes"]),
                "missing_theme_ids": [row["theme_id"] for row in distribution["missing_required_themes"]],
            }
        )
    min_theme_controls = int(requirements.get("minimum_positive_controls_per_required_theme") or 0)
    if min_theme_controls > 0:
        shallow_themes = [
            row
            for row in distribution["required_theme_coverage"]
            if int(row.get("matched_positive_count") or 0) < min_theme_controls
        ]
        if shallow_themes:
            failed_checks.append(
                {
                    "code": "positive_required_theme_depth_low",
                    "count": len(shallow_themes),
                    "minimum": min_theme_controls,
                    "theme_ids": [row["theme_id"] for row in shallow_themes],
                }
            )
    return {
        "audit_version": AUDIT_VERSION,
        "created_at_utc": utc_now(),
        "status": "passed" if not failed_checks else "failed",
        "positive_count": len(positives),
        "negative_control_count": len(negatives),
        "minimum_positive_conclusions": min_positive,
        "minimum_negative_controls": min_negative,
        "positive_field_failures": positive_field_failures,
        "negative_field_failures": negative_field_failures,
        "negative_trap_categories": category_rows,
        "missing_negative_trap_categories": missing_categories,
        "positive_distribution": distribution,
        "failed_checks": failed_checks,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit gold-standard conclusion quality.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--gold-path", default=DEFAULT_GOLD_PATH)
    parser.add_argument("--requirements-path", default=DEFAULT_REQUIREMENTS_PATH)
    parser.add_argument("--output", default="")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    workspace = Path(args.workspace).resolve()
    report = audit_quality(read_json(workspace / args.gold_path), read_json(workspace / args.requirements_path))
    if args.output:
        output = workspace / args.output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
