"""Audit gold-standard conclusion coverage against release benchmark targets."""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DEFAULT_GOLD_PATH = "config/gold_standard_conclusions.json"
DEFAULT_REQUIREMENTS_PATH = "config/gold_standard_coverage_requirements.json"
AUDIT_VERSION = "gold_standard_coverage_audit.v1"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def normalize_text(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").casefold()).strip()


def values_text(row: dict[str, Any], keys: list[str]) -> str:
    values: list[str] = []
    for key in keys:
        value = row.get(key)
        if isinstance(value, list):
            values.extend(str(item) for item in value)
        elif value not in {None, ""}:
            values.append(str(value))
    return normalize_text(" ".join(values))


def term_hit(text: str, terms: list[str]) -> bool:
    return any(normalize_text(term) and normalize_text(term) in text for term in terms)


def positive_gold_rows(gold: dict[str, Any]) -> list[dict[str, Any]]:
    rows = gold.get("conclusions", []) or gold.get("entries", []) or []
    result = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        polarity = normalize_text(row.get("polarity") or "positive")
        is_negative = bool(row.get("is_negative_control") is True or polarity == "negative")
        if not is_negative and polarity in {"positive", ""}:
            result.append(row)
    return result


def negative_control_count(gold: dict[str, Any]) -> int:
    rows = gold.get("conclusions", []) or gold.get("entries", []) or []
    count = 0
    for row in rows:
        if not isinstance(row, dict):
            continue
        polarity = normalize_text(row.get("polarity") or "")
        if row.get("is_negative_control") is True or polarity == "negative":
            count += 1
    return count


def row_matches_cancer(row: dict[str, Any], cancer_id: str, cancer_terms: list[str]) -> bool:
    if cancer_id == "pan_cancer":
        gold_id = normalize_text(row.get("gold_id") or row.get("id") or "")
        cancer_values = row.get("cancer_type") or row.get("disease") or []
        if not isinstance(cancer_values, list):
            cancer_values = [cancer_values]
        normalized_values = {normalize_text(value) for value in cancer_values if normalize_text(value)}
        return gold_id.startswith("pan_cancer_") or bool(normalized_values and normalized_values <= {normalize_text(term) for term in cancer_terms})
    cancer_text = values_text(row, ["cancer_type", "disease", "tissue", "cell_type", "comparison", "context_terms"])
    return term_hit(cancer_text, cancer_terms)


def explicit_cell_matches(row: dict[str, Any], cancer_id: str, theme_id: str) -> bool:
    cells = row.get("benchmark_coverage_cells")
    if cells is None:
        return True
    if not isinstance(cells, list):
        return False
    expected = f"{cancer_id}:{theme_id}"
    for item in cells:
        if isinstance(item, dict):
            item_cancer = str(item.get("cancer_id") or "")
            item_theme = str(item.get("theme_id") or "")
            if item_cancer == cancer_id and item_theme == theme_id:
                return True
        elif str(item) == expected:
            return True
    return False


def matched_gold_ids(
    rows: list[dict[str, Any]],
    cancer_id: str,
    theme_id: str,
    cancer_terms: list[str],
    theme_terms: list[str],
) -> list[str]:
    matches: list[str] = []
    for row in rows:
        theme_text = values_text(row, ["mechanism_axis", "mechanism", "pathway", "theme"])
        if (
            explicit_cell_matches(row, cancer_id, theme_id)
            and row_matches_cancer(row, cancer_id, cancer_terms)
            and term_hit(theme_text, theme_terms)
        ):
            matches.append(str(row.get("gold_id") or row.get("id") or ""))
    return sorted(set(matches))


def specificity_report(cells: list[dict[str, Any]], requirements: dict[str, Any]) -> dict[str, Any]:
    recommended_max = int(requirements.get("recommended_max_cells_per_gold_id") or 2)
    gold_to_cells: dict[str, list[dict[str, str]]] = {}
    for cell in cells:
        if cell.get("status") != "covered":
            continue
        for gold_id in cell.get("matched_gold_ids") or []:
            gold_to_cells.setdefault(str(gold_id), []).append(
                {
                    "cancer_id": str(cell.get("cancer_id") or ""),
                    "theme_id": str(cell.get("theme_id") or ""),
                }
            )
    multi_cell_matches = []
    over_recommended = []
    for gold_id, matched_cells in sorted(gold_to_cells.items()):
        row = {
            "gold_id": gold_id,
            "matched_cell_count": len(matched_cells),
            "matched_cells": matched_cells,
        }
        if len(matched_cells) > 1:
            multi_cell_matches.append(row)
        if len(matched_cells) > recommended_max:
            over_recommended.append(row)
    single_cell_count = sum(1 for matched_cells in gold_to_cells.values() if len(matched_cells) == 1)
    matched_gold_count = len(gold_to_cells)
    return {
        "recommended_max_cells_per_gold_id": recommended_max,
        "matched_gold_count": matched_gold_count,
        "single_cell_gold_count": single_cell_count,
        "single_cell_gold_fraction": round(single_cell_count / matched_gold_count, 6) if matched_gold_count else None,
        "multi_cell_gold_match_count": len(multi_cell_matches),
        "max_cells_per_gold_id_observed": max((len(value) for value in gold_to_cells.values()), default=0),
        "over_recommended_count": len(over_recommended),
        "over_recommended_gold_ids": [row["gold_id"] for row in over_recommended],
        "over_recommended_matches": over_recommended,
        "multi_cell_matches": multi_cell_matches,
        "status": "review_recommended" if over_recommended else "passed",
        "interpretation": "A broad gold conclusion may cover multiple benchmark cells, but excessive reuse weakens calibration specificity; split overbroad rows into narrower disease-theme conclusions during curation.",
    }


def audit_coverage(gold: dict[str, Any], requirements: dict[str, Any]) -> dict[str, Any]:
    cancers = {
        str(row.get("cancer_id") or ""): list(row.get("terms") or [])
        for row in ((requirements.get("dimensions") or {}).get("cancers") or [])
        if row.get("cancer_id")
    }
    themes = {
        str(row.get("theme_id") or ""): list(row.get("terms") or [])
        for row in ((requirements.get("dimensions") or {}).get("themes") or [])
        if row.get("theme_id")
    }
    positives = positive_gold_rows(gold)
    cells = []
    for cell in requirements.get("required_cells", []) or []:
        cancer_id = str(cell.get("cancer_id") or "")
        theme_id = str(cell.get("theme_id") or "")
        min_controls = int(cell.get("min_positive_controls") or 1)
        gold_ids = matched_gold_ids(positives, cancer_id, theme_id, cancers.get(cancer_id, []), themes.get(theme_id, []))
        status = "covered" if len(gold_ids) >= min_controls else "missing"
        cells.append(
            {
                "cancer_id": cancer_id,
                "theme_id": theme_id,
                "min_positive_controls": min_controls,
                "matched_positive_count": len(gold_ids),
                "matched_gold_ids": gold_ids,
                "status": status,
            }
        )
    missing = [row for row in cells if row["status"] != "covered"]
    minimum_negative = int(requirements.get("minimum_negative_controls") or 0)
    negative_count = negative_control_count(gold)
    specificity = specificity_report(cells, requirements)
    return {
        "audit_version": AUDIT_VERSION,
        "created_at_utc": utc_now(),
        "gold_positive_count": len(positives),
        "gold_negative_control_count": negative_count,
        "required_cell_count": len(cells),
        "covered_cell_count": len(cells) - len(missing),
        "missing_cell_count": len(missing),
        "coverage_fraction": round((len(cells) - len(missing)) / len(cells), 6) if cells else None,
        "negative_control_requirement": {
            "minimum": minimum_negative,
            "observed": negative_count,
            "status": "covered" if negative_count >= minimum_negative else "missing",
        },
        "coverage_specificity": specificity,
        "cells": cells,
        "missing_cells": missing,
        "status": "covered" if not missing and negative_count >= minimum_negative else "gaps_present",
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit gold-standard conclusion coverage.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--gold-path", default=DEFAULT_GOLD_PATH)
    parser.add_argument("--requirements-path", default=DEFAULT_REQUIREMENTS_PATH)
    parser.add_argument("--output", default="")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    workspace = Path(args.workspace).resolve()
    report = audit_coverage(read_json(workspace / args.gold_path), read_json(workspace / args.requirements_path))
    if args.output:
        output = workspace / args.output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if report["status"] == "covered" else 1


if __name__ == "__main__":
    raise SystemExit(main())
