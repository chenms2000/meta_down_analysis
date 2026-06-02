"""Audit registered real-world benchmark fixture files."""

from __future__ import annotations

import argparse
import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AUDIT_VERSION = "real_world_fixture_audit.v1"
DEFAULT_REGISTRY = "config/real_world_fixture_registry.json"
SKIP_SEARCH_DIRS = {".git", ".pytest_cache", "__pycache__", "node_modules"}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def parse_optional_float(value: Any) -> float | None:
    try:
        text = str(value or "").strip()
        return float(text) if text else None
    except (TypeError, ValueError):
        return None


def column_lookup(columns: list[str]) -> dict[str, str]:
    return {str(column).strip().casefold(): str(column) for column in columns if str(column).strip()}


def first_present_column(columns: list[str], candidates: list[str]) -> str:
    lookup = column_lookup(columns)
    for candidate in candidates:
        column = lookup.get(str(candidate).casefold())
        if column:
            return column
    return ""


def infer_input_mode(columns: list[str]) -> str:
    lookup = column_lookup(columns)
    if "trait" in lookup or "reported_trait" in lookup or "trait_name" in lookup:
        return "trait_score"
    if "hmdb" in lookup or "chebi" in lookup:
        return "stable_identifier"
    if "mz" in lookup:
        return "lcms_feature"
    if any(column in lookup for column in ("name", "metabolite", "compound")):
        return "compound_name"
    return "unknown"


def expected_schema(fixture: dict[str, Any]) -> dict[str, Any]:
    return {
        "identity_columns_any": [str(column) for column in fixture.get("identity_columns_any", []) or []],
        "effect_columns_any": [str(column) for column in fixture.get("effect_columns_any", []) or []],
        "expected_input_mode": str(fixture.get("expected_input_mode") or ""),
        "min_rows": int(fixture.get("min_rows") or 1),
    }


def recovery_rerun_plan(fixture: dict[str, Any]) -> dict[str, Any]:
    case_id = str(fixture.get("case_id") or "")
    path = str(fixture.get("path") or "")
    return {
        "plan_version": "fixture_recovery_rerun_plan.v1",
        "case_id": case_id,
        "required_for_strict_release": bool(fixture.get("required_for_strict_release")),
        "preconditions": [
            {
                "check_id": "restore_registered_path_or_update_registry",
                "description": "Restore the registered fixture file path or update the registry to a schema-valid candidate path.",
                "expected_path": path,
            },
            {
                "check_id": "schema_validates_before_analysis",
                "description": "The fixture must pass identity/effect column, numeric effect, min-row, and expected input-mode checks before analysis.",
                "expected_schema": expected_schema(fixture),
            },
        ],
        "validation_sequence": [
            {
                "step_id": "fixture_schema_audit",
                "purpose": "Validate restored fixture schema before running the full analysis path.",
                "command": [
                    "python",
                    "scripts/audit_real_world_fixtures.py",
                    "--workspace",
                    ".",
                    "--registry",
                    "config/real_world_fixture_registry.json",
                ],
            },
            {
                "step_id": "strict_batch_benchmark",
                "purpose": "Rerun strict batch benchmark across multi-cancer, multi-theme, multi-input release scenarios with held-out reruns.",
                "command": [
                    "python",
                    "scripts/run_release_batch_benchmark.py",
                    "--workspace",
                    ".",
                    "--release-id",
                    "<release_id>",
                    "--strict-release",
                    "--heldout-top-n",
                    "200",
                    "--max-paths",
                    "10",
                    "--max-hops",
                    "3",
                ],
            },
            {
                "step_id": "post_batch_release_audits",
                "purpose": "Refresh readiness, metric snapshot, parameter optimization, objective status, and strict workflow gates.",
                "command": [
                    "python",
                    "scripts/run_strict_release_workflow.py",
                    "--workspace",
                    ".",
                    "--release-id",
                    "<release_id>",
                    "--heldout-top-n",
                    "200",
                    "--max-paths",
                    "10",
                    "--max-hops",
                    "3",
                ],
            },
        ],
        "success_criteria": [
            "fixture audit status passed",
            "strict batch missing_case_count is 0",
            "benchmark matrix executed gaps are empty",
            "scenario traceability contract passed",
            "negative trap hits remain 0",
            "objective status no longer blocked by this fixture",
        ],
    }


def candidate_tokens(fixture: dict[str, Any]) -> list[str]:
    raw = " ".join(
        [
            str(fixture.get("case_id") or ""),
            str(fixture.get("path") or ""),
            str(fixture.get("fixture_type") or ""),
        ]
    )
    tokens = {
        token.casefold()
        for token in re.split(r"[^A-Za-z0-9]+", raw)
        if len(token) >= 3 and token.casefold() not in {"csv", "trait", "score", "fixture"}
    }
    return sorted(tokens)


def iter_workspace_files(workspace: Path):
    stack = [workspace]
    while stack:
        current = stack.pop()
        try:
            children = list(current.iterdir())
        except (OSError, PermissionError):
            continue
        for child in children:
            if child.is_dir():
                if child.name not in SKIP_SEARCH_DIRS:
                    stack.append(child)
            elif child.is_file():
                yield child


def find_candidate_paths(workspace: Path, fixture: dict[str, Any], limit: int = 8) -> list[dict[str, Any]]:
    tokens = candidate_tokens(fixture)
    if not tokens:
        return []
    candidates: list[dict[str, Any]] = []
    for path in iter_workspace_files(workspace):
        if path.suffix.casefold() not in {".csv", ".tsv", ".txt", ".xlsx", ".xls"}:
            continue
        relative = str(path.relative_to(workspace))
        haystack = relative.casefold()
        matched = [token for token in tokens if token in haystack]
        if not matched:
            continue
        candidates.append(
            {
                "path": str(path),
                "relative_path": relative,
                "matched_tokens": matched,
                "score": len(matched),
            }
        )
    candidates.sort(key=lambda row: (-int(row.get("score") or 0), str(row.get("relative_path") or "")))
    return candidates[:limit]


def audit_csv_fixture(path: Path, fixture: dict[str, Any]) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        columns = list(reader.fieldnames or [])
    identity_column = first_present_column(columns, [str(column) for column in fixture.get("identity_columns_any", []) or []])
    effect_column = first_present_column(columns, [str(column) for column in fixture.get("effect_columns_any", []) or []])
    numeric_effect_rows = 0
    if effect_column:
        numeric_effect_rows = sum(1 for row in rows if parse_optional_float(row.get(effect_column)) is not None)
    inferred_mode = infer_input_mode(columns)
    errors: list[dict[str, Any]] = []
    min_rows = int(fixture.get("min_rows") or 1)
    if len(rows) < min_rows:
        errors.append({"code": "row_count_below_minimum", "message": f"Observed {len(rows)} rows, required {min_rows}."})
    if not identity_column:
        errors.append({"code": "missing_identity_column", "message": "No registered identity column was found."})
    if not effect_column:
        errors.append({"code": "missing_effect_column", "message": "No registered numeric effect/ranking column was found."})
    elif numeric_effect_rows <= 0:
        errors.append({"code": "effect_column_not_numeric", "message": f"Effect column {effect_column} has no numeric values."})
    expected_mode = str(fixture.get("expected_input_mode") or "")
    if expected_mode and inferred_mode != expected_mode:
        errors.append({"code": "input_mode_mismatch", "message": f"Expected {expected_mode}, observed {inferred_mode}."})
    return {
        "path": str(path),
        "exists": True,
        "row_count": len(rows),
        "columns": columns,
        "identity_column": identity_column,
        "effect_column": effect_column,
        "numeric_effect_rows": numeric_effect_rows,
        "expected_input_mode": expected_mode,
        "observed_input_mode": inferred_mode,
        "status": "passed" if not errors else "failed",
        "errors": errors,
    }


def audit_fixture(workspace: Path, fixture: dict[str, Any]) -> dict[str, Any]:
    case_id = str(fixture.get("case_id") or "")
    path = workspace / str(fixture.get("path") or "")
    if not path.exists():
        return {
            "case_id": case_id,
            "path": str(path),
            "fixture_type": str(fixture.get("fixture_type") or ""),
            "required_for_strict_release": bool(fixture.get("required_for_strict_release")),
            "exists": False,
            "status": "missing",
            "expected_schema": expected_schema(fixture),
            "candidate_paths": find_candidate_paths(workspace, fixture),
            "recovery_action": "Restore the registered fixture path or update the registry to one of the candidate paths after schema validation.",
            "recovery_rerun_plan": recovery_rerun_plan(fixture),
            "errors": [{"code": "fixture_missing", "message": "Registered real-world fixture file is missing."}],
        }
    if str(fixture.get("fixture_type") or "").startswith("csv"):
        result = audit_csv_fixture(path, fixture)
    else:
        result = {"path": str(path), "exists": True, "status": "failed", "errors": [{"code": "unsupported_fixture_type"}]}
    result.update(
        {
            "case_id": case_id,
            "fixture_type": str(fixture.get("fixture_type") or ""),
            "required_for_strict_release": bool(fixture.get("required_for_strict_release")),
        }
    )
    return result


def audit_registry(workspace: Path, registry: dict[str, Any]) -> dict[str, Any]:
    fixtures = [row for row in registry.get("fixtures", []) or [] if isinstance(row, dict)]
    results = [audit_fixture(workspace, row) for row in fixtures]
    required = [row for row in results if row.get("required_for_strict_release")]
    failed_required = [row for row in required if row.get("status") != "passed"]
    failed_any = [row for row in results if row.get("status") not in {"passed"}]
    return {
        "audit_version": AUDIT_VERSION,
        "created_at_utc": utc_now(),
        "schema_version": registry.get("schema_version", ""),
        "fixture_count": len(results),
        "required_fixture_count": len(required),
        "passed_count": sum(1 for row in results if row.get("status") == "passed"),
        "missing_count": sum(1 for row in results if row.get("status") == "missing"),
        "failed_count": sum(1 for row in results if row.get("status") == "failed"),
        "required_failed_count": len(failed_required),
        "status": "passed" if not failed_required else "failed",
        "fixtures": results,
        "failed_required": failed_required,
        "failed_any": failed_any,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit registered real-world benchmark fixtures.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--registry", default=DEFAULT_REGISTRY)
    parser.add_argument("--output", default="")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    workspace = Path(args.workspace).resolve()
    report = audit_registry(workspace, read_json(workspace / args.registry))
    if args.output:
        write_json(workspace / args.output, report)
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
