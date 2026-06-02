"""Prepare a formal three-tier workflow for differential TraitScore tables.

The workflow keeps the full allDEGs table as the audit/background source while
promoting only selected high-signal rows into graph explanation runs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


WORKFLOW_VERSION = "formal_differential_workflow.v1"
DEFAULT_RELEASE_ID = "mvp_20260513T002254"
DEFAULT_OUTPUT_ROOT = "validation_reports"
DEFAULT_CORE_MAX_RECORDS = 120
DEFAULT_CORE_TOP_PER_GROUP = 60
DEFAULT_EXPLORATORY_MAX_RECORDS = 500
DEFAULT_EXPLORATORY_TOP_PER_GROUP = 250
Q_COLUMNS = ("q_wilcoxon", "q_ttest", "padj", "qvalue", "fdr")
P_COLUMNS = ("p_wilcoxon", "p_ttest", "pvalue", "p")
EFFECT_COLUMNS = (
    "mean_diff",
    "z_wilcoxon",
    "cohen_d",
    "pseudo_log2FC_shifted",
    "signed_log2FC_absmean",
    "log2FC",
    "median_diff",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_table(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [dict(row) for row in reader]
        return rows, list(reader.fieldnames or [])


def write_table(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def parse_float(value: Any) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        value_float = float(text)
    except ValueError:
        return None
    if not math.isfinite(value_float):
        return None
    return value_float


def first_present(fieldnames: list[str], candidates: tuple[str, ...]) -> str:
    lowered = {name.lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return ""


def row_effect(row: dict[str, str], effect_column: str) -> float:
    return parse_float(row.get(effect_column)) or 0.0


def row_q_value(row: dict[str, str], q_column: str) -> float | None:
    return parse_float(row.get(q_column)) if q_column else None


def row_direction(row: dict[str, str]) -> str:
    value = str(row.get("direction") or row.get("comparison_direction") or "unknown").strip()
    return value or "unknown"


def row_rank_key(row: dict[str, str], q_column: str, effect_column: str) -> tuple[float, float, str]:
    q_value = row_q_value(row, q_column)
    q_rank = q_value if q_value is not None else float("inf")
    effect_rank = abs(row_effect(row, effect_column))
    trait = str(row.get("trait") or row.get("GCST") or row.get("accession_id") or "")
    return (q_rank, -effect_rank, trait)


def select_rows(
    rows: list[dict[str, str]],
    q_column: str,
    effect_column: str,
    max_records: int,
    top_per_group: int,
    q_threshold: float | None,
    min_abs_effect: float,
) -> list[dict[str, str]]:
    filtered = []
    for row in rows:
        q_value = row_q_value(row, q_column)
        effect = abs(row_effect(row, effect_column))
        if q_threshold is not None and (q_value is None or q_value > q_threshold):
            continue
        if effect < min_abs_effect:
            continue
        filtered.append(row)

    grouped: dict[str, list[dict[str, str]]] = {}
    for row in filtered:
        grouped.setdefault(row_direction(row), []).append(row)

    selected: list[dict[str, str]] = []
    for direction in sorted(grouped):
        ranked = sorted(grouped[direction], key=lambda row: row_rank_key(row, q_column, effect_column))
        selected.extend(ranked[:top_per_group])
    return sorted(selected, key=lambda row: row_rank_key(row, q_column, effect_column))[:max_records]


def audit_rows(rows: list[dict[str, str]], fieldnames: list[str], q_column: str, effect_column: str) -> dict[str, Any]:
    q_values = [value for row in rows if (value := row_q_value(row, q_column)) is not None]
    effect_values = [abs(row_effect(row, effect_column)) for row in rows]
    direction_counts = Counter(row_direction(row) for row in rows)
    trait_examples = [
        str(row.get("trait") or row.get("GCST") or row.get("accession_id") or "")
        for row in sorted(rows, key=lambda row: row_rank_key(row, q_column, effect_column))[:10]
    ]
    required_columns = {"trait", "group1", "group2", "direction"}
    missing_required = sorted(column for column in required_columns if column not in fieldnames)
    return {
        "row_count": len(rows),
        "columns": fieldnames,
        "q_column": q_column,
        "effect_column": effect_column,
        "missing_required_columns": missing_required,
        "direction_counts": dict(sorted(direction_counts.items())),
        "q_threshold_counts": {
            "q_le_0.05": sum(1 for value in q_values if value <= 0.05),
            "q_le_0.10": sum(1 for value in q_values if value <= 0.10),
            "q_available": len(q_values),
        },
        "effect_summary": {
            "available": len(effect_values),
            "max_abs": max(effect_values) if effect_values else 0.0,
            "mean_abs": sum(effect_values) / len(effect_values) if effect_values else 0.0,
        },
        "top_trait_examples": trait_examples,
    }


def service_command(
    workspace: Path,
    release_id: str,
    input_path: Path,
    max_records: int,
    top_per_group: int,
    max_paths: int,
    max_hops: int,
) -> list[str]:
    return [
        sys.executable,
        str(workspace / "scripts" / "metabo_service.py"),
        "--workspace",
        str(workspace),
        "--release-id",
        release_id,
        "analyze-differential-table",
        str(input_path),
        "--max-records",
        str(max_records),
        "--top-per-group",
        str(top_per_group),
        "--max-paths",
        str(max_paths),
        "--max-hops",
        str(max_hops),
    ]


def powershell_command(command: list[str]) -> str:
    return " ".join(f'"{part}"' if " " in part else part for part in command)


def run_command(command: list[str], workspace: Path, output_path: Path) -> dict[str, Any]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(command, cwd=workspace, capture_output=True, text=True, encoding="utf-8")
    output_path.write_text(completed.stdout, encoding="utf-8")
    return {
        "command": command,
        "output_path": str(output_path),
        "returncode": completed.returncode,
        "passed": completed.returncode == 0,
        "stderr_tail": completed.stderr[-4000:],
    }


def workflow_markdown(report: dict[str, Any]) -> str:
    tiers = report.get("tiers") or {}
    lines = [
        "# Formal Differential Workflow",
        "",
        f"Status: `{report.get('status', '')}`",
        f"Workflow version: `{report.get('workflow_version', '')}`",
        f"Release ID: `{report.get('release_id', '')}`",
        "",
        "## Full-Table Audit",
        "",
        f"- Rows: `{(report.get('audit') or {}).get('row_count', 0)}`",
        f"- q column: `{(report.get('audit') or {}).get('q_column', '')}`",
        f"- effect column: `{(report.get('audit') or {}).get('effect_column', '')}`",
        "",
        "## Tiers",
        "",
    ]
    for tier_id in ("core_report", "exploratory_appendix", "full_audit"):
        tier = tiers.get(tier_id) or {}
        lines.extend(
            [
                f"### {tier.get('title', tier_id)}",
                "",
                f"- Role: {tier.get('role', '')}",
                f"- Input rows: `{tier.get('input_rows', 0)}`",
                f"- Output: `{tier.get('input_path', '')}`",
                f"- Interpretation boundary: {tier.get('boundary', '')}",
                "",
            ]
        )
        if tier.get("command"):
            lines.append(f"```powershell\n{powershell_command(tier['command'])}\n```")
            lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def find_default_significant_path(all_degs_path: Path) -> Path | None:
    candidate_name = all_degs_path.name.replace(".csv", "_significant_q0.05.csv")
    candidate = all_degs_path.parents[1] / "significant" / candidate_name if len(all_degs_path.parents) > 1 else all_degs_path
    return candidate if candidate.exists() else None


def run_workflow(args: argparse.Namespace) -> dict[str, Any]:
    workspace = Path(args.workspace).resolve()
    all_degs_path = Path(args.all_degs).resolve()
    rows, fieldnames = read_table(all_degs_path)
    q_column = first_present(fieldnames, Q_COLUMNS)
    effect_column = first_present(fieldnames, EFFECT_COLUMNS)
    if not q_column:
        raise ValueError(f"No q-value column found. Expected one of: {', '.join(Q_COLUMNS)}")
    if not effect_column:
        raise ValueError(f"No effect column found. Expected one of: {', '.join(EFFECT_COLUMNS)}")

    significant_path = Path(args.significant).resolve() if args.significant else find_default_significant_path(all_degs_path)
    core_source_path = significant_path if significant_path and significant_path.exists() else all_degs_path
    core_rows_source, core_fieldnames = read_table(core_source_path)
    core_q_column = first_present(core_fieldnames, Q_COLUMNS)
    core_effect_column = first_present(core_fieldnames, EFFECT_COLUMNS)

    output_dir = (
        workspace
        / args.output_root
        / args.release_id
        / "formal_differential_workflow"
        / all_degs_path.stem
    )
    core_input = output_dir / "core_input.csv"
    exploratory_input = output_dir / "exploratory_input.csv"
    commands_path = output_dir / "run_formal_analysis.ps1"
    report_json = output_dir / "formal_differential_workflow_report.json"
    report_md = output_dir / "formal_differential_workflow_report.md"

    core_rows = select_rows(
        core_rows_source,
        core_q_column,
        core_effect_column,
        args.core_max_records,
        args.core_top_per_group,
        args.core_q_threshold,
        args.min_abs_effect,
    )
    exploratory_rows = select_rows(
        rows,
        q_column,
        effect_column,
        args.exploratory_max_records,
        args.exploratory_top_per_group,
        args.exploratory_q_threshold,
        args.min_abs_effect,
    )
    write_table(core_input, core_rows, core_fieldnames)
    write_table(exploratory_input, exploratory_rows, fieldnames)

    core_command = service_command(
        workspace,
        args.release_id,
        core_input,
        args.core_max_records,
        args.core_top_per_group,
        args.core_max_paths,
        args.core_max_hops,
    )
    exploratory_command = service_command(
        workspace,
        args.release_id,
        exploratory_input,
        args.exploratory_max_records,
        args.exploratory_top_per_group,
        args.exploratory_max_paths,
        args.exploratory_max_hops,
    )
    commands_path.write_text(
        "\n".join(
            [
                "# Generated formal differential workflow commands.",
                "# Run core first; exploratory is appendix-only and may be slower.",
                powershell_command(core_command) + " | Set-Content -Encoding UTF8 core_analysis.json",
                powershell_command(exploratory_command) + " | Set-Content -Encoding UTF8 exploratory_appendix.json",
                "",
            ]
        ),
        encoding="utf-8",
    )

    report: dict[str, Any] = {
        "status": "prepared",
        "workflow_version": WORKFLOW_VERSION,
        "created_at": utc_now(),
        "release_id": args.release_id,
        "inputs": {
            "all_degs": str(all_degs_path),
            "significant": str(significant_path) if significant_path and significant_path.exists() else "",
            "core_source": str(core_source_path),
        },
        "audit": audit_rows(rows, fieldnames, q_column, effect_column),
        "tiers": {
            "core_report": {
                "title": "Core report",
                "role": "Primary research interpretation from high-signal rows.",
                "input_rows": len(core_rows),
                "input_path": str(core_input),
                "source_path": str(core_source_path),
                "command": core_command,
                "boundary": "Use for core conclusions only when entity resolution, evidence refs, and confidence gates pass.",
            },
            "exploratory_appendix": {
                "title": "Exploratory appendix",
                "role": "Broader candidate discovery and lower-confidence themes.",
                "input_rows": len(exploratory_rows),
                "input_path": str(exploratory_input),
                "source_path": str(all_degs_path),
                "command": exploratory_command,
                "boundary": "Appendix-only; graph/model-only or ambiguous results remain hypotheses.",
            },
            "full_audit": {
                "title": "Full audit",
                "role": "Whole-table row counts, direction balance, significance distribution, and review burden.",
                "input_rows": len(rows),
                "input_path": str(all_degs_path),
                "boundary": "Do not promote all full-table rows into strong conclusions.",
            },
        },
        "outputs": {
            "core_input_csv": str(core_input),
            "exploratory_input_csv": str(exploratory_input),
            "commands_ps1": str(commands_path),
            "workflow_json": str(report_json),
            "workflow_markdown": str(report_md),
        },
        "execution": {},
    }

    if args.execute_core:
        report["execution"]["core_report"] = run_command(core_command, workspace, output_dir / "core_analysis.json")
    if args.execute_exploratory:
        report["execution"]["exploratory_appendix"] = run_command(
            exploratory_command,
            workspace,
            output_dir / "exploratory_appendix.json",
        )

    write_json(report_json, report)
    report_md.write_text(workflow_markdown(report), encoding="utf-8")
    return report


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare a formal differential-table workflow.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", default=DEFAULT_RELEASE_ID)
    parser.add_argument("--all-degs", required=True, help="Full allDEGs differential table.")
    parser.add_argument("--significant", default="", help="Optional significant q0.05 subset. Auto-detected when omitted.")
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--core-max-records", type=int, default=DEFAULT_CORE_MAX_RECORDS)
    parser.add_argument("--core-top-per-group", type=int, default=DEFAULT_CORE_TOP_PER_GROUP)
    parser.add_argument("--core-q-threshold", type=float, default=0.05)
    parser.add_argument("--core-max-paths", type=int, default=25)
    parser.add_argument("--core-max-hops", type=int, default=4)
    parser.add_argument("--exploratory-max-records", type=int, default=DEFAULT_EXPLORATORY_MAX_RECORDS)
    parser.add_argument("--exploratory-top-per-group", type=int, default=DEFAULT_EXPLORATORY_TOP_PER_GROUP)
    parser.add_argument("--exploratory-q-threshold", type=float, default=None)
    parser.add_argument("--exploratory-max-paths", type=int, default=15)
    parser.add_argument("--exploratory-max-hops", type=int, default=3)
    parser.add_argument("--min-abs-effect", type=float, default=0.0)
    parser.add_argument("--execute-core", action="store_true")
    parser.add_argument("--execute-exploratory", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    report = run_workflow(parse_args(argv))
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
