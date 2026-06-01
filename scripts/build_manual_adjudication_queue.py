#!/usr/bin/env python3
"""Build a sampled manual adjudication queue from compound-name triage rows."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Callable

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from triage_compound_name_audit import TRIAGE_VERSION, policy_for_row  # noqa: E402


QUEUE_VERSION = "manual_adjudication_queue.20260520"


def intish(value: Any) -> int:
    try:
        return int(float(value or 0))
    except (TypeError, ValueError):
        return 0


def reasons(row: dict[str, str]) -> set[str]:
    return {reason for reason in str(row.get("risk_reasons") or "").split(";") if reason}


def has_reason(reason: str) -> Callable[[dict[str, str], str], bool]:
    return lambda row, policy: reason in reasons(row)


BUCKETS: list[dict[str, Any]] = [
    {
        "bucket": "highest_risk_name_only",
        "predicate": lambda row, policy: policy == "auto_abstain_name_only",
        "question": "Should this name-only surface remain abstained, or can a stable ID/structure be assigned?",
    },
    {
        "bucket": "lipid_identifier_required",
        "predicate": lambda row, policy: policy == "requires_identifier_or_orthogonal_feature"
        or "lipid_shorthand_or_chain_notation" in reasons(row),
        "question": "Which stable lipid identifier, exact structure, RT, or MS2 evidence would be required before exact identity?",
    },
    {
        "bucket": "multi_entity_surface",
        "predicate": has_reason("name_maps_to_multiple_metabolites"),
        "question": "Does this surface represent one entity, a metabolite pool, or an unresolved multi-map name?",
    },
    {
        "bucket": "stereo_or_isomer_surface",
        "predicate": has_reason("stereo_or_isomer_under_specified"),
        "question": "Is stereochemistry/isomer identity specified enough for exact resolution?",
    },
    {
        "bucket": "delimiter_sensitive_format",
        "predicate": lambda row, policy: policy == "format_warning",
        "question": "Is this a formatting/serialization problem rather than a biological identity problem?",
    },
    {
        "bucket": "low_priority_synonym_monitor",
        "predicate": lambda row, policy: policy in {"monitor", "monitor_or_downweight"}
        or "synonym_only_low_rank_surface" in reasons(row),
        "question": "Should this synonym surface stay monitored/downweighted, or is it safe under current thresholds?",
    },
]


OUTPUT_FIELDS = [
    "review_bucket",
    "review_priority",
    "triage_policy",
    "review_question",
    "adjudication_status",
    "reviewer_decision",
    "reviewer_entity_id",
    "reviewer_notes",
    "risk_score",
    "name_key",
    "recommended_action",
    "recommendation",
    "risk_reasons",
    "uid_count",
    "row_count",
    "canonical_examples",
    "authority_prefixes",
    "match_fields",
    "raw_examples",
]


def row_sort_key(row: dict[str, str]) -> tuple[int, int, int, str]:
    return (-intish(row.get("risk_score")), -intish(row.get("row_count")), -intish(row.get("uid_count")), str(row.get("name_key") or ""))


def read_rows(input_csv: Path) -> list[dict[str, str]]:
    with input_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def build_queue(rows: list[dict[str, str]], sample_per_bucket: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    prepared = []
    policy_counts: Counter[str] = Counter()
    for row in rows:
        policy = policy_for_row(row)
        policy_counts[policy] += 1
        prepared.append((row, policy))

    selected: list[dict[str, Any]] = []
    selected_keys: set[str] = set()
    bucket_counts: Counter[str] = Counter()
    for bucket in BUCKETS:
        matches = [
            (row, policy)
            for row, policy in prepared
            if row.get("name_key", "") not in selected_keys and bucket["predicate"](row, policy)
        ]
        for rank, (row, policy) in enumerate(sorted(matches, key=lambda item: row_sort_key(item[0]))[:sample_per_bucket], start=1):
            selected_keys.add(row.get("name_key", ""))
            bucket_counts[bucket["bucket"]] += 1
            selected.append(
                {
                    "review_bucket": bucket["bucket"],
                    "review_priority": rank,
                    "triage_policy": policy,
                    "review_question": bucket["question"],
                    "adjudication_status": "pending",
                    "reviewer_decision": "",
                    "reviewer_entity_id": "",
                    "reviewer_notes": "",
                    **{field: row.get(field, "") for field in OUTPUT_FIELDS if field not in {
                        "review_bucket",
                        "review_priority",
                        "triage_policy",
                        "review_question",
                        "adjudication_status",
                        "reviewer_decision",
                        "reviewer_entity_id",
                        "reviewer_notes",
                    }},
                }
            )

    summary = {
        "queue_version": QUEUE_VERSION,
        "triage_version": TRIAGE_VERSION,
        "total_input_rows": len(rows),
        "selected_rows": len(selected),
        "sample_per_bucket": sample_per_bucket,
        "bucket_counts": dict(sorted(bucket_counts.items())),
        "policy_counts": dict(sorted(policy_counts.items())),
        "review_fields": ["adjudication_status", "reviewer_decision", "reviewer_entity_id", "reviewer_notes"],
        "interpretation": "This queue is a sampled manual follow-up target, not a completed adjudication result.",
    }
    return selected, summary


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(item).replace("\n", " ") for item in row) + " |")
    return "\n".join(lines)


def render_markdown(summary: dict[str, Any], queue_csv: Path) -> str:
    return "\n".join(
        [
            "# Manual Adjudication Queue",
            "",
            f"Queue CSV: `{queue_csv}`",
            f"Total input risk rows: `{summary['total_input_rows']}`",
            f"Selected review rows: `{summary['selected_rows']}`",
            f"Sample per bucket: `{summary['sample_per_bucket']}`",
            "",
            "## Review Buckets",
            "",
            markdown_table(["bucket", "selected"], [[key, value] for key, value in summary["bucket_counts"].items()]),
            "",
            "## Triage Policy Counts",
            "",
            markdown_table(["policy", "input_rows"], [[key, value] for key, value in summary["policy_counts"].items()]),
            "",
            "## Interpretation",
            "",
            "- This is a sampled queue for future manual adjudication, not a completed gold-standard label set.",
            "- Rows are selected from the highest-risk deterministic triage surfaces first.",
            "- Reviewer fields are intentionally blank except for `adjudication_status=pending`.",
            "",
        ]
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a sampled manual adjudication queue from compound-name audit candidates.")
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--output-csv", default="")
    parser.add_argument("--output-json", default="")
    parser.add_argument("--output-md", default="")
    parser.add_argument("--sample-per-bucket", type=int, default=25)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    input_csv = Path(args.input_csv).resolve()
    output_csv = Path(args.output_csv).resolve() if args.output_csv else input_csv.parent / "manual_adjudication_queue.csv"
    output_json = Path(args.output_json).resolve() if args.output_json else input_csv.parent / "manual_adjudication_queue_summary.json"
    output_md = Path(args.output_md).resolve() if args.output_md else input_csv.parent / "manual_adjudication_queue_summary.md"
    rows = read_rows(input_csv)
    selected, summary = build_queue(rows, max(0, args.sample_per_bucket))
    write_csv(output_csv, selected)
    payload = {**summary, "input_csv": str(input_csv), "output_csv": str(output_csv), "output_md": str(output_md)}
    output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    output_md.write_text(render_markdown(payload, output_csv), encoding="utf-8")
    print(json.dumps({"output_csv": str(output_csv), "output_json": str(output_json), "output_md": str(output_md), "selected_rows": len(selected)}, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
