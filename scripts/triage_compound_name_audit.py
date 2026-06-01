#!/usr/bin/env python3
"""Triage compound-name audit rows into deterministic release policies."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


TRIAGE_VERSION = "compound_name_triage.20260520"


STRICT_ABSTAIN_ACTIONS = {
    "require_lipid_identifier",
    "require_specific_isomer_or_identifier",
    "manual_review_required",
    "avoid_broad_surface",
}
FORMAT_ACTIONS = {"quote_delimiter_sensitive_name"}
IDENTIFIER_PREFERRED_ACTIONS = {"prefer_lipid_identifier"}


def policy_for_row(row: dict[str, str]) -> str:
    action = str(row.get("recommended_action") or "")
    reasons = set(str(row.get("risk_reasons") or "").split(";"))
    uid_count = int(float(row.get("uid_count") or 0))
    risk_score = int(float(row.get("risk_score") or 0))
    if action in STRICT_ABSTAIN_ACTIONS:
        return "auto_abstain_name_only"
    if "name_maps_to_multiple_metabolites" in reasons and uid_count > 1:
        return "auto_abstain_name_only"
    if action in IDENTIFIER_PREFERRED_ACTIONS or "lipid_shorthand_or_chain_notation" in reasons:
        return "requires_identifier_or_orthogonal_feature"
    if action in FORMAT_ACTIONS:
        return "format_warning"
    if risk_score >= 70:
        return "auto_abstain_name_only"
    if risk_score >= 35:
        return "monitor_or_downweight"
    return "monitor"


def read_triage_rows(input_csv: Path, sample_per_policy: int) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    policy_counts: Counter[str] = Counter()
    action_counts: Counter[str] = Counter()
    reason_counts: Counter[str] = Counter()
    policy_risk_score_sum: defaultdict[str, int] = defaultdict(int)
    samples: dict[str, list[dict[str, Any]]] = defaultdict(list)
    total = 0
    with input_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            total += 1
            policy = policy_for_row(row)
            action = str(row.get("recommended_action") or "")
            risk_score = int(float(row.get("risk_score") or 0))
            policy_counts[policy] += 1
            action_counts[action] += 1
            policy_risk_score_sum[policy] += risk_score
            for reason in str(row.get("risk_reasons") or "").split(";"):
                if reason:
                    reason_counts[reason] += 1
            if len(samples[policy]) < sample_per_policy:
                samples[policy].append(
                    {
                        "name_key": row.get("name_key", ""),
                        "risk_score": risk_score,
                        "recommended_action": action,
                        "risk_reasons": row.get("risk_reasons", ""),
                        "uid_count": int(float(row.get("uid_count") or 0)),
                        "canonical_examples": row.get("canonical_examples", ""),
                        "recommendation": row.get("recommendation", ""),
                    }
                )
    summary = {
        "triage_version": TRIAGE_VERSION,
        "input_csv": str(input_csv),
        "total_risk_rows": total,
        "policy_counts": dict(sorted(policy_counts.items())),
        "recommended_action_counts": dict(sorted(action_counts.items())),
        "risk_reason_counts": dict(sorted(reason_counts.items())),
        "policy_mean_risk_score": {
            policy: round(policy_risk_score_sum[policy] / count, 3)
            for policy, count in sorted(policy_counts.items())
            if count
        },
        "release_policy": {
            "auto_abstain_name_only": "Do not auto-resolve from this surface alone; require stable ID, exact structure, RT/MS2, or other orthogonal evidence.",
            "requires_identifier_or_orthogonal_feature": "Use only with stable database IDs or orthogonal assay features; otherwise keep as ambiguous/pool-level.",
            "format_warning": "Do not split delimiter-sensitive names; prefer JSON or quoted CSV input.",
            "monitor_or_downweight": "Allow only under existing resolver thresholds and confidence margins; keep warning surface in audit reports.",
            "monitor": "Track as low-priority audit surface.",
        },
    }
    return summary, samples


def markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(item).replace("\n", " ") for item in row) + " |")
    return "\n".join(lines)


def render_markdown(summary: dict[str, Any], samples: dict[str, list[dict[str, Any]]]) -> str:
    policy_counts = summary.get("policy_counts") or {}
    action_counts = summary.get("recommended_action_counts") or {}
    reason_counts = summary.get("risk_reason_counts") or {}
    lines = [
        "# Compound Name Audit Triage",
        "",
        f"Input CSV: `{summary.get('input_csv', '')}`",
        f"Total risk rows: `{summary.get('total_risk_rows', 0)}`",
        "",
        "## Policy Counts",
        "",
        markdown_table(
            ["policy", "count", "mean_risk_score", "release handling"],
            [
                [
                    policy,
                    count,
                    (summary.get("policy_mean_risk_score") or {}).get(policy, ""),
                    (summary.get("release_policy") or {}).get(policy, ""),
                ]
                for policy, count in sorted(policy_counts.items())
            ],
        ),
        "",
        "## Recommended Actions",
        "",
        markdown_table(["action", "count"], [[key, value] for key, value in sorted(action_counts.items())]),
        "",
        "## Top Risk Reasons",
        "",
        markdown_table(
            ["reason", "count"],
            [[key, value] for key, value in sorted(reason_counts.items(), key=lambda item: (-item[1], item[0]))[:25]],
        ),
        "",
        "## Sample Rows By Policy",
        "",
    ]
    for policy, rows in sorted(samples.items()):
        lines.extend(
            [
                f"### {policy}",
                "",
                markdown_table(
                    ["name_key", "risk_score", "action", "uid_count", "reasons", "recommendation"],
                    [
                        [
                            row.get("name_key", ""),
                            row.get("risk_score", ""),
                            row.get("recommended_action", ""),
                            row.get("uid_count", ""),
                            row.get("risk_reasons", ""),
                            row.get("recommendation", ""),
                        ]
                        for row in rows
                    ],
                ),
                "",
            ]
        )
    lines.extend(
        [
            "## Interpretation",
            "",
            "- These rows are not manually adjudicated one by one.",
            "- The release policy is conservative: high-risk name-only surfaces abstain or require stronger evidence.",
            "- Manual adjudication remains a future optimization step for sampled rows or high-impact use cases.",
            "",
        ]
    )
    return "\n".join(lines)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Triage compound-name audit candidates into release policies.")
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--output-json", default="")
    parser.add_argument("--output-md", default="")
    parser.add_argument("--sample-per-policy", type=int, default=20)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    input_csv = Path(args.input_csv).resolve()
    output_json = Path(args.output_json).resolve() if args.output_json else input_csv.parent / "compound_name_triage_summary.json"
    output_md = Path(args.output_md).resolve() if args.output_md else input_csv.parent / "compound_name_triage_summary.md"
    summary, samples = read_triage_rows(input_csv, max(0, args.sample_per_policy))
    payload = {**summary, "samples_by_policy": samples}
    output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    output_md.write_text(render_markdown(summary, samples), encoding="utf-8")
    print(json.dumps({"output_json": str(output_json), "output_md": str(output_md), "total_risk_rows": summary["total_risk_rows"]}, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
