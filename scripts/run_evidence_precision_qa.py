#!/usr/bin/env python3
"""Generate deterministic Evidence Precision QA notes for a frozen release."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

try:
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover
    pq = None

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from build_literature_evidence import (  # noqa: E402
    DEFAULT_OUTPUT_ROOT,
    DEFAULT_PRECISION_FILTER_CONFIG,
    GENERIC_SURFACES,
    content_hash,
    load_precision_filters,
    normalized_surface_key,
    read_json,
    rule_matches,
)


DEFAULT_OUTPUT_PATH = "docs/evidence_precision_notes.md"
DEFAULT_MIN_SURFACE_COUNT = 250
DEFAULT_SUPPORTED_EDGE_RETENTION = 0.80


@dataclass
class SurfaceAgg:
    entity_type: str
    normalized_surface: str
    matched_field: str
    mention_count: int = 0
    entity_uids: set[str] = field(default_factory=set)
    display_names: Counter[str] = field(default_factory=Counter)
    pmids: set[str] = field(default_factory=set)


def pct_delta(current: int | float, baseline: int | float) -> str:
    if baseline == 0:
        return "n/a"
    return f"{((current - baseline) / baseline) * 100.0:.2f}%"


def metric_value(manifest: dict[str, Any], key: str) -> int:
    return int((manifest.get("metrics") or {}).get(key, 0) or 0)


def metric_diff_rows(current: dict[str, Any], baseline: dict[str, Any]) -> list[dict[str, Any]]:
    keys = [
        "sentence_mention_count",
        "relation_candidate_count",
        "edge_support_count",
        "supported_existing_edge_count",
        "novel_candidate_count",
        "conflict_candidate_count",
    ]
    rows = []
    for key in keys:
        cur = metric_value(current, key)
        base = metric_value(baseline, key)
        rows.append({"metric": key, "baseline": base, "current": cur, "delta": cur - base, "delta_pct": pct_delta(cur, base)})
    return rows


def filter_surface_sets(filters: dict[str, Any]) -> tuple[set[tuple[str, str, str]], set[tuple[str, str, str]]]:
    blocked = set()
    downweighted = set()
    for rule in filters.get("surface_blocklist", []) or []:
        surface = normalized_surface_key(rule.get("surface", ""))
        for entity_type in rule.get("entity_types") or ["*"]:
            matched_fields = rule.get("matched_fields") or ["*"]
            for matched_field in matched_fields:
                blocked.add((str(entity_type), str(matched_field), surface))
    for rule in filters.get("surface_downweight", []) or []:
        surface = normalized_surface_key(rule.get("surface", ""))
        for entity_type in rule.get("entity_types") or ["*"]:
            matched_fields = rule.get("matched_fields") or ["*"]
            for matched_field in matched_fields:
                downweighted.add((str(entity_type), str(matched_field), surface))
    return blocked, downweighted


def matches_rule_tuple(rule_set: set[tuple[str, str, str]], entity_type: str, matched_field: str, surface: str) -> bool:
    return any(
        rule_entity in {"*", entity_type}
        and rule_field in {"*", matched_field}
        and rule_surface == surface
        for rule_entity, rule_field, rule_surface in rule_set
    )


def surface_risk_flags(agg: SurfaceAgg, filters: dict[str, Any], blocked: set[tuple[str, str, str]], downweighted: set[tuple[str, str, str]]) -> list[str]:
    flags = []
    surface = agg.normalized_surface
    tokens = surface.split()
    if matches_rule_tuple(blocked, agg.entity_type, agg.matched_field, surface):
        flags.append("blocklisted")
    if matches_rule_tuple(downweighted, agg.entity_type, agg.matched_field, surface):
        flags.append("surface_downweighted")
    if any(
        rule_matches(rule, "", agg.entity_type, agg.matched_field)
        for rule in filters.get("matched_field_downweight", []) or []
    ):
        flags.append("field_downweighted")
    if agg.matched_field in {"alias", "synonym"} and agg.mention_count >= DEFAULT_MIN_SURFACE_COUNT:
        flags.append("high_frequency_alias_or_synonym")
    if surface in GENERIC_SURFACES or surface in {"tumour", "tumours", "protein kinase", "kinase", "metabolic process", "signal transduction"}:
        flags.append("generic_surface")
    if len(tokens) == 1 and len(tokens[0]) <= 3 and agg.mention_count >= DEFAULT_MIN_SURFACE_COUNT:
        flags.append("short_high_frequency_surface")
    if any(term in surface for term in ("protein kinase", "kinase activity", "metabolic process", "signal transduction")):
        flags.append("generic_pathway_or_family_phrase")
    return flags


def load_surface_aggs(evidence_dir: Path, batch_size: int) -> list[SurfaceAgg]:
    if pq is None:
        raise RuntimeError("pyarrow is required for evidence precision QA.")
    path = evidence_dir / "sentence_mentions.parquet"
    aggs: dict[tuple[str, str, str], SurfaceAgg] = {}
    columns = ["entity_type", "normalized_surface", "matched_field", "entity_uid", "display_name", "pmid"]
    for batch in pq.ParquetFile(path).iter_batches(columns=columns, batch_size=batch_size):
        for row in batch.to_pylist():
            key = (str(row.get("entity_type") or ""), str(row.get("normalized_surface") or ""), str(row.get("matched_field") or ""))
            agg = aggs.setdefault(key, SurfaceAgg(*key))
            agg.mention_count += 1
            if row.get("entity_uid"):
                agg.entity_uids.add(str(row["entity_uid"]))
            if row.get("display_name"):
                agg.display_names[str(row["display_name"])] += 1
            if row.get("pmid"):
                agg.pmids.add(str(row["pmid"]))
    return list(aggs.values())


def high_risk_surfaces(evidence_dir: Path, filters: dict[str, Any], batch_size: int, limit: int) -> list[dict[str, Any]]:
    blocked, downweighted = filter_surface_sets(filters)
    rows = []
    for agg in load_surface_aggs(evidence_dir, batch_size):
        flags = surface_risk_flags(agg, filters, blocked, downweighted)
        if not flags:
            continue
        action = "block" if "blocklisted" in flags else "downweight" if ("surface_downweighted" in flags or "field_downweighted" in flags) else "review"
        rows.append(
            {
                "entity_type": agg.entity_type,
                "surface": agg.normalized_surface,
                "matched_field": agg.matched_field,
                "mention_count": agg.mention_count,
                "distinct_entity_count": len(agg.entity_uids),
                "sample_entities": [name for name, _count in agg.display_names.most_common(3)],
                "pmid_count": len(agg.pmids),
                "flags": flags,
                "action": action,
            }
        )
    rows.sort(key=lambda row: (0 if row["action"] == "block" else 1 if row["action"] == "downweight" else 2, -row["mention_count"], row["entity_type"], row["surface"]))
    return rows[:limit]


def markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    output = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    for row in rows:
        output.append("| " + " | ".join(str(item).replace("\n", " ") for item in row) + " |")
    return "\n".join(output)


def rule_rows(rules: list[dict[str, Any]], include_multiplier: bool = False) -> list[list[Any]]:
    rows = []
    for rule in rules:
        row = [
            normalized_surface_key(rule.get("surface", "")),
            ", ".join(rule.get("entity_types") or ["*"]),
            ", ".join(rule.get("matched_fields") or ["*"]),
        ]
        if include_multiplier:
            row.append(rule.get("confidence_multiplier", ""))
        row.append(rule.get("reason", ""))
        rows.append(row)
    return rows


def build_notes(args: argparse.Namespace) -> str:
    workspace = Path(args.workspace).resolve()
    evidence_dir = (workspace / args.evidence_root / args.release_id).resolve()
    current_manifest = read_json(evidence_dir / "literature_evidence_manifest.json")
    baseline_manifest = read_json(Path(args.baseline_manifest).resolve()) if args.baseline_manifest else current_manifest
    filters, filter_path = load_precision_filters(workspace, args.precision_filter_config)
    filter_hash = content_hash(filters)[:16]
    risk_rows = high_risk_surfaces(evidence_dir, filters, args.batch_size, args.limit)
    diff_rows = metric_diff_rows(current_manifest, baseline_manifest)
    baseline_supported = metric_value(baseline_manifest, "supported_existing_edge_count")
    current_supported = metric_value(current_manifest, "supported_existing_edge_count")
    retention = 1.0 if baseline_supported == 0 else current_supported / baseline_supported
    baseline_noise = metric_value(baseline_manifest, "novel_candidate_count") + metric_value(baseline_manifest, "conflict_candidate_count")
    current_noise = metric_value(current_manifest, "novel_candidate_count") + metric_value(current_manifest, "conflict_candidate_count")
    noise_delta = current_noise - baseline_noise
    precision_status = "passed" if retention >= args.min_supported_edge_retention else "review_required"

    lines = [
        "# Evidence Precision Notes",
        "",
        f"Release: `{args.release_id}`",
        "",
        f"Precision status: `{precision_status}`",
        "",
        f"Filter config: `{filter_path}`",
        f"Filter hash: `{filter_hash}`",
        "",
        "## Metrics Diff",
        "",
        markdown_table(
            ["metric", "baseline", "current", "delta", "delta_pct"],
            [[row["metric"], row["baseline"], row["current"], row["delta"], row["delta_pct"]] for row in diff_rows],
        ),
        "",
        "## Acceptance Signals",
        "",
        markdown_table(
            ["signal", "value"],
            [
                ["supported_existing_edge_retention", f"{retention:.3f}"],
                ["supported_existing_edge_retention_min", f"{args.min_supported_edge_retention:.3f}"],
                ["novel_plus_conflict_delta", noise_delta],
                ["evidence_validation_required", "phase15 validation must remain pass_with_known_blocks or passed"],
            ],
        ),
        "",
        "## Precision Filter Hits",
        "",
        "Manifest `precision_filters` captures exact block/downweight hit counts after rebuild.",
        "",
        "## Configured Blocklist",
        "",
        markdown_table(
            ["surface", "entity_types", "matched_fields", "reason"],
            rule_rows(filters.get("surface_blocklist", []) or []),
        ),
        "",
        "## Configured Downweight List",
        "",
        markdown_table(
            ["surface", "entity_types", "matched_fields", "multiplier", "reason"],
            rule_rows(filters.get("surface_downweight", []) or [], include_multiplier=True),
        ),
        "",
        "## High-Risk Surface QA",
        "",
        markdown_table(
            ["action", "entity_type", "surface", "field", "mentions", "pmids", "flags", "sample_entities"],
            [
                [
                    row["action"],
                    row["entity_type"],
                    row["surface"],
                    row["matched_field"],
                    row["mention_count"],
                    row["pmid_count"],
                    ", ".join(row["flags"]),
                    "; ".join(row["sample_entities"]),
                ]
                for row in risk_rows[: args.limit]
            ],
        ),
        "",
        "## Notes",
        "",
        "- Blocklist/downweight rules affect only sentence evidence lexicon construction.",
        "- Novel/conflict candidates remain in the evidence layer and do not enter graph scoring.",
        "- C3 acceptance requires comparable manifest metrics, no abnormal supported-edge collapse, and passing evidence validation.",
        "",
    ]
    return "\n".join(lines)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Evidence Precision QA notes.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--evidence-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--precision-filter-config", default=DEFAULT_PRECISION_FILTER_CONFIG)
    parser.add_argument("--baseline-manifest", default="")
    parser.add_argument("--output", default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--limit", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=100_000)
    parser.add_argument("--min-supported-edge-retention", type=float, default=DEFAULT_SUPPORTED_EDGE_RETENTION)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    notes = build_notes(args)
    output = Path(args.output)
    if not output.is_absolute():
        output = Path(args.workspace).resolve() / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(notes, encoding="utf-8")
    print(json.dumps({"output": str(output), "release_id": args.release_id}, ensure_ascii=False, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
