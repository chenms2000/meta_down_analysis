"""Export a fillable European trait identity curation template.

The template is intentionally conservative: it helps reviewers add audited
chemical identity evidence without forcing ratio, pool, lipid shorthand, or
platform-code traits into exact metabolite identities.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Any


DEFAULT_OUTPUT = "manual_sources/european_trait_identity/european_trait_identity_review.csv"

FIELDNAMES = [
    "priority",
    "accession_id",
    "reported_trait",
    "source_resolution_status",
    "review_bucket",
    "suggested_identity_scope",
    "suggested_action",
    "component_index",
    "component_name",
    "relationship_to_trait",
    "candidate_names",
    "mapped_names",
    "pubchem_cids",
    "pubchem_titles",
    "review_status",
    "identity_scope",
    "reviewed_name",
    "reviewed_pubchem_cid",
    "reviewed_hmdb_id",
    "reviewed_chebi_id",
    "reviewed_inchikey",
    "evidence_source",
    "evidence_url",
    "notes",
    "summary_statistics_url",
    "pubmed_id",
    "paper_title",
]


def split_values(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    delimiter = "|" if "|" in text else ";"
    return [part.strip() for part in text.split(delimiter) if part.strip()]


def clean_text(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip())


def looks_ratio_or_composite(text: str) -> bool:
    return bool(re.search(r"\s+\+\s+|/|\bto\b|\bratio\b", text, flags=re.IGNORECASE))


def looks_pool_or_lipid(text: str) -> bool:
    lowered = text.casefold()
    return bool(
        re.search(
            r"\b(?:gpc|gpe|gpi|gps|gpg)\b|\d{1,2}:\d|sphingo|ceramide|stearoyl|oleoyl|"
            r"linoleoyl|palmitoyl|arachidon|acylcarnitine|carnitine|\bC\d",
            lowered,
        )
    )


def review_bucket_for_annotation(row: dict[str, Any], unresolved_review: dict[str, dict[str, Any]]) -> str:
    accession = row.get("accession_id", "")
    if accession in unresolved_review:
        return unresolved_review[accession].get("review_bucket", "")
    label = " ".join(
        value
        for value in [
            row.get("reported_trait", ""),
            row.get("candidate_names", ""),
            row.get("mapped_names", ""),
        ]
        if value
    )
    if re.search(r"\bX[- ]?\d+\b", label, flags=re.IGNORECASE):
        return "platform_x_code"
    if looks_ratio_or_composite(label):
        return "composite_or_ratio"
    if looks_pool_or_lipid(label):
        return "lipid_or_pool"
    if row.get("resolution_status") == "pubchem_only":
        return "pubchem_only"
    if row.get("resolution_status") == "local_name":
        return "local_name_needs_identity_confirmation"
    return row.get("resolution_status", "") or "manual_review"


def recommendation(row: dict[str, Any], bucket: str) -> tuple[int, str, str, str]:
    status = row.get("resolution_status", "")
    trait = row.get("reported_trait", "")
    label = f"{trait} {row.get('candidate_names', '')} {row.get('mapped_names', '')}"
    if bucket == "composite_or_ratio" or looks_ratio_or_composite(label):
        return (
            40,
            "ratio_component",
            "Keep the GCST trait as a ratio/composite; review each component separately and do not mark the whole trait as strict_identity.",
            "ratio_component",
        )
    if bucket in {"lipid_shorthand", "lipid_or_pool", "acylcarnitine_shorthand"} or looks_pool_or_lipid(label):
        return (
            50,
            "class_or_pool_unless_inchikey_confirms",
            "Use class_or_pool unless a source supplement, LIPID MAPS/HMDB entry, or full InChIKey confirms one exact structure.",
            "pool_or_isomer_candidate",
        )
    if status == "pubchem_only":
        return (
            10,
            "strict_identity_if_exact",
            "Open the PubChem CID/title; if it is the exact compound and not a ratio/isomer/pool, set review_status=reviewed and identity_scope=strict_identity.",
            "candidate_compound",
        )
    if bucket == "alias_variant_missing":
        return (
            20,
            "strict_identity_if_external_id_confirms",
            "Find HMDB/PubChem/ChEBI/InChIKey for the exact synonym/acid/salt form; fill one stable ID and evidence_source.",
            "candidate_compound",
        )
    if bucket in {"biochemical_abbrev", "local_name_needs_identity_confirmation"}:
        return (
            30,
            "strict_identity_if_unique",
            "Confirm the local/common name against HMDB/PubChem/ChEBI or InChIKey; fill a stable ID only when unique.",
            "candidate_compound",
        )
    if bucket == "platform_x_code":
        return (
            90,
            "unresolved_until_vendor_mapping",
            "Find the original study supplement/vendor mapping for the X-code; do not infer identity from context.",
            "platform_code",
        )
    return (
        80,
        "manual_review",
        "Leave review_status=needs_review until an external identifier or curated source confirms the identity.",
        "candidate_compound",
    )


def load_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def load_unresolved_review(path: Path) -> dict[str, dict[str, str]]:
    return {row.get("accession_id", ""): row for row in load_csv(path) if row.get("accession_id")}


def template_rows(annotation_rows: list[dict[str, str]], unresolved_review: dict[str, dict[str, str]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in annotation_rows:
        status = row.get("resolution_status", "")
        if status == "pubchem_cid_in_local_index":
            continue
        bucket = review_bucket_for_annotation(row, unresolved_review)
        priority, suggested_scope, action, relationship = recommendation(row, bucket)
        candidates = split_values(row.get("candidate_names")) or split_values(row.get("mapped_names")) or [row.get("reported_trait", "")]
        if suggested_scope == "ratio_component":
            component_names = candidates[:8]
        else:
            component_names = [row.get("pubchem_titles", "") or candidates[0]]
        pubchem_cids = split_values(row.get("pubchem_cids"))
        pubchem_titles = split_values(row.get("pubchem_titles"))
        for index, component in enumerate(component_names, start=1):
            reviewed_pubchem_cid = pubchem_cids[0] if status == "pubchem_only" and pubchem_cids else ""
            reviewed_name = pubchem_titles[0] if status == "pubchem_only" and pubchem_titles else clean_text(component)
            rows.append(
                {
                    "priority": priority,
                    "accession_id": row.get("accession_id", ""),
                    "reported_trait": row.get("reported_trait", ""),
                    "source_resolution_status": status,
                    "review_bucket": bucket,
                    "suggested_identity_scope": suggested_scope,
                    "suggested_action": action,
                    "component_index": index if len(component_names) > 1 else "",
                    "component_name": clean_text(component),
                    "relationship_to_trait": relationship,
                    "candidate_names": row.get("candidate_names", ""),
                    "mapped_names": row.get("mapped_names", ""),
                    "pubchem_cids": row.get("pubchem_cids", ""),
                    "pubchem_titles": row.get("pubchem_titles", ""),
                    "review_status": "needs_review",
                    "identity_scope": "ratio_component" if suggested_scope == "ratio_component" else "",
                    "reviewed_name": reviewed_name,
                    "reviewed_pubchem_cid": reviewed_pubchem_cid,
                    "reviewed_hmdb_id": "",
                    "reviewed_chebi_id": "",
                    "reviewed_inchikey": "",
                    "evidence_source": "",
                    "evidence_url": row.get("summary_statistics_url", ""),
                    "notes": "",
                    "summary_statistics_url": row.get("summary_statistics_url", ""),
                    "pubmed_id": row.get("pubmed_id", ""),
                    "paper_title": row.get("paper_title", ""),
                }
            )
    rows.sort(key=lambda item: (int(item["priority"]), item["accession_id"], str(item["component_index"])))
    return rows


def export_template(args: argparse.Namespace) -> Path:
    workspace = Path(args.workspace).resolve()
    annotations = load_csv(workspace / args.annotations)
    unresolved = load_unresolved_review(workspace / args.unresolved_review)
    rows = template_rows(annotations, unresolved)
    if args.limit and args.limit > 0:
        rows = rows[: args.limit]
    output = workspace / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    return output


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export a fillable European trait identity review template.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--annotations", default="raw_lake/European/European_trait_annotations.csv")
    parser.add_argument("--unresolved-review", default="raw_lake/European/European_unresolved_review.csv")
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--limit", type=int, default=0, help="Optional row limit; 0 exports all actionable rows.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    output = export_template(parse_args(argv))
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
