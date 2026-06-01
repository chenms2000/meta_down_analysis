"""Audit compound name surfaces for ambiguity and high-risk normalization drift.

The audit is intentionally conservative:
- Local mode scans frozen compound indexes for many-to-one / one-to-many names,
  delimiter-sensitive labels, lipid shorthand, and under-specified isomers.
- Optional online mode checks only the highest-risk names against official
  PubChem PUG REST and EBI OLS ChEBI search endpoints.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pyarrow.dataset as ds


DEFAULT_RELEASE_ID = "mvp_20260513T002254"
AUTHORITY_PREFIXES = ("CHEBI", "HMDB", "PUBCHEM", "CID", "KEGG", "LIPIDMAPS")
STEREO_OR_ISOMER_REQUIRED = {
    "aconitate": "Specify cis-aconitic acid or trans-aconitic acid.",
    "2-hydroxyglutarate": "Specify D/R or L/S 2-hydroxyglutaric acid when the assay supports it.",
    "2 hydroxyglutarate": "Specify D/R or L/S 2-hydroxyglutaric acid when the assay supports it.",
    "xanthine": "Specify a tautomer such as 7H-xanthine when the annotation source does so.",
    "sorbitol": "Prefer D-glucitol/D-sorbitol when the assay annotation is stereospecific.",
    "butyrylcarnitine": "Prefer O-butanoyl-L-carnitine / butyryl-L-carnitine.",
    "free carnitine": "Prefer L-carnitine or a stable identifier.",
    "asymmetric dimethylarginine": "Prefer N(omega),N(omega)-dimethyl-L-arginine or a stable identifier.",
}
BROAD_OR_SHORT_NAMES = {
    "acid",
    "base",
    "salt",
    "unknown",
    "lipid",
    "ceramide",
    "hexose",
    "pentose",
    "sugar",
    "protein",
    "peptide",
}
LIPID_SHORTHAND_RE = re.compile(
    r"(?i)(?:\b(?:[a-z]{1,4}\s*)?\d{1,2}:\d(?:/\d{1,2}:\d)?\b|ceramide\s+d\d{1,2}:\d/\d{1,2}:\d|d\d{1,2}:\d/\d{1,2}:\d)"
)


def normalize_key(value: Any) -> str:
    text = str(value or "").strip().casefold()
    text = re.sub(r"[\u2010-\u2015]", "-", text)
    text = re.sub(r"\s+", " ", text)
    return text


def safe_int(value: Any) -> int:
    try:
        return int(value)
    except Exception:
        return 0


def risk_reasons(raw_value: str, name_key: str, match_fields: set[str], min_rank: float, uid_count: int) -> list[str]:
    reasons: list[str] = []
    compact_key = normalize_key(name_key).replace("-", " ")
    if uid_count > 1:
        reasons.append("name_maps_to_multiple_metabolites")
    if "," in raw_value:
        reasons.append("delimiter_sensitive_name")
    if LIPID_SHORTHAND_RE.search(raw_value):
        reasons.append("lipid_shorthand_or_chain_notation")
    if compact_key in STEREO_OR_ISOMER_REQUIRED:
        reasons.append("stereo_or_isomer_under_specified")
    if compact_key in BROAD_OR_SHORT_NAMES or len(re.sub(r"[^a-z0-9]", "", compact_key)) <= 3:
        reasons.append("broad_or_short_surface")
    if "canonical_name" not in match_fields and min_rank <= 70:
        reasons.append("synonym_only_low_rank_surface")
    return reasons


def risk_score(reasons: list[str], uid_count: int, missing_authority: bool = False) -> int:
    score = 0
    weights = {
        "name_maps_to_multiple_metabolites": 35,
        "stereo_or_isomer_under_specified": 35,
        "lipid_shorthand_or_chain_notation": 30,
        "delimiter_sensitive_name": 20,
        "broad_or_short_surface": 20,
        "synonym_only_low_rank_surface": 10,
    }
    for reason in reasons:
        score += weights.get(reason, 5)
    if uid_count > 2:
        score += min(30, uid_count * 3)
    if missing_authority:
        score += 15
    return score


def scan_name_index(compound_dir: Path) -> dict[str, dict[str, Any]]:
    name_path = compound_dir / "compound_name_index.parquet"
    if not name_path.exists():
        raise FileNotFoundError(f"Missing {name_path}")

    stats: dict[str, dict[str, Any]] = {}
    dataset = ds.dataset(name_path, format="parquet")
    scanner = dataset.scanner(columns=["name_key", "raw_value", "metabolite_uid", "match_field", "rank"], batch_size=250_000)
    for batch in scanner.to_batches():
        for row in batch.to_pylist():
            name_key = str(row.get("name_key") or "").strip()
            if not name_key:
                continue
            stat = stats.get(name_key)
            uid = str(row.get("metabolite_uid") or "")
            raw_value = str(row.get("raw_value") or name_key)
            rank = float(row.get("rank") or 0.0)
            if stat is None:
                stat = {
                    "name_key": name_key,
                    "raw_examples": [],
                    "match_fields": set(),
                    "min_rank": rank,
                    "max_rank": rank,
                    "first_uid": uid,
                    "uid_set": None,
                    "row_count": 0,
                }
                stats[name_key] = stat
            stat["row_count"] += 1
            if raw_value and raw_value not in stat["raw_examples"] and len(stat["raw_examples"]) < 5:
                stat["raw_examples"].append(raw_value)
            stat["match_fields"].add(str(row.get("match_field") or ""))
            stat["min_rank"] = min(float(stat["min_rank"]), rank)
            stat["max_rank"] = max(float(stat["max_rank"]), rank)
            if uid and uid != stat["first_uid"]:
                if stat["uid_set"] is None:
                    stat["uid_set"] = {stat["first_uid"], uid}
                else:
                    stat["uid_set"].add(uid)
    return stats


def action_recommendation(reasons: list[str], name_key: str) -> dict[str, str]:
    normalized = normalize_key(name_key).replace("-", " ")
    if normalized in STEREO_OR_ISOMER_REQUIRED:
        return {
            "action": "require_specific_isomer_or_identifier",
            "recommendation": STEREO_OR_ISOMER_REQUIRED[normalized],
        }
    if "name_maps_to_multiple_metabolites" in reasons and "lipid_shorthand_or_chain_notation" in reasons:
        return {
            "action": "require_lipid_identifier",
            "recommendation": "Do not auto-map this lipid shorthand; provide LIPID MAPS/HMDB/ChEBI/InChIKey or exact structural name.",
        }
    if "name_maps_to_multiple_metabolites" in reasons:
        return {
            "action": "manual_review_required",
            "recommendation": "Surface maps to multiple frozen metabolites; keep in review queue unless an identifier or orthogonal feature resolves it.",
        }
    if "delimiter_sensitive_name" in reasons:
        return {
            "action": "quote_delimiter_sensitive_name",
            "recommendation": "Quote this name in CSV/TSV exports or use JSON records to avoid delimiter splitting.",
        }
    if "lipid_shorthand_or_chain_notation" in reasons:
        return {
            "action": "prefer_lipid_identifier",
            "recommendation": "Prefer LIPID MAPS/HMDB/InChIKey for chain-resolved lipid notation.",
        }
    if "broad_or_short_surface" in reasons:
        return {
            "action": "avoid_broad_surface",
            "recommendation": "Use a more specific metabolite name or stable database identifier.",
        }
    return {
        "action": "review_low_precision_surface",
        "recommendation": "Review before promoting this surface to a deterministic alias.",
    }


def collect_candidate_rows(stats: dict[str, dict[str, Any]], limit: int | None) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for stat in stats.values():
        uid_set = stat["uid_set"]
        uid_count = len(uid_set) if uid_set is not None else 1
        reasons = risk_reasons(
            stat["raw_examples"][0] if stat["raw_examples"] else stat["name_key"],
            stat["name_key"],
            stat["match_fields"],
            float(stat["min_rank"]),
            uid_count,
        )
        if not reasons:
            continue
        action = action_recommendation(reasons, stat["name_key"])
        rows.append(
            {
                "name_key": stat["name_key"],
                "raw_examples": stat["raw_examples"],
                "uid_count": uid_count,
                "uids": sorted(uid_set)[:20] if uid_set else [stat["first_uid"]],
                "row_count": int(stat["row_count"]),
                "match_fields": sorted(stat["match_fields"]),
                "min_rank": round(float(stat["min_rank"]), 3),
                "max_rank": round(float(stat["max_rank"]), 3),
                "risk_reasons": reasons,
                "risk_score": risk_score(reasons, uid_count),
                "recommended_action": action["action"],
                "recommendation": action["recommendation"],
            }
        )
    rows.sort(key=lambda row: (-safe_int(row["risk_score"]), -safe_int(row["uid_count"]), row["name_key"]))
    if limit is None:
        return rows
    return rows[:limit]


def annotate_authorities(rows: list[dict[str, Any]], normalized_dir: Path) -> None:
    wanted = {uid for row in rows for uid in row.get("uids", [])}
    if not wanted:
        return
    dataset = ds.dataset(normalized_dir / "metabolites.parquet", format="parquet")
    uid_to_xrefs: dict[str, list[str]] = {}
    uid_to_name: dict[str, str] = {}
    for batch in dataset.scanner(columns=["metabolite_uid", "canonical_name", "external_xrefs"], batch_size=100_000).to_batches():
        for row in batch.to_pylist():
            uid = str(row.get("metabolite_uid") or "")
            if uid not in wanted:
                continue
            uid_to_xrefs[uid] = [str(x) for x in (row.get("external_xrefs") or [])]
            uid_to_name[uid] = str(row.get("canonical_name") or "")
    for row in rows:
        authorities: dict[str, int] = defaultdict(int)
        canonical_examples: list[str] = []
        for uid in row.get("uids", []):
            if uid_to_name.get(uid) and uid_to_name[uid] not in canonical_examples:
                canonical_examples.append(uid_to_name[uid])
            for xref in uid_to_xrefs.get(uid, []):
                prefix = xref.split(":", 1)[0].upper()
                if prefix in AUTHORITY_PREFIXES or prefix.startswith("CHEBI") or prefix.startswith("KEGG") or prefix.startswith("PUBCHEM"):
                    authorities[prefix] += 1
        row["authority_prefixes"] = dict(sorted(authorities.items()))
        row["canonical_examples"] = canonical_examples[:5]
        if not authorities:
            row["risk_reasons"].append("no_local_authority_xref")
            row["risk_score"] = risk_score(row["risk_reasons"], safe_int(row["uid_count"]), missing_authority=True)


def fetch_json(url: str, timeout: int) -> tuple[dict[str, Any] | None, str]:
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "metabo-name-audit/1.0"})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8")), ""
    except urllib.error.HTTPError as exc:
        return None, f"http_{exc.code}"
    except Exception as exc:
        return None, type(exc).__name__


def online_verify_name(name: str, timeout: int = 15) -> dict[str, Any]:
    quoted = urllib.parse.quote(name)
    pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quoted}/cids/JSON"
    pubchem_payload, pubchem_error = fetch_json(pubchem_url, timeout)
    cids = []
    if pubchem_payload:
        cids = pubchem_payload.get("IdentifierList", {}).get("CID", [])[:10]

    ols_url = f"https://www.ebi.ac.uk/ols4/api/search?q={quoted}&ontology=chebi&rows=5"
    ols_payload, ols_error = fetch_json(ols_url, timeout)
    chebi_hits = []
    if ols_payload:
        for doc in ols_payload.get("response", {}).get("docs", [])[:5]:
            chebi_hits.append(
                {
                    "label": doc.get("label", ""),
                    "obo_id": doc.get("obo_id", ""),
                    "iri": doc.get("iri", ""),
                }
            )
    return {
        "query": name,
        "pubchem": {"url": pubchem_url, "cids": cids, "error": pubchem_error},
        "chebi_ols": {"url": ols_url, "hits": chebi_hits, "error": ols_error},
    }


def write_outputs(rows: list[dict[str, Any]], output_dir: Path, release_id: str, online: list[dict[str, Any]]) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "audit_version": "compound_name_audit.v1",
        "release_id": release_id,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "risk_row_count": len(rows),
        "online_check_count": len(online),
        "online_sources": [
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON",
            "https://www.ebi.ac.uk/ols4/api/search?q={name}&ontology=chebi",
        ],
        "summary": {
            "by_recommended_action": dict(
                sorted(
                    {
                        action: sum(1 for row in rows if row.get("recommended_action") == action)
                        for action in sorted({str(row.get("recommended_action") or "") for row in rows})
                    }.items()
                )
            ),
            "by_reason": dict(
                sorted(
                    {
                        reason: sum(1 for row in rows if reason in row.get("risk_reasons", []))
                        for reason in sorted({reason for row in rows for reason in row.get("risk_reasons", [])})
                    }.items()
                )
            ),
            "top_risk_names": rows[:25],
        },
        "online_checks": online,
    }
    report_path = output_dir / "compound_name_audit_report.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")

    csv_path = output_dir / "compound_name_audit_candidates.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "risk_score",
                "name_key",
                "recommended_action",
                "recommendation",
                "raw_examples",
                "risk_reasons",
                "uid_count",
                "canonical_examples",
                "authority_prefixes",
                "match_fields",
                "min_rank",
                "max_rank",
                "row_count",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "risk_score": row.get("risk_score", 0),
                    "name_key": row.get("name_key", ""),
                    "recommended_action": row.get("recommended_action", ""),
                    "recommendation": row.get("recommendation", ""),
                    "raw_examples": " | ".join(row.get("raw_examples", [])),
                    "risk_reasons": ";".join(row.get("risk_reasons", [])),
                    "uid_count": row.get("uid_count", 0),
                    "canonical_examples": " | ".join(row.get("canonical_examples", [])),
                    "authority_prefixes": json.dumps(row.get("authority_prefixes", {}), ensure_ascii=False, sort_keys=True),
                    "match_fields": ";".join(row.get("match_fields", [])),
                    "min_rank": row.get("min_rank", ""),
                    "max_rank": row.get("max_rank", ""),
                    "row_count": row.get("row_count", ""),
                }
            )
    suggestions_path = output_dir / "compound_name_cleaning_suggestions.csv"
    with suggestions_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["name_key", "recommended_action", "recommendation", "risk_reasons", "risk_score", "canonical_examples"],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "name_key": row.get("name_key", ""),
                    "recommended_action": row.get("recommended_action", ""),
                    "recommendation": row.get("recommendation", ""),
                    "risk_reasons": ";".join(row.get("risk_reasons", [])),
                    "risk_score": row.get("risk_score", 0),
                    "canonical_examples": " | ".join(row.get("canonical_examples", [])),
                }
            )
    return report_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit compound names for ambiguity and high-risk normalization surfaces.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", default=DEFAULT_RELEASE_ID)
    parser.add_argument("--limit", type=int, default=1000, help="Maximum local risk rows to write. Use 0 for all risk rows.")
    parser.add_argument("--online-limit", type=int, default=0, help="Verify this many top-risk names online.")
    parser.add_argument("--online-timeout", type=int, default=15)
    parser.add_argument("--output-dir", default="")
    args = parser.parse_args()

    workspace = Path(args.workspace).resolve()
    release_id = args.release_id
    compound_dir = workspace / "compound_match_index" / release_id
    normalized_dir = workspace / "normalized_store" / release_id
    output_dir = Path(args.output_dir).resolve() if args.output_dir else workspace / "validation_reports" / release_id / "name_audit"

    stats = scan_name_index(compound_dir)
    row_limit = None if args.limit == 0 else max(1, args.limit)
    rows = collect_candidate_rows(stats, row_limit)
    annotate_authorities(rows, normalized_dir)
    rows.sort(key=lambda row: (-safe_int(row["risk_score"]), -safe_int(row["uid_count"]), row["name_key"]))

    online_checks = []
    for row in rows[: max(0, args.online_limit)]:
        query = row.get("raw_examples", [row["name_key"]])[0]
        online_checks.append(online_verify_name(query, timeout=args.online_timeout))
        time.sleep(0.2)

    report_path = write_outputs(rows, output_dir, release_id, online_checks)
    print(
        json.dumps(
            {
                "release_id": release_id,
                "risk_row_count": len(rows),
                "unique_name_key_count": len(stats),
                "online_check_count": len(online_checks),
                "report_path": str(report_path),
                "csv_path": str(output_dir / "compound_name_audit_candidates.csv"),
                "suggestions_path": str(output_dir / "compound_name_cleaning_suggestions.csv"),
                "top_risk_names": [row["name_key"] for row in rows[:10]],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
