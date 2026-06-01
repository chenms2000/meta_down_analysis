"""Autofill high-confidence European trait identity review rows.

This script gathers conservative PubChem evidence and fills rows that can be
confirmed automatically. It marks safe exact-name/synonym matches as
`auto_reviewed` + `strict_identity`, and marks unsafe ratio/pool/X-code rows as
`auto_triaged` without promoting them to strict identity.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

try:
    import pyarrow.dataset as ds
except Exception:  # pragma: no cover
    ds = None

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from build_graph_projection import normalize_id_token  # noqa: E402


DEFAULT_REVIEW = "manual_sources/european_trait_identity/european_trait_identity_review.csv"
DEFAULT_CACHE = "raw_lake/European/pubchem_identity_autofill_cache.json"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
AUTO_REVIEWED = {"reviewed", "auto_reviewed", "accepted", "confirmed"}


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader], list(reader.fieldnames or [])


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_values(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    delimiter = "|" if "|" in text else ";"
    return [part.strip() for part in text.split(delimiter) if part.strip()]


def dedupe(values: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        text = re.sub(r"\s+", " ", str(value or "").strip())
        if not text:
            continue
        key = text.casefold()
        if key in seen:
            continue
        seen.add(key)
        out.append(text)
    return out


def normalize_name(value: Any) -> str:
    text = str(value or "").casefold()
    text = re.sub(r"\b(?:levels?|measurements?)\b", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def risky_identity_text(value: Any) -> bool:
    text = str(value or "")
    lowered = text.casefold()
    if re.search(r"\bX[- ]?\d+\b", text, flags=re.IGNORECASE):
        return True
    if re.search(r"\bratio\b|\s+to\s+|\s+\+\s+|/", text, flags=re.IGNORECASE):
        return True
    if re.search(r"\bor\b|\[[^\]]+\]|\(\d+\)", text, flags=re.IGNORECASE):
        return True
    return bool(
        re.search(
            r"\b(?:gpc|gpe|gpi|gps|gpg)\b|\d{1,2}:\d|sphingo|ceramide|stearoyl|oleoyl|"
            r"linoleoyl|palmitoyl|arachidon|acylcarnitine|carnitine|\bC\d",
            lowered,
        )
    )


def load_cache(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def save_cache(path: Path, cache: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cache, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def fetch_json(url: str, timeout: int) -> tuple[dict[str, Any] | None, str]:
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "metabo-european-identity-autofill/1.0"})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8")), ""
    except urllib.error.HTTPError as exc:
        return None, f"http_{exc.code}"
    except Exception as exc:
        return None, type(exc).__name__


def pubchem_name_cids(
    name: str,
    cache: dict[str, Any],
    timeout: int,
    sleep_seconds: float,
    offline: bool = False,
) -> dict[str, Any]:
    key = f"name::{name.casefold()}"
    if offline and key not in cache:
        return {"query": name, "url": "", "cids": [], "error": "offline_cache_miss"}
    if key not in cache:
        quoted = urllib.parse.quote(name)
        url = f"{PUBCHEM_BASE}/name/{quoted}/cids/JSON"
        payload, error = fetch_json(url, timeout)
        cids = [str(cid) for cid in (payload or {}).get("IdentifierList", {}).get("CID", [])[:5]]
        cache[key] = {"query": name, "url": url, "cids": cids, "error": error}
        if sleep_seconds:
            time.sleep(sleep_seconds)
    return cache[key]


def pubchem_cid_bundle(
    cid: str,
    cache: dict[str, Any],
    timeout: int,
    sleep_seconds: float,
    fetch_synonyms: bool = False,
    offline: bool = False,
) -> dict[str, Any]:
    cid = normalize_id_token(cid).removeprefix("pubchem:").removeprefix("cid:")
    key = f"cid_bundle::{cid}::synonyms_{int(bool(fetch_synonyms))}"
    if key in cache:
        return cache[key]
    if offline:
        return {
            "cid": cid,
            "property_url": "",
            "synonym_url": "",
            "properties": {},
            "synonyms": [],
            "property_error": "offline_cache_miss",
            "synonym_error": "offline_cache_miss",
        }
    prop_url = f"{PUBCHEM_BASE}/cid/{urllib.parse.quote(cid)}/property/Title,MolecularFormula,InChIKey,IUPACName/JSON"
    prop_payload, prop_error = fetch_json(prop_url, timeout)
    props = {}
    if prop_payload:
        rows = prop_payload.get("PropertyTable", {}).get("Properties", [])
        props = rows[0] if rows else {}
    synonyms = []
    syn_url = ""
    syn_error = "not_requested"
    if fetch_synonyms:
        syn_url = f"{PUBCHEM_BASE}/cid/{urllib.parse.quote(cid)}/synonyms/JSON"
        syn_payload, syn_error = fetch_json(syn_url, timeout)
        if syn_payload:
            rows = syn_payload.get("InformationList", {}).get("Information", [])
            synonyms = [str(value) for row in rows for value in row.get("Synonym", [])[:200]]
    bundle = {
        "cid": cid,
        "property_url": prop_url,
        "synonym_url": syn_url,
        "properties": props,
        "synonyms": synonyms,
        "property_error": prop_error,
        "synonym_error": syn_error,
    }
    cache[key] = bundle
    if sleep_seconds:
        time.sleep(sleep_seconds)
    return bundle


def local_pubchem_cids(compound_dir: Path) -> set[str]:
    path = compound_dir / "compound_identifier_index.parquet"
    if ds is None or not path.exists():
        return set()
    table = ds.dataset(str(path), format="parquet").to_table(columns=["namespace", "lookup_key"])
    out: set[str] = set()
    for namespace, lookup_key in zip(table.column("namespace").to_pylist(), table.column("lookup_key").to_pylist()):
        if str(namespace or "").upper() == "CID":
            out.add(normalize_id_token(lookup_key))
    return out


def candidate_names_for_row(row: dict[str, Any]) -> list[str]:
    return dedupe(
        [
            row.get("reviewed_name", ""),
            row.get("component_name", ""),
            *split_values(row.get("candidate_names", "")),
            *split_values(row.get("mapped_names", "")),
            row.get("reported_trait", ""),
        ]
    )


def pubchem_confirms_names(bundle: dict[str, Any], names: list[str]) -> tuple[bool, str]:
    props = bundle.get("properties") or {}
    reference_names = dedupe(
        [
            str(props.get("Title") or ""),
            str(props.get("IUPACName") or ""),
            *[str(value) for value in bundle.get("synonyms", [])],
        ]
    )
    reference_norms = {normalize_name(value) for value in reference_names if normalize_name(value)}
    for name in names:
        normalized = normalize_name(name)
        if normalized and normalized in reference_norms:
            return True, name
    return False, ""


def append_note(row: dict[str, Any], note: str) -> None:
    existing = str(row.get("notes") or "").strip()
    row["notes"] = f"{existing} | {note}" if existing else note


def triage_non_strict(row: dict[str, Any], scope: str, decision: str) -> None:
    if str(row.get("review_status") or "").strip().casefold() in AUTO_REVIEWED:
        return
    row["review_status"] = row.get("review_status") or "auto_triaged"
    if row["review_status"] == "needs_review":
        row["review_status"] = "auto_triaged"
    row["identity_scope"] = row.get("identity_scope") or scope
    row["autofill_decision"] = decision
    row["autofill_confidence"] = "triage"
    append_note(row, f"autofill: {decision}")


def maybe_autofill_strict(
    row: dict[str, Any],
    cache: dict[str, Any],
    local_cids: set[str],
    timeout: int,
    sleep_seconds: float,
    fetch_synonyms: bool,
    offline: bool,
) -> None:
    label = " ".join(
        str(row.get(field) or "")
        for field in ("reported_trait", "component_name", "candidate_names", "mapped_names", "suggested_identity_scope")
    )
    if risky_identity_text(label):
        triage_non_strict(row, row.get("identity_scope") or "class_or_pool", "risk_pattern_not_promoted")
        return
    names = candidate_names_for_row(row)
    cid = split_values(row.get("reviewed_pubchem_cid"))[:1]
    if not cid:
        for name in names[:5]:
            if not name or risky_identity_text(name):
                continue
            hit = pubchem_name_cids(name, cache, timeout, sleep_seconds, offline=offline)
            if hit.get("cids"):
                cid = [hit["cids"][0]]
                break
    if not cid:
        row["autofill_decision"] = "no_pubchem_cid_found"
        row["autofill_confidence"] = "none"
        return
    bundle = pubchem_cid_bundle(
        cid[0],
        cache,
        timeout,
        sleep_seconds,
        fetch_synonyms=fetch_synonyms,
        offline=offline,
    )
    props = bundle.get("properties") or {}
    confirmed, matched_name = pubchem_confirms_names(bundle, names)
    inchikey = str(props.get("InChIKey") or "").strip()
    title = str(props.get("Title") or "").strip()
    if not confirmed or not inchikey:
        row["autofill_decision"] = "pubchem_not_exact_enough"
        row["autofill_confidence"] = "low"
        append_note(row, "autofill: PubChem candidate found but synonym/InChIKey confirmation was insufficient")
        return
    row["review_status"] = "auto_reviewed"
    row["identity_scope"] = "strict_identity"
    row["reviewed_name"] = row.get("reviewed_name") or title or matched_name
    row["reviewed_pubchem_cid"] = cid[0]
    row["reviewed_inchikey"] = row.get("reviewed_inchikey") or inchikey
    row["evidence_source"] = (
        "PubChem PUG REST exact title/synonym match; auto rule strict_identity_v1"
    )
    row["evidence_url"] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid[0]}"
    row["autofill_decision"] = "auto_strict_identity"
    row["autofill_confidence"] = "high"
    row["local_cid_status"] = "local_cid_present" if normalize_id_token(cid[0]) in local_cids else "local_cid_absent"
    append_note(row, f"autofill: matched PubChem synonym/name '{matched_name}' to CID {cid[0]}")


def autofill_rows(args: argparse.Namespace) -> dict[str, Any]:
    workspace = Path(args.workspace).resolve()
    input_path = workspace / args.input
    output_path = workspace / (args.output or args.input)
    cache_path = workspace / args.cache
    rows, fieldnames = read_csv(input_path)
    cache = load_cache(cache_path)
    local_cids = local_pubchem_cids(workspace / args.compound_root / args.release_id)
    extra_fields = ["autofill_decision", "autofill_confidence", "local_cid_status"]
    fieldnames = [*fieldnames, *[field for field in extra_fields if field not in fieldnames]]
    processed = 0
    skipped_existing_decision = 0
    eligible_rows = [
        row
        for row in rows
        if not args.priority_max
        or int(str(row.get("priority") or "999").strip() or "999") <= args.priority_max
    ]
    for row in eligible_rows:
        if args.limit and processed >= args.limit:
            break
        status = str(row.get("review_status") or "").strip().casefold()
        if status in AUTO_REVIEWED and not args.force:
            continue
        if row.get("autofill_decision") and not args.force:
            skipped_existing_decision += 1
            continue
        suggested_scope = str(row.get("suggested_identity_scope") or "")
        bucket = str(row.get("review_bucket") or "")
        if suggested_scope == "ratio_component":
            triage_non_strict(row, "ratio_component", "ratio_component_not_strict")
        elif "class_or_pool" in suggested_scope or bucket in {"lipid_shorthand", "lipid_or_pool", "acylcarnitine_shorthand"}:
            triage_non_strict(row, "class_or_pool", "class_or_pool_not_strict")
        elif "unresolved" in suggested_scope or bucket == "platform_x_code":
            triage_non_strict(row, "unresolved", "platform_or_unresolved_not_strict")
        else:
            maybe_autofill_strict(
                row,
                cache,
                local_cids,
                args.timeout,
                args.sleep_seconds,
                args.fetch_synonyms,
                args.offline,
            )
        processed += 1
        if args.progress_every and processed % args.progress_every == 0:
            print(f"[autofill] processed {processed}/{len(eligible_rows)}", file=sys.stderr, flush=True)
        if args.cache_every and processed % args.cache_every == 0:
            save_cache(cache_path, cache)
            write_csv(output_path, rows, fieldnames)
    write_csv(output_path, rows, fieldnames)
    save_cache(cache_path, cache)
    return {
        "input": str(input_path),
        "output": str(output_path),
        "rows": len(rows),
        "eligible_rows": len(eligible_rows),
        "processed": processed,
        "skipped_existing_decision": skipped_existing_decision,
        "auto_strict_identity": sum(1 for row in rows if row.get("autofill_decision") == "auto_strict_identity"),
        "auto_triaged": sum(1 for row in rows if str(row.get("review_status") or "") == "auto_triaged"),
        "local_cid_present": sum(1 for row in rows if row.get("local_cid_status") == "local_cid_present"),
        "local_cid_absent": sum(1 for row in rows if row.get("local_cid_status") == "local_cid_absent"),
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Autofill safe European trait identity review rows.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--input", default=DEFAULT_REVIEW)
    parser.add_argument("--output", default="", help="Output path. Defaults to updating --input in place.")
    parser.add_argument("--cache", default=DEFAULT_CACHE)
    parser.add_argument("--compound-root", default="compound_match_index")
    parser.add_argument("--release-id", default="mvp_20260513T002254")
    parser.add_argument("--timeout", type=int, default=20)
    parser.add_argument("--sleep-seconds", type=float, default=0.08)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--priority-max", type=int, default=20, help="Only process rows with priority <= this value. Use 0 for all rows.")
    parser.add_argument("--progress-every", type=int, default=25)
    parser.add_argument("--cache-every", type=int, default=25)
    parser.add_argument("--fetch-synonyms", action="store_true", help="Fetch PubChem synonyms for deeper but slower confirmation.")
    parser.add_argument("--offline", action="store_true", help="Use cached PubChem responses only; never make network requests.")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    print(json.dumps(autofill_rows(parse_args(argv)), ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
