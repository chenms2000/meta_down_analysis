"""Build prediction overlay CSVs from external drug and cell-context sources.

The overlay remains separate from canonical graph facts. Rows produced here are
research-prioritization evidence used by metabo_service prediction_pack.v1.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any

try:
    import pyarrow.dataset as ds
except Exception:  # pragma: no cover - handled at runtime.
    ds = None

from metabo_service import (  # noqa: E402
    DEFAULT_GRAPH_ROOT,
    DEFAULT_PREDICTION_OVERLAY_ROOT,
    MetaboService,
    content_hash,
    load_context_arg,
    load_records_arg,
    normalize_prediction_text,
    prediction_terms_from_value,
)


CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
DRUGCENTRAL_TARGET_URL = "https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz"
DGIDB_GRAPHQL_URL = "https://dgidb.org/api/graphql"
OPENTARGETS_2603_BASE = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output"
CELLOSAURUS_TXT_URL = "https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt"
CELLOSAURUS_API_SEARCH_URL = ""
DEFAULT_HTTP_TIMEOUT_SECONDS = 20
DEFAULT_HTTP_RETRIES = 1

FIELDNAMES_DRUG = [
    "drug_id",
    "drug_name",
    "target_uid",
    "target_symbol",
    "target_name",
    "target_id",
    "mechanism",
    "confidence",
    "context_terms",
    "cell_state",
    "source_name",
    "source_record_id",
    "evidence_level",
    "license_id",
    "source_url",
    "source_version",
]

FIELDNAMES_CELL = [
    "cell_type_id",
    "cell_type_name",
    "cell_line",
    "target_uids",
    "target_symbols",
    "pathway_uids",
    "pathway_terms",
    "metabolite_terms",
    "context_terms",
    "cell_state",
    "confidence",
    "source_name",
    "source_record_id",
    "evidence_level",
    "license_id",
    "source_url",
    "source_version",
]

GENE_SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9._-]{1,19}$")
GENE_COLUMN_ID_RE = re.compile(r"\s*\([^)]*\)\s*$")
MODEL_ID_COLUMNS = ("ModelID", "DepMap_ID", "model_id", "model", "id", "cell_line_id")
CELL_NAME_COLUMNS = ("CellLineName", "stripped_cell_line_name", "CCLEName", "cell_line", "cell_line_name", "name")
MARKER_CELL_ID_COLUMNS = ("cell_type_id", "cell_id", "cluster_id", "annotation_id", "cluster", "leiden")
MARKER_CELL_NAME_COLUMNS = ("cell_type_name", "cell_type", "annotation", "cluster_name", "name")
MARKER_GENE_COLUMNS = ("target_symbol", "target_symbols", "marker_gene", "marker_genes", "gene", "genes", "gene_symbol", "gene_symbols")
MARKER_PATHWAY_COLUMNS = ("pathway_terms", "pathways", "pathway")
MARKER_METABOLITE_COLUMNS = ("metabolite_terms", "metabolites", "metabolite")
MARKER_CONTEXT_COLUMNS = (
    "context_terms",
    "tissue",
    "tissue_general",
    "disease",
    "organ",
    "assay",
    "dataset",
    "dataset_id",
    "cell_state",
    "state",
    "phenotype",
)


def clamp_unit(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        number = default
    if math.isnan(number) or math.isinf(number):
        number = default
    return min(1.0, max(0.0, number))


def request_json(url: str, timeout: int = DEFAULT_HTTP_TIMEOUT_SECONDS, retries: int = DEFAULT_HTTP_RETRIES) -> dict[str, Any]:
    headers = {"User-Agent": "metabo-data-end/0.1 prediction-overlay-builder"}
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as exc:  # pragma: no cover - network variability.
            last_error = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Failed to fetch JSON from {url}: {last_error}")


def request_json_post(
    url: str,
    payload: dict[str, Any],
    timeout: int = DEFAULT_HTTP_TIMEOUT_SECONDS,
    retries: int = DEFAULT_HTTP_RETRIES,
) -> dict[str, Any]:
    headers = {
        "User-Agent": "metabo-data-end/0.1 prediction-overlay-builder",
        "Content-Type": "application/json",
    }
    data = json.dumps(payload).encode("utf-8")
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                parsed = json.loads(response.read().decode("utf-8"))
            if parsed.get("errors"):
                raise RuntimeError(json.dumps(parsed["errors"], ensure_ascii=False)[:1000])
            return parsed
        except Exception as exc:  # pragma: no cover - network variability.
            last_error = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Failed to POST JSON to {url}: {last_error}")


def download_file(url: str, path: Path, timeout: int = 180, retries: int = 3) -> Path:
    if path.exists() and path.stat().st_size > 0:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".part")
    headers = {"User-Agent": "metabo-data-end/0.1 prediction-overlay-builder"}
    last_error: Exception | None = None
    for attempt in range(retries):
        if tmp_path.exists():
            tmp_path.unlink()
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response, tmp_path.open("wb") as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
            tmp_path.replace(path)
            return path
        except Exception as exc:
            last_error = exc
            if tmp_path.exists():
                tmp_path.unlink()
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"Failed to download {url}: {last_error}")
    return path


def list_index_files(url: str, include_regex: str) -> list[str]:
    payload = urllib.request.urlopen(url, timeout=DEFAULT_HTTP_TIMEOUT_SECONDS).read().decode("utf-8", errors="replace")
    pattern = re.compile(include_regex)
    files = []
    for match in re.finditer(r'href="([^"]+)"', payload):
        href = urllib.parse.unquote(match.group(1))
        name = href.rstrip("/")
        if "/" in name:
            name = name.rsplit("/", 1)[-1]
        if pattern.search(name):
            files.append(name)
    return sorted(set(files))


def download_index_files(url: str, output_dir: Path, include_regex: str, timeout: int = 240) -> list[Path]:
    paths: list[Path] = []
    for name in list_index_files(url, include_regex):
        if name == "_SUCCESS":
            continue
        paths.append(download_file(urllib.parse.urljoin(url, name), output_dir / name, timeout=timeout, retries=3))
    return paths


def file_sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def source_file_record(source: str, path: Path, url: str) -> dict[str, Any] | None:
    if not path.exists() or not path.is_file():
        return None
    return {
        "source": source,
        "path": str(path),
        "source_url": url,
        "bytes": path.stat().st_size,
        "sha256": file_sha256(path),
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def merge_terms(*values: Any) -> str:
    terms: list[str] = []
    seen: set[str] = set()
    for value in values:
        for term in prediction_terms_from_value(value):
            if term not in seen:
                terms.append(term)
                seen.add(term)
    return ";".join(terms)


def first_value(row: dict[str, Any], names: tuple[str, ...] | list[str], default: str = "") -> str:
    for name in names:
        value = row.get(name)
        if value is not None and str(value).strip() != "":
            return str(value).strip()
    lower_lookup = {str(key).casefold(): value for key, value in row.items()}
    for name in names:
        value = lower_lookup.get(str(name).casefold())
        if value is not None and str(value).strip() != "":
            return str(value).strip()
    return default


def normalize_gene_symbol(value: Any) -> str:
    text = GENE_COLUMN_ID_RE.sub("", str(value or "").strip())
    text = text.split("|", 1)[0].strip()
    text = text.replace(" ", "")
    if not text:
        return ""
    text = text.upper()
    return text if GENE_SYMBOL_RE.match(text) else ""


def split_gene_symbols(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        pieces = [str(item) for item in value]
    else:
        pieces = re.split(r"[;|,\s]+", str(value))
    symbols: list[str] = []
    seen: set[str] = set()
    for piece in pieces:
        symbol = normalize_gene_symbol(piece)
        if symbol and symbol not in seen:
            symbols.append(symbol)
            seen.add(symbol)
    return symbols


def target_uid_lookup_from_targets(targets: list[dict[str, Any]]) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for target in targets:
        uid = str(target.get("node_uid") or target.get("target_uid") or "")
        if not uid:
            continue
        values = [
            target.get("canonical_name"),
            target.get("display_name"),
            target.get("primary_external_id"),
            *extract_xrefs(target, "HGNC"),
            *extract_xrefs(target, "UNIPROT_SWISSPROT"),
        ]
        for value in values:
            symbol = normalize_gene_symbol(value)
            if symbol:
                lookup.setdefault(symbol, uid)
    return lookup


def target_uids_for_symbols(symbols: list[str], target_uid_by_symbol: dict[str, str]) -> str:
    uids: list[str] = []
    seen: set[str] = set()
    for symbol in symbols:
        uid = target_uid_by_symbol.get(normalize_gene_symbol(symbol), "")
        if uid and uid not in seen:
            uids.append(uid)
            seen.add(uid)
    return ";".join(uids)


def parse_score(value: Any) -> float | None:
    try:
        score = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(score) or math.isinf(score):
        return None
    return score


def select_depmap_gene_effect_symbols(
    row: dict[str, Any],
    top_n: int,
    threshold: float,
) -> tuple[list[str], float]:
    scored: list[tuple[float, str]] = []
    ignored = {name.casefold() for name in [*MODEL_ID_COLUMNS, *CELL_NAME_COLUMNS]}
    for column, value in row.items():
        if str(column).casefold() in ignored:
            continue
        symbol = normalize_gene_symbol(column)
        if not symbol:
            continue
        score = parse_score(value)
        if score is None:
            continue
        if score <= threshold:
            scored.append((score, symbol))
    scored.sort(key=lambda item: (item[0], item[1]))
    selected = [symbol for _score, symbol in scored[: max(0, top_n)]]
    strongest = abs(scored[0][0]) if scored else 0.0
    return selected, strongest


def source_file_record_local(source: str, path: Path | None, source_url: str = "") -> dict[str, Any] | None:
    if path is None:
        return None
    return source_file_record(source, path, source_url or str(path))


def extract_xrefs(node: dict[str, Any], prefix: str) -> list[str]:
    values = []
    wanted = prefix.upper() + ":"
    for xref in node.get("external_xrefs") or []:
        text = str(xref)
        if text.upper().startswith(wanted):
            values.append(text.split(":", 1)[1])
    return sorted(set(values))


def load_graph_targets(graph_root: Path, release_id: str) -> list[dict[str, Any]]:
    if ds is None:
        raise RuntimeError("pyarrow is required to read graph target nodes.")
    path = graph_root / release_id / "nodes.parquet"
    dataset = ds.dataset(path, format="parquet")
    return dataset.to_table(filter=ds.field("node_type") == "target").to_pylist()


def select_targets(args: argparse.Namespace) -> list[dict[str, Any]]:
    workspace = Path(args.workspace).resolve()
    release_id = args.release_id
    if args.analysis_records:
        service = MetaboService(workspace, release_id=release_id)
        analyzed = service.analyze_metabolites(
            load_records_arg(args.analysis_records),
            max_paths=args.max_paths,
            max_hops=args.max_hops,
            context=load_context_arg(args.context_json),
            context_mode=args.context_mode,
        )
        rows = analyzed.get("rankings", {}).get("targets", [])[: args.target_limit]
        targets = []
        for row in rows:
            node = service.get_node(row.get("target_uid", ""))
            if node:
                targets.append(node)
        return targets

    graph_root = workspace / args.graph_root
    targets = load_graph_targets(graph_root, release_id)
    symbols = {normalize_prediction_text(item) for item in re.split(r"[,;\s]+", args.target_symbols or "") if item.strip()}
    if symbols:
        targets = [
            row
            for row in targets
            if normalize_prediction_text(row.get("canonical_name")) in symbols
            or normalize_prediction_text(row.get("display_name")) in symbols
        ]
    else:
        targets = [row for row in targets if extract_xrefs(row, "CHEMBL") or extract_xrefs(row, "UNIPROT_SWISSPROT")]
    targets.sort(key=lambda row: (row.get("canonical_name", ""), row.get("node_uid", "")))
    return targets[: args.target_limit]


def molecule_pref_name(molecule_id: str, cache: dict[str, str]) -> str:
    if molecule_id in cache:
        return cache[molecule_id]
    url = f"{CHEMBL_BASE}/molecule/{urllib.parse.quote(molecule_id)}.json"
    try:
        payload = request_json(url, timeout=40, retries=2)
        name = str(payload.get("pref_name") or payload.get("molecule_chembl_id") or molecule_id)
    except Exception:
        name = molecule_id
    cache[molecule_id] = name
    return name


def fetch_paginated_chembl(resource: str, cache_path: Path, limit: int = 1000) -> list[dict[str, Any]]:
    if cache_path.exists() and cache_path.stat().st_size > 0:
        return [json.loads(line) for line in cache_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    offset = 0
    while True:
        url = f"{CHEMBL_BASE}/{resource}.json?limit={int(limit)}&offset={int(offset)}"
        payload = request_json(url, timeout=60, retries=3)
        chunk = payload.get(resource + "s", [])
        rows.extend(chunk)
        page_meta = payload.get("page_meta") or {}
        next_url = page_meta.get("next")
        if not next_url:
            break
        offset += int(page_meta.get("limit") or limit)
    cache_path.write_text("\n".join(json.dumps(row, ensure_ascii=False, sort_keys=True) for row in rows) + "\n", encoding="utf-8")
    return rows


def fetch_chembl_molecule_names(molecule_ids: set[str], cache_path: Path, batch_size: int = 100) -> dict[str, str]:
    if cache_path.exists() and cache_path.stat().st_size > 0:
        return json.loads(cache_path.read_text(encoding="utf-8"))
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    result: dict[str, str] = {}
    values = sorted(molecule_ids)
    for idx in range(0, len(values), batch_size):
        batch = values[idx : idx + batch_size]
        query = urllib.parse.urlencode({"molecule_chembl_id__in": ",".join(batch), "limit": len(batch)})
        try:
            payload = request_json(f"{CHEMBL_BASE}/molecule.json?{query}", timeout=60, retries=2)
        except Exception:
            continue
        for molecule in payload.get("molecules", []):
            molecule_id = str(molecule.get("molecule_chembl_id") or "")
            name = str(molecule.get("pref_name") or molecule.get("molecule_chembl_id") or molecule_id)
            if molecule_id:
                result[molecule_id] = name
    for molecule_id in values:
        result.setdefault(molecule_id, molecule_id)
    cache_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def build_chembl_rows(
    targets: list[dict[str, Any]],
    include_activities: bool,
    activity_limit: int,
    raw_dir: Path | None = None,
) -> list[dict[str, Any]]:
    rows = []
    molecule_cache: dict[str, str] = {}
    seen: set[tuple[str, str, str]] = set()
    target_by_chembl: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for target in targets:
        for chembl_target_id in extract_xrefs(target, "CHEMBL"):
            target_by_chembl[chembl_target_id].append(target)
    if raw_dir is not None:
        raw_dir.mkdir(parents=True, exist_ok=True)
        mechanisms = fetch_paginated_chembl("mechanism", raw_dir / "chembl_mechanisms.jsonl")
        molecule_ids = {
            str(item.get("parent_molecule_chembl_id") or item.get("molecule_chembl_id") or "")
            for item in mechanisms
            if str(item.get("target_chembl_id") or "") in target_by_chembl
        }
        molecule_cache = fetch_chembl_molecule_names({value for value in molecule_ids if value}, raw_dir / "chembl_molecule_names.json")
        for mechanism in mechanisms:
            chembl_target_id = str(mechanism.get("target_chembl_id") or "")
            if chembl_target_id not in target_by_chembl:
                continue
            molecule_id = str(mechanism.get("parent_molecule_chembl_id") or mechanism.get("molecule_chembl_id") or "")
            if not molecule_id:
                continue
            for target in target_by_chembl[chembl_target_id]:
                target_symbol = str(target.get("canonical_name") or "")
                key = ("chembl_mechanism", str(target.get("node_uid")), molecule_id)
                if key in seen:
                    continue
                seen.add(key)
                max_phase = clamp_unit(float(mechanism.get("max_phase") or 0.0) / 4.0)
                confidence = clamp_unit(
                    0.55
                    + (0.1 if mechanism.get("direct_interaction") else 0.0)
                    + (0.1 if mechanism.get("molecular_mechanism") else 0.0)
                    + (0.15 * max_phase)
                )
                rows.append(
                    {
                        "drug_id": molecule_id,
                        "drug_name": molecule_cache.get(molecule_id, molecule_id),
                        "target_uid": target.get("node_uid", ""),
                        "target_symbol": target_symbol,
                        "target_name": target.get("display_name", ""),
                        "target_id": chembl_target_id,
                        "mechanism": mechanism.get("action_type") or mechanism.get("mechanism_of_action") or "unknown",
                        "confidence": round(confidence, 6),
                        "context_terms": merge_terms(target_symbol, target.get("display_name"), target.get("primary_external_id")),
                        "source_name": "ChEMBL",
                        "source_record_id": f"{chembl_target_id}|{molecule_id}|mechanism:{mechanism.get('mec_id', '')}",
                        "evidence_level": "curated_mechanism",
                        "license_id": "open_core:chembl",
                        "source_url": f"{CHEMBL_BASE}/mechanism.json",
                        "source_version": "chembl_webresource_current",
                    }
                )
    if rows and not include_activities:
        rows.sort(key=lambda row: (-float(row["confidence"]), row["target_symbol"], row["drug_name"], row["source_record_id"]))
        return rows

    for target in targets:
        target_symbol = str(target.get("canonical_name") or "")
        context_terms = merge_terms(target_symbol, target.get("display_name"), target.get("primary_external_id"))
        for chembl_target_id in extract_xrefs(target, "CHEMBL"):
            url = f"{CHEMBL_BASE}/mechanism.json?target_chembl_id={urllib.parse.quote(chembl_target_id)}&limit=1000"
            try:
                mechanisms = request_json(url).get("mechanisms", [])
            except Exception:
                mechanisms = []
            for mechanism in mechanisms:
                molecule_id = str(mechanism.get("parent_molecule_chembl_id") or mechanism.get("molecule_chembl_id") or "")
                if not molecule_id:
                    continue
                key = ("chembl_mechanism", str(target.get("node_uid")), molecule_id)
                if key in seen:
                    continue
                seen.add(key)
                max_phase = clamp_unit(float(mechanism.get("max_phase") or 0.0) / 4.0)
                confidence = clamp_unit(
                    0.55
                    + (0.1 if mechanism.get("direct_interaction") else 0.0)
                    + (0.1 if mechanism.get("molecular_mechanism") else 0.0)
                    + (0.15 * max_phase)
                )
                rows.append(
                    {
                        "drug_id": molecule_id,
                        "drug_name": molecule_pref_name(molecule_id, molecule_cache),
                        "target_uid": target.get("node_uid", ""),
                        "target_symbol": target_symbol,
                        "target_name": target.get("display_name", ""),
                        "target_id": chembl_target_id,
                        "mechanism": mechanism.get("action_type") or mechanism.get("mechanism_of_action") or "unknown",
                        "confidence": round(confidence, 6),
                        "context_terms": context_terms,
                        "source_name": "ChEMBL",
                        "source_record_id": f"{chembl_target_id}|{molecule_id}|mechanism:{mechanism.get('mec_id', '')}",
                        "evidence_level": "curated_mechanism",
                        "license_id": "open_core:chembl",
                        "source_url": url,
                        "source_version": "chembl_webresource_current",
                    }
                )
            if not include_activities:
                continue
            activity_url = (
                f"{CHEMBL_BASE}/activity.json?target_chembl_id={urllib.parse.quote(chembl_target_id)}"
                f"&pchembl_value__gte=7&standard_flag=1&limit={int(activity_limit)}"
            )
            try:
                activities = request_json(activity_url).get("activities", [])
            except Exception:
                activities = []
            for activity in activities[:activity_limit]:
                molecule_id = str(activity.get("parent_molecule_chembl_id") or activity.get("molecule_chembl_id") or "")
                if not molecule_id:
                    continue
                key = ("chembl_activity", str(target.get("node_uid")), molecule_id)
                if key in seen:
                    continue
                seen.add(key)
                pchembl = clamp_unit(float(activity.get("pchembl_value") or 0.0) / 10.0)
                rows.append(
                    {
                        "drug_id": molecule_id,
                        "drug_name": activity.get("molecule_pref_name") or molecule_pref_name(molecule_id, molecule_cache),
                        "target_uid": target.get("node_uid", ""),
                        "target_symbol": target_symbol,
                        "target_name": target.get("display_name", ""),
                        "target_id": chembl_target_id,
                        "mechanism": activity.get("standard_type") or "bioactivity",
                        "confidence": round(0.65 * pchembl, 6),
                        "context_terms": context_terms,
                        "source_name": "ChEMBL",
                        "source_record_id": f"{chembl_target_id}|{molecule_id}|activity:{activity.get('activity_id', '')}",
                        "evidence_level": "bioactivity_assay",
                        "license_id": "open_core:chembl",
                        "source_url": activity_url,
                        "source_version": "chembl_webresource_current",
                    }
                )
    rows.sort(key=lambda row: (-float(row["confidence"]), row["target_symbol"], row["drug_name"], row["source_record_id"]))
    return rows


def build_drugcentral_rows(targets: list[dict[str, Any]], raw_dir: Path) -> list[dict[str, Any]]:
    path = download_file(DRUGCENTRAL_TARGET_URL, raw_dir / "drug.target.interaction.tsv.gz", timeout=240)
    by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    by_uniprot: dict[str, list[dict[str, Any]]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if str(row.get("ORGANISM", "")).casefold() != "homo sapiens":
                continue
            if row.get("GENE"):
                by_gene[normalize_prediction_text(row["GENE"])].append(row)
            if row.get("ACCESSION"):
                by_uniprot[normalize_prediction_text(row["ACCESSION"])].append(row)
    rows = []
    seen: set[tuple[str, str]] = set()
    for target in targets:
        target_symbol = str(target.get("canonical_name") or "")
        uniprots = [xref.split(":", 1)[1] for xref in target.get("external_xrefs", []) if str(xref).startswith("UNIPROT_SWISSPROT:")]
        candidates = list(by_gene.get(normalize_prediction_text(target_symbol), []))
        for accession in uniprots:
            candidates.extend(by_uniprot.get(normalize_prediction_text(accession), []))
        for row in candidates:
            drug_id = str(row.get("STRUCT_ID") or row.get("DRUG_NAME") or "")
            if not drug_id:
                continue
            key = (str(target.get("node_uid")), drug_id)
            if key in seen:
                continue
            seen.add(key)
            act_value = clamp_unit(float(row.get("ACT_VALUE") or 0.0) / 10.0)
            confidence = 0.35 + (0.25 if row.get("MOA") == "1" else 0.0) + (0.15 if row.get("ACTION_TYPE") else 0.0) + (0.2 * act_value)
            rows.append(
                {
                    "drug_id": f"DRUGCENTRAL:{drug_id}",
                    "drug_name": row.get("DRUG_NAME", ""),
                    "target_uid": target.get("node_uid", ""),
                    "target_symbol": target_symbol,
                    "target_name": target.get("display_name", ""),
                    "target_id": row.get("ACCESSION") or target.get("primary_external_id", ""),
                    "mechanism": row.get("ACTION_TYPE") or row.get("ACT_TYPE") or "unknown",
                    "confidence": round(clamp_unit(confidence), 6),
                    "context_terms": merge_terms(target_symbol, target.get("display_name"), row.get("TARGET_CLASS")),
                    "source_name": "DrugCentral",
                    "source_record_id": f"{drug_id}|{row.get('GENE', '')}|{row.get('ACCESSION', '')}",
                    "evidence_level": "drug_target_interaction",
                    "license_id": "open_core:drugcentral",
                    "source_url": DRUGCENTRAL_TARGET_URL,
                    "source_version": "DrugCentral 2021_09_01 download",
                }
            )
    rows.sort(key=lambda row: (-float(row["confidence"]), row["target_symbol"], row["drug_name"]))
    return rows


def phase_confidence(value: Any) -> float:
    text = str(value or "").upper()
    if "APPROVAL" in text:
        return 0.95
    if "PHASE_4" in text:
        return 0.92
    if "PHASE_3" in text:
        return 0.86
    if "PHASE_2" in text:
        return 0.74
    if "PHASE_1" in text:
        return 0.62
    return 0.50


def build_opentargets_rows(targets: list[dict[str, Any]], raw_dir: Path) -> list[dict[str, Any]]:
    if ds is None:
        raise RuntimeError("pyarrow is required to read Open Targets clinical target parquet.")
    clinical_url = f"{OPENTARGETS_2603_BASE}/clinical_target/"
    drug_url = f"{OPENTARGETS_2603_BASE}/drug_molecule/"
    clinical_dir = raw_dir / "clinical_target"
    drug_dir = raw_dir / "drug_molecule"
    clinical_files = download_index_files(clinical_url, clinical_dir, r"\.parquet$")
    drug_files = download_index_files(drug_url, drug_dir, r"\.parquet$")
    if not clinical_files:
        return []

    target_by_ensembl: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for target in targets:
        for ensembl_id in extract_xrefs(target, "ENSEMBL"):
            target_by_ensembl[ensembl_id].append(target)
    if not target_by_ensembl:
        return []

    drug_names: dict[str, dict[str, str]] = {}
    if drug_files:
        drug_table = ds.dataset(drug_dir, format="parquet").to_table(columns=["id", "name", "drugType", "maximumClinicalStage"])
        for row in drug_table.to_pylist():
            drug_id = str(row.get("id") or "")
            if not drug_id:
                continue
            drug_names[drug_id] = {
                "name": str(row.get("name") or drug_id),
                "drugType": str(row.get("drugType") or ""),
                "maximumClinicalStage": str(row.get("maximumClinicalStage") or ""),
            }

    rows = []
    seen: set[tuple[str, str, str]] = set()
    table = ds.dataset(clinical_dir, format="parquet").to_table(
        columns=["id", "drugId", "targetId", "diseases", "clinicalReportIds", "maxClinicalStage"]
    )
    for item in table.to_pylist():
        target_id = str(item.get("targetId") or "")
        if target_id not in target_by_ensembl:
            continue
        drug_id = str(item.get("drugId") or "")
        if not drug_id:
            continue
        drug_meta = drug_names.get(drug_id, {})
        disease_terms = []
        for disease in item.get("diseases") or []:
            if isinstance(disease, dict):
                disease_terms.append(disease.get("diseaseFromSource") or disease.get("diseaseId") or "")
        report_count = len(item.get("clinicalReportIds") or [])
        max_stage = str(item.get("maxClinicalStage") or drug_meta.get("maximumClinicalStage") or "")
        confidence = clamp_unit(phase_confidence(max_stage) + min(0.05, 0.01 * report_count))
        for target in target_by_ensembl[target_id]:
            key = (drug_id, str(target.get("node_uid")), str(item.get("id") or ""))
            if key in seen:
                continue
            seen.add(key)
            target_symbol = str(target.get("canonical_name") or "")
            rows.append(
                {
                    "drug_id": drug_id,
                    "drug_name": drug_meta.get("name") or drug_id,
                    "target_uid": target.get("node_uid", ""),
                    "target_symbol": target_symbol,
                    "target_name": target.get("display_name", ""),
                    "target_id": target_id,
                    "mechanism": f"clinical_target:{max_stage or 'UNKNOWN'}",
                    "confidence": round(confidence, 6),
                    "context_terms": merge_terms(target_symbol, target.get("display_name"), drug_meta.get("drugType"), disease_terms[:20]),
                    "source_name": "Open Targets",
                    "source_record_id": str(item.get("id") or f"{drug_id}|{target_id}"),
                    "evidence_level": "clinical_target",
                    "license_id": "open_core:opentargets_platform",
                    "source_url": clinical_url,
                    "source_version": "Open Targets Platform 26.03 clinical_target",
                }
            )
    rows.sort(key=lambda row: (-float(row["confidence"]), row["target_symbol"], row["drug_name"], row["source_record_id"]))
    return rows


def dgidb_query(names: list[str]) -> dict[str, Any]:
    query = """
    query($names:[String!]){
      genes(names:$names, first:1000){
        nodes {
          name
          conceptId
          interactions {
            id
            interactionScore
            evidenceScore
            drug { name conceptId approved antiNeoplastic }
            interactionTypes { type directionality }
            sources { sourceDbName sourceDbVersion license licenseLink sourceTrustLevel { level } }
            publications { pmid }
          }
        }
      }
    }
    """
    return request_json_post(DGIDB_GRAPHQL_URL, {"query": query, "variables": {"names": names}}, timeout=60, retries=3)


def fetch_dgidb_interactions(cache_path: Path, page_size: int = 100) -> list[dict[str, Any]]:
    if cache_path.exists() and cache_path.stat().st_size > 0:
        return [json.loads(line) for line in cache_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    query = """
    query($first:Int,$after:String){
      interactions(first:$first, after:$after){
        pageInfo { hasNextPage endCursor }
        nodes {
          id
          interactionScore
          evidenceScore
          gene { name conceptId }
          drug { name conceptId approved antiNeoplastic }
          interactionTypes { type directionality }
          sources { sourceDbName sourceDbVersion license licenseLink sourceTrustLevel { level } }
          publications { pmid }
        }
      }
    }
    """
    rows: list[dict[str, Any]] = []
    after: str | None = None
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w", encoding="utf-8") as handle:
        page = 0
        while True:
            payload = request_json_post(
                DGIDB_GRAPHQL_URL,
                {"query": query, "variables": {"first": int(page_size), "after": after}},
                timeout=60,
                retries=3,
            )
            interactions = ((payload.get("data") or {}).get("interactions") or {})
            chunk = interactions.get("nodes") or []
            rows.extend(chunk)
            for row in chunk:
                handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
            handle.flush()
            page += 1
            if page % 25 == 0:
                print(f"DGIdb cached {len(rows)} interaction rows...", flush=True)
            page_info = interactions.get("pageInfo") or {}
            if not page_info.get("hasNextPage"):
                break
            after = page_info.get("endCursor")
            if not after:
                break
    return rows


def build_dgidb_rows(targets: list[dict[str, Any]], batch_size: int, raw_dir: Path | None = None) -> list[dict[str, Any]]:
    target_by_symbol: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for target in targets:
        symbol = normalize_prediction_text(target.get("canonical_name"))
        if symbol:
            target_by_symbol[symbol].append(target)
    rows = []
    seen: set[tuple[str, str, str]] = set()
    if raw_dir is not None:
        interactions = fetch_dgidb_interactions(raw_dir / "dgidb_interactions.jsonl", page_size=max(1, batch_size))
        gene_items = []
        for interaction in interactions:
            gene = interaction.get("gene") or {}
            if normalize_prediction_text(gene.get("name")) in target_by_symbol:
                gene_items.append({"name": gene.get("name"), "conceptId": gene.get("conceptId"), "interactions": [interaction]})
    else:
        symbols = sorted(target_by_symbol)
        gene_items = []
        for idx in range(0, len(symbols), max(1, batch_size)):
            batch_symbols = symbols[idx : idx + max(1, batch_size)]
            try:
                payload = dgidb_query(batch_symbols)
            except Exception:
                continue
            gene_items.extend((((payload.get("data") or {}).get("genes") or {}).get("nodes") or []))

    for gene in gene_items:
        symbol = normalize_prediction_text(gene.get("name"))
        for target in target_by_symbol.get(symbol, []):
            target_symbol = str(target.get("canonical_name") or "")
            for interaction in gene.get("interactions") or []:
                drug = interaction.get("drug") or {}
                drug_id = str(drug.get("conceptId") or drug.get("name") or "")
                if not drug_id:
                    continue
                interaction_id = str(interaction.get("id") or "")
                key = (drug_id, str(target.get("node_uid")), interaction_id)
                if key in seen:
                    continue
                seen.add(key)
                interaction_types = interaction.get("interactionTypes") or []
                sources = interaction.get("sources") or []
                publications = interaction.get("publications") or []
                source_names = sorted({str(src.get("sourceDbName") or "") for src in sources if src.get("sourceDbName")})
                source_versions = sorted({str(src.get("sourceDbVersion") or "") for src in sources if src.get("sourceDbVersion")})
                mechanism = ";".join(
                    sorted(
                        {
                            str(item.get("type") or item.get("directionality") or "")
                            for item in interaction_types
                            if item.get("type") or item.get("directionality")
                        }
                    )
                ) or "drug_gene_interaction"
                evidence_score = clamp_unit(float(interaction.get("evidenceScore") or 0.0) / 10.0)
                interaction_score = clamp_unit(interaction.get("interactionScore"))
                confidence = clamp_unit(
                    0.42
                    + 0.25 * interaction_score
                    + 0.15 * evidence_score
                    + (0.08 if drug.get("approved") else 0.0)
                    + (0.03 if drug.get("antiNeoplastic") else 0.0)
                )
                rows.append(
                    {
                        "drug_id": f"DGIDB:{drug_id}",
                        "drug_name": drug.get("name") or drug_id,
                        "target_uid": target.get("node_uid", ""),
                        "target_symbol": target_symbol,
                        "target_name": target.get("display_name", ""),
                        "target_id": gene.get("conceptId") or target.get("primary_external_id", ""),
                        "mechanism": mechanism,
                        "confidence": round(confidence, 6),
                        "context_terms": merge_terms(
                            target_symbol,
                            target.get("display_name"),
                            mechanism,
                            source_names[:12],
                            [pub.get("pmid") for pub in publications[:12]],
                        ),
                        "source_name": "DGIdb",
                        "source_record_id": interaction_id or f"{drug_id}|{target_symbol}",
                        "evidence_level": "aggregated_drug_gene_interaction",
                        "license_id": "aggregated:dgidb_source_specific",
                        "source_url": DGIDB_GRAPHQL_URL,
                        "source_version": ";".join(source_versions[:8]) or "DGIdb GraphQL current",
                    }
                )
    rows.sort(key=lambda row: (-float(row["confidence"]), row["target_symbol"], row["drug_name"], row["source_record_id"]))
    return rows


def iter_cellosaurus_entries(path: Path):
    entry: dict[str, list[str]] = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("//"):
                if entry:
                    yield entry
                entry = defaultdict(list)
                continue
            if len(line) < 5:
                continue
            key = line[:2]
            value = line[5:].strip()
            if value:
                entry[key].append(value)


def cellosaurus_entry_terms(entry: dict[str, list[str]]) -> list[str]:
    values: list[str] = []
    values.extend(entry.get("ID", []))
    values.extend(entry.get("SY", []))
    values.extend(entry.get("DI", []))
    values.extend(entry.get("DR", []))
    values.extend(entry.get("CC", []))
    return prediction_terms_from_value(values)


def build_cellosaurus_rows(raw_dir: Path, context_terms: list[str], max_rows: int) -> list[dict[str, Any]]:
    path = download_file(CELLOSAURUS_TXT_URL, raw_dir / "cellosaurus.txt", timeout=240)
    wanted = set(context_terms)
    rows = []
    for entry in iter_cellosaurus_entries(path):
        ox = " ".join(entry.get("OX", [])).casefold()
        if "ncbi_taxid=9606" not in ox and "homo sapiens" not in ox:
            continue
        terms = cellosaurus_entry_terms(entry)
        if wanted and not (wanted & set(terms)):
            continue
        accession = (entry.get("AC") or [""])[0].rstrip(";")
        name = (entry.get("ID") or [accession])[0]
        if not accession or not name:
            continue
        aliases = ";".join(entry.get("SY", [])[:5])
        rows.append(
            {
                "cell_type_id": f"Cellosaurus:{accession}",
                "cell_type_name": name,
                "cell_line": aliases,
                "target_symbols": "",
                "pathway_terms": "",
                "metabolite_terms": "",
                "context_terms": ";".join(terms[:50]),
                "confidence": 0.35,
                "source_name": "Cellosaurus",
                "source_record_id": accession,
                "evidence_level": "cell_line_context",
                "license_id": "open_core:cellosaurus",
                "source_url": CELLOSAURUS_TXT_URL,
                "source_version": "cellosaurus_current_ftp",
            }
        )
        if max_rows and len(rows) >= max_rows:
            break
    return rows


def _first_cellosaurus_value(values: Any, kind: str | None = None) -> str:
    if not isinstance(values, list):
        return ""
    for item in values:
        if not isinstance(item, dict):
            continue
        if kind is None or str(item.get("type", "")).casefold() == kind.casefold():
            value = str(item.get("value") or "").strip()
            if value:
                return value
    return ""


def build_cellosaurus_api_rows(path: Path, source_url: str, context_terms: list[str], max_rows: int) -> list[dict[str, Any]]:
    if not path.exists() or not path.is_file():
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    records = ((payload.get("Cellosaurus") or {}).get("cell-line-list") or [])
    rows = []
    for record in records:
        accession = _first_cellosaurus_value(record.get("accession-list"), "primary")
        name = _first_cellosaurus_value(record.get("name-list"), "identifier") or accession
        aliases = [
            str(item.get("value")).strip()
            for item in record.get("name-list", [])
            if isinstance(item, dict) and str(item.get("type", "")).casefold() == "synonym" and item.get("value")
        ][:5]
        terms = merge_terms(name, aliases, context_terms, "Homo sapiens")
        if not accession or not name:
            continue
        rows.append(
            {
                "cell_type_id": f"Cellosaurus:{accession}",
                "cell_type_name": name,
                "cell_line": ";".join(aliases),
                "target_symbols": "",
                "pathway_terms": "",
                "metabolite_terms": "",
                "context_terms": terms,
                "confidence": 0.45,
                "source_name": "Cellosaurus",
                "source_record_id": accession,
                "evidence_level": "cell_line_context_api",
                "license_id": "open_core:cellosaurus",
                "source_url": source_url,
                "source_version": "cellosaurus_api_current",
            }
        )
        if max_rows and len(rows) >= max_rows:
            break
    return rows


def load_depmap_model_context(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None or not path.exists() or not path.is_file():
        return {}
    contexts: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            model_id = first_value(row, MODEL_ID_COLUMNS)
            name = first_value(row, CELL_NAME_COLUMNS, model_id)
            state = first_value(row, ("cell_state", "CellState", "lineage_subtype", "LineageSubtype"))
            context = merge_terms(
                row.get("OncotreeLineage"),
                row.get("OncotreePrimaryDisease"),
                row.get("OncotreeSubtype"),
                row.get("Tissue"),
                row.get("lineage"),
                row.get("primary_disease"),
                state,
                name,
            )
            if model_id and name:
                contexts[model_id] = {"name": name, "state": state, "context": context}
    return contexts


def build_depmap_rows(
    model_path: Path | None,
    gene_effect_path: Path | None = None,
    target_uid_by_symbol: dict[str, str] | None = None,
    top_genes_per_model: int = 50,
    dependency_threshold: float = -0.5,
) -> list[dict[str, Any]]:
    contexts = load_depmap_model_context(model_path)
    target_uid_by_symbol = target_uid_by_symbol or {}
    rows: list[dict[str, Any]] = []
    row_by_model: dict[str, dict[str, Any]] = {}
    if gene_effect_path is not None and gene_effect_path.exists() and gene_effect_path.is_file():
        with gene_effect_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                model_id = first_value(row, MODEL_ID_COLUMNS)
                if not model_id:
                    continue
                context = contexts.get(model_id, {})
                name = context.get("name") or first_value(row, CELL_NAME_COLUMNS, model_id)
                symbols, strongest = select_depmap_gene_effect_symbols(row, top_genes_per_model, dependency_threshold)
                if not symbols:
                    continue
                confidence = clamp_unit(0.55 + min(strongest, 1.5) * 0.18 + min(len(symbols), 50) * 0.002)
                item = {
                    "cell_type_id": f"DepMap:{model_id}",
                    "cell_type_name": name,
                    "cell_line": name,
                    "target_uids": target_uids_for_symbols(symbols, target_uid_by_symbol),
                    "target_symbols": ";".join(symbols),
                    "pathway_terms": "",
                    "metabolite_terms": "",
                    "context_terms": context.get("context", merge_terms(name)),
                    "cell_state": context.get("state", ""),
                    "confidence": round(confidence, 6),
                    "source_name": "DepMap",
                    "source_record_id": model_id,
                    "evidence_level": "gene_dependency_signature",
                    "license_id": "public_research:depmap",
                    "source_url": str(gene_effect_path),
                    "source_version": "local_depmap_gene_effect_csv",
                }
                rows.append(item)
                row_by_model[model_id] = item

    for model_id, context in contexts.items():
        if model_id in row_by_model:
            continue
        rows.append(
            {
                "cell_type_id": f"DepMap:{model_id}",
                "cell_type_name": context.get("name", model_id),
                "cell_line": context.get("name", model_id),
                "target_symbols": "",
                "pathway_terms": "",
                "metabolite_terms": "",
                "context_terms": context.get("context", ""),
                "cell_state": context.get("state", ""),
                "confidence": 0.45,
                "source_name": "DepMap",
                "source_record_id": model_id,
                "evidence_level": "cell_line_context",
                "license_id": "public_research:depmap",
                "source_url": str(model_path),
                "source_version": "local_depmap_model_csv",
            }
        )
    return rows


def marker_row_confidence(rows: list[dict[str, Any]], default: float) -> float:
    explicit = [clamp_unit(row.get("confidence"), default=-1.0) for row in rows if str(row.get("confidence", "")).strip()]
    explicit = [value for value in explicit if value >= 0.0]
    if explicit:
        return round(max(explicit), 6)
    scores = []
    for row in rows:
        for key in ("avg_log2FC", "avg_log2fc", "logfoldchanges", "score"):
            score = parse_score(row.get(key))
            if score is not None:
                scores.append(abs(score))
                break
    if scores:
        return round(clamp_unit(default + min(max(scores), 4.0) * 0.05), 6)
    return default


def build_marker_signature_rows(
    path: Path | None,
    source_name: str,
    target_uid_by_symbol: dict[str, str] | None = None,
    default_confidence: float = 0.55,
    source_version: str = "local_marker_signature_csv",
    license_id: str = "public_research:marker_signature",
) -> list[dict[str, Any]]:
    if path is None or not path.exists() or not path.is_file():
        return []
    target_uid_by_symbol = target_uid_by_symbol or {}
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            cell_type = first_value(row, MARKER_CELL_NAME_COLUMNS)
            if not cell_type:
                continue
            cell_type_id = first_value(row, MARKER_CELL_ID_COLUMNS, f"{source_name}:{cell_type}")
            groups[(cell_type_id, cell_type)].append(row)

    rows: list[dict[str, Any]] = []
    for (cell_type_id, cell_type), group in groups.items():
        symbols: list[str] = []
        pathway_terms: list[str] = []
        metabolite_terms: list[str] = []
        contexts: list[str] = []
        states: list[str] = []
        source_ids: list[str] = []
        for row in group:
            for column in MARKER_GENE_COLUMNS:
                symbols.extend(split_gene_symbols(row.get(column)))
            for column in MARKER_PATHWAY_COLUMNS:
                pathway_terms.extend(prediction_terms_from_value(row.get(column)))
            for column in MARKER_METABOLITE_COLUMNS:
                metabolite_terms.extend(prediction_terms_from_value(row.get(column)))
            contexts.append(merge_terms(*[row.get(column) for column in MARKER_CONTEXT_COLUMNS], cell_type))
            state = first_value(row, ("cell_state", "state", "phenotype", "annotation_state"))
            if state:
                states.append(state)
            source_id = first_value(row, ("source_record_id", "dataset_id", "study_id", "sample_id"))
            if source_id:
                source_ids.append(source_id)

        symbols = sorted(set(symbols))
        pathway_terms = sorted({term for term in pathway_terms if term})
        metabolite_terms = sorted({term for term in metabolite_terms if term})
        context_terms = merge_terms(contexts)
        if not any([symbols, pathway_terms, metabolite_terms, context_terms]):
            continue
        rows.append(
            {
                "cell_type_id": cell_type_id if ":" in cell_type_id else f"{source_name}:{cell_type_id}",
                "cell_type_name": cell_type,
                "cell_line": "",
                "target_uids": target_uids_for_symbols(symbols, target_uid_by_symbol),
                "target_symbols": ";".join(symbols),
                "pathway_terms": ";".join(pathway_terms),
                "metabolite_terms": ";".join(metabolite_terms),
                "context_terms": context_terms,
                "cell_state": ";".join(sorted(set(states))[:5]),
                "confidence": marker_row_confidence(group, default_confidence),
                "source_name": source_name,
                "source_record_id": ";".join(sorted(set(source_ids))[:5]) or cell_type_id,
                "evidence_level": "cell_type_signature",
                "license_id": first_value(group[0], ("license_id", "license"), license_id),
                "source_url": str(path),
                "source_version": first_value(group[0], ("source_version", "version"), source_version),
            }
        )
    return rows


def build_cellxgene_rows(path: Path | None, target_uid_by_symbol: dict[str, str] | None = None) -> list[dict[str, Any]]:
    return build_marker_signature_rows(
        path,
        source_name="CELLxGENE",
        target_uid_by_symbol=target_uid_by_symbol,
        default_confidence=0.55,
        source_version="local_cellxgene_signature_csv",
        license_id="public_research:cellxgene",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build drug/cell prediction overlay CSVs for a frozen release.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--graph-root", default=DEFAULT_GRAPH_ROOT)
    parser.add_argument("--output-root", default=DEFAULT_PREDICTION_OVERLAY_ROOT)
    parser.add_argument("--raw-root", default="raw_lake")
    parser.add_argument("--analysis-records", default="")
    parser.add_argument("--context-json", default="")
    parser.add_argument("--context-mode", choices=["soft", "hard"], default="soft")
    parser.add_argument("--target-symbols", default="")
    parser.add_argument("--target-limit", type=int, default=25)
    parser.add_argument("--max-paths", type=int, default=25)
    parser.add_argument("--max-hops", type=int, default=4)
    parser.add_argument("--include-chembl-activities", action="store_true")
    parser.add_argument("--chembl-activity-limit", type=int, default=25)
    parser.add_argument("--skip-chembl", action="store_true")
    parser.add_argument("--skip-drugcentral", action="store_true")
    parser.add_argument("--skip-dgidb", action="store_true")
    parser.add_argument("--skip-opentargets", action="store_true")
    parser.add_argument("--dgidb-batch-size", type=int, default=100)
    parser.add_argument("--skip-cellosaurus", action="store_true")
    parser.add_argument("--cellosaurus-max-rows", type=int, default=5000)
    parser.add_argument("--cellosaurus-api-json", default="")
    parser.add_argument("--cellosaurus-source-url", default="")
    parser.add_argument("--depmap-model-csv", default="")
    parser.add_argument("--depmap-gene-effect-csv", default="")
    parser.add_argument("--depmap-top-genes-per-model", type=int, default=50)
    parser.add_argument("--depmap-dependency-threshold", type=float, default=-0.5)
    parser.add_argument("--cellxgene-signature-csv", default="")
    parser.add_argument("--single-cell-marker-csv", default="")
    parser.add_argument(
        "--manifest-target-sample-limit",
        type=int,
        default=0,
        help="Number of selected targets to include in the manifest sample. Use 0 to record only counts.",
    )
    args = parser.parse_args(argv)

    workspace = Path(args.workspace).resolve()
    output_dir = workspace / args.output_root / args.release_id
    raw_dir = workspace / args.raw_root
    targets = select_targets(args)
    target_uid_by_symbol = target_uid_lookup_from_targets(targets)
    context = load_context_arg(args.context_json)
    context_terms = prediction_terms_from_value(context)

    drug_rows = []
    source_status = []
    if not args.skip_chembl:
        try:
            chembl_rows = build_chembl_rows(
                targets,
                include_activities=args.include_chembl_activities,
                activity_limit=args.chembl_activity_limit,
                raw_dir=raw_dir / "chembl" / args.release_id,
            )
            drug_rows.extend(chembl_rows)
            source_status.append(
                {
                    "source": "ChEMBL",
                    "rows": len(chembl_rows),
                    "status": "ok" if chembl_rows else "no_rows_or_unavailable",
                    "message": "" if chembl_rows else "No ChEMBL rows were fetched for the selected targets; this can happen when the endpoint is unavailable or targets have no matching mechanism/activity rows.",
                }
            )
        except Exception as exc:
            source_status.append({"source": "ChEMBL", "rows": 0, "status": "error", "message": str(exc)})
    else:
        source_status.append({"source": "ChEMBL", "rows": 0, "status": "skipped", "message": "Skipped for this run."})
    if not args.skip_drugcentral:
        try:
            dc_rows = build_drugcentral_rows(targets, raw_dir / "drugcentral" / args.release_id)
            drug_rows.extend(dc_rows)
            source_status.append({"source": "DrugCentral", "rows": len(dc_rows), "status": "ok"})
        except Exception as exc:
            source_status.append({"source": "DrugCentral", "rows": 0, "status": "error", "message": str(exc)})
    else:
        source_status.append({"source": "DrugCentral", "rows": 0, "status": "skipped", "message": "Skipped for this run."})
    if not args.skip_dgidb:
        try:
            dgidb_rows = build_dgidb_rows(
                targets,
                batch_size=args.dgidb_batch_size,
                raw_dir=raw_dir / "dgidb" / args.release_id,
            )
            drug_rows.extend(dgidb_rows)
            source_status.append({"source": "DGIdb", "rows": len(dgidb_rows), "status": "ok" if dgidb_rows else "no_rows_or_unavailable"})
        except Exception as exc:
            source_status.append({"source": "DGIdb", "rows": 0, "status": "error", "message": str(exc)})
    else:
        source_status.append({"source": "DGIdb", "rows": 0, "status": "skipped", "message": "Skipped for this run."})
    if not args.skip_opentargets:
        try:
            ot_rows = build_opentargets_rows(targets, raw_dir / "opentargets_drug_evidence" / args.release_id)
            drug_rows.extend(ot_rows)
            source_status.append({"source": "Open Targets", "rows": len(ot_rows), "status": "ok" if ot_rows else "no_rows_or_unavailable"})
        except Exception as exc:
            source_status.append({"source": "Open Targets", "rows": 0, "status": "error", "message": str(exc)})
    else:
        source_status.append({"source": "Open Targets", "rows": 0, "status": "skipped", "message": "Skipped for this run."})

    cell_rows = []
    if not args.skip_cellosaurus:
        try:
            if args.cellosaurus_api_json:
                cs_rows = build_cellosaurus_api_rows(
                    Path(args.cellosaurus_api_json),
                    args.cellosaurus_source_url or "local_cellosaurus_api_json",
                    context_terms,
                    args.cellosaurus_max_rows,
                )
            else:
                cs_rows = build_cellosaurus_rows(raw_dir / "cellosaurus" / args.release_id, context_terms, args.cellosaurus_max_rows)
            cell_rows.extend(cs_rows)
            source_status.append({"source": "Cellosaurus", "rows": len(cs_rows), "status": "ok"})
        except Exception as exc:
            source_status.append({"source": "Cellosaurus", "rows": 0, "status": "error", "message": str(exc)})
    depmap_path = Path(args.depmap_model_csv) if args.depmap_model_csv else None
    depmap_gene_effect_path = Path(args.depmap_gene_effect_csv) if args.depmap_gene_effect_csv else None
    depmap_rows = build_depmap_rows(
        depmap_path,
        gene_effect_path=depmap_gene_effect_path,
        target_uid_by_symbol=target_uid_by_symbol,
        top_genes_per_model=args.depmap_top_genes_per_model,
        dependency_threshold=args.depmap_dependency_threshold,
    )
    cell_rows.extend(depmap_rows)
    source_status.append({"source": "DepMap", "rows": len(depmap_rows), "status": "ok" if depmap_rows else "not_provided"})
    cellxgene_path = Path(args.cellxgene_signature_csv) if args.cellxgene_signature_csv else None
    cxg_rows = build_cellxgene_rows(cellxgene_path, target_uid_by_symbol=target_uid_by_symbol)
    cell_rows.extend(cxg_rows)
    source_status.append({"source": "CELLxGENE", "rows": len(cxg_rows), "status": "ok" if cxg_rows else "not_provided"})
    marker_path = Path(args.single_cell_marker_csv) if args.single_cell_marker_csv else None
    marker_rows = build_marker_signature_rows(
        marker_path,
        source_name="single_cell_marker",
        target_uid_by_symbol=target_uid_by_symbol,
        default_confidence=0.5,
        source_version="local_single_cell_marker_csv",
        license_id="public_research:single_cell_marker",
    )
    cell_rows.extend(marker_rows)
    source_status.append({"source": "single_cell_marker", "rows": len(marker_rows), "status": "ok" if marker_rows else "not_provided"})

    # Keep highest-confidence duplicate drug-target rows first.
    unique_drugs: dict[tuple[str, str, str], dict[str, Any]] = {}
    for row in sorted(drug_rows, key=lambda item: -float(item.get("confidence") or 0.0)):
        key = (row.get("drug_id", ""), row.get("target_uid", ""), row.get("mechanism", ""))
        unique_drugs.setdefault(key, row)
    drug_rows = sorted(unique_drugs.values(), key=lambda row: (-float(row.get("confidence") or 0.0), row.get("target_symbol", ""), row.get("drug_name", "")))

    write_csv(output_dir / "drug_targets.csv", drug_rows, FIELDNAMES_DRUG)
    write_csv(output_dir / "cell_type_signatures.csv", cell_rows, FIELDNAMES_CELL)
    successful_sources = {row.get("source") for row in source_status if row.get("status") == "ok"}
    source_files = []
    if "DrugCentral" in successful_sources:
        source_files.append(
            source_file_record(
                "DrugCentral",
                raw_dir / "drugcentral" / args.release_id / "drug.target.interaction.tsv.gz",
                DRUGCENTRAL_TARGET_URL,
            )
        )
    if "ChEMBL" in successful_sources:
        source_files.append(
            source_file_record(
                "ChEMBL",
                raw_dir / "chembl" / args.release_id / "chembl_mechanisms.jsonl",
                f"{CHEMBL_BASE}/mechanism.json",
            )
        )
    if "DGIdb" in successful_sources:
        source_files.append(
            source_file_record(
                "DGIdb",
                raw_dir / "dgidb" / args.release_id / "dgidb_interactions.jsonl",
                DGIDB_GRAPHQL_URL,
            )
        )
    if "Open Targets" in successful_sources:
        for source_name, source_path, source_url in (
            (
                "Open Targets clinical_target",
                raw_dir / "opentargets_drug_evidence" / args.release_id / "clinical_target" / "clinical_target.parquet",
                f"{OPENTARGETS_2603_BASE}/clinical_target/clinical_target.parquet",
            ),
        ):
            source_files.append(source_file_record(source_name, source_path, source_url))
    if "Cellosaurus" in successful_sources:
        cellosaurus_source_path = (
            Path(args.cellosaurus_api_json)
            if args.cellosaurus_api_json
            else raw_dir / "cellosaurus" / args.release_id / "cellosaurus.txt"
        )
        source_files.append(
            source_file_record(
                "Cellosaurus",
                cellosaurus_source_path,
                args.cellosaurus_source_url if args.cellosaurus_api_json else CELLOSAURUS_TXT_URL,
            )
        )
    if depmap_path is not None:
        source_files.append(source_file_record_local("DepMap model metadata", depmap_path))
    if depmap_gene_effect_path is not None:
        source_files.append(source_file_record_local("DepMap gene effect", depmap_gene_effect_path))
    if cellxgene_path is not None:
        source_files.append(source_file_record_local("CELLxGENE signatures", cellxgene_path))
    if marker_path is not None:
        source_files.append(source_file_record_local("single-cell marker signatures", marker_path))

    target_sample_limit = max(0, int(args.manifest_target_sample_limit))
    target_sample = [
        {
            "target_uid": row.get("node_uid", ""),
            "target_symbol": row.get("canonical_name", ""),
            "target_name": row.get("display_name", ""),
            "chembl_ids": extract_xrefs(row, "CHEMBL"),
            "uniprot_ids": extract_xrefs(row, "UNIPROT_SWISSPROT"),
        }
        for row in targets[:target_sample_limit]
    ]

    manifest = {
        "release_id": args.release_id,
        "created_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "output_dir": str(output_dir),
        "target_count": len(targets),
        "target_sample": target_sample,
        "target_sample_limit": target_sample_limit,
        "context_terms": context_terms,
        "source_status": source_status,
        "source_files": [item for item in source_files if item is not None],
        "tables": [
            {"table": "drug_targets", "path": str(output_dir / "drug_targets.csv"), "rows": len(drug_rows)},
            {"table": "cell_type_signatures", "path": str(output_dir / "cell_type_signatures.csv"), "rows": len(cell_rows)},
        ],
        "notes": [
            "Rows are prediction overlays, not canonical graph facts.",
            "ChEMBL rows are fetched from the ChEMBL webresource mechanism table for selected ranked/declared targets.",
            "DrugCentral rows use the public drug.target.interaction TSV snapshot.",
            "DGIdb rows use the public GraphQL API and retain aggregated source-specific license boundaries in license_id.",
            "Open Targets rows use the 26.03 clinical_target dataset joined to drug_molecule metadata.",
            "Cellosaurus rows provide context-only cell-line candidates from either the API search result or the full text export unless richer signatures are supplied through DepMap/CELLxGENE files.",
            "DepMap gene-effect and single-cell marker inputs are normalized into target_symbols/target_uids signature rows while retaining local file hashes in source_files.",
        ],
    }
    manifest["manifest_hash"] = content_hash(manifest)
    (output_dir / "prediction_overlay_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"release_id": args.release_id, "output_dir": str(output_dir), "tables": manifest["tables"], "source_status": source_status, "manifest_hash": manifest["manifest_hash"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
