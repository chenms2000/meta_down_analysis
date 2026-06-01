"""Prepare local cell signature source CSVs for prediction overlays.

Inputs are downloaded public source tables:
- DepMap files are already renamed by the download step.
- CellMarker Excel tables are converted into a normalized long marker table.
- CELLxGENE CellGuide marker-gene presence is converted into top human markers
  per cell type.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import time
import urllib.request
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd


CELLGUIDE_BASE = "https://cellguide.cellxgene.cziscience.com"
CELLGUIDE_MARKER_URL_SUFFIX = "computational_marker_genes/marker_gene_presence.json.gz"
STANDARD_COLUMNS = [
    "cell_type_id",
    "cell_type_name",
    "gene",
    "tissue",
    "disease",
    "dataset_id",
    "confidence",
    "avg_log2FC",
    "source_record_id",
    "license_id",
    "source_version",
    "source_url",
]


def utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def file_sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def csv_data_row_count(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("rb") as handle:
        line_count = sum(1 for _line in handle)
    return max(0, line_count - 1)


def download(url: str, path: Path, timeout: int = 180) -> Path:
    if path.exists() and path.stat().st_size > 0:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "metabo-data-end/0.1 cell-signature-prep"})
    with urllib.request.urlopen(request, timeout=timeout) as response, path.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    return path


def as_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if isinstance(value, float) and math.isnan(value):
            return ""
    except TypeError:
        pass
    return str(value).strip()


def normalize_cell_ontology_id(value: Any, fallback: str) -> str:
    text = as_text(value)
    if text.startswith("CL_"):
        return "CL:" + text.split("_", 1)[1]
    return text or fallback


def preserve_and_format_depmap_files(source_dir: Path, raw_dir: Path) -> None:
    raw_depmap_dir = raw_dir / "depmap_26q1"
    raw_depmap_dir.mkdir(parents=True, exist_ok=True)
    model = source_dir / "depmap_model.csv"
    effect = source_dir / "depmap_gene_effect.csv"
    if model.exists() and not (raw_depmap_dir / "Model.csv").exists():
        shutil.copy2(model, raw_depmap_dir / "Model.csv")
    if effect.exists() and not (raw_depmap_dir / "CRISPRGeneEffect.csv").exists():
        shutil.copy2(effect, raw_depmap_dir / "CRISPRGeneEffect.csv")
    if not effect.exists():
        return
    with effect.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline()
    first_column = header.split(",", 1)[0].strip().lstrip("\ufeff")
    if first_column == "ModelID":
        return
    tmp_path = effect.with_suffix(".csv.tmp")
    with effect.open("r", encoding="utf-8", newline="") as source, tmp_path.open("w", encoding="utf-8", newline="") as target:
        first = source.readline()
        target.write("ModelID," + first.split(",", 1)[1])
        shutil.copyfileobj(source, target, length=1024 * 1024)
    tmp_path.replace(effect)


def confidence_from_cellmarker(row: dict[str, Any], seq_file: bool = False) -> float:
    source = as_text(row.get("marker_source")).casefold()
    technology = as_text(row.get("technology_seq")).casefold()
    confidence = 0.58
    if "experiment" in source:
        confidence = 0.72
    elif "review" in source:
        confidence = 0.62
    elif "company" in source:
        confidence = 0.5
    if seq_file or "single-cell" in source or technology:
        confidence = max(confidence, 0.68)
    return confidence


def normalize_cellmarker_rows(paths: list[Path], output: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in paths:
        if not path.exists():
            continue
        frame = pd.read_excel(path, dtype=str).fillna("")
        seq_file = "seq" in path.name.casefold()
        for item in frame.to_dict(orient="records"):
            if as_text(item.get("species")).casefold() != "human":
                continue
            gene = as_text(item.get("Symbol") or item.get("marker")).upper()
            cell_name = as_text(item.get("cell_name") or item.get("cell_type"))
            if not gene or not cell_name:
                continue
            cell_type_id = normalize_cell_ontology_id(item.get("cellontology_id"), f"CellMarker:{cell_name}")
            pmid = as_text(item.get("PMID"))
            rows.append(
                {
                    "cell_type_id": cell_type_id,
                    "cell_type_name": cell_name,
                    "gene": gene,
                    "tissue": as_text(item.get("tissue_type") or item.get("tissue_class")),
                    "disease": as_text(item.get("cancer_type")),
                    "dataset_id": f"PMID:{pmid}" if pmid else "CellMarker",
                    "confidence": confidence_from_cellmarker(item, seq_file=seq_file),
                    "avg_log2FC": "",
                    "source_record_id": "|".join(
                        part
                        for part in [
                            "CellMarker2.0",
                            as_text(item.get("cellontology_id")),
                            cell_name,
                            gene,
                            pmid,
                        ]
                        if part
                    ),
                    "license_id": "public_research:cellmarker",
                    "source_version": "CellMarker_2.0_download_2026-05-19",
                    "source_url": "http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download.html",
                }
            )
    result = pd.DataFrame(rows, columns=STANDARD_COLUMNS)
    if not result.empty:
        result = result.drop_duplicates(["cell_type_id", "gene", "tissue", "disease", "dataset_id"])
        result = result.sort_values(["cell_type_name", "gene", "tissue", "dataset_id"])
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False, encoding="utf-8")
    return result


def load_cellguide_metadata(raw_dir: Path) -> tuple[str, dict[str, Any], dict[str, Any], Path]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    snapshot_path = raw_dir / "latest_snapshot_identifier.txt"
    if snapshot_path.exists() and snapshot_path.read_text(encoding="utf-8").strip():
        snapshot = snapshot_path.read_text(encoding="utf-8").strip()
    else:
        snapshot = urllib.request.urlopen(f"{CELLGUIDE_BASE}/latest_snapshot_identifier", timeout=60).read().decode("utf-8").strip()
        snapshot_path.write_text(snapshot + "\n", encoding="utf-8")
    celltype_path = download(f"{CELLGUIDE_BASE}/{snapshot}/celltype_metadata.json", raw_dir / "celltype_metadata.json")
    tissue_path = download(f"{CELLGUIDE_BASE}/{snapshot}/tissue_metadata.json", raw_dir / "tissue_metadata.json")
    marker_path = download(
        f"{CELLGUIDE_BASE}/{snapshot}/{CELLGUIDE_MARKER_URL_SUFFIX}",
        raw_dir / "marker_gene_presence.json.gz",
        timeout=300,
    )
    return (
        snapshot,
        json.loads(celltype_path.read_text(encoding="utf-8")),
        json.loads(tissue_path.read_text(encoding="utf-8")),
        marker_path,
    )


def normalize_cellxgene_rows(raw_dir: Path, output: Path, top_genes_per_cell_type: int, min_marker_score: float) -> pd.DataFrame:
    snapshot, celltypes, _tissues, marker_path = load_cellguide_metadata(raw_dir)
    payload = json.loads(gzip.decompress(marker_path.read_bytes()))
    best: dict[tuple[str, str], dict[str, Any]] = {}
    tissue_terms: dict[tuple[str, str], set[str]] = defaultdict(set)
    for gene, species_payload in payload.items():
        human_payload = species_payload.get("Homo sapiens")
        if not human_payload:
            continue
        symbol = as_text(gene).upper()
        for tissue, entries in human_payload.items():
            for entry in entries:
                score = float(entry.get("marker_score") or 0.0)
                if score < min_marker_score:
                    continue
                cell_type_id = as_text(entry.get("cell_type_id"))
                if not cell_type_id:
                    continue
                key = (cell_type_id, symbol)
                existing = best.get(key)
                if existing is None or score > float(existing.get("marker_score") or 0.0):
                    best[key] = {**entry, "gene": symbol, "marker_score": score, "tissue": tissue}
                if tissue and tissue != "All Tissues":
                    tissue_terms[key].add(tissue)

    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for (cell_type_id, symbol), entry in best.items():
        grouped[cell_type_id].append({**entry, "gene": symbol, "tissue_terms": sorted(tissue_terms.get((cell_type_id, symbol), set()))})

    rows: list[dict[str, Any]] = []
    for cell_type_id, entries in grouped.items():
        entries.sort(key=lambda row: (-float(row.get("marker_score") or 0.0), row.get("gene", "")))
        meta = celltypes.get(cell_type_id, {})
        cell_name = as_text(meta.get("name")) or cell_type_id
        for entry in entries[:top_genes_per_cell_type]:
            gene = as_text(entry.get("gene"))
            score = float(entry.get("marker_score") or 0.0)
            tissue = ";".join(entry.get("tissue_terms") or []) or as_text(entry.get("tissue"))
            rows.append(
                {
                    "cell_type_id": cell_type_id,
                    "cell_type_name": cell_name,
                    "gene": gene,
                    "tissue": tissue,
                    "disease": "",
                    "dataset_id": f"CELLxGENE_CellGuide:{snapshot}",
                    "confidence": round(min(0.9, 0.5 + min(score, 4.0) * 0.08), 6),
                    "avg_log2FC": round(score, 6),
                    "source_record_id": f"{snapshot}|{cell_type_id}|{gene}",
                    "license_id": "public_research:cellxgene_cellguide",
                    "source_version": f"CELLxGENE_CellGuide_snapshot_{snapshot}",
                    "source_url": f"{CELLGUIDE_BASE}/{snapshot}/{CELLGUIDE_MARKER_URL_SUFFIX}",
                }
            )
    result = pd.DataFrame(rows, columns=STANDARD_COLUMNS)
    if not result.empty:
        result = result.sort_values(["cell_type_name", "avg_log2FC", "gene"], ascending=[True, False, True])
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False, encoding="utf-8")
    return result


def build_manifest(paths: dict[str, Path], row_counts: dict[str, int], output: Path) -> None:
    manifest = {
        "created_at_utc": utc_now(),
        "purpose": "Prepared local source tables for prediction overlay cell signatures.",
        "outputs": {
            name: {
                "path": str(path),
                "rows": int(row_counts.get(name, 0)),
                "bytes": path.stat().st_size if path.exists() else 0,
                "sha256": file_sha256(path) if path.exists() else "",
            }
            for name, path in paths.items()
        },
        "notes": [
            "DepMap files are direct 26Q1 downloads renamed to the local pipeline names.",
            "CELLxGENE markers come from official CellGuide computational marker-gene presence data, limited to human cell types.",
            "CellMarker rows combine Human and Seq downloads and keep PMID provenance when present.",
        ],
    }
    output.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare cell signature CSVs for prediction overlay ingestion.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--top-cellxgene-genes-per-cell-type", type=int, default=50)
    parser.add_argument("--min-cellxgene-marker-score", type=float, default=0.75)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    source_dir = workspace / "manual_sources" / "cell_signatures"
    raw_dir = workspace / "raw_lake" / "cell_signatures"
    preserve_and_format_depmap_files(source_dir, raw_dir)
    cellmarker_paths = [
        raw_dir / "cellmarker" / "Cell_marker_Human.xlsx",
        raw_dir / "cellmarker" / "Cell_marker_Seq.xlsx",
    ]
    single_cell = normalize_cellmarker_rows(cellmarker_paths, source_dir / "single_cell_markers.csv")
    cellxgene = normalize_cellxgene_rows(
        raw_dir / "cellxgene",
        source_dir / "cellxgene_markers.csv",
        top_genes_per_cell_type=args.top_cellxgene_genes_per_cell_type,
        min_marker_score=args.min_cellxgene_marker_score,
    )
    outputs = {
        "depmap_model": source_dir / "depmap_model.csv",
        "depmap_gene_effect": source_dir / "depmap_gene_effect.csv",
        "cellxgene_markers": source_dir / "cellxgene_markers.csv",
        "single_cell_markers": source_dir / "single_cell_markers.csv",
    }
    build_manifest(
        outputs,
        {
            "depmap_model": csv_data_row_count(outputs["depmap_model"]),
            "depmap_gene_effect": csv_data_row_count(outputs["depmap_gene_effect"]),
            "cellxgene_markers": len(cellxgene),
            "single_cell_markers": len(single_cell),
        },
        source_dir / "cell_signature_source_manifest.json",
    )
    print(
        json.dumps(
            {
                "cellxgene_markers": len(cellxgene),
                "single_cell_markers": len(single_cell),
                "manifest": str(source_dir / "cell_signature_source_manifest.json"),
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
