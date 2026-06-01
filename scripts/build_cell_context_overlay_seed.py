"""Build a small local seed cell-context overlay for prediction research.

The rows produced here are not canonical biology facts. They are transparent
research-prioritization overlays that connect broad cell contexts to marker
targets and metabolic themes already handled by the prediction engine.
"""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq


DEFAULT_RELEASE_ID = "mvp_20260513T002254"
DEFAULT_OUTPUT_ROOT = "manual_sources/prediction_overlays"

FIELDNAMES = [
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

SEED_SIGNATURES = [
    {
        "cell_type_id": "seed:epithelial_barrier_context",
        "cell_type_name": "Epithelial barrier-like context",
        "target_symbols": ["EPCAM", "KRT8", "KRT18", "CDH1", "CLDN1", "OCLN", "MUC1"],
        "pathway_terms": ["epithelial", "barrier", "tight junction", "glycosylation", "phospholipid", "membrane repair"],
        "metabolite_terms": ["UDP-glucose", "UDP-N-acetylglucosamine", "phosphatidylcholine", "taurine"],
        "context_terms": ["epithelial", "barrier", "membrane", "glycosylation"],
        "cell_state": "barrier_membrane_repair",
        "confidence": 0.58,
    },
    {
        "cell_type_id": "seed:renal_proximal_tubular_context",
        "cell_type_name": "Renal proximal tubular epithelial context",
        "target_symbols": ["SLC22A6", "SLC22A8", "SLC22A12", "SLC5A2", "LRP2", "CUBN", "AQP1"],
        "pathway_terms": ["organic anion", "uremic", "renal transport", "osmolyte", "fatty acid oxidation", "mitochondrial"],
        "metabolite_terms": ["indoxyl sulfate", "creatinine", "urate", "taurine", "myo-inositol", "carnitine"],
        "context_terms": ["kidney", "renal", "proximal tubular epithelial", "organic anion", "uremic toxin"],
        "cell_state": "solute_transport_pressure",
        "confidence": 0.66,
    },
    {
        "cell_type_id": "seed:renal_osmolyte_tubular_context",
        "cell_type_name": "Renal tubular osmolyte homeostasis context",
        "target_symbols": ["SLC12A1", "SLC12A3", "AQP1", "AQP2", "UMOD", "SLC6A6"],
        "pathway_terms": ["osmolyte", "osmoregulation", "taurine", "inositol", "solute transport", "kidney"],
        "metabolite_terms": ["taurine", "myo-inositol", "betaine", "creatinine"],
        "context_terms": ["kidney", "renal tubular", "osmolyte", "solute transport"],
        "cell_state": "osmolyte_homeostasis",
        "confidence": 0.6,
    },
    {
        "cell_type_id": "seed:glioma_glial_context",
        "cell_type_name": "Glioma/glial glutamate-GABA context",
        "target_symbols": ["GFAP", "SLC1A2", "SLC1A3", "GLUL", "GAD1", "GAD2", "IDH1", "IDH2", "AQP4"],
        "pathway_terms": ["glutamate", "glutamine", "GABA", "glial", "astrocyte", "IDH", "2-hydroxyglutarate"],
        "metabolite_terms": ["glutamate", "glutamine", "GABA", "2-hydroxyglutarate", "alpha-ketoglutarate", "myo-inositol"],
        "context_terms": ["glioma", "brain", "glial", "astrocyte", "IDH"],
        "cell_state": "glial_neurotransmitter_metabolism",
        "confidence": 0.66,
    },
    {
        "cell_type_id": "seed:hypoxia_glycolytic_tumor_context",
        "cell_type_name": "Hypoxia glycolytic tumor-like context",
        "target_symbols": ["HIF1A", "SLC2A1", "HK1", "HK2", "LDHA", "CA9", "PDK1"],
        "pathway_terms": ["hypoxia", "glycolysis", "glucose", "lactate", "pyruvate", "hexokinase"],
        "metabolite_terms": ["glucose", "lactate", "pyruvate"],
        "context_terms": ["hypoxia", "tumor", "glycolytic", "cancer"],
        "cell_state": "hypoxia_glycolytic_shift",
        "confidence": 0.62,
    },
    {
        "cell_type_id": "seed:myeloid_inflammasome_context",
        "cell_type_name": "Myeloid inflammasome-like context",
        "target_symbols": ["NLRP3", "IL1B", "CASP1", "PYCARD", "CD68", "ITGAM", "TNF"],
        "pathway_terms": ["inflammasome", "NLRP3", "pyroptosis", "tryptophan", "kynurenine", "redox"],
        "metabolite_terms": ["kynurenine", "tryptophan", "succinate", "glutathione"],
        "context_terms": ["myeloid", "macrophage", "inflammatory", "inflammasome"],
        "cell_state": "inflammatory_activation",
        "confidence": 0.56,
    },
]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def load_target_uid_lookup(graph_dir: Path) -> dict[str, str]:
    path = graph_dir / "nodes.parquet"
    nodes = pq.read_table(path, columns=["node_uid", "node_type", "canonical_name"]).to_pandas()
    targets = nodes[nodes["node_type"].astype(str) == "target"].copy()
    lookup: dict[str, str] = {}
    for row in targets.to_dict(orient="records"):
        symbol = str(row.get("canonical_name") or "").strip().upper()
        uid = str(row.get("node_uid") or "").strip()
        if symbol and uid:
            lookup.setdefault(symbol, uid)
    return lookup


def semicolon(values: list[str]) -> str:
    return ";".join(str(value) for value in values if str(value or ""))


def build_rows(lookup: dict[str, str]) -> tuple[list[dict[str, Any]], dict[str, list[str]]]:
    rows: list[dict[str, Any]] = []
    unresolved: dict[str, list[str]] = {}
    for seed in SEED_SIGNATURES:
        symbols = [str(symbol).upper() for symbol in seed["target_symbols"]]
        target_uids = [lookup[symbol] for symbol in symbols if symbol in lookup]
        missing = [symbol for symbol in symbols if symbol not in lookup]
        if missing:
            unresolved[str(seed["cell_type_id"])] = missing
        rows.append(
            {
                "cell_type_id": seed["cell_type_id"],
                "cell_type_name": seed["cell_type_name"],
                "cell_line": "",
                "target_uids": semicolon(target_uids),
                "target_symbols": semicolon(symbols),
                "pathway_uids": "",
                "pathway_terms": semicolon(seed["pathway_terms"]),
                "metabolite_terms": semicolon(seed["metabolite_terms"]),
                "context_terms": semicolon(seed["context_terms"]),
                "cell_state": seed["cell_state"],
                "confidence": seed["confidence"],
                "source_name": "local_seed_cell_context_overlay",
                "source_record_id": seed["cell_type_id"],
                "evidence_level": "research_context_seed",
                "license_id": "local_research_seed:cell_context",
                "source_url": "local_seed",
                "source_version": "cell_context_seed_v1_20260519",
            }
        )
    return rows, unresolved


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build local seed cell-context overlay CSV.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--release-id", default=DEFAULT_RELEASE_ID)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--graph-root", default="graph_projection")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    graph_dir = workspace / args.graph_root / args.release_id
    output_dir = workspace / args.output_root / args.release_id
    target_lookup = load_target_uid_lookup(graph_dir)
    rows, unresolved = build_rows(target_lookup)
    output_path = output_dir / "cell_type_signatures.csv"
    write_csv(output_path, rows)
    report = {
        "created_at_utc": utc_now(),
        "release_id": args.release_id,
        "scope": "local_seed_cell_context_overlay",
        "rows": len(rows),
        "resolved_target_uid_count": sum(len(str(row.get("target_uids") or "").split(";")) for row in rows if row.get("target_uids")),
        "unresolved_symbols": unresolved,
        "outputs": {"cell_type_signatures": str(output_path)},
        "notes": [
            "Rows are research-prioritization overlays, not canonical graph facts.",
            "Target UIDs are resolved by exact target canonical symbol from graph_projection nodes.",
            "Use external DepMap/CELLxGENE signatures to replace or extend these local seeds when available.",
        ],
    }
    (output_dir / "cell_context_overlay_seed_manifest.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
