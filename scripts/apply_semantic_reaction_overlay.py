#!/usr/bin/env python
"""Apply RDF/BioPAX/SBML semantic reaction sides to an existing v2 store.

This is the lightweight path for improving mechanism fact coverage without
rebuilding the full entity/source store. It parses cached semantic sources,
appends only stable-ID-aligned side participants, rebuilds analysis_view, and
refreshes the manifest.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyarrow.parquet as pq

import database_accuracy_v2 as db


def refresh_manifest(output_dir: Path, release_id: str, normalized_dir: Path, literature_dir: Path) -> Path:
    audits = []
    for parquet_path in sorted(output_dir.glob("*/*.parquet")):
        table_name = parquet_path.stem
        try:
            schema = db.table_schema(table_name)
            schema_hash = db.short_hash(str(schema), 16)
        except Exception:
            schema_hash = ""
        audits.append(
            {
                "table": table_name,
                "path": str(parquet_path),
                "rows": pq.read_table(parquet_path).num_rows,
                "schema_hash": schema_hash,
            }
        )
    manifest = {
        "contract_version": "database_accuracy_store.v2",
        "release_id": release_id,
        "parser_name": "database_accuracy_v2",
        "parser_hash": db.short_hash(Path(db.__file__).read_text(encoding="utf-8"), 16),
        "table_count": len(audits),
        "total_rows": sum(int(row.get("rows") or 0) for row in audits),
        "source": {
            "normalized_dir": str(normalized_dir),
            "literature_dir": str(literature_dir),
        },
        "tables": audits,
        "notes": [
            "Semantic reaction overlay applied incrementally; v1 analysis_pack fields remain additive-compatible.",
            "RDF/BioPAX/SBML explicit side participants take precedence over Rhea SMILES and Reactome PE participation.",
        ],
    }
    manifest["manifest_hash"] = db.short_hash(manifest, 16)
    manifest_path = output_dir / "database_accuracy_store_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest_path


def apply_semantic_overlay(
    workspace: Path,
    release_id: str,
    normalized_root: str = db.DEFAULT_NORMALIZED_ROOT,
    literature_root: str = db.DEFAULT_LITERATURE_ROOT,
    output_root: str = db.DEFAULT_OUTPUT_ROOT,
) -> Path:
    normalized_dir = workspace / normalized_root / release_id
    literature_dir = workspace / literature_root / release_id
    output_dir = workspace / output_root / release_id
    if not output_dir.exists():
        raise FileNotFoundError(f"database_accuracy_store release not found: {output_dir}")

    relation_root = output_dir / "relation_store"
    entity_root = output_dir / "entity_store"
    reactions = db.read_rows(normalized_dir / "reactions.parquet")
    reaction_by_external_base = {db.strip_reactome_version(row.get("primary_external_id")): row for row in reactions}
    xrefs = db.read_rows(entity_root / "chemical_xrefs.parquet")
    chebi_to_chemical: dict[str, str] = {}
    for row in xrefs:
        if str(row.get("xref_source") or "").casefold() == "chebi":
            key = db.chebi_key(row.get("xref_id") or row.get("xref_key"))
            if key:
                chebi_to_chemical.setdefault(key, str(row.get("chemical_uid") or ""))

    participants_v2 = db.read_rows(relation_root / "reaction_participants_v2.parquet")
    reaction_side_participants = db.read_rows(relation_root / "reaction_side_participants.parquet")
    reaction_equations = db.read_rows(relation_root / "reaction_equations.parquet")
    reaction_xrefs = db.read_rows(relation_root / "reaction_xrefs.parquet")

    participant_seen = {
        (str(row.get("reaction_uid") or ""), str(row.get("chemical_uid") or row.get("physical_entity_id") or ""), str(row.get("participant_role") or ""), str(row.get("source_record_uid") or ""))
        for row in participants_v2
    }
    side_seen = {
        (str(row.get("reaction_uid") or ""), str(row.get("chemical_uid") or ""), str(row.get("side") or ""), str(row.get("source_record_uid") or ""))
        for row in reaction_side_participants
    }
    equation_seen = {
        (str(row.get("reaction_uid") or ""), str(row.get("source_reaction_id") or row.get("source_record_uid") or ""))
        for row in reaction_equations
    }
    xref_seen = {
        (str(row.get("reaction_uid") or ""), str(row.get("xref_source") or ""), str(row.get("xref_id") or ""))
        for row in reaction_xrefs
    }

    rhea_tsv_path = workspace / "raw_lake" / "rhea" / release_id / "rhea-tsv.tar.gz"
    direction_rows = db.tsv_records(db.read_tsv_from_tar(rhea_tsv_path, "tsv/rhea-directions.tsv"))
    direction_by_master = {
        str(row.get("RHEA_ID_MASTER") or ""): {
            "lr": str(row.get("RHEA_ID_LR") or ""),
            "rl": str(row.get("RHEA_ID_RL") or ""),
            "bi": str(row.get("RHEA_ID_BI") or ""),
        }
        for row in direction_rows
        if row.get("RHEA_ID_MASTER")
    }
    direction_by_rhea: dict[str, dict[str, str]] = {}
    for master, mapping in direction_by_master.items():
        if master:
            direction_by_rhea[master] = {**mapping, "direction": "UN", "master": master}
        for direction, rhea_id in (("LR", mapping.get("lr", "")), ("RL", mapping.get("rl", "")), ("BI", mapping.get("bi", ""))):
            if rhea_id:
                direction_by_rhea[rhea_id] = {**mapping, "direction": direction, "master": master}

    rhea_to_reaction_uids: dict[str, list[str]] = {}
    for row in db.tsv_records(db.read_tsv_from_tar(rhea_tsv_path, "tsv/rhea2reactome.tsv")):
        rhea_id = str(row.get("RHEA_ID") or "").strip()
        reactome_id = db.strip_reactome_version(row.get("ID"))
        reaction_uid = str((reaction_by_external_base.get(reactome_id) or {}).get("reaction_uid") or "")
        if rhea_id and reaction_uid:
            rhea_to_reaction_uids.setdefault(rhea_id, []).append(reaction_uid)

    paths = db.semantic_source_paths(workspace, release_id)
    for semantic_path in paths["rhea_biopax"]:
        if semantic_path.exists():
            if "rhea.rdf" in semantic_path.name.casefold():
                db.parse_rhea_rdf_semantic_sides(
                    semantic_path,
                    release_id=release_id,
                    chebi_to_chemical=chebi_to_chemical,
                    rhea_to_reaction_uids=rhea_to_reaction_uids,
                    direction_by_rhea=direction_by_rhea,
                    reaction_equations=reaction_equations,
                    reaction_side_participants=reaction_side_participants,
                    participants_v2=participants_v2,
                    reaction_xrefs=reaction_xrefs,
                    equation_seen=equation_seen,
                    side_seen=side_seen,
                    participant_seen=participant_seen,
                    xref_seen=xref_seen,
                )
            else:
                db.parse_biopax_semantic_sides(
                    semantic_path,
                    release_id=release_id,
                    chebi_to_chemical=chebi_to_chemical,
                    reaction_by_external_base=reaction_by_external_base,
                    rhea_to_reaction_uids=rhea_to_reaction_uids,
                    direction_by_rhea=direction_by_rhea,
                    reaction_equations=reaction_equations,
                    reaction_side_participants=reaction_side_participants,
                    participants_v2=participants_v2,
                    reaction_xrefs=reaction_xrefs,
                    equation_seen=equation_seen,
                    side_seen=side_seen,
                    participant_seen=participant_seen,
                    xref_seen=xref_seen,
                )
            break

    audits = [
        db.write_table(relation_root / "reaction_participants_v2.parquet", participants_v2, "reaction_participants_v2"),
        db.write_table(relation_root / "reaction_side_participants.parquet", reaction_side_participants, "reaction_side_participants"),
        db.write_table(relation_root / "reaction_equations.parquet", reaction_equations, "reaction_equations"),
        db.write_table(relation_root / "reaction_xrefs.parquet", reaction_xrefs, "reaction_xrefs"),
    ]
    audits.extend(db.build_analysis_view(output_dir))
    return refresh_manifest(output_dir, release_id, normalized_dir, literature_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply semantic RDF/BioPAX/SBML reaction overlay to an existing v2 store.")
    parser.add_argument("--workspace", type=Path, default=Path("."))
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--normalized-root", default=db.DEFAULT_NORMALIZED_ROOT)
    parser.add_argument("--literature-root", default=db.DEFAULT_LITERATURE_ROOT)
    parser.add_argument("--output-root", default=db.DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    manifest = apply_semantic_overlay(
        args.workspace.resolve(),
        args.release_id,
        normalized_root=args.normalized_root,
        literature_root=args.literature_root,
        output_root=args.output_root,
    )
    print(json.dumps({"database_accuracy_store_manifest": str(manifest)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
