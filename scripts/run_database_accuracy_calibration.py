from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

from metabo_service import MetaboService


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def decision_status(precheck: dict[str, Any], bucket: str = "ambiguous") -> str:
    rows = precheck.get(bucket) or []
    if not rows:
        return ""
    return str(((rows[0].get("identity_resolution_v2") or {}).get("decision") or {}).get("decision_status") or "")


def run_calibration(workspace: Path, release_id: str) -> dict[str, Any]:
    service = MetaboService(workspace, release_id=release_id)
    checks: list[dict[str, Any]] = []

    def add_check(check_id: str, passed: bool, details: dict[str, Any]) -> None:
        checks.append({"check_id": check_id, "passed": bool(passed), "details": details})

    ratio = service.precheck_metabolites(
        [{"trait": "GCST90201000", "reported_trait": "Tryptophan to Pyruvate ratio", "log2FC": 1.0}]
    )
    ratio_decisions = [
        ((row.get("identity_resolution_v2") or {}).get("decision") or {}).get("decision_status")
        for bucket in ("matched", "ambiguous", "unmatched", "invalid")
        for row in ratio.get(bucket, []) or []
    ]
    add_check(
        "ratio_trap_not_exact_abundance",
        "accepted_exact" not in ratio_decisions,
        {"decision_statuses": ratio_decisions, "summary": ratio.get("summary", {})},
    )

    lipid_class = service.precheck_metabolites([{"name": "Sphingomyelin"}])
    class_decisions = [
        ((row.get("identity_resolution_v2") or {}).get("decision") or {}).get("decision_status")
        for bucket in ("matched", "ambiguous", "unmatched", "invalid")
        for row in lipid_class.get(bucket, []) or []
    ]
    add_check(
        "class_trap_not_exact_chemical",
        "accepted_exact" not in class_decisions,
        {"decision_statuses": class_decisions, "summary": lipid_class.get("summary", {})},
    )

    zero_score_rows = [
        row
        for bucket in ("matched", "ambiguous", "unmatched", "invalid")
        for row in lipid_class.get(bucket, []) or []
        for candidate in ((row.get("identity_resolution_v2") or {}).get("candidates") or [])
        if float(candidate.get("match_score") or 0.0) <= 0.0
    ]
    add_check(
        "zero_score_candidates_blocked",
        all(
            candidate.get("blocking_reason") == "zero_score_candidate"
            for row in zero_score_rows
            for candidate in ((row.get("identity_resolution_v2") or {}).get("candidates") or [])
            if float(candidate.get("match_score") or 0.0) <= 0.0
        ),
        {"zero_score_input_count": len(zero_score_rows)},
    )

    analyzed = service.analyze_metabolites([{"name": "Glucose", "log2FC": 1.0, "padj": 0.01}], max_paths=5, max_hops=2)
    pack = analyzed.get("analysis_pack") or {}
    database_accuracy = pack.get("database_accuracy") or {}
    facts = database_accuracy.get("mechanism_ready_facts") or []
    add_check(
        "database_accuracy_pack_present",
        database_accuracy.get("contract_version") == "database_accuracy_store.v2",
        {"store_available": database_accuracy.get("store_available"), "identity_decision_summary": database_accuracy.get("identity_decision_summary", {})},
    )
    add_check(
        "mechanism_ready_facts_traceable_when_present",
        all(fact.get("source_record_uid") or fact.get("evidence_assertion_uid") for fact in facts),
        {"fact_count": len(facts)},
    )

    store_root = workspace / "database_accuracy_store" / release_id
    facts_path = store_root / "analysis_view" / "mechanism_ready_facts.parquet"
    candidates_path = store_root / "analysis_view" / "fact_candidates.parquet"
    store_facts = pq.read_table(facts_path).to_pylist() if facts_path.exists() else []
    store_candidates = pq.read_table(candidates_path).to_pylist() if candidates_path.exists() else []
    relation_counts = {}
    for table_name in (
        "reaction_equations",
        "reaction_side_participants",
        "reaction_xrefs",
        "enzyme_reaction_links",
        "reaction_publication_links",
    ):
        path = store_root / "relation_store" / f"{table_name}.parquet"
        relation_counts[table_name] = pq.read_table(path).num_rows if path.exists() else 0
    add_check(
        "equation_level_relation_tables_nonempty",
        relation_counts["reaction_equations"] > 0
        and relation_counts["reaction_side_participants"] > 0
        and relation_counts["enzyme_reaction_links"] > 0,
        relation_counts,
    )
    add_check(
        "store_mechanism_ready_facts_nonempty",
        len(store_facts) > 0,
        {"fact_count": len(store_facts), "candidate_count": len(store_candidates)},
    )
    allowed_scopes = {
        "directional_reaction_fact",
        "enzyme_reaction_fact",
        "transporter_reaction_fact",
        "bidirectional_reaction_fact",
        "exact_reaction_fact",
        "role_unknown_reaction_fact",
    }
    add_check(
        "store_mechanism_ready_facts_have_required_traces",
        all(
            row.get("source_record_uids")
            and row.get("identity_decision_uids")
            and row.get("allowed_claim_scope") in allowed_scopes
            for row in store_facts[:10000]
        ),
        {"checked_fact_count": min(len(store_facts), 10000)},
    )
    directional_count = sum(1 for row in store_facts if row.get("allowed_claim_scope") == "directional_reaction_fact")
    add_check(
        "directional_reaction_facts_present",
        directional_count > 0,
        {"directional_reaction_fact_count": directional_count},
    )
    semantic_sources = {
        "rhea_rdf_explicit_side",
        "rhea_biopax_explicit_side",
        "reactome_sbml_species_reference",
        "reactome_biopax_conversion",
    }
    side_rows_path = store_root / "relation_store" / "reaction_side_participants.parquet"
    side_rows = pq.read_table(side_rows_path).to_pylist() if side_rows_path.exists() else []
    semantic_side_count = sum(1 for row in side_rows if row.get("relation_source") in semantic_sources)
    add_check(
        "semantic_side_participants_present",
        semantic_side_count > 0,
        {"semantic_side_participants": semantic_side_count},
    )
    semantic_fact_count = 0
    semantic_fact_traceable = True
    for row in store_facts:
        try:
            metadata = json.loads(row.get("metadata_json") or "{}")
        except Exception:
            metadata = {}
        if metadata.get("relation_source") in semantic_sources:
            semantic_fact_count += 1
            if not (metadata.get("semantic_source_uri") or row.get("source_record_uids")):
                semantic_fact_traceable = False
    add_check(
        "semantic_mechanism_facts_trace_to_source_uri",
        semantic_fact_count > 0 and semantic_fact_traceable,
        {"semantic_mechanism_fact_count": semantic_fact_count, "semantic_fact_traceable": semantic_fact_traceable},
    )
    add_check(
        "pathway_module_candidates_not_promoted",
        all(row.get("fact_type") != "module_membership" for row in store_facts),
        {
            "module_candidate_count": sum(1 for row in store_candidates if row.get("candidate_fact_type") == "module_membership"),
            "promoted_module_fact_count": sum(1 for row in store_facts if row.get("fact_type") == "module_membership"),
        },
    )
    chemical_path = store_root / "entity_store" / "chemical_entities.parquet"
    chemicals = pq.read_table(chemical_path, columns=["chemical_uid", "canonical_name", "identity_status"]).to_pylist() if chemical_path.exists() else []
    fact_subjects = {str(row.get("subject_uid") or "") for row in store_facts}
    gold_terms = ["glucose", "lactate", "pyruvate", "glutamine", "glutamate", "succinate", "citrate", "cysteine", "glutathione", "carnitine"]
    recovered = {
        term: any(
            term in str(row.get("canonical_name") or "").casefold()
            and row.get("identity_status") == "canonical"
            and row.get("chemical_uid") in fact_subjects
            for row in chemicals
        )
        for term in gold_terms
    }
    add_check(
        "gold_metabolism_terms_recover_reaction_facts",
        sum(1 for ok in recovered.values() if ok) >= 8,
        {"recovered": recovered, "recovered_count": sum(1 for ok in recovered.values() if ok), "gold_count": len(gold_terms)},
    )

    by_status = Counter("passed" if row["passed"] else "failed" for row in checks)
    return {
        "report_version": "database_accuracy_calibration_report.v1",
        "workspace": str(workspace),
        "release_id": release_id,
        "passed": by_status.get("failed", 0) == 0,
        "summary": dict(sorted(by_status.items())),
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run proxy calibration checks for database_accuracy_store.v2.")
    parser.add_argument("--workspace", type=Path, default=Path("."))
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--output", type=Path, default=Path("validation_reports/database_accuracy_calibration_report.json"))
    args = parser.parse_args()
    report = run_calibration(args.workspace.resolve(), args.release_id)
    output = args.output if args.output.is_absolute() else args.workspace.resolve() / args.output
    write_json(output, report)
    print(json.dumps({"calibration_report": str(output), "passed": report["passed"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
