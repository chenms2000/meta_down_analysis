import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_phase15_validation.py"
SPEC = importlib.util.spec_from_file_location("run_phase15_validation", SCRIPT)
runner = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


class FakeConfig:
    def as_dict(self):
        return {"fake": True}


class FakeService:
    def __init__(self, root: Path):
        self.release_id = "mvp_fake"
        self.config = FakeConfig()
        self.compound_dir = root / "compound_match_index" / self.release_id
        self.compound_dir.mkdir(parents=True, exist_ok=True)

    def release_meta(self):
        return {"release_id": self.release_id, "config_hash": "fake"}

    def envelope(self, endpoint, request, payload):
        body = {
            "api_version": "mvp.v1",
            "endpoint": endpoint,
            "release": self.release_meta(),
            "determinism": {"request_hash": runner.content_hash(request)},
            **payload,
        }
        body["determinism"]["response_hash"] = runner.content_hash({k: v for k, v in body.items() if k != "determinism"})
        return body

    def resolve(self, query, entity_type=None):
        candidate = {
            "entity_uid": "met_glucose",
            "entity_type": "metabolite",
            "display_name": "D-glucose",
            "primary_external_id": "CHEBI:17234",
            "external_xrefs": ["HMDB:HMDB0000122", "WIKIPEDIA.EN:Glucose"],
            "score": 100.0,
            "score_components": {"identifier": 100.0},
        }
        return self.envelope(
            "/resolve",
            {"query": query, "entity_type": entity_type},
            {"status": "matched", "top_score": 100.0, "top_margin": 100.0, "candidates": [candidate]},
        )

    def precheck_metabolites(self, records):
        record = records[0]
        if isinstance(record, dict) and ("formula" in record or "mz" in record):
            row = {
                "input_id": "row_0",
                "record": record,
                "resolution": {
                    "status": "ambiguous",
                    "top_score": 55.0,
                    "top_margin": 0.0,
                    "candidates": [
                        {
                            "entity_uid": "met_glucose",
                            "entity_type": "metabolite",
                            "display_name": "D-glucose",
                            "score": 55.0,
                            "score_components": {"mass": 50.0, "formula": 20.0},
                        }
                    ],
                },
            }
            payload = {"summary": {"matched": 0, "ambiguous": 1, "unmatched": 0, "invalid": 0}, "matched": [], "ambiguous": [row], "unmatched": [], "invalid": []}
        else:
            row = {
                "input_id": "row_0",
                "metabolite_uid": "met_glucose",
                "record": record,
                "resolution": {
                    "status": "matched",
                    "top_score": 100.0,
                    "top_margin": 100.0,
                    "candidates": [
                        {
                            "entity_uid": "met_glucose",
                            "entity_type": "metabolite",
                            "display_name": "D-glucose",
                            "external_xrefs": ["HMDB:HMDB0000122", "WIKIPEDIA.EN:Glucose"],
                            "score": 100.0,
                            "score_components": {"identifier": 100.0},
                        }
                    ],
                },
            }
            payload = {"summary": {"matched": 1, "ambiguous": 0, "unmatched": 0, "invalid": 0}, "matched": [row], "ambiguous": [], "unmatched": [], "invalid": []}
        return self.envelope("/precheck/metabolites", {"records": records}, payload)

    def analyze_metabolites(self, records, max_paths=None, max_hops=None):
        evidence_ref = {
            "ref_type": "literature_support",
            "edge_uid": "edge_fake",
            "support_uid": "litsup_fake",
            "support_class": "confirm",
            "p_literature": 0.4,
            "pmids": ["123"],
            "pmcids": ["PMC123"],
            "relation_uids": ["litrel_fake"],
            "sentence_uids": ["sent_fake"],
            "license_id": "local_articles:test",
        }
        claim_refs = {
            "path_id": "path_1",
            "edge_uids": ["edge_fake"],
            "evidence_ref_uids": ["litsup_fake"],
            "source_records": ["Reactome:edge_fake"],
            "traceability_passed": True,
        }
        analysis_pack = {
            "contract_version": "analysis_pack.v1",
            "input_summary": {
                "analysis_mode": "metabolite_table",
                "input_count": len(records),
                "matched_count": 2,
                "ambiguous_count": 0,
                "unmatched_count": 0,
                "invalid_count": 0,
                "unique_matched_metabolite_count": 2,
                "duplicate_matched_metabolite_uids": [],
                "feature_summary": {
                    "matched_records": 2,
                    "directions": {"up": 2, "down": 0, "unchanged": 0, "unknown": 0},
                    "p_user": {"max": 0.5, "mean": 0.5},
                },
                "seed_weights": {"met_glucose": 1.5, "met_pyruvate": 1.5},
            },
            "matched": [{"input_id": "row_0", "metabolite_uid": "met_glucose", "resolution_status": "matched"}],
            "ambiguous": [],
            "unmatched": [],
            "biological_entity_pools": [],
            "quality_warnings": [],
            "pathway_rankings": [
                {
                    "rank": 1,
                    "pathway_uid": "path_glycolysis",
                    "display_name": "Glycolysis",
                    "score": 3.0,
                    "score_components": {"enrichment_score": 2.0},
                    "prediction_task": "pathway_prediction",
                    "confidence_tier": "medium",
                    "calibrated_confidence": 0.5,
                    "calibration_status": "calibrated_in_scope",
                    "boundary": "Medium-confidence research-prioritization row for validation fixtures.",
                    "claim_refs": claim_refs,
                    "evidence_refs": [evidence_ref],
                }
            ],
            "target_rankings": [
                {
                    "rank": 1,
                    "target_uid": "target_hk1",
                    "display_name": "HK1",
                    "score": 0.2,
                    "score_components": {"propagation_score": 0.2},
                    "prediction_task": "target_prediction",
                    "confidence_tier": "exploratory",
                    "calibrated_confidence": 0.2,
                    "calibration_status": "exploratory_in_scope",
                    "boundary": "Exploratory target-prioritization row for validation fixtures.",
                    "claim_refs": claim_refs,
                    "evidence_refs": [evidence_ref],
                }
            ],
            "disease_rankings": [
                {
                    "rank": 1,
                    "disease_uid": "disease_cancer",
                    "display_name": "Cancer",
                    "score": 0.1,
                    "score_components": {"propagation_score": 0.1},
                    "prediction_task": "disease_prediction",
                    "confidence_tier": "low",
                    "calibrated_confidence": 0.1,
                    "calibration_status": "appendix_low",
                    "boundary": "Low-confidence disease-prioritization row for validation fixtures.",
                    "claim_refs": claim_refs,
                    "evidence_refs": [evidence_ref],
                }
            ],
            "prediction_model": {
                "assessment": {"overall_calibration": {"grade": "B", "status": "usable_with_review"}},
                "predictions": [],
            },
            "top_explanation_paths": [
                {
                    "path_id": "path_1",
                    "cost": 1.0,
                    "path_confidence": 0.5,
                    "terminal_node_uid": "target_hk1",
                    "claim_refs": claim_refs,
                    "evidence_refs": [evidence_ref],
                }
            ],
            "literature_evidence_pack": {
                "support_count": 1,
                "max_p_literature": 0.4,
                "supported_pmids": ["123"],
                "support_classes": ["confirm"],
                "evidence_refs": [evidence_ref],
            },
            "evidence_refs": [evidence_ref],
            "release": self.release_meta(),
            "blocked_reasons": [],
        }
        analysis_pack["determinism"] = {
            "input_hash": runner.content_hash({"records": records, "max_paths": max_paths, "max_hops": max_hops}),
            "release_id": self.release_id,
            "config_hash": "fake",
            "contract_hash": "fake_contract",
        }
        analysis_pack["determinism"]["analysis_pack_hash"] = runner.content_hash({key: value for key, value in analysis_pack.items() if key != "determinism"})
        payload = {
            "analysis_pack": analysis_pack,
            "precheck": {"summary": {"matched": 2, "ambiguous": 0, "unmatched": 0, "invalid": 0}},
            "analysis_features": {
                "summary": {
                    "matched_records": 2,
                    "directions": {"up": 2, "down": 0, "unchanged": 0, "unknown": 0},
                    "p_user": {"max": 0.5, "mean": 0.5},
                }
            },
            "rankings": {
                "pathways": [{"pathway_uid": "path_glycolysis", "name": "Glycolysis", "score": 3.0, "score_components": {"enrichment_score": 2.0}}],
                "targets": [{"target_uid": "target_hk1", "display_name": "HK1", "score": 0.2, "score_components": {"propagation_score": 0.2}}],
                "diseases": [{"disease_uid": "disease_cancer", "display_name": "Cancer", "score": 0.1, "score_components": {"propagation_score": 0.1}}],
            },
            "directional_enrichment": {
                "summary": {"up": 2, "down": 0, "significant_up": 1, "significant_down": 0},
                "up": [{"pathway_uid": "path_glycolysis", "name": "Glycolysis", "score": 3.0}],
                "down": [],
            },
            "propagation": {
                "formula": "pi=(1-alpha)y+alpha*W*pi",
                "alpha": 0.85,
                "iterations": 3,
                "converged": True,
                "graph": {
                    "node_count": 4,
                    "edge_count": 3,
                    "truncated": {"nodes": False, "edges": False},
                    "compression": {"compressed": True, "retained_mass": 0.99},
                    "compressed_not_truncated": True,
                },
            },
            "explanation_paths": [{"path_id": "path_1"}],
            "compressed_explanation_paths": {
                "mode": "path_signature_compression",
                "input_path_count": 1,
                "signature_count": 1,
                "signatures": [{"signature_id": "sig_1", "path_count": 1}],
            },
        }
        return self.envelope("/analyze/metabolites", {"records": records, "max_paths": max_paths, "max_hops": max_hops}, payload)

    def evidence(
        self,
        edge_uid="",
        subject_uid="",
        predicate="",
        object_uid="",
        entity_uid="",
        sentence_uid="",
        relation_uid="",
        limit=25,
    ):
        request = {
            "edge_uid": edge_uid,
            "subject_uid": subject_uid,
            "predicate": predicate,
            "object_uid": object_uid,
            "entity_uid": entity_uid,
            "sentence_uid": sentence_uid,
            "relation_uid": relation_uid,
            "limit": limit,
        }
        payload = {
            "status": "found",
            "manifest": {"manifest_hash": "lit_hash", "metrics": {"evidence_candidate_count": 7}},
            "support": [
                {
                    "support_uid": "litsup_fake",
                    "subject_uid": "target_hk1",
                    "predicate": "target_associated_with_disease",
                    "object_uid": "disease_cancer",
                    "support_class": "confirm",
                    "supported_existing_edge_uids": ["edge_fake"],
                    "evidence_relation_uids": ["litrel_fake"],
                    "sentence_uids": ["sent_fake"],
                    "pmids": ["123"],
                    "pmcids": ["PMC123"],
                    "p_literature": 0.4,
                    "score_components": {
                        "p_literature_formula": "1-prod(1-p_sentence_i)",
                        "sentence_probabilities": [0.4],
                    },
                    "license_id": "local_articles:test",
                    "parser_hash": "lit_parser",
                    "config_hash": "lit_config",
                }
            ],
            "relation_candidates": [
                {
                    "relation_uid": "litrel_fake",
                    "subject_uid": "target_hk1",
                    "predicate": "target_associated_with_disease",
                    "object_uid": "disease_cancer",
                    "sentence_uid": "sent_fake",
                    "pmid": "123",
                    "pmcid": "PMC123",
                    "calibrated_prob": 0.4,
                    "license_id": "local_articles:test",
                    "parser_hash": "lit_parser",
                    "config_hash": "lit_config",
                }
            ],
            "sentences": [{"sentence_uid": "sent_fake", "pmid": "123", "pmcid": "PMC123", "sentence_text": "HK1 is associated with cancer."}],
            "mentions": [{"mention_uid": "mention_fake", "sentence_uid": "sent_fake", "entity_uid": "target_hk1"}],
            "score_notes": {
                "p_literature": "p_literature(edge)=1-prod(1-p_sentence_i)",
                "p_final": "p_final(e)=1-(1-p_curated)(1-p_literature)(1-p_topology)(1-p_user)",
                "evidence_role": "Literature evidence is a validation and calibration overlay; it does not overwrite curated graph facts.",
            },
        }
        return self.envelope("/evidence", request, payload)


class Phase15ValidationRunnerTests(unittest.TestCase):
    def test_builds_deterministic_report_with_blocked_manual_features(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture_dir = root / "fixtures"
            fixture_dir.mkdir()
            (fixture_dir / "common.json").write_text(
                json.dumps(
                    {
                        "fixture_version": "phase1.5.20260513",
                        "release_id": "mvp_fake",
                        "purpose": "common resolver",
                        "cases": [
                            {
                                "case_id": "glucose_hmdb",
                                "input": {"HMDB": "HMDB0000122"},
                                "expected": {"status": "matched", "entity_type": "metabolite", "label": "glucose"},
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            (fixture_dir / "analyze.json").write_text(
                json.dumps(
                    {
                        "fixture_version": "phase1.5.20260513",
                        "release_id": "mvp_fake",
                        "purpose": "analysis",
                        "cases": [
                            {
                                "case_id": "glycolysis",
                                "input": [{"HMDB": "HMDB0000122", "log2FC": 1.0}, {"name": "pyruvate", "log2FC": 1.0}],
                                "expected": {"top_pathway_label": "glycolysis", "directional_enrichment": "up", "requires_explanation_paths": True},
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            (fixture_dir / "rt_ms2.json").write_text(
                json.dumps(
                    {
                        "fixture_version": "phase1.5.20260513",
                        "release_id": "mvp_fake",
                        "purpose": "manual features",
                        "manual_feature_tables": ["compound_rt.csv", "compound_ms2.csv"],
                        "cases": [{"case_id": "rt_ms2", "input": {"formula": "C6H12O6"}, "expected": {"status": "matched"}}],
                    }
                ),
                encoding="utf-8",
            )
            (fixture_dir / "literature_evidence_cases.json").write_text(
                json.dumps(
                    {
                        "fixture_version": "phase2_lite.20260513",
                        "release_id": "mvp_fake",
                        "endpoint": "/evidence",
                        "cases": [
                            {
                                "case_id": "edge_trace",
                                "input": {"edge_uid": "edge_fake", "limit": 5},
                                "expected": {
                                    "status": "found",
                                    "support_classes": ["confirm"],
                                    "requires_sentence_trace": True,
                                    "requires_mentions": True,
                                    "requires_relation_candidate_provenance": True,
                                    "requires_p_literature_contract": True,
                                },
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )

            service = FakeService(root)
            evidence_dir = root / "literature_evidence" / service.release_id
            evidence_dir.mkdir(parents=True)
            (evidence_dir / "literature_evidence_manifest.json").write_text(
                json.dumps(
                    {
                        "manifest_hash": "lit_hash",
                        "metrics": {
                            "evidence_candidate_count": 7,
                            "normalized_mention_rate": 1.0,
                            "supported_existing_edge_count": 2,
                            "novel_candidate_count": 4,
                            "conflict_candidate_count": 1,
                            "sentence_mention_count": 11,
                            "edge_support_count": 6,
                        },
                    }
                ),
                encoding="utf-8",
            )
            fixtures = runner.load_fixtures(fixture_dir)
            first = runner.build_report(service, fixtures, fixture_dir, max_paths=5, max_hops=4, literature_root=root / "literature_evidence")
            second = runner.build_report(service, fixtures, fixture_dir, max_paths=5, max_hops=4, literature_root=root / "literature_evidence")

            self.assertEqual(first["report_hash"], second["report_hash"])
            self.assertEqual(first["literature_evidence"]["evidence_candidate_count"], 7)
            self.assertEqual(first["literature_evidence"]["supported_existing_edge_count"], 2)
            self.assertEqual(first["evidence_validation"]["status"], "passed")
            self.assertTrue(first["evidence_validation"]["contract_passed"])
            self.assertEqual(first["optimization_summary"]["entity_resolution"]["status"], "passed")
            self.assertEqual(first["optimization_summary"]["ranking_calibration"]["status"], "passed")
            self.assertTrue(first["optimization_summary"]["evidence_precision"]["checks"]["p_literature_formula_contract"])
            self.assertEqual(first["summary"]["passed"], 3)
            self.assertEqual(first["summary"]["blocked"], 1)
            self.assertEqual(first["summary"]["failed"], 0)
            self.assertEqual(first["summary"]["gate_status"], "pass_with_known_blocks")
            self.assertTrue(first["summary"]["gate_passed"])

            report_path = runner.write_report(first, root / "validation_reports", service.release_id)
            self.assertTrue(report_path.exists())
            saved = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(saved["report_hash"], first["report_hash"])


if __name__ == "__main__":
    unittest.main()
