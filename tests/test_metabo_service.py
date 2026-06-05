import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover
    pa = None
    pq = None


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "metabo_service.py"
SPEC = importlib.util.spec_from_file_location("metabo_service", SCRIPT)
metabo_service = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = metabo_service
SPEC.loader.exec_module(metabo_service)

DB_SCRIPT = ROOT / "scripts" / "database_accuracy_v2.py"
DB_SPEC = importlib.util.spec_from_file_location("database_accuracy_v2", DB_SCRIPT)
database_accuracy_v2 = importlib.util.module_from_spec(DB_SPEC)
assert DB_SPEC.loader is not None
sys.modules[DB_SPEC.name] = database_accuracy_v2
DB_SPEC.loader.exec_module(database_accuracy_v2)


def write_parquet(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


@unittest.skipIf(pa is None or pq is None, "pyarrow is required")
class MetaboServiceTests(unittest.TestCase):
    def make_service(self, root: Path):
        release = "mvp_20260101T000000"
        graph = root / "graph_projection" / release
        normalized = root / "normalized_store" / release
        pubchem = root / "pubchem_cid_cache" / release
        common = {"source_release": release, "license_id": "open_core:test"}

        nodes = [
            {
                "node_uid": "met_1",
                "node_idx": 0,
                "node_type": "metabolite",
                "canonical_name": "Glucose",
                "display_name": "Glucose",
                "primary_external_id": "CHEBI:17234",
                "external_xrefs": ["CHEBI:17234", "HMDB:HMDB0000122", "PUBCHEM.COMPOUND:5793"],
                "source_priority": "ChEBI",
                **common,
            },
            {
                "node_uid": "met_2",
                "node_idx": 1,
                "node_type": "metabolite",
                "canonical_name": "Shared",
                "display_name": "Shared",
                "primary_external_id": "CHEBI:1",
                "external_xrefs": ["CHEBI:1"],
                "source_priority": "ChEBI",
                **common,
            },
            {
                "node_uid": "met_3",
                "node_idx": 2,
                "node_type": "metabolite",
                "canonical_name": "Shared",
                "display_name": "Shared",
                "primary_external_id": "CHEBI:2",
                "external_xrefs": ["CHEBI:2"],
                "source_priority": "ChEBI",
                **common,
            },
            {
                "node_uid": "pathway_1",
                "node_idx": 3,
                "node_type": "pathway",
                "canonical_name": "Glycolysis",
                "display_name": "Glycolysis",
                "primary_external_id": "R-HSA-70171",
                "external_xrefs": ["REACTOME:R-HSA-70171"],
                "source_priority": "Reactome",
                **common,
            },
            {
                "node_uid": "reaction_1",
                "node_idx": 4,
                "node_type": "reaction",
                "canonical_name": "Glucose phosphorylation",
                "display_name": "Glucose phosphorylation",
                "primary_external_id": "R-HSA-1",
                "external_xrefs": ["REACTOME:R-HSA-1"],
                "source_priority": "Reactome",
                **common,
            },
            {
                "node_uid": "gene_1",
                "node_idx": 5,
                "node_type": "gene",
                "canonical_name": "HK1",
                "display_name": "HK1",
                "primary_external_id": "ENSEMBL:ENSG1",
                "external_xrefs": ["ENSEMBL:ENSG1"],
                "source_priority": "protein_coding",
                **common,
            },
            {
                "node_uid": "target_1",
                "node_idx": 6,
                "node_type": "target",
                "canonical_name": "HK1",
                "display_name": "hexokinase 1",
                "primary_external_id": "OPENTARGETS:ENSG1",
                "external_xrefs": ["OPENTARGETS:ENSG1"],
                "source_priority": "protein_coding",
                **common,
            },
            {
                "node_uid": "disease_1",
                "node_idx": 7,
                "node_type": "disease",
                "canonical_name": "Cancer",
                "display_name": "Cancer",
                "primary_external_id": "MONDO:0004992",
                "external_xrefs": ["MONDO:0004992"],
                "source_priority": "MONDO",
                **common,
            },
        ]
        write_parquet(graph / "nodes.parquet", nodes)

        resolver_rows = [
            {
                "lookup_uid": "lk_1",
                "entity_uid": "met_1",
                "entity_type": "metabolite",
                "namespace": "HMDB",
                "lookup_key": "hmdb0000122",
                "raw_value": "HMDB0000122",
                "match_field": "xref",
                "rank": 95.0,
                "source_table": "metabolite_xrefs",
                **common,
            },
            {
                "lookup_uid": "lk_2",
                "entity_uid": "met_1",
                "entity_type": "metabolite",
                "namespace": "TEXT",
                "lookup_key": "glucose",
                "raw_value": "Glucose",
                "match_field": "canonical_name",
                "rank": 100.0,
                "source_table": "metabolites",
                **common,
            },
            {
                "lookup_uid": "lk_3",
                "entity_uid": "met_2",
                "entity_type": "metabolite",
                "namespace": "TEXT",
                "lookup_key": "shared",
                "raw_value": "Shared",
                "match_field": "canonical_name",
                "rank": 100.0,
                "source_table": "metabolites",
                **common,
            },
            {
                "lookup_uid": "lk_4",
                "entity_uid": "met_3",
                "entity_type": "metabolite",
                "namespace": "TEXT",
                "lookup_key": "shared",
                "raw_value": "Shared",
                "match_field": "canonical_name",
                "rank": 100.0,
                "source_table": "metabolites",
                **common,
            },
        ]
        write_parquet(graph / "resolver_index.parquet", resolver_rows)

        edge_common = {
            "source_name": "Reactome",
            "source_record_id": "src",
            "evidence_level": "curated",
            "source_release": release,
            "license_id": "open_core:test",
            "parser_hash": "parser",
            "metadata_json": "{}",
        }
        edges = [
            {
                "edge_uid": "edge_mp",
                "edge_type": "metabolite_participates_in_pathway",
                "subject_uid": "met_1",
                "subject_type": "metabolite",
                "predicate": "participates_in",
                "object_uid": "pathway_1",
                "object_type": "pathway",
                "weight": 1.0,
                "source_table": "metabolite_pathway_edges",
                **edge_common,
            },
            {
                "edge_uid": "edge_mr",
                "edge_type": "metabolite_participates_in_reaction",
                "subject_uid": "met_1",
                "subject_type": "metabolite",
                "predicate": "participates_in_reaction",
                "object_uid": "reaction_1",
                "object_type": "reaction",
                "weight": 1.0,
                "source_table": "reaction_participants",
                **edge_common,
            },
            {
                "edge_uid": "edge_rp",
                "edge_type": "reaction_in_pathway",
                "subject_uid": "reaction_1",
                "subject_type": "reaction",
                "predicate": "in_pathway",
                "object_uid": "pathway_1",
                "object_type": "pathway",
                "weight": 1.0,
                "source_table": "reactions",
                **edge_common,
            },
            {
                "edge_uid": "edge_gp",
                "edge_type": "gene_involved_in_pathway",
                "subject_uid": "gene_1",
                "subject_type": "gene",
                "predicate": "involved_in",
                "object_uid": "pathway_1",
                "object_type": "pathway",
                "weight": 1.0,
                "source_table": "gene_pathway_edges",
                **edge_common,
            },
            {
                "edge_uid": "edge_tg",
                "edge_type": "target_maps_to_gene",
                "subject_uid": "target_1",
                "subject_type": "target",
                "predicate": "targets",
                "object_uid": "gene_1",
                "object_type": "gene",
                "weight": 1.0,
                "source_table": "target_gene_edges",
                **edge_common,
            },
            {
                "edge_uid": "edge_td",
                "edge_type": "target_associated_with_disease",
                "subject_uid": "target_1",
                "subject_type": "target",
                "predicate": "associated_with",
                "object_uid": "disease_1",
                "object_type": "disease",
                "weight": 0.5,
                "source_table": "target_disease_edges",
                **{**edge_common, "source_name": "Open Targets", "metadata_json": "{\"score_components_json\":\"{\\\"associationScore\\\":0.5}\"}"},
            },
        ]
        write_parquet(graph / "edges.parquet", edges)
        (graph / "graph_manifest.json").write_text(
            json.dumps({"release_id": release, "node_count": len(nodes), "tables": [{"table": "edges", "rows": len(edges)}]}),
            encoding="utf-8",
        )

        write_parquet(
            normalized / "metabolites.parquet",
            [
                {
                    "metabolite_uid": "met_1",
                    "canonical_name": "Glucose",
                    "synonyms": ["D-Glucose"],
                    "formula": "C6H12O6",
                    "exact_mass": 180.063388,
                    "charge": 0,
                    "inchikey": "WQZGKKKJIJFFOK-UHFFFAOYSA-N",
                    "smiles": "",
                    "external_xrefs": ["CHEBI:17234", "HMDB:HMDB0000122", "PUBCHEM.COMPOUND:5793"],
                    "source_priority": "ChEBI",
                    "checksum": "",
                    **common,
                },
                {
                    "metabolite_uid": "met_sphingomyelin",
                    "canonical_name": "Sphingomyelin",
                    "synonyms": ["sphingomyelin"],
                    "formula": "",
                    "exact_mass": None,
                    "charge": None,
                    "inchikey": "",
                    "smiles": "",
                    "external_xrefs": [],
                    "source_priority": "name_only",
                    "checksum": "",
                    **common,
                },
            ],
        )
        write_parquet(
            normalized / "metabolite_xrefs.parquet",
            [
                {
                    "xref_uid": "xref_1",
                    "metabolite_uid": "met_1",
                    "xref_source": "HMDB",
                    "xref_id": "HMDB0000122",
                    "xref_key": "HMDB:HMDB0000122",
                    "source_name": "HMDB",
                    **common,
                }
            ],
        )
        write_parquet(
            normalized / "sentences.parquet",
            [
                {
                    "sentence_uid": "sent_lit_1",
                    "article_uid": "article_lit_1",
                    "pmid": "123",
                    "pmcid": "PMC123",
                    "section": "abstract",
                    "sentence_text": "HK1 is associated with Cancer in metabolic studies.",
                    "text_hash": "hash",
                    "start_offset": 0,
                    "end_offset": 52,
                    "language": "en",
                    "source_release": release,
                    "parser_hash": "parser",
                }
            ],
        )
        write_parquet(
            normalized / "metabolite_pathway_edges.parquet",
            [
                {
                    "edge_uid": "edge_mp",
                    "metabolite_uid": "met_1",
                    "pathway_uid": "pathway_1",
                    "source_name": "Reactome",
                    "source_record_id": "CHEBI:17234|R-HSA-70171",
                    "evidence_level": "curated",
                    **common,
                },
                {
                    "edge_uid": "edge_mp_bg",
                    "metabolite_uid": "met_2",
                    "pathway_uid": "pathway_2",
                    "source_name": "Reactome",
                    "source_record_id": "CHEBI:1|R-HSA-2",
                    "evidence_level": "curated",
                    **common,
                },
            ],
        )
        write_parquet(
            normalized / "pathways.parquet",
            [
                {
                    "pathway_uid": "pathway_1",
                    "name": "Glycolysis",
                    "primary_external_id": "R-HSA-70171",
                    "source_name": "Reactome",
                    "species": "Homo sapiens",
                },
                {
                    "pathway_uid": "pathway_2",
                    "name": "Background",
                    "primary_external_id": "R-HSA-2",
                    "source_name": "Reactome",
                    "species": "Homo sapiens",
                },
            ],
        )
        write_parquet(normalized / "reactions.parquet", [{"reaction_uid": "reaction_1", "pathway_uid": "pathway_1"}])
        write_parquet(
            normalized / "reaction_participants.parquet",
            [
                {
                    "edge_uid": "edge_mr",
                    "reaction_uid": "reaction_1",
                    "participant_uid": "met_1",
                    "participant_type": "metabolite",
                    "source_name": "Reactome",
                    "source_record_id": "R-HSA-1|CHEBI:17234",
                    **common,
                }
            ],
        )
        (normalized / "normalized_manifest.json").write_text(
            json.dumps({"release_id": release, "tables": [{"table": "metabolite_pathway_edges", "rows": 2}]}),
            encoding="utf-8",
        )

        write_parquet(
            pubchem / "cid_properties.parquet",
            [
                {
                    "pubchem_cid": "5793",
                    "metabolite_uids": ["met_1"],
                    "title": "Glucose",
                    "formula": "C6H12O6",
                    "monoisotopic_mass": 180.063388,
                    "exact_mass": 180.063388,
                    "inchi": "",
                    "inchikey": "WQZGKKKJIJFFOK-UHFFFAOYSA-N",
                    "isomeric_smiles": "",
                    "synonyms": ["D-Glucose"],
                    "synonym_count": 1,
                    "fields_present": ["title", "formula", "inchikey"],
                    "source_release": release,
                    "license_id": "open_core:pubchem",
                    "parser_hash": "parser",
                }
            ],
        )
        write_parquet(
            pubchem / "cid_lookup_index.parquet",
            [
                {
                    "lookup_uid": "pc_1",
                    "pubchem_cid": "5793",
                    "metabolite_uids": ["met_1"],
                    "namespace": "CID",
                    "lookup_key": "5793",
                    "raw_value": "5793",
                    "match_field": "pubchem_cid",
                    "rank": 100.0,
                    "source_release": release,
                    "license_id": "open_core:pubchem",
                    "parser_hash": "parser",
                }
            ],
        )
        write_parquet(pubchem / "cid_missing.parquet", [{"pubchem_cid": "999", "metabolite_uids": ["met_x"], "reason": "not_in_cache"}])
        (pubchem / "pubchem_cid_cache_manifest.json").write_text(
            json.dumps({"release_id": release, "target_cid_count": 1, "matched_cid_count": 1}),
            encoding="utf-8",
        )

        literature = root / "literature_evidence" / release
        write_parquet(
            literature / "literature_edge_support.parquet",
            [
                {
                    "support_uid": "litsup_1",
                    "subject_uid": "target_1",
                    "subject_type": "target",
                    "predicate": "target_associated_with_disease",
                    "object_uid": "disease_1",
                    "object_type": "disease",
                    "polarity_set": ["association"],
                    "support_class": "confirm",
                    "supported_existing_edge_uids": ["edge_td"],
                    "evidence_relation_uids": ["litrel_1"],
                    "sentence_uids": ["sent_lit_1"],
                    "pmids": ["123"],
                    "pmcids": ["PMC123"],
                    "evidence_sentence_count": 1,
                    "distinct_article_count": 1,
                    "p_literature": 0.4,
                    "raw_score_max": 0.7,
                    "calibrated_prob_max": 0.4,
                    "score_components_json": "{\"p_literature_formula\":\"1-prod(1-p_sentence_i)\"}",
                    "license_id": "local_articles:test",
                    "source_release": release,
                    "parser_hash": "lit_parser",
                    "config_hash": "lit_config",
                },
                {
                    "support_uid": "litsup_overlay",
                    "subject_uid": "met_1",
                    "subject_type": "metabolite",
                    "predicate": "metabolite_regulates_gene",
                    "object_uid": "target_1",
                    "object_type": "target",
                    "polarity_set": ["activation"],
                    "support_class": "novel_candidate",
                    "supported_existing_edge_uids": [],
                    "evidence_relation_uids": ["litrel_overlay"],
                    "sentence_uids": ["sent_lit_1"],
                    "pmids": ["123"],
                    "pmcids": ["PMC123"],
                    "evidence_sentence_count": 1,
                    "distinct_article_count": 1,
                    "p_literature": 0.8,
                    "raw_score_max": 0.74,
                    "calibrated_prob_max": 0.8,
                    "score_components_json": "{\"p_literature_formula\":\"1-prod(1-p_sentence_i)\"}",
                    "license_id": "local_articles:test",
                    "source_release": release,
                    "parser_hash": "lit_parser",
                    "config_hash": "lit_config",
                }
            ],
        )
        write_parquet(
            literature / "relation_candidates.parquet",
            [
                {
                    "relation_uid": "litrel_1",
                    "subject_uid": "target_1",
                    "subject_type": "target",
                    "predicate": "target_associated_with_disease",
                    "object_uid": "disease_1",
                    "object_type": "disease",
                    "polarity": "association",
                    "cancer_context": "Cancer",
                    "sentence_uid": "sent_lit_1",
                    "pmid": "123",
                    "pmcid": "PMC123",
                    "section": "abstract",
                    "trigger_phrase": "associated with",
                    "extraction_rule_id": "association.explicit",
                    "raw_score": 0.7,
                    "calibrated_prob": 0.4,
                    "score_components_json": "{\"base_rule_confidence\":0.7}",
                    "license_id": "local_articles:test",
                    "source_release": release,
                    "parser_hash": "lit_parser",
                    "config_hash": "lit_config",
                }
            ],
        )
        write_parquet(
            literature / "sentence_mentions.parquet",
            [
                {
                    "mention_uid": "mention_1",
                    "sentence_uid": "sent_lit_1",
                    "article_uid": "article_lit_1",
                    "pmid": "123",
                    "pmcid": "PMC123",
                    "section": "abstract",
                    "entity_uid": "target_1",
                    "entity_type": "target",
                    "display_name": "HK1",
                    "surface": "HK1",
                    "normalized_surface": "hk1",
                    "start_offset": 0,
                    "end_offset": 3,
                    "matched_field": "approved_symbol",
                    "resolution_confidence": 0.94,
                    "source_release": release,
                    "license_id": "local_articles:test",
                    "parser_hash": "lit_parser",
                    "config_hash": "lit_config",
                },
                {
                    "mention_uid": "mention_2",
                    "sentence_uid": "sent_lit_1",
                    "article_uid": "article_lit_1",
                    "pmid": "123",
                    "pmcid": "PMC123",
                    "section": "abstract",
                    "entity_uid": "disease_1",
                    "entity_type": "disease",
                    "display_name": "Cancer",
                    "surface": "Cancer",
                    "normalized_surface": "cancer",
                    "start_offset": 23,
                    "end_offset": 29,
                    "matched_field": "name",
                    "resolution_confidence": 0.91,
                    "source_release": release,
                    "license_id": "local_articles:test",
                    "parser_hash": "lit_parser",
                    "config_hash": "lit_config",
                },
            ],
        )
        (literature / "literature_evidence_manifest.json").write_text(
            json.dumps(
                {
                    "release_id": release,
                    "manifest_hash": "lit_hash",
                    "builder_version": "test",
                    "metrics": {"evidence_candidate_count": 1, "supported_existing_edge_count": 1},
                    "tables": [
                        {"table": "sentence_mentions", "rows": 2},
                        {"table": "relation_candidates", "rows": 1},
                        {"table": "literature_edge_support", "rows": 1},
                    ],
                }
            ),
            encoding="utf-8",
        )

        return metabo_service.MetaboService(root, release_id=release)

    def test_resolve_precheck_and_pubchem_contracts(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))

            resolved = service.resolve("HMDB:HMDB0000122", entity_type="metabolite")
            self.assertEqual(resolved["status"], "matched")
            self.assertEqual(resolved["candidates"][0]["entity_uid"], "met_1")

            ambiguous = service.resolve("Shared", entity_type="metabolite")
            self.assertEqual(ambiguous["status"], "ambiguous")

            precheck = service.precheck_metabolites(
                [
                    {"HMDB": "HMDB0000122"},
                    {"name": "Shared"},
                    {"name": "No such metabolite"},
                    {"mz": 180.063},
                ]
            )
            self.assertEqual(precheck["summary"], {"matched": 1, "ambiguous": 1, "unmatched": 1, "invalid": 1})

            pubchem = service.pubchem("PubChem:5793")
            self.assertEqual(pubchem["status"], "found")
            self.assertEqual(pubchem["properties"]["formula"], "C6H12O6")

    def test_ambiguous_rows_create_low_weight_expanded_candidates(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))

            records = [{"name": "Shared", "log2FC": 1.0, "padj": 0.01, "direction": "up"}]
            precheck = service.precheck_metabolites(records)

            self.assertEqual(precheck["summary"]["matched"], 0)
            self.assertEqual(precheck["summary"]["ambiguous"], 1)
            self.assertEqual(precheck["expanded_summary"]["expanded_candidate_count"], 2)
            self.assertEqual(precheck["expanded_summary"]["expanded_by_class"], {"soft_identity": 2})
            self.assertEqual(precheck["expanded_summary"]["analysis_seed_input_count"], 1)
            self.assertAlmostEqual(precheck["expanded_candidates"][0]["expanded_candidate"]["weight_multiplier"], 0.25)

            analyzed = service.analyze_metabolites(records, max_paths=5, max_hops=2)
            pack = analyzed["analysis_pack"]
            self.assertEqual(pack["input_summary"]["matched_count"], 0)
            self.assertEqual(pack["input_summary"]["expanded_candidate_count"], 2)
            self.assertEqual(len(pack["matched"]), 0)
            self.assertEqual(len(pack["expanded_candidates"]), 2)
            self.assertTrue(all(row["seed_track"] == "expanded" for row in pack["expanded_candidates"]))
            self.assertTrue(all(row["seed_class"] == "soft_identity" for row in pack["expanded_candidates"]))
            self.assertTrue(all(feature["seed_track"] == "expanded" for feature in analyzed["analysis_features"]["matched"]))
            self.assertTrue(all(0.0 < weight < 1.0 for weight in analyzed["analysis_features"]["seed_weights"].values()))

    def test_european_ratio_traits_become_low_weight_ratio_component_seeds(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            annotation_dir = root / "raw_lake" / "European"
            annotation_dir.mkdir(parents=True, exist_ok=True)
            (annotation_dir / "European_trait_annotations.csv").write_text(
                "\n".join(
                    [
                        "accession_id,reported_trait,summary_statistics_url,pubmed_id,paper_title,resolution_status,mapped_names,local_name_hits,pubchem_cids,local_pubchem_cid_hits,pubchem_titles,pubchem_formulas,pubchem_inchikeys,candidate_names,pubchem_query,source_name",
                        "GCST90201000,Shared to Glucose ratio,http://example.org/GCST90201000,123456,European trait paper,local_name,Shared|Glucose,Shared|Glucose,,,,,,Shared|Glucose,,test",
                    ]
                ),
                encoding="utf-8",
            )
            service = self.make_service(root)

            precheck = service.precheck_metabolites([{"trait": "GCST90201000", "log2FC": 1.0, "padj": 0.01}])

            self.assertEqual(precheck["summary"]["matched"], 0)
            self.assertEqual(precheck["summary"]["ambiguous"], 1)
            self.assertEqual(precheck["expanded_summary"]["expanded_by_class"], {"ratio_component": 2})
            self.assertEqual(precheck["expanded_summary"]["expanded_candidate_count"], 2)
            self.assertAlmostEqual(precheck["expanded_summary"]["weighted_seed_mass"], 0.25)
            self.assertEqual(precheck["ambiguous"][0]["resolution"]["identity_review_reasons"], ["ratio_or_composite_trait_not_strict_identity"])

            analyzed = service.analyze_metabolites([{"trait": "GCST90201000", "log2FC": 1.0, "padj": 0.01}], max_paths=5, max_hops=2)
            evidence = analyzed["analysis_pack"]["ambiguous"][0]["source_trait_evidence"]
            self.assertEqual(evidence["summary_statistics_url"], "http://example.org/GCST90201000")
            self.assertEqual(evidence["summary_statistics_url_type"], "gwas_summary_statistics")
            self.assertEqual(evidence["pubmed_id"], "123456")
            exposures = analyzed["analysis_pack"]["genetic_exposures"]
            self.assertEqual(len(exposures), 1)
            self.assertEqual(exposures[0]["exposure_id"], "GCST90201000")
            self.assertEqual(exposures[0]["genetic_exposure_status"], "available")
            self.assertEqual(exposures[0]["chemical_identity_track"], "expanded_mechanism_seed")
            self.assertEqual(exposures[0]["chemical_seed_classes"], {"ratio_component": 2})
            self.assertEqual(analyzed["analysis_pack"]["input_summary"]["genetic_exposure_count"], 1)
            self.assertEqual(analyzed["analysis_pack"]["four_track_summary"]["genetic_exposure_count"], 1)
            self.assertIn("genetic_overlap", exposures[0]["allowed_uses"])

    def test_precheck_and_analyze_include_database_accuracy_v2_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            service = self.make_service(root)
            annotation_dir = root / "raw_lake" / "European"
            annotation_dir.mkdir(parents=True, exist_ok=True)
            (annotation_dir / "European_trait_annotations.csv").write_text(
                "\n".join(
                    [
                        "accession_id,reported_trait,summary_statistics_url,pubmed_id,paper_title,resolution_status,mapped_names,candidate_names",
                        "GCST90201000,Tryptophan to Pyruvate ratio,http://example.org/GCST90201000,123,Paper,local_name,Tryptophan|Pyruvate,Tryptophan|Pyruvate",
                    ]
                ),
                encoding="utf-8",
            )
            database_accuracy_v2.build_database_accuracy_store(root, service.release_id)
            service = metabo_service.MetaboService(root, release_id=service.release_id)

            precheck = service.precheck_metabolites(
                [
                    {"HMDB": "HMDB0000122", "log2FC": 1.0, "padj": 0.01},
                    {"trait": "GCST90201000", "reported_trait": "Tryptophan to Pyruvate ratio"},
                    {"name": "Sphingomyelin"},
                ]
            )

            self.assertTrue(precheck["identity_resolution_v2_summary"]["store_available"])
            self.assertEqual(precheck["matched"][0]["identity_resolution_v2"]["decision"]["decision_status"], "accepted_exact")
            all_rows = [*precheck["matched"], *precheck["ambiguous"], *precheck["unmatched"], *precheck["invalid"]]
            ratio_row = next(row for row in all_rows if row["identity_resolution_v2"]["input_feature"]["feature_type"] == "ratio")
            self.assertEqual(ratio_row["identity_resolution_v2"]["decision"]["decision_status"], "accepted_trait")
            self.assertNotEqual(ratio_row["identity_resolution_v2"]["decision"]["accepted_entity_type"], "chemical")
            class_row = next(row for row in all_rows if row["identity_resolution_v2"]["input_feature"]["feature_label"] == "Sphingomyelin")
            self.assertEqual(class_row["identity_resolution_v2"]["decision"]["decision_status"], "accepted_class")

            analyzed = service.analyze_metabolites([{"HMDB": "HMDB0000122", "log2FC": 1.0, "padj": 0.01}], max_paths=5, max_hops=2)
            accuracy = analyzed["analysis_pack"]["database_accuracy"]
            self.assertEqual(accuracy["contract_version"], "database_accuracy_store.v2")
            self.assertTrue(accuracy["store_available"])
            self.assertEqual(accuracy["identity_decision_summary"]["accepted_exact"], 1)
            self.assertTrue(accuracy["mechanism_ready_facts"])
            self.assertTrue(analyzed["interpretation_report"]["database_accuracy"]["store_available"])

    def test_database_accuracy_store_missing_falls_back_with_warning(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            analyzed = service.analyze_metabolites([{"HMDB": "HMDB0000122"}], max_paths=5, max_hops=2)
            accuracy = analyzed["analysis_pack"]["database_accuracy"]
            self.assertFalse(accuracy["store_available"])
            self.assertIn(
                "database_accuracy_store_missing",
                {warning["code"] for warning in analyzed["analysis_pack"]["quality_warnings"]},
            )
            self.assertEqual(analyzed["analysis_pack"]["matched"][0]["identity_resolution_v2"]["decision"]["decision_status"], "accepted_exact")

    def test_pathway_context_rerank_normalizes_legacy_query_notes(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            components = {
                "identifier": 0.0,
                "structure": 0.0,
                "name": 80.0,
                "mass": 0.0,
                "formula": 0.0,
                "rt": 0.0,
                "ms2": 0.0,
                "context": 0.0,
            }
            row_contracts = [
                {
                    "resolution": {
                        "resolver_version": "compound_resolver_v2",
                        "status": "matched",
                        "top_margin": 20.0,
                        "query_notes": ["legacy_name_note"],
                        "candidates": [
                            {
                                "entity_uid": "met_1",
                                "display_name": "Glucose",
                                "score": 80.0,
                                "score_components": dict(components),
                                "matches": [],
                            }
                        ],
                    }
                },
                {
                    "resolution": {
                        "resolver_version": "compound_resolver_v2",
                        "status": "matched",
                        "top_margin": 20.0,
                        "query_notes": [],
                        "candidates": [
                            {
                                "entity_uid": "met_2",
                                "display_name": "Shared",
                                "score": 80.0,
                                "score_components": dict(components),
                                "matches": [],
                            }
                        ],
                    }
                },
            ]

            service.apply_pathway_context_rerank(row_contracts)

            self.assertIsInstance(row_contracts[0]["resolution"]["query_notes"][0], dict)
            self.assertEqual(row_contracts[0]["resolution"]["query_notes"][0]["raw_note"], "legacy_name_note")
            self.assertTrue(any(note.get("mode") == "pathway_context" for note in row_contracts[0]["resolution"]["query_notes"]))

    def test_conclusion_evaluation_matches_gold_and_reports_negative_traps(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_dir = root / "config"
            config_dir.mkdir(parents=True, exist_ok=True)
            (config_dir / "gold_standard_conclusions.json").write_text(
                json.dumps(
                    {
                        "conclusions": [
                            {
                                "gold_id": "gold_glycolysis",
                                "mechanism_axis": "glycolysis",
                                "expected_entities": ["glucose"],
                                "polarity": "positive",
                                "source_refs": ["internal:test"],
                            },
                            {
                                "gold_id": "gold_unrelated_disease_trap",
                                "mechanism_axis": "cancer",
                                "expected_entities": ["cancer"],
                                "negative_trap_terms": ["cancer"],
                                "polarity": "negative",
                                "is_negative_control": True,
                                "source_refs": ["internal:test"],
                            },
                        ]
                    }
                ),
                encoding="utf-8",
            )
            service = self.make_service(root)

            analyzed = service.analyze_metabolites([{"HMDB": "HMDB0000122", "log2FC": 1.0, "padj": 0.01}], max_paths=5, max_hops=2)
            evaluation = analyzed["analysis_pack"]["conclusion_evaluation"]

            self.assertEqual(evaluation["contract_version"], "conclusion_evaluation.v1")
            self.assertEqual(evaluation["gold_standard"]["entry_count"], 2)
            by_id = {row["gold_id"]: row for row in evaluation["per_gold"]}
            self.assertIn(by_id["gold_glycolysis"]["gold_match_status"], {"partial_match", "exact_match"})
            self.assertEqual(by_id["gold_unrelated_disease_trap"]["gold_match_status"], "negative_control_not_triggered")
            self.assertEqual(evaluation["metrics"]["negative_trap_hits"], 0)
            self.assertTrue(evaluation["release_gate"]["passed"])
            self.assertEqual(
                analyzed["interpretation_report"]["conclusion_evaluation"]["contract_version"],
                "conclusion_evaluation.v1",
            )

            fake_pack = {
                "pathway_rankings": [
                    {
                        "pathway_uid": "pathway_bad",
                        "label": "Cancer",
                        "evidence_refs": [{"ref_type": "edge", "edge_uid": "edge_bad"}],
                        "appendix": False,
                    }
                ],
                "database_accuracy": {},
            }
            trap_eval = service.build_conclusion_evaluation(fake_pack, context={})
            trap_by_id = {row["gold_id"]: row for row in trap_eval["per_gold"]}
            self.assertEqual(trap_by_id["gold_unrelated_disease_trap"]["gold_match_status"], "false_positive_trap")
            self.assertFalse(trap_eval["release_gate"]["passed"])

            context_mismatch_pack = {
                "pathway_rankings": [
                    {
                        "pathway_uid": "pathway_context_mismatch",
                        "label": "Cancer",
                        "evidence_refs": [{"ref_type": "edge", "edge_uid": "edge_context_mismatch"}],
                        "appendix": False,
                        "context_mismatch": True,
                    }
                ],
                "database_accuracy": {},
            }
            context_eval = service.build_conclusion_evaluation(context_mismatch_pack, context={})
            context_by_id = {row["gold_id"]: row for row in context_eval["per_gold"]}
            self.assertEqual(
                context_by_id["gold_unrelated_disease_trap"]["gold_match_status"],
                "negative_control_appendix_only",
            )
            self.assertTrue(context_eval["release_gate"]["passed"])

            observation_pack = {
                "database_accuracy": {
                    "trait_observations": [
                        {
                            "feature_label": "Cancer",
                            "feature_type": "trait_score",
                            "decision_status": "accepted_trait",
                        }
                    ]
                }
            }
            observation_eval = service.build_conclusion_evaluation(observation_pack, context={})
            observation_by_id = {row["gold_id"]: row for row in observation_eval["per_gold"]}
            self.assertEqual(
                observation_by_id["gold_unrelated_disease_trap"]["gold_match_status"],
                "negative_control_appendix_only",
            )
            self.assertTrue(observation_eval["release_gate"]["passed"])

            positive_observation_pack = {
                "database_accuracy": {
                    "trait_observations": [
                        {
                            "feature_label": "glycolysis glucose ratio",
                            "feature_type": "ratio",
                            "decision_status": "accepted_trait",
                        }
                    ]
                }
            }
            positive_observation_eval = service.build_conclusion_evaluation(positive_observation_pack, context={})
            positive_observation_by_id = {row["gold_id"]: row for row in positive_observation_eval["per_gold"]}
            self.assertEqual(positive_observation_by_id["gold_glycolysis"]["gold_match_status"], "appendix_only")
            self.assertEqual(positive_observation_eval["metrics"]["gold_positive_matched"], 0)
            self.assertEqual(positive_observation_eval["metrics"]["unsupported_top_claim_rate"], None)

    def test_exact_input_observations_can_support_gold_entity_matches(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_dir = root / "config"
            config_dir.mkdir(parents=True, exist_ok=True)
            (config_dir / "gold_standard_conclusions.json").write_text(
                json.dumps(
                    {
                        "conclusions": [
                            {
                                "gold_id": "gold_serine",
                                "mechanism_axis": "serine one-carbon metabolism",
                                "expected_entities": ["serine", "glycine"],
                                "polarity": "positive",
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            service = self.make_service(root)
            pack = {
                "database_accuracy": {
                    "exact_input_observations": [
                        {
                            "feature_label": "Serine",
                            "decision_status": "accepted_exact",
                            "allowed_claim_scope": "input_exact_observation",
                            "input_trace": True,
                        }
                    ]
                }
            }

            evaluation = service.build_conclusion_evaluation(pack, context={})
            by_id = {row["gold_id"]: row for row in evaluation["per_gold"]}

            self.assertEqual(by_id["gold_serine"]["gold_match_status"], "partial_match")
            self.assertEqual(evaluation["metrics"]["unsupported_top_claim_rate"], 0.0)

    def test_gold_entity_matching_uses_common_metabolite_aliases(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_dir = root / "config"
            config_dir.mkdir(parents=True, exist_ok=True)
            (config_dir / "gold_standard_conclusions.json").write_text(
                json.dumps(
                    {
                        "conclusions": [
                            {
                                "gold_id": "gold_glucose",
                                "mechanism_axis": "glycolysis",
                                "expected_entities": ["glucose"],
                                "polarity": "positive",
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            service = self.make_service(root)
            pack = {
                "database_accuracy": {
                    "exact_input_observations": [
                        {
                            "feature_label": "D-glucopyranose",
                            "decision_status": "accepted_exact",
                            "allowed_claim_scope": "input_exact_observation",
                            "input_trace": True,
                        }
                    ]
                }
            }

            evaluation = service.build_conclusion_evaluation(pack, context={})
            by_id = {row["gold_id"]: row for row in evaluation["per_gold"]}

            self.assertEqual(by_id["gold_glucose"]["gold_match_status"], "partial_match")
            self.assertEqual(evaluation["metrics"]["precision_at_1"], 1.0)

    def test_generic_gold_terms_do_not_trigger_fatty_acid_false_matches(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_dir = root / "config"
            config_dir.mkdir(parents=True, exist_ok=True)
            (config_dir / "gold_standard_conclusions.json").write_text(
                json.dumps(
                    {
                        "conclusions": [
                            {
                                "gold_id": "gold_fatty_acid",
                                "mechanism_axis": "fatty acid metabolism",
                                "expected_entities": ["fatty acid"],
                                "polarity": "positive",
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            service = self.make_service(root)
            pack = {
                "pathway_rankings": [
                    {
                        "pathway_uid": "amino_acid_pathway",
                        "label": "Amino acid metabolism",
                        "evidence_refs": [{"ref_type": "edge", "edge_uid": "edge_amino_acid"}],
                    }
                ],
                "database_accuracy": {},
            }

            evaluation = service.build_conclusion_evaluation(pack, context={})
            by_id = {row["gold_id"]: row for row in evaluation["per_gold"]}

            self.assertEqual(by_id["gold_fatty_acid"]["gold_match_status"], "unsupported")
            self.assertEqual(by_id["gold_fatty_acid"]["top_matches"], [])

    def test_class_observations_can_support_class_level_gold_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_dir = root / "config"
            config_dir.mkdir(parents=True, exist_ok=True)
            (config_dir / "gold_standard_conclusions.json").write_text(
                json.dumps(
                    {
                        "conclusions": [
                            {
                                "gold_id": "gold_lipid_class",
                                "mechanism_axis": "lipid class remodeling",
                                "expected_entities": ["fatty acid"],
                                "required_evidence_type": "class_or_reaction_evidence",
                                "polarity": "positive",
                            },
                            {
                                "gold_id": "gold_serine_exact",
                                "mechanism_axis": "serine one-carbon metabolism",
                                "expected_entities": ["serine"],
                                "required_evidence_type": "direct_reaction_evidence",
                                "polarity": "positive",
                            },
                        ]
                    }
                ),
                encoding="utf-8",
            )
            service = self.make_service(root)
            pack = {
                "database_accuracy": {
                    "class_observations": [
                        {
                            "feature_label": "Fatty acid",
                            "decision_status": "accepted_class",
                            "allowed_claim_scope": "input_class_observation",
                            "input_trace": True,
                        },
                        {
                            "feature_label": "Serine",
                            "decision_status": "accepted_class",
                            "allowed_claim_scope": "input_class_observation",
                            "input_trace": True,
                        },
                    ]
                }
            }

            evaluation = service.build_conclusion_evaluation(pack, context={})
            by_id = {row["gold_id"]: row for row in evaluation["per_gold"]}

            self.assertEqual(by_id["gold_lipid_class"]["gold_match_status"], "partial_match")
            self.assertEqual(by_id["gold_serine_exact"]["gold_match_status"], "appendix_only")
            self.assertEqual(evaluation["metrics"]["gold_positive_matched"], 1)
            self.assertEqual(evaluation["metrics"]["precision_at_1"], 1.0)

    def test_exact_input_observation_label_uses_accepted_candidate_name_for_identifier_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            service = self.make_service(root)
            database_accuracy_v2.build_database_accuracy_store(root, service.release_id)

            analyzed = service.analyze_metabolites(
                [{"HMDB": "HMDB0000122", "log2FC": 1.0, "padj": 0.01}],
                max_paths=5,
                max_hops=2,
            )
            observations = analyzed["analysis_pack"]["database_accuracy"]["exact_input_observations"]

            self.assertTrue(observations)
            self.assertNotEqual(observations[0]["feature_label"], "row_0")
            self.assertTrue(observations[0]["accepted_entity_name"])
            self.assertEqual(observations[0]["feature_label"], observations[0]["accepted_entity_name"])

    def test_mechanism_ready_facts_apply_context_mismatch_appendix(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            service = self.make_service(root)
            database_accuracy_v2.build_database_accuracy_store(root, service.release_id)
            facts_path = root / "database_accuracy_store" / service.release_id / "analysis_view" / "mechanism_ready_facts.parquet"
            table = pq.read_table(facts_path)
            rows = table.to_pylist()
            rows[0]["object_name"] = "Sweet compound binds sweet taste receptor"
            rows[0]["metadata_json"] = json.dumps({"reaction_name": "Sweet compound binds sweet taste receptor"})
            pq.write_table(pa.Table.from_pylist(rows, schema=table.schema), facts_path)
            service = metabo_service.MetaboService(root, release_id=service.release_id)

            analyzed = service.analyze_metabolites(
                [{"HMDB": "HMDB0000122", "log2FC": 1.0, "padj": 0.01}],
                max_paths=5,
                max_hops=2,
                context={"cell_type": "epithelial", "context_terms": ["tumor", "epithelial"]},
            )
            facts = analyzed["analysis_pack"]["database_accuracy"]["mechanism_ready_facts"]
            taste_like = [
                fact
                for fact in facts
                if "taste" in str(fact.get("object_name") or "").casefold()
                or "taste" in str((fact.get("metadata") or {}).get("reaction_name") or "").casefold()
            ]
            self.assertTrue(taste_like)
            self.assertTrue(all(fact.get("appendix") for fact in taste_like))
            self.assertTrue(all(fact.get("context_mismatch") for fact in taste_like))

    def test_conclusion_evaluation_counts_only_context_applicable_positive_gold(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_dir = root / "config"
            config_dir.mkdir(parents=True, exist_ok=True)
            (config_dir / "gold_standard_conclusions.json").write_text(
                json.dumps(
                    {
                        "conclusions": [
                            {
                                "gold_id": "gold_brca",
                                "cancer_type": ["breast cancer", "brca"],
                                "mechanism_axis": "glycolysis",
                                "expected_entities": ["glucose"],
                                "polarity": "positive",
                            },
                            {
                                "gold_id": "gold_hcc",
                                "cancer_type": ["hepatocellular carcinoma", "hcc"],
                                "mechanism_axis": "bile acid",
                                "expected_entities": ["bile acid"],
                                "polarity": "positive",
                            },
                        ]
                    }
                ),
                encoding="utf-8",
            )
            service = self.make_service(root)
            pack = {
                "pathway_rankings": [
                    {
                        "pathway_uid": "glycolysis",
                        "label": "glycolysis glucose",
                        "evidence_refs": [{"ref_type": "edge", "edge_uid": "edge_glucose"}],
                    }
                ],
                "database_accuracy": {},
            }

            evaluation = service.build_conclusion_evaluation(pack, context={"terms": ["BRCA", "breast cancer"]})
            by_id = {row["gold_id"]: row for row in evaluation["per_gold"]}

            self.assertEqual(evaluation["metrics"]["gold_positive_count"], 1)
            self.assertEqual(evaluation["metrics"]["gold_positive_matched"], 1)
            self.assertEqual(evaluation["metrics"]["precision_at_1"], 1.0)
            self.assertEqual(evaluation["metrics"]["recall_at_1"], 1.0)
            self.assertEqual(by_id["gold_brca"]["gold_match_status"], "exact_match")
            self.assertEqual(by_id["gold_hcc"]["gold_match_status"], "not_applicable_context")

    def test_conclusion_evaluation_reports_precision_and_recall_at_k(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_dir = root / "config"
            config_dir.mkdir(parents=True, exist_ok=True)
            (config_dir / "gold_standard_conclusions.json").write_text(
                json.dumps(
                    {
                        "conclusions": [
                            {
                                "gold_id": "gold_glycolysis",
                                "mechanism_axis": "glycolysis",
                                "expected_entities": ["glucose"],
                                "polarity": "positive",
                            },
                            {
                                "gold_id": "gold_bile_acid",
                                "mechanism_axis": "bile acid",
                                "expected_entities": ["cholesterol"],
                                "polarity": "positive",
                            },
                            {
                                "gold_id": "gold_serine",
                                "mechanism_axis": "serine one-carbon",
                                "expected_entities": ["glycine"],
                                "polarity": "positive",
                            },
                        ]
                    }
                ),
                encoding="utf-8",
            )
            service = self.make_service(root)
            pack = {
                "pathway_rankings": [
                    {
                        "pathway_uid": "glycolysis",
                        "label": "glycolysis glucose",
                        "evidence_refs": [{"ref_type": "edge", "edge_uid": "edge_glucose"}],
                    },
                    {
                        "pathway_uid": "unrelated",
                        "label": "taste perception",
                        "evidence_refs": [{"ref_type": "edge", "edge_uid": "edge_unrelated"}],
                    },
                    {
                        "pathway_uid": "bile_acid",
                        "label": "bile acid cholesterol",
                        "evidence_refs": [{"ref_type": "edge", "edge_uid": "edge_bile"}],
                    },
                ],
                "database_accuracy": {},
            }

            evaluation = service.build_conclusion_evaluation(pack, context={})
            metrics = evaluation["metrics"]

            self.assertEqual(metrics["ranked_candidate_count"], 3)
            self.assertEqual(metrics["precision_at_1"], 1.0)
            self.assertEqual(metrics["recall_at_1"], 0.333333)
            self.assertEqual(metrics["precision_at_3"], 0.666667)
            self.assertEqual(metrics["recall_at_3"], 0.666667)

    def test_conclusion_evaluation_topk_prioritizes_mechanism_facts_over_auxiliary_rankings(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_dir = root / "config"
            config_dir.mkdir(parents=True, exist_ok=True)
            (config_dir / "gold_standard_conclusions.json").write_text(
                json.dumps(
                    {
                        "conclusions": [
                            {
                                "gold_id": "gold_glycolysis",
                                "mechanism_axis": "glycolysis",
                                "expected_entities": ["glucose"],
                                "polarity": "positive",
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            service = self.make_service(root)
            pack = {
                "pathway_rankings": [
                    {
                        "pathway_uid": "unrelated_pathway",
                        "label": "taste perception",
                        "evidence_refs": [{"ref_type": "edge", "edge_uid": "edge_taste"}],
                    }
                ],
                "database_accuracy": {
                    "mechanism_ready_facts": [
                        {
                            "fact_uid": "fact_glucose",
                            "subject_name": "glucose",
                            "predicate": "participates_in",
                            "object_name": "glycolysis",
                            "supporting_assertion_uids": ["assertion_glucose"],
                        }
                    ]
                },
            }

            evaluation = service.build_conclusion_evaluation(pack, context={})

            self.assertEqual(evaluation["top_ranked_candidates"][0]["source"], "mechanism_ready_fact")
            self.assertEqual(evaluation["metrics"]["precision_at_1"], 1.0)
            self.assertEqual(evaluation["metrics"]["recall_at_1"], 1.0)

    def test_expanded_seed_coverage_is_exploratory_usable_not_input_failed(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            input_summary = {
                "input_count": 100,
                "matched_count": 30,
                "ambiguous_count": 50,
                "unmatched_count": 20,
                "invalid_count": 0,
                "expanded_candidate_input_count": 50,
                "analysis_seed_input_count": 80,
            }
            self.assertEqual(service.match_confidence_weight(input_summary), 0.55)

            precheck = {
                "input_count": 100,
                "summary": {"matched": 30, "ambiguous": 50, "unmatched": 20, "invalid": 0},
                "expanded_summary": {"expanded_candidate_input_count": 50, "analysis_seed_input_count": 80},
            }
            assessment = service.prediction_model_assessment(
                precheck,
                [{"confidence_tier": "medium", "result_type": "metabolic_theme"}],
                [],
            )
            self.assertEqual(assessment["grade"], "C")
            self.assertEqual(assessment["level"], "exploratory_usable")
            self.assertGreaterEqual(assessment["analysis_seed_input_ratio"], 0.8)

    def test_analyze_metabolites_returns_pathways_and_explanations(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            records = [{"HMDB": "HMDB0000122", "log2FC": 2.0, "padj": 0.001, "direction": "up"}]
            analyzed = service.analyze_metabolites(records, max_paths=10, max_hops=5)

            self.assertEqual(analyzed["precheck"]["summary"]["matched"], 1)
            self.assertEqual(analyzed["analysis_features"]["summary"]["directions"]["up"], 1)
            self.assertGreater(analyzed["analysis_features"]["matched"][0]["p_user"], 0.0)
            pathways = analyzed["rankings"]["pathways"]
            self.assertEqual(pathways[0]["pathway_uid"], "pathway_1")
            self.assertEqual(pathways[0]["overlap_count"], 1)
            self.assertEqual(analyzed["directional_enrichment"]["summary"]["up"], 1)
            self.assertEqual(analyzed["directional_enrichment"]["up"][0]["pathway_uid"], "pathway_1")
            self.assertEqual(analyzed["propagation"]["formula"], "pi=(1-alpha)y+alpha*W*pi")
            self.assertIn("compression", analyzed["propagation"]["graph"])
            self.assertEqual(analyzed["propagation"]["graph"]["compression"]["max_adaptive_edges_per_node"], 160)
            self.assertEqual(analyzed["compressed_explanation_paths"]["mode"], "path_signature_compression")
            self.assertEqual(analyzed["rankings"]["targets"][0]["target_uid"], "target_1")
            self.assertEqual(analyzed["rankings"]["diseases"][0]["disease_uid"], "disease_1")
            self.assertIn("propagation_score", analyzed["rankings"]["targets"][0]["score_components"])
            overlay_edges = [edge for edge in service.incident_edges("met_1") if edge.get("evidence_level") == "literature_overlay"]
            self.assertFalse(overlay_edges)
            self.assertTrue(
                any(ref.get("support_uid") == "litsup_overlay" for ref in analyzed["rankings"]["targets"][0]["evidence_refs"])
            )
            for ranking_name in ("pathways", "targets", "diseases"):
                row = analyzed["rankings"][ranking_name][0]
                self.assertIn("literature_support_count", row)
                self.assertIn("max_p_literature", row)
                self.assertIn("supported_pmids", row)
                self.assertIn("support_classes", row)
                self.assertIn("evidence_refs", row)
                self.assertIn("literature_support_count", row["score_components"])
                self.assertIn("max_p_literature", row["score_components"])
            self.assertTrue(analyzed["explanation_paths"])
            self.assertTrue(any(path["nodes"][-1]["node_type"] == "target" for path in analyzed["explanation_paths"]))
            pack = analyzed["analysis_pack"]
            self.assertEqual(pack["contract_version"], "analysis_pack.v1")
            self.assertEqual(pack["input_summary"]["input_count"], 1)
            self.assertEqual(pack["input_summary"]["matched_count"], 1)
            self.assertEqual(pack["matched"][0]["metabolite_uid"], "met_1")
            self.assertFalse(pack["blocked_reasons"])
            self.assertTrue(pack["pathway_rankings"])
            self.assertTrue(pack["target_rankings"])
            self.assertTrue(pack["disease_rankings"])
            self.assertEqual(
                pack["top_explanation_paths"][0]["path_id"],
                pack["pathway_rankings"][0]["best_path_id"],
            )
            self.assertEqual(pack["prediction_model"]["contract_version"], "metabolic_prediction_model.v1")
            self.assertEqual(pack["prediction_model"]["assessment"]["grade"], "A")
            self.assertEqual(pack["pathway_rankings"][0]["prediction_task"], "pathway_prediction")
            self.assertEqual(pack["pathway_rankings"][0]["confidence_tier"], "medium")
            self.assertEqual(pack["pathway_rankings"][0]["calibration_status"], "calibrated_in_scope")
            self.assertFalse(pack["pathway_rankings"][0]["is_extrapolated"])
            self.assertEqual(pack["target_rankings"][0]["prediction_task"], "target_prediction")
            self.assertTrue(pack["target_rankings"][0]["is_extrapolated"])
            self.assertTrue(pack["target_rankings"][0]["needs_validation"])
            for ranking_name in ("pathway_rankings", "target_rankings", "disease_rankings"):
                self.assertTrue(pack[ranking_name][0]["claim_refs"]["traceability_passed"])
                self.assertTrue(pack[ranking_name][0]["claim_refs"]["edge_uids"] or pack[ranking_name][0]["claim_refs"]["evidence_ref_uids"])
            self.assertGreaterEqual(pack["literature_evidence_pack"]["support_count"], 1)

            repeated = service.analyze_metabolites(records, max_paths=10, max_hops=5)
            self.assertEqual(
                analyzed["determinism"]["response_hash"],
                repeated["determinism"]["response_hash"],
            )
            self.assertEqual(
                analyzed["analysis_pack"]["determinism"]["analysis_pack_hash"],
                repeated["analysis_pack"]["determinism"]["analysis_pack_hash"],
            )

            subgraph = service.subgraph(["met_1"], max_hops=2)
            self.assertTrue(any(node["node_uid"] == "pathway_1" for node in subgraph["nodes"]))
            self.assertTrue(any(edge["edge_uid"] == "edge_mp" for edge in subgraph["edges"]))

            evidence_subgraph = service.subgraph(["target_1"], max_hops=1)
            disease_edge = next(edge for edge in evidence_subgraph["edges"] if edge["edge_uid"] == "edge_td")
            self.assertEqual(disease_edge["score_components"]["p_literature"], 0.4)
            self.assertTrue(disease_edge["literature_support"])

            evidence = service.evidence(edge_uid="edge_td")
            self.assertEqual(evidence["status"], "found")
            self.assertEqual(evidence["support"][0]["support_uid"], "litsup_1")
            self.assertEqual(evidence["relation_candidates"][0]["relation_uid"], "litrel_1")
            self.assertEqual(evidence["sentences"][0]["sentence_uid"], "sent_lit_1")
            self.assertEqual(len(evidence["mentions"]), 2)

            releases = service.releases()
            self.assertEqual(releases["releases"][0]["release_id"], "mvp_20260101T000000")
            self.assertTrue(releases["releases"][0]["has_literature_evidence"])

    def test_analysis_pack_adds_fallback_chains_when_graph_paths_are_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            records = [{"HMDB": "HMDB0000122", "log2FC": 2.0, "padj": 0.001, "direction": "up"}]
            precheck = service.precheck_metabolites(records)
            analysis_features, features_by_uid = service.matched_analysis_features(precheck["matched"])
            seed_weights = service.seed_weights_from_features(features_by_uid)
            seed_presence_weights = service.seed_presence_weights_from_features(features_by_uid)
            matched_uids = sorted(features_by_uid)
            propagation = service.propagate_scores(seed_weights)
            pathways = service.merge_pathway_rankings(
                service.pathway_enrichment(matched_uids, seed_presence_weights),
                propagation,
                {},
                features_by_uid,
                input_records=records,
            )
            targets = service.propagation_ranking(propagation, "target", "target_uid", {})
            diseases = service.propagation_ranking(propagation, "disease", "disease_uid", {})

            pack = service.build_analysis_pack(
                {"records": records},
                precheck,
                analysis_features,
                features_by_uid,
                seed_weights,
                pathways[:5],
                targets[:5],
                diseases[:5],
                [],
                propagation,
            )

            self.assertFalse(pack["top_explanation_paths"])
            self.assertTrue(pack["fallback_explanation_chains"])
            chain = pack["fallback_explanation_chains"][0]
            self.assertEqual(chain["chain_type"], "ranking_evidence_fallback")
            self.assertTrue(chain["is_fallback_chain"])
            self.assertLessEqual(chain["path_confidence"], 0.25)
            self.assertTrue(chain["claim_refs"]["traceability_passed"])
            self.assertIn("not a stable stepwise graph path", chain["boundary"])
            warning = next(row for row in pack["quality_warnings"] if row["code"] == "no_explanation_paths")
            self.assertIn("fallback_chain_policy", warning["details"])

    def test_context_fit_rerank_moves_low_fit_pathways_to_appendix(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            context = service.normalize_prediction_context(
                "ICC, intrahepatic cholangiocarcinoma, liver, epithelial, tumor vs adjacent"
            )
            rows = [
                {
                    "rank": 1,
                    "pathway_uid": "taste_pathway",
                    "display_name": "Sensory perception of sweet, bitter, and umami taste",
                    "score": 6.0,
                    "overlap_count": 2,
                    "matched_metabolite_uids": ["met_1"],
                    "score_components": {"input_theme_boost": 4.0, "weighted_overlap": 1.0},
                    "directional_support": {"input_theme_hits": []},
                },
                {
                    "rank": 2,
                    "pathway_uid": "glycolysis_pathway",
                    "display_name": "Glycolysis and pyruvate metabolism",
                    "score": 4.0,
                    "overlap_count": 2,
                    "matched_metabolite_uids": ["met_1"],
                    "score_components": {"weighted_overlap": 1.0, "literature_support_count": 2},
                    "literature_support_count": 2,
                    "directional_support": {"input_theme_hits": []},
                },
            ]
            features_by_uid = {"met_1": [{"seed_class": "strict_identity"}]}

            reranked = service.apply_context_fit_rerank(rows, features_by_uid, context)

            self.assertEqual(reranked[0]["pathway_uid"], "glycolysis_pathway")
            self.assertEqual(reranked[0]["original_rank"], 2)
            self.assertGreater(reranked[0]["context_fit_score"], reranked[1]["context_fit_score"])
            self.assertEqual(reranked[1]["pathway_uid"], "taste_pathway")
            self.assertEqual(reranked[1]["context_fit_tier"], "appendix_fit")
            self.assertTrue(reranked[1]["context_fit_appendix"])
            self.assertTrue(any("low_context_fit_terms" in item for item in reranked[1]["context_fit_penalties"]))

    def test_context_fit_rerank_prioritizes_matching_disease_context(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            context = service.normalize_prediction_context(
                "OV, ovarian cancer, ovary, epithelial, tumor vs adjacent"
            )
            rows = [
                {
                    "rank": 1,
                    "disease_uid": "disease_melanoma",
                    "display_name": "Melanoma",
                    "score": 8.0,
                    "calibrated_confidence": 0.4,
                    "input_support_count": 2,
                    "literature_support_count": 2,
                    "evidence_refs": [],
                },
                {
                    "rank": 2,
                    "disease_uid": "disease_ovarian",
                    "display_name": "Epithelial ovarian cancer",
                    "score": 5.0,
                    "calibrated_confidence": 0.35,
                    "input_support_count": 2,
                    "literature_support_count": 1,
                    "evidence_refs": [],
                },
                {
                    "rank": 3,
                    "disease_uid": "disease_prostate",
                    "display_name": "Prostate cancer",
                    "score": 4.0,
                    "calibrated_confidence": 0.3,
                    "input_support_count": 2,
                    "literature_support_count": 1,
                    "evidence_refs": [],
                },
            ]

            reranked = service.apply_context_fit_rerank(rows, {}, context, id_key="disease_uid")

            self.assertEqual(reranked[0]["disease_uid"], "disease_ovarian")
            self.assertIn("cancer_context_group_matches", reranked[0]["context_fit_reasons"])
            mismatch_rows = {row["disease_uid"]: row for row in reranked[1:]}
            self.assertTrue(mismatch_rows["disease_melanoma"]["context_fit_appendix"])
            self.assertTrue(any("different_cancer_context" in item for item in mismatch_rows["disease_prostate"]["context_fit_penalties"]))

    def test_biomedbert_context_aware_literature_refs_rank_matching_sentences_first(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            service = self.make_service(root)
            write_parquet(
                service.normalized_dir / "sentences.parquet",
                [
                    {
                        "sentence_uid": "sent_ov",
                        "article_uid": "article_1",
                        "pmid": "1",
                        "pmcid": "",
                        "section": "abstract",
                        "sentence_text": "Ovarian cancer epithelial tumor metabolism shows glycolytic remodeling.",
                        "text_hash": "hash_1",
                        "start_offset": 0,
                        "end_offset": 80,
                        "source_release": service.release_id,
                        "parser_hash": "parser",
                    },
                    {
                        "sentence_uid": "sent_pan",
                        "article_uid": "article_2",
                        "pmid": "2",
                        "pmcid": "",
                        "section": "abstract",
                        "sentence_text": "Cancer cells can alter metabolism in many tumor types.",
                        "text_hash": "hash_2",
                        "start_offset": 0,
                        "end_offset": 55,
                        "source_release": service.release_id,
                        "parser_hash": "parser",
                    },
                ],
            )
            write_parquet(
                service.literature_dir / "sentence_relevance.parquet",
                [
                    {"sentence_uid": "sent_ov", "relevance_score": 0.85, "decision": "evidence_candidate", "model_name": "biomedbert", "config_hash": "cfg"},
                    {"sentence_uid": "sent_pan", "relevance_score": 0.9, "decision": "evidence_candidate", "model_name": "biomedbert", "config_hash": "cfg"},
                ],
            )
            context = service.normalize_prediction_context("ovarian cancer, epithelial, tumor metabolism")
            refs = [
                {"ref_type": "literature_support", "support_uid": "support_pan", "p_literature": 0.9, "sentence_uids": ["sent_pan"]},
                {"ref_type": "literature_support", "support_uid": "support_ov", "p_literature": 0.7, "sentence_uids": ["sent_ov"]},
            ]

            reranked = service.apply_context_aware_literature_rerank(refs, context)

            self.assertEqual(reranked[0]["support_uid"], "support_ov")
            self.assertEqual(reranked[0]["context_relevance_tier"], "context_high")
            self.assertIn("ovarian cancer", reranked[0]["context_relevance_context_hits"])
            self.assertGreater(reranked[0]["context_relevance_score"], reranked[1]["context_relevance_score"])

    def test_structured_required_context_fields_are_normalized(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            context = service.normalize_prediction_context(
                {
                    "cancer_type": "ovarian cancer",
                    "tissue": "ovary",
                    "cell_type": "epithelial",
                    "comparison": "tumor vs adjacent",
                    "analysis_goal": "tumor metabolism",
                }
            )

            self.assertTrue(context["has_context"])
            self.assertEqual(set(context["fields"]["cancer_type"]), {"cancer", "ovarian cancer", "ovarian"})
            self.assertIn("ovary", context["fields"]["tissue"])
            self.assertIn("epithelial", context["fields"]["cell_type"])
            self.assertIn("tumor vs adjacent", context["fields"]["comparison"])
            self.assertIn("tumor metabolism", context["fields"]["analysis_goal"])
            self.assertIn("tumor metabolism", context["terms"])

    def test_head_neck_context_group_flags_wrong_cancer_diseases(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            context = service.normalize_prediction_context(
                {
                    "cancer_type": "head and neck squamous cell carcinoma",
                    "tissue": "head and neck",
                    "cell_type": "epithelial",
                    "comparison": "tumor vs adjacent",
                    "analysis_goal": "tumor metabolism",
                }
            )
            hnsc = service.ranking_context_fit_profile({"display_name": "head and neck squamous cell carcinoma"}, "disease_uid", context)
            oral = service.ranking_context_fit_profile({"display_name": "oral cavity squamous cell carcinoma"}, "disease_uid", context)
            esophageal = service.ranking_context_fit_profile({"display_name": "esophageal squamous cell carcinoma"}, "disease_uid", context)
            melanoma = service.ranking_context_fit_profile({"display_name": "melanoma"}, "disease_uid", context)

            self.assertIn("head_neck", hnsc["context_fit_display_groups"])
            self.assertIn("cancer_context_group_matches", hnsc["context_fit_reasons"])
            self.assertIn("head_neck", oral["context_fit_display_groups"])
            self.assertTrue(esophageal["context_fit_appendix"])
            self.assertNotIn("cancer_context_group_matches", esophageal["context_fit_reasons"])
            self.assertTrue(melanoma["context_fit_appendix"])
            self.assertTrue(any("different_cancer_context" in item for item in melanoma["context_fit_penalties"]))

    def test_context_groups_cover_traitscore_batch_cancers(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            cases = {
                "cervical cancer": "cervical",
                "colorectal cancer": "colorectal",
                "cutaneous squamous cell carcinoma": "skin_squamous",
                "esophageal cancer": "esophageal",
                "hepatocellular carcinoma": "liver",
                "head and neck squamous cell carcinoma": "head_neck",
                "intrahepatic cholangiocarcinoma": "cholangiocarcinoma",
                "chromophobe renal cell carcinoma": "renal",
                "clear cell renal cell carcinoma": "renal",
                "lung adenocarcinoma": "lung",
                "ovarian cancer": "ovarian",
                "gastric cancer": "gastric",
                "thyroid cancer": "thyroid",
            }
            for text, expected_group in cases.items():
                self.assertIn(expected_group, service.context_cancer_groups(text), text)

    def test_prediction_pack_uses_overlay_targets_and_context(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            service = self.make_service(root)
            overlay = root / "manual_sources" / "prediction_overlays" / "mvp_20260101T000000"
            overlay.mkdir(parents=True, exist_ok=True)
            (overlay / "drug_targets.csv").write_text(
                "\n".join(
                    [
                        "drug_id,drug_name,target_symbol,mechanism,confidence,context_terms,source_name,source_record_id,evidence_level,license_id",
                        "drug_hk1,HK1 test inhibitor,HK1,inhibitor,0.8,cancer,manual_seed,drug_hk1|HK1,predicted,local:test",
                        "drug_other,Other drug,NOPE,inhibitor,0.9,cancer,manual_seed,drug_other|NOPE,predicted,local:test",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            (overlay / "cell_type_signatures.csv").write_text(
                "\n".join(
                    [
                        "cell_type_id,cell_type_name,target_symbols,pathway_terms,metabolite_terms,context_terms,cell_state,confidence,source_name,source_record_id,evidence_level,license_id",
                        "cell_glycolytic,Glycolytic tumor-like cell,HK1,glycolysis,glucose,cancer,glycolytic,0.7,manual_seed,cell_glycolytic,predicted,local:test",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            records = [{"HMDB": "HMDB0000122", "log2FC": 2.0, "padj": 0.001, "direction": "up"}]
            analyzed = service.analyze_metabolites(records, max_paths=10, max_hops=5, context={"cancer_type": "Cancer", "cell_state": "glycolytic"})
            predictions = analyzed["predictions"]
            self.assertEqual(predictions["contract_version"], "prediction_pack.v1")
            self.assertEqual(predictions["context"]["fields"]["cancer_type"], ["cancer"])
            self.assertEqual(predictions["context"]["fields"]["cell_state"], ["glycolytic"])
            self.assertEqual(predictions["drug_rankings"][0]["drug_id"], "drug_hk1")
            self.assertEqual(predictions["drug_rankings"][0]["prediction_task"], "drug_hypothesis_prediction")
            self.assertEqual(predictions["drug_rankings"][0]["calibration_status"], "overlay_informed")
            self.assertTrue(predictions["drug_rankings"][0]["needs_validation"])
            self.assertEqual(predictions["drug_rankings"][0]["matched_targets"][0]["target_uid"], "target_1")
            self.assertEqual(predictions["cell_type_rankings"][0]["cell_type_id"], "cell_glycolytic")
            self.assertEqual(predictions["cell_type_rankings"][0]["prediction_task"], "cell_context_prediction")
            self.assertTrue(predictions["drug_rankings"][0]["evidence_refs"])

            hard_filtered = service.analyze_metabolites(
                records,
                max_paths=10,
                max_hops=5,
                context={"cancer_type": "melanoma"},
                context_mode="hard",
            )
            self.assertFalse(hard_filtered["predictions"]["drug_rankings"])
            blocked_codes = {row["code"] for row in hard_filtered["predictions"]["blocked_reasons"]}
            self.assertIn("no_drug_target_overlap", blocked_codes)

    def test_input_theme_boost_lifts_biochemical_pathways(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            features_by_uid = {
                "met_bile_1": [
                    {
                        "input_name": "taurocholic acid",
                        "display_name": "taurocholic acid",
                        "direction": "up",
                        "significant": True,
                        "p_user": 0.6,
                    }
                ],
                "met_bile_2": [
                    {
                        "input_name": "glycochenodeoxycholic acid",
                        "display_name": "glycochenodeoxycholic acid",
                        "direction": "up",
                        "significant": True,
                        "p_user": 0.5,
                    }
                ],
            }
            pathways = [
                {
                    "pathway_uid": "path_bile",
                    "name": "Bile acid and bile salt metabolism",
                    "overlap_count": 0,
                    "matched_metabolite_uids": [],
                    "score": 0.0,
                    "score_components": {},
                },
                {
                    "pathway_uid": "path_generic",
                    "name": "Regulation of PTEN stability and activity",
                    "overlap_count": 1,
                    "matched_metabolite_uids": ["met_bile_1"],
                    "score": 0.4,
                    "score_components": {"overlap_count": 1},
                },
            ]

            ranked = sorted(
                service.attach_pathway_user_support(pathways, features_by_uid),
                key=lambda row: -float(row.get("score") or 0.0),
            )

            self.assertEqual(ranked[0]["pathway_uid"], "path_bile")
            self.assertGreater(ranked[0]["score_components"]["input_theme_boost"], 0.0)
            self.assertEqual(ranked[0]["directional_support"]["input_theme_hits"][0]["theme_id"], "bile_acid")
            self.assertGreater(ranked[1]["score_components"]["generic_pathway_penalty"], 0.0)

    def test_unresolved_input_names_can_lift_generic_pathway_themes(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            pathways = [
                {
                    "pathway_uid": "path_succinate",
                    "name": "Succinate metabolism",
                    "overlap_count": 0,
                    "matched_metabolite_uids": [],
                    "score": 0.0,
                    "score_components": {},
                }
            ]

            ranked = service.attach_pathway_user_support(
                pathways,
                {},
                input_records=[
                    {"metabolite": "succinate"},
                    {"metabolite": "succinate semialdehyde"},
                ],
            )

            self.assertGreater(ranked[0]["score_components"]["input_lexical_boost"], 0.0)
            self.assertEqual(ranked[0]["directional_support"]["input_lexical_hits"][0]["term"], "succinate")

    def test_explain_analysis_wires_safe_adapter_as_read_only_narrator(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            records = [{"HMDB": "HMDB0000122", "log2FC": 2.0, "padj": 0.001, "direction": "up"}]

            input_pack = service.build_llm_adapter_input_pack(
                records,
                question="Summarize evidence for the analysis.",
                max_paths=10,
                max_hops=5,
                evidence_limit=5,
                subgraph_max_hops=1,
            )
            self.assertEqual(input_pack["contract_version"], "llm_safe_adapter.input.v1")
            self.assertEqual(
                {row["endpoint"] for row in input_pack["source_endpoints"]},
                {"/analyze/metabolites", "/evidence", "/subgraph", "/releases"},
            )
            self.assertNotIn("/resolve", {row["endpoint"] for row in input_pack["source_endpoints"]})
            self.assertEqual(input_pack["analysis_pack"]["contract_version"], "analysis_pack.v1")
            self.assertIn("structured_prediction", input_pack["analysis_pack"])
            self.assertEqual(input_pack["analysis_pack"]["structured_prediction"]["mode"], "generalized")
            self.assertEqual(input_pack["evidence"]["support"][0]["support_uid"], "litsup_1")
            self.assertTrue(input_pack["subgraph"]["nodes"])
            self.assertEqual(input_pack["release"]["release_id"], "mvp_20260101T000000")

            explained = service.explain_analysis(
                records,
                question="Summarize evidence for the analysis.",
                max_paths=10,
                max_hops=5,
                evidence_limit=5,
                subgraph_max_hops=1,
                include_input_pack=True,
            )
            self.assertEqual(explained["endpoint"], "/explain")
            self.assertEqual(explained["contract_version"], "explain_endpoint.v1")
            self.assertEqual(explained["status"], "ok")
            self.assertTrue(explained["guard_passed"])
            self.assertEqual(explained["explanation"]["contract_version"], "llm_safe_adapter.output.v1")
            self.assertEqual(explained["explanation"]["status"], "ok")
            self.assertIn("narrative_summary", explained["explanation"])
            first_section = explained["explanation"]["narrative_summary"]["sections"][0]
            self.assertTrue(first_section["source_refs"])
            self.assertEqual(explained["adapter_input_pack"]["contract_version"], "llm_safe_adapter.input.v1")

            repeated = service.explain_analysis(
                records,
                question="Summarize evidence for the analysis.",
                max_paths=10,
                max_hops=5,
                evidence_limit=5,
                subgraph_max_hops=1,
            )
            repeated_again = service.explain_analysis(
                records,
                question="Summarize evidence for the analysis.",
                max_paths=10,
                max_hops=5,
                evidence_limit=5,
                subgraph_max_hops=1,
            )
            self.assertEqual(
                repeated["determinism"]["response_hash"],
                repeated_again["determinism"]["response_hash"],
            )

            bad_output = {
                "contract_version": "llm_safe_adapter.output.v1",
                "adapter_version": "external_llm_bad",
                "status": "ok",
                "narrative_summary": {
                    "sections": [
                        {
                            "section_id": "bad",
                            "text": "I resolved metabolite_fake and updated p_final score to 0.99 for edge_fake.",
                            "source_refs": [
                                {"source_type": "evidence", "ref_id": "fake_ref", "path": "$.evidence.support[99]"}
                            ],
                        }
                    ]
                },
                "p_final": 0.99,
            }
            blocked = service.explain_analysis(records, adapter_output=bad_output)
            self.assertEqual(blocked["status"], "blocked_by_guard")
            self.assertFalse(blocked["guard_passed"])
            issue_codes = {issue["code"] for issue in blocked["explanation"]["guard"]["issues"]}
            self.assertIn("forbidden_output_field", issue_codes)
            self.assertIn("source_ref_not_in_input", issue_codes)
            self.assertIn("entity_resolution_decision", issue_codes)

    def test_chat_endpoint_normalizes_uploaded_table_and_explains(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            chat = service.chat(
                {
                    "table_text": "代谢物,log2FC,p值\nGlucose,2.0,0.001\n",
                    "question": "请解释这个结果。",
                    "adapter_backend": "local",
                    "max_paths": 10,
                    "max_hops": 5,
                    "evidence_limit": 5,
                }
            )

            self.assertEqual(chat["endpoint"], "/chat")
            self.assertEqual(chat["contract_version"], "chat_endpoint.v1")
            self.assertEqual(chat["status"], "ok")
            self.assertEqual(chat["input_normalization"]["record_count"], 1)
            self.assertEqual(chat["records"][0]["name"], "Glucose")
            self.assertEqual(chat["records"][0]["log2FC"], 2.0)
            self.assertEqual(chat["records"][0]["pvalue"], 0.001)
            self.assertEqual(chat["analysis_pack"]["input_summary"]["matched_count"], 1)
            self.assertEqual(chat["analysis_pack"]["pathway_rankings"][0]["pathway_uid"], "pathway_1")
            self.assertTrue(chat["explain_response"]["guard_passed"])
            self.assertIn("narrative_summary", chat["explanation"])

            missing = service.chat({"question": "这个结果说明什么？", "adapter_backend": "local"})
            self.assertEqual(missing["status"], "needs_records")
            self.assertEqual(missing["input_normalization"]["record_count"], 0)

    def test_chat_table_parser_repairs_unquoted_delimiter_in_name(self):
        records, metadata = metabo_service.parse_chat_table_text(
            "metabolite,log2FC,pvalue,padj,direction\n"
            "D-fructose 1,6-bisphosphate,0.64,0.0085,0.038,up\n"
        )

        self.assertEqual(records[0]["name"], "D-fructose 1,6-bisphosphate")
        self.assertEqual(records[0]["log2FC"], 0.64)
        self.assertEqual(records[0]["pvalue"], 0.0085)
        self.assertEqual(records[0]["padj"], 0.038)
        self.assertEqual(records[0]["direction"], "up")
        self.assertIn("row_width_repaired", {warning["code"] for warning in metadata["warnings"]})

    def test_cli_record_loader_accepts_csv_trait_score_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "trait_score.csv"
            path.write_text(
                "trait,group1,group2,pseudo_log2FC_shifted,p_ttest,q_ttest,direction\n"
                "GCST000001,Tumor,Adjacent,0.445,0,0,group1_higher\n",
                encoding="utf-8",
            )

            records = metabo_service.load_records_arg(str(path))
            converted, metadata = metabo_service.prepare_analysis_records(records)

            self.assertEqual(records[0]["trait"], "GCST000001")
            self.assertEqual(metadata["analysis_mode"], "two_group_trait_comparison_table")
            self.assertEqual(converted[0]["log2FC"], 0.445)
            self.assertEqual(converted[0]["direction"], "up")

    def test_compact_precheck_summary_preserves_mode_counts_and_warnings(self):
        precheck = {
            "endpoint": "/precheck/metabolites",
            "analysis_mode": "two_group_trait_comparison_table",
            "input_count": 2,
            "summary": {"matched": 1, "ambiguous": 1, "unmatched": 0, "invalid": 0},
            "expanded_summary": {
                "strict_matched_count": 1,
                "expanded_candidate_count": 2,
                "expanded_candidate_input_count": 1,
                "unresolved_input_count": 0,
                "analysis_seed_input_count": 2,
                "weighted_seed_mass": 1.5,
                "weight_policy": {"strict_identity": 1.0},
            },
            "input_normalization": {
                "source": "records",
                "record_count": 2,
                "analysis_mode": "two_group_trait_comparison_table",
                "input_format_detected": "two_group_trait_comparison_table",
                "converted_record_count": 2,
                "warnings": [{"code": "trait_score_comparison_not_direct_abundance"}],
            },
            "matched": [{"large": "row"}],
            "supported_identifier_columns": ["HMDB", "name"],
        }

        summary = metabo_service.compact_precheck_summary(precheck)

        self.assertEqual(summary["analysis_mode"], "two_group_trait_comparison_table")
        self.assertEqual(summary["summary"]["matched"], 1)
        self.assertEqual(summary["expanded_summary"]["analysis_seed_input_count"], 2)
        self.assertEqual(summary["warning_codes"], ["trait_score_comparison_not_direct_abundance"])
        self.assertNotIn("matched", summary)

    def test_chat_table_parser_preserves_european_trait_columns(self):
        records, metadata = metabo_service.parse_chat_table_text(
            "trait\tlog2FC\tp_ttest\tq_ttest\tdirection\n"
            "GCST000001\t1.25\t0.001\t0.01\tup\n"
        )

        self.assertEqual(metadata["record_count"], 1)
        self.assertEqual(records[0]["trait"], "GCST000001")
        self.assertEqual(records[0]["log2FC"], 1.25)
        self.assertEqual(records[0]["direction"], "up")

    def test_two_group_trait_table_converts_tumor_primary_direction(self):
        records, parse_metadata = metabo_service.parse_chat_table_text(
            "trait\tcelltype_l\tcelltype\tgroup_col\tgroup1\tgroup2\tlog2FC\tp_ttest\tq_ttest\tdirection\n"
            "GCST000001\tLUAD epithelial tissue\tEpi_tissue\tsample_type\tTumor\tAdjacent\t1.25\t0\t0.01\tgroup1_hi\n"
        )

        converted, metadata = metabo_service.prepare_analysis_records(records, parse_metadata)

        self.assertEqual(metadata["analysis_mode"], "two_group_trait_comparison_table")
        self.assertEqual(metadata["converted_record_count"], 1)
        self.assertEqual(converted[0]["accession_id"], "GCST000001")
        self.assertEqual(converted[0]["log2FC"], 1.25)
        self.assertEqual(converted[0]["pvalue"], 1e-300)
        self.assertEqual(converted[0]["padj"], 0.01)
        self.assertEqual(converted[0]["direction"], "up")
        self.assertEqual(converted[0]["comparison_primary_group"], "Tumor")
        self.assertEqual(converted[0]["cell_type"], "LUAD epithelial tissue")
        self.assertEqual(converted[0]["source_analysis_mode"], "two_group_trait_comparison_table")

    def test_two_group_trait_table_flips_adjacent_group1_against_tumor(self):
        records, parse_metadata = metabo_service.parse_chat_table_text(
            "trait\tcelltype\tgroup1\tgroup2\tlog2FC\tq_wilcoxon\tdirection\n"
            "GCST000001\tEpi_tissue\tAdjacent\tTumor\t1.25\t0.02\tgroup1_hi\n"
        )

        converted, metadata = metabo_service.prepare_analysis_records(records, parse_metadata)

        self.assertEqual(metadata["analysis_mode"], "two_group_trait_comparison_table")
        self.assertEqual(converted[0]["comparison_primary_group"], "Tumor")
        self.assertEqual(converted[0]["comparison_reference_group"], "Adjacent")
        self.assertEqual(converted[0]["log2FC"], -1.25)
        self.assertEqual(converted[0]["padj"], 0.02)
        self.assertEqual(converted[0]["direction"], "down")

    def test_two_group_trait_table_uses_pseudo_log2fc_when_log2fc_blank(self):
        records, parse_metadata = metabo_service.parse_chat_table_text(
            "trait,group1,group2,log2FC,pseudo_log2FC_shifted,p_ttest,q_ttest,direction\n"
            "GCST000001,Tumor,Adjacent,,0.445,0,0,group1_higher\n"
        )

        converted, metadata = metabo_service.prepare_analysis_records(records, parse_metadata)

        self.assertEqual(metadata["analysis_mode"], "two_group_trait_comparison_table")
        self.assertEqual(converted[0]["log2FC"], 0.445)
        self.assertEqual(converted[0]["direction"], "up")
        self.assertEqual(converted[0]["pvalue"], 1e-300)
        self.assertEqual(converted[0]["padj"], 1e-300)
        self.assertEqual(converted[0]["comparison_effect_source"], "pseudo_log2FC_shifted_oriented_by_group_direction")
        self.assertEqual(converted[0]["effect_label"], "pseudo_log2FC_shifted")
        self.assertEqual(converted[0]["effect_value"], 0.445)
        self.assertEqual(converted[0]["log2fc_semantics"], "directional_effect_surrogate_not_abundance_log2fc")
        self.assertTrue(any(row["code"] == "trait_score_comparison_not_direct_abundance" for row in metadata["warnings"]))

    def test_differential_table_selector_accepts_precomputed_effects_without_abundance(self):
        records, _metadata = metabo_service.parse_chat_table_text(
            "metabolite,mean_diff,q_wilcoxon_global,direction,celltype\n"
            "Glucose,0.62,0.01,group1_higher,LUAD_Epi\n"
            "Lactate,0.10,0.20,group1_higher,LUAD_Epi\n"
        )

        selected, selection = metabo_service.select_differential_table_records(
            records,
            q_threshold=0.05,
            min_abs_effect=0.2,
            max_records=10,
        )

        self.assertEqual(selection["contract_version"], "differential_table_selection.v1")
        self.assertEqual(selection["selected_count"], 1)
        self.assertEqual(selection["rejected_counts"], {"q_above_threshold": 1})
        self.assertEqual(selected[0]["name"], "Glucose")
        self.assertEqual(selected[0]["effect_label"], "mean_diff")
        self.assertEqual(selected[0]["log2FC"], 0.62)
        self.assertEqual(selected[0]["log2fc_semantics"], "directional_effect_surrogate_not_abundance_log2fc")
        self.assertEqual(selected[0]["source_analysis_mode"], "differential_table")

    def test_ratio_trait_descriptor_preserves_numerator_denominator_direction_policy(self):
        ratio = metabo_service.ratio_trait_descriptor("Glucose to Lactate ratio", "pseudo_log2FC_shifted", 0.7, "up")

        self.assertTrue(ratio["dual_component_ratio"])
        self.assertEqual([row["name"] for row in ratio["components"]], ["Glucose", "Lactate"])
        self.assertEqual(ratio["components"][0]["role"], "numerator")
        self.assertEqual(ratio["components"][0]["signed_ratio_direction"], "same_as_ratio")
        self.assertEqual(ratio["components"][1]["role"], "denominator")
        self.assertEqual(ratio["components"][1]["signed_ratio_direction"], "opposite_to_ratio")
        self.assertIn("does not prove numerator abundance", ratio["interpretation_boundary"])

    def test_ratio_component_features_flip_denominator_direction_without_abundance_claim(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            ratio = metabo_service.ratio_trait_descriptor("Glucose to Lactate ratio", "pseudo_log2FC_shifted", 0.7, "up")
            denominator = ratio["components"][1]
            rows = [
                {
                    "input_id": "ratio_row::expanded_denominator",
                    "metabolite_uid": "met_lactate",
                    "record": {
                        "reported_trait": "Glucose to Lactate ratio",
                        "name": "Lactate",
                        "log2FC": 0.7,
                        "padj": 0.01,
                        "direction": "up",
                        "ratio_trait": ratio,
                    },
                    "resolution": {"candidates": [{"display_name": "Lactate", "primary_external_id": "CHEBI:16651"}]},
                    "expanded_candidate": {
                        "track": "expanded",
                        "seed_class": "ratio_component",
                        "weight_multiplier": 0.25,
                        "ratio_trait": {key: value for key, value in ratio.items() if key != "components"},
                        "ratio_component": {
                            "name": denominator["name"],
                            "role": denominator["role"],
                            "signed_ratio_direction": denominator["signed_ratio_direction"],
                            "component_weight": denominator["component_weight"],
                        },
                    },
                }
            ]

            features, by_uid = service.matched_analysis_features(rows)

            self.assertEqual(features[0]["direction"], "down")
            self.assertEqual(features[0]["log2_fold_change"], -0.7)
            self.assertEqual(features[0]["ratio_original_direction"], "up")
            self.assertEqual(features[0]["ratio_component"]["role"], "denominator")
            self.assertEqual(features[0]["ratio_effect_semantics"], "relative_ratio_consistent_direction_not_measured_abundance")
            self.assertEqual(by_uid["met_lactate"][0]["direction"], "down")

    def test_class_or_pool_candidate_gets_class_seed_constraints(self):
        descriptor = metabo_service.class_seed_descriptor(
            {"reported_trait": "Sphingomyelin (d18:2/21:0, d16:2/23:0) levels", "hmdb": "HMDB0000001"},
            {"display_name": "SM(d18:2/21:0)", "external_xrefs": ["HMDB:HMDB0000001"]},
        )

        self.assertEqual(descriptor["contract_version"], "class_seed.v1")
        self.assertIn("sphingomyelin", descriptor["class_ids"])
        self.assertIn("d18:2/21:0", descriptor["chain_constraints"])
        self.assertIn("hmdb", descriptor["stable_identifier_fields"])
        self.assertEqual(descriptor["allowed_claim_level"], "candidate_lipid_species_with_chain_constraints")

    def test_soft_identity_candidates_are_clustered_by_shared_identifier(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            row = {
                "input_id": "soft_row",
                "record": {"name": "Alpha-hydroxyisovalerate", "log2FC": 0.4, "padj": 0.01, "direction": "up"},
                "resolution": {
                    "status": "ambiguous",
                    "top_score": 100.0,
                    "top_margin": 0.0,
                    "candidates": [
                        {
                            "entity_uid": "met_a",
                            "display_name": "2-hydroxy-3-methylbutyrate",
                            "primary_external_id": "HMDB:HMDB0000407",
                            "score": 100.0,
                            "external_xrefs": ["HMDB:HMDB0000407"],
                            "matches": [{"namespace": "HMDB", "raw_value": "HMDB0000407"}],
                        },
                        {
                            "entity_uid": "met_b",
                            "display_name": "2-hydroxy-3-methylbutyric acid",
                            "primary_external_id": "HMDB:HMDB0000407",
                            "score": 100.0,
                            "external_xrefs": ["HMDB:HMDB0000407"],
                            "matches": [{"namespace": "HMDB", "raw_value": "HMDB0000407"}],
                        },
                    ],
                },
            }

            expanded = service.expanded_rows_for_resolution_row(row)
            features, _by_uid = service.matched_analysis_features(expanded)

            self.assertEqual(len(expanded), 1)
            cluster = expanded[0]["expanded_candidate"]["identity_cluster"]
            self.assertEqual(cluster["cluster_basis"], "hmdb")
            self.assertEqual(cluster["candidate_count"], 2)
            self.assertEqual(set(cluster["candidate_uids"]), {"met_a", "met_b"})
            self.assertEqual(features[0]["identity_effect_semantics"], "identity_cluster_or_consensus_only")

    def test_analyze_differential_table_uses_selected_difference_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            records = [
                {"metabolite": "Glucose", "mean_diff": 0.62, "q_wilcoxon_global": 0.01, "direction": "up"},
                {"metabolite": "Shared", "mean_diff": 0.2, "q_wilcoxon_global": 0.20, "direction": "up"},
            ]

            analyzed = service.analyze_differential_table(records, max_paths=5, max_hops=2, q_threshold=0.05)

            self.assertEqual(analyzed["endpoint"], "/analyze/differential-table")
            self.assertEqual(analyzed["differential_table_selection"]["selected_count"], 1)
            self.assertEqual(analyzed["precheck"]["analysis_mode"], "differential_table")
            self.assertEqual(analyzed["analysis_pack"]["input_summary"]["matched_count"], 1)
            self.assertEqual(analyzed["selected_records"][0]["effect_label"], "mean_diff")
            self.assertIn("raw abundance", analyzed["differential_table_selection"]["selection_rule"])
            self.assertIn("explanation_paths", analyzed)
            self.assertIn("compressed_explanation_paths", analyzed)
            self.assertTrue(analyzed["analysis_pack"]["top_explanation_paths"])
            self.assertLessEqual(
                len(analyzed["analysis_pack"]["top_explanation_paths"]),
                min(len(analyzed["explanation_paths"]), metabo_service.ANALYSIS_PACK_PATH_LIMIT),
            )

    def test_differential_chat_mode_wraps_selected_analysis_for_explanation(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            chat = service.chat(
                {
                    "input_mode": "trait_score_differential",
                    "records": [
                        {"metabolite": "Glucose", "mean_diff": 0.62, "q_wilcoxon_global": 0.01, "direction": "up"},
                        {"metabolite": "Shared", "mean_diff": 0.2, "q_wilcoxon_global": 0.20, "direction": "up"},
                    ],
                    "question": "Explain this differential table.",
                    "adapter_backend": "local",
                    "max_paths": 5,
                    "max_hops": 2,
                    "q_threshold": 0.05,
                    "evidence_limit": 5,
                }
            )

            self.assertEqual(chat["endpoint"], "/chat")
            self.assertEqual(chat["status"], "ok")
            self.assertEqual(chat["input_mode"], "trait_score_differential")
            self.assertEqual(chat["differential_table_selection"]["selected_count"], 1)
            self.assertEqual(chat["records"][0]["name"], "Glucose")
            self.assertEqual(chat["records"][0]["effect_label"], "mean_diff")
            self.assertEqual(chat["analysis_pack"]["input_summary"]["matched_count"], 1)
            self.assertEqual(chat["analysis_pack"]["differential_table_selection"]["selected_count"], 1)
            self.assertIn("narrative_summary", chat["explanation"])
            self.assertTrue(chat["explain_response"]["guard_passed"])
            self.assertEqual(
                chat["explain_response"]["adapter_input_source_endpoints"][0]["endpoint"],
                "/analyze/differential-table",
            )

    def test_precheck_enriches_trait_from_european_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            european = root / "raw_lake" / "European"
            european.mkdir(parents=True)
            (european / "European.csv").write_text(
                "accessionId,reportedTrait,accessionId\n"
                "GCST000001,Unknown PubChemBacked,GCST000001\n",
                encoding="utf-8",
            )
            (european / "European_trait_annotations.csv").write_text(
                "accession_id,reported_trait,resolution_status,mapped_names,local_name_hits,pubchem_cids,local_pubchem_cid_hits,pubchem_titles,pubchem_formulas,pubchem_inchikeys,candidate_names,pubchem_query,source_name\n"
                "GCST000001,Unknown PubChemBacked,pubchem_cid_in_local_index,Unknown PubChemBacked,,5793,5793,Glucose,C6H12O6,,,,test\n",
                encoding="utf-8",
            )
            service = self.make_service(root)

            precheck = service.precheck_metabolites([{"trait": "CGST000001", "log2FC": 1.25, "direction": "up"}])

            self.assertEqual(precheck["summary"]["matched"], 1)
            record = precheck["matched"][0]["record"]
            self.assertEqual(record["accession_id"], "GCST000001")
            self.assertEqual(record["reported_trait"], "Unknown PubChemBacked")
            self.assertEqual(record["pubchem_cid"], "5793")

    def test_chat_endpoint_extracts_records_from_loose_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            chat = service.chat(
                {
                    "free_text": "请解释：Glucose 上调 log2FC=2.0 p=0.001；Shared 下调 p=0.02",
                    "adapter_backend": "local",
                    "max_paths": 10,
                    "max_hops": 5,
                }
            )

            self.assertEqual(chat["input_normalization"]["source"], "free_text")
            self.assertEqual(chat["input_normalization"]["record_count"], 2)
            self.assertEqual(chat["records"][0]["name"], "Glucose")
            self.assertEqual(chat["records"][0]["direction"], "up")
            self.assertEqual(chat["records"][0]["log2FC"], 2.0)
            self.assertEqual(chat["records"][0]["pvalue"], 0.001)
            self.assertEqual(chat["records"][1]["name"], "Shared")
            self.assertEqual(chat["records"][1]["direction"], "down")
            self.assertEqual(chat["analysis_pack"]["input_summary"]["matched_count"], 1)
            self.assertEqual(chat["analysis_pack"]["input_summary"]["ambiguous_count"], 1)

    def test_explain_analysis_external_llm_backend_is_opt_in_and_guarded(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            records = [{"HMDB": "HMDB0000122", "log2FC": 2.0, "padj": 0.001, "direction": "up"}]

            disabled = service.explain_analysis(records, adapter_backend="external_llm")
            self.assertEqual(disabled["status"], "blocked_by_guard")
            disabled_codes = {issue["code"] for issue in disabled["explanation"]["guard"]["issues"]}
            self.assertIn("external_llm_disabled", disabled_codes)

            def safe_transport(config, payload):
                self.assertEqual(config.model, "test-model")
                self.assertEqual(payload["temperature"], 0)
                self.assertNotIn("response_format", payload)
                return {
                    "choices": [
                        {
                            "message": {
                                "content": "The external narrator rewrites only the frozen analysis and preserves matching, ranking, evidence, and research-only boundaries."
                            }
                        }
                    ]
                }

            service.llm_config = metabo_service.ExternalLLMConfig(
                enabled=True,
                provider="openai_compatible",
                endpoint="https://llm.example.test/v1/chat/completions",
                model="test-model",
                api_key="secret",
            )
            service.llm_transport = safe_transport
            explained = service.explain_analysis(records, adapter_backend="external_llm")
            self.assertEqual(explained["status"], "ok")
            self.assertTrue(explained["guard_passed"])
            self.assertEqual(explained["adapter_backend"], "external_llm")
            audit = explained["explanation"]["determinism"]["backend_audit"]
            self.assertEqual(audit["backend"], "external_llm")
            self.assertEqual(audit["model"], "test-model")
            self.assertIn("text_fixed_contract", audit["adapter_version"])
            self.assertIn("prompt_hash", audit)

            def unsafe_transport(_config, _payload):
                return {
                    "choices": [
                        {
                            "message": {
                                "content": "I resolved metabolite_fake and updated p_final score to 0.99 for edge_fake."
                            }
                        }
                    ]
                }

            service.llm_transport = unsafe_transport
            blocked = service.explain_analysis(records, adapter_backend="external_llm")
            self.assertEqual(blocked["status"], "blocked_by_guard")
            self.assertFalse(blocked["guard_passed"])
            blocked_codes = {issue["code"] for issue in blocked["explanation"]["guard"]["issues"]}
            self.assertIn("entity_resolution_decision", blocked_codes)
            self.assertIn("mutates_scores_or_graph", blocked_codes)

    def test_explain_analysis_external_text_backend_accepts_plain_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            records = [{"HMDB": "HMDB0000122", "log2FC": 2.0, "padj": 0.001, "direction": "up"}]

            def text_transport(config, payload):
                self.assertEqual(config.model, "test-model")
                self.assertNotIn("response_format", payload)
                return {
                    "choices": [
                        {
                            "message": {
                                "content": "这份结果应先看匹配、通路排序和证据来源；它不是诊断或治疗建议。"
                            }
                        }
                    ]
                }

            service.llm_config = metabo_service.ExternalLLMConfig(
                enabled=True,
                provider="openai_compatible",
                endpoint="https://llm.example.test/v1/chat/completions",
                model="test-model",
                api_key="secret",
            )
            service.llm_transport = text_transport
            explained = service.explain_analysis(records, adapter_backend="external_text")
            self.assertEqual(explained["status"], "ok")
            self.assertTrue(explained["guard_passed"])
            self.assertEqual(explained["adapter_backend"], "external_text")
            self.assertIn("narrative_summary", explained["explanation"])
            self.assertEqual(explained["explanation"]["determinism"]["backend_audit"]["backend"], "external_text")

    def test_explain_analysis_external_text_chunked_backend_guards_each_chunk(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            records = [{"HMDB": "HMDB0000122", "log2FC": 2.0, "padj": 0.001, "direction": "up"}]
            calls = []

            def chunked_transport(config, payload):
                self.assertEqual(config.model, "test-model")
                self.assertNotIn("response_format", payload)
                user_text = payload["messages"][1]["content"]
                calls.append(user_text)
                if "总述" in user_text:
                    content = "总述：这些通过校验的分块只支持研究解释和后续验证优先级。"
                else:
                    content = "本分块仅解释已有结构化结果，保留置信度和边界，不作诊断或治疗建议。"
                return {"choices": [{"message": {"content": content}}]}

            service.llm_config = metabo_service.ExternalLLMConfig(
                enabled=True,
                provider="openai_compatible",
                endpoint="https://llm.example.test/v1/chat/completions",
                model="test-model",
                api_key="secret",
            )
            service.llm_transport = chunked_transport
            explained = service.explain_analysis(records, adapter_backend="external_text_chunked")

            self.assertEqual(explained["status"], "ok")
            self.assertTrue(explained["guard_passed"])
            self.assertGreaterEqual(len(calls), 2)
            audit = explained["explanation"]["determinism"]["backend_audit"]
            self.assertEqual(audit["backend"], "external_text_chunked")
            self.assertGreaterEqual(audit["chunk_count"], 1)
            self.assertEqual(audit["passed_chunk_count"], audit["chunk_count"])
            sections = explained["explanation"]["narrative_summary"]["sections"]
            self.assertEqual(sections[0]["section_id"], "chunk_final_summary")
            self.assertTrue(all(section.get("source_refs") for section in sections))

    def test_chat_accepts_request_scoped_external_llm_config_without_key_leak(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))

            def text_transport(config, payload):
                self.assertTrue(config.enabled)
                self.assertEqual(config.endpoint, "https://llm.example.test/v1/chat/completions")
                self.assertEqual(config.model, "request-model")
                self.assertEqual(config.api_key, "request-secret")
                self.assertNotIn("response_format", payload)
                return {
                    "choices": [
                        {
                            "message": {
                                "content": "外部模型只按已有证据改写解释，不能修改匹配、评分或图谱。"
                            }
                        }
                    ]
                }

            service.llm_transport = text_transport
            chat = service.chat(
                {
                    "records": [{"HMDB": "HMDB0000122", "log2FC": 2.0, "padj": 0.001, "direction": "up"}],
                    "question": "请解释这个结果。",
                    "adapter_backend": "external_text",
                    "llm_config": {
                        "endpoint": "https://llm.example.test/v1/chat/completions",
                        "model": "request-model",
                        "api_key": "request-secret",
                    },
                }
            )

            self.assertEqual(chat["status"], "ok")
            self.assertTrue(chat["explain_response"]["guard_passed"])
            audit = chat["explanation"]["determinism"]["backend_audit"]
            self.assertEqual(audit["backend"], "external_text")
            self.assertEqual(audit["model"], "request-model")
            self.assertNotIn("request-secret", json.dumps(chat, ensure_ascii=False))

    def test_llm_test_accepts_request_scoped_config_without_key_leak(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            calls = []

            def ok_transport(config, payload):
                calls.append((config, payload))
                self.assertEqual(config.model, "request-model")
                self.assertEqual(config.api_key, "request-secret")
                self.assertEqual(config.proxy_url, "http://127.0.0.1:7890")
                self.assertLessEqual(payload["max_tokens"], 32)
                return {"choices": [{"message": {"content": "OK"}}], "usage": {"total_tokens": 3}}

            service.llm_transport = ok_transport
            result = service.llm_test(
                {
                    "llm_config": {
                        "endpoint": "https://llm.example.test/v1/chat/completions",
                        "model": "request-model",
                        "api_key": "request-secret",
                        "proxy_url": "http://127.0.0.1:7890",
                        "max_output_tokens": 1600,
                    }
                }
            )

            self.assertEqual(result["status"], "ok")
            self.assertTrue(result["llm_test"]["ready"])
            self.assertEqual(result["llm_test"]["code"], "ok")
            self.assertEqual(len(calls), 1)
            serialized = json.dumps(result, ensure_ascii=False)
            self.assertNotIn("request-secret", serialized)
            self.assertIn("api_key_configured", serialized)

    def test_llm_test_reports_configuration_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))

            result = service.llm_test(
                {
                    "llm_config": {
                        "endpoint": "https://llm.example.test/v1/chat/completions",
                        "api_key": "request-secret",
                    }
                }
            )

            self.assertEqual(result["status"], "blocked")
            self.assertFalse(result["llm_test"]["ready"])
            self.assertEqual(result["llm_test"]["phase"], "config")
            self.assertEqual(result["llm_test"]["code"], "missing_llm_model")
            self.assertIn("配置", result["llm_test"]["diagnostic"])
            self.assertNotIn("request-secret", json.dumps(result, ensure_ascii=False))

    def test_llm_test_reports_http_auth_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))

            def failing_transport(_config, _payload):
                raise metabo_service.LLMAdapterError(
                    "external_llm_http_error",
                    "External LLM request failed.",
                    {"status": 401, "body": "{\"error\":\"Unauthorized\"}"},
                )

            service.llm_transport = failing_transport
            result = service.llm_test(
                {
                    "llm_config": {
                        "endpoint": "https://llm.example.test/v1/chat/completions",
                        "model": "request-model",
                        "api_key": "request-secret",
                    }
                }
            )

            self.assertEqual(result["status"], "failed")
            self.assertFalse(result["llm_test"]["ready"])
            self.assertEqual(result["llm_test"]["code"], "external_llm_http_error")
            self.assertEqual(result["llm_test"]["http_status"], 401)
            self.assertIn("鉴权失败", result["llm_test"]["diagnostic"])
            self.assertNotIn("request-secret", json.dumps(result, ensure_ascii=False))

    def test_external_llm_transport_failure_falls_back_to_local_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))

            def failing_transport(_config, _payload):
                raise ConnectionResetError("connection reset for test")

            service.llm_transport = failing_transport
            chat = service.chat(
                {
                    "records": [{"HMDB": "HMDB0000122", "log2FC": 2.0, "padj": 0.001, "direction": "up"}],
                    "question": "请解释这个结果。",
                    "adapter_backend": "external_text",
                    "llm_config": {
                        "endpoint": "https://llm.example.test/v1/chat/completions",
                        "model": "request-model",
                        "api_key": "request-secret",
                    },
                }
            )

            self.assertEqual(chat["status"], "ok")
            self.assertTrue(chat["explain_response"]["guard_passed"])
            self.assertEqual(chat["explain_response"]["external_fallback_issue"]["code"], "external_llm_unhandled_error")
            audit = chat["explanation"]["determinism"]["backend_audit"]
            self.assertEqual(audit["fallback_backend"], "local")
            self.assertEqual(audit["backend"], "external_text")

    def test_biological_entity_pools_and_theme_coverage_keep_ambiguous_signals(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            records = [
                {"HMDB": "HMDB0000122", "metabolite": "Glucose", "log2FC": 1.1, "padj": 0.01},
                {"metabolite": "UDP-glucose", "log2FC": 0.8, "padj": 0.03},
                {"metabolite": "UDP-N-acetylglucosamine", "log2FC": 0.7, "padj": 0.04},
                {"metabolite": "phosphatidylcholine", "log2FC": -0.6, "padj": 0.05},
                {"metabolite": "taurine", "log2FC": 0.5, "padj": 0.02},
            ]

            analyzed = service.analyze_metabolites(records, max_paths=10, max_hops=5)
            pack = analyzed["analysis_pack"]
            pool_ids = {row["pool_id"] for row in pack["biological_entity_pools"]}
            self.assertIn("udp_sugar_pool", pool_ids)
            self.assertIn("phospholipid_membrane_pool", pool_ids)

            themes = {row["theme_id"]: row for row in pack["prediction_model"]["metabolic_themes"] if row.get("theme_id")}
            glycosylation = themes["nucleotide_sugar_glycosylation"]
            self.assertEqual(glycosylation["coverage"]["related_input_count"], 2)
            self.assertEqual(glycosylation["coverage"]["unmatched_support_count"], 2)
            self.assertEqual(glycosylation["confidence_tier"], "exploratory")
            self.assertIn("input_signal_without_stable_exact_match", glycosylation["coverage"]["downgrade_reason"])
            self.assertIn("no_matched_support", glycosylation["coverage"]["downgrade_reason"])
            self.assertIn(glycosylation, pack["prediction_model"]["degraded_covered_metabolic_themes"])

    def test_common_biochemical_exact_name_rescue_reduces_meaningless_ambiguity(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            cases = [
                ("UDP-N-acetylglucosamine", "UDP-N-acetyl-alpha-D-glucosamine"),
                ("oxidized glutathione", "glutathione disulfide"),
                ("glucose 6-phosphate", "D-glucopyranose 6-phosphate"),
                ("succinate", "succinate(2-)"),
                ("NAD+", "NAD(+)"),
                ("palmitoylcarnitine", "O-palmitoyl-L-carnitine"),
                ("sphingosine-1-phosphate", "sphingosine 1-phosphate"),
            ]
            notes = []
            for query, candidate_name in cases:
                candidates = [
                    {
                        "entity_uid": f"met_{query}",
                        "score": 55.0,
                        "display_name": candidate_name,
                        "score_components": {"name": 55.0},
                        "matches": [],
                    }
                ]
                status = service.maybe_rescue_common_biochemical_match(
                    "ambiguous",
                    candidates,
                    {"names": [query]},
                    notes,
                )

                self.assertEqual(status, "matched")
                self.assertEqual(candidates[0]["matches"][-1]["component"], "biological_exact_rescue")
            self.assertEqual(len([note for note in notes if note.get("mode") == "biological_exact_rescue"]), len(cases))

    def test_common_biochemical_zero_score_rescue_handles_conservative_name_abstention(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            notes = []
            candidates = [
                {
                    "entity_uid": "met_residue",
                    "score": 0.0,
                    "display_name": "L-glutamine residue",
                    "score_components": {"name": 0.0},
                    "matches": [
                        {
                            "component": "name",
                            "raw_value": "L-glutamine",
                            "triage_policy": "auto_abstain_name_only",
                            "match_field": "synonym",
                        }
                    ],
                },
                {
                    "entity_uid": "met_glutamine",
                    "score": 0.0,
                    "display_name": "L-glutamine",
                    "score_components": {"name": 0.0},
                    "matches": [
                        {
                            "component": "name",
                            "raw_value": "L-glutamine",
                            "triage_policy": "auto_abstain_name_only",
                            "match_field": "canonical_name",
                        }
                    ],
                },
            ]

            status = service.maybe_rescue_common_biochemical_match(
                "ambiguous",
                candidates,
                {"names": ["Glutamine"]},
                notes,
            )

            self.assertEqual(status, "matched")
            self.assertEqual(candidates[0]["entity_uid"], "met_glutamine")
            self.assertGreaterEqual(candidates[0]["score"], service.config.min_match_score)
            self.assertEqual(candidates[0]["matches"][-1]["component"], "common_biochemical_zero_score_rescue")
            self.assertEqual(notes[-1]["mode"], "common_biochemical_zero_score_rescue")

    def test_biochemical_rescue_query_classifier_covers_common_metabolism_not_sample_only(self):
        common_queries = [
            "glucose-6-phosphate",
            "fructose-1,6-bisphosphate",
            "alpha-ketoglutarate",
            "fumarate",
            "L-serine",
            "S-adenosylmethionine",
            "glutathione disulfide",
            "NADPH",
            "ATP",
            "UDP-glucuronic acid",
            "CDP-choline",
            "acetyl-CoA",
            "propionylcarnitine",
            "lysophosphatidylcholine",
            "sphingomyelin",
            "prostaglandin D2",
            "myo-inositol",
        ]
        non_metabolic_queries = [
            "random protein family",
            "immune cell receptor",
            "taste perception pathway",
        ]

        for query in common_queries:
            self.assertTrue(metabo_service.is_common_biochemical_rescue_query(query), query)
        for query in non_metabolic_queries:
            self.assertFalse(metabo_service.is_common_biochemical_rescue_query(query), query)

    def test_epithelial_context_mismatch_downgrades_but_keeps_prediction(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            packed = {
                "rank": 1,
                "pathway_uid": "pathway_1",
                "display_name": "Neuronal taste perception pathway",
                "score": 3.0,
                "score_components": {},
                "claim_refs": {"traceability_passed": True},
            }
            source_row = {
                "pathway_uid": "pathway_1",
                "name": "Neuronal taste perception pathway",
                "overlap_count": 2,
                "matched_metabolite_uids": ["met_1"],
                "score_components": {},
            }
            features_by_uid = {
                "met_1": [
                    {
                        "direction": "up",
                        "significant": True,
                        "p_user": 0.7,
                        "input_name": "glucose",
                        "display_name": "Glucose",
                    }
                ]
            }
            input_summary = {"input_count": 1, "matched_count": 1, "ambiguous_count": 0, "unmatched_count": 0, "invalid_count": 0}

            calibrated = service.calibrated_prediction_for_ranking_row(
                packed,
                source_row,
                "pathway_uid",
                None,
                features_by_uid,
                input_summary,
                {"fields": {"cell_type": ["epithelial cell"]}, "terms": ["epithelial", "cell"], "has_context": True},
            )

            self.assertTrue(calibrated["context_mismatch"])
            self.assertEqual(calibrated["calibration_status"], "context_mismatch")
            self.assertEqual(calibrated["confidence_tier"], "exploratory")
            self.assertIn("neuronal", calibrated["context_mismatch_groups"])

    def test_disease_cancer_context_mismatch_downgrades_wrong_cancer(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            context = service.normalize_prediction_context(
                {"cancer_type": "colorectal cancer", "tissue": "colon", "context_terms": ["CRC"]}
            )

            calibrated = service.calibrated_prediction_for_ranking_row(
                {
                    "rank": 1,
                    "disease_uid": "disease_melanoma",
                    "display_name": "melanoma",
                    "score": 5.0,
                    "score_components": {"propagation_score": 0.1},
                    "claim_refs": {"traceability_passed": True},
                },
                {
                    "disease_uid": "disease_melanoma",
                    "name": "melanoma",
                    "score_components": {"propagation_score": 0.1},
                },
                "disease_uid",
                None,
                {},
                {"input_count": 10, "matched_count": 8, "ambiguous_count": 0, "unmatched_count": 0, "invalid_count": 0},
                context,
            )

            self.assertTrue(calibrated["context_mismatch"])
            self.assertEqual(calibrated["calibration_status"], "context_mismatch")
            self.assertIn(calibrated["confidence_tier"], {"exploratory", "low"})
            self.assertLessEqual(calibrated["calibrated_confidence"], 0.25)
            self.assertIn("melanoma", calibrated["context_mismatch_groups"])
            self.assertIn("disease_context_mismatch_for_current_cancer_background", calibrated["downgrade_reason"])

    def test_trait_score_pathway_high_confidence_requires_stronger_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            calibrated = service.calibrated_prediction_for_ranking_row(
                {
                    "rank": 1,
                    "pathway_uid": "pathway_1",
                    "display_name": "Sphingolipid catabolism",
                    "score": 10.0,
                    "score_components": {"input_theme_boost": 4.0, "propagation_score": 0.01},
                    "claim_refs": {"traceability_passed": True},
                    "literature_support_count": 0,
                },
                {
                    "pathway_uid": "pathway_1",
                    "name": "Sphingolipid catabolism",
                    "overlap_count": 1,
                    "matched_metabolite_uids": ["met_1"],
                    "directional_support": {"up": 1, "down": 0, "significant": 1},
                    "score_components": {"input_theme_boost": 4.0, "propagation_score": 0.01},
                },
                "pathway_uid",
                None,
                {
                    "met_1": [
                        {
                            "direction": "up",
                            "significant": True,
                            "p_user": 0.9,
                            "input_name": "sphingomyelin",
                            "display_name": "sphingomyelin",
                        }
                    ]
                },
                {
                    "analysis_mode": "two_group_trait_comparison_table",
                    "input_count": 80,
                    "matched_count": 60,
                    "ambiguous_count": 5,
                    "unmatched_count": 0,
                    "invalid_count": 0,
                },
            )

            self.assertEqual(calibrated["confidence_tier"], "medium")
            self.assertLessEqual(calibrated["calibrated_confidence"], 0.74)
            self.assertIn("trait_score_requires_multi_seed_or_literature_support", calibrated["downgrade_reason"])

    def test_drug_predictions_inherit_low_target_confidence_and_appendix_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            service = self.make_service(root)
            overlay = root / "manual_sources" / "prediction_overlays" / "mvp_20260101T000000"
            overlay.mkdir(parents=True, exist_ok=True)
            (overlay / "drug_targets.csv").write_text(
                "\n".join(
                    [
                        "drug_id,drug_name,target_symbol,mechanism,confidence,context_terms,source_name,source_record_id,evidence_level,license_id",
                        "drug_hk1,HK1 test inhibitor,HK1,inhibitor,0.95,cancer,manual_seed,drug_hk1|HK1,predicted,local:test",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            drugs = service.build_drug_predictions(
                [
                    {
                        "target_uid": "target_1",
                        "display_name": "HK1",
                        "primary_external_id": "OPENTARGETS:ENSG1",
                        "score": 10.0,
                        "confidence_tier": "low",
                        "calibrated_confidence": 0.05,
                        "score_components": {"propagation_score": 1.0},
                    }
                ],
                {"fields": {"cancer_type": ["cancer"]}, "terms": ["cancer"], "has_context": True},
                "soft",
            )

            self.assertEqual(drugs[0]["confidence_tier"], "low")
            self.assertEqual(drugs[0]["display_confidence"], "low")
            self.assertEqual(drugs[0]["upstream_confidence_tier"], "low")
            self.assertGreater(drugs[0]["raw_overlay_score"], 0.0)
            self.assertTrue(drugs[0]["appendix"])
            self.assertIn("drug_confidence_inherits_target_confidence", drugs[0]["downgrade_reason"])

    def test_theme_confidence_gates_enforce_minimum_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))

            score, tier, reasons = service.apply_theme_confidence_gates(
                0.9,
                {
                    "matched_support_count": 0,
                    "significance_support_count": 2,
                    "direction_consistency": 1.0,
                },
            )
            self.assertEqual(tier, "exploratory")
            self.assertLessEqual(score, 0.44)
            self.assertIn("no_matched_support", reasons)

            score, tier, reasons = service.apply_theme_confidence_gates(
                0.9,
                {
                    "matched_support_count": 3,
                    "significance_support_count": 0,
                    "direction_consistency": 1.0,
                },
            )
            self.assertEqual(tier, "medium")
            self.assertLessEqual(score, 0.74)
            self.assertIn("no_significant_support", reasons)

            score, tier, reasons = service.apply_theme_confidence_gates(
                0.9,
                {
                    "matched_support_count": 3,
                    "significance_support_count": 2,
                    "direction_consistency": 0.0,
                },
            )
            self.assertEqual(tier, "exploratory")
            self.assertLessEqual(score, 0.44)
            self.assertIn("no_direction_consistency", reasons)

    def test_graph_loss_warning_prevents_overall_a_calibration(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            assessment = service.prediction_model_assessment(
                {"input_count": 10, "summary": {"matched": 10, "ambiguous": 0, "unmatched": 0, "invalid": 0}},
                [
                    {
                        "result_type": "metabolic_theme",
                        "confidence_tier": "high",
                        "calibrated_confidence": 0.9,
                    }
                ],
                [{"code": "propagation_compression_lossy", "severity": "warning"}],
            )

            self.assertEqual(assessment["theme_calibration"]["grade"], "A")
            self.assertEqual(assessment["pathway_graph_calibration"]["grade"], "C")
            self.assertEqual(assessment["overall_calibration"]["grade"], "B")
            self.assertEqual(assessment["grade"], "B")
            self.assertEqual(assessment["level"], "usable_with_caution")

    def test_disease_model_specific_pathway_is_appendix_exploratory(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            calibrated = service.calibrated_prediction_for_ranking_row(
                {
                    "rank": 1,
                    "pathway_uid": "pathway_1",
                    "display_name": "Glucose metabolism in triple negative breast cancer cells",
                    "score": 100.0,
                    "score_components": {},
                    "claim_refs": {"traceability_passed": True},
                },
                {
                    "pathway_uid": "pathway_1",
                    "name": "Glucose metabolism in triple negative breast cancer cells",
                    "overlap_count": 3,
                    "matched_metabolite_uids": ["met_1"],
                    "directional_support": {"up": 1, "down": 0, "unchanged": 0, "unknown": 0, "significant": 1},
                    "score_components": {},
                },
                "pathway_uid",
                None,
                {
                    "met_1": [
                        {
                            "direction": "up",
                            "significant": True,
                            "p_user": 0.9,
                            "input_name": "glucose",
                            "display_name": "Glucose",
                        }
                    ]
                },
                {"input_count": 1, "matched_count": 1, "ambiguous_count": 0, "unmatched_count": 0, "invalid_count": 0},
                {"fields": {}, "terms": [], "has_context": False},
            )

            self.assertEqual(calibrated["confidence_tier"], "exploratory")
            self.assertEqual(calibrated["calibration_status"], "appendix_low")
            self.assertTrue(calibrated["appendix"])
            self.assertIn("disease_or_model_specific_result", calibrated["downgrade_reason"])

    def test_ranking_high_confidence_requires_direction_and_significance(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            calibrated = service.calibrated_prediction_for_ranking_row(
                {
                    "rank": 1,
                    "pathway_uid": "pathway_1",
                    "display_name": "Carnitine fatty acid oxidation",
                    "score": 100.0,
                    "score_components": {},
                    "claim_refs": {"traceability_passed": True},
                },
                {
                    "pathway_uid": "pathway_1",
                    "name": "Carnitine fatty acid oxidation",
                    "overlap_count": 2,
                    "matched_metabolite_uids": ["met_1"],
                    "directional_support": {"up": 0, "down": 0, "unchanged": 1, "unknown": 0, "significant": 0},
                    "score_components": {},
                },
                "pathway_uid",
                None,
                {"met_1": [{"direction": "unchanged", "significant": False, "p_user": 0.0}]},
                {"input_count": 1, "matched_count": 1, "ambiguous_count": 0, "unmatched_count": 0, "invalid_count": 0},
                {"fields": {}, "terms": [], "has_context": False},
            )

            self.assertEqual(calibrated["confidence_tier"], "exploratory")
            self.assertIn("no_significant_support", calibrated["downgrade_reason"])
            self.assertIn("no_direction_consistency", calibrated["downgrade_reason"])

    def test_context_encoder_and_priority_themes_cover_requested_modes(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))

            def row(index, name, direction="up", padj=0.01):
                fc = 1.0 if direction == "up" else -1.0
                return {
                    "input_id": f"row_{index}",
                    "query": name,
                    "record": {"metabolite": name, "log2FC": fc, "padj": padj},
                    "metabolite_uid": f"met_extra_{index}",
                    "resolution": {
                        "top_score": 95,
                        "top_margin": 40,
                        "candidates": [{"entity_uid": f"met_extra_{index}", "display_name": name, "score": 95}],
                    },
                }

            renal_precheck = {
                "input_count": 8,
                "summary": {"matched": 8, "ambiguous": 0, "unmatched": 0, "invalid": 0},
                "matched": [
                    row(1, "indoxyl sulfate"),
                    row(2, "creatinine"),
                    row(3, "taurine"),
                    row(4, "palmitoylcarnitine"),
                    row(5, "carnitine", "down"),
                    row(6, "ADMA"),
                    row(7, "arginine", "down"),
                    row(8, "hypoxanthine"),
                ],
            }
            renal_context = service.normalize_prediction_context({"organ": "kidney", "cell_type": "renal proximal tubular epithelial cell"})
            renal_themes = {
                item["theme_id"]: item
                for item in service.build_metabolic_theme_predictions(
                    {},
                    [item["record"] for item in renal_precheck["matched"]],
                    {"input_count": 8, "matched_count": 8, "ambiguous_count": 0, "unmatched_count": 0, "invalid_count": 0},
                    renal_precheck,
                    renal_context,
                )
            }
            for theme_id in (
                "renal_organic_anion_uremic_toxin_handling",
                "renal_osmolyte_homeostasis",
                "acylcarnitine_fatty_acid_oxidation_pressure",
                "arginine_no_metabolism",
                "purine_degradation_nucleotide_stress",
            ):
                self.assertIn(theme_id, renal_themes)
                self.assertIn(renal_themes[theme_id]["context_fit"], {"context_boosted", "context_neutral"})

            glioma_precheck = {
                "input_count": 7,
                "summary": {"matched": 7, "ambiguous": 0, "unmatched": 0, "invalid": 0},
                "matched": [
                    row(11, "lactate"),
                    row(12, "glutamate"),
                    row(13, "GABA", "down"),
                    row(14, "2-hydroxyglutarate"),
                    row(15, "alpha-ketoglutarate", "down"),
                    row(16, "oxidized glutathione"),
                    row(17, "kynurenine"),
                ],
            }
            glioma_context = service.normalize_prediction_context({"disease": "glioma", "organ": "brain"})
            glioma_themes = {
                item["theme_id"]: item
                for item in service.build_metabolic_theme_predictions(
                    {},
                    [item["record"] for item in glioma_precheck["matched"]],
                    {"input_count": 7, "matched_count": 7, "ambiguous_count": 0, "unmatched_count": 0, "invalid_count": 0},
                    glioma_precheck,
                    glioma_context,
                )
            }
            for theme_id in (
                "glycolysis_glucose_utilization",
                "glial_glutamate_gaba_metabolism",
                "idh_like_metabolic_pressure",
                "glutathione_redox_stress",
                "tryptophan_kynurenine",
            ):
                self.assertIn(theme_id, glioma_themes)

            epithelial_context = service.normalize_prediction_context({"cell_type": "epithelial cell"})
            tags = {row["tag"] for row in epithelial_context["context_tags"]}
            self.assertIn("epithelial_barrier_membrane_repair", tags)
            self.assertIn("nucleotide_sugar_glycosylation", tags)
            self.assertIn("phospholipid_membrane_remodeling", tags)

    def test_special_direction_patterns_are_not_marked_mixed(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))

            def matched(index, name, direction):
                return {
                    "input_id": f"p_{index}",
                    "query": name,
                    "record": {"metabolite": name, "log2FC": 1.0 if direction == "up" else -1.0, "padj": 0.01},
                    "metabolite_uid": f"met_p_{index}",
                    "resolution": {"candidates": [{"entity_uid": f"met_p_{index}", "display_name": name, "score": 99}]},
                }

            precheck = {
                "input_count": 5,
                "summary": {"matched": 5, "ambiguous": 0, "unmatched": 0, "invalid": 0},
                "matched": [
                    matched(1, "palmitoylcarnitine", "up"),
                    matched(2, "stearoylcarnitine", "up"),
                    matched(3, "carnitine", "down"),
                    matched(4, "ADMA", "up"),
                    matched(5, "arginine", "down"),
                ],
            }
            coverage = service.build_theme_coverage(precheck, {"input_count": 5, "matched_count": 5, "ambiguous_count": 0, "unmatched_count": 0, "invalid_count": 0})
            fao = coverage["acylcarnitine_fatty_acid_oxidation_pressure"]["coverage"]
            self.assertEqual(fao["direction_consistency"], 1.0)
            self.assertEqual(fao["biological_pattern_consistency"], "fatty_acid_oxidation_pressure_consistent")

    def test_structured_prediction_json_contract_and_generalized_appendix(self):
        with tempfile.TemporaryDirectory() as tmp:
            service = self.make_service(Path(tmp))
            analyzed = service.analyze_metabolites([{"metabolite": "Glucose", "log2FC": 1.0, "padj": 0.01}], max_paths=5, max_hops=3)
            structured = analyzed["structured_prediction"]
            interpretation = analyzed["interpretation_report"]
            self.assertEqual(interpretation["report_version"], "interpretation_report.v1")
            self.assertTrue(interpretation["executive_summary"])
            self.assertTrue(interpretation["conclusion_chains"])
            self.assertTrue(interpretation["node_cards"])
            first_claim = interpretation["conclusion_chains"][0]
            first_card = interpretation["node_cards"][0]
            self.assertIn("headline", first_claim)
            self.assertIn("confidence_reasons", first_claim)
            self.assertIn("boundary", first_claim)
            self.assertTrue(first_claim["confidence_reasons"]["positive_factors"])
            self.assertIn("entity_summary", first_card)
            self.assertIn("label", first_card["entity_summary"])
            self.assertIn("node_uid", first_card["entity_summary"])
            self.assertIn("evidence_ref_count", first_card["entity_summary"])
            self.assertEqual(structured["mode"], "generalized")
            for key in (
                "context",
                "input_quality",
                "calibration",
                "high_confidence_themes",
                "medium_confidence_themes",
                "exploratory_themes",
                "downgraded_but_supported_themes",
                "pathway_evidence",
                "target_predictions",
                "disease_predictions",
                "drug_hypotheses",
                "context_mismatch_results",
                "product_layers",
                "review_queue",
                "warnings",
                "evidence_refs",
            ):
                self.assertIn(key, structured)
            product_layers = structured["product_layers"]
            for key in (
                "primary_research_candidates",
                "mechanistic_support",
                "appendix_overlay_only",
                "context_mismatch",
                "clinical_warning_or_research_only",
                "excluded_from_primary_reason",
            ):
                self.assertIn(key, product_layers)
            self.assertFalse(any(row["result_type"] in {"drug", "disease"} for row in product_layers["primary_research_candidates"]))
            for row in structured["disease_predictions"]:
                self.assertTrue(row["appendix"])
                self.assertIn("confidence_tier", row)
                self.assertIn("evidence_refs", row)
                self.assertIn("claim_refs", row)
                self.assertIn("excluded_from_primary_reason", row)
            for row in structured["drug_hypotheses"]:
                self.assertTrue(row["research_only"])
                self.assertTrue(row["appendix"])
                self.assertIn("excluded_from_primary_reason", row)
            for row in [
                *structured["high_confidence_themes"],
                *structured["medium_confidence_themes"],
                *structured["exploratory_themes"],
                *structured["downgraded_but_supported_themes"],
            ]:
                self.assertLessEqual(row["matched_support_count"], row["related_input_count"])

if __name__ == "__main__":
    unittest.main()
