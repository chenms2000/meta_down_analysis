import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover
    pa = None
    pq = None


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_literature_evidence.py"
SPEC = importlib.util.spec_from_file_location("build_literature_evidence", SCRIPT)
builder = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = builder
SPEC.loader.exec_module(builder)


def write_parquet(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


@unittest.skipIf(pa is None or pq is None, "pyarrow is required")
class LiteratureEvidenceBuilderTests(unittest.TestCase):
    def make_normalized_store(self, root: Path) -> str:
        release = "mvp_20260101T000000"
        normalized = root / "normalized_store" / release
        common = {"source_release": release, "license_id": "open_core:test", "parser_hash": "parser"}

        write_parquet(
            normalized / "sentences.parquet",
            [
                {
                    "sentence_uid": "sent_1",
                    "article_uid": "article_1",
                    "pmid": "1",
                    "pmcid": "PMC1",
                    "section": "abstract",
                    "sentence_text": "Lactate is elevated in hepatocellular carcinoma and is associated with poor prognosis.",
                    "source_release": release,
                    "parser_hash": "parser",
                },
                {
                    "sentence_uid": "sent_2",
                    "article_uid": "article_1",
                    "pmid": "1",
                    "pmcid": "PMC1",
                    "section": "abstract",
                    "sentence_text": "MYC is associated with hepatocellular carcinoma and MYC activates glycolysis.",
                    "source_release": release,
                    "parser_hash": "parser",
                },
                {
                    "sentence_uid": "sent_3",
                    "article_uid": "article_2",
                    "pmid": "2",
                    "pmcid": "PMC2",
                    "section": "abstract",
                    "sentence_text": "Water was measured in cell culture medium.",
                    "source_release": release,
                    "parser_hash": "parser",
                },
                {
                    "sentence_uid": "sent_4",
                    "article_uid": "article_2",
                    "pmid": "2",
                    "pmcid": "PMC2",
                    "section": "abstract",
                    "sentence_text": "Protein kinase is associated with hepatocellular carcinoma.",
                    "source_release": release,
                    "parser_hash": "parser",
                },
            ],
        )
        write_parquet(
            normalized / "articles.parquet",
            [
                {
                    "article_uid": "article_1",
                    "pmid": "1",
                    "pmcid": "PMC1",
                    "doi": "",
                    "title": "Metabolism in HCC",
                    "abstract": "",
                    "journal": "",
                    "pub_date": "",
                    "article_type": "",
                    "source_path": "",
                    "source_bytes": 0,
                    "checksum": "",
                    **common,
                }
            ],
        )
        write_parquet(
            normalized / "metabolites.parquet",
            [
                {
                    "metabolite_uid": "met_lactate",
                    "canonical_name": "L-lactate",
                    "synonyms": ["lactate"],
                    "formula": "",
                    "exact_mass": None,
                    "charge": 0,
                    "inchikey": "",
                    "smiles": "",
                    "external_xrefs": [],
                    "source_priority": "test",
                    "checksum": "",
                    **common,
                }
            ],
        )
        write_parquet(
            normalized / "genes.parquet",
            [
                {
                    "gene_uid": "gene_myc",
                    "ensembl_gene_id": "ENSG_MYC",
                    "entrez_gene_id": "",
                    "symbol": "MYC",
                    "aliases": ["c-Myc"],
                    "description": "",
                    "taxon": "9606",
                    "biotype": "protein_coding",
                    "chromosome": "",
                    "start": 0,
                    "end": 0,
                    "strand": "",
                    "uniprot_ids": [],
                    "external_xrefs": [],
                    "checksum": "",
                    **common,
                },
                {
                    "gene_uid": "gene_pk",
                    "ensembl_gene_id": "ENSG_PK",
                    "entrez_gene_id": "",
                    "symbol": "PKX1",
                    "aliases": ["protein kinase"],
                    "description": "",
                    "taxon": "9606",
                    "biotype": "protein_coding",
                    "chromosome": "",
                    "start": 0,
                    "end": 0,
                    "strand": "",
                    "uniprot_ids": [],
                    "external_xrefs": [],
                    "checksum": "",
                    **common,
                },
            ],
        )
        write_parquet(
            normalized / "targets.parquet",
            [
                {
                    "target_uid": "target_myc",
                    "target_external_id": "ENSG_MYC",
                    "preferred_name": "MYC proto-oncogene",
                    "approved_symbol": "MYC",
                    "target_type": "protein_coding",
                    "gene_uid": "gene_myc",
                    "tractability_flags_json": "{}",
                    "external_xrefs": [],
                    "checksum": "",
                    **common,
                }
            ],
        )
        write_parquet(
            normalized / "diseases.parquet",
            [
                {
                    "disease_uid": "disease_hcc",
                    "primary_external_id": "MONDO:HCC",
                    "name": "hepatocellular carcinoma",
                    "aliases": ["HCC"],
                    "description": "",
                    "parents": [],
                    "external_xrefs": [],
                    "source_priority": "test",
                    "checksum": "",
                    **common,
                }
            ],
        )
        write_parquet(
            normalized / "pathways.parquet",
            [
                {
                    "pathway_uid": "path_glycolysis",
                    "name": "glycolysis",
                    "species": "Homo sapiens",
                    "source_name": "Reactome",
                    "primary_external_id": "R-HSA-70171",
                    "hierarchy_path": "",
                    "external_xrefs": [],
                    "checksum": "",
                    **common,
                }
            ],
        )
        write_parquet(
            normalized / "metabolite_pathway_edges.parquet",
            [{"edge_uid": "edge_mp", "metabolite_uid": "met_lactate", "pathway_uid": "path_glycolysis", **common}],
        )
        write_parquet(
            normalized / "gene_pathway_edges.parquet",
            [
                {"edge_uid": "edge_gp", "gene_uid": "gene_myc", "pathway_uid": "path_glycolysis", **common},
                {"edge_uid": "edge_pk", "gene_uid": "gene_pk", "pathway_uid": "path_glycolysis", **common},
            ],
        )
        write_parquet(
            normalized / "target_disease_edges.parquet",
            [{"edge_uid": "edge_td", "target_uid": "target_myc", "disease_uid": "disease_hcc", **common}],
        )
        (normalized / "normalized_manifest.json").write_text(json.dumps({"release_id": release}), encoding="utf-8")
        (root / "articles_collect").mkdir()
        (root / "articles_collect" / "PMC1.txt").write_text("PMCID: PMC1\nPMID: 1\n", encoding="utf-8")
        return release

    def test_builds_mentions_candidates_and_literature_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            release = self.make_normalized_store(root)
            args = builder.parse_args(["--workspace", str(root), "--release-id", release])
            manifest = builder.build_evidence(args)
            evidence_dir = root / "literature_evidence" / release

            self.assertEqual(manifest["metrics"]["sentence_count_scanned"], 4)
            self.assertGreaterEqual(manifest["metrics"]["sentence_mention_count"], 5)
            self.assertEqual(manifest["metrics"]["relevance_scored_sentence_count"], 0)
            self.assertGreaterEqual(manifest["metrics"]["evidence_candidate_count"], 3)
            self.assertGreaterEqual(manifest["metrics"]["supported_existing_edge_count"], 1)
            self.assertGreaterEqual(manifest["metrics"]["novel_candidate_count"], 1)
            self.assertGreaterEqual(manifest["precision_filters"]["precision_blocked_surface_count"], 1)
            self.assertIn("gene:alias:protein kinase", manifest["precision_filters"]["blocklist_hits_by_surface"])

            candidates = pq.read_table(evidence_dir / "relation_candidates.parquet").to_pylist()
            predicates = {row["predicate"] for row in candidates}
            self.assertIn("metabolite_changed_in_cancer", predicates)
            self.assertIn("target_associated_with_disease", predicates)
            self.assertIn("gene_regulates_metabolic_process", predicates)
            mentions = pq.read_table(evidence_dir / "sentence_mentions.parquet").to_pylist()
            self.assertNotIn("protein kinase", {row["normalized_surface"] for row in mentions})
            relevance = pq.read_table(evidence_dir / "sentence_relevance.parquet").to_pylist()
            self.assertEqual(relevance, [])

            support = pq.read_table(evidence_dir / "literature_edge_support.parquet").to_pylist()
            target_support = [
                row
                for row in support
                if row["subject_uid"] == "target_myc"
                and row["predicate"] == "target_associated_with_disease"
                and row["object_uid"] == "disease_hcc"
            ]
            self.assertTrue(target_support)
            self.assertIn("edge_td", target_support[0]["supported_existing_edge_uids"])
            self.assertGreater(target_support[0]["p_literature"], 0.0)
            self.assertEqual(target_support[0]["support_status_set"], ["support"])

    def test_sentence_assertion_status_separates_support_background_uncertain_and_contradict(self):
        self.assertEqual(
            builder.sentence_assertion_status("Lactate is elevated in hepatocellular carcinoma.", "abstract"),
            "support",
        )
        self.assertEqual(
            builder.sentence_assertion_status("Lactate was elevated compared with adjacent normal tissue.", "abstract"),
            "support",
        )
        self.assertEqual(
            builder.sentence_assertion_status("This review summarizes lactate metabolism in cancer.", "abstract"),
            "background",
        )
        self.assertEqual(
            builder.sentence_assertion_status("Lactate may be associated with hepatocellular carcinoma.", "abstract"),
            "uncertain",
        )
        self.assertEqual(
            builder.sentence_assertion_status("Although lactate is associated with cancer, causality remains to be determined.", "abstract"),
            "uncertain",
        )
        self.assertEqual(
            builder.sentence_assertion_status("Lactate was not associated with hepatocellular carcinoma.", "abstract"),
            "contradict",
        )
        self.assertIn(
            "causal_weakening",
            builder.semantic_cues("Although lactate is associated with cancer, causality remains to be determined.", "abstract"),
        )
        self.assertIn(
            "comparison_context",
            builder.semantic_cues("Lactate was elevated compared with adjacent normal tissue.", "abstract"),
        )

    def test_complex_semantic_cues_do_not_promote_weak_or_negated_claims(self):
        cases = [
            (
                "Lactate showed no significant difference between tumor and adjacent tissue.",
                "abstract",
                "contradict",
                {"null_result", "negation", "comparison_context"},
            ),
            (
                "Lactate trended higher in tumors, but the difference was not statistically significant.",
                "results",
                "uncertain",
                {"weak_observation", "concession"},
            ),
            (
                "LDHA expression was not necessarily associated with lactate accumulation in all tumors.",
                "discussion",
                "uncertain",
                {"hedged_negation"},
            ),
            (
                "Glucose was not only increased but also quantified by LC-MS in tumor tissue.",
                "results",
                "support",
                {"non_negating_negation", "direct_assay"},
            ),
            (
                "Previous studies have implicated glutamine metabolism in cancer.",
                "introduction",
                "background",
                {"background"},
            ),
            (
                "Although lactate was measured by LC-MS, causality remains to be determined.",
                "results",
                "uncertain",
                {"causal_weakening", "direct_assay"},
            ),
            (
                "Reports of glutamine dependence in melanoma are conflicting across cohorts.",
                "discussion",
                "uncertain",
                {"conflict"},
            ),
            (
                "LDHA increased lactate only in cell lines but not patient tumors.",
                "results",
                "uncertain",
                {"context_boundary"},
            ),
            (
                "The association between serine abundance and survival was no longer significant after adjustment for proliferation.",
                "results",
                "uncertain",
                {"causal_weakening"},
            ),
            (
                "Glutamine depletion had no detectable effect on tumor growth.",
                "results",
                "contradict",
                {"null_result", "negation"},
            ),
            (
                "Hypoxia rather than lactate explained the observed survival association.",
                "discussion",
                "uncertain",
                {"alternative_explanation", "causal_weakening"},
            ),
            (
                "LDHA inhibition reduced lactate in cell lines; however, this effect was absent in patient tumors.",
                "results",
                "uncertain",
                {"context_boundary", "concession"},
            ),
            (
                "Adding lactate failed to rescue proliferation after LDHA inhibition.",
                "results",
                "contradict",
                {"null_result", "negation"},
            ),
            (
                "Lactate levels were variable and did not consistently differ across cohorts.",
                "results",
                "uncertain",
                {"hedged_negation", "negation"},
            ),
            (
                "The survival association was not attributable to lactate after adjustment for hypoxia.",
                "results",
                "uncertain",
                {"alternative_explanation", "causal_weakening", "negation"},
            ),
            (
                "Despite increased lactate, isotope tracing did not demonstrate increased glycolytic flux.",
                "results",
                "contradict",
                {"concession", "direct_assay", "negation", "null_result"},
            ),
            (
                "Serine abundance increased, but this was insufficient to drive nucleotide synthesis without folate availability.",
                "discussion",
                "uncertain",
                {"causal_weakening", "concession", "negation"},
            ),
        ]

        for text, section, expected_status, expected_cues in cases:
            with self.subTest(text=text):
                semantics = builder.assertion_semantics(text, section)
                self.assertEqual(semantics["support_status"], expected_status)
                self.assertTrue(expected_cues <= set(semantics["semantic_cues"]))

    def test_assertion_semantics_records_method_cues(self):
        semantics = builder.assertion_semantics(
            "Patient tumor tissue lactate was quantified by LC-MS after LDHA inhibitor treatment.",
            "results",
        )

        self.assertEqual(semantics["support_status"], "support")
        self.assertIn("patient_sample", semantics["method_cues"])
        self.assertIn("measurement_assay", semantics["method_cues"])
        self.assertIn("intervention_assay", semantics["method_cues"])

    def test_assertion_semantics_records_expanded_method_cues(self):
        cases = [
            ("Reports were conflicting across cohorts.", {"patient_sample"}),
            ("Lactate was not significantly higher in tumors than adjacent tissue.", {"patient_sample"}),
            ("Spatial metabolomics quantified lactate gradients in tumor regions compared with matched normal tissue.", {"patient_sample", "measurement_assay"}),
            ("Despite increased lactate, isotope tracing did not demonstrate increased glycolytic flux.", {"measurement_assay"}),
            ("Serine withdrawal impaired nucleotide synthesis except in tumors with folate rescue.", {"intervention_assay"}),
            ("Glutamine deprivation reduced proliferation in vitro but not in vivo.", {"cell_line_model", "animal_model", "intervention_assay"}),
        ]

        for text, expected_cues in cases:
            with self.subTest(text=text):
                semantics = builder.assertion_semantics(text, "results")
                self.assertTrue(expected_cues <= set(semantics["method_cues"]))

    def test_background_relations_do_not_confirm_existing_edges(self):
        with tempfile.TemporaryDirectory() as tmp:
            normalized = Path(tmp) / "normalized_store" / "mvp_20260101T000000"
            normalized.mkdir(parents=True)
            write_parquet(
                normalized / "target_disease_edges.parquet",
                [{"edge_uid": "edge_td", "target_uid": "target_myc", "disease_uid": "disease_hcc"}],
            )
            relation_rows = [
                {
                    "relation_uid": "rel_background",
                    "subject_uid": "target_myc",
                    "subject_type": "target",
                    "predicate": "target_associated_with_disease",
                    "object_uid": "disease_hcc",
                    "object_type": "disease",
                    "polarity": "association",
                    "support_status": "background",
                    "sentence_uid": "sent_1",
                    "pmid": "100",
                    "pmcid": "PMC100",
                    "calibrated_prob": 0.7,
                    "raw_score": 0.7,
                    "license_id": "local:test",
                    "source_release": "test",
                }
            ]
            support = builder.aggregate_support(relation_rows, normalized, 0.25, "parser", "config")
            self.assertEqual(len(support), 1)
            self.assertEqual(support[0]["support_class"], "background_candidate")
            self.assertEqual(support[0]["support_status_set"], ["background"])

    def test_support_probability_clusters_duplicate_sentences_by_article(self):
        with tempfile.TemporaryDirectory() as tmp:
            normalized = Path(tmp) / "normalized_store" / "mvp_20260101T000000"
            normalized.mkdir(parents=True)
            relation_rows = [
                {
                    "relation_uid": "rel_1",
                    "subject_uid": "met_lactate",
                    "subject_type": "metabolite",
                    "predicate": "metabolite_changed_in_cancer",
                    "object_uid": "disease_hcc",
                    "object_type": "disease",
                    "polarity": "increase",
                    "sentence_uid": "sent_1",
                    "pmid": "100",
                    "pmcid": "PMC100",
                    "calibrated_prob": 0.4,
                    "raw_score": 0.7,
                    "license_id": "local:test",
                    "source_release": "test",
                },
                {
                    "relation_uid": "rel_2",
                    "subject_uid": "met_lactate",
                    "subject_type": "metabolite",
                    "predicate": "metabolite_changed_in_cancer",
                    "object_uid": "disease_hcc",
                    "object_type": "disease",
                    "polarity": "increase",
                    "sentence_uid": "sent_2",
                    "pmid": "100",
                    "pmcid": "PMC100",
                    "calibrated_prob": 0.4,
                    "raw_score": 0.7,
                    "license_id": "local:test",
                    "source_release": "test",
                },
                {
                    "relation_uid": "rel_3",
                    "subject_uid": "met_lactate",
                    "subject_type": "metabolite",
                    "predicate": "metabolite_changed_in_cancer",
                    "object_uid": "disease_hcc",
                    "object_type": "disease",
                    "polarity": "increase",
                    "sentence_uid": "sent_3",
                    "pmid": "200",
                    "pmcid": "PMC200",
                    "calibrated_prob": 0.4,
                    "raw_score": 0.7,
                    "license_id": "local:test",
                    "source_release": "test",
                },
            ]
            support = builder.aggregate_support(relation_rows, normalized, 0.25, "parser", "config")
            self.assertEqual(len(support), 1)
            self.assertAlmostEqual(support[0]["p_literature"], 0.64)
            components = json.loads(support[0]["score_components_json"])
            self.assertEqual(components["article_count_for_probability"], 2)
            self.assertEqual(components["article_probabilities"], [0.4, 0.4])

    def test_low_precision_disease_filter_blocks_ncit_non_disease_surfaces(self):
        self.assertTrue(builder.is_low_precision_disease_surface("KRAS", "Oncogene K-Ras", "alias"))
        self.assertTrue(builder.is_low_precision_disease_surface("Targeting", "Targeted Therapy Agent", "alias"))
        self.assertTrue(builder.is_low_precision_disease_surface("Tumor Necrosis Factor", "Tumor Necrosis Factor", "name"))
        self.assertTrue(builder.is_low_precision_disease_surface("BRCA1", "BRCA1 Gene", "alias"))
        self.assertTrue(builder.is_low_precision_disease_surface("positive", "KIT Positive", "alias"))
        self.assertTrue(builder.is_low_precision_disease_surface("negative", "KIT Negative", "alias"))
        self.assertTrue(builder.is_low_precision_disease_surface("tumor suppressor", "regulation of cell cycle", "alias"))
        self.assertFalse(
            builder.is_low_precision_disease_surface(
                "TNBC",
                "Triple-Negative Breast Cancer Finding",
                "alias",
            )
        )
        self.assertFalse(
            builder.is_low_precision_disease_surface(
                "non-small cell lung cancer",
                "non-small cell lung carcinoma",
                "name",
            )
        )

    def test_relevance_mode_dependency_error_is_explicit_and_opt_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            release = self.make_normalized_store(root)
            args = builder.parse_args(["--workspace", str(root), "--release-id", release])
            builder.build_evidence(args)
            pubmedbert_args = builder.parse_args(["--workspace", str(root), "--release-id", release, "--relevance-mode", "pubmedbert"])
            with mock.patch.object(builder, "PubMedBertRelevanceScorer", side_effect=RuntimeError("missing torch")):
                with self.assertRaisesRegex(RuntimeError, "missing torch"):
                    builder.build_evidence(pubmedbert_args)

    def test_relevance_weight_only_downweights_existing_relation_probability(self):
        row = {
            "relation_uid": "rel_1",
            "calibrated_prob": 0.8,
            "score_components_json": builder.stable_json({"formula": "rule"}),
        }
        weighted = builder.apply_relevance_weight(row, 0.5, 0.2)
        self.assertIsNotNone(weighted)
        assert weighted is not None
        self.assertAlmostEqual(weighted["calibrated_prob"], 0.6)
        components = json.loads(weighted["score_components_json"])
        self.assertEqual(components["relevance_weight"]["formula"], "rule_relation_prob*(0.5+0.5*relevance_score)")
        self.assertIsNone(builder.apply_relevance_weight(row, 0.0, 0.5))


if __name__ == "__main__":
    unittest.main()
