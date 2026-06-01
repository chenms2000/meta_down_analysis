import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover
    pa = None
    pq = None


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))


def load_script(name: str):
    script = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


build_learning_views = load_script("build_learning_views")
build_weak_labels = load_script("build_weak_labels")
train_priority_ranker = load_script("train_priority_ranker")
build_cell_context_overlay_seed = load_script("build_cell_context_overlay_seed")
build_prediction_overlays = load_script("build_prediction_overlays")


def write_parquet(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


@unittest.skipIf(pa is None or pq is None, "pyarrow is required")
class LearningTemporalHoldoutTests(unittest.TestCase):
    def test_publication_years_are_joined_from_article_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            normalized = Path(tmp) / "normalized_store" / "mvp_test"
            write_parquet(
                normalized / "articles.parquet",
                [
                    {"article_uid": "a1", "pmid": "100", "pmcid": "PMC100", "pub_date": "2023 Mar"},
                    {"article_uid": "a2", "pmid": "200", "pmcid": "PMC200", "pub_date": "2025 Jan"},
                ],
            )
            lookup = build_learning_views.load_article_year_lookup(normalized)
            relations = pd.DataFrame(
                [
                    {"support_uid": "s1", "pmids": ["100", "200"], "pmcids": ["PMC100", "PMC200"]},
                    {"support_uid": "s2", "pmids": ["missing"], "pmcids": []},
                ]
            )

            enriched = build_learning_views.add_publication_years(relations, lookup)

            self.assertEqual(enriched.loc[0, "publication_years"], [2023, 2025])
            self.assertEqual(enriched.loc[0, "min_publication_year"], 2023)
            self.assertEqual(enriched.loc[0, "max_publication_year"], 2025)
            self.assertEqual(enriched.loc[0, "publication_year"], 2025)
            self.assertEqual(enriched.loc[1, "publication_year_source"], "")

    def test_weak_labels_use_year_and_keep_negative_split_with_source_paper(self):
        entity_features = {
            "target_a": {"display_name": "A"},
            "target_b": {"display_name": "B"},
            "met_1": {"display_name": "Metabolite"},
        }
        positive = build_weak_labels.feature_row(
            candidate_uid="pos_1",
            subject_uid="met_1",
            subject_type="metabolite",
            object_uid="target_a",
            object_type="target",
            predicate="metabolite_regulates_gene",
            source_kind="literature_relation",
            entity_features=entity_features,
            literature={"pmids": ["200"], "publication_year": 2025},
            weak_label="weak_positive",
            label_score=0.58,
            label_reason="novel_literature_candidate",
        )

        negatives = build_weak_labels.build_negative_candidates(
            positives=[positive],
            by_type={"target": ["target_b", "target_c"]},
            entity_features=entity_features,
            negatives_per_positive=1.0,
            seed=13,
        )

        self.assertEqual(positive["paper_id"], "200")
        self.assertEqual(positive["split"], "temporal_holdout")
        self.assertTrue(positive["temporal_holdout"])
        self.assertEqual(negatives[0]["paper_id"], "200")
        self.assertEqual(negatives[0]["split"], "temporal_holdout")
        self.assertTrue(negatives[0]["temporal_holdout"])

    def test_pathway_appendix_keeps_raw_score_and_caps_display_score(self):
        predictions = pd.DataFrame(
            [
                {
                    "candidate_uid": "c1",
                    "subject_uid": "met_1",
                    "subject_type": "metabolite",
                    "subject_name": "Metabolite",
                    "predicate": "associated_with",
                    "object_uid": "path_melanoma",
                    "object_type": "pathway",
                    "object_name": "Melanoma",
                    "priority_score": 0.95,
                    "priority_score_raw": 0.99,
                    "label_score": 1.0,
                    "source_kind": "literature_relation",
                    "label_reason": "curated_and_literature",
                }
            ]
        )

        pathways = train_priority_ranker.aggregate_priorities(predictions, "pathway")
        row = pathways.iloc[0]

        self.assertTrue(row["appendix"])
        self.assertEqual(row["display_tier"], "appendix")
        self.assertEqual(row["downgrade_reason"], "disease_or_model_specific_pathway_appendix_in_generalized_report")
        self.assertEqual(row["priority_score"], 0.25)
        self.assertEqual(row["priority_score_raw"], 0.99)

    def test_no_leakage_feature_mode_excludes_source_shortcuts(self):
        frame = pd.DataFrame(
            [
                {
                    "candidate_uid": "c1",
                    "subject_type": "drug",
                    "object_type": "target",
                    "predicate": "targets",
                    "source_kind": "drug_target_overlay",
                    "label_reason": "curated_drug_target_overlay",
                    "graph_degree_sum_log": 1.0,
                    "entity_semantic_overlap": 0.25,
                    "metabolic_theme_overlap": 1.0,
                },
                {
                    "candidate_uid": "c2",
                    "subject_type": "metabolite",
                    "object_type": "disease",
                    "predicate": "metabolite_changed_in_cancer",
                    "source_kind": "literature_relation",
                    "label_reason": "literature_support",
                    "graph_degree_sum_log": 2.0,
                    "entity_semantic_overlap": 0.0,
                    "metabolic_theme_overlap": 0.0,
                },
            ]
        )

        full_matrix, full_columns = train_priority_ranker.build_feature_matrix(frame, feature_mode="full")
        no_leakage_matrix, no_leakage_columns = train_priority_ranker.build_feature_matrix(frame, feature_mode="no_leakage")

        self.assertEqual(len(full_matrix), 2)
        self.assertEqual(len(no_leakage_matrix), 2)
        self.assertTrue(any(column.startswith("source_kind_") for column in full_columns))
        self.assertTrue(any(column.startswith("label_reason_") for column in full_columns))
        self.assertFalse(any(column.startswith("source_kind_") for column in no_leakage_columns))
        self.assertFalse(any(column.startswith("label_reason_") for column in no_leakage_columns))
        self.assertIn("graph_degree_sum_log", no_leakage_columns)
        self.assertIn("entity_semantic_overlap", no_leakage_columns)
        self.assertIn("metabolic_theme_overlap", no_leakage_columns)
        self.assertIn("shared_pathway_count", no_leakage_columns)
        self.assertIn("source_independent_prior_score", no_leakage_columns)

    def test_feature_row_adds_source_free_generalization_features(self):
        entity_features = {
            "met_lactate": {
                "display_name": "Lactate",
                "canonical_name": "lactate",
                "document": "lactate glycolysis cell cancer",
                "graph_degree": 4,
                "canonical_edge_count": 2,
                "pathway_context_count": 1,
                "pathway_context_uids": "path_gly",
            },
            "path_gly": {
                "display_name": "Glycolysis pathway",
                "canonical_name": "glycolysis",
                "document": "glycolysis lactate pyruvate cancer cell",
                "graph_degree": 6,
                "canonical_edge_count": 3,
                "pathway_context_count": 1,
                "pathway_context_uids": "path_gly",
            },
        }

        row = build_weak_labels.feature_row(
            candidate_uid="c1",
            subject_uid="met_lactate",
            subject_type="metabolite",
            object_uid="path_gly",
            object_type="pathway",
            predicate="participates_in",
            source_kind="literature_relation",
            entity_features=entity_features,
            literature={"supported_existing_edge_count": 1, "context_terms": "cancer cell"},
            weak_label="weak_positive",
            label_score=0.58,
            label_reason="literature_support",
        )

        self.assertGreater(row["entity_semantic_overlap"], 0.0)
        self.assertGreater(row["metabolic_theme_overlap"], 0.0)
        self.assertGreater(row["cell_context_overlap"], 0.0)
        self.assertEqual(row["shared_pathway_count"], 1.0)
        self.assertGreater(row["source_independent_prior_score"], 0.0)
        self.assertEqual(row["direct_graph_support_indicator"], 1.0)

    def test_source_balanced_weights_reduce_large_group_dominance(self):
        rows = []
        for index in range(100):
            rows.append({"candidate_uid": f"drug_{index}", "source_kind": "drug_target_overlay", "weak_label": "weak_positive"})
        for index in range(4):
            rows.append({"candidate_uid": f"lit_{index}", "source_kind": "literature_relation", "weak_label": "strong_positive"})
        frame = pd.DataFrame(rows)

        weights, diagnostics = train_priority_ranker.training_sample_weights(frame, source_balanced=True)
        frame["weight"] = weights
        large_group_weight = frame[frame["source_kind"] == "drug_target_overlay"]["weight"].mean()
        small_group_weight = frame[frame["source_kind"] == "literature_relation"]["weight"].mean()

        self.assertTrue(diagnostics["source_balanced"])
        self.assertGreater(small_group_weight, large_group_weight)

    def test_cell_context_seed_rows_include_target_uids_and_research_provenance(self):
        rows, unresolved = build_cell_context_overlay_seed.build_rows(
            {
                "EPCAM": "target_epcam",
                "KRT8": "target_krt8",
                "SLC22A6": "target_slc22a6",
                "SLC22A8": "target_slc22a8",
                "IDH1": "target_idh1",
                "GLUL": "target_glul",
            }
        )

        self.assertEqual(len(rows), 6)
        self.assertTrue(any(row["target_uids"] for row in rows))
        self.assertTrue(all(row["source_name"] == "local_seed_cell_context_overlay" for row in rows))
        self.assertTrue(all(row["evidence_level"] == "research_context_seed" for row in rows))
        self.assertIn("seed:renal_proximal_tubular_context", unresolved)

    def test_real_signature_imports_normalize_depmap_and_single_cell_markers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            depmap_model = root / "depmap_models.csv"
            depmap_effect = root / "depmap_gene_effect.csv"
            markers = root / "markers.csv"
            depmap_model.write_text(
                "ModelID,CellLineName,OncotreeLineage,Tissue,CellState\nACH-1,TestLine,Colon,large intestine,epithelial\n",
                encoding="utf-8",
            )
            depmap_effect.write_text(
                "ModelID,HK1 (3098),LDHA (3939),RPLP0 (6175)\nACH-1,-0.82,-0.61,-0.1\n",
                encoding="utf-8",
            )
            markers.write_text(
                "\n".join(
                    [
                        "cell_type_id,cell_type_name,gene,tissue,disease,dataset_id,avg_log2FC",
                        "ct1,Epithelial cell,EPCAM,colon,cancer,ds1,2.1",
                        "ct1,Epithelial cell,KRT8,colon,cancer,ds1,1.7",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            target_map = {"HK1": "target_hk1", "LDHA": "target_ldha", "EPCAM": "target_epcam"}

            depmap_rows = build_prediction_overlays.build_depmap_rows(
                depmap_model,
                gene_effect_path=depmap_effect,
                target_uid_by_symbol=target_map,
                top_genes_per_model=5,
                dependency_threshold=-0.5,
            )
            marker_rows = build_prediction_overlays.build_marker_signature_rows(
                markers,
                source_name="single_cell_marker",
                target_uid_by_symbol=target_map,
            )

            self.assertEqual(depmap_rows[0]["evidence_level"], "gene_dependency_signature")
            self.assertEqual(depmap_rows[0]["target_symbols"], "HK1;LDHA")
            self.assertEqual(depmap_rows[0]["target_uids"], "target_hk1;target_ldha")
            self.assertEqual(marker_rows[0]["target_symbols"], "EPCAM;KRT8")
            self.assertIn("target_epcam", marker_rows[0]["target_uids"])
            self.assertEqual(marker_rows[0]["source_name"], "single_cell_marker")


if __name__ == "__main__":
    unittest.main()
