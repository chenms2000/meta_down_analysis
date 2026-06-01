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


def load_script(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


build_graph_projection = load_script("build_graph_projection", ROOT / "scripts" / "build_graph_projection.py")
resolve_entity = load_script("resolve_entity", ROOT / "scripts" / "resolve_entity.py")


def write_rows(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


@unittest.skipIf(pa is None or pq is None, "pyarrow is required")
class GraphProjectionTests(unittest.TestCase):
    def test_projection_builds_sparse_graph_and_resolver(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            release = "mvp_20260101T000000"
            normalized = root / "normalized_store" / release
            output = root / "graph_projection"

            common = {"source_release": release, "license_id": "open_core:test"}
            write_rows(
                normalized / "metabolites.parquet",
                [
                    {
                        "metabolite_uid": "metabolite_1",
                        "canonical_name": "Glucose",
                        "synonyms": ["D-Glucose"],
                        "inchikey": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
                        "external_xrefs": ["CHEBI:17234", "PUBCHEM.COMPOUND:5793", "HMDB:HMDB0000122", "KEGG:C00031"],
                        "source_priority": "chebi",
                        **common,
                    }
                ],
            )
            write_rows(
                normalized / "genes.parquet",
                [
                    {
                        "gene_uid": "gene_1",
                        "ensembl_gene_id": "ENSG00000156515",
                        "entrez_gene_id": "3098",
                        "symbol": "HK1",
                        "aliases": ["Hexokinase 1"],
                        "biotype": "protein_coding",
                        "uniprot_ids": ["P19367"],
                        "external_xrefs": ["ENSEMBL:ENSG00000156515", "NCBI.GENE:3098"],
                        **common,
                    }
                ],
            )
            write_rows(
                normalized / "pathways.parquet",
                [
                    {
                        "pathway_uid": "pathway_1",
                        "name": "Glycolysis",
                        "source_name": "Reactome",
                        "primary_external_id": "R-HSA-70171",
                        "external_xrefs": ["REACTOME:R-HSA-70171"],
                        **common,
                    }
                ],
            )
            write_rows(
                normalized / "reactions.parquet",
                [
                    {
                        "reaction_uid": "reaction_1",
                        "name": "Glucose phosphorylation",
                        "species": "Homo sapiens",
                        "primary_external_id": "R-HSA-1",
                        "pathway_uid": "pathway_1",
                        "pathway_external_id": "R-HSA-70171",
                        "source_name": "Reactome",
                        "external_xrefs": ["REACTOME:R-HSA-1"],
                        **common,
                    }
                ],
            )
            write_rows(
                normalized / "diseases.parquet",
                [
                    {
                        "disease_uid": "disease_1",
                        "primary_external_id": "MONDO:0004992",
                        "name": "Cancer",
                        "aliases": ["malignant neoplasm"],
                        "external_xrefs": ["MONDO:0004992", "DOID:162"],
                        "source_priority": "mondo",
                        **common,
                    }
                ],
            )
            write_rows(
                normalized / "cell_types.parquet",
                [
                    {
                        "cell_type_uid": "cell_type_1",
                        "primary_external_id": "CL:0000066",
                        "name": "epithelial cell",
                        "aliases": ["epitheliocyte"],
                        "external_xrefs": ["CL:0000066"],
                        "source_priority": "Cell Ontology",
                        **common,
                    }
                ],
            )
            write_rows(
                normalized / "cell_states.parquet",
                [
                    {
                        "cell_state_uid": "cell_state_1",
                        "primary_external_id": "METABO_STATE:hypoxic",
                        "name": "hypoxic",
                        "aliases": ["hypoxia-associated"],
                        "external_xrefs": ["METABO_STATE:hypoxic"],
                        "state_category": "microenvironment",
                        **common,
                    }
                ],
            )
            write_rows(
                normalized / "tissues.parquet",
                [
                    {
                        "tissue_uid": "tissue_1",
                        "primary_external_id": "UBERON:0002048",
                        "name": "lung",
                        "aliases": ["pulmo"],
                        "external_xrefs": ["UBERON:0002048"],
                        "source_priority": "UBERON",
                        **common,
                    }
                ],
            )
            write_rows(
                normalized / "targets.parquet",
                [
                    {
                        "target_uid": "target_1",
                        "target_external_id": "ENSG00000156515",
                        "approved_symbol": "HK1",
                        "preferred_name": "hexokinase 1",
                        "target_type": "protein_coding",
                        "external_xrefs": ["OPENTARGETS:ENSG00000156515"],
                        **common,
                    }
                ],
            )
            for table, uid, source, xid in [
                ("metabolite_xrefs", "metabolite_1", "HMDB", "HMDB0000122"),
                ("gene_xrefs", "gene_1", "NCBI.GENE", "3098"),
                ("disease_xrefs", "disease_1", "MESH", "D009369"),
                ("cell_type_xrefs", "cell_type_1", "CL", "0000066"),
                ("cell_state_xrefs", "cell_state_1", "METABO_STATE", "hypoxic"),
                ("tissue_xrefs", "tissue_1", "UBERON", "0002048"),
                ("target_xrefs", "target_1", "ENSEMBL", "ENSG00000156515"),
            ]:
                uid_col = table.replace("_xrefs", "_uid")
                write_rows(
                    normalized / f"{table}.parquet",
                    [
                        {
                            uid_col: uid,
                            "xref_source": source,
                            "xref_id": xid,
                            "xref_key": f"{source}:{xid}",
                            **common,
                        }
                    ],
                )
            edge_common = {**common, "parser_hash": "parser", "source_name": "test"}
            write_rows(
                normalized / "metabolite_pathway_edges.parquet",
                [
                    {
                        "edge_uid": "edge_mp",
                        "metabolite_uid": "metabolite_1",
                        "pathway_uid": "pathway_1",
                        "predicate": "participates_in",
                        "source_record_id": "m|p",
                        "evidence_level": "curated",
                        "metabolite_external_id": "CHEBI:17234",
                        "pathway_external_id": "R-HSA-70171",
                        "evidence_code": "TAS",
                        "species": "Homo sapiens",
                        **edge_common,
                    }
                ],
            )
            write_rows(
                normalized / "metabolite_identity_edges.parquet",
                [
                    {
                        "edge_uid": "edge_same_as",
                        "subject_metabolite_uid": "metabolite_1",
                        "object_metabolite_uid": "metabolite_1",
                        "subject_external_id": "PUBCHEM.COMPOUND:5793",
                        "object_external_id": "CHEBI:17234",
                        "predicate": "same_as",
                        "source_record_id": "cid|chebi",
                        "evidence_level": "reviewed_identity_bridge",
                        "match_basis": "inchikey_full",
                        "linked_name": "Glucose",
                        **edge_common,
                    }
                ],
            )
            write_rows(
                normalized / "gene_pathway_edges.parquet",
                [
                    {
                        "edge_uid": "edge_gp",
                        "gene_uid": "gene_1",
                        "pathway_uid": "pathway_1",
                        "predicate": "involved_in",
                        "source_record_id": "g|p",
                        "evidence_level": "curated",
                        "gene_external_id": "ENSEMBL:ENSG00000156515",
                        "pathway_external_id": "R-HSA-70171",
                        "evidence_code": "TAS",
                        "species": "Homo sapiens",
                        **edge_common,
                    }
                ],
            )
            write_rows(
                normalized / "target_gene_edges.parquet",
                [
                    {
                        "edge_uid": "edge_tg",
                        "target_uid": "target_1",
                        "gene_uid": "gene_1",
                        "predicate": "targets",
                        "target_external_id": "ENSG00000156515",
                        "gene_external_id": "ENSEMBL:ENSG00000156515",
                        **edge_common,
                    }
                ],
            )
            write_rows(
                normalized / "target_disease_edges.parquet",
                [
                    {
                        "edge_uid": "edge_td",
                        "target_uid": "target_1",
                        "disease_uid": "disease_1",
                        "predicate": "associated_with",
                        "association_score": 0.9,
                        "evidence_count": 5,
                        "aggregation_type": "overall",
                        "aggregation_value": "None",
                        "current_novelty": 0.1,
                        "score_components_json": "{}",
                        "source_record_id": "t|d",
                        "evidence_level": "curated/inferred",
                        "target_external_id": "ENSG00000156515",
                        "disease_external_id": "MONDO:0004992",
                        **edge_common,
                    }
                ],
            )
            write_rows(
                normalized / "disease_parent_edges.parquet",
                [
                    {
                        "edge_uid": "edge_dp",
                        "child_disease_uid": "disease_1",
                        "parent_disease_uid": "disease_1",
                        "predicate": "is_a",
                        "child_external_id": "MONDO:0004992",
                        "parent_external_id": "MONDO:0004992",
                        **edge_common,
                    }
                ],
            )
            write_rows(
                normalized / "cell_type_parent_edges.parquet",
                [
                    {
                        "edge_uid": "edge_ctp",
                        "child_cell_type_uid": "cell_type_1",
                        "parent_cell_type_uid": "cell_type_1",
                        "predicate": "is_a",
                        "child_external_id": "CL:0000066",
                        "parent_external_id": "CL:0000066",
                        **edge_common,
                    }
                ],
            )
            write_rows(
                normalized / "cell_state_parent_edges.parquet",
                [
                    {
                        "edge_uid": "edge_csp",
                        "child_cell_state_uid": "cell_state_1",
                        "parent_cell_state_uid": "cell_state_1",
                        "predicate": "is_a",
                        "child_external_id": "METABO_STATE:hypoxic",
                        "parent_external_id": "METABO_STATE:hypoxic",
                        **edge_common,
                    }
                ],
            )
            write_rows(
                normalized / "tissue_parent_edges.parquet",
                [
                    {
                        "edge_uid": "edge_tp",
                        "child_tissue_uid": "tissue_1",
                        "parent_tissue_uid": "tissue_1",
                        "predicate": "is_a",
                        "child_external_id": "UBERON:0002048",
                        "parent_external_id": "UBERON:0002048",
                        **edge_common,
                    }
                ],
            )
            write_rows(
                normalized / "pathway_hierarchy_edges.parquet",
                [
                    {
                        "edge_uid": "edge_ph",
                        "parent_pathway_uid": "pathway_1",
                        "child_pathway_uid": "pathway_1",
                        "predicate": "parent_of",
                        "parent_external_id": "R-HSA-70171",
                        "child_external_id": "R-HSA-70171",
                        **edge_common,
                    }
                ],
            )
            write_rows(
                normalized / "reaction_participants.parquet",
                [
                    {
                        "edge_uid": "edge_rp",
                        "reaction_uid": "reaction_1",
                        "participant_uid": "metabolite_1",
                        "participant_external_id": "CHEBI:17234",
                        "participant_type": "metabolite",
                        "role": "input",
                        "source_record_id": "r|m",
                        **edge_common,
                    }
                ],
            )

            manifest_path = build_graph_projection.build_graph_projection(root / "normalized_store", output, release)
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["node_count"], 9)
            self.assertGreaterEqual(manifest["edge_type_count"], 12)

            nodes = pq.read_table(output / release / "nodes.parquet").to_pylist()
            sparse = pq.read_table(output / release / "sparse_edges.parquet").to_pylist()
            self.assertEqual(len(nodes), 9)
            self.assertEqual(len(sparse), 12)

            resolver_path = output / release / "resolver_index.parquet"
            hmdb = resolve_entity.resolve_entities(resolver_path, "HMDB0000122", entity_type="metabolite")
            self.assertEqual(hmdb[0]["entity_uid"], "metabolite_1")
            gene_alias = resolve_entity.resolve_entities(resolver_path, "Hexokinase 1", entity_type="gene")
            self.assertEqual(gene_alias[0]["entity_uid"], "gene_1")
            cell_type = resolve_entity.resolve_entities(resolver_path, "CL:0000066", entity_type="cell_type")
            self.assertEqual(cell_type[0]["entity_uid"], "cell_type_1")
            tissue = resolve_entity.resolve_entities(resolver_path, "UBERON:0002048", entity_type="tissue")
            self.assertEqual(tissue[0]["entity_uid"], "tissue_1")
            cell_state = resolve_entity.resolve_entities(resolver_path, "hypoxia-associated", entity_type="cell_state")
            self.assertEqual(cell_state[0]["entity_uid"], "cell_state_1")

    def test_candidate_queries_include_pubchem_and_disease_variants(self):
        candidates = resolve_entity.candidate_queries("PubChem:5793")
        self.assertIn(("PUBCHEM.COMPOUND", "5793", 120), candidates)
        candidates = resolve_entity.candidate_queries("EFO_0000311")
        self.assertIn(("OPENTARGETS_DISEASE", "0000311", 120), candidates)
        candidates = resolve_entity.candidate_queries("NCIT:C3512")
        self.assertIn(("OPENTARGETS_DISEASE", "c3512", 120), candidates)
        candidates = resolve_entity.candidate_queries("ONCOTREE:LUAD")
        self.assertIn(("OPENTARGETS_DISEASE", "luad", 120), candidates)
        candidates = resolve_entity.candidate_queries("P05231")
        self.assertIn(("UNIPROT", "p05231", 100), candidates)


if __name__ == "__main__":
    unittest.main()
