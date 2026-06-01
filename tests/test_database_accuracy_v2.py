import importlib.util
import io
import json
import sys
import tarfile
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
SCRIPT = ROOT / "scripts" / "database_accuracy_v2.py"
SPEC = importlib.util.spec_from_file_location("database_accuracy_v2", SCRIPT)
database_accuracy_v2 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = database_accuracy_v2
SPEC.loader.exec_module(database_accuracy_v2)


def write_parquet(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


def write_tsv_tar(path: Path, files: dict[str, str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as archive:
        for name, content in files.items():
            data = content.encode("utf-8")
            info = tarfile.TarInfo(name)
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))


@unittest.skipIf(pa is None or pq is None, "pyarrow is required")
class DatabaseAccuracyV2Tests(unittest.TestCase):
    def test_v2_store_separates_traits_classes_and_reaction_facts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            release = "mvp_20260101T000000"
            normalized = root / "normalized_store" / release
            common = {"source_release": release, "license_id": "open_core:test", "parser_hash": "parser"}
            write_parquet(
                normalized / "metabolites.parquet",
                [
                    {
                        "metabolite_uid": "met_glucose",
                        "canonical_name": "Glucose",
                        "synonyms": ["D-Glucose"],
                        "formula": "C6H12O6",
                        "exact_mass": 180.063,
                        "charge": 0,
                        "inchikey": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
                        "smiles": "",
                        "external_xrefs": ["CHEBI:17234"],
                        "source_priority": "ChEBI",
                        "checksum": "",
                        **common,
                    },
                    {
                        "metabolite_uid": "met_pyruvate",
                        "canonical_name": "Pyruvate",
                        "synonyms": ["pyruvic acid"],
                        "formula": "C3H3O3",
                        "exact_mass": 87.008,
                        "charge": -1,
                        "inchikey": "LCTONWCANYUPML-UHFFFAOYSA-M",
                        "smiles": "CC(=O)[O-]",
                        "external_xrefs": ["CHEBI:15361"],
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
                        "metabolite_uid": "met_glucose",
                        "xref_source": "CHEBI",
                        "xref_id": "17234",
                        "xref_key": "CHEBI:17234",
                        "source_name": "ChEBI",
                        **common,
                    },
                    {
                        "xref_uid": "xref_2",
                        "metabolite_uid": "met_pyruvate",
                        "xref_source": "CHEBI",
                        "xref_id": "15361",
                        "xref_key": "CHEBI:15361",
                        "source_name": "ChEBI",
                        **common,
                    }
                ],
            )
            write_parquet(
                normalized / "pathways.parquet",
                [
                    {
                        "pathway_uid": "path_glycolysis",
                        "name": "Glycolysis",
                        "species": "Homo sapiens",
                        "source_name": "Reactome",
                        "primary_external_id": "R-HSA-70171",
                        "hierarchy_path": "",
                        "external_xrefs": ["REACTOME:R-HSA-70171"],
                        "checksum": "",
                        **common,
                    }
                ],
            )
            write_parquet(
                normalized / "reactions.parquet",
                [
                    {
                        "reaction_uid": "rxn_glucose",
                        "name": "Glucose phosphorylation",
                        "species": "Homo sapiens",
                        "primary_external_id": "R-HSA-1",
                        "pathway_uid": "path_glycolysis",
                        "pathway_external_id": "R-HSA-70171",
                        "source_name": "Reactome",
                        "external_xrefs": ["REACTOME:R-HSA-1"],
                        "checksum": "",
                        **common,
                    },
                    {
                        "reaction_uid": "rxn_pyruvate_reverse",
                        "name": "Pyruvate reverse fixture",
                        "species": "Homo sapiens",
                        "primary_external_id": "R-HSA-2",
                        "pathway_uid": "path_glycolysis",
                        "pathway_external_id": "R-HSA-70171",
                        "source_name": "Reactome",
                        "external_xrefs": ["REACTOME:R-HSA-2"],
                        "checksum": "",
                        **common,
                    },
                    {
                        "reaction_uid": "rxn_bidirectional",
                        "name": "Bidirectional fixture",
                        "species": "Homo sapiens",
                        "primary_external_id": "R-HSA-3",
                        "pathway_uid": "path_glycolysis",
                        "pathway_external_id": "R-HSA-70171",
                        "source_name": "Reactome",
                        "external_xrefs": ["REACTOME:R-HSA-3"],
                        "checksum": "",
                        **common,
                    }
                ],
            )
            write_parquet(
                normalized / "reaction_participants.parquet",
                [
                    {
                        "edge_uid": "edge_rp_1",
                        "reaction_uid": "rxn_glucose",
                        "participant_uid": "met_glucose",
                        "participant_external_id": "CHEBI:17234",
                        "participant_type": "metabolite",
                        "role": "input",
                        "source_name": "Reactome",
                        "source_record_id": "R-HSA-1|CHEBI:17234",
                        **common,
                    }
                ],
            )
            write_parquet(
                normalized / "sentences.parquet",
                [
                    {
                        "sentence_uid": "sent_1",
                        "article_uid": "art_1",
                        "pmid": "123",
                        "pmcid": "",
                        "section": "abstract",
                        "sentence_text": "Glucose metabolism is altered in tumors.",
                        "text_hash": "hash",
                        "start_offset": 0,
                        "end_offset": 40,
                        "language": "en",
                        "source_release": release,
                        "parser_hash": "parser",
                    }
                ],
            )
            european = root / "raw_lake" / "European"
            european.mkdir(parents=True)
            (european / "European_trait_annotations.csv").write_text(
                "\n".join(
                    [
                        "accession_id,reported_trait,summary_statistics_url,pubmed_id,paper_title,resolution_status,mapped_names,candidate_names",
                        "GCST90201000,Tryptophan to Pyruvate ratio,http://example.org/GCST90201000,123,Paper,local_name,Tryptophan|Pyruvate,Tryptophan|Pyruvate",
                    ]
                ),
                encoding="utf-8",
            )
            reactome = root / "raw_lake" / "reactome" / release
            reactome.mkdir(parents=True)
            (reactome / "ChEBI2Reactome_PE_Reactions.txt").write_text(
                "17234\tR-ALL-1\tGlucose [cytosol]\tR-HSA-1\thttps://reactome.org\tGlucose phosphorylation\tTAS\tHomo sapiens\n",
                encoding="utf-8",
            )
            (reactome / "reactome_reaction_exporter.txt").write_text(
                "pathway_id\treaction_id\treaction_name\tuniprot_acc\trole_in_reaction\n"
                "R-HSA-70171\tR-HSA-1\tGlucose phosphorylation\tP12345\tcatalyst\n",
                encoding="utf-8",
            )
            (reactome / "ReactionPMIDS.txt").write_text("R-HSA-1\t123456\n", encoding="utf-8")
            (reactome / "Homo_sapiens.sbml").write_text(
                """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core">
  <model id="reactome_fixture">
    <listOfSpecies>
      <species id="s_glucose" name="Glucose">
        <annotation>
          <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
            <rdf:Description rdf:about="#s_glucose">
              <rdf:li rdf:resource="http://identifiers.org/chebi/CHEBI:17234"/>
            </rdf:Description>
          </rdf:RDF>
        </annotation>
      </species>
      <species id="s_pyruvate" name="Pyruvate">
        <annotation>
          <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
            <rdf:Description rdf:about="#s_pyruvate">
              <rdf:li rdf:resource="http://identifiers.org/chebi/CHEBI:15361"/>
            </rdf:Description>
          </rdf:RDF>
        </annotation>
      </species>
      <species id="s_name_only" name="Name only fixture"/>
    </listOfSpecies>
    <listOfReactions>
      <reaction id="reaction_R-HSA-1" name="R-HSA-1 semantic fixture" reversible="false">
        <listOfReactants>
          <speciesReference species="s_glucose" stoichiometry="1"/>
        </listOfReactants>
        <listOfProducts>
          <speciesReference species="s_pyruvate" stoichiometry="1"/>
          <speciesReference species="s_name_only" stoichiometry="1"/>
        </listOfProducts>
      </reaction>
    </listOfReactions>
  </model>
</sbml>
""",
                encoding="utf-8",
            )
            rhea_dir = root / "raw_lake" / "rhea" / release
            rhea_dir.mkdir(parents=True, exist_ok=True)
            (rhea_dir / "rhea-biopax.owl").write_text(
                """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
         xmlns:bp="http://www.biopax.org/release/biopax-level3.owl#">
  <bp:UnificationXref rdf:ID="xref_rhea_10001">
    <bp:db>Rhea</bp:db>
    <bp:id>10001</bp:id>
  </bp:UnificationXref>
  <bp:UnificationXref rdf:ID="xref_chebi_17234">
    <bp:db>ChEBI</bp:db>
    <bp:id>CHEBI:17234</bp:id>
  </bp:UnificationXref>
  <bp:UnificationXref rdf:ID="xref_chebi_15361">
    <bp:db>ChEBI</bp:db>
    <bp:id>CHEBI:15361</bp:id>
  </bp:UnificationXref>
  <bp:SmallMolecule rdf:ID="sm_glucose">
    <bp:displayName>Glucose semantic</bp:displayName>
    <bp:xref rdf:resource="#xref_chebi_17234"/>
  </bp:SmallMolecule>
  <bp:SmallMolecule rdf:ID="sm_pyruvate">
    <bp:displayName>Pyruvate semantic</bp:displayName>
    <bp:xref rdf:resource="#xref_chebi_15361"/>
  </bp:SmallMolecule>
  <bp:BiochemicalReaction rdf:ID="rhea_reaction_10001">
    <bp:xref rdf:resource="#xref_rhea_10001"/>
    <bp:conversionDirection>LEFT-TO-RIGHT</bp:conversionDirection>
    <bp:left rdf:resource="#sm_glucose"/>
    <bp:right rdf:resource="#sm_pyruvate"/>
  </bp:BiochemicalReaction>
</rdf:RDF>
""",
                encoding="utf-8",
            )
            write_tsv_tar(
                root / "raw_lake" / "rhea" / release / "rhea-tsv.tar.gz",
                {
                    "tsv/rhea-directions.tsv": "RHEA_ID_MASTER\tRHEA_ID_LR\tRHEA_ID_RL\tRHEA_ID_BI\n10000\t10001\t10002\t10003\n",
                    "tsv/rhea2reactome.tsv": (
                        "RHEA_ID\tDIRECTION\tMASTER_ID\tID\n"
                        "10001\tLR\t10000\tR-HSA-1.1\n"
                        "10002\tRL\t10000\tR-HSA-2.1\n"
                        "10003\tBI\t10000\tR-HSA-3.1\n"
                    ),
                    "tsv/rhea-reaction-smiles.tsv": (
                        "10001\tC(C1C(C(C(C(O1)O)O)O)O)O>>CC(=O)[O-]\n"
                        "10002\tCC(=O)[O-]>>C(C1C(C(C(C(O1)O)O)O)O)O\n"
                        "10003\tC(C1C(C(C(C(O1)O)O)O)O)O>>CC(=O)[O-]\n"
                    ),
                    "tsv/rhea-chebi-smiles.tsv": (
                        "CHEBI:17234\tC(C1C(C(C(C(O1)O)O)O)O)O\n"
                        "CHEBI:15361\tCC(=O)[O-]\n"
                    ),
                    "tsv/chebiId_name.tsv": "CHEBI:17234\tglucose\nCHEBI:15361\tpyruvate\n",
                    "tsv/rhea2uniprot_sprot.tsv": "RHEA_ID\tDIRECTION\tMASTER_ID\tID\n10001\tLR\t10000\tP12345\n",
                },
            )

            manifest = database_accuracy_v2.build_database_accuracy_store(root, release)

            self.assertTrue(manifest.exists())
            out = root / "database_accuracy_store" / release
            traits = pq.read_table(out / "entity_store" / "trait_entities.parquet").to_pylist()
            self.assertEqual(traits[0]["trait_type"], "metabolite_ratio")
            self.assertEqual(traits[0]["identity_scope"], "ratio")
            components = pq.read_table(out / "entity_store" / "trait_components.parquet").to_pylist()
            self.assertEqual([row["component_role"] for row in components], ["numerator", "denominator"])
            self.assertEqual([row["direction_semantics"] for row in components], ["same_direction", "inverse_direction"])

            chemicals = {row["chemical_uid"]: row for row in pq.read_table(out / "entity_store" / "chemical_entities.parquet").to_pylist()}
            self.assertEqual(chemicals["met_glucose"]["entity_granularity"], "exact_compound")
            self.assertEqual(chemicals["met_sphingomyelin"]["identity_status"], "review")
            names = pq.read_table(out / "entity_store" / "chemical_names.parquet").to_pylist()
            sphingo_names = [row for row in names if row["chemical_uid"] == "met_sphingomyelin"]
            self.assertTrue(all(row["name_risk"] == "high" for row in sphingo_names))

            participants = pq.read_table(out / "relation_store" / "reaction_participants_v2.parquet").to_pylist()
            self.assertEqual(participants[0]["participant_role"], "substrate")
            self.assertTrue(any(row["relation_source"] == "reactome_pe_participation" for row in participants))
            self.assertTrue(any(row["relation_source"] == "rhea_biopax_explicit_side" for row in participants))
            self.assertTrue(any(row["relation_source"] == "reactome_sbml_species_reference" for row in participants))
            equations = pq.read_table(out / "relation_store" / "reaction_equations.parquet").to_pylist()
            self.assertTrue({"LR", "RL", "BI"}.issubset({row["directionality"] for row in equations}))
            side_rows = pq.read_table(out / "relation_store" / "reaction_side_participants.parquet").to_pylist()
            self.assertTrue(any(row["relation_source"] == "rhea_biopax_explicit_side" and row["semantic_source_uri"] for row in side_rows))
            self.assertTrue(any(row["relation_source"] == "reactome_sbml_species_reference" and row["semantic_source_uri"] for row in side_rows))
            lr_glucose = next(row for row in side_rows if row["reaction_uid"] == "rxn_glucose" and row["chemical_uid"] == "met_glucose")
            lr_pyruvate = next(row for row in side_rows if row["reaction_uid"] == "rxn_glucose" and row["chemical_uid"] == "met_pyruvate")
            self.assertEqual(lr_glucose["participant_role"], "substrate")
            self.assertEqual(lr_pyruvate["participant_role"], "product")
            rl_pyruvate = next(row for row in side_rows if row["reaction_uid"] == "rxn_pyruvate_reverse" and row["chemical_uid"] == "met_pyruvate")
            rl_glucose = next(row for row in side_rows if row["reaction_uid"] == "rxn_pyruvate_reverse" and row["chemical_uid"] == "met_glucose")
            self.assertEqual(rl_pyruvate["participant_role"], "substrate")
            self.assertEqual(rl_glucose["participant_role"], "product")
            xrefs = pq.read_table(out / "relation_store" / "reaction_xrefs.parquet").to_pylist()
            self.assertTrue(any(row["reaction_uid"] == "rxn_glucose" and row["xref_id"] == "10001" for row in xrefs))
            enzymes = pq.read_table(out / "relation_store" / "enzyme_reaction_links.parquet").to_pylist()
            self.assertTrue(any(row["reaction_uid"] == "rxn_glucose" and row["protein_source_id"] == "P12345" for row in enzymes))
            publications = pq.read_table(out / "relation_store" / "reaction_publication_links.parquet").to_pylist()
            self.assertTrue(any(row["reaction_uid"] == "rxn_glucose" and row["pmid"] == "123456" for row in publications))
            candidates = pq.read_table(out / "analysis_view" / "fact_candidates.parquet").to_pylist()
            self.assertTrue(any(row["candidate_fact_type"] == "reaction_participation" for row in candidates))
            self.assertTrue(any(row["candidate_fact_type"] == "reaction_side_participant" for row in candidates))
            self.assertTrue(any("unresolved_semantic_participant" in (row["blocking_reasons"] or []) for row in candidates))
            self.assertTrue(any(row["candidate_fact_type"] == "enzyme_reaction_link" for row in candidates))
            self.assertTrue(any(row["candidate_fact_type"] == "module_membership" for row in candidates))
            facts = pq.read_table(out / "analysis_view" / "mechanism_ready_facts.parquet").to_pylist()
            self.assertTrue(any(row["allowed_claim_scope"] == "directional_reaction_fact" for row in facts))
            self.assertTrue(
                any(
                    json.loads(row["metadata_json"]).get("relation_source") in {"rhea_biopax_explicit_side", "reactome_sbml_species_reference"}
                    and json.loads(row["metadata_json"]).get("semantic_source_uri")
                    for row in facts
                )
            )
            self.assertTrue(any(row["allowed_claim_scope"] == "bidirectional_reaction_fact" for row in facts))
            self.assertTrue(any(row["allowed_claim_scope"] == "role_unknown_reaction_fact" for row in facts))
            self.assertTrue(facts[0]["source_record_uids"])
            self.assertTrue(facts[0]["identity_decision_uids"])

            input_features = pq.read_table(out / "identity_resolution_store" / "input_features.parquet").to_pylist()
            identity_candidates = pq.read_table(out / "identity_resolution_store" / "identity_candidates.parquet").to_pylist()
            identity_decisions = pq.read_table(out / "identity_resolution_store" / "identity_decisions.parquet").to_pylist()
            self.assertTrue(input_features)
            self.assertTrue(identity_candidates)
            self.assertTrue(identity_decisions)

            by_feature = {row["input_feature_uid"]: row for row in input_features}
            decisions_by_feature = {row["input_feature_uid"]: row for row in identity_decisions}
            ratio_feature = next(row for row in input_features if row["feature_type"] == "ratio")
            ratio_decision = decisions_by_feature[ratio_feature["input_feature_uid"]]
            self.assertEqual(ratio_decision["decision_status"], "accepted_trait")
            self.assertEqual(ratio_decision["accepted_entity_type"], "trait")
            ratio_candidates = [row for row in identity_candidates if row["input_feature_uid"] == ratio_feature["input_feature_uid"]]
            component_candidates = [
                row
                for row in ratio_candidates
                if row["candidate_entity_type"] in {"chemical", "unknown", "class"}
            ]
            self.assertTrue(component_candidates)
            self.assertTrue(all(row["blocking_reason"] == "ratio_component_not_abundance_seed" for row in component_candidates))

            class_feature = next(row for row in input_features if row["feature_label"] == "Sphingomyelin")
            class_decision = decisions_by_feature[class_feature["input_feature_uid"]]
            self.assertEqual(class_decision["decision_status"], "accepted_class")
            self.assertEqual(class_decision["accepted_entity_type"], "class")

            sphingo_feature = next(
                row
                for row in input_features
                if row["input_row_id"] == "met_sphingomyelin" and row["feature_type"] == "direct_metabolite"
            )
            sphingo_decision = decisions_by_feature[sphingo_feature["input_feature_uid"]]
            self.assertEqual(sphingo_decision["decision_status"], "rejected")
            self.assertNotEqual(sphingo_decision["accepted_entity_type"], "chemical")

            manifest_data = json.loads((out / "database_accuracy_store_manifest.json").read_text(encoding="utf-8"))
            self.assertGreaterEqual(manifest_data["table_count"], 26)
            self.assertGreater(manifest_data["total_rows"], 0)
            self.assertTrue(all(row.get("schema_hash") for row in manifest_data["tables"]))


if __name__ == "__main__":
    unittest.main()
