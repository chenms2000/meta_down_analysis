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


build_compound_match_index = load_script(
    "build_compound_match_index",
    ROOT / "scripts" / "build_compound_match_index.py",
)
metabo_service = load_script("metabo_service_for_compound_tests", ROOT / "scripts" / "metabo_service.py")


def write_parquet(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


@unittest.skipIf(pa is None or pq is None, "pyarrow is required")
class CompoundMatchIndexTests(unittest.TestCase):
    def test_builds_compound_index_and_resolves_feature_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            release = "mvp_20260101T000000"
            normalized = root / "normalized_store" / release
            graph = root / "graph_projection" / release
            pubchem = root / "pubchem_cid_cache" / release
            compound = root / "compound_match_index"
            common = {"source_release": release, "license_id": "open_core:test"}

            write_parquet(
                normalized / "metabolites.parquet",
                [
                    {
                        "metabolite_uid": "met_glucose",
                        "canonical_name": "Glucose",
                        "synonyms": ["D-Glucose", "Shared ambiguous"],
                        "formula": "C6H12O6",
                        "exact_mass": 180.063388,
                        "inchikey": "WQZGKKKJIJFFOK-UHFFFAOYSA-N",
                        "external_xrefs": ["HMDB:HMDB0000122", "CHEBI:17234", "PUBCHEM.COMPOUND:5793"],
                        **common,
                    },
                    {
                        "metabolite_uid": "met_l_glucose",
                        "canonical_name": "L-Glucose",
                        "synonyms": ["Shared ambiguous"],
                        "formula": "H12C6O6",
                        "exact_mass": 180.063388,
                        "inchikey": "WQZGKKKJIJFFOK-DVKNGEFBSA-N",
                        "external_xrefs": ["CHEBI:37631", "PUBCHEM.COMPOUND:79025"],
                        **common,
                    },
                    {
                        "metabolite_uid": "met_anchor",
                        "canonical_name": "Pyruvate",
                        "synonyms": [],
                        "formula": "C3H4O3",
                        "exact_mass": 88.016044,
                        "inchikey": "LCTONWCANYUPML-UHFFFAOYSA-N",
                        "external_xrefs": ["HMDB:HMDB0000243"],
                        **common,
                    },
                ],
            )
            write_parquet(
                normalized / "metabolite_xrefs.parquet",
                [
                    {
                        "metabolite_uid": "met_glucose",
                        "xref_source": "HMDB",
                        "xref_id": "HMDB0000122",
                        **common,
                    },
                    {
                        "metabolite_uid": "met_anchor",
                        "xref_source": "HMDB",
                        "xref_id": "HMDB0000243",
                        **common,
                    }
                ],
            )
            write_parquet(
                pubchem / "cid_properties.parquet",
                [
                    {
                        "pubchem_cid": "5793",
                        "metabolite_uids": ["met_glucose"],
                        "title": "D-Glucose",
                        "formula": "C6H12O6",
                        "exact_mass": 180.063388,
                        "inchikey": "WQZGKKKJIJFFOK-UHFFFAOYSA-N",
                        "synonyms": ["glucose"],
                        **common,
                    },
                    {
                        "pubchem_cid": "79025",
                        "metabolite_uids": ["met_l_glucose"],
                        "title": "L-Glucose",
                        "formula": "C6H12O6",
                        "exact_mass": 180.063388,
                        "inchikey": "WQZGKKKJIJFFOK-DVKNGEFBSA-N",
                        "synonyms": [],
                        **common,
                    },
                ],
            )
            write_parquet(
                graph / "nodes.parquet",
                [
                    {
                        "node_uid": "met_glucose",
                        "node_idx": 0,
                        "node_type": "metabolite",
                        "canonical_name": "Glucose",
                        "display_name": "Glucose",
                        "primary_external_id": "CHEBI:17234",
                        "external_xrefs": ["HMDB:HMDB0000122", "PUBCHEM.COMPOUND:5793"],
                        "source_priority": "ChEBI",
                        **common,
                    },
                    {
                        "node_uid": "met_l_glucose",
                        "node_idx": 1,
                        "node_type": "metabolite",
                        "canonical_name": "L-Glucose",
                        "display_name": "L-Glucose",
                        "primary_external_id": "CHEBI:37631",
                        "external_xrefs": ["PUBCHEM.COMPOUND:79025"],
                        "source_priority": "ChEBI",
                        **common,
                    },
                    {
                        "node_uid": "met_anchor",
                        "node_idx": 2,
                        "node_type": "metabolite",
                        "canonical_name": "Pyruvate",
                        "display_name": "Pyruvate",
                        "primary_external_id": "HMDB:HMDB0000243",
                        "external_xrefs": ["HMDB:HMDB0000243"],
                        "source_priority": "HMDB",
                        **common,
                    },
                ],
            )
            write_parquet(
                graph / "resolver_index.parquet",
                [
                    {
                        "lookup_uid": "lk",
                        "entity_uid": "met_glucose",
                        "entity_type": "metabolite",
                        "namespace": "TEXT",
                        "lookup_key": "glucose",
                        "raw_value": "Glucose",
                        "match_field": "canonical_name",
                        "rank": 100.0,
                        "source_table": "metabolites",
                        **common,
                    }
                ],
            )

            manifest_path = build_compound_match_index.build_compound_match_index(
                root / "normalized_store",
                root / "pubchem_cid_cache",
                compound,
                release,
            )
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["release_id"], release)
            self.assertTrue((compound / release / "compound_formula_mass_index.parquet").exists())
            self.assertTrue((compound / release / "compound_rt_index.parquet").exists())

            write_parquet(
                compound / release / "compound_rt_index.parquet",
                [
                    {
                        "lookup_uid": "rt_1",
                        "metabolite_uid": "met_glucose",
                        "rt": 5.1,
                        "rt_unit": "min",
                        "method_key": "hilic_pos",
                        "source_name": "fixture",
                        "source_record_id": "rt_glucose",
                        **common,
                        "parser_hash": "parser",
                    },
                    {
                        "lookup_uid": "rt_2",
                        "metabolite_uid": "met_l_glucose",
                        "rt": 7.0,
                        "rt_unit": "min",
                        "method_key": "hilic_pos",
                        "source_name": "fixture",
                        "source_record_id": "rt_l_glucose",
                        **common,
                        "parser_hash": "parser",
                    },
                ],
            )
            write_parquet(
                compound / release / "compound_ms2_index.parquet",
                [
                    {
                        "lookup_uid": "ms2_1",
                        "metabolite_uid": "met_glucose",
                        "precursor_mz": 181.070664,
                        "adduct": "[M+H]+",
                        "ion_mode": "positive",
                        "fragment_mz": 85.029,
                        "fragment_intensity": 100.0,
                        "source_name": "fixture",
                        "source_record_id": "spec_glucose",
                        **common,
                        "parser_hash": "parser",
                    },
                    {
                        "lookup_uid": "ms2_2",
                        "metabolite_uid": "met_glucose",
                        "precursor_mz": 181.070664,
                        "adduct": "[M+H]+",
                        "ion_mode": "positive",
                        "fragment_mz": 97.029,
                        "fragment_intensity": 50.0,
                        "source_name": "fixture",
                        "source_record_id": "spec_glucose",
                        **common,
                        "parser_hash": "parser",
                    },
                    {
                        "lookup_uid": "ms2_3",
                        "metabolite_uid": "met_l_glucose",
                        "precursor_mz": 181.070664,
                        "adduct": "[M+H]+",
                        "ion_mode": "positive",
                        "fragment_mz": 111.044,
                        "fragment_intensity": 100.0,
                        "source_name": "fixture",
                        "source_record_id": "spec_l_glucose",
                        **common,
                        "parser_hash": "parser",
                    },
                ],
            )
            write_parquet(
                normalized / "metabolite_pathway_edges.parquet",
                [
                    {
                        "edge_uid": "mp_glucose",
                        "metabolite_uid": "met_glucose",
                        "pathway_uid": "pathway_glycolysis",
                        "source_name": "Reactome",
                        "source_record_id": "glucose|glycolysis",
                        "evidence_level": "curated",
                        **common,
                    },
                    {
                        "edge_uid": "mp_anchor",
                        "metabolite_uid": "met_anchor",
                        "pathway_uid": "pathway_glycolysis",
                        "source_name": "Reactome",
                        "source_record_id": "pyruvate|glycolysis",
                        "evidence_level": "curated",
                        **common,
                    },
                    {
                        "edge_uid": "mp_l_glucose",
                        "metabolite_uid": "met_l_glucose",
                        "pathway_uid": "pathway_other",
                        "source_name": "Reactome",
                        "source_record_id": "l_glucose|other",
                        "evidence_level": "curated",
                        **common,
                    },
                ],
            )
            write_parquet(
                normalized / "pathways.parquet",
                [
                    {
                        "pathway_uid": "pathway_glycolysis",
                        "name": "Glycolysis",
                        "primary_external_id": "R-HSA-70171",
                        "source_name": "Reactome",
                        "species": "Homo sapiens",
                    },
                    {
                        "pathway_uid": "pathway_other",
                        "name": "Other",
                        "primary_external_id": "R-HSA-1",
                        "source_name": "Reactome",
                        "species": "Homo sapiens",
                    },
                ],
            )

            service = metabo_service.MetaboService(root, release_id=release)
            exact_id = service.precheck_metabolites([{"HMDB": "HMDB0000122"}])
            self.assertEqual(exact_id["summary"]["matched"], 1)
            self.assertEqual(exact_id["matched"][0]["metabolite_uid"], "met_glucose")

            ambiguous_name_only = service.precheck_metabolites([{"name": "Shared ambiguous"}])
            self.assertEqual(ambiguous_name_only["summary"]["ambiguous"], 1)
            name_note = ambiguous_name_only["ambiguous"][0]["resolution"]["query_notes"][0]
            self.assertEqual(name_note["triage_policy"], "auto_abstain_name_only")
            self.assertIn("name_maps_to_multiple_metabolites", name_note["triage_reasons"])
            candidates = ambiguous_name_only["ambiguous"][0]["resolution"]["candidates"]
            self.assertGreaterEqual(len(candidates), 2)
            self.assertTrue(all(candidate["score_components"]["name"] == 0.0 for candidate in candidates))

            feature = service.precheck_metabolites(
                [
                    {
                        "name": "glucose",
                        "formula": "C6H12O6",
                        "mz": 181.070664466812,
                        "adduct": "[M+H]+",
                        "ppm_tolerance": 5,
                    }
                ]
            )
            self.assertEqual(feature["summary"]["matched"], 1)
            components = feature["matched"][0]["resolution"]["candidates"][0]["score_components"]
            self.assertGreaterEqual(components["mass"], 49.0)
            self.assertEqual(components["formula"], 20.0)

            mass_only = service.precheck_metabolites([{"mz": 181.070664466812, "adduct": "[M+H]+", "ppm_tolerance": 5}])
            self.assertEqual(mass_only["summary"]["ambiguous"], 1)

            rt_ms2 = service.precheck_metabolites(
                [
                    {
                        "formula": "C6H12O6",
                        "mz": 181.070664466812,
                        "adduct": "[M+H]+",
                        "ppm_tolerance": 5,
                        "rt": 5.11,
                        "rt_method": "hilic_pos",
                        "rt_tolerance": 0.2,
                        "ms2_peaks": [[85.03, 100], [97.028, 50]],
                        "ms2_mz_tolerance": 0.02,
                        "ion_mode": "positive",
                    }
                ]
            )
            self.assertEqual(rt_ms2["summary"]["matched"], 1)
            rt_ms2_components = rt_ms2["matched"][0]["resolution"]["candidates"][0]["score_components"]
            self.assertGreater(rt_ms2_components["rt"], 10.0)
            self.assertGreater(rt_ms2_components["ms2"], 15.0)

            context_reranked = service.precheck_metabolites(
                [
                    {"HMDB": "HMDB0000243"},
                    {"formula": "C6H12O6", "mz": 181.070664466812, "adduct": "[M+H]+", "ppm_tolerance": 5},
                ]
            )
            self.assertEqual(context_reranked["summary"]["matched"], 1)
            self.assertEqual(context_reranked["summary"]["ambiguous"], 1)
            top_context_candidate = context_reranked["ambiguous"][0]["resolution"]["candidates"][0]
            self.assertEqual(top_context_candidate["entity_uid"], "met_glucose")
            self.assertGreater(top_context_candidate["score_components"]["context"], 0)

            connectivity_only = service.precheck_metabolites(
                [{"inchikey": "WQZGKKKJIJFFOK-AAAAAAAAAA-N"}]
            )
            self.assertEqual(connectivity_only["summary"]["ambiguous"], 1)

    def test_imports_reviewed_compound_identity_overlay(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            release = "mvp_20260101T000000"
            normalized = root / "normalized_store" / release
            compound = root / "compound_match_index"
            manual = root / "manual_sources" / "compound_identity_overlays" / release
            common = {"source_release": release, "license_id": "open_core:test"}

            write_parquet(
                normalized / "metabolites.parquet",
                [
                    {
                        "metabolite_uid": "met_glucose",
                        "canonical_name": "Glucose",
                        "synonyms": [],
                        "formula": "C6H12O6",
                        "exact_mass": 180.063388,
                        "inchikey": "WQZGKKKJIJFFOK-UHFFFAOYSA-N",
                        "external_xrefs": ["PUBCHEM.COMPOUND:5793"],
                        **common,
                    },
                    {
                        "metabolite_uid": "met_overlay_structure_duplicate",
                        "canonical_name": "Overlay structure duplicate",
                        "synonyms": [],
                        "formula": "C2H6O",
                        "exact_mass": 46.041865,
                        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                        "external_xrefs": ["CHEBI:999999"],
                        **common,
                    }
                ],
            )
            manual.mkdir(parents=True, exist_ok=True)
            (manual / "compound_identity_overlay.csv").write_text(
                "\n".join(
                    [
                        "metabolite_uid,canonical_name,synonyms,pubchem_cid,inchikey,formula,exact_mass,iupac_name,source_release,license_id",
                        f"manual_pubchem_cid_999999,Overlayine,Overlay acid|Overlay metabolite,999999,AAAAAAAAAAAAAA-BBBBBBBBBB-C,C2H6O,,overlay iupac,{release},local:test",
                        f"manual_pubchem_cid_5793,Duplicate glucose,duplicate,5793,WQZGKKKJIJFFOK-UHFFFAOYSA-N,C6H12O6,,,{release},local:test",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            manifest_path = build_compound_match_index.build_compound_match_index(
                root / "normalized_store",
                root / "pubchem_cid_cache",
                compound,
                release,
                manual_compound_root=root / "manual_sources" / "compound_identity_overlays",
            )
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["identity_overlay_import"]["imported_rows"], 1)
            self.assertEqual(manifest["identity_overlay_import"]["skipped_existing_cid"], 1)

            id_rows = pq.read_table(compound / release / "compound_identifier_index.parquet").to_pylist()
            name_rows = pq.read_table(compound / release / "compound_name_index.parquet").to_pylist()
            inchikey_rows = pq.read_table(compound / release / "compound_inchikey_index.parquet").to_pylist()
            overlay_ids = [
                row
                for row in id_rows
                if row["metabolite_uid"] == "manual_pubchem_cid_999999" and row["namespace"] == "CID"
            ]
            self.assertEqual({row["lookup_key"] for row in overlay_ids}, {"999999"})
            self.assertIn("Overlay acid", {row["raw_value"] for row in name_rows})
            self.assertIn("overlay iupac", {row["raw_value"] for row in name_rows})
            self.assertEqual(
                [row["inchikey"] for row in inchikey_rows if row["metabolite_uid"] == "manual_pubchem_cid_999999"],
                ["AAAAAAAAAAAAAA-BBBBBBBBBB-C"],
            )
            service = metabo_service.MetaboService(
                workspace=root,
                normalized_root=root / "normalized_store",
                graph_root=root / "graph_projection",
                pubchem_root=root / "pubchem_cid_cache",
                compound_root=compound,
                release_id=release,
            )
            exact_overlay = service.precheck_metabolites(
                [
                    {
                        "name": "Overlayine",
                        "pubchem_cid": "999999",
                        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                    }
                ]
            )
            self.assertEqual(exact_overlay["summary"]["matched"], 1)
            self.assertEqual(
                exact_overlay["matched"][0]["resolution"]["candidates"][0]["entity_uid"],
                "manual_pubchem_cid_999999",
            )

    def test_imports_rt_ms2_jsonl_and_msp_with_external_id_resolution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            release = "mvp_20260101T000001"
            normalized = root / "normalized_store" / release
            pubchem = root / "pubchem_cid_cache" / release
            compound = root / "compound_match_index"
            manual = root / "manual_sources" / "compound_features" / release
            common = {"source_release": release, "license_id": "open_core:test"}

            write_parquet(
                normalized / "metabolites.parquet",
                [
                    {
                        "metabolite_uid": "met_glucose",
                        "canonical_name": "Glucose",
                        "synonyms": ["D-Glucose"],
                        "formula": "C6H12O6",
                        "exact_mass": 180.063388,
                        "inchikey": "WQZGKKKJIJFFOK-UHFFFAOYSA-N",
                        "external_xrefs": ["HMDB:HMDB0000122", "CHEBI:17234"],
                        **common,
                    }
                ],
            )
            write_parquet(
                normalized / "metabolite_xrefs.parquet",
                [{"metabolite_uid": "met_glucose", "xref_source": "HMDB", "xref_id": "HMDB0000122", **common}],
            )
            manual.mkdir(parents=True)
            (manual / "compound_rt.jsonl").write_text(
                json.dumps(
                    {
                        "hmdb_id": "HMDB:HMDB0000122",
                        "rt": 5.12,
                        "method": "hilic_pos",
                        "source_name": "fixture_rt",
                        "source_record_id": "rt_hmdb_glucose",
                        "license_id": "open_core:fixture",
                    }
                )
                + "\n",
                encoding="utf-8",
            )
            (manual / "compound_ms2.msp").write_text(
                "\n".join(
                    [
                        "Name: Glucose",
                        "HMDB: HMDB0000122",
                        "PrecursorMZ: 181.070664",
                        "Precursor_type: [M+H]+",
                        "Ion_mode: positive",
                        "Num Peaks: 2",
                        "85.030 100",
                        "97.028 50",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            manifest_path = build_compound_match_index.build_compound_match_index(
                root / "normalized_store",
                root / "pubchem_cid_cache",
                compound,
                release,
                manual_feature_root=root / "manual_sources" / "compound_features",
            )
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["feature_import"]["rt_imported_rows"], 1)
            self.assertEqual(manifest["feature_import"]["ms2_imported_peaks"], 2)
            self.assertIn("compound_rt.jsonl", manifest["feature_import"]["formats"])
            self.assertIn("compound_ms2.msp", manifest["feature_import"]["formats"])

            rt_rows = pq.read_table(compound / release / "compound_rt_index.parquet").to_pylist()
            ms2_rows = pq.read_table(compound / release / "compound_ms2_index.parquet").to_pylist()
            self.assertEqual(rt_rows[0]["metabolite_uid"], "met_glucose")
            self.assertEqual({row["metabolite_uid"] for row in ms2_rows}, {"met_glucose"})
            self.assertEqual(sorted(row["fragment_mz"] for row in ms2_rows), [85.03, 97.028])


if __name__ == "__main__":
    unittest.main()
