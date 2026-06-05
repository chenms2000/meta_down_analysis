import importlib.util
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_normalized_store.py"
SPEC = importlib.util.spec_from_file_location("build_normalized_store", SCRIPT)
build_normalized_store = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = build_normalized_store
SPEC.loader.exec_module(build_normalized_store)


class NormalizedBuildTests(unittest.TestCase):
    def test_stable_uid_does_not_embed_external_identifier(self):
        uid = build_normalized_store.stable_uid("metabolite", "CHEBI:15377")
        self.assertTrue(uid.startswith("metabolite_"))
        self.assertNotIn("CHEBI", uid)
        self.assertEqual(uid, build_normalized_store.stable_uid("metabolite", "CHEBI:15377"))

    def test_normalize_xref_handles_obo_suffixes(self):
        parsed = build_normalized_store.normalize_xref('pubchem.compound:5288826 {source="pubchem.compound"}')
        self.assertEqual(parsed, ("PUBCHEM.COMPOUND", "5288826", "PUBCHEM.COMPOUND:5288826"))
        parsed = build_normalized_store.normalize_xref("MESH:D006394 ! hemangiosarcoma")
        self.assertEqual(parsed, ("MESH", "D006394", "MESH:D006394"))

    def test_bridgedb_xref_normalizes_system_codes(self):
        self.assertEqual(build_normalized_store.bridgedb_xref("Ce", "CHEBI:15377"), "CHEBI:15377")
        self.assertEqual(build_normalized_store.bridgedb_xref("Ce", "15377"), "CHEBI:15377")
        self.assertEqual(build_normalized_store.bridgedb_xref("Ch", "HMDB0000001"), "HMDB:HMDB0000001")
        self.assertEqual(build_normalized_store.bridgedb_xref("Cpc", "5793"), "PUBCHEM.COMPOUND:5793")
        self.assertIsNone(build_normalized_store.bridgedb_xref("Cpc", "100.044.355"))
        self.assertEqual(build_normalized_store.bridgedb_xref("Ik", "WQZGKKKJIJFFOK-GASJEMHNSA-N"), "INCHIKEY:WQZGKKKJIJFFOK-GASJEMHNSA-N")

    def test_iter_obo_terms_extracts_term_blocks(self):
        content = """format-version: 1.2

[Term]
id: MONDO:0000001
name: disease A
synonym: "A disease" EXACT []
xref: DOID:123
is_a: MONDO:0000000 ! root

[Typedef]
id: part_of
"""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "mini.obo"
            path.write_text(content, encoding="utf-8")
            terms = list(build_normalized_store.iter_obo_terms(path))
        self.assertEqual(len(terms), 1)
        self.assertEqual(terms[0]["id"], ["MONDO:0000001"])
        self.assertEqual(terms[0]["xref"], ["DOID:123"])

    def test_latest_release_id_uses_mvp_manifest_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            older = root / "mvp_20260101T000000.manifest.json"
            newer = root / "mvp_20260102T000000.manifest.json"
            older.write_text("{}", encoding="utf-8")
            newer.write_text("{}", encoding="utf-8")
            self.assertEqual(build_normalized_store.latest_release_id(root), "mvp_20260102T000000")

    def test_load_catalog_accepts_windows_style_manifest_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = root / "config"
            config.mkdir()
            (config / "source_catalog.toml").write_text('[[sources]]\nid = "test"\n', encoding="utf-8")
            catalog = build_normalized_store.load_catalog(root, {"catalog": "config\\source_catalog.toml"})
        self.assertEqual(catalog["sources"][0]["id"], "test")

    def test_parse_article_metadata_extracts_core_fields(self):
        text = """
JOURNAL INFORMATION
==============================
NLM Title Abbreviation: Cells

ARTICLE INFORMATION
==============================
PMCID: PMC123
PMID: 456
DOI: 10.1/example
Subjects: Article

Metabolism in Cancer Cells

License: CC BY

==============================
This is the abstract sentence. It mentions cancer metabolism.

1. Introduction

Body text.
"""
        meta = build_normalized_store.parse_article_metadata(text, Path("PMC123.1.txt"))
        self.assertEqual(meta["pmcid"], "PMC123")
        self.assertEqual(meta["pmid"], "456")
        self.assertEqual(meta["title"], "Metabolism in Cancer Cells")
        self.assertIn("abstract sentence", meta["abstract"])

    def test_hybrid_sentence_parser_falls_back_when_scispacy_unavailable(self):
        if build_normalized_store.pa is None:
            self.skipTest("pyarrow is required")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            articles = root / "articles_collect"
            articles.mkdir()
            (articles / "PMC1.txt").write_text(
                """
PMCID: PMC1
PMID: 1
Subjects: Article

Cancer metabolism

==============================
Lactate is elevated in cancer. Glucose was measured by LC-MS.

1. Introduction

Body text.
""",
                encoding="utf-8",
            )
            ctx = build_normalized_store.BuildContext(
                workspace=root,
                raw_root=root / "raw_lake",
                manifest_dir=root / "manifests",
                output_root=root / "normalized_store",
                release_id="mvp_20260101T000000",
                manifest={},
                catalog={},
                parser_hash="parser",
            )
            with mock.patch.object(build_normalized_store, "load_scispacy_pipeline", side_effect=RuntimeError("missing scispacy")):
                build_normalized_store.build_articles(ctx, articles, 500, "abstract", "hybrid", "en_core_sci_sm")
            notes = " ".join(row["message"] for row in ctx.notes)
            self.assertIn("falling back to rules", notes)
            table_names = {row["table"] for row in ctx.table_audit}
            self.assertIn("sentence_entity_candidates", table_names)
            sentence_path = ctx.output_dir / "sentences.parquet"
            rows = build_normalized_store.pq.read_table(sentence_path).to_pylist()
            self.assertTrue(rows)
            self.assertTrue(all(row["segmenter"] == "rules" for row in rows))

    def test_hmdb_metabolite_payload_extracts_xrefs_and_annotations(self):
        xml = """<hmdb xmlns="http://www.hmdb.ca">
  <metabolite>
    <accession>HMDB0000001</accession>
    <name>1-Methylhistidine</name>
    <status>Detected and quantified</status>
    <chemical_formula>C7H11N3O2</chemical_formula>
    <monisotopic_molecular_weight>169.085126611</monisotopic_molecular_weight>
    <inchikey>BRMWTNUJHUMWMS-LURJTMIESA-N</inchikey>
    <chebi_id>50599</chebi_id>
    <pubchem_compound_id>92105</pubchem_compound_id>
    <kegg_id>C01152</kegg_id>
    <synonyms><synonym>Pi-methylhistidine</synonym></synonyms>
    <biospecimen_locations><biospecimen>Blood</biospecimen></biospecimen_locations>
    <pathways><pathway><name>Histidine metabolism</name><smpdb_id>SMP00001</smpdb_id></pathway></pathways>
    <protein_associations><protein><name>Example enzyme</name><uniprot_id>P12345</uniprot_id></protein></protein_associations>
    <diseases><disease><name>Example disease</name><omim_id>123456</omim_id></disease></diseases>
  </metabolite>
</hmdb>"""
        root = build_normalized_store.ET.fromstring(xml)
        elem = build_normalized_store.child_by_name(root, "metabolite")
        self.assertIsNotNone(elem)
        payload = build_normalized_store.hmdb_metabolite_payload(elem)
        self.assertEqual(payload["accession"], "HMDB0000001")
        self.assertIn("HMDB:HMDB0000001", payload["external_xrefs"])
        self.assertIn("CHEBI:50599", payload["external_xrefs"])
        self.assertIn("PUBCHEM.COMPOUND:92105", payload["external_xrefs"])
        self.assertEqual(payload["biospecimen_locations"], ["Blood"])
        self.assertIn("Histidine metabolism|SMPDB:SMP00001", payload["pathway_refs"])

    def test_iter_hmdb_metabolite_elements_reads_zip_xml(self):
        xml = """<hmdb><metabolite><accession>HMDB0000001</accession><name>A</name></metabolite></hmdb>"""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hmdb_metabolites.zip"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("hmdb_metabolites.xml", xml)
            values = [
                build_normalized_store.child_text_local(elem, "accession")
                for elem in build_normalized_store.iter_hmdb_metabolite_elements(path)
            ]
        self.assertEqual(values, ["HMDB0000001"])

    def test_parse_ncit_cancer_terms_filters_non_oncology_terms(self):
        content = """format-version: 1.2

[Term]
id: NCIT:C2926
name: Lung Carcinoma
synonym: "Pulmonary carcinoma" EXACT []
is_a: NCIT:C3262 ! Neoplasm

[Term]
id: NCIT:C12468
name: Lung
is_a: NCIT:C12219 ! Organ
"""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ncit.obo"
            path.write_text(content, encoding="utf-8")
            rows = build_normalized_store.parse_ncit_cancer_terms(path, None)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["primary_external_id"], "NCIT:C2926")
        self.assertIn("NCIT:C3262", rows[0]["parents"])

    def test_parse_ncit_cancer_terms_filters_non_disease_concepts(self):
        content = """format-version: 1.2

[Term]
id: NCIT:C2926
name: Lung Carcinoma
synonym: "Pulmonary carcinoma" EXACT []
is_a: NCIT:C3262 ! Neoplasm

[Term]
id: NCIT:C20535
name: Tumor Necrosis Factor
def: "Tumor necrosis factor is encoded by the human TNF gene." []
synonym: "TNF" EXACT []

[Term]
id: NCIT:C163758
name: Targeted Therapy Agent
def: "An agent that targets cells of the tumor microenvironment." []
synonym: "Targeting" EXACT []

[Term]
id: NCIT:C202420
name: c-Myc Measurement
def: "The determination of the amount of c-Myc in a sample." []
synonym: "Tumor Protein c-myc" EXACT []

[Term]
id: NCIT:C17965
name: BRCA1 Gene
synonym: "BRCA1" EXACT []

[Term]
id: NCIT:C99999
name: KIT Positive
synonym: "Positive" EXACT []
"""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ncit.obo"
            path.write_text(content, encoding="utf-8")
            rows = build_normalized_store.parse_ncit_cancer_terms(path, None)
        self.assertEqual([row["primary_external_id"] for row in rows], ["NCIT:C2926"])

    def test_parse_oncotree_tumor_types_maps_nci_xrefs(self):
        payload = [
            {
                "code": "LUAD",
                "name": "Lung Adenocarcinoma",
                "mainType": "Non-Small Cell Lung Cancer",
                "tissue": "Lung",
                "parent": "NSCLC",
                "externalReferences": {"NCI": ["C3512"], "UMLS": ["C0152013"]},
                "children": {},
            }
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "oncotree.json"
            path.write_text(build_normalized_store.json.dumps(payload), encoding="utf-8")
            rows = build_normalized_store.parse_oncotree_tumor_types(path)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["primary_external_id"], "ONCOTREE:LUAD")
        self.assertIn("NCIT:C3512", rows[0]["external_xrefs"])
        self.assertIn("ONCOTREE:NSCLC", rows[0]["parents"])

    def test_parse_context_axis_obo_terms(self):
        content = """format-version: 1.2

[Term]
id: CL:0000066
name: epithelial cell
synonym: "epitheliocyte" EXACT []
is_a: CL:0000003 ! native cell
xref: UMLS:C0014597

[Term]
id: UBERON:0002048
name: lung
synonym: "pulmo" EXACT []
is_a: UBERON:0000062 ! organ
"""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "axis.obo"
            path.write_text(content, encoding="utf-8")
            cell_rows = build_normalized_store.parse_obo_axis_terms(path, {"CL"}, "Cell Ontology")
            tissue_rows = build_normalized_store.parse_obo_axis_terms(path, {"UBERON"}, "UBERON")
        self.assertEqual(cell_rows[0]["primary_external_id"], "CL:0000066")
        self.assertIn("epitheliocyte", cell_rows[0]["aliases"])
        self.assertIn("CL:0000003", cell_rows[0]["parents"])
        self.assertEqual(tissue_rows[0]["primary_external_id"], "UBERON:0002048")

    def test_cell_state_seed_axis_contains_metabolic_state(self):
        rows = [row for row, source_id in build_normalized_store.cell_state_source_rows()]
        ids = {row["primary_external_id"] for row in rows}
        self.assertIn("METABO_STATE:glycolytic", ids)
        self.assertTrue(all(row["external_xrefs"] for row in rows))


if __name__ == "__main__":
    unittest.main()
