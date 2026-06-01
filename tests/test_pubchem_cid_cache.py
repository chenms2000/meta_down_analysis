import gzip
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
SCRIPT = ROOT / "scripts" / "build_pubchem_cid_cache.py"
SPEC = importlib.util.spec_from_file_location("build_pubchem_cid_cache", SCRIPT)
build_pubchem_cid_cache = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = build_pubchem_cid_cache
SPEC.loader.exec_module(build_pubchem_cid_cache)


def write_parquet(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


def write_gzip(path: Path, lines: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("".join(lines))


@unittest.skipIf(pa is None or pq is None, "pyarrow is required")
class PubChemCidCacheTests(unittest.TestCase):
    def test_builds_cache_for_referenced_cids_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            release = "mvp_20260101T000000"
            normalized = root / "normalized_store" / release
            raw = root / "raw_lake" / "pubchem_compound_extras" / release
            output = root / "pubchem_cid_cache"
            write_parquet(
                normalized / "metabolite_xrefs.parquet",
                [
                    {
                        "xref_uid": "xref_1",
                        "metabolite_uid": "met_1",
                        "xref_source": "PUBCHEM.COMPOUND",
                        "xref_id": "123",
                        "xref_key": "PUBCHEM.COMPOUND:123",
                        "source_name": "ChEBI",
                        "source_release": release,
                        "license_id": "open_core:chebi",
                        "parser_hash": "parser",
                    },
                    {
                        "xref_uid": "xref_2",
                        "metabolite_uid": "met_2",
                        "xref_source": "PUBCHEM.COMPOUND",
                        "xref_id": "0",
                        "xref_key": "PUBCHEM.COMPOUND:0",
                        "source_name": "HMDB",
                        "source_release": "HMDB 5.0",
                        "license_id": "licensed_enhancement:hmdb",
                        "parser_hash": "parser",
                    },
                    {
                        "xref_uid": "xref_3",
                        "metabolite_uid": "gene_like",
                        "xref_source": "NCBI.GENE",
                        "xref_id": "123",
                        "xref_key": "NCBI.GENE:123",
                        "source_name": "NCBI",
                        "source_release": release,
                        "license_id": "open_core:ncbi_gene_human",
                        "parser_hash": "parser",
                    },
                ],
            )
            write_parquet(
                normalized / "raw_file_manifest.parquet",
                [
                    {
                        "source_id": "pubchem_compound_extras",
                        "path": f"raw_lake/pubchem_compound_extras/{release}/CID-Title.gz",
                        "bytes": 10,
                        "sha256": "abc",
                    }
                ],
            )
            write_gzip(raw / "CID-Title.gz", ["123\tGlucose\n", "999\tNot used\n"])
            write_gzip(raw / "CID-InChI-Key.gz", ["123\tInChI=1S/example\tABCDEFGHIJKLMN-ABCDEFGHIJ-A\n"])
            write_gzip(raw / "CID-SMILES.gz", ["123\tC(C1C(C(C(C(O1)O)O)O)O)O\n"])
            write_gzip(raw / "CID-Mass.gz", ["123\tC6H12O6\t180.063388\t180.063388\n"])
            write_gzip(raw / "CID-Synonym-filtered.gz", ["123\tD-Glucose\n", "123\tDextrose\n", "999\tNope\n"])

            manifest_path = build_pubchem_cid_cache.build_pubchem_cid_cache(
                root / "normalized_store",
                root / "raw_lake",
                output,
                release,
                max_synonyms_per_cid=1,
                include_synonyms=True,
            )

            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["target_cid_count"], 1)
            self.assertEqual(manifest["matched_cid_count"], 1)
            props = pq.read_table(output / release / "cid_properties.parquet").to_pylist()
            self.assertEqual(len(props), 1)
            self.assertEqual(props[0]["title"], "Glucose")
            self.assertEqual(props[0]["synonyms"], ["D-Glucose"])
            self.assertEqual(props[0]["synonym_count"], 2)
            links = pq.read_table(output / release / "cid_metabolite_links.parquet").to_pylist()
            self.assertEqual(len(links), 2)
            lookup = pq.read_table(output / release / "cid_lookup_index.parquet").to_pylist()
            self.assertTrue(any(row["namespace"] == "INCHIKEY" for row in lookup))
            self.assertTrue(any(row["raw_value"] == "Glucose" for row in lookup))


if __name__ == "__main__":
    unittest.main()
