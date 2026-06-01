import importlib.util
import csv
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_european_trait_annotations.py"
SPEC = importlib.util.spec_from_file_location("build_european_trait_annotations", SCRIPT)
build_european_trait_annotations = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = build_european_trait_annotations
SPEC.loader.exec_module(build_european_trait_annotations)


class EuropeanTraitAnnotationTests(unittest.TestCase):
    def test_expands_composite_hode_trait(self):
        candidates = build_european_trait_annotations.expanded_trait_candidates("13-HODE + 9-HODE levels")

        self.assertIn("13-HODE", candidates)
        self.assertIn("9-HODE", candidates)
        self.assertIn("13-hydroxyoctadecadienoic acid", candidates)
        self.assertIn("9-hydroxyoctadecadienoic acid", candidates)

    def test_expands_gpg_lipid_abbreviation(self):
        candidates = build_european_trait_annotations.expanded_trait_candidates("1-stearoyl-GPG (18:0) levels")

        self.assertIn("1-stearoyl-glycerophosphoglycerol (18:0)", candidates)

    def test_stereochemistry_relaxation_does_not_strip_plain_names(self):
        candidates = build_european_trait_annotations.expanded_trait_candidates("Serine levels")

        self.assertIn("Serine", candidates)
        self.assertNotIn("erine", candidates)

    def test_buckets_platform_codes_for_manual_mapping(self):
        bucket, action, note = build_european_trait_annotations.unresolved_review_bucket("X-21319 levels")

        self.assertEqual(bucket, "platform_x_code")
        self.assertEqual(action, "manual_source_mapping_required")
        self.assertIn("Platform-internal", note)

    def test_unresolved_review_preserves_manual_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "review.csv"
            path.write_text(
                "accession_id,review_status,reviewed_name,reviewed_pubchem_cid,reviewed_hmdb_id,reviewed_inchikey,evidence_source\n"
                "GCST000001,reviewed,Manual Name,12345,,,supplement\n",
                encoding="utf-8",
            )

            existing = build_european_trait_annotations.load_existing_review(path)
            build_european_trait_annotations.write_unresolved_review(
                path,
                [
                    {
                        "accession_id": "GCST000001",
                        "reported_trait": "X-1 levels",
                        "resolution_status": "unresolved",
                        "candidate_names": "X-1",
                        "mapped_names": "X-1",
                        "pubchem_query": "",
                    }
                ],
                existing,
            )

            row = next(csv.DictReader(path.open(encoding="utf-8-sig")))
            self.assertEqual(row["review_status"], "reviewed")
            self.assertEqual(row["reviewed_name"], "Manual Name")
            self.assertEqual(row["reviewed_pubchem_cid"], "12345")

    def test_build_annotations_carries_european_source_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            european_dir = root / "raw_lake" / "European"
            european_dir.mkdir(parents=True, exist_ok=True)
            european_csv = european_dir / "European.csv"
            european_csv.write_text(
                "firstAuthor,publicationDate,journal,title,efoTraits,bgTraits,summaryStatistics,pubmedId,initialSampleDescription,discoverySampleAncestry,accessionId,reportedTrait\n"
                "Smith,2024,Nature,Metabolite GWAS,metabolite measurement,biological process,http://example.org/GCST000001,123456,100 cases,European,GCST000001,Serine levels\n",
                encoding="utf-8",
            )

            args = type(
                "Args",
                (),
                {
                    "workspace": str(root),
                    "release_id": "mvp_test",
                    "compound_root": Path("compound_match_index"),
                    "european_csv": "raw_lake/European/European.csv",
                    "output": "raw_lake/European/European_trait_annotations.csv",
                    "review_output": "raw_lake/European/European_unresolved_review.csv",
                    "identity_review": "manual_sources/european_trait_identity/european_trait_identity_review.csv",
                    "cache": "raw_lake/European/pubchem_name_cache.json",
                    "online": False,
                    "timeout": 1,
                    "sleep_seconds": 0,
                    "max_queries_per_trait": 1,
                    "max_cids_per_trait": 1,
                    "cache_every": 0,
                },
            )()

            build_european_trait_annotations.build_annotations(args)

            row = next(csv.DictReader((european_dir / "European_trait_annotations.csv").open(encoding="utf-8-sig")))
            self.assertEqual(row["summary_statistics_url"], "http://example.org/GCST000001")
            self.assertEqual(row["pubmed_id"], "123456")
            self.assertEqual(row["paper_title"], "Metabolite GWAS")

    def test_build_annotations_imports_reviewed_identity_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            european_dir = root / "raw_lake" / "European"
            identity_dir = root / "manual_sources" / "european_trait_identity"
            european_dir.mkdir(parents=True, exist_ok=True)
            identity_dir.mkdir(parents=True, exist_ok=True)
            (european_dir / "European.csv").write_text(
                "summaryStatistics,pubmedId,title,accessionId,reportedTrait\n"
                "http://example.org/GCST000002,123456,Metabolite GWAS,GCST000002,Manual compound levels\n",
                encoding="utf-8",
            )
            (identity_dir / "european_trait_identity_review.csv").write_text(
                "accession_id,review_status,identity_scope,reviewed_name,reviewed_pubchem_cid,reviewed_hmdb_id,reviewed_chebi_id,reviewed_inchikey,evidence_source,evidence_url\n"
                "GCST000002,reviewed,strict_identity,Manual compound,12345,HMDB0000001,CHEBI:1,ABCDEFGHIJKLMNOPQRSTUV-N,source supplement,http://example.org/supp\n",
                encoding="utf-8",
            )

            args = type(
                "Args",
                (),
                {
                    "workspace": str(root),
                    "release_id": "mvp_test",
                    "compound_root": Path("compound_match_index"),
                    "european_csv": "raw_lake/European/European.csv",
                    "output": "raw_lake/European/European_trait_annotations.csv",
                    "review_output": "raw_lake/European/European_unresolved_review.csv",
                    "identity_review": "manual_sources/european_trait_identity/european_trait_identity_review.csv",
                    "cache": "raw_lake/European/pubchem_name_cache.json",
                    "online": False,
                    "timeout": 1,
                    "sleep_seconds": 0,
                    "max_queries_per_trait": 1,
                    "max_cids_per_trait": 1,
                    "cache_every": 0,
                },
            )()

            summary = build_european_trait_annotations.build_annotations(args)

            row = next(csv.DictReader((european_dir / "European_trait_annotations.csv").open(encoding="utf-8-sig")))
            self.assertEqual(row["manual_identity_scope"], "strict_identity")
            self.assertEqual(row["reviewed_pubchem_cid"], "12345")
            self.assertEqual(row["reviewed_hmdb_id"], "HMDB0000001")
            self.assertEqual(row["reviewed_chebi_id"], "CHEBI:1")
            self.assertEqual(row["reviewed_inchikey"], "ABCDEFGHIJKLMNOPQRSTUV-N")
            self.assertEqual(row["evidence_source"], "source supplement")
            self.assertEqual(summary["manual_strict_identity_count"], 1)


if __name__ == "__main__":
    unittest.main()
