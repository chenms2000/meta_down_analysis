import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_real_world_fixtures.py"
SPEC = importlib.util.spec_from_file_location("audit_real_world_fixtures", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)


def registry(path: str):
    return {
        "schema_version": "test",
        "fixtures": [
            {
                "case_id": "cscc_trait_score_csv",
                "path": path,
                "fixture_type": "csv_trait_score",
                "required_for_strict_release": True,
                "min_rows": 1,
                "identity_columns_any": ["trait", "reported_trait"],
                "effect_columns_any": ["cohen_d", "log2FC"],
                "expected_input_mode": "trait_score",
            }
        ],
    }


class RealWorldFixtureAuditTests(unittest.TestCase):
    def test_registry_reports_missing_required_fixture(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = audit.audit_registry(Path(tmp), registry("missing.csv"))

        self.assertEqual(report["status"], "failed")
        self.assertEqual(report["missing_count"], 1)
        self.assertEqual(report["failed_required"][0]["case_id"], "cscc_trait_score_csv")
        self.assertEqual(report["failed_required"][0]["expected_schema"]["expected_input_mode"], "trait_score")

    def test_missing_fixture_reports_candidate_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = root / "backup" / "trait_score_diff_cSCC_Epi_cSCC_Tumor_vs_Adjacent.csv"
            candidate.parent.mkdir(parents=True, exist_ok=True)
            candidate.write_text("trait,cohen_d\nGCST1,1.0\n", encoding="utf-8")

            report = audit.audit_registry(root, registry("registered_missing.csv"))

        fixture = report["failed_required"][0]
        self.assertEqual(fixture["status"], "missing")
        self.assertTrue(fixture["candidate_paths"])
        self.assertIn("cscc", fixture["candidate_paths"][0]["matched_tokens"])
        self.assertIn("recovery_action", fixture)
        self.assertEqual(fixture["recovery_rerun_plan"]["case_id"], "cscc_trait_score_csv")
        self.assertEqual(len(fixture["recovery_rerun_plan"]["validation_sequence"]), 3)
        self.assertIn("strict batch", fixture["recovery_rerun_plan"]["validation_sequence"][1]["purpose"].lower())

    def test_registry_accepts_valid_trait_score_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fixture.csv"
            path.write_text(
                "trait,reported_trait,cohen_d\n"
                "GCST1,Glutamine trait,1.2\n",
                encoding="utf-8",
            )

            report = audit.audit_registry(Path(tmp), registry("fixture.csv"))

        fixture = report["fixtures"][0]
        self.assertEqual(report["status"], "passed")
        self.assertEqual(fixture["status"], "passed")
        self.assertEqual(fixture["observed_input_mode"], "trait_score")
        self.assertEqual(fixture["effect_column"], "cohen_d")

    def test_registry_rejects_wrong_input_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fixture.csv"
            path.write_text(
                "name,cohen_d\n"
                "Glutamine,1.2\n",
                encoding="utf-8",
            )

            report = audit.audit_registry(Path(tmp), registry("fixture.csv"))

        fixture = report["fixtures"][0]
        self.assertEqual(report["status"], "failed")
        self.assertIn("input_mode_mismatch", {row["code"] for row in fixture["errors"]})


if __name__ == "__main__":
    unittest.main()
