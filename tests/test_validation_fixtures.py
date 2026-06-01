import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = ROOT / "tests" / "fixtures" / "validation"


class ValidationFixtureTests(unittest.TestCase):
    def test_phase15_validation_fixtures_have_contract_fields(self):
        paths = sorted(FIXTURE_DIR.glob("*.json"))
        self.assertGreaterEqual(len(paths), 4)
        for path in paths:
            payload = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(payload["fixture_version"], "phase1.5.20260513")
            self.assertTrue(payload["release_id"].startswith("mvp_"))
            self.assertTrue(payload["purpose"])
            self.assertTrue(payload["cases"])
            for case in payload["cases"]:
                self.assertIn("case_id", case)
                self.assertIn("input", case)
                self.assertIn("expected", case)


if __name__ == "__main__":
    unittest.main()
