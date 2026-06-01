import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "triage_compound_name_audit.py"
SPEC = importlib.util.spec_from_file_location("triage_compound_name_audit", SCRIPT)
triage = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = triage
SPEC.loader.exec_module(triage)


class CompoundNameTriageTests(unittest.TestCase):
    def test_triages_risk_rows_into_release_policies(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_csv = root / "compound_name_audit_candidates.csv"
            rows = [
                {
                    "risk_score": "90",
                    "name_key": "ambiguous lipid",
                    "recommended_action": "require_lipid_identifier",
                    "recommendation": "provide a stable ID",
                    "risk_reasons": "name_maps_to_multiple_metabolites;lipid_shorthand_or_chain_notation",
                    "uid_count": "3",
                    "canonical_examples": "x",
                },
                {
                    "risk_score": "30",
                    "name_key": "comma,name",
                    "recommended_action": "quote_delimiter_sensitive_name",
                    "recommendation": "quote it",
                    "risk_reasons": "delimiter_sensitive_name",
                    "uid_count": "1",
                    "canonical_examples": "y",
                },
                {
                    "risk_score": "10",
                    "name_key": "low rank synonym",
                    "recommended_action": "review_low_precision_surface",
                    "recommendation": "monitor",
                    "risk_reasons": "synonym_only_low_rank_surface",
                    "uid_count": "1",
                    "canonical_examples": "z",
                },
            ]
            with input_csv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)

            exit_code = triage.main(["--input-csv", str(input_csv), "--sample-per-policy", "1"])
            self.assertEqual(exit_code, 0)

            output_json = input_csv.parent / "compound_name_triage_summary.json"
            output_md = input_csv.parent / "compound_name_triage_summary.md"
            payload = json.loads(output_json.read_text(encoding="utf-8"))
            self.assertEqual(payload["total_risk_rows"], 3)
            self.assertEqual(payload["policy_counts"]["auto_abstain_name_only"], 1)
            self.assertEqual(payload["policy_counts"]["format_warning"], 1)
            self.assertEqual(payload["policy_counts"]["monitor"], 1)
            self.assertIn("Manual adjudication remains a future optimization step", output_md.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
