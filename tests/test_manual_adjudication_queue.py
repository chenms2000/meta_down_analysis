import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_manual_adjudication_queue.py"
SPEC = importlib.util.spec_from_file_location("build_manual_adjudication_queue", SCRIPT)
queue = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = queue
SPEC.loader.exec_module(queue)


class ManualAdjudicationQueueTests(unittest.TestCase):
    def test_builds_sampled_pending_review_queue(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_csv = root / "compound_name_audit_candidates.csv"
            rows = [
                {
                    "risk_score": "90",
                    "name_key": "ambiguous lipid",
                    "recommended_action": "require_lipid_identifier",
                    "recommendation": "provide stable ID",
                    "risk_reasons": "name_maps_to_multiple_metabolites;lipid_shorthand_or_chain_notation",
                    "uid_count": "3",
                    "row_count": "5",
                    "canonical_examples": "x",
                    "authority_prefixes": "HMDB",
                    "match_fields": "synonym",
                    "raw_examples": "raw x",
                },
                {
                    "risk_score": "88",
                    "name_key": "d/l isomer",
                    "recommended_action": "require_specific_isomer_or_identifier",
                    "recommendation": "provide isomer",
                    "risk_reasons": "stereo_or_isomer_under_specified",
                    "uid_count": "2",
                    "row_count": "2",
                    "canonical_examples": "y",
                    "authority_prefixes": "CHEBI",
                    "match_fields": "name",
                    "raw_examples": "raw y",
                },
                {
                    "risk_score": "30",
                    "name_key": "comma,name",
                    "recommended_action": "quote_delimiter_sensitive_name",
                    "recommendation": "quote it",
                    "risk_reasons": "delimiter_sensitive_name",
                    "uid_count": "1",
                    "row_count": "1",
                    "canonical_examples": "z",
                    "authority_prefixes": "HMDB",
                    "match_fields": "name",
                    "raw_examples": "raw z",
                },
            ]
            with input_csv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)

            exit_code = queue.main(["--input-csv", str(input_csv), "--sample-per-bucket", "2"])
            self.assertEqual(exit_code, 0)

            output_csv = input_csv.parent / "manual_adjudication_queue.csv"
            output_json = input_csv.parent / "manual_adjudication_queue_summary.json"
            payload = json.loads(output_json.read_text(encoding="utf-8"))
            self.assertEqual(payload["selected_rows"], 3)
            self.assertEqual(payload["interpretation"], "This queue is a sampled manual follow-up target, not a completed adjudication result.")
            review_rows = list(csv.DictReader(output_csv.open("r", encoding="utf-8")))
            self.assertTrue(all(row["adjudication_status"] == "pending" for row in review_rows))
            self.assertIn("reviewer_entity_id", review_rows[0])


if __name__ == "__main__":
    unittest.main()
