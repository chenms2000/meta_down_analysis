import importlib.util
import sys
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_source_heldout_validation.py"
SPEC = importlib.util.spec_from_file_location("run_source_heldout_validation", SCRIPT)
runner = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


class SourceHeldoutValidationTests(unittest.TestCase):
    def test_infer_source_specs_partitions_unique_literature_family(self):
        labels = pd.DataFrame(
            {
                "candidate_uid": [f"lit{i}" for i in range(12)] + [f"drug{i}" for i in range(4)],
                "source_kind": ["literature_relation"] * 12 + ["drug_target_overlay"] * 4,
                "source_name": ["literature:local"] * 12 + ["ChEMBL"] * 4,
                "source_partition": ["literature:local:partition_0"] * 6
                + ["literature:local:partition_1"] * 6
                + ["ChEMBL"] * 4,
                "weak_label": ["weak_positive"] * 12 + ["weak_positive"] * 4,
            }
        )

        specs = runner.infer_source_specs(labels, min_positive_rows=2)

        self.assertIn({"field": "source_partition", "value": "literature:local:partition_0"}, specs)
        self.assertIn({"field": "source_partition", "value": "literature:local:partition_1"}, specs)
        self.assertNotIn({"field": "source_name", "value": "literature:local"}, specs)

    def test_infer_source_specs_keeps_source_name_when_family_remains(self):
        labels = pd.DataFrame(
            {
                "candidate_uid": [f"drug{i}" for i in range(8)],
                "source_kind": ["drug_target_overlay"] * 8,
                "source_name": ["ChEMBL"] * 4 + ["DGIdb"] * 4,
                "source_partition": ["ChEMBL"] * 4 + ["DGIdb"] * 4,
                "weak_label": ["weak_positive"] * 8,
            }
        )

        specs = runner.infer_source_specs(labels, min_positive_rows=2)

        self.assertIn({"field": "source_name", "value": "ChEMBL"}, specs)
        self.assertIn({"field": "source_name", "value": "DGIdb"}, specs)

    def test_source_family_diagnostics_detects_missing_heldout_family(self):
        train = pd.DataFrame(
            {
                "source_kind": ["drug_target_overlay", "shuffled_negative", "drug_target_overlay"],
                "candidate_uid": ["a", "b", "c"],
            }
        )
        heldout = pd.DataFrame(
            {
                "source_kind": ["literature_relation", "literature_relation"],
                "candidate_uid": ["d", "e"],
            }
        )

        diagnostics = runner.source_family_diagnostics(train, heldout)

        self.assertFalse(diagnostics["same_family_available"])
        self.assertEqual(diagnostics["same_family_train_rows"], 0)
        self.assertEqual(diagnostics["missing_heldout_source_kinds_in_train"], ["literature_relation"])

    def test_feature_transfer_diagnostics_reports_shared_variation(self):
        train = pd.DataFrame(
            {
                "candidate_uid": ["a", "b", "c"],
                "p_literature": [0.0, 0.5, 1.0],
                "subject_type": ["gene", "gene", "drug"],
                "object_type": ["disease", "pathway", "pathway"],
                "predicate": ["a", "a", "b"],
                "relation_family": ["x", "x", "y"],
                "subject_object_family": ["g:d", "g:p", "d:p"],
            }
        )
        heldout = pd.DataFrame(
            {
                "candidate_uid": ["d", "e", "f"],
                "p_literature": [0.2, 0.4, 0.6],
                "subject_type": ["gene", "metabolite", "metabolite"],
                "object_type": ["disease", "disease", "pathway"],
                "predicate": ["a", "c", "c"],
                "relation_family": ["x", "z", "z"],
                "subject_object_family": ["g:d", "m:d", "m:p"],
            }
        )
        feature_columns = [
            "p_literature",
            "subject_type_gene",
            "subject_type_drug",
            "object_type_disease",
            "predicate_a",
        ]

        diagnostics = runner.feature_transfer_diagnostics(train, heldout, feature_columns)

        self.assertGreaterEqual(diagnostics["shared_nonconstant_feature_count"], 1)
        self.assertIn("p_literature", diagnostics["shared_nonconstant_features_sample"])


if __name__ == "__main__":
    unittest.main()
