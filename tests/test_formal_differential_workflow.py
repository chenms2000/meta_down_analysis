import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_formal_differential_workflow.py"
SPEC = importlib.util.spec_from_file_location("run_formal_differential_workflow", SCRIPT)
workflow = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = workflow
SPEC.loader.exec_module(workflow)


def write_trait_table(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {"trait": "GCST1", "group1": "Tumor", "group2": "Adjacent", "direction": "group1_higher", "q_wilcoxon": "0.001", "mean_diff": "0.8"},
        {"trait": "GCST2", "group1": "Tumor", "group2": "Adjacent", "direction": "group1_higher", "q_wilcoxon": "0.020", "mean_diff": "0.6"},
        {"trait": "GCST3", "group1": "Tumor", "group2": "Adjacent", "direction": "group2_higher", "q_wilcoxon": "0.030", "mean_diff": "-0.5"},
        {"trait": "GCST4", "group1": "Tumor", "group2": "Adjacent", "direction": "group2_higher", "q_wilcoxon": "0.200", "mean_diff": "-1.2"},
        {"trait": "GCST5", "group1": "Tumor", "group2": "Adjacent", "direction": "unknown", "q_wilcoxon": "0.500", "mean_diff": "0.1"},
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class FormalDifferentialWorkflowTests(unittest.TestCase):
    def test_workflow_prepares_core_exploratory_and_audit_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            all_degs = root / "TraitScore_GroupDiff_allCelltypes" / "allDEGs" / "trait_score_diff_LUAD_Epi_LUAD_Tumor_vs_Adjacent.csv"
            write_trait_table(all_degs)

            args = workflow.parse_args(
                [
                    "--workspace",
                    str(root),
                    "--release-id",
                    "test_release",
                    "--all-degs",
                    str(all_degs),
                    "--core-max-records",
                    "3",
                    "--core-top-per-group",
                    "2",
                    "--exploratory-max-records",
                    "4",
                    "--exploratory-top-per-group",
                    "2",
                ]
            )
            report = workflow.run_workflow(args)

            self.assertEqual(report["status"], "prepared")
            self.assertEqual(report["workflow_version"], "formal_differential_workflow.v1")
            self.assertEqual(report["audit"]["row_count"], 5)
            self.assertEqual(report["audit"]["q_threshold_counts"]["q_le_0.05"], 3)
            self.assertEqual(report["tiers"]["core_report"]["input_rows"], 3)
            self.assertEqual(report["tiers"]["exploratory_appendix"]["input_rows"], 4)
            self.assertTrue(Path(report["outputs"]["core_input_csv"]).exists())
            self.assertTrue(Path(report["outputs"]["exploratory_input_csv"]).exists())
            self.assertTrue(Path(report["outputs"]["commands_ps1"]).exists())
            self.assertTrue(Path(report["outputs"]["workflow_json"]).exists())
            self.assertTrue(Path(report["outputs"]["workflow_markdown"]).exists())
            markdown = Path(report["outputs"]["workflow_markdown"]).read_text(encoding="utf-8")
            self.assertIn("Core report", markdown)
            self.assertIn("Exploratory appendix", markdown)
            self.assertIn("Full audit", markdown)

    def test_default_significant_subset_is_preferred_for_core(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = root / "TraitScore_GroupDiff_allCelltypes"
            all_degs = base / "allDEGs" / "trait_score_diff_LUAD_Epi_LUAD_Tumor_vs_Adjacent.csv"
            significant = base / "significant" / "trait_score_diff_LUAD_Epi_LUAD_Tumor_vs_Adjacent_significant_q0.05.csv"
            write_trait_table(all_degs)
            write_trait_table(significant)

            args = workflow.parse_args(["--workspace", str(root), "--all-degs", str(all_degs)])
            report = workflow.run_workflow(args)

            self.assertEqual(report["inputs"]["core_source"], str(significant.resolve()))
            self.assertEqual(report["inputs"]["significant"], str(significant.resolve()))


if __name__ == "__main__":
    unittest.main()
