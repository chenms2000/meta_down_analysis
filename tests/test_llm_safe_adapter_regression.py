import importlib.util
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER_SCRIPT = ROOT / "scripts" / "run_llm_safe_adapter_regression.py"
FIXTURE_PATH = ROOT / "tests" / "fixtures" / "llm_safe_adapter" / "llm_safe_adapter_cases.json"

SPEC = importlib.util.spec_from_file_location("run_llm_safe_adapter_regression", RUNNER_SCRIPT)
runner = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


class LLMSafeAdapterRegressionTests(unittest.TestCase):
    def make_fixture(self):
        base = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))
        input_pack = base["cases"][0]["input_pack"]
        return {
            "fixture_version": "llm_safe_adapter_regression.test",
            "release_id": "mvp_test",
            "defaults": {},
            "cases": [
                {
                    "case_id": "fixture_input_pack",
                    "question": "Summarize safely.",
                    "input_pack": input_pack,
                    "expected": {
                        "required_output_fields": ["narrative_summary", "evidence_digest", "question_to_spec"],
                        "min_language_blocks": 2,
                        "guard_passed": True,
                    },
                }
            ],
            "adversarial_probes": [
                {
                    "probe_id": "fake_score",
                    "output": {
                        "contract_version": "llm_safe_adapter.output.v1",
                        "adapter_version": "bad",
                        "status": "ok",
                        "narrative_summary": {
                            "sections": [
                                {
                                    "section_id": "bad",
                                    "text": "I resolved metabolite_fake and updated p_final score to 0.99 for edge_fake.",
                                    "source_refs": [
                                        {"source_type": "evidence", "ref_id": "fake_ref", "path": "$.evidence.support[99]"}
                                    ],
                                }
                            ]
                        },
                        "p_final": 0.99,
                    },
                    "expected_issue_codes": ["forbidden_output_field", "source_ref_not_in_input", "entity_resolution_decision"],
                }
            ],
        }

    def test_regression_report_is_deterministic_without_external_llm(self):
        fixture = self.make_fixture()
        config = runner.ExternalLLMConfig(enabled=False, provider="openai_compatible", model="test", api_key="")
        first = runner.build_report(None, fixture, Path("fixture.json"), config, run_external=False)
        second = runner.build_report(None, fixture, Path("fixture.json"), config, run_external=False)
        self.assertEqual(first["report_hash"], second["report_hash"])
        self.assertEqual(first["summary"]["gate_status"], "passed")
        self.assertEqual(first["summary"]["local_passed"], 1)
        self.assertEqual(first["summary"]["adversarial_probe_passed"], 1)
        self.assertEqual(first["cases"][0]["external"]["status"], "not_run")

    def test_probe_requires_expected_issue_codes(self):
        fixture = self.make_fixture()
        input_pack = fixture["cases"][0]["input_pack"]
        probe = fixture["adversarial_probes"][0]
        result = runner.evaluate_probe(input_pack, probe)
        self.assertEqual(result["status"], "passed")
        self.assertFalse(result["guard_passed"])
        self.assertIn("forbidden_output_field", result["guard_issue_codes"])


if __name__ == "__main__":
    unittest.main()
