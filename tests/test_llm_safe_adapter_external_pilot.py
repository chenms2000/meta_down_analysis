import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_llm_safe_adapter_external_pilot.py"

SPEC = importlib.util.spec_from_file_location("run_llm_safe_adapter_external_pilot", SCRIPT)
pilot = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = pilot
SPEC.loader.exec_module(pilot)


class LLMSafeAdapterExternalPilotTests(unittest.TestCase):
    def test_readiness_requires_explicit_external_config(self):
        config = pilot.ExternalLLMConfig(enabled=False, provider="openai_compatible", model="", api_key="")
        ready = pilot.readiness(config, run_external=True)
        self.assertFalse(ready["ready"])
        self.assertIn("external_llm_not_enabled", ready["missing"])
        self.assertIn("missing_model", ready["missing"])
        self.assertIn("missing_api_key", ready["missing"])

    def test_score_not_run_is_non_failing_readiness_report(self):
        regression = {
            "summary": {"adversarial_probe_failed": 0},
            "cases": [{"external": {"status": "not_run"}}, {"external": {"status": "not_run"}}],
        }
        score = pilot.score_pilot(regression, {"run_external_requested": False, "ready": False}, min_external_pass_rate=1.0)
        self.assertEqual(score["gate_status"], "not_run")
        self.assertTrue(score["gate_passed"])
        self.assertEqual(score["external_status_counts"], {"not_run": 2})

    def test_score_live_external_requires_pass_rate(self):
        regression = {
            "summary": {"adversarial_probe_failed": 0},
            "cases": [{"external": {"status": "passed"}}, {"external": {"status": "blocked_by_guard", "guard_issue_codes": ["overstrong_claim"]}}],
        }
        score = pilot.score_pilot(regression, {"run_external_requested": True, "ready": True}, min_external_pass_rate=1.0)
        self.assertEqual(score["gate_status"], "needs_prompt_iteration")
        self.assertFalse(score["gate_passed"])
        self.assertEqual(score["external_guard_issue_counts"], {"overstrong_claim": 1})


if __name__ == "__main__":
    unittest.main()
