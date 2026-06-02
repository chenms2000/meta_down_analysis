import importlib.util
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ADAPTER_SCRIPT = ROOT / "scripts" / "llm_safe_adapter.py"
RUNNER_SCRIPT = ROOT / "scripts" / "run_llm_safe_adapter_validation.py"
FIXTURE_DIR = ROOT / "tests" / "fixtures" / "llm_safe_adapter"
SCHEMA_PATH = ROOT / "config" / "llm_safe_adapter_schema.json"

ADAPTER_SPEC = importlib.util.spec_from_file_location("llm_safe_adapter", ADAPTER_SCRIPT)
adapter = importlib.util.module_from_spec(ADAPTER_SPEC)
assert ADAPTER_SPEC.loader is not None
sys.modules[ADAPTER_SPEC.name] = adapter
ADAPTER_SPEC.loader.exec_module(adapter)

RUNNER_SPEC = importlib.util.spec_from_file_location("run_llm_safe_adapter_validation", RUNNER_SCRIPT)
runner = importlib.util.module_from_spec(RUNNER_SPEC)
assert RUNNER_SPEC.loader is not None
sys.modules[RUNNER_SPEC.name] = runner
RUNNER_SPEC.loader.exec_module(runner)


class LLMSafeAdapterTests(unittest.TestCase):
    def load_fixture(self):
        path = FIXTURE_DIR / "llm_safe_adapter_cases.json"
        return json.loads(path.read_text(encoding="utf-8"))

    def test_local_adapter_outputs_are_source_bound_and_guarded(self):
        fixture = self.load_fixture()
        local_cases = [case for case in fixture["cases"] if "adapter_output" not in case]
        self.assertGreaterEqual(len(local_cases), 4)
        for case in local_cases:
            with self.subTest(case_id=case["case_id"]):
                output = adapter.build_local_adapter_output(case["input_pack"])
                self.assertTrue(output["guard"]["passed"], output["guard"].get("issues"))
                self.assertEqual(output["status"], case["expected"]["status"])
                for field in case["expected"].get("required_outputs", []):
                    self.assertIn(field, output)
                self.assertGreater(output["guard"]["language_block_count"], 0)
                if output["status"] == "ok":
                    section_ids = {row.get("section_id") for row in output["narrative_summary"]["sections"]}
                    self.assertIn("calibration_boundaries", section_ids)
                    if case["input_pack"].get("analysis_pack", {}).get("input_summary"):
                        self.assertIn("precision_input_parsing", section_ids)
                self.assertEqual(
                    output["determinism"]["guard_policy_hash"],
                    adapter.guard_policy_hash(),
                )

    def test_guard_blocks_adversarial_mutation_output(self):
        fixture = self.load_fixture()
        case = next(row for row in fixture["cases"] if row["case_id"] == "adversarial_prompt_mutation_attempt")
        guard = adapter.guard_adapter_output(case["input_pack"], case["adapter_output"])
        codes = {issue["code"] for issue in guard["issues"]}
        self.assertFalse(guard["passed"])
        self.assertTrue(set(case["expected"]["required_issue_codes"]).issubset(codes))

    def test_validation_runner_is_deterministic(self):
        fixtures = runner.load_fixtures(FIXTURE_DIR)
        first = runner.build_report(fixtures, FIXTURE_DIR)
        second = runner.build_report(fixtures, FIXTURE_DIR)
        self.assertEqual(first["report_hash"], second["report_hash"])
        self.assertEqual(first["summary"]["gate_status"], "passed")
        self.assertEqual(first["summary"]["failed"], 0)
        self.assertEqual(first["summary"]["passed"], 5)

    def test_schema_matches_adapter_constants(self):
        schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
        self.assertEqual(schema["input_contract_version"], adapter.INPUT_CONTRACT_VERSION)
        self.assertEqual(schema["output_contract_version"], adapter.OUTPUT_CONTRACT_VERSION)
        self.assertEqual(set(schema["allowed_read_endpoints"]), adapter.ALLOWED_READ_ENDPOINTS)
        self.assertEqual(set(schema["source_types"]), adapter.SOURCE_TYPES)
        self.assertTrue(set(schema["allowed_output_fields"]).issubset(adapter.EXPLANATION_FIELDS))

    def test_external_llm_backend_requires_explicit_enablement(self):
        fixture = self.load_fixture()
        input_pack = fixture["cases"][0]["input_pack"]
        backend = adapter.ExternalLLMNarratorBackend(
            adapter.ExternalLLMConfig(
                enabled=False,
                provider="openai_compatible",
                endpoint="https://llm.example.test/v1/chat/completions",
                model="test-model",
                api_key="secret",
            ),
            transport=lambda _config, _payload: {},
        )
        with self.assertRaises(adapter.LLMAdapterError) as raised:
            backend.generate(input_pack)
        self.assertEqual(raised.exception.code, "external_llm_disabled")

    def test_external_llm_backend_wraps_text_and_keeps_audit_sanitized(self):
        fixture = self.load_fixture()
        input_pack = fixture["cases"][0]["input_pack"]

        def fake_transport(config, payload):
            self.assertEqual(config.model, "test-model")
            self.assertEqual(payload["temperature"], 0)
            self.assertNotIn("response_format", payload)
            return {
                "id": "chatcmpl_fixture",
                "choices": [
                    {
                        "message": {
                            "content": "The external narrator summarizes only the frozen analysis pack and preserves validation boundaries."
                        }
                    }
                ],
            }

        backend = adapter.ExternalLLMNarratorBackend(
            adapter.ExternalLLMConfig(
                enabled=True,
                provider="openai_compatible",
                endpoint="https://llm.example.test/v1/chat/completions",
                model="test-model",
                api_key="secret",
            ),
            transport=fake_transport,
        )
        result = backend.generate(input_pack)
        guard = adapter.guard_adapter_output(input_pack, result.output)
        self.assertTrue(guard["passed"], guard.get("issues"))
        self.assertEqual(result.audit["backend"], "external_llm")
        self.assertEqual(result.audit["model"], "test-model")
        self.assertIn("text_fixed_contract", result.audit["adapter_version"])
        self.assertEqual(result.output["narrative_summary"]["sections"][0]["section_id"], "external_text_summary")
        self.assertIn("prompt_hash", result.audit)
        self.assertNotIn("secret", json.dumps(result.audit, sort_keys=True))

    def test_prediction_claims_must_bind_to_structured_prediction_rows(self):
        fixture = self.load_fixture()
        input_pack = fixture["cases"][0]["input_pack"]
        input_pack["analysis_pack"]["structured_prediction"] = {
            "mode": "generalized",
            "target_predictions": [
                {
                    "prediction_id": "target_hk1_prediction",
                    "result_id": "target_hk1",
                    "display_name": "HK1",
                    "confidence_tier": "medium",
                    "evidence_refs": [{"support_uid": "litsup_hk1_cancer"}],
                    "claim_refs": {"traceability_passed": True, "evidence_ref_uids": ["litsup_hk1_cancer"]},
                    "appendix": False,
                    "research_only": False,
                }
            ],
        }
        output = {
            "contract_version": adapter.OUTPUT_CONTRACT_VERSION,
            "adapter_version": "fixture",
            "status": "ok",
            "narrative_summary": {
                "sections": [
                    {
                        "section_id": "target_claim",
                        "reported_display_names": ["HK1"],
                        "reported_prediction_ids": ["target_hk1_prediction"],
                        "prediction_refs": [
                            {
                                "prediction_id": "target_hk1_prediction",
                                "confidence_tier": "medium",
                                "evidence_refs": [{"support_uid": "litsup_hk1_cancer"}],
                                "claim_refs": {"traceability_passed": True, "evidence_ref_uids": ["litsup_hk1_cancer"]},
                                "appendix": False,
                                "research_only": False,
                            }
                        ],
                        "text": "HK1 is a medium-confidence target prediction in the structured prediction JSON.",
                        "source_refs": [
                            {"source_type": "analysis_pack", "ref_id": "analysis_pack:root", "path": "$.analysis_pack"}
                        ],
                    }
                ]
            },
        }
        guard = adapter.guard_adapter_output(input_pack, output)
        self.assertTrue(guard["passed"], guard.get("issues"))

    def test_unsupported_prediction_claim_is_rejected(self):
        fixture = self.load_fixture()
        input_pack = fixture["cases"][0]["input_pack"]
        input_pack["analysis_pack"]["structured_prediction"] = {
            "mode": "generalized",
            "target_predictions": [
                {
                    "prediction_id": "target_hk1_prediction",
                    "display_name": "HK1",
                    "confidence_tier": "medium",
                    "evidence_refs": [{"support_uid": "litsup_hk1_cancer"}],
                    "claim_refs": {"traceability_passed": True},
                    "appendix": False,
                    "research_only": False,
                }
            ],
        }
        output = {
            "contract_version": adapter.OUTPUT_CONTRACT_VERSION,
            "adapter_version": "fixture",
            "status": "ok",
            "narrative_summary": {
                "sections": [
                    {
                        "section_id": "unsupported_target_claim",
                        "reported_display_names": ["GLS"],
                        "reported_prediction_ids": ["target_gls_prediction"],
                        "claim": "GLS is a high-confidence target prediction.",
                        "source_refs": [
                            {"source_type": "analysis_pack", "ref_id": "analysis_pack:root", "path": "$.analysis_pack"}
                        ],
                    }
                ]
            },
        }
        guard = adapter.guard_adapter_output(input_pack, output)
        codes = {issue["code"] for issue in guard["issues"]}
        self.assertFalse(guard["passed"])
        self.assertIn("prediction_claim_without_structured_prediction_refs", codes)

    def test_external_llm_parser_extracts_json_from_wrapped_text(self):
        wrapped = {
            "choices": [
                {
                    "message": {
                        "content": (
                            "Here is the JSON:\n```json\n"
                            + json.dumps(
                                {
                                    "contract_version": "llm_safe_adapter.output.v1",
                                    "adapter_version": "fake_external_llm",
                                    "status": "insufficient_evidence",
                                    "insufficient_evidence": {
                                        "reason": "No supported evidence was present.",
                                        "source_refs": [
                                            {
                                                "source_type": "release",
                                                "ref_id": "release:mvp_20260513T002254",
                                                "path": "$.release.release_id",
                                            }
                                        ],
                                    },
                                }
                            )
                            + "\n```"
                        )
                    }
                }
            ]
        }
        parsed = adapter.parse_external_llm_output(wrapped)
        self.assertEqual(parsed["status"], "insufficient_evidence")
        self.assertIn("insufficient_evidence", parsed)

    def test_external_llm_empty_content_reports_shape(self):
        response = {
            "choices": [
                {
                    "finish_reason": "stop",
                    "message": {"role": "assistant", "content": ""},
                }
            ]
        }
        with self.assertRaises(adapter.LLMAdapterError) as raised:
            adapter.parse_external_llm_output(response)
        self.assertEqual(raised.exception.code, "external_llm_empty_content")
        self.assertIn("message_keys", raised.exception.detail)

    def test_external_text_backend_wraps_plain_language_output(self):
        fixture = self.load_fixture()
        input_pack = fixture["cases"][0]["input_pack"]

        def fake_transport(config, payload):
            self.assertEqual(config.model, "test-model")
            self.assertNotIn("response_format", payload)
            return {
                "choices": [
                    {
                        "message": {
                            "content": "这份结果可以先从匹配情况和可追溯证据两方面解读；不应把它当作诊断结论。"
                        }
                    }
                ]
            }

        backend = adapter.ExternalLLMTextNarratorBackend(
            adapter.ExternalLLMConfig(
                enabled=True,
                provider="openai_compatible",
                endpoint="https://llm.example.test/v1/chat/completions",
                model="test-model",
                api_key="secret",
            ),
            transport=fake_transport,
        )
        result = backend.generate(input_pack)
        self.assertEqual(result.output["status"], "ok")
        self.assertIn("narrative_summary", result.output)
        guard = adapter.guard_adapter_output(input_pack, result.output)
        self.assertTrue(guard["passed"], guard.get("issues"))
        self.assertEqual(result.audit["backend"], "external_text")

    def test_default_transport_uses_urlopen_without_proxy(self):
        request_payload = {"model": "test-model", "messages": [], "temperature": 0}
        config = adapter.ExternalLLMConfig(
            enabled=True,
            provider="openai_compatible",
            endpoint="https://llm.example.test/v1/chat/completions",
            model="test-model",
            api_key="secret",
        )

        class FakeResponse:
            def __enter__(self):
                return self

            def __exit__(self, *_args):
                return False

            def read(self):
                return b'{"choices":[{"message":{"content":"ok"}}]}'

        original_urlopen = adapter.urllib_request.urlopen
        try:
            calls = []

            def fake_urlopen(req, timeout):
                calls.append((req.full_url, timeout))
                return FakeResponse()

            adapter.urllib_request.urlopen = fake_urlopen
            response = adapter.default_openai_compatible_transport(config, request_payload)
        finally:
            adapter.urllib_request.urlopen = original_urlopen

        self.assertEqual(response["choices"][0]["message"]["content"], "ok")
        self.assertEqual(calls, [("https://llm.example.test/v1/chat/completions", 30.0)])


if __name__ == "__main__":
    unittest.main()
