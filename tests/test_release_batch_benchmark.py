import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_release_batch_benchmark.py"
SPEC = importlib.util.spec_from_file_location("run_release_batch_benchmark", SCRIPT)
runner = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


class ReleaseBatchBenchmarkTests(unittest.TestCase):
    def test_optimization_signals_classify_error_modes(self):
        scenarios = [
            {
                "case_id": "ambiguous_case",
                "status": "passed",
                "benchmark_tags": {
                    "cancer_ids": ["brca"],
                    "theme_ids": ["glycolysis_lactate"],
                    "input_modes": ["compound_name"],
                },
                "input_summary": {"input_count": 2, "ambiguous_count": 2, "unmatched_count": 0},
                "expected_identity_behavior": {},
                "scenario_expected_metrics": {"expected_count": 2, "expected_matched": 1},
                "conclusion_evaluation": {
                    "metrics": {
                        "gold_positive_count": 4,
                        "gold_positive_matched": 1,
                        "negative_trap_hits": 0,
                        "unsupported_top_claim_rate": 0.5,
                        "context_accuracy": 0.75,
                        "evidence_usefulness": 0.0,
                    },
                    "matched_gold": [{"gold_id": "gold_a", "status": "partial_match"}],
                    "per_gold_status_counts": {"context_mismatch": 1},
                    "bug_queue": [],
                },
            },
            {
                "case_id": "missing_fixture",
                "status": "missing",
                "path": "missing.csv",
            },
            {
                "case_id": "trap_case",
                "status": "failed",
                "input_summary": {"input_count": 1, "ambiguous_count": 0, "unmatched_count": 0},
                "conclusion_evaluation": {
                    "metrics": {
                        "gold_positive_count": 1,
                        "gold_positive_matched": 1,
                        "negative_trap_hits": 1,
                    },
                    "bug_queue": [{"code": "negative_trap_hit"}],
                },
            },
        ]

        summary = runner.aggregate_optimization_signals(scenarios)
        by_code = {row["code"]: row for row in summary["by_code"]}

        self.assertIn("identity_ambiguity_review", by_code)
        self.assertIn("false_negative_proxy", by_code)
        self.assertIn("unsupported_claim_risk", by_code)
        self.assertIn("context_drift_review", by_code)
        self.assertIn("gold_context_mismatch", by_code)
        self.assertIn("evidence_not_effective", by_code)
        self.assertIn("fixture_missing", by_code)
        self.assertIn("false_positive_trap", by_code)
        self.assertIn("negative_trap_hit", by_code)
        self.assertIn("resolver", by_code["identity_ambiguity_review"]["optimization_targets"])
        self.assertIn("evidence_scoring", by_code["unsupported_claim_risk"]["optimization_targets"])
        self.assertEqual(by_code["fixture_missing"]["cases"], ["missing_fixture"])
        by_target = {row["bucket"]: row for row in summary["by_target"]}
        self.assertIn("resolver", by_target)
        self.assertIn("negative_trap_block_rule", by_target)
        by_dimension = {row["bucket"]: row for row in summary["by_dimension"]}
        self.assertIn("cancers:brca", by_dimension)
        self.assertIn("themes:glycolysis_lactate", by_dimension)
        self.assertIn("input_modes:compound_name", by_dimension)

    def test_optimization_signals_flag_prediction_quality_overclaims(self):
        case = {
            "case_id": "crc_trait_score_csv",
            "status": "passed",
            "benchmark_tags": {
                "cancer_ids": ["coad_read"],
                "theme_ids": ["sphingolipid"],
                "input_modes": ["trait_score"],
            },
            "input_summary": {
                "analysis_mode": "two_group_trait_comparison_table",
                "input_count": 80,
                "matched_count": 25,
                "ambiguous_count": 30,
                "unmatched_count": 25,
            },
            "prediction_quality": {
                "high_confidence_trait_score_weak_support_count": 1,
                "context_mismatch_prediction_count": 1,
                "high_or_medium_context_mismatch_prediction_count": 0,
            },
            "conclusion_evaluation": {
                "metrics": {
                    "negative_trap_hits": 0,
                    "unsupported_top_claim_rate": 0.0,
                    "context_accuracy": 1.0,
                },
                "bug_queue": [],
                "per_gold_status_counts": {},
            },
        }

        summary = runner.aggregate_optimization_signals([case])
        by_code = {row["code"]: row for row in summary["by_code"]}

        self.assertIn("trait_score_confidence_overclaim", by_code)
        self.assertEqual(by_code["trait_score_confidence_overclaim"]["severity"], "error")
        self.assertIn("ranking_calibration", by_code["trait_score_confidence_overclaim"]["optimization_targets"])
        self.assertIn("context_drift_review", by_code)
        by_dimension = {row["bucket"]: row for row in summary["by_dimension"]}
        self.assertIn("input_modes:trait_score", by_dimension)

    def test_scenario_timeout_enters_error_signal_targets(self):
        case = {
            "case_id": "slow_case",
            "status": "failed",
            "benchmark_tags": {
                "cancer_ids": ["brca"],
                "theme_ids": ["glycolysis_lactate"],
                "input_modes": ["compound_name"],
            },
            "conclusion_evaluation": {
                "metrics": {"negative_trap_hits": 0},
                "bug_queue": [{"code": "scenario_timeout"}],
            },
        }

        signals = runner.aggregate_optimization_signals([case])
        by_code = {row["code"]: row for row in signals["by_code"]}

        self.assertIn("scenario_timeout", by_code)
        self.assertEqual(by_code["scenario_timeout"]["severity"], "error")
        self.assertIn("analysis_runtime", by_code["scenario_timeout"]["optimization_targets"])

    def test_scenario_timeout_result_preserves_traceability_failure(self):
        result = runner.scenario_timeout_result(
            {"case_id": "slow_case", "records": [{"name": "Glucose"}], "expected_gold_ids": ["gold_a"]},
            1,
            1.25,
        )

        self.assertEqual(result["status"], "failed")
        self.assertTrue(result["timed_out"])
        self.assertEqual(result["traceability"]["evidence_assertion_count"], 0)
        self.assertEqual(result["conclusion_evaluation"]["bug_queue"][0]["code"], "scenario_timeout")

    def test_markdown_includes_optimization_target_and_dimension_matrix(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report.md"
            report = {
                "release_id": "release",
                "created_at_utc": "2026-05-23T00:00:00Z",
                "status": {"passed": False},
                "release_profile": {"profile": "strict"},
                "progress_log": "progress.jsonl",
                "progress_events": [
                    {"event": "scenario_started", "case_id": "case_a"},
                    {"event": "scenario_completed", "case_id": "case_a", "status": "failed", "duration_seconds": 1.2, "error_type": "scenario_timeout"},
                ],
                "scenarios": [
                    {
                        "case_id": "crc_trait_score_csv",
                        "status": "passed",
                        "input_summary": {"matched_count": 25, "ambiguous_count": 30},
                        "prediction_quality": {
                            "analysis_mode": "two_group_trait_comparison_table",
                            "high_confidence_trait_score_weak_support_count": 1,
                            "context_mismatch_prediction_count": 1,
                            "high_or_medium_context_mismatch_prediction_count": 0,
                            "high_confidence_trait_score_weak_support": [
                                {
                                    "rank": 1,
                                    "display_name": "Sphingolipid catabolism",
                                    "overlap_count": 1,
                                    "significant_support_count": 1,
                                    "literature_support_count": 0,
                                }
                            ],
                            "context_mismatch_predictions": [
                                {
                                    "ranking": "disease_rankings",
                                    "rank": 1,
                                    "display_name": "melanoma",
                                    "confidence_tier": "low",
                                    "calibration_status": "context_mismatch",
                                }
                            ],
                        },
                    }
                ],
                "optimization_signals": {
                    "by_code": [
                        {
                            "code": "identity_ambiguity_review",
                            "severity": "warning",
                            "case_count": 1,
                            "cases": ["case_a"],
                            "optimization_targets": ["resolver", "input_parser"],
                            "recommended_action": "review resolver",
                        }
                    ],
                    "by_target": [
                        {
                            "bucket": "resolver",
                            "max_severity": "warning",
                            "codes": ["identity_ambiguity_review"],
                            "cases": ["case_a"],
                        }
                    ],
                    "by_dimension": [
                        {
                            "bucket": "input_modes:compound_name",
                            "max_severity": "warning",
                            "codes": ["identity_ambiguity_review"],
                            "optimization_targets": ["resolver"],
                            "cases": ["case_a"],
                        }
                    ],
                },
            }

            runner.write_markdown(path, report)

            text = path.read_text(encoding="utf-8")
        self.assertIn("Optimization targets:", text)
        self.assertIn("Progress log:", text)
        self.assertIn("scenario_completed", text)
        self.assertIn("scenario_timeout", text)
        self.assertIn("| resolver | warning | identity_ambiguity_review | 1 |", text)
        self.assertIn("Optimization dimensions:", text)
        self.assertIn("input_modes:compound_name", text)
        self.assertIn("Prediction Quality Review", text)
        self.assertIn("crc_trait_score_csv", text)
        self.assertIn("Sphingolipid catabolism", text)
        self.assertIn("melanoma", text)

    def test_expected_ambiguity_does_not_enter_warning_queue(self):
        scenarios = [
            {
                "case_id": "expected_ratio_abstention",
                "status": "passed",
                "input_summary": {"input_count": 1, "ambiguous_count": 1, "unmatched_count": 0},
                "expected_identity_behavior": {"allow_ambiguity": True, "reason": "ratio"},
                "scenario_expected_metrics": {"expected_count": 0, "expected_matched": 0},
                "conclusion_evaluation": {
                    "metrics": {"negative_trap_hits": 0},
                    "matched_gold": [],
                    "per_gold_status_counts": {},
                    "bug_queue": [],
                },
            }
        ]

        summary = runner.aggregate_optimization_signals(scenarios)

        self.assertNotIn("identity_ambiguity_review", {row["code"] for row in summary["by_code"]})

    def test_strict_release_profile_enables_hard_gates(self):
        args = runner.parse_args(["--strict-release"])

        applied = runner.apply_release_profile(args)
        summary = runner.release_profile_summary(applied)

        self.assertTrue(applied.require_real_fixtures)
        self.assertTrue(applied.require_heldout_rerun)
        self.assertTrue(applied.rerun_heldout)
        self.assertEqual(summary["profile"], "strict")
        self.assertTrue(summary["strict_release"])

    def test_exploratory_release_profile_preserves_explicit_flags(self):
        args = runner.parse_args(["--require-real-fixtures"])

        applied = runner.apply_release_profile(args)
        summary = runner.release_profile_summary(applied)

        self.assertEqual(summary["profile"], "exploratory")
        self.assertTrue(applied.require_real_fixtures)
        self.assertFalse(applied.require_heldout_rerun)
        self.assertFalse(applied.rerun_heldout)

    def test_gold_coverage_gap_merges_into_report_signals(self):
        scenario_signals = {
            "by_case": {},
            "by_code": [
                {
                    "code": "fixture_missing",
                    "severity": "warning",
                    "case_count": 1,
                    "cases": ["cscc_trait_score_csv"],
                    "recommended_action": "Restore fixture.",
                }
            ],
        }
        gold_coverage = {
            "status": "gaps_present",
            "missing_cell_count": 1,
            "missing_cells": [{"cancer_id": "pan_cancer", "theme_id": "nucleotide_metabolism"}],
            "negative_control_requirement": {"status": "covered"},
        }

        merged = runner.merge_report_signals(scenario_signals, [runner.gold_coverage_signal(gold_coverage)])
        by_code = {row["code"]: row for row in merged["by_code"]}

        self.assertIn("fixture_missing", by_code)
        self.assertIn("gold_coverage_gap", by_code)
        self.assertEqual(by_code["gold_coverage_gap"]["cases"], ["pan_cancer:nucleotide_metabolism"])
        by_target = {row["bucket"]: row for row in merged["by_target"]}
        self.assertIn("gold_standard", by_target)

    def test_gold_quality_gap_merges_into_report_signals(self):
        scenario_signals = {"by_case": {}, "by_code": []}
        gold_quality = {
            "status": "failed",
            "failed_checks": [{"code": "positive_gold_count_low"}, {"code": "negative_trap_category_gaps"}],
        }

        merged = runner.merge_report_signals(scenario_signals, [runner.gold_quality_signal(gold_quality)])
        by_code = {row["code"]: row for row in merged["by_code"]}

        self.assertIn("gold_quality_gap", by_code)
        self.assertEqual(by_code["gold_quality_gap"]["severity"], "error")
        self.assertEqual(by_code["gold_quality_gap"]["cases"], ["positive_gold_count_low", "negative_trap_category_gaps"])

    def test_literature_semantics_failure_merges_into_report_signals(self):
        scenario_signals = {"by_case": {}, "by_code": []}
        literature_semantics = {
            "status": "failed",
            "metrics": {"failed_case_count": 2},
            "failed_cases": [{"case_id": "negation_case"}, {"case_id": "causal_case"}],
        }

        merged = runner.merge_report_signals(
            scenario_signals,
            [runner.literature_semantics_signal(literature_semantics)],
        )
        by_code = {row["code"]: row for row in merged["by_code"]}

        self.assertIn("literature_semantics_regression", by_code)
        self.assertEqual(by_code["literature_semantics_regression"]["severity"], "error")
        self.assertEqual(by_code["literature_semantics_regression"]["cases"], ["negation_case", "causal_case"])

    def test_benchmark_status_fails_on_literature_semantics_regression(self):
        heldout = {"temporal_holdout": {"available": True}, "source_heldout": {"available": True}}

        status = runner.benchmark_status([], heldout, None, {"status": "failed"})

        self.assertFalse(status["passed"])
        self.assertFalse(status["literature_semantics_passed"])

    def test_benchmark_status_fails_on_matrix_gap(self):
        heldout = {"temporal_holdout": {"available": True}, "source_heldout": {"available": True}}

        status = runner.benchmark_status([], heldout, None, benchmark_matrix={"status": "failed"})

        self.assertFalse(status["passed"])
        self.assertFalse(status["benchmark_matrix_passed"])

    def test_summarize_analysis_records_traceability_contract_fields(self):
        analyzed = {
            "analysis_pack": {
                "input_summary": {"input_count": 1, "matched_count": 1, "ambiguous_count": 0, "unmatched_count": 0},
                "database_accuracy": {
                    "identity_decision_summary": {"accepted_exact": 1},
                    "mechanism_ready_facts": [{"fact_id": "f1"}],
                    "evidence_assertion_summary": {
                        "assertion_count": 3,
                        "by_support_status": {"support": 2, "background": 1},
                    },
                },
                "conclusion_evaluation": {
                    "metrics": {"negative_trap_hits": 0},
                    "release_gate": {"passed": True},
                    "top_ranked_candidates": [{"candidate_id": "pathway_ranking:p1"}],
                    "per_gold": [
                        {
                            "gold_id": "gold_a",
                            "gold_match_status": "exact_match",
                            "top_matches": [{"candidate_id": "pathway_ranking:p1"}],
                        }
                    ],
                },
                "determinism": {"analysis_pack_hash": "abc"},
            },
            "determinism": {"response_hash": "def"},
        }

        summary = runner.summarize_analysis({"case_id": "case", "expected_gold_ids": ["gold_a"]}, analyzed)
        trace = summary["traceability"]

        self.assertEqual(trace["contract_version"], runner.TRACEABILITY_CONTRACT_VERSION)
        self.assertEqual(trace["identity_decision_count"], 1)
        self.assertEqual(trace["mechanism_ready_fact_count"], 1)
        self.assertEqual(trace["evidence_assertion_count"], 3)
        self.assertEqual(trace["evidence_support_statuses"], ["background", "support"])
        metrics = summary["conclusion_evaluation"]["metrics"]
        self.assertTrue(metrics["expected_topk_eligible"])
        self.assertEqual(metrics["expected_precision_at_1"], 1.0)
        self.assertEqual(metrics["expected_recall_at_1"], 1.0)
        expected_topk = summary["conclusion_evaluation"]["expected_topk"]
        self.assertEqual(expected_topk["by_k"]["1"]["matched_gold_ids"], ["gold_a"])

    def test_expected_topk_prefers_gold_conclusion_ranking_when_available(self):
        payload = runner.scenario_expected_topk_metrics(
            {
                "top_ranked_gold_conclusions": [
                    {"gold_id": "gold_a"},
                    {"gold_id": "gold_b"},
                ],
                "top_ranked_candidates": [{"candidate_id": "raw_candidate:1"}],
                "per_gold": [
                    {
                        "gold_id": "gold_b",
                        "gold_match_status": "exact_match",
                        "top_matches": [{"candidate_id": "raw_candidate:1"}],
                    }
                ],
            },
            ["gold_b"],
        )

        self.assertEqual(payload["metrics"]["expected_recall_at_1"], 0.0)
        self.assertEqual(payload["metrics"]["expected_recall_at_3"], 1.0)
        self.assertEqual(payload["detail"]["by_k"]["1"]["ranking_basis"], "top_ranked_gold_conclusions")
        self.assertEqual(payload["detail"]["by_k"]["3"]["matched_gold_ids"], ["gold_b"])

    def test_expected_topk_metrics_skip_abstention_scenarios_without_expected_gold(self):
        payload = runner.scenario_expected_topk_metrics(
            {
                "top_ranked_candidates": [],
                "per_gold": [
                    {
                        "gold_id": "gold_a",
                        "gold_match_status": "unsupported",
                        "top_matches": [],
                    }
                ],
            },
            [],
        )
        metrics = payload["metrics"]

        self.assertFalse(metrics["expected_topk_eligible"])
        self.assertIsNone(metrics["expected_precision_at_1"])
        self.assertIsNone(metrics["expected_recall_at_1"])
        self.assertFalse(payload["detail"]["eligible"])

    def test_expected_topk_recall_gap_enters_optimization_signals(self):
        signals = runner.case_optimization_signals(
            {
                "case_id": "low_topk",
                "status": "passed",
                "input_summary": {"input_count": 1, "ambiguous_count": 0, "unmatched_count": 0},
                "scenario_expected_metrics": {"expected_count": 2, "expected_matched": 2},
                "conclusion_evaluation": {
                    "metrics": {
                        "expected_topk_eligible": True,
                        "expected_recall_at_5": 0.5,
                        "negative_trap_hits": 0,
                    },
                    "per_gold_status_counts": {},
                    "bug_queue": [],
                },
            }
        )

        self.assertIn("expected_topk_recall_gap", {row["code"] for row in signals})

    def test_traceability_contract_blocks_unstructured_scenario_outputs(self):
        scenarios = [
            {
                "case_id": "structured",
                "status": "passed",
                "input_summary": {"input_count": 1},
                "traceability": {
                    "input_summary_present": True,
                    "identity_decision_count": 1,
                    "mechanism_ready_fact_count": 1,
                    "evidence_assertion_count": 2,
                    "conclusion_evaluation_present": True,
                    "evaluation_metrics_present": True,
                    "release_gate_present": True,
                },
            },
            {"case_id": "missing_fixture", "status": "missing"},
            {
                "case_id": "unstructured",
                "status": "passed",
                "input_summary": {"input_count": 1},
                "traceability": {
                    "input_summary_present": True,
                    "identity_decision_count": 0,
                    "mechanism_ready_fact_count": 0,
                    "evidence_assertion_count": 0,
                    "conclusion_evaluation_present": True,
                    "evaluation_metrics_present": False,
                    "release_gate_present": False,
                },
            },
        ]

        audit = runner.audit_scenario_traceability(scenarios)
        signal = runner.traceability_contract_signal(audit)

        self.assertEqual(audit["status"], "failed")
        self.assertEqual(audit["checked_count"], 2)
        self.assertEqual(audit["skipped_missing_count"], 1)
        self.assertIn("identity_decisions_missing", audit["failure_codes"])
        self.assertIn("evidence_assertions_missing", audit["failure_codes"])
        self.assertEqual(signal["code"], "scenario_traceability_contract_failed")

    def test_benchmark_status_fails_on_traceability_contract_gap(self):
        heldout = {"temporal_holdout": {"available": True}, "source_heldout": {"available": True}}

        status = runner.benchmark_status([], heldout, None, traceability_contract={"status": "failed"})

        self.assertFalse(status["passed"])
        self.assertFalse(status["traceability_contract_passed"])

    def test_benchmark_matrix_allows_registered_but_missing_fixture_in_exploratory(self):
        scenarios = [
            {
                "case_id": "direct_hmdb",
                "status": "passed",
                "benchmark_tags": {
                    "cancer_ids": ["pan_cancer"],
                    "theme_ids": ["glycolysis_lactate"],
                    "input_modes": ["stable_identifier"],
                },
            },
            {
                "case_id": "cscc_trait_score_csv",
                "status": "missing",
                "benchmark_tags": {
                    "cancer_ids": ["cscc"],
                    "theme_ids": ["amino_acid"],
                    "input_modes": ["trait_score"],
                },
            },
        ]

        audit = runner.audit_benchmark_matrix(
            scenarios,
            ["pan_cancer", "cscc"],
            ["glycolysis_lactate", "amino_acid"],
            ["stable_identifier", "trait_score"],
            False,
        )

        self.assertEqual(audit["status"], "passed")
        self.assertEqual(audit["registered_gap_count"], 0)
        self.assertEqual(audit["executed_gap_count"], 3)
        self.assertIn("cancers:cscc", audit["executed_gaps"])

    def test_benchmark_matrix_blocks_missing_executed_coverage_in_strict_mode(self):
        scenarios = [
            {
                "case_id": "direct_hmdb",
                "status": "passed",
                "benchmark_tags": {
                    "cancer_ids": ["pan_cancer"],
                    "theme_ids": ["glycolysis_lactate"],
                    "input_modes": ["stable_identifier"],
                },
            },
            {
                "case_id": "cscc_trait_score_csv",
                "status": "missing",
                "benchmark_tags": {
                    "cancer_ids": ["cscc"],
                    "theme_ids": ["amino_acid"],
                    "input_modes": ["trait_score"],
                },
            },
        ]

        audit = runner.audit_benchmark_matrix(
            scenarios,
            ["pan_cancer", "cscc"],
            ["glycolysis_lactate", "amino_acid"],
            ["stable_identifier", "trait_score"],
            True,
        )

        self.assertEqual(audit["status"], "failed")
        self.assertEqual(audit["registered_gap_count"], 0)
        self.assertIn("input_modes:trait_score", audit["executed_gaps"])

    def test_default_scenarios_cover_required_benchmark_matrix_when_registered_fixtures_count(self):
        scenarios = runner.default_scenarios()
        scenarios.append(
            {
                "case_id": "cscc_trait_score_csv",
                "status": "missing",
                "benchmark_tags": runner.DEFAULT_CSSC_BENCHMARK_TAGS,
            }
        )

        audit = runner.audit_benchmark_matrix(
            scenarios,
            runner.DEFAULT_REQUIRED_BENCHMARK_CANCERS,
            runner.DEFAULT_REQUIRED_BENCHMARK_THEMES,
            runner.DEFAULT_REQUIRED_BENCHMARK_INPUT_MODES,
            False,
        )

        self.assertEqual(audit["status"], "passed")
        self.assertEqual(audit["registered_gap_count"], 0)

    def test_heldout_audit_requires_rerun_when_requested(self):
        heldout = {
            "temporal_holdout": {
                "available": True,
                "created_at_utc": runner.utc_now(),
                "validation_metrics": {"heldout_rows": 500, "heldout_positive_rows": 10},
            },
            "source_heldout": {
                "available": True,
                "created_at_utc": runner.utc_now(),
                "source_count": 2,
                "summary": [{"source": "A", "heldout_rows": 500}, {"source": "B", "heldout_rows": 500}],
            },
        }

        audit = runner.audit_heldout_validation(heldout, None, 168, True, 100, 2, 100)

        self.assertEqual(audit["status"], "failed")
        self.assertIn("heldout_rerun_required", {row["code"] for row in audit["failed_checks"]})

    def test_heldout_audit_requires_both_temporal_and_source_reruns(self):
        heldout = {
            "temporal_holdout": {
                "available": True,
                "created_at_utc": runner.utc_now(),
                "validation_metrics": {"heldout_rows": 500, "heldout_positive_rows": 10},
            },
            "source_heldout": {
                "available": True,
                "created_at_utc": runner.utc_now(),
                "source_count": 2,
                "summary": [{"source": "A", "heldout_rows": 500}, {"source": "B", "heldout_rows": 500}],
            },
        }
        rerun = {
            "temporal": {
                "command": ["python", "run_temporal_holdout_validation.py"],
                "started_at_utc": "2026-05-23T00:00:00Z",
                "completed_at_utc": "2026-05-23T00:00:03Z",
                "duration_seconds": 3.0,
                "returncode": 0,
                "passed": True,
            },
            "source": "passed",
        }

        audit = runner.audit_heldout_validation(heldout, rerun, 168, True, 100, 2, 100)
        codes = {row["code"] for row in audit["failed_checks"]}

        self.assertEqual(audit["status"], "failed")
        self.assertFalse(audit["rerun_performed"])
        self.assertIn("heldout_rerun_required", codes)
        self.assertIn("heldout_rerun_invalid:source", codes)

    def test_heldout_audit_records_rerun_provenance_and_freshness(self):
        heldout = {
            "temporal_holdout": {
                "available": True,
                "created_at_utc": runner.utc_now(),
                "validation_metrics": {"heldout_rows": 500, "heldout_positive_rows": 10},
            },
            "source_heldout": {
                "available": True,
                "created_at_utc": runner.utc_now(),
                "source_count": 2,
                "summary": [{"source": "A", "heldout_rows": 500}, {"source": "B", "heldout_rows": 500}],
            },
        }
        rerun = {
            "temporal": {
                "command": ["python", "run_temporal_holdout_validation.py"],
                "started_at_utc": "2026-05-23T00:00:00Z",
                "completed_at_utc": "2026-05-23T00:00:03Z",
                "duration_seconds": 3.0,
                "returncode": 0,
                "passed": True,
                "report_path": "learning_runs/run/reports/temporal_holdout/temporal_holdout_validation_report.json",
            },
            "source": {
                "command": ["python", "run_source_heldout_validation.py"],
                "started_at_utc": "2026-05-23T00:00:03Z",
                "completed_at_utc": "2026-05-23T00:00:06Z",
                "duration_seconds": 3.0,
                "returncode": 0,
                "passed": True,
                "report_path": "learning_runs/run/reports/source_heldout/source_heldout_validation_report.json",
            },
        }

        audit = runner.audit_heldout_validation(heldout, rerun, 168, True, 100, 2, 100)

        self.assertEqual(audit["status"], "passed")
        self.assertTrue(audit["rerun_performed"])
        self.assertEqual(audit["rerun_details"]["temporal"]["returncode"], 0)
        self.assertEqual(audit["rerun_details"]["source"]["duration_seconds"], 3.0)
        self.assertEqual(audit["report_freshness"]["temporal"]["status"], "fresh")
        self.assertEqual(audit["report_freshness"]["source"]["status"], "fresh")

    def test_heldout_audit_detects_quality_collapse(self):
        heldout = {
            "temporal_holdout": {
                "available": True,
                "created_at_utc": runner.utc_now(),
                "validation_metrics": {
                    "heldout_rows": 500,
                    "heldout_positive_rows": 10,
                    "heldout_nonpositive_rows": 490,
                    "heldout_lift_at_100": 1.2,
                    "heldout_roc_auc": 0.6,
                    "heldout_unique_transfer_scores": 10,
                },
            },
            "source_heldout": {
                "available": True,
                "created_at_utc": runner.utc_now(),
                "source_count": 2,
                "summary": [
                    {
                        "source": "A",
                        "heldout_rows": 500,
                        "validation_metrics": {
                            "heldout_rows": 500,
                            "heldout_positive_rows": 50,
                            "heldout_nonpositive_rows": 450,
                            "heldout_lift_at_100": 0.9,
                            "heldout_roc_auc": 0.5,
                            "heldout_unique_transfer_scores": 1,
                            "heldout_top_score_tie_fraction": 1.0,
                        },
                        "source_family_diagnostics": {
                            "missing_heldout_source_kinds_in_train": ["literature_relation"],
                            "same_family_train_rows": 0,
                        },
                        "feature_transfer_diagnostics": {
                            "shared_nonconstant_feature_count": 0,
                            "heldout_only_nonconstant_feature_count": 4,
                        },
                    },
                    {"source": "B", "heldout_rows": 500, "validation_metrics": {"heldout_rows": 500, "heldout_positive_rows": 500}},
                ],
            },
        }

        audit = runner.audit_heldout_validation(heldout, None, 168, False, 100, 2, 100)
        codes = {row["code"] for row in audit["failed_checks"]}

        self.assertEqual(audit["status"], "failed")
        self.assertIn("source_heldout_no_topk_lift", codes)
        self.assertIn("source_heldout_random_auc", codes)
        self.assertIn("source_heldout_score_collapse", codes)
        self.assertIn("source_heldout_score_tie_collapse", codes)
        lift_finding = next(row for row in audit["failed_checks"] if row["code"] == "source_heldout_no_topk_lift")
        self.assertEqual(
            lift_finding["metrics"]["diagnostics"]["source_family"]["missing_heldout_source_kinds_in_train"],
            ["literature_relation"],
        )

    def test_heldout_audit_warns_on_low_auc_when_topk_lift_is_useful(self):
        heldout = {
            "temporal_holdout": {
                "available": True,
                "created_at_utc": runner.utc_now(),
                "validation_metrics": {
                    "heldout_rows": 500,
                    "heldout_positive_rows": 10,
                    "heldout_nonpositive_rows": 490,
                    "heldout_lift_at_100": 1.2,
                    "heldout_roc_auc": 0.6,
                    "heldout_unique_transfer_scores": 10,
                },
            },
            "source_heldout": {
                "available": True,
                "created_at_utc": runner.utc_now(),
                "source_count": 2,
                "summary": [
                    {
                        "source": "literature_partition",
                        "heldout_rows": 500,
                        "validation_metrics": {
                            "heldout_rows": 500,
                            "heldout_positive_rows": 50,
                            "heldout_nonpositive_rows": 450,
                            "heldout_lift_at_100": 4.0,
                            "heldout_roc_auc": 0.45,
                            "heldout_unique_transfer_scores": 100,
                            "heldout_top_score_tie_fraction": 0.01,
                        },
                    },
                    {"source": "B", "heldout_rows": 500, "validation_metrics": {"heldout_rows": 500, "heldout_positive_rows": 500}},
                ],
            },
        }

        audit = runner.audit_heldout_validation(heldout, None, 168, False, 100, 2, 100)
        failed_codes = {row["code"] for row in audit["failed_checks"]}
        warning_codes = {row["code"] for row in audit["warnings"]}

        self.assertEqual(audit["status"], "passed")
        self.assertNotIn("source_heldout_random_auc", failed_codes)
        self.assertIn("source_heldout_random_auc", warning_codes)

    def test_heldout_audit_detects_stale_or_undersized_reports(self):
        heldout = {
            "temporal_holdout": {
                "available": True,
                "created_at_utc": "2000-01-01T00:00:00Z",
                "validation_metrics": {"heldout_rows": 2, "heldout_positive_rows": 0},
            },
            "source_heldout": {
                "available": True,
                "created_at_utc": "2000-01-01T00:00:00Z",
                "source_count": 1,
                "summary": [{"source": "A", "heldout_rows": 2}],
            },
        }

        audit = runner.audit_heldout_validation(heldout, None, 168, False, 100, 2, 100)
        codes = {row["code"] for row in audit["failed_checks"]}

        self.assertEqual(audit["status"], "failed")
        self.assertIn("temporal_heldout_stale", codes)
        self.assertIn("source_heldout_stale", codes)
        self.assertIn("temporal_heldout_too_small", codes)
        self.assertIn("temporal_heldout_no_positives", codes)
        self.assertIn("source_heldout_source_count_low", codes)

    def test_heldout_validation_failure_merges_into_report_signals(self):
        scenario_signals = {"by_case": {}, "by_code": []}
        heldout_validation = {
            "status": "failed",
            "failed_checks": [{"code": "temporal_heldout_stale"}, {"code": "source_heldout_missing"}],
        }

        merged = runner.merge_report_signals(scenario_signals, [runner.heldout_validation_signal(heldout_validation)])
        by_code = {row["code"]: row for row in merged["by_code"]}

        self.assertIn("heldout_validation_gap", by_code)
        self.assertEqual(by_code["heldout_validation_gap"]["severity"], "error")
        self.assertEqual(by_code["heldout_validation_gap"]["cases"], ["temporal_heldout_stale", "source_heldout_missing"])

    def test_fixture_audit_warns_or_fails_based_on_strict_mode(self):
        scenarios = [
            {
                "case_id": "cscc_trait_score_csv",
                "status": "missing",
                "path": "missing.csv",
                "message": "missing",
                "expected_schema": {"expected_input_mode": "trait_score"},
                "candidate_paths": [{"relative_path": "backup/cscc.csv"}],
                "recovery_action": "restore it",
                "recovery_rerun_plan": {"case_id": "cscc_trait_score_csv", "validation_sequence": [{"step_id": "fixture_schema_audit"}]},
            }
        ]

        exploratory = runner.audit_fixture_availability(scenarios, ["cscc_trait_score_csv"], False)
        strict = runner.audit_fixture_availability(scenarios, ["cscc_trait_score_csv"], True)

        self.assertEqual(exploratory["status"], "passed")
        self.assertEqual(strict["status"], "failed")
        self.assertEqual(strict["missing_required_count"], 1)
        self.assertEqual(strict["missing_required"][0]["expected_schema"]["expected_input_mode"], "trait_score")
        self.assertEqual(strict["missing_required"][0]["candidate_paths"][0]["relative_path"], "backup/cscc.csv")
        self.assertEqual(strict["missing_required"][0]["recovery_rerun_plan"]["case_id"], "cscc_trait_score_csv")

    def test_load_cssc_scenario_validates_trait_score_csv_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cscc_fixture_candidate.csv"
            path.write_text(
                "trait,reported_trait,cohen_d,padj\n"
                "GCST1,Glutamine trait,0.5,0.01\n"
                "GCST2,Arachidonate trait,-1.2,0.02\n",
                encoding="utf-8",
            )

            scenario = runner.load_cssc_scenario(path, 1)

            self.assertEqual(scenario["status"], "ready")
            self.assertEqual(scenario["selected_rows"], 1)
            self.assertEqual(scenario["records"][0]["trait"], "GCST2")
            self.assertEqual(scenario["fixture_schema"]["input_mode"], "trait_score")
            self.assertEqual(scenario["fixture_schema"]["effect_column"], "cohen_d")
            self.assertEqual(scenario["expected_gold_ids"], ["cscc_amino_acid_metabolism_axis"])
            self.assertEqual(scenario["benchmark_tags"]["theme_ids"], ["amino_acid"])

    def test_registered_fixture_scenario_loads_valid_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fixture.csv"
            path.write_text(
                "trait,reported_trait,cohen_d\n"
                "GCST1,Glutamine trait,0.1\n"
                "GCST2,Arachidonate trait,-2.0\n",
                encoding="utf-8",
            )
            fixture = {
                "case_id": "cscc_trait_score_csv",
                "path": "fixture.csv",
                "fixture_type": "csv_trait_score",
                "required_for_strict_release": True,
                "identity_columns_any": ["trait", "reported_trait"],
                "effect_columns_any": ["cohen_d"],
                "expected_input_mode": "trait_score",
                "top_n": 1,
                "context": {"cancer_type": "cutaneous squamous cell carcinoma"},
                "expected_gold_ids": ["cscc_amino_acid_metabolism_axis"],
            }

            scenario = runner.registered_fixture_scenario(Path(tmp), fixture)

            self.assertEqual(scenario["status"], "ready")
            self.assertEqual(scenario["records"][0]["trait"], "GCST2")
            self.assertEqual(scenario["context"]["cancer_type"], "cutaneous squamous cell carcinoma")
            self.assertEqual(scenario["expected_gold_ids"], ["cscc_amino_acid_metabolism_axis"])

    def test_registered_fixture_scenario_reports_missing_file(self):
        fixture = {
            "case_id": "cscc_trait_score_csv",
            "path": "missing.csv",
            "fixture_type": "csv_trait_score",
            "required_for_strict_release": True,
            "identity_columns_any": ["trait"],
            "effect_columns_any": ["cohen_d"],
            "expected_input_mode": "trait_score",
        }

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = root / "backup" / "cscc_fixture_candidate.csv"
            candidate.parent.mkdir(parents=True, exist_ok=True)
            candidate.write_text("trait,cohen_d\nGCST1,1.0\n", encoding="utf-8")
            scenario = runner.registered_fixture_scenario(root, fixture)

        self.assertEqual(scenario["status"], "missing")
        self.assertEqual(scenario["case_id"], "cscc_trait_score_csv")
        self.assertTrue(scenario["fixture_registry"]["required_for_strict_release"])
        self.assertEqual(scenario["expected_schema"]["expected_input_mode"], "trait_score")
        self.assertTrue(scenario["candidate_paths"])
        self.assertEqual(scenario["recovery_rerun_plan"]["case_id"], "cscc_trait_score_csv")
        self.assertEqual(len(scenario["recovery_rerun_plan"]["validation_sequence"]), 3)

    def test_registry_required_cases_dedupes_with_cli_defaults(self):
        registry = {
            "fixtures": [
                {"case_id": "cscc_trait_score_csv", "required_for_strict_release": True},
                {"case_id": "other_fixture", "required_for_strict_release": True},
            ]
        }

        combined = runner.dedupe_preserve(["cscc_trait_score_csv"] + runner.registry_required_cases(registry))

        self.assertEqual(combined, ["cscc_trait_score_csv", "other_fixture"])

    def test_load_cssc_scenario_rejects_missing_effect_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.csv"
            path.write_text("trait,reported_trait\nGCST1,Glutamine trait\n", encoding="utf-8")

            scenario = runner.load_cssc_scenario(path, 10)

            self.assertEqual(scenario["status"], "failed")
            self.assertEqual(scenario["error_type"], "invalid_fixture_schema")
            self.assertIn("missing_effect_column", {row["code"] for row in scenario["fixture_schema"]["errors"]})

    def test_load_cssc_scenario_rejects_non_numeric_effect_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.csv"
            path.write_text("trait,reported_trait,cohen_d\nGCST1,Glutamine trait,not-a-number\n", encoding="utf-8")

            scenario = runner.load_cssc_scenario(path, 10)

            self.assertEqual(scenario["status"], "failed")
            self.assertEqual(scenario["error_type"], "invalid_fixture_schema")
            self.assertIn("effect_column_not_numeric", {row["code"] for row in scenario["fixture_schema"]["errors"]})

    def test_invalid_fixture_schema_enters_optimization_signals(self):
        signals = runner.case_optimization_signals({"status": "failed", "error_type": "invalid_fixture_schema"})

        self.assertEqual(signals[0]["code"], "fixture_schema_invalid")
        self.assertEqual(signals[0]["severity"], "error")

    def test_required_fixture_missing_merges_into_report_signals(self):
        scenario_signals = {"by_case": {}, "by_code": []}
        fixture_audit = {
            "status": "failed",
            "missing_required": [{"case_id": "cscc_trait_score_csv"}],
        }

        merged = runner.merge_report_signals(scenario_signals, [runner.fixture_audit_signal(fixture_audit)])
        by_code = {row["code"]: row for row in merged["by_code"]}

        self.assertIn("required_fixture_missing", by_code)
        self.assertEqual(by_code["required_fixture_missing"]["severity"], "error")
        self.assertEqual(by_code["required_fixture_missing"]["cases"], ["cscc_trait_score_csv"])
        by_target = {row["bucket"]: row for row in merged["by_target"]}
        self.assertIn("real_world_fixture_registry", by_target)
        self.assertIn("database_schema", by_target)

    def test_report_signal_merge_preserves_matrix_dimension_targets(self):
        scenario_signals = {"by_case": {}, "by_code": [], "by_dimension": [], "by_target": []}
        matrix = {"status": "failed", "registered_gaps": [], "executed_gaps": ["cancers:cscc", "themes:arachidonate"]}

        merged = runner.merge_report_signals(scenario_signals, [runner.benchmark_matrix_signal(matrix)])

        by_target = {row["bucket"]: row for row in merged["by_target"]}
        self.assertIn("benchmark_runner", by_target)
        self.assertIn("real_world_fixture_registry", by_target)
        by_dimension = {row["bucket"]: row for row in merged["by_dimension"]}
        self.assertIn("cancers:cscc", by_dimension)
        self.assertIn("themes:arachidonate", by_dimension)

    def test_benchmark_status_fails_on_required_fixture_gap(self):
        heldout = {"temporal_holdout": {"available": True}, "source_heldout": {"available": True}}

        status = runner.benchmark_status(
            [],
            heldout,
            None,
            {"status": "passed"},
            {"status": "passed"},
            {"status": "failed"},
        )

        self.assertFalse(status["passed"])
        self.assertFalse(status["fixture_audit_passed"])

    def test_benchmark_status_fails_on_heldout_audit_regression(self):
        heldout = {"temporal_holdout": {"available": True}, "source_heldout": {"available": True}}

        status = runner.benchmark_status([], heldout, None, {"status": "passed"}, {"status": "failed"})

        self.assertFalse(status["passed"])
        self.assertFalse(status["heldout_validation_passed"])

    def test_scenario_payload_records_expected_gold_ids(self):
        payload = runner.scenario_payload(
            "case",
            [{"name": "Glucose"}],
            {"context_terms": ["tumor"]},
            ["pan_cancer_glycolysis_lactate_axis"],
            {"allow_ambiguity": True},
        )

        self.assertEqual(payload["expected_gold_ids"], ["pan_cancer_glycolysis_lactate_axis"])
        self.assertTrue(payload["expected_identity_behavior"]["allow_ambiguity"])
        self.assertEqual(payload["status"], "ready")


if __name__ == "__main__":
    unittest.main()


