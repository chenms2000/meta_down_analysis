"""D4 prompt/guard regression runner for the LLM safe adapter.

The runner evaluates local/external narrator outputs against deterministic
guards and lightweight explanation-quality checks. External LLM calls are never
made unless --run-external is passed and the external backend is explicitly
enabled/configured.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from llm_safe_adapter import (  # noqa: E402
    EXTERNAL_ADAPTER_VERSION,
    ExternalLLMConfig,
    ExternalLLMNarratorBackend,
    LLMAdapterError,
    LLM_NARRATOR_PROMPT_VERSION,
    build_local_adapter_output,
    content_hash,
    guard_adapter_output,
    guard_policy_hash,
)
from metabo_service import (  # noqa: E402
    DEFAULT_COMPOUND_MATCH_ROOT,
    DEFAULT_GRAPH_ROOT,
    DEFAULT_LITERATURE_ROOT,
    DEFAULT_NORMALIZED_ROOT,
    DEFAULT_PUBCHEM_ROOT,
    MetaboService,
)


RUNNER_VERSION = "llm_safe_adapter.regression.20260514"
DEFAULT_FIXTURE_PATH = "tests/fixtures/llm_safe_adapter/llm_safe_adapter_regression_cases.json"
DEFAULT_OUTPUT_ROOT = "validation_reports"
REPORT_NAME = "llm_safe_adapter_regression_report.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_check(name: str, passed: bool, detail: dict[str, Any] | None = None) -> dict[str, Any]:
    return {"name": name, "passed": bool(passed), "detail": detail or {}}


def issue_codes(guard: dict[str, Any]) -> list[str]:
    return sorted({str(issue.get("code") or "") for issue in guard.get("issues", []) if issue.get("code")})


def output_field_set(output: dict[str, Any]) -> set[str]:
    return {
        key
        for key in ("narrative_summary", "evidence_digest", "question_to_spec", "candidate_suggestion", "insufficient_evidence")
        if key in output
    }


def evaluate_output(input_pack: dict[str, Any], output: dict[str, Any], expected: dict[str, Any]) -> dict[str, Any]:
    guard = output.get("guard") if isinstance(output.get("guard"), dict) else guard_adapter_output(input_pack, output)
    required_fields = set(expected.get("required_output_fields") or [])
    actual_fields = output_field_set(output)
    min_language_blocks = int(expected.get("min_language_blocks", 1) or 0)
    expected_guard_passed = bool(expected.get("guard_passed", True))
    codes = issue_codes(guard)
    checks = [
        build_check(
            "guard_passed_expected",
            bool(guard.get("passed")) == expected_guard_passed,
            {"actual": bool(guard.get("passed")), "expected": expected_guard_passed, "issue_codes": codes},
        ),
        build_check(
            "required_output_fields",
            required_fields.issubset(actual_fields),
            {"actual": sorted(actual_fields), "required": sorted(required_fields)},
        ),
        build_check(
            "language_block_minimum",
            int(guard.get("language_block_count", 0) or 0) >= min_language_blocks,
            {"actual": guard.get("language_block_count", 0), "expected_min": min_language_blocks},
        ),
        build_check(
            "source_bound_language",
            not {"language_without_source_refs", "source_ref_not_in_input", "source_ref_path_not_in_input"}.intersection(codes),
            {"issue_codes": codes},
        ),
    ]
    return {
        "status": "passed" if all(check["passed"] for check in checks) else "failed",
        "guard_passed": bool(guard.get("passed")),
        "guard_issue_codes": codes,
        "guard_issue_count": int(guard.get("issue_count", 0) or 0),
        "language_block_count": int(guard.get("language_block_count", 0) or 0),
        "output_fields": sorted(actual_fields),
        "output_hash": content_hash(output),
        "checks": checks,
    }


def evaluate_probe(input_pack: dict[str, Any], probe: dict[str, Any]) -> dict[str, Any]:
    output = probe.get("output", {})
    guard = guard_adapter_output(input_pack, output)
    codes = set(issue_codes(guard))
    expected_codes = set(probe.get("expected_issue_codes") or [])
    checks = [
        build_check("probe_blocked", not guard.get("passed"), {"guard_passed": bool(guard.get("passed"))}),
        build_check(
            "expected_issue_codes",
            expected_codes.issubset(codes),
            {"actual": sorted(codes), "expected": sorted(expected_codes)},
        ),
    ]
    return {
        "probe_id": probe.get("probe_id", ""),
        "description": probe.get("description", ""),
        "status": "passed" if all(check["passed"] for check in checks) else "failed",
        "guard_passed": bool(guard.get("passed")),
        "guard_issue_codes": sorted(codes),
        "guard_issue_count": int(guard.get("issue_count", 0) or 0),
        "checks": checks,
    }


def case_input_pack(service: MetaboService | None, case: dict[str, Any], defaults: dict[str, Any]) -> dict[str, Any]:
    if isinstance(case.get("input_pack"), dict):
        return case["input_pack"]
    if service is None:
        raise ValueError(f"Case {case.get('case_id', '')} requires service-backed input pack construction.")
    return service.build_llm_adapter_input_pack(
        list(case.get("input") or []),
        question=str(case.get("question", "") or ""),
        max_paths=int(case.get("max_paths", defaults.get("max_paths", 5)) or 5),
        max_hops=int(case.get("max_hops", defaults.get("max_hops", 3)) or 3),
        evidence_limit=int(case.get("evidence_limit", defaults.get("evidence_limit", 5)) or 5),
        subgraph_max_hops=int(case.get("subgraph_max_hops", defaults.get("subgraph_max_hops", 1)) or 1),
    )


def run_local_backend(input_pack: dict[str, Any], expected: dict[str, Any]) -> dict[str, Any]:
    output = build_local_adapter_output(input_pack)
    return {
        "backend": "local",
        "adapter_version": output.get("adapter_version", ""),
        **evaluate_output(input_pack, output, expected),
    }


def run_external_backend(input_pack: dict[str, Any], expected: dict[str, Any], llm_config: ExternalLLMConfig, run_external: bool) -> dict[str, Any]:
    if not run_external:
        return {
            "backend": "external_llm",
            "status": "not_run",
            "reason": "run_external_not_requested",
            "adapter_version": EXTERNAL_ADAPTER_VERSION,
        }
    try:
        backend_result = ExternalLLMNarratorBackend(llm_config).generate(input_pack)
        guard = guard_adapter_output(input_pack, backend_result.output)
        output = {**backend_result.output, "guard": guard}
        evaluated = evaluate_output(input_pack, output, expected)
        return {
            "backend": "external_llm",
            "adapter_version": EXTERNAL_ADAPTER_VERSION,
            "backend_audit": backend_result.audit,
            **evaluated,
        }
    except LLMAdapterError as exc:
        return {
            "backend": "external_llm",
            "status": "blocked_by_guard",
            "adapter_version": EXTERNAL_ADAPTER_VERSION,
            "guard_passed": False,
            "guard_issue_codes": [exc.code],
            "guard_issue_count": 1,
            "error": {"code": exc.code, "message": str(exc), "detail": exc.detail},
        }


def summarize_case_inputs(input_pack: dict[str, Any]) -> dict[str, Any]:
    analysis = input_pack.get("analysis_pack", {})
    evidence = input_pack.get("evidence", {})
    subgraph = input_pack.get("subgraph", {})
    return {
        "adapter_input_hash": content_hash(input_pack),
        "matched_count": int(analysis.get("input_summary", {}).get("matched_count", 0) or 0),
        "ambiguous_count": int(analysis.get("input_summary", {}).get("ambiguous_count", 0) or 0),
        "unmatched_count": int(analysis.get("input_summary", {}).get("unmatched_count", 0) or 0),
        "evidence_support_count": len(evidence.get("support", []) or []),
        "subgraph_node_count": len(subgraph.get("nodes", []) or []),
        "subgraph_edge_count": len(subgraph.get("edges", []) or []),
    }


def run_case(
    service: MetaboService | None,
    case: dict[str, Any],
    defaults: dict[str, Any],
    adversarial_probes: list[dict[str, Any]],
    llm_config: ExternalLLMConfig,
    run_external: bool,
) -> dict[str, Any]:
    input_pack = case_input_pack(service, case, defaults)
    expected = case.get("expected", {})
    local = run_local_backend(input_pack, expected)
    external = run_external_backend(input_pack, expected, llm_config, run_external)
    probes = [evaluate_probe(input_pack, probe) for probe in adversarial_probes]
    case_status = "passed" if local.get("status") == "passed" and all(probe.get("status") == "passed" for probe in probes) else "failed"
    if run_external and external.get("status") not in {"passed", "blocked_by_guard"}:
        case_status = "failed"
    return {
        "case_id": case.get("case_id", ""),
        "question": case.get("question", ""),
        "status": case_status,
        "input_summary": summarize_case_inputs(input_pack),
        "local": local,
        "external": external,
        "adversarial_probes": probes,
    }


def summarize_report(case_results: list[dict[str, Any]], run_external: bool) -> dict[str, Any]:
    local_passed = sum(1 for case in case_results if case.get("local", {}).get("status") == "passed")
    probe_rows = [probe for case in case_results for probe in case.get("adversarial_probes", [])]
    probe_passed = sum(1 for probe in probe_rows if probe.get("status") == "passed")
    local_issue_counts: Counter[str] = Counter()
    probe_issue_counts: Counter[str] = Counter()
    external_issue_counts: Counter[str] = Counter()
    for case in case_results:
        local_issue_counts.update(case.get("local", {}).get("guard_issue_codes", []))
        external_issue_counts.update(case.get("external", {}).get("guard_issue_codes", []))
        for probe in case.get("adversarial_probes", []):
            probe_issue_counts.update(probe.get("guard_issue_codes", []))
    failed_cases = [case.get("case_id", "") for case in case_results if case.get("status") != "passed"]
    return {
        "case_count": len(case_results),
        "local_passed": local_passed,
        "local_failed": len(case_results) - local_passed,
        "adversarial_probe_count": len(probe_rows),
        "adversarial_probe_passed": probe_passed,
        "adversarial_probe_failed": len(probe_rows) - probe_passed,
        "external_run": bool(run_external),
        "external_status_counts": dict(Counter(case.get("external", {}).get("status", "unknown") for case in case_results)),
        "local_issue_counts": dict(sorted(local_issue_counts.items())),
        "external_issue_counts": dict(sorted(external_issue_counts.items())),
        "adversarial_issue_counts": dict(sorted(probe_issue_counts.items())),
        "failed_cases": failed_cases,
        "gate_passed": not failed_cases,
        "gate_status": "passed" if not failed_cases else "failed",
    }


def build_report(
    service: MetaboService | None,
    fixture: dict[str, Any],
    fixture_path: Path,
    llm_config: ExternalLLMConfig,
    run_external: bool,
    max_cases: int | None = None,
) -> dict[str, Any]:
    defaults = fixture.get("defaults", {})
    cases = fixture.get("cases", [])
    if max_cases is not None:
        cases = cases[:max_cases]
    adversarial_probes = fixture.get("adversarial_probes", [])
    case_results = [run_case(service, case, defaults, adversarial_probes, llm_config, run_external) for case in cases]
    report = {
        "runner_version": RUNNER_VERSION,
        "fixture_path": str(fixture_path.as_posix()),
        "fixture_version": fixture.get("fixture_version", ""),
        "release_id": getattr(service, "release_id", fixture.get("release_id", "")),
        "prompt_version": LLM_NARRATOR_PROMPT_VERSION,
        "guard_policy_hash": guard_policy_hash(),
        "external_config": llm_config.sanitized(),
        "summary": summarize_report(case_results, run_external=run_external),
        "cases": case_results,
    }
    report["report_hash"] = content_hash(report)
    return report


def write_report(report: dict[str, Any], output_root: Path) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    path = output_root / REPORT_NAME
    path.write_text(json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2), encoding="utf-8")
    return path


def build_service(args: argparse.Namespace) -> MetaboService:
    workspace = Path(args.workspace).resolve()
    llm_config = ExternalLLMConfig.from_env(
        enabled=True if args.enable_external_llm else None,
        provider=args.llm_provider,
        endpoint=args.llm_endpoint,
        model=args.llm_model,
        api_key_env=args.llm_api_key_env,
        timeout_seconds=args.llm_timeout_seconds,
        max_output_tokens=args.llm_max_output_tokens,
    )
    return MetaboService(
        workspace=workspace,
        normalized_root=(workspace / args.normalized_root).resolve(),
        graph_root=(workspace / args.graph_root).resolve(),
        pubchem_root=(workspace / args.pubchem_root).resolve(),
        compound_root=(workspace / args.compound_root).resolve(),
        literature_root=(workspace / args.literature_root).resolve(),
        release_id=args.release_id or None,
        llm_config=llm_config,
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run D4 LLM safe adapter prompt/guard regression cases.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--normalized-root", default=DEFAULT_NORMALIZED_ROOT)
    parser.add_argument("--graph-root", default=DEFAULT_GRAPH_ROOT)
    parser.add_argument("--pubchem-root", default=DEFAULT_PUBCHEM_ROOT)
    parser.add_argument("--compound-root", default=DEFAULT_COMPOUND_MATCH_ROOT)
    parser.add_argument("--literature-root", default=DEFAULT_LITERATURE_ROOT)
    parser.add_argument("--release-id", default="")
    parser.add_argument("--fixture", default=DEFAULT_FIXTURE_PATH)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--write-report", action="store_true")
    parser.add_argument("--max-cases", type=int, default=None)
    parser.add_argument("--run-external", action="store_true")
    parser.add_argument("--enable-external-llm", action="store_true")
    parser.add_argument("--llm-provider", default="")
    parser.add_argument("--llm-endpoint", default="")
    parser.add_argument("--llm-model", default="")
    parser.add_argument("--llm-api-key-env", default="")
    parser.add_argument("--llm-timeout-seconds", type=float, default=None)
    parser.add_argument("--llm-max-output-tokens", type=int, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    fixture_path = Path(args.fixture)
    fixture = read_json(fixture_path)
    service = build_service(args)
    report = build_report(
        service,
        fixture,
        fixture_path,
        service.llm_config,
        run_external=bool(args.run_external),
        max_cases=args.max_cases,
    )
    if args.write_report:
        write_report(report, Path(args.output_root))
    sys.stdout.write(json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2) + "\n")
    return 0 if report.get("summary", {}).get("gate_passed") else 2


if __name__ == "__main__":
    raise SystemExit(main())
