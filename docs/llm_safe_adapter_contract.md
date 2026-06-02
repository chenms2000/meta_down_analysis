# LLM Safe Adapter Contract

Contract versions:

- Input: `llm_safe_adapter.input.v1`
- Output: `llm_safe_adapter.output.v1`

The LLM safe adapter is an explanation boundary. AI is a narrator, not a
judge, not a resolver, and not a curator. It may rewrite structured analysis
into source-bound language, but it must not create facts or modify the graph.

## Read Boundary

The adapter may read only frozen, read-only service contracts:

- `/analyze/metabolites`
- `/analyze/differential-table`
- `/evidence`
- `/subgraph`
- `/releases`

The adapter input pack contains four payload slots: `analysis_pack`,
`evidence`, `subgraph`, and `release`. It must not receive `/resolve` output as
a decision surface, resolver internals, write APIs, or mutable graph state.

## Write Boundary

The adapter may output only these explanation objects:

- `narrative_summary`
- `evidence_digest`
- `question_to_spec`
- `candidate_suggestion`
- `insufficient_evidence`

Envelope metadata such as `contract_version`, `adapter_version`, `status`,
`guard`, and `determinism` is allowed for validation and reproducibility.

## Required Source Binding

Every natural-language field named `text`, `reason`, `suggestion`, `digest`,
`narrative`, or `claim` must include non-empty `source_refs`.

Each source ref must point back into the adapter input:

```json
{
  "source_type": "evidence",
  "ref_id": "sent_hk1_cancer",
  "path": "$.evidence.sentences[0].sentence_uid"
}
```

Allowed `source_type` values are `analysis_pack`, `evidence`, `subgraph`, and
`release`. PMID, PMCID, sentence, support, relation, edge, node, and release
refs are valid only when present in the input pack.

Entity and pathway display labels can contain words that look like biomedical
claims. The local adapter may return such labels as structured data fields, but
must not weave them into free-text claims unless the wording is explicitly
source-bound and safe.

## Prediction Model Boundary

`analysis_pack.prediction_model` is the structured prediction and calibration
surface. The service generates it before any local or external narrator is
called.

Each pathway, target, disease, drug overlay, or context overlay result carries
prediction and calibration fields such as:

```json
{
  "prediction_task": "pathway_prediction",
  "result_type": "metabolic",
  "confidence_tier": "medium",
  "calibrated_confidence": 0.48,
  "calibration_status": "calibrated_in_scope",
  "input_support_count": 1,
  "graph_distance": 1,
  "is_directly_supported": true,
  "is_extrapolated": false,
  "needs_validation": true,
  "boundary": "Medium-confidence prediction; useful as a candidate mechanism and should be reviewed with context."
}
```

Confidence tiers are `high`, `medium`, `exploratory`, and `low`. The narrator
may describe all predictions, but it must preserve the tier and boundary rather
than upgrading exploratory predictions into strong claims.

The global `prediction_model.assessment` assigns an A-D prediction quality
grade:

- `A / well_calibrated`: prediction quality is strong enough to prioritize
  high-confidence mechanisms.
- `B / usable_with_review`: predictions are usable for prioritization, with
  input and context review.
- `C / exploratory`: predictions are mainly exploratory and need additional
  evidence.
- `D / input_limited`: fix input matching before interpretation.

Ranking rows include `calibrated_confidence`, computed from raw score, input
support, match confidence, direction consistency, result type, graph distance,
and node genericity. This confidence is for prediction prioritization; it does
not mutate graph probabilities or curated facts.

## Explicit Prohibitions

The adapter must not:

- make entity resolution decisions
- choose canonical IDs
- change or create `p_final`
- emit new scores or score components
- write curated graph facts
- add entities or edges
- turn novel or conflict literature candidates into curated truth
- claim causality, proof, or definitive validation unless that exact claim is
  already represented in the structured input and source refs
- override `analysis_pack.prediction_model`, invent predictions outside it, or
  ignore `confidence_tier` and `boundary`

Forbidden output fields include `p_final`, `score`, `score_components`,
`calibrated_prob`, `resolved_entity_uid`, `entity_resolution_decision`,
`new_facts`, `new_entities`, `new_edges`, `curated_graph_write`, and
`graph_mutation`.

## Insufficient Evidence

If the input pack has no analysis rows, evidence support, or subgraph slice, the
adapter must return:

```json
{
  "status": "insufficient_evidence",
  "insufficient_evidence": {
    "reason": "...",
    "source_refs": [...]
  }
}
```

It must not fill the gap with plausible biomedical statements.

## Local Guard

`scripts/llm_safe_adapter.py` implements the deterministic guard. It checks:

- input endpoints are within the read boundary
- output top-level fields are authorized
- forbidden decision or score fields are absent
- every natural-language segment has source refs
- every cited ref and cited path exists in the input pack
- entity, edge, sentence, PMID, and PMCID tokens in text exist in the input
- numeric claims in text are present in the input
- mutating, resolver-like, curation, and overstrong claim patterns are blocked

The local adapter is deliberately small and deterministic. It is useful for
fixtures and safety tests before any external LLM is connected.

## Service Wiring

`scripts/metabo_service.py` exposes the D2 service integration through
`POST /explain` and the `explain` CLI subcommand.

The endpoint accepts the same metabolite records used by `/analyze/metabolites`
plus optional explanation controls:

```json
{
  "records": [{"HMDB": "HMDB0000122", "log2FC": 1.2}],
  "question": "Summarize the evidence.",
  "max_paths": 25,
  "max_hops": 4,
  "evidence_limit": 25,
  "subgraph_max_hops": 1,
  "adapter_backend": "local"
}
```

The service builds an `llm_safe_adapter.input.v1` pack only from these read-only
calls:

- `/analyze/metabolites` for `analysis_pack`
- `/evidence` for sentence and literature support
- `/subgraph` for a bounded graph slice
- `/releases` for release metadata

It then calls the local deterministic adapter and returns an
`explain_endpoint.v1` envelope:

```json
{
  "contract_version": "explain_endpoint.v1",
  "status": "ok",
  "adapter_backend": "local",
  "guard_passed": true,
  "explanation": {
    "contract_version": "llm_safe_adapter.output.v1",
    "status": "ok"
  }
}
```

If a supplied or future external adapter output fails guard validation, the
endpoint returns `status = "blocked_by_guard"` and withholds the unsafe
narrative. Guard issues remain visible for audit.

## External LLM Backend

D3 adds an optional external narrator backend behind the same `/explain`
contract. It is disabled by default.

External LLM modes accept plain assistant text and wrap it into the guarded
output contract with local deterministic code. `external_text_chunked` is the large-input path: it sends
compact, section-specific payloads for input parsing, ratio/class/identity
signals, rankings, evidence, and low-confidence appendices; each chunk is
guarded before a final summary is generated from the approved chunks.

Enable it only with explicit configuration:

```powershell
$env:LLM_SAFE_ADAPTER_ENABLE_EXTERNAL = "1"
$env:LLM_SAFE_ADAPTER_MODEL = "your-narrator-model"
$env:LLM_SAFE_ADAPTER_API_KEY = "<secret>"
python scripts\metabo_service.py --workspace . --enable-external-llm --llm-model your-narrator-model explain records.json --adapter-backend external_llm
```

The external backend uses an OpenAI-compatible chat-completions transport by
default. The endpoint can be configured with `LLM_SAFE_ADAPTER_ENDPOINT` or
`--llm-endpoint`. The API key is never written into audit metadata.

External backend rules:

- input is still only `llm_safe_adapter.input.v1`
- structured output is built locally as `llm_safe_adapter.output.v1`
- the external model only writes narrative text inside fixed local fields
- temperature is fixed at `0`
- JSON response mode is not required for external narration
- no tools, browsing, resolver access, parquet access, graph writes, or release
  writes are exposed
- output always passes through the local deterministic guard
- guard failure returns `blocked_by_guard` and withholds the unsafe narrative
- audit metadata records model, provider, endpoint, prompt hash, request hash,
  response hash, raw output hash, and adapter input hash

If the backend is not explicitly enabled, missing model/API key, emits unsafe
text, or raises a transport error, `/explain` returns `blocked_by_guard`
with a guard issue such as `external_llm_disabled`,
`missing_llm_model`, `missing_llm_api_key`, `overstrong_claim`, or
`external_llm_transport_error`.

## Validation Fixtures

Fixtures live in `tests/fixtures/llm_safe_adapter/` and cover:

- normal analysis pack summary
- ambiguous and unmatched review explanation
- evidence digest with PMID, PMCID, and sentence refs
- no-evidence `insufficient_evidence` fallback
- adversarial prompt that tries to resolve entities, add facts, and change
  scores

Run the dedicated validation runner with:

```powershell
python scripts\run_llm_safe_adapter_validation.py
```

This runner is separate from Phase 1.5 release validation so the existing
release gate and its known RT/MS2 block remain unchanged.

Run the D2 service wiring tests with:

```powershell
python -m unittest tests.test_metabo_service tests.test_llm_safe_adapter
```

These tests also cover the D3 external backend using a fake transport: external
LLM is opt-in, safe text can pass, and adversarial text is blocked.

## D4 Regression Pack

D4 adds prompt/guard regression over realistic analysis requests:

- `tests/fixtures/llm_safe_adapter/llm_safe_adapter_regression_cases.json`
- `scripts/run_llm_safe_adapter_regression.py`
- `validation_reports/llm_safe_adapter_regression_report.json`

The regression runner builds real `llm_safe_adapter.input.v1` packs from the
service, evaluates the local narrator, optionally evaluates the external
narrator, and runs adversarial supplied-output probes against the same input
pack.

Default run, with no external LLM call:

```powershell
python scripts\run_llm_safe_adapter_regression.py --workspace . --release-id mvp_20260513T002254 --write-report
```

External LLM run, explicit opt-in:

```powershell
$env:LLM_SAFE_ADAPTER_ENABLE_EXTERNAL = "1"
$env:LLM_SAFE_ADAPTER_MODEL = "your-narrator-model"
$env:LLM_SAFE_ADAPTER_API_KEY = "<secret>"
python scripts\run_llm_safe_adapter_regression.py --workspace . --release-id mvp_20260513T002254 --run-external --enable-external-llm --llm-model your-narrator-model --write-report
```

Reported metrics include:

- local pass/fail count
- external status counts
- adversarial probe pass/fail count
- local, external, and adversarial guard issue-code counts
- prompt version
- guard policy hash
- per-case adapter input hash and output hash

The default D4 gate passes only when every local narrator output passes quality
checks and every adversarial probe is blocked with its expected issue codes.
External LLM results are reported only when explicitly requested.
