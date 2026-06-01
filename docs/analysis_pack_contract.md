# Analysis Pack Contract

Contract version: `analysis_pack.v1`

The analysis pack is a deterministic, audit-oriented summary returned at
`/analyze/metabolites` under the top-level `analysis_pack` key. It is derived
only from the request, selected frozen release, service config, graph edges, and
sentence-level evidence artifacts.

Required top-level fields:

- `input_summary`: analysis mode, input counts, match counts, feature summary,
  duplicate matched metabolite UIDs, strict/expanded seed counts, genetic
  exposure counts, weighted seed mass, and seed weights.
- `genetic_exposures`: GCST-level exposure records for rows with GWAS summary
  statistics provenance. These records may be used for genetic overlap, locus
  mapping, colocalization, and MR-style exposure analyses. They explicitly carry
  `identity_boundary` because a GWAS summary-statistics URL is a phenotype
  evidence source, not an HMDB/PubChem/ChEBI chemical identity identifier.
- `four_track_summary`: count and interpretation summary for the four analysis
  tracks: GCST genetic exposure, reportedTrait-derived mechanism seeds, strict
  chemical identities, and review-only unresolved traits.
- `matched`: input rows that resolved to a metabolite UID, with resolver scores
  and the selected candidate. Rows may also include
  `biological_interpretation_entities`, which preserve the exact entity while
  exposing biological pools such as succinate pool, glutathione redox pair, or
  UDP-sugar pool for interpretation. European GCST rows may include
  `source_trait_evidence` with GWAS summary-statistics URL, PMID/title, and
  annotation provenance; these fields support traceability but are not chemical
  identity identifiers.
- `expanded_candidates`: low-weight analysis seeds generated from ambiguous or
  unmatched rows when there is usable candidate evidence. These rows are marked
  with `seed_track = "expanded"`, `seed_class` (`soft_identity`,
  `analog_candidate`, or `class_or_pool`), and `seed_weight_multiplier`. Ratio
  traits are review-only in the current contract and must not be decomposed into
  ordinary metabolite abundance seeds. Expanded candidates may support
  pathway/theme/propagation hypotheses but must not be narrated as exact
  chemical identities.
- `expanded_summary`: counts and policy for the dual-track resolver, including
  `strict_matched_count`, `expanded_candidate_count`, `expanded_input_count`,
  `expanded_candidate_input_count`, `unresolved_input_count`,
  `analysis_seed_input_count`, `analysis_seed_feature_count`,
  `weighted_seed_mass`, `expanded_by_class`, and `weight_policy`.
- `ambiguous`: input rows withheld from scoring because resolver confidence or
  margin was insufficient. Ambiguous rows can still contribute to metabolic
  theme coverage through biological interpretation entities.
- `unmatched`: input rows not matched to the release, plus invalid rows with
  `resolution_status = "invalid"`.
- `biological_entity_pools`: aggregate biological pool coverage over matched,
  ambiguous, and unmatched input rows. This is an interpretation layer, not a
  replacement for exact chemical identity.
- `quality_warnings`: ordered warning objects with `code`, `severity`,
  `message`, and machine-readable `details`.
- `pathway_rankings`, `target_rankings`, `disease_rankings`: compact ranking
  rows. Every row includes `claim_refs`, `evidence_refs`, and prediction
  calibration fields such as `prediction_task`, `confidence_tier`,
  `calibrated_confidence`, `calibration_status`, `graph_distance`, and
  `boundary`.
- `prediction_model`: the calibrated prediction surface with A-D prediction
  quality assessment, high/medium/exploratory/low confidence buckets,
  `metabolic_themes`, high/medium metabolic theme buckets,
  `degraded_covered_metabolic_themes`, `context_mismatch`, and
  `low_confidence_appendix`.
- `top_explanation_paths`: compact path rows with terminal node, edge UIDs,
  evidence refs, and path confidence.
- `literature_evidence_pack`: deduplicated literature support refs, support
  classes, PMIDs, and max literature probability.
- `evidence_refs`: deduplicated edge and literature refs used anywhere in the
  pack.
- `determinism`: request hash, release ID, config hash, contract hash, and
  `analysis_pack_hash`.
- `release`: release metadata copied from the service envelope.
- `blocked_reasons`: explicit reasons why scoring could not proceed, empty when
  analysis is usable.

Traceability rule:

Every ranking row must set `claim_refs.traceability_passed = true` by linking to
at least one graph edge UID, literature support UID, or source record. Ambiguous
and unmatched rows are returned as review queues and must not be silently
promoted into scoring.

Dual-track resolver rule:

`matched` remains strict: a row is matched only when it passes the configured
top score and margin thresholds. Expanded candidates are a separate low-weight
track. Strict identities use `1.00x`; soft identities use `0.50x`; analog
candidates use `0.30x`; class/pool candidates use `0.20x`; ratio components use
`0.00x` and remain review-only until a ratio-aware trait model is introduced;
unresolved rows use `0`. If one input row yields multiple expanded candidates,
that class weight is split across those candidates.

Prediction calibration rule:

The service may emit high, medium, exploratory, and low-confidence predictions.
Narration must preserve `confidence_tier` and `boundary`; it may not turn a
target, disease, drug, overlay, or long-distance graph prediction into a stronger
claim than the calibrated prediction layer provides.

Metabolic theme coverage rule:

Metabolic themes are predicted before pathway interpretation. Each theme should
carry a `coverage` object with `related_input_count`, `matched_support_count`,
`ambiguous_support_count`, `unmatched_support_count`, `direction_consistency`,
`significance_support_count`, and `downgrade_reason`. A theme with strong input
coverage but unstable entity matching must be shown as covered-but-downgraded,
not silently omitted.

Context and inheritance rule:

Predictions that semantically mismatch an explicit context, such as neuronal,
astrocytic, taste-perception, immune-cell-specific, or hematopoietic pathways in
an epithelial-cell background, must be retained with `context_mismatch = true`
and downgraded to exploratory when needed. Drug and disease outputs inherit
upstream confidence: overlay drugs remain appendix-style hypotheses, and a drug
cannot exceed the confidence tier of its supporting target.

Determinism rule:

For the same input, release, and config, both the endpoint
`determinism.response_hash` and `analysis_pack.determinism.analysis_pack_hash`
must remain identical.
