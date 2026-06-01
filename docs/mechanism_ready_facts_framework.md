# mechanism_ready_facts Framework

## Purpose

`mechanism_ready_facts` is the final accuracy gate before the report is allowed
to write a mechanism-style conclusion. It should not be a prettier pathway
ranking table. It should be a compact set of source-traceable factual units that
survived identity, relation, context, evidence, and direction checks.

The target flow is:

```text
input_feature
-> identity_decision
-> fact_candidates
-> readiness_gates
-> mechanism_ready_facts
-> interpretation_report core conclusions
```

## What Counts As Mechanism-Ready

A fact can enter `mechanism_ready_facts` only when it answers four questions:

1. What exact input entity is involved?
2. What biological relation is asserted?
3. What source record or evidence supports that relation?
4. What can and cannot be concluded from it?

Recommended fact types:

- `reaction_participation`: exact chemical participates in a curated reaction.
- `reaction_directional_change`: input direction can be mapped to substrate,
  product, or reversible/unknown reaction role without inversion mistakes.
- `enzyme_or_transporter_link`: reaction/module is connected to enzyme,
  transporter, or protein complex with curated source provenance.
- `module_coverage`: multiple exact chemicals hit a small curated mechanism
  module with consistent direction and source trace.
- `context_supported_assertion`: literature or omics evidence supports the
  relation in a compatible tissue/cancer/cell context.
- `trait_observation`: GCST or trait score is retained as trait-level evidence,
  not exact abundance.
- `class_observation`: lipid class or pool is retained as class-level evidence,
  not exact compound.

Facts that should not enter this table:

- pathway membership alone;
- disease/target/drug graph propagation;
- novel/background/conflict literature;
- ratio components treated as abundance;
- class/pool names promoted to exact compounds;
- input rows with ambiguous or unmatched identity decisions;
- context-mismatched facts unless explicitly labeled as appendix-only.

## Proposed Tables

### `fact_candidates`

Intermediate table. This is where broad coverage can live.

Fields:

- `candidate_fact_uid`
- `candidate_fact_type`
- `subject_uid`
- `subject_type`
- `predicate`
- `object_uid`
- `object_type`
- `relation_source`
- `source_record_uid`
- `context_uid`
- `direction`
- `role`
- `evidence_assertion_uids`
- `blocking_reasons`
- `readiness_score`

This table may contain pathway-level, module-level, reaction-level, literature,
and trait/class candidates. It is allowed to be noisy because it is not exposed
as a core mechanism conclusion.

### `mechanism_ready_facts`

Final table consumed by the service.

Fields:

- `fact_uid`
- `fact_type`
- `subject_uid`
- `subject_type`
- `predicate`
- `object_uid`
- `object_type`
- `context_uid`
- `direction`
- `role`
- `readiness_tier`: `high`, `medium`, `trait_level`, `class_level`
- `allowed_claim_scope`
- `blocking_reasons`
- `supporting_assertion_uids`
- `source_record_uids`
- `identity_decision_uids`
- `confidence_components_json`
- `boundary_text`

`allowed_claim_scope` should be machine-readable:

- `exact_reaction_fact`
- `module_level_mechanism`
- `trait_level_observation`
- `class_level_observation`
- `appendix_only`

## Readiness Gates

### Gate 1. Identity Gate

Accepted:

- `accepted_exact` for exact reaction and enzyme/transporter facts.
- `accepted_trait` only for trait-level observations.
- `accepted_class` only for class-level observations.

Rejected:

- `ambiguous`
- `unmatched`
- `rejected`
- score `0` candidate
- ratio component as chemical abundance
- class name as exact compound

### Gate 2. Relation Gate

Accepted relation sources:

- curated reaction participant with resolved chemical UID;
- curated catalyst/enzyme/transporter relation;
- curated small mechanism module;
- manually reviewed or high-precision structured evidence assertion.

Rejected or appendix-only:

- pathway membership alone;
- graph propagation path;
- disease/target/drug overlay;
- relation without source record;
- relation with only background text.

### Gate 3. Direction Gate

The system should not infer activation from membership. It can only say:

- input chemical increased/decreased;
- substrate/product role is known or unknown;
- module has direction consistency or mixed direction;
- mechanism direction is supported, unsupported, or blocked.

Ratio traits must keep numerator/denominator semantics separate:

```text
Tryptophan / pyruvate ratio up
!= tryptophan up
!= pyruvate down
```

### Gate 4. Context Gate

Context levels:

- `exact_context`: same cancer/tissue/cell type/comparison.
- `compatible_context`: same broad disease or epithelial/tumor context.
- `unspecified_context`: usable for generic biochemical facts only.
- `mismatch_context`: appendix-only.

Context mismatch should never be repaired by high pathway score.

### Gate 5. Evidence Gate

Evidence classes:

- `curated_database`: can support factual relation.
- `direct_assay`: can support mechanism wording when context compatible.
- `omics_association`: can support association wording.
- `literature_background`: appendix or explanation only.
- `novel_candidate`: overlay only.
- `conflict`: blocks strong wording.
- `model_prediction`: hypothesis only.

Only `curated_database`, `direct_assay`, and selected `omics_association` should
feed core mechanism facts. Background, novel, and conflict evidence should not.

## Current Audit

Local release checked:

```text
mvp_20260513T002254
```

Observed facts:

- `reaction_participants.parquet`: 150,077 rows.
- All observed `participant_type` values are `protein`.
- `participant_uid` is empty for all checked participants.
- `reaction_participants_v2.parquet`: 0 rows.
- `mechanism_ready_facts.parquet`: 0 rows.
- `evidence_assertions.parquet`: 51,774 rows.
- Evidence includes many `background` and `contradict` rows, so evidence
  filtering must remain strict.

Conclusion:

The current mechanism-ready layer is structurally present but effectively empty
because the reaction participant source is not resolving chemical participants.
The next accuracy improvement should focus on reaction/metabolite relation
ingestion and fact candidate gating, not on report wording.

## Implementation Status After First Repair

The first structural repair has been applied:

- Reactome `ChEBI2Reactome_PE_Reactions.txt` is parsed as a chemical reaction
  source.
- ChEBI IDs are mapped through `metabolite_xrefs.parquet`; names are not used
  as fallback for exact reaction facts.
- `reaction_participants_v2` now includes chemical source IDs, physical entity
  IDs/names, relation source, and source record UID.
- `fact_candidates.parquet` is generated before promotion.
- `mechanism_ready_facts.parquet` now includes claim scope, role, direction,
  blocking reasons, identity decision trace, confidence components, and
  boundary text.

Current release rebuild:

```text
reaction_participants_v2: 34,826 rows
fact_candidates: 121,426 rows
mechanism_ready_facts: 33,870 rows
```

Self-checks passed:

- ratio trap remains trait-only;
- class trap remains class-only;
- zero-score candidates remain blocked;
- pathway-projected module candidates are not promoted;
- mechanism-ready facts have source record and identity decision traces;
- gold metabolism terms recover reaction facts:
  glucose, lactate, pyruvate, glutamine, glutamate, succinate, citrate,
  cysteine, glutathione, carnitine.

Remaining limitation:

Reactome PE mapping gives strong chemical-to-reaction participation, but role
and direction are often still `participant` / `unknown`. These facts can support
reaction-level participation, but cannot by themselves prove reaction direction,
flux, pathway activation, or enzyme regulation.

## Implementation Status After Equation Repair

The second structural repair adds equation-level parsing and direction-aware
facts:

- Rhea `rhea-directions.tsv`, `rhea2reactome.tsv`, `rhea-reaction-smiles.tsv`,
  `rhea-chebi-smiles.tsv`, and `rhea2uniprot_sprot.tsv` are parsed from the
  local Rhea TSV archive.
- Reactome `reactome_reaction_exporter.txt` and `ReactionPMIDS.txt` are parsed
  for protein roles and publication provenance.
- New relation tables are emitted:
  - `reaction_equations`
  - `reaction_side_participants`
  - `reaction_xrefs`
  - `enzyme_reaction_links`
  - `reaction_publication_links`
- Rhea side participants are promoted only when the side component maps through
  a ChEBI-backed exact chemical identity.
- `interpretation_report` receives directional facts before role-unknown facts.

Current release rebuild:

```text
reaction_equations: 836 rows
reaction_side_participants: 1,880 rows
reaction_xrefs: 1,596 rows
enzyme_reaction_links: 207,457 rows
reaction_publication_links: 42,328 rows
fact_candidates: 330,763 rows
mechanism_ready_facts: 35,703 rows
directional_reaction_fact: 1,833 rows
role_unknown_reaction_fact: 33,870 rows
```

Calibration now checks:

- equation-level relation tables are non-empty;
- directional reaction facts are present;
- all checked mechanism facts carry source and identity traces;
- pathway-projected modules remain candidate-only;
- gold metabolism terms recover reaction facts.

Remaining limitation:

Rhea SMILES side parsing is intentionally exact-string based against Rhea ChEBI
SMILES. This avoids name fallback but may miss participants when equivalent
structures differ in serialization. A future optional pass can use structure
normalization or Rhea RDF/BioPAX to improve side coverage without weakening
identity precision.

## Optimizations

### Keep

- `identity_resolution_v2` as the first gate.
- `evidence_assertions` with support status separation.
- `mechanism_ready_facts` as the only source for core mechanism claims.
- Legacy pathway/target/disease rankings as appendix/supporting layer.

### Optimize

1. Add a `fact_candidates` table before `mechanism_ready_facts`.
2. Fix reaction participant ingestion so chemical participants are resolved from
   Reactome/Rhea/MetaCyc or equivalent sources.
3. Add catalyst/enzyme/transporter links, not just chemical-to-reaction links.
4. Replace pathway-projected modules with curated mechanism modules for the
   high-value metabolism themes.
5. Add `allowed_claim_scope` and `boundary_text` to each mechanism-ready fact.
6. Add context compatibility scoring.
7. Add direction consistency scoring across module hits.
8. Add a release gate requiring `mechanism_ready_facts` coverage for known gold
   metabolites such as glucose, lactate, pyruvate, glutamine, glutamate,
   succinate, citrate, cysteine, glutathione, carnitine, and acylcarnitines.

### Remove Or Downgrade

- Remove pathway membership-only rows from `mechanism_ready_facts`.
- Downgrade pathway-projected `pathway_modules` to candidate-only until curated
  module definitions exist.
- Keep text-mined `novel_candidate`, `background`, and `conflict` evidence out
  of mechanism-ready facts.
- Do not use disease/target/drug propagation to create mechanism facts.
- Do not use ratio components or lipid class members as exact chemicals.

## Minimal Implementation Slices

### Slice 1. Fact Candidate Layer

Add:

```text
analysis_view/fact_candidates.parquet
```

Populate from:

- exact identity decisions;
- curated metabolite-pathway edges as low-scope candidates;
- resolved reaction participants when available;
- evidence assertions with support/background/conflict separated.

No report behavior change yet.

### Slice 2. Chemical Reaction Coverage

Repair or extend ingestion so curated reactions contain chemical participants.
If the current Reactome parser only emits protein rows, add a dedicated
chemical participant extraction path or ingest Rhea/MetaCyc for biochemical
reactions.

Success condition:

```text
reaction_participants_v2 rows > 0
mechanism_ready_facts rows > 0
known gold metabolites recover reaction facts
```

### Slice 3. Readiness Gate

Promote only candidates that pass:

```text
identity accepted
source trace exists
relation source accepted
evidence not background/novel/conflict
context not mismatch
direction semantics not invalid
```

### Slice 4. Report Contract

The report should display:

- exact mechanism facts first;
- trait/class observations second;
- blocked candidates third;
- legacy rankings last.

### Slice 5. RDF/BioPAX/SBML Semantic Overlay

Add explicit reaction-side coverage without rebuilding the full release:

```text
fetch_semantic_reaction_sources.py
-> apply_semantic_reaction_overlay.py
-> reaction_side_participants
-> mechanism_ready_facts
```

Implementation rules:

- prefer Rhea RDF `substrates/products/contains/compound/accession` over
  string-level SMILES parsing when both exist;
- keep BioPAX/SBML parsers namespace-agnostic and source URI preserving;
- align participants only by stable ChEBI to local `chemical_xrefs`;
- never promote name-only, unresolved, or ambiguous semantic participants;
- use `rhea_smiles_exact_component` and `reactome_pe_participation` only as
  lower-priority fallbacks.

Current release status after overlay:

```text
reaction_side_participants: 5,674
semantic side participants: 3,794
mechanism_ready_facts: 39,479
directional_reaction_fact: 5,609
semantic mechanism facts: 3,776
```

## Calibration Traps

Add or keep fixtures for:

- ratio trap: ratio up does not mean numerator up or denominator down.
- class trap: class name does not become exact compound.
- pathway trap: membership alone does not become mechanism.
- evidence trap: background/novel/conflict does not become direct support.
- context trap: neuronal/astrocytic/taste modules do not become epithelial core.
- direction trap: substrate/product direction is not inferred when reaction
  directionality is unknown.
- coverage trap: glucose/lactate/pyruvate should recover curated biochemical
  reaction candidates once chemical reaction ingestion is fixed.
