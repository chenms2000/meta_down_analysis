# Database Accuracy Store v2 Schema

Contract version: `database_accuracy_store.v2`

This store is built by:

```text
python .\scripts\build_database_accuracy_store_v2.py --workspace . --release-id <release_id>
```

Optional semantic source overlay:

```text
python .\scripts\fetch_semantic_reaction_sources.py --workspace . --release-id <release_id>
python .\scripts\apply_semantic_reaction_overlay.py --workspace . --release-id <release_id>
```

The overlay path is the preferred way to add RDF/BioPAX/SBML reaction-side
coverage to an existing release without rebuilding the full source/entity store.

Output root:

```text
database_accuracy_store/<release_id>/
```

## Purpose

The v2 store separates source records, entities, identity semantics, relations,
evidence assertions, and analysis views. It is not a replacement for the
existing release in one step; it is an accuracy-oriented fact layer that the
service can migrate onto incrementally.

## Directory Layout

```text
database_accuracy_store/<release_id>/
  database_accuracy_store_manifest.json
  entity_store/
    source_records.parquet
    chemical_entities.parquet
    chemical_xrefs.parquet
    chemical_names.parquet
    metabolite_classes.parquet
    chemical_class_members.parquet
    trait_entities.parquet
    trait_components.parquet
  relation_store/
    reaction_entities.parquet
    reaction_participants_v2.parquet
    reaction_catalysts.parquet
    reaction_equations.parquet
    reaction_side_participants.parquet
    reaction_xrefs.parquet
    enzyme_reaction_links.parquet
    reaction_publication_links.parquet
    pathway_modules.parquet
    module_members.parquet
  evidence_store/
    evidence_contexts.parquet
    evidence_sentences.parquet
    evidence_assertions.parquet
  analysis_view/
    mechanism_ready_facts.parquet
```

## Entity Store

### `chemical_entities`

Exact or near-exact chemicals only. The table is derived from
`normalized_store/metabolites.parquet`.

Important fields:

- `chemical_uid`
- `entity_granularity`: `exact_compound`, `lipid_species`, or `unknown`
- `canonical_name`
- `formula`
- `monoisotopic_mass`
- `inchikey`
- `inchikey14`
- `identity_status`: `canonical` or `review`

Class-like or identifier-poor records are not promoted as strong exact chemical
facts.

### `chemical_names`

Names are no longer treated as identity facts by themselves.

Important fields:

- `name`
- `name_type`: `canonical`, `synonym`, or `class_name`
- `name_risk`: `low`, `medium`, or `high`

High-risk names, such as generic lipid classes, should create candidates or
class-level observations, not exact compound decisions.

### `metabolite_classes`

Class, pool, and family entities separated from exact chemicals.

Examples:

- sphingomyelin
- bile acid
- glutathione redox pair
- acylcarnitine pool

### `trait_entities`

GCST and other score/signature traits are first-class entities.

Important fields:

- `accession_id`
- `trait_type`: `metabolite_level`, `metabolite_ratio`, `class_trait`, or
  `unknown_trait`
- `identity_scope`: `exact_chemical`, `ratio`, `class_level`, or `unknown`

### `trait_components`

Ratio and composite traits keep their components without turning those
components into ordinary abundance seeds.

Important fields:

- `component_role`: `numerator`, `denominator`, or `component`
- `component_entity_type`: `chemical`, `class`, or `unknown`
- `direction_semantics`: `same_direction`, `inverse_direction`, or `undefined`
- `curation_status`: `rule_inferred` or `review_required`

## Relation Store

### `reaction_entities`

Reaction-level facts derived from `normalized_store/reactions.parquet`.

### `reaction_participants_v2`

Chemical participation in reactions with explicit role:

- `substrate`
- `product`
- other source roles preserved when role cannot be normalized

Important provenance fields:

- `relation_source`: e.g. `rhea_rdf_explicit_side`,
  `rhea_smiles_exact_component`, or `reactome_pe_participation`
- `semantic_source_uri`: RDF/BioPAX/SBML object URI when available

This is the layer that should eventually replace pathway-membership-first
mechanism interpretation.

### `reaction_equations`

Equation-level reaction records parsed primarily from Rhea and mapped to
Reactome reactions where cross-references exist. Includes directionality,
reaction SMILES when used, directional Rhea IDs, and source record provenance.

### `reaction_side_participants`

Substrate/product side assignments from equation sources. Preferred source
order:

1. Rhea RDF explicit side participants
2. Rhea/BioPAX or Reactome/SBML semantic side participants
3. Rhea directed reaction SMILES exact components
4. Reactome PE participation with unknown role

Only stable ChEBI-aligned participants enter this table. Name-only or
unresolved semantic participants remain candidate-only with
`unresolved_semantic_participant`.

### `reaction_xrefs`

Cross-references between local reaction UIDs and external reaction IDs such as
Rhea. Direction is preserved.

### `enzyme_reaction_links`

Reaction-to-protein/enzyme links from Rhea UniProt mappings and Reactome protein
role exports. These enrich mechanism facts but do not override chemical identity
or reaction-side gates.

### `reaction_publication_links`

Reaction-to-PMID provenance copied from Reactome reaction publication mappings.

### `pathway_modules`

Initial modules are projected from curated pathways. They are intentionally
named modules rather than final mechanisms; curated mechanism modules can later
replace or refine them.

### `module_members`

Connects modules to reactions and supporting chemicals.

## Evidence Store

### `evidence_contexts`

Context is a first-class object. Current builder emits `context_unspecified`
where source evidence lacks curated context. Context-specific mechanisms should
not be inferred from unspecified context alone.

### `evidence_sentences`

Sentence provenance copied from normalized article sentences.

### `evidence_assertions`

Literature support is converted into assertion rows:

- `subject_uid`
- `predicate`
- `object_uid`
- `evidence_class`
- `support_status`: `support`, `contradict`, `uncertain`, or `background`
- `sentence_uid`
- `confidence_components_json`

This prevents background or novel text-mined evidence from being treated the
same as direct curated support.

## Analysis View

### `mechanism_ready_facts`

Derived facts suitable for future mechanism scoring. This is not a final
mechanism result. It packages exact chemicals, reaction participation, available
assertion links, and source provenance into a stable view for the service layer.

Allowed claim scopes include:

- `directional_reaction_fact`
- `bidirectional_reaction_fact`
- `exact_reaction_fact`
- `role_unknown_reaction_fact`

Semantic directional facts are ranked before SMILES-derived facts in the service
output. Pathway/module membership alone is never promoted here.

## Migration Rule

The existing `graph_projection` remains usable for browsing and legacy ranking,
but high-accuracy mechanism analysis should migrate to:

```text
input_features -> identity_decisions -> reaction/module facts -> evidence_assertions -> mechanism_ready_facts
```

not:

```text
name match -> metabolite node -> pathway ranking -> conclusion
```
