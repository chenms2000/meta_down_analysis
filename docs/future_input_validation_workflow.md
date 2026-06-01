# Future Input Validation Workflow

This workflow keeps future user-uploaded files usable without turning fixtures
into an input whitelist.

## Intake

1. Run schema detection before analysis.
   - Detect `two_group_trait_comparison_table`, direct metabolite table,
     LC-MS feature table, lipid class input, and ratio trait input separately.
   - Report detected mode, row count, required columns, effect column, and
     missing/ambiguous fields to the user.
   - CLI users may run `precheck-metabolites --summary` directly on `.csv`,
     `.tsv`, `.json`, or `.jsonl` records before launching full analysis.
2. Capture context explicitly.
   - Prefer user-selected or file-derived `cancer_type`, `tissue`, `cell_type`,
     `comparison`, and `species`.
   - Show the normalized context back to the user before interpreting results.
3. Preserve raw input provenance.
   - Keep original row identifiers, effect columns, p/q fields, group labels,
     and fallback effect source such as `pseudo_log2FC_shifted`.

## Analysis Rules

1. Keep input modes separate.
   - Trait-score comparisons are not direct metabolite abundance evidence.
   - LC-MS mass-only or lipid class rows remain lower-confidence unless stable
     identifiers or orthogonal features resolve ambiguity.
2. Use conservative seed weights.
   - Strict identities can receive full seed weight.
   - Soft identity and class/pool candidates remain fractional.
   - Unresolved and ratio-observation rows must not force exact chemical claims.
3. Calibrate confidence after ranking.
   - High-confidence pathway results need multiple stable input signals or
     supporting literature, especially for trait-score tables.
   - Target, disease, drug, long-distance graph, and overlay-only predictions
     are research-priority hypotheses unless directly supported.
4. Enforce context checks.
   - Results from another explicit cancer background should be marked
     `context_mismatch` and moved to appendix/low-confidence surfaces.
   - Cell-type-specific results that conflict with the supplied cell context
     are retained for audit but downgraded.

## Fixture Policy

Fixtures are regression tests, not allowed-input lists.

- Add a fixture when a file represents a recurring real-world format,
  a high-risk ambiguity class, or a release-critical scenario.
- Do not add one fixture per user file by default.
- Keep a small matrix of representative fixtures:
  - cSCC epithelial trait-score tumor-vs-adjacent.
  - CRC epithelial trait-score tumor-vs-adjacent.
  - direct metabolite abundance with stable identifiers.
  - ambiguous lipid class/name input.
  - mass-only LC-MS input.
  - ratio trait or genetic exposure input.

## Output Review

Every report should expose:

- detected input mode and context;
- strict, expanded, unresolved, and invalid input counts;
- confidence tier, downgrade reasons, and interpretation boundary;
- evidence refs or graph/source refs for ranked claims;
- context-mismatch and appendix-only sections;
- validation priorities, not clinical or treatment recommendations.

## Release Gate Additions

Before release, verify:

- real fixture registry loads all required fixtures as `ready`;
- trait-score high-confidence pathways are not driven by a single input theme
  boost without multi-seed or literature support;
- wrong-cancer disease predictions are marked `context_mismatch`;
- LLM-safe adapter preserves confidence, boundaries, and evidence refs;
- manual adjudication remains listed as a future step when only proxy
  validation is available.
