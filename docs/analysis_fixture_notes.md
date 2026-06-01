# Analysis Fixture Notes

Scope: C2 real-analysis fixtures for `mvp_20260513T002254`.

Fixture file:

- `tests/fixtures/validation/real_analysis_fixture_cases.json`

Cases:

- `pure_metabolite_list_energy_core`: a metabolite-only list with no fold
  change or p-value columns. Expected behavior is nonempty pathway, target, and
  disease rankings, evidence refs, and a `metabolite_list_mode` warning.
- `differential_metabolomics_table_directional`: a clean differential table
  with `log2FC`, `pvalue`, `padj`, and `direction`. Expected behavior is stable
  idempotency, directional support, nonempty rankings, and traceable evidence.
- `messy_table_duplicate_unmatched_ambiguous_name`: a realistic dirty upload
  with duplicate pyruvate names, an ambiguous tryptophan name, an ambiguous
  hexose mass feature, and an unmatched fake name. Expected behavior is explicit
  review queues and quality warnings, with only matched metabolites entering
  scoring.
- `oncology_lung_adenocarcinoma_small_table`: a small cancer-context
  differential table using glucose, lactate, tryptophan, glutamine, and
  glutathione. Expected behavior is nonempty rankings and evidence traceability
  without adding oncology-specific facts outside the frozen release.

Validation additions:

- Minimum matched, ambiguous, and unmatched counts can be asserted per case.
- Direction counts and presence of fold-change or p-value columns can be
  asserted from `analysis_pack.input_summary.feature_summary`.
- Ranking sections can be required to be nonempty.
- Quality warning codes can be required as a subset.
- Analysis pack evidence refs and literature evidence packs can be required.

Current gate result after C2:

- `gate_status = pass_with_known_blocks`
- `passed = 14`
- `blocked = 1`
- `failed = 0`

Known block:

- RT/MS2 manual feature tables remain intentionally blocked.
