# Formal Differential Workflow

This workflow turns a full allDEGs TraitScore comparison table into three
auditable layers:

1. **Core report**: high-signal rows only. These rows may enter the primary
   interpretation report if entity resolution, evidence refs, and confidence
   gates pass.
2. **Exploratory appendix**: broader candidates for lower-confidence pathway,
   disease, target, or mechanism hints.
3. **Full audit**: the whole allDEGs table is used for row counts,
   significance/effect summaries, direction balance, and review burden. It is
   not promoted directly into strong conclusions.

The workflow is intentionally conservative: it keeps weak, ambiguous,
graph-only, or model-only signals out of primary conclusions unless they are
clearly labeled as exploratory.

## Prepare The Workflow

```powershell
python .\scripts\run_formal_differential_workflow.py `
  --workspace . `
  --release-id mvp_20260513T002254 `
  --all-degs "D:\Desktop\TraitScore_GroupDiff_allCelltypes\allDEGs\trait_score_diff_LUAD_Epi_LUAD_Tumor_vs_Adjacent.csv"
```

If a matching `significant/*_significant_q0.05.csv` file exists beside the
allDEGs folder, it is used as the core report source automatically. Otherwise
the script selects core rows from allDEGs by q-value and effect size.

Outputs are written under:

```text
validation_reports/<release_id>/formal_differential_workflow/<input_stem>/
```

Key files:

- `core_input.csv`: selected input for primary interpretation.
- `exploratory_input.csv`: broader selected input for appendix-only analysis.
- `formal_differential_workflow_report.json`: machine-readable audit report.
- `formal_differential_workflow_report.md`: human-readable workflow summary.
- `run_formal_analysis.ps1`: generated commands for graph explanation runs.

## Default Tiers

| Tier | Default rows | Purpose |
| --- | ---: | --- |
| Core report | 120 max, 60 per direction group | Primary research interpretation |
| Exploratory appendix | 500 max, 250 per direction group | Candidate discovery and lower-confidence signals |
| Full audit | all rows | Matching burden, ambiguity, direction and significance summaries |

## Execute Core Analysis

After reviewing `core_input.csv`, run the generated command in
`run_formal_analysis.ps1`, or execute the core tier directly:

```powershell
python .\scripts\run_formal_differential_workflow.py `
  --workspace . `
  --release-id mvp_20260513T002254 `
  --all-degs "D:\Desktop\TraitScore_GroupDiff_allCelltypes\allDEGs\trait_score_diff_LUAD_Epi_LUAD_Tumor_vs_Adjacent.csv" `
  --execute-core
```

Use `--execute-exploratory` only when you want the appendix analysis and accept
the longer runtime.

## Interpretation Boundary

TraitScore or two-group differential statistics are precomputed signals. When
ordinary `log2FC` is absent, fields such as `mean_diff`, `z_wilcoxon`,
`cohen_d`, or `pseudo_log2FC_shifted` are seed-weight surrogates, not direct
abundance measurements. Report wording must preserve this distinction.
