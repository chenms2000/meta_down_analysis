# Learning Pipeline MVP

This workflow builds two research layers from a frozen release:

- unsupervised discovery embeddings for clustering and similarity search
- weak-supervision priority rankings for relations, pathways, targets, and drugs

It does not mutate canonical graph or normalized release tables.

For cloud-only training, run the same commands on the training server after
preparing the frozen release inputs there. Large learning outputs should stay
in `learning_runs/` and remain outside Git; commit the code, commands, and
provenance notes instead.

Do not publish scripts or defaults that point to a private training server.
Other users should clone the repository, prepare their own local or cloud data
workspace, and run the training commands themselves.

## Run

```powershell
$release = "mvp_20260513T002254"
$run = "learn_20260518_general_mvp"

python .\scripts\build_learning_views.py `
  --workspace . `
  --release-id $release `
  --run-id $run

python .\scripts\train_unsupervised_embeddings.py `
  --workspace . `
  --run-id $run `
  --max-entities 60000 `
  --max-features 50000 `
  --components 48 `
  --neighbors 10 `
  --neighbor-query-limit 10000 `
  --clusters 60

python .\scripts\build_weak_labels.py `
  --workspace . `
  --run-id $run `
  --negatives-per-positive 1.0

python .\scripts\train_priority_ranker.py `
  --workspace . `
  --run-id $run
```

## Outputs

```text
learning_runs/<run>/
  views/entities.parquet
  views/literature_relations.parquet
  views/drug_target_overlay.parquet
  views/cell_context_overlay.parquet
  embeddings/entity_embeddings.parquet
  rankings/entity_neighbors.parquet
  rankings/entity_clusters.parquet
  weak_labels/weak_labels.parquet
  rankings/relation_priorities.parquet
  rankings/pathway_priorities.parquet
  rankings/target_priorities.parquet
  rankings/drug_priorities.parquet
  reports/priority_ranker_report.json
  models/priority_ranker.joblib
```

## General Report

The default report should be generalized and context-free:

```powershell
python .\scripts\build_general_learning_report.py `
  --workspace . `
  --run-id $run `
  --output-prefix global
```

Outputs:

```text
learning_runs/<run>/reports/global_learning_report.md
learning_runs/<run>/reports/global_top_pathways.csv
learning_runs/<run>/reports/global_top_targets.csv
learning_runs/<run>/reports/global_top_drugs.csv
learning_runs/<run>/reports/global_top_relations.csv
learning_runs/<run>/reports/global_cluster_summary.csv
learning_runs/<run>/reports/global_weak_label_summary.csv
```

This report does not apply any cancer, tissue, cell-line, or metabolite-table
filter. Use it as the primary generalized learning output.

## Optional Context Slices

Context-specific reports are optional case studies only. They are useful when
asking a specific biological question, but they should not be treated as the
generalized model output.

```powershell
python .\scripts\build_contextual_learning_report.py `
  --workspace . `
  --release-id $release `
  --run-id $run `
  --records-json .\path\to\records.json `
  --context-json .\path\to\context.json `
  --output-prefix context
```

## Interpretation

Priority scores are research-ranking signals, not truth probabilities. Strong
literature and curated-edge evidence can rank high. Context-only and abstain
rows are capped so they cannot become strong positives through model leakage.

Proxy validation uses weak labels and shuffled negatives, so high AUC means the
pipeline separated its own weak signals; it is not a biological gold-standard
accuracy estimate.

## Source-Held-Out Validation

Use source-held-out validation to test whether a temporary model can recover
relations from an evidence source that was removed from training. This is more
informative than row-level weak-label holdout when validating cross-source
generalization.

```powershell
python .\scripts\run_source_heldout_validation.py `
  --workspace . `
  --run-id $run `
  --heldout-source DrugCentral
```

Outputs:

```text
learning_runs/<run>/reports/source_heldout/
  source_label_distribution.csv
  source_heldout_validation_report.json
  source_name_DrugCentral_top_candidates.csv
  source_name_DrugCentral_all_heldout_scored.csv
  source_name_DrugCentral_temporary_model.joblib
```

The temporary held-out model is a validation artifact only. It does not replace
the frozen baseline model or rankings. Inspect `heldout_roc_auc`,
`heldout_positive_top_5pct_fraction`, and the top candidate CSV before treating
the source as recoverable from transferable graph/literature features.

For drug-source generalization, first materialize multiple drug evidence
sources into the release overlay:

```powershell
python .\scripts\build_prediction_overlays.py `
  --workspace . `
  --release-id $release `
  --target-limit 20000 `
  --depmap-model-csv .\manual_sources\cell_signatures\depmap_model.csv `
  --depmap-gene-effect-csv .\manual_sources\cell_signatures\depmap_gene_effect.csv `
  --cellxgene-signature-csv .\manual_sources\cell_signatures\cellxgene_markers.csv `
  --single-cell-marker-csv .\manual_sources\cell_signatures\single_cell_markers.csv `
  --skip-cellosaurus
```

Then rerun the learning pipeline under a new run id and validate each drug
source:

```powershell
python .\scripts\run_source_heldout_validation.py `
  --workspace . `
  --run-id $run `
  --heldout-source ChEMBL `
  --heldout-source DrugCentral `
  --heldout-source DGIdb `
  --heldout-source "Open Targets" `
  --heldout-source label_reason=curated_and_literature
```

The held-out validator uses an uncapped transfer model that excludes
`label_reason` from validation features. This avoids an easy proxy-label leak.
When a held-out source has no negative rows, AUC/AP are intentionally absent;
use global top-fraction, tie diagnostics, and the top candidate CSV instead.

Final release-gate validation:

```powershell
python .\scripts\run_release_gate_validation.py `
  --workspace . `
  --run-id $run `
  --analysis-json .\learning_runs\$run\reports\example_structured_analysis.json
```

The gate checks non-empty cell overlay provenance, readable temporal/source
held-out metrics, generalized-mode disease/drug downgrades, drug confidence caps
against upstream targets, context-mismatch downgrades, appendix pathway status,
and `matched_support_count <= related_input_count`.
