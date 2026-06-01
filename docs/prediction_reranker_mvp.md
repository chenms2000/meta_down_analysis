# Prediction Reranker MVP

This workflow turns `analysis_pack.prediction_model` into a trainable,
calibrated prediction-ranking loop. It does not mutate canonical graph tables or
change the online `analysis_pack` contract.

## Run

```powershell
$run = "prediction_reranker_mvp"

python .\scripts\build_prediction_training_set.py `
  --workspace . `
  --run-id $run

python .\scripts\train_prediction_reranker.py `
  --workspace . `
  --run-id $run

python .\scripts\evaluate_prediction_reranker.py `
  --workspace . `
  --run-id $run
```

By default, the exporter analyzes the local validation fixtures and writes:

```text
learning_runs/<run>/prediction_training/training_examples.parquet
learning_runs/<run>/prediction_training/feature_schema.json
learning_runs/<run>/prediction_training/export_manifest.json
learning_runs/<run>/prediction_training/scored_training_examples.parquet
learning_runs/<run>/models/prediction_reranker/<task>/model.joblib
learning_runs/<run>/reports/prediction_reranker_training_report.json
learning_runs/<run>/reports/prediction_reranker_eval_report.json
```

For batch outputs already containing `analysis_pack`, pass files or directories:

```powershell
python .\scripts\build_prediction_training_set.py `
  --workspace . `
  --run-id $run `
  --analysis-json .\path\to\analysis_outputs `
  --no-default-fixtures
```

## Labels

`weak_label` is a utility label, not a truth label:

```text
2 = should enter high-confidence or top-priority predictions
1 = medium-confidence or candidate mechanism
0 = low-value, noisy, distant, or under-supported prediction
```

If `human_label` is populated later, the trainer uses it over the weak label.

## Metrics

Reports include:

```text
Precision@5
nDCG@10
MRR
ECE
Brier score
overclaim rate
```

Early fixture-only runs may fall back to a constant baseline for tasks with too
little label diversity. That is expected until more analysis packs or
human-reviewed examples are exported.
