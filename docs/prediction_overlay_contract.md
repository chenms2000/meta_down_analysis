# Prediction Overlay Contract

Contract version: `prediction_pack.v1`

The prediction layer is a research prioritization overlay. It does not promote
drug, cell-type, or context associations into canonical graph facts, and it is
not clinical decision support.

## Inputs

`/analyze/metabolites`, `/explain`, `/chat`, and the CLI accept optional context:

```json
{
  "context": {
    "cancer_type": "<optional cancer type>",
    "tissue": "<optional tissue>",
    "cell_type": "<optional cell type>",
    "cell_state": "<optional cell state>",
    "cell_line": "<optional cell line>"
  },
  "context_mode": "soft"
}
```

`context_mode = soft` boosts matching overlay rows and downweights rows without
context. `context_mode = hard` keeps only overlay rows with explicit matching
context terms.

## Overlay Files

Place optional files under:

```text
manual_sources/prediction_overlays/<release_id>/
```

Supported formats are `.csv`, `.tsv`, `.json`, and `.jsonl`.

Drug overlay filenames:

- `drug_targets.*`
- `drug_target_edges.*`
- `drug_predictions.*`

Useful columns:

- `drug_id`, `drug_name`
- `target_uid`, `target_id`, `target_symbol`, `target_name`, `gene_symbol`
- `mechanism`
- `confidence`
- `context_terms`, `cancer_type`, `tissue`, `cell_type`, `cell_state`, `cell_line`
- `source_name`, `source_record_id`, `evidence_level`, `license_id`

Cell overlay filenames:

- `cell_type_signatures.*`
- `cell_signatures.*`
- `cell_context_signatures.*`

Useful columns:

- `cell_type_id`, `cell_type_name`
- `target_uids`, `target_symbols`, `gene_symbols`, `marker_genes`
- `pathway_uids`, `pathway_terms`
- `metabolite_terms`
- `confidence`
- `context_terms`, `cancer_type`, `tissue`, `cell_type`, `cell_state`, `cell_line`
- `source_name`, `source_record_id`, `evidence_level`, `license_id`

Multiple values can be separated by `;`, `|`, comma, or newline.

Drug target matching is exact for IDs and gene symbols. Do not rely on generic
tokens from names such as `kinase` or `pyruvate`; include a stable `target_uid`,
UniProt accession, ChEMBL target id, or gene symbol whenever possible.

## Building From Public Sources

Use `scripts/build_prediction_overlays.py` to materialize release-local overlay
files from public drug-target and cell-line context sources:

```powershell
python .\scripts\build_prediction_overlays.py `
  --workspace . `
  --release-id mvp_20260513T002254 `
  --analysis-records .\path\to\records.json `
  --context-json .\path\to\context.json `
  --context-mode soft `
  --target-limit 15
```

The builder can fetch ChEMBL mechanism rows, DrugCentral drug-target rows, DGIdb
aggregated drug-gene interactions, Open Targets 26.03 clinical target rows, and
Cellosaurus context rows. DepMap and CELLxGENE are accepted as local CSV inputs
via `--depmap-model-csv` and `--cellxgene-signature-csv`; these local inputs may
include `cell_state` or `state` columns. The manifest records source status and
SHA-256 checksums for downloaded/input files.

## CLI

```powershell
python .\scripts\metabo_service.py --workspace . --release-id mvp_20260513T002254 predict-metabolites .\records.json --context-json .\context.json
```

The output contains `drug_rankings`, `cell_type_rankings`, `overlay_status`,
`blocked_reasons`, and deterministic hashes.
