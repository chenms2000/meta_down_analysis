# High-Benefit Improvement Flow

Release focus: `mvp_20260513T002254`

This flow applies the highest-yield precision improvements first. Literature-derived evidence remains an overlay; it must not overwrite canonical graph facts without human review.

## 1. Compound Disambiguation

Goal: reduce false metabolite identities before graph scoring.

Inputs:

- `manual_sources/compound_features/<release_id>/compound_rt.csv`
- `manual_sources/compound_features/<release_id>/compound_rt.tsv`
- `manual_sources/compound_features/<release_id>/compound_rt.jsonl`
- `manual_sources/compound_features/<release_id>/compound_ms2.csv`
- `manual_sources/compound_features/<release_id>/compound_ms2.tsv`
- `manual_sources/compound_features/<release_id>/compound_ms2.jsonl`
- `manual_sources/compound_features/<release_id>/compound_ms2.msp`
- additional `*.msp` files in the same directory

Rows may use `metabolite_uid` directly or resolve through `HMDB`, `ChEBI`, `PubChem CID`, `KEGG`, or full `InChIKey`. Ambiguous external IDs are skipped rather than forced.

Build:

```powershell
python .\scripts\build_compound_match_index.py --workspace . --release-id mvp_20260513T002254
```

Acceptance:

- `compound_rt_index.parquet` has rows when RT evidence is provided.
- `compound_ms2_index.parquet` has rows when MS2 evidence is provided.
- `compound_match_index_manifest.json.feature_import.skipped_ambiguous` is reviewed.
- `tests/fixtures/validation/rt_ms2_disambiguation_set.json` passes once local RT/MS2 features exist.

## 2. Literature Precision

Goal: keep high recall while reducing co-mention and duplicate-sentence inflation.

Changes:

- relation probability now includes pair distance, mention density, and uncertainty/background/direct-assay context weights.
- literature support probability is clustered by article: multiple sentences from the same PMID/PMCID contribute only the strongest sentence probability.

Build:

```powershell
python .\scripts\build_literature_evidence.py --workspace . --release-id mvp_20260513T002254
python .\scripts\run_evidence_precision_qa.py --workspace . --release-id mvp_20260513T002254
```

Acceptance:

- supported existing edge retention remains at or above the configured gate.
- novel/conflict candidate volume does not rise unexpectedly.
- high-risk surfaces in `docs/evidence_precision_notes.md` are reviewed before promotion.

## 3. Propagation Stability

Goal: avoid dropping too much graph mass at high-degree nodes.

The propagation collector still starts with `max_edges_per_node = 40`, but can retain up to `propagation_max_adaptive_edges_per_node = 160` when the retained edge mass is below threshold.

Validate:

```powershell
python .\scripts\run_phase15_validation.py --workspace . --release-id mvp_20260513T002254
```

Acceptance:

- `propagation_compression_lossy` warnings should decrease or show higher retained mass.
- rankings remain deterministic for repeated requests.
- explanation path traceability remains true.

## 4. Data Source Expansion

New high-yield sources are registered in `config/source_catalog.toml`:

- `massbank`
- `gnps_spectral_libraries`
- `metabolomics_workbench_refmet`
- `metabolights_public`
- `rhea`
- `uniprot_human_reviewed`
- `ncit`
- `oncotree`
- `cellosaurus`

Use open/core sources directly when enabled. Manual/API-heavy sources should first be snapshotted into raw or manual feature folders with release, license, checksum, and query provenance.
