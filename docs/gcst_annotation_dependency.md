# GCST Annotation Dependency

This repository is a research demo for interpreting upstream GCST/trait
differential tables and ordinary metabolite tables. The code can parse the
GCST table format directly, but GCST-only inputs need a local annotation table
to map each GCST accession to reported traits, metabolite candidates, ratio
components, and optional stable compound IDs.

The annotation table is a data dependency and is not distributed in this code
repository. Place the local file at one of these paths:

```text
raw_lake/European/European_trait_annotations.csv
raw_lake/European_point/European_trait_annotations.csv
```

The service checks both locations. If the annotation file is absent, GCST-only
rows can still be parsed as trait-level records, but many rows will remain
unresolved or low-confidence because the system cannot infer the corresponding
metabolite, ratio component, or compound candidate from the accession alone.

## Required and Useful Columns

| Column | Required | Purpose |
| --- | --- | --- |
| `accession_id` | yes | GCST accession used by the upstream table. |
| `reported_trait` | recommended | Human-readable trait name used for name-based fallback parsing. |
| `mapped_names` | recommended | Pipe- or semicolon-separated metabolite/component names. |
| `manual_identity_scope` | optional | Identity confidence such as `strict_identity`, `ratio_component`, `class_or_pool`, or `unresolved`. |
| `reviewed_pubchem_cid` | optional | Reviewed PubChem CID for strict or high-confidence mapping. |
| `reviewed_hmdb_id` | optional | Reviewed HMDB identifier. |
| `reviewed_chebi_id` | optional | Reviewed ChEBI identifier. |
| `reviewed_inchikey` | optional | Reviewed InChIKey for chemical identity support. |
| `local_pubchem_cid_hits` | optional | PubChem CIDs found in the local frozen compound index. |
| `summary_statistics_url` | optional | Provenance link for the GCST accession. |
| `pubmed_id` | optional | Literature provenance. |
| `paper_title` | optional | Source paper title. |
| `annotation_note` | optional | Manual review notes and downgrade reasons. |

Use [config/european_trait_annotations.template.csv](../config/european_trait_annotations.template.csv)
as a schema template only. It is not a real annotation table.

## Build or Curate the Table

For local builds, use:

```powershell
python .\scripts\build_european_trait_annotations.py --input <local_gcst_metadata.csv> --output .\raw_lake\European\European_trait_annotations.csv
```

The resulting annotation file should be treated as a research data artifact. If
full reproducibility is required for a manuscript, distribute it separately as a
supplementary table or data archive, not as part of the code repository.

## Interpretation Boundary

The annotation table does not make every GCST accession a precise metabolite.
Rows can represent direct metabolite-like traits, ratios, lipid shorthand,
class-level traits, unknown platform features, or manually reviewed strict
identities. The analysis output should preserve those categories. Ambiguous
or ratio-derived mappings should remain lower-confidence research signals
unless independently reviewed.
