# Two-Group Trait Comparison Input

This input mode is separate from direct metabolite tables.

## Detection

Rows are treated as `two_group_trait_comparison_table` only when they contain:

- a trait identifier such as `trait`, `GCST`, `accession_id`, or `gwas_trait`;
- both `group1` and `group2`;
- at least one differential statistic such as `log2FC`, `p_ttest`,
  `p_wilcoxon`, `q_ttest`, `q_wilcoxon`, `direction`, `mean_diff`, or
  `median_diff`;
- no direct metabolite identifier such as `metabolite`, `name`, `HMDB`,
  `ChEBI`, `PubChem CID`, `KEGG`, `InChIKey`, `formula`, `mz`, or `ms2_peaks`.

Direct metabolite tables continue to use `metabolite_table`.

## Mapping

| Source column | Internal field |
| --- | --- |
| `trait`, `GCST`, `accession_id` | `accession_id` and `trait` |
| `log2FC` | oriented `log2FC`, expressed as primary group vs reference group |
| `pseudo_log2FC_shifted` | fallback effect size when `log2FC` is blank; sign is assigned from group direction |
| `signed_log2FC_absmean` | secondary fallback effect size when standard and pseudo `log2FC` are both blank |
| `p_ttest`, `p_wilcoxon` | `pvalue` |
| `q_ttest`, `q_wilcoxon`, `fdr`, `padj` | `padj` |
| `direction` | `up`, `down`, `unchanged`, or `unknown` |
| `celltype_l`, `cell_type_l`, `celltype` | `cell_type` context |
| `group1`, `group2`, `group_col` | retained as comparison provenance |

The original comparison direction and original `log2FC` are retained in
`comparison_original_direction` and `comparison_original_log2FC`.

## Direction Rule

The primary group is inferred from the group labels. Tumor-like labels
(`Tumor`, `Cancer`, `Malignant`, `Case`) are primary. Adjacent/normal/control
labels are reference. If neither side is recognizable, `group1` is primary.

Examples:

| `group1` | `group2` | source `direction` | source `log2FC` | internal `direction` | internal `log2FC` |
| --- | --- | --- | ---: | --- | ---: |
| `Tumor` | `Adjacent` | `group1_hi` | 1.25 | `up` | 1.25 |
| `Adjacent` | `Tumor` | `group1_hi` | 1.25 | `down` | -1.25 |
| `Adjacent` | `Tumor` | `group2_hi` | -1.25 | `up` | 1.25 |

After conversion, the usual analysis flow uses a dual-track resolver:

- `strict_identity`: unchanged strict `matched` rows. These require the resolver
  top score and margin thresholds and support exact entity detail, precise
  evidence lookup, and strong report wording.
- expanded candidates: ambiguous or unmatched rows with usable candidate
  evidence are retained separately as low-weight seeds. `soft_identity` uses
  `0.50x`, `analog_candidate` uses `0.30x`, `class_or_pool` uses `0.20x`, and
  `ratio_component` uses `0.25x`; if one trait yields multiple candidates, that
  class weight is split across them.

European GCST traits may use `raw_lake/European/European_trait_annotations.csv`
as a pre-annotation layer. Local PubChem CIDs in the frozen index remain strict;
local-name candidates can enter the expanded track; ratio traits are interpreted
as component seeds; lipid/isomer/shorthand traits can be mapped to class or pool
evidence; unresolved `X-...`-style traits stay in the review queue.

Direct metabolite tables remain separate and are analyzed directly from their
metabolite identifiers/names/statistics. The UI and report show strict matched,
expanded candidates, unresolved inputs, and weighted analysis seeds separately.
