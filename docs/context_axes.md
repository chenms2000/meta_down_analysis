# Context Axes

This project separates biological context into four independent axes. They can
be used together in requests, prediction overlays, resolver queries, and future
single-cell imports.

## Axes

| Axis | Meaning | Canonical source |
| --- | --- | --- |
| `cell_type` | Canonical cell identity, such as epithelial cell, glial cell, macrophage, fibroblast | `cell_types`, mapped primarily to Cell Ontology `CL:*` |
| `cell_state` | Transient or program-like state, such as malignant, proliferating, hypoxic, EMT-like, stem-like, glycolytic | `cell_states`, local controlled vocabulary first |
| `cancer_type` | Tumor or cancer disease context, such as LUAD, GBM, BRCA | existing `diseases` graph nodes from OncoTree, NCIt, MONDO, DOID, MeSH, Open Targets |
| `tissue` | Tissue or organ context, such as lung, brain, colon | `tissues`, mapped primarily to UBERON with EFO xrefs where available |

`cancer_type` and `cell_type` must remain separate. For example, lung
adenocarcinoma is a cancer context, while epithelial cell is a canonical cell
type. A malignant epithelial cluster would be represented as `cell_type =
epithelial cell`, `cell_state = malignant`, `cancer_type = LUAD`, and `tissue =
lung`.

## Normalized Tables

The normalized builder materializes these context tables when their raw sources
are present:

- `cell_types`, `cell_type_xrefs`, `cell_type_parent_edges`
- `cell_states`, `cell_state_xrefs`, `cell_state_parent_edges`
- `tissues`, `tissue_xrefs`, `tissue_parent_edges`

Graph projection treats these as optional tables for backward compatibility with
older releases. When present, they become `cell_type`, `cell_state`, and
`tissue` nodes with text and ID resolver entries.

## Context Input

Analysis and prediction APIs accept:

```json
{
  "context": {
    "cancer_type": "LUAD",
    "tissue": "lung",
    "cell_type": "epithelial cell",
    "cell_state": "glycolytic"
  },
  "context_mode": "soft"
}
```

Prediction overlays may use the same columns plus `context_terms`. `cell_line`
remains supported, but it is a model/sample axis rather than a replacement for
canonical `cell_type`.
