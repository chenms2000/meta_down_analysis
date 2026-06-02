# Supplementary Methods and Code Availability

## Supplementary Methods: Evidence-grounded tumor metabolism explanation workflow

We developed an open-source workflow, `meta_down_analysis`, to organize tumor-associated metabolite signals into auditable research explanations. The workflow was designed as a research interpretation layer rather than a clinical decision system. It links input metabolites or precomputed differential trait-score records to normalized metabolite entities, knowledge-graph relationships, literature evidence, ranking outputs, and confidence tiers.

Input records can include metabolite names, stable identifiers, LC-MS-like feature tables, or precomputed differential tables such as TraitScore group-comparison outputs. For precomputed differential tables, the workflow preserves the distinction between direct abundance measurements and derived statistical signals. When standard abundance fold-change fields are unavailable, fields such as `mean_diff`, `cohen_d`, `z_wilcoxon`, or pseudo-log2 fold-change values are treated as surrogate differential signals and are not interpreted as direct metabolite abundance changes.

The workflow first performs entity resolution against a frozen local release of normalized metabolite and graph resources. High-confidence matches are restricted to near-unique chemical identities supported by stable identifiers, names, or compound-index evidence. Ambiguous names, ratio traits, class-level features, and unresolved records are retained as lower-confidence or appendix-level signals rather than being promoted into primary conclusions. Ratio traits are interpreted as relative component-balance signals and are not treated as proof that individual numerator or denominator metabolites changed in absolute abundance.

Resolved and lower-confidence seed entities are propagated through a tumor metabolism knowledge graph containing metabolites, pathways, genes or targets, diseases, reactions, external cross-references, and literature evidence overlays. Candidate pathways, targets, disease-context links, and explanatory nodes are ranked using graph relationships, input support, evidence references, and context-aware downgrade rules. The output layer reports confidence tiers, traceability references, downgrade reasons, and interpretation boundaries for each conclusion.

Every conclusion is required to remain traceable to at least one structured support source, including input metabolites or trait rows, graph edges, literature evidence records, ranking fields, or explicit model outputs. The system marks weakly supported graph propagation, ambiguous identity matches, model-only rankings, and limited literature support as lower-confidence evidence. External large language models, when enabled, are used only as read-only narrative generators over frozen analysis packages; they do not create entities, alter graph scores, modify confidence tiers, or add unsupported claims.

For large differential tables, the repository provides a formal three-tier workflow. The core report uses the most statistically supported and interpretable rows for primary research interpretation. The exploratory appendix keeps broader candidate signals as lower-confidence hypotheses. The full audit layer records whole-table row counts, direction distributions, significance summaries, and review burden without promoting all rows into strong conclusions.

All outputs are intended for mechanistic research interpretation and validation prioritization. They should not be used to make clinical decisions, assign treatment recommendations, or infer patient-level actionability without independent experimental and clinical validation.

## Code Availability

The source code used for the tumor metabolism knowledge-base construction and evidence-grounded explanation workflow is available at:

```text
https://github.com/chenms2000/meta_down_analysis
```

The repository includes scripts for data-source manifest handling, normalized-store construction, compound matching, knowledge-graph projection, literature evidence integration, differential-table workflow preparation, confidence-tiered explanation generation, and regression testing. The current repository version verified for local use corresponds to commit:

```text
8640f52
```

The code is released under the MIT License. Third-party databases, literature-derived resources, external identifiers, and downloaded or rebuilt data resources remain subject to their original source licenses and citation requirements. The MIT License applies to the original code, scripts, and project documentation in the repository and does not relicense third-party data sources.

## Suggested Repository Citation

Before a Zenodo DOI is available, the repository can be cited as:

```text
chenms2000. meta_down_analysis: Evidence-grounded tumor metabolism knowledge base and research explanation workbench. GitHub; 2026. Available from: https://github.com/chenms2000/meta_down_analysis
```

After creating a GitHub Release and archiving it with Zenodo, replace the GitHub-only citation with the Zenodo DOI citation and update `CITATION.cff` and the manuscript data/code availability statement accordingly.

## Suggested Data and Code Availability Statement

```text
The source code for the tumor metabolism knowledge-base construction and evidence-grounded explanation workflow is available at https://github.com/chenms2000/meta_down_analysis. The version verified for this study corresponds to commit 8640f52. The repository is released under the MIT License. Third-party databases and literature-derived resources used or downloaded by the workflow remain subject to their original licenses and citation requirements.
```

