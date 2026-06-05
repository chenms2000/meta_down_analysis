# Evidence Precision Notes

Release: `mvp_20260513T002254`

Precision status: `passed`

Filter config: `D:\meta_down_analysis\config\evidence_precision_filters.json`
Filter hash: `5899029494ccb5a7`

## Metrics Diff

| metric | baseline | current | delta | delta_pct |
| --- | --- | --- | --- | --- |
| sentence_mention_count | 899685 | 899685 | 0 | 0.00% |
| relation_candidate_count | 54453 | 54453 | 0 | 0.00% |
| edge_support_count | 33656 | 33656 | 0 | 0.00% |
| supported_existing_edge_count | 5047 | 5047 | 0 | 0.00% |
| novel_candidate_count | 22119 | 22119 | 0 | 0.00% |
| conflict_candidate_count | 3727 | 3727 | 0 | 0.00% |

## Acceptance Signals

| signal | value |
| --- | --- |
| supported_existing_edge_retention | 1.000 |
| supported_existing_edge_retention_min | 0.800 |
| novel_plus_conflict_delta | 0 |
| evidence_validation_required | phase15 validation must remain pass_with_known_blocks or passed |

## Precision Filter Hits

Manifest `precision_filters` captures exact block/downweight hit counts after rebuild.

## Configured Blocklist

| surface | entity_types | matched_fields | reason |
| --- | --- | --- | --- |
| can | gene, target | alias | English modal verb; broad gene alias creates false co-mentions. |
| via | gene, target | alias | English preposition; broad gene alias creates false co-mentions. |
| hcc | gene, target | alias | Common hepatocellular carcinoma abbreviation; valid as disease alias, too broad as gene alias. |
| mass | gene, target | alias | Generic measurement word; broad gene alias. |
| age | gene, target | alias | Generic clinical covariate; broad gene alias. |
| large | gene, target | alias | Generic adjective; broad gene alias. |
| damage | gene, target | alias | Generic biology word; broad gene alias. |
| ros | gene, target | alias | Usually reactive oxygen species in oncology/metabolism text; too broad as gene alias. |
| protein kinase | gene, target, pathway | * | Generic protein family phrase, not a specific entity. |
| kinase | gene, target, pathway | * | Generic protein family word, not a specific entity. |
| proteins | metabolite | synonym | Generic biomolecule class, not a metabolite identity. |
| mrna | metabolite | synonym | Generic transcript molecule class, not a metabolite identity. |
| lung | disease | alias | Anatomical site alone is too broad for disease evidence. |
| tumour | disease | * | Generic oncology word; does not identify a specific disease node. |

## Configured Downweight List

| surface | entity_types | matched_fields | multiplier | reason |
| --- | --- | --- | --- | --- |
| p53 | gene, target | * | 0.92 | High-frequency oncology symbol retained but slightly downweighted. |
| akt | gene, target | * | 0.92 | High-frequency kinase family shorthand retained but slightly downweighted. |
| tnf | gene, target | * | 0.94 | High-frequency inflammatory symbol retained but slightly downweighted. |
| myc | gene, target | * | 0.94 | High-frequency oncogene symbol retained but slightly downweighted. |
| breast cancer | disease | alias | 0.96 | Common disease phrase retained, but alias form is less precise than canonical disease nodes. |

## High-Risk Surface QA

| action | entity_type | surface | field | mentions | pmids | flags | sample_entities |
| --- | --- | --- | --- | --- | --- | --- | --- |
| downweight | gene | p53 | alias | 4631 | 1213 | surface_downweighted, field_downweighted, high_frequency_alias_or_synonym, short_high_frequency_surface | TP53 |
| downweight | gene | akt | alias | 4266 | 1813 | surface_downweighted, field_downweighted, high_frequency_alias_or_synonym, short_high_frequency_surface | AKT1 |
| downweight | target | myc | approved_symbol | 3750 | 1093 | surface_downweighted, short_high_frequency_surface | MYC |
| downweight | gene | cd8 | alias | 3013 | 1253 | field_downweighted, high_frequency_alias_or_synonym, short_high_frequency_surface | CD8A |
| downweight | target | tnf | approved_symbol | 2851 | 1330 | surface_downweighted, short_high_frequency_surface | TNF |
| downweight | metabolite | lactate | synonym | 2847 | 1076 | field_downweighted, high_frequency_alias_or_synonym | (S)-lactate |
| downweight | gene | tnf | symbol | 2786 | 1296 | surface_downweighted, short_high_frequency_surface | TNF |
| downweight | gene | pd l1 | alias | 2635 | 735 | field_downweighted, high_frequency_alias_or_synonym | CD274 |
| downweight | metabolite | iron | synonym | 2559 | 835 | field_downweighted, high_frequency_alias_or_synonym | iron atom |
| downweight | gene | myc | symbol | 2230 | 636 | surface_downweighted, short_high_frequency_surface | MYC |
| downweight | metabolite | lead | synonym | 2023 | 1840 | field_downweighted, high_frequency_alias_or_synonym | lead(0) |
| downweight | gene | pca | alias | 1882 | 412 | field_downweighted, high_frequency_alias_or_synonym, short_high_frequency_surface | FLVCR1 |
| downweight | gene | led | alias | 1828 | 1673 | field_downweighted, high_frequency_alias_or_synonym, short_high_frequency_surface | SMIM10L2A |
| downweight | gene | atp | alias | 1768 | 1070 | field_downweighted, high_frequency_alias_or_synonym, short_high_frequency_surface | ATP8A2 |
| downweight | metabolite | oxygen | synonym | 1683 | 1003 | field_downweighted, high_frequency_alias_or_synonym | dioxygen |
| downweight | gene | her2 | alias | 1633 | 389 | field_downweighted, high_frequency_alias_or_synonym | ERBB2 |
| downweight | gene | find | alias | 1613 | 1440 | field_downweighted, high_frequency_alias_or_synonym | DCSTAMP |
| downweight | gene | hif 1 | alias | 1590 | 535 | field_downweighted, high_frequency_alias_or_synonym | SETD2 |
| downweight | metabolite | light | synonym | 1468 | 1139 | field_downweighted, high_frequency_alias_or_synonym | photon |
| downweight | gene | light | alias | 1441 | 1117 | field_downweighted, high_frequency_alias_or_synonym | TNFSF14 |
| downweight | gene | yap | alias | 1391 | 296 | field_downweighted, high_frequency_alias_or_synonym, short_high_frequency_surface | YAP1 |
| downweight | gene | c myc | alias | 1381 | 492 | field_downweighted, high_frequency_alias_or_synonym | MYC |
| downweight | gene | spatial | alias | 1367 | 767 | field_downweighted, high_frequency_alias_or_synonym | TBATA |
| downweight | gene | out | alias | 1349 | 1207 | field_downweighted, high_frequency_alias_or_synonym, short_high_frequency_surface | TCF23 |
| downweight | metabolite | fatty acids | synonym | 1335 | 903 | field_downweighted, high_frequency_alias_or_synonym | fatty acid |
| downweight | gene | tumor necrosis factor | alias | 1328 | 1243 | field_downweighted, high_frequency_alias_or_synonym | TNF |
| downweight | gene | il 6 | alias | 1234 | 752 | field_downweighted, high_frequency_alias_or_synonym | IL6 |
| downweight | gene | ifn | alias | 1228 | 522 | field_downweighted, high_frequency_alias_or_synonym, short_high_frequency_surface | IFNA1 |
| downweight | metabolite | alpha | synonym | 1198 | 973 | field_downweighted, high_frequency_alias_or_synonym | alpha-particle |
| downweight | metabolite | lipids | synonym | 1188 | 814 | field_downweighted, high_frequency_alias_or_synonym | lipid |

## Notes

- Blocklist/downweight rules affect only sentence evidence lexicon construction.
- Novel/conflict candidates remain in the evidence layer and do not enter graph scoring.
- C3 acceptance requires comparable manifest metrics, no abnormal supported-edge collapse, and passing evidence validation.
