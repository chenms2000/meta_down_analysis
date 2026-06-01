# European Trait Identity Curation

This file explains how to fill
`manual_sources/european_trait_identity/european_trait_identity_review.csv`.

The goal is not to make every GCST trait a strict metabolite. The goal is to
separate four evidence scopes:

- `strict_identity`: exact chemical identity is supported by a stable external
  identifier such as HMDB, PubChem CID, ChEBI, or full InChIKey.
- `ratio_component`: the GCST trait is a ratio or composite; components can
  support mechanism/pathway analysis, but the whole GCST is not one metabolite.
- `class_or_pool`: lipid shorthand, acylcarnitine pool, sphingolipid pool, or
  isomeric/pool-like traits without one exact structure.
- `unresolved`: platform X-code or insufficient evidence.

## Only Five Fields Are Usually Needed

For a row you can confirm, fill:

```text
review_status=reviewed
identity_scope=strict_identity | ratio_component | class_or_pool | unresolved
reviewed_name=<confirmed name>
reviewed_pubchem_cid or reviewed_hmdb_id or reviewed_chebi_id or reviewed_inchikey=<one stable ID>
evidence_source=<where the ID came from>
```

Use `review_status=defer` when you cannot confirm the identity.

## Easy Rows First

Start with rows where `priority` is 10 or 20:

- `pubchem_only`: PubChem found a candidate, but it is not yet a local strict
  identity. Open/check the CID. If it is the exact compound and not a ratio,
  isomer pool, or shorthand, set `identity_scope=strict_identity`.
- `alias_variant_missing`: likely synonym, acid/salt form, sulfate/glucuronide,
  or spelling variant. Confirm with HMDB/PubChem/ChEBI/InChIKey before strict.

## Ratio Rows

For traits like:

```text
Serine to threonine ratio
```

Do not mark the whole trait as `strict_identity`.

Use:

```text
identity_scope=ratio_component
reviewed_name=L-serine
reviewed_pubchem_cid=<CID if confirmed>
```

and a second row for threonine. These rows improve pathway/mechanism analysis
but remain low-weight ratio evidence.

## Lipids, Isomers, And Shorthand

For traits like:

```text
(16 or 17)-methylstearate (a19:0 or i19:0)
1-stearoyl-2-oleoyl-GPS (18:0/18:1)
N-palmitoyl-sphingosine (d18:1/16:0)
```

Use `class_or_pool` unless the source supplement, LIPID MAPS/HMDB, or full
InChIKey confirms one exact structure.

## X-Codes

For traits like:

```text
X-21319 levels
```

Use `unresolved` until the original study supplement or vendor mapping gives a
chemical ID. Do not infer identity from nearby traits or pathway context.

## Evidence Source Format

Good examples:

```text
PMID:36635386 supplementary table; PubChem CID:6992102
HMDB HMDB0000122 accessed from source supplement
Vendor Metabolon mapping table v2023; InChIKey=...
```

Weak examples:

```text
looks similar
probably glucose
ChatGPT
```

Weak evidence should stay `defer` or `unresolved`.

## Generate The Template

```powershell
python .\scripts\export_european_identity_review_template.py --workspace .
```

The output is:

```text
manual_sources/european_trait_identity/european_trait_identity_review.csv
```

For a smaller first batch:

```powershell
python .\scripts\export_european_identity_review_template.py --workspace . --limit 100
```

## Automatic Prefill

To let the system search PubChem and prefill high-confidence rows:

```powershell
python .\scripts\autofill_european_identity_review.py --workspace .
```

By default this processes the safest first pass, `priority <= 20`.
To continue into common local-name rows:

```powershell
python .\scripts\autofill_european_identity_review.py --workspace . --priority-max 30
```

To attempt every row, including rows that will mostly be triaged rather than
strictly identified:

```powershell
python .\scripts\autofill_european_identity_review.py --workspace . --priority-max 0
```

The autofill is conservative:

- exact PubChem title/synonym + InChIKey evidence becomes
  `review_status=auto_reviewed` and `identity_scope=strict_identity`;
- ratio traits become `auto_triaged` + `ratio_component`;
- lipid/isomer/pool-like traits become `auto_triaged` + `class_or_pool`;
- X-codes stay `unresolved`.

Rows marked `auto_reviewed` are treated as reviewed identity evidence by
`build_european_trait_annotations.py`. Rows marked `auto_triaged` document why
the trait should not be promoted to strict identity.

## Import Confirmed PubChem CIDs Into The Local Match Index

Some `auto_reviewed` PubChem/InChIKey identities are confirmed, but their CIDs
are not yet present in the local compound graph/index. Export them as a local
compound identity overlay:

```powershell
python .\scripts\export_reviewed_pubchem_compound_overlay.py --workspace . --release-id mvp_20260513T002254
```

Then rebuild the compound match index:

```powershell
python .\scripts\build_compound_match_index.py --workspace . --release-id mvp_20260513T002254
```

The overlay is imported from:

```text
manual_sources/compound_identity_overlays/<release_id>/compound_identity_overlay.csv
```

This promotes exact PubChem CID/InChIKey lookups to strict compound matching
without claiming full graph/pathway connectivity. If a CID later enters the
normalized metabolite store, the builder skips the overlay duplicate.
