# Lookup Bundle Audit

Date: 2026-04-26

This note records a quick internal audit of the bundled lookup resources in `data/lookup/` for traceability.

## Files Checked

- `lookup_table.npz`
- `annotations.json`
- `accessions.json`

## Internal Coherence Summary

- `annotations.json` entries: `124,397`
- `accessions.json` entries: `124,397`
- Shared keys between `annotations.json` and `accessions.json`: `124,397`
- Missing keys on either side: `0`

All model tables in `lookup_table.npz` were structurally coherent:

- embedding row count matched ID row count for every model
- no duplicate IDs were detected within any model table
- no model IDs were missing from `accessions.json`
- no model IDs were missing from `annotations.json`
- all embedding arrays were `float32`

## Per-Model Coverage

- `ESM`: `124,363` embeddings, dimension `1280`
- `Prost-T5`: `124,248` embeddings, dimension `1024`
- `Prot-T5`: `124,239` embeddings, dimension `1024`
- `Ankh3-Large`: `124,338` embeddings, dimension `1536`
- `ESM3c`: `124,397` embeddings, dimension `1152`

These counts show that the bundle is internally consistent, while also confirming that not every model space has full `124,397` coverage.

## Current-UniProt Spot Checks

Two simple comparisons were run against the current UniProt REST records on 2026-04-26:

### 10-entry check

- exact full-entry matches: `9/10`
- one mismatch looked consistent with normal UniProt release drift rather than extraction corruption

### 100-entry check

- sample size: `100`
- exact full-entry matches: `87/100`
- exact entry-level match rate: `87.0%`
- total local GO terms checked: `524`
- exact term matches: `472`
- exact term-level match rate: `90.08%`

Interpretation:

- the bundle behaves like a coherent November 2025 snapshot
- most sampled annotations still agree exactly with current UniProt
- the observed differences are consistent with expected UniProt updates over time

