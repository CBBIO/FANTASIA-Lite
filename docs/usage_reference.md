# FANTASIA Lite usage reference

This page contains details intentionally kept out of the main README.

## Lookup bundle

The reference is a UniProt November 2025 snapshot generated with PIS v3.1.0.
It includes EXP, IDA, IPI, IMP, IGI, IEP, TAS, and IC evidence codes.

| File | Purpose |
|---|---|
| `lookup_table.npz` | Reference embeddings |
| `annotations.json` | GO annotations and evidence |
| `accessions.json` | Reference and UniProt identifiers |

`annotations.json` and `accessions.json` each contain 124,397 entries, with
621,024 GO annotation rows. Embedding coverage is model-dependent:

| Reference space | Embeddings | Dimension |
|---|---:|---:|
| ESM-2 | 124,363 | 1,280 |
| ProstT5 | 124,248 | 1,024 |
| ProtT5 | 124,239 | 1,024 |
| Ankh3-Large | 124,338 | 1,536 |
| ESM3c | 124,397 | 1,152 |

The lookup bundle has five spaces, but the built-in Lite embedder currently
supports end-to-end ProtT5 and Ankh3-Large runs. Other spaces require compatible
externally generated query embeddings.

## Complete pipeline defaults

| Setting | Default |
|---|---|
| embedding model | `prot_t5` |
| device | automatic; CUDA preferred |
| lookup device | follows embedding device |
| distance | `cosine` |
| neighbors (`k`) | `1` |
| truncation | `0` (disabled) |
| FASTA chunk size | `500` |
| queue package | `100` |
| embedding batch | `8` |
| padded residues per batch | `12,000` |
| lookup batch | `256` |
| precision | float32 |
| raw output | off for `k=1`; automatic for `k>1` |
| TopGO | off |
| output | `outputs_<timestamp>/` |

The `minimal_pipeline.sh` helper changes embedding/model batch size to 4 and
lookup batch size to 1,024.

## Batch controls

- `--chunk-size`: FASTA records per temporary chunk.
- `--sequence-queue-package`: outer embedding package.
- `--embed-batch-size`: maximum sequences in a forward pass.
- `--model-batch-sizes MODEL=N`: per-model override.
- `--max-batch-residues`: upper bound on padded residues per batch.
- `--lookup-batch-size`: query vectors compared at once. It does not reduce the
  reference table being searched.

## Model revisions and run provenance

Lite passes the recorded revision to both `AutoTokenizer.from_pretrained` and
`AutoModel.from_pretrained`:

| Model | Repository | Revision |
|---|---|---|
| ProtT5 | `Rostlab/prot_t5_xl_uniref50` | `973be27c52ee6474de9c945952a8008aeb2a1a73` |
| Ankh3-Large | `ElnaggarLab/ankh3-large` | `2be091622e8a393f0ef21735070084123c874b6e` |

`model_provenance.yaml` is written before environment installation and updated
after installation with Python and package versions. It remains available even
when a later pipeline stage fails. `run_metadata.yaml` separately records run
status, command, Git commit, paths, parameters, coverage, and timings.

## Output columns

`results.csv` keeps the best supporting donor row for each
`(query, model, GO term, category)` combination. Raw output preserves every
selected donor/GO row before that consolidation.

| Column | Meaning |
|---|---|
| `query_accession` | Input FASTA identifier |
| `hit_rank` | Donor rank in raw output |
| `reference_id` | Internal donor identifier |
| `uniprot_accession` | Donor UniProt accession |
| `model_key` | Embedding space |
| `distance` | Query-reference embedding distance |
| `reliability_index` | Distance-derived score clipped to `[0,1]` |
| `distance_metric` | Cosine or Euclidean |
| `go_id` | Transferred GO term |
| `category` | F, P, or C namespace |
| `go_description` | GO description |
| `evidence_codes` | Donor evidence codes |

## Additional workflows

Run both built-in models sequentially:

```bash
python3 src/fantasia_pipeline.py \
  --serial-models --embed-models "prot_t5 ankh3" proteins.faa.gz
```

CPU fallback:

```bash
python3 src/fantasia_pipeline.py \
  --device cpu --no-gpu-lookup --embed-models prot_t5 proteins.faa.gz
```

Export selected ProtT5 layers without lookup:

```bash
python3 src/generate_embeddings.py \
  --fasta proteins.faa.gz --output layers.npz --models prot_t5 \
  --layer-indices 0 12 24
```

ProtT5 layer 24 and Ankh3 layer 48 are the final exported transformer layers.
Intermediate-layer export does not alter the standard last-layer lookup.

## Memory guidance

If CUDA runs out of memory:

1. process one model and enable `--serial-models`;
2. reduce `--model-batch-sizes prot_t5=2`;
3. reduce `--embed-batch-size`;
4. reduce `--max-batch-residues`;
5. only as an explicit methodological choice, set a positive `--length-filter`.

The default is uncapped. Proteins that still fail are reported without stopping
the remainder of the proteome.

## Revalidated scale benchmark

ProtT5, cosine, `k=1`, float32, uncapped runs on an RTX 3090 Ti (24 GB):

| Proteome | Proteins | Wall time | Coverage | Failed |
|---|---:|---:|---:|---:|
| *Paratomella rubra* | 20,223 | 22 min 39 s | 99.72% | 56 |
| *C. elegans* | 19,831 | 23 min 34 s | 99.75% | 49 |
| Mouse | 21,852 | 35 min 43 s | 99.31% | 150 |

These are scale estimates, not hardware-independent guarantees.
