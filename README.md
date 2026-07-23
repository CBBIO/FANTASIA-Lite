# FANTASIA Lite V1

FANTASIA Lite is a standalone protein-function annotation workflow. It embeds
protein sequences, searches a local reference embedding table, and transfers
experimentally supported Gene Ontology (GO) annotations from the nearest
reference proteins.

Lite uses flat files and does not require PostgreSQL or RabbitMQ. For
database-backed reference management, five built-in embedders, and larger
production workflows, use [FANTASIA Full](https://github.com/CBBIO/FANTASIA).

## Start here

The shortest complete Lite run uses ProtT5, cosine distance, `k=1`, full-length
query embeddings, and GPU acceleration by default when CUDA is available.
An NVIDIA GPU is strongly recommended; CPU mode is a supported but substantially
slower fallback.

### 1. Clone the repository

```bash
git clone https://github.com/CBBIO/FANTASIA-Lite.git
cd FANTASIA-Lite
```

The current default branch is `FANTASIA-Lite-V1`. The earlier implementation
remains available in `fantasia-lite-V0`.

### 2. Install the lookup bundle

`data/lookup/` is excluded from Git and is not present in a fresh clone. Create
it before installing the reference data:

```bash
mkdir -p data/lookup
```

Download `fantasia_lite_data_folder.tgz` from the
[FANTASIA Lite Zenodo record](https://doi.org/10.5281/zenodo.19742925), place it
in the repository root, and extract it:

```bash
tar -xzf fantasia_lite_data_folder.tgz
```

The archive already contains the `data/lookup/` path. Confirm the installation:

```bash
test -s data/lookup/lookup_table.npz
test -s data/lookup/annotations.json
test -s data/lookup/accessions.json
```

Alternatively, download those three files separately and place them directly
in `data/lookup/`.

### 3. Run the bundled test

```bash
./scripts/minimal_pipeline.sh fasta_test/test.fa
```

The first run creates `venv/` and installs the Python and PyTorch dependencies.
It therefore requires internet access and takes longer than subsequent runs.

### 4. Find the result

Each invocation creates a timestamped directory:

```text
outputs_<YYYYMMDD_HHMMSS>/
├── results.csv
├── query_embeddings.npz
├── failed_sequences.csv
├── fantasia_config.yaml
├── run_metadata.yaml
├── run.log
├── k.<N>.results.csv       # automatically written when k > 1
├── topgo/                  # only when --topgo is requested
└── tmp/
```

`results.csv` is the main final annotation table. Check run completion and
coverage with:

```bash
latest=$(find . -maxdepth 1 -type d -name 'outputs_*' | sort | tail -n 1)
head -n 5 "$latest/results.csv"
cat "$latest/run_metadata.yaml"
```

### 5. Run your own proteome

```bash
./scripts/minimal_pipeline.sh data/my_proteome.faa.gz
```

Plain and gzip-compressed FASTA files are supported. Gzipped FASTA input is
decompressed while reading; you do not need to unpack it first.

## Requirements

- Linux or macOS
- Python 3.9 or newer; Python 3.10+ is recommended
- Git
- internet access on the first run for Python packages and model weights
- the three-file Lite lookup bundle
- NVIDIA GPU strongly recommended for embedding and lookup; CPU is a slower fallback

The extracted lookup bundle is approximately 3.1 GB; the compressed archive is
approximately 1.8 GB. Allow at least 5 GB while retaining both. Model caches,
the virtual environment, query embeddings, and outputs require additional
space. For a first ProtT5 proteome run, approximately 12 GB of free disk space
is a practical starting point; large or multi-model projects need more.

There is no universal RAM/VRAM minimum because sequence lengths and batch sizes
matter. Start with at least 16 GB system RAM. The supplied tuned GPU profile was
tested on a 24 GB GPU; reduce `--embed-batch-size`,
`--model-batch-sizes`, and `--max-batch-residues` on smaller GPUs. CPU execution
is supported only as a substantially slower fallback.

## What the lookup bundle contains

The current reference is a UniProt November 2025 snapshot generated with PIS
v3.1.0. It contains proteins with experimental/curated evidence codes: EXP,
IDA, IPI, IMP, IGI, IEP, TAS, and IC.

| File | Purpose |
|---|---|
| `lookup_table.npz` | Reference embeddings searched by the selected model space |
| `annotations.json` | GO terms, namespaces, descriptions, and evidence codes |
| `accessions.json` | Mapping from internal reference IDs to UniProt accessions and metadata |

The identifier and annotation files each contain 124,397 reference entries and
621,024 GO annotation rows. Model-space coverage differs slightly:

| Reference model space | Embeddings | Dimension |
|---|---:|---:|
| ESM-2 | 124,363 | 1,280 |
| ProstT5 | 124,248 | 1,024 |
| ProtT5 | 124,239 | 1,024 |
| Ankh3-Large | 124,338 | 1,536 |
| ESM3c | 124,397 | 1,152 |

See the [lookup-bundle audit](docs/lookup_bundle_audit_2026-04-26.md) for
structural checks and release-drift spot checks.

## Models: embedder versus lookup table

The bundled reference table contains five embedding spaces, but the current
Lite embedder directly exposes two:

| Capability | Models |
|---|---|
| Built-in end-to-end Lite embedding and lookup | `prot_t5`, `ankh3` |
| Reference spaces present in the lookup bundle | ESM-2, ProstT5, ProtT5, Ankh3-Large, ESM3c |

ProtT5 is the recommended default. Ankh3-Large is available for model
comparison but is slower and generally requires more memory. The additional
reference spaces can be searched only when compatible externally generated
query embeddings are supplied in the expected archive format.

## Command-line defaults

The helper script makes a few tuned choices. The Python entry point itself uses
the following defaults:

| Setting | Python default | Meaning |
|---|---|---|
| embedding model | `prot_t5` | Built-in query embedder |
| device | automatic, GPU-first | CUDA when `nvidia-smi` is available; otherwise CPU fallback |
| lookup device | follows embedding device | Override with `--use-gpu-lookup` or `--no-gpu-lookup` |
| distance | `cosine` | Also accepts `euclidean` |
| neighbors (`k`) | `1` | `--limit-per-entry` |
| sequence truncation | `0` (disabled) | `--length-filter` / `--max-sequence-length` |
| FASTA chunk size | `500` | Sequences per temporary FASTA chunk |
| sequence queue package | `100` | Outer embedding package size |
| embedding batch size | `8` | Maximum sequences per forward batch |
| padded residues per batch | `12,000` | Memory-aware embedding limit |
| lookup batch size | `256` | Query embeddings per lookup batch |
| precision | `float32` | Mixed precision is not enabled by default |
| raw neighbor output | off for `k=1` | Automatically `k.<N>.results.csv` for `k>1` |
| TopGO export | off | Enable with `--topgo` |
| output directory | timestamped | `outputs_<YYYYMMDD_HHMMSS>/` |

`scripts/minimal_pipeline.sh` uses ProtT5, cosine, `k=1`, queue size 100,
embedding/model batch size 4, and lookup batch size 1,024. It does not truncate
sequences.

Show every option with:

```bash
python3 src/fantasia_pipeline.py --help
```

## Common workflows

### Explicit GPU run

```bash
python3 src/fantasia_pipeline.py \
  --device cuda \
  --use-gpu-lookup \
  --serial-models \
  --embed-models prot_t5 \
  --model-batch-sizes prot_t5=4 \
  --lookup-batch-size 1024 \
  data/my_proteome.faa.gz
```

### CPU fallback

```bash
python3 src/fantasia_pipeline.py \
  --device cpu \
  --no-gpu-lookup \
  --embed-models prot_t5 \
  data/my_proteome.faa.gz
```

### Retrieve five donors

```bash
python3 src/fantasia_pipeline.py \
  --embed-models prot_t5 \
  --distance-metric cosine \
  --limit-per-entry 5 \
  data/my_proteome.faa.gz
```

This writes consolidated annotations to `results.csv` and neighbor-level rows
to `k.5.results.csv`.

### Compare both built-in models

```bash
python3 src/fantasia_pipeline.py \
  --serial-models \
  --embed-models "prot_t5 ankh3" \
  data/my_proteome.faa.gz
```

Use `--serial-models` for multi-model runs to avoid loading both embedders at
the same time.

### Reuse query embeddings

```bash
python3 src/fantasia_pipeline.py \
  --reuse-embeddings \
  --embeddings-npz outputs_<timestamp>/query_embeddings.npz \
  --limit-per-entry 5 \
  data/my_proteome.faa.gz
```

The FASTA is still required for identifiers and run accounting. Reused
embeddings must match the selected lookup model space and accession order.

### Embedding only

The embedding-only entry point does not perform GO lookup:

```bash
python3 src/generate_embeddings.py \
  --fasta data/my_proteome.faa.gz \
  --output my_embeddings.npz \
  --models prot_t5 \
  --device cuda
```

Export selected ProtT5 layers:

```bash
python3 src/generate_embeddings.py \
  --fasta data/my_proteome.faa.gz \
  --output my_layers.npz \
  --models prot_t5 \
  --layer-indices 0 12 24
```

For embedding-only layer numbering, ProtT5 layer 24 and Ankh3 layer 48 are the
final exported transformer layers. Intermediate-layer export does not change
the standard last-layer lookup workflow.

## Output interpretation

### `results.csv`

This is the consolidated annotation table. It keeps the best supporting row
for each `(query, model, GO term, category)` combination. A protein can
therefore occupy multiple rows because it can receive multiple GO terms.

### `k.<N>.results.csv` or `--raw-results-csv`

This is the neighbor-level table before GO-term consolidation. It preserves
`hit_rank` and is the correct file for inspecting donor proteins or comparing
the selected neighborhood when `k>1`.

Important columns include:

| Column | Meaning |
|---|---|
| `query_accession` | Input FASTA identifier |
| `hit_rank` | Neighbor rank in raw output; 1 is nearest |
| `reference_id` | Internal lookup-table donor identifier |
| `uniprot_accession` | Donor UniProt accession |
| `model_key` | Embedding space used |
| `distance` | Query-reference embedding distance; lower is closer |
| `reliability_index` | Distance-derived score clipped to `[0,1]`; higher is better |
| `go_id` | Transferred GO identifier |
| `category` | `F` molecular function, `P` biological process, or `C` cellular component |
| `evidence_codes` | Evidence supporting the GO term in the donor |

### Other files

- `query_embeddings.npz`: reusable query embeddings.
- `failed_sequences.csv`: skipped sequences and failure reasons.
- `fantasia_config.yaml`: effective lookup configuration.
- `run_metadata.yaml`: command, Git branch/commit, resolved parameters, device,
  output paths, coverage, and stage timings.
- `run.log`: complete chronological console output.
- `topgo/`: optional TopGO tables generated only with `--topgo`.

`run_metadata.yaml` does not currently guarantee an immutable Hugging Face
model revision or weight checksum. For exact model-level reproducibility,
record the resolved cache snapshot and package environment with the run.

## Performance and memory tuning

These settings control different stages:

- `--chunk-size`: number of FASTA records in each temporary chunk.
- `--sequence-queue-package`: outer group passed to embedding generation.
- `--embed-batch-size`: sequences in a model forward pass.
- `--max-batch-residues`: cap on padded residues in a forward batch.
- `--lookup-batch-size`: query vectors compared at once; it does not reduce the
  number of reference proteins searched.

If CUDA runs out of memory:

1. use one model with `--serial-models`;
2. reduce `--model-batch-sizes prot_t5=2`;
3. reduce `--embed-batch-size`;
4. reduce `--max-batch-residues`;
5. only as an explicit methodological choice, set a positive
   `--length-filter`.

The default `--length-filter 0` is uncapped. Proteins that still cannot be
embedded are listed in `failed_sequences.csv`; the rest of the run continues.

## Revalidated benchmark

Current ProtT5, cosine, `k=1`, float32, full-length runs on an NVIDIA RTX 3090 Ti
(24 GB) produced:

| Proteome | Input proteins | Wall time | Coverage | Failed |
|---|---:|---:|---:|---:|
| *Paratomella rubra* | 20,223 | 22 min 39 s | 99.72% | 56 |
| *Caenorhabditis elegans* | 19,831 | 23 min 34 s | 99.75% | 49 |
| Mouse | 21,852 | 35 min 43 s | 99.31% | 150 |

Runtime depends strongly on sequence-length distribution, cache state, GPU,
software versions, and batch settings. Treat these values as scale estimates,
not hardware-independent guarantees.

## Repository layout

```text
FANTASIA-Lite/
├── README.md
├── data/lookup/                 # user-created; excluded from Git
├── fasta_test/                  # bundled test inputs
├── scripts/
│   ├── minimal_pipeline.sh
│   ├── default_last_layer.sh
│   └── embedding_only_with_layers.sh
├── src/
│   ├── fantasia_pipeline.py
│   ├── generate_embeddings.py
│   ├── fantasia_no_db.py
│   ├── fantasia_no_db_gpu.py
│   └── pipeline_timing_analyzer.py
└── outputs_<timestamp>/         # generated; excluded from Git
```

## Troubleshooting

### Lookup files are reported missing

Confirm that the files are directly under `data/lookup/`, not inside an extra
archive directory:

```bash
ls -lh data/lookup/{lookup_table.npz,annotations.json,accessions.json}
```

### CUDA is unavailable

Run `nvidia-smi`. If it is unavailable or incompatible with installed PyTorch,
use `--device cpu --no-gpu-lookup` or provide the appropriate PyTorch package
index with `--torch-index`.

### CUDA out of memory

Reduce model and residue batch sizes. Keep `--serial-models` enabled for
multi-model work. Inspect `failed_sequences.csv` for individual long proteins
that could not be embedded.

### A run appears to reinstall packages

Lite creates and manages `venv/`, then verifies/installs dependencies when the
pipeline starts. Delete `venv/` only when you deliberately want a clean
environment rebuild.

### The final table has multiple rows per protein

This is expected: `results.csv` has one row per retained query/model/GO/category
combination, not one row per query protein. Use `query_accession` to group rows.

## Citation and acknowledgements

Please cite the software repository, the version-specific Zenodo lookup record
used in the analysis, and the relevant FANTASIA publications:

- [Performance of protein language models in model organisms](https://doi.org/10.1093/nargab/lqae078)
- [Application of FANTASIA to functional annotation of dark proteomes](https://doi.org/10.1038/s42003-025-08651-2)
- [Protocol explaining FANTASIA](https://doi.org/10.1007/978-1-0716-4623-6_8)

FANTASIA Lite derives from
[FANTASIA Full](https://github.com/CBBIO/FANTASIA) and incorporates methods from
[GOPredSim](https://github.com/Rostlab/goPredSim). Protein language models are
distributed through [Hugging Face](https://huggingface.co/).
