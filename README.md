
# FANTASIA Lite V1.1

**FANTASIA Lite V1.1** is a streamlined, standalone version of the full [FANTASIA pipeline](https://github.com/CBBIO/FANTASIA), designed for fast and efficient Gene Ontology (GO) annotation of protein sequences from local FASTA files using embedding comparisons.

FANTASIA Lite generates deep learning embeddings and performs nearest-neighbor annotation transfer, while intentionally removing the service stack used by the full project. The bundled lookup table covers five reference embedding spaces, and the built-in Lite embedder can generate matching query embeddings for all five model families exposed through the local CLI.

The simplest way to run it is:

```bash
./scripts/minimal_pipeline.sh your_proteins.fa
```

Before any annotation run can work, create `data/lookup/`, then download the
required lookup data from [Zenodo](https://doi.org/10.5281/zenodo.19742925). These
directories are intentionally not tracked by Git and are not created by a
fresh clone:

```bash
mkdir -p data/lookup
```

Unlike the full FANTASIA pipeline, **FANTASIA Lite does not require PostgreSQL, RabbitMQ, or a database-backed orchestration layer**. It runs locally from flat files:
- `lookup_table.npz` for reference embeddings
- `annotations.json` for GO annotations
- `accessions.json` for accession mapping

## Warning

The tradeoff is that Lite is simpler to deploy but can be slower if embedding and lookup are not tuned well. This repository now includes GPU-aware embedding and lookup controls intended to recover as much performance as possible without reintroducing external services.

## Scope and Purpose

The main purpose of Lite V1.1is to provide a fast local annotator that can be dropped into other pipelines with minimal setup. The default fast path is intentionally simple: ProtT5 embeddings, cosine lookup, and `k=1` transfer. More advanced configuration is still available when you need multi-model runs, layered embedding export, larger `k`, or other custom settings.

This repository is ideal for users who want:
- Lightweight, local annotation of protein FASTA files
- No external database dependencies
- No PostgreSQL or RabbitMQ setup
- Simple setup and automated environment management
- High-quality functional annotation using experimental evidence from UniProt

For advanced features, large-scale annotation, or integration with external databases, see the full [FANTASIA repository](https://github.com/CBBIO/FANTASIA).

## What's New In V1.1

Compared with the earlier Lite V0 branch, Lite V1.1includes several practical improvements:

- Faster end-to-end execution through batched embedding inference and a merged one-pass lookup flow
- Lookup and merge stages run through the managed virtual environment, avoiding dependency leaks from the host Python interpreter
- Reuse mode for rerunning lookup from an existing embeddings archive with `--reuse-embeddings`
- Optional TopGO export, now disabled by default unless `--topgo` is requested
- Clearer GPU tuning controls for embedding and lookup batch sizes
- Full-length embeddings by default, with truncation disabled unless explicitly requested
- Better reporting of skipped sequences and end-of-run coverage summaries
- End-to-end local embedding support for `esm2`, `prost_t5`, `prot_t5`, `ankh3`, and `esmc`
- ESM-C installation avoids `torchvision` and `torchtext`, preventing fragile torch/torchvision wheel mismatches
- Embedding-only workflows via `src/generate_embeddings.py`, separate from lookup
- Optional export of selected layers or all available layers from the Lite embedder
- Helper scripts for a default last-layer end-to-end run and an embedding-only layered export run
- Updated documentation that explains the separation between embedding support and lookup-bundle coverage

## Installation Requirements
- Python 3.12 for full model support. ESM-C uses `esm==3.2.3`, which currently supports Python `>=3.12,<3.13`; all-five-model runs therefore require Python 3.12.
- Required lookup bundle (`lookup_table.npz`, `annotations.json`, `accessions.json`) from [Zenodo](https://doi.org/10.5281/zenodo.19742925) placed in `data/lookup/` before running the pipeline
- Internet connection for automatic dependency installation
- Sufficient disk space for outputs, embeddings, model downloads, and the managed virtual environment. CUDA environments and five-model runs can require several GB.
- Optional NVIDIA GPU with CUDA support. The default installer uses the current package index unless `--torch-index` or `TORCH_INDEX` is explicitly set.
- Git (for cloning the repository)
- `wget` or `curl` (for downloading the lookup bundle)


## Lookup Table Details


The FANTASIA Lite V1.1lookup table is built from the **UniProt November 2025 release** and includes only proteins with experimental evidence, ensuring high-quality functional annotations. All data was generated using PIS v3.1.0, the internal system used to extract and preprocess UniProt, PDB, and GOA data.

**Lookup bundle Zenodo concept DOI:** [10.5281/zenodo.19742925](https://doi.org/10.5281/zenodo.19742925)

Use this concept DOI to cite the lookup table or access its latest published
version. Individual releases also have version-specific DOIs.

### Core Statistics
- **Reference entries in `accessions.json`**: 124,397
- **Annotated reference entries in `annotations.json`**: 124,397
- **Total GO annotation rows in `annotations.json`**: 621,024

### Embedding Coverage
- **ESM-2**: 124,363 embeddings
- **ProstT5**: 124,248 embeddings
- **ProtT5-XL-UniRef50**: 124,239 embeddings
- **Ankh3-Large**: 124,338 embeddings
- **ESM3c**: 124,397 embeddings

These are the measured per-model embedding counts in the current lookup bundle. `annotations.json` and `accessions.json` both contain `124,397` entries, while the individual model spaces in `lookup_table.npz` can differ slightly in coverage.

### Package Contents
The Zenodo deposition provides the three essential lookup files both as the
single convenience archive `fantasia_lite_data_folder.tgz` and as separate
downloads. The archive contains:

1. **`lookup_table.npz`**
   - Precomputed protein embeddings for ESM-2, ProstT5, ProtT5, Ankh3-Large, and ESM3c
   - Last-layer compressed embeddings for all reference sequences
   - Enables fast nearest-neighbor search during annotation
   - Format: NumPy .npz archive

2. **`annotations.json`**
   - GO annotations of the reference proteins
   - Experimentally supported GO terms by category:
     - **F**: Molecular Function
     - **P**: Biological Process
     - **C**: Cellular Component
   - Format: JSON mapping from proteins to their GO terms

3. **`accessions.json`**
   - Mapping of internal indices to UniProt accessions
   - Contains UniProt ID, metadata, and sequence length
   - Allows the pipeline to retrieve source identifiers
   - Format: JSON list/dict

### How to Interpret the Lookup Bundle and Outputs

- `lookup_table.npz` is the reference embedding space. It does not contain GO terms by itself. It stores the precomputed vectors that query embeddings are compared against.
- `annotations.json` is the annotation layer. After the nearest reference proteins are found, this file is used to expand those hits into GO terms, GO categories, descriptions, and evidence codes.
- `accessions.json` is the identifier layer. It maps the internal lookup-table positions back to UniProt accessions and metadata.
- Lookup is always performed against the full reference table for the selected model space. `--lookup-batch-size` controls how many query embeddings are processed together, not how many reference proteins are searched.
- In Lite, `k` is controlled by `--limit-per-entry`. With `k=1`, the pipeline keeps only the top reference hit per query. With larger `k`, it keeps more nearest neighbors before GO-term consolidation.
- `results.csv` is the consolidated final annotation table. It keeps the best supporting row for each `(query, GO term, category)` combination.
- `raw_results.csv` is optional and preserves the neighbor-level rows before GO-term consolidation. It is especially useful when `k > 1`.
- The default fast Lite V1.1path is last-layer lookup with cosine distance, ProtT5 embeddings, and `k=1`. Layered embedding export is separate from lookup.
- Layer selection is intended for embedding-only workflows. The bundled lookup workflow uses the standard last-layer embeddings and does not switch lookup behavior based on `--layer-indices`.

### GO Evidence Codes (Experimental Only)
The lookup table includes only high-confidence experimental annotations:
- **EXP** — Inferred from Experiment
- **IDA** — Inferred from Direct Assay
- **IPI** — Inferred from Physical Interaction
- **IMP** — Inferred from Mutant Phenotype
- **IGI** — Inferred from Genetic Interaction
- **IEP** — Inferred from Expression Pattern
- **TAS** — Traceable Author Statement
- **IC** — Inferred by Curator

No database server or external dependencies are required.

## Quick Start

FANTASIA Lite V1.1is designed first as a fast local annotator that is easy to integrate into other pipelines, while still supporting more configurable research-style runs when needed.

The easiest default path is:

```bash
./scripts/minimal_pipeline.sh your_proteins.fa
```

This script keeps the intended Lite defaults explicit:
- `prot_t5`
- cosine lookup
- `k=1`
- full precision (`float32`)
- automatic CPU/GPU selection

In Lite, `k` is controlled by `--limit-per-entry`.

Beyond that, Lite V1.1provides two main Python entrypoints:

### 1. Standard Pipeline (`fantasia_pipeline.py`)
For processing protein sequences and obtaining GO annotations:

Lite V1.1keeps embedding and lookup in full precision (`float32`) for comparability with the full FANTASIA workflow. Mixed precision is not enabled by default.

```bash
# Basic usage - single model annotation
fantasia_pipeline your_proteins.fa

# Recommended general usage with a single model
fantasia_pipeline --serial-models --embed-models prot_t5 your_proteins.fa

# Recommended fast GPU usage without sequence truncation
fantasia_pipeline \
    --device cuda \
    --use-gpu-lookup \
    --serial-models \
    --embed-models prot_t5 \
    --sequence-queue-package 100 \
    --embed-batch-size 4 \
    --model-batch-sizes prot_t5=4 \
    --lookup-batch-size 1024 \
    your_proteins.fa

# Multiple models (slower but more comprehensive; keep --serial-models enabled)
python3.12 src/fantasia_pipeline.py \
    --serial-models \
    --embed-models "esm2 prost_t5 prot_t5 ankh3 esmc" \
    your_proteins.fa

# RTX 5090 / Blackwell CUDA example with conservative per-model batches
python3.12 src/fantasia_pipeline.py \
    --device cuda \
    --use-gpu-lookup \
    --serial-models \
    --embed-models "esm2 prost_t5 prot_t5 ankh3 esmc" \
    --model-batch-sizes esm2=1 prost_t5=1 prot_t5=4 ankh3=4 esmc=1 \
    your_proteins.fa

# Advanced configuration
fantasia_pipeline \
    --embed-models prot_t5 \
    --limit-per-entry 5 \
    --results-csv my_results.csv \
    your_proteins.fa
```

### Helper Scripts

Three ready-to-run helper scripts are included under `scripts/`:

These helper scripts use paths relative to the repository root. Run them from inside the cloned `FANTASIA-Lite` directory. If needed, you can override the Python interpreter with `PYTHON_BIN=/path/to/python`. Use a Python 3.12 interpreter for ESM-C or all-five-model runs.

- `scripts/minimal_pipeline.sh`
  This is the main Lite V1.1entrypoint for fast annotation and pipeline integration. It runs the standard end-to-end workflow with the intended Lite defaults:
  - ProtT5 only
  - cosine distance
  - `k=1`
  - automatic device detection
  - full precision
  - no truncation unless you explicitly request it in another command

  Example:

```bash
./scripts/minimal_pipeline.sh fasta_test/PRUB1_longiso.pep
```

- `scripts/default_last_layer.sh`
  Runs the standard end-to-end Lite workflow with the current recommended tuned settings for a 24 GB-class GPU:
  - ProtT5 only
  - last layer only
  - GPU lookup enabled
  - full-precision embeddings
  - TopGO disabled by default

  Example:

```bash
./scripts/default_last_layer.sh fasta_test/PRUB1_longiso.pep
```

- `scripts/embedding_only_with_layers.sh`
  Runs embedding generation only, with no lookup, and exports the default last-layer embeddings plus either:
  - all transformer layers when no extra layer arguments are given
  - selected layers when you pass them explicitly
  For ProtT5, layer numbering follows the model outputs directly: `0` is the earliest exported hidden-state representation, while `24` is the final encoder layer. The standard last-layer embedding is also always written separately as `Prot-T5_embeddings`.
  For Ankh3, the same ordering rule applies: `0` is the earliest exported hidden-state representation and `48` is the final exported layer. The standard last-layer embedding is also written separately as `Ankh3-Large_embeddings`.
  These layer options should be used in embedding-only mode, not as a replacement for the default lookup workflow.

  Examples:

```bash
# Export all available layers
./scripts/embedding_only_with_layers.sh fasta_test/PRUB1_longiso.pep outputs_layers.npz

# Export only selected layers
./scripts/embedding_only_with_layers.sh fasta_test/PRUB1_longiso.pep outputs_layers.npz 0 8 16 24
```

**Key Options:**
- `--embed-models`: Choose models (`esm2`, `prost_t5`, `prot_t5`, `ankh3`, `esmc`) - default: `prot_t5`. Aliases such as `ESM`, `Prost-T5`, `Prot-T5`, `Ankh3-Large`, and `ESM3c` are accepted.
- `--serial-models`: Process requested models one after another instead of as one combined model-group. For a single model, this usually makes no practical difference. For more than one model, it is the recommended and safer setting because it reduces GPU/CPU memory pressure and avoids loading multiple embedders at the same time.
- `--limit-per-entry N`: Return top N annotations per sequence (default: 1)
- `--raw-results-csv PATH`: Optional raw neighbor-level output before GO-term consolidation. If omitted, Lite writes only `results.csv` unless `--limit-per-entry > 1`, in which case a raw file is created automatically as `k.<N>.results.csv`
- `--topgo`: Optional. Generate TopGO files after lookup. TopGO export is disabled by default
- `--distance-metric {cosine,euclidean}`: Lookup metric. Lite V1.1defaults to `cosine`, but `euclidean` is also supported
- `--use-gpu-lookup`: Force GPU nearest-neighbor lookup when CUDA is available
- `--lookup-batch-size N`: Number of query embeddings compared per lookup batch
- `--sequence-queue-package N`: Outer packaging size before embedding forward passes
- `--embed-batch-size N`: Default forward-pass batch size for embeddings
- `--model-batch-sizes MODEL=N ...`: Per-model embedding batch size overrides, for example `esm2=1 prost_t5=1 prot_t5=4 ankh3=4 esmc=1`
- `--length-filter N`: Optional truncation before embedding; `0` disables truncation and is the default

## Getting Started

FANTASIA Lite V1.1is installed by cloning the repository, placing the lookup bundle in `data/lookup/`, and running a simple setup check from inside the cloned folder.

### Step 1: Clone the Repository
FANTASIA Lite annotates protein FASTA files by embedding each query, searching
a local reference table, and transferring experimentally supported Gene
Ontology (GO) terms from nearby reference proteins.

Lite runs from flat files without PostgreSQL or RabbitMQ. An NVIDIA GPU is
strongly recommended; CPU mode is available as a slower fallback. For
service-backed reference management and all five built-in embedders, use
[FANTASIA Full](https://github.com/CBBIO/FANTASIA).

## Run your first search

### 1. Clone Lite

```bash
git clone https://github.com/CBBIO/FANTASIA-Lite.git
cd FANTASIA-Lite
mkdir -p data/lookup
```

### 2. Install the reference bundle

Download `fantasia_lite_data_folder.tgz` from
[Zenodo](https://doi.org/10.5281/zenodo.19742925), place it in the repository
root, and extract it:

```bash
tar -xzf fantasia_lite_data_folder.tgz
test -s data/lookup/lookup_table.npz
test -s data/lookup/annotations.json
test -s data/lookup/accessions.json
```

The archive contains the `data/lookup/` path. Alternatively, download those
three files separately into `data/lookup/`.

### 3. Run the bundled example

```bash
./scripts/minimal_pipeline.sh fasta_test/test.fa
```

The default minimal script uses `prot_t5`. If your system default `python3` is not Python 3.12 and you plan to run ESM-C or all five models later, use an explicit interpreter for those runs:

```bash
PYTHON_BIN=/path/to/python3.12 ./scripts/minimal_pipeline.sh fasta_test/test.fa
```

The first run creates `venv/`, installs dependencies, and downloads ProtT5.
Lite automatically uses CUDA when available.

### 4. Run your proteome

```bash
./scripts/minimal_pipeline.sh data/my_proteome.faa.gz
```

Plain and gzip-compressed FASTA files are accepted. Gzip input is decompressed
while reading.

## Default behavior

The helper script runs ProtT5 with cosine distance, `k=1`, float32 embeddings,
no sequence truncation, GPU lookup when CUDA is available, and TopGO disabled.
Each invocation creates `outputs_<timestamp>/`.

| To change | Option |
|---|---|
| retrieve five donors | `--limit-per-entry 5` |
| force GPU | `--device cuda --use-gpu-lookup` |
| force CPU fallback | `--device cpu --no-gpu-lookup` |
| use Ankh3-Large | `--serial-models --embed-models ankh3` |
| reuse embeddings | `--reuse-embeddings --embeddings-npz FILE` |
| enable TopGO | `--topgo` |
| truncate explicitly | `--length-filter N` (`0` is uncapped) |

Show all options:

```bash
python3 src/fantasia_pipeline.py --help
```

## Results and provenance

```text
outputs_<timestamp>/
├── results.csv
├── query_embeddings.npz
├── failed_sequences.csv
├── fantasia_config.yaml
├── run_metadata.yaml
├── model_provenance.yaml
├── run.log
├── k.<N>.results.csv       # when k > 1
└── topgo/                  # only with --topgo
```

- `results.csv` contains final consolidated GO assignments. A protein can have
  multiple rows because it can receive multiple GO terms.
- `k.<N>.results.csv` is the donor-level table before GO consolidation. It
  preserves `hit_rank`, reference IDs, UniProt accessions, distances, and GO
  terms.
- `failed_sequences.csv` reports proteins that could not be embedded.
- `run_metadata.yaml` records the command, Git commit, parameters, device,
  coverage, and stage timings.
- `model_provenance.yaml` records model repositories, enforced Hugging Face
  revisions, lookup keys, device, Python, and package versions.

Lite pins both tokenizer and model downloads to these commits:

| Lite model | Repository | Enforced revision |
|---|---|---|
| ProtT5 | `Rostlab/prot_t5_xl_uniref50` | `973be27c52ee6474de9c945952a8008aeb2a1a73` |
| Ankh3-Large | `ElnaggarLab/ankh3-large` | `2be091622e8a393f0ef21735070084123c874b6e` |

Keep both YAML files with published results.

## Common commands

GPU run with tuned batches:

```bash
python3 src/fantasia_pipeline.py \
  --device cuda --use-gpu-lookup --serial-models \
  --embed-models prot_t5 --model-batch-sizes prot_t5=4 \
  --lookup-batch-size 1024 data/my_proteome.faa.gz
```

Retrieve five donors and retain the raw neighborhood:

```bash
python3 src/fantasia_pipeline.py \
  --embed-models prot_t5 --distance-metric cosine \
  --limit-per-entry 5 data/my_proteome.faa.gz
```

Reuse an existing query archive:

```bash
python3 src/fantasia_pipeline.py \
  --reuse-embeddings \
  --embeddings-npz outputs_<timestamp>/query_embeddings.npz \
  --limit-per-entry 5 data/my_proteome.faa.gz
```

Generate embeddings without GO lookup:

```bash
./scripts/default_last_layer.sh fasta_test/PRUB1_longiso.pep
./scripts/minimal_pipeline.sh fasta_test/UP000001940_6239.fasta
./scripts/minimal_pipeline.sh fasta_test/MUSM_10090.fasta
```

This corresponds to:
- `prot_t5`
- cosine lookup
- `k=1`
- full precision (`float32`)
- full-length embeddings by default
- GPU lookup enabled
- tuned batch sizes for a 24 GB-class GPU (`prot_t5=4`)

### NVIDIA GeForce RTX 3090 Ti (24 GB VRAM)

| Dataset | Sequences | Model | Runtime | Rate (seq/s) | Coverage | Failed |
|---------|-----------|-------|---------|--------------|----------|--------|
| **PRUB1_longiso.pep** | 20,223 | ProtT5 | 22m 38.60s | 14.89 | 99.72% | 56 |
| **UP000001940_6239.fasta** | 19,831 | ProtT5 | 23m 33.86s | 14.03 | 99.75% | 49 |
| **MUSM_10090.fasta** | 21,852 | ProtT5 | 35m 43.23s | 10.20 | 99.31% | 150 |

### PRUB1 `k=1` vs `k=5` on GPU

The following comparison keeps the device, model family, and full-length setting aligned, and changes only the lookup depth:

| Dataset | k | Runtime | Rate (seq/s) | Coverage | Failed |
|---------|---|---------|--------------|----------|--------|
| **PRUB1_longiso.pep** | 1 | 22m 38.60s | 14.89 | 99.72% | 56 |
| **PRUB1_longiso.pep** | 5 | 23m 19.46s | 14.45 | 99.83% | 34 |

For the revalidated `k=5` PRUB1 run, the recorded stage split was approximately 22m 41.35s embedding, 37.58s lookup, and 0.53s post-processing. In practice, increasing `k` from 1 to 5 added only a modest overhead because embedding remains the dominant cost.

### Five-Sequence Validation Tests

Small validation runs on the first 5 sequences of `fasta_test/test.fa` were used to confirm that Lite V1 behaves correctly in CPU-only mode, GPU `k=5` mode, and layered embedding export mode:

| Test | Device | Settings | Runtime | Notes |
|------|--------|----------|---------|-------|
| End-to-end lookup | CPU | ProtT5, cosine, `k=1` | 118.90s | Confirms Lite works without GPU |
| End-to-end lookup | GPU | ProtT5, cosine, `k=5` | 64.43s | Not directly comparable with the CPU row because the device changed |
| Embedding-only layered export | GPU | ProtT5 layers `0, 12, 24` | 12.97s | Wrote `Prot-T5_layer_0_embeddings`, `Prot-T5_layer_12_embeddings`, and `Prot-T5_layer_24_embeddings`, each with shape `(5, 1024)`; here `24` is the final encoder layer, not `0` |
| Embedding-only layered export | GPU | Ankh3 all layers | 70.58s | Wrote `Ankh3-Large_layer_0_embeddings` through `Ankh3-Large_layer_48_embeddings`; representative shapes confirmed for `layer_0`, `layer_24`, and `layer_48`, each `(5, 1536)` |

Notes:
- The skipped sequences were extreme long-protein CUDA OOM cases, not ordinary proteins.
- This benchmark uses the current Lite V1 merged-lookup flow rather than the older per-chunk lookup behavior.
- For the revalidated C. elegans run, the recorded stage split was approximately 22m 42s embedding, 46.46s lookup, and 1.09s post-processing.
- For the revalidated mouse run, the recorded stage split was approximately 35m 05s embedding, 30.79s lookup, and 0.86s post-processing.
- Additional benchmark tables for other datasets and hardware should be regenerated before being treated as representative of Lite V1.

## Advanced Usage

### Environment Management
The pipeline automatically manages Python virtual environments:
```bash
# Virtual environment is created automatically in venv/
# To clean up and force rebuild:
rm -rf venv/
./scripts/minimal_pipeline.sh your_file.fa  # Will recreate venv automatically
```

### Batch Processing
For processing multiple files systematically:
```bash
# Process multiple specific files
python3 src/pipeline_timing_analyzer.py \
    --files file1.fa file2.fa file3.fa \
    --model prot_t5

# Process all files in a directory
python3 src/pipeline_timing_analyzer.py \
    --fasta-dir my_proteomes/ \
    --report-csv batch_results.csv
```

The timing analyzer is useful for benchmarking and regression testing. In Lite V1 it now also reads per-run stage timings from `run_metadata.yaml`, so benchmark reports can separate embedding time, lookup time, and post-processing time.

### Memory Optimization
For large files or limited memory systems:
```bash
# Use a single model, smaller chunks, and explicit serial processing
fantasia_pipeline \
    --serial-models \
    --embed-models prot_t5 \
    --chunk-size 200 \
    large_proteome.fa
```

For the current Lite V1 fast path, `prot_t5` with `--model-batch-sizes prot_t5=4` is the best starting point on a 24 GB-class GPU. If you see CUDA OOM skips, reduce the model batch size further before changing lookup settings.

## Supported Models

### Embedder vs Lookup

In Lite V1, the embedding step and lookup step are intentionally separated:

- The **embedder** is the local model-inference component in `src/generate_embeddings.py`
- The **lookup** is the nearest-neighbor transfer step against the flat-file reference bundle in `data/lookup/`

These two layers are related, but they do not currently expose exactly the same model set through the Lite CLI.

### Built-in Lite Embedder

The built-in Lite embedder currently supports:
- **`prot_t5`**: Protein T5 model (recommended, good balance of speed and accuracy)
- **`ankh3`**: ANKH large protein language model (slower but potentially more accurate)

The Lite embedder can now also export:
- default last-layer embeddings
- selected intermediate layers with `--layer-indices`
- all available layers with `--all-layers`

This embedding-only mode is independent of lookup. In other words, Lite can generate embeddings locally for the models above, with optional layer export, even when you do not want to run annotation lookup.

### Lookup Bundle Coverage

The bundled Lite lookup table currently contains last-layer reference embeddings for:
- **ESM-2**
- **ProstT5**
- **ProtT5**
- **Ankh3-Large**
- **ESM3c**

So the lookup bundle is broader than the built-in Lite embedder. The lookup side is not limited to only two model spaces; the current local embedding CLI is the narrower component.

### Current End-to-End Lite Pipeline

The current built-in end-to-end Lite pipeline is intended for the last-layer models that the Lite embedder can generate directly today:
- **ProtT5**
- **Ankh3**

If you provide externally generated embeddings that match the lookup bundle keys and format, the lookup layer itself is separate and can in principle operate on the additional lookup-table models as well. So the practical distinction is:
- built-in Lite embedding CLI: currently `prot_t5` and `ankh3`, with optional layer export
- bundled lookup table: `ESM-2`, `ProstT5`, `ProtT5`, `Ankh3-Large`, and `ESM3c`, using last-layer reference embeddings

## File Format Support
- **Input**: FASTA files (`.fa`, `.faa`, `.fasta`) and gzip-compressed versions (`.fa.gz`, `.fasta.gz`)
- **Output**: CSV files for results, NPZ files for embeddings, TXT files for TopGO compatibility

## Troubleshooting

**Lookup files missing:** ensure the three files are directly under
`data/lookup/`, not in an additional nested directory.

**CUDA unavailable:** check `nvidia-smi`; otherwise use
`--device cpu --no-gpu-lookup` as a slower fallback.

**CUDA out of memory:** use one model with `--serial-models`, then reduce
`--model-batch-sizes`, `--embed-batch-size`, and `--max-batch-residues`.
Failures for individual proteins remain listed in `failed_sequences.csv`.

**Several rows per protein:** expected—`results.csv` has one row per retained
query/model/GO/category combination.

## Further documentation

- [Detailed usage reference](docs/usage_reference.md)
- [Lookup-bundle audit](docs/lookup_bundle_audit_2026-04-26.md)
- [FANTASIA Full](https://github.com/CBBIO/FANTASIA)

## Citation

Please cite this repository, the version-specific Zenodo lookup record used in
the analysis, and the relevant FANTASIA publications:

- [Performance of protein language models in model organisms](https://doi.org/10.1093/nargab/lqae078)
- [Application of FANTASIA to functional annotation of dark proteomes](https://doi.org/10.1038/s42003-025-08651-2)
- [Protocol explaining FANTASIA](https://doi.org/10.1007/978-1-0716-4623-6_8)
