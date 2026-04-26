
# FANTASIA Lite V1

**FANTASIA Lite V1** is a streamlined, standalone version of the full [FANTASIA pipeline](https://github.com/CBBIO/FANTASIA), designed for fast and efficient Gene Ontology (GO) annotation of protein sequences from local FASTA files using embedding comparisons.

FANTASIA Lite generates deep learning embeddings and perform nearest-neighbor annotation transfer, while intentionally removing the service stack used by the full project. The bundled lookup table covers multiple reference embedding spaces, while the built-in Lite embedder focuses on the models currently exposed directly through the local CLI.

The simplest way to run it is:

```bash
./scripts/minimal_pipeline.sh your_proteins.fa
```

Before any annotation run can work, the required lookup bundle must be downloaded from [Zenodo](https://zenodo.org/records/19742926) and placed in `data/lookup/`.

Unlike the full FANTASIA pipeline, **FANTASIA Lite does not require PostgreSQL, RabbitMQ, or a database-backed orchestration layer**. It runs locally from flat files:
- `lookup_table.npz` for reference embeddings
- `annotations.json` for GO annotations
- `accessions.json` for accession mapping

## Warning

The tradeoff is that Lite is simpler to deploy but can be slower if embedding and lookup are not tuned well. This repository now includes GPU-aware embedding and lookup controls intended to recover as much performance as possible without reintroducing external services.

## Scope and Purpose

The main purpose of Lite V1 is to provide a fast local annotator that can be dropped into other pipelines with minimal setup. The default fast path is intentionally simple: ProtT5 embeddings, cosine lookup, and `k=1` transfer. More advanced configuration is still available when you need multi-model runs, layered embedding export, larger `k`, or other custom settings.

This repository is ideal for users who want:
- Lightweight, local annotation of protein FASTA files
- No external database dependencies
- No PostgreSQL or RabbitMQ setup
- Simple setup and automated environment management
- High-quality functional annotation using experimental evidence from UniProt

For advanced features, large-scale annotation, or integration with external databases, see the full [FANTASIA repository](https://github.com/CBBIO/FANTASIA).

## What's New In V1

Compared with the earlier Lite V0 branch, Lite V1 includes several practical improvements:

- Faster end-to-end execution through batched embedding inference and a merged one-pass lookup flow
- In-process lookup execution instead of launching a separate lookup subprocess
- Reuse mode for rerunning lookup from an existing embeddings archive with `--reuse-embeddings`
- Optional TopGO export, now disabled by default unless `--topgo` is requested
- Clearer GPU tuning controls for embedding and lookup batch sizes
- Full-length embeddings by default, with truncation disabled unless explicitly requested
- Better reporting of skipped sequences and end-of-run coverage summaries
- Embedding-only workflows via `src/generate_embeddings.py`, separate from lookup
- Optional export of selected layers or all available layers from the Lite embedder
- Helper scripts for a default last-layer end-to-end run and an embedding-only layered export run
- Updated documentation that explains the separation between embedding support and lookup-bundle coverage

## Installation Requirements
- Python 3.10 or newer (the pipeline automatically creates and manages virtual environments)
- Required lookup bundle (`lookup_table.npz`, `annotations.json`, `accessions.json`) from [Zenodo: 19742926](https://zenodo.org/records/19742926) placed in `data/lookup/` before running the pipeline
- Internet connection for automatic dependency installation
- Sufficient disk space for outputs and embeddings (approximately 1-2 GB per run)
- Git (for cloning the repository)
- `wget` or `curl` (for downloading the lookup bundle)

## Lookup Table Details


The FANTASIA Lite V1 lookup table is built from the **UniProt November 2025 release** and includes only proteins with experimental evidence, ensuring high-quality functional annotations. All data was generated using PIS v3.1.0, the internal system used to extract and preprocess UniProt, PDB, and GOA data.

**Lookup bundle Zenodo DOI:** [10.5281/zenodo.19742926](https://doi.org/10.5281/zenodo.19742926)

Use this DOI to cite the lookup table or to access the official download page.

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
The lookup bundle (`fantasia_lite_data_folder.tgz`) contains three essential files:

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
- The default fast Lite V1 path is last-layer lookup with cosine distance, ProtT5 embeddings, and `k=1`. Layered embedding export is separate from lookup.
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

FANTASIA Lite V1 is designed first as a fast local annotator that is easy to integrate into other pipelines, while still supporting more configurable research-style runs when needed.

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

Beyond that, Lite V1 provides two main Python entrypoints:

### 1. Standard Pipeline (`fantasia_pipeline.py`)
For processing protein sequences and obtaining GO annotations:

Lite V1 keeps embedding and lookup in full precision (`float32`) for comparability with the full FANTASIA workflow. Mixed precision is not enabled by default.

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
fantasia_pipeline --serial-models --embed-models "prot_t5 ankh3" your_proteins.fa

# Advanced configuration
fantasia_pipeline \
    --embed-models prot_t5 \
    --limit-per-entry 5 \
    --results-csv my_results.csv \
    your_proteins.fa
```

### Helper Scripts

Three ready-to-run helper scripts are included under `scripts/`:

These helper scripts use paths relative to the repository root. Run them from inside the cloned `FANTASIA-Lite` directory. If needed, you can override the Python interpreter with `PYTHON_BIN=/path/to/python`.

- `scripts/minimal_pipeline.sh`
  This is the main Lite V1 entrypoint for fast annotation and pipeline integration. It runs the standard end-to-end workflow with the intended Lite defaults:
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
- `--embed-models`: Choose models (`prot_t5`, `ankh3`) - default: `prot_t5`
- `--serial-models`: Process requested models one after another instead of as one combined model-group. For a single model, this usually makes no practical difference. For more than one model, it is the recommended and safer setting because it reduces GPU/CPU memory pressure and avoids loading multiple embedders at the same time.
- `--limit-per-entry N`: Return top N annotations per sequence (default: 1)
- `--raw-results-csv PATH`: Optional raw neighbor-level output before GO-term consolidation. If omitted, Lite writes only `results.csv` unless `--limit-per-entry > 1`, in which case a raw file is created automatically as `k.<N>.results.csv`
- `--topgo`: Optional. Generate TopGO files after lookup. TopGO export is disabled by default
- `--distance-metric {cosine,euclidean}`: Lookup metric. Lite V1 defaults to `cosine`, but `euclidean` is also supported
- `--use-gpu-lookup`: Force GPU nearest-neighbor lookup when CUDA is available
- `--lookup-batch-size N`: Number of query embeddings compared per lookup batch
- `--sequence-queue-package N`: Outer packaging size before embedding forward passes
- `--embed-batch-size N`: Default forward-pass batch size for embeddings
- `--model-batch-sizes MODEL=N ...`: Per-model embedding batch size overrides
- `--length-filter N`: Optional truncation before embedding; `0` disables truncation and is the default

## Getting Started

FANTASIA Lite V1 is installed by cloning the repository, placing the lookup bundle in `data/lookup/`, and running a simple setup check from inside the cloned folder.

### Step 1: Clone the Repository
```bash
git clone https://github.com/CBBIO/FANTASIA-Lite.git
cd FANTASIA-Lite
```

To use THIS LATEST update of Lite V1 branch specifically:

```bash
git checkout FANTASIA-Lite-V1
```

### Step 2: Add the Lookup Bundle
Download the Lite lookup bundle from Zenodo and place these files in `data/lookup/`:

- `lookup_table.npz`
- `annotations.json`
- `accessions.json`

### Step 3: Run a Setup Check
The recommended validation step for Lite V1 is:

```bash
./scripts/minimal_pipeline.sh fasta_test/test.fa
```

This setup check:

1. creates or reuses the local virtual environment
2. installs the required dependencies automatically
3. runs the current Lite V1 default path on the small bundled test FASTA
4. writes results to a timestamped `outputs_YYYYMMDD_HHMMSS/` directory
5. confirms that embedding, lookup, and result writing are working correctly

**Expected output**: If everything works correctly, you'll see:
- Virtual environment creation (first run)
- Dependency installation progress
- Embedding generation progress bar
- Annotation results written to the output directory
- Success message

**Time estimate**: 
- First run: a few minutes to about 10-20 minutes, depending on dependency installation and model download state
- Subsequent runs: 1-2 minutes (only processes test file)

### Step 4: Process Your Own Data
The easiest way to run Lite V1 in its intended default mode is:

```bash
./scripts/minimal_pipeline.sh your_proteins.fa
```

This minimal runner is the recommended entrypoint for fast local annotation and for integration into larger workflows. It uses:
- automatic device detection
- `prot_t5`
- cosine lookup
- `k=1`
- full precision
- no sequence truncation unless you explicitly request it elsewhere

You can still call the full pipeline directly when you want more control:

```bash
fantasia_pipeline --serial-models --embed-models prot_t5 your_proteins.fa
```

## Device Selection

FANTASIA Lite automatically detects whether CUDA is available:
- If an NVIDIA GPU is detected, the pipeline defaults to `cuda`
- Otherwise it falls back to `cpu`

This detection is performed automatically by `fantasia_pipeline.py`, so most users do not need to set `--device` manually.

You can still override it explicitly:

```bash
# Force GPU
fantasia_pipeline --device cuda your_proteins.fa

# Force CPU
fantasia_pipeline --device cpu your_proteins.fa
```

## Performance Tuning

FANTASIA Lite is intentionally simpler than the full FANTASIA stack, but that simplicity means performance depends heavily on the embedding and lookup settings.

### Important Distinction: What Each Batch Size Means

- `--sequence-queue-package`
  Outer orchestration size. This controls how many input sequences are grouped into one work package before embedding. It does **not** change search completeness.

- `--embed-batch-size`
  Default PLM forward-pass batch size. This controls how many sequences are embedded together on the GPU or CPU.

- `--model-batch-sizes prot_t5=4 ankh3=4`
  Per-model overrides for embedding forward-pass batch size. These matter only for models you actually run. The best value depends on GPU memory and sequence lengths: if you see CUDA OOM skips in `failed_sequences.csv`, lower the relevant model batch size; if you see no skips and have spare memory, you can try increasing it.

- `--lookup-batch-size`
  Number of **query embeddings** processed together during nearest-neighbor search. Each lookup batch is still compared against the **full reference lookup table**. It does **not** search only the first `N` references and then stop.

### Full-Length Embeddings vs Truncation

By default, sequence truncation is **disabled**:

```bash
--length-filter 0
```

That means FANTASIA Lite will try to embed the full sequence. If the GPU can handle it, the full sequence is used. If you want to cap sequence length for safety or speed, set a positive value such as `--length-filter 2000`.

### Recommended Starting Profiles

For the fastest general GPU run without truncation:

```bash
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
```

This is the recommended starting profile for FANTASIA Lite on 24 GB-class GPUs when you want full-length embeddings and the fastest practical proteome-scale runtime. On a 20,223-sequence PRUB1 proteome, this profile completed end-to-end in about 24 minutes with GPU lookup enabled and processed 99.72% of proteins successfully.

At the end of each run, the pipeline prints a processed/skipped summary such as:

```text
Sequence summary: 20167/20223 processed (99.72%), 56/20223 skipped (0.28%).
```

Skipped proteins are reported in `failed_sequences.csv` and are typically extreme long-sequence CUDA OOM cases rather than ordinary proteins.

For smaller GPUs:

```bash
fantasia_pipeline \
    --device cuda \
    --use-gpu-lookup \
    --serial-models \
    --embed-models prot_t5 \
    --sequence-queue-package 50 \
    --embed-batch-size 4 \
    --model-batch-sizes prot_t5=4 \
    --lookup-batch-size 256 \
    your_proteins.fa
```

For CPU-only systems:

```bash
fantasia_pipeline \
    --device cpu \
    --serial-models \
    --embed-models prot_t5 \
    your_proteins.fa
```

### 2. Performance Analysis (`pipeline_timing_analyzer.py`)
For benchmarking, performance testing, and systematic analysis:

```bash
# Basic timing analysis - processes all test files with the default model
python3 src/pipeline_timing_analyzer.py

# Quick test with specific file and single model
python3 src/pipeline_timing_analyzer.py \
    --files fasta_test/test.fa \
    --model prot_t5 \
    --report-csv quick_benchmark.csv

# Compare models on specific files
python3 src/pipeline_timing_analyzer.py \
    --files fasta_test/test.fa fasta_test/UP000001940_6239.fasta \
    --model prot_t5 ankh3 \
    --report-csv model_comparison.csv

# Custom analysis with all options
python3 src/pipeline_timing_analyzer.py \
    --fasta-dir fasta_test \
    --model ankh3 \
    --files fasta_test/test.fa \
    --report-csv gpu_benchmark.csv
```

**Key Options:**
- `--files`: Specify individual FASTA files to process
- `--model`: Choose specific model(s) to test (default: `prot_t5`)
- `--report-csv`: Output file for timing results (default: `pipeline_timing_analysis.csv`)
- `--fasta-dir`: Directory containing FASTA files (default: `fasta_test`)

## Repository Structure
```
FANTASIA-Lite/
├── README.md                                # This documentation
├── LICENSE                                  # License information
├── .gitignore                               # Git ignore rules
├── data/
│   └── lookup/                              # Lookup database (download from Zenodo)
│       ├── accessions.json                  # Protein accession mappings
│       ├── annotations.json                 # GO annotation data
│       └── lookup_table.npz                 # Pre-computed embeddings database
├── fasta_test/                              # Test FASTA files for validation and benchmarking
│   ├── test.fa                              # Small test file (33 sequences)
│   ├── test_failure.fa                      # Test file with problematic sequences
│   ├── PRUB1_longiso.pep                    # Paratomella rubra proteome (non-model worm, not represented in standard databases)
│   ├── UP000001940_6239.fasta               # C. elegans proteome sample
│   └── MUSM_10090.fasta                     # Mouse proteome sample used for Lite V1 benchmarking
├── fantasia_pipeline.py                     # Main annotation pipeline
├── fantasia_no_db.py                        # Core lookup engine
├── generate_embeddings.py                   # Embedding generation module
└── pipeline_timing_analyzer.py              # Performance analysis and benchmarking tool
```

### Test Files (`fasta_test/`)
The repository includes comprehensive test files for validation and benchmarking:
- **`test.fa`**: Small test file with 33 valid sequences for quick validation
- **`test_failure.fa`**: Contains problematic sequences to test error handling
- **`PRUB1_longiso.pep`**: Proteome of the non-model worm *Paratomella rubra*, useful as a realistic dark-proteome style test case outside standard database-centric examples
- **`UP000001940_6239.fasta`**: Complete C. elegans proteome for realistic testing
- **`MUSM_10090.fasta`**: Mouse proteome sample used in the revalidated Lite V1 benchmark set

## Outputs

### Standard Pipeline Outputs (`outputs_YYYYMMDD_HHMMSS/`)
Each pipeline run creates a timestamped directory containing:
- **`results.csv`**: Main GO annotation results after GO-term consolidation and best-row selection
- **`raw_results.csv`** or **`k.<N>.results.csv`**: Optional raw neighbor-level results before GO-term consolidation
- **`query_embeddings.npz`**: Generated embeddings for input sequences  
- **`failed_sequences.csv`**: Sequences that failed processing with error details
- **`fantasia_config.yaml`**: Configuration used for the run
- **`run_metadata.yaml`**: Timestamped run metadata for traceability, including resolved parameters and output paths
- **`run.log`**: Timestamped pipeline log capturing console output from the run
- **`topgo/`**: Optional TopGO-compatible files for downstream analysis when `--topgo` is enabled
  - `<model>.topgo.<F|P|C>.txt`: GO terms by functional category (Function/Process/Component)

### Output File Types and How to Read Them

- **`results.csv`**: CSV table. This is the main file most users want. Each row is a final GO annotation kept after consolidation, so it is the best place to inspect the final biological interpretation of the run.
- **`raw_results.csv`** or **`k.<N>.results.csv`**: CSV table. This is the pre-consolidation lookup output. Use it when you want to inspect neighbor structure, `hit_rank`, or the effect of `k > 1`.
- **`query_embeddings.npz`**: NumPy archive. This stores the generated query embeddings and is mainly useful for reuse, benchmarking, or embedding-only downstream work rather than manual inspection.
- **`failed_sequences.csv`**: CSV table. This lists the sequences that were skipped or failed, together with the error reason. In full-length GPU runs, these are often extreme long-sequence CUDA OOM cases.
- **`fantasia_config.yaml`**: YAML text file. This captures the effective run configuration and is the main provenance file for reproducing a run.
- **`run_metadata.yaml`**: YAML text file. This records traceability metadata such as timestamps, resolved parameters, output paths, summary counts, and stage timings.
- **`run.log`**: Plain-text log. This is the chronological console transcript of the run and is the best file to inspect when debugging behavior or checking progress details after the fact.
- **`topgo/`**: Directory of plain-text TopGO helper files. These are only created when `--topgo` is enabled and are meant for downstream enrichment workflows, not for primary lookup interpretation.

**`results.csv` vs `raw_results.csv`**

- **`raw_results.csv`** preserves the direct lookup output. It keeps the selected nearest-neighbor hit structure, including `hit_rank`, before GO-term consolidation. If one reference hit carries multiple GO terms, or if `--limit-per-entry` is greater than `1`, this file can contain multiple rows for the same query and hit.
- **`results.csv`** is the cleaned final annotation table. It consolidates GO terms and keeps the best supporting row for each `(query, model, GO term, category)` combination.
- If you only need final annotations, `results.csv` is usually enough. If you want to inspect the underlying hit structure or retain multiple neighbors before consolidation, keep `raw_results.csv`.

**Column meanings**

- `query_accession`: input sequence identifier from the FASTA file
- `hit_rank`: nearest-neighbor rank in the raw output only (`1` = best hit, `2` = second hit, etc.)
- `reference_id`: internal lookup-bundle identifier for the matched reference entry
- `model_key`: embedding space used for the match, such as `Prot-T5` or `Ankh3-Large`
- `distance`: embedding-space distance between query and matched reference; lower is better
- `reliability_index`: normalized score derived from distance and clipped to `[0, 1]`; higher is better
- `distance_metric`: lookup metric used for the match, usually `cosine` in Lite V1
- `uniprot_accession`: UniProt accession of the matched reference protein
- `go_id`: transferred GO identifier
- `category`: GO namespace (`F` = Molecular Function, `P` = Biological Process, `C` = Cellular Component)
- `go_description`: human-readable GO term description
- `evidence_codes`: evidence code or merged evidence codes supporting that GO term in the reference protein

### Timing Analyzer Outputs
- **`pipeline_timing_analysis.csv`**: Comprehensive performance metrics including:
  - GPU model and memory specifications
  - Runtime and processing rates
  - GPU/CPU usage information  
  - Sequence processing statistics
  - Successfully processed vs failed sequences
  - GPU memory usage monitoring
  - Model comparison data
  - Timestamped pipeline output directory references

**Note**: Requires nvidia-smi for GPU monitoring (optional for CPU-only systems).

## Performance Analysis Features

The `pipeline_timing_analyzer.py` tool provides comprehensive benchmarking capabilities:

- **Hardware Comparison**: Compare GPU vs CPU performance across different systems
- **Model Evaluation**: Systematic comparison between prot_t5 and ankh3 models  
- **Scalability Testing**: Analyze performance across different file sizes and sequence counts
- **Regression Testing**: Track performance changes across pipeline versions
- **Resource Monitoring**: GPU memory usage and processing rate analysis

## Performance Benchmarks

The Lite V1 pipeline changed substantially, so older benchmark matrices from earlier Lite revisions are no longer directly representative. The benchmark section below only reports runs that have been revalidated after the current V1 optimization work.

### Currently Revalidated

The following end-to-end benchmarks reflect the current recommended fast profile:

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

### Common Issues
- **CUDA compatibility**: Set `TORCH_INDEX` environment variable for specific CUDA versions
- **Memory errors**: Use `--serial-models` and process one model at a time
- **Missing dependencies**: The pipeline automatically installs required packages
- **Lookup bundle missing**: Download from Zenodo and extract to `data/lookup/`
- **Out-of-memory errors**: Reduce `--embed-models` to a single model and keep `--serial-models` enabled

### Performance Optimization
- **GPU memory**: Use `--serial-models` to prevent multiple models loading simultaneously
- **Processing speed**: Start with `prot_t5` model for fastest results
- **Large files**: Increase `--chunk-size` for better memory management

### Important Notes
- Gzipped FASTA files are decompressed on the fly; no manual prep is required
- Sequences longer than the model limit are skipped and logged; other sequences continue
- Each pipeline run creates a timestamped directory (`outputs_YYYYMMDD_HHMMSS`)
- Parallel model execution is technically possible but rarely worth the memory cost

## FAQ
- **Can I use `.gz` FASTA files?** Yes. Compression is handled automatically.
- **What if a sequence is too long?** It is recorded in `outputs/failed_sequences.csv`; the rest of the batch continues.
- **Does the lookup bundle include ESM3c?** Yes. The current Lite lookup bundle includes `ESM3c`, but the built-in Lite embedder CLI is still focused on `prot_t5` and `ankh3`.


## Acknowledgements

FANTASIA Lite V1 is derived from the full [FANTASIA pipeline](https://github.com/CBBIO/FANTASIA) and incorporates methods from [GOPredSim](https://github.com/Rostlab/goPredSim). Transformer models are provided via [Hugging Face](https://huggingface.co/).

**Key Publications:**
- [Performance of protein language models in model organisms](https://doi.org/10.1093/nargab/lqae078)
- [Application of FANTASIA to functional annotation of dark proteomes](https://doi.org/10.1038/s42003-025-08651-2)
- [Protocol explaining FANTASIA](https://doi.org/10.1007/978-1-0716-4623-6_8)

**Citing FANTASIA**

If you use FANTASIA in your research, please cite the following publications:

- Martínez-Redondo, G. I., Barrios, I., Vázquez-Valls, M., Rojas, A. M., & Fernández, R. (2024).Illuminating the functional landscape of the dark proteome across the Animal Tree of Life.
    DOI: 10.1101/2024.02.28.582465

- Barrios-Núñez, I., Martínez-Redondo, G. I., Medina-Burgos, P., Cases, I., Fernández, R., & Rojas, A. M. (2024). Decoding proteome functional information in model organisms using protein language models.
    DOI: 10.1101/2024.02.14.580341


**Main Developers:**
- Ana M. Rojas: a.rojas.m@csic.es
- Àlex Domínguez Rodríguez: adomrod4@upo.es

**Project Team:**
- Ana M. Rojas: a.rojas.m@csic.es
- Rosa Fernández: rosa.fernandez@ibe.upf-csic.es
- Aureliano Bombarely: abombarely@ibmcp.upv.es
- Ildefonso Cases: icasdia@upo.es
- Àlex Domínguez Rodríguez: adomrod4@upo.es
- Gemma I. Martínez-Redondo: gemma.martinez@ibe.upf-csic.es
- Belén Carbonetto: belen.carbonetto.metazomics@gmail.com
- Iñigo de Martín: imartinagirre@gmail.com
- Sofía García Juan
---

**Version**: FANTASIA Lite V1  
**Last Updated**: November 2025
**Funded by** EOSC-OSCARS Fun4Biodiversity
