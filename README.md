
# FANTASIA Lite V0

**FANTASIA Lite V0** is a streamlined, standalone version of the full [FANTASIA pipeline](https://github.com/CBBIO/FANTASIA), designed for fast and efficient Gene Ontology (GO) annotation of protein sequences from local FASTA files.

FANTASIA Lite leverages state-of-the-art protein language models (ProtT5, Ankh3) to generate deep learning embeddings and perform nearest-neighbor annotation transfer—without requiring a database server or the full FANTASIA infrastructure.

This repository is ideal for users who want:
- Lightweight, local annotation of protein FASTA files
- No external database dependencies
- Simple setup and automated environment management
- High-quality functional annotation using experimental evidence from UniProt

For advanced features, large-scale annotation, or integration with external databases, see the full [FANTASIA repository](https://github.com/CBBIO/FANTASIA).

## What You Need
- Python 3.10 or newer (the pipeline automatically creates and manages virtual environments)
- Lookup bundle (`lookup_table.npz`, `annotations.json`, `accessions.json`) from [Zenodo: 17720428](https://zenodo.org/records/17720428) placed in `data/lookup/`
- Internet connection for automatic dependency installation
- Sufficient disk space for outputs and embeddings (approximately 1-2 GB per run)
- Git (for cloning the repository)
- `wget` or `curl` (for downloading the lookup bundle)

## Lookup Table Details


The FANTASIA Lite V0 lookup table is built from the **UniProt November 2025 release** and includes only proteins with experimental evidence, ensuring high-quality functional annotations. All data was generated using PIS v3.1.0, the internal system used to extract and preprocess UniProt, PDB, and GOA data.

**Lookup bundle Zenodo DOI:** [10.5281/zenodo.17720428](https://doi.org/10.5281/zenodo.17720428)

Use this DOI to cite the lookup table or to access the official download page.

### Core Statistics
- **Total accessions**: 127,546 UniProt entries
- **Total proteins**: 127,546 protein records
- **Total unique sequences**: 124,397 (3,149 proteins share identical sequences due to isoforms and redundancy)
- **Total GO annotations**: 627,932 experimental annotations
- **Sequence redundancy**: 2.47%

### Embedding Coverage
- **ProtT5-XL-UniRef50**: 124,397 embeddings (100% coverage)
- **Ankh3-Large**: 124,397 embeddings (100% coverage)

### Package Contents
The lookup bundle (`fantasia_lite_data_folder.tgz`) contains three essential files:

1. **`lookup_table.npz`**
   - Precomputed protein embeddings (ProtT5 and Ankh3)
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

## Getting Started

### Step 1: Clone the Repository
```bash
git clone https://github.com/CBBIO/FANTASIA-Lite.git
cd FANTASIA-Lite
```

### Step 2: Run Setup and Test Script
The easiest way to get started is to use the automated setup script:

```bash
./setup_and_test.sh
```

**What this script does:**

1. **Creates directory structure** (`data/lookup/`)
   - Sets up the necessary folders for lookup files

2. **Downloads lookup bundle** from Zenodo (~2-3 GB)
   - Fetches `fantasia_lite_data_folder.tgz` containing pre-computed embeddings
   - Skips download if the file already exists
   - Uses `wget` or `curl` depending on what's available on your system

3. **Extracts lookup files**
   - Automatically detects the archive structure
   - Places files in the correct location (`data/lookup/`)
   - Verifies that all three required files are present:
     - `lookup_table.npz` - Pre-computed protein embeddings database
     - `annotations.json` - GO annotation data
     - `accessions.json` - Protein accession mappings
   - Skips extraction if files already exist

4. **Runs validation test**
   - Executes: `python3 fantasia_pipeline.py --serial-models --embed-models prot_t5 fasta_test/test.fa`
   - Creates a Python virtual environment automatically (first run only)
   - Installs all required dependencies (PyTorch, transformers, etc.)
   - Processes 33 test sequences using the ProtT5 model
   - Generates output in a timestamped directory (`outputs_YYYYMMDD_HHMMSS/`)
   - Confirms the pipeline is working correctly

**Expected output**: If everything works correctly, you'll see:
- Download progress and extraction confirmation
- Virtual environment creation (first run)
- Dependency installation progress
- Embedding generation progress bar
- Annotation results written to the output directory
- Success message

**Time estimate**: 
- First run: 10-20 minutes (includes downloads and dependency installation)
- Subsequent runs: 1-2 minutes (only processes test file)

### Step 3: Process Your Own Data
Once the test completes successfully, you can annotate your own protein sequences:

```bash
python3 fantasia_pipeline.py --serial-models --embed-models prot_t5 your_proteins.fa
```

## Quick Start

FANTASIA Lite V0 provides two main ways to run protein annotation:

### 1. Standard Pipeline (`fantasia_pipeline.py`)
For processing protein sequences and obtaining GO annotations:

```bash
# Basic usage - single model annotation
python3 fantasia_pipeline.py your_proteins.fa

# Recommended usage with specific model
python3 fantasia_pipeline.py --serial-models --embed-models prot_t5 your_proteins.fa

# Multiple models (slower but more comprehensive)
python3 fantasia_pipeline.py --serial-models --embed-models "prot_t5 ankh3" your_proteins.fa

# Advanced configuration
python3 fantasia_pipeline.py \
    --embed-models prot_t5 \
    --limit-per-entry 5 \
    --results-csv my_results.csv \
    your_proteins.fa
```

**Key Options:**
- `--embed-models`: Choose models (`prot_t5`, `ankh3`) - default: `prot_t5`
- `--serial-models`: Process models sequentially (recommended for memory efficiency)
- `--limit-per-entry N`: Return top N annotations per sequence (default: 1)

### 2. Performance Analysis (`pipeline_timing_analyzer.py`)
For benchmarking, performance testing, and systematic analysis:

```bash
# Basic timing analysis - processes all test files with both models [prot_t5, ankh3]
python3 pipeline_timing_analyzer.py

# Quick test with specific file and single model
python3 pipeline_timing_analyzer.py \
    --files fasta_test/test.fa \
    --model prot_t5 \
    --report-csv quick_benchmark.csv

# Compare models on specific files
python3 pipeline_timing_analyzer.py \
    --files fasta_test/test.fa fasta_test/UP000001940_6239.fasta \
    --report-csv model_comparison.csv

# Custom analysis with all options
python3 pipeline_timing_analyzer.py \
    --fasta-dir fasta_test \
    --model ankh3 \
    --files fasta_test/test.fa \
    --report-csv gpu_benchmark.csv
```

**Key Options:**
- `--files`: Specify individual FASTA files to process
- `--model`: Choose specific model(s) to test (default: both `prot_t5` and `ankh3`)
- `--report-csv`: Output file for timing results (default: `pipeline_timing_analysis.csv`)
- `--fasta-dir`: Directory containing FASTA files (default: `fasta_test`)

## Repository Structure
```
FANTASIA-Lite/
├── README.md                                # This documentation
├── LICENSE                                  # License information
├── setup_and_test.sh                        # Automated setup and validation script
├── .gitignore                               # Git ignore rules
├── data/
│   └── lookup/                              # Lookup database (download from Zenodo)
│       ├── accessions.json                  # Protein accession mappings
│       ├── annotations.json                 # GO annotation data
│       └── lookup_table.npz                 # Pre-computed embeddings database
├── fasta_test/                              # Test FASTA files for validation and benchmarking
│   ├── test.fa                              # Small test file (33 sequences)
│   ├── test_failure.fa                      # Test file with problematic sequences
│   ├── UP000001940_6239.fasta               # C. elegans proteome sample
│   ├── UP000002311_2025_10_19.fasta         # Bacterial proteome sample
│   └── GCF_000001735.4_TAIR10.1_protein.faa # A. thaliana proteome sample
├── fantasia_pipeline.py                     # Main annotation pipeline
├── fantasia_no_db.py                        # Core lookup engine
├── generate_embeddings.py                   # Embedding generation module
└── pipeline_timing_analyzer.py              # Performance analysis and benchmarking tool
```

### Test Files (`fasta_test/`)
The repository includes comprehensive test files for validation and benchmarking:
- **`test.fa`**: Small test file with 33 valid sequences for quick validation
- **`test_failure.fa`**: Contains problematic sequences to test error handling
- **`UP000001940_6239.fasta`**: Complete C. elegans proteome for realistic testing
- **`UP000002311_2025_10_19.fasta`**: Bacterial proteome sample for diversity testing
- **`GCF_000001735.4_TAIR10.1_protein.faa`**: Arabidopsis thaliana proteome for plant protein testing

## Outputs

### Standard Pipeline Outputs (`outputs_YYYYMMDD_HHMMSS/`)
Each pipeline run creates a timestamped directory containing:
- **`results.csv`**: Main GO annotation results
- **`query_embeddings.npz`**: Generated embeddings for input sequences  
- **`failed_sequences.csv`**: Sequences that failed processing with error details
- **`fantasia_config.yaml`**: Configuration used for the run
- **`topgo/`**: TopGO-compatible files for downstream analysis
  - `<model>.topgo.<F|P|C>.txt`: GO terms by functional category (Function/Process/Component)

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

The following benchmarks provide realistic performance expectations across different GPU configurations:

### NVIDIA RTX A2000 12GB (12 GB VRAM)

| Dataset | Sequences | Model | Runtime | Rate (seq/s) | Annotations | Failed |
|---------|-----------|-------|---------|--------------|-------------|---------|
| **A. thaliana proteome** | 48,265 | ProtT5 | 4h 2m | 3.32 | 155,720 | 3,929 |
| (GCF_000001735.4) | | Ankh3 | 4h 52m | 2.76 | 177,643 | 35 |
| **C. elegans proteome** | 19,831 | ProtT5 | 1h 30m | 3.69 | 65,543 | 1,411 |
| (UP000001940_6239) | | Ankh3 | 1h 53m | 2.93 | 75,543 | 45 |
| **Bacterial proteome** | 6,067 | ProtT5 | 30m 34s | 3.31 | 27,380 | 705 |
| (UP000002311) | | Ankh3 | 41m 42s | 2.43 | 32,475 | 7 |
| **Small test file** | 33 | ProtT5 | 28s | 1.18 | 132 | 1 |
| (test.fa) | | Ankh3 | 18s | 1.84 | 142 | 0 |
| **Problematic sequences** | 176 | ProtT5 | 29s | 6.08 | 34 | 172 |
| (test_failure.fa) | | Ankh3 | 20s | 8.83 | 41 | 171 |

### Tesla T4 (15 GB VRAM)

| Dataset | Sequences | Model | Runtime | Rate (seq/s) | Annotations | Failed |
|---------|-----------|-------|---------|--------------|-------------|---------|
| **A. thaliana proteome** | 48,265 | ProtT5 | 6h 44m | 1.99 | 173,586 | 55 |
| (GCF_000001735.4) | | Ankh3 | 7h 17m | 1.84 | 177,830 | 3 |
| **C. elegans proteome** | 19,831 | ProtT5 | 2h 27m | 2.25 | 73,264 | 62 |
| (UP000001940_6239) | | Ankh3 | 2h 30m | 2.21 | 75,776 | 12 |
| **Bacterial proteome** | 6,067 | ProtT5 | 54m 8s | 1.87 | 32,187 | 8 |
| (UP000002311) | | Ankh3 | 54m 53s | 1.84 | 32,540 | 1 |
| **Small test file** | 33 | ProtT5 | 56s | 0.59 | 139 | 0 |
| (test.fa) | | Ankh3 | 45s | 0.77 | 142 | 0 |
| **Problematic sequences** | 176 | ProtT5 | 58s | 3.05 | 41 | 171 |
| (test_failure.fa) | | Ankh3 | 2m 33s | 1.23 | 44 | 170 |

### NVIDIA GeForce RTX 3090 Ti (24 GB VRAM)

| Dataset | Sequences | Model | Runtime | Rate (seq/s) | Annotations | Failed |
|---------|-----------|-------|---------|--------------|-------------|---------|
| **A. thaliana proteome** | 48,265 | ProtT5 | 1h 23m | 9.71 | 173,940 | 0 |
| (GCF_000001735.4) | | Ankh3 | 1h 28m | 9.15 | 177,848 | 0 |
| **C. elegans proteome** | 19,831 | ProtT5 | 33m 24s | 9.90 | 73,593 | 11 |
| (UP000001940_6239) | | Ankh3 | 35m 16s | 9.37 | 75,803 | 6 |
| **Bacterial proteome** | 6,067 | ProtT5 | 11m 32s | 8.77 | 32,265 | 0 |
| (UP000002311) | | Ankh3 | 12m 13s | 8.28 | 32,547 | 0 |
| **Small test file** | 33 | ProtT5 | 21s | 1.54 | 139 | 0 |
| (test.fa) | | Ankh3 | 12s | 2.70 | 142 | 0 |
| **Problematic sequences** | 176 | ProtT5 | 33s | 5.33 | 41 | 171 |
| (test_failure.fa) | | Ankh3 | 6m 6s | 0.48 | 550 | 71 |

### NVIDIA A100-SXM4-40GB (40 GB VRAM)

| Dataset | Sequences | Model | Runtime | Rate (seq/s) | Annotations | Failed |
|---------|-----------|-------|---------|--------------|-------------|---------|
| **A. thaliana proteome** | 48,265 | ProtT5 | 1h 50m | 7.30 | 173,940 | 0 |
| (GCF_000001735.4) | | Ankh3 | 1h 57m | 6.85 | 177,848 | 0 |
| **C. elegans proteome** | 19,831 | ProtT5 | 1h 41m | 3.27 | 73,636 | 3 |
| (UP000001940_6239) | | Ankh3 | 50m | 6.58 | 75,837 | 2 |
| **Bacterial proteome** | 6,067 | ProtT5 | 16m 33s | 6.11 | 32,265 | 0 |
| (UP000002311) | | Ankh3 | 15m 38s | 6.47 | 32,547 | 0 |
| **Small test file** | 33 | ProtT5 | 1m 38s | 0.34 | 139 | 0 |
| (test.fa) | | Ankh3 | 1m 3s | 0.52 | 142 | 0 |
| **Problematic sequences** | 176 | ProtT5 | 8m 31s | 0.34 | 360 | 65 |
| (test_failure.fa) | | Ankh3 | 11m 38s | 0.25 | 932 | 31 |

*All benchmarks performed with `--serial-models` flag (November 2025). The RTX 3090 Ti shows significantly faster processing rates, particularly beneficial for large-scale proteome annotation. The RTX A2000 12GB shows lower performance and higher failure rates due to memory constraints, but remains functional for most workloads. Tesla T4 values of the **Small test file** and **Problematic sequences** are averaged across 2 runs per configuration.*

## Advanced Usage

### Environment Management
The pipeline automatically manages Python virtual environments:
```bash
# Virtual environment is created automatically in venv/
# To clean up and force rebuild:
rm -rf venv/
python3 fantasia_pipeline.py your_file.fa  # Will recreate venv automatically
```

### Batch Processing
For processing multiple files systematically:
```bash
# Process multiple specific files
python3 pipeline_timing_analyzer.py \
    --files file1.fa file2.fa file3.fa \
    --model prot_t5

# Process all files in a directory
python3 pipeline_timing_analyzer.py \
    --fasta-dir my_proteomes/ \
    --report-csv batch_results.csv
```

### Memory Optimization
For large files or limited memory systems:
```bash
# Use serial processing and single model
python3 fantasia_pipeline.py \
    --serial-models \
    --embed-models prot_t5 \
    --chunk-size 200 \
    large_proteome.fa
```

## Supported Models

FANTASIA Lite V0 supports three embedding models:
- **`prot_t5`**: Protein T5 model (recommended, good balance of speed and accuracy)
- **`ankh3`**: ANKH large protein language model (slower but potentially more accurate)

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
- **Is ESM3c supported?** Not by default. Mirror the checkpoint locally and extend `MODEL_REGISTRY` if you need it.


## Acknowledgements

FANTASIA Lite V0 is derived from the full [FANTASIA pipeline](https://github.com/CBBIO/FANTASIA) and incorporates methods from [GOPredSim](https://github.com/Rostlab/goPredSim). Transformer models are provided via [Hugging Face](https://huggingface.co/).

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

**Version**: FANTASIA Lite V0  
**Last Updated**: November 2025
**Funded by** EOSC-OSCARS Fun4Biodiversity
