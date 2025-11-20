# FANTASIA Lite V0

FANTASIA Lite V0 is a streamlined version of the FANTASIA pipeline designed for efficient GO annotation transfer on local FASTA files. This lightweight Python implementation provides protein functional annotation through embedding-based nearest-neighbor lookup without requiring the full FANTASIA infrastructure.

## What You Need
- Python 3.10 or newer (the pipeline automatically creates and manages virtual environments)
- Lookup bundle (`lookup_table.npz`, `annotations.json`, `accessions.json`) from [Zenodo: 17393454](https://zenodo.org/records/17393454) placed in `data/lookup/`
- Internet connection for automatic dependency installation
- Sufficient disk space for outputs and embeddings (approximately 1-2 GB per run)
- Git (for cloning the repository)
- `wget` or `curl` (for downloading the lookup bundle)

## Lookup Table Details

The FANTASIA Lite V0 lookup table is built from the **UniProt November 2025 release** and includes only proteins with experimental evidence, ensuring high-quality functional annotations. All data was generated using PIS v3.1.0, the internal system used to extract and preprocess UniProt, PDB, and GOA data.

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
The lookup bundle (`fantasia_lookup_bundle.tgz`) contains three essential files:

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
   - Fetches `fantasia_lookup_bundle.tgz` containing pre-computed embeddings
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
- `--embed-models`: Choose models (`prot_t5`, `ankh3`, `esm2`) - default: `prot_t5`
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
FantasiaLiteV0/
├── README.md                          # This documentation
├── data/
│   └── lookup/                        # Lookup database (download from Zenodo)
│       ├── accessions.json           # Protein accession mappings
│       ├── annotations.json          # GO annotation data
│       └── lookup_table.npz          # Pre-computed embeddings database
├── fasta_test/                       # Test FASTA files for validation and benchmarking
│   ├── UP000001940_6239.fasta       # C. elegans proteome sample
│   ├── UP000002311_2025_10_19.fasta # Bacterial proteome sample
│   ├── test.fa                      # Small test file (33 sequences)
│   └── test_failure.fa              # Test file with problematic sequences
├── fantasia_pipeline.py             # Main annotation pipeline
├── fantasia_no_db.py                # Core lookup engine
├── generate_embeddings.py           # Embedding generation module
└── pipeline_timing_analyzer.py      # Performance analysis and benchmarking tool
```

### Test Files (`fasta_test/`)
The repository includes comprehensive test files for validation and benchmarking:
- **`test.fa`**: Small test file with 33 valid sequences for quick validation
- **`test_failure.fa`**: Contains problematic sequences to test error handling
- **`UP000001940_6239.fasta`**: Complete C. elegans proteome for realistic testing
- **`UP000002311_2025_10_19.fasta`**: Bacterial proteome sample for diversity testing

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
- **Model Evaluation**: Systematic comparison between prot_t5, ankh3, and esm2 models  
- **Scalability Testing**: Analyze performance across different file sizes and sequence counts
- **Regression Testing**: Track performance changes across pipeline versions
- **Resource Monitoring**: GPU memory usage and processing rate analysis

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
- **`esm2`**: ESM-2 evolutionary scale modeling (fast, good for large-scale analysis)

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

**Contact:**
- Ana M. Rojas Mendoza — a.rojas.m@csic.es
- Rosa M. Fernandez — rosa.fernandez@ibe.upf-csic.es  
- Àlex Dominguez Rodriguez — adomrod4@upo.es

---

**Version**: FANTASIA Lite V0  
**Last Updated**: November 2025
