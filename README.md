# FANTASIA Lite V1

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
python3 src/generate_embeddings.py \
  --fasta data/my_proteome.faa.gz \
  --output my_embeddings.npz --models prot_t5 --device cuda
```

## Requirements

- Linux or macOS with Python 3.9+; Python 3.10+ recommended
- NVIDIA GPU strongly recommended; the tuned profile was tested with 24 GB VRAM
- at least 16 GB system RAM as a practical starting point
- approximately 3.1 GB for the extracted lookup bundle
- approximately 12 GB free disk for a first ProtT5 proteome run, including
  environment, model cache, reference data, and outputs
- internet access for dependency and model installation

Smaller GPUs may work with reduced model batch and residue limits. See the
[detailed usage reference](docs/usage_reference.md) for every default, output
column, resource note, model/reference distinction, tuning advice, and
benchmark.

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
