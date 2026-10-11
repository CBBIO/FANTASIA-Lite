# Local Apptainer validation

Tested on 2026-10-11 on branch `feat/apptainer-offline-hpc`. Application files
under `src/` were unchanged.

## Environment and build

| Component | Observed version |
| --- | --- |
| Host | Ubuntu 26.04, x86_64 |
| Apptainer | 1.5.4 |
| GPU | NVIDIA GeForce RTX 5090, 32 GB |
| NVIDIA driver | 595.91.07 |
| Container Python | 3.12.3 |
| PyTorch | 2.7.1+cu128 |
| PyTorch CUDA | 12.8 |
| GPU compute capability | 12.0 |

The SIF was built successfully from `containers/fantasia-lite.def`:

```bash
mkdir -p dist/build-tmp dist/apptainer-cache
APPTAINER_TMPDIR="$PWD/dist/build-tmp" \
APPTAINER_CACHEDIR="$PWD/dist/apptainer-cache" \
apptainer build --fakeroot dist/fantasia-lite-cuda128.sif \
  containers/fantasia-lite.def > dist/build-cuda128.log 2>&1
```

The image is approximately 7.4 GiB. ESM's unused torchvision and torchtext
dependencies are deliberately omitted; the build dependency check allows only
those two omissions. The five embedding runtimes were exercised below.

## GPU and offline execution

From the repository root, the GPU smoke check passed:

```bash
apptainer exec --cleanenv --nv dist/fantasia-lite-cuda128.sif \
  python containers/smoke_gpu.py
```

The check imported the runtime classes, verified CUDA availability, and compared
a GPU matrix multiplication against its CPU result. PyTorch reported `sm_120`
among its compiled architectures.

For a stronger offline check, the local site configuration selected a test-only
Apptainer wrapper:

```sh
#!/bin/sh
action="$1"
shift
exec /usr/bin/apptainer "$action" --fakeroot --net --network none "$@"
```

Inside this container, `socket.if_nameindex()` reported only the loopback
interface, while PyTorch still detected the RTX 5090. This wrapper requires
site permission for network namespaces and is not needed on compute nodes
whose network policy already blocks external access.

Cached assets passed the launcher's `validate` command. The test used existing
local lookup data and the five pinned model snapshots, mounted read-only.
Container destinations were customized to `/validation/reference`,
`/validation/weights`, `/validation/input`, `/validation/output`, and
`/validation/scratch`.

```bash
scripts/fantasia-hpc --site-config dist/rtx5090-site.yaml validate
scripts/fantasia-hpc --site-config dist/rtx5090-site.yaml run \
  --input "$PWD/dist/smoke-input/query.fa" \
  --run-name rtx5090-all-models -- \
  --device cuda --use-gpu-lookup --serial-models \
  --embed-models 'esm2 prost_t5 prot_t5 ankh3 esmc' --embed-batch-size 1
```

The input contained one 17-residue sequence, `MKWVTFISLLFLFSSAYS`. All five
models generated embeddings and completed annotation lookup on CUDA. The job
exited successfully with one processed sequence and no skipped sequences.
`results.csv` contains all five model labels and 17 annotation rows.

| Measured stage | Seconds |
| --- | ---: |
| Embedding | 28.091 |
| Lookup | 17.468 |
| Postprocessing | 0.088 |
| Whole pipeline duration | 49.638 |

Local artifacts, excluded from Git, are available at:

- `dist/fantasia-lite-cuda128.sif`
- `dist/build-cuda128.log`
- `dist/rtx5090-all-models.log`
- `dist/smoke-output/rtx5090-all-models/results.csv`
- `dist/smoke-output/rtx5090-all-models/run_metadata.yaml`

## Automated checks and limits

The eight HPC tests and four container runtime adapter tests passed on the
host and subsequently inside the existing SIF with external networking disabled
while validating the GitHub Actions workflow. They cover site YAML, path mapping, bind modes, offline environment,
temporary lookup archive preparation/checksums, and dependency validation.

```bash
python3 -m unittest discover -s tests -p 'test_hpc_*.py'
python3 -m unittest discover -s tests -p 'test_container_runtime.py'
```

This is a functional smoke test, not a throughput or annotation-quality
benchmark. A100 hardware and Slurm submission were not available for this
test. The remote bundle/model download workflow was not rerun against a clean
cache; the GPU job used existing assets. The complete application test suite
was not run inside the SIF. Fully reproducible rebuilds still require transitive
dependency locks, immutable base-image references, and controlled OS packages.

## Future application changes

These improvements would require changes under `src/` and are intentionally
left for a separate change:

- Add a preinstalled-runtime mode that skips environment installation, removing
  the need for the container's validating `pip` adapter.
- Add first-class offline and cache-directory options to every model loader,
  with early errors naming missing files and the preparation command.
- Add common lookup/output/scratch directory options while retaining per-file
  overrides and clear CLI precedence.
- Record complete provenance for all five models and the application revision.
- Validate all required model files before embedding, with shared checksum
  verification and actionable errors for partial caches.
