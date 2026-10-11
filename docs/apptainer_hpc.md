# Apptainer on HPC systems with offline compute nodes

This integration leaves the programs under `src/` unchanged. The image contains
a preinstalled Python environment. Because the current pipeline retains its
dependency installer, the image replaces its `pip` executable with a validation
adapter. It checks that requested packages and versions are already installed,
and fails if they are missing: jobs do not install packages or access indexes.

## Build the image

Build on a system with Internet access. From the repository root:

```bash
mkdir -p dist
apptainer build --fakeroot dist/fantasia-lite-cuda128.sif containers/fantasia-lite.def
```

The image uses CUDA 12.8, Python 3.12, and PyTorch 2.7.1 with Blackwell support
(including RTX 5090), as well as Ampere/Hopper support. The built-in
test checks imports of PyTorch, Transformers, ESM-C, h5py, and PyYAML. A compute
node needs a compatible NVIDIA driver, and the container must be invoked with
`--nv`. The launcher adds this flag automatically when `runtime.use_nv: true`.

The RTX 5090 requires the newer Blackwell runtime; this changes the original
CUDA 12.4 target to CUDA 12.8. See the [local validation report](apptainer_validation.md)
for the tested hardware, commands, and remaining limitations.

Direct Python dependencies are pinned. Transitive dependencies are recorded in
`/opt/fantasia/containers/installed.txt` inside the SIF, but are not fully locked
as build inputs. The base image tag and Ubuntu package repositories can also
change. Reuse the tested SIF for consistent jobs; rebuilding is not yet
bit-for-bit reproducible.

## Configure the site

Copy [`config/hpc-site.example.yaml`](../config/hpc-site.example.yaml) to a
private site directory and replace every `/shared/...` value. Each resource
declares both its host filesystem path and its destination inside the
container. These destinations can differ entirely from `/lookup`, `/models`,
or `/output`.

The YAML also specifies the image, the Apptainer cache, the URL of a specific
Zenodo bundle release, and the archive's SHA-256. Do not rely solely on the
concept DOI: record a specific release URL and obtain or calculate its
checksum before publishing the configuration.

The launcher requires Python 3 with PyYAML on the host where it is invoked.
The host Python environment does not need PyTorch or Transformers.

## Prepare assets on node0

Run preparation exclusively on a login node with Internet access:

```bash
scripts/fantasia-hpc --site-config /private/path/hpc-site.yaml prepare
scripts/fantasia-hpc --site-config /private/path/hpc-site.yaml validate
```

`prepare` downloads the bundle only if any of the three lookup files is
missing. It downloads all five Hugging Face snapshots using the revisions
pinned by FANTASIA Lite, verifies the downloaded bundle and ESM-C weight
checksums, and writes `fantasia-assets-manifest.json` to `storage.staging.host`.
`validate` uses `local_files_only`; it downloads nothing and provides a useful
check before submitting a job.

To prepare a subset of models, for example:

```bash
scripts/fantasia-hpc --site-config /private/path/hpc-site.yaml prepare \
  --models prot_t5 esm2
```

## Submit jobs

The FASTA file must be under `storage.inputs.host`. The launcher creates a
subdirectory for each run under the output and scratch directories, mounts
lookup data, models, and inputs as read-only, and sets `HF_HUB_OFFLINE=1`,
`TRANSFORMERS_OFFLINE=1`, and `FANTASIA_OFFLINE=1`.

The Hugging Face variables disable model downloads. `FANTASIA_OFFLINE` is an
environment marker; the unchanged application does not consume it. The launcher
does not create a network namespace: network isolation is a cluster policy or
an additional Apptainer setting. The local validation used `--net --network none`
to verify execution without external connectivity. Apptainer's default host
mounts also remain active; the explicit data binds do not constitute filesystem
isolation.

```bash
scripts/fantasia-hpc --site-config /private/path/hpc-site.yaml run \
  --input /shared/fantasia/inputs/proteome.faa \
  --run-name proteome-001 \
  -- --serial-models --embed-models prot_t5
```

See [`scripts/fantasia-hpc.slurm.example`](../scripts/fantasia-hpc.slurm.example)
for a Slurm job example. The launcher's `--dry-run` mode prints the Apptainer
command and bind mounts without executing the container.

For example, from the repository root on node0:

```bash
export FANTASIA_ROOT="$PWD"
export SITE_CONFIG=/private/path/hpc-site.yaml
export INPUT_FASTA=/shared/fantasia/inputs/proteome.faa
sbatch scripts/fantasia-hpc.slurm.example
```

Adjust the partition and GPU request in the example for your cluster. The
explicit repository path is needed because Slurm executes a spooled copy of
the job script. Arguments after `--` are passed to the pipeline after the
YAML-derived defaults, so an explicit pipeline option overrides its default.
Paths in these extra options must use container paths.

A SIF built by Apptainer can also run with a compatible Singularity
installation. Set `runtime.apptainer: singularity` in the YAML if that is the
command provided by your site.
