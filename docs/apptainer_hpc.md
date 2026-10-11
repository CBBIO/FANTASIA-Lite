# Apptainer on HPC systems with offline compute nodes

This integration leaves the programs under `src/` unchanged. The image contains
a preinstalled Python environment. Because the current pipeline retains its
dependency installer, the image replaces its `pip` executable with a validation
adapter. It checks that requested packages and versions are already installed,
and fails if they are missing: jobs do not install packages or access indexes.

## Build the image

### Build and publish with GitHub Actions

The [Apptainer workflow](../.github/workflows/apptainer.yml) builds and tests a
SIF when a branch push changes the workflow, container, application, scripts,
tests, or `pyproject.toml`. Tags matching `apptainer-*` also trigger a build.
It publishes the tested image to GHCR using the repository's `GITHUB_TOKEN`
with `packages: write`; no additional publishing secret is required. The
organization must allow GitHub Actions to create packages. An existing GHCR
package must grant this repository Actions access.

After committing and pushing the workflow, open **Actions → Apptainer image**
in GitHub. The job summary provides the exact pull command and SIF SHA-256.
The package name is the lowercase repository name, and the tag is
`sha-<full-commit-sha>`; for this repository:

```bash
# On node0, replace FULL_COMMIT_SHA with the successful workflow's commit.
apptainer pull fantasia-lite-cuda128.sif \
  oras://ghcr.io/cbbio/fantasia-lite:sha-FULL_COMMIT_SHA
```

Successful publications from `feat/apptainer-offline-hpc` or the repository's
default branch also update `latest` to the same SIF. Other branches and tag
builds publish only their commit tag. To download the most recently published
image from either of those branches:

```bash
apptainer pull fantasia-lite-cuda128.sif \
  oras://ghcr.io/cbbio/fantasia-lite:latest
```

`latest` moves whenever an eligible build publishes successfully; it does not
identify a fixed release. Use the commit tag and recorded SIF checksum for
experiments that need a specific image.

Set `runtime.image` in the site YAML to the downloaded SIF's absolute path.
Download the `apptainer-build-...` artifact from that workflow run for the
checksum, dependency inventory, and build/test logs. Extract the checksum file
next to the image and verify it:

```bash
sha256sum -c fantasia-lite-cuda128.sif.sha256
```

Model weights and lookup data remain outside the SIF. Run `prepare` on node0
and `validate` before submitting jobs as described below.

GHCR packages are initially private, even for a public source repository.
For anonymous pulls, a package administrator can change the package visibility
to public. For private packages, authenticate interactively on node0 using a
GitHub personal access token (classic) with `read:packages` and access to the
package, entering the token at the password prompt:

```bash
apptainer registry login --username YOUR_GITHUB_USERNAME oras://ghcr.io
```

See the official [GHCR authentication and package access documentation](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry)
and [Apptainer ORAS support](https://apptainer.org/docs/user/latest/docker_and_oci.html).

The default runner is `ubuntu-24.04`. The job removes unrelated SDKs from this
disposable hosted runner and checks for at least 40 GiB free before building.
This is a starting capacity check, not a guarantee that future dependency
versions will fit. If disk space is insufficient, set the repository Actions
variable `APPTAINER_RUNNER` to an Ubuntu x86_64 larger-runner label with more
storage. A custom runner needs passwordless sudo and permission to install
Apptainer and create network namespaces. SDK cleanup is skipped on custom
runners. Do not use a shared production HPC login node as a build runner.

The job builds without a GPU, tests runtime imports and the launcher tests in
a container with external networking disabled, and publishes only after they
pass. It does not download model weights or perform GPU inference. Concurrent
pushes to the same branch cancel the older build. The workflow also supports
manual dispatch with `publish: false` for build/test-only runs once it is
available on the repository's default branch. Build-only runs retain the small
logs/metadata artifact, not the SIF.

Commit tags identify the source, but a rebuild of the same commit can overwrite
its tag: the dependency/base-image limitations below still apply. Save the SIF
checksum from the successful run when recording an experiment.

### Build locally

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

`storage.*.host`, `runtime.image`, and `runtime.apptainer_cache` expand exported
environment variables written as `$NAME` or `${NAME}`, followed by home-directory
expansion (`~`). For example:

```bash
export FANTASIA_DATA_ROOT="$HOME/fantasia"
```

```yaml
runtime:
  image: "${FANTASIA_DATA_ROOT}/fantasia-lite-cuda128.sif"
  apptainer_cache: "$HOME/.cache/apptainer"
  use_nv: true
storage:
  models:
    host: "${FANTASIA_DATA_ROOT}/models"
    container: /srv/fantasia/models
    mode: ro
  # Include the other storage entries from the example configuration.
```

Expansion happens on the host each time the launcher reads the YAML, before
Apptainer starts. Export the variables on node0 and in the Slurm job environment.
An undefined variable produces an error naming the field and variable; expanded
paths must still be absolute. Variables are substituted once, without shell
execution: `${NAME:-default}` and `$(command)` are not supported. Use `$$` for
a literal dollar sign. Container destinations are literal absolute paths and
do not undergo variable or home-directory expansion. This feature is in the
host launcher, so updating the repository is enough to use it with an existing
SIF.

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
