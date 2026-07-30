#!/usr/bin/env python3
"""
fantasia_pipeline.py
--------------------

Refactored orchestration for the FANTASIA lookup workflow.  This module
replaces the legacy shell pipeline with a Python-only implementation.

The pipeline performs the following high-level steps:
  1. Create (or reuse) a Python virtual environment and install dependencies.
  2. Validate that the lookup flat files (``lookup_table.npz``, ``annotations.json``,
     ``accessions.json``) are present.
  3. Chunk the input FASTA sequences to manageable batches.
  4. Generate embeddings for each chunk and run the lookup script.
  5. Merge intermediate CSV and NPZ artifacts and write a final config YAML.

Each stage is encapsulated in a dedicated function to keep the module
modular and reusable.
"""

from __future__ import annotations

import argparse
import dataclasses
import csv
import json
import os
import platform
import re
import shutil
import socket
from site import venv
import subprocess
import sys
import time
from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, TextIO

FANTASIA_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_VENV_DIR = FANTASIA_ROOT / "venv"
DEFAULT_EMBED_MODELS = "prot_t5"
SUPPORTED_EMBED_MODELS = {"esm2", "prost_t5", "prot_t5", "ankh3", "esmc"}
MODEL_ALIASES = {
    "esm": "esm2",
    "esm2": "esm2",
    "esm-2": "esm2",
    "prost-t5": "prost_t5",
    "prost_t5": "prost_t5",
    "prostt5": "prost_t5",
    "prot-t5": "prot_t5",
    "prot_t5": "prot_t5",
    "prott5": "prot_t5",
    "ankh3": "ankh3",
    "ankh3-large": "ankh3",
    "ankh3_large": "ankh3",
    "esmc": "esmc",
    "esm-c": "esmc",
    "esm_c": "esmc",
    "esm3c": "esmc",
}
ESMC_PYTHON_REQUIREMENT = ">=3.12,<3.13"
DEFAULT_CHUNK_SIZE = 500
DEFAULT_LIMIT_PER_ENTRY = 1
DEFAULT_DISTANCE_METRIC = "cosine"

DEFAULT_LOOKUP_DIR = FANTASIA_ROOT / "data" / "lookup"
# Generate timestamped output directory
_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
DEFAULT_OUTPUT_DIR = Path(f"outputs_{_timestamp}")
DEFAULT_TEMP_DIR = DEFAULT_OUTPUT_DIR / "tmp"
DEFAULT_FAILURE_REPORT = DEFAULT_OUTPUT_DIR / "failed_sequences.csv"
DEFAULT_CHUNK_FAILURE_DIR = DEFAULT_TEMP_DIR / "failures"
DEFAULT_TOPGO_DIR = DEFAULT_OUTPUT_DIR / "topgo"
DEFAULT_RUN_LOG = DEFAULT_OUTPUT_DIR / "run.log"
DEFAULT_RUN_METADATA = DEFAULT_OUTPUT_DIR / "run_metadata.yaml"

_ACTIVE_LOG_HANDLE: Optional[TextIO] = None


class PipelineError(RuntimeError):
    """Raised when a pipeline stage fails in a recoverable way."""


def parse_embed_model_names(models: str) -> List[str]:
    requested: List[str] = []
    unknown: List[str] = []
    for raw_model in models.strip().split():
        normalized = MODEL_ALIASES.get(raw_model.strip().lower())
        if normalized is None:
            unknown.append(raw_model)
            continue
        if normalized not in requested:
            requested.append(normalized)
    if unknown:
        supported = ", ".join(sorted(SUPPORTED_EMBED_MODELS))
        raise PipelineError(
            f"Unsupported embedding model(s): {', '.join(unknown)}. "
            f"Supported models: {supported}."
        )
    return requested


def validate_python_for_models(models: Sequence[str]) -> None:
    if "esmc" not in models:
        return
    if (3, 12) <= sys.version_info[:2] < (3, 13):
        return
    raise PipelineError(
        "ESM-C uses esm==3.2.3, which currently supports Python "
        f"{ESMC_PYTHON_REQUIREMENT}. Current interpreter is Python "
        f"{platform.python_version()} at {sys.executable}. Run the pipeline with "
        "Python 3.12 and recreate the virtual environment."
    )


@dataclasses.dataclass
class PipelineConfig:
    fasta_path: Path
    venv_dir: Path
    lookup_npz: Path
    annotations_json: Path
    accessions_json: Path
    embeddings_npz: Path
    config_yaml: Path
    results_csv: Path
    raw_results_csv: Optional[Path]
    topgo_dir: Path
    chunk_dir: Path
    chunk_embed_dir: Path
    chunk_results_dir: Path
    chunk_config_dir: Path
    chunk_size: int
    embed_models: str
    serial_models: bool
    device: Optional[str]
    use_gpu_lookup: Optional[bool]
    limit_per_entry: int
    distance_metric: str
    lookup_batch_size: int
    queue_batch_size: int
    max_sequence_length: int
    embed_batch_size: int
    max_batch_residues: int
    model_batch_sizes: str
    torch_index: Optional[str]
    chunk_failure_dir: Path
    failure_report: Path
    generate_topgo: bool
    reuse_embeddings: bool
    run_log: Path
    run_metadata_yaml: Path

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "PipelineConfig":
        fasta_path = Path(args.fasta_path).resolve()
        if not fasta_path.exists():
            raise PipelineError(f"FASTA file not found: {fasta_path}")

        env = os.environ
        embed_models = args.embed_models or env.get("EMBED_MODELS", DEFAULT_EMBED_MODELS)
        serial_models = args.serial_models or env.get("SERIAL_MODELS") == "1"

        device = args.device or env.get("PYTORCH_DEVICE")
        env_gpu_lookup = env.get("FANTASIA_GPU_LOOKUP")
        use_gpu_lookup: Optional[bool]
        if args.use_gpu_lookup:
            use_gpu_lookup = True
        elif args.no_gpu_lookup:
            use_gpu_lookup = False
        elif env_gpu_lookup is not None:
            use_gpu_lookup = env_gpu_lookup.strip().lower() in {"1", "true", "yes", "on"}
        else:
            use_gpu_lookup = None
        torch_index = args.torch_index or env.get("TORCH_INDEX")

        results_csv = Path(args.results_csv).resolve()
        raw_results_csv: Optional[Path]
        if args.raw_results_csv:
            raw_results_csv = Path(args.raw_results_csv).resolve()
        elif args.limit_per_entry > 1:
            raw_filename = f"k.{args.limit_per_entry}.results.csv"
            raw_results_csv = results_csv.parent / raw_filename
        else:
            raw_results_csv = None

        return cls(
            fasta_path=fasta_path,
            venv_dir=Path(args.venv_dir).resolve(),
            lookup_npz=Path(args.lookup_npz).resolve(),
            annotations_json=Path(args.annotations_json).resolve(),
            accessions_json=Path(args.accessions_json).resolve(),
            embeddings_npz=Path(args.embeddings_npz).resolve(),
            config_yaml=Path(args.config_yaml).resolve(),
            results_csv=results_csv,
            raw_results_csv=raw_results_csv,
            topgo_dir=Path(args.topgo_dir).resolve(),
            chunk_dir=Path(args.chunk_dir).resolve(),
            chunk_embed_dir=Path(args.chunk_embed_dir).resolve(),
            chunk_results_dir=Path(args.chunk_results_dir).resolve(),
            chunk_config_dir=Path(args.chunk_config_dir).resolve(),
            chunk_size=args.chunk_size,
            embed_models=embed_models,
            serial_models=serial_models,
            device=device,
            use_gpu_lookup=use_gpu_lookup,
            limit_per_entry=args.limit_per_entry,
            distance_metric=args.distance_metric,
            lookup_batch_size=args.lookup_batch_size,
            queue_batch_size=args.queue_batch_size,
            max_sequence_length=args.max_sequence_length,
            embed_batch_size=args.embed_batch_size,
            max_batch_residues=args.max_batch_residues,
            model_batch_sizes=args.model_batch_sizes,
            torch_index=torch_index,
            chunk_failure_dir=Path(args.chunk_failure_dir).resolve(),
            failure_report=Path(args.failure_report).resolve(),
            generate_topgo=args.topgo,
            reuse_embeddings=args.reuse_embeddings,
            run_log=Path(args.run_log).resolve(),
            run_metadata_yaml=Path(args.run_metadata_yaml).resolve(),
        )


@dataclasses.dataclass
class ChunkJob:
    chunk_path: Path
    chunk_embeddings: Path
    chunk_config: Path
    chunk_results: Path
    chunk_failures: Path
    chunk_raw_results: Optional[Path]


@dataclasses.dataclass
class PipelineRunStats:
    total_sequences: int
    skipped_sequences: int
    embedding_time_seconds: float = 0.0
    lookup_time_seconds: float = 0.0
    postprocess_time_seconds: float = 0.0

    @property
    def processed_sequences(self) -> int:
        return max(0, self.total_sequences - self.skipped_sequences)

    @property
    def processed_fraction(self) -> float:
        if self.total_sequences == 0:
            return 0.0
        return self.processed_sequences / self.total_sequences

    @property
    def skipped_fraction(self) -> float:
        if self.total_sequences == 0:
            return 0.0
        return self.skipped_sequences / self.total_sequences

    @property
    def measured_stage_time_seconds(self) -> float:
        return (
            self.embedding_time_seconds
            + self.lookup_time_seconds
            + self.postprocess_time_seconds
        )


class TeeStream:
    def __init__(self, primary: TextIO, secondary: TextIO):
        self.primary = primary
        self.secondary = secondary

    def write(self, data: str) -> int:
        self.primary.write(data)
        self.secondary.write(data)
        return len(data)

    def flush(self) -> None:
        self.primary.flush()
        self.secondary.flush()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    # Resolve the root directory of FANTASIA-Lite


    parser = argparse.ArgumentParser(
        description="Run the FANTASIA pipeline purely in Python."
    )
    parser.add_argument("fasta_path", help="Path to the FASTA file to annotate.")
    parser.add_argument(
        "--venv-dir",
        default=str(DEFAULT_VENV_DIR),
        help="Directory for the Python virtual environment (default: %(default)s).",
    )
    parser.add_argument(
        "--lookup-npz",
        default=str(DEFAULT_LOOKUP_DIR / "lookup_table.npz"),
        help="Output NPZ filename for the lookup table.",
    )
    parser.add_argument(
        "--annotations-json",
        default=str(DEFAULT_LOOKUP_DIR / "annotations.json"),
        help="Output JSON filename for annotations.",
    )
    parser.add_argument(
        "--accessions-json",
        default=str(DEFAULT_LOOKUP_DIR / "accessions.json"),
        help="Output JSON filename for accession mappings.",
    )
    parser.add_argument(
        "--embeddings-npz",
        default=str(DEFAULT_OUTPUT_DIR / "query_embeddings.npz"),
        help="Output NPZ filename for merged query embeddings.",
    )
    parser.add_argument(
        "--config-yaml",
        default=str(DEFAULT_OUTPUT_DIR / "fantasia_config.yaml"),
        help="Output YAML filename for the consolidated configuration.",
    )
    parser.add_argument(
        "--results-csv",
        default=str(DEFAULT_OUTPUT_DIR / "results.csv"),
        help="Output CSV filename for merged lookup results.",
    )
    parser.add_argument(
        "--raw-results-csv",
        default=None,
        help=(
            "Optional CSV filename for raw neighbour results before GO consolidation "
            "(defaults to k.<limit>.results.csv when --limit-per-entry > 1)."
        ),
    )
    parser.add_argument(
        "--topgo-dir",
        default=str(DEFAULT_TOPGO_DIR),
        help="Directory where TopGO-compatible tables are written.",
    )
    parser.add_argument(
        "--topgo",
        action="store_true",
        help="Generate TopGO-compatible outputs after the main lookup stage.",
    )
    parser.add_argument(
        "--chunk-dir",
        default=str(DEFAULT_TEMP_DIR / "fasta_chunks"),
        help="Directory used to store chunked FASTA files.",
    )
    parser.add_argument(
        "--chunk-embed-dir",
        default=str(DEFAULT_TEMP_DIR / "chunk_embeddings"),
        help="Directory used to store chunk embedding archives.",
    )
    parser.add_argument(
        "--chunk-results-dir",
        default=str(DEFAULT_TEMP_DIR / "chunk_results"),
        help="Directory used to store chunk-level lookup results.",
    )
    parser.add_argument(
        "--chunk-config-dir",
        default=str(DEFAULT_TEMP_DIR / "chunk_configs"),
        help="Directory used to store per-chunk lookup configurations.",
    )
    parser.add_argument(
        "--chunk-failure-dir",
        default=str(DEFAULT_CHUNK_FAILURE_DIR),
        help="Directory used to store per-chunk embedding failure reports.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help="Number of sequences per chunk (default: %(default)s).",
    )
    parser.add_argument(
        "--embed-models",
        default=None,
        help=(
            "Space-separated list of embedding models "
            "(choices: esm2, prost_t5, prot_t5, ankh3, esmc; default from env or 'prot_t5'). "
            "Full FANTASIA aliases such as ESM, Prost-T5, Prot-T5, Ankh3-Large, "
            "and ESM3c are accepted."
        ),
    )
    parser.add_argument(
        "--serial-models",
        action="store_true",
        help="Run each embedding model sequentially rather than as a group.",
    )
    parser.add_argument(
        "--device",
        default=None,
        help="Explicit device to use for PyTorch (overrides detection).",
    )
    parser.add_argument(
        "--use-gpu-lookup",
        action="store_true",
        help="Force the lookup stage to use the GPU lookup script.",
    )
    parser.add_argument(
        "--no-gpu-lookup",
        action="store_true",
        help="Force the lookup stage to use the CPU lookup script even on CUDA machines.",
    )
    parser.add_argument(
        "--limit-per-entry",
        type=int,
        default=DEFAULT_LIMIT_PER_ENTRY,
        help="Maximum number of annotations returned per entry.",
    )
    parser.add_argument(
        "--distance-metric",
        default=DEFAULT_DISTANCE_METRIC,
        help="Distance metric used during lookup (default: %(default)s).",
    )
    parser.add_argument(
        "--lookup-batch-size",
        type=int,
        default=256,
        help="Batch size used by the GPU lookup script (default: %(default)s).",
    )
    parser.add_argument(
        "--queue-batch-size",
        "--sequence-queue-package",
        type=int,
        default=100,
        help=(
            "Outer sequence package size before per-model forward batching "
            "(upstream embedding.queue_batch_size; default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--max-sequence-length",
        "--length-filter",
        type=int,
        default=0,
        help=(
            "Truncate sequences before embedding. "
            "0 disables truncation (upstream embedding.max_sequence_length)."
        ),
    )
    parser.add_argument(
        "--embed-batch-size",
        type=int,
        default=8,
        help="Maximum number of sequences embedded together per model (default: %(default)s).",
    )
    parser.add_argument(
        "--max-batch-residues",
        type=int,
        default=12000,
        help=(
            "Upper bound on padded residues per embedding batch. "
            "Larger values can be faster but use more memory (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--model-batch-sizes",
        nargs="*",
        default=(),
        help=(
            "Optional per-model forward batch overrides in the form "
            "'esm2=1 prost_t5=1 prot_t5=4 ankh3=4 esmc=1'."
        ),
    )
    parser.add_argument(
        "--torch-index",
        default=None,
        help="Custom package index for installing CUDA-enabled PyTorch wheels.",
    )
    parser.add_argument(
        "--failure-report",
        default=str(DEFAULT_FAILURE_REPORT),
        help="Path to a CSV summarising sequences skipped during embedding.",
    )
    parser.add_argument(
        "--run-log",
        default=str(DEFAULT_RUN_LOG),
        help="Path to a timestamped run log capturing pipeline output.",
    )
    parser.add_argument(
        "--run-metadata-yaml",
        default=str(DEFAULT_RUN_METADATA),
        help="Path to a YAML file recording run metadata and resolved parameters.",
    )
    parser.add_argument(
        "--reuse-embeddings",
        action="store_true",
        help=(
            "Skip embedding generation and reuse an existing --embeddings-npz archive "
            "for lookup-only reruns."
        ),
    )
    return parser.parse_args(argv)


def run_subprocess(
    cmd: Sequence[str],
    *,
    env: Optional[dict] = None,
    cwd: Optional[Path] = None,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """Run a subprocess command with optional environment overrides."""
    process_env = os.environ.copy()
    if env:
        process_env.update(env)
    stdout_chunks: List[str] = []
    process = subprocess.Popen(
        cmd,
        env=process_env,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert process.stdout is not None

    for line in process.stdout:
        sys.stdout.write(line)
        stdout_chunks.append(line)
        if _ACTIVE_LOG_HANDLE is not None:
            _ACTIVE_LOG_HANDLE.write(line)

    returncode = process.wait()
    result = subprocess.CompletedProcess(
        cmd,
        returncode,
        stdout="".join(stdout_chunks),
        stderr="",
    )
    if check and returncode != 0:
        raise subprocess.CalledProcessError(
            returncode,
            cmd,
            output=result.stdout,
            stderr=result.stderr,
        )
    return result


def create_virtualenv(venv_dir: Path) -> None:
    if venv_dir.exists():
        print(f"Reusing existing virtual environment at {venv_dir}")
        return
    print(f"Creating virtual environment at {venv_dir}")
    import venv

    builder = venv.EnvBuilder(with_pip=True, clear=False, upgrade=False)
    builder.create(str(venv_dir))


def venv_python_executable(venv_dir: Path) -> Path:
    scripts_dir = "Scripts" if platform.system() == "Windows" else "bin"
    python_path = venv_dir / scripts_dir / ("python.exe" if os.name == "nt" else "python")
    if not python_path.exists():
        raise PipelineError(f"Could not find python executable in venv: {python_path}")
    return python_path


def venv_pip_executable(venv_dir: Path) -> Path:
    scripts_dir = "Scripts" if platform.system() == "Windows" else "bin"
    pip_name = "pip.exe" if os.name == "nt" else "pip"
    pip_path = venv_dir / scripts_dir / pip_name
    if not pip_path.exists():
        raise PipelineError(f"Could not find pip executable in venv: {pip_path}")
    return pip_path


def count_fasta_sequences(fasta_path: Path) -> int:
    count = 0
    with fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def count_report_rows(report_path: Path) -> int:
    if not report_path.exists():
        return 0
    with report_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return sum(1 for row in reader if row)


def _safe_git_output(args: Sequence[str]) -> Optional[str]:
    try:
        completed = subprocess.run(
            list(args),
            cwd=str(FANTASIA_ROOT),
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return None
    text = completed.stdout.strip()
    return text or None


def build_run_metadata(
    config: PipelineConfig,
    *,
    argv: Sequence[str],
    started_at: datetime,
    finished_at: Optional[datetime] = None,
    device: Optional[str] = None,
    stats: Optional[PipelineRunStats] = None,
    status: str = "started",
) -> Dict[str, object]:
    metadata: Dict[str, object] = {
        "status": status,
        "started_at": started_at.isoformat(),
        "finished_at": finished_at.isoformat() if finished_at else None,
        "duration_seconds": (
            round((finished_at - started_at).total_seconds(), 3)
            if finished_at is not None
            else None
        ),
        "cwd": str(Path.cwd()),
        "hostname": socket.gethostname(),
        "python_executable": sys.executable,
        "python_version": sys.version.split()[0],
        "argv": list(argv),
        "git_branch": _safe_git_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        "git_commit": _safe_git_output(["git", "rev-parse", "HEAD"]),
        "resolved_device": device,
        "outputs": {
            "results_csv": str(config.results_csv),
            "raw_results_csv": str(config.raw_results_csv) if config.raw_results_csv else None,
            "embeddings_npz": str(config.embeddings_npz),
            "failure_report": str(config.failure_report),
            "config_yaml": str(config.config_yaml),
            "run_log": str(config.run_log),
            "run_metadata_yaml": str(config.run_metadata_yaml),
            "topgo_dir": str(config.topgo_dir),
        },
        "parameters": {},
    }
    parameters: Dict[str, object] = {}
    for field in dataclasses.fields(config):
        value = getattr(config, field.name)
        if isinstance(value, Path):
            parameters[field.name] = str(value)
        elif isinstance(value, (list, tuple)):
            parameters[field.name] = [str(item) for item in value]
        else:
            parameters[field.name] = value
    metadata["parameters"] = parameters
    if stats is not None:
        metadata["sequence_summary"] = {
            "total_sequences": stats.total_sequences,
            "processed_sequences": stats.processed_sequences,
            "processed_fraction": round(stats.processed_fraction, 6),
            "skipped_sequences": stats.skipped_sequences,
            "skipped_fraction": round(stats.skipped_fraction, 6),
        }
        metadata["stage_timing_seconds"] = {
            "embedding": round(stats.embedding_time_seconds, 3),
            "lookup": round(stats.lookup_time_seconds, 3),
            "postprocess": round(stats.postprocess_time_seconds, 3),
            "measured_total": round(stats.measured_stage_time_seconds, 3),
        }
    return metadata


def write_run_metadata(path: Path, metadata: Dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import yaml

        path.write_text(yaml.safe_dump(metadata, sort_keys=False), encoding="utf-8")
    except Exception:
        path.write_text(json.dumps(metadata, indent=2, default=str), encoding="utf-8")


def write_lookup_config(path: Path, conf: Dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import yaml

        path.write_text(yaml.safe_dump(conf, sort_keys=False), encoding="utf-8")
    except Exception:
        path.write_text(json.dumps(conf, indent=2, default=str), encoding="utf-8")


def run_lookup_subprocess(
    conf: Dict[str, object],
    *,
    use_gpu_lookup: bool,
    venv_python: Path,
    config_path: Path,
    device: str,
) -> None:
    write_lookup_config(config_path, conf)
    script_name = "fantasia_no_db_gpu.py" if use_gpu_lookup else "fantasia_no_db.py"
    run_subprocess(
        [str(venv_python), FANTASIA_ROOT / "src" / script_name, "--config", str(config_path)],
        env=torch_runtime_env(device),
    )


def detect_device(config_device: Optional[str]) -> str:
    if config_device:
        print(f"Using explicit PyTorch device setting: {config_device}")
        return config_device
    nvidia_smi = shutil.which("nvidia-smi")
    if nvidia_smi:
        print("Detected NVIDIA GPU via nvidia-smi; using CUDA.")
        return "cuda"
    print("No CUDA-capable GPU detected; defaulting to CPU.")
    return "cpu"


def validate_torch_device(venv_python: Path, device: str, *, explicit: bool) -> str:
    """Return a usable PyTorch device after a real allocation probe."""
    if not device.lower().startswith("cuda"):
        return device

    probe = (
        "import torch; "
        "assert torch.cuda.is_available(), 'torch.cuda.is_available() is false'; "
        "x = torch.empty((1,), device='cuda'); "
        "torch.cuda.synchronize(); "
        "print(torch.cuda.get_device_name(0))"
    )
    result = subprocess.run(
        [str(venv_python), "-c", probe],
        env={
            **os.environ,
            "PYTORCH_NO_CUDA_MEMORY_CACHING": os.environ.get(
                "PYTORCH_NO_CUDA_MEMORY_CACHING", "1"
            ),
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode == 0:
        print(f"Validated CUDA device with PyTorch: {result.stdout.strip()}")
        return device

    detail = (result.stderr or result.stdout).strip()
    message = (
        "PyTorch could not allocate on CUDA even though CUDA was requested/detected. "
        f"Probe error: {detail}"
    )
    if explicit:
        raise PipelineError(
            message
            + "\nUse --device cpu to continue without CUDA, or fix the CUDA/NVML driver "
            "environment before rerunning with --device cuda."
        )

    print(message)
    print("Falling back to CPU for embedding and lookup.")
    return "cpu"


def install_packages(config: PipelineConfig) -> str:
    device = detect_device(config.device)
    requested_models = parse_embed_model_names(
        str(getattr(config, "embed_models", DEFAULT_EMBED_MODELS))
    )
    install_esmc = "esmc" in requested_models
    validate_python_for_models(requested_models)
    pip_bin = venv_pip_executable(config.venv_dir)

    print("Upgrading pip…")
    run_subprocess([str(pip_bin), "install", "--upgrade", "pip"])

    torch_packages = ["torch"]
    if device.lower() == "cuda":
        if config.torch_index:
            print(f"Installing PyTorch with CUDA support from {config.torch_index}…")
            run_subprocess(
                [str(pip_bin), "install", "--index-url", config.torch_index, *torch_packages]
            )
        else:
            print("Installing PyTorch with CUDA support from the default package index…")
            run_subprocess([str(pip_bin), "install", *torch_packages])
    else:
        print("Installing CPU-only PyTorch packages…")
        run_subprocess([str(pip_bin), "install", *torch_packages])

    base_deps = [
        "numpy",
        "requests",
        "transformers<4.48.2",
        "biopython",
        "pandas",
        "scipy",
        "sentencepiece",
        "protobuf",
        "huggingface-hub",
        "httpx",
    ]
    esmc_runtime_deps = [
        "ipython",
        "einops",
        "biotite>=1.0.0",
        "msgpack-numpy",
        "scikit-learn",
        "brotli",
        "attrs",
        "cloudpathlib",
        "tenacity",
        "zstd",
        "ipywidgets",
        "py3dmol",
        "pydssp",
        "boto3",
        "pygtrie",
        "dna_features_viewer",
    ]

    system_name = platform.system()
    arch_name = platform.machine()
    if not (system_name == "Darwin" and arch_name == "arm64"):
        base_deps.append("h5py")
    else:
        print(
            "Skipping h5py installation on macOS arm64 to avoid HDF5 architecture issues."
        )

    print("Installing base dependencies…")
    run_subprocess([str(pip_bin), "install", *base_deps])
    if install_esmc:
        print("Installing ESM-C runtime dependencies without torchvision/torchtext…")
        run_subprocess([str(pip_bin), "install", *esmc_runtime_deps])
        # ESM 3.2.3 declares torchvision/torchtext, but the local ESM-C embedding
        # path does not use them. Installing the package without deps avoids
        # fragile torch/torchvision wheel mismatches.
        run_subprocess([str(pip_bin), "install", "--no-deps", "esm==3.2.3"])
    else:
        print("ESM-C model not requested; skipping ESM-C runtime dependencies.")

    return validate_torch_device(
        venv_python_executable(config.venv_dir),
        device,
        explicit=config.device is not None,
    )


def should_use_gpu_lookup(config: PipelineConfig, device: str) -> bool:
    """Decide whether the lookup stage should run via the GPU script."""
    if config.use_gpu_lookup is not None:
        if config.use_gpu_lookup and not device.lower().startswith("cuda"):
            raise PipelineError(
                "--use-gpu-lookup was requested, but the resolved PyTorch device is "
                f"{device!r}. Use --no-gpu-lookup or fix CUDA before forcing GPU lookup."
            )
        return config.use_gpu_lookup
    return device.lower().startswith("cuda")


def torch_runtime_env(device: str) -> dict[str, str]:
    """Environment overrides for PyTorch subprocesses."""
    if not device.lower().startswith("cuda"):
        return {}
    return {
        "PYTORCH_NO_CUDA_MEMORY_CACHING": os.environ.get(
            "PYTORCH_NO_CUDA_MEMORY_CACHING", "1"
        )
    }


def have_lookup_artifacts(config: PipelineConfig) -> bool:
    return (
        config.lookup_npz.exists()
        and config.annotations_json.exists()
        and config.accessions_json.exists()
    )


def ensure_lookup_artifacts(config: PipelineConfig) -> None:
    if have_lookup_artifacts(config):
        return
    missing = []
    for path in (config.lookup_npz, config.annotations_json, config.accessions_json):
        if not path.exists():
            missing.append(path)
    message = (
        "Lookup artifacts are missing. Provide the flat files from the pre-built lookup "
        "bundle (lookup_table.npz, annotations.json, accessions.json) before running the "
        "pipeline.\nMissing paths:\n"
        + "\n".join(f"  - {path}" for path in missing)
    )
    raise PipelineError(message)


def count_fasta_records(fasta_path: Path) -> int:
    opener = gzip_open if fasta_path.suffix in {".gz", ".gzip"} else Path.open
    count = 0
    with opener(fasta_path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def gzip_open(path: Path, mode: str):
    import gzip

    return gzip.open(path, mode)


def iter_fasta_records(fasta_path: Path) -> Iterable[tuple[str, str]]:
    opener = gzip_open if fasta_path.suffix in {".gz", ".gzip"} else Path.open
    with opener(fasta_path, "rt") as handle:
        header: Optional[str] = None
        sequence_lines: List[str] = []
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_lines)
                header = line
                sequence_lines = []
            else:
                sequence_lines.append(line)
        if header is not None:
            yield header, "".join(sequence_lines)


def chunk_fasta(config: PipelineConfig) -> List[Path]:
    total_records = count_fasta_records(config.fasta_path)
    if total_records == 0:
        return []
    if total_records <= config.chunk_size:
        return [config.fasta_path]

    if config.chunk_dir.exists():
        shutil.rmtree(config.chunk_dir)
    config.chunk_dir.mkdir(parents=True, exist_ok=True)

    chunk_paths: List[Path] = []
    chunk_index = 0
    records_in_chunk = 0
    chunk_file: Optional[TextIO] = None
    chunk_path: Optional[Path] = None

    def open_new_chunk() -> tuple[TextIO, Path]:
        nonlocal chunk_index
        chunk_index += 1
        stem = config.fasta_path.stem
        path = config.chunk_dir / f"{stem}_part{chunk_index:04d}.fa"
        file_handle = open(path, "w", encoding="utf-8")
        chunk_paths.append(path)
        return file_handle, path

    try:
        for header, sequence in iter_fasta_records(config.fasta_path):
            if records_in_chunk == 0:
                chunk_file, chunk_path = open_new_chunk()
            assert chunk_file is not None
            chunk_file.write(f"{header}\n")
            chunk_file.write("\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60)))
            chunk_file.write("\n")
            records_in_chunk += 1
            if records_in_chunk == config.chunk_size:
                chunk_file.close()
                chunk_file = None
                chunk_path = None
                records_in_chunk = 0
        if chunk_file is not None:
            chunk_file.close()
    finally:
        if chunk_file is not None and not chunk_file.closed:
            chunk_file.close()

    return chunk_paths


def ensure_output_directories(config: PipelineConfig) -> None:
    for directory in (
        config.chunk_dir,
        config.chunk_embed_dir,
        config.chunk_results_dir,
        config.chunk_config_dir,
        config.chunk_failure_dir,
    ):
        if directory.exists():
            shutil.rmtree(directory)
        directory.mkdir(parents=True, exist_ok=True)


def sanitize_tag(tag: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", tag)
    cleaned = cleaned.strip("_")
    return cleaned or "model"


def append_csv(src: Path, dest: Path, include_header: bool) -> None:
    with src.open("r", encoding="utf-8") as src_file, dest.open(
        "a", encoding="utf-8"
    ) as dest_file:
        first_line = src_file.readline()
        if include_header:
            dest_file.write(first_line)
        for line in src_file:
            dest_file.write(line)


def save_npz_fast(output_path: Path, **arrays: object) -> None:
    import numpy as np

    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(output_path, **arrays)


def merge_model_embeddings(output_path: Path, chunk_paths: Sequence[Path]) -> None:
    if not chunk_paths:
        raise PipelineError("No embedding chunks were generated for this model group.")
    import numpy as np

    merged_arrays: dict[str, "np.ndarray"] = {}
    merged_embeddings: dict[str, "np.ndarray"] = {}
    merged_layer_embeddings: dict[str, dict[str, "np.ndarray"]] = {}

    for idx, chunk_path in enumerate(chunk_paths):
        data = np.load(chunk_path, allow_pickle=True)
        for key in data.files:
            value = data[key]
            if key == "embeddings":
                emb_dict = value.item()
                if idx == 0:
                    merged_embeddings = {k: emb_dict[k].copy() for k in emb_dict}
                else:
                    for model_key, arr in emb_dict.items():
                        if model_key in merged_embeddings:
                            merged_embeddings[model_key] = np.concatenate(
                                (merged_embeddings[model_key], arr), axis=0
                            )
                        else:
                            merged_embeddings[model_key] = arr
            elif key == "layer_embeddings":
                layer_dict = value.item()
                if idx == 0:
                    merged_layer_embeddings = {
                        model_key: {layer_name: arr.copy() for layer_name, arr in model_layers.items()}
                        for model_key, model_layers in layer_dict.items()
                    }
                else:
                    for model_key, model_layers in layer_dict.items():
                        target_layers = merged_layer_embeddings.setdefault(model_key, {})
                        for layer_name, arr in model_layers.items():
                            if layer_name in target_layers:
                                target_layers[layer_name] = np.concatenate(
                                    (target_layers[layer_name], arr), axis=0
                                )
                            else:
                                target_layers[layer_name] = arr
            else:
                if idx == 0:
                    merged_arrays[key] = value.copy()
                else:
                    merged_arrays[key] = np.concatenate((merged_arrays[key], value), axis=0)

    merged_arrays["embeddings"] = np.array(merged_embeddings, dtype=object)
    if merged_layer_embeddings:
        merged_arrays["layer_embeddings"] = np.array(merged_layer_embeddings, dtype=object)
    save_npz_fast(output_path, **merged_arrays)


def merge_all_embeddings(output_path: Path, input_paths: Sequence[Path]) -> None:
    if not input_paths:
        raise PipelineError("No embedding archives provided for consolidation.")
    import numpy as np

    merged_arrays: dict[str, object] = {}
    merged_embeddings: dict[str, "np.ndarray"] = {}
    merged_layer_embeddings: dict[str, dict[str, "np.ndarray"]] = {}
    accessions_by_model: dict[str, "np.ndarray"] = {}
    sequences_by_model: dict[str, "np.ndarray"] = {}
    accession_sequence_map: dict[str, str] = {}
    base_accessions: Optional["np.ndarray"] = None
    base_sequences: Optional["np.ndarray"] = None
    aligned_accessions = True
    aligned_sequences = True

    for path in input_paths:
        data = np.load(path, allow_pickle=True)
        accessions = data["accessions"]
        sequences = data["sequences"]
        if base_accessions is None:
            base_accessions = accessions
            base_sequences = sequences
        else:
            if not np.array_equal(base_accessions, accessions):
                aligned_accessions = False
            if not np.array_equal(base_sequences, sequences):
                aligned_sequences = False

        for accession, sequence in zip(accessions, sequences):
            accession_sequence_map.setdefault(str(accession), str(sequence))

        model_prefixes = sorted(
            key[: -len("_ids")] for key in data.files if key.endswith("_ids")
        )
        for model_key in model_prefixes:
            accessions_by_model[model_key] = accessions
            sequences_by_model[model_key] = sequences
            merged_arrays[f"{model_key}_accessions"] = accessions
            merged_arrays[f"{model_key}_sequences"] = sequences

        for key in data.files:
            value = data[key]
            if key == "embeddings":
                for model_key, arr in value.item().items():
                    if model_key in merged_embeddings:
                        merged_embeddings[model_key] = np.concatenate(
                            (merged_embeddings[model_key], arr), axis=0
                        )
                    else:
                        merged_embeddings[model_key] = arr
            elif key == "layer_embeddings":
                for model_key, model_layers in value.item().items():
                    target_layers = merged_layer_embeddings.setdefault(model_key, {})
                    for layer_name, arr in model_layers.items():
                        if layer_name in target_layers:
                            target_layers[layer_name] = np.concatenate(
                                (target_layers[layer_name], arr), axis=0
                            )
                        else:
                            target_layers[layer_name] = arr
            elif key in {"accessions", "sequences"}:
                continue
            elif key in merged_arrays:
                merged_arrays[key] = np.concatenate((merged_arrays[key], value), axis=0)
            else:
                merged_arrays[key] = value

    if aligned_accessions and aligned_sequences and base_accessions is not None:
        merged_arrays["accessions"] = base_accessions
        merged_arrays["sequences"] = base_sequences
    else:
        merged_arrays["accessions"] = np.array(list(accession_sequence_map), dtype=object)
        merged_arrays["sequences"] = np.array(
            [accession_sequence_map[accession] for accession in accession_sequence_map],
            dtype=object,
        )
        merged_arrays["accessions_by_model"] = np.array(accessions_by_model, dtype=object)
        merged_arrays["sequences_by_model"] = np.array(sequences_by_model, dtype=object)
    merged_arrays["embeddings"] = np.array(merged_embeddings, dtype=object)
    if merged_layer_embeddings:
        merged_arrays["layer_embeddings"] = np.array(merged_layer_embeddings, dtype=object)
    save_npz_fast(output_path, **merged_arrays)


def run_merge_model_embeddings(
    venv_python: Path,
    output_path: Path,
    chunk_paths: Sequence[Path],
) -> None:
    source_dir = str(FANTASIA_ROOT / "src")
    code = (
        "import sys; "
        "from pathlib import Path; "
        f"sys.path.insert(0, {source_dir!r}); "
        "from fantasia_pipeline import merge_model_embeddings; "
        "merge_model_embeddings(Path(sys.argv[1]), [Path(p) for p in sys.argv[2:]])"
    )
    run_subprocess([str(venv_python), "-c", code, str(output_path), *map(str, chunk_paths)])


def run_merge_all_embeddings(
    venv_python: Path,
    output_path: Path,
    input_paths: Sequence[Path],
) -> None:
    source_dir = str(FANTASIA_ROOT / "src")
    code = (
        "import sys; "
        "from pathlib import Path; "
        f"sys.path.insert(0, {source_dir!r}); "
        "from fantasia_pipeline import merge_all_embeddings; "
        "merge_all_embeddings(Path(sys.argv[1]), [Path(p) for p in sys.argv[2:]])"
    )
    run_subprocess([str(venv_python), "-c", code, str(output_path), *map(str, input_paths)])


def merge_failure_reports(output_path: Path, report_paths: Sequence[Path]) -> int:
    """Combine per-chunk failure reports into a single CSV."""
    import csv

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sequence_id",
        "model_name",
        "error_type",
        "error_category",
        "error_message",
        "sequence_length",
        "fasta_path",
    ]

    rows_written = 0
    with output_path.open("w", newline="", encoding="utf-8") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames)
        writer.writeheader()
        for report in report_paths:
            if not report.exists():
                continue
            with report.open("r", encoding="utf-8") as in_handle:
                reader = csv.DictReader(in_handle)
                for row in reader:
                    if not row:
                        continue
                    writer.writerow(row)
                    rows_written += 1
    if rows_written == 0:
        print(f"No sequences were skipped; wrote empty report to {output_path}.")
    else:
        print(
            f"Wrote aggregate failure report with {rows_written} entr{'y' if rows_written == 1 else 'ies'} to {output_path}."
        )
    return rows_written


def write_config_yaml(config: PipelineConfig) -> None:
    lines = [
        f"lookup_table_path: {config.lookup_npz}",
        f"annotations_path: {config.annotations_json}",
        f"accession_path: {config.accessions_json}",
        f"embeddings_path: {config.embeddings_npz}",
        f"limit_per_entry: {config.limit_per_entry}",
        "embedding:",
        f"  distance_metric: {config.distance_metric}",
        f"results_path: {config.results_csv}",
    ]
    if config.raw_results_csv:
        lines.append(f"raw_results_path: {config.raw_results_csv}")
    config.config_yaml.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _topgo_filename(model_key: str, category: str, prefix: Optional[str] = None) -> str:
    safe_model = re.sub(r"[^A-Za-z0-9\-]+", "_", model_key).strip("_")
    safe_model = safe_model or "model"
    pieces = []
    if prefix:
        pieces.append(prefix)
    pieces.append(f"{safe_model}.topgo.{category.upper()}.txt")
    return ".".join(pieces)


def write_topgo_files(config: PipelineConfig) -> None:
    write_topgo_variant(config, Path(config.results_csv), variant_prefix=None)
    if config.raw_results_csv:
        raw_path = Path(config.raw_results_csv)
        if raw_path.exists():
            write_topgo_variant(config, raw_path, variant_prefix="k")
        else:
            print(f"Raw results file not found at {raw_path}; skipping k TopGO export.")


def write_topgo_variant(
    config: PipelineConfig,
    results_path: Path,
    variant_prefix: Optional[str],
) -> None:
    if not results_path.exists():
        print(f"Results file not found at {results_path}; skipping TopGO export.")
        return

    topgo_dir = Path(config.topgo_dir)
    topgo_dir.mkdir(parents=True, exist_ok=True)

    required_fields = {"query_accession", "model_key", "category", "go_id"}

    aggregated: dict[tuple[str, str], OrderedDict[str, List[str]]] = {}
    with results_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not required_fields.issubset(set(reader.fieldnames)):
            missing = required_fields - set(reader.fieldnames or [])
            raise PipelineError(
                f"results.csv is missing required columns for TopGO export: {', '.join(sorted(missing))}"
            )
        for row in reader:
            query = (row.get("query_accession") or "").strip()
            model = (row.get("model_key") or "").strip()
            category = (row.get("category") or "").strip().upper()
            go_id = (row.get("go_id") or "").strip()
            if not (query and model and category and go_id):
                continue
            key = (model, category)
            per_query = aggregated.setdefault(key, OrderedDict())
            go_list = per_query.setdefault(query, [])
            if go_id not in go_list:
                go_list.append(go_id)

    if not aggregated:
        print(f"No GO annotations found in {results_path}; skipping TopGO export.")
        return

    for (model, category), queries in aggregated.items():
        filename = _topgo_filename(model, category, prefix=variant_prefix)
        output_path = topgo_dir / filename
        with output_path.open("w", encoding="utf-8") as handle:
            for query, go_terms in queries.items():
                handle.write(f"{query}\t{','.join(go_terms)}\n")
        label = f"variant '{variant_prefix}'" if variant_prefix else "primary results"
        print(
            f"Wrote TopGO table for model '{model}' category '{category}' ({label}) to {output_path}"
        )


def build_model_groups(embed_models: str, serial_models: bool) -> List[str]:
    models = embed_models.strip()
    if not models:
        raise PipelineError("No embedding models configured.")
    requested = parse_embed_model_names(models)
    if serial_models:
        return requested
    return [" ".join(requested)]


def run_chunk_pipeline(
    config: PipelineConfig,
    venv_python: Path,
    device: str,
) -> PipelineRunStats:
    ensure_output_directories(config)
    use_gpu_lookup = should_use_gpu_lookup(config, device)
    embedding_time_seconds = 0.0
    lookup_time_seconds = 0.0
    postprocess_time_seconds = 0.0
    print(
        f"Lookup stage will use {'GPU' if use_gpu_lookup else 'CPU'} in-process lookup."
    )

    if config.reuse_embeddings:
        if not config.embeddings_npz.exists():
            raise PipelineError(
                f"--reuse-embeddings was requested but embeddings archive was not found: "
                f"{config.embeddings_npz}"
            )
        print(f"Reusing existing embeddings archive at {config.embeddings_npz}")
        lookup_conf: Dict[str, object] = {
            "lookup_table_path": str(config.lookup_npz),
            "annotations_path": str(config.annotations_json),
            "accession_path": str(config.accessions_json),
            "embeddings_path": str(config.embeddings_npz),
            "lookup_device": device,
            "lookup_batch_size": config.lookup_batch_size,
            "limit_per_entry": config.limit_per_entry,
            "embedding": {"distance_metric": config.distance_metric},
            "results_path": str(config.results_csv),
        }
        if config.raw_results_csv is not None:
            lookup_conf["raw_results_path"] = str(config.raw_results_csv)
        lookup_started = time.perf_counter()
        run_lookup_subprocess(
            lookup_conf,
            use_gpu_lookup=use_gpu_lookup,
            venv_python=venv_python,
            config_path=config.chunk_config_dir / "reuse_lookup_config.yaml",
            device=device,
        )
        lookup_time_seconds += time.perf_counter() - lookup_started
        skipped_sequences = count_report_rows(config.failure_report)
        total_sequences = count_fasta_sequences(config.fasta_path)
        return PipelineRunStats(
            total_sequences=total_sequences,
            skipped_sequences=skipped_sequences,
            embedding_time_seconds=embedding_time_seconds,
            lookup_time_seconds=lookup_time_seconds,
            postprocess_time_seconds=postprocess_time_seconds,
        )

    chunks = chunk_fasta(config)
    if not chunks:
        raise PipelineError(f"No sequences found in FASTA file: {config.fasta_path}")

    model_groups = build_model_groups(config.embed_models, config.serial_models)
    print(f"Processing {len(chunks)} FASTA chunk(s) with models: {model_groups}")

    master_results_tmp = config.results_csv.with_suffix(".csv.tmp")
    if master_results_tmp.exists():
        master_results_tmp.unlink()
    master_results_tmp.touch()
    master_header_written = False

    master_raw_tmp: Optional[Path] = None
    master_raw_header_written = False
    if config.raw_results_csv:
        master_raw_tmp = config.raw_results_csv.parent / f"{config.raw_results_csv.name}.tmp"
        if master_raw_tmp.exists():
            master_raw_tmp.unlink()
        master_raw_tmp.touch()

    master_embed_files: List[Path] = []
    failure_report_paths: List[Path] = []

    for model_group in model_groups:
        model_tag = sanitize_tag(model_group)
        print(f"--- Model group: {model_group} (tag: {model_tag}) ---")

        model_chunk_embeds: List[Path] = []
        chunk_jobs: List[ChunkJob] = []

        for chunk_path in chunks:
            chunk_name = chunk_path.stem
            chunk_embeddings = (
                config.chunk_embed_dir / f"{model_tag}_{chunk_name}_embeddings.npz"
            )
            chunk_config = config.chunk_config_dir / f"{model_tag}_{chunk_name}_config.yaml"
            chunk_results = config.chunk_results_dir / f"{model_tag}_{chunk_name}_results.csv"
            chunk_failures = config.chunk_failure_dir / f"{model_tag}_{chunk_name}_failures.csv"
            if chunk_failures.exists():
                chunk_failures.unlink()

            chunk_raw_results: Optional[Path] = None
            if master_raw_tmp is not None:
                chunk_raw_results = (
                    config.chunk_results_dir / f"{model_tag}_{chunk_name}_raw.csv"
                )

            chunk_jobs.append(
                ChunkJob(
                    chunk_path=chunk_path,
                    chunk_embeddings=chunk_embeddings,
                    chunk_config=chunk_config,
                    chunk_results=chunk_results,
                    chunk_failures=chunk_failures,
                    chunk_raw_results=chunk_raw_results,
                )
            )

        chunk_spec_path = config.chunk_embed_dir / f"{model_tag}_chunks.json"
        chunk_payload = [
            {
                "fasta": str(job.chunk_path),
                "output": str(job.chunk_embeddings),
                "failure_report": str(job.chunk_failures),
            }
            for job in chunk_jobs
        ]
        chunk_spec_path.write_text(json.dumps(chunk_payload), encoding="utf-8")

        print(
            f"Generating embeddings for {len(chunk_jobs)} chunk(s) with models: {model_group}"
        )
        embedding_cmd = [
            str(venv_python),
            FANTASIA_ROOT / "src" / "generate_embeddings.py",
            "--chunks-file",
            str(chunk_spec_path),
            "--device",
            device,
        ]
        embedding_cmd.extend(["--models", *model_group.split()])
        embedding_cmd.extend(
            [
                "--embed-batch-size",
                str(config.embed_batch_size),
                "--max-batch-residues",
                str(config.max_batch_residues),
                "--queue-batch-size",
                str(config.queue_batch_size),
                "--max-sequence-length",
                str(config.max_sequence_length),
            ]
        )
        if config.model_batch_sizes:
            embedding_cmd.extend(["--model-batch-sizes", *config.model_batch_sizes])
        embedding_started = time.perf_counter()
        run_subprocess(embedding_cmd, env=torch_runtime_env(device))
        embedding_time_seconds += time.perf_counter() - embedding_started

        for job in chunk_jobs:
            if job.chunk_failures.exists():
                failure_report_paths.append(job.chunk_failures)
            model_chunk_embeds.append(job.chunk_embeddings)

        model_merged_emb = config.chunk_embed_dir / f"{model_tag}_merged.npz"
        run_merge_model_embeddings(venv_python, model_merged_emb, model_chunk_embeds)

        model_merged_results = config.chunk_results_dir / f"{model_tag}_merged.csv"
        model_merged_raw: Optional[Path] = None
        lookup_conf: Dict[str, object] = {
            "lookup_table_path": str(config.lookup_npz),
            "annotations_path": str(config.annotations_json),
            "accession_path": str(config.accessions_json),
            "embeddings_path": str(model_merged_emb),
            "lookup_device": device,
            "lookup_batch_size": config.lookup_batch_size,
            "limit_per_entry": config.limit_per_entry,
            "embedding": {"distance_metric": config.distance_metric},
            "results_path": str(model_merged_results),
        }
        if master_raw_tmp is not None:
            model_merged_raw = config.chunk_results_dir / f"{model_tag}_raw.csv"
            lookup_conf["raw_results_path"] = str(model_merged_raw)

        print(f"Running annotation lookup once for merged model group {model_group}")
        lookup_started = time.perf_counter()
        run_lookup_subprocess(
            lookup_conf,
            use_gpu_lookup=use_gpu_lookup,
            venv_python=venv_python,
            config_path=config.chunk_config_dir / f"{model_tag}_lookup_config.yaml",
            device=device,
        )
        lookup_time_seconds += time.perf_counter() - lookup_started

        append_csv(
            model_merged_results,
            master_results_tmp,
            include_header=not master_header_written,
        )
        master_header_written = True
        master_embed_files.append(model_merged_emb)

        if master_raw_tmp is not None and model_merged_raw is not None:
            append_csv(
                model_merged_raw,
                master_raw_tmp,
                include_header=not master_raw_header_written,
            )
            master_raw_header_written = True

    postprocess_started = time.perf_counter()
    master_results_tmp.rename(config.results_csv)
    if master_raw_tmp is not None:
        master_raw_tmp.rename(config.raw_results_csv)
    run_merge_all_embeddings(venv_python, config.embeddings_npz, master_embed_files)
    skipped_sequences = merge_failure_reports(config.failure_report, failure_report_paths)
    postprocess_time_seconds += time.perf_counter() - postprocess_started
    total_sequences = count_fasta_sequences(config.fasta_path)
    return PipelineRunStats(
        total_sequences=total_sequences,
        skipped_sequences=skipped_sequences,
        embedding_time_seconds=embedding_time_seconds,
        lookup_time_seconds=lookup_time_seconds,
        postprocess_time_seconds=postprocess_time_seconds,
    )


def run_pipeline(config: PipelineConfig, *, argv: Sequence[str]) -> None:
    # Ensure base directories exist for lookup artifacts and outputs.
    for path in (
        config.lookup_npz,
        config.annotations_json,
        config.accessions_json,
        config.embeddings_npz,
        config.config_yaml,
        config.results_csv,
        config.raw_results_csv,
        config.failure_report,
        config.topgo_dir,
        config.run_log,
        config.run_metadata_yaml,
    ):
        if path:
            path.parent.mkdir(parents=True, exist_ok=True)
    started_at = datetime.now()
    global _ACTIVE_LOG_HANDLE
    with config.run_log.open("w", encoding="utf-8") as log_handle:
        _ACTIVE_LOG_HANDLE = log_handle
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        sys.stdout = TeeStream(original_stdout, log_handle)
        sys.stderr = TeeStream(original_stderr, log_handle)
        try:
            write_run_metadata(
                config.run_metadata_yaml,
                build_run_metadata(
                    config,
                    argv=argv,
                    started_at=started_at,
                    status="started",
                ),
            )
            device: Optional[str] = None
            stats: Optional[PipelineRunStats] = None
            try:
                requested_models = parse_embed_model_names(config.embed_models)
                validate_python_for_models(requested_models)
                create_virtualenv(config.venv_dir)
                device = install_packages(config)
                venv_python = venv_python_executable(config.venv_dir)
                ensure_lookup_artifacts(config)
                stats = run_chunk_pipeline(config, venv_python, device)
                write_config_yaml(config)
                if config.generate_topgo:
                    topgo_started = time.perf_counter()
                    write_topgo_files(config)
                    stats.postprocess_time_seconds += time.perf_counter() - topgo_started
                else:
                    print("TopGO export disabled; skipping TopGO file generation.")
                print(
                    "Stage timing (s): "
                    f"embedding={stats.embedding_time_seconds:.2f}, "
                    f"lookup={stats.lookup_time_seconds:.2f}, "
                    f"postprocess={stats.postprocess_time_seconds:.2f}, "
                    f"measured_total={stats.measured_stage_time_seconds:.2f}"
                )
                print(
                    "Sequence summary: "
                    f"{stats.processed_sequences}/{stats.total_sequences} processed "
                    f"({stats.processed_fraction * 100:.2f}%), "
                    f"{stats.skipped_sequences}/{stats.total_sequences} skipped "
                    f"({stats.skipped_fraction * 100:.2f}%)."
                )
                print(f"Pipeline complete. Results saved to {config.results_csv}")
                finished_at = datetime.now()
                write_run_metadata(
                    config.run_metadata_yaml,
                    build_run_metadata(
                        config,
                        argv=argv,
                        started_at=started_at,
                        finished_at=finished_at,
                        device=device,
                        stats=stats,
                        status="completed",
                    ),
                )
            except Exception:
                finished_at = datetime.now()
                write_run_metadata(
                    config.run_metadata_yaml,
                    build_run_metadata(
                        config,
                        argv=argv,
                        started_at=started_at,
                        finished_at=finished_at,
                        device=device,
                        stats=stats,
                        status="failed",
                    ),
                )
                raise
        finally:
            sys.stdout = original_stdout
            sys.stderr = original_stderr
            _ACTIVE_LOG_HANDLE = None


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        argv = list(argv) if argv is not None else sys.argv[1:]
        args = parse_args(argv)
        config = PipelineConfig.from_args(args)
        run_pipeline(config, argv=argv)
        return 0
    except PipelineError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as exc:
        cmd_text = " ".join(str(part) for part in exc.cmd)
        print(f"Command failed with exit code {exc.returncode}: {cmd_text}", file=sys.stderr)
        return exc.returncode or 1


if __name__ == "__main__":
    sys.exit(main())
