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
import subprocess
import sys
import time
from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, TextIO


DEFAULT_EMBED_MODELS = "prot_t5"
SUPPORTED_EMBED_MODELS = {"prot_t5", "esm2", "ankh3"}
DEFAULT_CHUNK_SIZE = 500
DEFAULT_LIMIT_PER_ENTRY = 1
DEFAULT_DISTANCE_METRIC = "cosine"

DEFAULT_LOOKUP_DIR = Path("data") / "lookup"
# Generate timestamped output directory
_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
DEFAULT_OUTPUT_DIR = Path(f"outputs_{_timestamp}")
DEFAULT_TEMP_DIR = DEFAULT_OUTPUT_DIR / "tmp"
DEFAULT_FAILURE_REPORT = DEFAULT_OUTPUT_DIR / "failed_sequences.csv"
DEFAULT_CHUNK_FAILURE_DIR = DEFAULT_TEMP_DIR / "failures"
DEFAULT_TOPGO_DIR = DEFAULT_OUTPUT_DIR / "topgo"


class PipelineError(RuntimeError):
    """Raised when a pipeline stage fails in a recoverable way."""


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
    limit_per_entry: int
    distance_metric: str
    torch_index: Optional[str]
    chunk_failure_dir: Path
    failure_report: Path

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "PipelineConfig":
        fasta_path = Path(args.fasta_path).resolve()
        if not fasta_path.exists():
            raise PipelineError(f"FASTA file not found: {fasta_path}")

        env = os.environ
        embed_models = args.embed_models or env.get("EMBED_MODELS", DEFAULT_EMBED_MODELS)
        serial_models = args.serial_models or env.get("SERIAL_MODELS") == "1"

        device = args.device or env.get("PYTORCH_DEVICE")
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
            limit_per_entry=args.limit_per_entry,
            distance_metric=args.distance_metric,
            torch_index=torch_index,
            chunk_failure_dir=Path(args.chunk_failure_dir).resolve(),
            failure_report=Path(args.failure_report).resolve(),
        )


@dataclasses.dataclass
class ChunkJob:
    chunk_path: Path
    chunk_embeddings: Path
    chunk_config: Path
    chunk_results: Path
    chunk_failures: Path
    chunk_raw_results: Optional[Path]


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the FANTASIA pipeline purely in Python."
    )
    parser.add_argument("fasta_path", help="Path to the FASTA file to annotate.")
    parser.add_argument(
        "--venv-dir",
        default="venv",
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
            "Space-separated list of embedding models (choices: prot_t5, esm2, ankh3; "
            "default from env or 'prot_t5')."
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
        "--torch-index",
        default=None,
        help="Custom extra-index for installing CUDA-enabled PyTorch wheels.",
    )
    parser.add_argument(
        "--failure-report",
        default=str(DEFAULT_FAILURE_REPORT),
        help="Path to a CSV summarising sequences skipped during embedding.",
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
    result = subprocess.run(
        cmd,
        env=process_env,
        cwd=str(cwd) if cwd else None,
        check=check,
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


def install_packages(config: PipelineConfig) -> str:
    device = detect_device(config.device)
    pip_bin = venv_pip_executable(config.venv_dir)

    print("Upgrading pip…")
    run_subprocess([str(pip_bin), "install", "--upgrade", "pip"])

    base_deps = [
        "numpy",
        "requests",
        "transformers",
        "biopython",
        "pandas",
        "scipy",
        "sentencepiece",
        "protobuf",
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

    torch_packages = ["torch", "torchvision", "torchaudio"]
    if device.lower() == "cuda":
        torch_index = config.torch_index or "https://download.pytorch.org/whl/cu124"
        print(f"Installing PyTorch with CUDA support from {torch_index}…")
        run_subprocess(
            [str(pip_bin), "install", "--extra-index-url", torch_index, *torch_packages]
        )
    else:
        print("Installing CPU-only PyTorch packages…")
        run_subprocess([str(pip_bin), "install", *torch_packages])

    return device


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


def merge_model_embeddings(output_path: Path, chunk_paths: Sequence[Path]) -> None:
    if not chunk_paths:
        raise PipelineError("No embedding chunks were generated for this model group.")
    import numpy as np

    merged_arrays: dict[str, "np.ndarray"] = {}
    merged_embeddings: dict[str, "np.ndarray"] = {}

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
            else:
                if idx == 0:
                    merged_arrays[key] = value.copy()
                else:
                    merged_arrays[key] = np.concatenate((merged_arrays[key], value), axis=0)

    merged_arrays["embeddings"] = np.array(merged_embeddings, dtype=object)
    np.savez_compressed(output_path, **merged_arrays)


def merge_all_embeddings(output_path: Path, input_paths: Sequence[Path]) -> None:
    if not input_paths:
        raise PipelineError("No embedding archives provided for consolidation.")
    import numpy as np

    base = np.load(input_paths[0], allow_pickle=True)
    accessions = base["accessions"]
    sequences = base["sequences"]
    merged_arrays = {k: base[k] for k in base.files if k not in {"embeddings", "accessions", "sequences"}}
    merged_embeddings = base["embeddings"].item()

    for path in input_paths[1:]:
        data = np.load(path, allow_pickle=True)
        if not np.array_equal(accessions, data["accessions"]):
            raise PipelineError("Accession lists do not match between embedding archives.")
        if not np.array_equal(sequences, data["sequences"]):
            raise PipelineError("Sequence lists do not match between embedding archives.")
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
            elif key in {"accessions", "sequences"}:
                continue
            elif key in merged_arrays:
                merged_arrays[key] = np.concatenate((merged_arrays[key], value), axis=0)
            else:
                merged_arrays[key] = value

    merged_arrays["accessions"] = accessions
    merged_arrays["sequences"] = sequences
    merged_arrays["embeddings"] = np.array(merged_embeddings, dtype=object)
    np.savez_compressed(output_path, **merged_arrays)


def merge_failure_reports(output_path: Path, report_paths: Sequence[Path]) -> None:
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
    requested = models.split()
    unknown = [model for model in requested if model not in SUPPORTED_EMBED_MODELS]
    if unknown:
        supported = ", ".join(sorted(SUPPORTED_EMBED_MODELS))
        raise PipelineError(
            f"Unsupported embedding model(s): {', '.join(unknown)}. Supported models: {supported}."
        )
    if serial_models:
        return requested
    return [" ".join(requested)]


def run_chunk_pipeline(
    config: PipelineConfig,
    venv_python: Path,
    device: str,
) -> None:
    ensure_output_directories(config)
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

        model_tmp_results = config.chunk_results_dir / f"{model_tag}_merged.csv.tmp"
        if model_tmp_results.exists():
            model_tmp_results.unlink()
        model_tmp_results.touch()
        model_header_written = False

        model_raw_tmp: Optional[Path] = None
        model_raw_header_written = False
        if master_raw_tmp is not None:
            model_raw_tmp = config.chunk_results_dir / f"{model_tag}_raw.csv.tmp"
            if model_raw_tmp.exists():
                model_raw_tmp.unlink()
            model_raw_tmp.touch()

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
            "generate_embeddings.py",
            "--chunks-file",
            str(chunk_spec_path),
            "--device",
            device,
        ]
        embedding_cmd.extend(["--models", *model_group.split()])
        run_subprocess(embedding_cmd)

        for job in chunk_jobs:
            chunk_path = job.chunk_path
            chunk_embeddings = job.chunk_embeddings
            chunk_config = job.chunk_config
            chunk_results = job.chunk_results
            chunk_failures = job.chunk_failures
            chunk_raw_results = job.chunk_raw_results

            if chunk_failures.exists():
                failure_report_paths.append(chunk_failures)

            chunk_config_lines = [
                f"lookup_table_path: {config.lookup_npz}",
                f"annotations_path: {config.annotations_json}",
                f"accession_path: {config.accessions_json}",
                f"embeddings_path: {chunk_embeddings}",
                f"limit_per_entry: {config.limit_per_entry}",
                "embedding:",
                f"  distance_metric: {config.distance_metric}",
                f"results_path: {chunk_results}",
            ]
            if chunk_raw_results is not None:
                chunk_config_lines.append(f"raw_results_path: {chunk_raw_results}")
            chunk_config_lines.append("")
            chunk_config.write_text("\n".join(chunk_config_lines), encoding="utf-8")

            print(f"Running annotation lookup for {chunk_path} (models: {model_group})")
            run_subprocess([str(venv_python), "fantasia_no_db.py", "--config", str(chunk_config)])

            append_csv(chunk_results, model_tmp_results, include_header=not model_header_written)
            model_header_written = True

            if chunk_raw_results is not None and chunk_raw_results.exists():
                if model_raw_tmp is None:
                    raise PipelineError("Internal error: raw results tmp path not prepared.")
                append_csv(
                    chunk_raw_results,
                    model_raw_tmp,
                    include_header=not model_raw_header_written,
                )
                model_raw_header_written = True

            model_chunk_embeds.append(chunk_embeddings)

        model_merged_results = config.chunk_results_dir / f"{model_tag}_merged.csv"
        model_tmp_results.rename(model_merged_results)

        model_merged_emb = config.chunk_embed_dir / f"{model_tag}_merged.npz"
        merge_model_embeddings(model_merged_emb, model_chunk_embeds)

        append_csv(
            model_merged_results,
            master_results_tmp,
            include_header=not master_header_written,
        )
        master_header_written = True
        master_embed_files.append(model_merged_emb)

        if master_raw_tmp is not None and model_raw_tmp is not None:
            model_merged_raw = config.chunk_results_dir / f"{model_tag}_raw.csv"
            model_raw_tmp.rename(model_merged_raw)
            append_csv(
                model_merged_raw,
                master_raw_tmp,
                include_header=not master_raw_header_written,
            )
            master_raw_header_written = True

    master_results_tmp.rename(config.results_csv)
    if master_raw_tmp is not None:
        master_raw_tmp.rename(config.raw_results_csv)
    merge_all_embeddings(config.embeddings_npz, master_embed_files)
    merge_failure_reports(config.failure_report, failure_report_paths)


def run_pipeline(config: PipelineConfig) -> None:
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
    ):
        if path:
            path.parent.mkdir(parents=True, exist_ok=True)

    create_virtualenv(config.venv_dir)
    device = install_packages(config)
    venv_python = venv_python_executable(config.venv_dir)
    ensure_lookup_artifacts(config)
    run_chunk_pipeline(config, venv_python, device)
    write_config_yaml(config)
    write_topgo_files(config)
    print(f"Pipeline complete. Results saved to {config.results_csv}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        args = parse_args(argv)
        config = PipelineConfig.from_args(args)
        run_pipeline(config)
        return 0
    except PipelineError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as exc:
        print(f"Command failed with exit code {exc.returncode}: {' '.join(exc.cmd)}", file=sys.stderr)
        return exc.returncode or 1


if __name__ == "__main__":
    sys.exit(main())
