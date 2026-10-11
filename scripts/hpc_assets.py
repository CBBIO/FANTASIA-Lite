#!/usr/bin/env python3
"""Prepare and validate FANTASIA assets outside isolated compute nodes.

This program runs inside the Apptainer image but is invoked only by the HPC
launcher on a connected login node. It never performs a download in
``validate`` mode.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import tarfile
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

LOOKUP_FILES = ("lookup_table.npz", "annotations.json", "accessions.json")


def load_model_registry() -> tuple[dict, object]:
    """Load the pipeline's pinned registry only inside the container runtime."""
    from generate_embeddings import MODEL_REGISTRY, sha256_file  # noqa: E402

    return MODEL_REGISTRY, sha256_file


def require_lookup(lookup_dir: Path) -> None:
    missing = [name for name in LOOKUP_FILES if not (lookup_dir / name).is_file()]
    if missing:
        raise RuntimeError(
            f"Lookup directory {lookup_dir} is incomplete; missing: {', '.join(missing)}"
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_extract(archive: Path, destination: Path) -> None:
    destination_resolved = destination.resolve()
    with tarfile.open(archive, "r:*") as tar:
        for member in tar.getmembers():
            target = (destination / member.name).resolve()
            if not target.is_relative_to(destination_resolved):
                raise RuntimeError(f"Unsafe archive member: {member.name}")
        tar.extractall(destination)


def download_lookup(url: str, expected_sha256: str, archive: Path, lookup_dir: Path) -> None:
    if not url or not expected_sha256 or len(expected_sha256) != 64:
        raise RuntimeError(
            "prepare requires assets.lookup.url and a 64-character assets.lookup.sha256."
        )
    archive.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url) as response, archive.open("wb") as output:
        shutil.copyfileobj(response, output)
    actual_sha256 = sha256(archive)
    if actual_sha256.lower() != expected_sha256.lower():
        raise RuntimeError(
            f"Lookup archive checksum mismatch: expected {expected_sha256}, got {actual_sha256}."
        )
    with tempfile.TemporaryDirectory(dir=archive.parent, prefix="fantasia-lookup-") as raw_temp:
        temp = Path(raw_temp)
        safe_extract(archive, temp)
        candidates = (temp / "data" / "lookup", temp / "lookup", temp)
        source = next(
            (candidate for candidate in candidates if all((candidate / item).is_file() for item in LOOKUP_FILES)),
            None,
        )
        if source is None:
            raise RuntimeError("The lookup archive does not contain the required lookup files.")
        lookup_dir.mkdir(parents=True, exist_ok=True)
        for item in LOOKUP_FILES:
            shutil.copy2(source / item, lookup_dir / item)
    require_lookup(lookup_dir)


def requested_models(names: Iterable[str], registry: dict) -> list[str]:
    selected = list(dict.fromkeys(names))
    unknown = [name for name in selected if name not in registry]
    if unknown:
        raise RuntimeError(f"Unknown model(s): {', '.join(unknown)}")
    return selected


def snapshot_for(
    model_name: str,
    cache_dir: Path,
    registry: dict,
    sha256_file_fn: object,
    *,
    local_files_only: bool,
) -> Path:
    from huggingface_hub import snapshot_download

    record = registry[model_name]
    try:
        path = Path(
            snapshot_download(
                repo_id=str(record["hf_id"]),
                revision=str(record["revision"]),
                cache_dir=str(cache_dir),
                local_files_only=local_files_only,
            )
        ).resolve()
    except Exception as exc:
        action = "prepare" if local_files_only else "download"
        raise RuntimeError(f"Could not {action} pinned model {model_name}: {exc}") from exc
    if path.name != str(record["revision"]):
        raise RuntimeError(
            f"Model {model_name} resolved revision {path.name}, expected {record['revision']}."
        )
    if model_name == "esmc":
        weight = path / "data" / "weights" / str(record["serialization"])
        if not weight.is_file():
            raise RuntimeError(f"Pinned ESM-C weights are missing: {weight}")
        actual = sha256_file_fn(weight)  # type: ignore[operator]
        if actual != record["weights_sha256"]:
            raise RuntimeError(f"ESM-C checksum mismatch: expected {record['weights_sha256']}, got {actual}")
    return path


def write_manifest(
    path: Path, lookup_dir: Path, cache_dir: Path, models: list[str], registry: dict
) -> None:
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "lookup": {name: sha256(lookup_dir / name) for name in LOOKUP_FILES},
        "models": {
            name: {
                "repository": registry[name]["hf_id"],
                "revision": registry[name]["revision"],
                "cache_dir": str(cache_dir),
            }
            for name in models
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    registry, sha256_file_fn = load_model_registry()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "validate"))
    parser.add_argument("--lookup-dir", type=Path, required=True)
    parser.add_argument("--model-cache-dir", type=Path, required=True)
    parser.add_argument("--staging-dir", type=Path, required=True)
    parser.add_argument("--lookup-url")
    parser.add_argument("--lookup-sha256")
    parser.add_argument("--lookup-archive-name", default="fantasia_lite_data_folder.tgz")
    parser.add_argument("--models", nargs="+", default=list(registry))
    args = parser.parse_args()

    models = requested_models(args.models, registry)
    os.environ["HF_HOME"] = str(args.model_cache_dir)
    args.model_cache_dir.mkdir(parents=True, exist_ok=True)
    if args.command == "prepare":
        if not all((args.lookup_dir / name).is_file() for name in LOOKUP_FILES):
            download_lookup(
                args.lookup_url or "",
                args.lookup_sha256 or "",
                args.staging_dir / args.lookup_archive_name,
                args.lookup_dir,
            )
        require_lookup(args.lookup_dir)
        for name in models:
            print(f"Downloading pinned model {name}…", flush=True)
            snapshot_for(name, args.model_cache_dir, registry, sha256_file_fn, local_files_only=False)
        write_manifest(
            args.staging_dir / "fantasia-assets-manifest.json",
            args.lookup_dir,
            args.model_cache_dir,
            models,
            registry,
        )
        print("Assets prepared successfully.")
    else:
        require_lookup(args.lookup_dir)
        for name in models:
            print(f"Validating cached model {name}…", flush=True)
            snapshot_for(name, args.model_cache_dir, registry, sha256_file_fn, local_files_only=True)
        print("Assets are ready for an offline compute job.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
