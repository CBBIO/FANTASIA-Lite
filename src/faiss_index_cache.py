"""Persistent exact FAISS indexes for FANTASIA Lite lookup tables."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any, Iterable

import numpy as np

CACHE_SCHEMA_VERSION = 1


def require_faiss() -> Any:
    """Import FAISS only when the optional lookup backend is selected."""
    try:
        import faiss  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "FAISS lookup was requested, but faiss-cpu is not installed. "
            "Run FANTASIA Lite with --lookup-backend faiss so its managed "
            "environment installs faiss-cpu, or install faiss-cpu yourself."
        ) from exc
    return faiss


def load_lookup_embeddings(
    lookup_path: str | Path, requested_models: Iterable[str] | None = None
) -> dict[str, np.ndarray]:
    """Read selected model matrices from a FANTASIA lookup NPZ archive."""
    requested = set(requested_models) if requested_models is not None else None
    with np.load(lookup_path, allow_pickle=True) as data:
        keys = set(data.files)
        prefixes = {
            key[: -len("_ids")]
            for key in keys
            if key.endswith("_ids")
        }
        prefixes.update(
            key[: -len("_embeddings")]
            for key in keys
            if key.endswith("_embeddings")
        )
        if not prefixes and {"ids", "embeddings"}.issubset(keys):
            prefixes = {"default"}

        unknown = sorted((requested or set()) - prefixes)
        if unknown:
            available = ", ".join(sorted(prefixes)) or "none"
            raise ValueError(
                f"Lookup models not found: {', '.join(unknown)}. Available models: {available}."
            )

        result: dict[str, np.ndarray] = {}
        for prefix in sorted(prefixes):
            if requested is not None and prefix not in requested:
                continue
            key = f"{prefix}_embeddings" if f"{prefix}_embeddings" in keys else "embeddings"
            embeddings = np.asarray(data[key], dtype=np.float32, order="C")
            if embeddings.ndim != 2:
                raise ValueError(
                    f"Embeddings for model '{prefix}' must be 2D; got {embeddings.shape}."
                )
            result[prefix] = embeddings
    return result


class FaissIndexCache:
    """Load or build exact CPU FAISS indexes with source-aware metadata."""

    def __init__(
        self,
        lookup_path: str | Path,
        cache_dir: str | Path,
        metric: str,
        *,
        rebuild: bool = False,
    ) -> None:
        if metric not in {"cosine", "euclidean"}:
            raise ValueError(f"Unsupported FAISS metric: {metric!r}.")
        self.lookup_path = Path(lookup_path).resolve()
        self.cache_dir = Path(cache_dir).resolve()
        self.metric = metric
        self.rebuild = rebuild

    @staticmethod
    def _safe_model_name(model_key: str) -> str:
        return "".join(
            char if char.isalnum() or char in {"-", "_", "."} else "_"
            for char in model_key
        )

    def _paths(self, model_key: str) -> tuple[Path, Path]:
        stem = f"{self._safe_model_name(model_key)}.{self.metric}"
        return self.cache_dir / f"{stem}.faiss", self.cache_dir / f"{stem}.json"

    def _manifest(self, model_key: str, embeddings: np.ndarray, faiss: Any) -> dict[str, Any]:
        source = self.lookup_path.stat()
        return {
            "schema_version": CACHE_SCHEMA_VERSION,
            "lookup_path": str(self.lookup_path),
            "lookup_size": source.st_size,
            "lookup_mtime_ns": source.st_mtime_ns,
            "model_key": model_key,
            "metric": self.metric,
            "rows": int(embeddings.shape[0]),
            "dimension": int(embeddings.shape[1]),
            "faiss_version": str(getattr(faiss, "__version__", "unknown")),
        }

    def _vectors_for_index(self, embeddings: np.ndarray, faiss: Any) -> np.ndarray:
        vectors = np.array(embeddings, dtype=np.float32, order="C", copy=True)
        if self.metric == "cosine":
            faiss.normalize_L2(vectors)
        return vectors

    def load_or_build(self, model_key: str, embeddings: np.ndarray) -> Any:
        """Return a valid flat index, rebuilding it when metadata does not match."""
        faiss = require_faiss()
        if embeddings.ndim != 2:
            raise ValueError(f"Embeddings for model '{model_key}' must be a 2D array.")
        index_path, manifest_path = self._paths(model_key)
        manifest = self._manifest(model_key, embeddings, faiss)
        if not self.rebuild:
            try:
                cached_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                if cached_manifest == manifest:
                    index = faiss.read_index(str(index_path))
                    if index.d == embeddings.shape[1] and index.ntotal == embeddings.shape[0]:
                        return index
            except (OSError, ValueError, json.JSONDecodeError, RuntimeError):
                pass

        self.cache_dir.mkdir(parents=True, exist_ok=True)
        vectors = self._vectors_for_index(embeddings, faiss)
        index = (
            faiss.IndexFlatIP(vectors.shape[1])
            if self.metric == "cosine"
            else faiss.IndexFlatL2(vectors.shape[1])
        )
        index.add(vectors)
        self._publish(index, manifest, index_path, manifest_path, faiss)
        return index

    def _publish(
        self,
        index: Any,
        manifest: dict[str, Any],
        index_path: Path,
        manifest_path: Path,
        faiss: Any,
    ) -> None:
        """Publish a complete index/manifest pair without exposing partial files."""
        fd, raw_index = tempfile.mkstemp(prefix=f".{index_path.name}.", suffix=".tmp", dir=self.cache_dir)
        os.close(fd)
        raw_manifest: str | None = None
        try:
            faiss.write_index(index, raw_index)
            fd, raw_manifest = tempfile.mkstemp(
                prefix=f".{manifest_path.name}.", suffix=".tmp", dir=self.cache_dir
            )
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(manifest, handle, indent=2, sort_keys=True)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(raw_index, index_path)
            os.replace(raw_manifest, manifest_path)
        finally:
            for raw_path in (raw_index, raw_manifest):
                if raw_path:
                    try:
                        os.unlink(raw_path)
                    except FileNotFoundError:
                        pass
