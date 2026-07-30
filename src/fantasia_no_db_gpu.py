"""
fantasia_no_db_gpu.py
=====================

GPU-enabled companion to ``fantasia_no_db.py`` for FANTASIA-Lite lookup.

This script preserves the same configuration style and output schema as the
CPU implementation, but runs the dense vector search with PyTorch tensors on a
CUDA device when available.  The main speedup comes from scoring batches of
query embeddings against an in-memory reference matrix resident on the GPU.
"""

from __future__ import annotations

import json
import math
import os
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
import pandas as pd

try:
    import h5py  # type: ignore
except ImportError:
    h5py = None  # type: ignore

try:
    import torch
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "PyTorch is required for fantasia_no_db_gpu.py. Install torch with CUDA support."
    ) from exc

try:
    from tqdm import tqdm  # type: ignore
except ImportError:  # pragma: no cover
    class _TqdmShim:
        def __init__(self, iterable, **_kwargs):
            self._iterable = iterable

        def __iter__(self):
            return iter(self._iterable)

        def close(self):
            pass

    def tqdm(iterable, **kwargs):  # type: ignore
        return _TqdmShim(iterable, **kwargs)

from fantasia_no_db import _load_yaml_config, extract_accession_id


class EmbeddingLookUpLocalGPU:
    """Perform database-free GO transfer with GPU-accelerated lookup."""

    def __init__(self, conf: Dict[str, Any]):
        self.conf = conf
        self.go_annotations: Dict[str, List[Dict[str, Any]]] = {}
        self.lookup_tables: Dict[str, Dict[str, Any]] = {}
        self.limit_per_entry = int(conf.get("limit_per_entry", 1))
        self.distance_metric = conf.get("embedding", {}).get(
            "distance_metric", "euclidean"
        ).lower()
        if self.distance_metric not in {"euclidean", "cosine"}:
            raise ValueError(
                f"Unsupported distance metric '{self.distance_metric}'. "
                "Choose 'euclidean' or 'cosine'."
            )

        self.device_name = str(
            conf.get("lookup_device")
            or conf.get("device")
            or conf.get("embedding", {}).get("device")
            or ("cuda" if torch.cuda.is_available() else "cpu")
        )
        self.device = torch.device(self.device_name)
        if self.device.type != "cuda":
            raise ValueError(
                "fantasia_no_db_gpu.py is intended for CUDA lookup. "
                f"Received device={self.device_name!r}."
            )
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA is not available in the current PyTorch installation.")

        self.lookup_batch_size = int(conf.get("lookup_batch_size", 256))
        if self.lookup_batch_size <= 0:
            raise ValueError("lookup_batch_size must be a positive integer.")

        self.results_path = conf.get(
            "results_path", os.path.join(os.getcwd(), "results.csv")
        )
        self.raw_results_path: Optional[str] = conf.get("raw_results_path")
        self.accession_path = conf.get("accession_path")
        self.accession_map: Dict[str, str] | None = None

    def load_lookup_table(self) -> None:
        path = self.conf.get("lookup_table_path")
        if not path:
            raise ValueError(
                "Configuration must include 'lookup_table_path' pointing to a .npz file."
            )
        if not os.path.exists(path):
            raise FileNotFoundError(f"Lookup table file not found: {path}")

        data = np.load(path, allow_pickle=True)
        keys = set(data.files)
        model_prefixes = set()
        for key in keys:
            if key.endswith("_ids"):
                model_prefixes.add(key[: -len("_ids")])
            elif key.endswith("_embeddings"):
                model_prefixes.add(key[: -len("_embeddings")])
        if not model_prefixes and "ids" in keys and "embeddings" in keys:
            model_prefixes = {"default"}

        for prefix in model_prefixes:
            ids_key = f"{prefix}_ids" if f"{prefix}_ids" in keys else "ids"
            embeds_key = (
                f"{prefix}_embeddings" if f"{prefix}_embeddings" in keys else "embeddings"
            )
            ids = data[ids_key]
            embeddings_np = np.asarray(data[embeds_key], dtype=np.float32, order="C")
            if embeddings_np.ndim != 2:
                raise ValueError(
                    f"Embeddings for model '{prefix}' must be a 2D array; got shape {embeddings_np.shape}."
                )
            if len(ids) != len(embeddings_np):
                raise ValueError(
                    f"Mismatched lengths for ids and embeddings in model '{prefix}'."
                )

            embeddings = torch.from_numpy(embeddings_np).to(self.device, non_blocking=True)
            table: Dict[str, Any] = {"ids": ids, "embeddings": embeddings}
            if self.distance_metric == "cosine":
                table["embeddings"] = torch.nn.functional.normalize(embeddings, dim=1)
            else:
                table["squared_norms"] = (embeddings * embeddings).sum(dim=1)
            self.lookup_tables[prefix] = table

    def load_annotations(self) -> None:
        path = self.conf.get("annotations_path")
        if not path:
            raise ValueError(
                "Configuration must include 'annotations_path' pointing to a JSON or CSV file."
            )
        if not os.path.exists(path):
            raise FileNotFoundError(f"Annotations file not found: {path}")

        ext = os.path.splitext(path)[1].lower()
        annotations: Dict[str, List[Dict[str, Any]]] = {}
        if ext == ".json":
            with open(path, "r", encoding="utf-8") as handle:
                raw = json.load(handle)
            if not isinstance(raw, dict):
                raise ValueError(
                    "JSON annotations file must contain an object at the top level."
                )
            for ref_id, recs in raw.items():
                if isinstance(recs, dict):
                    annotations.setdefault(str(ref_id), []).append(recs)
                elif isinstance(recs, list):
                    annotations[str(ref_id)] = list(recs)
                else:
                    raise ValueError(
                        "Annotation entries must be objects or lists of objects."
                    )
        else:
            raise ValueError(
                f"Unsupported annotation file extension '{ext}'. Use JSON for GPU lookup."
            )
        self.go_annotations = annotations

    def _load_accession_map(self) -> None:
        if not self.accession_path:
            return
        if not os.path.exists(self.accession_path):
            print(f"Warning: accession mapping file {self.accession_path} not found.")
            return
        try:
            with open(self.accession_path, "r", encoding="utf-8") as handle:
                self.accession_map = json.load(handle)
        except Exception as exc:
            print(
                f"Warning: failed to load accession mapping from {self.accession_path}: {exc}"
            )
            self.accession_map = None

    def _resolve_model_key(self, model_id_str: str) -> str | None:
        if model_id_str in self.lookup_tables:
            return model_id_str
        if model_id_str.isdigit():
            models_conf = self.conf.get("embedding", {}).get("models", {})
            for name, info in models_conf.items():
                if str(info.get("id")) == model_id_str and name in self.lookup_tables:
                    return name
        if len(self.lookup_tables) == 1:
            return next(iter(self.lookup_tables.keys()))
        models_conf = self.conf.get("embedding", {}).get("models", {})
        normalized = model_id_str.lower().replace("-", "").replace("_", "")
        for name in models_conf:
            if name.lower().replace("-", "").replace("_", "") == normalized:
                if name in self.lookup_tables:
                    return name
        return None

    def _load_query_entries(self) -> List[Dict[str, Any]]:
        h5_path = self.conf.get("embeddings_path")
        if not h5_path or not os.path.exists(h5_path):
            raise FileNotFoundError(
                f"Query embeddings file not found: {h5_path}. Set 'embeddings_path' in the configuration."
            )

        if h5_path.lower().endswith((".h5", ".hdf5")):
            if h5py is None:
                raise ImportError(
                    "h5py is required to read HDF5 files but is not installed."
                )
            return self._load_query_entries_h5(h5_path)
        return self._load_query_entries_npz(h5_path)

    def _load_query_entries_h5(self, path: str) -> List[Dict[str, Any]]:
        entries: List[Dict[str, Any]] = []
        with h5py.File(path, "r") as handle:
            for accession, group in handle.items():
                if "sequence" not in group:
                    continue
                seq_val = group["sequence"][()]
                if isinstance(seq_val, (bytes, bytearray, np.bytes_)):
                    sequence = seq_val.decode("utf-8")
                else:
                    sequence = str(seq_val)
                models: Dict[str, np.ndarray] = {}
                for item_name, item_group in group.items():
                    if not str(item_name).startswith("type_"):
                        continue
                    if "embedding" not in item_group:
                        continue
                    model_id_str = str(item_name).replace("type_", "")
                    model_key = self._resolve_model_key(model_id_str)
                    if model_key is None:
                        continue
                    models[model_key] = np.asarray(item_group["embedding"][:], dtype=np.float32)
                if models:
                    entries.append(
                        {
                            "query_accession": extract_accession_id(accession),
                            "query_sequence": sequence,
                            "models": models,
                        }
                    )
        return entries

    def _load_query_entries_npz(self, path: str) -> List[Dict[str, Any]]:
        npz = np.load(path, allow_pickle=True)
        accessions = npz.get("accessions")
        sequences = npz.get("sequences")
        embeddings_dict = npz.get("embeddings", {})
        if accessions is None or sequences is None:
            raise ValueError(
                "NPZ query file must contain 'accessions' and 'sequences' arrays."
            )

        index_by_model: Dict[str, np.ndarray] = {}
        if hasattr(embeddings_dict, "item"):
            for model_name, arr in embeddings_dict.item().items():
                model_key = self._resolve_model_key(str(model_name))
                if model_key is None:
                    continue
                index_by_model[model_key] = np.asarray(arr, dtype=np.float32)

        entries: List[Dict[str, Any]] = []
        for model_key, arr in index_by_model.items():
            model_accessions = npz.get(f"{model_key}_accessions", accessions)
            model_sequences = npz.get(f"{model_key}_sequences", sequences)
            for idx, (accession, sequence) in enumerate(
                zip(model_accessions, model_sequences, strict=False)
            ):
                if idx >= len(arr):
                    break
                entries.append(
                    {
                        "query_accession": extract_accession_id(str(accession)),
                        "query_sequence": str(sequence),
                        "models": {model_key: arr[idx]},
                    }
                )
        return entries

    def _compute_batch_distances(
        self, query_batch: torch.Tensor, ref_table: Dict[str, Any]
    ) -> torch.Tensor:
        ref_embeddings = ref_table["embeddings"]
        if self.distance_metric == "euclidean":
            query_sq = (query_batch * query_batch).sum(dim=1, keepdim=True)
            dot = query_batch @ ref_embeddings.T
            distances = query_sq + ref_table["squared_norms"].unsqueeze(0) - (2.0 * dot)
            return torch.sqrt(torch.clamp_min(distances, 0.0))

        query_unit = torch.nn.functional.normalize(query_batch, dim=1)
        similarities = query_unit @ ref_embeddings.T
        return torch.clamp_min(1.0 - similarities, 0.0)

    def _topk_indices(self, distances: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        n_refs = distances.shape[1]
        if self.limit_per_entry > 0:
            k = min(self.limit_per_entry, n_refs)
        else:
            k = n_refs
        values, indices = torch.topk(distances, k=k, dim=1, largest=False, sorted=True)
        return values, indices

    def _build_results(self, query_entries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        grouped_queries: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        for entry in query_entries:
            for model_key, embedding in entry["models"].items():
                grouped_queries[model_key].append(
                    {
                        "query_accession": entry["query_accession"],
                        "query_sequence": entry["query_sequence"],
                        "embedding": embedding,
                    }
                )

        results: List[Dict[str, Any]] = []
        for model_key, queries in grouped_queries.items():
            ref_table = self.lookup_tables.get(model_key)
            if ref_table is None:
                continue
            ref_ids = ref_table["ids"]
            batch_count = math.ceil(len(queries) / self.lookup_batch_size)
            for batch_index in tqdm(range(batch_count), desc=f"lookup[{model_key}]"):
                start = batch_index * self.lookup_batch_size
                end = min(start + self.lookup_batch_size, len(queries))
                batch = queries[start:end]
                batch_embeddings = np.stack([item["embedding"] for item in batch]).astype(np.float32)
                query_batch = torch.from_numpy(batch_embeddings).to(self.device, non_blocking=True)
                with torch.no_grad():
                    distances = self._compute_batch_distances(query_batch, ref_table)
                    top_values, top_indices = self._topk_indices(distances)
                top_values_cpu = top_values.detach().cpu().numpy()
                top_indices_cpu = top_indices.detach().cpu().numpy()

                for row_idx, query_item in enumerate(batch):
                    query_accession = query_item["query_accession"]
                    query_sequence = query_item["query_sequence"]
                    for rank_idx, ref_index in enumerate(top_indices_cpu[row_idx]):
                        ref_id = str(ref_ids[int(ref_index)])
                        distance = float(top_values_cpu[row_idx][rank_idx])
                        annotations = self.go_annotations.get(ref_id, [])
                        uniprot = None
                        if self.accession_map is not None:
                            uniprot = self.accession_map.get(ref_id)
                        if not annotations:
                            results.append(
                                {
                                    "query_accession": query_accession,
                                    "query_sequence": query_sequence,
                                    "reference_id": ref_id,
                                    "model_key": model_key,
                                    "distance": distance,
                                    "uniprot_accession": uniprot,
                                    "go_id": None,
                                    "category": None,
                                    "go_description": None,
                                    "evidence_code": None,
                                    "distance_metric": self.distance_metric,
                                }
                            )
                            continue
                        for anno in annotations:
                            results.append(
                                {
                                    "query_accession": query_accession,
                                    "query_sequence": query_sequence,
                                    "reference_id": ref_id,
                                    "model_key": model_key,
                                    "distance": distance,
                                    "uniprot_accession": uniprot,
                                    "go_id": anno.get("go_id"),
                                    "category": anno.get("category"),
                                    "go_description": anno.get("go_description"),
                                    "evidence_code": anno.get("evidence_code"),
                                    "distance_metric": self.distance_metric,
                                }
                            )
        return results

    def _write_results(self, results: List[Dict[str, Any]]) -> None:
        raw_columns = [
            "query_accession",
            "hit_rank",
            "reference_id",
            "model_key",
            "distance",
            "reliability_index",
            "distance_metric",
            "uniprot_accession",
            "go_id",
            "category",
            "go_description",
            "evidence_codes",
        ]
        if results:
            df = pd.DataFrame(results)
            if self.distance_metric == "cosine":
                df["reliability_index"] = 1.0 - df["distance"]
            else:
                df["reliability_index"] = 0.5 / (0.5 + df["distance"])
            df["reliability_index"] = df["reliability_index"].clip(lower=0.0, upper=1.0)

            if self.raw_results_path:
                raw_df = df.copy(deep=True)
                raw_df.rename(columns={"evidence_code": "evidence_codes"}, inplace=True)
                if "query_sequence" in raw_df.columns:
                    raw_df.drop(columns=["query_sequence"], inplace=True)
                raw_df["hit_rank"] = (
                    raw_df.groupby(["model_key", "query_accession"])["distance"]
                    .rank(method="first", ascending=True)
                    .astype(int)
                )
                raw_df.sort_values(
                    by=["model_key", "query_accession", "distance", "reference_id", "go_id"],
                    ascending=[True, True, True, True, True],
                    kind="mergesort",
                    inplace=True,
                )
                ordered_cols = [col for col in raw_columns if col in raw_df.columns]
                remaining_cols = [col for col in raw_df.columns if col not in ordered_cols]
                raw_df = raw_df[ordered_cols + remaining_cols]
                raw_df.to_csv(self.raw_results_path, index=False)

            per_reference_group = [
                "model_key",
                "query_accession",
                "go_id",
                "category",
                "reference_id",
                "uniprot_accession",
            ]
            first_columns = {
                "distance": "first",
                "reliability_index": "first",
                "distance_metric": "first",
                "go_description": "first",
            }

            def _combine_evidence(values: pd.Series) -> str | None:
                combined: List[str] = []
                for value in values:
                    if pd.isna(value):
                        continue
                    text = str(value).strip()
                    if text and text not in combined:
                        combined.append(text)
                return "|".join(combined) if combined else None

            df = (
                df.groupby(per_reference_group, dropna=False, sort=False)
                .agg({"evidence_code": _combine_evidence, **first_columns})
                .reset_index()
            )
            df.rename(columns={"evidence_code": "evidence_codes"}, inplace=True)
            selection_keys = ["model_key", "query_accession", "go_id", "category"]
            df["evidence_count"] = (
                df["evidence_codes"].fillna("").apply(lambda x: 0 if not x else len(x.split("|")))
            )
            df["_row_index"] = np.arange(len(df))
            df = df.sort_values(
                by=["distance", "evidence_count", "_row_index"],
                ascending=[True, False, True],
                kind="mergesort",
            )
            df = df.groupby(selection_keys, dropna=False, sort=False).first().reset_index()
            df.drop(columns=["evidence_count", "_row_index"], inplace=True)
            column_order = [
                "query_accession",
                "reference_id",
                "model_key",
                "distance",
                "reliability_index",
                "distance_metric",
                "uniprot_accession",
                "go_id",
                "category",
                "go_description",
                "evidence_codes",
            ]
            df = df[[col for col in column_order if col in df.columns]]
            df.to_csv(self.results_path, index=False)
            return

        with open(self.results_path, "w", encoding="utf-8") as handle:
            handle.write(
                "query_accession,reference_id,model_key,distance,reliability_index,distance_metric,uniprot_accession,go_id,category,go_description,evidence_codes\n"
            )
        if self.raw_results_path:
            pd.DataFrame(columns=raw_columns).to_csv(self.raw_results_path, index=False)

    def start(self) -> None:
        self.load_lookup_table()
        self.load_annotations()
        self._load_accession_map()
        query_entries = self._load_query_entries()
        results = self._build_results(query_entries)
        self._write_results(results)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Run FANTASIA GPU lookup without a database.")
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to a YAML configuration file describing the run.",
    )
    args = parser.parse_args()
    conf = _load_yaml_config(args.config)
    lookup = EmbeddingLookUpLocalGPU(conf)
    lookup.start()


if __name__ == "__main__":
    main()
