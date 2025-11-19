"""
fantasia_no_db.py
===================

This module provides a lightweight re‑implementation of the key lookup
functionality from the FANTASIA pipeline without relying on any
relational database.  In the official FANTASIA project, reference
embeddings and GO annotations are stored in a PostgreSQL database
augmented with the `pgvector` extension.  Query embeddings are
compared against those reference vectors in memory and the nearest
neighbours are used to transfer Gene Ontology (GO) terms.

The implementation below avoids any database dependency by loading
reference embeddings and their associated annotations from local
files (for example a NumPy ``.npz`` archive for the lookup table and
a JSON or CSV file for annotations).  It then performs an in‑memory
similarity search using NumPy/ SciPy and writes out a simple CSV
containing the transferred GO terms for each query.

This script can be run either as a library or directly from the
command line.  When executed as a script it accepts a YAML
configuration file specifying the locations of the lookup table,
annotations and query embeddings along with a handful of
parameters.  See the README for further details.
"""

from __future__ import annotations

import json
import os
import csv
import re
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
try:
    import h5py  # type: ignore
except ImportError:
    # h5py is an optional dependency.  It is only required when reading
    # HDF5 query embedding files.  If it is not installed the user
    # should either install ``h5py`` or convert their HDF5 files to
    # another supported format.
    h5py = None  # type: ignore
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


def extract_accession_id(fasta_header: str) -> str:
    """
    Extract just the accession ID from a FASTA header.
    
    For headers like: 'sp|P07648|RECC_ECOLI RecBCD enzyme subunit RecC OS=...'
    Returns just: 'sp|P07648|RECC_ECOLI'
    
    Parameters
    ----------
    fasta_header : str
        The full FASTA header line (without the '>' prefix)
    
    Returns
    -------
    str
        The cleaned accession ID
    """
    # Split on first space and take the first part 
    accession = fasta_header.split()[0] if fasta_header.strip() else fasta_header
    return accession.strip()


class EmbeddingLookUpLocal:
    """Perform embedding based GO transfer without a relational database.

    The class expects a configuration dictionary with the following
    keys:

    ``lookup_table_path``: path to a ``.npz`` file containing the
        reference embeddings.  Each key in the archive corresponds to a
        different model and should end with ``_ids`` or ``_embeddings``.
        For example, ``ESM_ids`` and ``ESM_embeddings`` would
        correspond to the same model.  If the file contains only
        ``ids`` and ``embeddings`` then those will be used as the sole
        lookup table.

    ``annotations_path``: path to a JSON or CSV file mapping reference
        sequence identifiers to lists of GO term records.  Each record
        should at least contain a ``go_id`` field; optional fields
        include ``category``, ``go_description`` and ``evidence_code``.

    ``embeddings_path``: path to an HDF5 file produced by the
        embedding step of FANTASIA (or any other tool).  The file
        should contain a group for each accession with datasets
        ``sequence`` and ``type_<id>/embedding`` similar to the
        official pipeline.

    ``embedding``: a dictionary controlling distance calculations.  It
        may contain ``distance_metric`` (either ``euclidean`` or
        ``cosine``) and a ``models`` subsection describing which
        models to use.  Each entry in the ``models`` dict should
        provide a ``distance_threshold`` and optionally ``batch_size``.

    ``limit_per_entry``: how many neighbours to retain per query.

    ``results_path``: where to write the resulting CSV file.  If not
        supplied it defaults to ``results.csv`` in the current working
        directory.

    This class operates entirely in memory and therefore may require
    substantial RAM depending on the size of your lookup table.
    """

    def __init__(self, conf: Dict):
        self.conf = conf
        self.lookup_tables: Dict[str, Dict[str, np.ndarray]] = {}
        self.go_annotations: Dict[str, List[Dict]] = {}

        # default parameters
        self.limit_per_entry = int(conf.get("limit_per_entry", 1))
        self.distance_metric = conf.get("embedding", {}).get(
            "distance_metric", "euclidean"
        ).lower()
        if self.distance_metric not in {"euclidean", "cosine"}:
            raise ValueError(
                f"Unsupported distance metric '{self.distance_metric}'. "
                "Choose 'euclidean' or 'cosine'."
            )

        # output paths
        self.results_path = conf.get(
            "results_path", os.path.join(os.getcwd(), "results.csv")
        )
        self.raw_results_path: Optional[str] = conf.get("raw_results_path")

        # optional path to an accession mapping.  If provided, the file should
        # contain a JSON object mapping reference sequence identifiers to
        # UniProt accessions.  When loaded, these values will be added to
        # the results under the 'uniprot_accession' column.
        self.accession_path = conf.get("accession_path")
        self.accession_map: Dict[str, str] | None = None

    def load_lookup_table(self) -> None:
        """Load reference embeddings from a NumPy ``.npz`` archive.

        The archive may contain multiple models.  Keys ending with
        ``_ids`` must have a matching key ending with ``_embeddings``.
        If only ``ids`` and ``embeddings`` are present, they are used
        under a default model name ``default``.  The loaded tables are
        stored in ``self.lookup_tables``.
        """
        path = self.conf.get("lookup_table_path")
        if not path:
            raise ValueError(
                "Configuration must include 'lookup_table_path' pointing to a .npz file."
            )
        if not os.path.exists(path):
            raise FileNotFoundError(f"Lookup table file not found: {path}")

        data = np.load(path, allow_pickle=True)
        keys = set(data.files)
        # Determine model prefixes
        model_prefixes = set()
        for k in keys:
            if k.endswith("_ids"):
                model_prefixes.add(k[: -len("_ids")])
            elif k.endswith("_embeddings"):
                model_prefixes.add(k[: -len("_embeddings")])
        # If no suffixes, default single model
        if not model_prefixes and "ids" in keys and "embeddings" in keys:
            model_prefixes = {"default"}

        metric = self.distance_metric
        for prefix in model_prefixes:
            ids_key = f"{prefix}_ids" if f"{prefix}_ids" in keys else "ids"
            embeds_key = (
                f"{prefix}_embeddings" if f"{prefix}_embeddings" in keys else "embeddings"
            )
            ids = data[ids_key]
            embeddings = np.asarray(data[embeds_key], dtype=np.float32, order="C")
            if embeddings.ndim != 2:
                raise ValueError(
                    f"Embeddings for model '{prefix}' must be a 2D array; got shape {embeddings.shape}."
                )
            if len(ids) != len(embeddings):
                raise ValueError(
                    f"Mismatched lengths for ids and embeddings in model '{prefix}'."
                )
            table: Dict[str, np.ndarray] = {"ids": ids, "embeddings": embeddings}

            if metric == "cosine":
                norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
                np.divide(
                    embeddings,
                    np.maximum(norms, 1e-12),
                    out=embeddings,
                )
            else:  # euclidean
                squared_norms = np.einsum("ij,ij->i", embeddings, embeddings, dtype=np.float32)
                table["squared_norms"] = squared_norms

            self.lookup_tables[prefix] = table

    def load_annotations(self) -> None:
        """Load GO annotations mapping from a JSON or CSV file.

        The annotations file should map reference identifiers to lists of
        annotation records.  For JSON the top level is expected to be
        a dictionary keyed by the identifier with a list of records as
        the value.  For CSV the file must contain at least a
        ``reference_id`` column and a ``go_id`` column; any additional
        columns will be preserved in the annotation record.
        """
        path = self.conf.get("annotations_path")
        if not path:
            raise ValueError(
                "Configuration must include 'annotations_path' pointing to a JSON or CSV file."
            )
        if not os.path.exists(path):
            raise FileNotFoundError(f"Annotations file not found: {path}")

        ext = os.path.splitext(path)[1].lower()
        annotations: Dict[str, List[Dict]] = {}
        if ext == ".json":
            with open(path, "r", encoding="utf-8") as f:
                raw = json.load(f)
            if not isinstance(raw, dict):
                raise ValueError(
                    "JSON annotations file must contain an object at the top level."
                )
            # ensure values are lists
            for ref_id, recs in raw.items():
                if isinstance(recs, dict):
                    annotations.setdefault(str(ref_id), []).append(recs)
                elif isinstance(recs, list):
                    annotations[str(ref_id)] = list(recs)
                else:
                    raise ValueError(
                        "Annotation entries must be objects or lists of objects."
                    )
        elif ext in {".csv", ".tsv"}:
            delimiter = "," if ext == ".csv" else "\t"
            with open(path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter=delimiter)
                if "reference_id" not in reader.fieldnames or "go_id" not in reader.fieldnames:
                    raise ValueError(
                        "CSV/TSV annotations require at least 'reference_id' and 'go_id' columns."
                    )
                for row in reader:
                    ref_id = str(row.pop("reference_id"))
                    annotations.setdefault(ref_id, []).append(row)
        else:
            raise ValueError(
                f"Unsupported annotation file extension '{ext}'. Use JSON or CSV/TSV."
            )
        self.go_annotations = annotations

    def _compute_distances(
        self, query_embedding: np.ndarray, ref_table: Dict[str, np.ndarray]
    ) -> np.ndarray:
        """Compute distances between a single query and all reference embeddings.

        Depending on ``self.distance_metric`` either Euclidean distances (smaller is
        more similar) or cosine distances (1 - cosine similarity) are returned.  The
        result is a 1D array of length equal to the number of reference vectors.
        """
        ref_embeddings = ref_table["embeddings"]
        query = np.asarray(query_embedding, dtype=np.float32)
        if self.distance_metric == "euclidean":
            squared_norms = ref_table.get("squared_norms")
            if squared_norms is None:
                squared_norms = np.einsum(
                    "ij,ij->i", ref_embeddings, ref_embeddings, dtype=np.float32
                )
                ref_table["squared_norms"] = squared_norms

            query_sq = float(np.dot(query, query))
            dot = ref_embeddings @ query
            distances = squared_norms + query_sq
            distances -= 2.0 * dot
            np.maximum(distances, 0.0, out=distances)
            np.sqrt(distances, out=distances)
            return distances

        # cosine metric
        query_norm = float(np.linalg.norm(query))
        if query_norm == 0.0:
            return np.ones(ref_embeddings.shape[0], dtype=np.float32)
        inv_norm = 1.0 / query_norm
        query_unit = query * inv_norm
        similarities = ref_embeddings @ query_unit
        distances = 1.0 - similarities
        return distances

    def start(self) -> None:
        """Run the lookup pipeline.

        This will load the lookup tables and annotations, iterate over
        the query embeddings stored in the provided HDF5 file and
        compute nearest neighbours.  A CSV file is written to
        ``self.results_path`` containing the transferred annotations.
        """
        self.load_lookup_table()
        self.load_annotations()

        # Load optional accession mapping if provided.  The mapping file is
        # expected to be a JSON dictionary keyed by reference sequence id
        # (as strings) with values equal to UniProt accessions.  If the path
        # is not set or the file does not exist, no accession information
        # will be included in the results.
        if self.accession_path:
            if os.path.exists(self.accession_path):
                try:
                    with open(self.accession_path, "r", encoding="utf-8") as f:
                        self.accession_map = json.load(f)
                except Exception as e:
                    print(f"Warning: failed to load accession mapping from {self.accession_path}: {e}")
                    self.accession_map = None
            else:
                print(f"Warning: accession mapping file {self.accession_path} not found.")
                self.accession_map = None

        h5_path = self.conf.get("embeddings_path")
        if not h5_path or not os.path.exists(h5_path):
            raise FileNotFoundError(
                f"Query embeddings file not found: {h5_path}. Set 'embeddings_path' in the configuration."
            )

        results: List[Dict] = []
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

        if h5_path.lower().endswith((".h5", ".hdf5")):
            if h5py is None:
                raise ImportError(
                    "h5py is required to read HDF5 files but is not installed. "
                    "Please install the 'h5py' package or convert your embeddings to a NumPy .npz file."
                )
            h5file_ctx = h5py.File(h5_path, "r")
        else:
            # allow query embeddings to be supplied as a .npz file for convenience
            h5file_ctx = None

        if h5file_ctx is not None:
            h5file = h5file_ctx
        else:
            # Load from npz: expect keys 'accessions', 'sequences', and 'embeddings'
            # where 'embeddings' is a nested dict keyed by model_key -> 2D array
            npz = np.load(h5_path, allow_pickle=True)
            h5file = {}
            accessions = npz.get("accessions")
            sequences = npz.get("sequences")
            embeddings_dict = npz.get("embeddings", {})
            if accessions is None or sequences is None:
                raise ValueError(
                    "NPZ query file must contain 'accessions' and 'sequences' arrays."
                )
            for accession, seq in zip(accessions, sequences):
                acc = str(accession)
                grp: Dict[str, Dict | str] = {}
                # Store sequence as a plain Python string.  In NumPy >=2.0 ``np.string_``
                # was removed, so avoid using it.  The downstream code handles
                # either strings or bytes.
                grp["sequence"] = seq
                # For each model key, store embedding under 'type_<key>/embedding'
                for model_key, arr in embeddings_dict.item().items():
                    g: Dict[str, np.ndarray] = {}
                    g["embedding"] = arr[accessions.tolist().index(accession)]
                    grp[f"type_{model_key}"] = g
                h5file[acc] = grp

        # iterate over entries
        for accession, group in h5file.items():
            # Clean the accession to extract only the ID part (e.g., sp|P07648|RECC_ECOLI)
            clean_accession = extract_accession_id(accession)
            # group may be a dict when reading from npz
            if "sequence" not in group:
                continue
            # decode sequences whether stored as bytes or string
            # Extract the amino acid sequence; supports HDF5 datasets and plain values.
            seq_dset = group["sequence"]
            # Try to obtain the underlying value.  HDF5 datasets require indexing
            try:
                val = seq_dset[()]
            except Exception:
                val = seq_dset
            # ``val`` may be a NumPy scalar, bytes, or python string
            if isinstance(val, (bytes, bytearray, np.bytes_)):
                sequence = val.decode("utf-8")
            else:
                sequence = str(val)
            for item_name, item_group in group.items():
                if not str(item_name).startswith("type_"):
                    continue
                # item_group could be a dict when reading from npz
                if isinstance(item_group, dict):
                    if "embedding" not in item_group:
                        continue
                    query_emb = np.asarray(item_group["embedding"], dtype=float)
                else:
                    if "embedding" not in item_group:
                        continue
                    query_emb = item_group["embedding"][:]
                model_id_str = str(item_name).replace("type_", "")
                model_key = None
                if model_id_str in self.lookup_tables:
                    model_key = model_id_str
                elif model_id_str.isdigit():
                    # numeric IDs might correspond to model names if specified in the config
                    models_conf = self.conf.get("embedding", {}).get("models", {})
                    for name in models_conf:
                        info = models_conf[name]
                        # check if config entry contains an id field equal to the numeric string
                        if str(info.get("id")) == model_id_str and name in self.lookup_tables:
                            model_key = name
                            break
                elif len(self.lookup_tables) == 1:
                    model_key = next(iter(self.lookup_tables.keys()))
                else:
                    # attempt to match by name ignoring punctuation
                    models_conf = self.conf.get("embedding", {}).get("models", {})
                    for name in models_conf:
                        if name.lower().replace("-", "").replace("_", "") == model_id_str.lower().replace("-", "").replace("_", ""):
                            if name in self.lookup_tables:
                                model_key = name
                                break
                if model_key is None:
                    continue
                ref_table = self.lookup_tables[model_key]
                ref_ids = ref_table["ids"]
                dists = self._compute_distances(query_emb, ref_table)
                if self.limit_per_entry > 0:
                    k = min(self.limit_per_entry, len(dists))
                    idx = np.argpartition(dists, k - 1)[:k]
                    top_idx = idx[np.argsort(dists[idx])]
                else:
                    top_idx = np.argsort(dists)
                for ref_index in top_idx:
                    ref_id = str(ref_ids[ref_index])
                    distance = float(dists[ref_index])
                    annotations = self.go_annotations.get(ref_id, [])
                    # Look up UniProt accession if a mapping is available
                    uniprot = None
                    if self.accession_map is not None:
                        uniprot = self.accession_map.get(ref_id)
                    if not annotations:
                        results.append(
                            {
                                "query_accession": clean_accession,
                                "query_sequence": sequence,
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
                    else:
                        for anno in annotations:
                            results.append(
                                {
                                    "query_accession": clean_accession,
                                    "query_sequence": sequence,
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

        # Close the h5file context if it was opened
        if h5file_ctx is not None:
            h5file_ctx.close()

        # Write out results
        if results:
            df = pd.DataFrame(results)
            # Compute reliability index per row based on the configured distance metric.
            metric = self.distance_metric
            if metric == "cosine":
                df["reliability_index"] = 1.0 - df["distance"]
            elif metric == "euclidean":
                df["reliability_index"] = 0.5 / (0.5 + df["distance"])
            else:
                df["reliability_index"] = np.nan
            # Clamp reliability to the [0, 1] interval where applicable.
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
                    if not text:
                        continue
                    if text not in combined:
                        combined.append(text)
                if not combined:
                    return None
                return "|".join(combined)

            df = (
                df.groupby(per_reference_group, dropna=False, sort=False)
                .agg({"evidence_code": _combine_evidence, **first_columns})
                .reset_index()
            )
            df.rename(columns={"evidence_code": "evidence_codes"}, inplace=True)
            # Select the best reference per model/query/GO/category.
            selection_keys = ["model_key", "query_accession", "go_id", "category"]
            df["evidence_count"] = (
                df["evidence_codes"]
                .fillna("")
                .apply(lambda x: 0 if not x else len(x.split("|")))
            )
            df["_row_index"] = np.arange(len(df))
            df = df.sort_values(
                by=["distance", "evidence_count", "_row_index"],
                ascending=[True, False, True],
                kind="mergesort",
            )
            df = (
                df.groupby(selection_keys, dropna=False, sort=False)
                .first()
                .reset_index()
            )
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
        else:
            # Ensure an empty file is created
            with open(self.results_path, "w", encoding="utf-8") as f:
                f.write(
                    "query_accession,reference_id,model_key,distance,reliability_index,distance_metric,uniprot_accession,go_id,category,go_description,evidence_codes\n"
                )
            if self.raw_results_path:
                empty_df = pd.DataFrame(columns=raw_columns)
                empty_df.to_csv(self.raw_results_path, index=False)


def _load_yaml_config(path: str) -> Dict:
    import yaml  # imported lazily to avoid dependency for library use
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main() -> None:
    """Entry point for command line execution."""
    import argparse

    parser = argparse.ArgumentParser(description="Run FANTASIA lookup without a database.")
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to a YAML configuration file describing the run.",
    )
    args = parser.parse_args()
    conf = _load_yaml_config(args.config)
    lookup = EmbeddingLookUpLocal(conf)
    lookup.start()


if __name__ == "__main__":
    main()
