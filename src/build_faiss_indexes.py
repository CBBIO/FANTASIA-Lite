#!/usr/bin/env python3
"""Prebuild exact CPU FAISS indexes for a FANTASIA Lite lookup bundle."""

from __future__ import annotations

import argparse
from pathlib import Path

from faiss_index_cache import FaissIndexCache, load_lookup_embeddings


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lookup-npz", type=Path, required=True, help="Reference lookup_table.npz.")
    parser.add_argument("--cache-dir", type=Path, default=None, help="FAISS cache directory.")
    parser.add_argument(
        "--models",
        nargs="+",
        required=True,
        help="One or more exact lookup-table model keys to index.",
    )
    parser.add_argument(
        "--distance-metric",
        choices=("cosine", "euclidean"),
        default="cosine",
        help="Metric represented by the generated indexes (default: %(default)s).",
    )
    parser.add_argument("--force", action="store_true", help="Rebuild matching cached indexes.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    lookup_npz = args.lookup_npz.resolve()
    if not lookup_npz.is_file():
        raise FileNotFoundError(f"Lookup table not found: {lookup_npz}")
    cache_dir = args.cache_dir.resolve() if args.cache_dir else lookup_npz.parent / "faiss"
    embeddings_by_model = load_lookup_embeddings(lookup_npz, args.models)
    cache = FaissIndexCache(
        lookup_npz, cache_dir, args.distance_metric, rebuild=args.force
    )
    for model_key, embeddings in embeddings_by_model.items():
        cache.load_or_build(model_key, embeddings)
        print(
            f"Ready: {model_key} ({args.distance_metric}, "
            f"{embeddings.shape[0]} x {embeddings.shape[1]}) in {cache_dir}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
