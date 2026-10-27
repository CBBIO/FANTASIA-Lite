from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
FAISS_AVAILABLE = all(
    importlib.util.find_spec(module) is not None for module in ("faiss", "pandas")
)
if FAISS_AVAILABLE:
    import pandas as pd

    import fantasia_pipeline
    from fantasia_no_db import EmbeddingLookUpLocal
    from faiss_index_cache import FaissIndexCache
ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(FAISS_AVAILABLE, "faiss-cpu and pandas are not installed")
class FaissLookupTests(unittest.TestCase):
    def make_bundle(self, root: Path) -> tuple[Path, Path, Path, Path]:
        lookup = root / "lookup_table.npz"
        np.savez(
            lookup,
            Toy_ids=np.array(["r1", "r2", "r3"]),
            Toy_embeddings=np.array([[1.0, 0.0], [0.0, 2.0], [-1.0, -1.0]], dtype=np.float32),
        )
        annotations = root / "annotations.json"
        annotations.write_text(
            json.dumps(
                {
                    "r1": [{"go_id": "GO:0001", "category": "F", "go_description": "one", "evidence_code": "EXP"}],
                    "r2": [{"go_id": "GO:0002", "category": "P", "go_description": "two", "evidence_code": "IDA"}],
                    "r3": [{"go_id": "GO:0003", "category": "C", "go_description": "three", "evidence_code": "IMP"}],
                }
            ),
            encoding="utf-8",
        )
        accessions = root / "accessions.json"
        accessions.write_text(json.dumps({"r1": "P1", "r2": "P2", "r3": "P3"}), encoding="utf-8")
        queries = root / "queries.npz"
        np.savez(
            queries,
            accessions=np.array(["q1", "q2"]),
            sequences=np.array(["AAAA", "BBBB"]),
            embeddings=np.array(
                {"Toy": np.array([[0.9, 0.1], [0.1, 1.8]], dtype=np.float32)}, dtype=object
            ),
        )
        return lookup, annotations, accessions, queries

    def run_lookup(
        self, root: Path, backend: str, metric: str, *, cache_dir: Path | None = None
    ) -> pd.DataFrame:
        lookup, annotations, accessions, queries = self.make_bundle(root)
        raw = root / f"{backend}-{metric}-raw.csv"
        conf = {
            "lookup_table_path": str(lookup),
            "annotations_path": str(annotations),
            "accession_path": str(accessions),
            "embeddings_path": str(queries),
            "lookup_backend": backend,
            "faiss_cache_dir": str(cache_dir or root / "faiss"),
            "lookup_batch_size": 1,
            "limit_per_entry": 2,
            "embedding": {"distance_metric": metric},
            "results_path": str(root / f"{backend}-{metric}.csv"),
            "raw_results_path": str(raw),
        }
        EmbeddingLookUpLocal(conf).start()
        return pd.read_csv(raw).sort_values(
            ["query_accession", "hit_rank", "reference_id"], kind="mergesort"
        ).reset_index(drop=True)

    def test_faiss_matches_numpy_for_both_exact_metrics(self):
        for metric in ("cosine", "euclidean"):
            with self.subTest(metric=metric), tempfile.TemporaryDirectory() as raw_root:
                root = Path(raw_root)
                expected = self.run_lookup(root, "numpy", metric)
                actual = self.run_lookup(root, "faiss", metric)
                self.assertEqual(
                    expected[["query_accession", "hit_rank", "reference_id"]].to_dict("records"),
                    actual[["query_accession", "hit_rank", "reference_id"]].to_dict("records"),
                )
                np.testing.assert_allclose(
                    expected["distance"].to_numpy(), actual["distance"].to_numpy(), rtol=1e-5, atol=1e-6
                )

    def test_cache_reuse_invalidation_and_corruption_recovery(self):
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            lookup, _, _, _ = self.make_bundle(root)
            vectors = np.load(lookup)["Toy_embeddings"]
            cache_dir = root / "cache"
            cache = FaissIndexCache(lookup, cache_dir, "cosine")
            first = cache.load_or_build("Toy", vectors)
            index_path = cache_dir / "Toy.cosine.faiss"
            manifest_path = cache_dir / "Toy.cosine.json"
            self.assertTrue(index_path.is_file())
            first_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            first_index_mtime = index_path.stat().st_mtime_ns

            second = FaissIndexCache(lookup, cache_dir, "cosine").load_or_build("Toy", vectors)
            self.assertEqual(first.ntotal, second.ntotal)
            self.assertEqual(first_index_mtime, index_path.stat().st_mtime_ns)

            stat = lookup.stat()
            os.utime(lookup, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1))
            FaissIndexCache(lookup, cache_dir, "cosine").load_or_build("Toy", vectors)
            self.assertNotEqual(
                first_manifest["lookup_mtime_ns"],
                json.loads(manifest_path.read_text(encoding="utf-8"))["lookup_mtime_ns"],
            )

            index_path.write_bytes(b"not a faiss index")
            rebuilt = FaissIndexCache(lookup, cache_dir, "cosine").load_or_build("Toy", vectors)
            self.assertEqual(rebuilt.ntotal, len(vectors))

    def test_prebuild_command_writes_requested_index(self):
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            lookup, _, _, _ = self.make_bundle(root)
            cache_dir = root / "prebuilt"
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "src" / "build_faiss_indexes.py"),
                    "--lookup-npz", str(lookup),
                    "--cache-dir", str(cache_dir),
                    "--models", "Toy",
                    "--distance-metric", "euclidean",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertIn("Ready: Toy", result.stdout)
            self.assertTrue((cache_dir / "Toy.euclidean.faiss").is_file())

    def test_pipeline_config_defaults_cache_and_rejects_gpu_lookup_conflict(self):
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            fasta = root / "query.fa"
            fasta.write_text(">q\nACDE\n", encoding="utf-8")
            lookup = root / "lookup" / "lookup_table.npz"
            args = fantasia_pipeline.parse_args(
                [
                    str(fasta), "--lookup-npz", str(lookup), "--lookup-backend", "faiss",
                    "--config-yaml", str(root / "config.yaml"),
                ]
            )
            config = fantasia_pipeline.PipelineConfig.from_args(args)
            self.assertEqual(config.lookup_backend, "faiss")
            self.assertEqual(config.faiss_cache_dir, lookup.parent.resolve() / "faiss")
            fantasia_pipeline.write_config_yaml(config)
            config_text = config.config_yaml.read_text(encoding="utf-8")
            self.assertIn("lookup_backend: faiss", config_text)
            self.assertIn(f"faiss_cache_dir: {config.faiss_cache_dir}", config_text)
            conflicting = fantasia_pipeline.parse_args(
                [
                    str(fasta), "--lookup-backend", "faiss", "--use-gpu-lookup",
                ]
            )
            with self.assertRaisesRegex(fantasia_pipeline.PipelineError, "cannot be combined"):
                fantasia_pipeline.PipelineConfig.from_args(conflicting)


if __name__ == "__main__":
    unittest.main()
