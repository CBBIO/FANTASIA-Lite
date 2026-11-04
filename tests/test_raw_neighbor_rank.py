import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

import fantasia_no_db

try:
    import torch

    import fantasia_no_db_gpu
except ImportError:  # pragma: no cover - torch is optional for this test
    torch = None
    fantasia_no_db_gpu = None


# Three donors at increasing distance from the query.  Donor "A" carries two GO
# terms and donor "B" none, so GO-row hit_rank (A=1,2 / B=3 / C=4) diverges
# from the true neighbour position (A=1 / B=2 / C=3).
REF_IDS = np.array(["A", "B", "C"], dtype=object)
REF_EMBEDDINGS = np.array([[1.0, 0.0], [2.0, 0.0], [4.0, 0.0]], dtype=np.float32)
ANNOTATIONS = {
    "A": [
        {"go_id": "GO:0000001", "category": "F", "go_description": "a1", "evidence_code": "IDA"},
        {"go_id": "GO:0000002", "category": "P", "go_description": "a2", "evidence_code": "IDA"},
    ],
    "C": [
        {"go_id": "GO:0000003", "category": "C", "go_description": "c1", "evidence_code": "IDA"},
    ],
}
QUERY = [
    {
        "query_accession": "q1",
        "query_sequence": "MSEQ",
        "models": {"Prot-T5": np.array([0.0, 0.0], dtype=np.float32)},
    }
]
EXPECTED_NEIGHBOR_RANK = {"A": 1, "B": 2, "C": 3}


def _configure(lookup, tmpdir: Path, table) -> None:
    lookup.lookup_tables = {"Prot-T5": table}
    lookup.go_annotations = ANNOTATIONS
    lookup.accession_map = None
    lookup.limit_per_entry = 3
    lookup.lookup_batch_size = 8
    lookup.distance_metric = "euclidean"
    lookup.results_path = str(tmpdir / "results.csv")
    lookup.raw_results_path = str(tmpdir / "k.3.results.csv")


class RawNeighborRankTests(unittest.TestCase):
    def _assert_outputs(self, tmpdir: Path) -> None:
        raw = pd.read_csv(tmpdir / "k.3.results.csv")
        self.assertEqual(list(raw.columns[:3]), ["query_accession", "hit_rank", "neighbor_rank"])
        per_donor = raw.groupby("reference_id")["neighbor_rank"].unique()
        for ref_id, expected in EXPECTED_NEIGHBOR_RANK.items():
            self.assertEqual(list(per_donor[ref_id]), [expected], ref_id)
        # hit_rank keeps its GO-row semantics for backward compatibility.
        self.assertEqual(int(raw.loc[raw.reference_id == "C", "hit_rank"].iloc[0]), 4)

        final = pd.read_csv(tmpdir / "results.csv")
        self.assertNotIn("neighbor_rank", final.columns)

    def test_cpu_raw_results_report_donor_neighbor_rank(self):
        lookup = fantasia_no_db.EmbeddingLookUpLocal.__new__(fantasia_no_db.EmbeddingLookUpLocal)
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            _configure(lookup, tmpdir, {"ids": REF_IDS, "embeddings": REF_EMBEDDINGS.copy()})
            lookup._write_results(lookup._build_results(QUERY))
            self._assert_outputs(tmpdir)

    @unittest.skipIf(torch is None, "torch not installed")
    def test_gpu_backend_raw_results_report_donor_neighbor_rank(self):
        cls = fantasia_no_db_gpu.EmbeddingLookUpLocalGPU
        lookup = cls.__new__(cls)
        embeddings = torch.from_numpy(REF_EMBEDDINGS.copy())
        table = {
            "ids": REF_IDS,
            "embeddings": embeddings,
            "squared_norms": (embeddings * embeddings).sum(dim=1),
        }
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            _configure(lookup, tmpdir, table)
            lookup.device = torch.device("cpu")
            lookup._write_results(lookup._build_results(QUERY))
            self._assert_outputs(tmpdir)


if __name__ == "__main__":
    unittest.main()
