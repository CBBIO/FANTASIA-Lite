import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import fantasia_pipeline
import generate_embeddings


class _Tokenizer:
    calls = []

    @classmethod
    def from_pretrained(cls, model_id, **kwargs):
        cls.calls.append((model_id, kwargs))
        return "tokenizer"


class _ModelInstance:
    def to(self, _device):
        return self

    def eval(self):
        return self


class _Model:
    calls = []

    @classmethod
    def from_pretrained(cls, model_id, **kwargs):
        cls.calls.append((model_id, kwargs))
        return _ModelInstance()


class LiteModelProvenanceTests(unittest.TestCase):
    def setUp(self):
        _Tokenizer.calls.clear()
        _Model.calls.clear()

    def test_registry_revisions_match_pipeline_provenance(self):
        for name, record in fantasia_pipeline.MODEL_PROVENANCE_DEFAULTS.items():
            self.assertEqual(generate_embeddings.MODEL_REGISTRY[name]["hf_id"], record["repository"])
            self.assertEqual(generate_embeddings.MODEL_REGISTRY[name]["revision"], record["revision"])
            self.assertEqual(len(record["revision"]), 40)

    def test_loader_enforces_same_revision_for_tokenizer_and_model(self):
        expected = generate_embeddings.MODEL_REGISTRY["prot_t5"]
        with patch.object(generate_embeddings, "AutoTokenizer", _Tokenizer), patch.object(
            generate_embeddings, "AutoModel", _Model
        ), patch.object(generate_embeddings, "infer_embedding_dim", return_value=1024):
            generate_embeddings.load_model_components("prot_t5", "cpu")
        self.assertEqual(_Tokenizer.calls[0][0], expected["hf_id"])
        self.assertEqual(_Tokenizer.calls[0][1]["revision"], expected["revision"])
        self.assertEqual(_Model.calls[0][0], expected["hf_id"])
        self.assertEqual(_Model.calls[0][1]["revision"], expected["revision"])

    def test_run_writes_model_provenance_before_other_stages(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "input.fa"
            fasta.write_text(">p1\nACDE\n")
            args = fantasia_pipeline.parse_args(
                [
                    str(fasta),
                    "--venv-dir", str(root / "venv"),
                    "--results-csv", str(root / "results.csv"),
                    "--embeddings-npz", str(root / "embeddings.npz"),
                    "--config-yaml", str(root / "config.yaml"),
                    "--failure-report", str(root / "failures.csv"),
                    "--run-log", str(root / "run.log"),
                    "--run-metadata-yaml", str(root / "run_metadata.yaml"),
                    "--model-provenance-yaml", str(root / "model_provenance.yaml"),
                    "--topgo-dir", str(root / "topgo"),
                ]
            )
            config = fantasia_pipeline.PipelineConfig.from_args(args)
            with patch.object(
                fantasia_pipeline, "create_virtualenv", side_effect=RuntimeError("stop")
            ):
                with self.assertRaises(RuntimeError):
                    fantasia_pipeline.run_pipeline(config, argv=[str(fasta)])
            provenance = (root / "model_provenance.yaml").read_text()
            self.assertIn("Rostlab/prot_t5_xl_uniref50", provenance)
            self.assertIn("973be27c52ee6474de9c945952a8008aeb2a1a73", provenance)
            self.assertIn("revision_enforced: true", provenance)


if __name__ == "__main__":
    unittest.main()
