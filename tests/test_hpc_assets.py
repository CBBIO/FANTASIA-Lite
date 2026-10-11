import hashlib
import importlib.util
import io
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "hpc_assets.py"
SPEC = importlib.util.spec_from_file_location("hpc_assets", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
hpc_assets = importlib.util.module_from_spec(SPEC)
sys.modules["hpc_assets"] = hpc_assets
SPEC.loader.exec_module(hpc_assets)


class HpcAssetsTests(unittest.TestCase):
    def test_download_lookup_verifies_checksum_and_installs_required_files(self):
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            source = root / "bundle.tgz"
            files = {
                "data/lookup/lookup_table.npz": b"npz",
                "data/lookup/annotations.json": b"{}",
                "data/lookup/accessions.json": b"{}",
            }
            with tarfile.open(source, "w:gz") as archive:
                for name, content in files.items():
                    metadata = tarfile.TarInfo(name)
                    metadata.size = len(content)
                    archive.addfile(metadata, io.BytesIO(content))
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            destination = root / "lookup"
            hpc_assets.download_lookup(
                source.as_uri(), digest, root / "staging" / "archive.tgz", destination
            )
            self.assertEqual((destination / "lookup_table.npz").read_bytes(), b"npz")
            hpc_assets.require_lookup(destination)

    def test_download_lookup_rejects_checksum_mismatch(self):
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            source = root / "bundle.tgz"
            with tarfile.open(source, "w:gz"):
                pass
            with self.assertRaises(RuntimeError):
                hpc_assets.download_lookup(
                    source.as_uri(), "0" * 64, root / "staging" / "archive.tgz", root / "lookup"
                )


if __name__ == "__main__":
    unittest.main()
