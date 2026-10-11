import subprocess
import tempfile
import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = ROOT / "scripts" / "fantasia-hpc"


class HpcLauncherTests(unittest.TestCase):
    def test_dry_run_uses_configured_mounts_and_offline_environment(self):
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            folders = {name: root / name for name in ("lookup", "models", "inputs", "outputs", "scratch", "staging")}
            for folder in folders.values():
                folder.mkdir()
            image = root / "fantasia.sif"
            image.touch()
            fasta = folders["inputs"] / "query.fa"
            fasta.write_text(">query\nACDE\n", encoding="utf-8")
            config = {
                "runtime": {
                    "image": str(image),
                    "apptainer": "apptainer",
                    "apptainer_cache": str(root / "apptainer-cache"),
                    "use_nv": True,
                },
                "storage": {
                    name: {
                        "host": str(folder),
                        "container": f"/site/{name}",
                        "mode": "ro" if name in {"lookup", "models", "inputs"} else "rw",
                    }
                    for name, folder in folders.items()
                },
            }
            config_path = root / "site.yaml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            result = subprocess.run(
                [
                    str(LAUNCHER),
                    "--site-config", str(config_path),
                    "--dry-run",
                    "run",
                    "--input", str(fasta),
                    "--run-name", "job-42",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertIn("--nv", result.stdout)
            self.assertIn(f"{folders['lookup']}:/site/lookup:ro", result.stdout)
            self.assertIn(f"{folders['outputs']}/job-42:/site/outputs:rw", result.stdout)
            self.assertIn("HF_HUB_OFFLINE=1", result.stdout)
            self.assertIn("HF_HUB_CACHE=/site/models", result.stdout)
            self.assertIn("/site/inputs/query.fa", result.stdout)


if __name__ == "__main__":
    unittest.main()
