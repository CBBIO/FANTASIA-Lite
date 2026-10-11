import importlib.util
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "hpc_site.py"
SPEC = importlib.util.spec_from_file_location("hpc_site", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
hpc_site = importlib.util.module_from_spec(SPEC)
import sys

sys.modules["hpc_site"] = hpc_site
SPEC.loader.exec_module(hpc_site)


def config_for(root: Path) -> dict:
    paths = {
        "lookup": root / "lookup",
        "models": root / "models",
        "inputs": root / "inputs",
        "outputs": root / "outputs",
        "scratch": root / "scratch",
        "staging": root / "staging",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    image = root / "fantasia.sif"
    image.touch()
    return {
        "runtime": {
            "image": str(image),
            "apptainer_cache": str(root / "apptainer-cache"),
            "use_nv": True,
        },
        "storage": {
            name: {
                "host": str(path),
                "container": f"/custom/{name}",
                "mode": "ro" if name in {"lookup", "models", "inputs"} else "rw",
            }
            for name, path in paths.items()
        },
        "assets": {
            "lookup": {
                "url": "https://example.invalid/lookup.tgz",
                "sha256": "a" * 64,
            }
        },
    }


class HpcSiteConfigTests(unittest.TestCase):
    def write_config(self, root: Path, raw: dict) -> Path:
        path = root / "site.yaml"
        path.write_text(yaml.safe_dump(raw), encoding="utf-8")
        return path

    def test_config_keeps_all_custom_container_paths(self):
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            site = hpc_site.load_site_config(self.write_config(root, config_for(root)))
            self.assertEqual(str(site.mount("lookup").container), "/custom/lookup")
            self.assertEqual(str(site.mount("models").container), "/custom/models")

    def test_expands_environment_in_all_host_paths(self):
        with tempfile.TemporaryDirectory(prefix="fantasia site ") as raw_root:
            root = Path(raw_root)
            raw = config_for(root)
            raw["runtime"]["image"] = "${FANTASIA_SITE_ROOT}/fantasia.sif"
            raw["runtime"]["apptainer_cache"] = "$FANTASIA_SITE_ROOT/apptainer-cache"
            for name in raw["storage"]:
                raw["storage"][name]["host"] = "${FANTASIA_SITE_ROOT}/" + name
            with patch.dict(os.environ, {"FANTASIA_SITE_ROOT": str(root)}):
                site = hpc_site.load_site_config(self.write_config(root, raw))
            self.assertEqual(site.image, root / "fantasia.sif")
            self.assertEqual(site.apptainer_cache, root / "apptainer-cache")
            for name in raw["storage"]:
                self.assertEqual(site.mount(name).host, root / name)
                self.assertEqual(str(site.mount(name).container), f"/custom/{name}")
            fasta = root / "inputs" / "query.fa"
            self.assertEqual(str(hpc_site.mapped_input(site, fasta)), "/custom/inputs/query.fa")

    def test_undefined_variable_reports_field_and_name(self):
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            for field in ("image", "apptainer_cache", *hpc_site.REQUIRED_MOUNTS):
                with self.subTest(field=field):
                    raw = config_for(root)
                    if field in {"image", "apptainer_cache"}:
                        raw["runtime"][field] = "/shared/${FANTASIA_MISSING}/asset"
                        location = "runtime." + field
                    else:
                        raw["storage"][field]["host"] = "$FANTASIA_MISSING/asset"
                        location = f"storage.{field}.host"
                    with patch.dict(os.environ, {}, clear=True):
                        with self.assertRaises(hpc_site.SiteConfigError) as caught:
                            hpc_site.load_site_config(self.write_config(root, raw))
                    self.assertIn(location, str(caught.exception))
                    self.assertIn("undefined environment variable 'FANTASIA_MISSING'", str(caught.exception))

    def test_tilde_expands_after_environment_substitution(self):
        with patch.dict(os.environ, {"FANTASIA_SITE_ROOT": "~/fantasia"}):
            path = hpc_site._absolute_path("${FANTASIA_SITE_ROOT}/models", "storage.models.host")
        self.assertEqual(path, Path.home() / "fantasia/models")
        self.assertEqual(hpc_site._absolute_path("~/models", "storage.models.host"), Path.home() / "models")

    def test_relative_expanded_path_is_rejected(self):
        with patch.dict(os.environ, {"FANTASIA_SITE_ROOT": "relative"}):
            with self.assertRaisesRegex(hpc_site.SiteConfigError, "must be absolute"):
                hpc_site._absolute_path("$FANTASIA_SITE_ROOT/models", "storage.models.host")

    def test_shell_expressions_are_not_supported(self):
        for value in ("/shared/${ROOT:-fallback}", "/shared/$(whoami)", "/shared/${ROOT"):
            with self.subTest(value=value), self.assertRaisesRegex(hpc_site.SiteConfigError, "invalid variable syntax"):
                hpc_site._absolute_path(value, "storage.models.host")

    def test_expansion_is_single_pass_and_supports_literal_dollars(self):
        with patch.dict(os.environ, {"FANTASIA_SITE_ROOT": "/shared/$literal"}):
            self.assertEqual(
                hpc_site._absolute_path("${FANTASIA_SITE_ROOT}/$$weights", "storage.models.host"),
                Path("/shared/$literal/$weights"),
            )

    def test_input_is_mapped_under_configured_mount(self):
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            raw = config_for(root)
            fasta = Path(raw["storage"]["inputs"]["host"]) / "nested" / "test.fa"
            fasta.parent.mkdir()
            fasta.write_text(">protein\nACDE\n", encoding="utf-8")
            site = hpc_site.load_site_config(self.write_config(root, raw))
            self.assertEqual(str(hpc_site.mapped_input(site, fasta)), "/custom/inputs/nested/test.fa")

    def test_compute_mounts_keep_assets_read_only_and_isolate_writes(self):
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            raw = config_for(root)
            site = hpc_site.load_site_config(self.write_config(root, raw))
            mounts = {mount.name: mount for mount in hpc_site.run_mounts(site, "12345")}
            self.assertEqual(mounts["lookup"].mode, "ro")
            self.assertEqual(mounts["models"].mode, "ro")
            self.assertEqual(mounts["outputs"].mode, "rw")
            self.assertTrue(str(mounts["outputs"].host).endswith("outputs/12345"))

    def test_rejects_compute_assets_configured_writable(self):
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            raw = config_for(root)
            raw["storage"]["models"]["mode"] = "rw"
            with self.assertRaises(hpc_site.SiteConfigError):
                hpc_site.load_site_config(self.write_config(root, raw))

    def test_requires_a_real_lookup_checksum_shape(self):
        with tempfile.TemporaryDirectory() as raw_root:
            root = Path(raw_root)
            raw = config_for(root)
            raw["assets"]["lookup"]["sha256"] = "short"
            with self.assertRaises(hpc_site.SiteConfigError):
                hpc_site.load_site_config(self.write_config(root, raw))


if __name__ == "__main__":
    unittest.main()
