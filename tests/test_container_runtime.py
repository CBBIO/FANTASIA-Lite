"""Regression checks for the immutable container's legacy pip adapter."""
import importlib.util
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "runtime_pip", Path(__file__).resolve().parents[1] / "containers/runtime_pip.py"
)
runtime = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runtime)


class RuntimeAdapterTests(unittest.TestCase):
    def test_legacy_installer_requests_are_checked(self):
        with patch.object(sys, "argv", ["pip", "install", "--no-deps", "esm==3.2.3"]), patch.object(
            runtime, "version", return_value="3.2.3"
        ) as installed:
            runtime.main()
            installed.assert_called_once_with("esm")

    def test_incompatible_version_fails(self):
        with patch.object(sys, "argv", ["pip", "install", "transformers<4.48.2"]), patch.object(
            runtime, "version", return_value="5.0.0"
        ), self.assertRaises(RuntimeError):
            runtime.main()

    def test_missing_package_fails(self):
        with patch.object(sys, "argv", ["pip", "install", "torch"]), patch.object(
            runtime, "version", side_effect=runtime.PackageNotFoundError("torch")
        ), self.assertRaises(runtime.PackageNotFoundError):
            runtime.main()

    def test_unknown_operation_fails(self):
        with patch.object(sys, "argv", ["pip", "uninstall", "torch"]), self.assertRaises(RuntimeError):
            runtime.main()
