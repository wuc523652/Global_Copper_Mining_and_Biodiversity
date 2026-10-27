"""Regression checks for bundled data, explicit overrides and safe outputs."""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


def load(relative: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


reproduce = load("reproduce.py", "layout_reproduce")
paths = load("analysis/_paths.py", "layout_analysis_paths")
plot_paths = load("plotting/data_paths.py", "layout_plot_paths")


class LocalDataLayoutTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name).resolve()
        self.package = self.base / "code"
        self.data = self.package / "data"
        self.data.mkdir(parents=True)
        for name in ("8Final_Results", "9Figure_Inputs", "10Raw_Inputs"):
            (self.data / name).mkdir()

    def invoke(self, *arguments):
        with patch.object(reproduce, "ROOT", self.package), \
                patch.object(sys, "argv", ["reproduce.py", *arguments]), \
                patch.object(reproduce, "run") as run:
            status = reproduce.main()
        return status, run.call_args_list

    def test_check_defaults_to_data_beside_script(self):
        previous = Path.cwd()
        try:
            os.chdir(self.base)
            status, calls = self.invoke("--mode", "check")
        finally:
            os.chdir(previous)
        self.assertEqual(status, 0)
        self.assertEqual(calls[0].args,
                         ("tools/data_manifest.py", "verify", "--data-root", str(self.data)))

    def test_explicit_data_root_overrides_bundled_data(self):
        external = self.base / "external-data"
        external.mkdir()
        _, calls = self.invoke("--data-root", str(external), "--mode", "check")
        self.assertEqual(calls[0].args[-1], str(external))

    def test_raw_mode_verifies_data_then_runs_real_pipeline(self):
        output = self.base / "raw-results"
        _, calls = self.invoke("--mode", "raw", "--output-root", str(output))
        self.assertEqual(len(calls), 2)
        self.assertEqual(calls[0].args[0], "tools/data_manifest.py")
        self.assertEqual(calls[1].args,
                         ("run_analysis.py", "--all", "--data-root", str(self.data),
                          "--work-root", str(output / "raw_analysis")))

    def test_generated_outputs_cannot_be_inside_data_or_code(self):
        for output in (self.data / "results", self.package / "analysis", self.base):
            with self.subTest(output=output), self.assertRaises(SystemExit) as error:
                self.invoke("--mode", "raw", "--output-root", str(output))
            self.assertEqual(error.exception.code, 2)

    def test_raw_script_default_points_to_actual_raw_inputs(self):
        with patch.object(paths, "PACKAGE_ROOT", self.package), patch.dict(os.environ, {}, clear=True):
            self.assertEqual(Path(paths.data_path("Global_GADM/GC_GADM.shp")),
                             self.data / "10Raw_Inputs/Global_GADM/GC_GADM.shp")

    def test_low_level_output_cannot_overwrite_other_bundled_data(self):
        with patch.object(paths, "PACKAGE_ROOT", self.package):
            with self.assertRaises(ValueError):
                paths.validate_roots(self.data / "10Raw_Inputs", self.data / "8Final_Results")

    def test_raw_script_explicit_environment_overrides_default(self):
        external = self.base / "external-raw"
        with patch.dict(os.environ, {"NC_DATA_ROOT": str(external)}):
            self.assertEqual(Path(paths.data_path("Global_GADM/GC_GADM.shp")),
                             external / "Global_GADM/GC_GADM.shp")

    def test_plotting_defaults_to_bundled_data(self):
        with patch.object(plot_paths, "PACKAGE_ROOT", self.package), patch.dict(os.environ, {}, clear=True):
            self.assertEqual(plot_paths.resolve_data_root(), self.data)

    def test_plotting_cli_then_environment_override_default(self):
        external = self.base / "external-data"
        for folder in ("8Final_Results", "9Figure_Inputs"):
            (external / folder).mkdir(parents=True)
        with patch.dict(os.environ, {"COPPER_DATA_ROOT": str(external)}):
            self.assertEqual(plot_paths.resolve_data_root(), external)
            self.assertEqual(plot_paths.resolve_data_root(self.data), self.data)

    def test_missing_explicit_plot_data_does_not_fall_back(self):
        with patch.object(plot_paths, "PACKAGE_ROOT", self.package), self.assertRaises(ValueError):
            plot_paths.resolve_data_root(self.base / "missing-data")


if __name__ == "__main__":
    unittest.main()
