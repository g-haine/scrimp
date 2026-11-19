"""User project starter contract (does not require PETSc, Docker or GetFEM).

Run from the SCRIMP source checkout:
    python -m unittest -v tests.test_project_starter
The Docker runtime deliberately omits project-template/, so these
source-distribution checks are skipped when run inside the container.
"""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "project-template"
CREATOR = ROOT / "scripts/create-scrimp-project.py"


@unittest.skipUnless(TEMPLATE.is_dir() and CREATOR.is_file(),
                     "project template is not bundled in the runtime Docker image")
class StarterProjectTest(unittest.TestCase):
    def test_devcontainer_contract(self):
        obj = json.loads(
            (TEMPLATE / ".devcontainer/devcontainer.json").read_text()
        )
        self.assertIn("SCRIMP_DOCKER_IMAGE", obj["image"])
        self.assertEqual(obj["workspaceFolder"], "/work")
        self.assertIn("type=bind", obj["workspaceMount"])
        self.assertIn("localWorkspaceFolder", obj["workspaceMount"])
        self.assertTrue(obj.get("updateRemoteUserUID"))
        self.assertEqual(
            obj["customizations"]["vscode"]["settings"]
               ["python.defaultInterpreterPath"],
            "/opt/conda/envs/scrimp/bin/python",
        )

    def test_vscode_configs_are_parseable(self):
        for filename in ("launch.json", "settings.json", "tasks.json"):
            cfg = json.loads((TEMPLATE / ".vscode" / filename).read_text())
            self.assertIsInstance(cfg, dict)

    def test_real_example_persists_outputs(self):
        example = (TEMPLATE / "simulations/first_simulation.py").read_text()
        self.assertIn("S.DPHS", example)
        self.assertIn("system.solve()", example)
        self.assertIn("outputs/oscillator.csv", example)
        self.assertNotIn("plt.show()", example)

    def test_generator_creates_isolated_project_and_never_overwrites(self):
        with tempfile.TemporaryDirectory() as temporary:
            dest = Path(temporary) / "MyResearch"
            result = subprocess.run(
                [sys.executable, str(CREATOR), str(dest),
                 "--image", "scrimp-mpi:v1.3.0-rc1"],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            config = json.loads(
                (dest / ".devcontainer/devcontainer.json").read_text()
            )
            self.assertEqual(config["image"], "scrimp-mpi:v1.3.0-rc1")
            self.assertTrue((dest / "simulations/first_simulation.py").is_file())
            self.assertTrue((dest / "outputs/.gitignore").is_file())
            self.assertFalse((dest / "sandbox").exists())
            self.assertFalse((dest / ".git").exists())
            marker = dest / "important_research_data.txt"
            marker.write_text("DO NOT OVERWRITE", encoding="utf-8")
            again = subprocess.run(
                [sys.executable, str(CREATOR), str(dest)],
                capture_output=True, text=True,
            )
            self.assertNotEqual(again.returncode, 0)
            self.assertEqual(marker.read_text(), "DO NOT OVERWRITE")

    def test_release_zip_excludes_user_files_and_matches_template(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "template.zip"
            result = subprocess.run(
                [sys.executable,
                 str(ROOT / "scripts/package-scrimp-project-template.py"),
                 "--output", str(output),
                 "--image", "ghcr.io/g-haine/scrimp:1.3.0"],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with ZipFile(output) as artifact:
                entries = artifact.namelist()
                self.assertIn(
                    "scrimp-project/.devcontainer/devcontainer.json", entries)
                cfg = json.loads(artifact.read(
                    "scrimp-project/.devcontainer/devcontainer.json"))
                self.assertEqual(
                    cfg["image"], "ghcr.io/g-haine/scrimp:1.3.0"
                )
                self.assertIn(
                    "scrimp-project/simulations/first_simulation.py", entries)
                self.assertNotIn("scrimp-project/outputs/oscillator.csv", entries)
                self.assertFalse(any(
                    "sandbox/" in x or ".git/" in x or ".build/" in x
                    for x in entries
                ))
                self.assertTrue(all(
                    x.startswith("scrimp-project/") for x in entries
                ))

    def test_generator_dry_run_does_not_create_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            dest = Path(temporary) / "preview"
            result = subprocess.run(
                [sys.executable, str(CREATOR), str(dest), "--dry-run"],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse(dest.exists())


if __name__ == "__main__":
    unittest.main()
