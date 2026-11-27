# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""Docker-first v1.3.0 release source contracts (no Docker/Conda needed)."""
from __future__ import annotations

import ast
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def setup_version() -> str:
    tree = ast.parse((ROOT / "setup.py").read_text(encoding="utf-8"))
    calls = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "setup"
    ]
    if len(calls) != 1:
        raise AssertionError("expected exactly one setuptools.setup()")
    values = [
        ast.literal_eval(kw.value)
        for kw in calls[0].keywords if kw.arg == "version"
    ]
    if len(values) != 1:
        raise AssertionError("expected exactly one package version")
    return values[0]


class DockerFirstReleaseSourceTests(unittest.TestCase):
    def test_versions_are_consistent_before_build(self):
        conf = (ROOT / "docs/source/conf.py").read_text(encoding="utf-8")
        self.assertEqual(setup_version(), "1.3.0")
        self.assertRegex(conf, r"(?m)^version = '1\.3\.0'$")
        self.assertRegex(conf, r"(?m)^release = '1\.3\.0'$")

    def test_docker_is_default_and_conda_is_alternative(self):
        for relative in ("README.md", "docs/source/install.rst"):
            with self.subTest(path=relative):
                text = (ROOT / relative).read_text(encoding="utf-8")
                docker_pos = text.index("Docker")
                conda_pos = text.index("Conda") if "Conda" in text else -1
                self.assertGreaterEqual(conda_pos, 0)
                self.assertLess(docker_pos, conda_pos)
                self.assertIn("ghcr.io/g-haine/scrimp:1.3.0", text)
                self.assertIn("not yet", text.lower())
                self.assertIn("independent", text.lower())
        self.assertIn("Docker", (ROOT / "docs/source/started.rst").read_text())

    def test_no_premature_release_instructions(self):
        text = (ROOT / "docs/V130_DOCKER_FIRST_RELEASE_GATE.md").read_text()
        self.assertIn("not yet", text.lower())
        self.assertIn("fresh-machine", text)
        self.assertIn("wheel", text.lower())
        self.assertIn("scrimp-project-template.zip", text)

    def test_distribution_gate_fail_closed_and_source_isolation(self):
        text = (ROOT / "scripts/run-v130-docker-distribution-gate.sh"
                ).read_text(encoding="utf-8")
        self.assertIn("SCRIMP_V130_EXPECT_REVISION", text)
        self.assertIn("SCRIMP_V130_DISTRIBUTION_IMAGE", text)
        self.assertIn("org.opencontainers.image.revision", text)
        self.assertIn('git -C "$ROOT" rev-parse HEAD', text)
        self.assertIn('status --porcelain --untracked-files=no', text)
        self.assertIn('version("scrimp")', text)
        self.assertIn('assert actual == "1.3.0"', text)
        self.assertIn("scrimp-project-template.zip", text)
        self.assertIn("outputs/oscillator.csv", text)
        self.assertIn("simulations/mpi_check.py", text)
        self.assertIn("testzip()", text)
        self.assertNotIn("--push", text)
        self.assertNotIn("docker pull", text)
        self.assertNotIn("docker system prune", text)
        self.assertIn('src=$PROJECT,dst=/work', text)
        self.assertNotIn('src=$ROOT,dst=/work', text)

    def test_docker_release_zip_is_separate_from_source_preview(self):
        import json
        config = json.loads((ROOT / "project-template/.devcontainer/devcontainer.json")
                            .read_text(encoding="utf-8"))
        self.assertEqual(config["image"], "scrimp-mpi:v1.3.0-rc1")
        creator = (ROOT / "scripts/create-scrimp-project.py").read_text()
        packager = (ROOT / "scripts/package-scrimp-project-template.py").read_text()
        self.assertIn("--image", creator)
        self.assertIn("--image", packager)
        self.assertIn('config["image"] = args.image', packager)


if __name__ == "__main__":
    unittest.main()
