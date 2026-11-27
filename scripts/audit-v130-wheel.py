#!/usr/bin/env python3
"""Build and audit a v1.3.0 wheel in an isolated directory.

This is a *package metadata and payload* acceptance test, not a native
GetFEM/PETSc runtime test. It uses no import of SCRIMP and no user data.
"""
from __future__ import annotations

import argparse
import json
from email.parser import Parser
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import venv
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="new wheel audit directory")
    args = parser.parse_args()
    (ROOT / "sandbox/outputs").mkdir(parents=True, exist_ok=True)
    output = (args.output or
              Path(tempfile.mkdtemp(prefix="v130-wheel-audit-",
                                     dir=ROOT / "sandbox/outputs"))).resolve()
    if args.output:
        if output.exists():
            parser.error(f"refusing existing audit directory: {output}")
        output.mkdir(parents=True)
    # Copy only Git-tracked Python/package metadata into an isolated
    # staging tree. Do not let setuptools scan untracked user files,
    # nor leave build/egg-info artifacts in the developer repository.
    tracked_changes = subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        cwd=ROOT, text=True)
    if tracked_changes.strip():
        parser.error("tracked source changes are present; commit first")
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    source = output / "source"
    source.mkdir()
    tracked = subprocess.check_output(
        ["git", "ls-files", "--cached", "-z", "--full-name"],
        cwd=ROOT).split(b"\0")
    needed = {"setup.py", "pyproject.toml", "README.md", "LICENSE"}
    copied = set()
    for raw in tracked:
        if not raw:
            continue
        relative = raw.decode("utf-8")
        if relative not in needed and not (
            relative.startswith("scrimp/") and relative.endswith(".py")
        ):
            continue
        incoming = ROOT / relative
        if incoming.is_symlink() or not incoming.is_file():
            parser.error(f"missing/unsafe tracked package source: {relative}")
        target = source / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(incoming, target)
        copied.add(relative)
    if not needed.issubset(copied) or "scrimp/__init__.py" not in copied:
        parser.error("missing required tracked wheel source files")
    wheel_dir = output / "wheel"
    wheel_dir.mkdir()
    print(f"Wheel audit evidence: {output}; source: {revision}", flush=True)
    subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-deps",
                    "--no-build-isolation", "--wheel-dir", str(wheel_dir),
                    str(source)], check=True)
    wheels = sorted(wheel_dir.glob("scrimp-1.3.0-*.whl"))
    assert len(wheels) == 1, f"expected exactly one SCRIMP 1.3.0 wheel: {wheels}"
    wheel = wheels[0]
    with ZipFile(wheel) as archive:
        assert archive.testzip() is None, "wheel CRC failure"
        names = archive.namelist()
        metadata_files = [x for x in names
                          if x.endswith(".dist-info/METADATA")]
        assert len(metadata_files) == 1, metadata_files
        meta = Parser().parsestr(
            archive.read(metadata_files[0]).decode("utf-8"))
        assert meta["Name"].lower() == "scrimp", meta["Name"]
        assert meta["Version"] == "1.3.0", meta["Version"]
        assert "scrimp/__init__.py" in names
        assert "scrimp/dphs.py" in names
        assert not any(x.startswith(("sandbox/", "outputs/", "tests/",
                                     "docs/", ".git/")) for x in names), (
            "unexpected release wheel payload"
        )
        for name in names:
            if name.startswith("scrimp/") and name.endswith(".py"):
                compile(archive.read(name), name, "exec")
    environment = output / "clean-venv"
    venv.EnvBuilder(with_pip=True).create(environment)
    interpreter = environment / "bin/python"
    if sys.platform == "win32":
        interpreter = environment / "Scripts/python.exe"
    subprocess.run([str(interpreter), "-m", "pip", "install",
                    "--no-index", "--no-deps", str(wheel)], check=True)
    # The clean environment deliberately does not contain native GetFEM,
    # PETSc or Gmsh. Query packaging metadata without importing SCRIMP.
    probe = ('import importlib.metadata as m; '
             'assert m.version("scrimp") == "1.3.0"; '
             'print("PASS: clean wheel metadata and installation")')
    subprocess.run([str(interpreter), "-c", probe], check=True)
    result = {
        "schema": "scrimp-v130-wheel-audit-v1",
        "status": "PASS_WHEEL_BUILD_PAYLOAD_CLEAN_METADATA_INSTALL",
        "wheel": str(wheel),
        "distribution_name": "scrimp",
        "distribution_version": "1.3.0",
        "source_revision": revision,
        "file_count": len(names),
        "scope": "wheel build, payload and isolated metadata installation",
        "not_tested": ["SCRIMP import without native GetFEM/PETSc",
                       "Docker runtime physics or MPI", "fresh machine"],
    }
    (output / "summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
