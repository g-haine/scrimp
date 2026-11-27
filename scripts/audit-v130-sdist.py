#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Audit the v1.3.0 source distribution without importing native SCRIMP.

Only Git-tracked Python package sources and package metadata are staged.
No user simulations, untracked data, GetFEM, PETSc, or build directories
are allowed into the source tarball. Rebuild a wheel *from the sdist* and
compare each Python package source byte-for-byte with the staged inputs.
This verifies packaging integrity, not native scientific runtime behavior.
"""
from __future__ import annotations

import argparse
from email.parser import Parser
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile
import tempfile
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
VERSION = "1.3.0"
PREFIX = f"scrimp-{VERSION}"


def run(*args: str, cwd: Path | None = None) -> str:
    return subprocess.check_output(args, cwd=cwd or ROOT, text=True).strip()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        help="new, nonexistent evidence directory")
    args = parser.parse_args()

    if run("git", "status", "--porcelain", "--untracked-files=no"):
        parser.error("tracked source changes: commit before source packaging")
    revision = run("git", "rev-parse", "HEAD")
    base = ROOT / "sandbox" / "outputs"
    base.mkdir(parents=True, exist_ok=True)
    out = (args.output or Path(tempfile.mkdtemp(
        prefix="v130-sdist-audit-", dir=base))).resolve()
    if args.output:
        if out.exists():
            parser.error("refusing existing audit directory")
        out.mkdir(parents=True)

    stage = out / "staged-source"
    stage.mkdir()
    mandatory = {"setup.py", "pyproject.toml", "README.md", "LICENSE"}
    staged: dict[str, bytes] = {}
    tracked = subprocess.check_output(
        ["git", "ls-files", "--cached", "-z", "--full-name"], cwd=ROOT
    ).split(b"\0")
    for raw in tracked:
        if not raw:
            continue
        rel = raw.decode("utf-8")
        if rel not in mandatory and not (
            rel.startswith("scrimp/") and rel.endswith(".py")
        ):
            continue
        src = ROOT / rel
        if src.is_symlink() or not src.is_file() or src.stat().st_size > 512 * 1024:
            raise RuntimeError(f"invalid or oversized Git-tracked package file: {rel}")
        data = src.read_bytes()
        target = stage / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        staged[rel] = data
    if not mandatory.issubset(staged) or not {
        "scrimp/__init__.py", "scrimp/dphs.py"
    }.issubset(staged):
        raise RuntimeError("required SCRIMP package files missing in stage")

    tar_dir = out / "sdist"
    tar_dir.mkdir()
    subprocess.run(
        [sys.executable, "-m", "build", "--sdist", "--no-isolation",
         "--outdir", str(tar_dir), str(stage)],
        check=True,
    )
    archives = sorted(tar_dir.glob(f"{PREFIX}.tar.gz"))
    if len(archives) != 1:
        raise AssertionError(f"expected one SCRIMP sdist, got {archives}")

    unpack = out / "unpacked"
    unpack.mkdir()
    files: dict[str, bytes] = {}
    with tarfile.open(archives[0], "r:gz") as archive:
        members = archive.getmembers()
        if not members or len(members) > 1000:
            raise AssertionError("unexpected sdist member count")
        for member in members:
            path = PurePosixPath(member.name)
            if (path.is_absolute() or ".." in path.parts
                or not path.parts or path.parts[0] != PREFIX
                or not (member.isfile() or member.isdir())
                or member.size > 5 * 1024 * 1024):
                raise AssertionError(f"unsafe or unexpected sdist entry: {member.name}")
            target = unpack.joinpath(*path.parts)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            stream = archive.extractfile(member)
            if stream is None:
                raise AssertionError(f"unreadable entry: {member.name}")
            contents = stream.read()
            rel = "/".join(path.parts[1:])
            if rel in files:
                raise AssertionError(f"duplicate file: {rel}")
            if rel.startswith(("sandbox/", "outputs/", "docs/", "tests/", ".git/")):
                raise AssertionError(f"unwanted source payload: {rel}")
            if rel.endswith((".png", ".zip", ".npz", ".msh", ".vtu", ".csv")):
                raise AssertionError(f"binary/output in source distribution: {rel}")
            files[rel] = contents
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(contents)

    for relative, expected in staged.items():
        if files.get(relative) != expected:
            raise AssertionError(f"sdist missing/altered source: {relative}")
    pkg_info = files.get("PKG-INFO")
    if pkg_info is None:
        raise AssertionError("PKG-INFO missing")
    metadata = Parser().parsestr(pkg_info.decode("utf-8"))
    if (metadata["Name"].lower(), metadata["Version"]) != ("scrimp", VERSION):
        raise AssertionError(f"incorrect sdist metadata: {metadata['Name']} {metadata['Version']}")

    wheels = out / "rebuilt-wheel"
    wheels.mkdir()
    subprocess.run([
        sys.executable, "-m", "pip", "wheel", "--no-deps",
        "--no-build-isolation", "--wheel-dir", str(wheels),
        str(unpack / PREFIX)], check=True)
    built = sorted(wheels.glob(f"scrimp-{VERSION}-*.whl"))
    if len(built) != 1:
        raise AssertionError(f"wheel from sdist not built: {built}")
    with ZipFile(built[0]) as zip_file:
        if zip_file.testzip() is not None:
            raise AssertionError("wheel from source has failed CRC")
        zipped = set(zip_file.namelist())
        for rel, data in staged.items():
            if rel.startswith("scrimp/") and rel.endswith(".py"):
                if rel not in zipped or zip_file.read(rel) != data:
                    raise AssertionError(f"wheel from sdist missing/altered: {rel}")
        metas = [p for p in zipped if p.endswith(".dist-info/METADATA")]
        if len(metas) != 1:
            raise AssertionError("missing wheel metadata")
        wheel_meta = Parser().parsestr(zip_file.read(metas[0]).decode("utf-8"))
        if (wheel_meta["Name"].lower(), wheel_meta["Version"]) != (
            "scrimp", VERSION
        ):
            raise AssertionError("reconstructed wheel metadata mismatch")

    result = {
        "schema": "scrimp-v130-sdist-audit-v1",
        "status": "PASS_SDIST_TRACKED_PAYLOAD_REBUILT_WHEEL",
        "source_revision": revision,
        "sdist": str(archives[0]),
        "wheel_from_sdist": str(built[0]),
        "source_file_count": len(staged),
        "sdist_file_count": len(files),
        "scope": "Git-tracked package source bytes, source archive, metadata, wheel rebuild",
        "not_tested": ["native GetFEM/PETSc import", "MPI physics",
                       "installation on non-Linux machines"],
    }
    (out / "summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
