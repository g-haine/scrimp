#!/usr/bin/env python3
"""Stage a physically isolated, Git-tracked SCRIMP Docker context.

The *real* repository may contain enormous untracked simulation outputs.
Never pass the repository root to Docker: BuildKit ignore negations may
re-expose whole directories. Instead stage only approved tracked files.
"""
from __future__ import annotations

import argparse
import collections
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
ALLOW = (
    "setup.py",
    "pyproject.toml",
    "README.md",
    "LICENSE",
    "docker/Dockerfile",
    "docker/environment.yml",
    "scrimp/*.py",
    "scrimp/GUI/*.py",
    "scrimp/utils/*.py",
    "examples/*.py",
    "examples/published/*/*.py",
    "sandbox/*.py",
    "sandbox/*.sh",
    "tests/*.py",
    "tests/perf/*.py",
    "scripts/*.py",
    "scripts/*.sh",
    "scripts/patches/getfem-5.5/*.patch",
)

# These specific entries must remain available or the image cannot build
# or run its minimum native smoke. The rest are controlled by ALLOW.
REQUIRED = (
    "setup.py",
    "pyproject.toml",
    "README.md",
    "LICENSE",
    "docker/Dockerfile",
    "docker/environment.yml",
    "scrimp/__init__.py",
    "scrimp/dphs.py",
    "sandbox/smoke_modern_stack.py",
    "scripts/install-petsc-python-bindings.sh",
    "scripts/build-getfem-mpi.sh",
    "scripts/patches/getfem-5.5/0000-mpi-assembly-timing.patch",
    "scripts/patches/getfem-5.5/0001-bdm-tetrahedron-unisolvent-dofs.patch",
    "scripts/patches/getfem-5.5/0002-python-interface-link-in-tree.patch",
    "scripts/patches/getfem-5.5/0003-nedelec-first-kind-simplices.patch",
    "scripts/patches/getfem-5.5/0004-nedelec-second-kind-simplices.patch",
)

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage-to", type=Path, metavar="EMPTY_DIRECTORY",
        help="copy only allowlisted Git-tracked sources into this "
             "physically isolated temporary Docker context")
    args = parser.parse_args()

    # A direct root-context build must be a FAIL-CLOSED operation.
    # Only the external staged directory can be used for a real build.
    ignore_lines = (ROOT / ".dockerignore").read_text().splitlines()
    rules = [line.strip() for line in ignore_lines
             if line.strip() and not line.lstrip().startswith("#")]
    assert rules == ["**"], (
        "Root .dockerignore must contain only deny-all '**' "
        "(no broad negation exceptions)"
    )
    dockerfile = (ROOT / "docker/Dockerfile").read_text()
    assert not re.search(r"(?m)^COPY\s+(?:--\S+\s+)*\.\s+", dockerfile), (
        "Broad Dockerfile COPY . is not permitted"
    )

    # Anchored, slash-aware pattern matching (never pathlib/glob of
    # the live repository: that can enumerate giant data directories).
    patterns = [
        re.compile("^" + re.escape(p).replace(r"\*", "[^/]*") + "$")
        for p in ALLOW
    ]
    tracked = subprocess.run(
        ["git", "ls-files", "--cached", "-z", "--full-name"],
        cwd=ROOT, check=True, stdout=subprocess.PIPE,
    ).stdout
    matched: dict[str, Path] = {}
    root_resolved = ROOT.resolve()
    for entry in tracked.split(b"\0"):
        if not entry:
            continue
        relative = os.fsdecode(entry)
        if not any(p.fullmatch(relative) for p in patterns):
            continue
        path = ROOT / relative
        assert path.resolve().is_relative_to(root_resolved), (
            f"Tracked source escapes repository: {relative!r}"
        )
        if any(parent.is_symlink() for parent in (path, *path.parents)):
            raise RuntimeError(f"Symlink path in Docker context: {relative}")
        if not path.is_file():
            raise RuntimeError(f"Missing Git-tracked source: {relative}")
        matched[relative] = path

    for relative in REQUIRED:
        assert relative in matched, (
            f"Required Git-tracked Docker source missing: {relative}"
        )

    # Check sizes of exactly the allowlisted sources, NOT user data
    # directories such as sandbox/outputs or hidden .build.
    cap = int(os.environ.get("SCRIMP_DOCKER_CONTEXT_MAX_MB", "64")) * 1024**2
    max_one_file = 16 * 1024**2
    totals = collections.Counter()
    count = collections.Counter()
    for relative, path in matched.items():
        size = path.stat().st_size
        if size > max_one_file:
            raise RuntimeError(f"Oversized allowed source ({size} bytes): {relative}")
        category = relative.split("/", 1)[0]
        totals[category] += size
        count[category] += 1
    overall = sum(totals.values())
    if overall > cap:
        raise RuntimeError(
            f"Allowlisted source {overall / 1024**2:.1f} MiB exceeds "
            f"{cap / 1024**2:.0f} MiB limit"
        )
    print("Docker staging: Git-tracked source ONLY; root context DENIED", flush=True)
    for category in sorted(totals):
        print(f"  {category:10s}: {count[category]:3d} files, "
              f"{totals[category] / 1024**2:.2f} MiB")
    print(f"  TOTAL: {len(matched)} files, {overall / 1024**2:.2f} MiB")
    print("  Excluded: ALL untracked files (including simulations), "
          "docs, notebooks, images, meshes, recordings, .build, "
          "sandbox/outputs and all other unlisted paths.")

    if args.stage_to is not None:
        dest = args.stage_to.resolve()
        if dest.is_relative_to(root_resolved):
            raise RuntimeError("Docker stage destination must be OUTSIDE repository")
        if not dest.is_dir() or any(dest.iterdir()):
            raise RuntimeError("Docker stage destination must be an empty directory")
        for relative, path in sorted(matched.items()):
            output = dest / relative
            output.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, output, follow_symlinks=False)
        # This is a physically isolated, already safe context. Do not
        # copy the root .dockerignore, which intentionally denies all.
        (dest / ".dockerignore").write_text(
            "# Audited, physically isolated Docker context.\n"
            "# No local untracked simulations can exist here.\n",
            encoding="utf-8",
        )
        staged_bytes = sum(
            (dest / relative).stat().st_size for relative in matched
        )
        assert staged_bytes == overall, (
            "Unexpected staged context size discrepancy"
        )
        print(f"STAGED: {len(matched)} files, "
              f"{staged_bytes / 1024**2:.2f} MiB in {dest}", flush=True)
        print("Docker MUST build from this directory, never the live root.")
    else:
        print("Run with --stage-to /tmp/EMPTY_DIRECTORY "
              "to create a safe Docker context.")
    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, RuntimeError, OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: Docker source context audit failed: {exc}", file=sys.stderr)
        sys.exit(2)
