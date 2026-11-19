#!/usr/bin/env python3
"""Preflight *source allowlist* for the SCRIMP v1.3 BuildKit context.

This is a fail-closed check for accidental broad COPY, symlinks,
unexpected source categories and a very large allowed source file.
The actual Docker context is filtered by the root .dockerignore;
the BuildKit progress log reports the transferred context size.
"""
from __future__ import annotations

import collections
import os
from pathlib import Path
import re
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
    ignore_lines = (ROOT / ".dockerignore").read_text().splitlines()
    ignore = [line.strip() for line in ignore_lines
              if line.strip() and not line.lstrip().startswith("#")]
    assert ignore and ignore[0] == "**", (
        "Docker context must start with a deny-all '**' rule"
    )
    # BuildKit pattern syntax differs from Python glob, so check that
    # *every* actual allow exception is exactly one audited pattern or
    # a parent of one. New exceptions must get an explicit review.
    allow_set = set(ALLOW)
    ancestor_rules = {"!" + str(Path(p)) for pattern in ALLOW
                      for p in Path(pattern).parents if str(p) != "."}
    expected = {"!" + p for p in ALLOW} | ancestor_rules
    actual = {x.rstrip("/") for x in ignore if x.startswith("!")}
    assert actual == expected, (
        f"Unexpected .dockerignore exceptions: "
        f"extra={sorted(actual - expected)}; "
        f"missing={sorted(expected - actual)}"
    )

    dockerfile = (ROOT / "docker/Dockerfile").read_text()
    assert not re.search(r"(?m)^COPY\s+(?:--\S+\s+)*\.\s+", dockerfile), (
        "Broad Dockerfile COPY . is not permitted"
    )
    for relative in REQUIRED:
        assert (ROOT / relative).is_file(), f"Required source missing: {relative}"

    matched: dict[str, Path] = {}
    for pattern in ALLOW:
        for path in ROOT.glob(pattern):
            if not path.is_file():
                continue
            relative = path.relative_to(ROOT).as_posix()
            if path.is_symlink():
                raise RuntimeError(f"Symlink in Docker source allowlist: {relative}")
            matched[relative] = path

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
    print("Docker input allowlist: DENY BY DEFAULT, no COPY .", flush=True)
    for category in sorted(totals):
        print(f"  {category:10s}: {count[category]:3d} files, "
              f"{totals[category] / 1024**2:.2f} MiB")
    print(f"  TOTAL: {len(matched)} files, {overall / 1024**2:.2f} MiB")
    print("  Excluded: docs, notebooks, images, meshes, recordings, "
          "sandbox/outputs, .build and all other unlisted paths.")
    print("  Docker BuildKit will report actual transferred context size.")
    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, RuntimeError, OSError, ValueError) as exc:
        print(f"ERROR: Docker source context audit failed: {exc}", file=sys.stderr)
        sys.exit(2)
