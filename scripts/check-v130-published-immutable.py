# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Fail-closed guard: examples/published is an IMMUTABLE scientific archive.

The baseline is anchored to the exact Git blob SHA-1 values at commit
b6fe4933415fb19eb02a8448a9fc587432196227 (23 Nov 2026).
Do not edit, reformat, reconfigure, rename, delete or add anything under
examples/published/, even if a modern SCRIMP API has since changed.

The file list, Git-index blobs AND working-tree bytes are independently
verified. A script outside the archive is required for any adapted run.
"""

import hashlib
from pathlib import Path
import subprocess
import sys


BASELINE_COMMIT = "b6fe4933415fb19eb02a8448a9fc587432196227"
ROOT = Path(__file__).resolve().parents[1]
PREFIX = "examples/published/"
EXPECTED_GIT_BLOBS = {
    "examples/published/COMPEL_2025/Maxwell_2D.py":
        "8b52fc9c70c6ac55295fda2b6ee140252cd251a1",
    "examples/published/COMPEL_2025/scrimp.yml":
        "e5b671606089ef88671e347d70afa11b79cc3ba0",
    "examples/published/LHMNC_2024/heat.py":
        "071c9115aa851895346706bd2628c340785bd584",
    "examples/published/MCSS_2024/PV_trace_experiment.py":
        "85eae75b951f4f2eaa7fd431b63b1b10b5fdafa3",
    "examples/published/MCSS_2024/launch_all_simulations.py":
        "771fa16d7b098315b03cb3f5676ad2050239fa4e",
    "examples/published/MCSS_2024/scrimp.yml":
        "601316f221b1a0095e8ba8fb4fab1699fda4efd2",
    "examples/published/MCSS_2024/shallow_water.py":
        "ce9e0452a6019e0563abe9196f64597bb0985f21",
    "examples/published/MCSS_2024/u_tangent.py":
        "e3599d1505f946a6d107a733ba178892ce69b428",
}


def git_blob_sha(data: bytes) -> str:
    """Git SHA-1 of the exact bytes of one file (not a text-normalized hash)."""
    return hashlib.sha1(
        b"blob " + str(len(data)).encode("ascii") + b"\0" + data
    ).hexdigest()


def validate(root: Path, index: dict, untracked=(), expected=None):
    """Validate tracked index AND original worktree against a pinned baseline.

    Args:
        index: {repository_relative_path: (git_mode, git_object_sha)}
        untracked: untracked (not ignored) paths inside protected subtree.
        expected: test-only pinned digest mapping (default: real snapshots).
    Returns:
        List of deterministic errors; no file is ever mutated.
    """
    if expected is None:
        expected = EXPECTED_GIT_BLOBS
    issues = []
    observed = set(index)
    baseline = set(expected)
    for name in sorted(baseline - observed):
        issues.append(f"DELETED archived source: {name}")
    for name in sorted(observed - baseline):
        issues.append(f"ADDED/RENAMED archived source: {name}")
    for name in sorted(untracked):
        issues.append(f"UNTRACKED new source inside immutable archive: {name}")
    for name in sorted(baseline & observed):
        mode, digest = index[name]
        if mode != "100644":
            issues.append(f"ARCHIVED FILE MODE changed: {name} ({mode})")
        if digest != expected[name]:
            issues.append(f"ARCHIVED GIT INDEX changed: {name}")
        file = root / name
        if file.is_symlink() or not file.is_file():
            issues.append(f"ARCHIVED WORKTREE FILE missing/symlinked: {name}")
            continue
        if git_blob_sha(file.read_bytes()) != expected[name]:
            issues.append(f"ARCHIVED WORKTREE BYTES changed: {name}")
    return issues


def read_repository_state(root: Path):
    """Read without staging, writing, running formatters or fetching any repo."""
    raw = subprocess.check_output(
        ["git", "-C", str(root), "ls-files", "--stage", "-z", "--", PREFIX]
    )
    index = {}
    for item in raw.split(b"\0"):
        if not item:
            continue
        metadata, path = item.split(b"\t", 1)
        mode, sha, stage = metadata.decode("ascii").split()
        if stage != "0":
            raise RuntimeError(
                "Unmerged conflict inside immutable published archive"
            )
        name = path.decode("utf-8", errors="strict")
        if name in index:
            raise RuntimeError(f"Duplicate index entry in archive: {name}")
        index[name] = (mode, sha)

    raw_untracked = subprocess.check_output(
        [
            "git", "-C", str(root), "ls-files", "--others",
            "--exclude-standard", "-z", "--", PREFIX
        ]
    )
    untracked = [
        p.decode("utf-8", errors="strict")
        for p in raw_untracked.split(b"\0") if p
    ]
    return index, untracked


def main():
    index, untracked = read_repository_state(ROOT)
    errors = validate(ROOT, index, untracked)
    if errors:
        print(
            "FAIL: immutable scientific snapshots differ from anchored "
            f"commit {BASELINE_COMMIT}:",
            file=sys.stderr,
        )
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(
        f"PASS: {len(EXPECTED_GIT_BLOBS)} exact, immutable "
        f"examples/published files; pinned to {BASELINE_COMMIT}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
