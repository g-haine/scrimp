#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Apply one strict unified-diff source patch without GNU patch.

Small, deliberately limited helper for the GetFEM v5.5 disposable Docker
smoke. Supports a single text file with ordinary context/add/remove hunks.
No fuzz, no offset search, no network access, no modifications on failure.
The source file must match the recorded patch context exactly.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re

HUNK = re.compile(r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@")


def apply_one(root: Path, patchfile: Path, dry_run: bool = False) -> Path:
    root = root.resolve()
    patch_lines = patchfile.read_text(encoding="utf-8").splitlines()

    if len(patch_lines) < 4 or not patch_lines[0].startswith("diff --git "):
        raise ValueError("Not a standard single-file unified diff")
    if not patch_lines[1].startswith("--- a/") or not patch_lines[2].startswith("+++ b/"):
        raise ValueError("Unexpected patch file path headers")
    old_name, new_name = patch_lines[1][6:], patch_lines[2][6:]
    if old_name != new_name or old_name.startswith("/"):
        raise ValueError("Expected one unchanged relative target path")
    target = (root / old_name).resolve()
    if not target.is_relative_to(root) or not target.is_file():
        raise ValueError(f"Refusing missing/outside source tree target: {target}")

    contents = target.read_text(encoding="utf-8")
    if not contents.endswith("\n"):
        raise ValueError(f"Source lacks trailing newline: {target}")
    original = contents.splitlines()
    output: list[str] = []
    cursor = 0
    i = 3
    hunks = 0

    while i < len(patch_lines):
        header = HUNK.match(patch_lines[i])
        if not header:
            raise ValueError(f"Unexpected patch header at line {i+1}")
        start, count, _, new_count = header.groups()
        start = int(start) - 1
        old_expected = int(count or 1)
        new_expected = int(new_count or 1)
        if start < cursor or start > len(original):
            raise ValueError(f"Overlapping/invalid hunk at source line {start+1}")
        output.extend(original[cursor:start])
        cursor = start
        old_seen = new_seen = 0
        i += 1
        while i < len(patch_lines) and not patch_lines[i].startswith("@@ "):
            line = patch_lines[i]
            if not line or line[0] not in " +-":
                raise ValueError(f"Unsupported patch operation at line {i+1}")
            op, value = line[0], line[1:]
            if op in " -":
                if cursor >= len(original) or original[cursor] != value:
                    raise ValueError(
                        f"Patch context mismatch at {target}:{cursor+1}; "
                        "source must match the GetFEM 5.5 baseline"
                    )
                cursor += 1
                old_seen += 1
            if op in " +":
                output.append(value)
                new_seen += 1
            i += 1
        if (old_seen, new_seen) != (old_expected, new_expected):
            raise ValueError(
                f"Incorrect hunk size: got {old_seen}/{new_seen}, "
                f"expected {old_expected}/{new_expected}"
            )
        hunks += 1
    if not hunks:
        raise ValueError("No patch hunks")
    output.extend(original[cursor:])
    patched = "\n".join(output) + "\n"
    if patched == contents:
        raise ValueError("Patch made no change")
    if not dry_run:
        target.write_text(patched, encoding="utf-8")
    print(f"{'CHECK' if dry_run else 'APPLIED'} {patchfile.name}: "
          f"{hunks} exact hunks in {target}")
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_tree", type=Path)
    parser.add_argument("patch", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    try:
        apply_one(args.source_tree, args.patch, args.dry_run)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
