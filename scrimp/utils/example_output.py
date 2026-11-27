# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# SPDX-License-Identifier: GPL-3.0-or-later
"""Deterministic, per-system user output layout (no GetFEM/PETSc dependency).

Examples should save to:
  <working_directory>/outputs/<system_name>/
      png/  svg/  gif/  mesh/  log/  npz/  csv/  json/

Directories are created only when actual files are written. A custom
output_dir is interpreted as the whole <system_name> root, so existing
users passing explicit output directories retain complete control.
No rank suffix is introduced for shared figures/VTK; MPI-only raw
evidence remains explicitly rank-labelled by the testing harnesses.
"""
from pathlib import Path
import re


KINDS = frozenset(("png", "svg", "gif", "mesh", "log", "npz", "csv", "json"))


def output_root(system_name, *, working_dir=None, output_dir=None):
    """Resolve the directory containing all outputs of one system.

    The name is a single safe path component. Output directory overrides
    are accepted without modifying or deleting their existing contents.
    """
    if output_dir is not None:
        return Path(output_dir).expanduser().resolve()
    name = str(system_name)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name) or name in (".", ".."):
        raise ValueError("system_name must be a nonempty safe folder name")
    work = Path.cwd() if working_dir is None else Path(working_dir)
    return (work / "outputs" / name).expanduser().resolve()


def output_folder(root, kind, *, create=True):
    """Return a typed output directory, optionally creating it."""
    if kind not in KINDS:
        raise ValueError("Unknown SCRIMP output kind: " + str(kind))
    folder = Path(root) / kind
    if create:
        folder.mkdir(parents=True, exist_ok=True)
    return folder


def output_file(root, name, kind, *, create=True):
    """Return a format-specific file path, rejecting path traversal."""
    name = str(name)
    if (not name or name in (".", "..") or Path(name).name != name
            or "/" in name or "\\" in name):
        raise ValueError("Output filename must be a simple basename")
    allowed = {
        "mesh": (".msh", ".pvd", ".vtu", ".pvtu"),
        "log": (".log", ".txt"),
        "json": (".json",),
        "png": (".png",),
        "svg": (".svg",),
        "gif": (".gif",),
        "npz": (".npz",),
        "csv": (".csv",),
    }
    if kind not in allowed:
        raise ValueError("Unknown SCRIMP output kind: " + str(kind))
    if not any(name.lower().endswith(suffix) for suffix in allowed[kind]):
        raise ValueError(
            f"Output filename {name!r} must use one of {allowed[kind]}"
        )
    return output_folder(root, kind, create=create) / name
