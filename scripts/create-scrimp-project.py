#!/usr/bin/env python3
"""Create an isolated, user-owned SCRIMP/VS Code project.

This script copies ONLY the small, reviewed project template. It does
not copy SCRIMP's sandbox, numerical datasets, build files, Git metadata
or any other contents of this repository.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

SOURCE = Path(__file__).resolve().parents[1] / "project-template"
TEMPLATE_FILES = (
    ".devcontainer/devcontainer.json",
    ".vscode/launch.json",
    ".vscode/settings.json",
    ".vscode/tasks.json",
    ".gitignore",
    "README.md",
    "geometries/README.md",
    "outputs/.gitignore",
    "simulations/first_simulation.py",
    "simulations/mpi_check.py",
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path, help="new project directory")
    parser.add_argument(
        "--image",
        help="optional Docker image tag (e.g. scrimp-mpi:v1.3.0-rc1). "
             "By default the template uses SCRIMP_DOCKER_IMAGE from VS "
             "Code's host environment.")
    parser.add_argument(
        "--dry-run", action="store_true",
        help="report files that would be created without writing anything")
    args = parser.parse_args()

    destination = args.destination.expanduser().resolve()
    if destination.exists():
        parser.error(f"project destination already exists: {destination}")
    if SOURCE.resolve() == destination or destination.is_relative_to(SOURCE.resolve()):
        parser.error("cannot create a project inside the template directory")
    if args.image is not None and (not args.image.strip()
                                   or any(c.isspace() for c in args.image)):
        parser.error("--image must be one nonempty image reference")

    for relative in TEMPLATE_FILES:
        path = SOURCE / relative
        if not path.is_file() or path.is_symlink():
            parser.error(f"missing or unsafe template file: {path}")
    print(f"SCRIMP project: {destination}")
    for relative in TEMPLATE_FILES:
        print(f"  {relative}")
    if args.dry_run:
        return 0

    destination.mkdir(parents=True, exist_ok=False)
    try:
        for relative in TEMPLATE_FILES:
            target = destination / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(SOURCE / relative, target)
        if args.image:
            config_path = destination / ".devcontainer/devcontainer.json"
            config = json.loads(config_path.read_text(encoding="utf-8"))
            config["image"] = args.image
            config_path.write_text(
                json.dumps(config, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
    except Exception:
        # We created the destination exclusively; removing our own
        # incomplete project is safe. Never remove a preexisting path.
        shutil.rmtree(destination)
        raise
    print("\nCreated. Open the directory in VS Code.")
    print("Docker path: Dev Containers → Reopen in Container.")
    print("Native Conda path: Python: Select Interpreter; then run a script.")
    print("User outputs will remain in its local outputs/ directory.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
