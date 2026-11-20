#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Package a standalone SCRIMP user project ZIP for a future GitHub release.

No Docker, conda, Git metadata, computed meshes or simulation outputs enter
this release artifact. This does not publish anything to GitHub/GHCR.
"""
from __future__ import annotations

import argparse
import json
from runpy import run_path
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "project-template"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "dist/scrimp-project-template.zip")
    parser.add_argument(
        "--image", help="embed a specific built/published Docker image "
                        "in devcontainer.json; omit for the pre-release default")
    args = parser.parse_args()
    if args.image is not None and (
        not args.image.strip() or any(c.isspace() for c in args.image)
    ):
        parser.error("--image must be one nonempty Docker image reference")
    output = args.output.expanduser().resolve()
    if output.is_relative_to(TEMPLATE.resolve()):
        parser.error("ZIP output must not be stored inside project-template/")
    # The project generator is the single authoritative allowlist.
    files = run_path(str(ROOT / "scripts/create-scrimp-project.py"))[
        "TEMPLATE_FILES"
    ]
    for path in files:
        source = TEMPLATE / path
        if not source.is_file() or source.is_symlink():
            parser.error(f"missing or unsafe release template file: {path}")
        if source.stat().st_size > 512 * 1024:
            parser.error(f"oversized template file: {path}")

    output.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(output, mode="w", compression=ZIP_DEFLATED,
                 compresslevel=9) as zip_file:
        for relative in sorted(files):
            filedata = (TEMPLATE / relative).read_bytes()
            if relative == ".devcontainer/devcontainer.json" and args.image:
                config = json.loads(filedata.decode("utf-8"))
                config["image"] = args.image
                filedata = (
                    json.dumps(config, indent=2, ensure_ascii=False) + "\n"
                ).encode("utf-8")
            member = ZipInfo(f"scrimp-project/{relative}",
                             date_time=(2026, 1, 1, 0, 0, 0))
            member.compress_type = ZIP_DEFLATED
            member.external_attr = 0o100644 << 16
            zip_file.writestr(member, filedata)
    print(f"SCRIMP user-project archive: {output}")
    print(f"Files: {len(files)}, bytes: {output.stat().st_size}")
    print("Not published: attach only to a tested GitHub release.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
