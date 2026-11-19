#!/usr/bin/env python3
"""Package a standalone SCRIMP user project ZIP for a future GitHub release.

No Docker, conda, Git metadata, computed meshes or simulation outputs enter
this release artifact. This does not publish anything to GitHub/GHCR.
"""
from __future__ import annotations

import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "project-template"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "dist/scrimp-project-template.zip")
    args = parser.parse_args()
    output = args.output.expanduser().resolve()
    if output.is_relative_to(TEMPLATE.resolve()):
        parser.error("ZIP output must not be stored inside project-template/")
    # The project generator is the single authoritative allowlist.
    generator = SourceFileLoader(
        "scrimp_project_generator", str(ROOT / "scripts/create-scrimp-project.py")
    ).load_module()
    files = generator.TEMPLATE_FILES
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
