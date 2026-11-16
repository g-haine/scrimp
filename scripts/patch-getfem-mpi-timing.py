#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Make GetFEM 5.5's unconditional MPI assembly print opt-in.

Upstream src/getfem_models.cc emits "Assembly time" on rank 0 for *every*
model assembly whenever GETFEM_PARA_LEVEL > 0, ignoring the supported
GetFEM trace level. Rather than capture/redirect process-wide stdout,
modify that exact *diagnostic* statement during the native build.

Enable the original timing messages only with
SCRIMP_GETFEM_ASSEMBLY_TIMING=1 (set on all ranks before launching MPI).
Fail on an unknown upstream source instead of silently patching the
wrong code. No numerical assembly statements are changed.
"""
import argparse
from pathlib import Path
import re


UPSTREAM_PRINT = re.compile(
    r'(?m)^(?P<indent>[ \t]*)if \(MPI_IS_MASTER\(\)\) '
    r'cout << "Assembly time " << MPI_Wtime\(\)-t_ref << endl;[ \t]*$'
)
ENV_NAME = "SCRIMP_GETFEM_ASSEMBLY_TIMING"
INCLUDE_ANCHOR = "#include <iomanip>"
INCLUDE_NEEDED = "#include <cstdlib>"


def patch_source(source: str) -> str:
    """Patch one upstream print, keeping the actual assembly untouched."""
    if ENV_NAME in source:
        if UPSTREAM_PRINT.search(source):
            raise ValueError("Mixed upstream and patched assembly print statements")
        # A second build invocation must not modify the native source.
        return source

    matches = list(UPSTREAM_PRINT.finditer(source))
    if len(matches) != 1:
        raise ValueError(
            "Expected exactly one unconditional GetFEM MPI assembly "
            f"timing print; found {len(matches)}. Audit upstream source."
        )
    if INCLUDE_NEEDED not in source:
        if source.count(INCLUDE_ANCHOR) != 1:
            raise ValueError("Cannot locate unique upstream C++ include anchor")
        source = source.replace(
            INCLUDE_ANCHOR, INCLUDE_ANCHOR + "\n" + INCLUDE_NEEDED, 1
        )
    # Use a second substitution because the include insertion changes
    # positions (match offsets) in the original string.
    def replacement(match):
        indent = match.group("indent")
        return (
            f"{indent}if (MPI_IS_MASTER()) {{\n"
            f'{indent}  const char *flag = std::getenv("{ENV_NAME}");\n'
            f"{indent}  if (flag && flag[0] == '1' && flag[1] == '\\0')\n"
            f'{indent}    cout << "Assembly time " << MPI_Wtime()-t_ref << endl;\n'
            f"{indent}}}"
        )

    return UPSTREAM_PRINT.sub(replacement, source, count=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="src/getfem_models.cc")
    args = parser.parse_args()
    original = args.source.read_text(encoding="utf-8")
    updated = patch_source(original)
    if updated != original:
        args.source.write_text(updated, encoding="utf-8")
        print(
            "GetFEM MPI assembly timing is opt-in: "
            f"set {ENV_NAME}=1 to restore native per-assembly messages."
        )
    else:
        print("GetFEM MPI assembly timing patch already applied.")


if __name__ == "__main__":
    main()
