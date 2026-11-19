# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Standalone sandbox compatibility entrypoint for installed SCRIMP writer.

Load the real production implementation by filename, without importing
`scrimp.__init__` (and thereby GetFEM/PETSc). This also lets source CI
exercise the I/O protocol with only NumPy installed.
"""
import importlib.util
from pathlib import Path

_SOURCE = (
    Path(__file__).resolve().parents[1]
    / "scrimp" / "utils" / "async_binary_trajectory.py"
)
_SPEC = importlib.util.spec_from_file_location(
    "_scrimp_binary_transport_impl", _SOURCE
)
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)

MAGIC = _MODULE.MAGIC
HEADER = _MODULE.HEADER
FRAME = _MODULE.FRAME
DATA = _MODULE.DATA
END = _MODULE.END
DTYPES = _MODULE.DTYPES
CODE = _MODULE.CODE
AsyncOutputBackpressure = _MODULE.AsyncOutputBackpressure
RankSender = _MODULE.RankSender
LazyTrajectory = _MODULE.LazyTrajectory
start_writer = _MODULE.start_writer
serve = _MODULE.serve
main = _MODULE.main

if __name__ == "__main__":
    main()
