#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# Build petsc4py / slepc4py against the EXACT conda-forge complex
# PETSc/SLEPc libraries already installed in CONDA_PREFIX. This
# deliberately avoids conda-forge petsc4py builds requiring NumPy<2.
set -euo pipefail
: "${CONDA_PREFIX:?Activate the target conda/micromamba environment}"
export PETSC_DIR="$CONDA_PREFIX"
export SLEPC_DIR="$CONDA_PREFIX"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

versions="$(python - <<'PY'
from pathlib import Path
import json
import os
root = Path(os.environ["CONDA_PREFIX"]) / "conda-meta"
for name in ("petsc", "slepc"):
    versions = [
        json.loads(path.read_text())["version"]
        for path in root.glob(f"{name}-*.json")
        if json.loads(path.read_text())["name"] == name
    ]
    if len(versions) != 1:
        raise RuntimeError(f"Expected exactly one native {name} package: {versions}")
    print(versions[0])
PY
)"
PETSC_VERSION="$(printf '%s\n' "$versions" | sed -n '1p')"
SLEPC_VERSION="$(printf '%s\n' "$versions" | sed -n '2p')"
echo "Native PETSc=$PETSC_VERSION, SLEPc=$SLEPC_VERSION; building Python bindings for NumPy=$(python -c 'import numpy; print(numpy.__version__)')"
python -m pip install --no-deps --no-build-isolation --no-cache-dir \
    "petsc4py==$PETSC_VERSION"
python -m pip install --no-deps --no-build-isolation --no-cache-dir \
    "slepc4py==$SLEPC_VERSION"
python - <<'PY'
import numpy as np
from petsc4py import PETSc
from slepc4py import SLEPc
assert int(np.__version__.split('.')[0]) >= 2
assert np.issubdtype(PETSc.ScalarType, np.complexfloating)
print("PETSc ABI OK:", PETSc.Sys.getVersion(), "complex", PETSc.ScalarType)
print("SLEPc ABI import OK:", SLEPc)
PY
