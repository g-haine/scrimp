#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# Opt-in production-path MPI root-only VTU/PVD writer qualification.
# All MPI ranks solve and compute physical powers. Only rank 0 writes
# the actual original full-mesh GetFEM visualization files, with a
# collective completion/error report. Does not change the default.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
export SCRIMP_V130_MPI_PV_ROOT_ONLY=1
unset SCRIMP_V130_MPI_PV_DIAGNOSTIC
# First native test also probes whether GetFEM's VTK writer can be
# entered by root only. Bound a potential hidden MPI collective hang
# before retrying with a longer user-selected deadline.
export SCRIMP_V130_MPI_CN_TIMEOUT="${SCRIMP_V130_MPI_CN_TIMEOUT:-90s}"
exec bash scripts/run-v130-auto-interface-mpi-cn.sh
