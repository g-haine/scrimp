#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# Final v1.3.0 release-branch default-path gate (no PILOT flags).
# This is *not* a tag, merge to master, Docker rebuild or published release.
# Reuses the real strict DAE/CN/Hamiltonian/VTK gates on the pinned image.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

export SCRIMP_V130_PRODUCTION_DEFAULT_SMOKE=1
export SCRIMP_V130_MPI_PV_ROOT_ONLY=1
unset SCRIMP_V130_MPI_PV_DIAGNOSTIC
unset SCRIMP_V130_AUTO_INTERFACE_PILOT
unset SCRIMP_V130_PV_ROOT_ONLY_PILOT
unset SCRIMP_AUTO_INTERFACE
unset SCRIMP_PV_EXPORT_MODE

exec bash scripts/run-v130-auto-interface-mpi-cn.sh
