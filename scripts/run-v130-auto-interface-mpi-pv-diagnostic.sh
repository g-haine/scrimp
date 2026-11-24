#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# RANK-ISOLATED original GetFEM VTU/PVD export diagnostic under real
# MPI1 / MPI2 and two physical region assignments. Does NOT qualify
# shared production ParaView output, concurrency or PNG plotting.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
export SCRIMP_V130_MPI_PV_DIAGNOSTIC=1
exec bash scripts/run-v130-auto-interface-mpi-cn.sh
