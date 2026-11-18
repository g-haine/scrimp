#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Standalone opt-in pH PCFIELDSPLIT experiment, 1/2 MPI ranks.
# Produces JSON, full logs; reports failed convergence honestly.
set -euo pipefail
mesh="${1:-0.08}"
x_pc="${2:-bjacobi}"
e_pc="${3:-lu}"
schur_pc="${4:-lu}"
if [[ $# -gt 4 || ! "$mesh" =~ ^0\.[0-9]+$ || ! "$x_pc" =~ ^(bjacobi|lu|schur-selfp)$ || ! "$e_pc" =~ ^(lu|bjacobi)$ || ! "$schur_pc" =~ ^(lu|bjacobi)$ ]]; then
    echo "Usage: bash sandbox/run_ph_fieldsplit_baseline.sh [mesh_h=0.08] [x_pc=bjacobi|lu|schur-selfp] [e_pc=lu|bjacobi] [schur_pc=lu|bjacobi]" >&2
    exit 2
fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
host_workspace="${SCRIMP_DOCKER_WORKSPACE:-$HOME/Simulations/workspace}"
mkdir -p "$host_workspace" sandbox/outputs
host_workspace="$(realpath "$host_workspace")"
[[ -w "$host_workspace" ]] || {
    echo "Not writable: $host_workspace" >&2
    exit 2
}
out="$(mktemp -d "$root/sandbox/outputs/ph-fieldsplit-${x_pc}-e${e_pc}-s${schur_pc}-h${mesh//./p}-$sha-XXXXXX")"
timeout_s="${SCRIMP_PH_FIELDSPLIT_TIMEOUT_S:-240}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] || (( timeout_s < 30 || timeout_s > 1800 )); then
    echo "SCRIMP_PH_FIELDSPLIT_TIMEOUT_S must be 30..1800" >&2
    exit 2
fi
echo "SCRIMP experimental FGMRES/PCFIELDSPLIT vs full MUMPS"
echo "Source commit: $sha; mesh h=$mesh; shifts=1000,2000; x_pc=$x_pc; e_pc=$e_pc; schur_pc=$schur_pc"
echo "No production solver changes; output: $out"
if ! git diff --quiet HEAD -- scrimp sandbox/probe_ph_fieldsplit_baseline.py; then
    echo "WARNING: source differs from HEAD; provenance is not clean."
    git status --short -- scrimp sandbox/probe_ph_fieldsplit_baseline.py
fi
sudo -v
fail=0
for ranks in 1 2; do
    cpus=0
    [[ "$ranks" == 2 ]] && cpus=0,2
    echo "RUN $ranks MPI ranks / CPUs=$cpus"
    if sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        --workdir /scrimp-work \
        -e HOME=/tmp -e MPLBACKEND=Agg \
        -e PYTHONPATH=/workspace/scrimp \
        -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
        -e SCRIMP_GETFEM_MATRIX_MODE=distributed \
        -e SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
        -e SCRIMP_MATRIX_MEMORY_MODE=compact \
        -e SCRIMP_MUMPS_RESET_POLICY=never \
        -e SCRIMP_BENCH_GIT_SHA="$sha" \
        --cpuset-cpus="$cpus" \
        -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
        -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
        -v "$root/examples:/workspace/scrimp/examples:ro" \
        -v "$host_workspace:/scrimp-work" \
        -v "$out:/artifacts" \
        scrimp-mpi:local \
        timeout -k 5s "${timeout_s}s" \
        mpiexec --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_ph_fieldsplit_baseline.py \
          --mesh "$mesh" --x-pc "$x_pc" --e-pc "$e_pc" --schur-pc "$schur_pc" --output "/artifacts/fieldsplit_${ranks}r.json" \
        > "$out/fieldsplit_${ranks}r.log" 2>&1; then
        grep -aE '^(PASS|FAIL) sigma=|^Saved:|^SCHUR-SELF-P FIELDS:|^SCHUR-Y_L COUPLINGS:' "$out/fieldsplit_${ranks}r.log" || true
    else
        echo "NON-PASS: $ranks MPI ranks; inspect $out/fieldsplit_${ranks}r.log" >&2
        tail -n 45 "$out/fieldsplit_${ranks}r.log" >&2
        fail=1
    fi
done
echo "Reports: $out"
if (( fail )); then
    echo "Experiment did not pass all MPI/accuracy gates" >&2
    exit 1
fi
echo "PASS: block-baseline manufactured solves match full MUMPS"
