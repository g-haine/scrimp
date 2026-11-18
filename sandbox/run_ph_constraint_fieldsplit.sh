#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Read-only MPI pH variable-block graph audit, 1/2 ranks.
# Produces JSON and full logs; no production solver changes.
set -euo pipefail
mesh="${1:-0.08}"
if [[ $# -gt 1 || ! "$mesh" =~ ^0\.[0-9]+$ ]]; then
    echo "Usage: bash sandbox/run_ph_constraint_fieldsplit.sh [mesh_h=0.08]" >&2
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
out="$(mktemp -d "$root/sandbox/outputs/ph-constraint-fieldsplit-h${mesh//./p}-$sha-XXXXXX")"
timeout_s="${SCRIMP_PH_FIELDSPLIT_TIMEOUT_S:-240}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] || (( timeout_s < 30 || timeout_s > 1800 )); then
    echo "SCRIMP_PH_GRAPH_TIMEOUT_S must be 30..1800" >&2
    exit 2
fi
echo "SCRIMP graph-selected constraint Schur FGMRES vs full MUMPS"
echo "Source commit: $sha; mesh h=$mesh; shifts=1000,2000"
echo "No production solver changes; output: $out"
if ! git diff --quiet HEAD -- scrimp sandbox/probe_ph_constraint_fieldsplit.py; then
    echo "WARNING: source differs from HEAD; provenance is not clean."
    git status --short -- scrimp sandbox/probe_ph_constraint_fieldsplit.py
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
        python -u /workspace/scrimp/sandbox/probe_ph_constraint_fieldsplit.py \
          --mesh "$mesh" --output "/artifacts/constraint_fieldsplit_${ranks}r.json" \
        > "$out/constraint_fieldsplit_${ranks}r.log" 2>&1; then
        grep -aE '^(PASS|REJECT) sigma=|^Saved:' "$out/constraint_fieldsplit_${ranks}r.log" || true
    else
        echo "NON-PASS: $ranks MPI ranks; inspect $out/constraint_fieldsplit_${ranks}r.log" >&2
        tail -n 45 "$out/constraint_fieldsplit_${ranks}r.log" >&2
        fail=1
    fi
done
echo "Reports: $out"
if (( fail )); then
    echo "Constraint fieldsplit failed on at least one MPI rank count" >&2
    exit 1
fi
echo "PASS: graph-selected Schur manufactured linear solves match MUMPS"
