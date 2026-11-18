#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Native MUMPS coenergy-block experiment, NOT a production solver.
# Checks exact non-condensed DAE Schur identity for two PETSc shifts
# and 1/2 MPI ranks, preserving all boundary variables.
#
#   bash sandbox/run_constitutive_mumps_schur.sh 0.08
set -euo pipefail

mesh="${1:-0.08}"
if [[ $# -gt 1 || ! "$mesh" =~ ^0\.[0-9]+$ ]]; then
    echo "Usage: bash sandbox/run_constitutive_mumps_schur.sh [mesh_h=0.08]" >&2
    exit 2
fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/ph-constitutive-schur-h${mesh//./p}-$sha-XXXXXX")"
timeout_s="${SCRIMP_PH_SCHUR_TIMEOUT_S:-180}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] ||
   (( timeout_s < 30 || timeout_s > 1800 )); then
    echo "SCRIMP_PH_SCHUR_TIMEOUT_S must be 30..1800 seconds" >&2
    exit 2
fi
echo "SCRIMP: constant-coefficient constitutive MUMPS and exact Schur audit"
echo "Mesh h=$mesh; 2 reference CN steps; PETSc shifts=1000/2000"
echo "Reference full-system MUMPS remains unchanged. No solver modification."
echo "Result: $out"
sudo -v
for ranks in 1 2; do
    cpus="0"
    [[ "$ranks" == 2 ]] && cpus="0,2"
    echo "RUN $ranks MPI ranks / cpus=$cpus"
    sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        -e HOME=/tmp -e MPLBACKEND=Agg \
        -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
        -e SCRIMP_GETFEM_MATRIX_MODE=distributed \
        -e SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
        -e SCRIMP_MUMPS_RESET_POLICY=never \
        -e SCRIMP_MATRIX_MEMORY_MODE=compact \
        -e SCRIMP_BENCH_GIT_SHA="$sha" \
        --cpuset-cpus="$cpus" \
        -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
        -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
        -v "$root/examples:/workspace/scrimp/examples:ro" \
        -v "$out:/workspace/scrimp/outputs" \
        -v "$out:/artifacts" \
        scrimp-mpi:local \
        timeout -k 5s "${timeout_s}s" \
        mpiexec --bind-to none -n "$ranks" \
        python -u sandbox/probe_constitutive_mumps_schur.py \
            --mesh "$mesh" --dt 0.002 \
            --shift1 1000 --shift2 2000 \
            --output "/artifacts/schur_${ranks}r.json" \
            --petsc-log \
        > "$out/schur_${ranks}r.log" 2>&1 || {
            echo "FAILED: $ranks MPI rank(s)" >&2
            tail -n 90 "$out/schur_${ranks}r.log" >&2
            exit 1
        }
    grep -aE '^PASS sigma=|^PASS: full pH-DAE' \
        "$out/schur_${ranks}r.log"
    echo "  PETSc numerical LU factor events (reference + Q):"
    grep -aE '^MatLUFactorSym|^MatLUFactorNum|^MatSolve[[:space:]]' \
        "$out/schur_${ranks}r.log" | head -n 8 || true
done
echo "CHECK: full uncondensed DAE Schur identity and inverse-Q parity"
sudo docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp -e MPLBACKEND=Agg \
    --cpuset-cpus="0,2" \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$out:/artifacts" \
    scrimp-mpi:local \
    python -u sandbox/compare_constitutive_mumps_schur.py \
        --directory /artifacts
echo "SUCCESS: native 1/2-MPI MUMPS/constitutive Schur identity gates."
echo "Output: $out"
