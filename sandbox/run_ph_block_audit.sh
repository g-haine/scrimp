#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Read-only phase-1 Schur PRECONDITIONER DESIGN audit:
#   * native GetFEM field offsets including boundary variables;
#   * PETSc M, K and last J actual variable-by-variable couplings;
#   * 1/2 MPI rank comparison.
#   * two CN steps JUST to assemble the actual reference PETSc matrices.
#
# This script changes NO solver, timestep algorithm, residual, physics
# or postprocessing and does NOT benchmark scalability.
#
#   bash sandbox/run_ph_block_audit.sh 0.08
set -euo pipefail

mesh="${1:-0.08}"
if [[ $# -gt 1 || ! "$mesh" =~ ^0\.[0-9]+$ ]]; then
    echo "Usage: bash sandbox/run_ph_block_audit.sh [mesh_h=0.08]" >&2
    exit 2
fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/ph-block-audit-h${mesh//./p}-$sha-XXXXXX")"
timeout_s="${SCRIMP_PH_BLOCK_AUDIT_TIMEOUT_S:-180}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] ||
   (( timeout_s < 30 || timeout_s > 1800 )); then
    echo "SCRIMP_PH_BLOCK_AUDIT_TIMEOUT_S must be 30..1800 seconds" >&2
    exit 2
fi
echo "SCRIMP pH-DAE field-block audit | source $sha | mesh h=$mesh"
echo "2 accepted CN steps, 1 and 2 MPI ranks, reference MUMPS unchanged"
echo "Output: $out"
sudo -v
for ranks in 1 2; do
    cpus="0"
    [[ "$ranks" == 2 ]] && cpus="0,2"
    echo "RUN rank=$ranks, cpus=$cpus"
    sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        -e HOME=/tmp -e MPLBACKEND=Agg \
        -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
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
        python -u sandbox/audit_ph_blocks.py \
            --mesh "$mesh" --dt 0.002 --steps 2 \
            --output "/artifacts/blocks_${ranks}r.json" \
        > "$out/blocks_${ranks}r.log" 2>&1 || {
            echo "FAILED: $ranks rank(s)" >&2
            tail -n 90 "$out/blocks_${ranks}r.log" >&2
            exit 1
        }
    tail -n 24 "$out/blocks_${ranks}r.log"
done
echo "CHECK: all field/block coefficients identical between 1 and 2 ranks"
sudo docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp -e MPLBACKEND=Agg \
    --cpuset-cpus="0,2" \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$out:/artifacts" \
    scrimp-mpi:local \
    python -u sandbox/compare_ph_block_audits.py \
        --directory /artifacts
echo "PASS: actual GetFEM/PETSc pH blocks classified consistently on 1/2 MPI ranks."
echo "JSON: $out/blocks_{1r,2r}.json"
echo "Output: $out"
