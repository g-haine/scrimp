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
#   bash sandbox/run_ph_block_audit.sh 0.12 wave_coenergy
set -euo pipefail

mesh="${1:-0.08}"
model="${2:-wave}"
if [[ $# -gt 2 || ! "$mesh" =~ ^0\.[0-9]+$ ||
      ( "$model" != "wave" && "$model" != "wave_coenergy" ) ]]; then
    echo "Usage: bash sandbox/run_ph_block_audit.sh [mesh_h=0.08] [wave|wave_coenergy]" >&2
    exit 2
fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
# SCRIMP creates auxiliary outputs under the Python cwd at import.
# Use the host user's own workspace as container cwd, not the
# image-owned source tree. Sources stay mounted read-only.
host_workspace="${SCRIMP_DOCKER_WORKSPACE:-$HOME/Simulations/workspace}"
mkdir -p "$host_workspace"
host_workspace="$(realpath "$host_workspace")"
if [[ ! -d "$host_workspace" || ! -w "$host_workspace" ]]; then
    echo "Docker host work directory is not writable: $host_workspace" >&2
    exit 2
fi
sha="$(git rev-parse --short=8 HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/ph-block-audit-${model}-h${mesh//./p}-$sha-XXXXXX")"
timeout_s="${SCRIMP_PH_BLOCK_AUDIT_TIMEOUT_S:-180}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] ||
   (( timeout_s < 30 || timeout_s > 1800 )); then
    echo "SCRIMP_PH_BLOCK_AUDIT_TIMEOUT_S must be 30..1800 seconds" >&2
    exit 2
fi
echo "SCRIMP pH-DAE field-block audit | source $sha | model=$model | mesh h=$mesh"
echo "2 accepted CN steps, 1 and 2 MPI ranks, reference MUMPS unchanged"
echo "User-owned host workdir: $host_workspace (container: /scrimp-work)"
echo "Output: $out"
# With this opt-in flag, snapshot the TRUE GetFEM matrix on each
# rank before it is transferred, plus final PETSc owned-row CSR.
# Distinguish a rank-dependent native assembly from a bridge error.
capture="${SCRIMP_PH_CAPTURE_MATRICES:-0}"
if [[ "$capture" != "0" && "$capture" != "1" ]]; then
    echo "SCRIMP_PH_CAPTURE_MATRICES must be 0 or 1" >&2
    exit 2
fi
capture_args=()
if [[ "$capture" == "1" ]]; then
    capture_args+=(--capture-matrices)
    echo "Raw GetFEM + PETSc sparse matrix capture enabled."
fi
sudo -v
for ranks in 1 2; do
    cpus="0"
    [[ "$ranks" == 2 ]] && cpus="0,2"
    echo "RUN rank=$ranks, cpus=$cpus"
    sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        --workdir /scrimp-work \
        -e HOME=/tmp -e MPLBACKEND=Agg \
        -e PYTHONPATH=/workspace/scrimp \
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
        -v "$host_workspace:/scrimp-work" \
        -v "$out:/artifacts" \
        scrimp-mpi:local \
        timeout -k 5s "${timeout_s}s" \
        mpiexec --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/audit_ph_blocks.py \
            --model "$model" --mesh "$mesh" --dt 0.002 --steps 2 \
            --output "/artifacts/blocks_${ranks}r.json" \
            "${capture_args[@]}" \
        > "$out/blocks_${ranks}r.log" 2>&1 || {
            echo "FAILED: $ranks rank(s)" >&2
            tail -n 90 "$out/blocks_${ranks}r.log" >&2
            exit 1
        }
    tail -n 24 "$out/blocks_${ranks}r.log"
done
if [[ "$capture" == "1" ]]; then
    echo "LOCALIZE: compare raw GetFEM matrices with final PETSc COO values"
    sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        --workdir /scrimp-work \
        -e HOME=/tmp -e MPLBACKEND=Agg \
        --cpuset-cpus="0,2" \
        -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
        -v "$host_workspace:/scrimp-work" \
        -v "$out:/artifacts" \
        scrimp-mpi:local \
        python -u /workspace/scrimp/sandbox/compare_raw_getfem_petsc.py \
            /artifacts
fi
echo "DIAGNOSE: per-variable block norm differences (offline; no extra solves)"
python3 sandbox/diagnose_ph_block_parity.py "$out"
echo "CHECK: 1/2-rank block Frobenius-norm parity (NOT entrywise equality)"
sudo docker run --rm \
    --user "$(id -u):$(id -g)" \
    --workdir /scrimp-work \
    -e HOME=/tmp -e MPLBACKEND=Agg \
    --cpuset-cpus="0,2" \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$host_workspace:/scrimp-work" \
    -v "$out:/artifacts" \
    scrimp-mpi:local \
    python -u /workspace/scrimp/sandbox/compare_ph_block_audits.py \
        --directory /artifacts
echo "PASS: GetFEM/PETSc pH block norm maps agree on 1/2 MPI ranks (not an entrywise operator proof)."
echo "JSON: $out/blocks_{1r,2r}.json"
echo "Output: $out"
