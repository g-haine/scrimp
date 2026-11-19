#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Short, UNMODIFIED PETSc/MUMPS baseline diagnostics for scaling audit.
# Count MatLUFactorNumeric vs MatSolve to discover whether fixed-CN
# acoustic waves are refactorized on every step despite MUMPS "never".
# Only FIVE CN steps; not a performance benchmark. Use matched pinning
# with the 96-run dataset and examine core topology separately.
#
#   bash sandbox/run_coo_factorization_probe.sh 0.02
#   bash sandbox/run_coo_factorization_probe.sh 0.01
set -euo pipefail

mesh="${1:-0.02}"
if [[ $# -gt 1 || ! "$mesh" =~ ^0\.[0-9]+$ ]]; then
    echo "Usage: bash sandbox/run_coo_factorization_probe.sh [mesh_h=0.02]" >&2
    exit 2
fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/mumps-factor-probe-h${mesh//./p}-$sha-XXXXXX")"
timeout_s="${SCRIMP_FACTOR_PROBE_TIMEOUT_S:-180}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] || ((timeout_s < 30 || timeout_s > 1800)); then
    echo "SCRIMP_FACTOR_PROBE_TIMEOUT_S must be in [30,1800]" >&2
    exit 2
fi
echo "SCRIMP MUMPS factorization profiling: COO/never, h=$mesh, five CN steps"
echo "Output: $out"
echo "This run only instruments the current solver; no physics changes."
echo "Inspect physical CPU topology before interpreting rank speedup:"
lscpu -e=CPU,CORE,SOCKET,NODE,ONLINE > "$out/cpu_topology.txt" || true
pcores=(0 2 4 6 8 10)
sudo -v
for n in 1 4 6; do
    cpus="$(IFS=,; echo "${pcores[*]:0:n}")"
    echo "PETSc Log: $n MPI rank(s), CPUs $cpus"
    sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        -e HOME=/tmp -e MPLBACKEND=Agg \
        -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
        -e SCRIMP_GETFEM_MATRIX_MODE=distributed \
        -e SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
        -e SCRIMP_MUMPS_RESET_POLICY=never \
        -e SCRIMP_MATRIX_MEMORY_MODE=compact \
        -e SCRIMP_PROFILE=1 \
        -e SCRIMP_BENCH_GIT_SHA="$sha" \
        --cpuset-cpus="$cpus" \
        -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
        -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
        -v "$root/examples:/workspace/scrimp/examples:ro" \
        -v "$out:/workspace/scrimp/outputs" \
        -v "$out:/artifacts" \
        scrimp-mpi:local \
        timeout -k 5s "${timeout_s}s" \
        mpiexec --bind-to none -n "$n" \
        python -u sandbox/bench_acoustic_2d.py \
        --mesh "$mesh" --dt 0.002 --final-time 0.01 \
        --save-dt 0.01 --series never \
        --profile --petsc-log \
        --output "/artifacts/petsc_${n}r.json" \
        > "$out/petsc_${n}r.log" 2>&1 || {
            echo "FAILED: inspect $out/petsc_${n}r.log" >&2
            tail -n 80 "$out/petsc_${n}r.log" >&2
            exit 1
        }
    echo "  PETSc events:"
    grep -aE '^MatLUFactorSym|^MatLUFactorNum|^MatSolve[[:space:]]|^PCSetUp[[:space:]]|^KSPSolve[[:space:]]|^MatAssemblyEnd' \
        "$out/petsc_${n}r.log" || {
            echo "PETSc log unavailable for $n ranks" >&2
            exit 1
        }
done
echo "PASS: 1/4/6-rank PETSc event logs and SCRIMP phase profiles created."
echo "Compare the NUMBER of numeric LU factorizations with 5 CN steps."
echo "Numeric LU count near 6 (initial + each main step) suggests repeated"
echo "factorization; count near 2 (initial + first main) suggests reuse."
echo "A combined event count does not by itself separate initial/main TS."
echo "Results: $out"
