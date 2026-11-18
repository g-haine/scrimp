#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Full 1-/2-rank physical parity of the production coenergy wave
# after the Concentric mixed-dimensional region fix.
#
# - unchanged source equations and global MUMPS;
# - fixed CN 5 steps, saves every full coupled DAE state;
# - independent GetFEM Hamiltonian integrals vs PETSc mass;
# - exact discrete midpoint loss = damping + boundary feedback;
# - full trajectory and physical parity across 1/2 ranks.
#
#   bash sandbox/run_coenergy_power_mpi.sh 0.12
set -euo pipefail
mesh="${1:-0.12}"
policy="${SCRIMP_COENERGY_MUMPS_POLICY:-never}"
if [[ "$policy" != "never" && "$policy" != "on_change" ]]; then
    echo "SCRIMP_COENERGY_MUMPS_POLICY must be never or on_change" >&2
    exit 2
fi
if [[ $# -gt 1 || ! "$mesh" =~ ^0\.[0-9]+$ ]]; then
    echo "Usage: bash sandbox/run_coenergy_power_mpi.sh [mesh_h=0.12]" >&2
    exit 2
fi

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
host_workspace="${SCRIMP_DOCKER_WORKSPACE:-$HOME/Simulations/workspace}"
mkdir -p "$host_workspace"
host_workspace="$(realpath "$host_workspace")"
[[ -w "$host_workspace" ]] || {
    echo "Docker work directory not writable: $host_workspace" >&2
    exit 2
}
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/coenergy-physics-h${mesh//./p}-$sha-XXXXXX")"
timeout_s="${SCRIMP_COENERGY_PHYSICS_TIMEOUT_S:-180}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] ||
   (( timeout_s < 30 || timeout_s > 1800 )); then
    echo "SCRIMP_COENERGY_PHYSICS_TIMEOUT_S must be 30..1800" >&2
    exit 2
fi
if ! git diff --quiet HEAD -- scrimp/utils/mesh.py; then
    echo "WARNING: mesh.py differs from Git HEAD $sha; results use mounted worktree."
    git status --short -- scrimp/utils/mesh.py
fi
echo "SCRIMP coenergy CN physical-power gate | commit=$sha | h=$mesh"
echo "5 fixed CN steps, same mesh, global PETSc MUMPS policy=$policy, 1/2 MPI ranks"
echo "Output: $out"
sudo -v
for ranks in 1 2; do
    cpus="0"
    [[ "$ranks" == 2 ]] && cpus="0,2"
    echo "RUN ranks=$ranks cpus=$cpus"
    sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        --workdir /scrimp-work \
        -e HOME=/tmp -e MPLBACKEND=Agg \
        -e PYTHONPATH=/workspace/scrimp \
        -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
        -e SCRIMP_GETFEM_MATRIX_MODE=distributed \
        -e SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
        -e SCRIMP_MATRIX_MEMORY_MODE=compact \
        -e SCRIMP_MUMPS_RESET_POLICY="$policy" \
        -e SCRIMP_BENCH_EXPECT_MUMPS_POLICY="$policy" \
        -e SCRIMP_TRAJECTORY_BACKEND=memory \
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
        python -u /workspace/scrimp/sandbox/verify_coenergy_power_mpi.py \
            --mesh "$mesh" --dt 0.002 --steps 5 \
            --output-dir /artifacts \
        > "$out/physics_${ranks}r.log" 2>&1 || {
            echo "FAILED: native physics gate on $ranks MPI ranks" >&2
            tail -n 90 "$out/physics_${ranks}r.log" >&2
            exit 1
        }
    grep -aE '^PASS coenergy physics' "$out/physics_${ranks}r.log"
done
echo "CHECK: complete 1/2 MPI states, Hamiltonian and port-power time series"
sudo docker run --rm \
    --user "$(id -u):$(id -g)" \
    --workdir /scrimp-work \
    -e HOME=/tmp -e MPLBACKEND=Agg \
    --cpuset-cpus="0,2" \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$host_workspace:/scrimp-work" \
    -v "$out:/artifacts" \
    scrimp-mpi:local \
    python -u /workspace/scrimp/sandbox/compare_coenergy_power_mpi.py \
        --directory /artifacts
echo "PASS: coenergy physical and MPI parity gates. Output: $out"
