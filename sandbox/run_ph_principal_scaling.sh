#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Strong scaling of the graph-selected principal DAE block H.
# Failures of exploratory distributed PCs are reported, not hidden.
set -euo pipefail
mesh="${1:-0.05}"
ranks_csv="${2:-1,2,4,6}"
repeat="${3:-4}"
if [[ ! "$mesh" =~ ^0\.[0-9]+$ || ! "$repeat" =~ ^[1-9][0-9]*$ || "$repeat" -gt 20 ]]; then
    echo "Usage: bash sandbox/run_ph_principal_scaling.sh [mesh=0.05] [ranks=1,2,4,6] [repeat=4]" >&2
    exit 2
fi
if [[ ! "$ranks_csv" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ]]; then
    echo "Ranks must be comma-separated positive integers" >&2
    exit 2
fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
host_workspace="${SCRIMP_DOCKER_WORKSPACE:-$HOME/Simulations/workspace}"
mkdir -p "$host_workspace" sandbox/outputs
host_workspace="$(realpath "$host_workspace")"
out="$(mktemp -d "$root/sandbox/outputs/ph-H-scaling-h${mesh//./p}-$sha-XXXXXX")"
timeout_s="${SCRIMP_H_SCALING_TIMEOUT_S:-360}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ || "$timeout_s" -lt 30 || "$timeout_s" -gt 3600 ]]; then
    echo "SCRIMP_H_SCALING_TIMEOUT_S must be 30..3600" >&2
    exit 2
fi

# Explicit physical CPU mapping. Override for other NUMA/affinity
# layouts: SCRIMP_SCALING_CPUS=0,2,4,6,8,10 (default).
cpu_csv="${SCRIMP_SCALING_CPUS:-0,2,4,6,8,10}"
IFS=, read -r -a cpus <<< "$cpu_csv"
IFS=, read -r -a ranks <<< "$ranks_csv"
echo "SCRIMP distributed principal H strong scaling"
echo "Source commit: $sha; mesh=$mesh; MPI ranks=$ranks_csv; repeats=$repeat"
echo "CPU affinity pool: $cpu_csv (override SCRIMP_SCALING_CPUS if needed)"
echo "Docker: scrimp-mpi:local; output: $out"
if ! git diff --quiet HEAD -- scrimp sandbox/probe_ph_principal_scaling.py; then
    echo "WARNING: working tree differs from HEAD"
fi
sudo -v
failed=0
for n in "${ranks[@]}"; do
    if (( n > ${#cpus[@]} )); then
        echo "SKIP ranks=$n: only ${#cpus[@]} CPU IDs configured" >&2
        failed=1
        continue
    fi
    affinity="$(IFS=,; echo "${cpus[*]:0:n}")"
    echo "RUN ranks=$n, pinned CPUs=$affinity"
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
        --cpuset-cpus="$affinity" \
        -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
        -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
        -v "$root/examples:/workspace/scrimp/examples:ro" \
        -v "$host_workspace:/scrimp-work" \
        -v "$out:/artifacts" \
        scrimp-mpi:local \
        timeout -k 5s "${timeout_s}s" \
        mpiexec --bind-to none -n "$n" \
        python -u /workspace/scrimp/sandbox/probe_ph_principal_scaling.py \
          --mesh "$mesh" --repeat "$repeat" \
          --output "/artifacts/H_scaling_${n}r.json" \
        > "$out/H_scaling_${n}r.log" 2>&1; then
        grep -aE '^(METHOD |Saved:)' "$out/H_scaling_${n}r.log" || true
    else
        echo "NON-PASS ranks=$n; see $out/H_scaling_${n}r.log" >&2
        tail -n 32 "$out/H_scaling_${n}r.log" >&2
        failed=1
    fi
done
echo "Reports: $out"
if (( failed )); then
    echo "At least one rank configuration failed. Preserve all JSON/logs." >&2
    exit 1
fi
echo "PASS: all requested MPI configurations produced H solver comparisons"
