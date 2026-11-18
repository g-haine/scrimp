#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.
#
# Experimental distributed Schur+ASM on complete manufactured pH-DAE.
# Never runs production iterative solver. MUMPS baseline remains reference.
set -euo pipefail
mesh="${1:-0.05}"
ranks_csv="${2:-1,2}"
repeat_rhs="${3:-2}"
schur_ksp="${4:-preonly}"
if [[ $# -gt 4 || ! "$schur_ksp" =~ ^(preonly|gmres)$ || ! "$mesh" =~ ^0\.[0-9]+$ ||
      ! "$ranks_csv" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ||
      ! "$repeat_rhs" =~ ^[1-9][0-9]*$ || "$repeat_rhs" -gt 20 ]]; then
    echo "Usage: bash sandbox/run_ph_distributed_dae.sh [mesh=0.05] [ranks=1,2] [repeat_rhs=1..20] [schur_ksp=preonly|gmres]" >&2
    exit 2
fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
host_workspace="${SCRIMP_DOCKER_WORKSPACE:-$HOME/Simulations/workspace}"
mkdir -p "$host_workspace" sandbox/outputs
host_workspace="$(realpath "$host_workspace")"
out="$(mktemp -d "$root/sandbox/outputs/ph-distributed-dae-${schur_ksp}-h${mesh//./p}-$sha-XXXXXX")"
timeout_s="${SCRIMP_DISTRIBUTED_DAE_TIMEOUT_S:-480}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ || "$timeout_s" -lt 30 || "$timeout_s" -gt 3600 ]]; then
    echo "SCRIMP_DISTRIBUTED_DAE_TIMEOUT_S must be 30..3600" >&2
    exit 2
fi
cpu_csv="${SCRIMP_SCALING_CPUS:-0,2,4,6,8,10}"
IFS=, read -r -a cpus <<< "$cpu_csv"
IFS=, read -r -a ranks <<< "$ranks_csv"
echo "SCRIMP full manufactured pH-DAE distributed Schur+ASM qualification"
echo "Source commit: $sha; mesh=$mesh; ranks=$ranks_csv; repeat_rhs=$repeat_rhs; schur_ksp=$schur_ksp"
echo "CPU pool: $cpu_csv; output: $out"
echo "No production solver change. Full MUMPS remains reference only."
sudo -v
failed=0
for n in "${ranks[@]}"; do
    if (( n > ${#cpus[@]} )); then
        echo "SKIP ranks=$n: CPU affinity list has ${#cpus[@]} entries" >&2
        failed=1
        continue
    fi
    affinity="$(IFS=,; echo "${cpus[*]:0:n}")"
    echo "RUN ranks=$n / CPUs=$affinity"
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
        python -u /workspace/scrimp/sandbox/probe_ph_distributed_dae.py \
        --mesh "$mesh" --schur-pre path2 --h-pc asm --schur-ksp "$schur_ksp" \
        --h-rtol "${SCRIMP_DISTRIBUTED_H_RTOL:-1e-10}" \
        --h-max-it "${SCRIMP_DISTRIBUTED_H_MAX_IT:-350}" \
        --repeat-rhs "$repeat_rhs" \
        --output "/artifacts/distributed_dae_${n}r.json" \
        > "$out/distributed_dae_${n}r.log" 2>&1; then
        grep -aE '^(PASS|FAIL|REJECT|NESTED|REUSE|MUMPS-REUSE) sigma=|^PATH2 built:|^Saved:' "$out/distributed_dae_${n}r.log" || true
    else
        echo "NON-PASS ranks=$n: $out/distributed_dae_${n}r.log" >&2
        tail -n 45 "$out/distributed_dae_${n}r.log" >&2
        failed=1
    fi
done
echo "Reports: $out"
if (( failed )); then
    echo "NON-PASS: inspect logs for convergence failures; do not relax accuracy gates" >&2
    exit 1
fi
echo "PASS: distributed full pH-DAE manufactured solves qualified on selected ranks"
