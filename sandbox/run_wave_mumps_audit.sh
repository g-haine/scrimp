#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.
# MUMPS factorization and triangular solve diagnostic, full wave Jacobian.
set -euo pipefail
mesh="${1:-0.025}"
ranks_csv="${2:-1,2,4,6}"
orderings_csv="${3:-default,sequential,parallel}"
repeat="${4:-4}"
if [[ ! "$mesh" =~ ^0\.[0-9]+$ ||
      ! "$ranks_csv" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ||
      ! "$orderings_csv" =~ ^(default|sequential|parallel)(,(default|sequential|parallel))*$ ||
      ! "$repeat" =~ ^[0-9]+$ ]]; then
    echo "Usage: bash sandbox/run_wave_mumps_audit.sh [h=0.025] [ranks=1,2,4,6] [orderings=default,sequential,parallel] [repeat=4]" >&2
    exit 2
fi
if (( repeat < 2 || repeat > 12 )); then
    echo "repeat must be 2..12" >&2
    exit 2
fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
workdir="${SCRIMP_DOCKER_WORKSPACE:-$HOME/Simulations/workspace}"
mkdir -p "$workdir" sandbox/outputs
workdir="$(realpath "$workdir")"
out="$(mktemp -d "$root/sandbox/outputs/wave-mumps-audit-h${mesh//./p}-$sha-XXXXXX")"
timeout_s="${SCRIMP_MUMPS_AUDIT_TIMEOUT_S:-600}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ || "$timeout_s" -lt 30 || "$timeout_s" -gt 7200 ]]; then
    echo "SCRIMP_MUMPS_AUDIT_TIMEOUT_S must be 30..7200" >&2
    exit 2
fi
cpu_csv="${SCRIMP_SCALING_CPUS:-0,2,4,6,8,10}"
IFS=, read -r -a cpu_pool <<< "$cpu_csv"
IFS=, read -r -a ranks <<< "$ranks_csv"
IFS=, read -r -a orderings <<< "$orderings_csv"
echo "SCRIMP wave full Jacobian MUMPS audit"
echo "Commit: $sha; mesh=$mesh; ranks=$ranks_csv; orderings=$orderings_csv; repeat=$repeat"
echo "CPU affinity pool=$cpu_csv; container=scrimp-mpi:local"
echo "Output: $out"
echo "Baseline operator: distributed COO; no changes to production solver"
sudo -v
failed=0
for n in "${ranks[@]}"; do
    if (( n > ${#cpu_pool[@]} )); then
        echo "SKIP ranks=$n: CPU pool insufficient" >&2
        failed=1
        continue
    fi
    cpus="$(IFS=,; echo "${cpu_pool[*]:0:n}")"
    for ordering in "${orderings[@]}"; do
        label="${ordering}_${n}r"
        echo "RUN ranks=$n ordering=$ordering CPUs=$cpus"
        if sudo docker run --rm \
            --user "$(id -u):$(id -g)" \
            --workdir /scrimp-work \
            --cpuset-cpus="$cpus" \
            -e HOME=/tmp -e MPLBACKEND=Agg \
            -e PYTHONPATH=/workspace/scrimp \
            -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
            -e SCRIMP_GETFEM_MATRIX_MODE=distributed \
            -e SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
            -e SCRIMP_MATRIX_MEMORY_MODE=compact \
            -e SCRIMP_MUMPS_RESET_POLICY=never \
            -e SCRIMP_BENCH_GIT_SHA="$sha" \
            -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
            -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
            -v "$root/examples:/workspace/scrimp/examples:ro" \
            -v "$workdir:/scrimp-work" \
            -v "$out:/artifacts" \
            scrimp-mpi:local \
            timeout -k 5s "${timeout_s}s" \
            mpiexec --bind-to none -n "$n" \
            python -u /workspace/scrimp/sandbox/probe_wave_mumps_audit.py \
              --mesh "$mesh" --ordering "$ordering" --repeat "$repeat" \
              --output "/artifacts/$label.json" \
              > "$out/$label.log" 2>&1; then
            grep -aE '^(MUMPS-AUDIT |Saved:)' "$out/$label.log" || true
        else
            echo "NON-PASS ranks=$n ordering=$ordering: $out/$label.log" >&2
            tail -n 40 "$out/$label.log" >&2
            failed=1
        fi
    done
done
echo "Reports: $out"
if (( failed )); then
    echo "Some configurations failed. Parallel ordering may be absent in MUMPS build." >&2
    exit 1
fi
echo "PASS: complete wave MUMPS factor/setup/solve audit finished"
