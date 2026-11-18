#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.
#
# DIAGNOSIS, not another MUMPS ordering sweep:
#   4 vs 6 MPI, ICNTL(28)=2, MUMPS reset NEVER, actual CN wave.py,
#   stage-separated PETSc log: symbolic, numeric, triangular, assembly.
# All changes are sandbox-only, no production modification.
set -euo pipefail

mesh="${1:-0.01}"
ranks_csv="${2:-4,6}"
repeat="${3:-1}"
if [[ $# -gt 3 || ! "$mesh" =~ ^0\.[0-9]+$ ||
      ! "$ranks_csv" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ||
      ! "$repeat" =~ ^[1-3]$ ]]; then
    echo "Usage: bash sandbox/run_wave_mumps_4v6_diag.sh [mesh=0.01] [ranks=4,6] [repeat=1..3]" >&2
    exit 2
fi

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/wave-mumps-parallel-petsc-diag-h${mesh//./p}-$sha-XXXXXX")"
cpu_csv="${SCRIMP_SCALING_CPUS:-0,2,4,6,8,10}"
IFS=, read -r -a cpu_pool <<< "$cpu_csv"
IFS=, read -r -a ranks <<< "$ranks_csv"
timeout_s="${SCRIMP_WAVE_CN_TIMEOUT_S:-900}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] || (( timeout_s < 60 || timeout_s > 7200 )); then
    echo "SCRIMP_WAVE_CN_TIMEOUT_S must be 60..7200" >&2; exit 2
fi
for cpu in "${cpu_pool[@]}"; do
    if [[ ! "$cpu" =~ ^[0-9]+$ ]]; then
        echo "Invalid SCRIMP_SCALING_CPUS=$cpu_csv" >&2
        exit 2
    fi
done

echo "SCRIMP full CN / PETSc-MUMPS targeted strong-scaling diagnosis"
echo "Commit=$sha mesh=$mesh ranks=$ranks_csv repeat=$repeat"
echo "MUMPS ordering=parallel ICNTL(28)=2, reset=never ONLY"
echo "CPU affinity pool=$cpu_csv; MPI binding=none (unchanged from reference)"
echo "Output: $out"

# Record physical topology of precisely the CPUs in the test.
# Bind-to-none means all MPI processes can migrate within the cpuset;
# recording socket/core/SMT placement is essential before changing it.
{
    echo "Host $(hostname) / $(date -Is)"
    echo "CPU pool=$cpu_csv; mpiexec --bind-to none"
    if command -v lscpu >/dev/null 2>&1; then
        echo "=== lscpu summary ==="
        lscpu | grep -E '^(Architecture:|CPU\(s\):|Thread\(s\) per core:|Core\(s\) per socket:|Socket\(s\):|NUMA node\(s\):|Model name:)' || true
    fi
    echo "=== selected logical CPU topology ==="
    echo "logical_cpu socket core thread_siblings numa"
    for cpu in "${cpu_pool[@]}"; do
        base="/sys/devices/system/cpu/cpu${cpu}"
        socket="$(cat "$base/topology/physical_package_id" 2>/dev/null || echo unknown)"
        core="$(cat "$base/topology/core_id" 2>/dev/null || echo unknown)"
        siblings="$(cat "$base/topology/thread_siblings_list" 2>/dev/null || echo unknown)"
        nodes=( "$base"/node[0-9]* )
        node="unknown"
        if [[ -e "${nodes[0]}" ]]; then node="${nodes[0]##*/}"; fi
        echo "$cpu $socket $core $siblings $node"
    done
    echo "=== cgroup CPU quota when available ==="
    cat /sys/fs/cgroup/cpu.max 2>/dev/null || true
} > "$out/cpu_topology.txt"
cat "$out/cpu_topology.txt"

sudo -v
for rep in $(seq 1 "$repeat"); do
    for n in "${ranks[@]}"; do
        if (( n > ${#cpu_pool[@]} )); then
            echo "CPU pool too short for ranks=$n" >&2; exit 2
        fi
        affinity="$(IFS=,; echo "${cpu_pool[*]:0:n}")"
        label="parallel_${n}r_rep${rep}"
        echo "RUN ranks=$n rep=$rep CPUs=$affinity (PETSc native log ON)"
        if sudo docker run --rm \
            --user "$(id -u):$(id -g)" \
            --cpuset-cpus="$affinity" \
            --workdir /workspace/scrimp \
            -e HOME=/tmp -e MPLBACKEND=Agg \
            -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
            -e SCRIMP_GETFEM_MATRIX_MODE=distributed \
            -e SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
            -e SCRIMP_MATRIX_MEMORY_MODE=compact \
            -e SCRIMP_MUMPS_RESET_POLICY=never \
            -e SCRIMP_BENCH_GIT_SHA="$sha" \
            -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
            -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
            -v "$root/examples:/workspace/scrimp/examples:ro" \
            -v "$out:/workspace/scrimp/outputs" \
            -v "$out:/artifacts" \
            scrimp-mpi:local \
            timeout -k 10s "${timeout_s}s" \
            mpiexec --bind-to none -n "$n" \
            python -u sandbox/bench_wave_mumps_ordering.py \
                --mesh "$mesh" --dt 0.002 --final-time 0.02 \
                --save-dt 0.01 --series never --profile --petsc-log \
                --ordering parallel \
                --output "/artifacts/$label.json" \
            > "$out/$label.log" 2>&1; then
            grep -aE '^MAIN-TS-MUMPS |^INITIAL-TS-MUMPS |^Elapsed time:|^Steps:' "$out/$label.log" || true
            echo "  Native PETSc events (stage-aware):"
            grep -aE '^MatLUFactorSym|^MatLUFactorNum|^MatSolve |^PCSetUp ' "$out/$label.log" | head -25 || true
        else
            echo "NON-PASS ranks=$n rep=$rep -- $out/$label.log" >&2
            tail -n 60 "$out/$label.log" >&2
            echo "STOP after first failed test; no other configurations started" >&2
            echo "Reports: $out" >&2
            exit 1
        fi
    done
done
echo "Reports: $out"
echo "PASS: parallel MUMPS actual CN + stage-resolved PETSc diagnostics captured"
