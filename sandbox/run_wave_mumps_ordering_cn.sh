#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.
#
# Full CN trajectory comparison using ONLY MUMPS RESET POLICY NEVER.
# Default MUMPS vs explicitly verified ICNTL(28)=2 inside main TS.
set -euo pipefail
mesh="${1:-0.005}"
ranks_csv="${2:-1,2,4,6}"
repeats="${3:-1}"
if [[ $# -gt 3 || ! "$mesh" =~ ^0\.[0-9]+$ ||
      ! "$ranks_csv" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ||
      ! "$repeats" =~ ^[0-9]+$ ]]; then
    echo "Usage: bash sandbox/run_wave_mumps_ordering_cn.sh [mesh=0.005] [ranks=1,2,4,6] [repeats=1..5]" >&2
    exit 2
fi
if (( repeats < 1 || repeats > 5 )); then echo "repeats must be 1..5" >&2; exit 2; fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/wave-cn-mumps-ordering-never-h${mesh//./p}-$sha-XXXXXX")"
cpu_csv="${SCRIMP_SCALING_CPUS:-0,2,4,6,8,10}"
IFS=, read -r -a cpus <<< "$cpu_csv"
IFS=, read -r -a ranks <<< "$ranks_csv"
timeout_s="${SCRIMP_WAVE_CN_TIMEOUT_S:-900}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] || (( timeout_s < 60 || timeout_s > 7200 )); then
    echo "SCRIMP_WAVE_CN_TIMEOUT_S must be 60..7200" >&2; exit 2
fi
echo "SCRIMP full CN wave MUMPS default vs parallel ordering"
echo "Commit=$sha mesh=$mesh ranks=$ranks_csv repeats=$repeats MUMPS reset policy=never ONLY"
echo "Output: $out"
sudo -v
failed=0
for rep in $(seq 1 "$repeats"); do
  for n in "${ranks[@]}"; do
    if (( n > ${#cpus[@]} )); then echo "Missing CPUs for $n ranks" >&2; exit 2; fi
    cpuset="$(IFS=,; echo "${cpus[*]:0:n}")"
    # Alternate ordering to reduce systematic thermal/cache ordering bias.
    if (( rep % 2 )); then orderings=(default parallel); else orderings=(parallel default); fi
    for ordering in "${orderings[@]}"; do
      label="wave_${ordering}_${n}r_rep${rep}"
      echo "RUN ranks=$n ordering=$ordering rep=$rep cpus=$cpuset policy=never"
      if sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        --cpuset-cpus="$cpuset" \
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
        -v "$out:/artifacts" \
        scrimp-mpi:local \
        timeout -k 10s "${timeout_s}s" \
        mpiexec --bind-to none -n "$n" \
        python -u sandbox/bench_wave_mumps_ordering.py \
            --mesh "$mesh" --dt 0.002 --final-time 0.02 \
            --save-dt 0.01 --series never --profile \
            --ordering "$ordering" \
            --output "/artifacts/$label.json" \
            > "$out/$label.log" 2>&1; then
          grep -aE '^MAIN-TS-MUMPS |^Elapsed time:|^Steps:' "$out/$label.log" || true
        else
          echo "NON-PASS $label -- $out/$label.log" >&2
          tail -n 50 "$out/$label.log" >&2
          failed=1
        fi
    done
  done
done
echo "Reports: $out"
if (( failed )); then echo "NON-PASS: read logs, do not relax numerical gates" >&2; exit 1; fi
echo "PASS: all full CN runs verified actual MUMPS ordering with policy NEVER"
