#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Current SCRIMP ONLY: full strong-scaling qualification at h=0.005
# (~3 million DOFs), 1/2/4/6 MPI ranks, 3 full repetitions.
# No pilot, no v1.2.0, no worktree, no release source anywhere.
# Physical wave.py, CN dt=.002 t_f=.02, global MUMPS and exact
# on_change operator-reuse policy unchanged. Human PNG/SVG figures.
set -euo pipefail

mesh="${1:-0.005}"
repetitions="${2:-3}"
ranks_csv="${SCRIMP_V130_LARGE_RANKS:-1,2,4,6}"
if [[ $# -gt 2 || ! "$mesh" =~ ^0\.[0-9]+$ ||
      ! "$repetitions" =~ ^[1-5]$ ||
      ! "$ranks_csv" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ]]; then
    echo "Usage: bash sandbox/run_v130_large_wave_scaling.sh [mesh=0.005] [repetitions=3]" >&2
    echo "Optional: SCRIMP_V130_LARGE_RANKS=1,2,4,6; SCRIMP_V130_LARGE_TIMEOUT_S=7200" >&2
    exit 2
fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
mkdir -p "$root/sandbox/outputs"
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-large-wave-h${mesh//./p}-${sha:0:8}-XXXXXX")"
IFS=, read -r -a ranks <<< "$ranks_csv"
IFS=, read -r -a cpu_pool <<< "${SCRIMP_SCALING_CPUS:-0,2,4,6,8,10}"
timeout_s="${SCRIMP_V130_LARGE_TIMEOUT_S:-7200}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] ||
   (( timeout_s < 300 || timeout_s > 14400 )); then
    echo "SCRIMP_V130_LARGE_TIMEOUT_S must be 300..14400 seconds" >&2
    exit 2
fi
for n in "${ranks[@]}"; do
    if (( n > ${#cpu_pool[@]} )); then
        echo "Insufficient CPU affinity entries for $n ranks" >&2
        exit 2
    fi
done
echo "SCRIMP v1.3.0 current-only large acoustic strong scaling"
echo "Commit=$sha mesh=$mesh MPI=$ranks_csv repetitions=$repetitions"
echo "Time=CN dt=0.002 t_f=0.02, 10 accepted steps"
echo "GetFEM distributed COO, compact matrices, on_change MUMPS, MPI ordering"
echo "NO release v1.2.0; NO pilot; direct complete campaign"
echo "Output: $out"
printf '{"schema":"scrimp-v130-current-large-wave-v1","current_sha":"%s","mesh":%s,"repetitions":%s,"ranks":"%s","dt":0.002,"final_time":0.02,"steps":10,"source":"examples/wave.py","variant":"current","mumps_policy":"on_change"}\n' \
  "$sha" "$mesh" "$repetitions" "$ranks_csv" > "$out/manifest.json"
sudo -v

# Interleave ranks between repetitions, as in the original scaling
# protocol, rather than running all 1-rank cases before all 6-rank.
for rep in $(seq 1 "$repetitions"); do
    for n in "${ranks[@]}"; do
        label="current_on_change_${n}r_rep${rep}"
        cpus="$(IFS=,; echo "${cpu_pool[*]:0:n}")"
        echo "RUN $label cpu_affinity=$cpus mesh=$mesh"
        if ! sudo docker run --rm \
            --user "$(id -u):$(id -g)" \
            --cpuset-cpus="$cpus" --workdir /workspace/scrimp \
            -e HOME=/tmp -e MPLBACKEND=Agg \
            -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
            -e SCRIMP_GETFEM_MATRIX_MODE=distributed \
            -e SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
            -e SCRIMP_MATRIX_MEMORY_MODE=compact \
            -e SCRIMP_TRAJECTORY_BACKEND=memory \
            -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
            -v "$root/examples:/workspace/scrimp/examples:ro" \
            -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
            -v "$out:/workspace/scrimp/outputs" \
            -v "$out:/artifacts" \
            scrimp-mpi:local timeout -k 10s "${timeout_s}s" \
            env -u SCRIMP_MUMPS_RESET_POLICY \
                -u SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN \
                -u SCRIMP_NONLINEAR_RESIDUAL_FORM \
            mpiexec --bind-to none -n "$n" \
            python -u sandbox/bench_wave_release_scaling.py \
                --variant current --policy on_change --mesh "$mesh" \
                --dt 0.002 --final-time 0.02 --save-dt 0.01 \
                --source-sha "$sha" --repetition "$rep" \
                --profile --petsc-log \
                --out "/artifacts/${label}.json" \
            > "$out/${label}.log" 2>&1; then
            echo "FAIL: $label (see $out/${label}.log)" >&2
            tail -n 90 "$out/${label}.log" >&2
            echo "Existing successful runs are preserved in $out" >&2
            exit 1
        fi
        if [[ ! -s "$out/${label}.json" ||
              ! -s "$out/${label}_state.npz" ||
              ! -s "$out/${label}.log" ]]; then
            echo "FAIL: Missing report, complete state or native PETSc log for $label" >&2
            exit 1
        fi
        grep -aE '^PASS current/on_change/' "$out/${label}.log" | tail -1 || true
    done
done

echo "ANALYZE: exact MPI full-state parity, PETSc LU events, 3-repeat distributions"
sudo docker run --rm \
    --user "$(id -u):$(id -g)" --workdir /workspace/scrimp \
    -e HOME=/tmp -e MPLBACKEND=Agg \
    -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$out:/workspace/scrimp/outputs" \
    -v "$out:/artifacts" \
    scrimp-mpi:local \
    python -u sandbox/analyze_v130_large_wave_scaling.py \
        --directory /artifacts --repetitions "$repetitions" \
        --ranks "$ranks_csv"
echo "PASS: current-only large wave native MPI qualification"
echo "PNG figures: $out/png"
echo "SVG vector figures: $out/svg"
echo "Reports: $out"
