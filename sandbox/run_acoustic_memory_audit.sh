#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Native MPI memory audit, identical 2D acoustics with two saving cadences.
# The original SCRIMP DPHS.monitor remains active; only measurements are
# added outside its normal behavior. No host Conda simulation.
#
#   bash sandbox/run_acoustic_memory_audit.sh
#   bash sandbox/run_acoustic_memory_audit.sh 0.04 0.002 0.5 2
#
# Default: h=.04, dt=.002, tf=.5 (250 CN steps), 2 MPI ranks,
# COO default, distributed GetFEM, MUMPS reset=never.
# 2 runs: initial/final-only vs every accepted step.
#
# Outputs are versioned by source SHA and physical simulation settings.

set -euo pipefail
mesh="${1:-0.04}"
dt="${2:-0.002}"
final="${3:-0.5}"
ranks="${4:-2}"

if [[ ! "$mesh" =~ ^0\.[0-9]+$ ]] ||
   [[ ! "$dt" =~ ^0\.[0-9]+$ ]] ||
   [[ ! "$final" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
   [[ "$ranks" != 1 && "$ranks" != 2 && "$ranks" != 4 && "$ranks" != 6 ]]; then
    echo "Usage: bash sandbox/run_acoustic_memory_audit.sh 0.04 0.002 0.5 [1|2|4|6]" >&2
    exit 2
fi

timeout_s="${SCRIMP_MEMORY_TIMEOUT_S:-900}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] ||
   (( timeout_s < 60 || timeout_s > 7200 )); then
    echo "SCRIMP_MEMORY_TIMEOUT_S must be 60..7200 seconds" >&2
    exit 2
fi

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
tag="h${mesh//./p}-dt${dt//./p}-T${final//./p}-${ranks}r"
out="$root/sandbox/outputs/acoustic-memory-${tag}-${sha}"
mkdir -p "$out"

pcores=(0 2 4 6 8 10)
cpus="$(IFS=,; echo "${pcores[*]:0:ranks}")"
if ! sudo -v; then
    echo "sudo credentials unavailable" >&2
    exit 1
fi

echo "SCRIMP memory audit | commit=$sha | h=$mesh | dt=$dt | t_f=$final | MPI=$ranks | cpus=$cpus"
echo "Distributed GetFEM + default COO; MUMPS never; results: $out"

echo "PREFLIGHT: check SCRIMP import/write permissions and report analyzer"
if ! sudo docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp \
    -e MPLBACKEND=Agg \
    --cpuset-cpus="0,2" \
    -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$root/tests:/workspace/scrimp/tests:ro" \
    -v "$out:/workspace/scrimp/outputs" \
    -v "$out:/artifacts" \
    scrimp-mpi:local \
    bash -c 'python -c "import scrimp" && python -m unittest -v tests.test_acoustic_memory_analysis tests.test_mpi_algebra.TestMPIAlgebra.test_final_solution_filter_avoids_multiple_copies_with_exact_values' \
    >"$out/preflight.log" 2>&1; then
    echo "FAILED report-analyzer preflight; inspect $out/preflight.log" >&2
    tail -n 80 "$out/preflight.log" >&2
    exit 1
fi
echo "  PASS: SCRIMP import + writable output directory + analyzer checks"

for mode in sparse dense; do
    logfile="$out/${mode}.log"
    echo "RUN: $mode saved-state frequency"
    if ! sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        -e HOME=/tmp \
        -e MPLBACKEND=Agg \
        --cpuset-cpus="$cpus" \
        -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
        -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
        -v "$root/examples:/workspace/scrimp/examples:ro" \
        -v "$out:/workspace/scrimp/outputs" \
        -v "$out:/artifacts" \
        scrimp-mpi:local \
        timeout -k 10s "${timeout_s}s" \
        env -u SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN \
            -u SCRIMP_GETFEM_BRIDGE_DISTRIBUTED \
            SCRIMP_GETFEM_MATRIX_MODE=distributed \
            SCRIMP_MUMPS_RESET_POLICY=never \
            SCRIMP_BENCH_GIT_SHA="$sha" \
        mpiexec --bind-to none -n "$ranks" \
        python -u sandbox/probe_acoustic_memory.py \
            --mesh "$mesh" --dt "$dt" --final-time "$final" \
            --mode "$mode" --sample-every 10 \
            --output-dir /artifacts \
        >"$logfile" 2>&1; then
        echo "FAILED $mode: inspect $logfile" >&2
        tail -n 80 "$logfile" >&2
        exit 1
    fi
    for (( rank=0; rank<ranks; rank++ )); do
        for suffix in ".json" "_state.npz"; do
            f="$out/memory_${mode}_rank${rank}${suffix}"
            if [[ ! -s "$f" ]]; then
                echo "Missing required rank-local report: $f" >&2
                exit 1
            fi
        done
    done
    echo "  PASS: $mode completed; all $ranks MPI ranks wrote reports"
done

echo "ANALYZE: strict whole-state parity and memory traces"
if ! sudo docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp \
    -e MPLBACKEND=Agg \
    --cpuset-cpus="0,2" \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$out:/workspace/scrimp/outputs" \
    -v "$out:/artifacts" \
    scrimp-mpi:local \
    python -u sandbox/analyze_acoustic_memory.py \
        --directory /artifacts \
        --prefix /artifacts/acoustic_memory \
        >"$out/analyze.log" 2>&1; then
    echo "FAILED analysis; inspect $out/analyze.log" >&2
    tail -n 80 "$out/analyze.log" >&2
    exit 1
fi
cat "$out/analyze.log"
echo "Files:"
ls -lh "$out"/acoustic_memory_*.{png,csv,json}
echo "SUCCESS: memory audit and initial/final full-state parity passed."
echo "Results: $out"
