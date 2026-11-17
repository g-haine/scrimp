#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# ANALYSIS ONLY: reuse completed native runs from --mpi-only or a
# historical-state mismatch. No new PDE solves, no mpiexec, no
# fabricated release gains. Output only current MPI figures.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"

if (( $# > 1 )); then
    echo "Usage: bash sandbox/analyze_existing_release120_scaling.sh [output_directory]" >&2
    exit 2
fi
if (( $# == 1 )); then
    out="$(realpath "$1")"
else
    # Select the newest completed report bundle, irrespective of the
    # current Git commit (which changes as we update the plotting code).
    out=""
    shopt -s nullglob
    dirs=(
        "$root"/sandbox/outputs/release120-coo-*
        "$root"/sandbox/outputs/mpi-coo-*
    )
    shopt -u nullglob
    if (( ${#dirs[@]} == 0 )); then
        echo "No release120-coo output directories found" >&2
        exit 2
    fi
    while IFS= read -r candidate; do
        if [[ -s "$candidate/current_never_1r_rep1.json" &&
              -s "$candidate/current_on_change_6r_rep1.json" ]]; then
            out="$candidate"
            break
        fi
    done < <(ls -dt "${dirs[@]}")
fi
if [[ ! -d "$out" || -z "$out" ]]; then
    echo "No complete native benchmarking directory found" >&2
    exit 2
fi
if [[ ! -s "$out/current_on_change_6r_rep1.json" ]]; then
    echo "Missing complete MPI 6-rank COO results: $out" >&2
    exit 1
fi
shopt -s nullglob
files=( "$out"/current_never_1r_rep*.json )
shopt -u nullglob
reps="${#files[@]}"
if (( reps < 1 || reps > 10 )); then
    echo "Invalid 1..10 repetitions in $out" >&2
    exit 2
fi
echo "SCRIMP existing-data analysis ONLY; NO new physical simulation."
echo "Source results: $out"
echo "Number of repetitions: $reps"
sudo docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp -e MPLBACKEND=Agg \
    --cpuset-cpus="0,2" \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$out:/artifacts" \
    scrimp-mpi:local \
    bash -c '
        set -euo pipefail
        if [[ -s /artifacts/release_120_release_native_1r_rep1.json ]]; then
            echo "DIAGNOSE: v1.2.0 vs current raw/sorted full state"
            python -u sandbox/diagnose_release120_states.py \
                --directory /artifacts
        else
            echo "MPI-ONLY: no v1.2.0 data; historical diagnostics skipped"
        fi
        echo "VERIFY+PLOT: current COO/MPI 1/2/4/6 and never/on_change"
        python -u sandbox/plot_release120_coo_scaling.py \
            --directory /artifacts --repetitions "$1" --current-only
    ' _ "$reps"
if [[ -s "$out/png/05_steady_step_time_coo.png" &&
      -s "$out/png/08_first_step_overhead_coo.png" ]]; then
    echo "PASS: 7 current COO/MPI figures (global + steady-step metrics)."
else
    echo "PASS: 3 global COO/MPI figures; no step timing in legacy data."
fi
echo "NOTE: no cross-version gain is asserted in current-only figures."
echo "PNG figures: $out/png/"
echo "SVG figures: $out/svg/"
echo "Reports: $out/performance_summary.json (and optional release_state_diagnostic.json)"
