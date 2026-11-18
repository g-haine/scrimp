#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.
#
# v1.3.0 qualification: serial v1.2.0 vs current + 1/2/4/6 MPI
# with exact on_change symbolic/numeric reuse. Original wave.py,
# fixed Crank-Nicolson and full final-state parity, no production edits.
set -euo pipefail

mesh="${1:-0.04}"
repetitions="${2:-3}"
ranks_csv="${SCRIMP_V130_BENCH_RANKS:-1,2,4,6}"
if [[ $# -gt 2 || ! "$mesh" =~ ^0\.[0-9]+$ ||
      ! "$repetitions" =~ ^[1-5]$ ||
      ! "$ranks_csv" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ]]; then
    echo "Usage: bash sandbox/run_v130_wave_benchmark.sh [mesh=0.04] [reps=3]" >&2
    echo "Optional SCRIMP_V130_BENCH_RANKS=1,2,4,6; SCRIMP_BENCH_TIMEOUT_S=900" >&2
    exit 2
fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
current_sha="$(git rev-parse HEAD)"
git rev-parse -q --verify "refs/tags/v1.2.0^{commit}" >/dev/null || git fetch --tags origin
release_sha="$(git rev-parse "refs/tags/v1.2.0^{commit}")"
out="$(mktemp -d "$root/sandbox/outputs/v130-wave-h${mesh//./p}-${current_sha:0:8}-XXXXXX")"
worktree="$(mktemp -d "$root/sandbox/outputs/.v120-source-XXXXXX")"
rmdir "$worktree"
cleanup() {
    if [[ -e "$worktree/.git" ]]; then
        git worktree remove --force "$worktree" >/dev/null 2>&1 || true
    fi
    [[ ! -d "$worktree" ]] || rm -rf "$worktree"
}
trap cleanup EXIT
git worktree add --detach "$worktree" v1.2.0 >/dev/null
if ! cmp -s "$worktree/examples/wave.py" "$root/examples/wave.py"; then
    echo "FAIL: wave.py differs between release and current; no apples-to-oranges benchmark" >&2
    exit 1
fi
IFS=, read -r -a ranks <<< "$ranks_csv"
IFS=, read -r -a cpu_pool <<< "${SCRIMP_SCALING_CPUS:-0,2,4,6,8,10}"
timeout_s="${SCRIMP_BENCH_TIMEOUT_S:-900}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] || (( timeout_s < 60 || timeout_s > 7200 )); then
    echo "SCRIMP_BENCH_TIMEOUT_S must be between 60 and 7200 seconds" >&2; exit 2
fi
echo "SCRIMP v1.3.0 qualification: original v1.2.0 serial vs current serial/MPI"
echo "Release=$release_sha current=$current_sha mesh=$mesh reps=$repetitions"
echo "Time: CN dt=.002, final=.02 (10 steps); current MUMPS on_change default+MPI ordering"
echo "Output: $out"
printf '{"release_sha":"%s","current_sha":"%s","mesh":%s,"repetitions":%s,"ranks":"%s"}\n' \
    "$release_sha" "$current_sha" "$mesh" "$repetitions" "$ranks_csv" > "$out/manifest.json"
sudo -v

run_case() {
    local src="$1" variant="$2" policy="$3" p="$4" rep="$5" revision="$6"
    if (( p > ${#cpu_pool[@]} )); then echo "Insufficient CPUs for $p ranks" >&2; return 2; fi
    local cpus label
    cpus="$(IFS=,; echo "${cpu_pool[*]:0:p}")"
    label="${variant}_${policy}_${p}r_rep${rep}"
    echo "RUN $label cpus=$cpus"
    local launcher=(python -u)
    local config=()
    local switches=()
    if [[ "$variant" == "current" ]]; then
        launcher=(mpiexec --bind-to none -n "$p" python -u)
        config=(SCRIMP_GETFEM_MATRIX_MODE=distributed
                SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo
                SCRIMP_MATRIX_MEMORY_MODE=compact
                SCRIMP_MUMPS_RESET_POLICY=on_change)
        switches=(--profile --petsc-log)
    fi
    if ! sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        --cpuset-cpus="$cpus" --workdir /workspace/scrimp \
        -e HOME=/tmp -e MPLBACKEND=Agg \
        -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
        -v "$src/scrimp:/workspace/scrimp/scrimp:ro" \
        -v "$src/examples:/workspace/scrimp/examples:ro" \
        -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
        -v "$out:/workspace/scrimp/outputs" \
        -v "$out:/artifacts" \
        scrimp-mpi:local timeout -k 10s "${timeout_s}s" \
        env -u SCRIMP_MUMPS_RESET_POLICY -u SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN \
            -u SCRIMP_GETFEM_MATRIX_MODE -u SCRIMP_GETFEM_BRIDGE_DISTRIBUTED \
            -u SCRIMP_MATRIX_MEMORY_MODE \
            "${config[@]}" "${launcher[@]}" \
            sandbox/bench_wave_release_scaling.py \
            --variant "$variant" --policy "$policy" --mesh "$mesh" \
            --dt 0.002 --final-time 0.02 --save-dt 0.01 \
            --repetition "$rep" --source-sha "$revision" \
            --source-tag "$( [[ "$variant" == release_120 ]] && echo v1.2.0 || : )" \
            --out "/artifacts/$label.json" "${switches[@]}" \
        > "$out/$label.log" 2>&1; then
        echo "NON-PASS $label; stopping campaign" >&2
        tail -n 70 "$out/$label.log" >&2
        return 1
    fi
    grep -aE '^PASS (release_120|current)/' "$out/$label.log" | tail -1 || true
    if [[ ! -s "$out/$label.json" || ! -s "$out/${label}_state.npz" ]]; then
        echo "FAIL: report or full-state output missing for $label" >&2
        return 1
    fi
}

for rep in $(seq 1 "$repetitions"); do
    # Native release: only one sequential process, unmodified source.
    run_case "$worktree" release_120 release_native 1 "$rep" "$release_sha"
    for p in "${ranks[@]}"; do
        run_case "$root" current on_change "$p" "$rep" "$current_sha"
    done
done
sudo docker run --rm \
    --user "$(id -u):$(id -g)" --workdir /workspace/scrimp \
    -e HOME=/tmp -e MPLBACKEND=Agg \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$out:/artifacts" \
    scrimp-mpi:local python -u sandbox/analyze_v130_wave_benchmark.py \
        --directory /artifacts --repetitions "$repetitions" --ranks "$ranks_csv"
echo "PASS: v1.3.0 wave benchmark, full-state checks and stage report"
echo "Reports: $out"
