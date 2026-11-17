#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Reproducible comparison:
#   historical tag v1.2.0, 1 rank, original code (no retrofitted COO);
#   current source, distributed PETSc COO + MUMPS never/on_change,
#   1,2,4,6 ranks, exactly the same 2D acoustic example.
# No long physical simulations: 10 fixed Crank-Nicolson steps per run.
# Examples:
#   bash sandbox/run_release120_coo_scaling.sh 0.08 3
#   bash sandbox/run_release120_coo_scaling.sh 0.04 3 --mpi-only
#   bash sandbox/run_release120_coo_scaling.sh --mpi-only 0.025 2
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: bash sandbox/run_release120_coo_scaling.sh [--mpi-only] [mesh_h=0.08] [repetitions=2]

  --mpi-only  Only CURRENT SCRIMP, COO + MUMPS never/on_change, at
              1/2/4/6 MPI ranks. No v1.2.0 tag, worktree, run or
              cross-version speedup claim. Seven PNG + SVG figures
              (total solve + measured cost per step, startup excluded).
  -h, --help  Show this message.

Default (without --mpi-only): additionally benchmark original v1.2.0
on one process and require cross-version full-state parity before
producing eight figures, including per-step strong scaling.
USAGE
}

mpi_only=0
positional=()
for arg in "$@"; do
    case "$arg" in
        --mpi-only) mpi_only=1 ;;
        -h|--help) usage; exit 0 ;;
        --*) echo "Unknown option: $arg" >&2; usage >&2; exit 2 ;;
        *) positional+=("$arg") ;;
    esac
done
if (( ${#positional[@]} > 2 )); then
    usage >&2
    exit 2
fi
mesh="${positional[0]:-0.08}"
repetitions="${positional[1]:-2}"
if [[ ! "$mesh" =~ ^0\.[0-9]+$ ]] ||
   [[ ! "$repetitions" =~ ^[1-9][0-9]*$ ]] ||
   (( repetitions > 10 )); then
    usage >&2
    exit 2
fi

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
current_sha="$(git rev-parse HEAD)"
mkdir -p "$root/sandbox/outputs"

baseline_src=""
release_sha=""
if (( mpi_only )); then
    out="$(mktemp -d "$root/sandbox/outputs/mpi-coo-h${mesh//./p}-${current_sha:0:8}-XXXXXX")"
else
    # The historical release is ONLY checked out for cross-version
    # comparisons. --mpi-only never looks up the tag or creates a worktree.
    git rev-parse -q --verify "refs/tags/v1.2.0^{commit}" >/dev/null ||
        git fetch --tags origin
    release_sha="$(git rev-parse "refs/tags/v1.2.0^{commit}")"
    if [[ "$release_sha" == "$current_sha" ]]; then
        echo "Current HEAD is the original tag; cannot compare identical sources" >&2
        exit 2
    fi
    out="$(mktemp -d "$root/sandbox/outputs/release120-coo-h${mesh//./p}-${current_sha:0:8}-XXXXXX")"
    baseline_src="$(mktemp -d "$root/sandbox/outputs/.bench-v120-source-XXXXXX")"
    # git worktree add requires a non-existent destination.
    rmdir "$baseline_src"
    cleanup() {
        if [[ -f "$baseline_src/.git" || -d "$baseline_src/.git" ]]; then
            git worktree remove --force "$baseline_src" >/dev/null 2>&1 || true
        fi
        [[ ! -d "$baseline_src" ]] || rm -rf "$baseline_src"
    }
    trap cleanup EXIT
    git worktree add --detach "$baseline_src" v1.2.0 >/dev/null
    if ! cmp -s "$baseline_src/examples/wave.py" "$root/examples/wave.py"; then
        echo "STOP: wave.py differs between v1.2.0 and current source!" >&2
        exit 1
    fi
fi
if ! sudo -v; then
    echo "sudo authentication failed" >&2
    exit 1
fi

if (( mpi_only )); then
    printf 'SCRIMP CURRENT COO+MUMPS MPI ONLY | h=%s | reps=%s\n' \
        "$mesh" "$repetitions"
    printf 'Current commit: %s | policies: never, on_change | ranks: 1,2,4,6\n' \
        "${current_sha:0:10}"
    echo 'No v1.2.0 worktree or release run.'
else
    printf 'SCRIMP original v1.2.0 vs COO MPI scaling | h=%s | reps=%s\n' \
        "$mesh" "$repetitions"
    printf 'Baseline tag %s | Current commit %s\n' \
        "${release_sha:0:10}" "${current_sha:0:10}"
    echo 'Example file: identical at BOTH revisions.'
fi
printf 'Output: %s\n' "$out"
pcores=(0 2 4 6 8 10)
ranks=(1 2 4 6)
timeout_s="${SCRIMP_BENCH_TIMEOUT_S:-240}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] ||
   (( timeout_s < 30 || timeout_s > 3600 )); then
    echo "SCRIMP_BENCH_TIMEOUT_S must be 30..3600 seconds" >&2
    exit 2
fi

valid_json() {
    /usr/bin/python3 - "$1" "$2" "$3" "$4" "$5" <<'PY'
import json, math, sys
path, variant, policy, ranks, commit = (
    sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]),
    sys.argv[5],
)
try:
    r = json.load(open(path, encoding="utf-8"))
    assert r["variant"] == variant
    assert r["mumps_policy"] == policy
    assert r["mpi_ranks"] == ranks
    assert r["source_sha"] == commit
    assert r["steps"] == 10
    assert r["rejections"] == 0 and r["snes_failures"] == 0
    assert abs(r["time_final_monitor"] - 0.02) <= 1e-10
    assert math.isfinite(r["wall_solve_s"]) and r["wall_solve_s"] > 0
    assert r["step_warmup_excluded"] == 1
    assert r["steady_step_count"] == 9
    assert len(r["step_times_max_rank_s"]) == 10
    assert len(r["step_times_by_rank_s"]) == ranks
    assert all(len(times) == 10 for times in r["step_times_by_rank_s"])
    assert math.isfinite(r["first_step_wall_s"]) and r["first_step_wall_s"] > 0
    assert math.isfinite(r["steady_step_mean_s"]) and r["steady_step_mean_s"] > 0
except (AssertionError, OSError, ValueError, KeyError, TypeError):
    sys.exit(1)
PY
}

run_one() {
    local src="$1" variant="$2" policy="$3" n="$4" rep="$5" sha="$6"
    local cpus label report log
    cpus="$(IFS=,; echo "${pcores[*]:0:n}")"
    label="${variant}_${policy}_${n}r_rep${rep}"
    report="$out/$label.json"
    log="$out/$label.log"
    printf 'RUN %s source=%s cpus=%s\n' "$label" "${sha:0:8}" "$cpus"

    local extras=()
    # The historical control is a genuinely SEQUENTIAL Python process,
    # not a one-rank MPI launcher. PETSc itself may be MPI-linked.
    local launcher=(python -u)
    if [[ "$variant" == current ]]; then
        launcher=(mpiexec --bind-to none -n "$n" python -u)
        extras=(
            SCRIMP_GETFEM_MATRIX_MODE=distributed
            SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo
            SCRIMP_MUMPS_RESET_POLICY="$policy"
            SCRIMP_MATRIX_MEMORY_MODE=compact
        )
    fi
    if ! sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        -e HOME=/tmp -e MPLBACKEND=Agg \
        -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
        --cpuset-cpus="$cpus" \
        -v "$src/scrimp:/workspace/scrimp/scrimp:ro" \
        -v "$src/examples:/workspace/scrimp/examples:ro" \
        -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
        -v "$out:/workspace/scrimp/outputs" \
        -v "$out:/artifacts" \
        scrimp-mpi:local \
        timeout -k 10s "${timeout_s}s" \
        env -u SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN \
            -u SCRIMP_MUMPS_RESET_POLICY \
            -u SCRIMP_GETFEM_BRIDGE_DISTRIBUTED \
            -u SCRIMP_MATRIX_MEMORY_MODE \
            "${extras[@]}" \
        "${launcher[@]}" sandbox/bench_wave_release_scaling.py \
            --variant "$variant" --policy "$policy" \
            --mesh "$mesh" --dt 0.002 --final-time 0.02 \
            --save-dt 0.01 --source-sha "$sha" \
            --source-tag "$( [[ "$variant" == release_120 ]] && echo v1.2.0 )" \
            --repetition "$rep" \
            --out "/artifacts/$label.json" \
        >"$log" 2>&1; then
        echo "FAILED: $label; first/last 60 lines follow" >&2
        head -n 60 "$log" >&2
        tail -n 60 "$log" >&2
        return 1
    fi
    if ! valid_json "$report" "$variant" "$policy" "$n" "$sha"; then
        echo "INVALID benchmark output: $label" >&2
        tail -n 60 "$log" >&2
        return 1
    fi
    /usr/bin/python3 - "$report" <<'PY'
import json,sys
a=json.load(open(sys.argv[1],encoding="utf-8"))
print("  PASS: %s DoFs, 10 steps, solve %.3fs, first %.4fs, steady %.4fs, HWM %.1f MiB"
      %(a["dofs"],a["wall_solve_s"],a["first_step_wall_s"],
        a["steady_step_mean_s"],a["max_rank_hwm_kib"]/1024.))
PY
}

for rep in $(seq 1 "$repetitions"); do
    # The RELEASE control is run in its original source worktree,
    # always sequentially; never mislabel COO/on_change as v1.2.0.
    if (( !mpi_only )); then
        run_one "$baseline_src" release_120 release_native 1 "$rep" "$release_sha"
    fi
    for n in "${ranks[@]}"; do
        # Rotate policy order to reduce cache/thermal bias.
        if (( rep % 2 == 0 )); then
            order=(on_change never)
        else
            order=(never on_change)
        fi
        for policy in "${order[@]}"; do
            run_one "$root" current "$policy" "$n" "$rep" "$current_sha"
        done
    done
done

echo "VALIDATE + PLOT: all initial/final coefficients and all rank cases"
plot_flags=()
if (( mpi_only )); then
    plot_flags+=(--current-only)
fi
sudo docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp -e MPLBACKEND=Agg \
    --cpuset-cpus="0,2" \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$out:/artifacts" \
    scrimp-mpi:local \
    python -u sandbox/plot_release120_coo_scaling.py \
        --directory /artifacts --repetitions "$repetitions" \
        "${plot_flags[@]}"

if (( mpi_only )); then
    echo "SUCCESS. Seven COO/MPI-only PNG + SVG figures, CSV and full data:"
else
    echo "SUCCESS. Eight PNG + SVG figures, CSV and full data:"
fi
ls -lh "$out"/png/*.png "$out"/svg/*.svg "$out"/*.csv "$out"/*.json | tail -n 26
echo "PNG figures: $out/png/"
echo "SVG figures: $out/svg/"
echo "Data and reports: $out/"
