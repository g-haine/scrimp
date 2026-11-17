#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Strong scaling, fixed 2D acoustic mesh, historical CSR vs PETSc COO,
# crossed independently with explicit MUMPS reset policies.
#
#   bash sandbox/run_acoustic_mumps_scaling.sh 0.02 sweep both never
#   bash sandbox/run_acoustic_mumps_scaling.sh 0.04 pilot both never
#   bash sandbox/run_acoustic_mumps_scaling.sh 0.04 full both never,always
#   bash sandbox/run_acoustic_mumps_scaling.sh 0.04 full both never,on_change,always
#
# pilot = ranks 1,6 (1 repetition per variant);
# sweep = ranks 1,2,4,6 (1 repetition per variant);
# full  = ranks 1,2,4,6 (5 repetitions per variant; reuses valid pilot/sweep).
# Default policy set: never,always; on_change is opt-in for linear acoustics.
# All SIMULATIONS and PLOTTING run inside scrimp-mpi:local (NOT host Conda).
# Old scaling-mumps-h* results are left untouched; new results are
# versioned by source commit and have explicit (bridge,policy) labels.
# Run with bash, never source.

set -euo pipefail

mesh="${1:-0.04}"
mode="${2:-pilot}"
bridge_arg="${3:-both}"
policies_arg="${4:-never,always}"

if [[ ! "$mesh" =~ ^0\.[0-9]+$ ]] ||
   [[ "$mode" != pilot && "$mode" != sweep && "$mode" != full ]] ||
   [[ "$bridge_arg" != csr && "$bridge_arg" != coo && "$bridge_arg" != both ]]; then
    echo "Usage: bash sandbox/run_acoustic_mumps_scaling.sh 0.04 [pilot|sweep|full] [csr|coo|both] [never,always|never,on_change,always|...]" >&2
    exit 2
fi
case "$bridge_arg" in
    csr) bridges=(csr) ;;
    coo) bridges=(coo) ;;
    both) bridges=(csr coo) ;;
esac

IFS=, read -r -a policies <<< "$policies_arg"
[[ "${#policies[@]}" -ge 1 ]] || {
    echo "At least one MUMPS policy is required" >&2
    exit 2
}
declare -A seen_policy=()
for policy in "${policies[@]}"; do
    case "$policy" in
        never|on_change|always|auto) ;;
        *)
            echo "Unknown policy: $policy (never/on_change/always/auto)" >&2
            exit 2 ;;
    esac
    normalized="$policy"
    [[ "$policy" == auto ]] && normalized=never
    if [[ -n "${seen_policy[$normalized]:-}" ]]; then
        echo "Duplicate effective MUMPS policy: $policy" >&2
        exit 2
    fi
    seen_policy["$normalized"]=1
done

# Large h=0.02 pilots may require more than the previous 240 s/run.
# Explicit, bounded per-run timeout, independent of Docker/Conda state.
timeout_s="${SCRIMP_SCALING_TIMEOUT_S:-240}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] ||
   (( timeout_s < 30 || timeout_s > 3600 )); then
    echo "SCRIMP_SCALING_TIMEOUT_S must be 30..3600 seconds" >&2
    exit 2
fi

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
out="$root/sandbox/outputs/scaling-bridge-h${mesh//./p}-${sha}"
mkdir -p "$out"

pcores=(0 2 4 6 8 10)
if [[ "$mode" == pilot ]]; then
    ranks=(1 6)
    repetitions=1
elif [[ "$mode" == sweep ]]; then
    # Run the finer-mesh pilot from most to least parallel first.
    # The 1-rank historical CSR case can be the slowest / memory-heavy.
    # This does not affect the rank labels or the scaling analysis.
    ranks=(6 4 2 1)
    repetitions=1
else
    ranks=(1 2 4 6)
    repetitions=5
fi

# A report can only be reused for this exact benchmark configuration,
# source commit, method, policy, rank count and 20-step CN endpoint.
valid_report() {
    /usr/bin/python3 - "$1" "$2" "$3" "$4" "$mesh" "$sha" <<'PY'
import json
import math
import sys
from pathlib import Path

path, ranks, policy, bridge, h, sha = (
    sys.argv[1], int(sys.argv[2]), sys.argv[3],
    sys.argv[4], float(sys.argv[5]), sys.argv[6]
)
try:
    run = json.loads(Path(path).read_text(encoding="utf-8"))
    expected = "never" if policy == "auto" else policy
    assert run["series"] == policy
    assert run["bridge_method"] == bridge
    assert run["bridge_backend"] == bridge
    assert run["getfem_matrix_mode"] == "distributed"
    assert str(run["benchmark_source_commit"]) == sha
    assert run["mpi_size"] == ranks
    assert run["effective_mumps_reset_policy"] == expected
    assert math.isclose(float(run["rectangle_mesh_h"]), h, abs_tol=1e-12)
    assert run["degrees_of_freedom"] > 0
    assert run["ts_type"] == "cn"
    assert run["ts_adapt_type"] == "none"
    assert run["ts_converged_reason"] > 0
    assert run["ts_steps"] == 20
    assert run["ts_rejected_steps"] == 0
    assert run["ts_snes_failures"] == 0
    assert abs(float(run["requested_final_time"]) - 0.04) < 1e-12
    assert abs(float(run["actual_final_time"]) - 0.04) < 1e-12
    assert abs(float(run["last_saved_time"]) - 0.04) < 1e-12
    assert math.isfinite(float(run["wall_solve_seconds_max"]))
    assert run["wall_solve_seconds_max"] > 0
    norms = run["final_variable_norms"]
    assert norms and all(math.isfinite(float(x)) for x in norms.values())
except (OSError, ValueError, TypeError, KeyError, AssertionError):
    sys.exit(1)
PY
}

if ! sudo -v; then
    echo "sudo authentication failed" >&2
    exit 1
fi
echo "SCRIMP $sha | h=$mesh | policies=$policies_arg | bridge=$bridge_arg | $mode | timeout=${timeout_s}s/run"
echo "Output: $out"

for rep in $(seq 1 "$repetitions"); do
    for n in "${ranks[@]}"; do
        cpus="$(IFS=,; echo "${pcores[*]:0:n}")"
        # Alternate IJV/COO ordering with repetitions to limit
        # systematic thermal/cache/order bias, without random seeds.
        if (( rep % 2 == 0 )); then
            ordered_bridges=()
            for (( k=${#bridges[@]}-1; k>=0; k-- )); do
                ordered_bridges+=("${bridges[k]}")
            done
        else
            ordered_bridges=("${bridges[@]}")
        fi
        for policy in "${policies[@]}"; do
            for bridge in "${ordered_bridges[@]}"; do
                label="wave_${bridge}_${policy}_${n}ranks_${rep}"
                report="$out/$label.json"
                logfile="$out/$label.log"
                if valid_report "$report" "$n" "$policy" "$bridge"; then
                    printf 'SKIP valid: %s\n' "$label"
                    continue
                fi
                # Explicit MUMPS policy, not a falsely accepted 'auto'
                # environment value. Historical --series auto still
                # means no override (automatic never on linear bricks).
                policy_env=()
                if [[ "$policy" != auto ]]; then
                    policy_env=(SCRIMP_MUMPS_RESET_POLICY="$policy")
                fi
                printf 'RUN h=%s MPI=%s bridge=%s MUMPS=%s rep=%s cpus=%s\n' \
                    "$mesh" "$n" "$bridge" "$policy" "$rep" "$cpus"

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
                        -u SCRIMP_MUMPS_RESET_POLICY \
                        SCRIMP_GETFEM_MATRIX_MODE=distributed \
                        SCRIMP_GETFEM_BRIDGE_DISTRIBUTED="$bridge" \
                        SCRIMP_BENCH_GIT_SHA="$sha" \
                        "${policy_env[@]}" \
                    mpiexec --bind-to none -n "$n" \
                    python -u sandbox/bench_acoustic_2d.py \
                        --mesh "$mesh" --dt 0.002 --final-time 0.04 \
                        --save-dt 0.01 --series "$policy" --profile \
                        --output "/artifacts/$label.json" \
                    >"$logfile" 2>&1; then
                    echo "FAILED: $label" >&2
                    tail -n 45 "$logfile" >&2
                    exit 1
                fi

                if ! valid_report "$report" "$n" "$policy" "$bridge"; then
                    echo "INVALID report: $label" >&2
                    tail -n 45 "$logfile" >&2
                    exit 1
                fi
                /usr/bin/python3 - "$report" <<'PY'
import json
import sys
r = json.load(open(sys.argv[1], encoding="utf-8"))
print("  PASS: bridge={bridge_method}, MUMPS={effective_mumps_reset_policy}, DOFs={degrees_of_freedom}, solve={wall_solve_seconds_max:.3f}s, peak-rank-RSS={peak_rank_rss_kib} KiB".format(**r))
PY
            done
        done
    done
done

# Graphs are meaningful only if BOTH methods are sampled.
# The analysis validates the complete fixed mesh/policy/method/rank grid.
if [[ "$bridge_arg" == both ]]; then
    echo "Analyzing and plotting the CSR/COO x MUMPS comparison in Docker..."
    shopt -s nullglob
    all_reports=("$out"/wave_*ranks_*.json)
    container_reports=()
    for report in "${all_reports[@]}"; do
        container_reports+=("/artifacts/${report##*/}")
    done
    if ! sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        -e HOME=/tmp \
        -e MPLBACKEND=Agg \
        --cpuset-cpus="0,2" \
        -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
        -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
        -v "$out:/workspace/scrimp/outputs" \
        -v "$out:/artifacts" \
        scrimp-mpi:local \
        python -u sandbox/analyze_acoustic_scaling.py \
            "${container_reports[@]}" \
            --plot-bridge-comparison /artifacts/comparaison_bridge \
            --bridge-policies "$policies_arg"; then
        echo "FAILED chart/report validation; raw JSON left intact" >&2
        exit 1
    fi
    echo "Graphs (speedup, solve time and COO gain):"
    ls -lh "$out"/comparaison_bridge_*.png
fi

if [[ "$mode" == pilot ]]; then
    echo "Pilot done. For 5 repetitions across 1/2/4/6 MPI ranks:"
    echo "  bash sandbox/run_acoustic_mumps_scaling.sh $mesh full both $policies_arg"
fi
echo "Completed. Reports in $out"
