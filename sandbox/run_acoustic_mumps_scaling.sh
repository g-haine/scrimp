#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Same fixed-mesh acoustic comparison as bench_acoustic_2d.py:
#   bash sandbox/run_acoustic_mumps_scaling.sh 0.04 pilot
#   bash sandbox/run_acoustic_mumps_scaling.sh 0.04 full
#
# Run with "bash", NOT "source". The script cannot close the caller's shell.
# Pilot: rank 1 and rank 6, once each, both policies (4 runs).
# Full:  ranks 1/2/4/6, 5 reps, both policies (40 runs); reuse pilot.
# Every mesh uses a separate directory and the existing analysis/plots.
# The CPU list uses one thread per physical core on port-haine.

set -uo pipefail

mesh="${1:-0.04}"
mode="${2:-pilot}"
if [[ ! "$mesh" =~ ^0\.[0-9]+$ ]] ||
   [[ "$mode" != pilot && "$mode" != full ]]; then
    echo "Usage: bash sandbox/run_acoustic_mumps_scaling.sh 0.04 [pilot|full]" >&2
    exit 2
fi

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
out="$root/sandbox/outputs/scaling-mumps-h${mesh//./p}"
mkdir -p "$out" || exit 1

pcores=(0 2 4 6 8 10)
if [[ "$mode" == pilot ]]; then
    ranks=(1 6)
    repetitions=1
else
    ranks=(1 2 4 6)
    repetitions=5
fi

# Strictly validate the JSON output, not just its existence.
valid_report() {
    python - "$1" "$2" "$3" "$mesh" <<'PY'
import json
import math
import sys
from pathlib import Path

path, ranks, policy, h = sys.argv[1], int(sys.argv[2]), sys.argv[3], float(sys.argv[4])
try:
    run = json.loads(Path(path).read_text(encoding="utf-8"))
    expected = "never" if policy == "auto" else "always"
    assert run["series"] == policy
    assert run["mpi_size"] == ranks
    assert run["effective_mumps_reset_policy"] == expected
    assert abs(float(run["rectangle_mesh_h"]) - h) <= 1e-12
    assert run["degrees_of_freedom"] > 12682
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
    assert run["final_variable_norms"]
    assert all(math.isfinite(float(v)) for v in run["final_variable_norms"].values())
except (OSError, ValueError, TypeError, KeyError, AssertionError):
    sys.exit(1)
PY
}

# Fetch sudo credentials before long-running measurements.
if ! sudo -v; then
    echo "sudo authentication failed" >&2
    exit 1
fi

failed=0
for rep in $(seq 1 "$repetitions"); do
    for n in "${ranks[@]}"; do
        cpus=$(IFS=,; echo "${pcores[*]:0:n}")
        for policy in auto always; do
            label="wave_${policy}_${n}ranks_${rep}"
            report="$out/$label.json"
            logfile="$out/$label.log"
            if valid_report "$report" "$n" "$policy"; then
                printf 'SKIP valid: %s\n' "$label"
                continue
            fi

            policy_env=()
            if [[ "$policy" == always ]]; then
                policy_env=(SCRIMP_MUMPS_RESET_POLICY=always)
            fi
            printf 'RUN h=%s MPI=%s policy=%s repeat=%s cpus=%s\n' \
                "$mesh" "$n" "$policy" "$rep" "$cpus"

            if ! sudo docker run --rm \
                --user "$(id -u):$(id -g)" \
                -e HOME=/tmp \
                --cpuset-cpus="$cpus" \
                -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
                -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
                -v "$root/examples:/workspace/scrimp/examples:ro" \
                -v "$out:/workspace/scrimp/outputs" \
                -v "$out:/artifacts" \
                scrimp-mpi:local \
                timeout -k 10s 240s \
                env -u SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN \
                    -u SCRIMP_MUMPS_RESET_POLICY "${policy_env[@]}" \
                mpiexec --bind-to none -n "$n" \
                python -u sandbox/bench_acoustic_2d.py \
                    --mesh "$mesh" --dt 0.002 --final-time 0.04 \
                    --save-dt 0.01 --series "$policy" \
                    --output "/artifacts/$label.json" \
                > "$logfile" 2>&1; then
                printf 'FAILED execution: %s\n' "$label" >&2
                tail -n 45 "$logfile"
                failed=1
                break 3
            fi
            if ! valid_report "$report" "$n" "$policy"; then
                printf 'FAILED report validation: %s\n' "$label" >&2
                tail -n 45 "$logfile"
                failed=1
                break 3
            fi
            python - "$report" <<'PY'
import json, sys
r = json.load(open(sys.argv[1], encoding="utf-8"))
print("  PASS: DOFs={degrees_of_freedom}, solve={wall_solve_seconds_max:.3f}s, policy={effective_mumps_reset_policy}, rank-RSS={peak_rank_rss_kib} KiB".format(**r))
PY
        done
    done
done

if [[ "$failed" -ne 0 ]]; then
    echo "Stopped safely. Existing validated JSON reports remain in $out" >&2
    exit 1
fi

if [[ "$mode" == pilot ]]; then
    echo
    echo "Pilot finished. Inspect times and memory before the 40-run series:"
    echo "    bash sandbox/run_acoustic_mumps_scaling.sh $mesh full"
else
    echo
    echo "Validating and plotting the complete fixed-mesh comparison..."
    python sandbox/analyze_acoustic_scaling.py \
        "$out"/wave_*ranks_*.json \
        --plot-comparison "$out/comparaison_mumps"
    echo "Graphs:"
    ls -lh "$out"/comparaison_mumps_*.png
fi
