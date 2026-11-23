#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.
#
# Qualification is PARITY (uninterrupted == paused/resumed), not
# mathematical energy conservation of nonlinear-mass TSCN.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
[[ "$#" == 0 ]] || {
    echo "Usage: $0 (full gamma=0.25, h=.2/.1, MPI source/target=1/2 campaign)" >&2
    exit 2
}
python3 -m py_compile sandbox/probe_v130_rotating_controlled_ic.py \
    sandbox/compare_v130_rotating_ic_parity.py
python3 tests/test_v130_rotating_controlled_ic.py
python3 tests/test_v130_rotating_ic_parity.py
mkdir -p sandbox/outputs
campaign="$(mktemp -d "$root/sandbox/outputs/v130-rotating-parity-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
[[ "$campaign" == "$root/sandbox/outputs/"* ]] || exit 2
printf 'case\thost_evidence\n' > "$campaign/cases.tsv"
printf '{"schema":"scrimp-v130-rotating-checkpoint-native-parity-campaign-v1","git_source":"%s","image":"%s","experiment":"controlled_rotating_swe_gamma_0.25","physical_model_changed":false,"energy_magnitude_not_a_gate":true,"meshes":[0.2,0.1],"mpi":[1,2],"source_time":0.02,"final_time":0.024}\n' \
    "$(git rev-parse HEAD)" "${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}" > "$campaign/manifest.json"

echo "SCRIMP v1.3 rotating SWE checkpoint native parity campaign"
echo "h=.2/.1; 0->.02 pause; 0->.024 reference, .02->.024 resumed; MPI 1/2"
echo "Only difference is stopping/restarting and MPI rank count."
echo "Original physical relative Hamiltonian defect RECORDED, not a gate."
echo "All evidence retained under $campaign; each native case stays immutable."

run_case() {
    local name="$1" phase="$2" mesh="$3" ranks="$4" tf="$5"
    local checkpoint="${6:-}"
    local logfile="$campaign/$name.host.log"
    local argv=( "$phase" "$mesh" "$ranks" "$tf" )
    if [[ "$phase" == resume ]]; then
        [[ -s "$checkpoint" ]] || {
            echo "Checkpoint missing: $checkpoint" >&2; return 1;
        }
        argv+=( "$checkpoint" )
    fi
    echo "RUN $name [physical model unchanged; parity mode]" >&2
    if ! SCRIMP_V130_ROTATING_IC_QUALIFICATION=restart_parity \
        bash scripts/run-v130-rotating-controlled-ic.sh "${argv[@]}" > "$logfile" 2>&1; then
        echo "FAIL native case $name; evidence retained: $campaign" >&2
        tail -n 120 "$logfile" >&2
        return 1
    fi
    local path
    path="$(sed -n 's/^Evidence on host: //p' "$logfile" | tail -n 1)"
    [[ -n "$path" && "$path" == "$root/sandbox/outputs/"* \
        && -s "$path/accepted_state.npz" \
        && -s "$path/initial_state_report.json" \
        && -s "$path/initial_state_trace.npz" ]] || {
        echo "FAIL: incomplete native original checkpoint evidence for $name" >&2
        return 1
    }
    printf '%s\t%s\n' "$name" "$path" >> "$campaign/cases.tsv"
    grep -aE '^ROTATING CONTROLLED IC ACCEPTED|^HOST EVIDENCE PASS' "$logfile" >&2
    printf '%s\n' "$path"
}
declare -A seeds
declare -A continuations

for mesh in 0.2 0.1; do
    label="h${mesh/./p}"
    for ranks in 1 2; do
        key="${label}_${ranks}"
        seeds["$key"]="$(run_case "${label}_seed_${ranks}r" seed "$mesh" "$ranks" 0.02)"
        continuations["$key"]="$(run_case "${label}_continuous_${ranks}r" seed "$mesh" "$ranks" 0.024)"
    done
    for source in 1 2; do
        for target in 1 2; do
            name="${label}_restart_${source}to${target}"
            source_dir="${seeds[${label}_${source}]}"
            reference_dir="${continuations[${label}_${target}]}"
            resumed_dir="$(run_case "$name" resume "$mesh" "$target" 0.024 \
                "$source_dir/accepted_state.npz")"
            python3 sandbox/compare_v130_rotating_ic_parity.py \
                --seed "$source_dir" \
                --continuous "$reference_dir" \
                --resume "$resumed_dir" \
                --report "$campaign/$name.parity.json" | tee "$campaign/$name.comparison.log"
        done
    done
done

python3 - "$campaign" <<'PY'
import json, sys
from pathlib import Path
out=Path(sys.argv[1])
results=sorted(out.glob("*.parity.json"))
assert len(results)==8, f"Expected 8 native physical comparisons, got {len(results)}"
rows=[json.loads(path.read_text(encoding="utf-8")) for path in results]
assert all(r["status"]=="PASS_RESTART_EQUIVALENT_TO_UNINTERRUPTED" for r in rows)
actual={(r["mesh"],r["source_ranks"],r["resume_ranks"]) for r in rows}
expected={(m,s,t) for m in (.2,.1) for s in (1,2) for t in (1,2)}
assert actual==expected, f"Missing original physical MPI combinations: {expected-actual}"
summary={
    "status":"PASS_COMPLETE_NATIVE_CHECKPOINT_PARITY_CAMPAIGN",
    "schema":"scrimp-v130-controlled-rotating-swe-checkpoint-parity-v1",
    "independent_native_comparisons":len(rows),
    "meshes":[.2,.1], "MPI_source_target":[[1,1],[1,2],[2,1],[2,2]],
    "max_full_trajectory_FE_L2":max(
        r["maximum_full_FE_trajectory_L2_difference"] for r in rows
    ),
    "max_H_plus_D_minus_W_difference":max(
        r["H_plus_D_minus_W_maximum_absolute_difference"] for r in rows
    ),
    "power_balance_integrator_accuracy_claimed":False,
    "original_physics_modified":False,
    "checkpoint_TS_history_restored":False,
}
(out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
print(json.dumps(summary,indent=2))
PY
echo "COMPLETE: full native SCRIMP checkpoint equivalence evidence in $campaign"
