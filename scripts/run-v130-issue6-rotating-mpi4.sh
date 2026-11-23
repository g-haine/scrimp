#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.
#
# Issue #6 native gamma=.25 compatible rotating SWE MPI1/2/4 powers.
# No model/weak form/time scheme/production checkpoint changes.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
[[ "$#" == 0 ]] || {
    echo "Usage: $0 (two meshes; 14 native cases; 2 cross-MPI4 restarts/mesh)" >&2
    exit 2
}
python3 -m py_compile sandbox/probe_v130_rotating_controlled_ic.py \
    sandbox/compare_v130_issue6_active_mpi.py \
    sandbox/compare_v130_rotating_ic_parity.py \
    tests/test_v130_issue6_active_mpi.py
python3 tests/test_v130_issue6_active_mpi.py
python3 tests/test_v130_rotating_controlled_ic.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-issue6-active-mpi4-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf 'name\tpath\n' > "$out/native_cases.tsv"
printf '{"schema":"scrimp-v130-issue6-compatible-rotating-mpi4-campaign-v1","source":"%s","model":"original controlled gamma=0.25 curved-P2 rotating SWE","meshes":[0.2,0.1],"mpi_ranks":[1,2,4],"initial_physical_time":0,"pause_time":0.02,"final_physical_time":0.024,"dt":0.002,"normal_feedback_not_active_before_75":true,"generic_port_power_not_physical_h_weighted_work":true,"scientific_status":"NATIVE_PENDING"}\n' \
    "$(git rev-parse HEAD)" > "$out/manifest.json"
echo "SCRIMP #6 original controlled rotating SWE MPI 1/2/4 qualified-host campaign"
echo "Controlled gamma=.25, P2 curved disk, actual wall spinup (not stationary rest)."
echo "Original normal feedback remains inactive before t=75."
echo "Full FE states, native h-weighted power and original generic MPI port powers."
echo "Per-rank JSON and root-only PNG outputs; original production public checkpoints."
echo "No model energy-accuracy gate: ONLY strict cross-rank and checkpoint equivalence."
echo "EVIDENCE: $out"

run_native() {
    local name="$1" phase="$2" mesh="$3" ranks="$4" final="$5"
    local checkpoint="${6:-}"
    local args=( "$phase" "$mesh" "$ranks" "$final" )
    if [[ "$phase" == resume ]]; then
        [[ -s "$checkpoint" ]] || {
            echo "FAIL missing source checkpoint for $name: $checkpoint" >&2
            return 1
        }
        args+=( "$checkpoint" )
    fi
    echo "RUN $name original physical SWE MPI=$ranks" >&2
    local log="$out/$name.host.log"
    if ! SCRIMP_V130_ROTATING_IC_QUALIFICATION=restart_parity \
        SCRIMP_V130_ISSUE6_AUDIT_PORTS=1 \
        bash scripts/run-v130-rotating-controlled-ic.sh "${args[@]}" > "$log" 2>&1; then
        echo "FAIL original native model case $name; preserved $log" >&2
        tail -n 120 "$log" >&2
        return 1
    fi
    local path
    path="$(sed -n 's/^Evidence on host: //p' "$log" | tail -n 1)"
    [[ -n "$path" && "$path" == "$root/sandbox/outputs/"* \
        && -s "$path/initial_state_report.json" \
        && -s "$path/accepted_state.npz" \
        && -s "$path/initial_state_trace.npz" \
        && -s "$path/mpi_physical_powers.png" ]] || {
        echo "FAIL: incomplete original active MPI evidence: $name" >&2
        return 1
    }
    printf '%s\t%s\n' "$name" "$path" >> "$out/native_cases.tsv"
    grep -aE '^ROTATING CONTROLLED IC ACCEPTED|^HOST EVIDENCE PASS|^MPI #6 PORT AUDIT PASS' "$log" >&2
    printf '%s\n' "$path"
}

for mesh in 0.2 0.1; do
    label="h${mesh/./p}"
    # One source at MPI1 and one at MPI4, both physically compatible.
    seed1="$(run_native "${label}_seed_mpi1" seed "$mesh" 1 0.02)"
    seed4="$(run_native "${label}_seed_mpi4" seed "$mesh" 4 0.02)"

    # The same original active controlled SWE on MPI1/2/4.
    cont1="$(run_native "${label}_continuous_mpi1" seed "$mesh" 1 0.024)"
    cont2="$(run_native "${label}_continuous_mpi2" seed "$mesh" 2 0.024)"
    cont4="$(run_native "${label}_continuous_mpi4" seed "$mesh" 4 0.024)"

    # Exercise a changed rank decomposition in both directions,
    # including the reassembled port powers on all target ranks.
    restart14="$(run_native "${label}_restart_1to4" resume "$mesh" 4 0.024 \
        "$seed1/accepted_state.npz")"
    restart41="$(run_native "${label}_restart_4to1" resume "$mesh" 1 0.024 \
        "$seed4/accepted_state.npz")"

    python3 sandbox/compare_v130_issue6_active_mpi.py \
        --seed1 "$seed1" --seed4 "$seed4" \
        --cont1 "$cont1" --cont2 "$cont2" --cont4 "$cont4" \
        --resume1to4 "$restart14" --resume4to1 "$restart41" \
        --report "$out/${label}_native_mpi4.json" |
        tee "$out/${label}_native_mpi4.comparison.log"
done

python3 - "$out" <<'PY'
import json, sys
from pathlib import Path
out=Path(sys.argv[1])
files=sorted(out.glob("h*_native_mpi4.json"))
assert len(files)==2, f"Need both independent original meshes: {files}"
rows=[json.loads(p.read_text(encoding="utf-8")) for p in files]
assert {x["mesh"] for x in rows}=={.1,.2}
assert all(x["status"]=="PASS_NATIVE_MPI_1_2_4_ACTIVE_ROTATING_SWE" for x in rows)
assert all(x["all_rank_local_port_files_checked"]==17 for x in rows)
manifest=json.loads((out/"manifest.json").read_text(encoding="utf-8"))
manifest["scientific_status"]="PASS_NATIVE_COMPATIBLE_ROTATING_SWE_MPI_1_2_4"
manifest["native_case_count"]=14
manifest["independent_mpi1_2_4_reference_comparisons"]=4
manifest["native_mpi4_cold_restart_comparisons"]=4
manifest["checked_rank_unique_original_power_files"]=34
manifest["max_full_FE_L2_cross_MPI"]=max(
    r["max_full_path_FE_L2_MPI"] for r in rows
)
manifest["max_original_h_weighted_physical_trace_difference"]=max(
    r["max_original_physical_power_observable_difference_MPI"] for r in rows
)
manifest["max_generic_port_power_difference"]=max(
    r["max_generic_port_power_difference_MPI"] for r in rows
)
(out/"summary.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
print(json.dumps(manifest,indent=2))
PY

echo "COMPLETE native MPI #6 active rotating evidence: $out"
