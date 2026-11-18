#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See the LICENSE file for license information.
#
# Diagnostic only. Replays legacy GetFEM state synchronization after
# accepted initial pseudo-TS steps, in the current SCRIMP WITHOUT
# changing its production solver. Compare against SAVED v1.2.0 and
# current reference states; do not rerun the 70s release baseline.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
historical="${1:-$root/sandbox/outputs/v130-wave-h0p04-97bb8585-mAO3U3}"
mesh="0.04"
if [[ $# -gt 1 ]]; then
    echo "Usage: bash sandbox/run_wave_initial_sync_probe.sh [path-to-v130-wave-h0p04-97bb8585-mAO3U3]" >&2
    exit 2
fi
historical="$(realpath "$historical")"
for required in \
    release_120_release_native_1r_rep1_state.npz \
    release_120_release_native_1r_rep1.json \
    current_on_change_1r_rep1_state.npz \
    current_on_change_1r_rep1.json; do
    if [[ ! -s "$historical/$required" ]]; then
        echo "Missing historical artifact: $historical/$required" >&2
        echo "The previous v130-wave archive is sufficient; no old run needed." >&2
        exit 2
    fi
done
mkdir -p "$root/sandbox/outputs"
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/wave-initial-sync-diagnosis-h0p04-${sha:0:8}-XXXXXX")"
echo "SCRIMP wave initial pseudo-TS sync diagnosis"
echo "Current=$sha old/current reference=$historical"
echo "Original v1.2.0 binaries are read-only; ONLY current SCRIMP is run"
echo "Output: $out"
sudo -v
if ! sudo docker run --rm \
    --user "$(id -u):$(id -g)" \
    --cpuset-cpus=0 --workdir /workspace/scrimp \
    -e HOME=/tmp -e MPLBACKEND=Agg \
    -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
    -e SCRIMP_GETFEM_MATRIX_MODE=distributed \
    -e SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
    -e SCRIMP_MATRIX_MEMORY_MODE=compact \
    -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
    -v "$root/examples:/workspace/scrimp/examples:ro" \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$out:/artifacts" \
    scrimp-mpi:local timeout -k 10s 180s \
    env -u SCRIMP_MUMPS_RESET_POLICY \
        -u SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN \
    mpiexec --bind-to none -n 1 \
    python -u sandbox/bench_wave_release_scaling.py \
        --variant current --policy on_change \
        --probe-pseudo-state-sync \
        --mesh "$mesh" --dt 0.002 --final-time 0.02 \
        --save-dt 0.01 --source-sha "$sha" --repetition 1 \
        --out /artifacts/pseudo_sync_probe.json \
    > "$out/pseudo_sync_probe.log" 2>&1; then
    echo "FAILED: current-only sync probe" >&2
    tail -n 70 "$out/pseudo_sync_probe.log" >&2
    exit 1
fi
grep -aE '^(WAVE-PSEUDO-SYNC |PASS current/)' "$out/pseudo_sync_probe.log" || true

sudo docker run --rm -i \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp -e MPLBACKEND=Agg \
    -v "$historical:/historical:ro" \
    -v "$out:/artifacts" \
    scrimp-mpi:local python -u - <<'PY'
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base = Path("/historical")
out = Path("/artifacts")
cases = {
    "v1.2.0": base / "release_120_release_native_1r_rep1_state.npz",
    "v1.3.0 current": base / "current_on_change_1r_rep1_state.npz",
    "v1.3.0 pseudo-sync probe": out / "pseudo_sync_probe_state.npz",
}
meta = json.loads((out / "pseudo_sync_probe.json").read_text())
if not meta.get("initial_sync_probe", {}).get("enabled"):
    raise SystemExit("FAIL: initial pseudo sync probe was not enabled")
events = meta["initial_sync_probe"]["accepted_step_events"]
if not any(event["step"] >= 1 for event in events):
    raise SystemExit("FAIL: no accepted pseudo-initialization monitor captured")
intervals = events[-1]["variables"]
arrays = {}
for name, file in cases.items():
    with np.load(file) as z:
        arrays[name] = {kind: np.array(z[kind]) for kind in ("initial", "final")}
    if any(len(v) != meta["dofs"] for v in arrays[name].values()):
        raise SystemExit("FAIL: mismatching global DOF count")

def distance(a, b):
    d = a - b
    return {
        "abs_l2": float(np.linalg.norm(d)),
        "scaled_l2": float(np.linalg.norm(d) / max(1., np.linalg.norm(a))),
        "scaled_linf": float(np.max(abs(d))/max(1.,np.max(abs(a)))),
    }

summary={"current_commit":meta["source_sha"],
         "historical":"v1.2.0 native rep1 and current rep1, same 0.04 mesh",
         "accepted_pseudo_events":events,
         "full_state_distances":{},
         "named_field_intervals":intervals,
         "named_field_norms":{}}
v120=arrays["v1.2.0"]
for k in ("initial", "final"):
    for name in ("v1.3.0 current", "v1.3.0 pseudo-sync probe"):
        stats=distance(v120[k],arrays[name][k])
        summary["full_state_distances"][f"{k}: v1.2.0 vs {name}"]=stats
        print("STATE-PARITY",k,name,stats)
for var, detail in intervals.items():
    lo,hi=detail["range"]
    summary["named_field_norms"][var]={}
    for kind in ("initial", "final"):
        summary["named_field_norms"][var][kind]={
           name:float(np.linalg.norm(vals[kind][lo:hi]))
           for name,vals in arrays.items()
        }
    print("FIELD-INITIAL",var,
          summary["named_field_norms"][var]["initial"])

summary_file=out/"wave_initial_sync_diagnosis.json"
summary_file.write_text(json.dumps(summary,indent=2)+"\n")

png=out/"png";svg=out/"svg";png.mkdir(exist_ok=True);svg.mkdir(exist_ok=True)
names=list(intervals)
x=np.arange(len(names)); width=0.25
fig,ax=plt.subplots(figsize=(10,5.5))
for i,(case,values) in enumerate(arrays.items()):
    ax.bar(x+(i-1)*width, [
       np.linalg.norm(values["initial"][
          slice(*intervals[var]["range"])
       ]) for var in names
    ], width, label=case)
ax.set_xticks(x,names)
ax.set_ylabel("Euclidean norm of initial state block")
ax.set_title("Wave.py: legacy vs current pseudo-initial-state synchronization")
ax.legend()
ax.grid(axis="y",alpha=.25)
fig.tight_layout()
fig.savefig(png/"wave_initial_named_blocks.png",dpi=190)
fig.savefig(svg/"wave_initial_named_blocks.svg")
plt.close(fig)
print("REPORT:",summary_file)
print("HUMAN PLOT: png/wave_initial_named_blocks.png")
print("Interpretation: reproducing old vectors does NOT establish physical correctness.")
PY
echo "Reports: $out"
