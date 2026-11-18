#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See the LICENSE file for license information.
#
# Native production regression: single current wave solve, archived
# unchanged v1.2.0 state, native MUMPS event counts and human-readable
# initial/final full-vector parity plot. NO legacy simulation rerun.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
history="${1:-$root/sandbox/outputs/v130-wave-h0p04-97bb8585-mAO3U3}"
if [[ $# -gt 1 ]]; then
    echo "Usage: bash sandbox/run_wave_initial_state_regression.sh [archived-v130-wave-output]" >&2
    exit 2
fi
history="$(realpath "$history")"
for name in \
    release_120_release_native_1r_rep1_state.npz \
    release_120_release_native_1r_rep1.json; do
    if [[ ! -s "$history/$name" ]]; then
        echo "Missing v1.2.0 archived reference: $history/$name" >&2
        exit 2
    fi
done
mkdir -p "$root/sandbox/outputs"
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/wave-pseudo-state-parity-${sha:0:8}-XXXXXX")"
echo "SCRIMP wave production pseudo-state parity — current only"
echo "Commit=$sha h=0.04 1 MPI; v1.2.0 reference=$history"
echo "Output: $out"
sudo -v

if ! sudo docker run --rm \
    --user "$(id -u):$(id -g)" --cpuset-cpus=0 \
    --workdir /workspace/scrimp \
    -e HOME=/tmp -e MPLBACKEND=Agg \
    -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
    -e SCRIMP_GETFEM_MATRIX_MODE=distributed \
    -e SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
    -e SCRIMP_MATRIX_MEMORY_MODE=compact \
    -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
    -v "$root/examples:/workspace/scrimp/examples:ro" \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$out:/workspace/scrimp/outputs" \
    -v "$out:/artifacts" \
    scrimp-mpi:local timeout -k 10s 180s \
    env -u SCRIMP_MUMPS_RESET_POLICY \
        -u SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN \
        -u SCRIMP_NONLINEAR_RESIDUAL_FORM \
    mpiexec --bind-to none -n 1 \
    python -u sandbox/bench_wave_release_scaling.py \
        --variant current --policy on_change \
        --mesh 0.04 --dt 0.002 --final-time 0.02 \
        --save-dt 0.01 --source-sha "$sha" --repetition 1 \
        --profile --petsc-log \
        --out /artifacts/current_production.json \
    > "$out/current_production.log" 2>&1; then
    echo "FAIL current production run" >&2
    tail -n 90 "$out/current_production.log" >&2
    echo "Reports: $out" >&2
    exit 1
fi
grep -a '^PASS current/' "$out/current_production.log" || true

sudo docker run --rm -i \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp -e MPLBACKEND=Agg \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$history:/historical:ro" \
    -v "$out:/artifacts" \
    scrimp-mpi:local python -u - <<'PY'
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/workspace/scrimp/sandbox")
from analyze_v130_wave_benchmark import (
    _expected_lu_updates, _petsc_events, _state_distance
)

d = Path("/artifacts")
h = Path("/historical")
data = json.loads((d / "current_production.json").read_text())
historical = json.loads((h / "release_120_release_native_1r_rep1.json").read_text())
if (data["variant"] != "current" or data["mpi_ranks"] != 1
        or data["steps"] != 10 or data["rejections"] or data["snes_failures"]
        or data["mumps_policy"] != "on_change"
        or data["main_cn_mumps_icntl28"] != 1
        or data["dofs"] != historical["dofs"]):
    raise SystemExit("FAIL: invalid current wave model/TS/solver metadata")
n, m = _expected_lu_updates(data)
events = _petsc_events(
    d / "current_production.log",
    expected_numeric=n,
    expected_symbolic=m,
)
with np.load(h / "release_120_release_native_1r_rep1_state.npz") as z:
    original = {key: np.asarray(z[key]) for key in ("initial", "final")}
with np.load(d / "current_production_state.npz") as z:
    candidate = {key: np.asarray(z[key]) for key in ("initial", "final")}

errors = {}
for key in ("initial", "final"):
    ref = original[key]
    cur = candidate[key]
    _state_distance(ref, cur, "current vs native v1.2.0 " + key)
    err = cur - ref
    abs_l2 = float(np.linalg.norm(err))
    scaled_l2 = abs_l2 / max(1.0, float(np.linalg.norm(ref)))
    scaled_linf = float(np.max(np.abs(err))) / max(1.0, float(np.max(np.abs(ref))))
    errors[key] = {
        "absolute_l2": abs_l2, "scaled_l2": scaled_l2,
        "scaled_linf": scaled_linf,
    }
    print(f"STATE-PARITY {key}: L2={abs_l2:.9g}, scaled_L2={scaled_l2:.9g}, "
          f"scaled_Linf={scaled_linf:.9g}")
print("PETSc MUMPS main TS:", events)
report = {
    "source_sha": data["source_sha"],
    "reference_source_sha": historical["source_sha"],
    "full_vector_tolerance": 5e-8,
    "full_initial_and_final_state_parity": True,
    "mumps_native_events_match_exact_operator_changes": True,
    "errors": errors,
    "main_lu_events": events,
    "current_solve_wall_s": data["wall_solve_s"],
    "historical_solve_wall_s": historical["wall_solve_s"],
    "warning": ("Full-state parity is necessary but not sufficient to "
                "validate the DAE initial constraints and Hamiltonian balance."),
}
(d / "wave_initial_state_regression.json").write_text(
    json.dumps(report, indent=2) + "\n", encoding="utf-8"
)
png = d / "png"; png.mkdir(exist_ok=True)
svg = d / "svg"; svg.mkdir(exist_ok=True)
fig, ax = plt.subplots(figsize=(8.5, 5.0))
bars = ax.bar(
    ("Initial state", "Final state"),
    [max(errors[x]["scaled_l2"], 1e-18) for x in ("initial", "final")],
    label="Scaled full-vector L2 error",
)
ax.axhline(5e-8, linestyle="--", label="Parity tolerance 5e-8")
ax.set_yscale("log")
ax.set_ylabel("Scaled L2 error (log scale; zeros plotted at 1e-18)")
ax.set_title("wave.py — v1.2.0 vs corrected current initial state")
ax.legend()
ax.grid(axis="y", alpha=.25)
fig.tight_layout()
fig.savefig(png / "wave_v120_current_full_state_parity.png", dpi=190)
fig.savefig(svg / "wave_v120_current_full_state_parity.svg")
plt.close(fig)
print("PASS complete state parity and native MUMPS reuse")
print("PNG: png/wave_v120_current_full_state_parity.png")
PY

echo "Reports: $out"
