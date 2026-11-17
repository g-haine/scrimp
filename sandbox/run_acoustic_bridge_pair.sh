#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Native Docker A/B test for the opt-in distributed GetFEM -> PETSc COO bridge.
# One fixed-step CN acoustic benchmark and one complete complex-trajectory
# validation for each mode (ijv and coo). NO host Conda simulation.
#
#   bash sandbox/run_acoustic_bridge_pair.sh 0.08 2
#   bash sandbox/run_acoustic_bridge_pair.sh 0.04 2
#
# Run using bash, never source. Requires scrimp-mpi:local and sudo docker.

set -euo pipefail

mesh="${1:-0.08}"
ranks="${2:-2}"
if [[ ! "$mesh" =~ ^0\.[0-9]+$ ]] || \
   [[ ! "$ranks" =~ ^(1|2|4|6)$ ]]; then
    echo "Usage: bash sandbox/run_acoustic_bridge_pair.sh 0.08 [1|2|4|6]" >&2
    exit 2
fi

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
short_sha="$(git rev-parse --short=8 HEAD)"
out="$root/sandbox/outputs/bridge-pair-h${mesh//./p}-${ranks}r-${short_sha}"
mkdir -p "$out"

pcores=(0 2 4 6 8 10)
cpus="$(IFS=,; echo "${pcores[*]:0:ranks}")"

if ! sudo -v; then
    echo "sudo authentication failed" >&2
    exit 1
fi

run_docker() {
    local method="$1"
    shift
    sudo docker run --rm \
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
            -u SCRIMP_MUMPS_RESET_POLICY \
            SCRIMP_GETFEM_MATRIX_MODE=distributed \
            SCRIMP_GETFEM_BRIDGE_DISTRIBUTED="$method" \
        mpiexec --bind-to none -n "$ranks" \
        python -u "$@"
}

echo "SCRIMP commit $short_sha | GetFEM 5.5 Docker | h=$mesh MPI=$ranks cpus=$cpus"
echo "Results: $out"

# The IJV baseline and candidate must use identical source and solver
# options. 'auto' labels SCRIMP's unchanged automatic MUMPS reset policy.
for method in ijv coo; do
    label="wave_${method}"
    logfile="$out/${label}.log"
    report="$out/${label}.json"
    echo "BENCH: $method"
    if ! run_docker "$method" sandbox/bench_acoustic_2d.py \
        --mesh "$mesh" --dt 0.002 --final-time 0.04 \
        --save-dt 0.01 --series auto --profile \
        --output "/artifacts/${label}.json" \
        >"$logfile" 2>&1; then
        echo "FAILED benchmark $method: $logfile" >&2
        tail -n 60 "$logfile" >&2
        exit 1
    fi
    test -s "$report" || {
        echo "Missing report $report" >&2
        exit 1
    }
done

# Compare the ENTIRE saved complex trajectory: times, all coefficients,
# rank-to-rank parity and reference-vs-candidate. No norms-only shortcut.
for method in ijv coo; do
    label="state_${method}"
    logfile="$out/${label}.log"
    echo "TRAJECTORY: $method"
    extra=()
    if [[ "$method" == coo ]]; then
        extra=(--reference /artifacts/state_ijv.npz)
    fi
    if ! run_docker "$method" sandbox/check_acoustic_state_parity.py \
        --mesh "$mesh" --dt 0.002 --final-time 0.04 \
        --output "/artifacts/${label}.npz" "${extra[@]}" \
        >"$logfile" 2>&1; then
        echo "FAILED full trajectory parity $method: $logfile" >&2
        tail -n 60 "$logfile" >&2
        exit 1
    fi
done

# Host /usr/bin/python3 only parses the JSON results from Docker. Its
# libraries are never used for simulations or GetFEM/PETSc.
if ! /usr/bin/python3 - "$out/wave_ijv.json" "$out/wave_coo.json" <<'PY'
import json
import math
import sys
from pathlib import Path

base, candidate = [
    json.loads(Path(path).read_text(encoding="utf-8"))
    for path in sys.argv[1:3]
]

for report, label in ((base, "ijv"), (candidate, "coo")):
    assert report["effective_mumps_reset_policy"] == "never", label
    assert report["ts_type"] == "cn", label
    assert report["ts_adapt_type"] == "none", label
    assert report["ts_converged_reason"] > 0, label
    assert report["ts_steps"] == 20, label
    assert report["ts_rejected_steps"] == 0, label
    assert report["ts_snes_failures"] == 0, label
    assert abs(report["actual_final_time"] - 0.04) < 1e-12, label
    assert abs(report["last_saved_time"] - 0.04) < 1e-12, label
    assert all(
        math.isfinite(float(v))
        for v in report["final_variable_norms"].values()
    ), label

assert base["degrees_of_freedom"] == candidate["degrees_of_freedom"]
assert base["mpi_size"] == candidate["mpi_size"]
assert base["rectangle_mesh_h"] == candidate["rectangle_mesh_h"]
assert base["final_variable_norms"].keys() == candidate["final_variable_norms"].keys()
for name, v in base["final_variable_norms"].items():
    w = candidate["final_variable_norms"][name]
    assert math.isclose(v, w, abs_tol=1e-11, rel_tol=1e-10), (name, v, w)

events = (
    "setup.bridge.mass",
    "setup.bridge.stiffness",
    "setup.bridge.mass.petsc.coo.prepare",
    "setup.bridge.mass.petsc.coo.preallocate",
    "setup.bridge.mass.petsc.coo.set_values",
    "setup.bridge.stiffness.petsc.coo.prepare",
    "setup.bridge.stiffness.petsc.coo.preallocate",
    "setup.bridge.stiffness.petsc.coo.set_values",
)
print("\n=== End-to-end GetFEM/PETSc A/B (one run each) ===")
print("DOF:", base["degrees_of_freedom"], "MPI:", base["mpi_size"])
for label, report in (("ijv", base), ("coo", candidate)):
    print(
        f"{label:3s} solve={report['wall_solve_seconds_max']:.6f}s "
        f"RSS/rank={report['peak_rank_rss_kib']/1024:.1f} MiB"
    )
    prof = report["profile"]["events"]
    for event in events:
        if event in prof:
            print(
                f"  {event:48s} "
                f"{prof[event]['seconds_rank_max']:.6f}s"
            )
print("Complete complex trajectory parity: PASS (checked in Docker)")
print("Final variable norms: PASS")
PY
then
    echo "FAILED result validation" >&2
    exit 1
fi

echo "SUCCESS: all JSON and complete complex trajectory checks passed."
echo "Logs and reports: $out"
