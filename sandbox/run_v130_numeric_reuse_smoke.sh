#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See the LICENSE file for license information.
#
# Short production on_change qualifier, current version ONLY.
# Avoids rerunning the expensive historical v1.2.0 baseline until
# numeric reuse and complete 1-rank/2-rank state parity have passed.
set -euo pipefail

mesh="${1:-0.04}"
ranks_csv="${2:-1,2}"
if [[ $# -gt 2 || ! "$mesh" =~ ^0\.[0-9]+$ ||
      ! "$ranks_csv" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ]]; then
    echo "Usage: bash sandbox/run_v130_numeric_reuse_smoke.sh [mesh=0.04] [ranks=1,2]" >&2
    exit 2
fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse HEAD)"
mkdir -p "$root/sandbox/outputs"
out="$(mktemp -d "$root/sandbox/outputs/v130-numeric-smoke-h${mesh//./p}-${sha:0:8}-XXXXXX")"
IFS=, read -r -a ranks <<< "$ranks_csv"
IFS=, read -r -a cpus <<< "${SCRIMP_SCALING_CPUS:-0,2,4,6,8,10}"
timeout_s="${SCRIMP_BENCH_TIMEOUT_S:-900}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] || (( timeout_s < 60 || timeout_s > 7200 )); then
    echo "SCRIMP_BENCH_TIMEOUT_S must be in 60..7200" >&2
    exit 2
fi
echo "SCRIMP v1.3.0 actual production on_change numeric LU smoke"
echo "Commit=$sha mesh=$mesh MPI ranks=$ranks_csv; no v1.2.0 execution"
echo "Output: $out"
sudo -v
for n in "${ranks[@]}"; do
    if (( n > ${#cpus[@]} )); then echo "Too few CPUs for $n ranks" >&2; exit 2; fi
    affinity="$(IFS=,; echo "${cpus[*]:0:n}")"
    label="current_on_change_${n}r"
    echo "RUN current on_change ranks=$n cpus=$affinity (full PETSc stages)"
    if ! sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        --cpuset-cpus="$affinity" --workdir /workspace/scrimp \
        -e HOME=/tmp -e MPLBACKEND=Agg \
        -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
        -e SCRIMP_GETFEM_MATRIX_MODE=distributed \
        -e SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
        -e SCRIMP_MATRIX_MEMORY_MODE=compact \
        -e SCRIMP_MUMPS_REUSE_DIAG=1 \
        -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
        -v "$root/examples:/workspace/scrimp/examples:ro" \
        -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
        -v "$out:/workspace/scrimp/outputs" \
        -v "$out:/artifacts" \
        scrimp-mpi:local timeout -k 10s "${timeout_s}s" \
        env -u SCRIMP_MUMPS_RESET_POLICY \
            -u SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN \
        mpiexec --bind-to none -n "$n" \
        python -u sandbox/bench_wave_release_scaling.py \
            --variant current --policy on_change --mesh "$mesh" \
            --dt 0.002 --final-time 0.02 --save-dt 0.01 \
            --source-sha "$sha" --repetition 1 \
            --profile --petsc-log \
            --out "/artifacts/${label}.json" \
        > "$out/${label}.log" 2>&1; then
        echo "NON-PASS: $label" >&2
        tail -n 70 "$out/${label}.log" >&2
        exit 1
    fi
    grep -aE '^PASS current/' "$out/${label}.log" || true
    echo "  Production on_change decisions (MAIN CN):"
    grep -a '^MUMPS-REUSE-AUDIT phase=main ' "$out/${label}.log" || true
    if ! grep -aq '^MatLUFactorNum' "$out/${label}.log"; then
        echo "NON-PASS: missing PETSc numeric LU events ($label)" >&2
        exit 1
    fi
done
sudo docker run --rm -i \
    --user "$(id -u):$(id -g)" -e HOME=/tmp -e MPLBACKEND=Agg \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$out:/artifacts" \
    scrimp-mpi:local python -u - "$ranks_csv" <<'PY'
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "/workspace/scrimp/sandbox")
from analyze_v130_wave_benchmark import _petsc_events, _state_distance, publish_smoke_figures
root = Path("/artifacts")
ref = None
for n in [int(x) for x in sys.argv[1].split(",")]:
    label = f"current_on_change_{n}r"
    data = json.loads((root / f"{label}.json").read_text(encoding="utf-8"))
    if (data["mpi_ranks"] != n or data["steps"] != 10
            or data["rejections"] or data["snes_failures"]
            or data["mumps_policy"] != "on_change"
            or data["main_cn_mumps_icntl28"] != (2 if n>1 else 1)):
        raise SystemExit(f"FAIL {label}: invalid solver settings / TS steps")
    native = _petsc_events(root / f"{label}.log")
    state = np.load(root / f"{label}_state.npz")
    initial, final = state["initial"], state["final"]
    if ref is None:
        ref = (initial, final)
    else:
        _state_distance(ref[0], initial, f"{label}: initial state")
        _state_distance(ref[1], final, f"{label}: final state")
    print(f"PASS {label}: main CN numeric LU calls="
          f"{native['main.MatLUFactorNum']['calls']}, "
          f"symbolic={native['main.MatLUFactorSym']['calls']}, "
          f"MatSolve={native['main.MatSolve']['calls']}, "
          f"solve={data['wall_solve_s']:.5f}s")
print("PASS production on_change exact numeric reuse and complete MPI state parity")
plots = publish_smoke_figures(root, [int(x) for x in sys.argv[1].split(",")])
print("PNG figures:", ", ".join("png/" + name + ".png" for name in plots))
PY
echo "Human-readable figures: $out/png (SVG: $out/svg)"
echo "Reports: $out"
