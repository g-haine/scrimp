#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.
#
# PETSc TS / MUMPS numeric-LU reuse: full wave, baseline vs MAIN CN static.
# No edits to production and NEVER reset policy exclusively.
set -euo pipefail
mesh="${1:-0.04}"
ranks="${2:-2}"
if [[ $# -gt 2 || ! "$mesh" =~ ^0\.[0-9]+$ ||
      ! "$ranks" =~ ^[1-9][0-9]*$ ]]; then
    echo "Usage: bash sandbox/run_wave_mumps_numeric_reuse_probe.sh [mesh=0.04] [ranks=2]" >&2
    exit 2
fi
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/wave-mumps-numeric-reuse-h${mesh//./p}-$sha-XXXXXX")"
IFS=, read -r -a cpu_pool <<< "${SCRIMP_SCALING_CPUS:-0,2,4,6,8,10}"
if (( ranks > ${#cpu_pool[@]} )); then
    echo "Insufficient CPU affinity pool for MPI=$ranks" >&2; exit 2
fi
cpus="$(IFS=,; echo "${cpu_pool[*]:0:ranks}")"
timeout_s="${SCRIMP_WAVE_CN_TIMEOUT_S:-900}"
if [[ ! "$timeout_s" =~ ^[0-9]+$ ]] || (( timeout_s < 60 || timeout_s > 7200 )); then
    echo "Invalid SCRIMP_WAVE_CN_TIMEOUT_S" >&2; exit 2
fi
echo "SCRIMP numeric MUMPS reuse probe: full CN wave with strict accuracy gates"
echo "Commit=$sha mesh=$mesh MPI=$ranks cpus=$cpus MUMPS ordering=parallel reset=never ONLY"
echo "Reference=baseline (10 numeric factorizations observed earlier), trial=static_cn"
echo "Output: $out"
sudo -v
for mode in baseline static_cn; do
    echo "RUN mode=$mode MPI=$ranks"
    if sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        --cpuset-cpus="$cpus" \
        --workdir /workspace/scrimp \
        -e HOME=/tmp -e MPLBACKEND=Agg \
        -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 \
        -e SCRIMP_GETFEM_MATRIX_MODE=distributed \
        -e SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
        -e SCRIMP_MATRIX_MEMORY_MODE=compact \
        -e SCRIMP_MUMPS_RESET_POLICY=never \
        -e SCRIMP_BENCH_GIT_SHA="$sha" \
        -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
        -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
        -v "$root/examples:/workspace/scrimp/examples:ro" \
        -v "$out:/workspace/scrimp/outputs" \
        -v "$out:/artifacts" \
        scrimp-mpi:local \
        timeout -k 10s "${timeout_s}s" \
        mpiexec --bind-to none -n "$ranks" \
        python -u sandbox/bench_wave_mumps_ordering.py \
            --mesh "$mesh" --dt 0.002 --final-time 0.02 \
            --save-dt 0.01 --series never --profile --petsc-log \
            --ordering parallel --factor-reuse "$mode" \
            --output "/artifacts/$mode.json" \
        > "$out/$mode.log" 2>&1; then
        grep -aE '^MAIN-TS-MUMPS |^INITIAL-TS-MUMPS |^FACTOR-REUSE |^Elapsed time:|^Steps:' "$out/$mode.log" || true
        echo "  PETSc LU event summary (initial stage then main CN):"
        grep -aE '^MatLUFactorSym|^MatLUFactorNum|^MatSolve ' "$out/$mode.log" || true
    else
        echo "NON-PASS mode=$mode" >&2
        tail -n 60 "$out/$mode.log" >&2
        echo "STOP: examine $out" >&2
        exit 1
    fi
done
/usr/bin/python3 - "$out/baseline.json" "$out/static_cn.json" "$out/baseline.log" "$out/static_cn.log" <<'PY'
import json, math, sys
base, static = (json.load(open(p, encoding="utf-8")) for p in sys.argv[1:3])
for report in (base, static):
    if (report["effective_mumps_reset_policy"] != "never"
            or report["ts_steps"] != 10
            or report["ts_rejected_steps"] or report["ts_snes_failures"]
            or not report["mumps_icntl28_verified"]):
        raise SystemExit("Invalid TS/PC baseline metadata")
if not static["numeric_factor_reuse_main_ts_options"]:
    raise SystemExit("static_cn options were not applied")
norms1, norms2 = base["final_variable_norms"], static["final_variable_norms"]
if set(norms1) != set(norms2):
    raise SystemExit("Mismatched physical variable sets")
reldiff = {k: abs(norms2[k]-norms1[k]) / max(1.0, abs(norms1[k]))
           for k in norms1}
print("COMPARISON same ranks, same CN: normalized absolute norm diffs", reldiff)
if any(not math.isfinite(v) or v > 5e-8 for v in reldiff.values()):
    raise SystemExit("FAIL: state norm gate (full vector/power budget not yet compared)")
print("  baseline solve:", base["wall_solve_seconds_max"], "seconds")
print("  static_cn solve:", static["wall_solve_seconds_max"], "seconds")
def main_numeric_count(log):
    stage = None
    found = []
    for line in open(log, encoding="utf-8", errors="replace"):
        if "Event Stage" in line and ("WAVE_MAIN_CN" in line
                                      or "WAVE_INITIAL_CONSISTENCY" in line
                                      or "Main Stage" in line):
            stage = "main" if "WAVE_MAIN_CN" in line else "other"
            continue
        parts = line.strip().split()
        if len(parts) > 1 and parts[0] == "MatLUFactorNum" and stage == "main":
            found.append(int(parts[1]))
    if len(found) != 1:
        raise SystemExit(
            f"FAIL: expected exactly one MAIN CN MatLUFactorNum event "
            f"in {log}, got {found}; inspect the PETSc log stage names"
        )
    return found[0]
baseline_count = main_numeric_count(sys.argv[3])
reuse_count = main_numeric_count(sys.argv[4])
print("Main CN MatLUFactorNum counts: baseline=", baseline_count,
      " static_cn=", reuse_count)
if baseline_count != 10 or reuse_count != 1:
    raise SystemExit(
        "FAIL: actual numeric LU reuse not confirmed (expected 10 -> 1)"
    )
print("NOTE: norm agreement does not replace full-state and power-balance checks.")
PY
echo "Reports: $out"
echo "PASS: main CN numerical LU reuse (10 -> 1) confirmed by PETSc event counts"
