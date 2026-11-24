#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.
#
# Snapshot EXACT P=K+shift*M before the first genuine CN factorization,
# compare MUMPS/SuperLU_DIST, analyze off-line with NON-REGULARIZING
# SciPy SuperLU. This is NOT a solve, and cannot qualify SuperLU_DIST.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
  echo "Unqualified native GetFEM/PETSc image: $revision" >&2
  exit 2
}
python3 -m py_compile \
  sandbox/probe_v130_wave2d_factor_matrix.py \
  sandbox/analyze_v130_wave2d_factor_matrix.py \
  tests/test_v130_wave2d_factor_matrix.py \
  examples/wave.py
python3 tests/test_v130_wave2d_factor_matrix.py
python3 scripts/check-v130-published-immutable.py
python3 -c 'import scipy, numpy'
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-wave2d-first-cn-jacobian-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-wave2d-first-CN-sparse-operator-audit-v1","commit":"%s","image":"%s","image_revision":"%s","shift_expected":250,"MPI":1,"model":"examples.wave","not_a_solver_acceptance":true,"previous_zero_pivot_evidence":"v130-wave2d-superlu-dist-ordering-8cb2a51f-422NJxVt.zip"}\n' \
  "$(git rev-parse HEAD)" "$image" "$revision" > "$out/manifest.json"
echo "FIRST CN Jacobian: capture actual PETSc PC P and original K, M"
echo "Same 2230-DOF 2D anisotropic wave, shift 250, dt 0.008"
echo "MUMPS control + failing SuperLU_DIST (MPI1 only); no parameter changes"
echo "No artificial pivots, matrix shifts, T changes or LU tolerances"
echo "EVIDENCE: $out"
for solver in mumps superlu_dist; do
  echo "RUN first original $solver CN step"
  mkdir -p "$out/$solver" "$out/work_$solver"
  rc=0
  "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp --env MPLBACKEND=Agg \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
    --env SCRIMP_MPI_TRACE=1 \
    --env SCRIMP_TS_FAILURE_DIAGNOSTICS=1 \
    --env SCRIMP_V130_WAVE2D_CAPTURE_DIR=/artifacts \
    --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
    --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
    --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
    --env SCRIMP_MATRIX_MEMORY_MODE=compact \
    --env SCRIMP_MUMPS_RESET_POLICY=on_change \
    --env SCRIMP_GETFEM_RHS_MODE=auto \
    --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
    --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
    --volume "$root/examples/wave.py:/workspace/scrimp/examples/wave.py:ro" \
    --volume "$root/sandbox/probe_v130_wave2d_factor_matrix.py:/workspace/scrimp/sandbox/probe_v130_wave2d_factor_matrix.py:ro" \
    --volume "$out/$solver:/artifacts:rw" \
    --volume "$out/work_$solver:/scrimp-work:rw" \
    --workdir /scrimp-work \
    "$image" timeout -k 5s "${SCRIMP_V130_WAVE2D_FACTOR_TIMEOUT:-90s}" \
    mpiexec --oversubscribe --bind-to none -n 1 \
    python -u /workspace/scrimp/sandbox/probe_v130_wave2d_factor_matrix.py \
      --solver "$solver" --mesh-size 0.2 \
      --final-time 0.008 --time-step 0.008 \
      --check --output /artifacts \
    > "$out/$solver.host.log" 2>&1 || rc=$?
  printf '%s\n' "$rc" > "$out/$solver.exitcode"
  echo "EXIT $solver=$rc captured=$([[ -r "$out/$solver/first_main_CN_factor.npz" ]] && echo yes || echo no)"
  if [[ "$rc" != 0 ]]; then
    grep -aE 'MAIN TS FAILURE DIAGNOSTIC|CAPTURED_FIRST_REAL_CN_PC_MATRIX|Traceback|zero_diag_ids' \
      "$out/$solver.host.log" | tail -n 10 >&2 || true
  fi
done
analyze_rc=0
python3 sandbox/analyze_v130_wave2d_factor_matrix.py \
  --mumps "$out/mumps" --superlu-dist "$out/superlu_dist" \
  --output "$out/sparse_matrix_analysis.json" \
  > "$out/analyzer.host.log" 2>&1 || analyze_rc=$?
printf '%s\n' "$analyze_rc" > "$out/analyzer.exitcode"
python3 - "$out" <<'PY'
import json
from pathlib import Path
import sys
import zipfile
folder=Path(sys.argv[1])
exits={name:int((folder/f"{name}.exitcode").read_text())
       for name in ("mumps","superlu_dist")}
captured={
    name: all((folder/name/f"first_main_CN_factor.{ext}").is_file()
              for ext in ("npz","json"))
    for name in exits
}
analysis_rc=int((folder/"analyzer.exitcode").read_text())
analysis=None
if analysis_rc==0:
    analysis=json.loads((folder/"sparse_matrix_analysis.json").read_text())
status=(
    "PASS_EXACT_FIRST_CN_JACOBIAN_CAPTURE_ANALYSIS_NOT_SOLVER_ACCEPTANCE"
    if captured=={"mumps":True,"superlu_dist":True}
       and exits["mumps"]==0 and analysis_rc==0
    else "FAIL_OR_INCOMPLETE_FIRST_CN_JACOBIAN_DIAGNOSTIC"
)
summary={
    "status":status, "native_exit_codes":exits,
    "all_exact_factor_matrices_saved":captured,
    "scipy_analysis_exit":analysis_rc,
    "scipy_analysis_status":analysis["status"] if analysis else "NOT_AVAILABLE",
    "SuperLU_DIST_backend_passed":False,
    "three_original_native_LU_passes_preserved":True,
    "full_five_backend_gate_still": "FAIL",
}
(folder/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
archive=Path(str(folder)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for path in sorted(folder.rglob("*")):
        if path.is_file():
            z.write(path,arcname=str(path.relative_to(folder)))
print(json.dumps({**summary,"archive":str(archive)},indent=2),flush=True)
if status!="PASS_EXACT_FIRST_CN_JACOBIAN_CAPTURE_ANALYSIS_NOT_SOLVER_ACCEPTANCE":
    raise SystemExit(1)
PY
echo "FIRST CN MATRIX ANALYSIS DONE — physical LU 3/5 unchanged; evidence $out"
