#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
rev="$("${docker_cmd[@]}" image inspect --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$rev" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
 echo "Unexpected image revision $rev" >&2; exit 2;
}
python3 tests/test_v130_heat_wave_superlu_step_diagnostic.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-hw-superlu-diagnostic-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"source":"%s","image_revision":"%s","schema":"scrimp-v130-superlu-step-diagnostic-v1"}\n' "$(git rev-parse HEAD)" "$rev" > "$out/manifest.json"
echo "Native SuperLU_DIST first-step SNES/KSP/PC diagnostics, MPI1"
echo "No relaxation, no changes to published Heat-Wave source"
echo "EVIDENCE: $out"
for candidate in clamp constrained_l2; do
 mkdir "$out/$candidate" "$out/work_$candidate"
 echo "RUN SuperLU_DIST $candidate"
 if ! "${docker_cmd[@]}" run --rm --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env MPLBACKEND=Agg \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
 --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
 --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
 --env SCRIMP_MATRIX_MEMORY_MODE=compact \
 --env SCRIMP_MUMPS_RESET_POLICY=on_change \
 --env SCRIMP_GETFEM_RHS_MODE=auto \
 --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
 --env SCRIMP_TS_FAILURE_DIAGNOSTICS=1 \
 --env SCRIMP_HW_TS_DIAGNOSTIC_MONITORS=1 \
 --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
 --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
 --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
 --volume "$root/sandbox/probe_v130_main_2d_examples.py:/workspace/scrimp/sandbox/probe_v130_main_2d_examples.py:ro" \
 --volume "$root/sandbox/measure_v130_heat_wave_physics.py:/workspace/scrimp/sandbox/measure_v130_heat_wave_physics.py:ro" \
 --volume "$root/sandbox/probe_v130_heat_wave_boundary_compatibility.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_boundary_compatibility.py:ro" \
 --volume "$root/sandbox/probe_v130_heat_wave_compatible_cn.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_compatible_cn.py:ro" \
 --volume "$out/$candidate:/artifacts:rw" \
 --volume "$out/work_$candidate:/scrimp-work:rw" --workdir /scrimp-work \
 "$image" timeout -k 5s "${SCRIMP_V130_HW_SUPERLU_DIAG_TIMEOUT:-300s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /workspace/scrimp/sandbox/probe_v130_heat_wave_compatible_cn.py \
 --candidate "$candidate" --backend superlu_dist \
 --capture-ts-failure --out /artifacts \
 > "$out/$candidate.host.log" 2>&1; then
  echo "UNEXPECTED native error for $candidate: $out" >&2
  tail -n 150 "$out/$candidate.host.log" >&2; exit 1
 fi
done
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
cases={}
for c in ("clamp","constrained_l2"):
 path=root/c/"compatible_cn_failure.json"
 completed=root/c/"compatible_cn_report.json"
 log=(root/f"{c}.host.log").read_text(errors="replace")
 if path.is_file():
  r=json.loads(path.read_text())
  assert r["status"]=="DIAGNOSTIC_TS_REJECTED_NOT_PHYSICAL_PASS"
  assert r["strict_DAE_pass"] and not r["physics_qualified"]
  assert r["initial_algebraic_residual"]<r["initial_algebraic_tolerance"]
  assert "MAIN TS FAILURE DIAGNOSTIC:" in log
  cases[c]={"status":r["status"],"PETSc_error":r["PETSc_error"],
            "DAE_residual":r["initial_algebraic_residual"]}
 elif completed.is_file():
  r=json.loads(completed.read_text())
  assert r["original_numerical_report"]["steps"]==4
  cases[c]={"status":r["status"],
            "physics_gate_passed":r["physical_balance_gate_passed"]}
 else: raise RuntimeError(f"Missing actual PETSc result for {c}")
 interesting=[line for line in log.splitlines()
             if any(s in line for s in (
               "MAIN TS FAILURE DIAGNOSTIC:",
               "DIVERGED_", "converged reason", "SNES Function norm",
               "KSP reason", "snes_monitor"
             ))]
 print(f"--- {c} PETSc solver diagnostics ---")
 print("\n".join(interesting[-25:]) or "(see native host log)")
status=("DIAGNOSTIC_ONLY_TS_FAILURE_REPRODUCED" if any(
 x["status"]=="DIAGNOSTIC_TS_REJECTED_NOT_PHYSICAL_PASS"
 for x in cases.values()) else "NO_TS_FAILURE_IN_DIAGNOSTIC_RUN")
summary={"schema":"scrimp-v130-superlu-step-diagnostic-v1",
         "status":status,"cases":cases,
         "physics_qualified":False,
         "original_model_modified":False,
         "solver_tolerances_modified":False}
(root/"comparison.json").write_text(json.dumps(summary,indent=2)+"\n")
with zipfile.ZipFile(str(root)+".zip","x",compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(root.rglob("*")):
  if f.is_file():z.write(f,arcname=str(f.relative_to(root)))
print(json.dumps(summary,indent=2))
print("ORIGINAL SHAREABLE SUPERLU STEP DIAGNOSTICS:",str(root)+".zip")
PY
echo "COMPLETE native SuperLU_DIST diagnostic: $out"
