#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Single native original Heat-Wave strict-DAE CN run, inspect whether
# small signed interface power is weak coupling or cancellation.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2;exit 2; }
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
 --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
 echo "Unexpected GetFEM/PETSc image revision: $revision" >&2;exit 2
}
python3 -m py_compile sandbox/probe_v130_hw_interface_weak_coupling.py
python3 tests/test_v130_hw_interface_weak_coupling.py
python3 tests/test_v130_hw_interface_weak_coupling.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-hw-interface-audit-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-original-HW-interface-audit-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
mkdir "$out/physical" "$out/work"
echo "Original Heat-Wave weak interface observation, MUMPS/clamp MPI1"
echo "True GetFEM region=10 traces, individual port norms and signed power-density squared"
echo "Actual original sparse stiffness trace/feedback coupling blocks"
echo "No altered published sources, no forced excitation, no physical coupling PASS"
echo "EVIDENCE: $out"
if ! "${docker_cmd[@]}" run --rm \
 --user "$(id -u):$(id -g)" --env HOME=/tmp --env MPLBACKEND=Agg \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
 --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
 --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
 --env SCRIMP_MATRIX_MEMORY_MODE=compact \
 --env SCRIMP_MUMPS_RESET_POLICY=on_change \
 --env SCRIMP_GETFEM_RHS_MODE=auto \
 --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
 --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
 --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
 --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
 --volume "$root/sandbox/probe_v130_main_2d_examples.py:/workspace/scrimp/sandbox/probe_v130_main_2d_examples.py:ro" \
 --volume "$root/sandbox/measure_v130_heat_wave_physics.py:/workspace/scrimp/sandbox/measure_v130_heat_wave_physics.py:ro" \
 --volume "$root/sandbox/probe_v130_heat_wave_boundary_compatibility.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_boundary_compatibility.py:ro" \
 --volume "$root/sandbox/probe_v130_heat_wave_compatible_cn.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_compatible_cn.py:ro" \
 --volume "$root/sandbox/probe_v130_hw_interface_weak_coupling.py:/workspace/scrimp/sandbox/probe_v130_hw_interface_weak_coupling.py:ro" \
 --volume "$out/physical:/artifacts:rw" \
 --volume "$out/work:/scrimp-work:rw" --workdir /scrimp-work \
 "$image" timeout -k 5s "${SCRIMP_V130_HW_INTERFACE_AUDIT_TIMEOUT:-300s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /workspace/scrimp/sandbox/probe_v130_hw_interface_weak_coupling.py \
 --candidate clamp --out /artifacts \
 > "$out/native.host.log" 2>&1; then
 echo "FAIL original native weak interface audit; evidence retained: $out" >&2
 tail -n 175 "$out/native.host.log" >&2;exit 1
fi
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
report=json.loads((root/"physical/interface_weak_coupling_diagnostic.json").read_text())
cn=json.loads((root/"physical/compatible_cn_report.json").read_text())
assert report["status"]=="DIAGNOSTIC_ONLY_INTERFACE_COUPLING_PHYSICS_OPEN"
assert cn["full_original_CN_TS_executed"] is True
assert cn["physical_balance_gate_passed"] is True
assert cn["DAE_correction"]["DAE_initializer"]["status"]=="PASS"
assert report["matrix"]["full_state_DOFS"]==7689
assert report["raw_trace_mismatch_is_not_proof_of_weak_FEM_constraint_failure"]
assert not report["full_multiphysics_coupling_qualified"]
for name in ("U_T_CN","Y_w_CN","U_w_CN","Y_T_CN","T_CN","p_CN"):
 assert name in report["true_region10_squared_integrals"]
for block in ("p_to_Y_w_trace","e_Q_to_Y_T_flux",
              "Y_w_to_U_T_feedback","Y_T_to_U_w_feedback"):
 assert block in report["matrix"]["coupling_blocks"]
compact={
 "status":"PASS_DIAGNOSTIC_CAPTURED_NOT_COUPLING_QUALIFIED",
 "original_signed_CN_power_max":report["original_max_signed_CN_interface_power"],
 "original_CN_power_exchange_defect":report["original_max_CN_power_exchange_defect"],
 "surface_norms_CN":{
  k:v for k,v in report["true_region10_squared_integrals"].items()
  if k.endswith("_CN")
 },
 "coupling_blocks":report["matrix"]["coupling_blocks"],
 "optional_trace_errors":report["optional_GWFL_trace_errors"],
 "full_multiphysics_coupling_qualified":False,
}
(root/"comparison.json").write_text(json.dumps(compact,indent=2)+"\n",encoding="utf-8")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(root.rglob("*")):
  if f.is_file(): z.write(f,arcname=str(f.relative_to(root)))
print(json.dumps(compact,indent=2))
print("SHAREABLE ORIGINAL INTERFACE WEAK-COUPLING DIAGNOSTICS:",archive)
PY
echo "COMPLETE original Heat-Wave interface audit: $out"
