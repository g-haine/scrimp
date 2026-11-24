#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Real MPI1 PETSc 4-step Crank-Nicolson of ORIGINAL Heat-Wave
# with correct one-sided port integration, strict DAE, MUMPS.
# Physical balance and actual interface excitation are SEPARATE gates.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2;exit 2; }
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker);fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
 --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
 echo "Pinned image revision required: $revision" >&2;exit 2
}
python3 -m py_compile \
 sandbox/measure_v130_hw_dual_facet_physics.py \
 sandbox/probe_v130_hw_dual_sided_cn.py \
 sandbox/v130_hw_interface_geometry.py
python3 tests/test_v130_hw_dual_sided_cn.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-hw-dual-facet-CN-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-HW-original-corrected-CN-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Genuine original Heat-Wave 4-step CN with corrected heat/wave oriented faces"
echo "Strict PETSc DAE, original Gaussian projected L2, clamp and constrained-L2 p0"
echo "Real MUMPS factorization and all original postprocess files"
echo "CRUCIAL: physical port powers use region 30 and 31, never doubled 10"
echo "Separate gates for CN convergence, energy identities, actual port excitation"
echo "No production policy changes; original example source files read-only"
echo "EVIDENCE: $out"

for candidate in clamp constrained_l2; do
 mkdir "$out/$candidate" "$out/work_$candidate"
 echo "RUN real original four-step CN corrected faces candidate=$candidate"
 rc=0
 "${docker_cmd[@]}" run --rm --user "$(id -u):$(id -g)" \
  --env HOME=/tmp --env MPLBACKEND=Agg \
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
  --volume "$root/sandbox/probe_v130_heat_wave_boundary_compatibility.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_boundary_compatibility.py:ro" \
  --volume "$root/sandbox/probe_v130_heat_wave_compatible_cn.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_compatible_cn.py:ro" \
  --volume "$root/sandbox/measure_v130_heat_wave_physics.py:/workspace/scrimp/sandbox/measure_v130_heat_wave_physics.py:ro" \
  --volume "$root/sandbox/measure_v130_hw_dual_facet_physics.py:/workspace/scrimp/sandbox/measure_v130_hw_dual_facet_physics.py:ro" \
  --volume "$root/sandbox/v130_hw_interface_geometry.py:/workspace/scrimp/sandbox/v130_hw_interface_geometry.py:ro" \
  --volume "$root/sandbox/probe_v130_hw_dual_sided_cn.py:/workspace/scrimp/sandbox/probe_v130_hw_dual_sided_cn.py:ro" \
  --volume "$out/$candidate:/artifacts:rw" \
  --volume "$out/work_$candidate:/scrimp-work:rw" --workdir /scrimp-work \
  "$image" timeout -k 5s "${SCRIMP_V130_HW_DUAL_CN_TIMEOUT:-300s}" \
  mpiexec --oversubscribe --bind-to none -n 1 \
  python -u /workspace/scrimp/sandbox/probe_v130_hw_dual_sided_cn.py \
  --candidate "$candidate" --out /artifacts \
  > "$out/$candidate.host.log" 2>&1 || rc=$?
 printf '%s\n' "$rc" > "$out/$candidate.exitcode"
 if [[ "$rc" -ne 0 ]]; then
   echo "FAIL native corrected CN $candidate exit=$rc; preserving logs" >&2
   tail -n 110 "$out/$candidate.host.log" >&2
 fi
done

python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
cases={}
failures={}
for candidate in ("clamp","constrained_l2"):
 code=int((root/f"{candidate}.exitcode").read_text())
 path=root/candidate/"dual_facet_cn_audit.json"
 if code or not path.exists():
  cases[candidate]={"status":"UNEXPECTED_NATIVE_TS_FAILURE_NOT_PHYSICAL_RESULT",
    "exit_code":code}
  failures[candidate]=(root/f"{candidate}.host.log").read_text(
    errors="replace"
  ).splitlines()[-35:]
  continue
 report=json.loads(path.read_text())
 assert report["actual_real_4_step_PETSc_CN_PASSED"] is True
 assert report["strict_DAE_initialization_PASSED"] is True
 assert report["interface_heat_region"]==30
 assert report["interface_wave_region"]==31
 assert not report["example_src_modified"]
 assert not report["MPI2_or_parallel_export_qualified"]
 assert not report["production_solver_policy_changed"]
 physics=report["actual_energy_power_diagnostics"]
 assert physics["old_doubled_interface_region_excluded_from_physical_power"]
 cases[candidate]={
  "status":report["status"],
  "max_CN_total_balance_defect":
    physics["max_CN_step_H_balance_defect"],
  "max_CN_heat_balance_defect":
    physics["max_CN_step_heat_balance_defect"],
  "max_CN_wave_balance_defect":
    physics["max_CN_step_wave_balance_defect"],
  "max_CN_interface_power_defect":
    physics["max_CN_interface_power_defect"],
  "max_individual_interface_power":
    physics["max_individual_CN_interface_power"],
  "interface_excitation_relative_to_heat_dissipation":
    physics["interface_exchange_ratio_to_heat_dissipation"],
  "physical_balance_gate_passed":
    report["CN_discrete_midpoint_energy_gate_PASSED"],
  "interface_exchange_nontrivial":
    report["CN_interface_exchange_nontrivial"],
 }
summary={
 "schema":"scrimp-v130-HW-original-dual-facet-real-CN-comparison-v1",
 "status":(
    "DIAGNOSTIC_NATIVE_FAILURE_OR_TS_INCOMPLETE" if failures
    else "BOTH_DUAL_FACET_CN_CONVERGED_MODEL_PHYSICS_REQUIRES_REVIEW"
 ),
 "cases":cases,
 "native_failures":failures,
 "original_public_example_sources_unchanged":True,
 "actual_corrected_one_sided_interfaces_used":True,
 "actual_main_TS_full_CN_qualified_only_where_recorded":True,
 "MPI2_qualified":False,
 "SuperLU_DIST_qualified_for_corrected_Jacobian":False,
}
(root/"comparison.json").write_text(
 json.dumps(summary,indent=2)+"\n",encoding="utf-8"
)
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
 for f in sorted(root.rglob("*")):
  if f.is_file():
   z.write(f,arcname=str(f.relative_to(root)))
print(json.dumps(summary,indent=2))
print("SHAREABLE NATIVE PHYSICALLY CORRECTED HEAT-WAVE CN EVIDENCE:",archive)
if failures:raise SystemExit(
  "One or more actual PETSc corrected CN trajectories failed; see archive"
)
PY
echo "COMPLETE corrected original Heat-Wave CN MPI1 audit: $out"
