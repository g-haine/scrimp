#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Strict native DAE feasibility of separately oriented original Heat-Wave
# interface (one baseline + two candidate p0 treatments). NO CN/TS.
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
 echo "Pinned image revision required: $revision" >&2;exit 2
}
python3 -m py_compile \
 sandbox/probe_v130_hw_dual_sided_dae.py \
 sandbox/v130_hw_interface_geometry.py
python3 tests/test_v130_hw_dual_sided_dae.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-hw-dual-interface-DAE-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-HW-dual-interface-strict-DAE-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Native PETSc strict DAE on one original and two corrected GetFEM systems"
echo "Baseline clamp; dual-sided clamp; dual-sided constrained-L2"
echo "Every run independent MPI1; source original unchanged; no CN/TS"
echo "Incompatibility is evidence, not a reason to relax DAE or change physical states"
echo "EVIDENCE: $out"
for spec in baseline:clamp dual_sided:clamp dual_sided:constrained_l2; do
 IFS=: read -r variant candidate <<<"$spec"
 case_name="${variant}_${candidate}"
 mkdir "$out/$case_name" "$out/work_$case_name"
 echo "RUN original assembly=$variant DAE-candidate=$candidate"
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
  --volume "$root/sandbox/v130_hw_interface_geometry.py:/workspace/scrimp/sandbox/v130_hw_interface_geometry.py:ro" \
  --volume "$root/sandbox/probe_v130_heat_wave_boundary_compatibility.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_boundary_compatibility.py:ro" \
  --volume "$root/sandbox/probe_v130_hw_dual_sided_dae.py:/workspace/scrimp/sandbox/probe_v130_hw_dual_sided_dae.py:ro" \
  --volume "$out/$case_name:/artifacts:rw" \
  --volume "$out/work_$case_name:/scrimp-work:rw" --workdir /scrimp-work \
  "$image" timeout -k 5s "${SCRIMP_V130_HW_DUAL_DAE_TIMEOUT:-240s}" \
  mpiexec --oversubscribe --bind-to none -n 1 \
  python -u /workspace/scrimp/sandbox/probe_v130_hw_dual_sided_dae.py \
  --variant "$variant" --candidate "$candidate" \
  --out /artifacts/strict_DAE.json \
  > "$out/$case_name.host.log" 2>&1 || rc=$?
 printf '%s\n' "$rc" > "$out/$case_name.exitcode"
 if [[ "$rc" -ne 0 ]]; then
  echo "UNEXPECTED infrastructure/native failure $case_name exit=$rc" >&2
  tail -n 90 "$out/$case_name.host.log" >&2
 fi
done
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
runs={}
unexpected={}
for variant,candidate in (
 ("baseline","clamp"),
 ("dual_sided","clamp"),
 ("dual_sided","constrained_l2"),
):
 name=variant+"_"+candidate
 code=int((root/f"{name}.exitcode").read_text())
 nativefile=root/name/"strict_DAE.json"
 if code or not nativefile.is_file():
  unexpected[name]={"exit_code":code,
   "last_log_lines":(root/f"{name}.host.log").read_text(
    errors="replace"
   ).splitlines()[-16:]}
  runs[name]={"status":"UNEXPECTED_NATIVE_PROCESS_FAILURE"}
  continue
 report=json.loads(nativefile.read_text())
 assert report["variant"]==variant and report["candidate"]==candidate
 assert report["assembly"]["DOFs"]==7689
 assert report["strict_initializer_called"] is True
 assert report["no_pseudo_TS_or_tolerance_relaxation"] is True
 assert report["main_CN_TS_executed"] is False
 assert report["physical_power_balance_qualified"] is False
 native=report["native_initialization"]
 assert native["differential_initial_exact"] is True
 assert native["pseudo_ts_executed"] is False
 assert report["p_boundary_after_l2"]==0.
 status=report["status"]
 if status not in (
  "PASS_NATIVE_STRICT_DAE_FEASIBILITY_ONLY",
  "INCOMPATIBLE_INITIAL_DAE_UNALTERED_PHYSICS",
 ):
  raise RuntimeError("Missing explicit native DAE verdict: "+name)
 runs[name]={
  "status":status,
  "native_algebraic_initial_l2":native["native_initial_residual_l2"],
  "native_algebraic_final_l2":native["native_final_residual_l2"],
  "native_feasibility_tolerance":native["feasibility_tolerance"],
  "accuracy_target_met":native["accuracy_target_met"],
  "dominant_incompatible_variables":
    native.get("dominant_algebraic_variables"),
  "assembly_nnz":report["assembly"]["stored_nnz"],
  "heat_flux_interface_block":
    report["assembly"]["blocks"]["heat_normal_flux_trace"]["Frobenius_norm"],
  "wave_trace_block":
    report["assembly"]["blocks"]["wave_velocity_trace"]["Frobenius_norm"],
 }
baseline_ok=runs["baseline_clamp"]["status"]=="PASS_NATIVE_STRICT_DAE_FEASIBILITY_ONLY"
dual_ok=any(
 runs[name]["status"]=="PASS_NATIVE_STRICT_DAE_FEASIBILITY_ONLY"
 for name in ("dual_sided_clamp","dual_sided_constrained_l2")
)
summary={
 "schema":"scrimp-v130-HW-original-dual-sided-initial-DAE-v1",
 "status":(
    "DIAGNOSTIC_INCOMPLETE_UNEXPECTED_NATIVE_ERROR" if unexpected
    else "FAIL_BASELINE_STRICT_DAE_REGRESSION" if not baseline_ok
    else "PASS_DUAL_FACET_STRICT_DAE_INITIAL_FEASIBILITY_ONLY" if dual_ok
    else "DUAL_FACET_DAE_INITIAL_INCOMPATIBILITY_DIAGNOSED"
 ),
 "baseline_clamp_is_strict_DAE_feasible":baseline_ok,
 "dual_sided_any_strict_DAE_feasible":dual_ok,
 "cases":runs,
 "unexpected_native_failures":unexpected,
 "original_public_example_bytes_unchanged":True,
 "main_TS_or_CN_tested":False,
 "physical_energy_balance_qualified":False,
 "production_solver_changed":False,
}
(root/"comparison.json").write_text(
 json.dumps(summary,indent=2)+"\n",encoding="utf-8"
)
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(root.rglob("*")):
  if f.is_file():z.write(f,arcname=str(f.relative_to(root)))
print(json.dumps(summary,indent=2))
print("SHAREABLE ORIGINAL DUAL-SIDED STRICT DAE EVIDENCE:",archive)
if unexpected or not baseline_ok:
 raise SystemExit("Unexpected native failure or control DAE regression; evidence archived")
PY
echo "COMPLETE native original vs dual-sided strict DAE feasibility: $out"
