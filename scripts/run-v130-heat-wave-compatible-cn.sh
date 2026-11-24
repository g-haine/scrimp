#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Native original Heat-Wave 4-step CN physical boundary candidates.
# Both explicit corrections, MUMPS vs SuperLU_DIST, serial ParaView.
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
 echo "Pinned GetFEM/PETSc revision required: $revision" >&2; exit 2;
}
python3 -m py_compile sandbox/probe_v130_heat_wave_compatible_cn.py
python3 tests/test_v130_heat_wave_compatible_cn.py
python3 tests/test_v130_heat_wave_compatible_cn.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-hw-compatible-cn-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-heat-wave-compatible-CN-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Native original Heat-Wave CN with corrected initial p exterior trace"
echo "Clamp and constrained L2, each with MUMPS/SuperLU_DIST, 4 CN steps, MPI1"
echo "Full Hamiltonian partition, heat dissipation, signed interface/outer powers"
echo "Source examples read-only, 8 immutable originals, no automatic DAE relaxation"
echo "Separate physics gate; numerical LU parity is NOT physical acceptance"
echo "EVIDENCE: $out"
for candidate in clamp constrained_l2; do
 for backend in mumps superlu_dist; do
  name="${candidate}_${backend}"
  folder="$out/$name"
  mkdir "$folder" "$out/work_$name"
  echo "RUN original Heat-Wave compatible CN candidate=$candidate LU=$backend"
  if ! "${docker_cmd[@]}" run --rm \
   --user "$(id -u):$(id -g)" --env HOME=/tmp --env MPLBACKEND=Agg \
   --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
   --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
   --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
   --env SCRIMP_MATRIX_MEMORY_MODE=compact \
   --env SCRIMP_MUMPS_RESET_POLICY=on_change \
   --env SCRIMP_GETFEM_RHS_MODE=auto \
   --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
   --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
   --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
   --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
   --volume "$root/sandbox/probe_v130_main_2d_examples.py:/workspace/scrimp/sandbox/probe_v130_main_2d_examples.py:ro" \
   --volume "$root/sandbox/measure_v130_heat_wave_physics.py:/workspace/scrimp/sandbox/measure_v130_heat_wave_physics.py:ro" \
   --volume "$root/sandbox/probe_v130_heat_wave_boundary_compatibility.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_boundary_compatibility.py:ro" \
   --volume "$root/sandbox/probe_v130_heat_wave_compatible_cn.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_compatible_cn.py:ro" \
   --volume "$folder:/artifacts:rw" \
   --volume "$out/work_$name:/scrimp-work:rw" --workdir /scrimp-work \
   "$image" timeout -k 5s "${SCRIMP_V130_HW_COMPATIBLE_CN_TIMEOUT:-300s}" \
   mpiexec --oversubscribe --bind-to none -n 1 \
   python -u /workspace/scrimp/sandbox/probe_v130_heat_wave_compatible_cn.py \
    --candidate "$candidate" --backend "$backend" --out /artifacts \
   > "$out/$name.host.log" 2>&1; then
   echo "FAIL original Heat-Wave physical CN candidate=$candidate solver=$backend" >&2
   echo "Original evidence retained at: $out" >&2
   tail -n 180 "$out/$name.host.log" >&2
   exit 1
  fi
 done
done
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
import numpy as np
root=Path(sys.argv[1])
summary={}
physics_gates={}
for candidate in ("clamp","constrained_l2"):
 reports={}
 for backend in ("mumps","superlu_dist"):
  folder=root/f"{candidate}_{backend}"
  r=json.loads((folder/"compatible_cn_report.json").read_text())
  assert r["candidate"]==candidate and r["backend"]==backend
  assert r["full_original_CN_TS_executed"] is True
  assert r["DAE_correction"]["DAE_initializer"]["status"]=="PASS"
  assert r["DAE_correction"]["T_q_unchanged"] is True
  assert r["time_steps"]==4
  assert r["original_example_source_modified"] is False
  with np.load(folder/"original_2d_example_trace.npz",allow_pickle=False) as a:
   traces={key:a[key] for key in a.files}
  assert traces["times"].shape==(5,)
  assert traces["states"].shape==(5,7689)
  assert np.allclose(traces["times"],np.arange(5)*0.001,atol=1e-11,rtol=0.)
  assert np.all(np.isfinite(traces["states"]))
  assert "phys_CN_per_step_H_balance_error" in traces
  assert "phys_CN_P_interface_physical" in traces
  assert "phys_CN_D_heat_flux" in traces
  reports[backend]=(r,traces)
 left,right=(reports["mumps"][1],reports["superlu_dist"][1])
 assert set(left)==set(right),(candidate,"different fields")
 differences={}
 for key in sorted(left):
  x,y=left[key],right[key]
  assert x.shape==y.shape,(candidate,key,x.shape,y.shape)
  if key=="states":
   err=float(np.max(np.linalg.norm(x-y,axis=1)))
  else:
   err=float(np.max(np.abs(x-y)))
  assert np.isfinite(err) and err<1e-8,(candidate,key,err)
  differences[key]=err
 passed=all(r[0]["physical_balance_gate_passed"] for r in reports.values())
 physics_gates[candidate]=passed
 summary[candidate]={
  "status":"PASS_NUMERICAL_CN_BOTH_LU_PHYSICAL_GATES_"+(
    "PASS" if passed else "OPEN"),
  "max_MUMPS_vs_SuperLU_dist_observable_differences":differences,
  "MUMPS_physical_report":reports["mumps"][0]["original_physics_report"],
  "initial_DAE_final_residual":reports["mumps"][0][
   "DAE_correction"]["DAE_initializer"]["native_final_residual_l2"],
  "physical_balance_gate_passed_both_LU":passed,
 }
status=("PASS_COMPATIBLE_CN_PARITY_AND_SCIENTIFIC_THRESHOLDS"
        if all(physics_gates.values()) else
        "PASS_CN_BACKEND_PARITY_PHYSICAL_SCIENCE_OPEN")
result={
 "schema":"scrimp-v130-heat-wave-compatible-CN-MPI1-v1",
 "status":status,"candidates":summary,
 "published_examples_modified":False,"MPI1_only":True,
 "parallel_export_qualified":False,
 "heat_wave_physical_acceptance":all(physics_gates.values()),
}
(root/"comparison.json").write_text(json.dumps(result,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
 for path in sorted(root.rglob("*")):
  if path.is_file():z.write(path,arcname=str(path.relative_to(root)))
print(json.dumps(result,indent=2))
print("SHAREABLE ORIGINAL COMPATIBLE CN HEAT-WAVE EVIDENCE:",archive)
PY
echo "COMPLETE original Heat-Wave corrected-initial-data CN audit: $out"
