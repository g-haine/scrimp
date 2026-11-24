#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# One native unchanged original Heat-Wave, then postsolve-only
# GetFEM one-sided 2D-face normal diagnostics.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
rev="$("${docker_cmd[@]}" image inspect \
 --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$rev" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
 echo "Unqualified original PETSc/GetFEM image: $rev" >&2; exit 2;
}
python3 -m py_compile sandbox/probe_v130_hw_one_sided_facets.py
python3 tests/test_v130_hw_one_sided_facets.py
python3 tests/test_v130_hw_one_sided_facets.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-hw-one-sided-facet-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-HW-one-sided-facet-diagnostic-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$rev" > "$out/manifest.json"
mkdir "$out/physical" "$out/work"
echo "Original Heat-Wave native MUMPS/clamp 4-step CN, strict DAE, MPI1"
echo "Probe GetFEM region 10 versus heat 2D faces 30 and wave 2D faces 31"
echo "Additional facet regions created ONLY AFTER original solve and physical collector"
echo "Read-only original source; no changed initial state or matrix"
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
 --volume "$root/sandbox/probe_v130_hw_one_sided_facets.py:/workspace/scrimp/sandbox/probe_v130_hw_one_sided_facets.py:ro" \
 --volume "$out/physical:/artifacts:rw" \
 --volume "$out/work:/scrimp-work:rw" --workdir /scrimp-work \
 "$image" timeout -k 5s "${SCRIMP_V130_HW_FACET_TIMEOUT:-300s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /workspace/scrimp/sandbox/probe_v130_hw_one_sided_facets.py \
 --out /artifacts \
 > "$out/native.host.log" 2>&1; then
 echo "FAIL native one-sided facet diagnosis, evidence retained: $out" >&2
 tail -n 140 "$out/native.host.log" >&2;exit 1
fi
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
diag=json.loads((root/"physical/one_sided_facet_normal_diagnostic.json").read_text())
physics=json.loads((root/"physical/compatible_cn_report.json").read_text())
assert diag["original_source_and_weak_forms_unmodified"]
assert diag["diagnostic_regions_added_after_actual_TS_solve"]
assert physics["full_original_CN_TS_executed"] is True
assert physics["DAE_correction"]["DAE_initializer"]["status"]=="PASS"
assert physics["physical_balance_gate_passed"] is True
assert diag["main_real_CN_steps"]==4
for reg in ("original_region10","heat_faces_region30","wave_faces_region31"):
 assert reg in diag["integrals_CN"]
compact={
 "schema":"scrimp-v130-HW-one-sided-facet-comparison-v1",
 "status":"PASS_FACET_DIAGNOSTIC_CAPTURED_PHYSICS_REQUIRES_FOLLOWUP",
 "topology":diag["topology"],
 "normals":{
  r:{k:v for k,v in data.items() if
     k.endswith("_sq") or "power" in k}
  for r,data in diag["integrals_CN"].items()
 },
 "unavailable_traces":diag["unavailable_cross_FEM_traces"],
 "original_strict_DAE_and_CN_pass":True,
 "physical_interface_coupling_qualified":False
}
(root/"comparison.json").write_text(json.dumps(compact,indent=2)+"\n",encoding="utf-8")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(root.rglob("*")):
  if f.is_file():z.write(f,arcname=str(f.relative_to(root)))
print(json.dumps(compact,indent=2))
print("SHAREABLE ORIGINAL HEAT-WAVE ONE-SIDED FACET EVIDENCE:",archive)
PY
echo "COMPLETE native original facet-normal diagnosis: $out"
