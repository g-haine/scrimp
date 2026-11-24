#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Stage 1 only: audit actual original Heat-Wave stiffness K, vs a
# sandbox-only dual-one-sided-face interface. Do NOT run PETSc TS/DAE.
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
 echo "Pinned original GetFEM/PETSc container required: $rev" >&2;exit 2
}
python3 -m py_compile sandbox/probe_v130_hw_dual_sided_assembly.py
python3 tests/test_v130_hw_dual_sided_assembly.py
python3 tests/test_v130_hw_dual_sided_assembly.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-hw-dual-interface-K-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-HW-dual-face-assembly-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$rev" > "$out/manifest.json"
echo "Original Heat-Wave K GetFEM/PETSc: two-sided original vs dual one-sided facets"
echo "Stage 1: actual mass/stiffness assembly ONLY; NO DAE, NO CN, NO PHYSICS PASS"
echo "Original source files pinned and immutable, MPI1"
echo "EVIDENCE: $out"
for variant in baseline dual_sided; do
 mkdir "$out/$variant" "$out/work_$variant"
 echo "ASSEMBLE variant=$variant"
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
  --volume "$root/sandbox/probe_v130_hw_one_sided_facets.py:/workspace/scrimp/sandbox/probe_v130_hw_one_sided_facets.py:ro" \
  --volume "$root/sandbox/probe_v130_hw_dual_sided_assembly.py:/workspace/scrimp/sandbox/probe_v130_hw_dual_sided_assembly.py:ro" \
  --volume "$out/$variant:/artifacts:rw" \
  --volume "$out/work_$variant:/scrimp-work:rw" --workdir /scrimp-work \
  "$image" timeout -k 5s "${SCRIMP_V130_HW_DUAL_ASSEMBLY_TIMEOUT:-240s}" \
  mpiexec --oversubscribe --bind-to none -n 1 \
  python -u /workspace/scrimp/sandbox/probe_v130_hw_dual_sided_assembly.py \
    --variant "$variant" --out /artifacts/stiffness_audit.json \
  > "$out/$variant.host.log" 2>&1; then
   echo "UNEXPECTED assembly diagnostic error in $variant: $out" >&2
   tail -n 130 "$out/$variant.host.log" >&2
   exit 1
 fi
done
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
records={
 k:json.loads((root/k/"stiffness_audit.json").read_text())
 for k in ("baseline","dual_sided")
}
a,b=records["baseline"],records["dual_sided"]
assert a["variant"]=="baseline" and b["variant"]=="dual_sided"
assert a["status"]==b["status"]=="DIAGNOSTIC_ASSEMBLY_ONLY_NO_DAE_NO_TS"
assert not a["time_integration_executed"] and not b["DAE_initialized"]
assert a["original_example_source_modified"] is False
assert b["original_example_source_modified"] is False
assert b["HEAT_on_disk_side_region30"] and b["WAVE_on_annulus_side_region31"]
keys=("heat_flux_into_interface_output","heat_boundary_input_into_flux",
      "wave_velocity_into_interface_output","wave_boundary_input_into_velocity",
      "feedback_wave_output_into_heat_input",
      "feedback_heat_output_into_wave_input")
blocks={}
for key in keys:
 before=float(a["blocks"][key]["Frobenius_norm"])
 after=float(b["blocks"][key]["Frobenius_norm"])
 blocks[key]={"original_region10":before,"dual_one_sided":after,
              "ratio":after/max(before,1e-30)}
heat_both=all(b["blocks"][k]["Frobenius_norm"]>1e-4 for k in keys[:2])
wave_both=all(b["blocks"][k]["Frobenius_norm"]>1e-4 for k in keys[2:4])
feedback_both=all(b["blocks"][k]["Frobenius_norm"]>1e-4 for k in keys[4:])
dofs_same=a["stiffness_matrix_full_DOFS"]==b["stiffness_matrix_full_DOFS"]
ready=heat_both and wave_both and feedback_both and dofs_same
result={
 "schema":"scrimp-v130-HW-original-vs-dual-interface-K-v1",
 "status":"PASS_EXPECTED_SPARSE_INTERFACE_BLOCKS_RESTORED_ASSEMBLY_ONLY"
     if ready else "ASSEMBLY_ONLY_INTERFACE_BLOCKS_STILL_UNQUALIFIED",
 "actual_baseline_dofs":a["stiffness_matrix_full_DOFS"],
 "actual_dual_sided_dofs":b["stiffness_matrix_full_DOFS"],
 "same_original_dynamic_and_algebraic_dimension":dofs_same,
 "coupling_blocks":blocks,
 "nonzero_heat_normal_blocks":heat_both,
 "nonzero_wave_trace_blocks":wave_both,
 "nonzero_bidirectional_feedback_blocks":feedback_both,
 "no_original_weak_form_text_changes":True,
 "DAE_initialization_executed":False,
 "CN_trajectory_executed":False,
 "energy_or_physical_coupling_qualified":False,
}
(root/"comparison.json").write_text(json.dumps(result,indent=2)+"\n")
with zipfile.ZipFile(str(root)+".zip","x",compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(root.rglob("*")):
  if f.is_file():z.write(f,arcname=str(f.relative_to(root)))
print(json.dumps(result,indent=2))
print("SHAREABLE ORIGINAL HEAT-WAVE DUAL-INTERFACE ASSEMBLY EVIDENCE:",str(root)+".zip")
PY
echo "COMPLETE original vs dual-sided assembly comparison: $out"
