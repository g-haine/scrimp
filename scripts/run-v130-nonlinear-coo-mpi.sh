#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Native GetFEM/PETSc COO dynamic-pattern gate and original nonlinear
# shallow-water mass/stiffness IJacobian FD gate (MPI1/2).
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
 echo "Unqualified original GetFEM/PETSc image $revision" >&2; exit 2;
}
python3 -m py_compile sandbox/probe_v130_nonlinear_coo_dynamic_graph.py
python3 tests/test_v130_nonlinear_coo_assembly.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-nonlinear-coo-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-nonlinear-coo-gate-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Original nonlinear GetFEM assembly -> PETSc COO -> AIJ, MPI1/2"
echo "Dynamic sparse graph re-preallocation plus actual shallow-water IJacobian FD"
echo "No article-era sources changed; no main TS/physical trajectory"
echo "EVIDENCE: $out"
for ranks in 1 2; do
 mkdir "$out/${ranks}r" "$out/work_${ranks}"
 echo "RUN nonlinear dynamic COO graph MPI=$ranks"
 if ! "${docker_cmd[@]}" run --rm \
   --user "$(id -u):$(id -g)" --env HOME=/tmp \
   --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
   --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
   --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
   --env SCRIMP_MPI_TRACE=1 \
   --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
   --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
   --volume "$root/sandbox/probe_v130_nonlinear_coo_dynamic_graph.py:/workspace/scrimp/sandbox/probe_v130_nonlinear_coo_dynamic_graph.py:ro" \
   --volume "$out/${ranks}r:/artifacts:rw" \
   --workdir /artifacts \
   "$image" timeout -k 5s "${SCRIMP_V130_NL_COO_TIMEOUT:-240s}" \
   mpiexec --oversubscribe --bind-to none -n "$ranks" \
   python -u /workspace/scrimp/sandbox/probe_v130_nonlinear_coo_dynamic_graph.py \
   --out /artifacts/dynamic.json \
   > "$out/${ranks}r/dynamic.host.log" 2>&1; then
  echo "FAIL dynamic COO graph MPI=$ranks; evidence $out" >&2
  tail -n 150 "$out/${ranks}r/dynamic.host.log" >&2; exit 1
 fi
 echo "RUN original nonlinear shallow-water FD IJacobian MPI=$ranks"
 if ! "${docker_cmd[@]}" run --rm \
   --user "$(id -u):$(id -g)" --env HOME=/tmp \
   --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
   --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
   --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
   --env SCRIMP_MATRIX_MEMORY_MODE=compact \
   --env SCRIMP_TRAJECTORY_BACKEND=memory \
   --env SCRIMP_MUMPS_RESET_POLICY=on_change \
   --env SCRIMP_GETFEM_RHS_MODE=auto \
   --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
   --env SCRIMP_MPI_TRACE=1 \
   --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
   --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
   --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
   --volume "$root/sandbox/check_v130_shallow_water_mass_jacobian.py:/workspace/scrimp/sandbox/check_v130_shallow_water_mass_jacobian.py:ro" \
   --volume "$out/${ranks}r:/artifacts:rw" \
   --volume "$out/work_${ranks}:/scrimp-work:rw" \
   --workdir /scrimp-work \
   "$image" timeout -k 5s "${SCRIMP_V130_NL_COO_TIMEOUT:-240s}" \
   mpiexec --oversubscribe --bind-to none -n "$ranks" \
   python -u /workspace/scrimp/sandbox/check_v130_shallow_water_mass_jacobian.py \
     --mesh .2 --epsilon 1e-6 --state-perturb .01 \
     --mode default --output-dir /artifacts \
   > "$out/${ranks}r/shallow.host.log" 2>&1; then
  echo "FAIL original shallow-water nonlinear Jacobian MPI=$ranks; evidence $out" >&2
  tail -n 150 "$out/${ranks}r/shallow.host.log" >&2; exit 1
 fi
 echo "PASS original nonlinear SWE directional test, MPI=$ranks"
done
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
import numpy as np
root=Path(sys.argv[1])
summary={}
global_states={}
for ranks in (1,2):
 folder=root/f"{ranks}r"
 # Rank-local graph source must be a full two- or one-file set.
 records=[json.loads((folder/f"dynamic_rank_{r}.json").read_text())
          for r in range(ranks)]
 assert all(r["status"]=="PASS_NATIVE_DYNAMIC_COO_REPREALLOCATION"
            for r in records)
 assert all(r["same_petsc_mat_reused"] for r in records)
 snapshots=[]
 for step in range(4):
  dense=np.zeros((records[0]["n"],records[0]["n"]),dtype=float)
  rows_owned=0
  for r in records:
   row=r["steps"][step]
   lo,hi=row["ownership"]
   dense[lo:hi]=np.asarray(row["owned_dense_rows"],dtype=float)
   rows_owned+=hi-lo
  assert rows_owned==dense.shape[0],(ranks,step,rows_owned)
  snapshots.append(dense)
 assert np.allclose(snapshots[0],snapshots[2],rtol=1e-12,atol=1e-12)
 assert np.allclose(snapshots[1],snapshots[3],rtol=1e-12,atol=1e-12)
 assert np.linalg.norm(snapshots[0]-snapshots[1])>1e-9
 global_states[ranks]=snapshots
 shallow=json.loads((folder/f"shallow_default_{ranks}r.json").read_text())
 errs=[float(a["relative_error"]) for a in shallow["outcomes"]]
 assert len(errs)==4 and max(errs)<1e-5,(ranks,errs)
 assert not shallow["main_TS_executed"]
 assert shallow["frozen_explicit_RHS_drift"]==0.
 for key in ("dynamic.host.log","shallow.host.log"):
  log=(folder/key).read_text(errors="replace")
  assert "GetFEM -> PETSc effective=coo" in log,(ranks,key)
  assert "GetFEM -> PETSc effective=ijv" not in log,(ranks,key)
 summary[ranks]={
  "dynamic_pattern_graph_counts": [
   records[0]["steps"][i]["global_stored_nnz"] for i in range(4)],
  "shallow_max_directional_Jacobian_relative_error":max(errs),
  "effective_coo_traced_on_all_GetFEM_bridges":True,
 }
for step in range(4):
 rel=float(np.linalg.norm(global_states[2][step]-global_states[1][step])/
     max(np.linalg.norm(global_states[1][step]),1e-30))
 assert rel<1e-10,(step,rel)
summary["status"]="PASS_NONLINEAR_COO_DYNAMIC_GRAPH_AND_SWE_JACOBIAN_MPI1_2"
summary["DAE_or_time_trajectory_qualified"]=False
(root/"comparison.json").write_text(json.dumps(summary,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as z:
 for file in sorted(root.rglob("*")):
  if file.is_file(): z.write(file,arcname=str(file.relative_to(root)))
print(json.dumps(summary,indent=2))
print("SHAREABLE NONLINEAR COO MPI EVIDENCE:",archive)
PY
echo "COMPLETE nonlinear COO GetFEM/PETSc MPI1/2 audit: $out"
