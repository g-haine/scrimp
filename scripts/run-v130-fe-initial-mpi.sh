#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Native MPI1/2 universal sparse FE projection and lumped numeric gate.
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
 echo "Qualified GetFEM/PETSc image required: $revision" >&2
 exit 2
}
python3 -m py_compile sandbox/probe_v130_fe_initial_mpi.py
python3 tests/test_v130_fe_initial_projection.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-fe-default-mpi-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-default-fe-l2-mpi-gate-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Original GetFEM MPI level 2, PETSc sparse MUMPS, DEFAULT L2 projection"
echo "GetFEM -> PETSc MATRIX MODE: distributed; reference IJV/ADD_VALUES bridge"
echo "CG2 scalar, CG1 vector, RT0, BDM1, Nedelec, Argyris, Gaussian + lumped"
echo "MPI 1/2, rank agreement, cross-MPI global coefficient comparison"
echo "Actual DAE feasibility, energy/power and PDE convergence NOT qualified"
echo "EVIDENCE: $out"
for ranks in 1 2; do
 mkdir "$out/${ranks}r"
 echo "RUN original DEFAULT L2 MPI ${ranks} ranks"
 if ! "${docker_cmd[@]}" run --rm \
   --user "$(id -u):$(id -g)" \
   --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp \
   --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
   --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=csr \
   --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
   --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
   --volume "$root/sandbox/probe_v130_fe_initial_mpi.py:/pilot.py:ro" \
   --volume "$root/sandbox/probe_v130_universal_fem_initialization.py:/original_pilot.py:ro" \
   --volume "$out/${ranks}r:/artifacts:rw" \
   --workdir /artifacts \
   "$image" timeout -k 5s "${SCRIMP_V130_FE_MPI_TIMEOUT:-150s}" \
   mpiexec --oversubscribe --bind-to none -n "$ranks" \
   python -u /pilot.py --outdir /artifacts \
   > "$out/${ranks}r/native.host.log" 2>&1; then
   echo "FAIL native DEFAULT L2 MPI${ranks}; original logs retained: $out" >&2
   tail -n 150 "$out/${ranks}r/native.host.log" >&2
   exit 1
 fi
done
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
import numpy as np
root=Path(sys.argv[1])
s=json.loads((root/"1r"/"rank-0.json").read_text())
r=[json.loads((root/"2r"/f"rank-{i}.json").read_text())
   for i in range(2)]
assert s["status"]=="PASS_RANK" and all(x["status"]=="PASS_RANK" for x in r)
assert s["mpi_ranks"]==1 and all(x["mpi_ranks"]==2 for x in r)
assert len(s["cases"])==8
results={}
for name,v in s["cases"].items():
 assert v["status"].startswith("PASS_")
 differences=[]
 integral_ratios=[]
 for q in r:
  w=q["cases"][name]
  assert w["status"].startswith("PASS_"),(name,w)
  if name=="lumped":
   assert w["value"]==v["value"]==2.5
   continue
  x=np.asarray(v["global_coefficients"],dtype=float)
  y=np.asarray(w["global_coefficients"],dtype=float)
  assert x.shape==y.shape,(name,x.shape,y.shape)
  rel=float(np.linalg.norm(y-x)/max(np.linalg.norm(x),1e-30))
  differences.append(rel)
  assert rel<1e-9,(name,q["rank"],rel)
  integral_ratios.append(w["getfem_raw_target_L2_squared"]/
                         v["getfem_raw_target_L2_squared"])
  if name!="gaussian_cg1_scalar":
   assert w["physical_relative_L2_error"]<5e-8,(name,q["rank"])
 results[name]={
   "max_MPI2_vs_MPI1_relative_coefficient_error":
    max(differences) if differences else 0.0,
   "raw_getfem_target_integral_MPI2_over_MPI1":
    integral_ratios,
 }
assert abs(s["getfem_rhs_multiplicative_scale"]-1.)<1e-12
assert all(abs(x["getfem_rhs_multiplicative_scale"]-1.)<1e-12
           or abs(x["getfem_rhs_multiplicative_scale"]-0.5)<1e-12
           for x in r)
summary={
 "schema":"scrimp-v130-default-fe-l2-mpi-comparison-v1",
 "status":"PASS_DEFAULT_L2_MPI1_MPI2_COEFFICIENTS_AND_FIELDS",
 "mpi2_rhs_multipliers":[x["getfem_rhs_multiplicative_scale"] for x in r],
 "cases":results,
 "lumped_numeric_default_passed":True,
 "DAE_and_full_physics_qualified":False,
 "published_sources_unchanged":True,
}
(root/"comparison.json").write_text(json.dumps(summary,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as z:
 for path in sorted(root.rglob("*")):
  if path.is_file(): z.write(path,arcname=str(path.relative_to(root)))
print(json.dumps(summary,indent=2))
print("SHAREABLE ORIGINAL DEFAULT FE MPI1/2 PROJECTION EVIDENCE:",archive)
PY
echo "COMPLETE DEFAULT FEM L2 MPI1/2 native comparison: $out"
