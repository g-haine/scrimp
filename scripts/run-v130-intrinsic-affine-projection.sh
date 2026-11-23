#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Independent small-mesh exact affine FE projection (MPI1); research only.
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
python3 -m py_compile sandbox/probe_v130_intrinsic_affine_projection.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-intrinsic-affine-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-intrinsic-affine-pilot-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Six exact affine intrinsic vector FE fields (RTk, BDMk, Nedelec)"
echo "Independent native GetFEM mass and reconstructed-field L2 validation"
echo "EVIDENCE: $out"
if ! "${docker_cmd[@]}" run --rm \
 --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
 --volume "$root/sandbox/probe_v130_intrinsic_affine_projection.py:/pilot.py:ro" \
 --volume "$out:/artifacts:rw" \
 --workdir /artifacts \
 "$image" timeout -k 5s "${SCRIMP_V130_FE_INIT_TIMEOUT:-90s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /pilot.py --out /artifacts/results.json \
 > "$out/native.host.log" 2>&1; then
 echo "FAIL original GetFEM affine FE pilot; logs: $out" >&2
 tail -n 120 "$out/native.host.log" >&2
 exit 1
fi
grep '^PASS NATIVE SIX INTRINSIC AFFINE PHYSICAL L2$' "$out/native.host.log"
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
p=Path(sys.argv[1])
report=json.loads((p/"results.json").read_text())
assert report["status"]=="PASS_NATIVE_AFFINE_SIX_FAMILY"
assert len(report["cases"])==6
for v in report["cases"].values():
 assert v["status"]=="PASS_NATIVE_AFFINE_FIELD_L2"
 assert v["field_relative_L2_error"]<5e-8
 assert v["mass_relative_residual"]<1e-10
assert report["MPI_qualified"] is False
assert report["production_initialization_changed"] is False
archive=Path(str(p)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(p.glob("*")):
  if f.is_file(): z.write(f,arcname=f.name)
print(json.dumps(report,indent=2))
print("SHAREABLE NATIVE INTRINSIC AFFINE L2 EVIDENCE:",archive)
PY
echo "COMPLETE native intrinsic affine FE qualification: $out"
