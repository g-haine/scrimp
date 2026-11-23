#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Deliberately separate native GetFEM L2-projection experiment.
# Does not alter SCRIMP production initialization or published sources.
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
 echo "Wrong GetFEM/PETSc image revision: $revision" >&2
 exit 2
}
python3 -m py_compile sandbox/probe_v130_intrinsic_l2_projection.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-intrinsic-l2-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-intrinsic-l2-pilot-v1","source":"%s","image_revision":"%s","native_status":"PENDING"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "DIAGNOSTIC ONLY: original GetFEM RTk/BDMk/Nedelec L2 projection"
echo "Captures native MeshFem.eval shape and current DPHS.set_initial_value rejection"
echo "Small, MPI1, isolated native GetFEM mass projection. NOT production qualification"
echo "EVIDENCE: $out"
if ! "${docker_cmd[@]}" run --rm \
 --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
 --volume "$root/sandbox/probe_v130_intrinsic_l2_projection.py:/pilot.py:ro" \
 --volume "$out:/artifacts:rw" \
 --workdir /artifacts \
 "$image" timeout -k 5s "${SCRIMP_V130_FE_INIT_TIMEOUT:-90s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /pilot.py --out /artifacts/results.json \
 > "$out/native.host.log" 2>&1; then
 echo "FAIL original native intrinsic projection diagnostic; evidence: $out" >&2
 tail -n 120 "$out/native.host.log" >&2
 exit 1
fi
grep '^PASS GETFEM INTRINSIC L2 DIAGNOSTIC ONLY$' "$out/native.host.log"
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
p=Path(sys.argv[1])
d=json.loads((p/"results.json").read_text())
assert d["status"]=="PASS_NATIVE_GETFEM_MASS_PROJECTION_DIAGNOSTIC_ONLY"
assert d["production_initialization_qualified"] is False
assert d["DAE_and_MPI_qualified"] is False
archive=Path(str(p)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(p.glob("*")):
  if f.is_file(): z.write(f,arcname=f.name)
print("SHAREABLE ORIGINAL GETFEM PROJECTION DIAGNOSTIC:",archive)
PY
echo "COMPLETE NATIVE GETFEM L2 DIAGNOSTIC ONLY: $out"
