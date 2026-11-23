#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Native FE-wide scalar, intrinsic vector and derivative DOF mass projection.
# Diagnostics only: production DPHS.set_initial_value is unmodified.
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
 echo "Unqualified GetFEM/PETSc image revision: $revision" >&2
 exit 2
}
python3 -m py_compile sandbox/probe_v130_universal_fem_initialization.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-universal-fem-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-universal-fem-pilot-v1","source":"%s","image_revision":"%s","native_status":"PENDING"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Generic GetFEM L2 projection across CG/DG scalar/vector, Hermite, Argyris, Morley, HCT, RT, BDM, Nedelec"
echo "Compare native physical-field L2 reconstruction with unchanged actual production initializer"
echo "Small dense MPI1 pilot: NOT sparse PETSc, canonical derivative interpolation, MPI or DAE qualification"
echo "EVIDENCE: $out"
if ! "${docker_cmd[@]}" run --rm \
 --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
 --volume "$root/sandbox/probe_v130_universal_fem_initialization.py:/pilot.py:ro" \
 --volume "$out:/artifacts:rw" \
 --workdir /artifacts \
 "$image" timeout -k 5s "${SCRIMP_V130_FE_INIT_TIMEOUT:-180s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /pilot.py --out /artifacts/results.json \
 > "$out/native.host.log" 2>&1; then
 echo "FAIL native universal FE pilot; logs and any JSON: $out" >&2
 tail -n 110 "$out/native.host.log" >&2
 exit 1
fi
grep '^PASS NATIVE UNIVERSAL SCALAR VECTOR HIGH REGULARITY PROJECTION$' "$out/native.host.log"
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
p=Path(sys.argv[1])
obj=json.loads((p/"results.json").read_text())
assert obj["status"]=="PASS_REQUIRED_AND_ALL_AVAILABLE_FEM"
assert len(obj["mandatory_cases"])==9
assert obj["production_changed"] is False
assert obj["MPI_qualified"] is False
assert obj["canonical_derivative_moment_interpolation_qualified"] is False
archive=Path(str(p)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(p.glob("*")):
  if f.is_file(): z.write(f,arcname=f.name)
print(json.dumps(obj,indent=2))
print("SHAREABLE ORIGINAL UNIVERSAL FE PILOT:",archive)
PY
echo "COMPLETE native universal FE diagnostic: $out"
