#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Native production-default FE L2 projection: fifteen spaces and Gaussian.
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
 echo "Qualified rc3 image revision required: $revision" >&2
 exit 2
}
python3 -m py_compile sandbox/probe_v130_production_fe_initialization.py
python3 tests/test_v130_fe_initial_projection.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-production-fe-l2-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-native-production-fe-l2-v1","source":"%s","image_revision":"%s","native_status":"PENDING"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "ACTUAL DPHS.set_initial_value() default L2, PETSc sparse MUMPS, MPI1"
echo "15 scalar/vector/derivative FEM spaces plus Gaussian nonpolynomial CG1"
echo "New production L2 default, explicit legacy interpolation, no image rebuild"
echo "EVIDENCE: $out"
if ! "${docker_cmd[@]}" run --rm \
 --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
 --volume "$root/sandbox/probe_v130_production_fe_initialization.py:/pilot.py:ro" \
 --volume "$root/sandbox/probe_v130_universal_fem_initialization.py:/original_pilot.py:ro" \
 --volume "$out:/artifacts:rw" \
 --workdir /artifacts \
 "$image" timeout -k 5s "${SCRIMP_V130_FE_INIT_TIMEOUT:-180s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /pilot.py --out /artifacts/results.json \
 > "$out/native.host.log" 2>&1; then
 echo "FAIL native production FE L2 test; original logs: $out" >&2
 tail -n 120 "$out/native.host.log" >&2
 exit 1
fi
grep '^PASS ACTUAL PRODUCTION SPARSE FEM L2 15 PLUS GAUSSIAN$' "$out/native.host.log"
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
p=Path(sys.argv[1])
r=json.loads((p/"results.json").read_text())
assert r["status"]=="PASS_NATIVE_SPARSE_15_PLUS_NONPOLYNOMIAL"
assert len(r["cases"])==16
assert r["MPI2_qualified"] is False
assert r["DAE_physics_qualified"] is False
assert r["default_initialization_method"]=="l2"
assert r["legacy_interpolation_explicit"] is True
archive=Path(str(p)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(p.glob("*")):
  if f.is_file(): z.write(f,arcname=f.name)
print(json.dumps(r,indent=2))
print("SHAREABLE ORIGINAL NATIVE PRODUCTION FE L2:",archive)
PY
echo "COMPLETE actual production opt-in sparse FE projection evidence: $out"
