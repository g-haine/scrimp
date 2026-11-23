#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Standalone actual GetFEM intrinsic vector FE coefficient compatibility.
# Does NOT certify higher-order or second-kind Nedelec availability,
# PDE well-posedness, MPI stability, or H(curl)/H(div) trace identities.
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
    echo "Refusing unqualified original GetFEM/PETSc image revision: $revision" >&2
    exit 2
}
python3 -m py_compile sandbox/probe_v130_intrinsic_vector_fem_initialization.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-intrinsic-fe-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-intrinsic-vector-fe-pilot-v1","source":"%s","image_revision":"%s","published_unchanged":true,"native_status":"PENDING"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Original GetFEM MeshFem RT0,RT1,RT2,BDM1,BDM2,Nedelec-lowest"
echo "Actual production DPHS.set_initial_value, filtered FE coefficient comparisons"
echo "No Docker rebuild, no published source/output mounts, MPI1 only"
echo "EVIDENCE: $out"
if ! "${docker_cmd[@]}" run --rm \
 --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
 --volume "$root/sandbox/probe_v130_intrinsic_vector_fem_initialization.py:/pilot.py:ro" \
 --volume "$out:/artifacts:rw" \
 --workdir /artifacts \
 "$image" timeout -k 5s "${SCRIMP_V130_FE_INIT_TIMEOUT:-90s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /pilot.py --out /artifacts/results.json \
 > "$out/native.host.log" 2>&1; then
    echo "FAIL native intrinsic FE test; original logs retained in $out" >&2
    tail -n 110 "$out/native.host.log" >&2
    exit 1
fi
grep '^FEM PASS' "$out/native.host.log" || {
    echo "Missing original native FE-family pass" >&2; exit 1
}
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
p=Path(sys.argv[1])
obj=json.loads((p/"results.json").read_text())
assert obj["native_status"]=="PASS_AVAILABLE_GETFEM_INTRINSIC_FE_COEFFICIENTS"
assert obj["cases"]["rt0"]["status"]=="PASS_ACTUAL_DPHS_SET_INITIAL_VALUE_ON_GETFEM"
assert obj["Nedelec_second_kind_supported_or_tested"] is False
assert obj["higher_order_Nedelec_supported_or_tested"] is False
archive=Path(str(p)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as f:
 for x in sorted(p.glob("*")):
  if x.is_file(): f.write(x,arcname=x.name)
print(json.dumps(obj,indent=2))
print("SHAREABLE ORIGINAL GETFEM INTRINSIC FE ZIP:",archive)
PY
echo "COMPLETE original native FE family coefficients evidence: $out"
