#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# OPTIONAL one-solver long visualisation run (not a five-solver parity test).
# Genuine GetFEM/PETSc fixed CN, actual q/p VTU/PVD, PNG/SVG/GIF.
set -euo pipefail
root="$(cd "$(dirname "\${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
docker_cmd=(sudo docker)
if [[ "\${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="\${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("\${docker_cmd[@]}" image inspect \
 --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
expected="1d6d0c90115213298047b28237a7f38087563a35"
[[ "$revision" == "$expected" ]] || {
  echo "Cannot render with unqualified pinned GetFEM/PETSc image" >&2
  exit 2
}
python3 -m py_compile scripts/render-v130-wave1d-doc.py examples/wave_1d.py
bash -n "$0"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-wave1d-long-doc-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
mkdir -p "$out/native" "$out/work"
printf '{"schema":"scrimp-v130-wave1d-native-documentation-v1","source":"%s","pinned_image":"%s","image_revision":"%s","solver":"mumps","MPI":1,"time_final":1.5,"dt":0.01,"NOT_five_backend_acceptance":true}\n' \
 "$(git rev-parse HEAD)" "$image" "$revision" > "$out/native/manifest.json"
echo "OPTIONAL real native wave 1D longer documentation illustration"
echo "150 fixed CN steps, 151 accepted states, true GetFEM VTU/PVD"
echo "EVIDENCE: $out"
status=0
"\${docker_cmd[@]}" run --rm \
 --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env MPLBACKEND=Agg \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
 --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
 --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
 --env SCRIMP_MATRIX_MEMORY_MODE=compact \
 --env SCRIMP_GETFEM_RHS_MODE=auto \
 --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
 --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
 --volume "$root/examples/wave_1d.py:/workspace/scrimp/examples/wave_1d.py:ro" \
 --volume "$out/native:/artifacts:rw" \
 --volume "$out/work:/scrimp-work:rw" \
 --workdir /scrimp-work \
 "$image" timeout -k 5s "\${SCRIMP_V130_WAVE1D_LONG_TIMEOUT:-300s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u -m examples.wave_1d --solver mumps \
 --mesh-size 0.05 --time-step 0.01 --final-time 1.5 \
 --check --export --output /artifacts \
 > "$out/native.host.log" 2>&1 || status=$?
if [[ "$status" != 0 ]]; then
  echo "Long native tutorial failed (exit=$status)" >&2
  tail -n 80 "$out/native.host.log" >&2
else
  echo "Rendering actual accepted GetFEM fields and CN observations"
  python3 scripts/render-v130-wave1d-doc.py \
    --input "$out/native" --output "$out/figures" \
    --case mumps_mpi1 --frame-ms 100 || status=$?
fi
python3 - "$out" <<'PY'
from pathlib import Path
import sys
import zipfile
root=Path(sys.argv[1])
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for f in sorted(root.rglob("*")):
        if f.is_file():
            z.write(f,arcname=str(f.relative_to(root)))
print("LONG WAVE1D PHYSICS AND FIGURES ARCHIVE:",archive)
PY
if [[ "$status" != 0 ]]; then
  echo "FAIL_OR_UNQUALIFIED_LONG_WAVE_1D_DOC_FIGURES" >&2
  exit "$status"
fi
echo "PASS_REAL_NATIVE_LONG_WAVE_1D_DOC_FIGURES"
echo "Long GIF pending human visual review: $out/figures/wave_1d_string.gif"
