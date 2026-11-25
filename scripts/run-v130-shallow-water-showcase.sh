#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# ORIGINAL rectangular inviscid Saint-Venant 2D h CG2 / p DG1.
# True GetFEM native VTU, original height and velocity PNG/GIF,
# authentic independent GetFEM Hamiltonian and mass integrals.
# Do NOT run viscous rotating research or change published examples.
# Source image pinned; no fake FE fields and no exact BDF2 energy claim.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
docker_cmd=(sudo docker)
[[ "${SCRIMP_DOCKER_NO_SUDO:-0}" != 1 ]] || docker_cmd=(docker)
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == 1d6d0c90115213298047b28237a7f38087563a35 ]] || {
  echo "Wrong pinned scientific image revision $revision" >&2; exit 2;
}

tf="${SCRIMP_V130_SW_SHOWCASE_FINAL_TIME:-0.12}"
dt="${SCRIMP_V130_SW_SHOWCASE_DT:-0.0001}"
save="${SCRIMP_V130_SW_SHOWCASE_SAVE_STEP:-0.002}"
python3 - "$tf" "$dt" "$save" <<'PY'
import math,sys
tf,dt,ds=map(float,sys.argv[1:])
steps=tf/dt
if not (all(math.isfinite(x) and x>0 for x in (tf,dt,ds))
        and 5<=round(steps)<=6000
        and math.isclose(steps,round(steps),abs_tol=1e-8)
        and ds>=dt and ds<=tf
        and math.isclose(ds/dt,round(ds/dt),abs_tol=1e-8)):
    raise SystemExit("Choose 5..6000 fixed BDF2 steps, dt<=save<=tf and integral save/dt")
PY
python3 -m py_compile \
  examples/shallow_water.py \
  scripts/verify-v130-shallow-water-showcase.py \
  tests/test_v130_shallow_water_visualization_source.py
bash -n "$0"
python3 tests/test_v130_shallow_water_visualization_source.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-SWE-native-visuals-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
mkdir "$out/native" "$out/work"
printf '{"schema":"v130-original-SWE-showcase-v1","source":"%s","image":"%s","image_revision":"%s","model":"examples/shallow_water.py","MPI":1,"initialization":"nonlinear_affine","time_integrator":"BDF2 lagged original gyroscopic G","final_time":%s,"dt":%s,"save_step":%s,"original_PDE_unchanged":true,"rotating_viscous_research_included":false,"exact_H_conservation_claimed":false}\n' \
  "$(git rev-parse HEAD)" "$image" "$revision" "$tf" "$dt" "$save" > "$out/manifest.json"
echo "INVISCID SHALLOW WATER 2D NATIVE VISUALS: original h CG2 and p DG1"
echo "BDF2 nonlinear_affine showcase only; t_f=$tf, dt=$dt, save_step=$save"
echo "3D native height, 2D height map, velocity p/rho and original H/mass."
echo "No synthetic waves, no rotating viscous research, no exact BDF2 H claim."
echo "EVIDENCE: $out"

rc=0
"${docker_cmd[@]}" run --rm \
  --user "$(id -u):$(id -g)" --env HOME=/tmp \
  --env MPLBACKEND=Agg --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env PYTHONPATH=/workspace/scrimp \
  --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
  --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
  --env SCRIMP_MATRIX_MEMORY_MODE=compact \
  --env SCRIMP_MUMPS_RESET_POLICY=on_change \
  --env SCRIMP_GETFEM_RHS_MODE=auto \
  --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
  --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
  --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
  --volume "$out/native:/artifacts:rw" \
  --volume "$out/work:/scrimp-work:rw" \
  --workdir /scrimp-work "$image" \
  timeout -k 10s "${SCRIMP_V130_SW_SHOWCASE_TIMEOUT:-3600s}" \
  mpiexec --oversubscribe --bind-to none -n 1 \
  python -u -m examples.shallow_water \
    --showcase-output /artifacts \
    --final-time "$tf" --time-step "$dt" --save-step "$save" \
    --factor-solver mumps \
  >"$out/native.host.log" 2>&1 || rc=$?
printf '%s\n' "$rc" > "$out/native.exitcode"
if [[ "$rc" -ne 0 ]]; then
  echo "UNQUALIFIED original BDF2 native visual run exit=$rc" >&2
  tail -n 85 "$out/native.host.log" >&2
fi

verify_rc=0
python3 scripts/verify-v130-shallow-water-showcase.py "$out" \
  >"$out/verify.host.log" 2>&1 || verify_rc=$?
if [[ "$verify_rc" -ne 0 ]]; then
  echo "UNQUALIFIED visual evidence; original saved native fields retained" >&2
  tail -n 65 "$out/verify.host.log" >&2
fi
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
report=json.loads((root/"summary.json").read_text())
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for f in sorted(root.rglob("*")):
        if f.is_file():z.write(f,arcname=str(f.relative_to(root)))
print(json.dumps({**report,"archive":str(archive),
                  "evidence":str(root)},indent=2),flush=True)
PY
exit "$verify_rc"
