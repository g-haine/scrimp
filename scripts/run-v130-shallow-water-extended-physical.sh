#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# Genuine canonical inviscid Saint-Venant example, fixed T=0.002:
# BDF2 dt=0.0001 and 0.00005, MPI1+MPI2, original FEM and data.
# No exact BDF2 energy conservation or pseudo-TS qualification claim.
# Native logs, full original FE rank trajectories and independent host
# gates are archived even after rejection of one physical run.
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
  echo "Incorrect pinned GetFEM/PETSc image revision $revision" >&2; exit 2;
}
python3 -m py_compile \
  examples/shallow_water.py \
  sandbox/probe_v130_issue28_shallow_short_physical.py \
  scripts/verify-v130-shallow-water-extended.py \
  tests/test_v130_shallow_water_extended_source.py
bash -n "$0"
python3 tests/test_v130_shallow_water_extended_source.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-shallow-original-extended-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"v130-shallow-original-extended-physical-v1","source":"%s","image_revision":"%s","model":"examples/shallow_water.py","fixed_final_time":0.002,"steps":[20,40],"MPI":[1,2],"source_published_immutable":true,"initialization":"nonlinear_affine","integrator":"BDF2","gyro":"explicit_previous_accepted_state","exact_BDF2_energy_conservation_claimed":false,"rotating_viscous_SWE_qualified":false}\n' \
  "$(git rev-parse HEAD)" "$revision" >"$out/manifest.json"

echo "INVISCID SAINT-VENANT 2D: authentic canonical nonlinear mass BDF2 physical gate"
echo "Same T=.002, dt=.0001 (20 steps) and .00005 (40 steps), MUMPS MPI1/2"
echo "Strict original nonlinear_affine, FE full states each rank and actual port work."
echo "NO exact energy-preservation claim for lagged explicit gyro + BDF2."
echo "EVIDENCE: $out"

for sample in "20:.0001" "40:.00005"; do
  steps="${sample%%:*}"
  dt="${sample#*:}"
  for ranks in 1 2; do
    label="steps${steps}_mpi${ranks}"
    mkdir "$out/$label" "$out/work_$label"
    echo "RUN $label dt=$dt final_time=0.002"
    rc=0
    "${docker_cmd[@]}" run --rm \
      --user "$(id -u):$(id -g)" \
      --env HOME=/tmp --env MPLBACKEND=Agg \
      --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
      --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
      --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
      --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
      --env SCRIMP_MATRIX_MEMORY_MODE=compact \
      --env SCRIMP_MUMPS_RESET_POLICY=on_change \
      --env SCRIMP_GETFEM_RHS_MODE=auto \
      --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
      --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
      --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
      --volume "$root/sandbox/probe_v130_issue28_shallow_short_physical.py:/workspace/scrimp/sandbox/probe_v130_issue28_shallow_short_physical.py:ro" \
      --volume "$out/$label:/artifacts:rw" \
      --volume "$out/work_$label:/scrimp-work:rw" \
      --workdir /scrimp-work \
      "$image" timeout -k 5s "${SCRIMP_V130_SW_EXTENDED_TIMEOUT:-420s}" \
      mpiexec --oversubscribe --bind-to none -n "$ranks" \
      python -u /workspace/scrimp/sandbox/probe_v130_issue28_shallow_short_physical.py \
      --out /artifacts --dt "$dt" --final-time 0.002 --capture-full-fe \
      >"$out/$label.host.log" 2>&1 || rc=$?
    printf '%s\n' "$rc" >"$out/$label.exitcode"
    if [[ "$rc" -ne 0 ]]; then
      echo "UNQUALIFIED native $label exit=$rc; evidence retained" >&2
      tail -n 65 "$out/$label.host.log" >&2
    else
      grep -a "ISSUE28 SHALLOW SHORT PHYSICAL" "$out/$label.host.log" || true
    fi
  done
done
verify_rc=0
python3 scripts/verify-v130-shallow-water-extended.py "$out" \
  >"$out/verify.host.log" 2>&1 || verify_rc=$?
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
report=json.loads((root/"summary.json").read_text())
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for f in sorted(root.rglob("*")):
        if f.is_file():
            z.write(f,arcname=str(f.relative_to(root)))
print(json.dumps({**report,"evidence":str(root),
                  "archive":str(archive)},indent=2),flush=True)
PY
exit "$verify_rc"
