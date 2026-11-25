#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# Stage 1: genuine Heat-Wave MPI1 vs MPI2 four-step DAE/CN physics, full FE.
# Stage 2: genuine Heat-Wave MUMPS MPI1 long-horizon EVERY CN STEP.
# Fail closed and archive every run, even a failed native execution.
# No parallel plots, no unsafe shared ParaView writers, no fake CN midpoints.
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
  echo "Wrong pinned Docker scientific revision: $revision" >&2; exit 2;
}
long_profile="${SCRIMP_V130_HW_LONG_PROFILE:-gaussian}"
case "$long_profile" in gaussian|gallery_packet) ;; *)
  echo "Invalid long profile (gaussian or gallery_packet)" >&2; exit 2;;
esac
long_steps="${SCRIMP_V130_HW_LONG_STEPS:-3000}"
long_dt="${SCRIMP_V130_HW_LONG_DT:-0.005}"
python3 - "$long_steps" "$long_dt" <<'PY'
import math,sys
steps=int(sys.argv[1]); dt=float(sys.argv[2])
if not (4<=steps<=10000 and math.isfinite(dt) and 0.<dt<=.01):
    raise SystemExit("Long Heat-Wave steps/dt out of safety range")
PY
python3 -m py_compile \
  scrimp/dphs.py \
  scrimp/utils/trajectory_sampling.py \
  examples/heat_wave.py \
  sandbox/measure_v130_hw_dual_facet_physics.py \
  sandbox/probe_v130_heat_wave_mpi_extended.py \
  scripts/verify-v130-heat-wave-mpi-long-physics.py \
  tests/test_v130_trajectory_sampling.py \
  tests/test_v130_heat_wave_mpi_long_source.py
bash -n "$0"
python3 tests/test_v130_trajectory_sampling.py
python3 tests/test_v130_heat_wave_mpi_long_source.py
python3 tests/test_v130_heat_wave_boundary_initial_data.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-heat-wave-MPI2-longCN-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-heat-wave-MPI2-longCN-v1","source":"%s","image":"%s","image_revision":"%s","long_dt":%s,"long_steps":%s,"long_profile":"%s","short_stage":{"MPI":[1,2],"dt":0.001,"steps":4},"long_stage":{"MPI":1,"save_every_accepted_CN_step":true},"MPI2_ParaView_qualified":false,"SuperLU_DIST_qualified":false}\n' \
  "$(git rev-parse HEAD)" "$image" "$revision" "$long_dt" "$long_steps" "$long_profile" \
  > "$out/manifest.json"
echo "HEAT-WAVE MPI2 + long CN strict physical acceptance, original GetFEM/PETSc"
echo "SHORT: 4 CN steps dt=.001 MPI1 and MPI2, full FE and physical parities"
echo "LONG: $long_steps CN steps dt=$long_dt, profile=$long_profile MPI1, ALL accepted states"
echo "No parallel visuals or manufactured/under-sampled CN energy identities."
echo "EVIDENCE: $out"

run_case() {
  local label="$1" nranks="$2" dt="$3" steps="$4" profile="$5"
  mkdir -p "$out/$label" "$out/work_$label"
  echo "RUN $label MPI=$nranks steps=$steps dt=$dt profile=$profile"
  local rc=0
  local timeout_seconds="${SCRIMP_V130_HW_LONG_PHYSICAL_TIMEOUT:-7200s}"
  if [[ "$label" == mpi1_long ]]; then
    echo "Full 15s stage may use longer CPU for 3000 GetFEM CN power integrals"
  else
    timeout_seconds="${SCRIMP_V130_HW_MPI_PHYSICAL_TIMEOUT:-600s}"
  fi
  local additional=()
  if [[ "$label" != mpi1_long ]]; then
    additional=(--store-full-states)
  fi
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
    --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
    --volume "$root/sandbox/measure_v130_hw_dual_facet_physics.py:/workspace/scrimp/sandbox/measure_v130_hw_dual_facet_physics.py:ro" \
    --volume "$root/sandbox/probe_v130_heat_wave_mpi_extended.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_mpi_extended.py:ro" \
    --volume "$out/$label:/artifacts:rw" \
    --volume "$out/work_$label:/scrimp-work:rw" \
    --workdir /scrimp-work "$image" \
    timeout -k 10s "$timeout_seconds" \
    mpiexec --oversubscribe --bind-to none -n "$nranks" \
    python -u /workspace/scrimp/sandbox/probe_v130_heat_wave_mpi_extended.py \
      --out /artifacts --dt "$dt" --steps "$steps" --profile "$profile" \
      "${additional[@]}" \
    >"$out/$label.host.log" 2>&1 || rc=$?
  printf '%s\n' "$rc" > "$out/$label.exitcode"
  if [[ "$rc" -ne 0 ]]; then
    echo "UNQUALIFIED native case: $label exit=$rc" >&2
    tail -n 75 "$out/$label.host.log" >&2
  fi
}

run_case mpi1_short 1 .001 4 gaussian
run_case mpi2_short 2 .001 4 gaussian

short_rc=0
python3 scripts/verify-v130-heat-wave-mpi-long-physics.py \
  "$out" --phase short >"$out/short-verification.host.log" 2>&1 || short_rc=$?
if [[ "$short_rc" == 0 ]]; then
  echo "PASS native MPI1/MPI2 four-step full FE + physical parity"
  run_case mpi1_long 1 "$long_dt" "$long_steps" "$long_profile"
else
  echo "UNQUALIFIED short native MPI1/MPI2; long stage NOT RUN" >&2
  tail -n 80 "$out/short-verification.host.log" >&2
fi

final_rc=0
python3 scripts/verify-v130-heat-wave-mpi-long-physics.py \
  "$out" --phase complete >"$out/final-verification.host.log" 2>&1 || final_rc=$?
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
rep=json.loads((root/"summary.json").read_text())
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for file in sorted(root.rglob("*")):
        if file.is_file():
            z.write(file,arcname=str(file.relative_to(root)))
print(json.dumps({**rep,"archive":str(archive),"evidence":str(root)},
                 indent=2),flush=True)
PY
exit "$final_rc"
