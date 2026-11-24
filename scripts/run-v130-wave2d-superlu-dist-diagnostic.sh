#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.
#
# Diagnose the failed FIRST TSCN step of the real 2D wave with SuperLU_DIST.
# No solver retries, factor tuning, physics/threshold changes or legacy fallbacks.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0 (no positional arguments)" >&2; exit 2; }
ranks="${SCRIMP_V130_WAVE2D_DIAG_RANKS:-1}"
case "$ranks" in
  1|2) ;;
  *) echo "Set SCRIMP_V130_WAVE2D_DIAG_RANKS=1 or 2" >&2; exit 2 ;;
esac
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
  echo "Unqualified native Docker image revision: $revision" >&2
  exit 2
}
python3 -m py_compile examples/wave.py tests/test_v130_wave2d_conservative_pilot.py
python3 tests/test_v130_wave2d_conservative_pilot.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-wave2d-superlu-dist-first-step-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
mkdir -p "$out/native" "$out/work"
printf '{"schema":"wave2d-superlu-dist-first-CN-step-diagnostic-v1","source":"%s","image":"%s","image_revision":"%s","mpi":%s,"not_full_acceptance":true,"baseline_evidence":"v130-wave2d-conservative-five-lu-99d99f2d-tEOjqgHy.zip"}\n' \
  "$(git rev-parse HEAD)" "$image" "$revision" "$ranks" > "$out/manifest.json"
echo "DIAG: SuperLU_DIST MPI=$ranks, same 2230-DOF FE, dt=0.008, one step"
echo "DIAG: existing TS/SNES/KSP/PC failure introspection (serial only)"
echo "EVIDENCE: $out"
rc=0
"${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp --env MPLBACKEND=Agg \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
    --env SCRIMP_TS_FAILURE_DIAGNOSTICS=1 \
    --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
    --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
    --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
    --env SCRIMP_MATRIX_MEMORY_MODE=compact \
    --env SCRIMP_MUMPS_RESET_POLICY=on_change \
    --env SCRIMP_GETFEM_RHS_MODE=auto \
    --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
    --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
    --volume "$root/examples/wave.py:/workspace/scrimp/examples/wave.py:ro" \
    --volume "$out/native:/artifacts:rw" \
    --volume "$out/work:/scrimp-work:rw" \
    --workdir /scrimp-work \
    "$image" timeout -k 5s "${SCRIMP_V130_WAVE2D_DIAG_TIMEOUT:-120s}" \
    mpiexec --oversubscribe --bind-to none -n "$ranks" \
    python -u -m examples.wave \
      --solver superlu_dist --mesh-size 0.2 \
      --final-time 0.008 --time-step 0.008 \
      --check --output /artifacts \
    > "$out/native.host.log" 2>&1 || rc=$?
printf '%s\n' "$rc" > "$out/native.exitcode"
python3 - "$out" "$rc" <<'PY'
from pathlib import Path
import json
import sys
import zipfile
folder = Path(sys.argv[1])
exit_code = int(sys.argv[2])
log = (folder/"native.host.log").read_text(errors="replace")
messages = [line for line in log.splitlines() if any(
    needle in line for needle in (
        "MAIN TS FAILURE DIAGNOSTIC:",
        "MAIN TS FAILURE diagnostic unavailable:",
        "MAIN TS PETSc linear-solve failure:",
        "DIVERGED_STEP_REJECTED",
        "KSP_DIVERGED",
        "SNES_DIVERGED",
    )
)]
summary = {
    "status": (
        "FIRST_STEP_NATIVE_PASS_NOT_FIVE_LU_ACCEPTANCE" if not exit_code
        else "SUPERLU_DIST_FIRST_STEP_NATIVE_FAILURE_RECORDED"
    ),
    "exit_code": exit_code, "diagnostic_lines": messages,
    "missing_serial_diagnostics": bool(exit_code and not any(
        "MAIN TS FAILURE DIAGNOSTIC:" in line for line in messages
    )),
    "full_five_backend_gate_passed": False,
    "original_three_passes_preserved": True,
}
(folder/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
archive = Path(str(folder)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     allowZip64=True) as z:
    for file in sorted(folder.rglob("*")):
        if file.is_file(): z.write(file,arcname=str(file.relative_to(folder)))
print(json.dumps({**summary,"archive":str(archive)},indent=2),flush=True)
PY
if [[ "$rc" != 0 ]]; then
  echo "SUPERLU_DIST_UNQUALIFIED — error and solver reasons preserved" >&2
  exit "$rc"
fi
echo "FIRST_STEP_ONLY_PASS — full 10-step/five-backend acceptance still pending"
