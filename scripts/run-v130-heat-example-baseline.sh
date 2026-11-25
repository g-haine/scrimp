#!/usr/bin/env bash
# SCRIMP v1.3.0: reproduce the canonical 2D heat example, unchanged.
# First gate: native execution and its own Hamiltonian PNG; physics is separate.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ "$#" == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
  echo "Unqualified Docker image revision: $revision" >&2; exit 2;
}
python3 -m py_compile examples/heat.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-heat-original-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
mkdir -p "$out/native"
echo "Canonical 2D heat example, unchanged: original t_f=5 and BDF2"
echo "Docker: $image; MPI=1; GetFEM to PETSc: COO; MUMPS default"
echo "EVIDENCE: $out"

rc=0
"${docker_cmd[@]}" run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp --env MPLBACKEND=Agg \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
  --env PYTHONPATH=/workspace/scrimp \
  --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
  --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
  --env SCRIMP_MATRIX_MEMORY_MODE=compact \
  --env SCRIMP_MUMPS_RESET_POLICY=on_change \
  --env SCRIMP_GETFEM_RHS_MODE=auto \
  --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
  --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
  --volume "$root/examples/heat.py:/workspace/scrimp/examples/heat.py:ro" \
  --volume "$out/native:/artifacts:rw" \
  --workdir /artifacts \
  "$image" timeout -k 10s "${SCRIMP_V130_HEAT_TIMEOUT:-1800s}" \
  mpiexec --oversubscribe --bind-to none -n 1 \
  python -u /workspace/scrimp/examples/heat.py \
  >"$out/native.host.log" 2>&1 || rc=$?
printf '%s\n' "$rc" > "$out/native.exitcode"

final_rc=0
python3 - "$out" "$rc" "$revision" <<'PY' || final_rc=$?
import json
import re
import sys
import zipfile
from pathlib import Path

out = Path(sys.argv[1])
rc = int(sys.argv[2])
revision = sys.argv[3]
log = (out / "native.host.log").read_text(encoding="utf-8", errors="replace")
png = out / "native/outputs/png/Hamiltonian_Heat_2D.png"
step_matches = list(re.finditer(
    r"Steps:\s*(\d+)\s*\(\s*(\d+)\s+rejected,\s*"
    r"(\d+)\s+Nonlinear solver failures\s*\)", log
))
steps = int(step_matches[-1].group(1)) if step_matches else None
rejects = int(step_matches[-1].group(2)) if step_matches else None
snes_failures = int(step_matches[-1].group(3)) if step_matches else None
errors = []
if rc:
    errors.append(f"Native Docker command exit code {rc}")
if not png.is_file() or png.stat().st_size < 1000:
    errors.append("Canonical Hamiltonian_Heat_2D.png absent or empty")
if steps is None or steps < 1:
    errors.append("Missing successful PETSc time-step summary")
if snes_failures is not None and snes_failures:
    errors.append(f"Nonlinear PETSc solve failures: {snes_failures}")
status = ("PASS_NATIVE_HEAT_EXECUTION_AND_PNG_ONLY"
          if not errors else "FAIL_OR_UNQUALIFIED_NATIVE_HEAT")
summary = {
    "schema": "scrimp-v130-heat-canonical-baseline-v1",
    "status": status,
    "errors": errors,
    "original_example": "examples/heat.py",
    "original_model_and_time_configuration_unchanged": True,
    "docker_revision": revision,
    "MPI": 1,
    "steps": steps,
    "rejected_steps": rejects,
    "nonlinear_solver_failures": snes_failures,
    "native_exit_code": rc,
    "canonical_png": str(png),
    "canonical_png_bytes": png.stat().st_size if png.is_file() else None,
    "energy_power_balance_qualified": False,
    "MPI_parity_qualified": False,
    "not_a_scientific_physics_PASS": True,
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
archive = Path(str(out) + ".zip")
with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4, allowZip64=True) as zf:
    for file in sorted(out.rglob("*")):
        if file.is_file():
            zf.write(file, arcname=str(file.relative_to(out)))
print(json.dumps({**summary, "evidence": str(out),
                  "archive": str(archive)}, indent=2), flush=True)
raise SystemExit(bool(errors))
PY
final_rc=$?
if [[ "$final_rc" != 0 ]]; then
  echo "Heat example unqualified; last native log lines:" >&2
  tail -n 70 "$out/native.host.log" >&2
fi
exit "$final_rc"
