#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Capture original SuperLU first CN Jacobian, then factor the SAME
# immutable PETSc binary matrix with MUMPS and SuperLU_DIST.
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
 echo "Qualified GetFEM/PETSc Docker image required: $revision" >&2; exit 2;
}
python3 -m py_compile sandbox/probe_v130_hw_cn_identical_lu.py
python3 tests/test_v130_hw_cn_identical_lu.py
python3 tests/test_v130_hw_cn_identical_lu.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-hw-cn-identical-lu-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-hw-identical-LU-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
mkdir "$out/physical" "$out/work"
echo "Original Heat-Wave, clamp boundary, true DAE, SuperLU_DIST first CN stage"
echo "Capture actual production P via pre-existing SCRIMP_JACOBIAN_SNAPSHOT_PREFIX"
echo "Replay EXACT SAME PETSc matrix with MUMPS and SuperLU_DIST and inspect structure"
echo "No physical-TS PASS, no numerical tolerance changes, no article-era source edits"
echo "EVIDENCE: $out"
if ! "${docker_cmd[@]}" run --rm --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env MPLBACKEND=Agg \
 --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
 --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
 --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
 --env SCRIMP_MATRIX_MEMORY_MODE=compact \
 --env SCRIMP_MUMPS_RESET_POLICY=on_change \
 --env SCRIMP_GETFEM_RHS_MODE=auto \
 --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
 --env SCRIMP_TS_FAILURE_DIAGNOSTICS=1 \
 --env SCRIMP_JACOBIAN_SNAPSHOT_PREFIX=/artifacts/first_cn_J \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
 --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
 --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
 --volume "$root/sandbox/probe_v130_main_2d_examples.py:/workspace/scrimp/sandbox/probe_v130_main_2d_examples.py:ro" \
 --volume "$root/sandbox/measure_v130_heat_wave_physics.py:/workspace/scrimp/sandbox/measure_v130_heat_wave_physics.py:ro" \
 --volume "$root/sandbox/probe_v130_heat_wave_boundary_compatibility.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_boundary_compatibility.py:ro" \
 --volume "$root/sandbox/probe_v130_heat_wave_compatible_cn.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_compatible_cn.py:ro" \
 --volume "$out/physical:/artifacts:rw" \
 --volume "$out/work:/scrimp-work:rw" --workdir /scrimp-work \
 "$image" timeout -k 5s "${SCRIMP_V130_HW_LU_REPLAY_TIMEOUT:-240s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /workspace/scrimp/sandbox/probe_v130_heat_wave_compatible_cn.py \
  --candidate clamp --backend superlu_dist \
  --capture-ts-failure --out /artifacts \
 > "$out/physical.host.log" 2>&1; then
 echo "Original first-CN Jacobian capture failed; evidence retained: $out" >&2
 tail -n 130 "$out/physical.host.log" >&2; exit 1
fi
snapshot="$out/physical/first_cn_J_001.petscbin"
[[ -s "$snapshot" ]] || {
 echo "Missing actual FIRST CN assembled PETSc matrix: $snapshot" >&2
 exit 1
}
[[ -s "$out/physical/compatible_cn_failure.json" ]] || {
 echo "Missing original observed SuperLU rejection: $out" >&2
 exit 1
}
echo "RUN independent PETSc LU/MUMPS and LU/SuperLU_DIST on first original CN snapshot"
if ! "${docker_cmd[@]}" run --rm --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --volume "$root/sandbox/probe_v130_hw_cn_identical_lu.py:/replay.py:ro" \
 --volume "$out/physical:/artifacts:rw" --workdir /artifacts \
 "$image" timeout -k 5s "${SCRIMP_V130_HW_LU_REPLAY_TIMEOUT:-240s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /replay.py \
 --matrix /artifacts/first_cn_J_001.petscbin \
 --output /artifacts/replay_same_matrix.json \
 > "$out/replay.host.log" 2>&1; then
 echo "Matrix replay infrastructure failed; evidence retained: $out" >&2
 tail -n 170 "$out/replay.host.log" >&2; exit 1
fi
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
failure=json.loads((root/"physical/compatible_cn_failure.json").read_text())
replay=json.loads((root/"physical/replay_same_matrix.json").read_text())
assert failure["status"]=="DIAGNOSTIC_TS_REJECTED_NOT_PHYSICAL_PASS"
assert failure["strict_DAE_pass"] is True
assert replay["status"]=="DIAGNOSTIC_ONLY_NO_PHYSICAL_TS_CLAIM"
assert replay["matrix"]["size"]==7689
assert replay["same_original_matrix_for_both_factorizers"] is True
assert replay["full_Crank_Nicolson_TS_completed"] is False
assert replay["scientific_energy_power_qualified"] is False
results=replay["factorizations"]
assert set(results)=={"mumps","superlu_dist"}
assert results["mumps"]["operator_CSR_before"]==results[
 "superlu_dist"]["operator_CSR_before"]
summary={
 "schema":"scrimp-v130-original-CN-identical-LU-audit-v1",
 "status":"PASS_DIAGNOSTIC_EVIDENCE_CAPTURED_NOT_PHYSICAL_TS",
 "original_SuperLU_TS_rejected":True,
 "matrix":replay["matrix"],
 "factorizations":results,
 "model_or_solver_tol_modified":False,
 "physical_TS_or_power_balance_qualified":False,
}
(root/"comparison.json").write_text(json.dumps(summary,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
 for f in sorted(root.rglob("*")):
  if f.is_file():z.write(f,arcname=str(f.relative_to(root)))
print(json.dumps(summary,indent=2))
print("SHAREABLE ORIGINAL CN JACOBIAN REPLAY EVIDENCE:",archive)
PY
echo "COMPLETE actual original CN Jacobian LU replay: $out"
