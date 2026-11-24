#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Real PETSc PCPYTHON (P A)^(-1) P on original Heat-Wave first CN and TS.
# Stage 1 snapshot KSP must pass before stage 2 physical 4-step CN.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 1 ]] || {
 echo "Usage: $0 sandbox/outputs/.../physical/first_cn_J_001.petscbin" >&2
 exit 2
}
matrix="$(realpath "$1")"
[[ -s "$matrix" && "$(basename "$matrix")" == "first_cn_J_001.petscbin" ]] || {
 echo "Require original 7689x7689 CN Jacobian binary" >&2;exit 2
}
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
 --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
 echo "Pinned GetFEM/PETSc image required: $revision" >&2; exit 2
}
python3 -m py_compile \
 sandbox/probe_v130_hw_matching_pc.py \
 sandbox/probe_v130_hw_matched_pc_snapshot.py \
 sandbox/probe_v130_hw_matched_pc_cn.py
python3 tests/test_v130_hw_matched_pc.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-hw-matched-PC-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-HW-original-matched-PC-v1","source":"%s","image_revision":"%s","original_binary":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" "$matrix" > "$out/manifest.json"
mkdir "$out/snapshot" "$out/physical" "$out/work"
echo "Original Heat-Wave PETSc Python PC: exact (P A)^(-1) P b"
echo "Stage 1 identical CN snapshot; stage 2 actual 4-step TS CN, clamp"
echo "MPI1 ONLY; normal SCRIMP DPHS unchanged; source examples read-only"
echo "No regularization, no change to physical state or energy/power monitor"
echo "EVIDENCE: $out"
if ! "${docker_cmd[@]}" run --rm --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --volume "$(dirname "$matrix"):/original:ro" \
 --volume "$root/sandbox/probe_v130_hw_cn_identical_lu.py:/replay_base.py:ro" \
 --volume "$root/sandbox/probe_v130_hw_cn_exact_row_matching.py:/probe_v130_hw_cn_exact_row_matching.py:ro" \
 --volume "$root/sandbox/probe_v130_hw_matching_pc.py:/probe_v130_hw_matching_pc.py:ro" \
 --volume "$root/sandbox/probe_v130_hw_matched_pc_snapshot.py:/snapshot.py:ro" \
 --volume "$out/snapshot:/artifacts:rw" --workdir /artifacts \
 "$image" timeout -k 5s "${SCRIMP_V130_MATCHED_PC_TIMEOUT:-240s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /snapshot.py \
 --matrix /original/first_cn_J_001.petscbin \
 --out /artifacts/snapshot.json \
 > "$out/snapshot.host.log" 2>&1; then
 echo "FAIL matched Python PC on identical original CN snapshot: $out" >&2
 tail -n 140 "$out/snapshot.host.log" >&2;exit 1
fi
echo "PASS matched PC original snapshot; RUN real original Heat-Wave CN clamp"
if ! "${docker_cmd[@]}" run --rm --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env MPLBACKEND=Agg \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
 --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
 --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
 --env SCRIMP_MATRIX_MEMORY_MODE=compact \
 --env SCRIMP_MUMPS_RESET_POLICY=on_change \
 --env SCRIMP_GETFEM_RHS_MODE=auto \
 --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
 --env SCRIMP_TS_FAILURE_DIAGNOSTICS=1 \
 --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
 --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
 --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
 --volume "$root/sandbox/probe_v130_main_2d_examples.py:/workspace/scrimp/sandbox/probe_v130_main_2d_examples.py:ro" \
 --volume "$root/sandbox/measure_v130_heat_wave_physics.py:/workspace/scrimp/sandbox/measure_v130_heat_wave_physics.py:ro" \
 --volume "$root/sandbox/probe_v130_heat_wave_boundary_compatibility.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_boundary_compatibility.py:ro" \
 --volume "$root/sandbox/probe_v130_heat_wave_compatible_cn.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_compatible_cn.py:ro" \
 --volume "$root/sandbox/probe_v130_hw_cn_identical_lu.py:/workspace/scrimp/sandbox/replay_base.py:ro" \
 --volume "$root/sandbox/probe_v130_hw_cn_exact_row_matching.py:/workspace/scrimp/sandbox/probe_v130_hw_cn_exact_row_matching.py:ro" \
 --volume "$root/sandbox/probe_v130_hw_matching_pc.py:/workspace/scrimp/sandbox/probe_v130_hw_matching_pc.py:ro" \
 --volume "$root/sandbox/probe_v130_hw_matched_pc_cn.py:/workspace/scrimp/sandbox/probe_v130_hw_matched_pc_cn.py:ro" \
 --volume "$out/physical:/artifacts:rw" \
 --volume "$out/work:/scrimp-work:rw" --workdir /scrimp-work \
 "$image" timeout -k 5s "${SCRIMP_V130_MATCHED_PC_TIMEOUT:-240s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /workspace/scrimp/sandbox/probe_v130_hw_matched_pc_cn.py \
 --candidate clamp --out /artifacts \
 > "$out/physical.host.log" 2>&1; then
 echo "FAIL original physical CN under experimental matched PC: $out" >&2
 tail -n 160 "$out/physical.host.log" >&2; exit 1
fi
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
s=json.loads((root/"snapshot/snapshot.json").read_text())
p=json.loads((root/"physical/rowmatched_cn_audit.json").read_text())
numerical=json.loads((root/"physical/original_2d_example_report.json").read_text())
assert s["status"]=="PASS_SNAPSHOT_PYTHON_PC_ORIGINAL_EQUATIONS_ONLY"
assert s["matched_PC"]["setup_count"]>=1 and s["matched_PC"]["apply_count"]>=1
assert p["native_rowmatched_PC"]["setup_count"]>=1
assert p["native_rowmatched_PC"]["apply_count"]>=4
assert p["true_original_Ts_4steps"] and numerical["steps"]==4
assert not p["automatically_enabled_in_production"]
assert not p["MPI2_qualified"]
assert not p["original_example_modified"]
summary={
 "schema":"scrimp-v130-original-heat-wave-matched-PC-CN-v1",
 "status": "PASS_REAL_CN_NUMERICAL_AND_PHYSICAL_THRESHOLDS_MPI1"
           if p["scientific_energy_power_gate_passed"] else
           "PASS_REAL_CN_NUMERICAL_SCIENTIFIC_BALANCES_OPEN_MPI1",
 "snapshot_original_relative_residual":s["residual_original_equations"],
 "snapshot_known_state_error":s["known_state_relative_error"],
 "physical_initial_candidate":p["candidate"],
 "strict_DAE_status":numerical["initialization_report"]["status"] if "initialization_report" in numerical else "see compatible_cn_report.json",
 "four_accepted_steps":numerical["steps"],
 "physical_thresholds":p["actual_original_physics_thresholds"],
 "native_rowmatched_PC":p["native_rowmatched_PC"],
 "no_production_change":True,
 "no_MPI2_or_parallel_export_qualification":True,
}
(root/"comparison.json").write_text(json.dumps(summary,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,allowZip64=True) as z:
 for f in sorted(root.rglob("*")):
  if f.is_file():z.write(f,arcname=str(f.relative_to(root)))
print(json.dumps(summary,indent=2))
print("SHAREABLE MATCHED ORIGINAL CN/TS EVIDENCE:",archive)
PY
echo "COMPLETE original Heat-Wave real TS matched-PC audit: $out"
