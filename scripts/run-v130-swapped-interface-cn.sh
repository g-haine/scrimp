#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# Genuine strict-Dae, four-step PETSc CN and sided pH energy for
# heat=annulus(2) / wave=disk(1) in the CURRENT Heat-Wave example.
# No legacy geometric/physical sign patch. MPI1 only.
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
  echo "Pinned image mismatch: $revision" >&2; exit 2
}

python3 -m py_compile \
  scrimp/interface.py scrimp/interface_deferred.py \
  scrimp/interface_discovery.py scrimp/port.py scrimp/brick.py \
  scrimp/dphs.py examples/heat_wave.py \
  sandbox/measure_v130_hw_dual_facet_physics.py \
  sandbox/probe_v130_swapped_interface_cn.py
python3 tests/test_v130_interface_topology.py
python3 tests/test_v130_interface_ownership.py
python3 tests/test_v130_interface_assembly_plan.py
python3 tests/test_v130_interface_deferred.py
python3 tests/test_v130_hw_facet_collector_regions.py
python3 tests/test_v130_heat_wave_outward_signs.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-swapped-interface-CN-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-swapped-real-HW-CN-v1","source":"%s","image_revision":"%s"}\n' \
  "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
mkdir "$out/run" "$out/work"
echo "REAL CURRENT HEAT-WAVE heat_region=2 wave_region=1"
echo "Experimental outward-normal compiler; actual L2/strict DAE/CN/energy"
echo "No corrected initial conditions; MPI1; original published source unaffected"
echo "EVIDENCE: $out"

rc=0
"${docker_cmd[@]}" run --rm \
  --user "$(id -u):$(id -g)" --env HOME=/tmp --env MPLBACKEND=Agg \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
  --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
  --env SCRIMP_V130_AUTO_INTERFACE_PILOT=1 \
  --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
  --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
  --env SCRIMP_MATRIX_MEMORY_MODE=compact \
  --env SCRIMP_MUMPS_RESET_POLICY=on_change \
  --env SCRIMP_GETFEM_RHS_MODE=auto \
  --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
  --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
  --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
  --volume "$root/sandbox/measure_v130_hw_dual_facet_physics.py:/workspace/scrimp/sandbox/measure_v130_hw_dual_facet_physics.py:ro" \
  --volume "$root/sandbox/probe_v130_swapped_interface_cn.py:/workspace/scrimp/sandbox/probe_v130_swapped_interface_cn.py:ro" \
  --volume "$out/run:/artifacts:rw" \
  --volume "$out/work:/scrimp-work:rw" --workdir /scrimp-work \
  "$image" timeout -k 5s "${SCRIMP_V130_SWAPPED_CN_TIMEOUT:-360s}" \
  mpiexec --oversubscribe --bind-to none -n 1 \
  python -u /workspace/scrimp/sandbox/probe_v130_swapped_interface_cn.py \
    --out /artifacts \
  > "$out/run.host.log" 2>&1 || rc=$?
printf '%s\n' "$rc" > "$out/run.exitcode"
if [[ "$rc" -ne 0 ]]; then
  echo "Swapped strict DAE/CN/energy native gate failed, rc=$rc" >&2
  tail -n 100 "$out/run.host.log" >&2
fi

python3 - "$out" <<'PY'
import json, sys, zipfile
from pathlib import Path
root = Path(sys.argv[1])
exit_code = int((root / "run.exitcode").read_text())
report_path = root / "run" / "swapped_auto_cn_audit.json"
reported = json.loads(report_path.read_text()) if report_path.is_file() else {}
accepted = (
    exit_code == 0
    and reported.get("status") == "PASS_SWAPPED_AUTO_INTERFACE_STRICT_DAE_CN_ENERGY_MPI1"
    and reported.get("heat_bulk") == 2
    and reported.get("wave_bulk") == 1
    and reported.get("actual_4_step_CN")
    and all(reported.get("physical_gates", {}).values())
)
summary = {
    "schema": "scrimp-v130-swapped-real-HW-CN-check-v1",
    "status": (
        "PASS_SWAPPED_AUTO_INTERFACE_DAE_CN_ENERGY_MPI1"
        if accepted else "FAIL_OR_UNQUALIFIED_SWAPPED_AUTO_INTERFACE_CN"
    ),
    "native_exit_code": exit_code,
    "raw_original_initial_condition_used": True,
    "original_live_example_not_rewritten_at_runtime": True,
    "published_archival_sources_unchanged": True,
    "physical_report": reported if report_path.is_file() else None,
    "MPI2_qualified": False,
}
(root / "comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
with zipfile.ZipFile(str(root) + ".zip", "x", compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4, allowZip64=True) as z:
    for item in sorted(root.rglob("*")):
        if item.is_file():
            z.write(item, arcname=str(item.relative_to(root)))
print(json.dumps({
    "status": summary["status"], "native_exit_code": exit_code,
    "raw_evidence": str(root), "archive": str(root) + ".zip",
}, indent=2))
if not accepted:
    raise SystemExit(1)
PY
echo "COMPLETE swapped native auto interface CN: $out"
