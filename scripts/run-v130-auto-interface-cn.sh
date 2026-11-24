#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# Genuine original Heat-Wave strict DAE / TSCN with manually assigned
# reference facets versus the experimental automatic interface compiler.
# Neither physical example source nor physical weak form is edited.
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
 echo "Pinned GetFEM/PETSc image required, got: $revision" >&2; exit 2
}
python3 -m py_compile scrimp/dphs.py scrimp/port.py scrimp/brick.py \
 scrimp/interface.py scrimp/interface_deferred.py scrimp/interface_discovery.py \
 sandbox/probe_v130_auto_interface_cn.py \
 sandbox/probe_v130_hw_dual_sided_cn.py \
 sandbox/measure_v130_hw_dual_facet_physics.py
python3 tests/test_v130_interface_topology.py
python3 tests/test_v130_interface_ownership.py
python3 tests/test_v130_interface_assembly_plan.py
python3 tests/test_v130_interface_deferred.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-auto-interface-CN-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-auto-interface-real-CN-comparison-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "EXPERIMENTAL original Heat-Wave automatic vs reference four-step CN"
echo "Strict DAE, two boundary-compatible initializations, MUMPS, MPI1"
echo "Ports use dynamically resolved interface faces (never doubled region 10)"
echo "EVIDENCE: $out"

for candidate in clamp constrained_l2; do
 for variant in reference auto; do
  mkdir "$out/${candidate}_${variant}" "$out/work_${candidate}_${variant}"
  if [[ "$variant" == reference ]]; then
   entry=/workspace/scrimp/sandbox/probe_v130_hw_dual_sided_cn.py
   toggle=0
  else
   entry=/workspace/scrimp/sandbox/probe_v130_auto_interface_cn.py
   toggle=1
  fi
  echo "RUN candidate=$candidate variant=$variant"
  rc=0
  "${docker_cmd[@]}" run --rm --user "$(id -u):$(id -g)" \
   --env HOME=/tmp --env MPLBACKEND=Agg \
   --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
   --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
   --env SCRIMP_V130_AUTO_INTERFACE_PILOT="$toggle" \
   --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
   --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
   --env SCRIMP_MATRIX_MEMORY_MODE=compact \
   --env SCRIMP_MUMPS_RESET_POLICY=on_change \
   --env SCRIMP_GETFEM_RHS_MODE=auto \
   --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
   --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
   --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
   --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
   --volume "$root/sandbox/probe_v130_main_2d_examples.py:/workspace/scrimp/sandbox/probe_v130_main_2d_examples.py:ro" \
   --volume "$root/sandbox/probe_v130_heat_wave_boundary_compatibility.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_boundary_compatibility.py:ro" \
   --volume "$root/sandbox/probe_v130_heat_wave_compatible_cn.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_compatible_cn.py:ro" \
   --volume "$root/sandbox/measure_v130_heat_wave_physics.py:/workspace/scrimp/sandbox/measure_v130_heat_wave_physics.py:ro" \
   --volume "$root/sandbox/measure_v130_hw_dual_facet_physics.py:/workspace/scrimp/sandbox/measure_v130_hw_dual_facet_physics.py:ro" \
   --volume "$root/sandbox/v130_hw_interface_geometry.py:/workspace/scrimp/sandbox/v130_hw_interface_geometry.py:ro" \
   --volume "$root/sandbox/probe_v130_hw_dual_sided_cn.py:/workspace/scrimp/sandbox/probe_v130_hw_dual_sided_cn.py:ro" \
   --volume "$root/sandbox/probe_v130_auto_interface_cn.py:/workspace/scrimp/sandbox/probe_v130_auto_interface_cn.py:ro" \
   --volume "$out/${candidate}_${variant}:/artifacts:rw" \
   --volume "$out/work_${candidate}_${variant}:/scrimp-work:rw" \
   --workdir /scrimp-work \
   "$image" timeout -k 5s "${SCRIMP_V130_AUTO_CN_TIMEOUT:-360s}" \
   mpiexec --oversubscribe --bind-to none -n 1 \
   python -u "$entry" --candidate "$candidate" --out /artifacts \
   > "$out/${candidate}_${variant}.host.log" 2>&1 || rc=$?
  printf '%s\n' "$rc" > "$out/${candidate}_${variant}.exitcode"
  if [[ "$rc" -ne 0 ]]; then
   echo "FAIL real CN candidate=$candidate variant=$variant exit=$rc" >&2
   tail -n 100 "$out/${candidate}_${variant}.host.log" >&2
  fi
 done
done

python3 - "$out" <<'PY'
import json, sys, zipfile
from pathlib import Path

import numpy as np

root = Path(sys.argv[1])
results = {}
for candidate in ("clamp", "constrained_l2"):
    base = root / (candidate + "_reference")
    auto = root / (candidate + "_auto")
    rc_ref = int((root / (candidate + "_reference.exitcode")).read_text())
    rc_auto = int((root / (candidate + "_auto.exitcode")).read_text())
    entry = {"reference_exit": rc_ref, "auto_exit": rc_auto}
    results[candidate] = entry
    if rc_ref or rc_auto:
        entry["status"] = "FAIL_REAL_NATIVE_CN_EXECUTION"
        continue
    ref_report = json.loads((base / "dual_facet_cn_audit.json").read_text())
    auto_report = json.loads((auto / "auto_interface_cn_audit.json").read_text())
    ref_physics = ref_report["actual_energy_power_diagnostics"]
    auto_physics = auto_report["actual_physical_diagnostics"]
    ref_trace = np.load(base / "original_2d_example_trace.npz")
    auto_trace = np.load(auto / "original_2d_example_trace.npz")

    keys = ("times", "states", "H", "T", "p", "q")
    differences = {}
    for key in keys:
        a, b = np.asarray(ref_trace[key]), np.asarray(auto_trace[key])
        if a.shape != b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            raise RuntimeError(f"Invalid full trace {candidate}:{key}")
        differences[key] = {
            "shape": list(a.shape),
            "max_abs": float(np.max(np.abs(a - b))),
            "equal_with_tolerance": bool(np.allclose(a, b, rtol=2e-9, atol=2e-10)),
        }
    physical_keys = (
        "max_CN_step_H_balance_defect",
        "max_CN_step_heat_balance_defect",
        "max_CN_step_wave_balance_defect",
        "max_CN_interface_power_defect",
        "max_individual_CN_interface_power",
        "max_CN_heat_dissipation",
        "interface_exchange_ratio_to_heat_dissipation",
    )
    physics_difference = {}
    for key in physical_keys:
        a, b = float(ref_physics[key]), float(auto_physics[key])
        physics_difference[key] = {
            "reference": a, "automatic": b,
            "abs_difference": abs(a - b),
            "close": bool(np.isclose(a, b, rtol=2e-8, atol=2e-10)),
        }

    # Separate accepted-step/numerical parity from the *actual*
    # energy and nontrivial-interface physical gates.
    numerical = all(v["equal_with_tolerance"] for v in differences.values())
    diagnostics = all(v["close"] for v in physics_difference.values())
    science = bool(
        ref_report["strict_DAE_initialization_PASSED"]
        and ref_report["actual_real_4_step_PETSc_CN_PASSED"]
        and ref_report["CN_discrete_midpoint_energy_gate_PASSED"]
        and ref_report["CN_interface_exchange_nontrivial"]
        and auto_report["actual_true_strict_DAE_initialization_passed"]
        and auto_report["actual_real_4_step_PETSc_CN_passed"]
        and auto_report["discrete_CN_energy_gate_passed"]
        and auto_report["nontrivial_interface_exchange_passed"]
        and auto_physics["old_doubled_interface_region_excluded_from_physical_power"]
        and auto_physics["interface_heat_region"] != 10
        and auto_physics["interface_wave_region"] != 10
    )
    status = (
        "PASS_AUTO_INTERFACE_REAL_CN_DAE_ENERGY_REFERENCE_PARITY"
        if numerical and diagnostics and science
        else "FAIL_OR_UNQUALIFIED_AUTO_INTERFACE_REAL_CN_PARITY"
    )
    entry.update({
        "status": status,
        "trajectory": differences,
        "physical_diagnostics": physics_difference,
        "actual_DAE_and_4_step_CN_and_energy_gate": science,
        "original_example_source_unchanged": True,
        "no_manual_constructor_remapping_in_auto_arm": True,
        "auto_resolved_regions": [
            auto_report["automatic_interface_manifest"]["one_sided_heat_region"],
            auto_report["automatic_interface_manifest"]["one_sided_wave_region"],
        ],
    })

all_pass = all(
    value["status"] == "PASS_AUTO_INTERFACE_REAL_CN_DAE_ENERGY_REFERENCE_PARITY"
    for value in results.values()
)
summary = {
    "schema": "scrimp-v130-auto-interface-real-CN-physical-parity-v1",
    "status": (
        "PASS_BOTH_AUTO_INTERFACE_CANDIDATES_DAE_CN_ENERGY_PARITY"
        if all_pass else "FAIL_OR_UNQUALIFIED_AUTO_INTERFACE_DAE_CN_PHYSICS"
    ),
    "cases": results,
    "production_default_auto_interface_enabled": False,
    "MPI2_qualified": False,
    "swapped_subdomain_physical_sign_convention_qualified": False,
    "original_published_example_unchanged": True,
}
(root / "comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
archive = Path(str(root) + ".zip")
with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4, allowZip64=True) as z:
    for file in sorted(root.rglob("*")):
        if file.is_file():
            z.write(file, arcname=str(file.relative_to(root)))
print(json.dumps(summary, indent=2))
print("SHAREABLE ORIGINAL AUTO INTERFACE CN EVIDENCE:", archive)
if not all_pass:
    raise SystemExit(1)
PY

echo "COMPLETE original Heat-Wave automatic CN vs reference: $out"
