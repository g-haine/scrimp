#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# Automatic internal Heat-Wave, real MPI1 vs MPI2:
# strict DAE + true PETSc CN + physical sided energy/power comparison.
# Parallel ParaView and plotting intentionally excluded from this gate.
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
 echo "Pinned original image revision required: $revision" >&2
 exit 2
}

python3 -m py_compile \
 scrimp/dphs.py scrimp/port.py scrimp/interface_deferred.py \
 sandbox/probe_v130_auto_interface_mpi_cn.py \
 sandbox/measure_v130_hw_dual_facet_physics.py \
 tests/test_v130_auto_interface_mpi_cn_comparator.py
python3 tests/test_v130_interface_topology.py
python3 tests/test_v130_interface_ownership.py
python3 tests/test_v130_interface_assembly_plan.py
python3 tests/test_v130_interface_deferred.py
python3 tests/test_v130_hw_facet_collector_regions.py
python3 tests/test_v130_heat_wave_outward_signs.py
python3 tests/test_v130_interface_mpi_csr_parity.py
python3 tests/test_v130_auto_interface_mpi_cn_comparator.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-auto-interface-MPI-CN-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-MPI1-MPI2-auto-interface-DAE-CN-v1","source":"%s","image_revision":"%s"}\n' \
  "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "REAL GetFEM/PETSc strict DAE / CN / sided pH energy: MPI1 vs MPI2"
echo "Both 2D Heat-Wave orientations; original model and physical weak forms"
echo "Plot/ParaView disabled only in sandbox; native physical assembly on ALL ranks"
echo "EVIDENCE: $out"

for heat in 1 2; do
 for ranks in 1 2; do
  label="heat${heat}_mpi${ranks}"
  mkdir "$out/$label" "$out/work_$label"
  echo "RUN heat=$heat wave=$((3-heat)) MPI=$ranks initial=$([[ $heat == 1 ]] && echo clamp || echo original)"
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
    --volume "$root/sandbox/probe_v130_heat_wave_boundary_compatibility.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_boundary_compatibility.py:ro" \
    --volume "$root/sandbox/probe_v130_auto_interface_mpi_cn.py:/workspace/scrimp/sandbox/probe_v130_auto_interface_mpi_cn.py:ro" \
    --volume "$out/$label:/artifacts:rw" \
    --volume "$out/work_$label:/scrimp-work:rw" --workdir /scrimp-work \
    "$image" timeout -k 5s "${SCRIMP_V130_MPI_CN_TIMEOUT:-360s}" \
    mpiexec --oversubscribe --bind-to none -n "$ranks" \
    python -u /workspace/scrimp/sandbox/probe_v130_auto_interface_mpi_cn.py \
      --heat "$heat" --out /artifacts \
    > "$out/$label.host.log" 2>&1 || rc=$?
  printf '%s\n' "$rc" > "$out/$label.exitcode"
  if [[ "$rc" -ne 0 ]]; then
   echo "Native MPI DAE/CN/energy not qualified: heat=$heat MPI=$ranks exit=$rc" >&2
   tail -n 95 "$out/$label.host.log" >&2
  fi
 done
done

python3 - "$out" <<'PY'
import json, sys, zipfile
from pathlib import Path
import numpy as np

root = Path(sys.argv[1])
results = {}
any_failure = False
host_error = None  # success path must also serialize a complete report
try:
    for heat in (1, 2):
        run = {}
        for size in (1, 2):
            label = f"heat{heat}_mpi{size}"
            code = int((root/f"{label}.exitcode").read_text())
            record = {"exit_code": code, "rank_reports": [], "rank_traces": []}
            run[size] = record
            if code:
                record["status"] = "FAIL_NATIVE_MPI_DAE_CN_OR_PHYSICS"
                any_failure = True
                continue
            for rank in range(size):
                path = root/label
                raw = json.loads((path/f"rank{rank}.json").read_text())
                if (raw["MPI_ranks"] != size or raw["heat_region"] != heat
                        or raw["rank"] != rank
                        or raw["status"] != "PASS_RANK_MPI_CN_AND_PHYSICAL_GATES"
                        or not raw["rank_consensus"]
                        or raw["real_PETSc_CN_steps"] != 4
                        or not all(raw["physical_gates"].values())
                        or not raw["plot_and_ParaView_export_skipped_on_all_ranks"]
                        or not raw["physical_GetFEM_assembly_executed_on_all_ranks"]):
                    raise RuntimeError(f"Native rank report gate failed: {label}/{rank}")
                record["rank_reports"].append(raw)
                record["rank_traces"].append(str(path/f"rank{rank}.npz"))
            if any(
                item["rank_fingerprints"] != record["rank_reports"][0]["rank_fingerprints"]
                or item["automatic_interface"] != record["rank_reports"][0]["automatic_interface"]
                for item in record["rank_reports"]
            ):
                raise RuntimeError(f"Cross-rank automatic-interface mismatch: {label}")
            record["status"] = "PASS_NATIVE_MPI_RANK_DAE_CN_PHYSICS"
        candidate = {}
        if not any(s["exit_code"] != 0 for s in run.values()):
            reference = np.load(run[1]["rank_traces"][0])
            parallel = np.load(run[2]["rank_traces"][0])
            if set(reference.files) != set(parallel.files):
                raise RuntimeError(
                    f"Different MPI1/MPI2 physical trace fields for heat={heat}"
                )
            passed = True
            for name in reference.files:
                a = np.asarray(reference[name], dtype=float)
                b = np.asarray(parallel[name], dtype=float)
                if (a.shape != b.shape or not np.all(np.isfinite(a))
                        or not np.all(np.isfinite(b))):
                    raise RuntimeError(
                        f"Nonfinite/mismatched MPI trajectory: heat={heat}, {name}"
                    )
                max_abs = float(np.max(np.abs(a-b))) if a.size else 0.
                equal = bool(np.allclose(a, b, rtol=2e-8, atol=5e-10))
                candidate[name] = {
                    "shape": list(a.shape),
                    "max_abs_difference": max_abs,
                    "within_explicit_tolerance": equal,
                }
                passed = passed and equal
            # MPI2 rank1 and rank0 must produce exactly matching full FE
            # history and physical traces, independently of MPI1/2 tolerances.
            with np.load(run[2]["rank_traces"][1]) as second:
                with np.load(run[2]["rank_traces"][0]) as first:
                    rank_equal = all(np.array_equal(first[k], second[k])
                                     for k in first.files)
            passed = passed and rank_equal
        else:
            passed = False
            rank_equal = None
        results[f"heat{heat}_wave{3-heat}"] = {
            "status": (
                "PASS_REAL_MPI2_DAE_CN_AND_PHYSICAL_PARITY"
                if passed else "FAIL_OR_UNQUALIFIED_MPI2_DAE_CN_PARITY"
            ),
            "MPI1": {k:v for k,v in run[1].items() if k != "rank_traces"},
            "MPI2": {k:v for k,v in run[2].items() if k != "rank_traces"},
            "trace_comparison": candidate,
            "MPI2_rank0_rank1_full_trace_bitwise_identical": rank_equal,
        }
        any_failure = any_failure or not passed
    
except Exception as exc:
    any_failure = True
    host_error = f"{type(exc).__name__}: {exc}"

summary = {
    "schema": "scrimp-v130-genuine-auto-interface-MPI2-DAE-CN-physics-v1",
    "status": (
        "PASS_BOTH_AUTO_INTERFACE_ORIENTATIONS_MPI2_DAE_CN_ENERGY_PARITY"
        if not any_failure else
        "FAIL_OR_UNQUALIFIED_AUTO_INTERFACE_MPI2_DAE_CN_PHYSICS"
    ),
    "cases": results,
    "host_comparison_failure": host_error,
    "current_original_example_source_unchanged": True,
    "all_GetFEM_physical_integrals_on_all_ranks": True,
    "parallel_PNG_and_ParaView_qualified": False,
    "production_auto_compiler_enabled_by_default": False,
}
(root/"comparison.json").write_text(json.dumps(summary,indent=2)+"\n")
with zipfile.ZipFile(str(root)+".zip", "x", compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for p in sorted(root.rglob("*")):
        if p.is_file():
            z.write(p,arcname=str(p.relative_to(root)))
print(json.dumps({
    "status": summary["status"],
    "case_statuses": {k:v["status"] for k,v in results.items()},
    "archive": str(root)+".zip",
}, indent=2))
if any_failure:
    raise SystemExit(1)
PY
echo "COMPLETE actual auto-interface MPI1/MPI2 strict DAE CN/energy: $out"
