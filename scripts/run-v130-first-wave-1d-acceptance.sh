#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# FIRST DOCUMENTED EXAMPLE — true GetFEM/PETSc source-install gate.
# One real port-Hamiltonian string model, five supported LU/MPI cases.
# Fail closed on physical power, differential DAE, actual accepted CN,
# complete FE/port parity, geometry and optional q/p VTK export.
# This is NOT a solver performance benchmark or published paper replay.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
expected="1d6d0c90115213298047b28237a7f38087563a35"
[[ "$revision" == "$expected" ]] || {
  echo "Expected qualified immutable GetFEM/PETSc image $expected, found $revision" >&2
  exit 2
}

python3 -m py_compile examples/wave_1d.py tests/test_v130_wave_1d_tutorial.py
python3 tests/test_v130_wave_1d_tutorial.py
python3 scripts/check-v130-published-immutable.py
bash -n "$0"

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-first-tutorial-wave1d-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-wave1d-install-acceptance-v1","source":"%s","pinned_image":"%s","image_revision":"%s","native_status":"PENDING","mathematical_model":"1D PFEM conservative string, outward Neumann forces","publication_connection":"Haine-Matignon-Serhani 2023; Ferraro-Fournie-Haine 2024; NOT paper numerical reproduction","solver_cases":["mumps_mpi1","mumps_mpi2","superlu_mpi1","superlu_dist_mpi1","superlu_dist_mpi2"],"superlu_MPI2_deliberately_unsupported":true}\n' \
  "$(git rev-parse HEAD)" "$image" "$revision" > "$out/manifest.json"

echo "FIRST USER TUTORIAL: original 1D vibrating string, compatible initial q,p"
echo "Five genuine GetFEM/PETSc direct LU cases: MUMPS MPI1/2; SuperLU MPI1; SuperLU_DIST MPI1/2"
echo "Strict DAE, ten accepted fixed CN steps, signed boundary work and 1D GetFEM ParaView"
echo "Host-owned immutable evidence: $out"
echo "No published historical sources/outputs changed; no PETSc speed claims"

cases=(mumps:1 mumps:2 superlu:1 superlu_dist:1 superlu_dist:2)
for item in "${cases[@]}"; do
  solver="${item%%:*}"
  ranks="${item##*:}"
  name="${solver}_mpi${ranks}"
  mkdir -p "$out/$name" "$out/work_$name"
  echo "RUN $name"
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
    --volume "$root/examples/wave_1d.py:/workspace/scrimp/examples/wave_1d.py:ro" \
    --volume "$out/$name:/artifacts:rw" \
    --volume "$out/work_$name:/scrimp-work:rw" \
    --workdir /scrimp-work \
    "$image" timeout -k 5s "${SCRIMP_V130_WAVE1D_TIMEOUT:-150s}" \
    mpiexec --oversubscribe --bind-to none -n "$ranks" \
    python -u -m examples.wave_1d --solver "$solver" \
      --mesh-size 0.05 --time-step 0.004 --final-time 0.04 \
      --check --export --output /artifacts \
    > "$out/$name.host.log" 2>&1 || rc=$?
  printf '%s\n' "$rc" > "$out/$name.exitcode"
  if [[ "$rc" != 0 ]]; then
    echo "NOT QUALIFIED original first tutorial $name exit=$rc; logging retained" >&2
    tail -n 55 "$out/$name.host.log" >&2
  fi
done

python3 - "$out" <<'PY'
import json, sys, zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np

root = Path(sys.argv[1])
cases = [("mumps",1),("mumps",2),("superlu",1),
         ("superlu_dist",1),("superlu_dist",2)]
report={"schema":"scrimp-v130-first-tutorial-real-native-acceptance-v1",
        "status":"FAIL_OR_UNQUALIFIED_WAVE_1D_BACKENDS_PHYSICS_MPI_VTK",
        "cases":{}, "comparisons":{},
        "original_published_snapshots_unmodified":True,
        "model_is_actual_1D_PFEM":True,
        "speed_benchmark":False}
errors=[]
traces={}
try:
    for backend,ranks in cases:
        label=f"{backend}_mpi{ranks}"
        folder=root/label
        item={"exit_code":int((root/f"{label}.exitcode").read_text())}
        report["cases"][label]=item
        if item["exit_code"]!=0:
            errors.append(f"{label} native exit={item['exit_code']}")
            continue
        all_ranks=[]
        for rank in range(ranks):
            info=json.loads((folder/f"rank{rank}.json").read_text())
            with np.load(folder/f"rank{rank}.npz",
                         allow_pickle=False) as loaded:
                arrays={k:np.asarray(loaded[k]).copy() for k in loaded.files}
            if (info["rank"]!=rank or info["MPI"]!=ranks
                    or info["requested_solver"]!=backend
                    or info["selected_PETSc_option"]!=backend
                    or info["accepted_steps"]!=10
                    or not all(info["physical_gates"].values())
                    or "strict_DAE_init_without_changing_dynamic_data"
                       not in info["physical_gates"]):
                raise RuntimeError(f"Solver/DAE/CN/physics rank gate: {label}/{rank}")
            t,z,H=arrays["times"],arrays["states"],arrays["H"]
            if (t.shape!=(11,) or z.shape[0]!=11
                    or H.shape!=(11,) or z.shape[1]<=10
                    or not np.all(np.isfinite(z))
                    or not np.all(np.isfinite(H))
                    or not np.allclose(t,np.arange(11)*.004,rtol=0,atol=1e-11)):
                raise RuntimeError(f"Wrong accepted FE trajectory {label}/{rank}")
            if (info["max_abs_step_energy_defect"]>1e-8
                    or info["max_abs_boundary_power"]<=1e-9):
                raise RuntimeError(f"Missing nontrivial pH power balance: {label}")
            all_ranks.append((info,arrays))
        ref_info,ref=all_ranks[0]
        for rank,(info,trace) in enumerate(all_ranks[1:],start=1):
            for key in ref:
                if (trace[key].shape!=ref[key].shape
                        or not np.allclose(trace[key],ref[key],
                                           rtol=2e-9,atol=1e-10)):
                    raise RuntimeError(f"MPI rank full-data disagreement {label}/{rank}/{key}")
        visual={}
        for field in ("q","p"):
            pvd=folder/"pv"/field/f"{field}.pvd"
            entries=ET.parse(pvd).findall(".//DataSet")
            if len(entries)!=11:
                raise RuntimeError(f"Expected 11 genuine VTK frames in {pvd}")
            valid=[]
            for k,entry in enumerate(entries):
                tick=float(entry.attrib["timestep"])
                if abs(tick-k*.004)>1e-11:
                    raise RuntimeError(f"Wrong VTK time in {pvd}")
                filename=entry.attrib["file"]
                if Path(filename).name!=filename:
                    raise RuntimeError(f"Unsafe VTK filename in {pvd}")
                vtupath=pvd.parent/filename
                payload=vtupath.read_bytes()
                piece=ET.fromstring(payload).find(".//Piece")
                if piece is None or not payload:
                    raise RuntimeError(f"Invalid native GetFEM VTK output {vtupath}")
                valid.append({
                    "t":tick, "bytes":len(payload),
                    "points":int(piece.attrib["NumberOfPoints"]),
                    "cells":int(piece.attrib["NumberOfCells"]),
                })
            visual[field]=valid
        item.update({
            "status":"PASS_NATIVE_REAL_DAE_CN_ENERGY_AND_VTK",
            "DOFs":int(ref["states"].shape[1]),
            "max_abs_energy_defect":
                float(ref_info["max_abs_step_energy_defect"]),
            "max_abs_boundary_power":
                float(ref_info["max_abs_boundary_power"]),
            "MPI_rank_traces_agree":True, "VTK":visual,
        })
        traces[label]=ref
    if not errors and len(traces)==len(cases):
        anchor=traces["mumps_mpi1"]
        for label,run in traces.items():
            differences={}
            for key in ("times","states","H","power_left","power_right",
                        "signed_work","energy_defect"):
                a,b=anchor[key],run[key]
                if a.shape!=b.shape:
                    raise RuntimeError(f"Different source 1D FE model for {label}: {key}")
                delta=float(np.max(np.abs(a-b)))
                if not np.isfinite(delta) or not np.allclose(
                    a,b,rtol=2e-8,atol=5e-10
                ):
                    raise RuntimeError(f"MUMPS serial vs {label}/{key} FAIL: {delta}")
                differences[key]=delta
            report["comparisons"][label]=differences
        report["status"]="PASS_FIRST_WAVE_1D_ALL_AVAILABLE_LU_MPI_PHYSICS_VTK"
except Exception as exc:
    errors.append(f"{type(exc).__name__}: {exc}")

report["errors"]=errors
(root/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as zipped:
    for path in sorted(root.rglob("*")):
        if path.is_file():
            zipped.write(path,arcname=str(path.relative_to(root)))
print(json.dumps({
    "status":report["status"],
    "cases":{k:v.get("status","FAIL_OR_UNQUALIFIED")
             for k,v in report["cases"].items()},
    "errors":errors, "archive":str(archive),
},indent=2))
if report["status"]!="PASS_FIRST_WAVE_1D_ALL_AVAILABLE_LU_MPI_PHYSICS_VTK":
    raise SystemExit(1)
PY
echo "COMPLETE first-user-tutorial native acceptance: $out"
