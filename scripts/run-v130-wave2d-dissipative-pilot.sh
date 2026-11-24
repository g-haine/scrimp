#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE.
#
# Real 2B: anisotropic dissipative wave on the same 2A rectangle,
# identical 2A 4-boundary excitation, new CG2 f_r/e_r volume port.
# MUMPS MPI1 and MPI2, 10 genuine CN steps; signed work MINUS
# strictly positive resistive loss; two independent GetFEM quadratures.
# This is NOT a SuperLU_DIST release qualification or a performance run.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == "1" ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
  echo "Wrong pinned original GetFEM/PETSc image: $revision" >&2
  exit 2
}
python3 -m py_compile \
  examples/wave.py \
  examples/wave_dissipative.py \
  tests/test_v130_wave2d_dissipative_source.py
python tests/test_v130_wave2d_dissipative_source.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-wave2d-dissipative-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-dissipative-wave2d-native-physical-pilot-v1","source":"%s","image":"%s","image_revision":"%s","original_model":"examples/wave_dissipative.py","base_model":"examples/wave.py","solver":"mumps","MPI":[1,2],"requested_CN_steps":10,"damping_strength":1.0,"not_superlu_dist_qualification":true}\n' \
  "$(git rev-parse HEAD)" "$image" "$revision" > "$out/manifest.json"
echo "2B REAL DISSIPATIVE ANISOTROPIC 2D WAVE — strict pH CN audit"
echo "Same 2A tensor, density, left velocity and 3 other boundary ports"
echo "New resistive port: f_r=-e_p, e_r=nu*f_r, nu=0.5*(2-x)>=0"
echo "Validate H(t)-H(0)=W_four_boundaries(t)-integral_0^t D(s)ds"
echo "Both native Damping power and independent GetFEM nu*f_r^2 at CN mean"
echo "MUMPS MPI1/MPI2; 10 CN steps; 6 original FE variables x 11 VTU"
echo "EVIDENCE: $out"
for ranks in 1 2; do
  name="mumps_mpi${ranks}"
  echo "RUN $name"
  mkdir -p "$out/$name" "$out/work_$name"
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
    --volume "$root/examples/wave_dissipative.py:/workspace/scrimp/examples/wave_dissipative.py:ro" \
    --volume "$out/$name:/artifacts:rw" \
    --volume "$out/work_$name:/scrimp-work:rw" \
    --workdir /scrimp-work \
    "$image" timeout -k 5s "${SCRIMP_V130_WAVE2D_DISSIPATIVE_TIMEOUT:-180s}" \
    mpiexec --oversubscribe --bind-to none -n "$ranks" \
    python -u -m examples.wave_dissipative \
      --solver mumps --mesh-size 0.2 \
      --final-time 0.08 --time-step 0.008 --damping-strength 1 \
      --check --export --output /artifacts \
    > "$out/$name.host.log" 2>&1 || rc=$?
  echo "$rc" > "$out/$name.exitcode"
  if [[ "$rc" != 0 ]]; then
    echo "UNQUALIFIED $name exit=$rc; real original logs saved" >&2
    tail -n 35 "$out/$name.host.log" >&2
  fi
done
python3 - "$out" <<'PY'
import json
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

folder=Path(sys.argv[1])
report={
    "schema":"v130-dissipative-wave2d-MUMPS-MPI-physical-v1",
    "status":"FAIL_OR_UNQUALIFIED_NATIVE_2B_WAVE_DISSIPATION",
    "cases":{}, "MPI1_vs_MPI2":{}, "errors":[],
    "superlu_dist_qualified":False,
    "canonical_2A_source_unchanged":True,
    "published_sources_unchanged":True,
}
keys=("t","z","H","power_bottom","power_right","power_top",
      "power_left","signed_work","cumulative_loss","dissipated_power",
      "independently_assembled_volume_loss","resistive_flow_constraint",
      "resistive_effort_constraint","energy_step_defects")
fields=("q","p","e_q","e_p","f_r","e_r")
cases=[(1,"mumps_mpi1"),(2,"mumps_mpi2")]
for ranks,name in cases:
    report["cases"][name]={"exit_code":int((folder/f"{name}.exitcode").read_text()),
                           "MPI":ranks}
try:
    states={}
    for ranks,name in cases:
        record=report["cases"][name]
        if record["exit_code"]:
            raise RuntimeError(f"Native 2B {name} failed: {record['exit_code']}")
        work=folder/name
        values=[]
        for rank in range(ranks):
            meta=json.loads((work/f"rank{rank}.json").read_text())
            with np.load(work/f"rank{rank}.npz",allow_pickle=False) as inp:
                data={key:inp[key].copy() for key in inp.files}
            gates=meta.get("physical_gates",{})
            if (meta.get("status")!="PASS_NATIVE_2B_DISSIPATIVE_FE_PHYSICS"
                or meta.get("rank")!=rank or meta.get("MPI")!=ranks
                or meta.get("requested_LU")!="mumps"
                or meta.get("selected_LU")!="mumps"
                or meta.get("steps")!=10
                or meta.get("DOFs",0)<=2230
                or meta.get("damping_strength")!=1.0
                or meta.get("TS_reason",0)<=0
                or meta.get("DAE_initialization",{}).get("status")!="PASS"
                or not meta.get("DAE_initialization",{}).get("differential_initial_exact")
                or len(gates)!=10 or not all(gates.values())
                or any(key not in data for key in keys)):
                raise RuntimeError(f"Incomplete physical 2B MUMPS evidence: {name}/r{rank}")
            t=data["t"]; z=data["z"]; H=data["H"]
            if (t.shape!=(11,) or z.shape!=(11,meta["DOFs"])
                or H.shape!=t.shape
                or not np.allclose(t,np.arange(11)*.008,atol=1e-10,rtol=0)
                or any(not np.all(np.isfinite(v)) for v in data.values())
                or any(data[k].shape!=(11,) for k in keys if k not in ("z","energy_step_defects"))
                or data["energy_step_defects"].shape!=(10,)):
                raise RuntimeError(f"Invalid original 2B CN/FE trajectory: {name}/r{rank}")
            P=data["dissipated_power"]; Pv=data["independently_assembled_volume_loss"]
            loss=data["cumulative_loss"]; W=data["signed_work"]
            defects=np.diff(H)-np.diff(t)*(
                data["power_bottom"][1:]+data["power_right"][1:]
                +data["power_top"][1:]+data["power_left"][1:]) + np.diff(t)*P[1:]
            if (np.min(P[1:]) < -1e-11 or np.max(P[1:])<=1e-10
                or loss[-1] <= 1e-12
                or np.min(np.diff(loss)) < -1e-12
                or not np.allclose(np.diff(loss),np.diff(t)*P[1:],atol=1e-11,rtol=1e-7)
                or np.max(np.abs(Pv-P)) >= 1e-9
                or not np.allclose(defects,data["energy_step_defects"],atol=1e-10,rtol=0)
                or np.max(np.abs(defects)) >= 1e-8
                or np.max(np.abs(H-H[0]-W+loss)) >= 1e-8
                or abs(float(loss[-1])-float(meta["total_dissipated_energy"])) >= 1e-10):
                raise RuntimeError(f"True H=W-D or independent damping quadrature failed: {name}/r{rank}")
            values.append(data)
        reference=values[0]
        for rank,other in enumerate(values[1:],start=1):
            for key in keys:
                if (other[key].shape!=reference[key].shape
                    or not np.allclose(other[key],reference[key],rtol=2e-8,atol=5e-9)):
                    raise RuntimeError(f"Original MPI rank field mismatch: {name}/r{rank}/{key}")
        VTK={}
        for field in fields:
            path=work/"pv"/field/f"{field}.pvd"
            frames=ET.parse(path).findall(".//DataSet")
            if len(frames)!=11:
                raise RuntimeError(f"Missing original GetFEM {field} frames in {name}")
            grids=[]
            for i,frame in enumerate(frames):
                if abs(float(frame.attrib["timestep"])-reference["t"][i]) > 1e-10:
                    raise RuntimeError(f"Invalid native time in {path}")
                filename=Path(frame.attrib["file"])
                if filename.name!=str(filename):
                    raise RuntimeError(f"Unsafe relative VTK link {filename}")
                piece=ET.fromstring((path.parent/filename).read_bytes()).find(".//Piece")
                if (piece is None or int(piece.attrib["NumberOfPoints"])<2
                    or int(piece.attrib["NumberOfCells"])<1
                    or piece.find(f"./PointData/DataArray[@Name='{field}']") is None):
                    raise RuntimeError(f"Missing native FE field {field} in {name}")
                grids.append((int(piece.attrib["NumberOfPoints"]),
                              int(piece.attrib["NumberOfCells"])))
            if len(set(grids))!=1:
                raise RuntimeError(f"Nonstationary VTK geometry in {field}")
            VTK[field]=len(frames)
        record.update({
            "status":"PASS_NATIVE_DISSIPATIVE_WAVE2D_MUMPS_CN_POWER_VTK",
            "DOFs":int(reference["z"].shape[1]),
            "steps":10,
            "max_abs_step_defect":float(np.max(np.abs(reference["energy_step_defects"]))),
            "max_abs_cumulative_defect":
                float(np.max(np.abs(reference["H"]-reference["H"][0]
                                    -reference["signed_work"]+reference["cumulative_loss"]))),
            "max_abs_dissipation_quadrature_difference":
                float(np.max(np.abs(reference["dissipated_power"]-
                                    reference["independently_assembled_volume_loss"]))),
            "total_viscous_loss":float(reference["cumulative_loss"][-1]),
            "max_abs_left_power":float(np.max(np.abs(reference["power_left"]))),
            "rank_traces_agree":True,
            "VTK_fields":VTK,
        })
        states[name]=reference
    for key in keys:
        left,right=states["mumps_mpi1"][key],states["mumps_mpi2"][key]
        delta=float(np.max(np.abs(left-right)))
        report["MPI1_vs_MPI2"][key]=delta
        if left.shape!=right.shape or not np.allclose(left,right,atol=5e-9,rtol=2e-8):
            raise RuntimeError(f"Different original 2B FE trajectories/powers: {key}")
    report["status"]="PASS_WAVE2D_DISSIPATIVE_MUMPS_MPI1_MPI2_PHYSICAL_VTK"
except Exception as exc:
    report["errors"].append(f"{type(exc).__name__}: {exc}")
for ranks,name in cases:
    record=report["cases"][name]
    if record["exit_code"]!=0:
        message=f"Native {name} exit={record['exit_code']}"
        if message not in report["errors"]:report["errors"].append(message)
(folder/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
archive=Path(str(folder)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for path in sorted(folder.rglob("*")):
        if path.is_file():z.write(path,arcname=str(path.relative_to(folder)))
print(json.dumps({**report,"archive":str(archive)},indent=2),flush=True)
if report["status"]!="PASS_WAVE2D_DISSIPATIVE_MUMPS_MPI1_MPI2_PHYSICAL_VTK":
    raise SystemExit(1)
PY
echo "COMPLETE native 2B dissipative original GetFEM/PETSc pilot: $out"
