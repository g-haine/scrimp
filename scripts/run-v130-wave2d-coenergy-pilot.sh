#!/usr/bin/env bash
# SCRIMP v1.3.0: original 2C coenergy, native MUMPS MPI1/2 physical pilot.
# Real GetFEM CN midpoint energy, independent inner damping and rim feedback.
# All outputs remain in sandbox/outputs (never examples/published).
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
  echo "Unexpected pinned GetFEM/PETSc native image: $revision" >&2; exit 2;
}
python3 -m py_compile examples/wave_coenergy.py \
  tests/test_v130_wave2d_coenergy_source.py \
  scripts/render-v130-wave2d-native.py
bash -n "$0"
python3 tests/test_v130_wave2d_coenergy_source.py
python3 tests/test_v130_interface_deferred.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-wave2d-coenergy-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-coenergy-wave2d-MUMPS-MPI-pilot-v1","source":"%s","image":"%s","image_revision":"%s","canonical":"examples/wave_coenergy.py","physical_scenario":"true compliance T inverse, two-region disk, positive inner damping and negative outer boundary feedback","steps":10,"not_superlu_dist_qualification":true}\n' \
  "$(git rev-parse HEAD)" "$image" "$revision" > "$out/manifest.json"
echo "COENERGY 2C: true T^-1, internal disk resistance and outer passive feedback"
echo "10 authentic CN steps, MUMPS MPI1/2, FE q/p/e_r, physical H+Li+Lb"
echo "EVIDENCE: $out"
for ranks in 1 2; do
  name="mumps_mpi$ranks"
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
    --volume "$root/examples/wave_coenergy.py:/workspace/scrimp/examples/wave_coenergy.py:ro" \
    --volume "$out/$name:/artifacts:rw" \
    --volume "$out/work_$name:/scrimp-work:rw" \
    --workdir /scrimp-work \
    "$image" timeout -k 5s "${SCRIMP_V130_WAVE2D_COENERGY_TIMEOUT:-180s}" \
    mpiexec --oversubscribe --bind-to none -n "$ranks" \
    python -u -m examples.wave_coenergy \
      --solver mumps --mesh-size 0.2 --final-time 0.08 --time-step 0.008 \
      --check --export --output /artifacts \
    >"$out/$name.host.log" 2>&1 || rc=$?
  echo "$rc" >"$out/$name.exitcode"
  if [[ "$rc" != 0 ]]; then
    echo "UNQUALIFIED $name exit=$rc; saved native traces" >&2
    tail -n 35 "$out/$name.host.log" >&2
  fi
done

# Do not hide native failures. Preserve evidence and zip even on failure.
python3 - "$out" <<'PY'
import json
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

folder=Path(sys.argv[1])
names=("mumps_mpi1","mumps_mpi2")
fields=("q","p","e_r")
required=("t","z","H","power_boundary","dissipation_internal",
          "dissipation_internal_weak_check",
          "dissipation_internal_pointwise_diagnostic",
          "dissipation_feedback","feedback_constraint",
          "resistive_weak_residual_L2","cumulative_internal_loss",
          "cumulative_feedback_loss","energy_step_defects",
          "energy_cumulative_defects")
summary={
    "schema":"v130-wave2d-coenergy-real-MUMPS-MPI-v1",
    "status":"FAIL_OR_UNQUALIFIED_2C_COENERGY_NATIVE_PHYSICS",
    "cases":{}, "MPI1_vs_MPI2":{}, "errors":[],
    "superlu_dist_qualified":False,
    "published_sources_unchanged":True,
}
for ranks,name in enumerate(names,start=1):
    summary["cases"][name]={"MPI":ranks,
        "exit_code":int((folder/f"{name}.exitcode").read_text())}
try:
    trajectories={}
    for ranks,name in enumerate(names,start=1):
        rec=summary["cases"][name]
        if rec["exit_code"]!=0:
            raise RuntimeError(f"{name}: native exit {rec['exit_code']}")
        traces=[]
        work=folder/name
        for rank in range(ranks):
            meta=json.loads((work/f"rank{rank}.json").read_text())
            with np.load(work/f"rank{rank}.npz",allow_pickle=False) as npz:
                data={k:npz[k].copy() for k in npz.files}
            gates=meta.get("physical_gates",{})
            if (meta.get("status")!="PASS_NATIVE_WAVE2D_COENERGY_PHYSICS"
                or meta.get("MPI")!=ranks or meta.get("rank")!=rank
                or meta.get("steps")!=10 or meta.get("TS_reason",0)<=0
                or meta.get("selected_LU")!="mumps"
                or meta.get("DAE_initialization",{}).get("status")!="PASS"
                or len(gates)!=10 or not all(gates.values())
                or any(key not in data for key in required)):
                raise RuntimeError(f"{name}/rank{rank} native physics gate FAIL")
            t=data["t"];z=data["z"];H=data["H"]
            if (t.shape!=(11,) or z.shape!=(11,meta["DOFs"])
                or H.shape!=t.shape
                or not np.allclose(t,np.arange(11)*.008,atol=1e-10,rtol=0)
                or any(not np.all(np.isfinite(x)) for x in data.values())
                or any(data[k].shape!=(11,) for k in required
                       if k not in ("z","energy_step_defects"))
                or data["energy_step_defects"].shape!=(10,)):
                raise RuntimeError(f"{name}/rank{rank} invalid true FE samples")
            dt=np.diff(t)
            P=data["power_boundary"]
            Di=data["dissipation_internal"]
            Li=data["cumulative_internal_loss"]
            Lb=data["cumulative_feedback_loss"]
            defect=np.diff(H)-dt*P[1:]+dt*Di[1:]
            if (np.min(Di[1:])<-1e-10 or Li[-1]<=1e-8
                or np.max(P[1:])>1e-9 or Lb[-1]<=1e-8
                or np.max(np.abs(defect))>1e-8
                or np.max(np.abs(H-H[0]+Li+Lb))>1e-8
                or not np.allclose(defect,data["energy_step_defects"],
                                   rtol=0,atol=1e-10)
                or np.max(np.abs(Di-data["dissipation_internal_weak_check"]))>1e-8
                or np.max(np.abs(data["dissipation_feedback"]+P))>1e-8):
                raise RuntimeError(f"{name}/rank{rank} independent energy FAIL")
            traces.append(data)
        first=traces[0]
        for rank,other in enumerate(traces[1:],start=1):
            for key in required:
                if (other[key].shape!=first[key].shape
                    or not np.allclose(other[key],first[key],atol=5e-9,rtol=2e-8)):
                    raise RuntimeError(f"{name}/rank{rank} MPI replica mismatch {key}")
        vtk={}
        for field in fields:
            path=work/"pv"/field/f"{field}.pvd"
            frames=ET.parse(path).findall(".//DataSet")
            if len(frames)!=11:
                raise RuntimeError(f"Missing native {field} 11 VTK frames in {name}")
            geometries=set()
            for i,frame in enumerate(frames):
                if abs(float(frame.attrib["timestep"])-first["t"][i])>1e-10:
                    raise RuntimeError(f"Bad VTK time {field}/{i}")
                filename=Path(frame.attrib["file"])
                if filename.name!=str(filename):
                    raise RuntimeError("Invalid relative VTK path")
                piece=ET.fromstring((path.parent/filename).read_bytes()).find(".//Piece")
                if (piece is None or int(piece.attrib["NumberOfPoints"])<2
                    or int(piece.attrib["NumberOfCells"])<1
                    or piece.find(f"./PointData/DataArray[@Name='{field}']") is None):
                    raise RuntimeError(f"Missing actual GetFEM field {field}/{i}")
                geometries.add((piece.attrib["NumberOfPoints"],
                                piece.attrib["NumberOfCells"]))
            if len(geometries)!=1:
                raise RuntimeError(f"VTK geometry changes in {field}")
            vtk[field]=len(frames)
        rec.update({
            "status":"PASS_NATIVE_2C_COENERGY_MUMPS_CN_POWER_VTK",
            "DOFs":int(first["z"].shape[1]),
            "steps":10,
            "total_internal_loss":float(first["cumulative_internal_loss"][-1]),
            "total_feedback_loss":float(first["cumulative_feedback_loss"][-1]),
            "max_abs_step_defect":float(np.max(np.abs(first["energy_step_defects"]))),
            "max_abs_cumulative_defect":float(np.max(np.abs(first["energy_cumulative_defects"]))),
            "rank_traces_agree":True,
            "VTK_fields":vtk,
        })
        trajectories[name]=first
    for key in required:
        a,b=trajectories[names[0]][key],trajectories[names[1]][key]
        summary["MPI1_vs_MPI2"][key]=float(np.max(np.abs(a-b)))
        if a.shape!=b.shape or not np.allclose(a,b,atol=5e-9,rtol=2e-8):
            raise RuntimeError(f"Different real FE trajectories MPI1/MPI2: {key}")
    summary["status"]="PASS_WAVE2D_COENERGY_MUMPS_MPI1_MPI2_PHYSICAL_VTK"
except Exception as exc:
    summary["errors"].append(f"{type(exc).__name__}: {exc}")

(folder/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2),flush=True)
PY

if python3 - "$out" <<'PY'
import json,sys
from pathlib import Path
s=json.loads((Path(sys.argv[1])/"summary.json").read_text())
raise SystemExit(0 if s["status"].startswith("PASS_") else 1)
PY
then
  echo "Rendering genuine accepted FE / VTK (no reconstructed physics)"
  if ! python3 scripts/render-v130-wave2d-native.py \
    --input "$out" --output "$out/visuals" --case mumps_mpi1; then
    python3 - "$out" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1])/"summary.json"
s=json.loads(p.read_text())
s["errors"].append("Native physics passed, but required visual output failed")
s["status"]="FAIL_WAVE2D_COENERGY_REQUIRED_VISUALS"
p.write_text(json.dumps(s,indent=2)+"\n")
PY
  fi
fi
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
folder=Path(sys.argv[1]); archive=Path(str(folder)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for file in sorted(folder.rglob("*")):
        if file.is_file():
            z.write(file,arcname=str(file.relative_to(folder)))
summary=json.loads((folder/"summary.json").read_text())
print(json.dumps({"status":summary["status"],"archive":str(archive),
                  "evidence":str(folder),"errors":summary["errors"]},indent=2))
if summary["status"]!="PASS_WAVE2D_COENERGY_MUMPS_MPI1_MPI2_PHYSICAL_VTK":
    raise SystemExit(1)
PY
