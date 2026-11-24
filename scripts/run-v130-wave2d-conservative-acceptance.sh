#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# SCIENTIFIC INSTALL ACCEPTANCE 2A — canonical 2D conservative wave.
# Five PETSc LU/MPI configurations on the same anisotropic mixed-port PFEM.
# Do NOT confuse successful LU factorization with pH energy acceptance.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
    docker_cmd=(docker)
fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
    --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
expected="1d6d0c90115213298047b28237a7f38087563a35"
[[ "$revision" == "$expected" ]] || {
    echo "Expected qualified GetFEM/PETSc image rev=$expected; got $revision" >&2
    exit 2
}

python3 -m py_compile \
    examples/wave.py \
    tests/test_v130_wave2d_conservative_pilot.py \
    tests/test_v130_interface_pv_production_policy.py
python3 tests/test_v130_wave2d_conservative_pilot.py
python3 tests/test_v130_interface_pv_production_policy.py
python3 scripts/check-v130-published-immutable.py
bash -n "$0"

mkdir -p sandbox/outputs
output="$(mktemp -d "$root/sandbox/outputs/v130-wave2d-conservative-five-lu-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-wave2d-conservative-five-lu-acceptance-v1","source":"%s","image":"%s","image_revision":"%s","case_status":"PENDING_NATIVE","canonical_example":"examples/wave.py","scenario":"genuine anisotropic membrane T(x,y), 3 zero Neumann tractions, 1 prescribed velocity on left, zero compatible initial state","not_a_published_paper_reproduction":true}\n' \
 "$(git rev-parse HEAD)" "$image" "$revision" > "$output/manifest.json"

echo "CONSERVATIVE 2D anisotropic pH wave, original q.T.q = q dot T dot q"
echo "MUMPS MPI1/2; SuperLU MPI1; SuperLU_DIST MPI1/2; 10 CN steps"
echo "Four measured boundary powers and physical H-H0-signed-work"
echo "q,p,e_q,e_p ParaView; wave2d and real MPI VTK validator regressions; published examples immutable"
echo "EVIDENCE: $output"

for item in mumps:1 mumps:2 superlu:1 superlu_dist:1 superlu_dist:2; do
    solver="${item%%:*}"
    ranks="${item##*:}"
    name="${solver}_mpi${ranks}"
    mkdir -p "$output/$name" "$output/work_$name"
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
        --volume "$root/examples/wave.py:/workspace/scrimp/examples/wave.py:ro" \
        --volume "$output/$name:/artifacts:rw" \
        --volume "$output/work_$name:/scrimp-work:rw" \
        --workdir /scrimp-work \
        "$image" timeout -k 5s "${SCRIMP_V130_WAVE2D_TIMEOUT:-150s}" \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u -m examples.wave \
            --solver "$solver" --mesh-size 0.2 \
            --final-time 0.08 --time-step 0.008 \
            --check --export --output /artifacts \
        > "$output/$name.host.log" 2>&1 || rc=$?
    printf '%s\n' "$rc" > "$output/$name.exitcode"
    if [[ "$rc" != 0 ]]; then
        echo "NOT QUALIFIED $name exit=$rc (all evidence retained)" >&2
        tail -n 65 "$output/$name.host.log" >&2
    fi
done

python3 - "$output" <<'PY'
import json
from pathlib import Path
import sys
from xml.etree import ElementTree as ET
import zipfile
import numpy as np

folder = Path(sys.argv[1])
report={
  "status": "FAIL_OR_UNQUALIFIED_WAVE2D_FIVE_LU_MPI_CN",
  "cases": {},
  "comparisons_to_mumps_mpi1": {},
  "errors": [],
  "historical_published_source_modified": False,
  "not_a_paper_reproduction": True,
}
trace_by_case={}
cases=(("mumps",1),("mumps",2),("superlu",1),
       ("superlu_dist",1),("superlu_dist",2))
# All five actual simulations have already finished before this
# reporter starts. Record EVERY exit status, even when the first
# failed case prevents qualifying the following cases physically.
for solver,ranks in cases:
    name=f"{solver}_mpi{ranks}"
    report["cases"][name]={
        "exit_code":int((folder/f"{name}.exitcode").read_text())
    }
try:
    grids={}
    for solver,ranks in cases:
        name=f"{solver}_mpi{ranks}"
        case=folder/name
        result=report["cases"][name]
        if result["exit_code"] != 0:
            raise RuntimeError(f"Native simulation failed: {name} exit={result['exit_code']}")
        rank_data=[]
        for rank in range(ranks):
            meta=json.loads((case/f"rank{rank}.json").read_text())
            with np.load(case/f"rank{rank}.npz",allow_pickle=False) as f:
                arrays={k:f[k].copy() for k in f.files}
            if (meta["status"]!="PASS_ORIGINAL_FE_PHYSICS"
                    or meta["rank"]!=rank
                    or meta["MPI"]!=ranks
                    or meta["requested_LU"]!=solver
                    or meta["selected_LU"]!=solver
                    or meta["steps"]!=10
                    or not all(meta["physical_gates"].values())):
                raise RuntimeError(f"Real physical gate failed: {name} rank{rank}")
            t=arrays["t"]; z=arrays["z"]
            if (t.shape!=(11,) or z.ndim!=2 or z.shape[0]!=len(t)
                    or z.shape[1]<100 or arrays["H"].shape!=(11,)
                    or not np.all(np.isfinite(z))
                    or not np.allclose(t,np.arange(11)*.008,
                                       atol=1e-10,rtol=0)):
                raise RuntimeError(f"Wrong original complete FE series {name} rank{rank}")
            if (meta["max_abs_step_defect"]>1e-8
                    or meta["max_abs_cumulative_defect"]>1e-8
                    or meta["max_abs_left_power"]<=1e-8):
                raise RuntimeError(f"Real port-work energy proof failed {name} rank{rank}")
            rank_data.append((meta,arrays))
        ref_meta,ref=rank_data[0]
        for rank,(_,a) in enumerate(rank_data[1:],start=1):
            for key in ref:
                if (a[key].shape!=ref[key].shape or not np.allclose(
                    a[key],ref[key],rtol=5e-9,atol=2e-9
                )):
                    raise RuntimeError(f"MPI rank replication mismatch: {name} rank{rank} {key}")
        vtk={}
        grid_signature={}
        for variable in ("q","p","e_q","e_p"):
            pvd=case/"pv"/variable/f"{variable}.pvd"
            entries=ET.parse(pvd).findall(".//DataSet")
            if len(entries)!=11:
                raise RuntimeError(f"Missing full VTK collection: {pvd}")
            geometries=[]
            for k,e in enumerate(entries):
                if abs(float(e.attrib["timestep"])-ref["t"][k])>1e-10:
                    raise RuntimeError(f"Wrong physical VTK sampling: {pvd}")
                filename=Path(e.attrib["file"])
                if filename.name!=str(filename):
                    raise RuntimeError(f"Unsafe shared VTK relative path {pvd}")
                vtu=pvd.parent/filename
                piece=ET.fromstring(vtu.read_bytes()).find(".//Piece")
                if piece is None or int(piece.attrib["NumberOfPoints"])<=1:
                    raise RuntimeError(f"Invalid genuine GetFEM VTU: {vtu}")
                if int(piece.attrib["NumberOfCells"]) < 1:
                    raise RuntimeError(f"Missing GetFEM VTK cells: {vtu}")
                value=piece.find(f"./PointData/DataArray[@Name='{variable}']")
                if value is None:
                    raise RuntimeError(f"Missing original field in VTU: {vtu}")
                geometry=(
                    piece.attrib["NumberOfPoints"],
                    piece.attrib["NumberOfCells"],
                    (piece.find("./Points/DataArray").text or "").strip(),
                    tuple((field.attrib["Name"], (field.text or "").strip())
                          for field in piece.findall("./Cells/DataArray"))
                )
                geometries.append(geometry)
            if any(geometry != geometries[0] for geometry in geometries):
                raise RuntimeError(f"Changing true VTK mesh: {name}/{variable}")
            grid_signature[variable]=geometries[0]
            vtk[variable]=len(entries)
        result.update({
            "status":"PASS_REAL_ANISOTROPIC_MIXED_CN_ENERGY_AND_VTK",
            "dofs":int(ref["z"].shape[1]),
            "max_abs_step_defect":float(ref_meta["max_abs_step_defect"]),
            "max_abs_cumulative_defect":float(ref_meta["max_abs_cumulative_defect"]),
            "max_abs_left_power":float(ref_meta["max_abs_left_power"]),
            "rank_traces_agree":True,
            "VTK_frames":vtk,
        })
        trace_by_case[name]=ref
        grids[name]=grid_signature
    for name, run in trace_by_case.items():
        errors={}
        for key in ("t","z","H","power_left","power_bottom","power_right",
                    "power_top","signed_work","energy_step_defects"):
            a,b=trace_by_case["mumps_mpi1"][key],run[key]
            if a.shape!=b.shape:
                raise RuntimeError(f"Different FE DOFs/time shape {name}/{key}")
            delta=float(np.max(np.abs(a-b)))
            errors[key]=delta
            if not np.allclose(a,b,rtol=2e-8,atol=5e-9):
                raise RuntimeError(f"Actual full-state LU/MPI parity FAIL {name}/{key}: {delta}")
        if grids[name]!=grids["mumps_mpi1"]:
            raise RuntimeError(f"Different original GetFEM VTK mesh: {name}")
        report["comparisons_to_mumps_mpi1"][name]=errors
    report["status"]="PASS_WAVE2D_CONSERVATIVE_FIVE_LU_MPI_CN_PORTS_VTK"
except Exception as exc:
    report["errors"].append(f"{type(exc).__name__}: {exc}")
# Report all failed native configurations, not only the first one
# encountered by the sequential scientific data verifier.
for name,case in report["cases"].items():
    if case["exit_code"] != 0:
        message=f"Native {name} ended with exit={case['exit_code']}"
        if message not in report["errors"]:
            report["errors"].append(message)

(folder/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
archive=Path(str(folder)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for f in sorted(folder.rglob("*")):
        if f.is_file():
            z.write(f,arcname=str(f.relative_to(folder)))
print(json.dumps({
    "status":report["status"],
    "cases":report["cases"],
    "comparisons_to_mumps_mpi1":report["comparisons_to_mumps_mpi1"],
    "errors":report["errors"],
    "archive":str(archive),
},indent=2))
if report["status"]!="PASS_WAVE2D_CONSERVATIVE_FIVE_LU_MPI_CN_PORTS_VTK":
    raise SystemExit(1)
PY
echo "COMPLETE canonical 2D anisotropic five-backend physical installation acceptance: $output"
