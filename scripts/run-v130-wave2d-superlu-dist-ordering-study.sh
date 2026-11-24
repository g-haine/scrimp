#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.
#
# Research-only SuperLU_DIST numerical zero-pivot discrimination.
# TWO strictly non-regularizing permutations, on the unchanged wave.py,
# with FULL original CN/DAE/energy/four-port/VTK acceptance. MPI1 only.
# Neither a success here nor a failure alters the five-solver release gate.
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
    echo "Unqualified GetFEM/PETSc native image revision: $revision" >&2
    exit 2
}
baseline="$root/sandbox/outputs/v130-wave2d-conservative-five-lu-99d99f2d-tEOjqgHy.zip"
[[ -r "$baseline" ]] || {
    echo "Missing verified ORIGINAL MUMPS baseline $baseline" >&2
    echo "Do not manufacture reference data; preserve the previous archive." >&2
    exit 2
}
python3 -m py_compile \
    sandbox/probe_v130_wave2d_superlu_dist_ordering.py \
    tests/test_v130_wave2d_conservative_pilot.py
python3 tests/test_v130_wave2d_conservative_pilot.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-wave2d-superlu-dist-ordering-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"wave2d-superlu-dist-permutation-study-v1","source":"%s","image":"%s","image_revision":"%s","original_reference":"%s","baseline_sha256":"%s","not_full_five_lu_acceptance":true,"no_regularization":true}\n' \
    "$(git rev-parse HEAD)" "$image" "$revision" "$(basename "$baseline")" \
    "$(sha256sum "$baseline" | cut -d' ' -f1)" > "$out/manifest.json"
echo "SUPERLU_DIST exact pH-DAE, 10 CN steps MPI1, original 4 signed ports"
echo "Test ONLY distinct row/column permutations; NO pivot replacement"
echo "Original MUMPS baseline kept; compare full 2230 FE states and VTK"
echo "EVIDENCE: $out"
for selection in rowperm_awpm colperm_mmd; do
    echo "RUN $selection"
    mkdir -p "$out/$selection" "$out/work_$selection"
    rc=0
    "${docker_cmd[@]}" run --rm \
        --user "$(id -u):$(id -g)" \
        --env HOME=/tmp --env MPLBACKEND=Agg \
        --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
        --env SCRIMP_TS_FAILURE_DIAGNOSTICS=1 \
        --env SCRIMP_MPI_TRACE=1 \
        --env SCRIMP_V130_SUPERLU_DIST_ORDERING="$selection" \
        --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
        --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
        --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
        --env SCRIMP_MATRIX_MEMORY_MODE=compact \
        --env SCRIMP_MUMPS_RESET_POLICY=on_change \
        --env SCRIMP_GETFEM_RHS_MODE=auto \
        --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
        --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
        --volume "$root/examples/wave.py:/workspace/scrimp/examples/wave.py:ro" \
        --volume "$root/sandbox/probe_v130_wave2d_superlu_dist_ordering.py:/workspace/scrimp/sandbox/probe_v130_wave2d_superlu_dist_ordering.py:ro" \
        --volume "$out/$selection:/artifacts:rw" \
        --volume "$out/work_$selection:/scrimp-work:rw" \
        --workdir /scrimp-work \
        "$image" timeout -k 5s "${SCRIMP_V130_WAVE2D_ORDER_TIMEOUT:-120s}" \
        mpiexec --oversubscribe --bind-to none -n 1 \
        python -u /workspace/scrimp/sandbox/probe_v130_wave2d_superlu_dist_ordering.py \
            --solver superlu_dist --mesh-size 0.2 \
            --final-time 0.08 --time-step 0.008 \
            --check --export --output /artifacts \
        > "$out/$selection.host.log" 2>&1 || rc=$?
    echo "$rc" > "$out/$selection.exitcode"
    if [[ "$rc" != 0 ]]; then
        echo "UNQUALIFIED $selection exit=$rc (full native log preserved)" >&2
        grep -aE 'MAIN TS FAILURE DIAGNOSTIC|IJacobian: operator shape|ZERO|zeropivot|Unknown|Traceback' \
            "$out/$selection.host.log" | tail -n 12 >&2 || true
    fi
done

python3 - "$out" "$baseline" <<'PY'
import io
import json
import sys
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile, ZIP_DEFLATED
import numpy as np

folder, reference_file = map(Path, sys.argv[1:3])
with ZipFile(reference_file) as z:
    with np.load(io.BytesIO(z.read("mumps_mpi1/rank0.npz")),
                 allow_pickle=False) as data:
        baseline = {k:data[k].copy() for k in data.files}
report = {
    "schema":"wave2d-true-superlu-dist-permutation-mpi1-v1",
    "status":"FAIL_OR_UNQUALIFIED_SUPERLU_DIST_ORDERING_EXPERIMENT",
    "cases":{},"baseline":"MUMPS_MPI1_original_99d99f2d",
    "full_five_solver_gate_passed":False,
    "original_three_solver_PASSes_preserved":True,
}
passes=0
for selection in ("rowperm_awpm","colperm_mmd"):
    path=folder/selection
    rc=int((folder/f"{selection}.exitcode").read_text())
    outcome={"exit_code":rc,"status":"UNQUALIFIED"}
    report["cases"][selection]=outcome
    log=(folder/f"{selection}.host.log").read_text(errors="replace")
    outcome["zero_row_traces"]=[
        row for row in log.splitlines()
        if "IJacobian: operator shape" in row
    ][:8]
    outcome["diagnostics"]=[
        row for row in log.splitlines()
        if "MAIN TS FAILURE DIAGNOSTIC:" in row
    ][:4]
    outcome["ordering_confirmed"]=(
        "V130_SUPERLU_DIST_ORDERING_APPLIED "+selection in log
    )
    option_name = (
        "mat_superlu_dist_rowperm" if selection == "rowperm_awpm"
        else "mat_superlu_dist_colperm"
    )
    # Do not mistake an inserted-but-ignored PETSc option for a
    # successful ordering experiment: inspect PETSc's unused list.
    unused_lines = [line for line in log.splitlines()
                    if "Option left:" in line and option_name in line]
    outcome["ordering_option_consumed"]=not bool(unused_lines)
    outcome["unused_ordering_lines"]=unused_lines
    if rc:
        outcome["status"]="FAIL_FIRST_CN_OR_SOLVER"
        continue
    try:
        if not outcome["ordering_confirmed"]:
            raise RuntimeError("Configured PETSc permutation not confirmed")
        if not outcome["ordering_option_consumed"]:
            raise RuntimeError(
                "PETSc did not consume intended ordering option "
                + option_name
            )
        meta=json.loads((path/"rank0.json").read_text())
        if (meta["status"]!="PASS_ORIGINAL_FE_PHYSICS"
            or meta["MPI"]!=1 or meta["requested_LU"]!="superlu_dist"
            or meta["selected_LU"]!="superlu_dist"
            or meta["DOFs"]!=2230 or meta["steps"]!=10
            or not all(meta["physical_gates"].values())
            or meta["max_abs_step_defect"]>1e-10
            or meta["max_abs_cumulative_defect"]>1e-10
            or meta["max_abs_left_power"]<=1e-8):
            raise RuntimeError("Genuine native TSCN/DAE/energy/port gates fail")
        with np.load(path/"rank0.npz",allow_pickle=False) as data:
            actual={key:data[key].copy() for key in data.files}
        for key in ("t","z","H","power_left","power_bottom",
                    "power_right","power_top","signed_work",
                    "energy_step_defects"):
            a,b=baseline[key],actual[key]
            if a.shape!=b.shape or not np.all(np.isfinite(b)):
                raise RuntimeError(f"Invalid original state/power {key}")
            delta=float(np.max(np.abs(a-b)))
            outcome.setdefault("vs_original_MUMPS_max_abs",{})[key]=delta
            if not np.allclose(a,b,rtol=2e-8,atol=5e-9):
                raise RuntimeError(f"Fail full FE parity for {key}: {delta}")
        vtk={}
        for name in ("q","p","e_q","e_p"):
            pvd=path/"pv"/name/f"{name}.pvd"
            frames=ET.parse(pvd).findall(".//DataSet")
            if len(frames)!=11:
                raise RuntimeError(f"Missing physical VTK frames for {name}")
            for j,frame in enumerate(frames):
                if abs(float(frame.attrib["timestep"])-actual["t"][j])>1e-10:
                    raise RuntimeError(f"Incorrect actual time in {pvd}")
                file=Path(frame.attrib["file"])
                if file.name!=str(file):
                    raise RuntimeError(f"Unsafe VTK path in {pvd}")
                piece=ET.fromstring((pvd.parent/file).read_bytes()).find(".//Piece")
                if piece is None or int(piece.attrib["NumberOfPoints"])<=1:
                    raise RuntimeError(f"Malformed GetFEM VTK {pvd.parent/file}")
                if piece.find(f"./PointData/DataArray[@Name='{name}']") is None:
                    raise RuntimeError(f"No actual field {name} in VTU")
            vtk[name]=len(frames)
        outcome["VTK_frames"]=vtk
        outcome["step_defect"]=meta["max_abs_step_defect"]
        outcome["status"]="PASS_NATIVE_10_CN_PHYSICAL_ORDERING_CANDIDATE_MPI1"
        passes+=1
    except Exception as exc:
        outcome["status"]="UNQUALIFIED_POSTPROCESS"
        outcome["error"]=f"{type(exc).__name__}: {exc}"
if passes:
    report["status"]="PASS_AT_LEAST_ONE_NONREGULARIZING_ORDERING_CANDIDATE_MPI1_ONLY"
(folder/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
archive=Path(str(folder)+".zip")
with ZipFile(archive,"x",compression=ZIP_DEFLATED,compresslevel=4,
             allowZip64=True) as zipped:
    for file in sorted(folder.rglob("*")):
        if file.is_file():
            zipped.write(file,arcname=str(file.relative_to(folder)))
print(json.dumps({**report,"archive":str(archive)},indent=2),flush=True)
if passes==0:
    raise SystemExit(1)
PY
echo "CANDIDATE ONLY — five-default-backend installation acceptance still FAIL"
