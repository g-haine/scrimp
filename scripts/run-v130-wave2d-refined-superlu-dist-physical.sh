#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE.
#
# NATIVE SCIENTIFIC GATE (CANDIDATE ONLY):
# Original examples.wave.py, SuperLU_DIST MPI1/2, 10 CN steps,
# opt-in ReplaceTinyPivot+IterRefine, strict true four-port H-work
# balances, exact global FE history, rank parity, four genuine VTK.
# Compare against already-qualified ORIGINAL MUMPS MPI1 source archive.
# Never change the canonical example, PETSc defaults or accepted gates.
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
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
    echo "Wrong pinned original GetFEM/PETSc image revision: $revision" >&2
    exit 2
}
baseline="$root/sandbox/outputs/v130-wave2d-conservative-five-lu-99d99f2d-tEOjqgHy.zip"
[[ -r "$baseline" ]] || {
    echo "Missing real MUMPS 1-rank canonical wave evidence: $baseline" >&2
    exit 2
}
python3 -m py_compile \
    examples/wave.py \
    sandbox/probe_v130_wave2d_refined_superlu_dist.py \
    tests/test_v130_wave2d_refined_dist_native.py
python3 tests/test_v130_wave2d_refined_dist_native.py
python3 tests/test_v130_wave2d_conservative_pilot.py
python3 tests/test_v130_interface_pv_production_policy.py
python3 scripts/check-v130-published-immutable.py
bash -n "$0"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-wave2d-refined-dist-physical-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-wave2d-optin-refined-superlu-dist-native-physics-v1","commit":"%s","image":"%s","image_revision":"%s","baseline":"%s","not_a_default_backend_qualification":true,"physical_model":"examples.wave:wave_eq","explicit_optin":true,"getfem_unchanged":true}\n' \
    "$(git rev-parse HEAD)" "$image" "$revision" "$(basename "$baseline")" \
    > "$out/manifest.json"

echo "NATIVE 2D ANISOTROPIC WAVE — optional refined SuperLU_DIST"
echo "MUMPS 10-step CN reference loaded from original archive (NOT rerun)"
echo "SuperLU_DIST tiny-pivot + internal iterative-refinement MPI1/2"
echo "Same q,p,e_q,e_p + four physical ports; true H-H0-work check"
echo "11 original FE states, complete MPI rank parity, 4 x 11 VTK per case"
echo "NO change to SCRIMP's canonical model/backend defaults"
echo "EVIDENCE: $out"

for ranks in 1 2; do
    name="refined_dist_mpi${ranks}"
    mkdir -p "$out/$name" "$out/work_$name"
    echo "RUN $name"
    rc=0
    "${docker_cmd[@]}" run --rm \
        --user "$(id -u):$(id -g)" \
        --env HOME=/tmp --env MPLBACKEND=Agg \
        --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
        --env SCRIMP_V130_WAVE2D_ENABLE_REFINED_DIST=1 \
        --env SCRIMP_TS_FAILURE_DIAGNOSTICS=1 \
        --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
        --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
        --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
        --env SCRIMP_MATRIX_MEMORY_MODE=compact \
        --env SCRIMP_MUMPS_RESET_POLICY=on_change \
        --env SCRIMP_GETFEM_RHS_MODE=auto \
        --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
        --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
        --volume "$root/examples/wave.py:/workspace/scrimp/examples/wave.py:ro" \
        --volume "$root/sandbox/probe_v130_wave2d_refined_superlu_dist.py:/workspace/scrimp/sandbox/probe_v130_wave2d_refined_superlu_dist.py:ro" \
        --volume "$out/$name:/artifacts:rw" \
        --volume "$out/work_$name:/scrimp-work:rw" \
        --workdir /scrimp-work \
        "$image" timeout -k 5s "${SCRIMP_V130_WAVE2D_REFINED_TIMEOUT:-150s}" \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_wave2d_refined_superlu_dist.py \
            --solver superlu_dist --mesh-size 0.2 \
            --final-time 0.08 --time-step 0.008 \
            --check --export --output /artifacts \
        > "$out/$name.host.log" 2>&1 || rc=$?
    printf '%s\n' "$rc" > "$out/$name.exitcode"
    if [[ "$rc" != 0 ]]; then
        echo "NOT QUALIFIED $name exit=$rc; native logs retained" >&2
        tail -n 35 "$out/$name.host.log" >&2
    fi
done

python3 - "$out" "$baseline" <<'PY'
import hashlib
import io
import json
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

root, baseline_file = map(Path, sys.argv[1:])
report = {
    "schema": "scrimp-v130-wave2d-refined-superlu-dist-full-physical-v1",
    "status": "FAIL_OR_UNQUALIFIED_FULL_PHYSICAL_REFINED_SUPERLU_DIST",
    "original_unmodified_model": "examples/wave.py",
    "experimental_solver_options": {
        "mat_superlu_dist_replacetinypivot": True,
        "mat_superlu_dist_iterrefine": True,
    },
    "cases": {}, "baseline": {}, "errors": [],
    "superlu_dist_default_backend_qualified": False,
    "previous_three_genuine_wave_backend_passes_preserved": True,
    "original_five_default_backends_gate": "3/5",
}

keys = ("t", "z", "H", "power_bottom", "power_right", "power_top",
        "power_left", "signed_work", "energy_step_defects")
cases = ((1, "refined_dist_mpi1"), (2, "refined_dist_mpi2"))
for ranks, name in cases:
    report["cases"][name] = {
        "exit_code": int((root/f"{name}.exitcode").read_text()),
        "MPI": ranks,
    }

def read_numpy(blob):
    with np.load(io.BytesIO(blob), allow_pickle=False) as source:
        return {key: source[key].copy() for key in source.files}

def geometry(vtu_blob, field):
    piece = ET.fromstring(vtu_blob).find(".//Piece")
    if piece is None:
        raise RuntimeError(f"Missing original GetFEM VTK Piece: {field}")
    n_points = int(piece.attrib["NumberOfPoints"])
    n_cells = int(piece.attrib["NumberOfCells"])
    if n_points < 2 or n_cells < 1:
        raise RuntimeError(f"Empty physical GetFEM VTK mesh: {field}")
    if piece.find(f"./PointData/DataArray[@Name='{field}']") is None:
        raise RuntimeError(f"Missing actual GetFEM PointData {field}")
    points = piece.find("./Points/DataArray")
    cells = piece.findall("./Cells/DataArray")
    if points is None or not cells:
        raise RuntimeError(f"Malformed original GetFEM grid: {field}")
    hasher = hashlib.sha256()
    hasher.update(str((n_points,n_cells)).encode())
    hasher.update((points.text or "").strip().encode())
    for cell in cells:
        hasher.update((cell.attrib.get("Name","")).encode())
        hasher.update((cell.text or "").strip().encode())
    return hasher.hexdigest()

def native_case_pvd(case, field, time):
    folder = case/"pv"/field
    pvd = folder/f"{field}.pvd"
    entries = ET.parse(pvd).findall(".//DataSet")
    if len(entries) != 11:
        raise RuntimeError(f"Wrong genuine VTK time-frame count: {pvd}")
    signatures = []
    for k, frame in enumerate(entries):
        if abs(float(frame.attrib["timestep"])-float(time[k])) > 1e-10:
            raise RuntimeError(f"Incorrect VTK timestamp: {pvd}")
        name = Path(frame.attrib["file"])
        if name.name != str(name):
            raise RuntimeError(f"Unsafe VTK reference: {pvd} / {name}")
        signatures.append(geometry((folder/name).read_bytes(),field))
    if len(set(signatures)) != 1:
        raise RuntimeError(f"Nonstationary VTK geometry: {pvd}")
    return signatures[0]

def baseline_pvd(arch, field, time):
    prefix = f"mumps_mpi1/pv/{field}/"
    entries = ET.fromstring(arch.read(prefix+f"{field}.pvd")).findall(".//DataSet")
    if len(entries) != 11:
        raise RuntimeError(f"Wrong MUMPS reference VTK series: {field}")
    signatures=[]
    for k, frame in enumerate(entries):
        if abs(float(frame.attrib["timestep"])-float(time[k])) > 1e-10:
            raise RuntimeError(f"Wrong MUMPS reference timestamp: {field}")
        filename = Path(frame.attrib["file"])
        if filename.name != str(filename):
            raise RuntimeError(f"Unsafe reference VTK name: {field}")
        signatures.append(geometry(arch.read(prefix+filename.name),field))
    if len(set(signatures)) != 1:
        raise RuntimeError(f"MUMPS reference geometry varies: {field}")
    return signatures[0]

try:
    with zipfile.ZipFile(baseline_file) as arch:
        baseline = read_numpy(arch.read("mumps_mpi1/rank0.npz"))
        base_meta = json.loads(arch.read("mumps_mpi1/rank0.json"))
        if (base_meta.get("status") != "PASS_ORIGINAL_FE_PHYSICS"
                or base_meta.get("requested_LU") != "mumps"
                or base_meta.get("selected_LU") != "mumps"
                or base_meta.get("steps") != 10
                or base_meta.get("DOFs") != 2230
                or not all(base_meta.get("physical_gates",{}).values())):
            raise RuntimeError("Original archived MUMPS baseline not qualified")
        if baseline["z"].shape != (11,2230):
            raise RuntimeError("Wrong original full FE reference dimensions")
        base_geometry = {
            field: baseline_pvd(arch,field,baseline["t"])
            for field in ("q","p","e_q","e_p")
        }
    report["baseline"] = {
        "verified_real_native_MUMPS_MPI1": True,
        "state_SHA256": base_meta.get("state_SHA256"),
        "DOFs": 2230, "steps": 10,
        "VTK_original_fields": list(base_geometry),
    }

    results={}
    for ranks,name in cases:
        case=root/name
        recorded=report["cases"][name]
        log=(root/f"{name}.host.log").read_text(errors="replace")
        marker="V130_WAVE2D_REFINED_SUPERLU_DIST_ENABLED"
        count = log.count(marker)
        recorded["options_injection_markers"] = count
        relevant_unused = [line for line in log.splitlines()
            if "Option left:" in line and any(x in line for x in (
                "mat_superlu_dist_replacetinypivot",
                "mat_superlu_dist_iterrefine",
            ))]
        recorded["unused_relevant_PETSc_options"] = relevant_unused
        if recorded["exit_code"] != 0:
            raise RuntimeError(f"Full physical example exited nonzero: {name}")
        if count != ranks or relevant_unused:
            raise RuntimeError(f"PETSc factor options not verified in all ranks: {name}")
        rank_data=[]
        for rank in range(ranks):
            meta=json.loads((case/f"rank{rank}.json").read_text())
            with np.load(case/f"rank{rank}.npz",allow_pickle=False) as f:
                arrays={key:f[key].copy() for key in f.files}
            if (
                meta.get("status") != "PASS_ORIGINAL_FE_PHYSICS"
                or meta.get("rank") != rank
                or meta.get("MPI") != ranks
                or meta.get("requested_LU") != "superlu_dist"
                or meta.get("selected_LU") != "superlu_dist"
                or meta.get("steps") != 10
                or meta.get("DOFs") != 2230
                or meta.get("DAE_initialization",{}).get("status") != "PASS"
                or not meta.get("DAE_initialization",{}).get("differential_initial_exact")
                or not all(meta.get("physical_gates",{}).values())
                or any(key not in arrays for key in keys)
            ):
                raise RuntimeError(f"Native SNES/TS/DAE/4-port gates failed: {name}/rank{rank}")
            t=arrays["t"]; z=arrays["z"]; H=arrays["H"]
            if (
                t.shape != (11,) or z.shape != (11,2230)
                or H.shape != (11,)
                or not np.all(np.isfinite(z))
                or not np.all(np.isfinite(H))
                or not np.allclose(t,np.arange(11)*.008,rtol=0,atol=1e-10)
                or any(arrays[key].shape != baseline[key].shape for key in keys)
            ):
                raise RuntimeError(f"Unqualified full original FE series: {name}/rank{rank}")
            if (
                meta["max_abs_step_defect"] >= 1e-10
                or meta["max_abs_cumulative_defect"] >= 1e-10
                or meta["max_abs_left_power"] <= 1e-8
            ):
                raise RuntimeError(f"Insufficient real CN work/Hamiltonian accuracy: {name}/rank{rank}")
            digest=hashlib.sha256(np.ascontiguousarray(z).tobytes()).hexdigest()
            if digest != meta["state_SHA256"]:
                raise RuntimeError(f"Original accepted-state SHA256 mismatch: {name}/rank{rank}")
            for key in keys:
                if not np.all(np.isfinite(arrays[key])):
                    raise RuntimeError(f"Nonfinite original FE/power values: {name}/{key}")
            rank_data.append(arrays)
        ref=rank_data[0]
        for rank,other in enumerate(rank_data[1:],1):
            for key in keys:
                if (other[key].shape != ref[key].shape or
                    not np.allclose(other[key],ref[key],rtol=5e-9,atol=2e-9)):
                    raise RuntimeError(f"MPI rank-replication mismatch {name}/rank{rank}/{key}")
        diffs={}
        for key in keys:
            delta=float(np.max(np.abs(ref[key]-baseline[key])))
            diffs[key]=delta
            if not np.allclose(ref[key],baseline[key],rtol=2e-8,atol=5e-9):
                raise RuntimeError(f"Full original MUMPS/PETSc FE parity fails: {name}/{key}: {delta}")
        grids={field:native_case_pvd(case,field,ref["t"])
               for field in ("q","p","e_q","e_p")}
        if grids != base_geometry:
            raise RuntimeError(f"Wrong original GetFEM VTK geometry vs MUMPS: {name}")
        recorded.update({
            "status":"PASS_OPTIN_REFINED_SUPERLU_DIST_FULL_CN_PHYSICS",
            "rank_trajectories_agree":True,
            "max_abs_step_energy_defect":float(np.max(np.abs(ref["energy_step_defects"]))),
            "max_abs_energy_work_defect":
                float(np.max(np.abs(ref["H"]-ref["H"][0]-ref["signed_work"]))),
            "max_abs_power_left":float(np.max(np.abs(ref["power_left"]))),
            "max_abs_vs_MUMPS_reference":diffs,
            "VTK_fields":{name:11 for name in grids},
            "original_VTK_geometry_matches_MUMPS":True,
        })
        results[name]=ref
    for key in keys:
        diff=float(np.max(np.abs(results["refined_dist_mpi1"][key]-
                                 results["refined_dist_mpi2"][key])))
        if not np.allclose(results["refined_dist_mpi1"][key],
                           results["refined_dist_mpi2"][key],
                           rtol=2e-8,atol=5e-9):
            raise RuntimeError(f"Opt-in SuperLU_DIST full-MPI state/power parity FAIL {key}: {diff}")
    report["status"]="PASS_OPTIN_REFINED_SUPERLU_DIST_FULL_WAVE2D_MPI1_MPI2_CANDIDATE_ONLY"
except Exception as exc:
    report["errors"].append(f"{type(exc).__name__}: {exc}")
for ranks,name in cases:
    case=report["cases"][name]
    if case["exit_code"] != 0:
        message=f"Original {name} exited with code {case['exit_code']}"
        if message not in report["errors"]:
            report["errors"].append(message)
(root/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for file in sorted(root.rglob("*")):
        if file.is_file():
            z.write(file,arcname=str(file.relative_to(root)))
print(json.dumps({
    "status":report["status"],
    "cases":report["cases"],
    "errors":report["errors"],
    "archive":str(archive),
    "default_superlu_dist_backend_still_unqualified":True,
    "original_five_backend_gate_still":"3/5",
},indent=2))
if report["status"]!="PASS_OPTIN_REFINED_SUPERLU_DIST_FULL_WAVE2D_MPI1_MPI2_CANDIDATE_ONLY":
    raise SystemExit(1)
PY
echo "COMPLETE opt-in refined SuperLU_DIST REAL 2D wave physics pilot: $out"
