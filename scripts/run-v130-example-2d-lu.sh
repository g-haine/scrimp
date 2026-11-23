#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Actual original 2D example P0 acceptance: MUMPS vs SuperLU_DIST.
# MPI=1 ONLY until public ParaView export is verified MPI-safe.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
    docker_exec=(docker)
else
    docker_exec=(sudo docker)
fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_exec[@]}" image inspect \
    --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == 1d6d0c90115213298047b28237a7f38087563a35 ]] || {
    echo "Qualified GetFEM/PETSc image required, found revision=$revision" >&2
    exit 2
}
python3 -m py_compile \
    sandbox/probe_v130_main_2d_examples.py \
    examples/heat_wave.py examples/shallow_water.py \
    tests/test_v130_main_example_solver_policy.py
python3 tests/test_v130_main_example_solver_policy.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-main-2d-lu-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-main-2d-lu-pilot-campaign-v1","source":"%s","image_revision":"%s","status":"NATIVE_PENDING","models":["heat_wave","shallow_water"],"solvers":["mumps","superlu_dist"],"MPI":[1],"dt":0.001,"t_final":0.004,"parallel_ParaView_export_not_qualified":true,"not_a_benchmark":true}\n' \
  "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Original SCRIMP 2-D Heat-Wave and inviscid Saint-Venant pilot"
echo "MUMPS vs SuperLU_DIST, 1 rank only: first native results before MPI I/O audit"
echo "Original weak forms, original public Hamiltonian PNG and ParaView PVD/VTU"
echo "t=0..0.004; dt=0.001; 4 accepted BDF2 steps expected"
echo "No changes or mounts in examples/published"
echo "EVIDENCE: $out"

for example in heat_wave shallow_water; do
    for backend in mumps superlu_dist; do
        name="${example}_${backend}_mpi1"
        mkdir -p "$out/$name" "$out/work_$name"
        echo "RUN $name original physical 2D example"
        if ! "${docker_exec[@]}" run --rm \
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
            --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
            --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
            --volume "$root/sandbox/probe_v130_main_2d_examples.py:/probe-2d.py:ro" \
            --volume "$out/$name:/artifacts:rw" \
            --volume "$out/work_$name:/scrimp-work:rw" \
            --workdir /scrimp-work \
            "$image" timeout -k 5s "${SCRIMP_V130_2D_PILOT_TIMEOUT:-240s}" \
            mpiexec --oversubscribe --bind-to none -n 1 \
            python -u /probe-2d.py \
              --example "$example" --backend "$backend" --out /artifacts \
            > "$out/$name.host.log" 2>&1; then
            echo "FAIL original native $name; all evidence retained at $out" >&2
            tail -n 130 "$out/$name.host.log" >&2
            exit 1
        fi
        grep -a '^PASS ORIGINAL 2D EXAMPLE' "$out/$name.host.log" || {
            echo "FAIL missing physical example proof in $name" >&2
            exit 1
        }
        [[ -s "$out/$name/original_2d_example_report.json" \
            && -s "$out/$name/original_2d_example_trace.npz" ]] || {
            echo "FAIL incomplete native evidence $name" >&2
            exit 1
        }
    done
done

python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
import numpy as np
root=Path(sys.argv[1])
summaries={}
for example in ("heat_wave","shallow_water"):
    all_runs={}
    for solver in ("mumps","superlu_dist"):
        folder=root/f"{example}_{solver}_mpi1"
        report=json.loads((folder/"original_2d_example_report.json").read_text())
        assert report["status"]=="PASS_NATIVE_ORIGINAL_2D_EXAMPLE_SINGLE_MPI"
        assert report["example"]==example and report["factor_solver"]==solver
        assert report["steps"]==4 and report["MPI"]==1
        assert report["rejected_steps"]==0 and report["SNES_failures"]==0
        with np.load(folder/"original_2d_example_trace.npz",
                     allow_pickle=False) as a:
            trace={k:a[k] for k in a.files}
        assert trace["times"].shape==(5,)
        assert np.all(np.isfinite(trace["states"]))
        assert trace["states"].shape[0]==5
        assert trace["H"].shape==(5,)
        assert report["full_state_DOFs"]==trace["states"].shape[1]
        all_runs[solver]=(report,trace)
    left,right=all_runs["mumps"][1],all_runs["superlu_dist"][1]
    if set(left)!=set(right):
        raise RuntimeError(f"Different original physical observables: {example}")
    if not np.array_equal(left["times"],right["times"]):
        raise RuntimeError(f"Different original time grid: {example}")
    if left["states"].shape!=right["states"].shape:
        raise RuntimeError(f"Different original finite-element model: {example}")
    differences={}
    for key in sorted(left):
        a,b=left[key],right[key]
        if a.shape!=b.shape:
            raise RuntimeError(f"Different original observable shape: {example} {key}")
        if key=="states":
            delta=float(np.max(np.linalg.norm(a-b,axis=1)))
        else:
            delta=float(np.max(abs(a-b)))
        if not np.isfinite(delta) or delta>=1e-8:
            raise RuntimeError(
                f"Original {example} MUMPS/SuperLU_DIST {key} parity FAIL: {delta}"
            )
        differences[key]=delta
    summaries[example]={
        "full_state_DOFs":int(left["states"].shape[1]),
        "max_discrepancies":differences,
        "original_water_height_mass_change_mumps":(
            None if example!="shallow_water" else
            float(left["integral_h"][-1]-left["integral_h"][0])
        ),
        "original_model_energy_change_mumps":float(left["H"][-1]-left["H"][0]),
        "original_paraview_and_H_png_both_backends":True,
    }

summary={
    "schema":"scrimp-v130-main-2d-lu-native-campaign-v1",
    "status":"PASS_2D_MAIN_EXAMPLES_MUMPS_SUPERLU_DIST_MPI1",
    "actual_native_runs":4,
    "models":summaries,
    "full_shipped_time_horizon_not_qualified":True,
    "parallel_ParaView_export_not_qualified":True,
    "physically_correct_model_not_inferred_from_solver_parity":True,
    "published_snapshots_or_outputs_modified":False,
    "speed_benchmark_claimed":False,
}
(root/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for f in sorted(root.rglob("*")):
        if f.is_file():
            z.write(f,arcname=str(f.relative_to(root)))
print("SHAREABLE ORIGINAL 2D EVIDENCE ZIP:",archive)
PY
echo "COMPLETE physical 2D main examples native evidence: $out"
