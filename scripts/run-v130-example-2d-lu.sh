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
[[ $# == 0 ]] || { echo "Usage: $0 (set SCRIMP_V130_2D_SCOPE=heat_wave|shallow_water|all)" >&2; exit 2; }
scope="${SCRIMP_V130_2D_SCOPE:-heat_wave}"
case "$scope" in
    heat_wave) examples=(heat_wave) ;;
    shallow_water) examples=(shallow_water) ;;
    all) examples=(heat_wave shallow_water) ;;
    *) echo "Unknown SCRIMP_V130_2D_SCOPE=$scope" >&2; exit 2 ;;
esac

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
    sandbox/measure_v130_heat_wave_physics.py \
    examples/heat_wave.py examples/shallow_water.py \
    tests/test_v130_main_example_solver_policy.py
python3 tests/test_v130_main_example_solver_policy.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-main-2d-lu-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-main-2d-lu-pilot-campaign-v2","source":"%s","image_revision":"%s","status":"NATIVE_PENDING","scope":"%s","solvers":["mumps","superlu_dist"],"MPI":[1],"requested_initial_dt":0.001,"t_final":0.004,"heat_wave_fixed_CN":true,"shallow_water_BDF_adaptive":true,"coupled_heat_wave_model_physics_not_yet_qualified":true,"parallel_ParaView_export_not_qualified":true,"not_a_benchmark":true}\n' \
  "$(git rev-parse HEAD)" "$revision" "$scope" > "$out/manifest.json"
echo "SCRIMP 2-D Heat-Wave PHYSICAL REVISION and inviscid Saint-Venant pilot"
echo "MUMPS vs SuperLU_DIST, 1 rank only, prior to parallel ParaView I/O audit"
echo "Heat-Wave: mixed DG0/RT0/CG1, corrected interface and exterior, fixed CN"
echo "Shallow water: unchanged original adaptive BDF2 numerical model"
echo "t=0..0.004; heat-wave CN dt=0.001, 4 accepted steps"
echo "Per-subsystem Hamiltonians, CN midpoint dissipation and oriented interface work"
echo "NO heat-wave physical PASS without independent scientific audit"
echo "No changes or mounts in examples/published"
echo "EVIDENCE: $out"

for example in "${examples[@]}"; do
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
            --volume "$root/sandbox/measure_v130_heat_wave_physics.py:/workspace/scrimp/sandbox/measure_v130_heat_wave_physics.py:ro" \
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

python3 - "$out" "${examples[@]}" <<'PY'
import json,sys,zipfile
from pathlib import Path
import numpy as np
root=Path(sys.argv[1])
examples=tuple(sys.argv[2:])
assert examples and set(examples)<={"heat_wave","shallow_water"}
summaries={}
for example in examples:
    all_runs={}
    for solver in ("mumps","superlu_dist"):
        folder=root/f"{example}_{solver}_mpi1"
        report=json.loads((folder/"original_2d_example_report.json").read_text())
        assert report["status"]==(
            "PASS_NATIVE_INTEGRATION_HEAT_WAVE_PHYSICAL_AUDIT_PENDING"
            if example=="heat_wave" else
            "PASS_NATIVE_ORIGINAL_2D_EXAMPLE_SINGLE_MPI"
        )
        assert report["example"]==example and report["factor_solver"]==solver
        assert report["steps"]>=1 and report["MPI"]==1
        assert report["rejected_steps"]>=0 and report["SNES_failures"]==0
        assert report["PETSc_reason"]>0
        if example=="heat_wave":
            assert report["heat_wave_fixed_CN"] is True
            assert report["steps"]==4 and report["rejected_steps"]==0
        else:
            assert report["adaptive_integrator_original_behavior_preserved"] is True
        with np.load(folder/"original_2d_example_trace.npz",
                     allow_pickle=False) as a:
            trace={k:a[k] for k in a.files}
        times=trace["times"]
        assert times.ndim==1 and times.size>=2
        assert np.all(np.isfinite(times)) and np.all(np.diff(times)>0)
        assert abs(float(times[0]))<1e-12
        assert abs(float(times[-1])-0.004)<1e-11
        assert np.all(np.isfinite(trace["states"]))
        assert trace["states"].shape[0]==len(times)
        assert trace["H"].shape==(len(times),)
        if example=="heat_wave":
            diagnostics=json.loads(
                (folder/"heat_wave_physics_diagnostic.json").read_text()
            )
            assert diagnostics["scientific_status"]==(
                "DIAGNOSTIC_ONLY_MODEL_NOT_YET_QUALIFIED"
            )
            assert times.shape==(5,)
            assert np.allclose(times,np.linspace(0.,.004,5),atol=1e-11,rtol=0)
            assert "phys_D_heat_flux" in trace
            assert "phys_P_interface_sum" in trace
            assert "phys_CN_D_heat_flux" in trace
            assert "phys_CN_per_step_H_balance_error" in trace
            assert "phys_CN_per_step_heat_balance_error" in trace
            assert "phys_CN_per_step_wave_balance_error" in trace
            assert "phys_CN_cumulative_total_balance_error" in trace
            assert trace["phys_CN_D_heat_flux"].shape==times.shape
            assert diagnostics["dissipation_min"]>=-1e-10
            assert diagnostics["CN_step_powers_are_assembled_at_mean_FE_state"] is True
        assert report["full_state_DOFs"]==trace["states"].shape[1]
        all_runs[solver]=(report,trace)
    left,right=all_runs["mumps"][1],all_runs["superlu_dist"][1]
    if set(left)!=set(right):
        raise RuntimeError(f"Different original physical observables: {example}")
    common_clock=(
        left["times"].shape==right["times"].shape
        and np.allclose(left["times"],right["times"],
                        atol=1e-12,rtol=0.)
    )
    if example=="heat_wave" and not common_clock:
        raise RuntimeError(
            "Fixed CN must yield the identical accepted time grid "
            "for MUMPS and SuperLU_DIST"
        )
    # Shallow-water adaptive BDF2 may have different accepted
    # time grids. Only compare its original endpoints when they differ;
    # do not invent whole-path parity via interpolation.
    comparison_scope=(
        "WHOLE_SAVED_TRAJECTORY_FIXED_CN"
        if example=="heat_wave" else
        "WHOLE_SAVED_TRAJECTORY_AT_SAME_ADAPTIVE_CLOCK"
        if common_clock else "ENDPOINTS_ONLY_DIFFERENT_ADAPTIVE_CLOCKS"
    )
    if left["states"].shape!=right["states"].shape:
        raise RuntimeError(f"Different original finite-element model: {example}")
    differences={}
    for key in sorted(left):
        a,b=left[key],right[key]
        if not common_clock:
            if a.ndim<1 or b.ndim<1 or a.shape[1:]!=b.shape[1:]:
                raise RuntimeError(
                    f"Adaptive endpoint FE/observable mismatch: {example} {key}"
                )
            a=np.stack((a[0],a[-1]))
            b=np.stack((b[0],b[-1]))
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
        "numerical_comparison_scope":comparison_scope,
        "saved_time_count_mumps":len(left["times"]),
        "saved_time_count_superlu_dist":len(right["times"]),
        "original_heat_wave_physical_model_qualified":False,
        "original_heat_wave_model_revised_for_CN_RT_physics":
            example=="heat_wave",
        "original_water_height_mass_change_mumps":(
            None if example!="shallow_water" else
            float(left["integral_h"][-1]-left["integral_h"][0])
        ),
        "original_model_energy_change_mumps":float(left["H"][-1]-left["H"][0]),
        "original_paraview_and_H_png_both_backends":True,
    }

summary={
    "schema":"scrimp-v130-main-2d-lu-native-campaign-v1",
    "status":"PASS_BACKEND_NUMERICAL_PARITY_PHYSICAL_AUDIT_OPEN",
    "actual_native_runs":len(examples)*2,
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
