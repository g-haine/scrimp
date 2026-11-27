#!/usr/bin/env bash
# Native canonical Heat-Wave fixed-mesh CN physical refinement diagnostic.
# Original GetFEM Concentric h=.05 is unchanged; no mesh refinement.
# Four real full pH-DAE solutions at common final time .2:
# 20/40/80/160 fixed-CN steps, dt=.01/.005/.0025/.00125.
# Exact all-accepted GetFEM coefficient vectors and physical power
# checks are archived at each dt. User-facing original example is
# still the sole canonical physics; this is a scientific test runner.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
 --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
if [[ "$revision" != 1d6d0c90115213298047b28237a7f38087563a35 ]]; then
    echo "Unqualified Docker GetFEM/PETSc image revision: $revision" >&2
    exit 2
fi

python3 -m py_compile examples/heat_wave.py \
  sandbox/probe_v130_heat_wave_mpi_extended.py \
  sandbox/compare_v130_heat_wave_CN_dt_refinement.py \
  tests/test_v130_heat_wave_CN_dt_refinement_source.py
python3 tests/test_v130_heat_wave_CN_dt_refinement_source.py
python3 tests/test_v130_heat_wave_boundary_initial_data.py
python3 scripts/check-v130-published-immutable.py
bash -n "$0"

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-heat-wave-CN-dt-refinement-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
python3 - "$out" "$image" "$revision" <<'PY'
import json,subprocess,sys
from pathlib import Path
root=Path(sys.argv[1])
(root/"manifest.json").write_text(json.dumps({
    "schema":"scrimp-v130-original-heatwave-CN-dt-convergence-v1",
    "source":subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
    "image":sys.argv[2],
    "image_revision":sys.argv[3],
    "canonical_example":"examples/heat_wave.py",
    "original_model_unchanged":True,
    "same_mesh_and_original_FEM_weak_forms":True,
    "space_mesh_size_h":.05,
    "original_FE_DOF":7689,
    "profile":"gaussian",
    "final_time":.2,
    "solver":"mumps",
    "MPI":1,
    "dts":[.01,.005,.0025,.00125],
    "steps":[20,40,80,160],
    "output_every_real_accepted_CN_step":True,
    "no_spatial_convergence_claim":True,
    "no_MPI_parity_retest_required":True,
    "no_gallery_eigenmode_claim":True,
},indent=2)+"\n")
PY
echo "SCRIMP v1.3.0 original Heat-Wave full GetFEM FE CN physical refinement"
echo "ONE fixed original h=.05 mesh; MUMPS MPI1; gaussian physically compatible state"
echo "dt=.01, .005, .0025, .00125; t_f=.2, 20/40/80/160 real CN steps"
echo "All 7689 original GetFEM coefficients / saved state, no reconstruction"
echo "Physical heat/wave/interface power must PASS at every native dt"
echo "Original full FE time-refinement ratios REPORTED only (no order-two claim)"
echo "Evidence: $out"
overall=0
for record in \
    "dt_0p01:0.01:20" \
    "dt_0p005:0.005:40" \
    "dt_0p0025:0.0025:80" \
    "dt_0p00125:0.00125:160"
do
    IFS=: read -r label dt steps <<< "$record"
    mkdir -p "$out/$label" "$out/work_$label"
    if [[ "$overall" != 0 ]]; then
        echo "SKIPPED_AFTER_PREVIOUS_FAILURE" > "$out/$label.exitcode"
        echo "SKIPPED $label after an unqualified native run"
        continue
    fi
    echo "RUN original CN native $label dt=$dt steps=$steps MPI1"
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
        --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
        --volume "$root/sandbox/measure_v130_hw_dual_facet_physics.py:/workspace/scrimp/sandbox/measure_v130_hw_dual_facet_physics.py:ro" \
        --volume "$root/sandbox/probe_v130_heat_wave_mpi_extended.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_mpi_extended.py:ro" \
        --volume "$out/$label:/artifacts:rw" \
        --volume "$out/work_$label:/scrimp-work:rw" \
        --workdir /scrimp-work "$image" \
        timeout -k 10s "${SCRIMP_V130_HW_DT_REFINEMENT_TIMEOUT:-1800s}" \
        mpiexec --oversubscribe --bind-to none -n 1 \
        python -u /workspace/scrimp/sandbox/probe_v130_heat_wave_mpi_extended.py \
            --out /artifacts --dt "$dt" --steps "$steps" \
            --profile gaussian --store-full-states \
        >"$out/$label.host.log" 2>&1 || rc=$?
    printf '%s\n' "$rc" > "$out/$label.exitcode"
    if [[ "$rc" == 0 ]]; then
        echo "ORIGINAL GetFEM/PETSc CN complete: $label"
    else
        echo "UNQUALIFIED: $label returned exit $rc" >&2
        tail -n 90 "$out/$label.host.log" >&2
        overall=1
    fi
done

compare_rc=0
python3 sandbox/compare_v130_heat_wave_CN_dt_refinement.py "$out" \
    || compare_rc=$?
if [[ "$compare_rc" != 0 ]]; then overall=1; fi
if [[ "$overall" == 0 ]]; then
    echo "PASS_HEATWAVE_CN_FIXED_MESH_PHYSICAL_REFINEMENT_DIAGNOSTIC"
else
    echo "FAIL_OR_UNQUALIFIED_HEATWAVE_CN_PHYSICAL_REFINEMENT — archived" >&2
fi
echo "ARCHIVE: $out.zip"
exit "$overall"
