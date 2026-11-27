#!/usr/bin/env bash
# Canonical nonlinear finite-dimensional full-sine PETSc TSCN scientific gate.
# Two physical cases (free, forced); original 12s/1200-step native
# MUMPS MPI1/MPI2; dt=.005 MPI1 refinement; source-generated PNG/SVG/GIF.
# Trapezoidal CN is not an exact discrete gradient of nonlinear H.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == "1" ]]; then
    docker_cmd=(docker)
fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
    --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
if [[ "$revision" != 1d6d0c90115213298047b28237a7f38087563a35 ]]; then
    echo "Unqualified GetFEM/PETSc Docker image revision: $revision" >&2
    exit 2
fi

python3 -m py_compile \
    examples/nonlinear_pendulum.py \
    sandbox/probe_v130_nonlinear_pendulum_physics.py \
    sandbox/compare_v130_nonlinear_pendulum_CN_full.py \
    tests/test_v130_nonlinear_pendulum_source.py
python3 tests/test_v130_nonlinear_pendulum_source.py
python3 tests/test_v130_finite_dimensional_source.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-pendulum-CN-12s-showcase-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
python3 - "$out" "$image" "$revision" <<'PY'
import json,subprocess,sys
from pathlib import Path
root=Path(sys.argv[1])
(root/"manifest.json").write_text(json.dumps({
    "schema":"scrimp-v130-pendulum-CN-12s-canonical-mpi-refinement-v1",
    "source":subprocess.check_output(
        ["git","rev-parse","HEAD"],text=True
    ).strip(),
    "image":sys.argv[2],
    "image_revision":sys.argv[3],
    "canonical_example":"examples/nonlinear_pendulum.py",
    "original_full_nonlinear_sine_model_unchanged":True,
    "only_physical_time_scheme":"PETSc_TSCN_endpoint_trapezoidal",
    "exact_initial_theta_p_coenergies_without_pseudo":True,
    "original_DAE_DOFs":4,
    "torque_cases":[0.,.1],
    "default_t_f":12.,
    "coarse_dt":.01,
    "refined_dt":.005,
    "preview_steps":8,
    "MPI_coarse":[1,2],
    "MPI_refined":[1],
    "scientific_native_gate_at_start":"PENDING",
    "original_example_produces_figures_not_checker":True,
},indent=2)+"\n")
PY
echo "SCRIMP 1.3.0 original nonlinear PENDULUM CN full 12s qualification"
echo "MUMPS on actual original GetFEM four-DOF nonlinear sin(theta) model."
echo "First: source-generated PNG/SVG/GIF short preflights, FREE and FORCED."
echo "Then: 1200 CN steps MPI1+MPI2, dt=.01, t_f=12, both physical cases."
echo "Finally: 2400 CN steps MPI1, dt=.005, order-two refinement vs RK4."
echo "Inspect complete state/coenergy/energy/power records; no BDF2 claims."
echo "EVIDENCE: $out"

overall=0
for record in \
    "preview_free:free:0:1:0.08:0.01:1:0.003" \
    "preview_forced:forced:0.1:1:0.08:0.01:1:0.003" \
    "free_mpi1:free:0:1:12:0.01:1:0.02" \
    "free_mpi2:free:0:2:12:0.01:0:0.02" \
    "free_refined:free:0:1:12:0.005:0:0.006" \
    "forced_mpi1:forced:0.1:1:12:0.01:1:0.02" \
    "forced_mpi2:forced:0.1:2:12:0.01:0:0.02" \
    "forced_refined:forced:0.1:1:12:0.005:0:0.006"
do
    IFS=: read -r label variant torque mpi tf dt showcase reference_atol <<< "$record"
    mkdir -p "$out/native_$label" "$out/work_$label"
    if [[ "$overall" != 0 ]]; then
        echo "SKIPPED_AFTER_FAILED_NATIVE_RUN" > "$out/$label.exitcode"
        echo "SKIP $label: a prerequisite run failed"
        continue
    fi
    options=()
    if [[ "$showcase" == 1 ]]; then options=(--showcase); fi
    timeout="${SCRIMP_V130_PENDULUM_CN_FULL_TIMEOUT:-2400s}"
    if [[ "$label" == preview_* ]]; then
        timeout="${SCRIMP_V130_PENDULUM_CN_PREVIEW_TIMEOUT:-420s}"
    fi
    echo "RUN $label / original CN: MPI=$mpi torque=$torque t_f=$tf dt=$dt graphics=$showcase"
    rc=0
    "${docker_cmd[@]}" run --rm \
        --user "$(id -u):$(id -g)" \
        --env HOME=/tmp --env MPLBACKEND=Agg \
        --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
        --env PYTHONPATH=/workspace/scrimp \
        --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
        --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
        --env SCRIMP_MATRIX_MEMORY_MODE=compact \
        --env SCRIMP_MUMPS_RESET_POLICY=on_change \
        --env SCRIMP_GETFEM_RHS_MODE=auto \
        --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
        --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
        --volume "$root/examples/nonlinear_pendulum.py:/workspace/scrimp/examples/nonlinear_pendulum.py:ro" \
        --volume "$root/sandbox/probe_v130_nonlinear_pendulum_physics.py:/probe-pend.py:ro" \
        --volume "$out/native_$label:/artifacts:rw" \
        --volume "$out/work_$label:/scrimp-work:rw" \
        --workdir /scrimp-work \
        "$image" timeout -k 10s "$timeout" \
        mpiexec --oversubscribe --bind-to none -n "$mpi" \
        python -u /probe-pend.py \
            --final-time "$tf" --dt "$dt" --angle 1.2 \
            --torque "$torque" --solver mumps \
            --reference-atol "$reference_atol" \
            "${options[@]}" --out /artifacts \
        > "$out/$label.host.log" 2>&1 || rc=$?
    printf '%s\n' "$rc" > "$out/$label.exitcode"
    if [[ "$rc" != 0 ]]; then
        overall=1
        echo "UNQUALIFIED $label original CN exit=$rc" >&2
        tail -n 100 "$out/$label.host.log" >&2
    else
        echo "ORIGINAL NATIVE CN COMPLETED: $label"
    fi

    # Fail before the long campaigns if the canonical user-facing
    # example does not really produce its required scientific graphics.
    if [[ "$overall" == 0 && "$label" == preview_* ]]; then
        if ! python3 - "$out" "$label" "$variant" <<'PY'
import json,sys
from pathlib import Path
from PIL import Image
root=Path(sys.argv[1])
label,variant=sys.argv[2:]
folder=root/f"work_{label}"/"outputs"/f"nonlinear_pendulum_{variant}"
manifest=json.loads((folder/"json"/"manifest.json").read_text())
assert manifest["time_stepper"]=="cn"
assert manifest["original_saved_GetFEM_DAE_states"]==9
assert manifest["GIF_genuine_accepted_frames"]==9
assert len(manifest["figures"])==7
for name in manifest["figures"]:
    f=folder/name.rsplit(".",1)[-1]/name
    assert f.is_file() and f.stat().st_size>1200,str(f)
with Image.open(folder/"gif/pendulum_motion.gif") as image:
    assert image.n_frames==9
print("PASS original canonical nonlinear CN "+variant+" PNG/SVG/GIF preview")
PY
        then
            overall=1
            echo "FAIL: original CN canonical user figure preflight $label" >&2
        fi
    fi
done

compare_rc=0
python3 sandbox/compare_v130_nonlinear_pendulum_CN_full.py "$out" \
    || compare_rc=$?
if [[ "$compare_rc" != 0 ]]; then overall=1; fi
if [[ "$overall" == 0 ]]; then
    echo "PASS_ORIGINAL_PENDULUM_CN_12S_MPI1_MPI2_REFINEMENT_AND_VISUALS"
    echo "FULL FREE FIGURES: $out/work_free_mpi1/outputs/nonlinear_pendulum_free/{png,svg,gif,npz,json}"
    echo "FULL FORCED FIGURES: $out/work_forced_mpi1/outputs/nonlinear_pendulum_forced/{png,svg,gif,npz,json}"
else
    echo "FAIL_OR_UNQUALIFIED_PENDULUM_CN_12S — all native evidence retained" >&2
fi
echo "ARCHIVE: $out.zip"
exit "$overall"
