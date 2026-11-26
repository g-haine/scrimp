#!/usr/bin/env bash
# SCRIMP v1.3.0 original mass-spring-damper full 8-second qualification.
# This calls the CANONICAL example for both physics and visualization.
# No analysis script draws figures. Native MUMPS, exact q0/p0/e0,
# 4 true lumped GetFEM degrees of freedom, 801 saved states, MPI1/2.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
    --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == 1d6d0c90115213298047b28237a7f38087563a35 ]] || {
    echo "Unqualified GetFEM/PETSc Docker revision: $revision" >&2
    exit 2
}
python3 -m py_compile \
    examples/mass_spring_damper.py \
    sandbox/probe_v130_mass_spring_physics.py \
    sandbox/compare_v130_mass_spring_full_showcase.py
python3 tests/test_v130_finite_dimensional_source.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-mass-spring-8s-showcase-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
python3 - "$out" "$image" "$revision" <<'PY'
import json,subprocess,sys
from pathlib import Path
root=Path(sys.argv[1])
revision=subprocess.check_output(
    ["git","rev-parse","HEAD"],text=True
).strip()
(root/"manifest.json").write_text(json.dumps({
    "schema":"scrimp-v130-mass-spring-full-canonical-v1",
    "source":revision,
    "image":sys.argv[2],
    "image_revision":sys.argv[3],
    "canonical_entrypoint":"examples.mass_spring_damper.mass_spring_damper",
    "original_model_m_k_c_force_bricks_preserved":True,
    "original_BDF2_dt":.01,
    "original_final_time":8.,
    "initial_q0":.4,
    "initial_p0":0.,
    "initial_e_q0":2.,
    "initial_e_p0":0.,
    "initialization":"exact_without_pseudo_TS",
    "original_MUMPS_MPI":[1,2],
    "staged_not_yet_a_native_PASS":True,
    "runs":[
        {"label":"preview","MPI":1,"t_f":.08,"showcase":True},
        {"label":"mpi1","MPI":1,"t_f":8.,"showcase":True},
        {"label":"mpi2","MPI":2,"t_f":8.,"showcase":False},
    ]
},indent=2)+"\n")
PY
echo "SCRIMP v1.3 original MASS-SPRING: full 800-step BDF2 qualification"
echo "Original exact initial state; 4 authentic lumped GetFEM DOFs"
echo "Short visual preflight -> 8s MPI1 canonical figures -> 8s MPI2 DAE parity"
echo "Independent ODE RK4, constitutive laws, exact BDF2 identity for k>=2"
echo "Canonical example itself generates the PNG/SVG/GIF, not the probe."
echo "EVIDENCE: $out"

overall=0
for run in "preview:1:0.08:1" "mpi1:1:8:1" "mpi2:2:8:0"; do
    IFS=: read -r label mpi tf showcase <<< "$run"
    mkdir -p "$out/native_$label" "$out/work_$label"
    if [[ "$overall" != 0 ]]; then
        echo SKIPPED_DUE_TO_PREVIOUS_FAILURE > "$out/$label.exitcode"
        continue
    fi
    echo "RUN ORIGINAL BDF2 $label MPI=$mpi t_f=$tf canonical_plots=$showcase"
    options=()
    if [[ "$showcase" == 1 ]]; then options=(--showcase); fi
    rc=0
    timeout="${SCRIMP_V130_MSD_FULL_TIMEOUT:-1200s}"
    if [[ "$label" == preview ]]; then
        timeout="${SCRIMP_V130_MSD_PREVIEW_TIMEOUT:-300s}"
    fi
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
        --volume "$root/examples/mass_spring_damper.py:/workspace/scrimp/examples/mass_spring_damper.py:ro" \
        --volume "$root/sandbox/probe_v130_mass_spring_physics.py:/probe-msd.py:ro" \
        --volume "$out/native_$label:/artifacts:rw" \
        --volume "$out/work_$label:/scrimp-work:rw" \
        --workdir /scrimp-work \
        "$image" timeout -k 10s "$timeout" \
        mpiexec --oversubscribe --bind-to none -n "$mpi" \
        python -u /probe-msd.py \
            --final-time "$tf" --time-step .01 \
            --solver mumps \
            "${options[@]}" --out /artifacts \
        > "$out/$label.host.log" 2>&1 || rc=$?
    printf '%s\n' "$rc" > "$out/$label.exitcode"
    if [[ "$rc" != 0 ]]; then
        overall=1
        echo "UNQUALIFIED native original $label, exit $rc" >&2
        tail -n 90 "$out/$label.host.log" >&2
    else
        echo "ORIGINAL NATIVE COMPLETE: $label"
    fi

    # Catch broken canonical user plots on four steps BEFORE paying
    # for two full time integrations. Raw field traces still retained.
    if [[ "$label" == preview && "$overall" == 0 ]]; then
        if ! python3 - "$out" <<'PY'
import json,sys
from pathlib import Path
from PIL import Image
root=Path(sys.argv[1])
fig=root/"work_preview/outputs/showcase/mass_spring_damper/figures"
manifest=fig.parent/"manifest.json"
assert manifest.is_file(),"Canonical short preview manifest absent"
info=json.loads(manifest.read_text())
assert info["native_GetFEM_saved_states"]==9
assert info["original_accepted_states_used_for_GIF"]==9
assert len(info["figures"])==7
for name in info["figures"]:
    file=fig/name
    assert file.is_file() and file.stat().st_size>1200, str(file)
with Image.open(fig/"mass_spring_motion.gif") as gif:
    assert gif.n_frames==9
print("PASS: original canonical mass-spring short PNG/SVG/GIF preflight")
PY
        then
            overall=1
            echo "FAIL: canonical 4-DOF lumped example visual preflight" >&2
        fi
    fi
done

compare_rc=0
python3 sandbox/compare_v130_mass_spring_full_showcase.py "$out" || compare_rc=$?
if [[ "$compare_rc" != 0 ]]; then overall=1; fi
if [[ "$overall" == 0 ]]; then
    echo "PASS_NATIVE_ORIGINAL_MASS_SPRING_800STEP_MPI1_MPI2_AND_VISUALS"
    echo "TRUE EXAMPLE FIGURES: $out/work_mpi1/outputs/showcase/mass_spring_damper/figures"
else
    echo "FAIL_OR_UNQUALIFIED_MASS_SPRING_8S — all evidence retained" >&2
fi
echo "Original lumped GetFEM states and source-generated figures:"
echo "ARCHIVE: $out.zip"
exit "$overall"
