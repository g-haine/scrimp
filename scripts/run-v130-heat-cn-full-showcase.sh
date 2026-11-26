#!/usr/bin/env bash
# SCRIMP v1.3.0: original improved canonical 2-D heat example.
# One short true-VTK renderer preflight, then the full t=5 physical
# MPI1 showcase and MPI2 all-state parity. No manufactured data.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ "$#" == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
    docker_cmd=(docker)
fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
    --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
    echo "Unqualified native GetFEM/PETSc image revision: $revision" >&2
    exit 2
}

python3 -m py_compile examples/heat.py sandbox/compare_v130_heat_cn_full_showcase.py
python3 tests/test_v130_heat_cn_source.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-heat-CN-full-showcase-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
python3 - "$out" "$revision" <<'PY'
import json,subprocess,sys
from pathlib import Path
root=Path(sys.argv[1])
sha=subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip()
(root/"manifest.json").write_text(json.dumps({
    "schema":"scrimp-v130-heat-CN-strict-full-qualification-v1",
    "revision":sha,
    "image_revision":sys.argv[2],
    "canonical_entrypoint":"python -m examples.heat",
    "original_model_bricks_FEM_controls_initial_data":True,
    "user_CLI_init_profile_argument_omitted":True,
    "expected_CN_default_init_profile":"strict",
    "runs":[
        {"label":"preview","MPI":1,"t_f":0.04,"showcase":True},
        {"label":"mpi1","MPI":1,"t_f":5.0,"showcase":True},
        {"label":"mpi2","MPI":2,"t_f":5.0,"showcase":False},
    ],
    "source_only_not_a_native_PASS":True,
},indent=2)+"\n")
PY

echo "SCRIMP 2-D HEAT — full original 500-step CN+strict-algebraic physical qualification"
echo "Default canonical python -m examples.heat, NO override of initializer profile"
echo "Native GetFEM VTU/PNG/SVG/GIF: preflight then MPI1 full; full MPI2 FE/power parity"
echo "Original published examples remain bytewise immutable"
echo "EVIDENCE: $out"

failure=0
for entry in "preview:1:0.04:showcase" "mpi1:1:5:showcase" "mpi2:2:5:no-showcase"; do
    IFS=: read -r label mpi final_time visibility <<< "$entry"
    mkdir -p "$out/native_$label" "$out/work_$label"
    if [[ "$failure" -ne 0 ]]; then
        printf '%s\n' "SKIPPED_AFTER_EARLIER_NATIVE_FAILURE" > "$out/$label.exitcode"
        continue
    fi

    echo "RUN canonical heat CN: $label MPI=$mpi t_f=$final_time visuals=$visibility"
    rc=0
    timeout="${SCRIMP_V130_HEAT_LONG_TIMEOUT:-3600s}"
    if [[ "$label" == "preview" ]]; then
        timeout="${SCRIMP_V130_HEAT_PREVIEW_TIMEOUT:-900s}"
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
        --volume "$root/examples/heat.py:/workspace/scrimp/examples/heat.py:ro" \
        --volume "$out/native_$label:/artifacts:rw" \
        --volume "$out/work_$label:/scrimp-work:rw" \
        --workdir /scrimp-work \
        "$image" timeout -k 10s "$timeout" \
        mpiexec --oversubscribe --bind-to none -n "$mpi" \
        python -u -m examples.heat \
            --time-stepper cn --initializer algebraic \
            --mesh-size .1 --time-step .01 --final-time "$final_time" \
            --solver mumps --check --"$visibility" --output /artifacts \
        > "$out/$label.host.log" 2>&1 || rc=$?
    printf '%s\n' "$rc" > "$out/$label.exitcode"
    if [[ "$rc" -ne 0 ]]; then
        failure=1
        echo "NATIVE CASE UNQUALIFIED: $label, exit=$rc" >&2
        tail -n 75 "$out/$label.host.log" >&2
    else
        echo "NATIVE FINISHED: $label"
    fi
done

compare_rc=0
python3 sandbox/compare_v130_heat_cn_full_showcase.py "$out" || compare_rc=$?
if [[ "$compare_rc" -ne 0 ]]; then failure=1; fi

if [[ "$failure" -eq 0 ]]; then
    echo "PASS_NATIVE_HEAT_CN_500_STEPS_MPI1_MPI2_AND_VISUALS"
    echo "Original user figures: $out/work_mpi1/outputs/showcase/heat/figures"
    echo "Original complete ParaView VTK/PVD: $out/work_mpi1/outputs/pv/T"
else
    echo "FAIL_OR_UNQUALIFIED_HEAT_CN_LONG — inspect original logs and archived traces" >&2
fi
echo "The ZIP holds original FE traces, figures, manifests, logs and HASHED selected raw VTUs."
echo "All 501 unmodified full-series VTUs remain in the host evidence folder."
echo "ARCHIVE: $out.zip"
exit "$failure"
