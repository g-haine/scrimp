#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# P0 native regression: exact shipped shallow_water weak forms,
# independently FD-check actual default vs historical mass Jacobian.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
    docker_cmd=(docker)
else
    docker_cmd=(sudo docker)
fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
expected="1d6d0c90115213298047b28237a7f38087563a35"
revision="$("${docker_cmd[@]}" image inspect --format \
    '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
if [[ "$revision" != "$expected" ]]; then
    echo "Refusing unqualified frozen runtime $revision" >&2
    exit 2
fi
for f in \
    scrimp/dphs.py \
    scrimp/utils/dae_initialization.py \
    examples/shallow_water.py \
    sandbox/check_v130_shallow_water_mass_jacobian.py \
    sandbox/compare_v130_shallow_water_mass_jacobian.py \
    sandbox/compare_v130_shallow_explicit_tangent.py; do
    [[ -s "$f" ]] || { echo "Missing source $f" >&2; exit 2; }
done
sha="$(git rev-parse HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-shallow-mass-J-${sha:0:8}-XXXXXXXX")"
work="$(mktemp -d "$out/work-XXXXXXXX")"
printf '{"image":"%s","image_source":"%s","diagnostic_source":"%s","model":"original examples/shallow_water.py"}\n' \
    "$image" "$revision" "$sha" > "$out/manifest.json"
echo "SCRIMP v1.3.0 default mass Jacobian native regression"
echo "Original shallow_water weak forms, no main TS, MPI1/2"
echo "A/B: production tangent vs original GetFEM explicit G included diagnostically" 
echo "Compare corrected mass default (env unset) vs legacy on MPI 1/2; never relax thresholds"
echo "Frozen image: $image @ $revision; core: $sha"
echo "Output: $out"
common=( \
    --user "$(id -u):$(id -g)" \
    --workdir /scrimp-work --env HOME=/tmp \
    --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
    --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
    --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
    --env SCRIMP_MATRIX_MEMORY_MODE=compact \
    --env SCRIMP_TRAJECTORY_BACKEND=memory \
    --env SCRIMP_MUMPS_RESET_POLICY=on_change \
    --env SCRIMP_GETFEM_RHS_MODE=auto \
    --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
    --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
    --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
    --volume "$root/sandbox/check_v130_shallow_water_mass_jacobian.py:/workspace/scrimp/sandbox/check_v130_shallow_water_mass_jacobian.py:ro" \
)
# Fail before any expensive GetFEM mesh/assembly if the overlay
# cannot import the *current* production Python core against rc3.
# rc3 itself predates scrimp/utils/dae_initialization.py; mounting
# dphs.py without its companion module is not a valid source overlay.
echo "PREFLIGHT: current SCRIMP core and DAE initialization import in frozen rc3"
if ! "${docker_cmd[@]}" run --rm "${common[@]}" "$image" \
    python -c 'from scrimp.dphs import DPHS; from scrimp.utils.dae_initialization import initialize_linear_algebraic; print("PASS: current DPHS and DAE helper import")' \
    > "$out/import_preflight.log" 2>&1; then
    echo "FAIL: frozen rc3 current source-overlay import" >&2
    tail -n 100 "$out/import_preflight.log" >&2
    echo "Evidence: $out" >&2
    exit 1
fi
cat "$out/import_preflight.log"
for ranks in 1 2; do
    for tangent in production diagnostic_full; do
        for mode in default legacy; do
            echo "RUN SWE MASS J mode=$mode tangent=$tangent MPI=$ranks"
            extra=()
            if [[ "$mode" == legacy ]]; then
                extra=(--env SCRIMP_NONLINEAR_MASS_JACOBIAN=legacy)
            fi
            if ! "${docker_cmd[@]}" run --rm "${common[@]}" "${extra[@]}" "$image" \
                timeout -k 5s 240s \
                mpiexec --oversubscribe --bind-to none -n "$ranks" \
                python -u /workspace/scrimp/sandbox/check_v130_shallow_water_mass_jacobian.py \
                    --mesh .2 --epsilon 1e-6 --state-perturb .01 \
                    --mode "$mode" --tangent-policy "$tangent" \
                    --output-dir /artifacts \
                > "$out/${tangent}_${mode}_${ranks}r.log" 2>&1; then
                echo "FAIL: shallow mass Jacobian mode=$mode tangent=$tangent MPI=$ranks" >&2
                tail -n 100 "$out/${tangent}_${mode}_${ranks}r.log" >&2
                echo "Evidence: $out" >&2
                exit 1
            fi
            grep -a "SHALLOW MASS JACOBIAN" "$out/${tangent}_${mode}_${ranks}r.log" || true
        done
    done
done
echo "COMPARE diagnostic A/B: effect of original explicitly marked gyroscopic G on the tangent"
if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" --workdir /scrimp-work \
    --env HOME=/tmp \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$root/sandbox/compare_v130_shallow_explicit_tangent.py:/workspace/scrimp/sandbox/compare_v130_shallow_explicit_tangent.py:ro" \
    "$image" \
    python -u /workspace/scrimp/sandbox/compare_v130_shallow_explicit_tangent.py \
        /artifacts > "$out/compare_tangent_hypothesis.log" 2>&1; then
    tail -n 110 "$out/compare_tangent_hypothesis.log" >&2
    echo "FAIL: independent explicit-gyro diagnostic, see $out" >&2
    exit 1
fi
cat "$out/compare_tangent_hypothesis.log"
echo "STRICT PRODUCTION GATE: no diagnostic G-tangent override (threshold unchanged)"
if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" --workdir /scrimp-work \
    --env HOME=/tmp \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$root/sandbox/compare_v130_shallow_water_mass_jacobian.py:/workspace/scrimp/sandbox/compare_v130_shallow_water_mass_jacobian.py:ro" \
    "$image" \
    python -u /workspace/scrimp/sandbox/compare_v130_shallow_water_mass_jacobian.py \
        /artifacts > "$out/compare.log" 2>&1; then
    echo "FAIL: strict nonlinear mass-Jacobian qualification" >&2
    tail -n 110 "$out/compare.log" >&2
    echo "Evidence: $out" >&2
    exit 1
fi
cat "$out/compare.log"
echo "COMPLETE: shipped original nonlinear SWE mass Jacobian gate; $out"
echo "This is NOT nonlinear DAE consistent-initialization qualification."
