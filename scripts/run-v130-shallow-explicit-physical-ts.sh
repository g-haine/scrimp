#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Short genuine nonlinear SWE PETSc TS; independent of previously
# PASSED directional Jacobian FD evidence. Original physical weak
# forms, accepted-step lagging, measured Hamiltonian, MPI parity.
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
    sandbox/probe_v130_shallow_explicit_ts.py \
    sandbox/compare_v130_shallow_explicit_ts.py; do
    [[ -s "$f" ]] || { echo "Missing source $f" >&2; exit 2; }
done
sha="$(git rev-parse HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-shallow-explicit-TS-${sha:0:8}-XXXXXXXX")"
work="$(mktemp -d "$out/work-XXXXXXXX")"
printf '{"image":"%s","image_source":"%s","source":"%s","test":"actual short PETSc TS with original SWE weak forms, explicit G, MPI1/2"}\n' \
    "$image" "$revision" "$sha" > "$out/manifest.json"
echo "SCRIMP v1.3.0 TRUE shallow-water physical TS/explicit G qualification"
echo "Independent of PASSED nonlinear mass-Jacobian directional FD gate."
echo "Original weak forms; published p0 and manufactured nonzero G, MPI 1/2"
echo "Image: $image @ $revision; source: $sha"
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
    --volume "$root/sandbox/probe_v130_shallow_explicit_ts.py:/workspace/scrimp/sandbox/probe_v130_shallow_explicit_ts.py:ro" \
)
echo "PREFLIGHT production SCRIMP imports"
if ! "${docker_cmd[@]}" run --rm "${common[@]}" "$image" \
    python -c 'from scrimp.dphs import DPHS; from scrimp.utils.dae_initialization import initialize_linear_algebraic; print("PASS: production imports")' \
    > "$out/preflight.log" 2>&1; then
    tail -n 70 "$out/preflight.log" >&2
    echo "FAIL: frozen image/source overlay; $out" >&2
    exit 1
fi
cat "$out/preflight.log"
for ranks in 1 2; do
    for case in published manufactured; do
        echo "RUN actual SWE physical TS MPI=$ranks initial=$case"
        if ! "${docker_cmd[@]}" run --rm "${common[@]}" "$image" \
            timeout -k 5s 240s \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /workspace/scrimp/sandbox/probe_v130_shallow_explicit_ts.py \
                --mesh 0.2 --dt 0.0001 --steps 4 --case "$case" \
                --output-dir /artifacts \
            > "$out/${case}_${ranks}r.log" 2>&1; then
            echo "FAIL: actual shallow-water TS MPI=$ranks case=$case" >&2
            tail -n 125 "$out/${case}_${ranks}r.log" >&2
            echo "Evidence: $out" >&2
            exit 1
        fi
        grep -aF "SHALLOW PHYSICAL TS" "$out/${case}_${ranks}r.log" || true
    done
done
echo "COMPARE full-state and Hamiltonian MPI1/2 parity, accepted-step G"
if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" --workdir /scrimp-work \
    --env HOME=/tmp \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$root/sandbox/compare_v130_shallow_explicit_ts.py:/workspace/scrimp/sandbox/compare_v130_shallow_explicit_ts.py:ro" \
    "$image" \
    python -u /workspace/scrimp/sandbox/compare_v130_shallow_explicit_ts.py \
        /artifacts > "$out/compare.log" 2>&1; then
    tail -n 100 "$out/compare.log" >&2
    echo "FAIL: actual physical TS MPI/explicit contract; $out" >&2
    exit 1
fi
cat "$out/compare.log"
echo "COMPLETE: actual shallow-water physical TS qualification; evidence $out"
echo "Hamiltonian drift RECORDED, not asserted to be 0 for first-order lagged G."
