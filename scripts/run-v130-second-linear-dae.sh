#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Second independent pH-DAE: dissipative 2D wave + boundary control,
# production DPHS.init_step(algebraic) in frozen qualified rc3 runtime.
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
    echo "Refusing unqualified image source $revision" >&2
    exit 2
fi
for file in \
    scrimp/dphs.py \
    scrimp/utils/dae_initialization.py \
    sandbox/probe_v130_dissipative_core_dae.py \
    sandbox/compare_v130_dissipative_core_dae.py \
    examples/wave_dissipative.py; do
    [[ -s "$file" ]] || { echo "Missing $file" >&2; exit 2; }
done
sha="$(git rev-parse HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-second-linear-dae-${sha:0:8}-XXXXXXXX")"
work="$(mktemp -d "$out/work-XXXXXXXX")"
echo "SCRIMP independent second linear DAE (dissipative wave)"
echo "Frozen native image $image / source $revision; core commit $sha"
echo "Original dissipative bricks and port, manufactured compatible p0"
echo "h=0.2, fixed CN dt=.002 through 10 steps, MPI=1/2"
echo "Check *both* signed boundary work and positive internal damping work."
echo "No Docker rebuild, no production example or default-mode change."
echo "Output: $out"
printf '{"schema":"scrimp-v130-second-linear-dae-run-v1","image_source":"%s","core_source":"%s","model":"examples/wave_dissipative.py","mpi":[1,2]}\n' \
    "$revision" "$sha" > "$out/manifest.json"
common=( \
    --user "$(id -u):$(id -g)" \
    --workdir /scrimp-work \
    --env HOME=/tmp --env MPLBACKEND=Agg \
    --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
    --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
    --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
    --env SCRIMP_MATRIX_MEMORY_MODE=compact \
    --env SCRIMP_TRAJECTORY_BACKEND=memory \
    --env SCRIMP_MUMPS_RESET_POLICY=on_change \
    --env SCRIMP_INIT_FACTOR_SOLVER_TYPE=mumps \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
    --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
    --volume "$root/examples/wave_dissipative.py:/workspace/scrimp/examples/wave_dissipative.py:ro" \
    --volume "$root/sandbox/probe_v130_dissipative_core_dae.py:/workspace/scrimp/sandbox/probe_v130_dissipative_core_dae.py:ro" \
)
for ranks in 1 2; do
    echo "RUN SECOND LINEAR DAE MPI=$ranks"
    if ! "${docker_cmd[@]}" run --rm "${common[@]}" "$image" \
        timeout -k 5s 240s \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_dissipative_core_dae.py \
            --mesh 0.2 --dt 0.002 --steps 10 \
            --output-dir /artifacts \
        > "$out/dissipative_${ranks}r.log" 2>&1; then
        echo "FAIL: independent dissipative DAE MPI=$ranks" >&2
        tail -n 105 "$out/dissipative_${ranks}r.log" >&2
        echo "Evidence: $out" >&2
        exit 1
    fi
    grep -aF "SECOND LINEAR DAE MPI=$ranks" "$out/dissipative_${ranks}r.log" || true
done
if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" \
    --workdir /scrimp-work --env HOME=/tmp \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$root/sandbox/compare_v130_dissipative_core_dae.py:/workspace/scrimp/sandbox/compare_v130_dissipative_core_dae.py:ro" \
    "$image" \
    python -u /workspace/scrimp/sandbox/compare_v130_dissipative_core_dae.py \
        /artifacts \
    > "$out/compare.log" 2>&1; then
    echo "FAIL: independent dissipative DAE comparator" >&2
    tail -n 105 "$out/compare.log" >&2
    echo "Evidence: $out" >&2
    exit 1
fi
cat "$out/compare.log"
echo "COMPLETE: second linear DAE native qualification; evidence $out"
echo "DEFAULT REMAINS pseudo until independent default migration gate."
