#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Qualify the TRUE implicit default DPHS.init_step(algebraic) on frozen rc3.
# No Docker image rebuild; mount ONLY changed core Python files and
# independent diagnostic harness into qualified image, read-only.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
    docker_cmd=(docker)
else
    docker_cmd=(sudo docker)
fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect --format \
    '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
expected="1d6d0c90115213298047b28237a7f38087563a35"
if [[ "$revision" != "$expected" ]]; then
    echo "Refusing unqualified native rc3 source: $revision" >&2
    exit 2
fi
reference="${SCRIMP_V130_DAE_REFERENCE_DIR:-$root/sandbox/outputs/v130-petsc-dae-mpi-h0p2-4d3b4d0c-gb0qNKB8}"
legacy_reference="${SCRIMP_V130_LEGACY_REFERENCE_DIR:-$root/sandbox/outputs/v130-wave-init-power-h0p2-62e31d9b-zDa95p5J}"
for p in \
    scrimp/dphs.py \
    scrimp/utils/dae_initialization.py \
    sandbox/probe_v130_core_dae_initial.py \
    sandbox/compare_v130_core_dae_initial.py \
    sandbox/audit_v130_wave_initializers.py \
    sandbox/probe_v130_wave_fixed_dynamics.py \
    sandbox/audit_v130_wave_dae_initial.py \
    sandbox/solve_v130_dae_algebraic_block.py \
    sandbox/check_acoustic_energy_mpi.py \
    sandbox/probe_v130_dissipative_core_dae.py \
    sandbox/compare_v130_dissipative_core_dae.py \
    sandbox/verify_v130_default_promotion.py \
    examples/wave_dissipative.py; do
    [[ -s "$p" ]] || { echo "Missing required source $p" >&2; exit 2; }
done
for path in \
    "$reference/distributed_polynomial_compatible_1r.npz" \
    "$reference/distributed_polynomial_compatible_2r.npz" \
    "$reference/distributed_polynomial_compatible_1r.json" \
    "$reference/distributed_polynomial_compatible_2r.json" \
    "$legacy_reference/wave_pseudo_1r.npz"; do
    [[ -s "$path" ]] || {
        echo "Missing independent frozen physical baseline: $path" >&2
        echo "Cannot claim default promotion without verifying reference." >&2
        exit 2
    }
done
head="$(git rev-parse HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-algebraic-default-h0p2-${head:0:8}-XXXXXXXX")"
work="$(mktemp -d "$out/work-XXXXXXXX")"
echo "SCRIMP v1.3.0 true implicit algebraic DAE DEFAULT qualification"
echo "Image (no rebuild): $image @ $revision"
echo "Source under test: $head"
echo "Frozen prototype: $reference"
echo "Frozen default-pseudo physics: $legacy_reference"
echo "No-option algebraic conservative/dissipative x MPI1/2 + Gaussian reject + explicit pseudo/BDF2"
echo "Output: $out"
printf '{"schema":"scrimp-v130-implicit-default-dae-native-v1","native_image":"%s","source":"%s"}\n' \
    "$revision" "$head" > "$out/manifest.json"
common=( \
    --user "$(id -u):$(id -g)" \
    --workdir /scrimp-work \
    --env HOME=/tmp --env MPLBACKEND=Agg \
    --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
    --env SCRIMP_IMAGE_REVISION="$revision" \
    --env SCRIMP_BENCH_GIT_SHA="$head" \
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
    --volume "$root/sandbox/probe_v130_core_dae_initial.py:/workspace/scrimp/sandbox/probe_v130_core_dae_initial.py:ro" \
    --volume "$root/sandbox/audit_v130_wave_initializers.py:/workspace/scrimp/sandbox/audit_v130_wave_initializers.py:ro" \
    --volume "$root/sandbox/probe_v130_wave_fixed_dynamics.py:/workspace/scrimp/sandbox/probe_v130_wave_fixed_dynamics.py:ro" \
    --volume "$root/sandbox/audit_v130_wave_dae_initial.py:/workspace/scrimp/sandbox/audit_v130_wave_dae_initial.py:ro" \
    --volume "$root/sandbox/solve_v130_dae_algebraic_block.py:/workspace/scrimp/sandbox/solve_v130_dae_algebraic_block.py:ro" \
    --volume "$root/sandbox/check_acoustic_energy_mpi.py:/workspace/scrimp/sandbox/check_acoustic_energy_mpi.py:ro" \
    --volume "$root/examples/wave_dissipative.py:/workspace/scrimp/examples/wave_dissipative.py:ro" \
    --volume "$root/sandbox/probe_v130_dissipative_core_dae.py:/workspace/scrimp/sandbox/probe_v130_dissipative_core_dae.py:ro" \
)
run_case() {
    local mode="$1" scenario="$2" ranks="$3"
    local -a extra_args=()
    if [[ "$mode" == "algebraic" ]]; then
        extra_args=(--use-default)
    fi
    echo "RUN DEFAULT CORE mode=$mode scenario=$scenario MPI=$ranks"
    if ! "${docker_cmd[@]}" run --rm "${common[@]}" "$image" \
        timeout -k 5s 240s \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_core_dae_initial.py \
            --mode "$mode" --scenario "$scenario" \
            --mesh 0.2 --dt 0.002 --steps 10 \
            --tiny-bdf-dt 0.000001 --output-dir /artifacts \
            "${extra_args[@]}" \
            > "$out/${mode}_${scenario}_${ranks}r.log" 2>&1; then
        echo "FAIL: core DAE mode=$mode scenario=$scenario MPI=$ranks" >&2
        tail -n 120 "$out/${mode}_${scenario}_${ranks}r.log" >&2
        echo "Evidence retained: $out" >&2
        exit 1
    fi
    grep -aF "CORE DAE $mode MPI=$ranks" \
        "$out/${mode}_${scenario}_${ranks}r.log" || true
}
for ranks in 1 2; do
    run_case algebraic polynomial_compatible "$ranks"
done
for ranks in 1 2; do
    run_case algebraic published "$ranks"
done
run_case pseudo polynomial_compatible 1
run_case bdf2 polynomial_compatible 1
for ranks in 1 2; do
    echo "RUN DEFAULT DISSIPATIVE MPI=$ranks"
    if ! "${docker_cmd[@]}" run --rm "${common[@]}" "$image" \
        timeout -k 5s 240s \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_dissipative_core_dae.py \
            --mesh 0.2 --dt 0.002 --steps 10 \
            --output-dir /artifacts --use-default \
        > "$out/dissipative_default_${ranks}r.log" 2>&1; then
        echo "FAIL: default dissipative MPI=$ranks" >&2
        tail -n 120 "$out/dissipative_default_${ranks}r.log" >&2
        echo "Host evidence retained: $out" >&2
        exit 1
    fi
    grep -aF "SECOND LINEAR DAE MPI=$ranks" \
        "$out/dissipative_default_${ranks}r.log" || true
done
echo "COMPARE real default core against BOTH frozen independent scientific references"
if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" \
    --workdir /scrimp-work \
    --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp/sandbox \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$reference:/frozen-dae:ro" \
    --volume "$legacy_reference:/frozen-legacy:ro" \
    --volume "$root/sandbox/compare_v130_core_dae_initial.py:/workspace/scrimp/sandbox/compare_v130_core_dae_initial.py:ro" \
    "$image" \
    python -u /workspace/scrimp/sandbox/compare_v130_core_dae_initial.py \
        --out /artifacts --reference /frozen-dae \
        --legacy-reference /frozen-legacy \
    > "$out/compare.log" 2>&1; then
    tail -n 120 "$out/compare.log" >&2
    echo "FAIL: strict opt-in core physics acceptance, see $out" >&2
    exit 1
fi
cat "$out/compare.log"
echo "COMPARE default dissipative CN parity and signed physical damping work"
if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" --workdir /scrimp-work \
    --env HOME=/tmp \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$root/sandbox/compare_v130_dissipative_core_dae.py:/workspace/scrimp/sandbox/compare_v130_dissipative_core_dae.py:ro" \
    "$image" \
    python -u /workspace/scrimp/sandbox/compare_v130_dissipative_core_dae.py \
        /artifacts > "$out/compare_dissipative.log" 2>&1; then
    tail -n 100 "$out/compare_dissipative.log" >&2
    echo "FAIL: default dissipative physical MPI gate, see $out" >&2
    exit 1
fi
cat "$out/compare_dissipative.log"
echo "VERIFY neither DAE model specified init_step_ts_type"
if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" --workdir /scrimp-work \
    --env HOME=/tmp \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$root/sandbox/verify_v130_default_promotion.py:/workspace/scrimp/sandbox/verify_v130_default_promotion.py:ro" \
    "$image" \
    python -u /workspace/scrimp/sandbox/verify_v130_default_promotion.py \
        /artifacts > "$out/verify_default.log" 2>&1; then
    tail -n 110 "$out/verify_default.log" >&2
    echo "FAIL: true implicit-default qualification, see $out" >&2
    exit 1
fi
cat "$out/verify_default.log"
echo "COMPLETE: staged implicit algebraic default native gate; evidence $out"
echo "Do not tag or rebuild release image until broader compatibility/memory gates are checked."
