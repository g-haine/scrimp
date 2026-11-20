#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Issue #6 bounded native MPI 1/2/4 smoke. Does not claim long-time physics.
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
    echo "FAIL: unqualified image $revision (expected $expected)" >&2
    exit 2
fi
for p in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py \
         sandbox/viscous_swe_2d_mpi.py \
         sandbox/compare_v130_issue6_viscous_mpi.py; do
    [[ -s "$p" ]] || { echo "Missing required file: $p" >&2; exit 2; }
done
profile="${1:-physical}"
case "$profile" in
    rest) cases=(open); rotation=0; initial_profile=rigid ;;
    physical) cases=(open closed); rotation=0.2; initial_profile=tangent ;;
    *) echo "Usage: $0 [rest|physical]" >&2; exit 2 ;;
esac
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue6-viscous-${profile}-${sha:0:8}-XXXXXXXX")"
tf=0.02
dt=0.002
printf '{"schema":"scrimp-v130-issue6-native-smoke-v1","profile":"%s","source":"%s","image":"%s","image_source":"%s","time_end":%s,"dt":%s}\n' \
    "$profile" "$sha" "$image" "$revision" "$tf" "$dt" > "$out/manifest.json"
echo "Issue #6: short viscous Saint-Venant true MPI 1/2/4: $profile"
echo "Source: $sha / immutable Docker $image at $revision"
echo "Host evidence: $out"
echo "Original physical weak forms. Only bounded time/mesh and ramp; no release PASS assumed."
base=(--user "$(id -u):$(id -g)" --workdir /scrimp-work
    --env HOME=/tmp --env MPLBACKEND=Agg
    --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
    --env SCRIMP_GETFEM_MATRIX_MODE=distributed
    --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo
    --env SCRIMP_MATRIX_MEMORY_MODE=compact
    --env SCRIMP_TRAJECTORY_BACKEND=memory
    --env SCRIMP_MUMPS_RESET_POLICY=on_change
    --env SCRIMP_GETFEM_RHS_MODE=auto
    --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem
    --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro"
    --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro"
    --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro"
    --volume "$root/sandbox/viscous_swe_2d_mpi.py:/workspace/scrimp/sandbox/viscous_swe_2d_mpi.py:ro")
for scenario in "${cases[@]}"; do
    for ranks in 1 2 4; do
        artifacts="$out/${scenario}_${ranks}r"
        work="$out/work_${scenario}_${ranks}r"
        mkdir -p "$artifacts" "$work"
        echo "RUN issue #6 scenario=$scenario MPI=$ranks rotation=$rotation"
        if ! "${docker_cmd[@]}" run --rm "${base[@]}" \
            --volume "$work:/scrimp-work:rw" \
            --volume "$artifacts:/artifacts:rw" "$image" \
            timeout -k 5s "${SCRIMP_ISSUE6_TIMEOUT:-300s}" \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /workspace/scrimp/sandbox/viscous_swe_2d_mpi.py \
                --case "$scenario" --mesh 0.2 --height-order 2 \
                --final-time "$tf" --save-dt "$dt" --max-step "$dt" \
                --ts-type bdf --adapt-type none \
                --rotation "$rotation" --initial-profile "$initial_profile" \
                --wall-velocity-profile facet_tangent \
                --ramp-start 0.003 --ramp-end 0.009 \
                --output-dir /artifacts \
            > "$out/${scenario}_${ranks}r.log" 2>&1; then
            echo "FAIL: native case $scenario MPI=$ranks; log $out/${scenario}_${ranks}r.log" >&2
            tail -n 120 "$out/${scenario}_${ranks}r.log" >&2
            echo "Retained all evidence: $out" >&2
            exit 1
        fi
        grep -aE '"maximum_relative_balance_error"|"maximum_balance_error_over_actual_exchange"|"maximum_relative_volume_balance_error"' \
            "$out/${scenario}_${ranks}r.log" || true
    done
done
compare=(--cases "${cases[@]}" --final-time "$tf")
if [[ "$profile" == "physical" ]]; then compare+=(--require-activity); fi
echo "COMPARE MPI1/2/4 full-state endpoints, h, signed physical powers and each-rank port powers"
if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" --workdir /scrimp-work \
    --env HOME=/tmp \
    --volume "$out:/artifacts:rw" \
    --volume "$root/sandbox/compare_v130_issue6_viscous_mpi.py:/compare.py:ro" \
    "$image" python -u /compare.py /artifacts "${compare[@]}" \
    > "$out/compare.log" 2>&1; then
    echo "FAIL: MPI parity or physical constraints; see $out/compare.log" >&2
    tail -n 120 "$out/compare.log" >&2
    echo "Retained evidence: $out" >&2
    exit 1
fi
cat "$out/compare.log"
echo "PASS: bounded issue #6 $profile smoke, not long-time/paper qualification. Evidence: $out"
