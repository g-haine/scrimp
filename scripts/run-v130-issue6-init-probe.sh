#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Diagnostic-only issue #6 rotating SWE initialization probes.
# Never bypass the physical guard or claim these as release qualification.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 baseline|dt1e-7|rigid-facet|warm-ep|unphysical-matched|from-rest" >&2
    exit 2
fi
probe="$1"
case "$probe" in
    baseline)
        set -- --initial-profile tangent --wall-velocity-profile facet_tangent --algebraic-initial-guess zero --init-dt 1e-9 ;;
    dt1e-7)
        set -- --initial-profile tangent --wall-velocity-profile facet_tangent --algebraic-initial-guess zero --init-dt 1e-7 ;;
    rigid-facet)
        set -- --initial-profile rigid --wall-velocity-profile facet_tangent --algebraic-initial-guess zero --init-dt 1e-9 ;;
    warm-ep)
        set -- --initial-profile rigid --wall-velocity-profile facet_tangent --algebraic-initial-guess ep --init-dt 1e-9 ;;
    unphysical-matched)
        set -- --initial-profile rigid --wall-velocity-profile rigid_trace --algebraic-initial-guess ep_wall --init-dt 1e-9 ;;
    from-rest)
        set -- --initial-profile rigid --wall-velocity-profile facet_tangent --algebraic-initial-guess zero --init-dt 1e-9 --start-from-rest-duration 0.004 ;;
    *)
        echo "Unknown issue #6 initialization diagnostic: $probe" >&2
        exit 2 ;;
esac
docker_exec() {
    if [[ "$(printenv SCRIMP_DOCKER_NO_SUDO || true)" == 1 ]]; then
        docker "$@"
    else
        sudo docker "$@"
    fi
}
image="$(printenv SCRIMP_V130_DAE_IMAGE || echo scrimp-mpi:v1.3.0-rc3)"
expected="1d6d0c90115213298047b28237a7f38087563a35"
revision="$(docker_exec image inspect --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
if [[ "$revision" != "$expected" ]]; then
    echo "Refusing unqualified image source: $revision" >&2
    exit 2
fi
for file in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py sandbox/viscous_swe_2d_mpi.py; do
    [[ -s "$file" ]] || { echo "Missing $file" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse --short=8 HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue6-initial-$probe-$sha-XXXXXXXX")"
mkdir -p "$out/work" "$out/artifacts"
printf '{"schema":"scrimp-v130-issue6-init-probe-v1","probe":"%s","source":"%s","image_revision":"%s","scientific_qualification":false}\n' "$probe" "$sha" "$revision" > "$out/manifest.json"
echo "Issue #6 diagnostic-only probe=$probe source=$sha"
echo "Evidence: $out"
echo "Unchanged nonlinear SWE weak forms, BDF2 solver and 100x guard."
status=0
if docker_exec run --rm --user "$(id -u):$(id -g)" --workdir /scrimp-work \
    --env HOME=/tmp --env MPLBACKEND=Agg \
    --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
    --env SCRIMP_INIT_DYNAMIC_STATE_AUDIT=1 \
    --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
    --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
    --env SCRIMP_MATRIX_MEMORY_MODE=compact \
    --env SCRIMP_TRAJECTORY_BACKEND=memory \
    --env SCRIMP_MUMPS_RESET_POLICY=on_change \
    --env SCRIMP_GETFEM_RHS_MODE=auto \
    --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
    --volume "$out/work:/scrimp-work:rw" \
    --volume "$out/artifacts:/artifacts:rw" \
    --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
    --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro" \
    --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
    --volume "$root/sandbox/viscous_swe_2d_mpi.py:/workspace/scrimp/sandbox/viscous_swe_2d_mpi.py:ro" \
    "$image" \
    timeout -k 5s "$(printenv SCRIMP_ISSUE6_INIT_TIMEOUT || echo 120s)" \
    mpiexec --oversubscribe --bind-to none -n 1 \
    python -u /workspace/scrimp/sandbox/viscous_swe_2d_mpi.py \
        --case open --mesh 0.2 --height-order 2 \
        --final-time 0.02 --save-dt 0.002 --max-step 0.002 \
        --ts-type bdf --adapt-type none --rotation 0.2 \
        --ramp-start 0.003 --ramp-end 0.009 \
        --newton-block-trace --output-dir /artifacts "$@" \
    > "$out/probe.log" 2>&1; then
    echo "DIAGNOSTIC completed; NOT original rotating-data release qualification."
else
    status=$?
    echo "DIAGNOSTIC stopped (native status=$status); guard remains active."
fi
echo "Relevant variable/initializer evidence:"
grep -aE 'Initial consistency:|SCRIMP Newton blocks|SCRIMP initial dynamic audit|SCRIMP algebraic guess|Initial BDF2 state diverged|Initial BDF2 consistency solve failed|RuntimeError:|TS_reason=' "$out/probe.log" | tail -n 65 || true
printf '{"probe":"%s","return_code":%d,"scientific_qualification":false}\n' "$probe" "$status" > "$out/result.json"
echo "Complete log: $out/probe.log"
exit "$status"
