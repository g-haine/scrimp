#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# #28: true physical RHS pilot for first-differentiated Y0 at unchanged q0.
# This is a READ-ONLY diagnostic, not a qualified rotating initializer.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
# If the host-only comparison fails, replay saved native evidence safely.
if [[ "${1:-}" == "--compare-only" ]]; then
    if [[ "$#" != 2 ]]; then
        echo "Usage: $0 --compare-only <saved-evidence-directory>" >&2
        exit 2
    fi
    python3 sandbox/compare_v130_issue28_physical_rhs_y0.py --evidence "$2"
    exit 0
fi
if [[ "$#" != 0 ]]; then
    echo "Usage: $0 [--compare-only <saved-evidence-directory>]" >&2
    exit 2
fi

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
    echo "Refusing unqualified immutable Docker image $revision" >&2
    exit 2
fi
for f in scrimp/dphs.py scrimp/port.py scrimp/utils/mesh.py \
    scrimp/utils/dae_initialization.py scrimp/utils/nonlinear_affine_dae.py \
    sandbox/viscous_swe_2d_mpi.py sandbox/probe_v130_issue28_structure.py \
    sandbox/probe_v130_issue28_physical_rhs_y0.py \
    sandbox/compare_v130_issue28_physical_rhs_y0.py; do
    [[ -s "$f" ]] || { echo "Missing original source $f" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-physical-rhs-Y0-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"v130-issue28-physical-rhs-Y0-pilot-v2","source":"%s","image_revision":"%s","geometry_order":2,"original_fixed_h0_p0":true,"physical_Y0_accepted":false,"no_physical_time_steps":true}\n' \
    "$sha" "$revision" > "$out/manifest.json"
echo "#28 original disk P2: physical RHS after regular algebraic projection"
echo "Source=$sha, pinned image=$image OCI=$revision"
echo "Rest and original omega=.2 rotating, MPI1/2, true native IFunction."
echo "IMPORTANT: computed Y0 is DIAGNOSTIC ONLY until Fa compatibility verified."
echo "Host-owned output: $out"

for case in rest rotating; do
    for ranks in 1 2; do
        mkdir -p "$out/work_${case}_${ranks}" "$out/result_${case}_${ranks}"
        echo "RUN physical RHS case=$case MPI=$ranks"
        if ! docker_exec run --rm \
            --user "$(id -u):$(id -g)" --workdir /scrimp-work \
            --env HOME=/tmp --env MPLBACKEND=Agg \
            --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
            --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
            --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
            --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
            --env SCRIMP_MATRIX_MEMORY_MODE=compact \
            --env SCRIMP_MUMPS_RESET_POLICY=on_change \
            --env SCRIMP_GETFEM_RHS_MODE=auto \
            --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
            --volume "$out/work_${case}_${ranks}:/scrimp-work:rw" \
            --volume "$out/result_${case}_${ranks}:/artifacts:rw" \
            --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
            --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro" \
            --volume "$root/scrimp/utils/mesh.py:/workspace/scrimp/scrimp/utils/mesh.py:ro" \
            --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
            --volume "$root/scrimp/utils/nonlinear_affine_dae.py:/workspace/scrimp/scrimp/utils/nonlinear_affine_dae.py:ro" \
            --volume "$root/sandbox/viscous_swe_2d_mpi.py:/workspace/scrimp/sandbox/viscous_swe_2d_mpi.py:ro" \
            --volume "$root/sandbox/probe_v130_issue28_structure.py:/workspace/scrimp/sandbox/probe_v130_issue28_structure.py:ro" \
            --volume "$root/sandbox/probe_v130_issue28_physical_rhs_y0.py:/workspace/scrimp/sandbox/probe_v130_issue28_physical_rhs_y0.py:ro" \
            "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_RHS_TIMEOUT || echo 420s)" \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /workspace/scrimp/sandbox/probe_v130_issue28_physical_rhs_y0.py \
                --case "$case" --out /artifacts \
            > "$out/${case}_mpi${ranks}.log" 2>&1; then
            echo "FAIL: physical RHS diagnostic case=$case ranks=$ranks" >&2
            echo "Evidence retained: $out" >&2
            tail -n 140 "$out/${case}_mpi${ranks}.log" >&2
            exit 1
        fi
        grep -a 'ISSUE28 TRUE PHYSICAL RHS' "$out/${case}_mpi${ranks}.log"
    done
done

python3 sandbox/compare_v130_issue28_physical_rhs_y0.py --evidence "$out"
echo "COMPLETE: native read-only evidence $out"
