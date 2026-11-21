#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# #28 genuine PETSc BDF2 physical TS following production
# nonlinear_affine DAE projection, original published shallow_water.py.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
if [[ "$(printenv SCRIMP_DOCKER_NO_SUDO || true)" == 1 ]]; then
    docker_cmd=(docker)
else
    docker_cmd=(sudo docker)
fi
if [[ "${1:-}" == "--compare-only" ]]; then
    if [[ "$#" != 2 ]]; then
        echo "Usage: $0 --compare-only <saved-evidence-directory>" >&2
        exit 2
    fi
    python3 sandbox/compare_v130_issue28_production_short_physical.py --evidence "$2"
    exit 0
fi
if [[ "$#" != 0 ]]; then
    echo "Usage: $0 [--compare-only <saved-evidence-directory>]" >&2
    exit 2
fi
image="$(printenv SCRIMP_V130_DAE_IMAGE || echo scrimp-mpi:v1.3.0-rc3)"
expected="1d6d0c90115213298047b28237a7f38087563a35"
revision="$("${docker_cmd[@]}" image inspect --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
if [[ "$revision" != "$expected" ]]; then
    echo "Refusing unqualified Docker image revision $revision" >&2
    exit 2
fi
for f in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py \
    scrimp/utils/nonlinear_affine_dae.py examples/shallow_water.py \
    sandbox/probe_v130_shallow_explicit_ts.py \
    sandbox/probe_v130_issue28_production_short_physical.py \
    sandbox/compare_v130_issue28_production_short_physical.py; do
    [[ -s "$f" ]] || { echo "Missing source $f" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-production-short-ts-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"issue28-production-short-physical-v1","source":"%s","image_revision":"%s","original_physical_model":true,"initialization_strategy":"nonlinear_affine","true_physical_BDF2_steps":4,"exact_power_budget_qualified":false}\n' "$sha" "$revision" > "$out/manifest.json"
echo "#28 production nonlinear_affine -> 4 genuine physical BDF2 steps"
echo "Source=$sha Docker=$image pinned revision=$revision"
echo "Original rectangle .1, published h0,p0, original zero boundary controls."
echo "Signed Hamiltonian drift is MEASURED, not assumed exactly conserved."
echo "Host evidence: $out"
for ranks in 1 2; do
    mkdir -p "$out/work_${ranks}" "$out/result_${ranks}"
    echo "RUN actual production initialized original shallow-water TS MPI=$ranks"
    if ! "${docker_cmd[@]}" run --rm \
        --user "$(id -u):$(id -g)" --workdir /scrimp-work \
        --env HOME=/tmp --env MPLBACKEND=Agg \
        --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
        --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
        --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
        --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
        --env SCRIMP_MATRIX_MEMORY_MODE=compact \
        --env SCRIMP_TRAJECTORY_BACKEND=memory \
        --env SCRIMP_MUMPS_RESET_POLICY=on_change \
        --env SCRIMP_GETFEM_RHS_MODE=auto \
        --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
        --volume "$out/work_${ranks}:/scrimp-work:rw" \
        --volume "$out/result_${ranks}:/artifacts:rw" \
        --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
        --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro" \
        --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
        --volume "$root/scrimp/utils/nonlinear_affine_dae.py:/workspace/scrimp/scrimp/utils/nonlinear_affine_dae.py:ro" \
        --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
        --volume "$root/sandbox/probe_v130_shallow_explicit_ts.py:/workspace/scrimp/sandbox/probe_v130_shallow_explicit_ts.py:ro" \
        --volume "$root/sandbox/probe_v130_issue28_production_short_physical.py:/workspace/scrimp/sandbox/probe_v130_issue28_production_short_physical.py:ro" \
        "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_TIMEOUT || echo 240s)" \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_issue28_production_short_physical.py \
            --out /artifacts --dt .0001 --steps 4 \
        > "$out/mpi_${ranks}.log" 2>&1; then
        echo "FAIL: native physical TS MPI=$ranks; evidence retained $out" >&2
        tail -n 150 "$out/mpi_${ranks}.log" >&2
        exit 1
    fi
    grep -a 'ISSUE28 PRODUCTION AFFINE REAL TS' "$out/mpi_${ranks}.log"
done
python3 sandbox/compare_v130_issue28_production_short_physical.py --evidence "$out"
echo "COMPLETE: native evidence retained $out"
