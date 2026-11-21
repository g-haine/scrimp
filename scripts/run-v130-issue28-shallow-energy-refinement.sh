#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Issue #28: original shallow SWE port-power defect vs dt, constant
# physical final time, original mesh/data, real production init and BDF2.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

if [[ "${1:-}" == "--compare-only" ]]; then
    if [[ "$#" != 2 ]]; then
        echo "Usage: $0 --compare-only <saved-evidence-directory>" >&2
        exit 2
    fi
    python3 sandbox/compare_v130_issue28_shallow_energy_refinement.py --evidence "$2"
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
    echo "Refusing wrong Docker OCI revision $revision" >&2
    exit 2
fi
for f in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py \
    scrimp/utils/nonlinear_affine_dae.py \
    examples/shallow_water.py \
    sandbox/probe_v130_issue28_shallow_short_physical.py \
    sandbox/compare_v130_issue28_shallow_energy_refinement.py; do
    [[ -s "$f" ]] || { echo "Missing original source $f" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-shallow-energy-refine-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"issue28-shallow-energy-time-refinement-v1","source":"%s","image_revision":"%s","original_physical_model":true,"spatial_mesh_step":0.1,"time_final":0.0004,"true_physical_step_counts":[4,8,16],"exact_discrete_power_conservation_claimed":false}\n' \
    "$sha" "$revision" > "$out/manifest.json"

echo "#28: original rectangular SWE true BDF2 physical energy defect vs dt"
echo "source=$sha immutable image=$image @$revision"
echo "h_mesh=.1, FINAL T=.0004 held constant; dt=.0001/.00005/.000025"
echo "Original h0,p0, original zero boundary power; production nonlinear_affine."
echo "These are SIX native runs (MPI 1/2 x 3 dt); no exact order presupposed."
echo "Evidence: $out"

for entry in "4:0.0001" "8:0.00005" "16:0.000025"; do
    steps="${entry%%:*}"
    dt="${entry#*:}"
    for ranks in 1 2; do
        mkdir -p "$out/work_steps${steps}_${ranks}" "$out/result_steps${steps}_${ranks}"
        echo "RUN original SWE dt=$dt T=.0004 steps=$steps MPI=$ranks"
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
            --volume "$out/work_steps${steps}_${ranks}:/scrimp-work:rw" \
            --volume "$out/result_steps${steps}_${ranks}:/artifacts:rw" \
            --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
            --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro" \
            --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
            --volume "$root/scrimp/utils/nonlinear_affine_dae.py:/workspace/scrimp/scrimp/utils/nonlinear_affine_dae.py:ro" \
            --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
            --volume "$root/sandbox/probe_v130_issue28_shallow_short_physical.py:/workspace/scrimp/sandbox/probe_v130_issue28_shallow_short_physical.py:ro" \
            "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_REFINEMENT_TIMEOUT || echo 420s)" \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /workspace/scrimp/sandbox/probe_v130_issue28_shallow_short_physical.py \
                --out /artifacts --dt "$dt" --final-time 0.0004 \
            > "$out/steps${steps}_mpi${ranks}.log" 2>&1; then
            echo "FAIL: original SWE time refinement steps=$steps MPI=$ranks" >&2
            echo "All logs/JSON already computed are retained: $out" >&2
            tail -n 125 "$out/steps${steps}_mpi${ranks}.log" >&2
            exit 1
        fi
        grep -a 'ISSUE28 SHALLOW SHORT PHYSICAL' "$out/steps${steps}_mpi${ranks}.log"
    done
done
python3 sandbox/compare_v130_issue28_shallow_energy_refinement.py --evidence "$out"
echo "COMPLETE: fixed-T native time-refinement evidence preserved $out"
