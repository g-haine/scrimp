#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# #28/#6: original P2 physics; actual PETSc endpoint TSCN at fixed t_f, dt refinement.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

if [[ "${1:-}" == "--compare-only" ]]; then
    [[ "$#" == 2 ]] || { echo "Usage: $0 --compare-only <native-evidence-directory>" >&2; exit 2; }
    python3 sandbox/compare_v130_issue28_full_spinup.py --evidence "$2"
    exit 0
fi
[[ "$#" == 0 ]] || { echo "Usage: $0 [--compare-only <native-evidence-directory>]" >&2; exit 2; }

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
[[ "$revision" == "$expected" ]] || {
    echo "Refusing different/untested Docker image OCI=$revision" >&2
    exit 2
}
for f in scrimp/dphs.py scrimp/port.py scrimp/utils/mesh.py \
    scrimp/utils/dae_initialization.py scrimp/utils/nonlinear_affine_dae.py \
    sandbox/viscous_swe_2d_mpi.py sandbox/probe_v130_issue28_structure.py \
    sandbox/probe_v130_issue28_weak_wall.py \
    sandbox/probe_v130_issue28_one_cn_spinup.py \
    sandbox/probe_v130_issue28_full_spinup.py \
    sandbox/compare_v130_issue28_full_spinup.py; do
    [[ -s "$f" ]] || { echo "Missing source $f" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-full-spinup-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"v130-issue28-full-spinup-v1","source":"%s","image_revision":"%s","original_production_sources_unchanged":true,"fixed_physical_final_time":0.02,"energy_convergence_claimed":false}\n' \
    "$sha" "$revision" > "$out/manifest.json"

echo "#28/#6 original curved P2 SWE complete C1 spinup TSCN dt-refinement and weighted power"
echo "Source=$sha, image=$image OCI=$revision"
echo "NEW p0=0, original physical facet wall; t_f=0.02; n=10/20/40, h=.2/.1, MPI1/2"
echo "Important: sandbox-only initial q0; PETSc endpoint nonlinear mass assumption not guaranteed; no general accuracy claim."
echo "Host-owned output: $out"

for h in 0.2 0.1; do
    label="h${h/./p}"
    for n in 10 20 40; do
        for ranks in 1 2; do
            mkdir -p "$out/work_${label}_n${n}_mpi${ranks}" \
                     "$out/result_${label}_n${n}_mpi${ranks}"
            echo "RUN one physical CN spinup h=$h steps=$n MPI=$ranks"
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
                --volume "$out/work_${label}_n${n}_mpi${ranks}:/scrimp-work:rw" \
                --volume "$out/result_${label}_n${n}_mpi${ranks}:/artifacts:rw" \
                --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
                --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro" \
                --volume "$root/scrimp/utils/mesh.py:/workspace/scrimp/scrimp/utils/mesh.py:ro" \
                --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
                --volume "$root/scrimp/utils/nonlinear_affine_dae.py:/workspace/scrimp/scrimp/utils/nonlinear_affine_dae.py:ro" \
                --volume "$root/sandbox/viscous_swe_2d_mpi.py:/workspace/scrimp/sandbox/viscous_swe_2d_mpi.py:ro" \
                --volume "$root/sandbox/probe_v130_issue28_structure.py:/workspace/scrimp/sandbox/probe_v130_issue28_structure.py:ro" \
                --volume "$root/sandbox/probe_v130_issue28_weak_wall.py:/workspace/scrimp/sandbox/probe_v130_issue28_weak_wall.py:ro" \
                --volume "$root/sandbox/probe_v130_issue28_one_cn_spinup.py:/workspace/scrimp/sandbox/probe_v130_issue28_one_cn_spinup.py:ro" \
                --volume "$root/sandbox/probe_v130_issue28_full_spinup.py:/workspace/scrimp/sandbox/probe_v130_issue28_full_spinup.py:ro" \
                "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_FULL_SPINUP_TIMEOUT || echo 420s)" \
                mpiexec --oversubscribe --bind-to none -n "$ranks" \
                python -u /workspace/scrimp/sandbox/probe_v130_issue28_full_spinup.py \
                    --mesh "$h" --steps "$n" --out /artifacts \
                > "$out/${label}_n${n}_mpi${ranks}.log" 2>&1; then
                echo "FAIL: native one-step CN spinup h=$h steps=$n MPI=$ranks" >&2
                echo "Evidence retained: $out" >&2
                tail -n 140 "$out/${label}_n${n}_mpi${ranks}.log" >&2
                exit 1
            fi
            grep -a 'ISSUE28 FULL SPINUP' \
                "$out/${label}_n${n}_mpi${ranks}.log"
        done
    done
done
python3 sandbox/compare_v130_issue28_full_spinup.py --evidence "$out"
echo "COMPLETE: archived native full C1 physical wall spinup evidence $out"
