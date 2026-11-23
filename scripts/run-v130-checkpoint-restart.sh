#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.
#
# #29: native actual GetFEM/PETSc complete-state cold checkpoint/restart.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

if [[ "${1:-}" == "--compare-only" ]]; then
    [[ "$#" == 2 ]] || { echo "Usage: $0 --compare-only <evidence-dir>" >&2; exit 2; }
    python3 sandbox/compare_v130_checkpoint_restart.py --evidence "$2"
    exit 0
fi
[[ "$#" == 0 ]] || { echo "Usage: $0 [--compare-only <evidence-dir>]" >&2; exit 2; }

# Catch Python grammar failures BEFORE any Docker image is launched.
python3 -m py_compile \
    sandbox/v130_checkpoint_state.py \
    sandbox/probe_v130_checkpoint_restart.py \
    sandbox/compare_v130_checkpoint_restart.py \
    tests/test_v130_checkpoint_state.py
python3 tests/test_v130_checkpoint_state.py
python3 tests/test_v130_checkpoint_restart.py

docker_exec() {
    if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
        docker "$@"
    else
        sudo docker "$@"
    fi
}
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
expected="1d6d0c90115213298047b28237a7f38087563a35"
revision="$(docker_exec image inspect --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "$expected" ]] || {
    echo "Refusing unqualified Docker image OCI $revision" >&2
    exit 2
}
files=(
    scrimp/dphs.py scrimp/port.py scrimp/utils/mesh.py
    scrimp/utils/dae_initialization.py scrimp/utils/nonlinear_affine_dae.py
    sandbox/viscous_swe_2d_mpi.py sandbox/probe_v130_issue28_structure.py
    sandbox/probe_v130_issue28_weak_wall.py
    sandbox/probe_v130_issue28_one_cn_spinup.py
    sandbox/v130_checkpoint_state.py sandbox/probe_v130_checkpoint_restart.py
    sandbox/compare_v130_checkpoint_restart.py
)
for f in "${files[@]}"; do
    [[ -s "$f" ]] || { echo "Missing source $f" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue29-cold-restart-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-issue29-runner-v1","source":"%s","image_revision":"%s","original_physics_untouched":true,"cold_state_restart_only":true,"source_spinup_time":0.02,"physical_final_time":0.024}\n' \
    "$sha" "$revision" > "$out/manifest.json"

echo "#29 TRUE native accepted-state GetFEM/PETSc cold restart; NO TS history"
echo "Source=$sha, image=$image OCI=$revision"
echo "P2 physical SWE C1 spinup to 0.02, continuous/restarted to 0.024."
echo "Mesh h=.2/.1; full rank permutations 1->1,1->2,2->1,2->2."
echo "Full host-owned evidence: $out"

run_case() {
    local mesh="$1" phase="$2" ranks="$3" case_dir="$4" checkpoint="$5"
    local work="$out/work_${case_dir}" dest="$out/${case_dir}" log="$out/${case_dir}.log"
    mkdir -p "$work" "$dest"
    echo "RUN mesh=$mesh phase=$phase MPI=$ranks case=$case_dir"
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
        --volume "$work:/scrimp-work:rw" \
        --volume "$dest:/artifacts:rw" \
        --volume "$out:/checkpoints:rw" \
        --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
        --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro" \
        --volume "$root/scrimp/utils/mesh.py:/workspace/scrimp/scrimp/utils/mesh.py:ro" \
        --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
        --volume "$root/scrimp/utils/nonlinear_affine_dae.py:/workspace/scrimp/scrimp/utils/nonlinear_affine_dae.py:ro" \
        --volume "$root/sandbox/viscous_swe_2d_mpi.py:/workspace/scrimp/sandbox/viscous_swe_2d_mpi.py:ro" \
        --volume "$root/sandbox/probe_v130_issue28_structure.py:/workspace/scrimp/sandbox/probe_v130_issue28_structure.py:ro" \
        --volume "$root/sandbox/probe_v130_issue28_weak_wall.py:/workspace/scrimp/sandbox/probe_v130_issue28_weak_wall.py:ro" \
        --volume "$root/sandbox/probe_v130_issue28_one_cn_spinup.py:/workspace/scrimp/sandbox/probe_v130_issue28_one_cn_spinup.py:ro" \
        --volume "$root/sandbox/v130_checkpoint_state.py:/workspace/scrimp/sandbox/v130_checkpoint_state.py:ro" \
        --volume "$root/sandbox/probe_v130_checkpoint_restart.py:/workspace/scrimp/sandbox/probe_v130_checkpoint_restart.py:ro" \
        "$image" timeout -k 5s "${SCRIMP_V130_CHECKPOINT_TIMEOUT:-420s}" \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_checkpoint_restart.py \
            --phase "$phase" --mesh "$mesh" \
            --checkpoint "$checkpoint" --out /artifacts \
        > "$log" 2>&1; then
        echo "FAIL #29 mesh=$mesh phase=$phase MPI=$ranks. Evidence retained: $out" >&2
        tail -n 130 "$log" >&2
        exit 1
    fi
    grep -aE 'CHECKPOINT (PHASE|RESTORE) ' "$log"
}
for mesh in 0.2 0.1; do
    label="h${mesh/./p}"
    for source in 1 2; do
        checkpoint="/checkpoints/${label}_source_mpi${source}/accepted_state.npz"
        mkdir -p "$out/${label}_source_mpi${source}"
        run_case "$mesh" seed "$source" "${label}_seed_mpi${source}" "$checkpoint"
    done
    for target in 1 2; do
        run_case "$mesh" continuous "$target" "${label}_continuous_mpi${target}" \
            "/checkpoints/not-used.npz"
    done
    for source in 1 2; do
        checkpoint="/checkpoints/${label}_source_mpi${source}/accepted_state.npz"
        for target in 1 2; do
            run_case "$mesh" resume "$target" "${label}_restart_from${source}_to${target}" \
                "$checkpoint"
        done
    done
done
python3 sandbox/compare_v130_checkpoint_restart.py --evidence "$out"
echo "COMPLETE: archived original native #29 cold-restart pilot evidence $out"
