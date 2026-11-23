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
resume_dir=""
if [[ "${1:-}" == "--resume" ]]; then
    [[ "$#" == 2 ]] || {
        echo "Usage: $0 --resume <existing-host-evidence-directory>" >&2
        exit 2
    }
    resume_dir="$2"
else
    [[ "$#" == 0 ]] || {
        echo "Usage: $0 [--resume <existing-evidence-dir> | --compare-only <evidence-dir>]" >&2
        exit 2
    }
fi

# Catch Python grammar failures BEFORE any Docker image is launched.
python3 -m py_compile \
    sandbox/v130_checkpoint_state.py \
    sandbox/probe_v130_checkpoint_restart.py \
    sandbox/compare_v130_checkpoint_restart.py \
    tests/test_v130_checkpoint_state.py \
    tests/test_v130_checkpoint_restart.py \
    sandbox/validate_v130_checkpoint_case.py
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
    sandbox/validate_v130_checkpoint_case.py
)
for f in "${files[@]}"; do
    [[ -s "$f" ]] || { echo "Missing source $f" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
if [[ -n "$resume_dir" ]]; then
    out="$(realpath -e -- "$resume_dir")"
    [[ "$out" == "$root/sandbox/outputs/"* ]] || {
        echo "Refusing --resume outside the host-owned sandbox/outputs directory" >&2
        exit 2
    }
    [[ -s "$out/manifest.json" ]] || {
        echo "Refusing --resume: missing original evidence manifest" >&2
        exit 2
    }
    # Source commit is recorded by the original native evidence. It
    # need not equal HEAD after this launcher is fixed, but the physics,
    # solver, checkpoint schema and native probe must be IDENTICAL.
    source_sha="$(python3 - "$out/manifest.json" "$revision" <<'PY'
import json, re, sys
from pathlib import Path
p = Path(sys.argv[1])
doc = json.loads(p.read_text(encoding="utf-8"))
assert doc.get("schema") == "scrimp-v130-issue29-runner-v1"
assert doc.get("original_physics_untouched") is True
assert doc.get("cold_state_restart_only") is True
assert doc.get("source_spinup_time") == 0.02
assert doc.get("physical_final_time") == 0.024
assert doc.get("image_revision") == sys.argv[2]
source = doc.get("source", "")
assert re.fullmatch(r"[a-f0-9]{40}", source)
print(source)
PY
)" || { echo "Invalid original evidence manifest" >&2; exit 2; }
    git merge-base --is-ancestor "$source_sha" "$sha" || {
        echo "Refusing --resume: original native source is not a HEAD ancestor" >&2
        exit 2
    }
    if ! git diff --quiet "$source_sha" "$sha" -- \
        scrimp/dphs.py scrimp/port.py scrimp/utils/mesh.py \
        scrimp/utils/dae_initialization.py scrimp/utils/nonlinear_affine_dae.py \
        sandbox/viscous_swe_2d_mpi.py \
        sandbox/probe_v130_issue28_structure.py \
        sandbox/probe_v130_issue28_weak_wall.py \
        sandbox/probe_v130_issue28_one_cn_spinup.py \
        sandbox/v130_checkpoint_state.py \
        sandbox/probe_v130_checkpoint_restart.py; then
        echo "Refusing --resume: production physics, solver or native probe changed" >&2
        exit 2
    fi
    echo "RESUMING immutable native evidence from source=$source_sha"
else
    out="$(mktemp -d "$root/sandbox/outputs/v130-issue29-cold-restart-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
    printf '{"schema":"scrimp-v130-issue29-runner-v1","source":"%s","image_revision":"%s","original_physics_untouched":true,"cold_state_restart_only":true,"source_spinup_time":0.02,"physical_final_time":0.024}\n' \
        "$sha" "$revision" > "$out/manifest.json"
fi

echo "#29 TRUE native accepted-state GetFEM/PETSc cold restart; NO TS history"
echo "Source=$sha, image=$image OCI=$revision"
echo "P2 physical SWE C1 spinup to 0.02, continuous/restarted to 0.024."
echo "Mesh h=.2/.1; full rank permutations 1->1,1->2,2->1,2->2."
echo "Full host-owned evidence: $out"

run_case() {
    local mesh="$1" phase="$2" ranks="$3" case_dir="$4" checkpoint="$5"
    local work="$out/work_${case_dir}" dest="$out/${case_dir}" log="$out/${case_dir}.log"
    local validate=(python3 sandbox/validate_v130_checkpoint_case.py
        --evidence "$out" --mesh "$mesh" --phase "$phase" --ranks "$ranks")
    if [[ "$phase" == resume ]]; then
        # h0p2_restart_from1_to2 -> source_rank=1
        local original_rank="${case_dir#*_restart_from}"
        original_rank="${original_rank%%_to*}"
        validate+=(--source-ranks "$original_rank")
    fi
    if [[ -e "$log" || -e "$dest/segment_report.json" ||
          -e "$dest/segment_states.npz" ]]; then
        if [[ -z "$resume_dir" ]]; then
            echo "Refusing to overwrite existing native case evidence: $case_dir" >&2
            exit 1
        fi
        echo "VERIFY existing native case before skipping $case_dir"
        "${validate[@]}" || {
            echo "FAIL: partial or invalid case cannot be skipped: $case_dir" >&2
            exit 1
        }
        echo "SKIP valid original native evidence $case_dir"
        return 0
    fi
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
    # Real native success is CHECKPOINT PHASE=seed (EQUALS).
    # The old regex required a SPACE after PHASE and silently killed
    # the whole runner under set -e, even after a successful TS step.
    if ! grep -aE 'CHECKPOINT (PHASE=|RESTORE )' "$log"; then
        echo "FAIL: native success marker missing from $log" >&2
        tail -n 60 "$log" >&2
        exit 1
    fi
    "${validate[@]}" || {
        echo "FAIL: native case evidence not valid: $case_dir" >&2
        exit 1
    }
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
