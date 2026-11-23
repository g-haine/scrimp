#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.
#
# #29 second independent native linear pH-DAE model: original wave example.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

if [[ "${1:-}" == "--compare-only" ]]; then
    [[ "$#" == 2 ]] || { echo "Usage: $0 --compare-only <evidence-dir>" >&2; exit 2; }
    python3 sandbox/compare_v130_checkpoint_linear_dae.py --evidence "$2"
    exit 0
fi
resume_dir=""
if [[ "${1:-}" == "--resume" ]]; then
    [[ "$#" == 2 ]] || { echo "Usage: $0 --resume <existing-evidence-dir>" >&2; exit 2; }
    resume_dir="$2"
else
    [[ "$#" == 0 ]] || {
        echo "Usage: $0 [--resume DIR | --compare-only DIR]" >&2
        exit 2
    }
fi

# Catch all Python syntax/host checkpoint failures before any Docker.
python3 -m py_compile \
    scrimp/dphs.py scrimp/utils/dae_checkpoint.py \
    sandbox/probe_v130_checkpoint_linear_dae.py \
    sandbox/compare_v130_checkpoint_linear_dae.py \
    tests/test_v130_checkpoint_state.py \
    tests/test_v130_checkpoint_public_api.py
python3 tests/test_v130_checkpoint_state.py
python3 tests/test_v130_checkpoint_public_api.py

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
    echo "Refusing image with unqualified original GetFEM/PETSc revision=$revision" >&2
    exit 2
}
files=(
    scrimp/dphs.py scrimp/port.py
    scrimp/utils/dae_checkpoint.py
    scrimp/utils/mesh.py
    scrimp/utils/dae_initialization.py
    scrimp/utils/nonlinear_affine_dae.py
    sandbox/probe_v130_checkpoint_linear_dae.py
    sandbox/compare_v130_checkpoint_linear_dae.py
)
for f in "${files[@]}"; do
    [[ -s "$f" ]] || { echo "Missing native source $f" >&2; exit 2; }
done
mkdir -p "$root/sandbox/outputs"
sha="$(git rev-parse HEAD)"
if [[ -n "$resume_dir" ]]; then
    out="$(realpath -e -- "$resume_dir")"
    [[ "$out" == "$root/sandbox/outputs/"* ]] || {
        echo "Refusing --resume outside original host-owned sandbox/outputs" >&2
        exit 2
    }
    source_sha="$(python3 - "$out/manifest.json" "$revision" <<'PY'
import json,re,sys
from pathlib import Path
doc=json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
assert doc.get("schema")=="scrimp-v130-issue29-linear-DAE-runner-v1"
assert doc.get("original_example")=="wave_2d_linear_DAE"
assert doc.get("cold_state_restart_only") is True
assert doc.get("physical_final_time")==.024
assert doc.get("image_revision")==sys.argv[2]
sha=doc.get("source","")
assert re.fullmatch(r"[a-f0-9]{40}",sha)
print(sha)
PY
)" || { echo "Refusing invalid original linear DAE source manifest" >&2; exit 2; }
    git merge-base --is-ancestor "$source_sha" "$sha" || {
        echo "Refusing --resume: unrelated source revisions" >&2; exit 2;
    }
    if ! git diff --quiet "$source_sha" "$sha" -- \
        scrimp/dphs.py scrimp/port.py \
        scrimp/utils/dae_checkpoint.py \
        scrimp/utils/mesh.py \
        scrimp/utils/dae_initialization.py \
        scrimp/utils/nonlinear_affine_dae.py \
        sandbox/probe_v130_checkpoint_linear_dae.py; then
        echo "Refusing --resume: native physical model/solver/checkpoint changed" >&2
        exit 2
    fi
    echo "RESUME original linear DAE evidence from source=$source_sha"
else
    out="$(mktemp -d "$root/sandbox/outputs/v130-issue29-linear-DAE-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
    printf '{"schema":"scrimp-v130-issue29-linear-DAE-runner-v1","source":"%s","image_revision":"%s","original_example":"wave_2d_linear_DAE","cold_state_restart_only":true,"physical_final_time":0.024}\n' \
        "$sha" "$revision" > "$out/manifest.json"
fi
echo "#29 original linear pH-DAE wave 2D; NO SWE, NO PETSc TS history"
echo "Source=$sha image=$image OCI=$revision"
echo "Continuous 0→0.024 vs full accepted checkpoint at .02 + restart .02→.024"
echo "Native MPI=1,2 sources and 1→1, 1→2, 2→1, 2→2 cold restarts"
echo "Host-owned evidence: $out"

run_case() {
    local phase="$1" ranks="$2" case_name="$3" cp="$4" source="${5:-}"
    local work="$out/work_${case_name}" dest="$out/${case_name}" log="$out/${case_name}.log"
    local verify=(python3 sandbox/compare_v130_checkpoint_linear_dae.py
        --evidence "$out" --verify-case "$case_name"
        --phase "$phase" --ranks "$ranks")
    if [[ "$phase" == resume ]]; then
        verify+=(--source-ranks "$source")
    fi
    if [[ -e "$log" || -e "$dest/linear_report.json" ||
          -e "$dest/linear_states.npz" ]]; then
        if [[ -z "$resume_dir" ]]; then
            echo "Refusing to overwrite previously produced linear DAE evidence" >&2
            exit 1
        fi
        echo "VERIFY original native linear DAE case $case_name"
        "${verify[@]}" || {
            echo "Existing linear DAE native evidence partial or invalid: $case_name" >&2
            exit 1
        }
        echo "SKIP original qualified case $case_name"
        return
    fi
    mkdir -p "$work" "$dest"
    echo "RUN linear DAE phase=$phase MPI=$ranks case=$case_name"
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
        --volume "$root/scrimp/utils/dae_checkpoint.py:/workspace/scrimp/scrimp/utils/dae_checkpoint.py:ro" \
        --volume "$root/scrimp/utils/mesh.py:/workspace/scrimp/scrimp/utils/mesh.py:ro" \
        --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
        --volume "$root/scrimp/utils/nonlinear_affine_dae.py:/workspace/scrimp/scrimp/utils/nonlinear_affine_dae.py:ro" \
        --volume "$root/sandbox/probe_v130_checkpoint_linear_dae.py:/workspace/scrimp/sandbox/probe_v130_checkpoint_linear_dae.py:ro" \
        "$image" timeout -k 5s "${SCRIMP_V130_CHECKPOINT_LINEAR_TIMEOUT:-360s}" \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_checkpoint_linear_dae.py \
            --phase "$phase" --checkpoint "$cp" --out /artifacts \
        > "$log" 2>&1; then
        echo "FAIL native linear DAE $case_name. Logs/results: $out" >&2
        tail -n 130 "$log" >&2
        exit 1
    fi
    if ! grep -aE 'LINEAR DAE PHASE=' "$log"; then
        echo "FAIL: expected real linear DAE success marker missing" >&2
        tail -n 90 "$log" >&2
        exit 1
    fi
    "${verify[@]}" || {
        echo "FAIL: native linear DAE case $case_name not physically valid" >&2
        exit 1
    }
}
for source in 1 2; do
    mkdir -p "$out/wave_source_mpi${source}"
    run_case seed "$source" "wave_seed_mpi${source}" \
        "/checkpoints/wave_source_mpi${source}/accepted_state.npz"
done
for target in 1 2; do
    run_case continuous "$target" "wave_continuous_mpi${target}" \
        "/checkpoints/unused.npz"
done
for source in 1 2; do
    for target in 1 2; do
        run_case resume "$target" "wave_restart_from${source}_to${target}" \
            "/checkpoints/wave_source_mpi${source}/accepted_state.npz" "$source"
    done
done
python3 sandbox/compare_v130_checkpoint_linear_dae.py --evidence "$out"
echo "COMPLETE: archived real native independent linear pH-DAE checkpoint API gate $out"
