#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.
#
# Standalone physical-generation pilot, not a general stationary solver.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

usage() {
    echo "Usage: $0 seed|resume 0.2|0.1 1|2|4 <absolute-final-time> [input-checkpoint.npz]" >&2
}
[[ $# -ge 4 && $# -le 5 ]] || { usage; exit 2; }
phase="$1"; mesh="$2"; ranks="$3"; final_time="$4"
[[ "$phase" == seed || "$phase" == resume ]] || { usage; exit 2; }
[[ "$mesh" == 0.2 || "$mesh" == 0.1 ]] || { usage; exit 2; }
[[ "$ranks" == 1 || "$ranks" == 2 || "$ranks" == 4 ]] || { usage; exit 2; }
[[ "$phase" == seed && $# == 4 || "$phase" == resume && $# == 5 ]] || {
    usage; exit 2;
}
checkpoint_dir=""
checkpoint_name=""
if [[ "$phase" == resume ]]; then
    source_checkpoint="$(realpath -e -- "$5")"
    [[ -f "$source_checkpoint" ]] || { echo "Checkpoint is not a file" >&2; exit 2; }
    checkpoint_dir="$(dirname "$source_checkpoint")"
    checkpoint_name="$(basename "$source_checkpoint")"
fi

python3 -m py_compile \
    sandbox/viscous_swe_2d_mpi.py \
    sandbox/probe_v130_issue28_one_cn_spinup.py \
    sandbox/probe_v130_rotating_controlled_ic.py \
    scrimp/utils/dae_checkpoint.py
python3 tests/test_v130_rotating_controlled_ic.py
bash -n "$0"

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
    echo "Refusing original-physics pilot on unqualified image revision=$revision" >&2
    exit 2
}
files=(
    scrimp/dphs.py scrimp/port.py
    scrimp/utils/mesh.py scrimp/utils/dae_checkpoint.py
    scrimp/utils/dae_initialization.py
    scrimp/utils/nonlinear_affine_dae.py
    sandbox/viscous_swe_2d_mpi.py
    sandbox/probe_v130_issue28_structure.py
    sandbox/probe_v130_issue28_weak_wall.py
    sandbox/probe_v130_issue28_one_cn_spinup.py
    sandbox/probe_v130_rotating_controlled_ic.py
)
for f in "${files[@]}"; do
    [[ -s "$f" ]] || { echo "Required source missing: $f" >&2; exit 2; }
done
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-rotating-controlled-ic-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
mkdir -p "$out/work"
source_sha="$(git rev-parse HEAD)"
printf '{"schema":"scrimp-v130-rotating-controlled-ic-runner-v1","source":"%s","image_revision":"%s","phase":"%s","mesh":%s,"MPI":%s,"absolute_final_time":%s,"normal_feedback":"original_decay","scientific_status":"NATIVE_NOT_YET_RUN"}\n' \
    "$source_sha" "$revision" "$phase" "$mesh" "$ranks" "$final_time" \
    > "$out/manifest.json"
echo "Physical controlled rotating SWE preparation: gamma=0.25, mesh=$mesh, MPI=$ranks"
echo "Original GetFEM/PETSc weak forms, C1 wall startup, original normal feedback envelope"
echo "Input phase=$phase, final absolute physical time=$final_time, dt=0.002"
echo "Evidence on host: $out"

docker_args=(
    run --rm
    --user "$(id -u):$(id -g)"
    --workdir /scrimp-work
    --env HOME=/tmp --env MPLBACKEND=Agg
    --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
    --env SCRIMP_GETFEM_MATRIX_MODE=distributed
    --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo
    --env SCRIMP_MATRIX_MEMORY_MODE=compact
    --env SCRIMP_MUMPS_RESET_POLICY=on_change
    --env SCRIMP_GETFEM_RHS_MODE=auto
    --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem
    --volume "$out/work:/scrimp-work:rw"
    --volume "$out:/artifacts:rw"
)
if [[ "$phase" == resume ]]; then
    docker_args+=(--volume "$checkpoint_dir:/source-checkpoint:ro")
fi
for f in "${files[@]}"; do
    docker_args+=(--volume "$root/$f:/workspace/scrimp/$f:ro")
done
native_args=(
    --phase "$phase" --mesh "$mesh" --gamma 0.25
    --dt 0.002 --final-time "$final_time"
    --checkpoint-out /artifacts/accepted_state.npz
    --out /artifacts
)
if [[ "$phase" == resume ]]; then
    native_args+=(--checkpoint-in "/source-checkpoint/$checkpoint_name")
fi
if ! docker_exec "${docker_args[@]}" "$image" \
    timeout -k 5s "${SCRIMP_V130_ROTATING_IC_TIMEOUT:-600s}" \
    mpiexec --oversubscribe --bind-to none -n "$ranks" \
    python -u /workspace/scrimp/sandbox/probe_v130_rotating_controlled_ic.py \
    "${native_args[@]}" > "$out/native.log" 2>&1; then
    echo "FAIL native rotating initial-state pilot; preserved: $out" >&2
    tail -n 130 "$out/native.log" >&2
    exit 1
fi
if ! grep -a 'ROTATING CONTROLLED IC ACCEPTED' "$out/native.log"; then
    echo "FAIL missing native accepted-state marker: $out" >&2
    tail -n 100 "$out/native.log" >&2
    exit 1
fi
python3 - "$out" "$mesh" "$ranks" "$final_time" <<'PY'
import json, math, sys
from pathlib import Path
out=Path(sys.argv[1])
r=json.loads((out/"initial_state_report.json").read_text(encoding="utf-8"))
assert r["status"]=="ACCEPTED_COMPATIBLE_ROTATING_DAE_INITIAL_STATE_NOT_STEADY"
assert r["gamma"]==0.25 and r["mesh"]==float(sys.argv[2])
assert r["MPI"]==int(sys.argv[3])
assert math.isclose(r["physical_time_final"],float(sys.argv[4]),abs_tol=1e-12)
assert r["original_TSCN_max_l2"]<1e-7
assert r["original_Fa_max_l2"]<1e-8
assert r["original_power_defect_over_exchange"]<1e-3
assert r["original_volume_balance_defect"]<1e-8
assert r["minimum_positive_height_FE_coefficient"]>0
assert r["stationary_rotating_state_certified"] is False
assert (out/"accepted_state.npz").stat().st_size>0
print("HOST EVIDENCE PASS: accepted full DAE state; equilibrium NOT claimed")
PY
echo "CHECKPOINT (host): $out/accepted_state.npz"
echo "REPORT (host): $out/initial_state_report.json"
echo "Future segments: $0 resume $mesh <MPI> <next-absolute-time> '$out/accepted_state.npz'"
