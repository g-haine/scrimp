#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# #28 native expected-rejection high-index rotating Y, original P2 geometry.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
docker_exec() {
    if [[ "$(printenv SCRIMP_DOCKER_NO_SUDO || true)" == 1 ]]; then
        docker "$@"
    else
        sudo docker "$@"
    fi
}
image="$(printenv SCRIMP_V130_DAE_IMAGE || echo scrimp-mpi:v1.3.0-rc3)"
revision="$(docker_exec image inspect --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
expected="1d6d0c90115213298047b28237a7f38087563a35"
if [[ "$revision" != "$expected" ]]; then
    echo "Refuse unqualified Docker image revision $revision" >&2
    exit 2
fi
for f in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py \
    scrimp/utils/nonlinear_affine_dae.py scrimp/utils/mesh.py \
    sandbox/viscous_swe_2d_mpi.py \
    sandbox/probe_v130_issue28_negative_high_index.py; do
    [[ -s "$f" ]] || { echo "Missing source $f" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-negative-high-index-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"issue28-native-negative-high-index-v1","source":"%s","image_revision":"%s","expected_rejection":true,"physical_time_steps":0}\n' "$sha" "$revision" > "$out/manifest.json"
echo "#28: production nonlinear_affine MUST REJECT original rotating viscous SWE"
echo "source=$sha image=$image at immutable revision=$revision"
echo "Original h0,p0, omega=.2, true P2 geometric disk, corrected GWFL wall"
echo "No pseudo time stepping or accepted physical TS. Host evidence: $out"
for ranks in 1 2; do
    mkdir -p "$out/work_${ranks}" "$out/result_${ranks}"
    echo "RUN expected rejection rotating P2 MPI=$ranks"
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
        --volume "$out/work_${ranks}:/scrimp-work:rw" \
        --volume "$out/result_${ranks}:/artifacts:rw" \
        --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
        --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro" \
        --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
        --volume "$root/scrimp/utils/nonlinear_affine_dae.py:/workspace/scrimp/scrimp/utils/nonlinear_affine_dae.py:ro" \
        --volume "$root/scrimp/utils/mesh.py:/workspace/scrimp/scrimp/utils/mesh.py:ro" \
        --volume "$root/sandbox/viscous_swe_2d_mpi.py:/workspace/scrimp/sandbox/viscous_swe_2d_mpi.py:ro" \
        --volume "$root/sandbox/probe_v130_issue28_negative_high_index.py:/workspace/scrimp/sandbox/probe_v130_issue28_negative_high_index.py:ro" \
        "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_TIMEOUT || echo 240s)" \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_issue28_negative_high_index.py \
            --out /artifacts \
        > "$out/mpi_${ranks}.log" 2>&1; then
        echo "FAIL: expected high-index rejection MPI=$ranks; evidence $out" >&2
        tail -n 140 "$out/mpi_${ranks}.log" >&2
        exit 1
    fi
    grep -a 'ISSUE28 NEGATIVE HIGH INDEX' "$out/mpi_${ranks}.log"
done
python3 - "$out" <<'PY'
import json
import pathlib
import sys

base=pathlib.Path(sys.argv[1])
reports=[
    json.loads((base/f"result_{n}/negative_high_index.json").read_text())
    for n in (1,2)
]
for n, d in zip((1,2), reports):
    assert d["schema"]=="scrimp-v130-issue28-negative-high-index-v1"
    assert d["status"]=="PASS_EXPECTED_REJECTION"
    assert d["mpi_ranks"]==n
    assert d["geometry_order"]==2 and d["omega"]==0.2
    assert d["detected_algebraic_null_columns"]>=d["algebraic_Y_dofs"]>0
    assert d["full_original_GetFEM_model_rolled_back"]
    assert d["original_fixed_h0_p0_preserved"]
    assert not d["physical_TS_executed"]
    assert not d["pseudo_TS_executed"]
    assert not d["rotating_high_index_Y0_qualified"]
    assert d["production_initialization_report"]["status"]=="REJECTED"
assert reports[0]["algebraic_Y_dofs"]==reports[1]["algebraic_Y_dofs"]
assert reports[0]["detected_algebraic_null_columns"]==reports[1]["detected_algebraic_null_columns"]
print("PASS: original rotating high-index Y structurally rejected in MPI 1/2.")
print("This is an expected NEGATIVE control, not a rotating DAE initialization.")
PY
echo "Retained evidence: $out"
