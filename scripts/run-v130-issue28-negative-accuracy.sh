#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# #28: numerical accuracy rejection AFTER computed algebraic correction;
# original rectangular shallow water, true MPI1/2 production DPHS.
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
expected="1d6d0c90115213298047b28237a7f38087563a35"
revision="$(docker_exec image inspect --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
if [[ "$revision" != "$expected" ]]; then
    echo "Refusing unqualified image revision $revision" >&2
    exit 2
fi
for f in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py \
    scrimp/utils/nonlinear_affine_dae.py \
    examples/shallow_water.py sandbox/check_v130_shallow_water_mass_jacobian.py \
    sandbox/probe_v130_issue28_negative_accuracy.py; do
    [[ -s "$f" ]] || { echo "Missing source $f" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-negative-accuracy-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"issue28-negative-accuracy-v1","source":"%s","image_revision":"%s","numerical_accuracy_only":"1e-16","physical_time_steps":0}\n' "$sha" "$revision" > "$out/manifest.json"
echo "#28 negative control: reject post-solve trial on strict accuracy 1e-16"
echo "Original shallow_water.py h0,p0 and boundary laws unchanged."
echo "Check rollback GetFEM + PETSc GlobalVec + explicit RHS cache, MPI 1/2."
echo "Host evidence: $out"
for ranks in 1 2; do
    mkdir -p "$out/work_${ranks}" "$out/result_${ranks}"
    echo "RUN expected numerical rejection after solved algebraic trial MPI=$ranks"
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
        --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
        --volume "$root/sandbox/check_v130_shallow_water_mass_jacobian.py:/workspace/scrimp/sandbox/check_v130_shallow_water_mass_jacobian.py:ro" \
        --volume "$root/sandbox/probe_v130_issue28_negative_accuracy.py:/workspace/scrimp/sandbox/probe_v130_issue28_negative_accuracy.py:ro" \
        "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_TIMEOUT || echo 240s)" \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_issue28_negative_accuracy.py \
            --out /artifacts \
        > "$out/mpi_${ranks}.log" 2>&1; then
        echo "FAIL: negative accuracy rollback MPI=$ranks; evidence $out" >&2
        tail -n 150 "$out/mpi_${ranks}.log" >&2
        exit 1
    fi
    grep -a 'ISSUE28 NEGATIVE ACCURACY' "$out/mpi_${ranks}.log"
done
python3 - "$out" <<'PY'
import json
import math
import pathlib
import sys
p=pathlib.Path(sys.argv[1])
reports=[
    json.loads((p/f"result_{n}/negative_accuracy.json").read_text())
    for n in (1,2)
]
for n, r in zip((1,2),reports):
    assert r["schema"]=="scrimp-v130-issue28-negative-accuracy-v1"
    assert r["status"]=="PASS_EXPECTED_NUMERICAL_REJECTION"
    assert r["mpi_ranks"]==n
    assert r["computed_candidate_Fa_l2"]>r["numerical_accuracy_target"]
    for key in (
        "getfem_bitwise_rolled_back",
        "petsc_GlobalVec_bitwise_rolled_back",
        "explicit_cache_ready_restored",
        "explicit_cache_time_restored",
        "explicit_cache_previous_contents_preserved",
        "algebraic_correction_computed",
    ):
        assert r[key] is True,key
    assert r["physical_initial_h0_p0_modified"] is False
    assert r["physically_incompatible_initial_data_demonstrated"] is False
    assert r["production_initialization_report"]["status"]=="REJECTED"
    assert not r["pseudo_TS_executed"] and not r["main_TS_executed"]
assert math.isclose(
    reports[0]["computed_candidate_Fa_l2"],
    reports[1]["computed_candidate_Fa_l2"],
    rel_tol=1e-3,abs_tol=1e-11,
), "MPI 1/2 candidate-residual mismatch"
print("PASS native MPI1/2 expected numerical-accuracy rejection after MUMPS solve")
print("Model state, PETSc GlobalVec and explicit-cache rollback verified.")
print("This does NOT demonstrate physically incompatible differential initial data.")
PY
echo "Evidence retained: $out"
