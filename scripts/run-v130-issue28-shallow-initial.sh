#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# #28 original shipped shallow_water.py fixed-h,p DAE diagnostic.
# No BDF pseudo initializer, no main TS. True GetFEM/PETSc MPI 1/2.
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
rev="$(docker_exec image inspect --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
if [[ "$rev" != "$expected" ]]; then
    echo "Unqualified immutable Docker image revision $rev" >&2
    exit 2
fi
for f in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py \
    examples/shallow_water.py sandbox/check_v130_shallow_water_mass_jacobian.py \
    sandbox/probe_v130_issue28_structure.py \
    sandbox/probe_v130_issue28_shallow_initial.py; do
    [[ -s "$f" ]] || { echo "Missing source $f" >&2; exit 2; }
done
mesh="$(printenv SCRIMP_ISSUE28_SHALLOW_MESH || echo 0.1)"
if [[ "$mesh" != 0.1 ]]; then
    echo "This qualified setup preserves the original example mesh=0.1" >&2
    exit 2
fi
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-shallow-initial-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"v130-issue28-shallow-initial-v1","source":"%s","image_revision":"%s","original_example":true,"original_mesh":0.1,"physical_time_steps":0,"production_initializer_qualified":false}\n' "$sha" "$rev" > "$out/manifest.json"
echo "Issue #28: native ORIGINAL examples/shallow_water.py nonlinear-DAE structural/one-step algebraic diagnostic"
echo "Source=$sha image=$image revision=$rev; exact published h0,p0 and mesh=0.1"
echo "Explicit gyroscopic brick G frozen at accepted initial state; no TS/no accepted result."
echo "Host evidence: $out"
for ranks in 1 2; do
    mkdir -p "$out/work_${ranks}" "$out/result_${ranks}"
    echo "RUN original shallow-water fixed dynamics MPI=$ranks"
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
        --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
        --volume "$root/sandbox/check_v130_shallow_water_mass_jacobian.py:/workspace/scrimp/sandbox/check_v130_shallow_water_mass_jacobian.py:ro" \
        --volume "$root/sandbox/probe_v130_issue28_structure.py:/workspace/scrimp/sandbox/probe_v130_issue28_structure.py:ro" \
        --volume "$root/sandbox/probe_v130_issue28_shallow_initial.py:/workspace/scrimp/sandbox/probe_v130_issue28_shallow_initial.py:ro" \
        "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_TIMEOUT || echo 240s)" \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_issue28_shallow_initial.py \
            --mesh .1 --out /artifacts \
        > "$out/mpi_${ranks}.log" 2>&1; then
        echo "FAIL: original shallow-water MPI=$ranks; evidence retained at $out" >&2
        tail -n 140 "$out/mpi_${ranks}.log" >&2
        exit 1
    fi
    grep -a 'ISSUE28 SHALLOW FIXED DYNAMICS PASS' "$out/mpi_${ranks}.log"
done
python3 - "$out" <<'PY'
import json, math, pathlib, sys
p=pathlib.Path(sys.argv[1])
a=json.loads((p/"result_1/shallow_fixed_dynamics.json").read_text())
b=json.loads((p/"result_2/shallow_fixed_dynamics.json").read_text())
for x in (a,b):
    assert x["status"]=="PASS_DIAGNOSTIC"
    assert x["nonlinear_mass"]
    assert x["explicit_gyroscopic_G_frozen_at_original_initial_state"]
    assert x["h0_p0_bitwise_preserved"]
    assert x["getfem_initial_state_rollback_bitwise"]
    assert x["max_abs_algebraic_mass_rows"]<=1e-12
    assert x["after_one_direct_algebraic_Newton_Fa_l2"]<=1e-9
    assert x["relative_algebraic_Jacobian_FD_error"]<=1e-7
    assert x["frozen_gyro_rhs_drift"]==0
    assert not x["main_TS_executed"] and not x["initializer_TS_executed"]
    assert not x["genuine_nonlinear_algebraic_Newton_qualified"]
    assert not x["production_nonlinear_initializer_qualified"]
for key in ("dofs","dynamic_dofs","algebraic_dofs","algebraic_names"):
    assert a[key]==b[key], (key,a[key],b[key])
for key in ("initial_Fa_l2","after_one_direct_algebraic_Newton_Fa_l2",
            "relative_algebraic_Jacobian_FD_error"):
    assert math.isclose(a[key],b[key],rel_tol=1e-4,abs_tol=2e-9), (key,a[key],b[key])
print("PASS: original nonlinear shallow-water fixed-state algebraic diagnostic, MPI 1/2 parity")
print("Production nonlinear DAE initialization and full physical TS remain UNQUALIFIED.")
PY
echo "COMPLETE: host-owned evidence $out"
