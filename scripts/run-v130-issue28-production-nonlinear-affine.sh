#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Native qualification of *actual production* opt-in nonlinear_affine
# DAE initial projection; original rectangular shallow-water, MPI 1/2.
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
    echo "Unqualified immutable Docker image revision $revision" >&2
    exit 2
fi
for f in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py \
    scrimp/utils/nonlinear_affine_dae.py \
    examples/shallow_water.py sandbox/check_v130_shallow_water_mass_jacobian.py \
    sandbox/probe_v130_issue28_production_nonlinear_affine.py; do
    [[ -s "$f" ]] || { echo "Missing source $f" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-production-na-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"issue28-production-nonlinear-affine-v1","source":"%s","immutable_revision":"%s","opt_in_only":true,"physical_steps":0}\n' "$sha" "$revision" > "$out/manifest.json"
echo "#28 actual production DPHS nonlinear_affine initial projection, MPI 1/2"
echo "Source=$sha image=$image @ $revision"
echo "Original rectangular shallow_water.py, mesh=.1; no TS/physical time stepping."
echo "Host evidence: $out"
for ranks in 1 2; do
    mkdir -p "$out/work_${ranks}" "$out/result_${ranks}"
    echo "RUN actual production init_step nonlinear_affine MPI=$ranks"
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
        --volume "$root/sandbox/probe_v130_issue28_production_nonlinear_affine.py:/workspace/scrimp/sandbox/probe_v130_issue28_production_nonlinear_affine.py:ro" \
        "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_TIMEOUT || echo 240s)" \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_issue28_production_nonlinear_affine.py \
            --out /artifacts \
        > "$out/mpi_${ranks}.log" 2>&1; then
        echo "FAIL: production nonlinear_affine MPI=$ranks, evidence retained $out" >&2
        tail -n 150 "$out/mpi_${ranks}.log" >&2
        exit 1
    fi
    grep -a 'ISSUE28 PRODUCTION NONLINEAR_AFFINE' "$out/mpi_${ranks}.log"
done
python3 - "$out" <<'PY'
import json, math, pathlib, sys
base=pathlib.Path(sys.argv[1])
a=json.loads((base/"result_1/production_nonlinear_affine.json").read_text())
b=json.loads((base/"result_2/production_nonlinear_affine.json").read_text())
for x in (a,b):
    assert x["status"]=="PASS_NATIVE_PRODUCTION_INITIALIZER"
    assert x["original_dynamic_h0_p0_bitwise"]
    assert x["published_algebraic_dofs_modified"]
    assert x["actual_published_IFunction_Fa_l2"]<=1e-11
    assert x["explicit_gyro_from_original_shallow_water"]
    assert x["explicit_gyro_cache_ready"]
    assert not x["pseudo_initial_TS_executed"]
    assert not x["main_physical_TS_executed"]
    assert not x["high_index_rotating_Y_initialization_qualified"]
    r=x["real_production_init_step_report"]
    assert r["status"]=="PASS" and r["differential_initial_exact"]
    assert r["native_final_algebraic_residual_l2"]<=r["accuracy_target"]
    assert r["algebraic_null_columns"]==0
for field in ("dofs",):
    assert a[field]==b[field],field
for field in ("actual_published_IFunction_Fa_l2",):
    assert math.isclose(a[field],b[field],rtol=1e-3,abs_tol=1e-11)
print("PASS native production nonlinear_affine initializer path, MPI 1/2.")
print("Not a general nonlinear Newton or rotating high-index DAE qualification.")
PY
echo "COMPLETE: native retained evidence $out"
