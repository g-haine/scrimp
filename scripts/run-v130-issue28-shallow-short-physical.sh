#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Issue #28: SHORT BDF2 physical time gate after opt-in nonlinear_affine
# production initial projection, unchanged ORIGINAL rectangular shallow SWE.
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
for p in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py \
    scrimp/utils/nonlinear_affine_dae.py \
    examples/shallow_water.py \
    sandbox/probe_v130_issue28_shallow_short_physical.py; do
    [[ -s "$p" ]] || { echo "Missing source $p" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-shallow-short-physical-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"issue28-shallow-short-physical-v1","source":"%s","image_revision":"%s","original_mesh":0.1,"time_integrator":"BDF2","dt":0.0001,"t_final":0.0003,"fem_and_initial_data_unchanged":true}\n' "$sha" "$revision" > "$out/manifest.json"
echo "#28: ORIGINAL shallow_water.py short PHYSICAL BDF2 after production nonlinear_affine"
echo "Source=$sha image=$image @ immutable OCI=$revision"
echo "Original rectangular FEM/h0/p0/zero wall/explicit gyro, 3 short steps."
echo "No claim of exact BDF2 energy conservation or rotating-Y0 qualification."
echo "Host evidence: $out"
for ranks in 1 2; do
    mkdir -p "$out/work_${ranks}" "$out/result_${ranks}"
    echo "RUN original shallow-water short physical MPI=$ranks"
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
        --volume "$root/sandbox/probe_v130_issue28_shallow_short_physical.py:/workspace/scrimp/sandbox/probe_v130_issue28_shallow_short_physical.py:ro" \
        "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_SHORT_TS_TIMEOUT || echo 300s)" \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_issue28_shallow_short_physical.py \
            --out /artifacts \
        > "$out/mpi_${ranks}.log" 2>&1; then
        echo "FAIL original shallow short PHYSICAL MPI=$ranks" >&2
        tail -n 145 "$out/mpi_${ranks}.log" >&2
        echo "All evidence retained: $out" >&2
        exit 1
    fi
    grep -a 'ISSUE28 SHALLOW SHORT PHYSICAL' "$out/mpi_${ranks}.log"
done
python3 - "$out" <<'PY'
import json
import math
import pathlib
import sys
root=pathlib.Path(sys.argv[1])
runs=[
  json.loads((root/f"result_{n}/shallow_short_physical.json").read_text())
  for n in (1,2)
]
for rank,r in zip((1,2),runs):
    assert r["schema"]=="scrimp-v130-issue28-original-shallow-short-physical-v1"
    assert r["status"]=="PASS_SHORT_PHYSICAL_TIME_GATE"
    assert r["mpi_ranks"]==rank
    assert r["original_prescribed_h0_p0_bitwise"]
    assert r["original_boundary_controls_unchanged_zero"]
    assert r["original_explicit_gyro_unchanged"]
    assert r["production_initializer_report"]["status"]=="PASS"
    assert r["production_initializer_report"]["differential_initial_exact"]
    assert 2 <= r["accepted_steps"] <= 6 and r["ts_converged_reason"]>0
    assert r["rejected_steps"]==r["nonlinear_solver_failures"]==0
    assert math.isclose(r["t_final_actual"],0.0003,abs_tol=1e-12)
    assert r["relative_mass_drift"]<=1e-5
    assert r["relative_hamiltonian_drift"]<=5e-3
    assert r["maximum_abs_boundary_power"]<=1e-7
    assert r["Hamiltonian_GetFEM_and_registered_agree"]
    assert not r["BDF2_exact_power_conservation_claimed"]
    assert not r["original_rotating_viscous_Y0_qualified"]
for key in ("accepted_steps","saved_times"):
    assert runs[0][key]==runs[1][key],f"MPI rank parity {key}"
for key in ("mass","hamiltonian","total_boundary_power"):
    a,b=runs[0][key],runs[1][key]
    assert len(a)==len(b)>1
    for x,y in zip(a,b):
        assert math.isclose(x,y,rel_tol=5e-6,abs_tol=1e-9), (key,x,y)
print("PASS: short BDF2 original shallow SWE physical TS after production nonlinear_affine, MPI 1/2.")
print("All Hamiltonian/power defects recorded; BDF2 exact power conservation NOT claimed.")
print("Original rotating viscous Y0 remains an open higher-index initializer problem.")
PY
echo "COMPLETE: retained native short-time evidence $out"
