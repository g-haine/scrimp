#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Issue #28.0: no initial or physical PETSc TS step is executed.
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
[[ "$revision" == "$expected" ]] || {
    echo "FAIL: wrong Docker image revision: $revision" >&2
    exit 2
}
for p in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py sandbox/viscous_swe_2d_mpi.py sandbox/probe_v130_issue28_structure.py; do
    [[ -s "$p" ]] || { echo "Missing source $p" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-structure-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"source":"%s","image":"%s","image_revision":"%s","main_TS_executed":false,"initialization_TS_executed":false}\n' "$sha" "$image" "$revision" > "$out/manifest.json"
echo "SCRIMP #28.0 native original SWE structure: $sha, immutable $image at $revision"
echo "Host evidence: $out"
for ranks in 1 2; do
    mkdir -p "$out/work_$ranks" "$out/result_$ranks"
    echo "AUDIT MPI=$ranks: original h0,p0, mass rows, Jacobian Y and IFunction FD"
    if ! docker_exec run --rm \
        --user "$(id -u):$(id -g)" --workdir /scrimp-work \
        --env HOME=/tmp --env MPLBACKEND=Agg \
        --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
        --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
        --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
        --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
        --env SCRIMP_MATRIX_MEMORY_MODE=compact \
        --env SCRIMP_TRAJECTORY_BACKEND=memory \
        --env SCRIMP_MUMPS_RESET_POLICY=on_change \
        --env SCRIMP_GETFEM_RHS_MODE=auto \
        --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
        --volume "$out/work_$ranks:/scrimp-work:rw" \
        --volume "$out/result_$ranks:/artifacts:rw" \
        --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
        --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro" \
        --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
        --volume "$root/sandbox/viscous_swe_2d_mpi.py:/workspace/scrimp/sandbox/viscous_swe_2d_mpi.py:ro" \
        --volume "$root/sandbox/probe_v130_issue28_structure.py:/workspace/scrimp/sandbox/probe_v130_issue28_structure.py:ro" \
        "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_TIMEOUT || echo 240s)" \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_issue28_structure.py \
            --mesh 0.2 --out /artifacts \
        > "$out/mpi_$ranks.log" 2>&1; then
        echo "FAIL: native structure MPI=$ranks; inspect $out/mpi_$ranks.log" >&2
        tail -n 120 "$out/mpi_$ranks.log" >&2
        exit 1
    fi
    grep -a "ISSUE28 STRUCTURE PASS" "$out/mpi_$ranks.log" || true
done
python3 - "$out" <<'PY'
import json,math,pathlib,sys
root=pathlib.Path(sys.argv[1])
a=json.loads((root/"result_1/structure.json").read_text())
b=json.loads((root/"result_2/structure.json").read_text())
for k in ("dofs","dynamic_dofs","algebraic_dofs","traction_Y_dofs","intervals"):
    assert a[k]==b[k], "MPI structural mismatch "+k
for k in ("max_abs_algebraic_mass_rows","max_abs_Jaa_Y_columns",
          "max_abs_Jda_Y_columns","Y_direction_Fa_FD_l2","Y_direction_Fd_FD_l2"):
    assert math.isclose(a[k],b[k],rel_tol=1e-6,abs_tol=1e-9), "MPI numerical mismatch "+k
assert a["differential_initial_bitwise_preserved"]
assert b["differential_initial_bitwise_preserved"]
print("PASS: original SWE structural MPI 1/2 parity; no time step/initialization performed")
PY
echo "COMPLETE: native structural evidence $out"
