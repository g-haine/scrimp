#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Issue #28 native differentiated-block diagnostic, no physical Y0 solve.
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
    echo "Refuse unqualified immutable Docker revision $revision" >&2
    exit 2
fi
for f in scrimp/dphs.py scrimp/port.py scrimp/utils/mesh.py \
    scrimp/utils/dae_initialization.py \
    scrimp/utils/nonlinear_affine_dae.py \
    sandbox/viscous_swe_2d_mpi.py sandbox/probe_v130_issue28_structure.py \
    sandbox/probe_v130_issue28_differentiated_rank.py; do
    [[ -s "$f" ]] || { echo "Missing original source: $f" >&2; exit 2; }
done

mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-differentiated-rank-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"v130-issue28-differentiated-rank-v1","source":"%s","image_revision":"%s","original_geometry_order":2,"original_mesh_h":0.2,"physical_initial_Y0_calculated":false,"index_classification_proven":false}\n' "$sha" "$revision" > "$out/manifest.json"

echo "#28 native FIRST-differentiated constraints, original viscous rotating disk"
echo "source=$sha image=$image @ immutable revision=$revision"
echo "Geometric P2 disk, original physical q0, corrected Rot in GWFL."
echo "Synthetic linear operator tests ONLY; no Y0, no main/pseudo TS."
echo "Evidence: $out"
for case in rest rotating; do
    for ranks in 1 2; do
        mkdir -p "$out/work_${case}_${ranks}" "$out/result_${case}_${ranks}"
        echo "RUN differentiated block audit case=$case MPI=$ranks"
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
            --volume "$out/work_${case}_${ranks}:/scrimp-work:rw" \
            --volume "$out/result_${case}_${ranks}:/artifacts:rw" \
            --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
            --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro" \
            --volume "$root/scrimp/utils/mesh.py:/workspace/scrimp/scrimp/utils/mesh.py:ro" \
            --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
            --volume "$root/scrimp/utils/nonlinear_affine_dae.py:/workspace/scrimp/scrimp/utils/nonlinear_affine_dae.py:ro" \
            --volume "$root/sandbox/viscous_swe_2d_mpi.py:/workspace/scrimp/sandbox/viscous_swe_2d_mpi.py:ro" \
            --volume "$root/sandbox/probe_v130_issue28_structure.py:/workspace/scrimp/sandbox/probe_v130_issue28_structure.py:ro" \
            --volume "$root/sandbox/probe_v130_issue28_differentiated_rank.py:/workspace/scrimp/sandbox/probe_v130_issue28_differentiated_rank.py:ro" \
            "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_RANK_TIMEOUT || echo 300s)" \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /workspace/scrimp/sandbox/probe_v130_issue28_differentiated_rank.py \
                --case "$case" --out /artifacts \
            > "$out/${case}_mpi${ranks}.log" 2>&1; then
            echo "FAIL: differentiated native probe case=$case MPI=$ranks" >&2
            echo "Host-owned evidence retained $out" >&2
            tail -n 150 "$out/${case}_mpi${ranks}.log" >&2
            exit 1
        fi
        grep -a 'ISSUE28 DIFFERENTIATED RANK' "$out/${case}_mpi${ranks}.log"
    done
done

python3 - "$out" <<'PY'
import json
import math
from pathlib import Path
import sys
base = Path(sys.argv[1])
for label in ("rest", "rotating"):
    result = []
    for rank in (1, 2):
        f = base/f"result_{label}_{rank}"/"differentiated_rank.json"
        p = json.loads(f.read_text(encoding="utf-8"))
        assert p["schema"] == "scrimp-v130-issue28-differentiated-rank-v1"
        assert p["initial_case"] == label and p["mpi_ranks"] == rank
        assert p["dofs"] == 3300 and p["dynamic_dofs"] == 1452
        assert p["algebraic_dofs"] == 1848
        assert p["traction_Y_dofs"] == 132
        assert p["regular_algebraic_dofs"] == 1716
        assert p["differentiated_bordered_shape"] == [3300, 3300]
        assert p["J_aa_Y_column_max_abs"] <= 1e-12
        assert p["J_dynamic_Y_column_max_abs"] > 1e-12
        assert p["original_q0_bitwise_preserved"]
        assert not p["physical_traction_Y0_calculated"]
        assert not p["rotating_algebraic_compatibility_established"]
        assert not p["differentiation_index_proven"]
        assert not p["real_initialization_TS_executed"]
        assert not p["real_physical_TS_executed"]
        assert p["status"] in (
            "NUMERIC_REGULAR_AT_LINEARIZATION",
            "NOT_NUMERICALLY_CERTIFIED",
        )
        if label == "rest":
            assert p["native_original_Fa_at_linearization_l2"] < 1e-11
            assert p["native_original_Fd_zdot_zero_at_linearization_l2"] < 1e-11
        if p["status"] == "NUMERIC_REGULAR_AT_LINEARIZATION":
            assert len(p["synthetic_solve_tests"]) == 2
            assert all(t["relative_known_synthetic_solution_error"]<=1e-8
                       and t["relative_synthetic_linear_residual"]<=1e-8
                       and t["ksp_converged_reason"]>0
                       for t in p["synthetic_solve_tests"])
        result.append(p)
    a,b = result
    assert a["status"] == b["status"], "MPI disagreement on bordered numerical status"
    assert math.isclose(
        a["native_original_Fa_at_linearization_l2"],
        b["native_original_Fa_at_linearization_l2"],
        rel_tol=1e-5, abs_tol=1e-10
    ), "MPI initial algebraic residual differs"
    print(f"{label} P2: native matrix audit MPI1/2 status={a['status']}, "
          f"Fa={a['native_original_Fa_at_linearization_l2']:.9g}, "
          f"Y={a['traction_Y_dofs']} columns")
print("AUDIT COMPLETE: matrix rank pilot, NOT proof of DAE index or rotating Y0.")
PY
echo "COMPLETE: retained raw native evidence $out"
