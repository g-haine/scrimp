#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Native #28 P2 geometry feasibility diagnostic: original wall law.
# No initialization, time step, physics modification, or silent pass.
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
    echo "Refusing unqualified Docker image: $revision" >&2
    exit 2
fi
for file in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py \
    scrimp/utils/mesh.py sandbox/viscous_swe_2d_mpi.py \
    sandbox/probe_v130_issue28_structure.py sandbox/probe_v130_issue28_feasibility.py; do
    [[ -s "$file" ]] || { echo "Missing source file $file" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-P2-geometry-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"issue28-P2-curved-geometry-v1","source":"%s","image_revision":"%s","geometry_order":2,"h_order":2,"p_order":2,"physical_weak_forms_modified":false,"time_step_executed":false}\n' "$sha" "$revision" > "$out/manifest.json"
echo "Issue #28 curved GEOMETRIC P2 (field FEM P2, unchanged), source=$sha"
echo "Host evidence: $out"
echo "Original physical U=omega*Rot*Normal, p0=tangent, no initial TS."
for scenario in rest original; do
    case "$scenario" in
        rest) omega=0 ;;
        original) omega=0.2 ;;
    esac
    for ranks in 1 2; do
        mkdir -p "$out/work_${scenario}_${ranks}" "$out/result_${scenario}_${ranks}"
        echo "RUN $scenario rotation=$omega geometry=P2 MPI=$ranks"
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
            --volume "$out/work_${scenario}_${ranks}:/scrimp-work:rw" \
            --volume "$out/result_${scenario}_${ranks}:/artifacts:rw" \
            --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
            --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro" \
            --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
            --volume "$root/scrimp/utils/mesh.py:/workspace/scrimp/scrimp/utils/mesh.py:ro" \
            --volume "$root/sandbox/viscous_swe_2d_mpi.py:/workspace/scrimp/sandbox/viscous_swe_2d_mpi.py:ro" \
            --volume "$root/sandbox/probe_v130_issue28_structure.py:/workspace/scrimp/sandbox/probe_v130_issue28_structure.py:ro" \
            --volume "$root/sandbox/probe_v130_issue28_feasibility.py:/workspace/scrimp/sandbox/probe_v130_issue28_feasibility.py:ro" \
            "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_TIMEOUT || echo 240s)" \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /workspace/scrimp/sandbox/probe_v130_issue28_feasibility.py \
                --rotation "$omega" --mesh .2 --geometry-order 2 --out /artifacts \
            > "$out/${scenario}_${ranks}.log" 2>&1; then
            echo "FAIL: native curved-geometry $scenario MPI=$ranks" >&2
            tail -n 120 "$out/${scenario}_${ranks}.log" >&2
            echo "Retained evidence: $out" >&2
            exit 1
        fi
        grep -a 'ISSUE28 FEASIBILITY:' "$out/${scenario}_${ranks}.log" || true
    done
done
python3 - "$out" <<'PY'
import json, math, pathlib, sys
p=pathlib.Path(sys.argv[1])
for scenario in ("rest","original"):
    one=json.loads((p/f"result_{scenario}_1/feasibility.json").read_text())
    two=json.loads((p/f"result_{scenario}_2/feasibility.json").read_text())
    for x in (one,two):
        assert x["geometry_order_requested"]==2
        assert x["getfem_import_has_nonlinear_geometric_transform"], (
            "P2 requested, but GetFEM did not preserve curved geometry",
            x["getfem_imported_geotransform_descriptors"])
        assert x["differential_h0_p0_bitwise_preserved"]
        assert x["entire_initial_state_rollback_bitwise"]
        assert x["linear_dependence_verified"]
        assert not x["main_TS_executed"]
        assert not x["initialization_TS_executed"]
    assert one["status"]==two["status"], ("MPI status mismatch",scenario)
    for k in ("dofs","algebraic_rows","adjustable_algebraic_dofs",
              "excluded_null_traction_Y_dofs","getfem_imported_geometry_nodes"):
        assert one[k]==two[k], ("MPI structure mismatch",scenario,k)
    for k in ("initial_Fa_l2","least_squares_final_Fa_l2",
              "normal_equation_AT_Fa_l2"):
        assert math.isclose(one[k],two[k],rel_tol=2e-5,abs_tol=2e-8), (
            "MPI numerical mismatch",scenario,k,one[k],two[k])
    print(f"{scenario}: status={one['status']}, "
          f"|Fa|={one['least_squares_final_Fa_l2']:.9g}, "
          f"GetFEM nodes={one['getfem_imported_geometry_nodes']}, "
          f"elements={one['getfem_imported_element_count']}")
    print(" GetFEM transformations:",one["getfem_imported_geotransform_descriptors"])
    if scenario=="original":
        reference=0.02656138  # prior P1 original, supplied native evidence
        ratio=one["least_squares_final_Fa_l2"]/reference
        print(f"Original P2/P1 residual ratio (prior P1 source 9082bbee): {ratio:.9g}")
        print("A smaller residual is evidence of improvement, not a "
              "certificate of exact initial compatibility.")
print("PASS diagnostic P2 import and MPI agreement; NOT a physical release qualification")
PY
echo "Evidence preserved: $out"
