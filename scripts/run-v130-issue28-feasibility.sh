#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Issue #28.1: fixed h,p native algebraic feasibility; NO time steps.
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
    echo "FAIL: unqualified Docker image revision $rev" >&2
    exit 2
fi
for path in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py sandbox/viscous_swe_2d_mpi.py sandbox/probe_v130_issue28_structure.py sandbox/probe_v130_issue28_feasibility.py; do
    [[ -s "$path" ]] || { echo "FAIL missing $path" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-feasibility-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"source":"%s","image":"%s","image_revision":"%s","no_initialization":true,"no_physical_TS":true,"read_only_diagnostic":true}\n' \
    "$sha" "$image" "$rev" > "$out/manifest.json"
echo "Issue #28.1: rotating initial p0 vs stationary control, native MPI 1/2."
echo "Source=$sha; image=$image; immutable revision=$rev"
echo "Host evidence: $out"
for scenario in rest original; do
    if [[ "$scenario" == rest ]]; then
        omega="0"
    else
        omega="0.2"
    fi
    for ranks in 1 2; do
        mkdir -p "$out/work_${scenario}_${ranks}" "$out/result_${scenario}_${ranks}"
        echo "AUDIT $scenario rotation=$omega MPI=$ranks (no time stepping)"
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
            --volume "$root/sandbox/viscous_swe_2d_mpi.py:/workspace/scrimp/sandbox/viscous_swe_2d_mpi.py:ro" \
            --volume "$root/sandbox/probe_v130_issue28_structure.py:/workspace/scrimp/sandbox/probe_v130_issue28_structure.py:ro" \
            --volume "$root/sandbox/probe_v130_issue28_feasibility.py:/workspace/scrimp/sandbox/probe_v130_issue28_feasibility.py:ro" \
            "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_TIMEOUT || echo 240s)" \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /workspace/scrimp/sandbox/probe_v130_issue28_feasibility.py \
                --rotation "$omega" --mesh .2 --out /artifacts \
            > "$out/${scenario}_${ranks}.log" 2>&1; then
            echo "FAIL: native $scenario MPI$ranks; evidence retained $out" >&2
            tail -n 140 "$out/${scenario}_${ranks}.log" >&2
            exit 1
        fi
        grep -a 'ISSUE28 FEASIBILITY:' "$out/${scenario}_${ranks}.log" || true
    done
done
python3 - "$out" <<'PY'
import json, math, pathlib, sys
p=pathlib.Path(sys.argv[1])
for scenario in ("rest","original"):
    a=json.loads((p/f"result_{scenario}_1/feasibility.json").read_text())
    b=json.loads((p/f"result_{scenario}_2/feasibility.json").read_text())
    assert a["status"]==b["status"], (scenario,a["status"],b["status"])
    assert a["status"]!="INCONCLUSIVE_LSQR", ("uncertified least squares",scenario)
    for k in ("dofs","algebraic_rows","adjustable_algebraic_dofs",
              "excluded_null_traction_Y_dofs","reduced_matrix_shape"):
        assert a[k]==b[k], (scenario,k)
    for k in ("initial_Fa_l2","least_squares_final_Fa_l2",
              "affine_model_error_l2","normal_equation_AT_Fa_l2"):
        assert math.isclose(a[k],b[k],rel_tol=2e-5,abs_tol=1e-8), (scenario,k,a[k],b[k])
    for x in (a,b):
        assert x["differential_h0_p0_bitwise_preserved"]
        assert x["entire_initial_state_rollback_bitwise"]
        assert not x["initialization_TS_executed"]
        assert not x["main_TS_executed"]
        assert x["linear_dependence_verified"]
    if scenario=="rest":
        assert a["status"]=="FEASIBLE_AT_TOLERANCE", "rest case not calibrated"
    print(f"PASS {scenario}: {a['status']} / "
          f"initial |Fa|={a['initial_Fa_l2']:.9g}, "
          f"best |Fa|={a['least_squares_final_Fa_l2']:.9g} "
          f"(MPI1/2 parity)")
print("PASS: feasibility DIAGNOSTICS, NOT an index or nonlinear initializer qualification")
PY
echo "Retained evidence: $out"
