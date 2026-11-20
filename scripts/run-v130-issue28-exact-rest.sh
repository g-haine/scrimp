#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Stage #28: exact native resting control (no LSQR, no TS).
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
    echo "Refusing unexpected immutable Docker source $revision" >&2
    exit 2
fi
for file in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py \
    scrimp/utils/mesh.py sandbox/viscous_swe_2d_mpi.py \
    sandbox/probe_v130_issue28_structure.py \
    sandbox/probe_v130_issue28_exact_rest.py; do
    [[ -s "$file" ]] || { echo "Missing source file $file" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-exact-rest-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"source":"%s","image_revision":"%s","real_getfem_IFunction":true,"LSQR_executed":false,"initialization_TS_executed":false,"main_TS_executed":false}\n' "$sha" "$revision" > "$out/manifest.json"
echo "#28 exact resting equilibrium vs LSQR false-negative, unchanged viscous SWE"
echo "Source=$sha image=$image @ $revision"
echo "Evidence=$out"
for order in 1 2; do
    for ranks in 1 2; do
        mkdir -p "$out/work_p${order}_${ranks}" "$out/result_p${order}_${ranks}"
        echo "RUN exact rest P$order MPI=$ranks"
        if ! docker_exec run --rm \
            --user "$(id -u):$(id -g)" --workdir /scrimp-work \
            --env HOME=/tmp --env MPLBACKEND=Agg \
            --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
            --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
            --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
            --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
            --env SCRIMP_MATRIX_MEMORY_MODE=compact \
            --env SCRIMP_GETFEM_RHS_MODE=auto \
            --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
            --volume "$out/work_p${order}_${ranks}:/scrimp-work:rw" \
            --volume "$out/result_p${order}_${ranks}:/artifacts:rw" \
            --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
            --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro" \
            --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
            --volume "$root/scrimp/utils/mesh.py:/workspace/scrimp/scrimp/utils/mesh.py:ro" \
            --volume "$root/sandbox/viscous_swe_2d_mpi.py:/workspace/scrimp/sandbox/viscous_swe_2d_mpi.py:ro" \
            --volume "$root/sandbox/probe_v130_issue28_structure.py:/workspace/scrimp/sandbox/probe_v130_issue28_structure.py:ro" \
            --volume "$root/sandbox/probe_v130_issue28_exact_rest.py:/workspace/scrimp/sandbox/probe_v130_issue28_exact_rest.py:ro" \
            "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_TIMEOUT || echo 240s)" \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /workspace/scrimp/sandbox/probe_v130_issue28_exact_rest.py \
                --geometry-order "$order" --out /artifacts \
            > "$out/p${order}_${ranks}.log" 2>&1; then
            echo "FAIL native exact-rest P$order MPI=$ranks; logs retained $out" >&2
            tail -n 110 "$out/p${order}_${ranks}.log" >&2
            exit 1
        fi
        grep -a 'ISSUE28 EXACT REST PASS' "$out/p${order}_${ranks}.log"
    done
done
python3 - "$out" <<'PY'
import json, math, pathlib, sys
base = pathlib.Path(sys.argv[1])
for order in (1, 2):
    a = json.loads((base/f"result_p{order}_1/exact_rest.json").read_text())
    b = json.loads((base/f"result_p{order}_2/exact_rest.json").read_text())
    for x in (a, b):
        assert x["status"] == "EXACT_REST_NATIVE_PASS"
        assert x["dynamic_h0_p0_bitwise_preserved"]
        assert x["real_IFunction_evaluated"]
        assert not x["least_squares_used"]
        assert not x["main_TS_executed"]
        assert not x["initializer_TS_executed"]
        assert x["F_a_reference_l2"] <= x["tolerance_absolute"]
        assert x["F_d_reference_l2"] <= x["tolerance_absolute"]
        assert x["F_a_wrong_Yn_sign_l2"] > 1e-4
    for name in ("F_a_reference_l2","F_d_reference_l2","F_a_wrong_Yn_sign_l2"):
        assert math.isclose(a[name], b[name], rel_tol=1e-7, abs_tol=1e-12), name
    print(f"P{order} MPI1/2 PASS: exact rest |Fa|={a['F_a_reference_l2']:.9g}, "
          f"|Fd|={a['F_d_reference_l2']:.9g}, "
          f"wrong-Yn-sign |Fa|={a['F_a_wrong_Yn_sign_l2']:.9g}")
print("PASS: resting weak-form equilibrium reference, NOT rotating DAE qualification")
PY
echo "COMPLETE: host-owned exact-rest evidence $out"
