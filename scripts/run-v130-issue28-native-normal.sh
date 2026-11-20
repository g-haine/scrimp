#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# #28: measure GetFEM's actual Normal sign on original rotating p0.
# Independent of DAE LSQR and without any initial/physical TS step.
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
for file in scrimp/dphs.py scrimp/port.py scrimp/utils/dae_initialization.py \
        scrimp/utils/mesh.py sandbox/viscous_swe_2d_mpi.py \
        sandbox/probe_v130_issue28_normal_orientation.py; do
    [[ -s "$file" ]] || { echo "Missing source $file" >&2; exit 2; }
done
mkdir -p sandbox/outputs
sha="$(git rev-parse HEAD)"
out="$(mktemp -d "$root/sandbox/outputs/v130-issue28-native-normal-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"source":"%s","image_revision":"%s","original_omega":0.2,"main_TS_executed":false,"initializer_executed":false}\n' "$sha" "$revision" > "$out/manifest.json"
echo "#28: GetFEM actual boundary Normal/Tangent sign for ORIGINAL rotating p0"
echo "source=$sha image=$image at $revision"
echo "No LSQR, no initial/physical TS, no modification to PDE or p0."
echo "Evidence: $out"
for order in 1 2; do
    for ranks in 1 2; do
        mkdir -p "$out/work_p${order}_${ranks}" "$out/result_p${order}_${ranks}"
        echo "RUN original GetFEM boundary-normal P$order MPI=$ranks"
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
            --volume "$root/sandbox/probe_v130_issue28_normal_orientation.py:/workspace/scrimp/sandbox/probe_v130_issue28_normal_orientation.py:ro" \
            "$image" timeout -k 5s "$(printenv SCRIMP_ISSUE28_TIMEOUT || echo 240s)" \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /workspace/scrimp/sandbox/probe_v130_issue28_normal_orientation.py \
                --geometry-order "$order" --out /artifacts \
            > "$out/p${order}_${ranks}.log" 2>&1; then
            echo "FAIL native boundary-normal diagnostic P$order MPI=$ranks" >&2
            tail -n 110 "$out/p${order}_${ranks}.log" >&2
            echo "Saved all evidence: $out" >&2
            exit 1
        fi
        grep -a "ISSUE28 NATIVE NORMAL" "$out/p${order}_${ranks}.log" || true
    done
done
python3 - "$out" <<'PY'
import json
import math
import pathlib
import sys
base=pathlib.Path(sys.argv[1])
for order in (1,2):
    one=json.loads((base/f"result_p{order}_1/normal_orientation.json").read_text())
    two=json.loads((base/f"result_p{order}_2/normal_orientation.json").read_text())
    for d in (one,two):
        assert d["geometry_order"]==order
        assert d["initial_state_bitwise_preserved"]
        assert not d["nonlinear_initializer_executed"]
        assert not d["physical_time_step_executed"]
        assert len(d["rank_local_measurements"])==d["mpi_ranks"]
    if one["normal_orientation_result"]!=two["normal_orientation_result"]:
        raise SystemExit("GetFEM orientation mismatch between MPI 1 and MPI 2")
    print(f"P{order}: GetFEM Native Normal sense="
          f"{one['normal_orientation_result']} "
          f"aligned mismatch={one['relative_aligned_trace_error']:.9g}, "
          f"reversed mismatch={one['relative_reversed_trace_error']:.9g}, "
          f"normal leakage={one['relative_normal_leakage']:.9g}")
    if one["normal_orientation_result"]!="ALIGNED":
        print("WARNING: original tangent sign is not aligned with "
              "actual GetFEM Normal. Investigate boundary region/mesh tags.")
print("COMPLETE: native GetFEM sign diagnosis; this is NOT a DAE compatibility certificate")
PY
echo "Evidence preserved: $out"
