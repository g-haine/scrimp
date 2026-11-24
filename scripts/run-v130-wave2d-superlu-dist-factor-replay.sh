#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE.
#
# REPRODUCE real PC failure on the archived first-CN PETSc matrix,
# then evaluate the sole optional numerical pivot replacement under
# precisely the same PETSc library, on MPI 1 and 2.
# NO wave/DAE reassembly, NO change to matrix, NO fake solver fallback.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
    docker_cmd=(docker)
fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
    --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
    echo "Wrong pinned GetFEM/PETSc image revision=$revision" >&2; exit 2;
}
source_zip="$root/sandbox/outputs/v130-wave2d-first-cn-jacobian-574d9c71-NKqIZ0BZ.zip"
[[ -r "$source_zip" ]] || {
    echo "Missing actual first-CN factor-matrix archive: $source_zip" >&2
    exit 2
}
python3 -m py_compile \
    sandbox/replay_v130_wave2d_superlu_dist_factor.py \
    tests/test_v130_wave2d_factor_replay.py
python3 tests/test_v130_wave2d_factor_replay.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-wave2d-superlu-dist-factor-replay-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
mkdir -p "$out/input"
python3 - "$source_zip" "$out/input/matrix.npz" <<'PY'
import hashlib, io, sys, zipfile
import numpy as np
z, target = sys.argv[1:]
with zipfile.ZipFile(z) as f:
    blob = f.read("mumps/first_main_CN_factor.npz")
with np.load(io.BytesIO(blob),allow_pickle=False) as f:
    data = {key:f[key].copy() for key in f.files}
digest=hashlib.sha256()
for key in sorted(data):
    digest.update(key.encode())
    digest.update(np.ascontiguousarray(data[key]).tobytes())
assert digest.hexdigest()=="24346856f358e4f81cf332ebc974e329dc084c60f5d6b98bfbcd788fdff37f54", (
    "Wrong native P/M/K matrix: no replay")
open(target,"wb").write(blob)
print("PINNED ORIGINAL 2230x2230 P/M/K CSR verified SHA256:",digest.hexdigest())
PY
printf '{"schema":"v130-wave2d-true-PETSc-SuperLU_DIST-LU-replay-v1","source":"%s","image":"%s","image_revision":"%s","real_input_archive":"%s","operator_SHA256":"24346856f358e4f81cf332ebc974e329dc084c60f5d6b98bfbcd788fdff37f54","candidate_only":true,"full_wave_acceptance":"3/5"}\n' \
    "$(git rev-parse HEAD)" "$image" "$revision" "$(basename "$source_zip")" \
    > "$out/manifest.json"
echo "Real GetFEM CN factor matrix, no source modification, separate LU replay"
echo "Default vs optional SuperLU_DIST tiny-pivot replacement, MPI=1,2"
echo "Compare PETSc-owned coefficients and true Ax-b against exact x*"
echo "EVIDENCE: $out"
for mode in default replace_tiny_pivot; do
    for ranks in 1 2; do
        name="${mode}_mpi${ranks}"
        echo "RUN $name"
        mkdir -p "$out/$name" "$out/work_$name"
        rc=0
        "${docker_cmd[@]}" run --rm \
            --user "$(id -u):$(id -g)" \
            --env HOME=/tmp --env OMP_NUM_THREADS=1 \
            --env OPENBLAS_NUM_THREADS=1 \
            --volume "$out/input:/input:ro" \
            --volume "$root/sandbox/replay_v130_wave2d_superlu_dist_factor.py:/workspace/scrimp/sandbox/replay_v130_wave2d_superlu_dist_factor.py:ro" \
            --volume "$out/$name:/artifacts:rw" \
            --volume "$out/work_$name:/scrimp-work:rw" \
            --workdir /scrimp-work \
            "$image" timeout -k 5s "${SCRIMP_V130_WAVE2D_REPLAY_TIMEOUT:-120s}" \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /workspace/scrimp/sandbox/replay_v130_wave2d_superlu_dist_factor.py \
                --matrix /input/matrix.npz \
                --mode "$mode" \
                --output /artifacts \
            > "$out/$name.host.log" 2>&1 || rc=$?
        printf '%s\n' "$rc" > "$out/$name.exitcode"
        if [[ "$rc" != 0 ]]; then
            echo "ORIGINAL_BACKEND_FACTOR_FAILED $name exit=$rc — log saved" >&2
            grep -aE 'KSP_reason|PC_failed|Traceback|petsc4py|error code|RuntimeError' \
                "$out/$name.host.log" | tail -n 14 >&2 || true
        fi
    done
done
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
cases={}
passed=0
for mode in ("default","replace_tiny_pivot"):
    for ranks in (1,2):
        label=f"{mode}_mpi{ranks}"
        folder=root/label
        rc=int((root/f"{label}.exitcode").read_text())
        summary={
            "exit_code":rc, "MPI":ranks,"mode":mode,
            "rank_json_present":all((folder/f"rank{i}.json").is_file()
                                    for i in range(ranks)),
        }
        if summary["rank_json_present"]:
            evidence=[json.loads((folder/f"rank{i}.json").read_text())
                      for i in range(ranks)]
            status=all(
                e.get("factorization_passed_all_ranks") and
                e.get("original_full_P_M_K_SHA256") ==
                    "24346856f358e4f81cf332ebc974e329dc084c60f5d6b98bfbcd788fdff37f54" and
                e.get("PETSc_owned_CSR_graph_and_values_identical_to_archive") and
                e.get("factor_solver")=="superlu_dist" and
                e.get("relative_backward_residual") is not None and
                e["relative_backward_residual"] < 1e-11 and
                e.get("relative_known_solution_error") is not None and
                e["relative_known_solution_error"] < 1e-7
                for e in evidence
            )
            summary["reported_rank_statuses"]=[e.get("status") for e in evidence]
            summary["KSP_reasons"]=[e.get("KSP_reason") for e in evidence]
            summary["PC_failed"]=[e.get("PC_failed") for e in evidence]
            summary["max_relative_backward_error"]=max(
                (e.get("relative_backward_residual") or 0)
                for e in evidence
            )
            summary["max_relative_solution_error"]=max(
                (e.get("relative_known_solution_error") or 0)
                for e in evidence
            )
            summary["native_replay_pass"]=bool(rc==0 and status)
        else:
            summary["native_replay_pass"]=False
        cases[label]=summary
        if summary["native_replay_pass"]:
            passed+=1
pivot_mpi_pass=all(
    cases[f"replace_tiny_pivot_mpi{ranks}"]["native_replay_pass"]
    for ranks in (1,2)
)
report={
    "status": (
        "PASS_OPTIONAL_PETSC_LU_REPLAY_MPI1_MPI2_NOT_WAVE_ACCEPTANCE"
        if pivot_mpi_pass else
        "FAIL_OR_INCOMPLETE_OPTIONAL_PIVOT_PETSC_LU_REPLAY"
    ),
    "cases":cases,
    "default_solver_still_not_qualified":True,
    "full_wave_superlu_dist_solver_qualified":False,
    "previous_three_genuine_wave_backend_passes_preserved":True,
    "wave_installation_gate":"3/5",
}
(root/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as output:
    for path in sorted(root.rglob("*")):
        if path.is_file():
            output.write(path,arcname=str(path.relative_to(root)))
print(json.dumps({**report,"archive":str(archive)},indent=2),flush=True)
if not pivot_mpi_pass:
    raise SystemExit(1)
PY
echo "PIVOT_REPLACEMENT_PETSC_FACTOR_REPLAY_CANDIDATE_ONLY — full CN tests pending"
