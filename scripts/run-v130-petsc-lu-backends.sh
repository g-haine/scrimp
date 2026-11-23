#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# QUALIFICATION PREFLIGHT, NOT A PERFORMANCE BENCHMARK.
# Do not touch/build/retag/prune the existing qualified Docker image.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" != 1 ]]; then
    docker_cmd=(sudo docker)
fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
    --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
expected="1d6d0c90115213298047b28237a7f38087563a35"
[[ "$revision" == "$expected" ]] || {
    echo "Refusing to use unqualified PETSc/GetFEM image: $revision" >&2
    exit 2
}
python3 -m py_compile scripts/probe-v130-petsc-lu-backends.py
bash -n "$0"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-petsc-lu-preflight-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-petsc-lu-backends-campaign-v1","image":"%s","revision":"%s","source":"%s","status":"PENDING","matrix":"distributed complex nonsymmetric AIJ","not_an_examples_gate":true,"not_a_benchmark":true}\n' \
    "$image" "$revision" "$(git rev-parse HEAD)" > "$out/manifest.json"
echo "Original qualified PETSc image: $image ($revision)"
echo "Compare availability with ACTUAL distributed LU accuracy; no rebuilt Docker."
echo "MUMPS and SuperLU_DIST: MPI 1,2,4. Sequential SuperLU: MPI 1 only."
echo "EVIDENCE: $out"
for backend in mumps superlu superlu_dist; do
    ranks_list=(1 2 4)
    if [[ "$backend" == superlu ]]; then ranks_list=(1); fi
    for ranks in "${ranks_list[@]}"; do
        case_name="${backend}_mpi${ranks}"
        echo "RUN $case_name"
        if ! "${docker_cmd[@]}" run --rm \
            --user "$(id -u):$(id -g)" \
            --env HOME=/tmp --env OMP_NUM_THREADS=1 \
            --env OPENBLAS_NUM_THREADS=1 \
            --volume "$root/scripts/probe-v130-petsc-lu-backends.py:/probe.py:ro" \
            --volume "$out:/results:rw" \
            --workdir /tmp \
            "$image" timeout -k 5s "${SCRIMP_V130_LU_PREFLIGHT_TIMEOUT:-90s}" \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /probe.py --backend "$backend" --size 384 \
            --out "/results/$case_name.json" \
            > "$out/$case_name.host.log" 2>&1; then
            echo "FAIL actual PETSc LU backend $case_name (logs kept: $out)" >&2
            tail -n 100 "$out/$case_name.host.log" >&2
            exit 1
        fi
        grep -a '^PASS PETSc LU backend=' "$out/$case_name.host.log" || {
            echo "Missing native factorization PASS for $case_name" >&2
            exit 1
        }
    done
done

python3 - "$out" <<'PY'
import json
import sys
import zipfile
from pathlib import Path
directory=Path(sys.argv[1])
cases={(backend,rank) for backend in ("mumps","superlu_dist")
       for rank in (1,2,4)}
cases.add(("superlu",1))
results=[]
for backend,ranks in sorted(cases):
    path=directory/f"{backend}_mpi{ranks}.json"
    row=json.loads(path.read_text(encoding="utf-8"))
    assert row["status"]=="PASS_ACTUAL_COMPLEX_SPARSE_LU"
    assert row["backend"]==backend and row["ranks"]==ranks
    assert row["selected_factor_solver"]==backend
    assert row["matrix_size"]==384 and row["KSP_reason"]>0
    assert row["matrix_type"]=="mpiaij" if ranks>1 else row["matrix_type"] in ("seqaij","mpiaij")
    assert row["relative_residual_l2"]<1e-11
    assert row["relative_solution_l2_error"]<1e-11
    results.append(row)
summary={
    "status":"PASS_ALL_AVAILABLE_PETSC_LU_BACKENDS",
    "cases":len(results),"mumps_rank_counts":[1,2,4],
    "superlu_dist_rank_counts":[1,2,4],
    "sequential_superlu_rank_counts":[1],
    "max_relative_residual":max(x["relative_residual_l2"] for x in results),
    "max_relative_solution_error":max(x["relative_solution_l2_error"] for x in results),
    "not_a_benchmark":True,
    "not_a_SCRIMP_examples_or_GetFEM_qualification":True,
}
(directory/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
archive=Path(str(directory)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as zipped:
    for file in sorted(directory.glob("*")):
        if file.is_file() and file.suffix in (".json",".log"):
            zipped.write(file,arcname=file.name)
print(json.dumps(summary,indent=2))
print("SHAREABLE ORIGINAL PETSC LU PREFLIGHT ZIP:",archive)
PY
echo "COMPLETE: PETSc native direct-solver preflight evidence in $out"
