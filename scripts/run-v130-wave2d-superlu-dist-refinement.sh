#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# See LICENSE for license information.
#
# Direct sparse PETSc refinement only. The original first-CN Jacobian
# never changes. FULL SCRIMP wave physics remains 3/5 until independently
# qualified with true TS/DAE, 4 port powers and original 4-field VTK.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0 (no positional args)" >&2; exit 2; }
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
  echo "Not the native pinned SCRIMP GetFEM/PETSc candidate: $revision" >&2
  exit 2
}
source_zip="$root/sandbox/outputs/v130-wave2d-first-cn-jacobian-574d9c71-NKqIZ0BZ.zip"
[[ -r "$source_zip" ]] || {
  echo "Missing ORIGINAL first-CN factor matrix archive: $source_zip" >&2
  exit 2
}
python3 -m py_compile \
  sandbox/replay_v130_wave2d_superlu_dist_refinement.py \
  tests/test_v130_wave2d_refinement.py
python3 tests/test_v130_wave2d_refinement.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-wave2d-superlu-dist-refinement-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
mkdir -p "$out/input"
python3 - "$source_zip" "$out/input/matrix.npz" <<'PY'
import io,hashlib,sys,zipfile
import numpy as np
z,target=sys.argv[1:]
with zipfile.ZipFile(z) as archive:
    blob=archive.read("mumps/first_main_CN_factor.npz")
with np.load(io.BytesIO(blob),allow_pickle=False) as opened:
    arrays={k:opened[k].copy() for k in opened.files}
hash_=hashlib.sha256()
for key in sorted(arrays):
    hash_.update(key.encode())
    hash_.update(np.ascontiguousarray(arrays[key]).tobytes())
assert hash_.hexdigest()=="24346856f358e4f81cf332ebc974e329dc084c60f5d6b98bfbcd788fdff37f54", (
    "UNQUALIFIED source Jacobian; refuse to use unverified P,K,M")
with open(target,"wb") as f: f.write(blob)
print("PINNED original 2230x2230 P,K,M:",hash_.hexdigest())
PY
printf '{"schema":"scrimp-v130-2d-superlu-dist-CN-PETSc-refinement-v1","source":"%s","image":"%s","revision":"%s","original_PMK_archive":"%s","original_PMK_SHA256":"24346856f358e4f81cf332ebc974e329dc084c60f5d6b98bfbcd788fdff37f54","candidate_only":true,"existing_five_backend_physics_gate":"3/5"}\n' \
  "$(git rev-parse HEAD)" "$image" "$revision" "$(basename "$source_zip")" > "$out/manifest.json"
echo "SuperLU_DIST original PETSc 2D wave CN matrix, same SHA256, MPI1/2"
echo "Internal iterative refinement vs RIGHT-preconditioned true-residual GMRES"
echo "BOTH require opt-in ReplaceTinyPivot. NO physical model modifications."
echo "STILL REQUIRE: true backward error <1e-11 and exact solution error <1e-7"
echo "EVIDENCE: $out"
for refinement in internal_refine right_gmres; do
  for ranks in 1 2; do
    name="${refinement}_mpi${ranks}"
    echo "RUN $name"
    mkdir -p "$out/$name" "$out/work_$name"
    rc=0
    "${docker_cmd[@]}" run --rm \
      --user "$(id -u):$(id -g)" \
      --env HOME=/tmp --env OPENBLAS_NUM_THREADS=1 --env OMP_NUM_THREADS=1 \
      --volume "$out/input:/input:ro" \
      --volume "$root/sandbox/replay_v130_wave2d_superlu_dist_refinement.py:/workspace/scrimp/sandbox/replay_v130_wave2d_superlu_dist_refinement.py:ro" \
      --volume "$out/$name:/artifacts:rw" \
      --volume "$out/work_$name:/scrimp-work:rw" \
      --workdir /scrimp-work \
      "$image" timeout -k 5s "${SCRIMP_V130_WAVE2D_REFINEMENT_TIMEOUT:-120s}" \
      mpiexec --oversubscribe --bind-to none -n "$ranks" \
      python -u /workspace/scrimp/sandbox/replay_v130_wave2d_superlu_dist_refinement.py \
        --matrix /input/matrix.npz \
        --refinement "$refinement" \
        --output /artifacts \
      > "$out/$name.host.log" 2>&1 || rc=$?
    echo "$rc" > "$out/$name.exitcode"
    if [[ "$rc" != 0 ]]; then
      echo "UNQUALIFIED $name exit=$rc, original evidence retained" >&2
      grep -aE 'KSP_reason|PC_failed|Traceback|error code|relative_' "$out/$name.host.log" \
        | tail -n 12 >&2 || true
    fi
  done
done
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
folder=Path(sys.argv[1])
cases={}
for refinement in ("internal_refine","right_gmres"):
    for mpi in (1,2):
        name=f"{refinement}_mpi{mpi}"
        case=folder/name
        exitcode=int((folder/f"{name}.exitcode").read_text())
        log=(folder/f"{name}.host.log").read_text(errors="replace")
        expected=["mat_superlu_dist_replacetinypivot"]
        if refinement=="internal_refine":
            expected.append("mat_superlu_dist_iterrefine")
        unused=[line for line in log.splitlines()
                if "Option left:" in line and any(v in line for v in expected)]
        got=all((case/f"rank{rank}.json").exists() for rank in range(mpi))
        record={
            "exit_code":exitcode,"rank_json_present":got,
            "unused_relevant_PETSc_options":unused,
            "native_refinement_pass":False,
        }
        if got:
            data=[json.loads((case/f"rank{rank}.json").read_text())
                  for rank in range(mpi)]
            actual_solver_type="preonly" if refinement=="internal_refine" else "gmres"
            accepted=all(
                entry.get("status")=="PASS_INDEPENDENT_ORIGINAL_PETSC_PC_REFINED_REPLAY_ONLY"
                and entry.get("MPI")==mpi and entry.get("rank")==rank
                and entry.get("mode")==refinement
                and entry.get("factorization_passed_all_ranks")
                and entry.get("original_full_P_M_K_SHA256") ==
                    "24346856f358e4f81cf332ebc974e329dc084c60f5d6b98bfbcd788fdff37f54"
                and entry.get("PETSc_owned_CSR_graph_and_values_identical_to_archive")
                and entry.get("original_J_not_modified")
                and entry.get("factor_solver")=="superlu_dist"
                and entry.get("ksp_type")==actual_solver_type
                and entry.get("replace_tiny_pivot_requested")
                and entry.get("internal_iterative_refinement_requested")==(refinement=="internal_refine")
                and entry.get("KSP_reason", -1)>0
                and entry.get("relative_backward_residual") is not None
                and entry["relative_backward_residual"]<1e-11
                and entry.get("relative_known_solution_error") is not None
                and entry["relative_known_solution_error"]<1e-7
                and (refinement!="right_gmres" or entry.get("ksp_iterations",0)>0)
                for rank,entry in enumerate(data)
            )
            record.update({
                "KSP_reasons":[x.get("KSP_reason") for x in data],
                "KSP_iterations":[x.get("ksp_iterations") for x in data],
                "PC_failed":[x.get("PC_failed") for x in data],
                "max_backward_error":max(
                    x.get("relative_backward_residual") or 0 for x in data),
                "max_known_solution_error":max(
                    x.get("relative_known_solution_error") or 0 for x in data),
                "KSP_PC_side":[x.get("ksp_pc_side") for x in data],
            })
            record["native_refinement_pass"]=bool(
                accepted and not unused and exitcode==0
            )
        cases[name]=record
qualified_methods={
    method:all(cases[f"{method}_mpi{r}"]["native_refinement_pass"]
               for r in (1,2))
    for method in ("internal_refine","right_gmres")
}
any_candidate=any(qualified_methods.values())
status=("PASS_AT_LEAST_ONE_STRICT_REFINEMENT_PETSC_MPI1_MPI2_ONLY"
        if any_candidate else
        "FAIL_ALL_STRICT_REFINEMENT_PETSC_CANDIDATES")
report={
    "status":status,
    "qualified_strict_methods_on_both_ranks":qualified_methods,
    "cases":cases,
    "default_superlu_dist_still_not_qualified":True,
    "full_ten_step_wave_2d_physics_qualified_for_superlu_dist":False,
    "original_three_backends_still_PASS":True,
    "wave_installation_gate":"3/5",
}
(folder/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
archive=Path(str(folder)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for f in sorted(folder.rglob("*")):
        if f.is_file():z.write(f,arcname=str(f.relative_to(folder)))
print(json.dumps({**report,"archive":str(archive)},indent=2),flush=True)
if not any_candidate:
    raise SystemExit(1)
PY
echo "CANDIDATE_ONLY: original physical 10-step CN/4-port/VTK gate still FAIL 3/5"
