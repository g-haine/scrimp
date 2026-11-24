#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Reorder original Heat-Wave CN equations via exact bipartite matching
# and replay with genuine MUMPS/SuperLU_DIST, no matrix modification.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 1 ]] || {
 echo "Usage: $0 sandbox/outputs/.../physical/first_cn_J_001.petscbin" >&2; exit 2
}
matrix="$(realpath "$1")"
[[ -s "$matrix" && "$(basename "$matrix")" == "first_cn_J_001.petscbin" ]] || {
 echo "Require pinned original Heat-Wave CN snapshot" >&2;exit 2
}
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
 --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
 echo "Original pinned GetFEM/PETSc revision required: $revision" >&2;exit 2
}
python3 -m py_compile sandbox/probe_v130_hw_cn_exact_row_matching.py
python3 tests/test_v130_hw_cn_exact_row_matching.py
python3 tests/test_v130_hw_cn_exact_row_matching.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-hw-cn-row-matching-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-cn-exact-row-matching-v1","source":"%s","image_revision":"%s","original_binary":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" "$matrix" > "$out/manifest.json"
echo "One identical original 7689x7689 CN matrix; structural matching on true nonzeros"
echo "Exact left row permutation PAx=Pb; no HSL, no CombBLAS, no pivot replacement"
echo "Independent native PETSc KSP/PC/MUMPS/SuperLU_DIST process per candidate, MPI1"
echo "EVIDENCE: $out"
cases=(
 "original:superlu_dist:NATURAL:original_natural"
 "row_matching:mumps:default:matching_mumps"
 "row_matching:superlu_dist:NATURAL:matching_natural"
 "row_matching:superlu_dist:MMD_AT_PLUS_A:matching_mmd"
 "row_matching:superlu_dist:METIS_AT_PLUS_A:matching_metis"
)
for item in "${cases[@]}"; do
 IFS=: read -r mode solver colperm label <<<"$item"
 echo "RUN $label"
 rc=0
 "${docker_cmd[@]}" run --rm --user "$(id -u):$(id -g)" \
  --env HOME=/tmp --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
  --volume "$(dirname "$matrix"):/original:ro" \
  --volume "$root/sandbox/probe_v130_hw_cn_identical_lu.py:/replay_base.py:ro" \
  --volume "$root/sandbox/probe_v130_hw_cn_exact_row_matching.py:/pilot.py:ro" \
  --volume "$out:/artifacts:rw" --workdir /artifacts \
  "$image" timeout -k 5s "${SCRIMP_V130_HW_MATCHING_TIMEOUT:-180s}" \
  mpiexec --oversubscribe --bind-to none -n 1 \
  python -u /pilot.py \
  --matrix /original/first_cn_J_001.petscbin \
  --mode "$mode" --solver "$solver" --colperm "$colperm" \
  --output "/artifacts/$label.json" \
  > "$out/$label.host.log" 2>&1 || rc=$?
 printf '%s\n' "$rc" > "$out/$label.exitcode"
 if [[ "$rc" -ne 0 ]]; then
  echo "UNEXPECTED native process failure $label exit=$rc; continuing for diagnostics" >&2
  tail -n 55 "$out/$label.host.log" >&2
 fi
done
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
cases=(
 ("original","superlu_dist","NATURAL","original_natural"),
 ("row_matching","mumps","default","matching_mumps"),
 ("row_matching","superlu_dist","NATURAL","matching_natural"),
 ("row_matching","superlu_dist","MMD_AT_PLUS_A","matching_mmd"),
 ("row_matching","superlu_dist","METIS_AT_PLUS_A","matching_metis"),
)
entries={}
unexpected={}
for mode,solver,colperm,label in cases:
 log=root/f"{label}.host.log"
 code=root/f"{label}.exitcode"
 if not log.is_file() or not code.is_file():
  raise RuntimeError("Missing provenance "+label)
 rc=int(code.read_text())
 file=root/f"{label}.json"
 if rc!=0 or not file.is_file():
  entries[label]={"status":"UNEXPECTED_NATIVE_PROCESS_FAILURE",
                   "exit_code":rc,
                   "log_tail":log.read_text(errors="replace").splitlines()[-10:]}
  unexpected[label]=rc
  continue
 data=json.loads(file.read_text())
 assert data["original_CSR_SHA256_before"]==data["original_CSR_SHA256_after"]
 assert data["original_CSR_SHA256_before"]=="c94c69c4e4307d7a5aa61d9eba186782c0495b49431cadcb7281d7a7b5bffaff"
 assert data["mode"]==mode and data["solver"]==solver
 assert data["colperm_requested"]==colperm
 assert data["tiny_pivot_replacement"] is False
 assert data["diagonal_shift_or_regularization"] is False
 assert data["relative_PA_action_error"]<1e-12
 assert data["matched_graph"]["bipartite_matching_complete"]
 assert data["matched_graph"]["matched_numeric_zero_diagonal"]==0
 assert data["main_CN_TS_completed"] is False
 entries[label]={
  "status":data["status"],
  "KSP_reason":data["KSP_convergence_reason"],
  "PC_failed":data["PC_failure_reason"],
  "original_residual":data["relative_residual_in_ORIGINAL_equations"],
  "known_state_error":data["relative_known_state_error"],
  "matched_graph":data["matched_graph"],
 }
pass_cases=[k for k,v in entries.items()
            if v["status"]=="PASS_PERMUTED_JACOBIAN_KNOWN_STATE_LU_ONLY"]
summary={
 "schema":"scrimp-v130-original-heat-wave-CN-row-matched-LU-v1",
 "status":"DIAGNOSTIC_PARTIAL_UNEXPECTED_PROCESS_FAILURE" if unexpected else
          "DIAGNOSTIC_COMPLETE_NO_FULL_TS_CLAIM",
 "successful_known_state_LU_cases":pass_cases,
 "cases":entries,
 "unexpected_process_exits":unexpected,
 "strict_original_matrix_preserved":True,
 "no_chemical_or_numerical_regularization":True,
 "full_CN_TS_and_energy_not_qualified":True,
}
(root/"comparison.json").write_text(json.dumps(summary,indent=2)+"\n")
with zipfile.ZipFile(str(root)+".zip","x",compression=zipfile.ZIP_DEFLATED) as z:
 for file in sorted(root.glob("*")):
  if file.is_file():z.write(file,arcname=file.name)
print(json.dumps(summary,indent=2))
print("ORIGINAL CN STRUCTURAL MATCHING LU EVIDENCE:",str(root)+".zip")
if unexpected:raise SystemExit("Some native process failures: logs preserved")
PY
echo "COMPLETE original Heat-Wave exact structural row-matching replay: $out"
