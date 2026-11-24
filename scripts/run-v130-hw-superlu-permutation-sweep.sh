#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Run different true PETSc SuperLU_DIST permutations on the SAME
# original CN matrix captured in the previous campaign, MPI1.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 1 ]] || {
 echo "Usage: $0 sandbox/outputs/.../physical/first_cn_J_001.petscbin" >&2
 exit 2
}
matrix="$(realpath "$1")"
[[ -s "$matrix" && "$(basename "$matrix")" == first_cn_J_001.petscbin ]] || {
 echo "Require original first CN PETSc binary snapshot: $matrix" >&2
 exit 2
}
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
 --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
 echo "Wrong original GetFEM/PETSc image $revision" >&2; exit 2
}
python3 -m py_compile sandbox/probe_v130_hw_superlu_permutation_sweep.py
python3 tests/test_v130_hw_superlu_permutation_sweep.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-hw-superlu-perm-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-superlu-original-CN-sweep-v1","source":"%s","original_binary":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$matrix" "$revision" > "$out/manifest.json"
echo "SAME original first Heat-Wave CN 7689x7689 matrix / immutable SHA"
echo "SuperLU_DIST row/column permutations only, fresh process per case"
echo "No tiny-pivot replacement, no regularization, no physical TS inference"
echo "EVIDENCE: $out"
# first is the default reproduction; subsequent cases vary only
# documented SuperLU_DIST row/column permutations.
cases=(
 "default:default:baseline"
 "LargeDiag_MC64:MMD_AT_PLUS_A:mc64_mmd"
 "LargeDiag_AWPM:MMD_AT_PLUS_A:awpm_mmd"
 "LargeDiag_AWPM:METIS_AT_PLUS_A:awpm_metis"
 "LargeDiag_MC64:MMD_ATA:mc64_mmd_ata"
 "NOROWPERM:MMD_AT_PLUS_A:norow_mmd"
)
for item in "${cases[@]}"; do
 IFS=: read -r rowperm colperm label <<<"$item"
 echo "RUN $label row=$rowperm col=$colperm"
 if ! "${docker_cmd[@]}" run --rm --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --volume "$(dirname "$matrix"):/original:ro" \
 --volume "$root/sandbox/probe_v130_hw_cn_identical_lu.py:/replay_base.py:ro" \
 --volume "$root/sandbox/probe_v130_hw_superlu_permutation_sweep.py:/probe.py:ro" \
 --volume "$out:/artifacts:rw" --workdir /artifacts \
 "$image" timeout -k 5s "${SCRIMP_V130_SUPERLU_PERM_TIMEOUT:-180s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /probe.py --matrix /original/first_cn_J_001.petscbin \
 --rowperm "$rowperm" --colperm "$colperm" \
 --output "/artifacts/$label.json" \
 > "$out/$label.host.log" 2>&1; then
  echo "UNEXPECTED infrastructure error in $label; logs retained in $out" >&2
  tail -n 100 "$out/$label.host.log" >&2; exit 1
 fi
done
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
out=Path(sys.argv[1])
items={}
for label in ("baseline","mc64_mmd","awpm_mmd","awpm_metis",
              "mc64_mmd_ata","norow_mmd"):
 data=json.loads((out/f"{label}.json").read_text())
 assert data["input_CSR_SHA256"]==data["output_CSR_SHA256"]
 assert data["input_CSR_SHA256"]=="c94c69c4e4307d7a5aa61d9eba186782c0495b49431cadcb7281d7a7b5bffaff"
 assert data["tiny_pivot_replacement_disabled"]
 items[label]={
  "rowperm":data["rowperm_requested"],
  "colperm":data["colperm_requested"],
  "status":data["status"],
  "KSP":data["KSP_reason"],"PC":data["PC_failed"],
  "relative_linear_residual":data["relative_linear_residual"],
  "relative_state_error":data["relative_state_error"],
  "selected_options_verified":(
   data["selected_rowperm_verified"] and data["selected_colperm_verified"]),
  "actual_factor_options":data.get("factor_options_excerpt"),
 }
passing=[label for label,x in items.items() if x["status"].startswith("PASS_")]
summary={
 "schema":"scrimp-v130-superlu-CN-permutation-sweep-v1",
 "status":"DIAGNOSTIC_COMPLETE_NOT_PHYSICAL_QUALIFICATION",
 "identical_original_CSR_SHA256":"c94c69c4e4307d7a5aa61d9eba186782c0495b49431cadcb7281d7a7b5bffaff",
 "successful_LU_orderings":passing,
 "cases":items,
 "tiny_pivot_replacement_or_regularization":False,
 "heat_wave_CN_trajectory_qualified":False,
 "original_sources_unchanged":True,
}
(out/"comparison.json").write_text(json.dumps(summary,indent=2)+"\n")
archive=Path(str(out)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(out.glob("*")):
  if f.is_file(): z.write(f,arcname=f.name)
print(json.dumps(summary,indent=2))
print("ORIGINAL PETSc SUPERLU DIST PERMUTATION EVIDENCE:",archive)
PY
echo "COMPLETE original SuperLU_DIST permutation sweep: $out"
