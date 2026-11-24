#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Run native PETSc SuperLU_DIST permutations on ONE original CN Mat.
# MC64 is unavailable in the pinned redistribution-safe build (HSL).
# AWPM also requires CombBLAS; skip build-disabled variants with evidence
# rather than aborting before testing usable row/column permutations.
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
# The pinned image provably lacks MC64/HSL. DO NOT invoke
# LargeDiag_MC64: it calls abort(6) in libsuperlu_dist.so.
# LargeDiag_AWPM is optional and needs the CombBLAS build feature.
# All cases run in separate processes so one unsupported backend
# cannot silently invalidate other factorization evidence.
cases=(
 "default:default:baseline"
 "default:MMD_AT_PLUS_A:default_mmd"
 "default:MMD_ATA:default_ata"
 "NOROWPERM:MMD_AT_PLUS_A:norow_mmd"
 "NOROWPERM:MMD_ATA:norow_ata"
 "NOROWPERM:METIS_AT_PLUS_A:norow_metis"
 "NOROWPERM:NATURAL:norow_natural"
 "LargeDiag_AWPM:MMD_AT_PLUS_A:awpm_mmd"
 "LargeDiag_AWPM:METIS_AT_PLUS_A:awpm_metis"
)
echo "SKIP rowperm=LargeDiag_MC64: HSL MC64 was not compiled into pinned SuperLU_DIST"
for item in "${cases[@]}"; do
 IFS=: read -r rowperm colperm label <<<"$item"
 echo "RUN $label row=$rowperm col=$colperm"
 rc=0
 "${docker_cmd[@]}" run --rm --user "$(id -u):$(id -g)" \
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
 > "$out/$label.host.log" 2>&1 || rc=$?
 printf '%s\n' "$rc" > "$out/$label.exitcode"
 if [[ "$rc" -ne 0 ]]; then
  if [[ "$rowperm" == "LargeDiag_AWPM" ]] && \
    grep -Eiq 'CombBLAS.*(not available|not enabled|disabled)|AWPM.*(unavailable|not available)' "$out/$label.host.log"; then
   echo "SKIP $label: this SuperLU_DIST build lacks CombBLAS/AWPM; original log retained"
  else
   echo "FAIL subprocess $label exit=$rc; preserving log and continuing other cases" >&2
   tail -n 35 "$out/$label.host.log" >&2
  fi
 fi
done
python3 - "$out" <<'PY'
import json
import re
import sys
import zipfile
from pathlib import Path

out = Path(sys.argv[1])
identity = "c94c69c4e4307d7a5aa61d9eba186782c0495b49431cadcb7281d7a7b5bffaff"
cases = (
    ("default", "default", "baseline"),
    ("default", "MMD_AT_PLUS_A", "default_mmd"),
    ("default", "MMD_ATA", "default_ata"),
    ("NOROWPERM", "MMD_AT_PLUS_A", "norow_mmd"),
    ("NOROWPERM", "MMD_ATA", "norow_ata"),
    ("NOROWPERM", "METIS_AT_PLUS_A", "norow_metis"),
    ("NOROWPERM", "NATURAL", "norow_natural"),
    ("LargeDiag_AWPM", "MMD_AT_PLUS_A", "awpm_mmd"),
    ("LargeDiag_AWPM", "METIS_AT_PLUS_A", "awpm_metis"),
)
items = {}
unavailable = {}
unexpected = {}
for row, col, label in cases:
    logfile = out / f"{label}.host.log"
    exitfile = out / f"{label}.exitcode"
    if not logfile.is_file() or not exitfile.is_file():
        raise RuntimeError(f"Missing native process provenance for {label}")
    log = logfile.read_text(encoding="utf-8", errors="replace")
    code = int(exitfile.read_text().strip())
    missing_combblas = bool(re.search(
        r"CombBLAS.*(?:not available|not enabled|disabled)"
        r"|AWPM.*(?:unavailable|not available)", log, re.I
    ))
    if missing_combblas and row == "LargeDiag_AWPM":
        items[label] = {
            "rowperm": row, "colperm": col, "exit_code": code,
            "status": "NOT_RUN_TO_COMPLETION_AWPM_NEEDS_COMBBLAS",
            "factorization_validated": False,
        }
        unavailable[label] = items[label]["status"]
        continue
    if code != 0:
        items[label] = {
            "rowperm": row, "colperm": col, "exit_code": code,
            "status": "UNEXPECTED_NATIVE_PROCESS_FAILURE_NOT_A_NUMERICAL_RESULT",
            "last_log_lines": log.splitlines()[-12:],
        }
        unexpected[label] = items[label]["status"]
        continue
    resultfile = out / f"{label}.json"
    if not resultfile.is_file():
        items[label] = {
            "rowperm": row, "colperm": col, "exit_code": code,
            "status": "UNEXPECTED_MISSING_JSON",
        }
        unexpected[label] = items[label]["status"]
        continue
    data = json.loads(resultfile.read_text())
    if (data["input_CSR_SHA256"] != identity
            or data["output_CSR_SHA256"] != identity
            or not data["tiny_pivot_replacement_disabled"]
            or data["rowperm_requested"] != row
            or data["colperm_requested"] != col):
        raise RuntimeError(f"Unqualified or altered CN matrix in {label}")
    items[label] = {
        "rowperm": row,
        "colperm": col,
        "exit_code": code,
        "status": data["status"],
        "KSP": data["KSP_reason"],
        "PC": data["PC_failed"],
        "relative_linear_residual": data["relative_linear_residual"],
        "relative_state_error": data["relative_state_error"],
        "selected_options_verified": bool(
            data["selected_rowperm_verified"]
            and data["selected_colperm_verified"]
        ),
        "actual_factor_options": data.get("factor_options_excerpt"),
    }

passing = [label for label, entry in items.items()
           if entry["status"] ==
           "PASS_IDENTICAL_CN_MATRIX_SOLVED_WITH_THIS_ORDERING"]
summary = {
    "schema": "scrimp-v130-superlu-CN-permutation-sweep-v2",
    "status": (
        "DIAGNOSTIC_PARTIAL_UNEXPECTED_PROCESS_FAILURE"
        if unexpected else
        "DIAGNOSTIC_COMPLETE_WITH_UNAVAILABLE_BUILD_FEATURES"
    ),
    "original_source_failure": "HSL MC64 unavailable, abort(6)",
    "MC64_cases": {
        "mc64_mmd": "NOT_ATTEMPTED_HSL_NONREDISTRIBUTABLE",
        "mc64_mmd_ata": "NOT_ATTEMPTED_HSL_NONREDISTRIBUTABLE",
    },
    "optional_AWPM_unavailable": unavailable,
    "unexpected_native_failures": unexpected,
    "identical_original_CSR_SHA256": identity,
    "successful_LU_orderings": passing,
    "cases": items,
    "tiny_pivot_replacement_or_regularization": False,
    "heat_wave_CN_trajectory_qualified": False,
    "original_sources_unchanged": True,
}
(out / "comparison.json").write_text(
    json.dumps(summary, indent=2) + "\n", encoding="utf-8"
)
archive = Path(str(out) + ".zip")
with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED) as z:
    for file in sorted(out.glob("*")):
        if file.is_file():
            z.write(file, arcname=file.name)
print(json.dumps(summary, indent=2))
print("ORIGINAL PETSc SUPERLU DIST PERMUTATION EVIDENCE:", archive)
if unexpected:
    raise SystemExit("Unexpected native process failures: see archived logs")
PY
echo "COMPLETE original SuperLU_DIST permutation sweep: $out"
