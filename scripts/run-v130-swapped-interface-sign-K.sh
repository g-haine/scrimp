#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# Diagnostic ONLY: mathematical outward-normal convention versus
# reconstructed pre-5042705 disk-based sign corrections for heat=2/wave=1.
# Real GetFEM/PETSc K, unchanged current example in both Docker runs.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
 --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
 echo "Pinned GetFEM/PETSc image required, got: $revision" >&2; exit 2
}
python3 -m py_compile examples/heat_wave.py \
 sandbox/probe_v130_swapped_interface_sign_K.py \
 tests/test_v130_heat_wave_outward_signs.py
python3 tests/test_v130_heat_wave_outward_signs.py
python3 tests/test_v130_interface_topology.py
python3 tests/test_v130_interface_ownership.py
python3 tests/test_v130_interface_assembly_plan.py
python3 tests/test_v130_interface_deferred.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-swapped-interface-sign-K-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-swapped-sign-contrast-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "GetFEM/PETSc heat=2 wave=1: outward normal vs reconstructed legacy signs"
echo "Diagnostic stiffness K ONLY; no DAE/CN/energy qualification"
echo "Evidence: $out"

for variant in published outward; do
 mkdir "$out/$variant" "$out/work_$variant"
 echo "ASSEMBLE swapped-sign variant=$variant"
 if ! "${docker_cmd[@]}" run --rm \
   --user "$(id -u):$(id -g)" --env HOME=/tmp --env MPLBACKEND=Agg \
   --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
   --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
   --env SCRIMP_V130_AUTO_INTERFACE_PILOT=1 \
   --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
   --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
   --env SCRIMP_MATRIX_MEMORY_MODE=compact \
   --env SCRIMP_MUMPS_RESET_POLICY=on_change \
   --env SCRIMP_GETFEM_RHS_MODE=auto \
   --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
   --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
   --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
   --volume "$root/sandbox/probe_v130_swapped_interface_sign_K.py:/workspace/scrimp/sandbox/probe_v130_swapped_interface_sign_K.py:ro" \
   --volume "$out/$variant:/artifacts:rw" \
   --volume "$out/work_$variant:/scrimp-work:rw" --workdir /scrimp-work \
   "$image" timeout -k 5s "${SCRIMP_V130_SWAPPED_SIGN_K_TIMEOUT:-240s}" \
   mpiexec --oversubscribe --bind-to none -n 1 \
   python -u /workspace/scrimp/sandbox/probe_v130_swapped_interface_sign_K.py \
     --variant "$variant" --out /artifacts/stiffness_audit.json \
   > "$out/$variant.host.log" 2>&1; then
  echo "FAIL native swapped assembly variant=$variant: $out" >&2
  tail -n 100 "$out/$variant.host.log" >&2
  exit 1
 fi
done

python3 - "$out" <<'PY'
import json, sys, zipfile
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz

root = Path(sys.argv[1])
info = {
    name: json.loads((root / name / "stiffness_audit.json").read_text())
    for name in ("published", "outward")
}
K_old = load_npz(root / "published" / "stiffness_audit.npz").tocsr()
K_new = load_npz(root / "outward" / "stiffness_audit.npz").tocsr()
for K in (K_old, K_new):
    K.sum_duplicates()
    K.sort_indices()
diff = (K_old - K_new).tocsr()
diff.eliminate_zeros()
intervals = info["outward"]["intervals"]

def block(matrix, row, col):
    rs, rn = intervals[row]
    cs, cn = intervals[col]
    return matrix[rs:rs+rn, cs:cs+cn]

blocks = (("e_Q", "U_T"), ("U_w", "Y_T"))
changes = {}
remainder = diff.tolil()
for row, col in blocks:
    old = block(K_old, row, col)
    new = block(K_new, row, col)
    changed = block(diff, row, col)
    doubled = (old + new).tocsr()
    doubled.eliminate_zeros()
    changes[row + "<-" + col] = {
        "old_nnz": int(old.nnz),
        "new_nnz": int(new.nnz),
        "changed_nnz": int(changed.nnz),
        "old_frobenius_norm": float(np.linalg.norm(old.data)),
        "new_frobenius_norm": float(np.linalg.norm(new.data)),
        "exact_opposite_sign": doubled.nnz == 0,
    }
    rs, rn = intervals[row]
    cs, cn = intervals[col]
    remainder[rs:rs+rn, cs:cs+cn] = 0
remainder = remainder.tocsr()
remainder.eliminate_zeros()
valid = (
    K_old.shape == K_new.shape
    and info["published"]["variant"] == "published"
    and info["outward"]["variant"] == "outward"
    and info["outward"]["runtime_only_physical_sign_overrides"] == []
    and all(x["changed_nnz"] > 0 and x["exact_opposite_sign"]
            for x in changes.values())
    and remainder.nnz == 0
)
result = {
    "schema": "scrimp-v130-swapped-heat-wave-outward-sign-assembly-v1",
    "status": (
        "PASS_SWAPPED_OUTWARD_SIGN_CHANGE_LOCALIZED_IN_TWO_INTERFACE_BLOCKS"
        if valid else "FAIL_SWAPPED_SIGN_CHANGE_NOT_AS_EXPECTED"
    ),
    "source_example_modified_at_runtime": False,
    "current_example": "outward",
    "historical_disk_signs_reconstructed_in_sandbox_only": True,
    "full_K_shape": list(K_new.shape),
    "comparison_difference_nnz": int(diff.nnz),
    "unexplained_difference_nnz": int(remainder.nnz),
    "changed_blocks": changes,
    "DAE_initialized": False,
    "energy_balance_qualified": False,
    "MPI2_qualified": False,
}
(root / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
with zipfile.ZipFile(str(root)+".zip", "x", compression=zipfile.ZIP_DEFLATED) as z:
    for file in sorted(root.rglob("*")):
        if file.is_file():
            z.write(file, arcname=str(file.relative_to(root)))
print(json.dumps(result, indent=2))
print("SHAREABLE SWAPPED HEAT-WAVE SIGN EVIDENCE:", str(root)+".zip")
if not valid:
    raise SystemExit(1)
PY
echo "COMPLETE swapped Heat-Wave sign assembly diagnostic: $out"
