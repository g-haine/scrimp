#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# Actual GetFEM/PETSc sparse stiffness parity, MPI1, strict pinned image.
# Same unchanged published Heat-Wave example and interface brick mapping.
# Only the oriented facet detection differs between the two runs.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
rev="$("${docker_cmd[@]}" image inspect \
 --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$rev" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
 echo "Pinned original GetFEM/PETSc image required, got: $rev" >&2
 exit 2
}

python3 -m py_compile \
 scrimp/interface.py \
 sandbox/v130_hw_interface_geometry.py \
 sandbox/probe_v130_hw_dual_sided_assembly.py \
 sandbox/probe_v130_generic_interface_matrix_parity.py
python3 tests/test_v130_interface_topology.py
python3 tests/test_v130_interface_ownership.py
python3 tests/test_v130_interface_assembly_plan.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-generic-interface-K-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-generic-facet-matrix-parity-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$rev" > "$out/manifest.json"
echo "Generic versus reference oriented facets: FULL Heat-Wave stiffness matrix"
echo "Stage: actual GetFEM/PETSc assembly only, not DAE or CN"
echo "EVIDENCE: $out"

for geometry in reference generic; do
 mkdir "$out/$geometry" "$out/work_$geometry"
 echo "ASSEMBLE geometry=$geometry"
 if ! "${docker_cmd[@]}" run --rm \
  --user "$(id -u):$(id -g)" --env HOME=/tmp --env MPLBACKEND=Agg \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
  --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
  --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
  --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
  --env SCRIMP_MATRIX_MEMORY_MODE=compact \
  --env SCRIMP_MUMPS_RESET_POLICY=on_change \
  --env SCRIMP_GETFEM_RHS_MODE=auto \
  --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
  --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
  --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
  --volume "$root/sandbox/v130_hw_interface_geometry.py:/workspace/scrimp/sandbox/v130_hw_interface_geometry.py:ro" \
  --volume "$root/sandbox/probe_v130_hw_dual_sided_assembly.py:/workspace/scrimp/sandbox/probe_v130_hw_dual_sided_assembly.py:ro" \
  --volume "$root/sandbox/probe_v130_generic_interface_matrix_parity.py:/workspace/scrimp/sandbox/probe_v130_generic_interface_matrix_parity.py:ro" \
  --volume "$out/$geometry:/artifacts:rw" \
  --volume "$out/work_$geometry:/scrimp-work:rw" --workdir /scrimp-work \
  "$image" timeout -k 5s "${SCRIMP_V130_GENERIC_K_TIMEOUT:-240s}" \
  mpiexec --oversubscribe --bind-to none -n 1 \
  python -u /workspace/scrimp/sandbox/probe_v130_generic_interface_matrix_parity.py \
  --geometry "$geometry" --out /artifacts/stiffness_audit.json \
  > "$out/$geometry.host.log" 2>&1; then
  echo "Assembly FAILED for $geometry; details at $out" >&2
  tail -n 100 "$out/$geometry.host.log" >&2
  exit 1
 fi
done

python3 - "$out" <<'PY'
import json, sys, zipfile
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz

out = Path(sys.argv[1])
reports = {
    name: json.loads((out / name / "stiffness_audit.json").read_text())
    for name in ("reference", "generic")
}
A = load_npz(out / "reference" / "stiffness_audit.npz").tocsr()
B = load_npz(out / "generic" / "stiffness_audit.npz").tocsr()
A.sum_duplicates()
B.sum_duplicates()
difference = (A - B).tocsr()
difference.eliminate_zeros()
max_abs = float(np.max(np.abs(difference.data))) if difference.nnz else 0.
reference_max = float(np.max(np.abs(A.data))) if A.nnz else 0.
abs_tol = 1.e-12
relative_tol = 1.e-13
threshold = abs_tol + relative_tol * reference_max
valid = (
    A.shape == B.shape
    and A.nnz == B.nnz
    and max_abs <= threshold
    and reports["reference"]["stiffness_matrix_full_DOFS"] == A.shape[0]
    and reports["generic"]["stiffness_matrix_full_DOFS"] == B.shape[0]
    and reports["reference"]["mass_matrix_assembled"]
    and reports["generic"]["mass_matrix_assembled"]
    and reports["generic"]["variant"] == "dual_sided"
)
result = {
    "schema": "scrimp-v130-full-native-generic-interface-K-parity-v1",
    "status": "PASS_NATIVE_GENERIC_AND_REFERENCE_ASSEMBLED_K_PARITY_ONLY"
        if valid else "FAIL_GENERIC_AND_REFERENCE_ASSEMBLED_K_DIFFER",
    "full_stiffness_shape": list(A.shape),
    "reference_nnz": int(A.nnz),
    "generic_nnz": int(B.nnz),
    "difference_nnz": int(difference.nnz),
    "max_absolute_difference": max_abs,
    "max_reference_absolute_entry": reference_max,
    "comparison_threshold": threshold,
    "same_original_weak_forms": True,
    "same_original_port_and_brick_region_mapping": True,
    "geometry_selection_only_varied": True,
    "production_auto_assembly_enabled": False,
    "DAE_initialized": False,
    "CN_or_physical_energy_qualified": False,
    "MPI2_qualified": False,
}
(out / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
with zipfile.ZipFile(str(out) + ".zip", "x", compression=zipfile.ZIP_DEFLATED) as z:
    for f in sorted(out.rglob("*")):
        if f.is_file():
            z.write(f, arcname=str(f.relative_to(out)))
print(json.dumps(result, indent=2))
print("SHAREABLE GENERIC FACET K ASSEMBLY EVIDENCE:", str(out) + ".zip")
if not valid:
    sys.exit(1)
PY
echo "COMPLETE full native K parity: $out"
