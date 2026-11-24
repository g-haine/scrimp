#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# EXPERIMENTAL auto interface assembly versus previously verified native K.
# The only user-facing difference is a test-only environment toggle.
# Published Heat-Wave source files remain immutable in both Docker runs.
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
 echo "Pinned original GetFEM/PETSc image required: $rev" >&2; exit 2
}
python3 -m py_compile \
 scrimp/interface.py scrimp/interface_discovery.py \
 scrimp/interface_deferred.py scrimp/interface_ownership.py \
 scrimp/interface_assembly_plan.py scrimp/port.py scrimp/brick.py \
 scrimp/dphs.py sandbox/probe_v130_auto_interface_K.py \
 sandbox/probe_v130_generic_interface_matrix_parity.py
python3 tests/test_v130_interface_topology.py
python3 tests/test_v130_interface_ownership.py
python3 tests/test_v130_interface_assembly_plan.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-auto-interface-K-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-experimental-auto-interface-K-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$rev" > "$out/manifest.json"
echo "Original Heat-Wave, manual reference K vs EXPERIMENTAL auto DPHS compiler"
echo "Actual MPI1 assembly only; no DAE initialization or CN/energy qualification"
echo "EVIDENCE: $out"

for variant in reference auto; do
 mkdir "$out/$variant" "$out/work_$variant"
 if [[ "$variant" == reference ]]; then
  entry=/workspace/scrimp/sandbox/probe_v130_generic_interface_matrix_parity.py
  args=(--geometry generic --out /artifacts/stiffness_audit.json)
  toggle=0
 else
  entry=/workspace/scrimp/sandbox/probe_v130_auto_interface_K.py
  args=(--out /artifacts/stiffness_audit.json)
  toggle=1
 fi
 echo "ASSEMBLE variant=$variant (experimental_compiler=$toggle)"
 if ! "${docker_cmd[@]}" run --rm \
  --user "$(id -u):$(id -g)" --env HOME=/tmp --env MPLBACKEND=Agg \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
  --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
  --env SCRIMP_V130_AUTO_INTERFACE_PILOT="$toggle" \
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
  --volume "$root/sandbox/probe_v130_auto_interface_K.py:/workspace/scrimp/sandbox/probe_v130_auto_interface_K.py:ro" \
  --volume "$out/$variant:/artifacts:rw" \
  --volume "$out/work_$variant:/scrimp-work:rw" --workdir /scrimp-work \
  "$image" timeout -k 5s "${SCRIMP_V130_AUTO_K_TIMEOUT:-240s}" \
  mpiexec --oversubscribe --bind-to none -n 1 \
  python -u "$entry" "${args[@]}" \
  > "$out/$variant.host.log" 2>&1; then
  echo "FAIL native assembly for $variant: $out" >&2
  tail -n 120 "$out/$variant.host.log" >&2
  exit 1
 fi
done

python3 - "$out" <<'PY'
import json, sys, zipfile
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz

path = Path(sys.argv[1])
ref = load_npz(path / "reference" / "stiffness_audit.npz").tocsr()
auto = load_npz(path / "auto" / "stiffness_audit.npz").tocsr()
for matrix in (ref, auto):
    matrix.sum_duplicates()
    matrix.sort_indices()
delta = (ref - auto).tocsr()
delta.eliminate_zeros()
max_diff = float(np.max(np.abs(delta.data))) if delta.nnz else 0.
ref_report = json.loads((path / "reference" / "stiffness_audit.json").read_text())
auto_report = json.loads((path / "auto" / "stiffness_audit.json").read_text())
strict_CSR_equal = (
    np.array_equal(ref.indptr, auto.indptr)
    and np.array_equal(ref.indices, auto.indices)
    and np.array_equal(ref.data, auto.data)
)
valid = (
    strict_CSR_equal
    and ref.shape == auto.shape == (7689, 7689)
    and ref.nnz == auto.nnz == 46935
    and auto_report["status"] == "DIAGNOSTIC_ACTUAL_AUTO_GETFEM_K_ONLY"
    and auto_report["compiler"]["status"] ==
        "EXPERIMENTAL_COMPILED_NOT_PHYSICALLY_QUALIFIED"
    and ref_report["variant"] == "dual_sided"
)
result = {
    "schema": "scrimp-v130-native-experimental-automatic-internal-K-parity-v1",
    "status": "PASS_EXPERIMENTAL_AUTOMATIC_INTERNAL_INTERFACE_K_IDENTICAL"
        if valid else "FAIL_EXPERIMENTAL_AUTO_INTERFACE_K_NOT_IDENTICAL",
    "original_example_source_unchanged": True,
    "auto_compiler_requires_test_only_environment_flag": True,
    "full_stiffness_shape": list(ref.shape),
    "reference_stored_nnz": int(ref.nnz),
    "auto_stored_nnz": int(auto.nnz),
    "strict_CSR_arrays_identical": bool(strict_CSR_equal),
    "numeric_difference_nnz": int(delta.nnz),
    "max_absolute_difference": max_diff,
    "DAE_initialized": False,
    "CN_energy_or_physical_exchange_qualified": False,
    "MPI2_qualified": False,
}
(path / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
with zipfile.ZipFile(str(path) + ".zip", "x", compression=zipfile.ZIP_DEFLATED) as z:
    for file in sorted(path.rglob("*")):
        if file.is_file():
            z.write(file, arcname=str(file.relative_to(path)))
print(json.dumps(result, indent=2))
print("SHAREABLE EXPERIMENTAL AUTO K EVIDENCE:", str(path)+".zip")
if not valid:
    sys.exit(1)
PY
echo "COMPLETE native experimental automatic interface K parity: $out"
