#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# Experimental automatic interface: compare ACTUAL PETSc-owned K rows,
# MPI1 versus MPI2, for both disk and annulus assignments.
# No GetFEM/PETSc fake, no rank0-only matrix capture, no serial emulation.
# Strict DAE/CN, parallel export and physical MPI2 power are later gates.
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
  echo "Pinned GetFEM/PETSc image required: $revision" >&2
  exit 2
}
python3 -m py_compile \
  scrimp/interface.py scrimp/interface_discovery.py scrimp/interface_deferred.py \
  scrimp/dphs.py scrimp/port.py scrimp/brick.py examples/heat_wave.py \
  sandbox/probe_v130_auto_interface_mpi_K.py
python3 tests/test_v130_interface_topology.py
python3 tests/test_v130_interface_ownership.py
python3 tests/test_v130_interface_assembly_plan.py
python3 tests/test_v130_interface_deferred.py
python3 tests/test_v130_hw_facet_collector_regions.py
python3 tests/test_v130_heat_wave_outward_signs.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-auto-interface-MPI-K-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-auto-interface-MPI-K-v1","source":"%s","image_revision":"%s"}\n' \
  "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "ACTUAL distributed PETSc sparse K: MPI1 vs MPI2; Heat=1 and Heat=2"
echo "MPI rank-owned CSR rows, full global reconstruction; no DAE/CN"
echo "EVIDENCE: $out"

for heat in 1 2; do
  for ranks in 1 2; do
    label="heat${heat}_mpi${ranks}"
    mkdir "$out/$label" "$out/work_$label"
    echo "ASSEMBLE heat=$heat wave=$((3-heat)) ranks=$ranks"
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
      --volume "$root/sandbox/probe_v130_auto_interface_mpi_K.py:/workspace/scrimp/sandbox/probe_v130_auto_interface_mpi_K.py:ro" \
      --volume "$out/$label:/artifacts:rw" \
      --volume "$out/work_$label:/scrimp-work:rw" --workdir /scrimp-work \
      "$image" timeout -k 5s "${SCRIMP_V130_AUTO_MPI_K_TIMEOUT:-300s}" \
      mpiexec --oversubscribe --bind-to none -n "$ranks" \
      python -u /workspace/scrimp/sandbox/probe_v130_auto_interface_mpi_K.py \
        --heat "$heat" --out /artifacts \
      > "$out/$label.host.log" 2>&1; then
      echo "FAIL native PETSc assembly heat=$heat MPI=$ranks: $out" >&2
      tail -n 120 "$out/$label.host.log" >&2
      printf '%s\n' "FAIL heat=$heat ranks=$ranks" > "$out/FAILED"
      break 2
    fi
  done
done

python3 - "$out" <<'PY'
import json, sys, zipfile
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz, vstack

root = Path(sys.argv[1])
summary = {"schema": "scrimp-v130-auto-interface-actual-MPI-K-parity-v1",
           "status": "FAIL_OR_UNQUALIFIED_AUTO_INTERFACE_MPI2_K",
           "cases": {}, "actual_MPI2_assembled": False,
           "DAE_initialized": False, "CN_energy_qualified_MPI2": False,
           "parallel_exports_qualified": False}
try:
    if (root / "FAILED").exists():
        raise RuntimeError((root / "FAILED").read_text().strip())

    def reconstruct(label, ranks):
        records = []
        for rank in range(ranks):
            path = root / label
            meta = json.loads((path / f"rank{rank}.json").read_text())
            local = load_npz(path / f"rank{rank}.npz").tocsr()
            if (meta["rank"] != rank or meta["mpi_size"] != ranks
                    or meta["local_shape"] != list(local.shape)):
                raise RuntimeError(f"Corrupt rank-owned CSR evidence: {label}:{rank}")
            records.append((meta, local))
        ranges = sorted(
            (int(m["PETSc_owned_rows"][0]),
             int(m["PETSc_owned_rows"][1]), m, K)
            for m, K in records
        )
        global_n = ranges[0][2]["global_shape"][0]
        if (ranges[0][0] != 0 or ranges[-1][1] != global_n
                or any(a[1] != b[0] for a, b in zip(ranges, ranges[1:]))
                or any(
                    x[2]["global_shape"] != [global_n, global_n]
                    or x[3].shape != (x[1]-x[0], global_n)
                    for x in ranges
                )):
            raise RuntimeError(f"Distributed row coverage failure: {label}")
        identity = [{
            "rank": int(x[2]["rank"]),
            "rows": [x[0], x[1]],
            "stored_nnz": int(x[3].nnz),
            "owner_regions": x[2]["oriented_port_regions"],
        } for x in ranges]
        for x in ranges[1:]:
            if x[2]["interface_compiler"] != ranges[0][2]["interface_compiler"]:
                raise RuntimeError(f"MPI ranks disagree on oriented compiler: {label}")
        K = vstack((x[3] for x in ranges), format="csr")
        K.sum_duplicates()
        K.sort_indices()
        return K, identity

    valid = True
    for heat in (1, 2):
        A, owned1 = reconstruct(f"heat{heat}_mpi1", 1)
        B, owned2 = reconstruct(f"heat{heat}_mpi2", 2)
        delta = (A - B).tocsr()
        delta.eliminate_zeros()
        abs_max = float(np.max(np.abs(delta.data))) if delta.nnz else 0.
        ref_max = float(np.max(np.abs(A.data))) if A.nnz else 0.
        threshold = 1e-12 + 1e-11 * ref_max
        strict_CSR = bool(
            np.array_equal(A.indptr, B.indptr)
            and np.array_equal(A.indices, B.indices)
            and np.array_equal(A.data, B.data)
        )
        expected_n = 7689 if heat == 1 else 8056
        passed = bool(
            A.shape == B.shape == (expected_n, expected_n)
            and A.nnz == B.nnz and abs_max <= threshold
            and len(owned2) == 2 and owned2[0]["rows"][1] ==
                owned2[1]["rows"][0]
        )
        summary["cases"][f"heat{heat}_wave{3-heat}"] = {
            "status": "PASS_MPI2_ASSEMBLED_K_PARITY" if passed
                      else "FAIL_MPI2_ASSEMBLED_K_NOT_EQUAL",
            "global_shape": list(A.shape),
            "MPI1_nnz": int(A.nnz),
            "MPI2_nnz": int(B.nnz),
            "difference_nnz": int(delta.nnz),
            "max_abs_difference": abs_max,
            "threshold": threshold,
            "bitwise_CSR_identical": strict_CSR,
            "MPI1_owned_rows": owned1,
            "MPI2_owned_rows": owned2,
        }
        valid = valid and passed
    summary["actual_MPI2_assembled"] = True
    if valid:
        summary["status"] = "PASS_BOTH_AUTO_INTERFACE_ORIENTATIONS_MPI2_K_PARITY"
except Exception as exc:
    summary["failure"] = f"{type(exc).__name__}: {exc}"

(root / "comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
archive = Path(str(root) + ".zip")
with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4, allowZip64=True) as z:
    for f in sorted(root.rglob("*")):
        if f.is_file():
            z.write(f, arcname=str(f.relative_to(root)))
print(json.dumps(summary, indent=2))
print("SHAREABLE ACTUAL MPI2 AUTO INTERFACE ASSEMBLY EVIDENCE:", archive)
if summary["status"] != "PASS_BOTH_AUTO_INTERFACE_ORIENTATIONS_MPI2_K_PARITY":
    raise SystemExit(1)
PY
echo "COMPLETE actual MPI1/MPI2 GetFEM/PETSc automatic interface K audit: $out"
