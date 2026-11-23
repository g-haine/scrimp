#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Native comparison of default FE L2 initialization: qualified
# distributed IJV vs explicitly experimental distributed PETSc COO.
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
 echo "Incorrect GetFEM/PETSc image revision: $revision" >&2; exit 2;
}
python3 -m py_compile sandbox/probe_v130_fe_initial_mpi.py
python3 tests/test_v130_fe_initial_projection.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-fe-initial-coo-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-fe-initial-coo-comparison-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Native GetFEM->PETSc distributed IJV reference vs opt-in COO"
echo "Same original FE L2 default, eight test cases and MPI 1/2"
echo "No release/tag changes; DAE/PDE physics not qualified"
echo "EVIDENCE: $out"
for bridge in ijv coo; do
 for ranks in 1 2; do
  folder="$out/${bridge}-${ranks}r"
  mkdir "$folder"
  if [[ "$bridge" == "coo" ]]; then experimental=1; else experimental=0; fi
  echo "RUN bridge=$bridge ranks=$ranks"
  if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp \
    --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
    --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED="$bridge" \
    --env SCRIMP_FE_INITIAL_COO="$experimental" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
    --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
    --volume "$root/sandbox/probe_v130_fe_initial_mpi.py:/pilot.py:ro" \
    --volume "$root/sandbox/probe_v130_universal_fem_initialization.py:/original_pilot.py:ro" \
    --volume "$folder:/artifacts:rw" \
    --workdir /artifacts \
    "$image" timeout -k 5s "${SCRIMP_V130_FE_MPI_TIMEOUT:-150s}" \
    mpiexec --oversubscribe --bind-to none -n "$ranks" \
    python -u /pilot.py --outdir /artifacts \
    > "$folder/native.host.log" 2>&1; then
      echo "FAIL $bridge MPI$ranks, original logs: $out" >&2
      tail -n 150 "$folder/native.host.log" >&2
      exit 1
  fi
 done
done
python3 - "$out" <<'PY'
import json
import sys
import zipfile
from pathlib import Path
import numpy as np

root = Path(sys.argv[1])
def load(mode, count, rank):
    return json.loads((root / f"{mode}-{count}r" /
                       f"rank-{rank}.json").read_text())

reference = load("ijv", 1, 0)
assert reference["status"] == "PASS_RANK"
assert len(reference["cases"]) == 8
comparison = {}
for mode in ("ijv", "coo"):
    for count in (1, 2):
        for rank in range(count):
            check = load(mode, count, rank)
            assert check["status"] == "PASS_RANK", (mode, count, rank)
            assert check["mpi_ranks"] == count
            assert check["getfem_matrix_mode"] == "distributed"
            assert check["getfem_distributed_bridge"] == mode
            assert check["fe_initial_experimental_coo"] == (mode == "coo")

for name, first in reference["cases"].items():
    discrepancies, physical = [], []
    for mode in ("ijv", "coo"):
        for count in (1, 2):
            for rank in range(count):
                test = load(mode, count, rank)["cases"][name]
                assert test["status"].startswith("PASS_"), (
                    name, mode, count, rank, test
                )
                if name == "lumped":
                    assert test["value"] == first["value"] == 2.5
                    continue
                baseline = np.asarray(first["global_coefficients"],
                                      dtype=np.float64)
                candidate = np.asarray(test["global_coefficients"],
                                       dtype=np.float64)
                assert baseline.shape == candidate.shape
                relative = float(
                    np.linalg.norm(candidate - baseline) /
                    max(float(np.linalg.norm(baseline)), 1e-30)
                )
                discrepancies.append(relative)
                assert relative < 1e-9, (
                    name, mode, count, rank, relative
                )
                l2 = float(test["physical_relative_L2_error"])
                physical.append(l2)
                if name != "gaussian_cg1_scalar":
                    assert l2 < 5e-8, (name, mode, count, rank, l2)
    comparison[name] = {
        "max_relative_coeff_difference_vs_ijv_MPI1":
            max(discrepancies) if discrepancies else 0.0,
        "max_physical_L2_relative_error":
            max(physical) if physical else None,
    }

summary = {
    "schema": "scrimp-v130-fe-initial-coo-comparison-v1",
    "status": "PASS_EXPERIMENTAL_COO_EQUALS_IJV_MPI1_MPI2",
    "cases": comparison,
    "production_coo_default_enabled": False,
    "DAE_and_full_physics_qualified": False,
    "published_sources_unchanged": True,
}
(root / "comparison.json").write_text(
    json.dumps(summary, indent=2) + "\n", encoding="utf-8"
)
archive = Path(str(root) + ".zip")
with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED) as z:
    for path in sorted(root.rglob("*")):
        if path.is_file():
            z.write(path, arcname=str(path.relative_to(root)))
print(json.dumps(summary, indent=2))
print("SHAREABLE ORIGINAL EXPERIMENTAL COO EVIDENCE:", archive)
PY
echo "COMPLETE native FE initial COO comparison: $out"
