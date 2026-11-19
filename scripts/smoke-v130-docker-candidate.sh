#!/usr/bin/env bash
# Minimal native ABI, GetFEM FEM and MPI smoke for an UNQUALIFIED candidate.
# Does not modify or delete scrimp-mpi:local or any local Conda environment.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${SCRIMP_DOCKER_CANDIDATE_IMAGE:-scrimp-mpi:v1.3.0-rc1}"
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
  DOCKER=(docker)
else
  DOCKER=(sudo docker)
fi

cd "$ROOT"
mkdir -p sandbox/outputs
LOG="$ROOT/sandbox/outputs/v130-docker-smoke.log"

{
  echo "SCRIMP v1.3.0 candidate native and MPI smoke: $IMAGE"
  echo "No changes to scrimp-mpi:local or host Conda."
  "${DOCKER[@]}" image inspect --format \
    '{{.Id}}  source-revision={{index .Config.Labels "org.opencontainers.image.revision"}}' \
    "$IMAGE"

  echo "=== Native GetFEM simplicial descriptors, tested k=1/2 ==="
  "${DOCKER[@]}" run --rm -i "$IMAGE" python - <<'PY'
import getfem as gf
from petsc4py import PETSc

assert gf.util_mpi_parallelism_level() == 2
assert PETSc.ScalarType(1j).imag == 1
counts = {
    "FEM_RTK(2,0)": 3,
    "FEM_RTK(3,0)": 4,
    "FEM_BDMK(2,1)": 6,
    "FEM_BDMK(3,1)": 12,
    "FEM_BDMK(3,2)": 30,
    "FEM_NEDELEC1(2,1)": 8,
    "FEM_NEDELEC1(3,1)": 20,
    "FEM_NEDELEC1(3,2)": 45,
    "FEM_NEDELEC2(2,1)": 6,
    "FEM_NEDELEC2(2,2)": 12,
    "FEM_NEDELEC2(3,1)": 12,
    "FEM_NEDELEC2(3,2)": 30,
}
for fem_name, expected in counts.items():
    observed = int(gf.Fem(fem_name).nbdof())
    print(f"{fem_name}: {observed} (expected {expected})", flush=True)
    assert observed == expected, (fem_name, observed, expected)

# Do not accidentally advertise or silently enable Nedelec II k>=3.
try:
    gf.Fem("FEM_NEDELEC2(3,3)")
except Exception:
    print("Nedelec-II k>=3 guard: present", flush=True)
else:
    raise AssertionError("Unexpected Nedelec-II k=3 support: not qualified")
print("PASS: candidate FEM descriptors and complex PETSc", flush=True)
PY

  echo "=== Existing stack/MUMPS smoke: 1 MPI rank ==="
  "${DOCKER[@]}" run --rm "$IMAGE" python sandbox/smoke_modern_stack.py

  echo "=== Existing stack/MUMPS smoke: 2 MPI ranks ==="
  "${DOCKER[@]}" run --rm "$IMAGE" \
    mpiexec --oversubscribe -n 2 python -u sandbox/smoke_modern_stack.py

  echo "=== Persist explicit Conda spec and approved GetFEM patch checksums ==="
  "${DOCKER[@]}" run --rm "$IMAGE" \
    cat docker-resolved-spec.txt > "$ROOT/sandbox/outputs/v130-docker-resolved-spec.txt"
  "${DOCKER[@]}" run --rm "$IMAGE" \
    sha256sum scripts/patches/getfem-5.5/*.patch \
    > "$ROOT/sandbox/outputs/v130-getfem-patch-sha256.txt"
  echo "PASS: candidate native, 1/2-rank MPI and FEM descriptor smoke."
} 2>&1 | tee "$LOG"
