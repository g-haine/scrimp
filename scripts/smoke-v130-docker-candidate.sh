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

  echo "=== Independent user workspace: serial simulation and host-persistent CSV ==="
  # Use only a new temporary test workspace, NEVER the owner's sandbox.
  # The runtime volume mount is distinct from the build context.
  USER_TEST_ROOT="$(mktemp -d -t scrimp-v130-user-smoke.XXXXXXXX)"
  trap 'rm -rf -- "$USER_TEST_ROOT"' EXIT
  python3 scripts/create-scrimp-project.py "$USER_TEST_ROOT/project" \
    --image "$IMAGE"
  "${DOCKER[@]}" run --rm --user "$(id -u):$(id -g)" \
    --env HOME=/tmp --env XDG_CACHE_HOME=/tmp/.cache \
    --volume "$USER_TEST_ROOT/project:/work:rw" --workdir /work \
    --entrypoint /opt/conda/envs/scrimp/bin/python \
    "$IMAGE" -c '
import runpy
import warnings
import numpy as np
# NumPy 2: do not accept a mere logged ComplexWarning.
warnings.simplefilter("error", np.exceptions.ComplexWarning)
runpy.run_path("simulations/first_simulation.py", run_name="__main__")
'
  test -s "$USER_TEST_ROOT/project/outputs/oscillator.csv" || {
    echo "ERROR: first simulation did not create a persistent HOST CSV." >&2
    exit 2
  }
  # Check the host CSV using only Python's standard library (CI host need
  # not have NumPy). Validate the actual equation, not just file existence.
  python3 - "$USER_TEST_ROOT/project/outputs/oscillator.csv" <<'PY'
import csv
import math
import sys

with open(sys.argv[1], newline="", encoding="utf-8") as stream:
    reader = csv.DictReader(stream)
    assert reader.fieldnames == ["time", "q", "p", "H"], reader.fieldnames
    rows = list(reader)
assert len(rows) == 11, f"Expected 11 saved samples, got {len(rows)}"
t = [float(row["time"]) for row in rows]
assert abs(t[0]) < 1e-12 and abs(t[-1] - 0.1) < 1e-10, t
assert all(b > a for a, b in zip(t, t[1:])), t
for row in rows:
    q, p, H = (float(row[k]) for k in ("q", "p", "H"))
    assert all(math.isfinite(x) for x in (q, p, H))
    assert math.isclose(H, 0.5 * 5.0 * q * q + p * p / (2.0 * 1.5),
                        rel_tol=1e-10, abs_tol=1e-12), row
print("PASS: 11 finite host CSV samples, t=0..0.1, H(q,p) consistency")
PY
  echo "PASS: no NumPy ComplexWarning; user output exists on host."

  echo "=== Independent user workspace: two-rank MPI startup ==="
  "${DOCKER[@]}" run --rm --user "$(id -u):$(id -g)" \
    --env HOME=/tmp --env XDG_CACHE_HOME=/tmp/.cache \
    --volume "$USER_TEST_ROOT/project:/work:rw" --workdir /work \
    --entrypoint /opt/conda/envs/scrimp/bin/mpiexec \
    "$IMAGE" --oversubscribe -n 2 \
      /opt/conda/envs/scrimp/bin/python simulations/mpi_check.py
  echo "=== PETSc optional diagnostics: no implicit -info under two MPI ranks ==="
  # Run an actual .py script, not python -c: PETSc sees Python -c as an
  # unused command-line option and emits a false-positive warning at exit.
  # This script is created only inside the disposable USER workspace and
  # mounted at RUNTIME; it never enters the Docker build context/image.
  cat > "$USER_TEST_ROOT/project/petsc_options_check.py" <<'PY'
import scrimp  # Exercise the real package import and linalg initialization.
from petsc4py import PETSc

comm = PETSc.COMM_WORLD
assert comm.getSize() == 2, f"Expected two MPI ranks, got {comm.getSize()}"
options = PETSc.Options()
assert not options.hasName("info"), (
    "SCRIMP must not implicitly enable PETSc -info", options.getAll()
)
assert not options.hasName("c"), (
    "PETSc received a stray Python -c option", options.getAll()
)
print(f"PASS: PETSc import options clean (rank {comm.getRank() + 1}/2)",
      flush=True)
comm.barrier()
PY
  PETSC_OPTIONS_LOG="$ROOT/sandbox/outputs/v130-petsc-options-smoke.log"
  "${DOCKER[@]}" run --rm --user "$(id -u):$(id -g)" \
    --env HOME=/tmp --env XDG_CACHE_HOME=/tmp/.cache \
    --volume "$USER_TEST_ROOT/project:/work:rw" --workdir /work \
    --entrypoint /opt/conda/envs/scrimp/bin/mpiexec \
    "$IMAGE" --oversubscribe -n 2 \
      /opt/conda/envs/scrimp/bin/python petsc_options_check.py \
    2>&1 | tee "$PETSC_OPTIONS_LOG"
  if grep -Eq "There are options you set that were not used|Option left:" \
      "$PETSC_OPTIONS_LOG"; then
    echo "ERROR: unused PETSc option warning in import-only MPI run." >&2
    exit 2
  fi
  echo "PASS: PETSc -info and -c absent, no unused-option warning."
  rm -rf -- "$USER_TEST_ROOT"
  trap - EXIT

  echo "=== Persist explicit Conda spec and approved GetFEM patch checksums ==="
  "${DOCKER[@]}" run --rm "$IMAGE" \
    cat docker-resolved-spec.txt > "$ROOT/sandbox/outputs/v130-docker-resolved-spec.txt"
  "${DOCKER[@]}" run --rm "$IMAGE" \
    sha256sum scripts/patches/getfem-5.5/*.patch \
    > "$ROOT/sandbox/outputs/v130-getfem-patch-sha256.txt"
  echo "PASS: native, MPI and independent user-workspace host-output smoke."
} 2>&1 | tee "$LOG"
