#!/usr/bin/env bash
# Test the staged GetFEM BDM tetrahedron patch in a DISPOSABLE container.
# This recompiles only as needed; it never changes scrimp-mpi:local or host conda.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${SCRIMP_GETFEM_TEST_IMAGE:-scrimp-mpi:local}"
PATCHES="$REPO_ROOT/scripts/patches/getfem-5.5"

echo "Image: $IMAGE (ephemeral run; existing image and point-0 results untouched)"
echo "Patch: $PATCHES/0001-bdm-tetrahedron-unisolvent-dofs.patch"

mkdir -p "$REPO_ROOT/sandbox/outputs"
LOG="$REPO_ROOT/sandbox/outputs/getfem55_bdm_patch_smoke.log"

sudo docker run --rm \
  -v "$PATCHES:/tmp/scrimp-getfem-patches:ro" \
  "$IMAGE" bash -lc '
    set -euo pipefail
    src=/workspace/scrimp/.build/getfem/getfem-5.5
    patchfile=/tmp/scrimp-getfem-patches/0001-bdm-tetrahedron-unisolvent-dofs.patch

    if [[ ! -f "$src/src/getfem_fem.cc" ]]; then
      echo "ERROR: GetFEM source tree not present at $src in this image" >&2
      exit 2
    fi
    if ! command -v patch >/dev/null 2>&1; then
      echo "ERROR: GNU patch is not installed in this preexisting image." >&2
      echo "The next Docker build includes the required dependency; do not alter the qualified image." >&2
      exit 2
    fi
    patch --batch --forward --dry-run -p1 -d "$src" < "$patchfile"
    patch --batch --forward -p1 -d "$src" < "$patchfile"
    make -C "$src" -j2
    make -C "$src" install "am__py_compile=:"
    python - <<'"'"'PY'"'"'
import getfem
print("GetFEM imported from:", getfem.__file__)
print("GetFEM MPI level:", getfem.util_mpi_parallelism_level())
assert int(getfem.util_mpi_parallelism_level()) == 2
tests = [
    ("FEM_BDMK(2,1)", 6),
    ("FEM_BDMK(2,2)", 12),
    ("FEM_BDMK(3,1)", 12),
    ("FEM_BDMK(3,2)", 30),
]
for descriptor, expected in tests:
    f = getfem.Fem(descriptor)
    actual = int(f.nbdof())
    print(f"{descriptor}: local DoFs={actual} (expected {expected})")
    assert actual == expected, descriptor
print("PASS: patched BDM local constructors/dimensions only.")
print("Still required: 3D shared-facet H(div) check and small assembly.")
PY  ' 2>&1 | tee "$LOG"
