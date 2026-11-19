#!/usr/bin/env bash
# Probe the approved simplicial Nedelec-I patch inside a disposable
# GetFEM 5.5 container. Does NOT mutate scrimp-mpi:local or host conda.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${SCRIMP_GETFEM_TEST_IMAGE:-scrimp-mpi:local}"
PATCHES="$ROOT/scripts/patches/getfem-5.5"
LOG="$ROOT/sandbox/outputs/getfem55_nedelec1_patch_smoke.log"
mkdir -p "$ROOT/sandbox/outputs"

echo "Nedelec-I experimental Docker smoke (triangles/tetrahedra, k<=2)"
echo "Image: $IMAGE (disposable container, no Docker image rebuild)"
echo "Patches: 0001 BDM, 0002 Python link, 0003 Nedelec-I (now approved)"
echo "Output: $LOG"

sudo docker run --rm \
  -v "$PATCHES:/tmp/scrimp-patches:ro" \
  -v "$ROOT/scripts/apply_unified_patch.py:/tmp/scrimp-apply-patch.py:ro" \
  -v "$ROOT/sandbox/audit_getfem55_nedelec_simplices.py:/tmp/scrimp-nedelec-smoke.py:ro" \
  "$IMAGE" bash -lc '
    set -euo pipefail
    src=/workspace/scrimp/.build/getfem/getfem-5.5
    test -f "$src/src/getfem_fem.cc"
    for patchfile in \
      /tmp/scrimp-patches/0001-bdm-tetrahedron-unisolvent-dofs.patch \
      /tmp/scrimp-patches/0002-python-interface-link-in-tree.patch \
      /tmp/scrimp-patches/0003-nedelec-first-kind-simplices.patch; do
      python /tmp/scrimp-apply-patch.py "$src" "$patchfile" --dry-run
      python /tmp/scrimp-apply-patch.py "$src" "$patchfile"
    done
    (cd "$src" && ./config.status interface/src/python/setup.py)
    make -C "$src" -j2
    make -C "$src" install "am__py_compile=:"
    python -u /tmp/scrimp-nedelec-smoke.py
  ' 2>&1 | tee "$LOG"
