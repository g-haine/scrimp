#!/usr/bin/env bash
# Test the staged GetFEM BDM tetrahedron patch in a DISPOSABLE container.
# This recompiles only as needed; it never changes scrimp-mpi:local or host conda.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${SCRIMP_GETFEM_TEST_IMAGE:-scrimp-mpi:local}"
PATCHES="$REPO_ROOT/scripts/patches/getfem-5.5"

echo "Image: $IMAGE (ephemeral run; existing image and point-0 results untouched)"
echo "Patches: 0001 BDM tetrahedra + 0002 Python extension in-tree archive linkage"
echo "Applier: Python standard library (no GNU patch installation)"

mkdir -p "$REPO_ROOT/sandbox/outputs"
LOG="$REPO_ROOT/sandbox/outputs/getfem55_bdm_patch_smoke.log"

sudo docker run --rm \
  -v "$PATCHES:/tmp/scrimp-getfem-patches:ro" \
  -v "$REPO_ROOT/scripts/apply_unified_patch.py:/tmp/scrimp-apply-unified-patch.py:ro" \
  "$IMAGE" bash -lc '
    set -euo pipefail
    src=/workspace/scrimp/.build/getfem/getfem-5.5

    if [[ ! -f "$src/src/getfem_fem.cc" ]]; then
      echo "ERROR: GetFEM source tree not present at $src in this image" >&2
      exit 2
    fi
    # The qualified image lacks GNU patch. Verify the exact upstream
    # contexts before applying both source fixes in this disposable container.
    for patchfile in \
      /tmp/scrimp-getfem-patches/0001-bdm-tetrahedron-unisolvent-dofs.patch \
      /tmp/scrimp-getfem-patches/0002-python-interface-link-in-tree.patch; do
      python /tmp/scrimp-apply-unified-patch.py "$src" "$patchfile" --dry-run
      python /tmp/scrimp-apply-unified-patch.py "$src" "$patchfile"
    done

    # setup.py is generated from setup.py.in by configure. Refresh it,
    # otherwise make keeps linking GetFEM from the old linker arguments.
    (cd "$src" && ./config.status interface/src/python/setup.py)
    grep -nF "extra_objects = [" "$src/interface/src/python/setup.py"
    make -C "$src" -j2
    make -C "$src" install "am__py_compile=:"
    echo "GetFEM source: post-patch BDM assertion location:"
    grep -nF "GMM_ASSERT2(ipoint == nddl" "$src/src/getfem_fem.cc"
    echo "GetFEM source object / installed extension provenance:"
    ls -l "$src/src/getfem_fem.o" "$src/src/.libs/libgetfem.a"
    python - <<'"'"'PY'"'"'
import getfem
import importlib
import pathlib
import subprocess
import sys
extension = importlib.import_module("getfem._getfem")
print("GetFEM package:", getfem.__file__, flush=True)
print("Loaded GetFEM native extension:", extension.__file__, flush=True)
print("Extension mtime:", pathlib.Path(extension.__file__).stat().st_mtime, flush=True)
marker = b"SCRIMP_BDM3D_PATCH_0001"
src = pathlib.Path("/workspace/scrimp/.build/getfem/getfem-5.5")
artifacts = [
    src / "src/getfem_fem.o",
    src / "src/.libs/libgetfem.a",
    pathlib.Path(extension.__file__),
]
for artifact in artifacts:
    blob = artifact.read_bytes()
    count = blob.count(marker)
    print(f"PATCH MARKER: {artifact}: count={count}", flush=True)
    if count == 0:
        raise RuntimeError(
            f"Stale GetFEM binary: SCRIMP patch marker absent from {artifact}. "
            "The freshly compiled source has not reached the loaded extension."
        )
ldd = subprocess.run(["ldd", extension.__file__], check=True, capture_output=True, text=True)
loaded = [row.strip() for row in ldd.stdout.splitlines() if "libgetfem" in row]
print("External shared GetFEM dependencies:", loaded or "(none: linked statically)", flush=True)
print("GetFEM MPI level:", getfem.util_mpi_parallelism_level(), flush=True)
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
print("PASS: patched BDM local constructors/dimensions.")
print("Checking two tetrahedra with a common triangular face:", flush=True)
import numpy as np

mesh = getfem.Mesh("empty", 3)
geo = getfem.GeoTrans("GT_PK(3,1)")
# Same shared triangle, opposite vertices, and reversed local face ordering.
tetra1 = np.array([[0., 1., 0., 0.],
                   [0., 0., 1., 0.],
                   [0., 0., 0., 1.]])
tetra2 = np.array([[0., 0., 1., 0.],
                   [0., 1., 0., 0.],
                   [0., 0., 0., -1.]])
mesh.add_convex(geo, tetra1)
mesh.add_convex(geo, tetra2)
mim = getfem.MeshIm(mesh, 6)

for k, local, face in ((1, 12, 3), (2, 30, 6)):
    mf = getfem.MeshFem(mesh, 1)
    mf.set_fem(getfem.Fem(f"FEM_BDMK(3,{k})"))
    expected = 2 * local - face
    n = int(mf.nbdof())
    print(f"BDM{k}, two tetrahedra: global DoFs={n}, expected={expected}",
          flush=True)
    assert n == expected, (
        f"Face DoFs are not shared for BDM{k}: got {n}, expected {expected}"
    )
    mat = getfem.asm("mass matrix", mim, mf)
    dense = mat.toarray() if hasattr(mat, "toarray") else np.asarray(mat)
    assert dense.shape == (expected, expected), dense.shape
    assert np.all(np.isfinite(dense)), f"BDM{k}: nonfinite mass matrix"
    sym = float(np.linalg.norm(dense - dense.T) /
                max(1., np.linalg.norm(dense)))
    min_eig = float(np.linalg.eigvalsh((dense + dense.T)/2.).min())
    print(f"BDM{k}: mass symmetry error={sym:.3e}, "
          f"lambda_min={min_eig:.3e}", flush=True)
    assert sym < 1e-10 and min_eig > 0, (
        f"BDM{k}: mass matrix not symmetric positive definite"
    )
print("PASS: BDM3D shared-face DoF count and positive mass assembly.")
print("NOTE: full pointwise normal-flux continuity is not yet checked.")
PY
  ' 2>&1 | tee "$LOG"
