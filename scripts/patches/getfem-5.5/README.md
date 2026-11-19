# SCRIMP-specific GetFEM 5.5 source patches

These are *small source patches* applied to the **official GetFEM 5.5 release**
by `scripts/build-getfem-mpi.sh` before the native GetFEM build. The goal is
a useful SCRIMP build first; offering these changes to GetFEM upstream may
happen later. They are not presented as an upstream-certified FEM library.

**Only triangles and tetrahedra are in scope.**

| Patch | Why | State |
|---|---|---|
| `0001-bdm-tetrahedron-unisolvent-dofs.patch` | Fixes BDM facet normal DoF lattice for tetrahedra; uses interior nullspace completion from k>=2, while preserving existing 2D behavior | **Staged, NOT compiled or qualified yet** |
| `experimental/0002-nedelec-first-kind-simplices.patch` | Prototype `FEM_NEDELEC1(N,K)` on triangles/tetrahedra, reusing `FEM_NEDELEC(N)` at k=0 | **EXPERIMENTAL, not applied by default, no native build yet** |

The original Docker GetFEM 5.5 MPI-2 image produced `FEM_BDMK(3,1)`
and `FEM_BDMK(3,2)` assertion failures at `getfem_fem.cc:2319`.
The staged patch fixes the identified *construction logic*, but no
claims of working H(div) conformity or high-order stability are made
until the build and lightweight smoke tests pass.

### Build wiring

The build script uses GNU `patch` (provided by `docker/environment.yml`)
and fails if a patch no longer applies cleanly. Source patches are only
applied for `GETFEM_VERSION=5.5`; git-source builds must use the
`v5.5` tag or the dereferenced release commit
`c242f7fe0989e9a90846c63eaccc6b495d774724`.

The existing, qualified `scrimp-mpi:local` image is **not modified**.
Build and test a *separately named* image before replacing it during
SCRIMP v1.3.0 Docker workstream.

### Test in a disposable container (without touching the qualified image)

From the SCRIMP repository root:

```bash
bash scripts/test-getfem55-bdm-patch-in-docker.sh
```

This runs on a *temporary* container created from `scrimp-mpi:local`,
patches its source, performs an incremental GetFEM build/install, and
checks BDM1/2 constructor dimensions in 2D/3D. The compile/test output
is saved by the host script to
`sandbox/outputs/getfem55_bdm_patch_smoke.log`.

The preexisting image needs `patch` installed for this convenience
smoke; the next Docker build includes it explicitly. The script fails
cleanly if the legacy image lacks the tool or GetFEM source tree.
Passing this local smoke does **not** constitute conforming H(div) or
MPI qualification.

### Minimal checks (before enabling for SCRIMP users)

1. Constructor and expected local dimensions: BDM 2D k=1,2 and 3D k=1,2.
2. A second, higher-order constructor smoke if LAPACK is detected.
3. One two-tetrahedron shared normal-trace check and a small mass/div
   operator assembly. Retain SCRIMP point-0 baseline.

### Nédélec

`scripts/check_nedelec_simplicial_basis.py` provides an independent,
parameterized first-kind reference-space and unisolvency prototype.
It is *not yet* a GetFEM patch. The higher-order first-kind and
second-kind source implementations will be added here only once their
local basis and shared-edge/shared-face DoFs are coherent.

The Nédélec 1st-kind prototype resides in `experimental/` intentionally:
the build script applies only `scripts/patches/getfem-5.5/*.patch`, not
nested experimental patches. Its basis interpolation matrix was checked
in exact arithmetic for 2D/3D at k=0..2 by
`scripts/check_nedelec_simplicial_basis.py`, but C++ compilation,
global DoF linking and 3D physical-face orientation are **not yet tested**.
Do not enable it for end users before these checks. The parameterized
2nd-kind construction is still to be developed.
