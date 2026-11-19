# SCRIMP-specific GetFEM 5.5 source patches

These are *small source patches* applied to the **official GetFEM 5.5 release**
by `scripts/build-getfem-mpi.sh` before the native GetFEM build. The goal is
a useful SCRIMP build first; offering these changes to GetFEM upstream may
happen later. They are not presented as an upstream-certified FEM library.

**Only triangles and tetrahedra are in scope.**

| Patch | Why | State |
|---|---|---|
| `0000-mpi-assembly-timing.patch` | Makes GetFEM's unconditional MPI `Assembly time` print opt-in with `SCRIMP_GETFEM_ASSEMBLY_TIMING=1`, without changing assembly | **Included in the next Docker build; matches original GetFEM v5.5 source; native rebuild pending** |
| `0001-bdm-tetrahedron-unisolvent-dofs.patch` | Fixes BDM facet normal DoF lattice for tetrahedra; uses interior nullspace completion from k>=2, while preserving existing 2D behavior | **Native Docker lean acceptance PASS: BDM1/2 3D constructors, shared-face DOFs, SPD mass and normal-trace samples** |
| `0002-python-interface-link-in-tree.patch` | Links the Python `_getfem` extension explicitly against the **fresh in-tree** `libgetfemint.a` and `libgetfem.a` instead of resolving potentially stale `-lgetfem` | **Native Docker relink PASS; patch marker present in loaded extension** |
| `experimental/0003-nedelec-first-kind-simplices.patch` | Prototype `FEM_NEDELEC1(N,K)` on triangles/tetrahedra, reusing `FEM_NEDELEC(N)` at k=0 | **EXPERIMENTAL, not applied by default, no native build yet** |

The original Docker GetFEM 5.5 MPI-2 image produced `FEM_BDMK(3,1)`
and `FEM_BDMK(3,2)` assertion failures at `getfem_fem.cc:2319`.
The patch now compiles and constructs BDM1/2 on triangles and tetrahedra
with their expected local dimensions in Docker (19 November 2026).
The 3D lean smoke also passed two-tetrahedron shared-face DoF counts,
positive-definite mass assembly and physical normal-trace checks at
±1e-7 (largest relative mismatch 4.529e-6). This is enough for the
v1.3.0 SCRIMP-specific BDM1/2 scope; higher orders, arbitrary curved
meshes and upstream certification are **not** claimed.

### Build wiring

The build script applies all top-level `*.patch` files in filename order (`0000`, `0001`, ...), using GNU `patch` (provided by `docker/environment.yml`), and fails if a patch no longer applies cleanly. Source patches are only
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

The preexisting, qualified image **does not have GNU `patch`**.
The disposable smoke therefore uses
`scripts/apply_unified_patch.py` (strict source-context matching,
Python standard library only) to apply the unchanged, standard
`.patch` files **inside the disposable container**. It applies
`0001` (BDM) and `0002` (Python linkage), then regenerates the
configured `interface/src/python/setup.py` using `config.status`
before compiling the binding. The resulting extension contained the BDM patch marker in the latest
native constructor run (PASS 4/4); the same test also passed global
face DoFs, mass-matrix assembly and near-face normal-trace samples.
This requires no Docker rebuild, download, root-package install,
or modification of the qualified image. The later v1.3 Docker build
still uses normal GNU `patch` for all approved top-level patches.

The script fails cleanly if the expected GetFEM source tree or exact
patch context is missing. Passing this local constructor smoke does
**not** constitute conforming H(div) or MPI qualification.

### Minimal checks (before enabling for SCRIMP users)

1. Constructor and expected local dimensions: BDM 2D k=1,2 and 3D k=1,2.
2. A second, higher-order constructor smoke if LAPACK is detected.
3. One two-tetrahedron shared normal-trace check and a small mass/div
   operator assembly. Retain SCRIMP point-0 baseline.

### Nédélec

`scripts/check_nedelec_simplicial_basis.py` provides an independent,
exact-arithmetic first-kind reference-space and unisolvency check.
The **experimental first-kind GetFEM source patch exists** at
`experimental/0003-nedelec-first-kind-simplices.patch`; the
second-kind patch is **not yet implemented**.

The Nédélec 1st-kind prototype resides in `experimental/` intentionally:
the build script applies only `scripts/patches/getfem-5.5/*.patch`, not
nested experimental patches. Its basis interpolation matrix was checked
in exact arithmetic for 2D/3D at k=0..2 by
`scripts/check_nedelec_simplicial_basis.py`, but C++ compilation,
global DoF linking and 3D physical-face orientation are **not yet tested**.
Do not enable it for end users before these checks. The parameterized
2nd-kind construction is still to be developed.

### MPI assembly timing patch

`0000-mpi-assembly-timing.patch` replaces the former custom
`scripts/patch-getfem-mpi-timing.py` rewrite. The unified diff changes
only `src/getfem_models.cc`: it adds `#include <cstdlib>` and checks
`SCRIMP_GETFEM_ASSEMBLY_TIMING=1` before the master-rank
`Assembly time` message. Nothing in the numerical assembly changes.
The script was removed so *every production GetFEM modification* now
uses the same source patch mechanism. The qualified existing image is
not changed. To get the new behavior, use the planned Docker rebuild.

### Why the Python-link patch was needed

A native disposable-container run on 19 November 2026 confirmed that
GetFEM's patched `getfem_fem.o` and `src/.libs/libgetfem.a` contained
the unique BDM diagnostic marker, whereas the installed Python
`getfem/_getfem*.so` did **not**. The original setuptools command
used `-lgetfem`, which does not unambiguously identify the just-built
archive. `0002` replaces the implicit GetFEM archive lookup with
**explicit paths to the two source-tree archives**, retaining the
same other link dependencies. This is a linking hypothesis pending
native confirmation by the next disposable smoke; it does not alter
GetFEM finite-element mathematics.

The disposable BDM smoke has now **fully passed for our lean scope**:
local constructors, 21/54 global DoFs, positive-definite mass matrices,
and four near-face normal-trace checks (BDM1/BDM2, two fields each).

### Experimental Nédélec-I native smoke

The experimental `0003` file has been rebased so its strict unified
diff applies **after** the approved BDM `0001` changes to
`src/getfem_fem.cc`. It is intentionally still excluded from the
standard `scripts/build-getfem-mpi.sh` Docker build.

To compile and probe it in a disposable container:

```bash
bash scripts/test-getfem55-nedelec1-patch-in-docker.sh
```

The script applies `0001`, `0002`, then the experimental `0003`
using the no-install Python patch helper. It regenerates the GetFEM
Python build configuration, compiles and installs **inside the
disposable container**, and runs
`sandbox/audit_getfem55_nedelec_simplices.py`.

The small smoke checks reference constructors for k=0/1/2, and k=1/2
two-simplex global DoF linking, SPD mass matrices, and physical
tangential trace continuity sampled from both cells. It covers
**triangles and tetrahedra only**. Output is saved in
`sandbox/outputs/getfem55_nedelec1_patch_smoke.log`.

**Status: prepared, not yet compiled or run by the user.** Do not
activate `0003` automatically until these checks pass. Nédélec-II
construction will then be addressed separately.
