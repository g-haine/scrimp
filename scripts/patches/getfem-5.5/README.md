# SCRIMP-maintained GetFEM 5.5 source patches

**Scope:** triangles and tetrahedra only. This is a lean,
SCRIMP-first implementation intended to help users run structure-preserving
finite-element models. The patches can subsequently serve as a starting
point for the GetFEM community; *upstream acceptance and exhaustive
certification are not SCRIMP v1.3.0 release requirements.*

`scripts/build-getfem-mpi.sh` applies the top-level `*.patch` files in
lexicographic order to the official GetFEM **v5.5** source tree before
compilation. `docker/environment.yml` supplies GNU `patch` for the
next production image build. It fails on incompatible GetFEM sources.

| Patch | Purpose | Status (19 Nov 2026) |
| --- | --- | --- |
| `0000-mpi-assembly-timing.patch` | Make MPI assembly-time prints opt-in with `SCRIMP_GETFEM_ASSEMBLY_TIMING=1` | Source diff checked; next Docker rebuild pending |
| `0001-bdm-tetrahedron-unisolvent-dofs.patch` | Fix tetrahedral BDM k=1,2 face DoFs and interior basis completion | **PASS** — native constructors, shared-face DoFs, SPD mass, normal traces |
| `0002-python-interface-link-in-tree.patch` | Link Python bindings against newly built GetFEM static archives explicitly | **PASS** — patched native marker present in loaded `_getfem` extension |
| `0003-nedelec-first-kind-simplices.patch` | Provide `FEM_NEDELEC1(N,K)`, triangles/tetrahedra, k=0 delegates to built-in Nédélec | **PASS** — native k=1/2 constructors, shared-edge/face DoFs, SPD mass, tangential traces |
| `experimental/0004-nedelec-second-kind-simplices.patch` | Provide `FEM_NEDELEC2(N,K)` for triangles/tetrahedra, currently k=1,2 | **EXPERIMENTAL: exact reference rank PASS; native smoke NOT RUN** |

The standard build loop deliberately ignores nested
`experimental/*.patch` files. Nédélec-II will only be promoted to
the top level after it compiles and passes the same small two-cell
smoke as Nédélec-I.

## Lean native evidence

The BDM Docker smoke passed k=1/2 in 2D/3D, with 21 and 54 global DoFs
on two tetrahedra, SPD mass matrices and normal-trace differences
below 4.529e-6 relative for fields evaluated 1e-7 to either side.

The Nédélec-I Docker smoke passed k=1/2 on two triangles and two
tetrahedra, with the expected local and global DoF counts, SPD mass
matrices and maximum tangential mismatch **2.385e-6 relative** for
fields evaluated 1e-7 to either side. The GetFEM build reports
`util_mpi_parallelism_level() == 2`, but these are sequential tests
inside an MPI-enabled image, **not MPI rank-parity tests**.

These are sufficient for SCRIMP's first release-oriented gate, not
proofs for all curved cells or higher degrees. More involved
Maxwell energy/power and 1/2-rank tests belong to the Maxwell
workstream, not to the GetFEM contribution base.

## Reproducible disposable Docker smoke commands

From the SCRIMP repository root:

```bash
bash scripts/test-getfem55-bdm-patch-in-docker.sh
bash scripts/test-getfem55-nedelec1-patch-in-docker.sh
bash scripts/test-getfem55-nedelec2-patch-in-docker.sh  # experimental
```

They run **inside disposable containers** created from the
qualified `scrimp-mpi:local` image. They do not change that image or
the host Conda environment and require no complete Docker rebuild.
The image lacks GNU `patch`, so the disposable probes instead use
`scripts/apply_unified_patch.py`, a strict Python-standard-library
unified-diff helper. The probes regenerate GetFEM's configured
`interface/src/python/setup.py` after applying the Python-link patch,
then rebuild the core and binding inside the container.

Logs remain in `sandbox/outputs/`:

- `getfem55_bdm_patch_smoke.log`
- `getfem55_nedelec1_patch_smoke.log`
- `getfem55_nedelec2_patch_smoke.log`

**Nédélec-II is not yet a production element.** Its reference scheme
uses `[P_k]^d`: k+1 tangential point samples on every edge, plus,
at k=2, three radial tangential point samples on every triangular
face (interior points in 2D, linkable face points in 3D). The face
direction is chosen by the unique dominant barycentric vertex, so
it does not depend on local simplex vertex numbering. The resulting
four reference interpolation matrices have exact full rank:
6/6 and 12/12 on triangles, 12/12 and 30/30 on tetrahedra
(`scripts/check_nedelec_simplicial_basis.py`).
These **nodal functionals** are an alternative to the canonical
integral moments used in
[DefElement's Nédélec-II definition](https://defelement.org/elements/nedelec2.html).
The native two-cell tangential-continuity and assembly probe will
determine whether the GetFEM implementation is appropriate for SCRIMP.

## GetFEM source provenance

The approved patches target official GetFEM v5.5 (release tag or
dereferenced commit
`c242f7fe0989e9a90846c63eaccc6b495d774724`).
The existing point-0 qualified Docker image is untouched. After
Nédélec-II validation, a **new separately named Docker image** must
be built and validated before replacing it.

The previous `scripts/patch-getfem-mpi-timing.py` rewriter was
replaced by the standard `0000-mpi-assembly-timing.patch`; no
special-case GetFEM source-rewriting scripts remain in production.
