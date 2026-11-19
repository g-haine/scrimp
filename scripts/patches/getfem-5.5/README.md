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
| `0004-nedelec-second-kind-simplices.patch` | Provide `FEM_NEDELEC2(N,K)` on triangles/tetrahedra; **explicitly only k=1,2** | **PASS** — constructors, shared DoFs, SPD mass, tangential traces; k>=3 unsupported |

All five approved top-level patches (`0000`–`0004`) are applied in
order by the next GetFEM 5.5 build. **The newly approved Nédélec-II
patch supports only k=1,2**; this release does not claim arbitrary
order polynomial families. For k>=3, the Nédélec-II DoFs need a
separate higher-order construction and validation.

## Lean native evidence

The BDM Docker smoke passed k=1/2 in 2D/3D, with 21 and 54 global DoFs
on two tetrahedra, SPD mass matrices and normal-trace differences
below 4.529e-6 relative for fields evaluated 1e-7 to either side.

The Nédélec-I Docker smoke passed k=1/2 on two triangles and two
tetrahedra, with the expected local and global DoF counts, SPD mass
matrices and maximum tangential mismatch **2.385e-6 relative** for
fields evaluated 1e-7 to either side.

The Nédélec-II Docker smoke passed the same targeted checks for
**k=1,2 on triangles and tetrahedra**: local DoFs 6/12 (2D),
12/30 (3D), global two-cell DoFs 10/21 (2D), 18/48 (3D),
SPD mass matrices, and maximum relative tangential mismatch
**1.446e-6** at offsets ±1e-7. All four cases passed. The GetFEM build reports
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
bash scripts/test-getfem55-nedelec2-patch-in-docker.sh  # k=1,2
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

**Nédélec-II is production-ready only for the lean k=1,2 SCRIMP
scope**, not for arbitrary k. Its reference scheme uses `[P_k]^d`:
k+1 tangential point samples on every edge, plus, at k=2, three
radial tangential point samples on every triangular face (interior
points in 2D, linkable face points in 3D). The face
direction is chosen by the unique dominant barycentric vertex, so
it does not depend on local simplex vertex numbering. The resulting
four reference interpolation matrices have exact full rank:
6/6 and 12/12 on triangles, 12/12 and 30/30 on tetrahedra
(`scripts/check_nedelec_simplicial_basis.py`).
These **nodal functionals** are an alternative to the canonical
integral moments used in
[DefElement's Nédélec-II definition](https://defelement.org/elements/nedelec2.html).
The native two-cell tangential-continuity and assembly smoke has
**passed for k=1,2**. Extending to k>=3 requires additional
face/interior functionals and remains an open enhancement.

## GetFEM source provenance

The approved patches target official GetFEM v5.5 (release tag or
dereferenced commit
`c242f7fe0989e9a90846c63eaccc6b495d774724`).
The existing point-0 qualified Docker image is untouched. With
k=1,2 Nédélec-II validation complete, a **new separately named
Docker image** must be built and validated before replacing the
qualified original.

The previous `scripts/patch-getfem-mpi-timing.py` rewriter was
replaced by the standard `0000-mpi-assembly-timing.patch`; no
special-case GetFEM source-rewriting scripts remain in production.

## Exact order support and evidence

| Family | In-source order range | Orders covered by SCRIMP native tests |
| --- | --- | --- |
| RT | GetFEM original supports parameterized k>=0 | No additional targeted tests in this workstream |
| BDM | GetFEM original family with tetrahedral patch, k>=1 | k=1,2 |
| Nédélec I | Parameterized k>=0 (k=0 aliases built-in) | Constructors k=0,1,2; two-cell tests k=1,2 |
| Nédélec II | **k=1,2 only** (explicit constructor restriction) | k=1,2 |

**Do not advertise the whole family for arbitrary polynomial order.**
In particular, k>=3 is not implemented for Nédélec II. A true
arbitrary-order second-kind implementation needs polynomial moments
or another complete set of face and interior DoFs with orientation
and unisolvency; simply lifting the restriction is invalid.

The precise **deferred higher-order implementation checklist**,
including per-entity DoF formulas, a concrete k=3 example and
small acceptance criteria for Nédélec I/II, BDM and RT, is in
[the higher-order GetFEM follow-up](../../docs/V130_GETFEM_HIGHER_ORDER_TODO.md).
