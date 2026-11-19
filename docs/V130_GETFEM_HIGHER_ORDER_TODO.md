# GetFEM 5.5 simplicial FEM — higher-order follow-up (deferred)

**Status (19 November 2026): deferred, not a SCRIMP v1.3.0 release gate.**
The supported native smoke scope is **triangles and tetrahedra,
Nédélec I/II and BDM at k=1,2**; Nédélec-I k=0 uses GetFEM's built-in
lowest-order element. Existing Raviart–Thomas (RT) implementation is
reused; this workstream did not independently requalify every RT order.
The complete evidence lives in
[V130_GETFEM_FEM_AUDIT.md](V130_GETFEM_FEM_AUDIT.md).

This note identifies exactly what remains **before claiming arbitrary
polynomial order**, without expanding the current release workload.
The base patches under `scripts/patches/getfem-5.5/` remain
SCRIMP-first contributions rather than GetFEM-upstream-ready packages.

## 1. Primary implementation gap: Nédélec second kind for k >= 3

**Code:** `scripts/patches/getfem-5.5/0004-nedelec-second-kind-simplices.patch`.

The current `FEM_NEDELEC2(N,K)` constructor **explicitly rejects
k >= 3**. It constructs the full vector polynomial space
`[P_k]^d` but its **pointwise tangential interpolation functionals**
cover only k=1 and k=2: k+1 samples on each edge, plus three
radial-tangent samples per triangle interior at k=2. The latter are
face-interior, globally linkable samples when the triangle is a
tetrahedral face; in 2D they are cell-interior functionals.

For arbitrary k, choose a complete set of functionals based on
integral tangential moments (the canonical family) **or** a
provably unisolvent, geometrically invariant nodal construction.
The present `k+2` barycentric lattice + one radial direction at
each face-interior point is **not** sufficient for k >= 3; merely
removing the constructor guard is incorrect.

Required per-entity DoF counts for `[P_k]^d` (k >= 1):

| Entity / count | Triangle (2D) | Tetrahedron (3D) |
| --- | ---: | ---: |
| Local space dimension | `(k+1)(k+2)` | `(k+1)(k+2)(k+3)/2` |
| Each edge | `k+1` tangential | `k+1` tangential |
| Each triangular face interior | n/a (cell interior) | `(k+1)(k-1)` |
| Cell interior | `(k+1)(k-1)` | `(k+1)(k-1)(k-2)/2` |
| Vertices | 0 | 0 |

For tetrahedra, multiply the edge row by 6 and the face-interior
row by 4 to recover the total number of boundary functionals.
The remaining volume-interior count is the last row.
The listed counts are algebraic **requirements**, not claims that
corresponding functionals already exist in the patch.

**Concrete first missing order, k=3:**

- Triangle: 20 local DoFs = 12 edge + **8 cell-interior**.
- Tetrahedron: 60 local DoFs = 24 edge + **32 face-interior
  (8 on each face)** + **4 volume-interior**.

The first higher-order implementation must therefore:

1. Define independent polynomial moments or nodal tangential
   functionals with **exactly** these counts on each entity; in
   particular generalize 2D cell-interior and 3D face/volume-interior
   functionals simultaneously.
2. Preserve **edge/face orientation independence** when two simplices
   order their vertices differently, including geometric
   covariant (H(curl)) transformation through GetFEM's `mat_trans`
   and correct global DoF identification.
3. Establish local unisolvence (exact symbolic rank initially for
   k=3, then conditioning of interpolation matrices numerically);
   avoid assuming that degree-k monomial sets and arbitrary point
   samples remain stable at higher degree.
4. Ensure the implementation and polynomial-degree convention agree
   with the public `FEM_NEDELEC2(N,K)` descriptor, and document its
   relation to conventional Nédélec-II moment definitions.

## 2. Other families: what is implemented versus what is tested

**Nédélec first kind**, `FEM_NEDELEC1(N,K)`,
`0003-nedelec-first-kind-simplices.patch`:

- The constructor is **parameterized** for k >= 0; k=0 delegates
  to `FEM_NEDELEC(N)`. Its local dimension formulas are
  `(k+1)(k+3)` (triangle) and
  `(k+1)(k+3)(k+4)/2` (tetrahedron).
- Native GetFEM tests currently cover **k=0,1,2 constructors**
  and **k=1,2** two-cell linkage, mass and tangential traces.
- Before claiming k>=3: test the local interpolation-matrix rank,
  numerical conditioning and the intended polynomial span; then
  test opposite edge/face orientations, tangential continuity and
  representative mass/curl assembly at k=3. Investigate any
  construction, geometry or conditioning failures rather than
  treating the parameterized constructor as proof of support.

**BDM**, `FEM_BDMK(N,K)`,
`0001-bdm-tetrahedron-unisolvent-dofs.patch`:

- The tetrahedral patch fixes **k=1,2** in the native Docker smoke.
  Existing higher orders are not newly qualified.
- For k>=3, validate the count and rank of normal-face functionals,
  completion of the interior nullspace (including the SVD threshold
  used by the patch), numerical conditioning, opposite-face
  orientation/normal-trace matching and small mass/divergence
  assembly. Retain original working 2D behavior.

**RT**, `FEM_RTK(N,K)`:

- Reuse the upstream GetFEM 5.5 parameterized family; no new
  implementation is requested in the current scope.
- If all-order support is later advertised explicitly by SCRIMP,
  add at least a representative higher-order simplex constructor,
  shared-normal-trace and operator-assembly smoke at k=3. Keep
  this separate from any future upstream certification effort.

## 3. Lean future acceptance criteria (only when work resumes)

For each newly supported order, and in both 2D and 3D:

1. **Reference:** expected dimension and full-rank interpolation
   functionals; record condition number for growing k.
2. **Interface:** two neighboring simplices with reversed edge/face
   ordering; expected global DoFs and physical normal/tangential
   traces sampled on both sides.
3. **Assembly:** finite, symmetric positive-definite mass matrix
   and one appropriate curl/divergence operator; exercise the
   correct mapped basis rather than relying only on cardinality.
4. **Delivery:** compile the patch in a **separately named disposable
   or candidate Docker image**, compare the expected GetFEM symbol
   version/patch set, and add a short corresponding acceptance entry
   to the FEM audit. Preserve `scrimp-mpi:local` as qualified baseline
   until the new image passes.

Do **not** add a large GetFEM-upstream test suite to SCRIMP solely
to satisfy this deferred task. Maxwell's equations, power/energy
balances and MPI rank-parity are separate integration gates in
the [release plan](V130_RELEASE_PLAN.md), not substitutes for
the element's H(curl)/H(div) continuity checks.

## 4. Current v1.3.0 boundary

**Accepted now:** k=1,2 in the small native GetFEM simplicial
smokes for BDM and Nédélec I/II, using their current patched source;
the original RT remains available.

**Not claimed:** Nédélec-II k>=3; native validation for Nédélec-I
or BDM k>=3; arbitrary geometries or polynomial degrees; MPI
rank-parity for these elements; a complete Maxwell model or upstream
acceptance.

**Next release-plan action:** point 2, rebuild a reproducible,
cache-friendly SCRIMP v1.3.0 Docker candidate with the approved
`0000`–`0004` patches. Test it under a **new image tag**, then
address installer/launcher and cross-platform documentation.
No part of this higher-order backlog blocks that next action.
