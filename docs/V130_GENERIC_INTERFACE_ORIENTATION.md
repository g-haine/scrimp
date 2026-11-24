# SCRIMP v1.3.0 — Generic conforming interface orientation (stage 1)

Status: **NATIVE GETFEM 2-D CONCENTRIC TOPOLOGY PASS (MPI1); NOT ACTIVE IN DPHS ASSEMBLY**.

## Contract

An interface region (for example `10`) is a *logical* physical boundary
shared by two adjacent bulk regions. The GetFEM mesh contains two outward
parent-convex facet representations. A variational form involving
`Normal` must use the side belonging to the appropriate physical
subdomain, never a union of both sides.

`scrimp.interface.resolve_conforming_interface(mesh, dim, interface,
(subdomain_a, subdomain_b))` resolves the two geometric sides. The
public `Domain.resolve_interface(interface, (a,b), mesh_id=0)` exposes
the same read-only information without changing the interface region or
allocating GetFEM regions.

This stage uses geometric-transform dimensions (rather than node counts),
boundary faces from each volumetric region, canonical facet point-ID
signatures, and a check that corresponding outward normals oppose. It
supports GetFEM logical interface regions stored as (parent convex,
local face) entries **or** imported whole (dim-1) convex entries.

**Fail closed**: a missing, extra, ambiguous, overlapping, or
non-conforming facet raises `InterfaceTopologyError`. This stage does
not implement nonmatching meshes, mortar coupling, periodic identification,
or curved high-order interface matching by coordinates. Exact point-ID
equality is required. It is intentionally a non-mutating geometric query;
it does not infer physical transmission laws or alter signs.

## Source-only checks

```bash
python tests/test_v130_interface_topology.py
python -m py_compile scrimp/interface.py tests/test_v130_interface_topology.py \
  sandbox/probe_v130_generic_interface_topology.py
```

Tests use a dependency-light simulated GetFEM mesh to exercise 1D, 2D,
3D, reordered subdomains, imported lower-dimensional facets, rejection of
inconsistent normals/topology, and defensive copies.

## Native qualification required

In the qualified GetFEM/PETSc environment:

```bash
python sandbox/probe_v130_generic_interface_topology.py
```

This checks both sides of the real `Concentric` 2D mesh, subdomain
permutation invariance, outward-normal cancellation, and preservation of
all original GetFEM regions. No example is modified.

**Do not interpret success as a Heat-Wave, DAE, CN, or MPI validation.**
The previous `sandbox/probe_v130_hw_dual_sided_cn.py` establishes
physical CN behavior only for its own reversible sandbox re-targeting,
which is not yet in core.

## Native evidence — 2026-11-24

User-run local GetFEM/PETSc environment on the 2-D Concentric mesh,
branch `release/1.3.0-preparation` at `4981299`:

- `python tests/test_v130_interface_topology.py`: **10 tests, OK**.
- `python -m py_compile`: passed.
- `python sandbox/probe_v130_generic_interface_topology.py`:
  `PASS_ORIENTED_TOPOLOGY_ONLY`.
- Logical interface region: **10**; subdomains: **1, 2**.
- Facets per side: **78**.
- Maximum opposite-normal *vector* defect:
  **2.482534153247273e-16**.
- Permutation invariant: **true**.
- Original logical interface region unchanged: **true**.
- No additional GetFEM regions allocated: **true**.

The genuine native run confirms the 2-D topology and the Python
`mesh.geotrans()` behavior (compressed GT catalog, singleton proxy).
These are **not** FE-assembly, DAE, discrete-energy, 3-D, or MPI2
qualifications. See the separate reversible Heat-Wave CN audit for
the numerical evidence on manually separated sides.

## Next steps (not implemented here)

1. Native test on real 1D/2D/3D imported meshes and curved/high-order
   elements; check different Gmsh physical group representations.
2. Define variable/port ownership of one-sided interface traces, including
   ambiguous ports and feedback across both subdomains.
3. Defer registration of affected GetFEM filtered variables and bricks
   until the side assignment is resolved. This requires a carefully
   staged DPHS compilation path; `add_FEM` and `add_brick` currently
   register immediately. Preserve legacy behavior outside conforming
   two-sided interfaces.
4. Bind automatic one-sided internal regions in the GetFEM model without
   duplicate quadrature/DOFs; validate sign conventions for all
   `Normal` terms and the port power accounting.
5. Validate original Heat-Wave both ways, actual assembled Jacobians,
   strictly initialized CN energy identities, physical nonzero exchange,
   MPI 1/2, and external boundary regressions. Only then remove the
   example's manual geometric sign adaptation.

**Do not silently auto-select an interface side from variable names,
brick names, region numbers, geometry radii, or the ordering of bulk
region IDs.** The resolved geometric normal and the physical interconnection
law are separate concepts.
