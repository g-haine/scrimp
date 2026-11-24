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

## Source-only ownership analysis (stage 2, no assembly change)

`scrimp/interface_ownership.py` now contains a conservative
`infer_interface_port_owners(...)` helper. For each logical interface,
it associates a port with a unique owning subdomain only if interface
weak forms explicitly reference a variable whose physical bulk region
is known. Interface-only masses and cross-feedback equations are *not*
used as ownership anchors: feedback such as
`Y_w*Test_U_T` intentionally crosses subsystem boundaries.

It is a **standalone audit only**: its string tokenization is NOT a full
GWFL syntax parser. In particular, nontrivial user macros, aliases,
multi-region forms and implicit operators must be handled or explicitly
rejected before the production compiler can rely on it. Source-only
unit tests verify positive Heat-Wave-style ownership and negative
ambiguous/contradictory cases:

```bash
python tests/test_v130_interface_ownership.py
```

No GetFEM region, FEM variable, brick, physical sign or power-report
policy is changed by this addition.

## DPHS declaration-level audit (stage 3, native check pending)

`DPHS.audit_interface_ownership(interface, (a,b), mesh_id=0)` is a
new **explicit, read-only** model inspection method. It combines the
native geometry of `Domain.resolve_interface` with the actual declared
SCRIMP `Port` and `Brick` objects. It derives bulk-variable ownership
from port regions and derives interface-port ownership only from weak
forms anchored by those bulk variables. It returns the owner for each
interface port and the verified geometric facet count.

This method is NOT invoked by `DPHS.solve`, `add_FEM`, `add_brick`
or `set_control`. It does not allocate GetFEM regions, change the weak
forms, alter control expressions, affect Jacobians or DAE initialization,
or imply the energy law is correct. It is deliberately not an implicit
feature flag or production orientation policy.

The native source-immutable audit of the existing Heat-Wave example is:

```bash
python sandbox/probe_v130_native_interface_ownership.py
```

The probe intercepts the final `DPHS.solve` call, inspects the complete
set of already declared ports/bricks, and exits before the TS solver.
For `heat_region=1, wave_region=2`, it must report
`"Interface Heat": 1` and `"Interface Wave": 2` and 78 matched faces.
No actual simulation is claimed.

## Geometric regression: two adjacent rectangles (native pending)

A two-sided interface is generally **not the full boundary** of either
subdomain. For two adjacent rectangles, the shared vertical segment is a
proper subset of BOTH rectangle perimeters. The topology resolver
constructs both full subdomain boundary maps, intersects their facet
keys, then checks the logical interface exactly matches that shared set.
Nonshared exterior facets are deliberately excluded. If the region
label fails to identify some truly shared facets, the resolver rejects
the mismatch rather than silently discarding them.

```bash
python sandbox/probe_v130_two_rectangles_interface.py
```

This independent **native** GetFEM probe constructs 8 Cartesian
quadrangles across two adjacent rectangles, identifies the two shared
interface facets, verifies that other exterior facets are untouched,
checks subdomain-permutation invariance, and rejects a deliberately
incomplete logical interface. It does not invoke a PETSc time solver,
mass-matrix projection or reorientation of existing DPHS weak forms.
The native probe must be executed before claiming validation.

## Native declared-model audit: projection avoidance (pending)

`sandbox/probe_v130_native_interface_ownership.py` now intercepts
the first `DPHS.set_initial_value` call, after all original Heat-Wave
ports, bricks and controls have been declared but **before** any
L2 initial-value projection. The preceding attempt intercepted
`DPHS.solve` too late and failed on a local petsc4py installation
without `Mat.setPreallocationCOO`. That failure is not evidence of
wrong interface ownership; nor is it solved by this audit-only
interception. Production PETSc COO compatibility remains a separate
environment qualification.

## Next steps (not implemented here)

1. Native test on real 1D/2D/3D imported meshes and curved/high-order
   elements; check different Gmsh physical group representations.
2. Extend the ownership audit to the actual declared DPHS model, a full
   GWFL expression/dependency resolver, and explicit rejection or handling
   of ambiguous ports and cross-subdomain feedback.
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
