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

The probe intercepts the first `DPHS.set_initial_value` call, inspects
all declared ports/bricks, and exits before projection or the TS solver.
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

## Read-only weak-form placement plan (stage 4; native pending)

`DPHS.audit_interface_brick_plan(interface, (a, b), mesh_id=0)`
combines the already native-qualified topology and declared-port
ownership with a **non-mutating placement proposal** for each brick
on the logical interface. The decisive criterion for a residual-form
brick is its `Test_...` variable. In particular, the cross-subdomain
feedback `Y_w*Test_U_T` is placed on the **heat** side (the owner
of `U_T`), not the wave side (the owner of `Y_w`).

The diagnostic now verifies both region permutations and the ten
declared original Heat-Wave interface bricks before initial-value
projection. Every actual region 10 GetFEM brick still assembles as
before: **the proposal does not modify any FEM, face region, Jacobian,
port power, or weak form**. This is not yet the generic compiler.
Non-test interface forms, ambiguous multi-owner residuals, unknown
macro aliases and unsupported nonconforming couplings must be
rejected or explicitly handled before enabling automatic assembly.

```bash
python tests/test_v130_interface_assembly_plan.py
python sandbox/probe_v130_native_interface_ownership.py
```

The tests and the updated native audit have not yet been executed in
the qualified GetFEM environment.

## Native qualification completed — 2026-11-24

The user reran the stage-4 source-only suites and original Heat-Wave
declaration-level audit at commit `da9cb8f`:

- Interface topology: **10/10 PASS**.
- Port ownership: **6/6 PASS**.
- Residual-side brick plan: **7/7 PASS**.
- Native original Heat-Wave declarative audit:
  `PASS_DECLARED_DPHS_INTERFACE_OWNERSHIP_AND_BRICK_PLAN_ONLY`.
- Both `(heat,wave)=(1,2)` and `(2,1)` assign the
  ten interface bricks to the intended side, including both feedback
  equations and masses.
- No GetFEM weak-form reorientation, DAE initialization, CN, or MPI2
  execution occurred. The PETSc `-isset` unused-options warning
  is incidental to early diagnostic exit.

The user also ran
`sandbox/probe_v130_two_rectangles_interface.py` successfully
at commit `c2d64b4`:
`PASS_PARTIAL_BOUNDARY_INTERFACE_TOPOLOGY_ONLY`,
with 8 quadrangles (4 in each rectangle), 2 shared facets out of 8
boundary facets on each side, and rejection of incomplete labeling.

## Reference-facet parity gate (native pending)

`sandbox/probe_v130_generic_interface_reference_parity.py` compares
the generic `Domain.resolve_interface(10, (1,2))` face sets against the
previous physically qualified, geometry-specific Heat-Wave sandbox
helper. Exact parent convex and local face IDs must match, not merely
their coordinates or numbers. This pilot creates temporary legacy
reference regions 30/31 *only after* generic facet resolution; the
original logical region 10 stays unchanged.

Run from the qualified environment with the sandbox module on the
Python search path:

```bash
PYTHONPATH="$PWD/sandbox:$PWD" \\
  python sandbox/probe_v130_generic_interface_reference_parity.py
```

This is geometric parity, not yet assembly parity. The next gate must
construct the actual GetFEM variables and bricks from the automatically
resolved face sets and compare assembled sparse matrices with the
previous native qualified two-sided pilot. No production assembly policy
is activated yet.

## Full native assembled stiffness parity (stage 5, pending)

`scripts/run-v130-generic-interface-K-parity.sh` runs two instances
of the **original unchanged published Heat-Wave example** inside the
pinned `scrimp-mpi:v1.3.0-rc3` image. Both runs use identical,
reference-qualified one-sided placement of interface control ports
and bricks. They vary **only** the source of the oriented
facets: (1) the historical geometry-specific Concentric helper,
versus (2) `scrimp.interface.resolve_conforming_interface`.

The sandbox wrapper
`sandbox/probe_v130_generic_interface_matrix_parity.py` reuses the
existing native assembly auditor, but captures the **entire**
PETSc stiffness matrix as a sparse NPZ on each run. The driver checks
the full matrix dimensions, stored nonzero counts, and maximum
entrywise absolute difference, with explicit absolute/relative
tolerances. It writes both full matrices, individual logs,
`comparison.json`, and a zip archive of the evidence.

```bash
SCRIMP_DOCKER_NO_SUDO=1 \\
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
bash scripts/run-v130-generic-interface-K-parity.sh
```

**This remains a sandbox test, not automatic assembly.**
The same manual constructor region remapping is intentionally
present in both arms, so any mismatch isolates the geometry-selection
algorithm. A PASS would establish MPI1 assembled-K parity only;
not mass-matrix identity, PETSc TS/CN conservation, DAE initialization,
production compilation, arbitrary FEM-family support or MPI2.

## Full native GetFEM/PETSc assembled-K parity — PASS (2026-11-24)

The user ran `scripts/run-v130-generic-interface-K-parity.sh`
in the pinned Docker image at commit `734b628` and supplied the entire
native evidence archive. Its two original-example stiffness matrices
(`generic/stiffness_audit.npz`,
`reference/stiffness_audit.npz`) were independently compared:

- Matrix dimension: **7689 × 7689**.
- Stored CSR entries: **46935** in each.
- **All CSR indptr, indices and data arrays are exactly identical**.
- Difference matrix nnz: **0**.
- Maximum absolute difference: **0.0**.
- The unchanged published Heat-Wave example and its manually
  selected reference port/brick regions were identical in both arms.
  Only the geometric facet discovery algorithm changed.
- No DAE initialization, CN energy test, or MPI2 result was claimed.

Status:
`PASS_NATIVE_GENERIC_AND_REFERENCE_ASSEMBLED_K_PARITY_ONLY`.
The archive also contains raw logs and both sparse matrices.

## Experimental deferred DPHS compilation (stage 6; NOT qualified)

The source now contains a **test-only, disabled-by-default**
deferred GetFEM compiler. Enable it **only for qualification** with:

```bash
SCRIMP_V130_AUTO_INTERFACE_PILOT=1
```

In this mode, `DPHS` discovers uniquely labeled conforming internal
interfaces from the `Domain` physical-region topology, defers the
GetFEM FEM-variable and brick registration, and compiles them before
initial values or `solve()`, once the original user-defined weak
forms have revealed their tested-variable ownership. Two oriented
GetFEM regions are allocated from currently unused numerical IDs.
The user-facing `Port.get_region()` and `Brick.get_regions()` keep
their original logical labels. Effective assembly regions are stored
privately. Variable and brick declaration order is replayed unchanged.
`Port.compute()` integrates power on its effective one-sided region.

For now, unsupported or ambiguous cases **fail closed**: e.g. global
lumped ports in the deferred model, unlabeled or multiply labeled
shared interfaces, multi-region internal bricks, interface prescribed
sources lacking a residual test owner, and unresolved GWFL macros.
This is not a general mortar or nonmatching-mesh implementation.
The baseline solver is unaffected when the environment switch is unset.

**The experimental compiler has NOT yet passed native assembly, CN,
energy or MPI2 tests. Do not use it for production simulations.**

The next qualification runs the *unchanged original published*
Heat-Wave source both ways: reference manually mapped one-sided
interfaces versus experimental automatic compilation. It compares the
entire GetFEM/PETSc assembled stiffness K, including exact original
DOFs and CSR structure, not just block norms:

```bash
python tests/test_v130_interface_deferred.py
SCRIMP_DOCKER_NO_SUDO=1 \\
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
bash scripts/run-v130-auto-interface-K.sh
```

A PASS would qualify the experimental **MPI1 assembled stiffness K**
only; the strict DAE/CN/physical power and MPI2 gates remain open.

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
