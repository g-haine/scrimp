# SCRIMP v1.3.0 — Generic conforming interface orientation (stage 1)

> **Mise à jour — 24 novembre 2026, sur `release/1.3.0-preparation`
> UNIQUEMENT.** Le compilateur des interfaces internes conformes et
> l'export ParaView MPI par le rang 0 sont désormais **des chemins de
> code activés par défaut sur cette branche** (pas de tag, de merge
> `master` ni de publication). Les mentions de « pilote désactivé »
> plus bas retracent l'état des phases historiques et ne décrivent
> **plus** la valeur par défaut du code de préparation. Les qualifications
> natives, tous les échecs/correctifs, les exclusions et les options
> de retour arrière sont consignés dans
> [la décision de production du 24 novembre 2026](
> V130_INTERFACE_PV_PRODUCTION_DECISION_20261124.md).
> Le PASS natif des chemins **sans drapeaux** est désormais
> **confirmé** (code `241f58fb`, 71 tests source et quatre
> simulations MPI1/2 DAE/CN/énergie/ParaView ; archive
> `v130-auto-interface-MPI-CN-241f58fb-YEJXJhSK.zip`).
> Se reporter à la décision de production pour les limites
> et l'anomalie historique de métadonnée du rapport corrigée ensuite.

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

## Experimental fully automatic MPI1 K parity — PASS (2026-11-24)

User-run pinned-Docker audit at commit `f0d3fa7`, with original
Heat-Wave source unchanged, compared manual one-sided reference
assembly to the **production-path DPHS** deferred compiler enabled
only by `SCRIMP_V130_AUTO_INTERFACE_PILOT=1`:

- 10 topology, 6 port ownership, 7 brick placement and 5 deferred
  compiler source-only tests passed (**28/28**).
- `PASS_EXPERIMENTAL_AUTOMATIC_INTERNAL_INTERFACE_K_IDENTICAL`.
- GetFEM/PETSc K shape **7689 × 7689**, nnz **46935** in both runs.
- **Strict bitwise CSR parity**: indptr, indices and data arrays exactly
  equal; difference nnz **0**; max absolute difference **0.0**.
- Compiler spontaneously chose oriented region IDs **21** (heat)
  and **22** (wave) from available labels, preserved logical region
  **10**, assigned all 10 interface bricks correctly, and replayed
  **7 FEM** plus **23 bricks** in original declaration order.
- MUMPS/Gmsh/ParaView source examples remained immutable; production
  automatic mode was disabled in the reference arm.
- Both the run's archive and complete numerical CSR arrays were
  independently examined. This is not DAE/CN/physics/MPI2 evidence.

## Genuine strict DAE / CN / interface-energy qualification (pending)

`sandbox/probe_v130_auto_interface_cn.py` exercises the existing
unchanged published original Heat-Wave example with **actual** automatic
DPHS compilation and **no** runtime port, brick or domain-constructor
region overrides. At postprocessing it obtains effective heat/wave
interface regions (e.g., 21/22) directly from `Port.get_assembly_region()`,
and measures CN surface powers on the correct oriented sides instead of
hardcoded reference labels 30/31 or doubled logical region 10.
The existing physical collector and its numerical balances/thresholds
are otherwise unchanged. The pilot also checks that the normal public
`Port.compute(..., sampling="cn_step")` computes those same sided
physical traces.

`scripts/run-v130-auto-interface-cn.sh` runs the strict native
DAE + four-step PETSc TSCN/MUMPS for **both** boundary-compatible
initialization candidates (`clamp`, `constrained_l2`), each with
the reference manually sided interface and the new automatic compiler.
It compares accepted times, **full FE trajectories**, Hamiltonians,
original state components, individual heat and wave power/energy
defects, nontrivial interface excitation and DAE initial consistency.
Separate gates ensure that physical balance PASS is not inferred from
solver convergence.

```bash
SCRIMP_DOCKER_NO_SUDO=1 \\
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
bash scripts/run-v130-auto-interface-cn.sh
```

Success target:
`PASS_BOTH_AUTO_INTERFACE_CANDIDATES_DAE_CN_ENERGY_PARITY`.

These new native gates have **not been run yet**. They validate only
the specific 2-D `heat_region=1, wave_region=2` Heat-Wave convention;
the reversed physical sign convention from the unchanged published
example still requires an explicit separate mathematical audit.
The baseline default mode remains unchanged, and no MPI2 or arbitrary
higher-order/nonconforming family is qualified.

## Native automatic CN run at a44a80d — solver PASS, collector BLOCKER (2026-11-24)

The user provided the full pinned-Docker evidence archive of
`scripts/run-v130-auto-interface-cn.sh` at commit `a44a80d`.
For both boundary-compatible initializations (`clamp` and
`constrained_l2`), original unchanged Heat-Wave with the experimental
automatic compiler completed the real strict DAE and full four-step
fixed-step PETSc CN with **zero rejected steps and zero SNES failures**.
The algebraic initialization residuals were approximately
`9.620818e-14` and `9.91492022e-14`, respectively. Original
Hamiltonian and public `Port.compute(cn_step)` processing completed.

**Neither automatic case passed the physics postprocessing stage**:
the reference `measure_v130_hw_dual_facet_physics.py` collector
validated the *literal* region labels 30/31 even though the compiler
generated 21/22, and raised `Original corrected 2D facet regions missing`.
The reference manual-facet arms were physically qualified in the
same archive (4 CN steps and nontrivial port exchange); for `clamp`
its total CN energy defect was about `5.46e-15`, and for
`constrained_l2` about `6.89e-15`.

Corrective change: the physical collector now accepts explicit
`heat_faces` and `wave_faces` keyword arguments with historical
defaults **30/31**, validating and using those IDs for **all** power
expressions and reporting. The automatic native pilot passes its
resolved one-sided region IDs as arguments instead of temporarily
overriding module globals. A dependency-light regression suite
`tests/test_v130_hw_facet_collector_regions.py` covers default
reference regions, dynamic 21/22 regions, the absence of global
state leakage, and missing/doubled region rejection.

This correction is **postprocessing only**: no DAE, finite element,
matrix assembly, PETSc CN or user example source was changed.
It does **not** convert the failed postprocessing run into a
physical PASS retroactively. The new complete native run is required:

```bash
SCRIMP_DOCKER_NO_SUDO=1 \\
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
bash scripts/run-v130-auto-interface-cn.sh
```

## Full automatic Heat-Wave DAE/CN energy equivalence — PASS (2026-11-24)

The user ran the pinned Docker physical comparison at commit `b001f2e`,
and supplied its complete 86-file evidence archive
`v130-auto-interface-CN-b001f2eb-7wIwM8qq.zip`. The run is a
**genuine, original, unmodified Heat-Wave** execution with
`SCRIMP_V130_AUTO_INTERFACE_PILOT=1`, comparing to separately mapped
historical one-sided interfaces.

- Status: **`PASS_BOTH_AUTO_INTERFACE_CANDIDATES_DAE_CN_ENERGY_PARITY`**.
- All 32 source-only unit tests passed; immutable example audit passed.
- Both `clamp` and `constrained_l2` passed original strict algebraic
  DAE initial projection, 4/4 fixed CN accepted steps, 0 rejects,
  0 nonlinear solver failures.
- Full accepted state vectors (**5 × 7689**), all saved times,
  Hamiltonians and T/p/q trajectories were **elementwise identical
  to the reference** (maximum absolute difference 0).
- All compared physical energy and interface-power diagnostics were
  **numerically identical** between automatic and reference.
- `clamp`: max total CN step balance defect
  **5.457440055423035e-15**, max interface-power defect
  **8.881784197001252e-16**.
- `constrained_l2`: max total balance defect
  **6.893791093531831e-15**, max interface-power defect
  **3.1086244689504383e-15**.
- Strongly nonvacuous exchange: max individual CN interface power
  about **3.56907**, relative to heat dissipation about **0.11139**.
- Compiled internal heat and wave side regions: **21**, **22**;
  original public logical interface region **10** unchanged.
- Both physical collector and public `Port.compute(cn_step)` used
  the one-sided effective regions.

**Scope:** This qualifies the opt-in automatic interface compiler for
original Heat-Wave **heat=1, wave=2, MPI1 only**, including strict
DAE/CN and nontrivial pH energy identities. The flag is still
disabled by default. It does **not** qualify a swapped heat/wave
physical sign convention, MPI2, arbitrary geometries with dynamic
solvers, mixed/nonmatching meshes or higher-order FEM spaces.

## Physical transmission signs: updated live Heat-Wave example

The numerical example `examples/heat_wave.py` has now been **explicitly
updated with user authorization** at commit `5042705` to follow
the mathematical model with **an outward normal on each bulk subdomain**.

For both heat=1/wave=2 and heat=2/wave=1, the constitutive transmission
uses the **same expressions**:

```python
S.Brick("B_T", "-U_T*Test_e_Q.Normal", [10], position="effort")
hw.set_control("Interface Heat", "Y_w")
hw.set_control("Interface Wave", "Y_T")
```

In the mixed heat/grad-grad wave weak form this represents the signed
interface power `-U_T*Y_T+U_w*Y_w = 0`. The former disk-based
`heat_normal_orientation` factor was part of a legacy **common**
interface normal, not a physical parameter of the paper model.
Now SCRIMP resolves each outward parent-face normal automatically;
the application must **not** apply the legacy factor a second time.

The source change **does not modify the eight scientific snapshots**
under `examples/published/`, which remain protected by
`scripts/check-v130-published-immutable.py`. The historical native
Heat-Wave acceptance evidence recorded above used the earlier live
example source before this user-authorized improvement. For
heat=1/wave=2, both old and new text specify identical forms.

**Important compatibility gate:** the corrected example requires
the oriented-interface compiler for physical validity. The
experimental compiler is still *off by default* in SCRIMP v1.3.0
release preparation; enable it with
`SCRIMP_V130_AUTO_INTERFACE_PILOT=1` for the corrected physics.
The legacy doubled logical region `10` is not an acceptable
substitute for separately outward-oriented face integration.

### Swapped geometric placement versus physical sign diagnostics

`tests/test_v130_heat_wave_outward_signs.py` guards the new
region-number-independent law at source level.
`sandbox/probe_v130_swapped_interface_sign_K.py` and
`scripts/run-v130-swapped-interface-sign-K.sh` compare real GetFEM
stiffness on the heat=2/wave=1 configuration, using current live
example as-is for the `outward` arm, and reconstructing the
pre-5042705 historical signs **only inside the sandbox** for the
`published` arm.

The test must show that the old/new changes are localized to the
`(e_Q,U_T)` and `(U_w,Y_T)` stiffness blocks, with a full sign
reversal there. This is an **assembly-only** contrast; a separate
strict DAE/CN energy qualification for the heat=2/wave=1 variant
is required. This diagnostic has not yet been executed natively.


## Native swapped Heat-Wave outward-sign K blocks — PASS (2026-11-24)

The user ran `scripts/run-v130-swapped-interface-sign-K.sh` at
commit `1c08c78`, provided the archive
`v130-swapped-interface-sign-K-1c08c78b-0cTy2vOb.zip`,
and confirmed:

- `PASS_SWAPPED_OUTWARD_SIGN_CHANGE_LOCALIZED_IN_TWO_INTERFACE_BLOCKS`.
- Heat in **annulus 2**, wave in **disk 1**; one-sided physical regions
  are assigned **22** (heat) and **21** (wave).
- Full GetFEM/PETSc assembled matrix: **8056 × 8056**;
  **390** nonzero difference entries between the reconstructed
  historical and current physical sign conventions.
- Exactly two sign-reversed blocks:
  `e_Q <- U_T` (**156 entries**) and
  `U_w <- Y_T` (**234 entries**).
- All other blocks identical: unexplained nonzeros **0**.
- Source-only suites: 4 outward-law, 10 topology, 6 port-owner,
  7 brick-plan, 5 deferred tests passed. The 8 historically published
  source snapshots remain untouched.

**Scope**: This verifies intended GetFEM *assembly* changes caused
by replacing the historical common-normal sign compensation with
the mathematical outward-normal transmission convention. It is not
proof of strict DAE initialization, energy preservation or MPI2
for the heat=2/wave=1 configuration.

## Swapped Heat-Wave actual strict DAE / CN physical qualification (pending)

`sandbox/measure_v130_hw_dual_facet_physics.py` now accepts
`heat_bulk` and `wave_bulk` explicitly, in addition to the
dynamically resolved `heat_faces` and `wave_faces`; the existing
default test remains heat=1/wave=2. The volume Hamiltonians and
heat dissipation are integrated on the correct physical bulk
subdomains. The signed outer-boundary power is **positive** when
the exterior belongs to the wave equation, **negative** when it
belongs to the heat equation; the published outer control is
homogeneous, and the new test verifies its work is negligible.

The new pinned-image diagnostic
`sandbox/probe_v130_swapped_interface_cn.py` runs the **current**
`examples/heat_wave.py` with heat=2, wave=1 and the opt-in
automatic compiler. Unlike the previous wave-annulus case,
the **outer boundary belongs to the heat equation**, so the
wave velocity field must not be artificially clamped at the
outer perimeter. The test initially uses the original
L2-projected Gaussian fields **without any compatibility patch**.
An incompatible projected initial state will be reported as a
real DAE failure; it will not be converted to an energy PASS.

The test verifies actual DAE, four accepted CN steps, full state
dimension (8056), separate heat/wave/total energy and interface
power balances, nontrivial interface excitation, zero prescribed
outer-boundary work, and public `Port.compute` sided CN power:

```bash
SCRIMP_DOCKER_NO_SUDO=1 \\
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
bash scripts/run-v130-swapped-interface-cn.sh
```

A genuine native **PASS** would qualify MPI1 swapped Heat-Wave
dynamics and physical balance for the revised outward-normal laws,
not yet MPI2 or generic higher-order/multiphysics interfaces.

## Swapped Heat-Wave strict DAE / CN / physical balance — PASS (2026-11-24)

The user executed `scripts/run-v130-swapped-interface-cn.sh` at
commit `3cbc738` and supplied the full 20-file archive
`v130-swapped-interface-CN-3cbc738a-djMlJFrk.zip`.
The live, outward-normal `examples/heat_wave.py` was run with
**heat=annulus 2, wave=disk 1** and automatic one-sided interface
compilation, without ANY initial-data adjustment or geometric/sign
runtime patches.

- **`PASS_SWAPPED_AUTO_INTERFACE_DAE_CN_ENERGY_MPI1`**.
- All source-only tests (10+6+7+5+6+4 = 38) PASS; 8 protected
  `examples/published/` files remain immutable.
- Strict algebraic DAE initial projection: **PASS**, final true residual
  **9.758748179398196e-14**, tolerance 8.218631422610572e-9;
  exact differential initial values preserved; no pseudo-TS.
- True fixed PETSc CN: **4 of 4 accepted steps**, no rejected steps or
  nonlinear failures; full accepted trajectories shape **5 × 8056**.
- Automatic effective heat side region **22** and wave side **21**
  (78 interface facets per side), public logical region **10** unchanged.
- Max total per-step discrete Hamiltonian/heat dissipation/port balance
  defect **1.6792123247455493e-15**;
  heat/wave separate defects **2.3141211169530607e-15** and
  **6.904199434387692e-16**.
- Max interface power continuity defect **8.881784197001252e-16**,
  nontrivial individual CN interface power **4.390457971145539**,
  heat dissipation max **39.50200106717261**, and no positive
  Hamiltonian step increments.
- Homogeneous thermal exterior boundary power max
  **4.663842928891199e-14**.
- The public `Port.compute(cn_step)` heat/wave port powers agree
  exactly with independently integrated GetFEM one-sided CN traces
  (zero maximum discrepancy).
- The native report explicitly records all ten physical gate booleans
  true, the unchanged live physical source and no fake boundary
  compatibility clamp.

**Conclusion:** both assignments of the 2-D Concentric coupled
Heat-Wave system are now scientifically qualified through strict DAE,
physical CN balance and nontrivial power exchange in **MPI1**.
Automatic compilation remains test-only/off by default; **MPI2 and
other geometric/FE configurations are NOT qualified by this result**.

## MPI2 first gate: actual distributed assembled stiffness (native pending)

`sandbox/probe_v130_auto_interface_mpi_K.py` and
`scripts/run-v130-auto-interface-mpi-K.sh` now test the **actual
automatic interface compilation** on the **current outward-normal**
Heat-Wave example for both bulk assignments, on **1 and 2 real MPI
ranks**. The pinned PETSc/GetFEM image and original FEM/weak forms are
retained. At the point immediately before strict algebraic
initialization, each PETSc rank writes ONLY its owned CSR rows and
records its PETSc ownership interval. The driver reconstructs the full
stiffness matrices, validates complete and disjoint row ownership,
checks that each rank sees the same oriented-interface manifest,
and compares both full MPI2 operators to their MPI1 counterparts.
No dense allgather, no rank-zero-only GetFEM assembly, and no fake
MPI2 implementation.

```bash
SCRIMP_DOCKER_NO_SUDO=1 \\
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
bash scripts/run-v130-auto-interface-mpi-K.sh
```

Success target:
`PASS_BOTH_AUTO_INTERFACE_ORIENTATIONS_MPI2_K_PARITY`.
The script archives raw per-rank matrices and metadata plus logs,
and fails rather than claiming a PASS if a rank fails, row coverage
is invalid or a numerical coefficient differs beyond tolerance.

**This gate has not yet been executed natively.** It is strictly
assembled-`K` qualification. Even a PASS will **not** qualify MPI2
initial DAE, Crank–Nicolson trajectories, energy/port-power identities,
parallel ParaView output, other geometries or arbitrary FEM families.
Those are separate later stages; the automatically oriented compiler
remains disabled by default.

## MPI1 vs genuine MPI2 automatic stiffness: source evidence PASS; comparator corrected (2026-11-24)

The user ran the original MPI1/MPI2 native assembly pilot at commit
`83cf20b` and supplied the full archived logs and **six per-rank CSR
matrices** for 4 runs (heat=1/2 × MPI1/2), archived as
`v130-auto-interface-MPI-K-83cf20b5-zgjERGVM.zip`.

The FOUR actual PETSc/GetFEM assemblies and per-rank owned-row
exports **completed** (rank ownership [0:3845,3845:7689] and
[0:4028,4028:8056] for the MPI2 cases). The original offline
comparison raised `TypeError: iteration over a 0-d array` from
passing a **generator**, instead of a sequence, to
`scipy.sparse.vstack()`. Therefore the original driver printed
`FAIL_OR_UNQUALIFIED_AUTO_INTERFACE_MPI2_K`, **not because any
native assembly failed**, but because its comparison script crashed.

Independent direct reanalysis of the original archived CSR matrices
found:

| Bulk assignment | Matrix shape | Raw MPI1/MPI2 stored nnz | Explicit stored zeros MPI1/MPI2 | Canonical numerical nnz, both | Difference |
|---|---:|---:|---:|---:|---:|
| Heat=1, Wave=2 | 7689 × 7689 | 46935 / 46953 | 5490 / 5508 | 41445 | 0 |
| Heat=2, Wave=1 | 8056 × 8056 | 47570 / 47594 | 4381 / 4405 | 43189 | 0 |

PETSc stores respectively **18 and 24 more explicit zeros** in MPI2.
After `sum_duplicates()`, `eliminate_zeros()`, and canonical
index sorting, the **entire indptr, indices, and data arrays are
bitwise identical** in MPI1 and MPI2, for both configurations.
Numerical maximum absolute difference is precisely **0**.
The compiler manifest matches between MPI1/2 and across MPI2 ranks,
with continuous, disjoint owned-row partitions. No numerical tolerance
was needed for this independent retrospective conclusion.

Thus the **original native assembly evidence itself establishes
MPI2 assembled-`K` parity**, but the original on-host script did
not print a PASS due to its postprocessor bug. The fix:

- Changes `vstack(generator)` into `vstack(list)`.
- Introduces `sandbox/v130_interface_csr_parity.py` with an
  explicit canonical nonzero comparison that reports **both** raw
  PETSc stored nnz/explicit-zero counts and canonical numerical nnz.
- Adds `tests/test_v130_interface_mpi_csr_parity.py` covering
  false negatives caused by stored zeros, genuine coefficient
  differences, changed numerical sparsity, wrong dimensions,
  nonfinite values and rank-row stacking.

Re-executing `scripts/run-v130-auto-interface-mpi-K.sh` with the
corrected comparison is still required to give the on-host test its
intended final status
`PASS_BOTH_AUTO_INTERFACE_ORIENTATIONS_MPI2_K_PARITY`.

**Scope:** MPI2 GetFEM/PETSc **assembled stiffness K only**.
This evidence does NOT qualify MPI2 DAE initialization, numerical
CN trajectories, physical port energy or parallel plotting/exports.
Automatic orientation remains disabled by default.

## True automatic-interface MPI2 assembled K — native PASS (2026-11-24)

The user reran the corrected pinned-Docker script at commit `d00052f`
and supplied `v130-auto-interface-MPI-K-d00052f5-AyVWkC3j.zip`,
containing all four genuine GetFEM/PETSc runs, per-rank CSR snapshots,
full logs and the passing comparison. Independently inspected archive
`comparison.json`:

- **`PASS_BOTH_AUTO_INTERFACE_ORIENTATIONS_MPI2_K_PARITY`**.
- Heat=1, Wave=2: full **7689×7689**, 41445 canonical numerical
  nonzeros in both MPI1 and MPI2; max absolute difference **0**;
  CSR arrays bitwise identical after eliminating explicit zeros.
  Owned MPI2 rows [0:3845] and [3845:7689]; port regions 21/22.
- Heat=2, Wave=1: full **8056×8056**, 43189 canonical numerical
  nonzeros in both MPI1 and MPI2; max absolute difference **0**;
  CSR arrays bitwise identical after eliminating explicit zeros.
  Owned MPI2 rows [0:4028] and [4028:8056]; port regions 22/21.
- Raw MPI2 matrices held 18 and 24 additional *explicit zero*
  storage entries; these are not numerical differences.
- 45 tests passed (10 topology, 6 ownership, 7 brick placement,
  5 deferred, 6 CN physical region, 4 outward normal laws,
  7 CSR canonicalization). The eight archived scientific files
  were unchanged.

This **qualifies automatic oriented interface GetFEM/PETSc stiffness
assembly for both original Heat-Wave configurations in MPI2**, but
cannot be extrapolated to full time integration or energy yet.

## MPI2 strict DAE, physical CN and sided port-power gate (native pending)

The new `sandbox/probe_v130_auto_interface_mpi_cn.py` executes the
**current** `examples/heat_wave.py` under the experimental automatic
interface compiler, for Heat=1/Wave=2 and Heat=2/Wave=1, at both
MPI1 and real MPI2. The physical weak forms and FEM definitions are
unchanged. To isolate GetFEM/PETSc time integration and physics from
unqualified parallel file I/O, only the calls to plot the Hamiltonian
and write ParaView files are suppressed in this **sandbox probe**.
All ranks then enter the same original collective
`compute_Hamiltonian()`, `compute_powers()`, and GetFEM
`get_quantity()` physically sided quadratures.

For Heat=1/Wave=2, the separately qualified `clamp` initial-data
candidate enforces zero exterior wave CG1 boundary coefficients
before invoking the *real* strict algebraic DAE initializer.
For Heat=2/Wave=1, the unmodified original Gaussian is already
compatible and no correction is applied. There are **no**
interface-region overrides or physical sign corrections.

The native driver `scripts/run-v130-auto-interface-mpi-cn.sh`
runs 4 cases (two region assignments × MPI1/MPI2) and saves,
**on every MPI rank**, full 5×N accepted FE histories, GetFEM
Hamiltonians, all numerical physical traces, original strict DAE
reports, compiler manifest, and SHA256 fingerprints of rank-global
outputs. It checks exact agreement between the two MPI2 ranks and
compares all MPI1/MPI2 traces with stated numerical tolerances
(`rtol=2e-8`, `atol=5e-10`). Independent energy-identity,
power-exchange, dissipation, homogeneous-exterior and strict DAE
gates are required. Even if native execution or host comparison
fails, raw logs and any available per-rank evidence are retained in
a shareable zip archive.

```bash
SCRIMP_DOCKER_NO_SUDO=1 \\
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
bash scripts/run-v130-auto-interface-mpi-cn.sh
```

Expected status:
`PASS_BOTH_AUTO_INTERFACE_ORIENTATIONS_MPI2_DAE_CN_ENERGY_PARITY`.

**Not yet executed in the native image.** MPI2 ParaView/plot output,
arbitrary conforming or nonconforming FE spaces, and enabling the
compiler by default remain independent gates.

## Native MPI2 strict DAE and CN: solver PASS, postprocessing blocker (2026-11-24)

The user ran `scripts/run-v130-auto-interface-mpi-cn.sh` at
commit `48dc01d` and supplied the complete native archive
`v130-auto-interface-MPI-CN-48dc01d7-wMhkOOwQ.zip`.
**All four real original Heat-Wave simulations completed their
strict algebraic DAE initializations and four accepted PETSc
Crank–Nicolson steps** with zero rejects / SNES failures:

| Heat / Wave | MPI | Real algebraic DAE residual | CN accepted |
|---|---:|---:|---:|
| 1 / 2 | 1 | 9.620818004658717e-14 | 4/4 |
| 1 / 2 | 2 | 9.782179042671012e-14 | 4/4 |
| 2 / 1 | 1 | 9.758748179398196e-14 | 4/4 |
| 2 / 1 | 2 | 9.853778732954403e-14 | 4/4 |

Both MPI2 ranks reached completion and public port powers were
computed successfully. The boundary-compatible heat=1/wave=2 case
used its qualified `clamp` initial state; the inverted case used
the unmodified original Gaussian.

**No MPI2 physical energy PASS can be claimed from this run.**
Each of four runs stopped in the sandbox postprocessing due to an
unconditional second `system.compute_Hamiltonian()`. The original
live Heat-Wave example had already computed each Hamiltonian term
once (five saved accepted instants). `Hamiltonian.compute()`
appends, rather than replacing old values, so the sandbox made
each term history length ten, and the subsequent
`DPHS.get_Hamiltonian()` raised
`ValueError: operands could not be broadcast together with shapes
(5,) (10,) (5,)`.

A second, independent *host-comparator-only* issue then raised
`NameError: name 'host_error' is not defined`, as `host_error`
had only been assigned in the `except` block of the inline
comparator. This masked the proper final failure summary/archive
report (the original native archive itself was nevertheless supplied).

Corrections are limited to sandbox/host test infrastructure:

- If the original example's Hamiltonian is already computed, the
  sandbox reuses it. It computes once only if absent, and checks
  explicitly that every term has exactly one sample per accepted
  time, before proceeding to real independent GetFEM quadratures.
- Initialize `host_error = None` before the comparator's `try`
  to guarantee valid success and failure report serialization.
- Add `tests/test_v130_auto_interface_mpi_cn_comparator.py`,
  which runs the **actual inline comparator Python** on generated
  synthetic MPI1/MPI2 reports and traces, checking success,
  numerical mismatch, native nonzero exit, corrupt metadata and
  archive creation. Synthetic data serve *only* to test the host
  comparator, not to claim GetFEM/PETSc qualification.

The actual auto interface compiler, native weak forms, PETSc
DAE/CN algorithms and original Heat-Wave example are untouched.
**A new full native run is mandatory for MPI2 energy and full-trajectory
parity**, with the original status target
`PASS_BOTH_AUTO_INTERFACE_ORIENTATIONS_MPI2_DAE_CN_ENERGY_PARITY`.

## Host comparator regression preflight: subprocess SIGKILL (2026-11-24)

The user reran `scripts/run-v130-auto-interface-mpi-cn.sh` at
commit `f81f6ff`. All previously qualified 45 source-only tests passed,
but the five newly added *host comparator* regression tests did not:
the test that should have accepted generated synthetic 1/2-rank
reports launched a second Python interpreter through
`subprocess.run([sys.executable, "-c", ...])`. That child process
terminated with **signal 9 (SIGKILL, returncode -9)**, so none of
the four synthetic comparator cases could produce
`comparison.json`. The exact external cause of SIGKILL is **not
established** by this trace; it is not evidence of a failed GetFEM
assembly or DAE/CN solve, because the shell driver stopped **before
starting Docker** at the failed preflight. No new native MPI2 physics
result was obtained from this run.

The source-only test now compiles and executes the **actual inline
Python comparator code in the same test process**, using
`unittest.mock.patch.object(sys, "argv", ...)` and redirected
stdout/stderr, and explicitly catches the comparator's
`SystemExit(1)` to preserve the same success/failure semantics.
This eliminates a second Python/NumPy interpreter while retaining
the intended tests of normal success, numerical mismatch, native
exit failure, malformed rank metadata, and archive generation.

The production MPI2/DAE/CN/energy pilot, the actual inline comparator
it uses, all original Heat-Wave physical forms, and SCRIMP's numerical
routines **are unchanged** by this corrective commit.
The new in-process tests and native MPI2 energy comparison must
still be rerun; no physics PASS is inferred from this patch.

## Automatic Heat-Wave MPI2 strict DAE / CN / physical energy — NATIVE PASS (2026-11-24)

The user executed `scripts/run-v130-auto-interface-mpi-cn.sh`
at commit `c22e274` and uploaded the full native evidence archive
`v130-auto-interface-MPI-CN-c22e2748-MTltxsjk.zip` (26 entries).
Its final status is
**`PASS_BOTH_AUTO_INTERFACE_ORIENTATIONS_MPI2_DAE_CN_ENERGY_PARITY`**.
The archive's 4 independent run logs, six per-rank JSON reports and
six complete NPZ trajectories were inspected independently of the
host summary. All **50** source-only tests and the 8 immutable
published-example source checks pass.

| Physical configuration | MPI2 DOFs | Strict DAE residual | Accepted fixed CN steps | Maximum total CN energy defect | Maximum closed-interface power defect |
| --- | ---: | ---: | ---: | ---: | ---: |
| Heat=1, Wave=2 | 7689 | 9.782179042671012e-14 | 4/4, no rejects/SNES failures | 9.610368056911511e-16 | 8.881784197001252e-16 |
| Heat=2, Wave=1 | 8056 | 9.853778732954403e-14 | 4/4, no rejects/SNES failures | 2.6853519408120974e-15 | 4.440892098500626e-16 |

- Both interface excitations are nonvacuous: peak individual
  CN power **3.5690759218838686** and **4.3904579711455245**;
  maximum heat dissipations **32.04166757122651** and
  **39.50200106717281**, respectively.
- The full **5 × 7689** / **5 × 8056** FE state histories
  agree MPI1 vs MPI2 with maximum absolute differences
  **2.2e-12** and **1.45e-11**; accepted times and
  Hamiltonian totals are identical.
- On each MPI2 run, rank0 and rank1 have **bitwise-identical full
  recorded trajectories**, including all physical GetFEM traces;
  rank SHA256 fingerprints agree.
- All 10 gates per case are true: strict algebraic DAE, accepted CN
  steps, nonzero thermal dissipation, nontrivial interface exchange,
  total/heat/wave energy balances, closed signed interface power,
  negligible homogeneous outer work, and decreasing Hamiltonian.
- Effective one-sided GetFEM regions are **21/22** for heat=1/wave=2
  and **22/21** in the opposite assignment; public interface region
  remains **10**. The live `examples/heat_wave.py` uses the paper's
  outward-normal transmission law; it was not rewritten by the pilot.
- The heat=1/wave=2 case uses its separately qualified clamped
  exterior wave initial condition. The reversed heat=2/wave=1 case
  uses the original, unmodified L2-projected Gaussian.

**Resulting qualification scope:** for this original, conforming 2D
Concentric Heat-Wave model, the opt-in automatic interface compiler is
qualified **end to end** for both region permutations through exact
geometry, GetFEM/PETSc K, genuine MPI1/MPI2 strict DAE/CN and
nontrivial port-Hamiltonian energy identities. There are **no**
parallel export/plot results in this audit, because the qualification
pilot intentionally disabled those user I/O operations. Neither
arbitrary FEM families/geometries nor default-on production use can
be inferred from this very specific model.

## Rank-isolated native ParaView diagnostic (pending)

Now that the original Heat-Wave physical simulation has passed MPI2
in both orientations, the next independent gate concerns the
**file-export path**. Existing `DPHS.export_to_pv()` creates and
removes its output directory only on rank zero, but currently runs
`open(...pvd)` and `MeshFem.export_to_vtu(...)` **on every rank
with the same filename**. Concurrent writes could race, and it
is still unknown whether each rank's exported VTK mesh represents
the complete model or merely its MPI partition. Therefore a
shared ParaView exporter is **not yet qualified**.

An opt-in *sandbox-only* diagnostic has been prepared:

```bash
SCRIMP_DOCKER_NO_SUDO=1 \\
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
bash scripts/run-v130-auto-interface-mpi-pv-diagnostic.sh
```

The driver reuses the same four *native* strict-DAE/CN/energy runs,
but calls the original `export_to_pv()` with a **different parent
path per MPI rank**. Each rank writes real `p` and `T` PVD/VTU
files independently, without collisions. The host postprocessor
validates five timesteps and XML VTK piece metadata, records
SHA256/file sizes/point and cell counts, and compares mesh
coverage to serial exports on each rank. The output's status is
deliberately labeled
`DIAGNOSTIC_MPI_RANK_ISOLATED_PV_FILES_COMPLETE_CN_ENERGY_PASSED`,
**not** a parallel shared-production export PASS.

A successful diagnostic would provide the evidence needed to
choose and validate a correct production policy (e.g. collectively
export rank-specific pieces with a global `.pvtu`, or a verified
root-only complete-mesh export). No policy has been changed in
SCRIMP production by this diagnostic. Serial and MPI2 numeric
solver behavior remains independent.

**This new export diagnostic has not yet been run natively.**

## Rank-isolated ParaView preflight: stale host failure label (2026-11-24)

The first user execution of
`scripts/run-v130-auto-interface-mpi-pv-diagnostic.sh` at commit
`3f1d10c` stopped in the **source-only host comparator preflight
before Docker**. The previous 45 source-only tests passed and four
of five comparator checks passed. One old regression fixture
expected the original failure string
`FAIL_OR_UNQUALIFIED_AUTO_INTERFACE_MPI2_DAE_CN_PHYSICS`,
but the extended host comparator used the new joint
`FAIL_OR_UNQUALIFIED_AUTO_INTERFACE_MPI2_DAE_CN_OR_PV_DIAGNOSTIC`
string even when PV diagnostics were **disabled**.

The host comparator now **preserves the original status when
PV diagnostics are disabled** and uses the extended string only
when `SCRIMP_V130_MPI_PV_DIAGNOSTIC=1`. This is a
compatibility fix to test infrastructure; no FEM, solver,
physics, production export or source example was changed.

Two additional dependency-light tests exercise the **actual
embedded comparator**, with *synthetically generated XML* files,
for isolated-PV success (6 per-rank reports, five PVD timesteps
per p/T and eight serial-vs-parallel mesh-comparison entries)
and for native-failure archiving with the PV-specific status.
All synthetic evidence is explicitly test-only and cannot
qualify the real GetFEM VTU exporter. The test fixtures write
valid JSON/line endings.

The expanded comparator suite now contains **7 tests** (52 total
source-only checks when combined with the other 45).
The corrective preflight and the new genuine four-run VTU export
diagnostic remain **not yet executed in the user's environment**.
No physical or parallel I/O conclusion is inferred from the
failed preflight.

## Actual native rank-isolated MPI2 ParaView: PASS (2026-11-24)

The user ran
`scripts/run-v130-auto-interface-mpi-pv-diagnostic.sh`
at commit `4bb09ed` and supplied
`v130-auto-interface-MPI-CN-4bb09edc-yY0sPh1V.zip`
(98 archived files). The **7** host-comparator regression tests,
**52** total source-only tests and the eight historically published
file invariants passed. The four actual DAE/CN/energy runs passed,
and the final diagnostic status was
`DIAGNOSTIC_MPI_RANK_ISOLATED_PV_FILES_COMPLETE_CN_ENERGY_PASSED`.

An independent audit of the actual native VTU/PVD bytes using
the Python VTK reader confirms, in *both* heat/wave region
assignments:

- Every serial/MPI2 output has five correctly indexed timesteps.
- Every `p` snapshot contains **1641 points / 3154 cells**;
  every `T` snapshot has **9462 points / 3154 cells**
  (the GetFEM export duplicates visualization vertices for T).
- On both MPI2 ranks the **entire mesh** is exported, **not**
  an incomplete per-rank partition.
- Rank0 and rank1 output bytes are **exactly identical** for all
  2 × 2 physical fields × 5 saved instants; this includes VTK
  geometry and field values, not merely counts.
- MPI1 and MPI2 mesh points/connectivity are bitwise identical.
  At t=0.004, decoded VTK Float32 scalar fields agree with
  maximum numerical differences **2.442e-17** (p, heat=1)
  and **1.616e-27** (T, heat=2); the other final fields
  match exactly. VTK reads these files successfully.

**Critical production issue:** The original
`DPHS.export_to_pv()` makes/removes its output directory on
rank0, then executes `open(...pvd)` and `export_to_vtu`
on **all MPI ranks with the same filenames**. Since every rank
has a full and identical output, concurrent writes are
unnecessary and create a race. The rank-isolated audit
demonstrated complete rank0 data, but deliberately never
exercised this shared-file race.

## Opt-in rank-0 shared ParaView export: new native gate (pending)

An explicitly gated, backwards-compatible candidate now exists
in `DPHS.export_to_pv()`:
`SCRIMP_V130_PV_ROOT_ONLY_PILOT=1`. The default path is
unchanged. In pilot mode, all ranks enter the method and its
MPI collectives, but only **rank0 writes** the complete
`.vtu` and `.pvd` files at the shared requested path;
root-only export errors are communicated to all ranks via
`MPI allgather`. This avoids concurrent writes, assuming
the native GetFEM VTK writer itself does not require
all MPI ranks to call it (the new pinned-image test will
verify that directly).

`scripts/run-v130-auto-interface-mpi-pv-root-only.sh`
reuses four native strict DAE/CN/energy runs, actually calls
the pilot production exporter for both p and T on every
MPI rank, and then checks:

- Exactly one complete shared output per field and timestep,
  and no shadow rank1 export directory;
- Full VTK mesh, coordinates and cell connectivity, not just
  PVD presence or mesh sizes;
- Decoded full Float32 scalar field data agree with MPI1
  field arrays within a stated precision-aware tolerance;
- All original MPI1/MPI2 accepted-state, DAE and energy/power
  gates still pass.

The helper `sandbox/v130_pv_root_only_audit.py` parses
the actual GetFEM binary VTK scalar arrays and checks
mesh hashes and scalar data. Two additional synthetic
host-comparator tests distinguish a valid root-only export
from a deliberately corrupted field, including failure
archiving. These tests **cannot substitute for the native run**.

```bash
SCRIMP_DOCKER_NO_SUDO=1 \\
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
bash scripts/run-v130-auto-interface-mpi-pv-root-only.sh
```

Expected status:
`PASS_MPI2_ROOT_ONLY_SHARED_VTU_PVD_COMPLETE_CN_ENERGY`.
**Not yet run natively.** A successful opt-in pilot would
qualify the proposed root-only VTK writer for this Heat-Wave
model. Enabling that policy by default for all supported
MPI cases requires an additional release decision and
regression checks.

### Root-only preflight hardening

Before the first native root-only run, the shared-PV driver received
additional protections:

- The exact host-side comparator is exercised by **9** dependency-light
  tests, including generated GetFEM-like base64 scalar PVD/VTU arrays,
  a valid MPI1/MPI2 shared-file comparison, and an intentionally
  corrupted MPI2 scalar which must fail with a preserved archive.
- The Python helper is loaded from the current checkout by absolute
  path, since the native driver uses Python standard input and the
  unit-test process has a different initial `sys.path`.
- Root-only writer mode is **off by default** in production; it is
  enabled only by the pinned native runner. The new runner also uses
  a bounded 90-second per-case timeout by default to detect a
  potential hidden MPI collective inside GetFEM's exporter without
  waiting for the general CN timeout.
- If this native test passes, it establishes safe root-only shared
  export for the two tested fully replicated GetFEM Heat-Wave cases
  only. Domain-decomposed meshes, vector fields, and truly
  rank-partitioned field histories still require other export
  strategies; no global default policy is being silently changed.

## Native shared root-only MPI2 ParaView: PASS (2026-11-24)

The user ran the actual `scripts/run-v130-auto-interface-mpi-pv-root-only.sh`
at commit `0385237` and uploaded the **74-entry evidence archive**
`v130-auto-interface-MPI-CN-03852372-55Ii0zEg.zip`.
The opt-in rank0 writer in production `DPHS.export_to_pv()` was
exercised on real pinned GetFEM/PETSc MPI1 and MPI2 for both
Heat-Wave permutations, with **four** complete DAE/CN simulations,
physical energy gates and the genuine original VTK writer.

- Final status: **`PASS_MPI2_ROOT_ONLY_SHARED_VTU_PVD_COMPLETE_CN_ENERGY`**.
- All **54** source-only checks (10 topology, 6 ownership,
  7 assembly plan, 5 deferred, 6 CN physical collector, 4 outward
  signs, 7 MPI CSR, 9 host comparator) PASS; eight protected
  historical publications unchanged.
- **Both MPI2 ranks** completed strict DAE initialization,
  four accepted fixed CN steps and real GetFEM energy/power
  measurements. There were no rejects, SNES failures, or
  MPI2 rank discrepancies. The accepted trajectories still
  agree closely with serial.
- Heat=1/wave=2 MPI2: DAE residual
  `9.782179042671012e-14`, total energy-step defect
  `9.610368056911511e-16`, nontrivial interface power
  `3.5690759218838686`.
- Heat=2/wave=1 MPI2: DAE residual
  `9.853778732954403e-14`, total energy-step defect
  `2.6853519408120974e-15`, nontrivial interface power
  `4.3904579711455245`.
- For each heat permutation, the shared output includes exactly
  **5 time-indexed VTUs for p and 5 for T**, plus two PVD collections.
  The `p` visualization has **1641 vertices / 3154 cells**;
  `T` has **9462 vertices / 3154 cells**. All four PVD
  time series are correctly indexed from 0 to 0.004.
- The independently decoded VTK binary **geometry/coordinates
  and cell connectivity agree exactly** between serial and MPI2;
  full scalar values agree within Float32 precision. Over
  all saved times and both physical orientations, the largest
  absolute p-field difference is `2.441726822362823e-17`
  and the largest T-field difference is
  `6.462348535570529e-27`. No shadow rank1 output directory
  was produced.
- Port powers from public `Port.compute()` still agree exactly
  with the independent oriented GetFEM CN traces on each MPI rank.
- The archive's comparison reports `host_comparison_failure=null`
  and `shared_production_MPI_ParaView_qualified=true`.

**Scientific qualification boundary:** the opt-in
`SCRIMP_V130_PV_ROOT_ONLY_PILOT=1` validates a
**race-free, shared directory, full-mesh ParaView export
from rank0** for *these two original Heat-Wave configurations*.
The writer is still disabled by default in SCRIMP production;
the `shared_production_MPI_ParaView_qualified` report field
must be interpreted only in this tested context, **not**
as a universal claim for truly partitioned meshes, non-replicated
distributed trajectories, vector/Tensor VTK output or a
native `.pvtu` decomposition.

Before making this route a production default for a larger class,
qualify `Init` and `Final` export paths, collective failure
propagation, behavior in external models and whether the root
actually owns complete visualization data; otherwise use a
partition-aware export strategy. The automatic two-sided
interface compiler remains opt-in and off by default.

## Remaining generality and production gates

1. Establish a race-free, geometrically complete MPI2 ParaView
   export policy. Validate both **field values and mesh coverage**
   independently of successful numerical solves; do not consider
   rank-isolated diagnostics equivalent to shared-file production I/O.
2. Run **actual assembled and physically balanced** coupled
   heat/elastic or Heat-Wave interfaces on at least one second
   geometry (including two adjacent rectangles with a **partial
   shared boundary**, already topology-qualified) without
   Concentric-specific heuristics.
3. Extend and test exact opposite-normal interface extraction on
   1D and 3D imported meshes, Gmsh physical groups, curved
   and higher-order facets.
4. Qualify RT(k), BDM, and Nédélec elements on GetFEM versions
   where the requested families actually exist. Validate
   intrinsic H(div)/H(curl) traces rather than assuming CG
   or one scalar DOF per mesh point.
5. Test ambiguous/multiple interface physical labels, multi-region
   interface forms, multiple meshes and supported lumped systems.
   Fail closed with actionable errors for unsupported configurations.
6. Run the full public baseline regression, multi-process,
   long-time/large-mesh and release-install checks. Only then
   decide whether the currently opt-in experimental compiler
   can be enabled by default for the subset of proven-safe
   conforming internal interfaces.

**Never infer an interface side from variable names, brick names,
mesh region numbers, geometric radii, or subdomain order.**
Geometrical facet orientation and the user's physical
interconnection signs remain separate responsibilities.
