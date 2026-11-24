# SCRIMP v1.3.0 — universal FE initialization diagnostic

**23 November 2026, experimental; scalar/high-regularity native gate pending.**

The `dbfa7f1a` original native pilot qualified exact **affine**
field reconstruction in six intrinsic vector GetFEM families:
RT0/1/2, BDM1/2, lowest Nedelec. Every actual reconstructed
field has relative L2 error between 1.65e-16 and 4.70e-14;
every mass-system residual is at roundoff. The initial
scalar/vector classification question remains open until direct
scalar and higher-regularity tests.

## Single projection algorithm, many finite element spaces

For any supported scalar or vector finite element space Vh, evaluate
the analytical initial field f by **GetFEM generic weak-form
assembly** and solve M c = b, where

```text
M_ij = integral(phi_j . phi_i)
b_i  = integral(f . phi_i)
```

For scalar qdim=1, GWFL uses multiplication `u*Test_u`;
for vector qdim=2 it uses dot product `u.Test_u`.
These are contractions of the same variational principle, NOT
separate scalar and vector *algorithms*. It is equally valid for
Lagrange, intrinsic H(div)/H(curl), and scalar derivative-DOF spaces
provided the FEM mass matrix is well-defined and the source is
representable in the accepted expression language.

The source-only implementation must NOT infer coefficient layouts by
counting nodal values. In particular, Hermite and Argyris use DOFs
related to derivatives and C1 continuity. A projection does NOT
promise preservation of those canonical derivative DOFs, even when
the reconstructed field is a good approximation. This distinction
is important for high-order interpolation and commuting elements.

## Native candidate matrix

`scripts/run-v130-universal-fem-initialization.sh` tests 15 FE
cases on a selected region, original qualified image, MPI1:

- scalar CG1/CG2, DG0/DG2;
- vector CG1/DG1;
- scalar Hermite(2), Argyris, Morley, HCT and reduced HCT;
- nonconforming scalar P1;
- intrinsic RT0, BDM1 and lowest Nedelec.

For each available family it assembles exactly representable
polynomial initial data (degree <= element space), solves a small
dense reference L2 projection, and measures **physical** integrated
GetFEM L2 reconstruction error. It also calls actual production
`DPHS.set_initial_value` on the same FE and compares the old
physical L2 error where the old implementation returns.

The pilot requires the six ordinary CG/DG and three intrinsic
families. Advanced derivative-DOF families may be reported
`UNAVAILABLE_IN_PINNED_GETFEM` **only** if the descriptor cannot be
constructed; a failure to assemble or reconstruct any *available*
space is a hard failure. For HCT-family composite FEMs use
`IM_HCT_COMPOSITE(IM_TRIANGLE(13))`. Other spaces use
`IM_TRIANGLE(13)`, sufficient for the test polynomial degrees.
Acceptance: relative mass algebraic residual < 1e-10 and physical
L2 error < 5e-8.

### Interpolation vs projection

On exact representable polynomials both should reconstruct the
same field, when the original initializer is correct. On general
non-polynomial initial data, however, nodal interpolation and
L2 projection are **mathematically different operators**.
Unifying the production initializer around the L2 projection
therefore requires an explicit user-facing semantic decision,
backward-compatibility testing, and a separate check of the
user's Python/NumPy expression syntax vs GetFEM GWFL syntax.

This **does not qualify** production changes, high-order Hermite/
Argyris canonical derivatives, large sparse PETSc scalability,
MPI2, physical PDE/DAE integration, or arbitrary analytic fields.
Publication-era example files must remain immutable.

## Failed first universal-pilot run and scientific interpretation (23 November 2026)

Original host evidence at source `fea5fb65`, folder
`v130-universal-fem-fea5fb65-M1kq1fqm`:

- Scalar CG1 and DG0 and vector RT0/BDM1/lowest Nédélec,
  as well as reduced HCT, nonconforming P1 and CG1/DG1
  vector spaces: projection itself succeeded and matched
  GetFEM-reconstructed fields to roundoff.
- CG2, DG2, Hermite, Argyris, Morley and full HCT failed
  to assemble **because the diagnostic wrote `^` as a
  polynomial power in GWFL**. The GetFEM parser rejected
  the expressions; these failures do **not** test whether
  their finite element spaces support projection.
- Production `DPHS.set_initial_value` returned for CG1/DG1
  vector fields but produced relative physical L2 errors
  **0.9281246565019555** and **0.7598154957809599**;
  projection errors were 4.15e-16 and 3.75e-16.
  BDM1 produced old relative physical error
  **0.9536667458578806** vs projection 1.78e-15.
- Reduced HCT old scalar interpolation returned with error
  0.07276178095132534 while projection had 6.77e-15.
  Its old behavior is recorded separately, not a projection
  failure or an expected scalar nodal baseline.

The diagnostic now uses explicit products `X(1)*X(1)`,
`X(2)*X(2)*X(2)` etc., and keeps its old-production error
as an independent, visible negative result. Only classical
scalar CG/DG elements are required to preserve the known
historical physical interpolation baseline in this pilot.
**No production method has been changed.**

## 23 November 2026 — universal 15/15 native gate and opt-in production candidate

The original native user ZIP `v130-universal-fem-22890c8e-JnYxY8FX.zip`
has been independently inspected. All 15 scalar/vector/high-regularity
spaces passed on the pinned GetFEM/PETSc image, including all six
advanced families; `unavailable` is empty. The worst relative
physical L2 error is 2.7467e-13 (Argyris), far below 5e-8.

Eight **silent** old-production physical errors are confirmed:
CG1/DG1 vector, Hermite, Argyris, Morley, HCT, reduced HCT and BDM1.
The first four scalar Lagrangian CG/DG positive controls retain
physically correct old interpolation. All original example/published
sources remain unchanged.

The first production candidate adds an **explicit-only**
`DPHS.set_initial_value(name, python_expression, method="l2")`.
This is one mathematical L2 projection for every supported FEM:
the only qdim distinction is between GWFL scalar multiplication and
vector dot contraction. It uses the same GetFEM sparse matrix bridge
and PETSc MUMPS backend as SCRIMP's linear solver. The default remains
`method="interpolation"` because L2 and nodal interpolation are not
equivalent for non-polynomial fields; silently changing the default
would alter scientific initial conditions. A restricted Python AST
translator supports x/y/z, NumPy/math common scalar calls,
arithmetic, powers and explicit vectors; it refuses unsupported
syntax before any GetFEM model mutation. For HCT composite spaces,
callers can provide `projection_mim` to override the integration
method, which must be valid for the FEM.

**Native production evidence PENDING.** The source-only checks
verify dispatch, expression translation and use of sparse PETSc,
not actual native mass projection, MPI collectives, nonlinear
DAE or pH physics. The probe `run-v130-production-fe-initialization.sh`
is the next acceptance gate before a default change.

## 23 November 2026 — initial native sparse production run failed before validation

Original host output `v130-production-fe-l2-d00f42aa-25gzYEa2`,
source `d00f42aa`: all 15 element cases and the Gaussian
rejected with the same native GetFEM message,
`Argument 4 should be a DOUBLE REAL data array`.
This is **not** a failure of the already-qualified GetFEM L2 projection,
nor evidence of a physical FE inconsistency. The output is consistent
with passing a PETSc complex-valued coefficient vector into a
`gf.Model("real")` setter. The source now explicitly converts
real-valued PETSc solutions to contiguous `float64`, validating that
any imaginary component is at roundoff scale (relative max tolerance
`1e-12`) and rejecting anything larger; it never silently discards
an actual complex physical state. The native runner prints complete
per-case tracebacks on any remaining failure and records
`PETSc_ScalarType` in the JSON to confirm the source of the mismatch.
**Native rerun pending; production method opt-in only.**

## 23 November 2026 — verified 16/16 native sparse production and default promotion

The user-provided original archive
`v130-production-fe-l2-693153c1-WB4iwQmI.zip`
has been checked directly: `results.json` reports
`PASS_NATIVE_SPARSE_15_PLUS_NONPOLYNOMIAL` for 15 FEM spaces and one
Gaussian. Maximum exact-polynomial physical relative L2 error
`1.0327239069599589e-14` (Argyris). Gaussian sparse PETSc MUMPS
versus original GetFEM dense reference relative coefficient discrepancy
`2.1230373608471312e-16`. Gaussian approximation error
`0.24622684405878315` is a discrete **field approximation** error,
not the numerical linear-solver error. Legacy nodal vs L2 coefficient
difference is `0.18364836397333328` (18.36%), so a default change
is a genuine initial-data change for nonpolynomial inputs.

**Authorized default change:** The v1.3.0 preparation branch now uses
`DPHS.set_initial_value(name, expression, method="l2")` by default.
Historical interpolation is still available by explicitly requesting
`method="interpolation"`, for reproducibility. Scalar lumped states
remain initialized numerically under either method; they do not
have a spatial FE projection. Unsupported Python expressions fail
explicitly, with no silent interpolation fallback.

The production native pilot now exercises the **no-keyword** default
for all 15 FE families and the Gaussian and separately calls explicit
legacy interpolation for the Gaussian comparator. This is a new
native gate and must be rerun following the default change.

**Evidence boundary:** the original 693153c1 archive proves the
explicit `method="l2"` path on MPI1. The new implicit-default dispatch,
full pH-DAE simulations, Heat-Wave, MPI2 and general PDE energy/power
acceptance are **not** qualified by that archive. No v1.3.0 release
tagging is authorized by the initial-state FE tests. All eight
published research sources remain byte-for-byte immutable.

## 23 November 2026 — MPI RHS scaling correction before parallel gate

The SCRIMP GetFEM calibration helper `_getfem_mpi_rhs_scale()`
returns a **multiplicative** factor: `1` for an already-correct
GetFEM RHS or `1/P` if the MPI-enabled GetFEM build duplicates
a fully reduced RHS `P` times. The initial `method="l2"`
implementation incorrectly **divided** its RHS by this factor,
which was invisible in all MPI1 tests (`P=1`). The projection
now applies this factor by **multiplication**, matching SCRIMP's
existing PDE RHS bridge. A source-only regression checks
`[2,4,-6] * 0.5 = [1,2,-3]` and verifies the production call site.

This corrects the direction of normalization, but the full MPI2
GetFEM/PETSc mass-assembly and physical reconstructed-field gate is
still required before claiming parallel initial-value support.

## 23 November 2026 — next native gate: MPI1 vs MPI2 and lumped default

The native source `9dbace0e` test succeeded on all 16
original GetFEM field cases using the **no-keyword** default
`set_initial_value(...)`, and confirmed that the purely
numeric `lumped_variables` initialization remains valid.
The archival report explicitly recorded `MPI2_qualified=false`.

Source `ba022483` corrected the projection's GetFEM MPI RHS
calibration to multiply the RHS by `_getfem_mpi_rhs_scale()`
(in the duplicated-RHS case, `1/P`), rather than dividing.
A new native original-image comparison gate
`scripts/run-v130-fe-initial-mpi.sh` runs six representative
finite element spaces, the Gaussian and a lumped state on
MPI1 and MPI2, and rejects coefficient mismatch (>1e-9),
inter-rank disagreement, wrong physical fields and mixed
MPI calibration. It retains each original rank's output and
records any differences in raw GetFEM scalar integrations
separately (not conflated with field coefficients).

**Not a DAE/power/energy gate:** even a pass on this test
will not prove coupled Heat-Wave algebraic feasibility,
physics of interface conditions, or multi-step time
integration. No research-published source or output is changed.

## 23 November 2026 — MPI2 failed, corrected incorrect replicated runner

The native pilot at `936ecb2a` passed MPI1 but all FE
mass projections raised PETSc `KSP_DIVERGED_PC_FAILED (-11)`
in MPI2. The runner mistakenly **forced**
`SCRIMP_GETFEM_MATRIX_MODE=replicated`, unlike the
previously qualified SCRIMP MPI distributed operators. In
GetFEM MPI level 2, an individual rank's matrix can contain
only a fraction of the globally assembled contributions;
treating those local fragments as complete replicated matrices
can leave incorrect sparse entries, leading to factorization
failure. This is a **hypothesis**, not proof of singularity.

The runner now explicitly selects `distributed` with
`SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=csr`: this intentionally
uses the established IJV/ADD_VALUES reference transfer for
an unqualified new call site. The native pilot refuses
`replicated`, and a source-only regression protects the
settings. We do **not** change PETSc MUMPS parameters or
modify the production mass bridge to paper over a failed test.
The same MPI1/2 coefficient and physical-field comparisons
must be rerun. Only after the distributed reference passes
should the faster COO path be tested separately.

## 23 November 2026 — distributed MPI1/2 native qualification

User native output `v130-fe-default-mpi-9b8e34c7-wndKnTX1`
reported `PASS_DEFAULT_L2_MPI1_MPI2_COEFFICIENTS_AND_FIELDS`
for CG2 scalar, CG1 vector, RT0, BDM1, Nedelec, Argyris, Gaussian
and a numeric `lumped_variable`. The largest cross-MPI relative
coefficient difference is `2.321088974812009e-11` (Argyris),
below acceptance 1e-9. Both MPI2 ranks report RHS scale 1.0,
and GetFEM integrated reference norms agree at roundoff.
This qualifies **distributed IJV/ADD_VALUES** at MPI1/2,
not the replicated bridge (which had failed with KSP -11).

For safety, the production FE project now refuses MPI>1
when `SCRIMP_GETFEM_MATRIX_MODE` is not `distributed`.
Its default IJV transfer remains native-qualified.
Experimental `SCRIMP_FE_INITIAL_COO=1` is only permitted
with `SCRIMP_GETFEM_MATRIX_MODE=distributed` and
`SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo`; this opts into
the existing PETSc COO mass assembly, which needs its own
MPI1/2 comparison before qualification.

Neither this test nor a future COO transfer pass establishes
coupled DAE feasibility, Heat-Wave energy/power, or
time-trajectory PDE acceptance.

## 23 November 2026 — unify FE projection with standard linear COO bridge

The user identified a conceptual mismatch in the pilot's naming:
COO is the **assembly insertion API**, whereas the resulting PETSc
matrix is AIJ/CSR-like sparse storage; they are not two competing
mathematical FE initialization methods. The previous pilot
temporarily selected the historical IJV/ADD_VALUES transfer, which
was useful as a stable MPI1/2 reference but does not certify the
standard linear-matrix COO path.

The FE L2 mass matrix now follows the **same
`distributed_coo_eligible=True` path as DPHS's initial linear mass
and stiffness** whenever `SCRIMP_GETFEM_MATRIX_MODE=distributed`
(the standard bridge setting defaults to
`SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo`). The extra
`SCRIMP_FE_INITIAL_COO` option was removed; no second FE-specific
switch is needed. Serial `replicated` remains compatible; MPI>1
`replicated` is still rejected because GetFEM's rank-local matrix
must not be mistaken for a fully replicated global one.

The standard MPI1/2 runner now checks COO directly and compares
coefficients, physical field reconstruction, and a numeric lumped
initial value. The historical IJV-vs-COO comparison runner remains
optional research instrumentation, not a production choice.
COO projection in this new call site remains **pending the native
rerun**. No pH-DAE, Heat-Wave trajectory or power balance is
qualified by these source-only changes.

Other non-linear Jacobian/dynamic GetFEM call sites still use a
historical IJV fallback as documented in `convert_gmm_to_petsc`.
This change does NOT assert that all SCRIMP matrix assembly call
sites were already converted to COO.

## 23 November 2026 — final production nonlinear GetFEM bridges to COO (native pending)

The user confirmed the original default FE projection is now MPI1/2
qualified on the standard **distributed COO insertion -> PETSc AIJ**
path, including RT0, BDM1, Nedelec, Argyris, Gaussian and lumped states:
`v130-fe-default-mpi-a78eed06-tYJZ7Chu.zip`, source `a78eed06`.

Repository audit exposed two missed production call sites:
`DPHS.assemble_nl_mass` and `DPHS.assemble_nl_stiffness`
had not passed `distributed_coo_eligible=True` to the GetFEM bridge.
Both are now explicitly COO-eligible, matching their **linear**
counterparts; the shared bridge default was also changed to
COO-eligible. This means nonlinear Newton tangents are assembled via
PETSc COO, not silently forced onto IJV when the global mode is
`distributed`. The PETSc matrices remain AIJ/CSR-like in storage.
The existing COO bridge recomputes `setPreallocationCOO` with
the current stored nonzero graph and owned diagonal at EACH assembly,
before `setValuesCOO`. Dynamic `MatAXPY` still uses
`DIFFERENT_NONZERO_PATTERN` and retains stable PETSc Mat handles.

An explicit `SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=csr` or `ijv`
remains for historical benchmark reproduction only, not normal
production assembly; removing it from all archival benchmarks would
destroy previously generated A/B comparisons. It is not an accepted
permanent performance fallback. Separate native gates on the actual
nonlinear shallow-water Jacobian in MPI1/2 and on an artificially
changing COO graph are REQUIRED before qualification, even though the
linear-mass/initial-state COO bridge was already qualified.

## 23 November 2026 — native nonlinear COO acceptance gate (source-only)

The next native gate `scripts/run-v130-nonlinear-coo-mpi.sh`
checks both 1- and 2-rank distributed COO assembly on the
original pinned Docker image, without changing any published sources.
First, the *same* PETSc AIJ matrix is reassembled with four
alternating GetFEM graphs (uncoupled -> coupled -> uncoupled ->
coupled) with nonzero graph count, stored CSR and numeric values
independently compared across MPI1/MPI2. Then the original
`examples/shallow_water.py` nonlinear mass/stiffness Jacobian is
checked against an independent directional finite-difference
derivative with nonzero time derivative and fully explicit/frozen G.
The native runner traces every actual GetFEM bridge operation as
`effective=coo` and fails if any production bridge uses IJV.

This is mandatory before claiming that nonlinear/dynamic PETSc
COO is qualified. It does NOT qualify the nonlinear DAE initializer,
complete physical time integration, or energy/power trajectories.
Explicit `csr/ijv` is retained only for historical A/B benchmark
reproduction, not as the production assembly choice.

## 23 November 2026 — distributed COO is the *implicit* default

The common `convert_gmm_to_petsc()` and the universal initial-state
projection both now select `SCRIMP_GETFEM_MATRIX_MODE=distributed`
when the environment variable is unset. The standard bridge also
defaults to `SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo`, and all four
production DPHS mass/stiffness call sites are COO-eligible.
PETSc still stores matrices as AIJ/CSR-like sparse arrays; COO
is the insertion/assembly method.

`replicated` remains an explicitly requested compatibility option
for serial scripts; the FE initial projector refuses it in MPI>1.
`csr/ijv` remains an **explicit** archival transfer override for
reproduction of previous benchmarks, not a production default.
This default change must be covered by the next original-image
nonlinear/native physical gate before general release qualification.

## 24 November 2026 — original nonlinear COO run: initial Python sign gap

User original native log `v130-nonlinear-coo-84365587-Lq5WZ2LM`
passed the synthetic MPI1 **dynamic COO graph** and stopped in the
next stage while constructing the original
`examples/shallow_water.py` physical model. The unmodified model
initial condition is
`3. - (np.sign(x-0.5)+1)/3.`. Universal default L2 initial
projection refused `np.sign` in its restricted Python -> GWFL
translator **before** assembling the original nonlinear Jacobian.
The failed stage is therefore neither a COO accuracy test nor
a nonlinear Jacobian numerical failure. MPI2 was not reached.

GetFEM GWFL has `sign(t)`, but GetFEM's native `sign(0)`
is **+1**, unlike NumPy's `np.sign(0)=0`. Its
`Heaviside(t)` returns 1 at zero, so an exact translation of
NumPy's scalar sign is `Heaviside(t)-Heaviside(-t)`.
This retains the original height's defined value exactly at
the discontinuity `x=0.5`. It does not pretend that a smooth
polynomial quadrature integrates a discontinuity without error;
the field is still L2 projected with the chosen GetFEM integration
rule. The same untouched original source expression is now
extracted by AST in a source-only regression test and its
sign-zero semantics checked explicitly.

**Required next gate:** rerun
`scripts/run-v130-nonlinear-coo-mpi.sh` on the native pinned
GetFEM/PETSc image to reach both physical shallow-water Jacobians
and complete MPI1/2 graph comparison. No published example, solver
MUMPS policy or legacy interpolation has been modified.
