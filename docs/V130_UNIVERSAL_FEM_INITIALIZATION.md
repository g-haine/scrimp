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
