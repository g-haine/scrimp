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
