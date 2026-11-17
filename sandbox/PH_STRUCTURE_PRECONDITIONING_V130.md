# SCRIMP v1.3.0 — scalable solvers from port-Hamiltonian structure

**Decision (17 November 2026):** The 50% MPI-efficiency goal targets
**fully linear systems with coefficients constant in physical time**,
with the 2D acoustic wave DAE as the primary benchmark. It must NOT
depend on a fixed numerical timestep: users commonly require adaptive
TS/BDF/CN timesteps. Nonlinear systems (such as dissipative
Saint-Venant) require separate performance tests without imposing the
same E_p>=0.5 acceptance gate. Preserve physical power balance.

This supersedes the proposed *numerical LU caching for one fixed CN
shift* in SCALABILITY_V130_AUDIT.md. That opt-in proof of concept was
withdrawn. The observations of six PETSc numeric factorizations in
one DAE-initialization + five-CN-step run remain valid and useful,
but **not** a reason to build a fixed-step-only production feature.

## 1. Structural model and its time-discrete Jacobian

Consider the linear pH formulation

```text
M dz/dt = (J-R) Q z + B u
M_y y  = -B^T Q z

J^T = -J,   R^T = R >= 0.
```

For the full physical energy to be a quadratic form, one additionally
needs the actual discrete energy pairing: for example, with constant
matrices and `W = Q^T M = M^T Q >= 0`,

```text
H(z) = 0.5 z^T W z.
```

This is a **separate assumption**, not implied by skewness of J and
positivity of R. M may be singular in a pH DAE and Q may only be
defined via finite-element mass/constitutive operators rather than a
standalone invertible assembled matrix.

PETSc's *implicit TS Jacobian* for the residual F(t,z,zdot) is

```text
A(sigma) = sigma M - (J-R)Q
         = sigma M + (R-J)Q,
```

where sigma is supplied by TS and varies with the timestep and
integration method. In a fixed-step CN picture, multiplying the
midpoint equation by dt gives the equivalent operator

```text
M + (dt/2)(R-J)Q.
```

The goal is **not** to reuse numeric LU when sigma changes; it is
to exploit the block structure and update a cheap preconditioner
when sigma changes.

If Q is invertible and the above energy pairing is symmetric, the
effort variable `e=Qz` leads to

```text
A_e(sigma) = sigma M Q^{-1} + R - J
           = H_sigma - J,  H_sigma= sigma M Q^{-1}+R.
```

`M Q^{-1}` is symmetric positive semidefinite **only under the
quadratic-energy condition** `Q^T M=(Q^T M)^T>=0`. If H_sigma is
positive definite on the actual (possibly reduced) unknown space,
it may be used as an SPD preconditioner and the transformed
operator has symmetric positive part plus skew part.
This *does not guarantee* bounded condition number or mesh-independent
GMRES convergence: the skew contribution can become increasingly
stiff under mesh refinement. More importantly, **the full SCRIMP
wave DAE contains zero-mass algebraic rows**: a global SPD
H_sigma must NOT be presumed.

## 2. Derived passivity-based numerical acceptance test

For Crank-Nicolson applied to this constant-coefficient model,
`z_m = (z_{n+1}+z_n)/2`, `e_m=Q z_m` and
`r_m = M(z_{n+1}-z_n)/dt - (J-R)e_m - B u_m`.
If the energy pairing and boundary conventions above hold,

```text
 (H_{n+1} - H_n)/dt
   = - e_m^T R e_m - u_m^T M_y y_m + e_m^T r_m.
```

The last term is the *discrete power defect caused by the
algebraic linear-system residual*. It directly connects the Krylov
stopping criteria to the port-Hamiltonian energy accuracy.

For other time schemes and non-constant coefficients, derive the
appropriate discrete balance; DO NOT claim CN's exact discrete
identity for an arbitrary adaptive BDF integration. Variable-step
CN can retain the midpoint quadratic-energy identity on each
accepted step when the model and stage values satisfy its
hypotheses. Boundary flow/control sign conventions matter.

Physical validation must include the **full state vector**,
algebraic constraints, output/port power and the energy balance.
Monitor the KSP *true residual* as well as the reported
preconditioned residual, and the power defect |e_m^T r_m|.

## 3. Priority prototype: fieldsplit / Schur / AMG

The concrete `examples/wave.py` model uses:

- dynamic states `q` (vector DG1), `p` (scalar CG2);
- co-energies `e_q` (DG1), `e_p` (CG2);
- `-M_e_q e_q + C_q q=0` (Hooke);
- `-M_e_p e_p + C_p p=0` (momentum constitutive law);
- mixed boundary control, including a Dirichlet-type
  Lagrange-multiplier output on the left boundary.

The DG mass for `e_q` admits **element-local** inverses.
The consistent CG mass for `e_p` is generally global:
never assert its inverse is local. A lumped or diagonal mass
**may serve as an approximate preconditioner only**, and
must not silently change the original discrete equations.

Study two (possibly nested) block decompositions:

1. Split dynamic variables `(q,p)` from algebraic
   co-energies/ports `(e_q,e_p,boundary)`.
   The block preconditioner must respect the saddle structure
   of the boundary multipliers and the zero-mass rows.
2. Or eliminate appropriate **invertible constitutive blocks
   in the preconditioner**, then construct the wave-equation
   Schur complement with the sign/transpose structure derived
   from `J^T=-J`. Schematically (not an implementation formula),
   the pressure-side reduced operator resembles

   ```text
   S_p(sigma) ≈ A_p(sigma)
            + D^T A_q(sigma)^{-1} D,
   ```

   with consistent mass and coefficients, and the time-shift
   factors dictated by the actual block ordering.

Test `KSPFGMRES` with `PCFIELDSPLIT` (PETSc `IS`
assembled from exact GetFEM variable intervals) using
Schur `full/lower/upper` factors, and AMG (`PCGAMG`
or hypre if built) on the *appropriate elliptic Schur
approximation*. For development, the `selfp` Schur
option uses an approximate diagonal inverse of A00.
This is only appropriate if that diagonal is a sensible
approximation (it is not universally so for saddle
blocks). Compare also a cheap block-diagonal energy
preconditioner as a scientific baseline.

Preconditioner construction can be shift-dependent:
`P(sigma)` is updated/reconfigured when PETSc's TS
shift changes, without assuming numeric LU reuse.
Reuse constant *geometric* information (GetFEM FEM
subspaces, sparsity pattern, indices and optionally
AMG interpolation when legitimate). Do not assume
numeric factors from sigma1 solve sigma2.

If an inner Krylov method or inexact application makes
the preconditioner nonlinear, PETSc documentation
recommends **FGMRES**, not ordinary GMRES; CG and MINRES
are not generally valid on `H_sigma-J`, which is
nonsymmetric.

PETSc references:

- https://petsc.org/release/manualpages/TS/TSSetIJacobian/
- https://petsc.org/release/manualpages/PC/PCFIELDSPLIT/
- https://petsc.org/release/manualpages/PC/PCFieldSplitSetIS/
- https://petsc.org/release/manualpages/PC/PCFieldSplitSetSchurPre/
- https://petsc.org/release/manualpages/PC/PCGAMG/
- https://petsc.org/release/manualpages/KSP/KSPFGMRES/

## 4. Required first gates before an invasive solver change

1. **Extract exact block sizes and signs** from the production
   `GF_MODEL.interval_of_variable` and PETSc `J`,
   including ownership ranges, co-energy laws and ports.
   Verify all relevant discrete transpose/adjoint identities
   under the actual FEM mass pairings. Do NOT assume that
   the continuum J skewness transfers automatically to
   the assembled Euclidean matrix after variable elimination.
2. **Independent block operator check** on 1 then 2 ranks.
   Verify `||J+J^T||`, `||R-R^T||`, quadratic
   energy positivity on the admissible subspace (when
   appropriate), and correspondence to original residual.
3. **Small fixed mesh, variable steps and shift** (at least
   two distinct accepted dt values) using the same operator
   coefficients, compared against the full MUMPS reference:
   exact saved vectors, convergence, constraints, output
   powers, Hamiltonian, and discrete power defect.
4. **Strong scaling** p=1,2,4,6 and h=0.02, 0.01, 0.005,
   with 3+ repetitions, fixed comparison workload and
   CPU topology already checked on the user's machine.
   Report BOTH global runtime and steady-step time.
   The design target E6>=0.5 applies to fully linear,
   time-constant-coefficient systems; separately measure
   nonlinear Saint-Venant but do not impose that gate.

**Do not retrofit an adaptive-unsafe special case as a
production feature.** Keep original MUMPS direct solver
available as the reference and fallback throughout.
