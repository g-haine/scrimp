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

## 1. General pH-DAE model: constitutive relations are implicit

**Correction of an earlier oversimplification.** In SCRIMP the
state/co-energy relation is generally an **algebraic variational
relation**, not the explicit pointwise map `e=Qz`. Abstractly,

```text
F dz/dt = (J-R)e + B u
C z     = E e
M_y y   = -B^T e

J^T=-J,  R=R^T >= 0.
```

We use distinct symbols to avoid overloading the original `M`:

- `F`: dynamic/flow mass or descriptor operator;
- `C` and `E`: the *two* assembled constitutive operators,
  corresponding to the user's generalized law **M z = Q e**
  (here `C=M`, `E=Q` in that law);
- `e`: independent co-energy unknowns, not automatically `Qz`;
- `J` and `R`: interconnection and dissipation operators **in
  power-conjugate effort/flow coordinates**;
- `B, M_y`: input/output-port operators, with the stated sign
  convention. SCRIMP's mixed controlled boundary formulation can
  introduce extra saddle-point rows and unknowns.

We **never assume** that `C`, `E` or `F` are invertible, nor
that `C` or `E` equal an identity matrix, nor that the
finite-element algebraic rows can be discarded. Here `F`
is often singular on the whole augmented DAE state due to
zero-time-derivative co-energy and boundary variables.

For PETSc implicit TS, let `sigma = d(zdot)/d(z)` as provided
by `IJacobian`. Ignoring prescribed `u` and separate output
rows, the **uncondensed** operator on independent `(z,e)`
unknowns is

```text
             [ sigma F     R-J  ]
A(sigma)  =   [                   ].
             [    C         -E   ]
```

Additional boundary/control rows are retained for the actual
acoustic benchmark. The shift `sigma` changes under adaptive
timestepping, even if all physical coefficients `F,C,E,J,R,B`
are constant. A useful preconditioner should be based on
`A(sigma)`, updating or approximating its shift-dependent parts
**without requiring repeated global MUMPS numeric LU**.

If, **and only if**, `E` is invertible on the relevant constrained
subspace, reorder the **same** two block equations by rows and
columns into constitutive/dynamic rows and `(e,z)` columns:

```text
                [   -E           C     ]
A_(e,z)(sigma) = [                       ].
                [  R-J        sigma F  ]
```

The Schur complement with respect to the `-E` block is

```text
e = E^{-1} C z

S_z(sigma) = sigma F + (R-J) E^{-1} C.
```

This is an *exact algebraic identity when E is invertible*.
A scalable `PCFIELDSPLIT` implementation would replace
the hard global inverse of `E` by a suitable approximation
*within the preconditioner* and keep the true mixed
Jacobian as the Krylov operator.

This is a derived *special case*, not the assumed input
formulation. In particular, `E^{-1}` may be nonlocal; an
approximate solve with `E` may be suitable INSIDE a
preconditioner, but the full DAE must remain the physical
operator used by PETSc.

The skew-adjoint `J` and symmetric positive semidefinite
`R` properties guide the **signs, adjoints and choice of
approximate Schur factors**, but they do not by themselves
make the assembled uncondensed DAE, its Schur complement
or its symmetric part SPD.

## 2. Dissipation, ports and the physically meaningful residual

Assume the actual discrete co-energy variables are power
conjugate to the dynamic flow; in particular, for the
quadratic Hamiltonian of a compatible linear discretization,

```text
dH/dt = e^T F dz/dt.
```

This identity is a **separate compatibility condition** to
verify on the concrete FE spaces and constitutive relation
`C z=E e`; neither the algebraic equality alone nor
`J^T=-J` implies it.

Then the pH equations yield

```text
dH/dt = - e^T R e + e^T B u
      = - e^T R e - u^T M_y y.
```

For a *quadratic* Hamiltonian and a CN midpoint step,
with `z_m=(z_{n+1}+z_n)/2` and an admissible
`e_m=(e_{n+1}+e_n)/2` satisfying the midpoint
constitutive law, define separate dynamic and
algebraic residuals

```text
r_f = F (z_{n+1}-z_n)/dt - (J-R)e_m - B u_m
r_c = C z_m - E e_m
r_y = M_y y_m + B^T e_m.
```

If constitutive and port equations are satisfied to their
required tolerances and the midpoint energy pairing holds,

```text
(H_{n+1}-H_n)/dt
   = -e_m^T R e_m - u_m^T M_y y_m + e_m^T r_f.
```

An inexact Krylov solve can also leave **nonzero `r_c`
and `r_y`**, which may cause additional energy/port
defects: `e_m^T r_f` alone does *not* certify the
full DAE solution. Record all three residuals and
the physically computed Hamiltonian/powers. For
variable-step CN, this midpoint identity can still
hold **step by step** under the same hypotheses; do
not assume it holds for arbitrary adaptive BDF
or nonlinear systems.

Consequently, use the port-Hamiltonian power balance
as a **solver-accuracy test**, alongside full-vector
comparison and all algebraic/port constraints, rather
than as a substitute for residual accuracy.

## 3. Priority prototype: fieldsplit / Schur / AMG

The concrete `examples/wave.py` model uses:

- dynamic states `q` (vector DG1), `p` (scalar CG2);
- co-energies `e_q` (DG1), `e_p` (CG2);
- `-M_e_q e_q + C_q q=0` (Hooke);
- `-M_e_p e_p + C_p p=0` (momentum constitutive law);
- mixed boundary control, including a Dirichlet-type
  Lagrange-multiplier output on the left boundary.

The uncondensed blocks above must first be checked against the
actual GetFEM stencil and test-function orientation: in particular,
SCRIMP's `-M_e_q e_q + C_q q=0` and
`-M_e_p e_p + C_p p=0` are **not** the identity
map `e=Qz`. PETSc can build a split from the actual variable
index sets and can apply approximate local/AMG solvers to
constitutive blocks without changing the physical residual.

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

## 4. First implementation step — 18 November 2026

The **read-only native block inspector** is now available. It does
not change the numerical time integrator, MUMPS, the Jacobian, or
any of the physical DAE equations.

Run:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash sandbox/run_ph_block_audit.sh 0.08
```

The launcher performs **two fixed CN steps only to assemble and
retain PETSc matrices** with the existing reference MUMPS solver
for 1 and 2 MPI ranks. It is a *structural audit*, **not**
a prototype that depends on fixed steps. A later Schur PC must
handle arbitrary positive TS shifts `sigma`.

The inspector uses the actual
`GetFEM Model.interval_of_variable()` values for every
registered dynamic state, co-energy and boundary variable.
It rejects missing DoF intervals and overlaps, rather than
assuming `(z,e)` are contiguous global index ranges. On
each rank it inspects **owned PETSc CSR rows** of the
assembled `mass`, `stiffness` and last TS Jacobian.
It aggregates globally only the small **block norm and
nonzero-count arrays**; it never gathers a full PETSc matrix.

Result directory:

```text
sandbox/outputs/ph-block-audit-h0p08-<commit>-<id>/
├── blocks_1r.json
├── blocks_2r.json
├── blocks_1r.log
├── blocks_2r.log
└── ph_block_parity.json
```

The 1-/2-rank comparator checks complete field coverage,
partition ownership, agreement of **all block Frobenius
norms**, dynamic/co-energy couplings and the **zero
time-derivative mass rows** on co-energy and boundary
variables for this specific acoustic model.

**Critical distinction:** Frobenius norms do not establish
entry-by-entry operator parity, skew-adjointness or a
spectral equivalence bound for the preconditioner. The
existing trajectory+Hamiltonian parity tests remain required;
the next gate will check adjoint signs and actually solve
`A(sigma)x=b` for two distinct values of `sigma`,
compared to the exact direct MUMPS reference.

**Initial proposed split** (subject to native measurements):

- block 0: independent co-energy unknowns `e_q, e_p`,
  representing the constitutive `-Q_c` rows;
- block 1: dynamic `q,p` **plus all boundary-control/
  multiplier unknowns**. These form a potentially
  indefinite saddle-point Schur block.

This split is a *hypothesis*, not an automatic solver
configuration: the actual port and boundary couplings
must be checked before attempting an AMG preconditioner.

Read-only helper: `scrimp/utils/ph_block_layout.py`.
Native audit: `sandbox/audit_ph_blocks.py`.
One/two-rank gate: `sandbox/compare_ph_block_audits.py`.
Regression tests: `tests/test_ph_block_layout.py`
and `tests/test_ph_block_comparison.py`.

## 5. Next native gate: MUMPS on constitutive Q and Schur action

Following the successful h=0.08 1-/2-rank audit (12,682 DoFs:
6,189 dynamic, 6,189 coenergies and 304 boundary unknowns),
we can now probe **an actual time-independent constituent block**.

The baseline acoustic `K=system.stiffness`, `M=system.mass`
and the full implicit linear operator `A(sigma)=K+sigma M`
are assembled by the *existing MUMPS-based SCRIMP solve*.
On the actual coenergy-DoF set `e`, test

```text
A_ee(sigma) == K_ee == -Q_c
```

at **two distinct PETSc shifts** (1000, 2000). There is NO
assumption that the physical time step is fixed: the two shifts
are explicitly varied after the reference run.

On the complementary set `x` of dynamic **and every boundary
port/multiplier** unknown, extract the blocks

```text
        [ -Q_c   C(sigma) ]
A    =  [                 ]
        [  D       G(sigma)]
```

and verify, for deterministic distributed test vectors,

```text
e=Q_c^{-1} C x
S(sigma) x = G(sigma)x + D Q_c^{-1} C(sigma)x

A(sigma)[e;x] = [0; S(sigma)x].
```

The **same PETSc KSPPREONLY/PCLU/MUMPS operator** is created for
`Q_c` and then solved for multiple right-hand sides and shifts,
with no `KSP.setOperators`/PC reset in between. Numeric factor
reuse must be **confirmed in the PETSc event log**; identical
operator configuration is not a measured factor count.

The gate also checks the sign, symmetry and positive diagonal
of `Q_c` for this acoustic case and verifies a manufactured
`Q_c e=b` solve. Symmetry + positive diagonal **do not prove SPD**,
nor does the Schur identity guarantee that a cheap scalable
preconditioner exists. Any failure is reported without modifying
the production solver.

Native command (one/two ranks, h=0.08, 2 CN reference steps):

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash sandbox/run_constitutive_mumps_schur.sh 0.08
```

**Docker ownership (revised after successful first native gate):**
the original `/workspace/scrimp` is a path **inside Docker**,
not the host root. Running as the host UID does not make its
image-owned parent writable. The launcher now explicitly mounts
the user's `~/Simulations/workspace` as **`/scrimp-work`**
and sets `--workdir /scrimp-work`. This is where SCRIMP's
import-time `outputs/{log,png,mesh,pv}` are created. The
production source mounts remain **read-only**, and reports and
PETSc logs still go to the unique
`sandbox/outputs/ph-constitutive-schur-*/` campaign folder.
Override the host working directory, if needed, with
`SCRIMP_DOCKER_WORKSPACE=/another/owned/path`.
The launcher sets `PYTHONPATH=/workspace/scrimp` so changing
the working directory does not accidentally load an installed
different version of SCRIMP.

**First native evidence (h=0.08, 12,682 DDL, 1/2 MPI ranks):**
both shifts `sigma=1000,2000` gave
`||Aee + Q||/||Q||=0`, `Q` solve residual ~1.1e-16,
Schur relative action error ~1.1e-16 (one rank)
and 1.5e-16 (two ranks), with all algebraic residuals
~1.4–1.5e-16. PETSc aggregated counters showed
`MatLUFactorSym=3`, `MatLUFactorNum=4`, `MatSolve=6`
on each run: consistent with one Q factorization and
three subsequent applications of the factors, in addition
to the untouched reference TS operations. This does not
yet show a scalable outer preconditioner, or actual
adaptive-TS trajectories. Source: user's native
`ph-constitutive-schur-h0p08-1a808652-iYqIFC` archive.



Output:

```text
sandbox/outputs/ph-constitutive-schur-h0p08-<commit>-<id>/
├── schur_1r.json
├── schur_2r.json
├── schur_1r.log
├── schur_2r.log
└── constitutive_schur_parity.json
```

Read-only numerical probe:
`sandbox/probe_constitutive_mumps_schur.py`. MPI gate:
`sandbox/compare_constitutive_mumps_schur.py`. The pure
NumPy unit test checks the Schur sign with an independent small
matrix: `tests/test_constitutive_schur.py`.

**Next step if native gate passes:** construct a configurable
FGMRES/PCFIELDSPLIT prototype with `-Q_c` solved by MUMPS
and a *shift-dependent approximate Schur block*, keeping
the original full MUMPS route as reference. Then check
full state, constitutive residuals, boundary powers and
Hamiltonians with two genuinely accepted adaptive TS
steps of different sizes before any speed claim.

## 6. Required first gates before an invasive solver change

1. **Extract exact block sizes and signs** from the production
   `GF_MODEL.interval_of_variable` and PETSc `J`,
   including ownership ranges, the separate flow (`F`) and
   constitutive (`C,E`) mass/operators, and every port
   constraint. Do not conflate PETSc's Jacobian object
   `self.J` with the **Dirac interconnection** matrix `J`.
   Verify the actual FE power-duality and discrete adjoints.
2. **Independent uncondensed block check**, first on 1 then
   2 MPI ranks. Inspect `F,C,E,J,R` and boundary blocks,
   verify their actual adjoint/symmetry conditions under
   discrete dualities, test constitutive and dynamic residuals,
   and check energy positivity on the admissible subspace
   where appropriate. The complete PETSc Jacobian is neither
   expected to be skew nor SPD.
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
