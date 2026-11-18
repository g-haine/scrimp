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

## 6. Important: the co-energy/substituted formulation

The other published example is
`examples/wave_coenergy.py` (underscore, not
`wave-coenergy.py`). This **is not the same numerical benchmark**
as `examples/wave.py`: it uses the Concentric domain with two
subdomains, a dissipative algebraic port `e_r`, and a
boundary feedback `U = 0.5*Y`. Direct performance/physics
comparisons across the two examples would therefore be invalid.

It declares:

```python
S.CoState("e_q", "Stress", states[0], substituted=True)
S.CoState("e_p", "Velocity", states[1], substituted=True)
```

There are **no independent GetFEM DoFs** for `e_q` or `e_p`.
The dynamic unknowns `q`, `p` already represent physical
coenergies (stress, velocity). The constitutive material
laws are incorporated into their **time-derivative mass
bricks**, e.g. `q.Tinv.Test_q`, `p*rho*Test_p`.

**A separate constitutive `-Q_c` cannot be extracted**
from this reduced system. The MUMPS-`Q_c` Schur path
developed for `wave.py` is **not applicable as written**.
It must not silently treat the absent `e_q/e_p` as zero
rows, artificially create duplicate DoFs, or form
`Q_c^{-1}`.

The reduced/co-energy formulation instead has the
uncondensed Jacobian

```text
A_reduced(sigma) = sigma M_coenergy + K_Dirac
                   + K_dissipation + K_boundary/feedback
```

where the feedback term may be assembled via the
*nonlinear GetFEM tangent mechanism* even when the
prescribed expression `0.5*Y` is mathematically
linear with constant coefficients. In particular,
`system.linear_stiffness=False` **does not prove that
this example is a nonlinear physical model**.
At the same time, a solver must not infer matrix
constancy purely from a formal linear expression:
runtime dependence must be checked on assembled
Jacobians, including any boundary feedback.

For this variant, investigate a *different* block
preconditioner using the `q,p` dynamic mass/Dirac
couplings, dissipative `e_r` and boundary `U,Y`.
The DG mass block for `q` may permit element-local
solves, while the CG pressure/velocity block and
boundary saddle couplings require more care.
**MUMPS remains the robust global reference** and may
also be used on suitable constant subblocks.

The same guiding principle covers mixed systems
where some co-states have `substituted=True`
and others remain independent. Such cases are
classified distinctly; no unverified Schur solver
is enabled automatically.

### First native co-energy MPI result — INCOMPLETE parity (18 Nov 2026)

User's h=0.12 `wave_coenergy` run at commit `8f498e0b`:

- 1 rank and 2 ranks each converge after two CN steps, with
  5,790 DoFs and the same five registered fields:
  `U` 108, `Y` 108, `e_r` 735, `p` 1,251,
  `q` 3,588. The substituted `e_q/e_p` have no DoFs.
- Both runs report the same **counts** of nonzero
  stiffness couplings between those physical roles.
- The *independent block norm comparison FAILED*:
  `flow_mass_F_and_algebraic_zero_rows/block_frobenius_norm`
  differed by up to **2.64246**, so physical MPI matrix
  parity has NOT been established and the prototype
  must not be enabled.
- Counting nonzero entries does **not** prove identical
  coefficients. This result does not establish whether the
  source is mesh geometry, GetFEM assembly, the COO bridge,
  or our CSR block attribution.
- The `logic_error exception caught` messages (two
  per rank) were consistent with querying GetFEM for
  the deliberately substituted co-states; the diagnostic
  now skips those known-absent interval lookups.
- The `Concentric` Gmsh mesh was *regenerated for each*
  1-/2-rank run; the old JSON did not record file
  identity. The audit now records SHA-256 of its original
  `Concentric.msh` in each rank report, and compares it
  before assessing matrix block norms.
- The next audit cross-checks native PETSc
  `Mat.norm()` with the norms derived from owned
  PETSc CSR rows. This distinguishes a matrix-wide
  difference from a possible indexing/block-read error.

**Update after inspecting the actual Docker environment and old
archive:** `scrimp-mpi:local` uses Python 3.12 from
`/opt/conda/envs/scrimp/bin/python`, GetFEM MPI level **2**,
and PETSc **3.26.1**. The source commit for the failed
coenergy audit was `8f498e0b`, provided by the host's
read-only SCRIMP source mount. The image uses the SAME
GetFEM executable/binary in the one- and two-rank runs,
so an alleged local-vs-Docker version difference cannot
directly explain a **within-image** MPI rank discrepancy.
Neither `gf.__file__` nor MPI level 2 identifies the
native GetFEM C++ build revision.

The original ZIP (`ph-block-audit-wave_coenergy-h0p12-8f498e0b-vHjj3O`)
contains only JSON reports and logs, not raw per-rank
GetFEM tangents or both original mesh files. Per-block
Frobenius ratios of ~0.505–0.507 (mass and q/p
couplings) and ~0.7071 (damping-related blocks) justify
investigating **multiple assembly semantics**, but
do NOT prove an RHS bug: the discrepant quantities are
matrices, not the GetFEM RHS vector. The calibrated
`SCRIMP_GETFEM_RHS_MODE` does not modify matrix values.

**Decisive native experiment, now instrumented:**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_PH_CAPTURE_MATRICES=1 \
  bash sandbox/run_ph_block_audit.sh 0.12 wave_coenergy
```

This leaves production SCRIMP untouched, reruns the same
short 2-step full-MUMPS reference on 1/2 MPI ranks and
captures each rank's **raw GetFEM global CSR tangent**
BEFORE `convert_gmm_to_petsc()`, plus the final PETSc
**owned-row** CSR AFTER the transfer. This is enabled
only when `SCRIMP_PH_CAPTURE_MATRICES=1`.
The driver then compares the **exact sparse coefficients**
using `sandbox/compare_raw_getfem_petsc.py` and writes
`native_getfem_vs_petsc_bridge.json`, before running
the existing (intentionally strict) 1/2-rank parity gate.

Its four central metrics for both `mass` and
`stiffness` are:

1. `native_sum_2r_vs_native_1r`: numerical MPI difference
   present **before PETSc**? If large, examine GetFEM
   assembly, partitioning and mesh/DoF numbering.
2. `petsc_2r_vs_native_sum_2r`: difference created by
   the distributed COO bridge in the two-rank run?
3. `petsc_1r_vs_native_1r`: serial GetFEM-to-PETSc
   value preservation?
4. `petsc_2r_vs_petsc_1r`: end-to-end 1/2-rank
   physical matrix equality?

Gmsh input mesh SHA256 is recorded; if the mesh files
differ, **do not** draw coefficient-parity conclusions.
The script does **not** compensate the results by
an unjustified factor two. A failed final parity gate
is EXPECTED until the root cause is resolved.

**No expensive rerun is required as the next step.**
First inspect existing reports with the standard-library
offline tool:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python3 sandbox/diagnose_ph_block_parity.py \
  sandbox/outputs/ph-block-audit-wave_coenergy-h0p12-8f498e0b-vHjj3O
```

It writes `ph_block_discrepancy.json` in that directory
and prints sorted per-variable norm differences. If the
problem is localized to `q<-q` or `p<-p`, check the
corresponding physical mass brick and mesh; if all
blocks scale together, investigate MPI duplication or
data transfer. These are hypotheses, not established
causes. After the original JSONs have been examined,
rerun the short native audit once to obtain mesh
fingerprints, and if needed do an exact sparse
entrywise matrix comparison.

### Resolved: Concentric's mixed-dimensional region merge (18 Nov 2026)

**Root cause confirmed by the user's two raw-matrix audits** with
*identical Gmsh .msh SHA256*
`ec213dfd79ad4ade52515820295587083d3d30937717792cc0e44dea1f2045b3`
at `h=0.12` and the same Docker GetFEM MPI-2 binary:

```python
# Incorrect (removed in scrimp/utils/mesh.py):
mesh.region_merge(1, 10)
mesh.region_merge(2, 10)
```

These lines merged the **lower-dimensional 1D interface** region
`10` into both **2D volume** regions `1` and `2`.
GetFEM's `gmsh_with_lower_dim_elt` import includes those
lower-dimensional elements, so the mixed region changes both
integration and FEM DoF restrictions. The resulting matrices
were inconsistent between one and two ranks. No generic
`1/mpi_ranks` compensation is justified.

The user edited only these geometry declarations locally
between successive runs. **Both runs printed source HEAD
`3dc1fff4`**, demonstrating that `git rev-parse HEAD` does
*not* identify an uncommitted bind-mounted source. The native
audit now records `scrimp_mesh_source_sha256` and warns when
the host `mesh.py` differs from HEAD.

| Observed quantity, native one/two MPI ranks | Old geometry | Corrected geometry |
|---|---:|---:|
| Total GetFEM DoFs | 5790 | 5724 |
| Dissipation `e_r` DoFs | 735 | 669 |
| Raw GetFEM mass 2r vs 1r rel Fro error | 0.494526 | 5.08731e-18 |
| Raw GetFEM stiffness 2r vs 1r rel Fro error | 0.707107 | 0 |
| Raw-to-PETSc COO error (both rank counts) | 0 | 0 |
| `||K[e_r,e_r]||_F` at one rank | 7.9304e15 | 0.0843037 |

The 66 artificial dissipative DoFs were associated with the
mixed-dimensional region inclusion. The huge dissipative
coefficient was consistent with
`nu=10*(0.36-x*x-y*y)=0` on the internal curve `r=0.6`,
while `M_r` uses `e_r/nu*Test_e_r`. The clean geometry
does not integrate this volumetric term on the 1D interface.
These were not PETSc or MUMPS defects.

**Permanent fix:** SCRIMP now keeps `region_merge(10,11)`,
`region_merge(10,12)` for internal-curve labels, and
`region_merge(20,21)`, `region_merge(20,22)` for the
external boundary, but does **NOT** merge either into
the volume regions. The separately accessible region
`Interface: 10` is unchanged. The regression
`tests/test_concentric_region_integrity.py` enforces
these source-level invariants without native dependencies.
The native 1-/2-rank `wave_coenergy` audit above is the
numerical evidence. No other built-in geometry was changed.

**Separate feedback check.** The original coenergy example
uses `U=0.5*Y`. Its final Jacobian (where the feedback
tangent is assembled) has `||J[U,Y]||_F=0.225193954`
on both rank counts and the linear control block has
`||K[U,U]||_F=0.450387908`. This ratio supports the
expected derivative but is **not a coefficientwise proof**.
The current optional raw-matrix audit additionally saves
the final PETSc TS Jacobian on each rank and checks
`J[U,Y] = -0.5*K[U,U]` and `J[U,U]=K[U,U]`
coefficientwise in `sandbox/compare_raw_getfem_petsc.py`.
The feedback test must be run on a clean source worktree
before the model is used for strong-scaling and power
preservation claims.

**Final native source-committed matrix gate (commit `33996977`,
user-provided ZIP `ph-block-audit-wave_coenergy-h0p12-33996977-Xp0Afk`):**
geometry uses the permanent fix with 5,724 DoFs
(`U=108,Y=108,e_r=669,p=1251,q=3588`).
Both native ranks used identical Gmsh SHA256.
An independent SciPy CSR difference of the full 5,724x5,724
PETSc matrices confirmed:

| Operator | Relative 2r vs 1r Fro error | Value mismatches |
|---|---:|---:|
| Mass `M` | 5.0873087e-18 | 8 floating roundoff |
| Linear `K` | 0 | 0 |
| Complete TS Jacobian `J` | 5.0196435e-18 | 7 floating roundoff |

The report `native_getfem_vs_petsc_bridge.json` also
records `feedback_validation=passed`. For both rank
counts `J[U,Y]=-0.5*K[U,U]`, `J[U,U]=K[U,U]`
have **zero coefficientwise difference**, and
`J[U,Y]` matches between one and two ranks.
This confirms the boundary-feedback tangent is correct.

**Next gate: actual trajectory and physical power parity.**
The scripts
`sandbox/verify_coenergy_power_mpi.py`,
`sandbox/compare_coenergy_power_mpi.py` and
`sandbox/run_coenergy_power_mpi.sh` perform a short
5-step CN test with full-system MUMPS unchanged.
They compare complete saved vectors `q,p,e_r,U,Y` and
two independent Hamiltonian calculations:
`H=.5*z^T*M*z` and GetFEM's original
`0.5*q.Tinv.q+0.5*p*p*rho` volume integrals.

They also check the exact implicit midpoint residual,
all algebraic rows, and the discrete energy balance
including **both** dissipative ports:

```text
P_r = e_r^T K[e_r,e_r] e_r
P_boundary = U^T (-K[U,U]) Y
H_(n+1) - H_n + dt*(P_r(mid) + P_boundary(mid)) = 0
```

This follows the sign convention of the actual assembled
weak form, with `K=system.tangent_stiffness` including
the feedback tangent. Because `e_r` is substituted,
SCRIMP's generic `compute_powers()` does not automatically
account for its dissipative power. The diagnostic explicitly
includes it instead of assuming `plot_Hamiltonian` alone
certifies the balance. We expect two independent model
constraints: `K[U,:]z_mid=0`,
`K[Y,:]z_mid=0`, and the damping constitutive row
`K[e_r,:]z_mid=0`.

Run on the cleaned, committed worktree:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash sandbox/run_coenergy_power_mpi.sh 0.12
```

The method checks only five **fixed CN steps** at this
stage. A pass does NOT establish correctness under actual
TS adaptivity or any speedup; those are subsequent gates.
Only after confirming the physical parity should a
pH-structured coenergy preconditioner be enabled.
Keep MUMPS global as baseline and fallback.

### Code and native diagnostic

`scrimp/utils/ph_block_layout.py` now exposes
`classify_ph_formulation()` and permits **zero
independent coenergy DoFs only if all registered
co-states are explicitly substituted**. Missing
non-substituted variables, duplicate/overlapping
GetFEM intervals and invented effort DoFs are
errors. Plans include:

| Classification | Mathematically valid candidate | Production today |
|---|---|---|
| `implicit_constitutive` | extract and factor `Q_c` with MUMPS, Schur complement | unchanged global MUMPS |
| `substituted_coenergy` | dynamic mass/Dirac + damping/boundary blocks | unchanged global MUMPS |
| `partially_substituted` | hybrid block structure, no blanket `Q_c` assumption | unchanged global MUMPS |
| `other_or_no_costates` | require structure examination | unchanged global MUMPS |

The read-only MPI audit supports both `examples/wave.py`
and `examples/wave_coenergy.py`:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash sandbox/run_ph_block_audit.sh 0.12 wave_coenergy
```

It uses the same short two-step native MUMPS solve,
GetFEM variable intervals and PETSc local CSR block
norm/parity checks on 1 and 2 MPI ranks, but
**changes the expected block-coupling gates according
to the actual formulation**. Reports include
`formulation` and `solver_candidate`; this is
a design diagnostic, not a solver switch or
performance benchmark. The original
`bash sandbox/run_ph_block_audit.sh 0.08 wave`
is unchanged.

**Acceptance to begin a reduced-form Schur prototype:**
native evidence of the actual mass/Dirac,
dissipation and port blocks, plus successful
1/2-rank parity. Then check two distinct TS shifts,
full DAE constraints and physical power/energy
balance before performance claims. Do not demand
the existence of a constitutive `Q_c` block in
a substituted formulation.

## 7. Required first gates before an invasive solver change

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
