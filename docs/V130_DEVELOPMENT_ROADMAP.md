# SCRIMP v1.3.0 — scope, limits, and further work

Status: **development branch, not yet a qualified release** (18 Nov 2026).
Branch: feature/getfem55-mpi-container-scalability.

## Confirmed results before the v1.3.0 integration

- The GetFEM-to-PETSc distributed COO bridge and compact matrix mode
  eliminated major sequential bottlenecks. Historical single-process
  improvement exceeded 30x on the previously profiled wave workload,
  but a new same-model v1.2.0 benchmark is required for a version claim.
- In 24 complete h=0.01, 744,850-DOF fixed-CN simulations (10 steps,
  3 repeats at 1/2/4/6 ranks, on the user's six physical P-cores),
  parallel MUMPS symbolic analysis ICNTL(28)=2 improved full solve
  time vs the prior sequential ordering by ~20.6% at 4 MPI ranks and
  ~23.3% at 6. Effective MUMPS settings and accepted CN steps were checked.
- On that workload 4 ranks were *faster* than 6 ranks: averages
  ~24.829 s vs ~25.712 s using parallel symbolic ordering.
  MPI scaling remains limited; causes include numerical factorization,
  triangular solve and synchronization/placement effects.
- Native PETSc events showed **one** MatLUFactorSym call and **ten**
  MatLUFactorNum calls for the ten main CN steps. The old 'never'
  reset policy reused symbolic analysis but **not** numerical LU.
- An isolated sandbox static-CN prototype, using persistent PETSc SNES
  lag only after initial consistency, changed main-CN numeric LU
  count from ten to one at h=0.04, 2 ranks, with physical block
  norm agreement near machine precision. The saved observed full
  solve wall time was 1.493 s -> 0.904 s. This prototype is not
  suitable as a universal default because it freezes the operator
  without testing whether coefficients or the TS shift have changed.

## First production on_change qualification: diagnosed and corrected

On 18 Nov 2026, the first complete wave h=0.04 qualification passed
native v1.2.0 1-rank (78.117 s) and current 1/2-rank time stepping
(2.033 s / 1.552 s), but the **native PETSc** audit correctly refused
to certify numeric LU reuse: the 1-rank main TS still performed
MatLUFactorNum 10 times for 10 accepted steps. No validated release
performance claim can be made from this failed gate.

Root cause: PETSc SNESComputeJacobian overwrites the KSP
reuse-preconditioner setting *after* the user IJacobian callback.
Calling PC.setReusePreconditioner(True) from IJacobian was ineffective.
The production fix now assembles into a separate reusable candidate
Mat, compares its complete stored CSR with the previous operator on
all ranks, and **does not touch the factored PETSc Mat at all if
nothing changed**. PCSetUp then observes an unchanged PETSc MatState.
For a real graph or coefficient change the candidate is copied into
the operator and PETSc re-factorizes normally; a graph change forces
a symbolic reset. The candidate and previous CSR snapshots are freed
after solve().

The revised mechanism is **not yet native qualified**. Its added
workspace memory must also be measured on refined meshes. Before
repeating the 78-second original release control, use only current
1/2-rank native smoke:

~~~bash
bash sandbox/run_v130_numeric_reuse_smoke.sh 0.04 1,2
~~~

The smoke enforces equality between main-stage MatLUFactorNum and
exact collective operator changes, and between MatLUFactorSym and
graph changes, plus MatSolve=10, TS convergence and full current MPI
initial/final parity. If the smoke fails, do not proceed to
cross-release performance tests.

## Exact final-step CN shift: two numeric LU factorizations justified

The diagnostic smoke on 18 Nov 2026 (h=0.04, 47,942 DOFs,
one MPI rank, 10 accepted CN steps) showed **two** native
main-TS `MatLUFactorNum` calls. Detailed trace:

- Steps 1--9: PETSc CN shift exactly `1000`; unchanged stored
  CSR values on every step; a single numeric factor is reused.
- Step 10: CN shift becomes `1000.0000000000008` owing to
  floating-point final-time alignment (MATCHSTEP);
  121,141 stored coefficients change, with maximal absolute
  coefficient delta `2.7755575615628914e-16`.
- Native PETSc numeric LU count stays 1 through step 9, then
  increases to 2 at step 10. Symbolic analysis stays at 1.

This is **not** an unnecessary factorization under the specified
exact `on_change` semantics. Do not introduce a hidden numeric
tolerance or freeze the Jacobian in order to force a 10-to-1 result:
the coefficient change is small but real, and other models may
exhibit physically meaningful coefficient changes.

The v1.3.0 acceptance gate now checks that the actual native PETSc
numeric/symbolic LU counts match the **collectively determined**
operator-value/graph updates, while preserving full-vector parity.
For the measured run the expected main TS counts are two numeric LU,
one symbolic LU, and ten triangular solves. Other meshes, ranks or
time grids may legitimately yield different counts, but must satisfy
the same exact update/event relationship. Numeric and symbolic
update counters are retained after solve without storing CSR data.

## v1.3.0 implementation and release gates

- [x] Switch default MUMPS policy to on_change rather than inferring
  permanence from linearity alone.
- [x] Distinguish exact stored CSR graph changes from coefficient
  changes. Invalidate symbolic LU for a graph change; refactor
  numerically for any coefficient change across all MPI ranks;
  reuse both only when values are unchanged.
- [x] Enable the experimentally validated MUMPS parallel analysis
  ordering by default on MPI>1, respecting explicit PETSc overrides.
- [x] Add isolated exact-CSR regression tests.
- [x] Prepare same-source v1.2.0-vs-current sequential benchmark and
  staged 1/2/4/6-rank MPI scaling benchmark with per-step timings.
- [ ] **Run** the new production on_change wave h=0.04 qualification
  (10 CN steps, 3 repeats, full initial/final state comparison).
  The main TS must show native MatLUFactorNum calls equal to
  exact collective operator updates (two for the observed final-step
  rounding case), and MatLUFactorSym equal graph changes.
- [ ] Verify the **numeric transition case** on a changing Jacobian:
  one fresh numeric LU for every actual numeric change, but no
  unnecessary symbolic reanalysis for an unchanged graph.
- [ ] Validate the actual full-state trajectory and discrete
  Hamiltonian/input-power identity (not merely final vector norms).
- [ ] Repeat physically meaningful MPI tests (linear driven wave and
  dissipative Saint-Venant) at 1/2/4 ranks. Ensure 'on_change' does
  not freeze nonlinear/variable-step tangents.
- [ ] Run all unit and MPI regression tests, and Sphinx HTML/LaTeX
  builds. Release tag/version updates come **after** these gates.
- [ ] Record reproducible host CPU topology, MPI/OpenMP/BLAS
  placement, PETSc/GetFEM versions, memory peaks and benchmarks.

## Deliberately deferred

### Why 6 MPI ranks do not beat 4

The h=0.01 PETSc log with the old refactor-every-step mode showed
4-rank vs 6-rank MAIN CN events:

| Event | 4 ranks | 6 ranks |
| --- | ---: | ---: |
| Symbolic LU (one call) | 2.987 s | 3.257 s |
| Numeric LU (ten calls) | 7.569 s | 8.380 s |
| Triangular solve (ten calls) | 0.923 s | 1.060 s |

Those numbers illustrate a plateau across the solver's native phases.
They are **not** an extrapolation to the forthcoming one-numeric-LU
v1.3.0 configuration. Rerun after numeric reuse to see which operations
still dominate. Consider MUMPS multifrontal load balance, etree,
MPI distribution and rank/core placement, and per-step collective
synchronization. Benchmark on >1 NUMA domain or additional machines
only if hardware is available.

### Alternative Schur/field-block solvers

Earlier pH-DAE graph discovery and sparse path-2 Schur prototypes
correctly identified constraints and achieved residuals near 1e-15,
but their distributed nested Krylov/ASM setup was roughly **100--120x
slower than reused global MUMPS** on the investigated small wave
system. This is an experimental negative result. Do not choose Schur
or ASM as the v1.3.0 production default. Revisit only with a
better Schur approximation and realistic large-scale benefit.

### Physics, memory and models

- General nonlinear/state-dependent mass and the viscous
  Saint-Venant power balance remain scientifically demanding; linear
  wave performance does not qualify those models.
- Validate long-time memory stability of both the CSR snapshot cache
  and the LU factor after the on_change integration.
- Native GetFEM MPI-2 generic RHS replication/over-summation remains
  an upstream/build configuration issue (#23) until root cause is known.
  Never divide arbitrary RHS contributions by MPI rank count without
  the existing build-specific calibration and physical parity proof.
- Keep independent v1.3 documentation modernization (#26), Docker
  cache/context optimization (#25), and MPI power/physics completion
  (#6) on their own acceptance gates.

## Repository cleanup policy

The Git history retains all removed experimental implementations and
their runner scripts; the negative results are summarized above.
Keep production tests, small numerical regressions, the v1.2.0
historical source comparator, the new v1.3 wave qualifier, MPI
physics checks and the source files used by active tests.
Do **not** remove large ignored data from the user's working tree,
delete existing benchmark reports, or rewrite published Git history.

## Open issue triage (provisional, not a closure decision)

- #6: MPI power-plot deadlock fix and various algebra/IO corrections
  are represented in the current branch, but the umbrella issue is
  not closed until its complete MPI and physical parity gates pass.
- #23: GetFEM MPI-2 assembly root cause unresolved: keep open.
- #25: Docker context/build cache/development ergonomics: keep open.
- #26: full documentation modernization: keep open; this focused
  solver guide is a first step only.
- #1 and #2 were closed in earlier releases, unrelated to v1.3.0.
