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
  The main TS must show exactly one native MatLUFactorNum call.
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
