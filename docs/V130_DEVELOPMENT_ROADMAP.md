# SCRIMP v1.3.0 — scope, limits, and further work

> **User scope decision — 23 November 2026:** local GetFEM MPI-2 issue
> #23: retire the investigation of the old miscompiled GetFEM build,
> without asserting an upstream defect; do not remove tested RHS
> safeguards. Schur/field-block/ASM solver alternatives: **retired,
> no future investigation**, use proven MUMPS. FE descriptors:
> Nédélec-I code parameterized k>=0 but qualified only k<=2;
> Nédélec-II explicitly k=1,2; RT upstream parameterized; BDM
> k>=3 not natively qualified. Windows/macOS: best-effort docs
> for real users, no full native acceptance gate. The new compatible
> controlled rotating SWE **initial-state-generation pilot** is
> staged, not yet native-qualified; see
> [V130_ROTATING_CONTROLLED_INITIAL_STATE.md](V130_ROTATING_CONTROLLED_INITIAL_STATE.md).
>
> **Updated 23 November 2026:** issue #29's **qualified cold checkpoint/restart** criterion has passed (production SWE and second linear pH-DAE; cross-MPI 1↔2; CN energy/work at roundoff). Refer to [V130_CHECKPOINT_ISSUE29_QUALIFICATION.md](V130_CHECKPOINT_ISSUE29_QUALIFICATION.md). This does not imply a general DAE initializer or historical TS state restoration; release v1.3.0 remains untagged pending other gates.

**Port-power issue #30 COMPLETE (23 November):** original
GetFEM/PETSc 2-D wave DAE CN and BDF2 × MPI1/2,
4/4 genuine native public `Port.compute` and
`DPHS.plot_powers` qualification PASS and original ZIP
independently recomputed. Max whole-path MPI FE discrepancy
6.7117e-13 (CN) and 2.6174e-14 (BDF2);
CN model-specific discrete energy step residual 2.10e-19.
`instantaneous` generic nodal sampling and `cn_step`
conventions explicit; BDF trapezoid is an approximation,
not an exact scheme-specific stage integral.
Evidence `v130-issue30-port-generic-native-883b428f-N4zJAtcX`;
see [V130_PORT_POWER_QUADRATURE_ISSUE30.md](V130_PORT_POWER_QUADRATURE_ISSUE30.md).

**Native MPI4 qualification (23 November):** 14/14 original rotating
Saint-Venant PETSc runs on h=.2/.1, MPI 1/2/4 original
full-trajectory/power parity and cold checkpoint 1↔4 **PASS**.
Full original ZIP independently audited; max full trajectory
L2 discrepancy 2.5751e-13, original physical trace discrepancy
2.3270e-13, generic port-power discrepancy 3.4001e-16.
See [V130_ISSUE6_COMPATIBLE_MPI4_QUALIFICATION.md](V130_ISSUE6_COMPATIBLE_MPI4_QUALIFICATION.md).
This does not validate t>75 active normal feedback nor long-time
control; broader #6 remains open.

**Release plan (workstreams 0–9 and acceptance gates):** [V130_RELEASE_PLAN.md](V130_RELEASE_PLAN.md). The present file records the prior technical investigations and measured results; the release plan governs current scope and priorities.

Status: **development / untagged release preparation** (updated 21 Nov 2026).
Current integration branch: `release/1.3.0-preparation`; the historical
`feature/getfem55-mpi-container-scalability` branch supplied earlier performance work.

## Binding priority update — 21 November 2026 (supersedes the earlier ordering)

**P0 / before tagging v1.3.0: issue [#29](https://github.com/g-haine/scrimp/issues/29), generic checkpoint and state-based restart for pH-DAEs.**
Required by future LHMNC simulations to avoid replaying expensive
Saint-Venant rotation preparation. Build this as a **general SCRIMP
feature**, not a Saint-Venant-only numerical trick. Implement a minimal
safe `save_checkpoint`/`load_checkpoint` or equivalent after
a sandbox proof, with native MPI acceptance and documentation.

**Checkpoint design constraints:** save the *complete accepted FE coefficient
state* (differential, coenergy, boundary/algebraic variables), physical
time `t_*`, source/model/parameters, exact mesh topology/geometry and
global DoF ordering, FEM layout, time-dependent control provenance,
TS method, scientific status and checksums in an **atomically committed
versioned format**. No serializing live PETSc TS/KSP/MUMPS factors.
On load, verify the identical physical model and compatible mesh, restore
GetFEM + distributed PETSc values, continue with **absolute time**
`t_0=t_*`, reconstruct operators and **skip** pseudo/BDF
consistency initialization for this already compatible state.
Check original algebraic rows `F_a(t_*,z_*,0)` collectively; the full
`F(t_*,z_*,0)` need *not* vanish for a nonstationary state.
Reject incompatibility, corrupt files and unsupported schemes collectively
without silently changing prescribed dynamic states or re-running spinup.
A checkpoint of only h/p is insufficient. This MVP is a **cold
state restart**, not an exact serialization of BDF/ARKIMEX history.

**Acceptance:** native same-mesh continuous-vs-split
`0→t_*→t_f` full-state/constraints/Hamiltonian/true weighted
boundary-work/dissipation/mass comparison and MPI
`1→1,1→2,2→1,2→2` restarts; negative tests for corrupted/incomplete
files, changed mesh/FEM/control/time and unqualified TS history; one
additional linear DAE case to establish generality. Save the actual
mesh/order or prove identical fingerprints; do not assume a regenerated
Gmsh mesh has the same global DoF numbering.

**Maxwell 3D workstream 3: P2 / optional before v1.3.0,
otherwise immediate post-release.** Do not remove its 3D H(curl)/H(div),
orientation, Poynting-power and MPI scientific acceptance criteria.
Only start the *new model/benchmarks* if #29 and all other
release-blocking work are green and time remains.
The previously qualified lean RT/BDM/Nédélec GetFEM simplex patches and
existing published Maxwell examples remain supported;
their compatibility tests stay in the release checks.
Any release without a **new** validated 3D Maxwell demonstration must
say so explicitly; do not silently claim that gate 3 passed.

**Full C1 Saint-Venant spinup new evidence (not checkpoint qualification):**
the user's archived native complete-ramp campaign
`sandbox/outputs/v130-issue28-full-spinup-c0b8f3d2-NSdC2dDW/`
(source `c0b8f3d2`, immutable rc3 Docker) contains all **12/12**
mesh h=.2/.1 × MPI 1/2 × N=10/20/40 physical TSCN runs to
`t_f=.02`, full original FE vectors and power traces, and strict
MPI/full-vector comparator. True endpoint TS/DAE checks PASS
(max `||F_TSCN||≈1.60e-10`, max `||F_a||≈3.25e-14`),
mass defect `<1e-13`, and empirical dynamic h,p order
**1.84153** (h=.2) / **1.96182** (h=.1).
Relative to exchanged energy, the final physical weighted
`ΔH + D - W` defect improves with N but remains
about **2.80e-5** / **6.94e-5** at N=40,
rather than the early t=.002 roundoff floor.
Do not call the state at `t=.02` a stationary rotating solution:
only the **wall velocity** reaches its plateau. Continued monitoring of
interior transients and energy accuracy remains a separate LHMNC
qualification. The original prescribed rotating-at-t0 data remain
incompatible on tested meshes, not repaired by a checkpoint.

This priority change supersedes the release-plan historical instruction
that *all* LHMNC-enabling work starts only after gates 0–9:
the **generic checkpoint** is now a pre-release gate; the full
LHMNC scientific programme remains post-release.

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

The revised mechanism has now **passed a native production smoke
at h=0.04 on one and two MPI ranks** (see results below). The
additional workspace memory must still be quantified on refined
meshes and during long trajectories. The short reproducibility gate
remains:

~~~bash
bash sandbox/run_v130_numeric_reuse_smoke.sh 0.04 1,2
~~~

The smoke enforces equality between main-stage MatLUFactorNum and
exact collective operator changes, and between MatLUFactorSym and
graph changes, plus MatSolve=10, TS convergence and full current MPI
initial/final parity. It passed on commit b7f3c7c6, archived under
`sandbox/outputs/v130-numeric-smoke-h0p04-b7f3c7c6-rrUp27/`.

| Production on_change, h=0.04 | 1 rank | 2 ranks |
|---|---:|---:|
| Accepted CN steps | 10 | 10 |
| Main CN IJacobian calls | 10 | 10 |
| Exact unchanged-operator reuse decisions | 8 | 8 |
| Required numerical factor updates | 2 | 2 |
| Observed MatLUFactorNum | 2 | 2 |
| Required symbolic factor updates | 1 | 1 |
| Observed MatLUFactorSym | 1 | 1 |
| Observed MatSolve | 10 | 10 |
| Complete solve time | 1.563 s | 1.264 s |
| Max process RSS | 588.0 MiB | 525.6 MiB |

Full initial PETSc vectors are bitwise identical; the full final
vectors differ by L2=2.463e-12, relative L2=1.659e-13 and
Linf=5.089e-13, consistent with the specified 5e-8
scaled parity threshold. Numeric LU counts correspond to the
last-step rounding-induced Jacobian change, not a solver failure.
`png/` and `svg/` contain human-readable factor-decision and
timing plots. These are single-run timings, **not** replicated
v1.3.0 performance estimates. They do not establish the power
balance, nonlinear-tangent correctness or long-run memory behavior.

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

## Qualified v1.2.0 → corrected candidate sequential comparison (18 Nov)

After explicit propagation of the accepted PETSc pseudo-TS state
to GetFEM variables, the complete native wave h=0.04, ten-step
Crank–Nicolson v1.2.0 (tag 88e2c144) vs corrected current
(commit 24e70e92) qualification **PASSED**.

The user-supplied complete performance archive was inspected
directly: full initial PETSc vector is **bit-for-bit identical**
(47,942 DOFs), and full final state differences are L2
`6.215e-30`, scaled L2 `3.924e-31`, L-inf `1.281e-30`.
One native source-matched serial run of each version gave:

| Measure | v1.2.0 | Corrected v1.3.0 candidate |
|---|---:|---:|
| Full wall solve incl. assembly | 73.268 s | 1.512 s |
| Initial CN step | 393.7 ms | 381.3 ms |
| Steps 2–10 mean | 131.1 ms | 53.6 ms |
| Maximum single-process RSS | 592.2 MiB | 582.6 MiB |
| Mass assembly + matrix transfer | 15.15 s | 0.0452 s |
| Stiffness assembly + matrix transfer | 56.08 s | 0.1693 s |

**Raw 48.45x now passes the full-state numerical regression gate**
for this configuration, but was measured **once**; do not label it
a three-run median. Stage timings show that the dominant benefit
is GetFEM/PETSc assembly and transfer, not a 48x speedup of
the underlying PETSc factorization. Main current TS: 1 symbolic
+ 2 numerical MUMPS LU events and 10 solves, matching exact
operator-change decisions. The final-step extra numerical LU
is from CN `MATCHSTEP` shift roundoff.

This resolves the previous ~0.376/0.397 cross-release full-state
mismatch. Its root cause was accidental loss of accepted pseudo-TS
state synchronization when prescribed RHS-only sources were
correctly excluded from nonlinear stiffness flags. Both source
classification and new exact MUMPS policy are retained.

**Remaining separate scientific gates:** the pseudo-TS step perturbs
`q(0)=0` to `||q|| ≈ 0.0047176`, so verify actual DAE
constraint satisfaction and physical initialization semantics.
Qualify full trajectories and power balance, corrected 1/2/4/6
MPI comparisons, larger h=0.01 and long-run memory before
declaring the optimization programme fully complete.

## Final corrected v1.3.0 strong scaling at 2,966,786 DOFs — PASS

The **19 Nov 2026** large-wave qualification at h=0.005 used the
corrected current-only source commit `25837de0`, fixed-CN dt=.002,
ten steps, 1/2/4/6 MPI ranks, and **five independent repetitions**
per configuration (20 native runs). **No v1.2.0 runs** were made.
The local archive is
`sandbox/outputs/v130-large-wave-h0p005-25837de0-CBukYK/`.

| MPI | Solve median (s) | Total median (s) | Speedup | Efficiency | Largest rank RSS HWM (GiB) | Sum rank RSS HWM (GiB) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 134.679 | 142.195 | 1.000x | 100.0% | 14.05 | 14.05 |
| 2 | 109.466 | 116.766 | 1.230x | 61.5% | 9.20 | 18.34 |
| **4** | **86.239** | **94.724** | **1.562x** | **39.0%** | **6.14** | **22.17** |
| 6 | 91.493 | 101.151 | 1.472x | 24.5% | 5.17 | 25.64 |

The full five-run 4-rank solve interval [85.655, 86.660] s
and 6-rank interval [91.063, 93.093] s do not overlap. **Four
MPI ranks are the measured optimum on this node**. Going from
4 to 6 ranks increases solve wall time by 6.1% and the sum
of independent per-rank RSS high-water marks by 15.7%.
The latter is a memory-capacity proxy, **not** a simultaneous
whole-job resident peak.

Native PETSc stage logs for **all 20 runs** were independently
parsed from the user-supplied reduced ZIP (states omitted).
Every main CN TS performs 1 symbolic LU, 2 numeric LU, 10
triangular solves; the native event counts agree with
exact collective `on_change` coefficient/graph decisions.
The original native analyzer compared all complete PETSc
initial/final vectors before the user removed the bulky state
files from the shared ZIP: maximum reported scaled L2 errors
**1.04e-14 initial** and **7.56e-13 final**, vs 5e-8 tolerance.
These state errors cannot independently be recalculated without
the removed NPZ; they remain supported by the archived native
validation report and per-run fingerprints.

Main reasons for the **4→6** rank regression (median-stage
observations; inclusive profiler clocks must not be summed):

| Phase / event | 4 ranks | 6 ranks |
|---|---:|---:|
| DAE initial consistency | 23.77 s | 26.40 s |
| Main PETSc TS solve | 52.48 s | 55.14 s |
| Main native MUMPS symbolic LU | 15.98 s | 18.08 s |
| Main native MUMPS numeric LU (2 calls) | 9.38 s | 9.07 s |
| Main native MatSolve (10 calls) | 3.73 s | 3.47 s |
| GetFEM stiffness phase | 2.93 s | 3.55 s |
| Stiffness GetFEM→PETSc bridge | 4.32 s | 3.70 s |

Therefore the six-rank plateau is associated with initial
consistency, symbolic ordering and other overheads, **not**
a failure of numeric LU scaling or reuse. We record this
limit rather than introduce speculative new solver changes
before v1.3.0. Twelve PNG/SVG scientific performance figures
were successfully generated, including scaling, per-step,
startup, PETSc phases, bridge, LU and memory.

On 1/2 ranks, the separate h=0.12 coenergy wave physical
test also **PASSED** with relative discrete power defects
`7.04e-16` and `3.02e-16`, full trajectory-state relative
MPI discrepancy `4.19e-16`, and independent GetFEM Hamiltonian
parity. That establishes the claimed discrete identity for
the **tested dissipative linear model** only.

**Point-0 performance qualification is complete for these test
cases**. Do not silently close separate open questions on
the DAE pseudo-time perturbation of prescribed `q(0)=0`,
nonlinear-tangent invalidation and long-duration memory
stability; these remain documented limitations/tests.

## v1.3.0 implementation and release gates

- [x] **P0 #29 — native-qualified generic checkpoint/restart** before
  release: full accepted FE state and physical-time/mesh/control
  provenance, atomic save, cold restart without algebraic reinitialization,
  strict original `F_a` compatibility, same/cross MPI native split-vs-
  continuous trajectory and physical power gates, negative controls
  and one additional DAE. Do not mark completed by saving an NPZ alone.
- [ ] **P2 — new Maxwell 3D model deferred conditionally**:
  before v1.3.0 only if all P0 gates pass and time permits; else
  immediately post-release, with the original 3D scientific tests
  retained. Existing Maxwell examples and lean FE patches must not regress.

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
- [x] **Pass** production on_change h=0.04 smoke with one and two
  MPI ranks: 10 CN steps, full initial/final MPI vector comparison,
  numeric and symbolic native PETSc LU counts equal exact collective
  matrix update decisions.
- [x] **Pass** matched original-tag v1.2.0 serial vs corrected
  v1.3.0 candidate sequential benchmark: identical full initial
  state, final-state relative L2 ~3.924e-31, valid measured
  single-run 48.45x full-wall-solve ratio and PNG/SVG reports.
- [x] Repeat on corrected source with 1/2/4/6 MPI ranks,
  5 native repetitions at 2,966,786 DOFs. All 20 native LU and
  complete initial/final state-parity gates PASS; four-rank
  performance optimum and RSS capacity cost documented.
- [ ] Verify the **numeric transition case** on a changing Jacobian:
  one fresh numeric LU for every actual numeric change, but no
  unnecessary symbolic reanalysis for an unchanged graph.
- [x] Validate a full-state 1/2-rank coenergy trajectory,
  discrete Hamiltonian and damping/boundary power balance at h=0.12
  (5,724 DoFs; five CN steps): ~1e-16 relative power defects.
  This is not proof of arbitrary nonlinear/time-varying behavior.
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
They are **not** an extrapolation to the newly validated
exact on_change v1.3.0 behavior (two numeric LU events on this
10-step test owing to final-time shift rounding). Rerun after numeric reuse to see which operations
still dominate. Consider MUMPS multifrontal load balance, etree,
MPI distribution and rank/core placement, and per-step collective
synchronization. Benchmark on >1 NUMA domain or additional machines
only if hardware is available.

### Alternative Schur/field-block solvers — retired, NOT planned

Earlier pH-DAE graph discovery and sparse path-2 Schur prototypes
correctly identified constraints and achieved residuals near 1e-15,
but their distributed nested Krylov/ASM setup was roughly **100--120x
slower than reused global MUMPS** on the investigated small wave
system. This experimental negative result is retained only as
historical evidence. **Decision 23 November:** no future research,
optimization, benchmark or deferred roadmap item for Schur/ASM;
this is outside the project scientific remit. MUMPS remains the
qualified production direct solver.

### Physics, memory and models

- General nonlinear/state-dependent mass and the viscous
  Saint-Venant power balance remain scientifically demanding; linear
  wave performance does not qualify those models.
- Validate long-time memory stability of both the CSR snapshot cache
  and the LU factor after the on_change integration.
- The historical GetFEM 5.4.2 MPI-2 raw RHS/partition anomaly
  (#23) arose in an untrusted locally miscompiled older build.
  **No upstream bug is established or further investigation planned**;
  issue closed as not planned, not fixed. Keep the actual build-specific
  RHS auto-calibration and physical MPI parity guards. Never divide
  arbitrary RHS contributions by rank count without validation.
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
- #23: old local GetFEM MPI-2 misbuild investigation **closed as not planned**; no upstream root-cause finding or blanket removal of calibrated safeguards.
- #25: Docker context/build cache/development ergonomics: keep open.
- #26: full documentation modernization: keep open; this focused
  solver guide is a first step only.
- #1 and #2 were closed in earlier releases, unrelated to v1.3.0.

## 23 November 2026 — #29 accepted-state checkpoint/restart: source-only sandbox pilot staged

**New source staged on `release/1.3.0-preparation`, NOT yet native-qualified**:
`sandbox/v130_checkpoint_state.py`,
`sandbox/probe_v130_checkpoint_restart.py`,
`sandbox/compare_v130_checkpoint_restart.py`,
`tests/test_v130_checkpoint_state.py`,
`tests/test_v130_checkpoint_restart.py`,
`scripts/run-v130-checkpoint-restart.sh`,
`.github/workflows/v130-issue29-checkpoint-host.yml`.
Production `scrimp/dphs.py`, native SWE weak bricks,
original wall velocity and original DAE initializers
have **not** changed.

**Checkpoint format pilot:** one versioned
`np.savez_compressed` file containing
the complete real IEEE64 GetFEM coefficient vector
and canonical SHA-256 JSON manifest.
Manifest records accepted physical time, PETSc
TSCN cold-state-only contract, full state checksum,
and original model signature. Signature hashes
actual GetFEM mesh point IDs, convex-order
geometry/connectivity, native MeshFem
`char()` element definition, per-convex
basic DOF ordering, ordered FE DOF positions,
full variable intervals, weak-form expressions,
parameter values, regions, port types
and physical boundary controls. Unlike checking
only h/p DOFs or mesh element counts,
this refuses mismatched algebraic variables or
a different FE basis. Source provenance is
compared on all MPI ranks.

**Atomic publishing:** rank 0 creates an
fsynced temporary NPZ in the target directory
and commits using `os.link` to avoid
overwriting an existing accepted checkpoint.
MPI rank state/time/model SHA parity is
checked collectively; rank 0 errors are
broadcast to all. No pickle, old PETSc
TS/KSP/MUMPS factor matrices, or BDF2
history are saved. Corruption, incomplete
payload, wrong physical time or different
model/control/mesh must fail closed.
This is a **cold restart**, not bitwise
restoration of an integrator's internal stage
history; only PETSc endpoint CN is
supported by this first pilot.

**Native pilot** (new, independent from the
already-qualified ramp evidence): same original
curved-P2 physical SWE `h0=5,p0=0`,
same `gamma=0`, `omega=.2`,
C1 rising ramp over .02, original
decreasing envelope after .03,
original actual full nonlinear mass.
For each mesh h=.2/.1:
- **seed**: exactly 10 accepted real
  `ts_type=cn` PETSc steps `0→.02`,
  actual original `Fa` and TSCN checks,
  save entire accepted original FE z(t*).
- **continuous**: exactly 12 accepted
  physical steps `0→.024`.
- **resume**: reconstruct a new model
  and fresh TS on ranks 1/2 from a
  verified checkpoint saved on ranks
  1/2, set absolute `t0=.02`,
  restore all algebraic and dynamic
  GetFEM and replicated PETSc GlobalVec
  coefficients **bitwise** in the
  sandbox-only initialization hook,
  check original `Fa(.02,z*,0)`
  before physical time integration,
  execute precisely two accepted
  steps to `.024`. No pseudo/BDF
  or algebraic re-initialization is
  called in this pilot.

The native runner schedules **16 actual
independent GetFEM/PETSc executions**
(2 meshes × [2 seeds + 2 continuous
+ 4 cross-rank resumes]), using the
existing immutable original
`scrimp-mpi:v1.3.0-rc3` image
OCI `1d6d0c90115213298047b28237a7f38087563a35`.
The comparator independently checks all
original FE state coefficients and
their continuity at the source checkpoint,
original per-step TSCN and DAE constraints,
MPI1→1,1→2,2→1,2→2, and the
total true original h-weighted
`ΔH+D−W` and mass balance of
continuous vs interrupted/restarted
physical simulations. Real PETSc
accepted step counts/time are asserted
on both the saved and resumed runs.
Data are saved in host-owned
`sandbox/outputs/`, not in a
disposable Docker filesystem.
These strict sandbox tolerances are
initial gates, subject to review if
real native results expose a scientific
or implementation problem; never
weaken them just to print PASS.

**Run preflight and actual first native pilot:**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python -m py_compile \
  sandbox/v130_checkpoint_state.py \
  sandbox/probe_v130_checkpoint_restart.py \
  sandbox/compare_v130_checkpoint_restart.py \
  tests/test_v130_checkpoint_state.py \
  tests/test_v130_checkpoint_restart.py
python tests/test_v130_checkpoint_state.py
python tests/test_v130_checkpoint_restart.py
bash -n scripts/run-v130-checkpoint-restart.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-checkpoint-restart.sh
```

**Evidence status:** synthetic host regression
source is staged but its execution has not
yet been independently observed; native
**0/16 executed for this new checkpoint gate**.
The already-qualified full C1 ramp 12/12
separately remains valid and is the
source scientific motivation, not proof
of checkpoint restart. No `DPHS`
production API has been implemented,
no generic linear second model proven
and no claim of exact BDF/IMEX history
restoration. The #29 P0 release gate
therefore remains **OPEN**.

### 23 Nov 2026 — #29 first host preflight: 25/26 tests passed, two fixes staged (native still not run)

At checkpoint pilot source `7715c538`, the user confirmed `py_compile`
of all new Python files and `bash -n` of the launcher succeeded.
`tests/test_v130_checkpoint_state.py` executed **16 tests, 15 passed,
1 failed**. The only failure was a synthetic fake MPI rank whose
`allgather()` returned malformed evidence flagged as *invalid source*;
the test incorrectly expected the *valid but divergent source* error
`MPI inconsistent`. The actual fail-closed rejection occurred correctly.
`tests/test_v130_checkpoint_restart.py` ran **9/9 PASS** (synthetic
host fixtures, not native GetFEM).

The checkpoint schema tests now separate two cases explicitly:
(1) two structurally valid ranks with different accepted data must
raise `MPI inconsistent`;
(2) one malformed/invalid MPI rank must raise
`MPI checkpoint source rejected`.
The corrupted-NPZ read path now explicitly owns the file descriptor
using `with open(..., 'rb')` around NumPy loading so a deliberately
truncated file raises without `ResourceWarning: unclosed file`.
The native launcher now executes **both** host suites, in addition
to compiling all Python source, BEFORE starting Docker. There are
now **17 checkpoint-format + 9 replay = 26 synthetic host tests**.

The fixes were committed to the source branch; **post-fix host tests
and 16 native cold-restart cases have not yet been rerun**. No
physics, PETSc/GetFEM or numerical tolerances changed. #29 remains
OPEN and the public DPHS checkpoint API is still unimplemented.

### 23 Nov 2026 — #29 first native pilot: original script silent exit traced to grep; checkpoint intact; verified evidence resume staged

**User second independent native attempt** at source `d3c269304cf0fe7fe6b46d45844f50b1120a84c2`,
immutable `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`,
archive `sandbox/outputs/v130-issue29-cold-restart-d3c26930-nRhUY0qA/`,
returned to the shell immediately after printing
`RUN mesh=0.2 phase=seed MPI=1` with no
`FAIL` or `COMPLETE`. The uploaded ZIP was
inspected directly: it **does** contain the
original complete native `h0p2_seed_mpi1.log`,
`segment_report.json`, full
`segment_states.npz` and an accepted full
`h0p2_source_mpi1/accepted_state.npz`.
The actual 10-step curved-P2 original PETSc
TSCN spinup `0→.02` completed successfully,
no TS rejections or failed SNES steps,
`max||F_TSCN||≈8.1411e-11`,
`max||Fa||≈3.2423e-14`, original physical
weighted power defect `≈3.4775e-7`.
The seed's 3300-vector is byte-identical to
the accepted checkpoint and matches
the reported SHA-256 and model fingerprint.
**No native restart is yet present in this archive.**

**Exact root cause: launcher-only grep.**
`scripts/run-v130-checkpoint-restart.sh` used
`grep -aE 'CHECKPOINT (PHASE|RESTORE) ' "$log"`
after successful Docker completion. Actual
success line is `CHECKPOINT PHASE=seed ...`
(**equals** after PHASE), which did not match.
Under `set -euo pipefail`, grep exit code
1 aborted the entire runner silently and
prevented the other 15 configurations.
This was **not** a GetFEM, PETSc, DAE,
checkpoint serialization, or boundary-power
failure. The original complete seed remains
qualified as **creation only**.

**Source fix staged**: correctly match
`CHECKPOINT (PHASE=|RESTORE )`, but
**never rely on an unvalidated grep as proof
of a successful physical case**.
New `sandbox/validate_v130_checkpoint_case.py`
independently checks each native JSON,
NPZ original complete FE vectors, per-step
TS/DAE/power metrics, mandatory success
log line, original source checkpoint
hash/time/physical control provenance
and FE fingerprint before a case may
be skipped or accepted. The runner now
also supports `--resume <existing-host-evidence-dir>`.
For resume, it **preserves original manifest,
source SHA, log and checkpoint** and
accepts the old evidence only if the
source commit is an ancestor and original
production physics + pilot checkpoint
schema/probe files are identical. An
existing partial/corrupted case FAILS
CLOSED (no overwrite). Missing cases are
then executed and independently validated.
The new host synthetic fixture covers all
16 phase-skipping scenarios and rejects
missing PHASE=, wrong checkpoint source and
incomplete seed evidence. CI compiles the
validator and runs the amended host tests.

**Next user run should REUSE the already
computed seed, not discard it:**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_checkpoint_state.py
python tests/test_v130_checkpoint_restart.py
bash -n scripts/run-v130-checkpoint-restart.sh

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-checkpoint-restart.sh --resume \
  sandbox/outputs/v130-issue29-cold-restart-d3c26930-nRhUY0qA
```

**All post-fix host cases and native
cross-MPI restarts are still PENDING** until
the user executes them. Do not close #29,
claim a production DPHS checkpoint API,
or claim a successful native MPI restart
from the existing single-seed evidence.

## 23 November 2026 — #29 original Saint-Venant native COLD checkpoint/restart pilot **16/16 PASSED**, production API staged separately

**Qualified native sandbox evidence**:
`sandbox/outputs/v130-issue29-cold-restart-d3c26930-nRhUY0qA/`,
source `d3c269304cf0fe7fe6b46d45844f50b1120a84c2`,
pinned immutable Docker `scrimp-mpi:v1.3.0-rc3`
OCI `1d6d0c90115213298047b28237a7f38087563a35`.
The user supplied actual final ZIP
`v130-issue29-cold-restart-d3c26930-nRhUY0qA-restart.zip`,
which was independently opened: **16 original
real PETSc/GetFEM phase reports, 16 full state+power
NPZ trajectories, 4 accepted full FE checkpoints,
one immutable source manifest and complete
`checkpoint_comparison.json`**. ZIP CRC intact,
comparison status
`NATIVE_COLD_CHECKPOINT_SPLIT_VS_CONTINUOUS_PILOT_PASSED`.
All 17/17 checkpoint-state + 13/13 continuous/restart
synthetic host regression controls PASS.

Experiment: original curved P2 dissipative
Saint-Venant with nonlinear mass, physical
facet-tangent wall C1 rest-to-rotation
from `h0=5,p0=0` (the previously
qualified compatible sandbox initial state),
`dt=.002`: real PETSc `0→.02` seed
(10 accepted), `0→.024` reference
(12 accepted), and `.02→.024`
fresh PETSc TSCN from the full saved state
(2 accepted). Both physical meshes
h=.2/.1 and **all 1→1,1→2,2→1,2→2**
source-to-target MPI ranks independently
tested (8 true cold restarts). No
PETSc TS/KSP/MUMPS history serialized;
no pseudo/BDF init/projection on restart.

**Bitwise checkpoint creation:** all 4
accepted-state full vectors match their
native seed's original `states[-1]`
exactly and SHA256 checksum passes:
3300 DOFs for nominal h=.2,
10272 DOFs for h=.1, at absolute
physical `t_*=0.02`. Each contains
full original dynamic and algebraic
GetFEM coefficients and immutable mesh/
FE/order/weak-form/control provenance.

**Observed genuine continuous-vs-split
FULL FE vector L2 differences at `t=.024`:**

| Nominal mesh | MPI 1→1 | 1→2 | 2→1 | 2→2 |
|---|---:|---:|---:|---:|
| h=.2 | 1.79515e-13 | 2.48915e-13 | 2.05195e-13 | 2.10372e-13 |
| h=.1 | 1.67651e-13 | 1.52415e-13 | 1.69267e-13 | 1.71498e-13 |

Hamiltonian **restart-time jump** max
`7.55e-15`. Independent original
h-weighted boundary+viscous+Hamiltonian
cumulative physical `ΔH+D−W`
difference between uninterrupted and
two-segment cold restart max
`3.02007e-14`; all algebraic rows
at resumed t* were checked to about
`1e-15`, and the two resumed
original endpoint TSCN steps and
volume balances were accepted
without failed SNES or rejected steps.
This establishes **native scientific
validity of the narrow sandbox
complete-accepted-state cold restart
for this P2 viscous SWE experiment**.

**Not established by this test:**
general arbitrary high-index DAE
initialization, long-term nonlinear-mass
TSCN power/accuracy, a stationary
interior rotating equilibrium, exact
BDF/IMEX history restart,
repartitioned/disparate FE numbering,
or the PRODUCTION `DPHS` API.
Do not close #29 on the sandbox pilot
alone or silently qualify a BDF restart.

**NEXT SOURCE ONLY: generic DPHS public
checkpoint API staged** (no native
results for this implementation yet):
- `scrimp/utils/dae_checkpoint.py`
  promotes the native-compatible
  versioned atomic full FE + model/
  mesh/FE DoF/weak-form fingerprint
  checkpoint reader/writer into the
  package.
- `scrimp/dphs.py` now exposes
  `save_checkpoint(path, extra=...)`
  after a fully accepted final TS
  and `load_checkpoint(path)`
  on a FRESH assembled model.
  The latter stages a complete state
  without mutating the physical initial
  differential variables until
  after the real PETSc/GetFEM
  assembly; before TS creation the
  original algebraic rows
  `F_a(t_*,z_*,0)` are rechecked.
  The solve then uses absolute t*
  and skips pseudo/BDF/algebraic
  consistency initialization.
  It is COLD state restart only,
  with no integrator history.
- Separate
  `sandbox/probe_v130_checkpoint_public_api.py`,
  `sandbox/compare_v130_checkpoint_public_api.py`,
  `sandbox/validate_v130_checkpoint_public_api_case.py`,
  `scripts/run-v130-checkpoint-public-api.sh`
  repeat true 16-case original
  physical SWE continuous/split and
  MPI cross-rank campaign using
  the actual **production DPHS
  save/load methods**. Native pending.
- `tests/test_v130_checkpoint_public_api.py`
  has host-only API-method AST,
  old/new file-schema interoperability
  and fail-closed physical-control
  regression checks. CI
  `v130-issue29-checkpoint-host.yml`
  now compiles the public API and
  runs the new suite.
- **Scientific gate still OPEN:**
  test new generic public API on a
  **second independent linear DAE
  model** and qualify its
  full native continuous/restart
  energy and algebraic constraints.

Preflight and **NEW API** native study:
```bash
cd ~/Simulations/scrimp
git pull --ff-only

python -m py_compile \
  scrimp/dphs.py scrimp/utils/dae_checkpoint.py \
  sandbox/probe_v130_checkpoint_public_api.py \
  sandbox/compare_v130_checkpoint_public_api.py \
  sandbox/validate_v130_checkpoint_public_api_case.py \
  tests/test_v130_checkpoint_public_api.py
python tests/test_v130_checkpoint_state.py
python tests/test_v130_checkpoint_restart.py
python tests/test_v130_checkpoint_public_api.py
bash -n scripts/run-v130-checkpoint-public-api.sh

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-checkpoint-public-api.sh
```

**Important:** public API is currently
source-only and may need further
corrections from first real
GETFEM/PETSc runs. The earlier
16/16 **sandbox** evidence is not
a substitute for native validation
of the new production entry points.

## 23 November 2026 — #29 production DPHS API physically native-qualified on SWE (16/16 PASS), independent linear DAE staged

**Qualified new production API evidence (NOT merely the earlier
sandbox hook):** user-supplied native archive
`v130-issue29-public-api-restart-260aeb55-nq2ksRwp.zip`,
original runner source `260aeb558ffec8783044e7c93c2f02504d75eb41`,
Docker `scrimp-mpi:v1.3.0-rc3` with immutable OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
I independently inspected the ZIP, verified its CRC,
**16 native phase reports and full FE/power trajectories**,
the complete four SHA-256-verified accepted-state
checkpoints (`t_*=.02`, respectively **3300/10272 DOFs**),
and `checkpoint_public_api_comparison.json`.
All four saved checkpoint vectors are byte-identical to
their own native PETSc TS accepted seed endpoints, and
**all eight resumed initial FE vectors** are bit-identical
to their source rank's original complete DAE state.

All **16/16 original curved-P2 viscous SWE** phases using
the REAL production `DPHS.save_checkpoint()/load_checkpoint()`
API PASS; both meshes h=.2/.1, MPI1/2, continuous
`0→.024` vs seed `0→.02` and fresh PETSc CN
`.02→.024` for source→target MPI
`1→1,1→2,2→1,2→2`.
Total host **35/35** (17+13+5).
Max full original FE vector L2 endpoint discrepancy
**2.2588229756e-13**.
Max Hamiltonian jump at t* **7.5495165675e-15**.
Max signed ORIGINAL h-weighted
`ΔH + D − W` difference continuous vs
split **3.0200668355e-14**.
At restored t* original `F_a` is O(1e-15).
The existing true PETSc original TSCN and
algebraic every-step native checks PASS with no
SNES/rejection failure. All four MPI change
combinations verified with actual new model/TS creation.
The archive has source fingerprint and full model
metadata; no general arbitrary index/high-index
initializer is claimed.

**Qualification status:** `DPHS` production checkpoint
API is native-qualified **for this nonlinear
viscous SWE experiment**; independent second-model
DAE still pending. Thus #29 remains OPEN and
generalized cold checkpoint is not yet
release-qualified. PETSc TS BDF/IMEX stage history
is NOT serialized; CN cold state restart only.

**Nonphysical audit cosmetic issue fixed:** the previous
host evidence checker imported
`scrimp.utils.dae_checkpoint` normally,
triggering `scrimp.__init__ → petsc4py.init(sys.argv)`
with pure-checker arguments `--evidence`, `--phase`,
`--mesh` and `--source-ranks`. PETSc printed
warnings about *unused options* after otherwise
successfully verified cases. The host-only
`sandbox/compare_v130_checkpoint_public_api.py`
and `sandbox/validate_v130_checkpoint_public_api_case.py`
now load the exact pure NumPy/JSON module by
`importlib.util.spec_from_file_location`,
without initializing PETSc. The comparator's prose
no longer incorrectly calls the real production
API result a sandbox-only experiment.
**No original SWE physics or PETSc TS numeric
source changed; do NOT rerun 16 SWE native cases
solely to eliminate a CLI-warning display.**

**NEXT native scientific release gate, second
genuinely independent linear pH-DAE**:
- `sandbox/probe_v130_checkpoint_linear_dae.py`
  copies the native linear 2D wave physical
  Dirac/constitutive/boundary bricks/FEM
  from `examples/wave.py` with a smaller
  `Rectangle` mesh h=.5; real states
  `q,p`, efforts `e_q,e_p`,
  four mixed boundary ports and real
  time-dependent left velocity input
  `U_L=.2*sin(25*t)*sin(pi*y)`.
  Exact compatible initial zero q,p,
  boundary input zero at t0,
  `init_step=False` (no pseudo or
  BDF initializer). The forcing grows,
  so this is a NONTRIVIAL linear DAE
  evolution, not a zero-state roundtrip.
  Same CN dt=.002 and physical
  `0→.02` seed, `0→.024` continuous,
  `.02→.024` fresh production
  `load_checkpoint` resume.
- `sandbox/compare_v130_checkpoint_linear_dae.py`
  validates all real TS accepted steps,
  independent true original full
  GetFEM IFunction endpoint CN equation,
  **all algebraic rows** at every
  state, complete FE state/checkpoint
  SHA256, source/target rank control
  provenance, original Hamiltonian,
  true GetFEM boundary-port powers,
  cumulative boundary work and
  continuous-vs-split matching with
  strict fail-closed gates.
- `scripts/run-v130-checkpoint-linear-dae.sh`:
  **8 new native cases** (2 MPI1/2
  seeds, 2 continuous, and 4 cross-rank
  resumes). Immutable original Docker
  image; host-owned logs, `--resume`
  verified partial-evidence reuse,
  `--compare-only` audit.
- `tests/test_v130_checkpoint_linear_dae.py`:
  7 host positive/negative synthetic
  controls; the CI compiles source,
  checks Bash and runs them.
  No host synthetic result is to be
  represented as native physical success.

**Native linear wave gate is entirely PENDING**.
Given its real physical/FE differences from SWE,
errors in its first actual GetFEM run must
be debugged on this second model, not hidden
through a weaker comparator. No release/tag
or #29 closure until original linear DAE
and negative controls are native/host qualified.

```bash
cd ~/Simulations/scrimp
git pull --ff-only

python -m py_compile \
  sandbox/probe_v130_checkpoint_linear_dae.py \
  sandbox/compare_v130_checkpoint_linear_dae.py \
  tests/test_v130_checkpoint_linear_dae.py
python tests/test_v130_checkpoint_linear_dae.py
bash -n scripts/run-v130-checkpoint-linear-dae.sh

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-checkpoint-linear-dae.sh
```

**Release priority:** checkpoint #29 still P0,
new Maxwell 3D P2 optional before v1.3.0,
otherwise post-release. Original already-rotating
incompatible Saint-Venant p0 remains unchanged;
restored saved state is physically compatible,
but reaching the C1 wall plateau does NOT imply
interior stationary rotation.

### 23 November 2026 — linear second-model host preflight: 6/7; diagnostic specificity corrected (native pending)

The user ran the staged second-independent-model
`tests/test_v130_checkpoint_linear_dae.py`
at source `5c9d924`:
Python `py_compile` for both linear DAE probes and
test file PASS, `bash -n` for the 8-case native
runner PASS; **7 synthetic host tests executed,
6 PASS, 1 FAIL**. The only failing test,
`test_rejects_wrong_model_fingerprint`,
correctly rejected a tampered resumed segment
report, but the comparator merged its model
SHA mismatch into a broad
`Actual DPHS linear DAE resume projected/changed
q0 or violated Fa` error, whereas the negative
test required the more accurate
`physically different model` reason.

**Fix staged on the source branch**:
`sandbox/compare_v130_checkpoint_linear_dae.py`
checks the accepted checkpoint's exact original
FE/weak-form model fingerprint against the
restored report *before* inspecting full state
byte equality or DAE constraints; seed source
fingerprint checks are likewise isolated.
This retains every previous strict
fail-closed safeguard. New explicit tests
separately cover (a) a **syntactically valid
but physically different 64-hex-digit model
hash** in restart, (b) wrong seed-model
fingerprint, and (c) a genuine false
`no_pseudo_BDF_or_projection` restart flag
with otherwise identical physical model.
There are now **9 synthetic host regression
tests** for the second model. No SWE/native
physics, PETSc solver, checkpoint format or
numerical acceptance threshold changed.

**Status**: new 9/9 host tests are NOT YET
confirmed by a rerun; **zero native second
linear DAE trials have yet been observed**.
Issue #29 remains OPEN. Actual production
DPHS checkpoint API's earlier original
Saint-Venant native **16/16 PASS** evidence
remains unchanged.

Revalidate before Docker:
```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_checkpoint_linear_dae.py
bash -n scripts/run-v130-checkpoint-linear-dae.sh

# Only if the host suite passes 9/9:
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-checkpoint-linear-dae.sh
```

## 23 November 2026 — #29 independent linear wave pH-DAE: 8/8 native TS PASS; corrected original CN port-work **without re-running** (final comparator revalidation pending)

**Immutable original native evidence (user-uploaded ZIP independently
inspected, ZIP CRC PASS):**
`sandbox/outputs/v130-issue29-linear-DAE-e795ac44-Bo62JD0X/`,
source `e795ac4443e68973052fdfcc4c131ca53ebe1f03`,
Docker `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Host preflight: **17/17 checkpoint-format,
5/5 public API, 10/10 linear-wave tests PASS**.
The eight actual native physical runs finished:
two 10-step seeds `0→.02`, two 12-step uninterrupted
`0→.024`, and four two-step source-to-target
MPI `1→1/1→2/2→1/2→2` cold restarts
`.02→.024`, using the true production
`DPHS.save_checkpoint/load_checkpoint`.
All have **426 complete original FE DOFs**,
positive PETSc convergence reason, zero
rejected steps and zero SNES failures;
original maximum `||F_a||=5.13e-17` and
`||F_TSCN||=1.47e-15`.
Both full original seed checkpoints are SHA-256 valid
and **bit-identical** to the accepted source state;
every resumed full original FE vector at t=.02
is bit-identical to its source checkpoint.

**Observed continuous-vs-restarted final
full 426-coefficient L2 errors:**
MPI `1→1: 1.10466e-13`,
`1→2: 5.82907e-13`,
`2→1: 6.42876e-13`,
`2→2: 4.18018e-14`.
Maximum Hamiltonian jump at checkpoint
`5.42e-20` (practically zero).

**Why the first comparator failed AFTER all native
runs:** the original wave-port `Port.compute()`
returns, for every `k≥1`, a **Crank–Nicolson
mean-state STEP power** `P_k`, even though
stored at endpoint `t_k`; only `P_0`
is instantaneous initial-state power.
The original wave's control-port orientation
yields the exact discrete balance
`H_k-H_{k-1}=-Δt_k P_k`, so
physical supplied work is
`W(t_n)=-Σ_{k=1}^n Δt_k P_k`.
The *v1 diagnostic alone* incorrectly
treated these as true endpoint powers and
computed positive-sign trapezoidal
`Σ Δt(P_{k-1}+P_k)/2`.
At checkpoint, prior step
`P_{10}=-0.02393131124488`
is **not equal** to restarted first
instantaneous sample
`P_0=+0.00017376987419`,
even though all 426 coefficients are
restored bitwise. Combining the wrong
v1 trapezoidal powers invented a
`2.411e-5` split-versus-continuous
energy-work mismatch, while the
correctly reconstructed source physical
Hamiltonian and port power show:

- **Maximum native CN per-step discrete
  Hamiltonian+port-work defect
  `2.44e-19`** over all 8 original
  segments, computed independently
  from immutable NPZ arrays.
- Maximum corrected split-versus-continuous
  **total incoming boundary work difference
  below `7e-19`**.
- Maximum corrected continuous-vs-split
  **signed `ΔH-W` discrepancy
  below `2e-19`**.
- The uninterrupted original physical
  Hamiltonian increases to
  `H(.024)=3.66098987026944e-4`,
  and correctly computed incoming
  discrete port work equals this
  within `3.3e-19`.
  Thus the original reported apparent
  `dH-W=7.03835494e-4` was a
  **diagnostic quadrature/sign error**,
  *not* actual energy non-conservation.

**Source correction isolated to diagnostic/qualification:**
- `sandbox/compare_v130_checkpoint_linear_dae.py`
  explicitly accepts both original **v1
  native reports** and future **v2
  reports**. It validates the archived
  v1 trapezoid as the actual recorded
  legacy convention, then independently
  recomputes correct work using original
  stored untouched `port_power[k≥1]`
  arrays. It verifies CN Hamiltonian
  identity **at every native step**, full
  state SHA256/DAE continuity and all
  four cross-MPI work comparisons.
  It creates NEW
  `linear_dae_comparison.json` v2
  with correct per-step physics and
  detailed original provenance, **without
  overwriting/altering any original
  native segment JSON, NPZ or checkpoint**.
- `sandbox/probe_v130_checkpoint_linear_dae.py`
  uses correct CN step-integrated signed
  work and explicitly saves a **v2**
  physical report for any future run.
  No public runtime `DPHS`/
  `Port.compute` or original weak
  forms/time-step operators changed.
- `tests/test_v130_checkpoint_linear_dae.py`
  adds regressions distinguishing
  a restored initial instantaneous
  `P_0` from the previous CN step
  `P_{10}`, validating immutable
  v1 re-auditing, native-v2 semantics,
  and deliberately corrupted
  step-power/energy identity.
  **13 host tests now need rerun**.
- A separate follow-up
  **[#30](https://github.com/g-haine/scrimp/issues/30)**
  records the analogous generic
  **`DPHS.plot_powers()` trapezoid
  over stage powers** problem. This
  needs a broader, time-scheme-aware
  solution for visualizations, but
  does not justify modifying the
  physically qualified #29
  checkpoint implementation.

**MOST IMPORTANT: REUSE THE ALREADY QUALIFIED
8 NATIVE RUNS. DO NOT RERUN DOCKER.**
The only remaining immediate action is
to run the new host regressions, then
the **offline** corrected comparator
against the existing evidence directory:

```bash
cd ~/Simulations/scrimp
git pull --ff-only

python -m py_compile \
  sandbox/probe_v130_checkpoint_linear_dae.py \
  sandbox/compare_v130_checkpoint_linear_dae.py \
  tests/test_v130_checkpoint_linear_dae.py
python tests/test_v130_checkpoint_linear_dae.py

bash scripts/run-v130-checkpoint-linear-dae.sh --compare-only \
  sandbox/outputs/v130-issue29-linear-DAE-e795ac44-Bo62JD0X
```

**Qualification reporting:** actual physical
8/8 native linear-wave simulations are
accepted and independent source NPZ reanalysis
supports all original CN energy/state/MPI gates.
The patched public comparator's final
`NATIVE_LINEAR_PHDae_DPHS_COLD_RESTART_PASSED`
output and 13/13 new host tests are **NOT
YET OBSERVED**; issue #29 remains OPEN
until their execution and evidence review.
No BDF2/IMEX integrator-history replay
or general inconsistent high-index
initializer is claimed.
