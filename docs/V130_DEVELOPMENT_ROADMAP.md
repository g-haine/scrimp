# SCRIMP v1.3.0 — scope, limits, and further work

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

- [ ] **P0 #29 — native-qualified generic checkpoint/restart** before
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
