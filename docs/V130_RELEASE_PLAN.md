# SCRIMP v1.3.0 — release plan and acceptance gates

**Scope agreed 18 November 2026.** This document is the authoritative
*v1.3.0 work plan*; [V130_DEVELOPMENT_ROADMAP.md](V130_DEVELOPMENT_ROADMAP.md)
retains the historical performance diagnostics, negative experiments
and outstanding technical limitations.

**Version principle:** v1.3.0 is a backward-compatible scientific
release. No general Saint-Venant 2D tangential-velocity boundary-control
development enters its critical path. That model is a separate,
**post-release** project for LHMNC 2027.

**DAE initialization (20 November 2026):
TRUE IMPLICIT DEFAULT — FULL NATIVE MPI 1/2 PHYSICAL PASS.
The scientific and native default-promotion gate is CLOSED.**
The production linear algebraic initializer
has passed native MPI 1/2 tests for both the
conservative acoustic wave and the internally
dissipative acoustic wave. The second model
returned `||F_a||~1e-13`, zero drift of q0/p0,
`D=0.0132790026>0`,
`W=-5.40249757e-5`, maximum signed
energy-defect `3.55e-15` and full
CN MPI1/2 trajectory discrepancy
`6.37e-14` (native host:
`sandbox/outputs/v130-second-linear-dae-7eec0172-Ostu9LKD/`).
The historical `pseudo` and tiny-step BDF2
methods remain available explicitly and
retain their measured nonzero initial-data drift.

**Development-branch source candidate:** the
implicit `init_step_ts_type` fallback has
now been changed to `algebraic` in both
`DPHS.set_time_scheme()` and
`DPHS.init_step()`. Original top-level
examples, notebooks, published Maxwell code,
GUI-generated scripts and selected sandbox
simulations explicitly request historical
`pseudo` to avoid accidental physics changes.

**First native no-option trial (20 November 2026) FAILED
during PETSc option preflight, before any MPI physics test.**
The previous core code still inserted
`init_step_ts_type="pseudo"` into
`PETSc.Options` when no option was supplied,
despite having changed both the fallback string
and `DPHS.init_step()` dispatch to `algebraic`.
The native preflight also incorrectly assumed
`PETSc.Options.hasName()` proved caller
explicitness. Both were fixed at the source:
register `algebraic` in PETSc and examine
`kwargs` for explicitness. Synthetic tests
now reject the previous silent pseudo-default
bug, and native conservative/dissipative probes
verify both absence of a caller keyword and
the actual effective PETSc strategy. Original
native failed output:
`sandbox/outputs/v130-algebraic-default-h0p2-973842a3-rV7dnIEX/`.
See section 25 of
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md).
**Native rerun of the default gate remains REQUIRED.**

**Second native retry (source 3c28e919, 20 November
2026):** PETSc.Options smoke **PASS** (actual
default `algebraic`, `init_step=False` and
explicit `pseudo` compatibility). The first
conservative wave diagnostic then failed
**before solver assembly** because its
pre-run benchmark guard still expected
the inherited strategy to be `pseudo`,
whereas the corrected library now
materializes `algebraic` before the
diagnostic reconfigures time integration.
This diagnostic-only guard now clears
either of the two known inherited choices
for the no-keyword test, and a pure unit
test checks this behavior. **The 1/2-rank
conservative and dissipative default
physics gates remain NOT EXECUTED in the
second attempt.** Native evidence:
`sandbox/outputs/v130-algebraic-default-h0p2-3c28e919-WPt6Z7wP/`.
See section 26 of
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md).
No release/acceptance gate was weakened.

**Updated final native result (20 November 2026):
implicit algebraic DEFAULT PASSED on both
conservative and dissipative pH-DAE waves,
MPI 1/2.** The actual no-keyword call
selected `algebraic`; all conservative
11-state CN snapshots reproduce the frozen
independent references (scaled errors
`0` and `5.30974918e-18`), true initial
algebraic residuals are ~`1e-13`,
dynamic initial conditions remain exact,
signed Hamiltonian/port power defects are
~`6e-15`, and the original incompatible
Gaussian is rejected before CN. The
dissipative full-power identity has positive
`D=0.0132790026`, nonzero signed boundary
work and defects ~`3.6e-15`, with
MPI 1/2 trajectory parity ~`6.37e-14`.
Explicit historical pseudo and tiny BDF2
both PASS. The independent strict default
comparator ends with five PASS lines.
Frozen native host evidence:
`sandbox/outputs/v130-algebraic-default-h0p2-3a17a65b-H38miGDO/`.
See section 27 of
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md).
**The library DEFAULT is therefore qualified
for these two linear model classes**;
broader nonlinear models, memory/exceptional
MPI cleanup, all examples and packaging
remain *separate release gates*.

**Accelerated release week 23–27 November 2026:**
the current dated priorities, issue-closure
conditions, P0 Jacobian gate on the *original*
`examples/shallow_water.py`, Maxwell 3D
schedule risk and the precise go/no-go
criterion are maintained in
[V130_RELEASE_SPRINT_2026W48.md](V130_RELEASE_SPRINT_2026W48.md).
For the actual Jacobian source code,
`SCRIMP_NONLINEAR_MASS_JACOBIAN=finite_difference`
is already the default, and its **full directional
FD native gate is now PASSED** on the original
shipped shallow-water model in MPI 1/2
(`1.06105781e-8` vs `3.84378444e-6`
historical legacy; independent nonzero explicit
gyroscopic term, no tolerance relaxation).
Frozen output:
`sandbox/outputs/v130-shallow-mass-J-fae066a6-VvvbTRp4/`.
The complete native acceptance gate for issue
[#27](https://github.com/g-haine/scrimp/issues/27)
is **PASS**: the independent directional-FD
qualification above and the actual short
physical PETSc BDF2 integration on the original
shallow-water weak forms, with fully explicit
accepted-state-frozen gyroscopic `G`, both passed
on MPI 1/2. The latter is frozen at source
`9ecd173983252795a3cba13b5879d29127b57638`,
output
`sandbox/outputs/v130-shallow-explicit-TS-9ecd1739-LZjvda5C/`.
Four accepted steps per case, no rejects/SNES
failures, exact frozen-cache drift zero and full
MPI state parity `6.39215556e-16`.
No exact discrete Hamiltonian balance, global order
two, nonlinear DAE initializer or MPI rank-4
qualification is implied. Do not rerun these
frozen gates without relevant core changes.
The general nonlinear DAE Newton extension
[#28](https://github.com/g-haine/scrimp/issues/28)
is a separate, explicitly guarded future
qualification, not a reason to make untested
default assumptions.
See [technical design](V130_DAE_ALGEBRAIC_INITIALIZER_DESIGN.md),
[default migration](V130_DAE_DEFAULT_MIGRATION.md),
and [chronological evidence](V130_DAE_INITIAL_VALIDATION.md).

## 0 — Finish solver optimization and establish scientific equivalence

**Already achieved (production wave, h=0.04, 47,942 unknowns):**
three independent complete runs with 10 accepted fixed-CN steps at
1/2/4/6 MPI ranks. All current MPI final and initial states agree
within the prescribed tolerance, MUMPS parallel symbolic analysis
is effective, and native symbolic/numeric factor counts equal
exact collectively detected CSR pattern/value changes.
Current solution-time medians: **1.501 / 1.224 / 0.945 / 0.926 s**
at 1/2/4/6 ranks, respectively; efficiency at six ranks is
**27.0%**, a documented strong-scaling limitation.
The 1/2-rank native factor-reuse smoke passes with 1 symbolic,
2 numeric and 10 triangular solve operations over 10 CN steps.

**Resolved cross-version discrepancy (18 November 2026).** The
previous unqualified three-repeat campaign on the branch before the
initialization fix found scaled L2 state differences 0.376/0.397.
The cause was missing GetFEM synchronization of the accepted pseudo-TS
state after correctly excluding prescribed source bricks from
nonlinear Jacobian assembly (see diagnostic below).

**Native v1.2.0 → corrected current release comparison PASSED** on
commit `24e70e92f3db357121e6eb41b1e0e38dc6c94913`, with the
original tagged v1.2.0 commit
`88e2c1442d6aca3c7971c856e049c26d9d1caf0a`.
The user supplied the complete native Docker archive for direct
independent inspection. Fixed wave.py, h=0.04, **47,942 DOFs**,
CN dt=0.002, t_f=0.02, 10 accepted steps:

| Same-model sequential comparison | v1.2.0 native | Corrected candidate |
|---|---:|---:|
| Full wall solve including assembly | 73.267577 s | 1.512175 s |
| Setup (benchmark wrapper) | 0.095716 s | 0.096254 s |
| Mass assembly + PETSc transfer (logged) | 15.15 s | 0.0452 s |
| Stiffness assembly + PETSc transfer (logged) | 56.08 s | 0.1693 s |
| First accepted CN step | 0.39369 s | 0.38129 s |
| Mean CN steps 2–10 | 0.13113 s | 0.05357 s |
| Peak single-process RSS | 592.2 MiB | 582.6 MiB |

**One-run, state-qualified end-to-end ratio: 48.45x**. Most of
this gain is faster linear assembly/GetFEM-to-PETSc transfer, *not*
a 48x faster MUMPS solve; mean steady-CN step is ~2.45x faster.
Do not present the one-run ratio as a 3-repetition median or
extrapolate to nonlinear systems.

Independent full-vector comparison of archived NumPy snapshots:
the complete initial PETSc states are **bit-for-bit identical**;
final state L2 difference is **6.2150211e-30** (relative L2
**3.923971e-31**, maximum absolute component delta
**1.2809612e-30**). Current PETSc main CN records exactly
**1 symbolic LU, 2 numeric LU, 10 triangular solves**, matching
10 Jacobian callbacks with 2 exact numeric updates and 8 reuses.

Reference archive created by
`SCRIMP_V130_BENCH_RANKS=1 bash sandbox/run_v130_wave_benchmark.sh 0.04 1`,
in local `sandbox/outputs/v130-wave-h0p04-24e70e92-F8rBHU/`
(the output directory is deliberately ignored by Git). Both
sequential versions were rerun once from original source tags;
neither program was altered during the timed run.

**Scientific scope:** matching full endpoints validates the numerical
regression, **not yet** the continuous-in-time power/Hamiltonian
balance or fidelity of pseudo-TS initialization to the physically
prescribed `q(0)=0`. The corrected branch's MPI scaling at
1/2/4/6 ranks and full trajectories still need requalification.

**Focused source-level diagnosis (pending native confirmation).**
In v1.2.0, `set_control()` declares prescribed input source bricks
with `linear=False`. The old `DPHS.set_linear_flags()` therefore
marked the entire stiffness as nonlinear and `DPHS.monitor()`
copied every accepted pseudo-TS state back into GetFEM variables.
The current implementation correctly excludes RHS-only source bricks
when classifying Jacobian nonlinearity, so that *implicit* state
synchronization no longer occurs in the linear acoustic example.
This is a plausible explanation for the 17,640-component initial-state
difference, **not yet proof** that restoring the legacy synchronization
is physically correct.

A second change affects the residual: the current default is
`SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem` rather than v1.2.0's
legacy reconstruction from the tangent matrix. A quick A/B probe
is available, using the already archived v1.2.0 vectors; it runs only
two short **current-version** simulations and does not touch the
production solver:

~~~bash
bash sandbox/run_wave_initial_sync_probe.sh
~~~

It reports named GetFEM variable ranges, the pseudo-TS accepted
state, full initial/final vector differences against archived v1.2.0,
and `png/wave_initial_named_blocks.png`:

- current with accepted pseudo-step synchronization, direct residual;
- current with the same synchronization and legacy residual.

**Decision gate:** identify whether v1.2.0 agreement arises from
an undocumented pseudo-initialization side effect, changed residual
semantics, or both; then check initial conditions and DAE constraints
mathematically. **Do not change the production model to match an
unvalidated legacy numerical artifact.**

### 18 Nov diagnostic: root cause of cross-version state mismatch confirmed

The two current-version experimental runs with replayed pseudo-TS
monitor synchronization reproduce the archived **unmodified v1.2.0**
complete state exactly at initialization, and to ~6.22e-30 absolute
L2 at t_f, for **both** `direct_getfem` and `legacy` residual forms.
Without replay, the initial and final L2 differences were ~5.82
and ~6.29. The accepted pseudo-TS state is thus the actual source of
the state mismatch, **not** an effect of residual formulation.

The accepted pseudo-state has field norms
`q=0.00471761, p=14.36784, e_q=0.0202897, e_p=5.76279`,
while before the pseudo step `q=e_q=e_p=0`. The user-prescribed
`q(0)=0` is therefore perturbed by pseudo-time advancement; the
initialization's *physical* DAE-consistency and fidelity to prescribed
differential initial conditions still require a separate check. Do
not claim scientific validity simply from equality with v1.2.0.

**Production fix implemented:** after the *accepted* pseudo-TS solve,
unconditionally scatter its global PETSc state and transfer it to
GetFEM variables, before constructing the main TS initial vector.
This is independent of `linear_mass`/`linear_stiffness` and preserves
the correctly classified prescribed source bricks and MUMPS reuse.
A direct unit regression tests that both linear and nonlinear flag
configurations transfer the accepted state, without aliasing memory.

**Native production qualification PASSED:** the corrected
current source and untouched original v1.2.0 source were both run
once, with equal full initial/final states (numeric details above).
The optional standalone regression
`sandbox/run_wave_initial_state_regression.sh` remains useful
when both reference JSON and NPZ are available; the full clean-run
benchmark above is the authoritative regenerated reference.
The cross-release PNG/SVG comparison was generated and the 
original diagnostic mismatch is resolved.

Next: native current 1/2/4/6-rank full-vector parity after the fix,
then an independent scientific DAE/port-power check.

Exit gates:

- [x] Explain historical cross-version state mismatch:
  initial pseudo-TS state was not propagated to GetFEM after fixing
  the source/nonlinearity classification. Both experimental replays
  reproduce the native v1.2.0 complete trajectory endpoints.
- [x] Verify the production pseudo-state sync patch with a
  fresh native v1.2.0 and corrected-current sequential comparison;
  complete initial states bitwise identical and final states agree
  to ~4e-31 scaled L2; exact native MUMPS reuse audit PASS.
- [ ] Verify physical DAE consistency and whether pseudo-time
  advancement respects the intended prescribed initial q(0)=0.
- [ ] Correct any regression exposed; keep both versions' numerical
  outputs and the build/solve timing categories distinct.
- [ ] Compare full *time trajectories*, not just initial/final
  vectors, and check discrete Hamiltonian, control-port work and
  dissipative contribution when applicable.
- [x] Current 1/2/4/6 rank full-state and exact LU-reuse
  performance gate at h=0.04, three repetitions.
- [x] Re-run a fresh original-source cross-release sequential
  comparison after the initial-state synchronization fix (one
  repetition, fully matched initial/final states, ×48.45 measured).
- [ ] Verify corrected source at 1/2/4/6 MPI ranks and repeat
  qualified sequential performance before presenting a robust median;
  then qualify larger h=0.01 wave runs or document remaining limits.
- [ ] Document the additional independent candidate Jacobian's
  memory cost on a larger mesh and during long integrations.
- [ ] Inspect remaining native MUMPS/TS performance without
  introducing another solver optimization into the v1.3 scope.

**Git decision at this gate:** keep working on the existing
`feature/getfem55-mpi-container-scalability` branch until the
equivalence issue is understood. Then cut a single long-lived
`release/1.3.0-preparation` integration branch. Do not prematurely
merge an unqualified solver change to `master`.

### Final point-0 large-wave strong-scaling qualification: h=0.005 — PASSED

**19 November 2026: full native production campaign PASS, five independent
repetitions at each of 1, 2, 4, 6 MPI ranks — 20 complete simulations.**
Not a pilot; no v1.2.0 process was run. Exact source commit:
`25837de018ae53f8056f41cedf2f66f7ddb3841c`; archived native
output `sandbox/outputs/v130-large-wave-h0p005-25837de0-CBukYK/`.
The user supplied the complete 4.3-MB reduced ZIP containing **all
20 PETSc staged logs and JSON reports, 12 PNG + 12 SVG plots,
two CSV files and the original analysis**; large NPZ state vectors
were omitted from the shared ZIP after native verification.

Fixed `examples/wave.py`, mesh **h=0.005**, actual global
**2,966,786 DoFs**, ten accepted Crank–Nicolson steps,
dt=0.002, final time=0.02. Distributed GetFEM→PETSc COO,
compact matrix mode, default exact MUMPS `on_change`, ICNTL(28)
parallel ordering on multi-rank solves. CPU affinity:
0; 0,2; 0,2,4,6; 0,2,4,6,8,10. No production source changed
during the twenty measured runs.

| MPI ranks | Full solve median (s) | Min–max (s) | Full solve + build median (s) | Speedup vs current 1-rank | Efficiency | Largest rank RSS HWM (GiB) | Sum of rank RSS HWM (GiB) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 134.679 | 134.292–136.065 | 142.195 | 1.000× | 100.0% | 14.05 | 14.05 |
| 2 | 109.466 | 109.335–109.734 | 116.766 | 1.230× | 61.5% | 9.20 | 18.34 |
| **4** | **86.239** | **85.655–86.660** | **94.724** | **1.562×** | **39.0%** | **6.14** | **22.17** |
| 6 | 91.493 | 91.063–93.093 | 101.151 | 1.472× | 24.5% | 5.17 | 25.64 |

**Strong-scaling optimum is four MPI ranks on the tested six-physical-
core machine**, not six. Five-run distributions do not overlap.
From 4→6 ranks: wall solve **+6.1%**, total incl. build **+6.8%**,
and summed process RSS lifetime high-water marks **+15.7%**.
Per-rank HWM falls, however, from 6.14 to 5.17 GiB, so six ranks
may be useful for per-process memory constraints. RSS sums are **not
an instantaneous synchronized total** and include MPI replication.
No claims about other node topologies.

**Native scientific/numerical gate:** all 20 PETSc stage logs were
independently parsed from the shared ZIP and agree exactly with
SCRIMP's recorded collective matrix-update decisions:
`MatLUFactorSym=1`, `MatLUFactorNum=2`,
`MatSolve=10` per main TS. The on_change policy reused unchanged
numeric LU on eight of ten CN Jacobian calls. The native analyzer
previously compared **full** initial/final vectors in all 20 runs
and reported PASS; from the reduced ZIP its audit results can be
checked, but the large vectors cannot be independently reconstructed.
Maximum reported scaled full-vector L2 discrepancies (across all
MPI ranks and repetitions): initial **1.04e-14**, final **7.56e-13**
(tolerance 5e-8).

**Why 6 ranks is slower than 4:** the measured native stages indicate
higher startup/communication/ordering cost rather than a failure of
numeric factorization scaling. Median initial DAE consistency rises
from **23.77 to 26.40 s**, main TS from **52.48 to 55.14 s**,
symbolic MUMPS LU from **15.98 to 18.08 s**; yet numeric LU gets
slightly faster (**9.38→9.07 s**) and triangular solves improve
(**3.73→3.47 s**). The GetFEM stiffness phase also grows
(**2.93→3.55 s**) despite faster stiffness bridge
(**4.32→3.70 s**). Profiling phases overlap: **do not sum
nested event timings**. Strong-scaling saturation is documented,
not a justification to retune the solver mid-release.

The full campaign created **12 human-readable figures in `png/` and
`svg/`**: strong scaling, speedup, efficiency, five-repeat spread,
first/warm CN and every CN step, SCRIMP phases, GetFEM/COO bridge,
MUMPS factorization breakdown/counts, memory and separate
first/warm speedup. The command to reproduce (default three reps,
here explicitly FIVE):

~~~bash
bash sandbox/run_v130_large_wave_scaling.sh 0.005 5
~~~

**Conclusion:** performance/scaling and MPI full-state/numeric-LU
qualification at approximately 3M DoFs **passed**. The separate
coenergy discrete power-identity test passed on 1/2 MPI ranks
(h=0.12, 5,724 DoFs, 5 CN steps; defects ~1e-16).
This campaign does **not** establish physically faithful `q(0)=0`
in the artificial pseudo-TS initial consistency step, universal
nonlinear solver correctness, or flat RSS over long integrations.
Retain those as documented scientific/memory limitations and
dedicated follow-up tests. The v1.3.0 release itself remains
gated by the remaining workstreams below.

### 20 November — physical DAE initial-data diagnostic: native result and follow-up

The Ubuntu Docker BuildKit/native/PETSc import qualification is
closed. Scientific qualification resumes with **two separate**
initialization obligations: (1) algebraic DAE constraint
consistency, (2) fidelity to prescribed differential initial data
`q(0)=0`, `p(0)=3**(-20*((x-0.5)^2+(y-0.5)^2))`.
The historical accepted pseudo-TS changed the initial q field
(norm ~4.71761e-3 in an earlier production wave diagnosis).
Matching v1.2.0 is therefore **not** a scientific initial-data PASS.

A source-level diagnostic and reproducible, small, **rc3 image**
native runner now exist. They measure GetFEM pre/post pseudo-TS
full states, physical-field q/p drift, owned algebraic rows of the
actual PETSc `IFunction`, zero algebraic mass-row entries,
the first physical CN state and **all six saved complete state
vectors** across MPI 1/2. They never modify production solvers,
the accepted pseudo-state or any Docker native image. The
diagnostic **reports** a physical fidelity FAIL separately
rather than hiding it behind a numerical MPI PASS.

Procedure and explicit acceptance criteria:
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md).
**Run status: pending native Ubuntu evidence.** The next
step is the small h=0.2, five-step 1/2-rank diagnostic; study the
algebraic residual and dynamic-field drift before implementing
any algebraic-only projection. After a mathematically defensible
initialization, use the existing driven-wave
`check_acoustic_energy_mpi.py` /`analyze_acoustic_power_balance.py`
and coenergy physical-power tests, then Maxwell 3D. The ~3M-DoF
strong-scaling campaign need not be repeated just to perform
this initial physics diagnosis.

**New native evidence, 20 November 2026, completed.**
The unchanged published acoustic wave (`h=0.2`, 2,230 DoFs,
fixed CN 5 steps) was run in the qualified rc3 image at MPI 1
and 2. Initial q was exactly zero but drifted by coefficient
L2 **0.00114424447**; initial p drifted by **0.00751282216**
(~0.248% of its original coefficient L2 norm).
The algebraic residual decreased from **0.010689895** to
**2.44609825e-18 / 2.31888379e-18** at MPI 1/2.
Every saved full CN state agrees at scaled L2 maximum
**1.78851051e-14**. Thus MPI parity and algebraic consistency
**PASS**, and physical prescribed differential initial data
**FAIL**. Native output directory:
`sandbox/outputs/v130-dae-wave-h0p2-b56ceedb-7TyrnhWS/`.
Full interpretation:
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md),
section 5. **No production initial-condition fix is justified yet.**

**Compatibility hypothesis and next executable gate:** the left
velocity input is `U_L(0)=0`, but the prescribed Gaussian p0
has a nonzero left-boundary trace; with `e_p=p/rho`, the
constraints could be incompatible if q/p are held fixed.
This is not yet proved for the actual discretization. A
diagnostic-only finite-element feasibility test
(`scripts/run-v130-fixed-dynamics-probe.sh`) now forms the
actual PETSc algebraic-only block, solves it in least squares
while freezing every q/p coefficient, and independently
reevaluates the production `IFunction` on the candidate.
**First native result: UNRESOLVED at both 1/2 MPI ranks.**
With all q/p coefficients fixed, the remaining algebraic
residual is **4.82244319e-6** on both ranks, versus
**2.4461e-18/2.2127e-18** after the original pseudo-TS.
The extracted block's linear prediction matches the
actual PETSc `IFunction` to ~1.9e-18; the LSMR condition
estimates differ (167.7/126.4), which is not a rank
certificate. Original evidence:
`sandbox/outputs/v130-algebraic-only-h0p2-de3b2c73-RSTntEzH/`.

**Follow-up native evidence (20 November 2026): published control
PASS, alternate control NOT TESTED due to diagnostic parser error.**
The enhanced original-control 1/2-rank run confirms that
`e_p` accounts for most of the `4.82244319e-6`
fixed-dynamics residual (`e_p: 4.82044425e-6`), with
left input/output residuals ~`9.817e-8`.
LSMR stopped with normal-equation residuals ~2e-18.
The bounded dual witness establishes a numerical
`4.82243693e-6`/`4.82243747e-6` lower bound
for corrections with `||delta_a|| <= 12.09304`,
**not** a global exact rank proof.

The `matched_left_control` case was stopped during
`GetFEM.add_source_term_brick` by the diagnostic-only
Python exponent `3**...` (GetFEM generic assembly
rejects `**`). No result for the alternative physics
or its boundary compatibility can be inferred from
that aborted run. The temporary source expression has
been replaced by the equivalent `exp(...)` analytic
Gaussian trace, and the runner now executes a dedicated
**native matched-control syntax preflight** before
the four A/B cases. A source-only regression checks the
expression's numerical boundary trace. Existing rc3
Docker image and production equations are untouched.
**Actual A/B native rerun completed:** preflight passed;
both 1/2-rank `published` and `matched_left_control`
simulations passed their numerical/MPI checks. Fixed-data
algebraic-only residuals were **4.82244319e-6** and
**9.18141097e-6**, respectively. Matching the **analytic**
left-boundary trace increased the residual by ~90.4%;
neither scenario is qualified as fully consistent with
prescribed q/p. Both scenarios are accurate linear
`IFunction` reconstructions, and q/p were unchanged by
the attempted algebraic-only projection. Host evidence:
`sandbox/outputs/v130-algebraic-only-h0p2-25f9bf27-DbVTy7aw/`.
This rejects the simple continuous-trace boundary cure,
without excluding weak/discrete compatibility effects.
Full scientific interpretation: section 9 of
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md).

**Native nonzero FE-compatible positive control COMPLETED
and PASS, 20 November 2026:** with manufactured
`p0=x*(3-x)`, `e_p=x`, `q0=0`,
`U_L(0)=0`, the algebraic-only correction passes the
original PETSc `IFunction` residual at
**7.06518655e-12 / 7.13966896e-12** (MPI 1/2),
starting from **0.162506819**, while **exactly
preserving every prescribed q/p coefficient**.
The original unchanged pseudo-TS reaches ~1.8e-16
algebraic residual **but changes q by
0.00195958812 and p by 1.30196443e-5**
at both ranks, despite the existence of a fully
compatible algebraic-only state.
Host evidence:
`sandbox/outputs/v130-dae-fe-positive-h0p2-1dfbc184-MEYxdCWf/`.
The issue is therefore not simply an inability of the
algebraic projection procedure to solve the actual
weak-form constraints; pseudo-time integration is not
a fixed-differential-data DAE projection.
It is still premature to enable a production initializer
or generalize to nonlinear/high-index DAE models.
See section 11 of
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md).

**Full CN and signed Hamiltonian/port-power A/B
native qualification PASSED on both 1 and 2 MPI ranks
(20 November 2026).** With exactly the same manufactured
wave and nonzero left drive `U_L(t)=0.08*sin(20*t)`,
10 fixed CN steps, the pseudo-TS changes
`q0=0` by **0.00195958812** and `p0` by
**1.30196443e-5**. The algebraic-only initializer
changes neither prescribed state. Both methods satisfy
the non-vacuous port-Hamiltonian balance to ~
**4e-15–6e-15** (actual defects), and reproduce
all 11 full states between MPI 1 and 2 to
max scaled L2 ~**1.76e-14** (pseudo) and
**2.05e-12** (algebraic). Measured absolute
port work is **5.41131226e-5** and
**5.39545892e-5**, respectively. This establishes
the desired physical distinction on the qualified
linear, compatible acoustic example, not a
generic nonlinear/Maxwell DAE proof.
Host result:
`sandbox/outputs/v130-wave-init-power-h0p2-62e31d9b-zDa95p5J/`.
See section 13 of
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md).

**Next engineering gate: scalable PETSc-distributed
algebraic projection — SOURCE READY, NATIVE PENDING.**
The validated diagnostic used rank-0 SciPy LSMR
and ran the old pseudo-TS before replacing its
state; it is **not** production-ready. A separate
PETSc `IS`/distributed `Mat.createSubMatrix` /
`KSP LSQR` prototype now performs a
fixed-dynamic algebraic correction **without
running pseudo-time TS at all**. Its native
control cases are the same compatible nonzero
polynomial wave (require real CN/power PASS)
and the published Gaussian with incompatible
constraints (require explicit t0 rejection),
on MPI 1 and 2, frozen rc3.
**First native distributed probe failed at the PETSc API
boundary before a physical solve:** the positive
polynomial MPI=1 test constructed the FE model,
created the PETSc algebraic submatrix and returned
from LSQR, but `Vec.setValues` rejected
NumPy `int64` indices on an rc3 build whose
`PetscInt` is `int32`. No DAE feasibility
result or Hamiltonian/port balance was measured
in that abort. The diagnostic now inserts
algebraic indices as `PETSc.IntType`, checks
integer range, and includes a source-only
regression. Re-run the same launcher
`scripts/run-v130-petsc-mpi-dae-initial.sh`
after pulling; do **not** qualify PETSc LSQR
physics until native positive/negative MPI 1/2
output is received. Host failed attempt:
`sandbox/outputs/v130-petsc-dae-mpi-h0p2-62d92b6a-3urxt6Di/`.
No production `dphs.py` change, Docker
rebuild, v1.2.0 or scaling rerun.
**Second native distributed PETSc trial (20 November):
positive polynomial FAIL at the *unchanged* physical
acceptance threshold, after passing PETSc index insertion.**
True algebraic residual measured
`1.88801542e-9` against
`1.72506819e-9` threshold, with original
`0.162506819`. This is ~9.45% above threshold;
LSQR termination reason/iterations were not
emitted by the old failure path. No main
CN step was taken, no 2-rank or Gaussian
result is claimed. Failed native host evidence:
`sandbox/outputs/v130-petsc-dae-mpi-h0p2-ce6e8edf-yuonmmPp/`.

The diagnostic is now extended with **up to
three bounded distributed LSQR defect-correction
passes on the actual native F_a**, preserving
q/p bitwise and keeping the original acceptance
threshold unchanged. Every pass records true
`IFunction` residual, KSP reason, iterations,
internal reported residual and correction norm.
Stagnation/divergence terminates refinement;
the best actual state is kept and verified.
Positive rejection now archives its JSON before
returning FAIL, rather than losing the decisive
KSP diagnostics on exception.
**Next action:** rerun only
`scripts/run-v130-petsc-mpi-dae-initial.sh`
after pulling. PETSc MPI projector validation
remains **PENDING**, not PASS. No Docker
rebuild or production change.
See sections 14–17 of
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md).

**Third native fully distributed PETSc DAE run, 20 November
2026 — FOUR CASES COMPLETE, strict MPI physics gate NOT YET PASS.**
The polynomial positive case passes its
`1.7250682e-9` DAE feasibility limit at MPI 1
(`F_a=9.8141843e-14`, two LSQR passes) and
MPI 2 (`F_a=1.7163317e-9`, **one** LSQR
pass). Both published Gaussian runs are
explicitly rejected at t0 with residual
`4.8224432e-6` versus threshold
`2.0689895e-10`; no pseudo-TS was
executed anywhere. Both positive runs complete
the 10-step physical CN and independently satisfy
the non-vacuous port-power Hamiltonian gate.
The final comparator nonetheless correctly
**FAILS** the unrelaxed MPI full-trajectory
gate: `6.20613748944739e-8 > 5e-8`.
Host: `sandbox/outputs/v130-petsc-dae-mpi-h0p2-6005b09f-8kuP5pGk/`.

**New controlled diagnostic — NATIVE PENDING:**
distinguish the unchanged DAE feasibility criterion
from a *stricter optional positive-control
refinement goal* `min(tolerance,1e-11)`.
This avoids stopping MPI 2 after one only-marginally
accurate correction and may make both projected
initial efforts comparable. No feasibility or
MPI parity tolerance is relaxed. The same small
launcher `scripts/run-v130-petsc-mpi-dae-initial.sh`
now saves t0 projection discrepancy and per-CN-step
errors to `distributed_dae_parity_diagnostics.json`
even on comparator FAIL, so a residual-accuracy
explanation can be tested against a possible
nullspace ambiguity. Do not promote code to
`scrimp/dphs.py` until the full strict physical
MPI gate passes. See section 18 of
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md).

**FINAL SCIENTIFIC DAE PROTOTYPE GATE: PASSED
(20 November 2026; qualified frozen rc3 image).**

All four native PETSc distributed DAE scenarios complete:
positive nonzero polynomial initial condition
`p0=x*(3-x)` gives true algebraic
`||F_a(t0)||_2 = 9.8141843e-14` (MPI 1) and
`9.6187779e-14` (MPI 2), each in **two** distributed
LSQR passes without any pseudo-TS and without
gathering `K_aa` on rank 0. The explicit original
Gaussian negative controls both produce
`||F_a||_2 = 4.8224432e-6` and are
**EXPECTED_REJECTION** before any physical CN step.

Both positive cases preserve prescribed q0/p0
**bitwise** and execute the unchanged original
10-step Crank–Nicolson integrator. Max
Hamiltonian + signed, **nonzero** boundary work
defects are `5.6085329e-15` and
`6.0526306e-15` for MPI 1/2.
Cross-rank projected initial-state discrepancy
is `9.21579389e-13`; maximum scaled
full-trajectory MPI discrepancy is
`9.22182758e-13`, versus the original,
unchanged `5e-8` acceptance limit.
All strict native scientific comparator gates
**PASS**. Host evidence:
`sandbox/outputs/v130-petsc-dae-mpi-h0p2-4d3b4d0c-gb0qNKB8/`.
Source Git SHA:
`4d3b4d0cdee4e9b66d49170f5d128466c74e5ad8`.
See section 19 of
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md).

**NEXT: opt-in production engineering, NOT yet qualified.**
The source-ready prototype uses a distributed
PETSc algebraic submatrix, but still replicates
the full initial state using SCRIMP's existing
`Scatter.toAll` and GetFEM
`from_variables`/`to_variables`; a distributed
matrix is not an end-to-end memory-scalability
proof. Design the opt-in
`init_step_ts_type=algebraic` path while
leaving the default `pseudo` untouched.
Require exact preservation of differential
initial DOFs, true `IFunction` feasibility,
a separate strict refinement target, bounded
iterations, explicit unsupported-model guards
and atomic rejection with no fallback to
pseudo-TS. Use the qualified native
wave 1/2 MPI JSON/NPZ as independent regression
fixtures; test production lifecycle, cleanup,
memory overhead and larger/other linear models
before broadening to nonlinear DAE or Maxwell.
**No SCRIMP production solver code was modified
by this scientific validation.**

**Next bounded engineering deliverable (separate from the already PASS
scientific gate):** add a production-owned, opt-in
`init_step_ts_type=algebraic` mode for certified linear DAE layouts.
The integration must use the *existing* pre-main-TS initialization
hook but bypass pseudo-time TS construction entirely; retain
`init_step_ts_type=pseudo` as the default. Implement the distributed
PETSc algebraic extraction and bounded LSQR correction in a
testable module rather than importing `sandbox` diagnostics.
Add guards for the complete state/port partition, algebraic
mass rows, unsupported lumped/substituted and nonlinear
models, and safe handling of incompatible data with native
residual, original-state restoration and collective MPI
cleanup. Keep the physical feasibility threshold separate
from an independent projection accuracy target.

**Acceptance for this NEW core change:**
verify legacy default pseudo and any supported BDF
initialization still behave identically; check
exact initial dynamic coefficients, native
IFunction constraints and a positive/negative
1/2-rank MPI test; reproduce **all** CN states,
signed boundary powers and Hamiltonian residuals
against the independent frozen diagnostic `rc3`
evidence. Then measure initialization wall time,
per-rank RSS and temporary PETSc object lifetimes,
and extend to a second independently compatible
linear model. Do not claim general nonlinear/high-index
or fully memory-scalable initialization until those
separate tests pass.

The method, current scope and completed scientific
measurements are documented in
[V130_DAE_ALGEBRAIC_INITIALIZER_DESIGN.md](V130_DAE_ALGEBRAIC_INITIALIZER_DESIGN.md).

**CORE OPT-IN PHYSICAL QUALIFICATION: PASSED
(20 November 2026).** The production entry point
`init_step_ts_type="algebraic"` uses the distributed
PETSc helper without pseudo-time TS and reproduces
all eleven independent frozen acoustic CN states
to a max scaled error **0** (MPI 1) /
**5.30974918e-18** (MPI 2). Prescribed q0/p0
are bitwise unchanged; true algebraic residuals
are ~**1e-13**; signed boundary-power energy
defects ~**6e-15**. Published Gaussian is rejected
on MPI 1/2 at t0 with
`4.82244319e-6` against `2.0689895e-10`.
The explicit legacy `pseudo` physical trajectory
matches the frozen reference exactly.
The negative-time BDF2 `init_step_dt=1e-6`
run has nonzero drifts
`1.9595912757338825e-5` (q) and
`8.680599894760661e-10` (p), despite
a similarly accurate physical power balance.
Native output:
`sandbox/outputs/v130-core-dae-h0p2-66fe02c2-WykmP41Y/`.

**SECOND MODEL NATIVE STATUS (20 November 2026):
MPI-1 PROJECTION AND CN PASS; COMPLETE PHYSICAL POWER/MPI
QUALIFICATION STILL PENDING.**
The original `examples/wave_dissipative.py` model
and its physical damping bricks were constructed on
`h=0.2`, 2,804 DOFs. The production algebraic
initializer reached true
`||F_a||=9.99416399e-14` in two passes
and the unchanged main CN TS accepted 10/10 steps
with no rejections or nonlinear failures.
Power diagnostics then failed at GetFEM
`gf.asm` argument 5 because the interior
damping port's `get_region()` is `None`.
The existing SCRIMP `Port.compute()` convention
maps whole-domain ports to `region=-1`;
the independent audit is now corrected and a
source regression was added.
**The actual signed Hamiltonian + boundary work +
dissipation balance and MPI 1/2 parity have
NOT YET BEEN MEASURED**, so this gate does
not yet PASS. Failed native evidence:
`sandbox/outputs/v130-second-linear-dae-29aa01dd-ULacuhnd/`.

**NEXT NATIVE RETRY:**

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-second-linear-dae.sh
~~~

Eight historical top-level examples have been
migrated **source-only** to explicit
`init_step_ts_type="pseudo"` without changing
their physics, initial data or time steps.
Global default promotion to `algebraic`
will follow only after the second model,
unsupported-model audit, cleanup/memory and
native **no-explicit-initializer** regression.
No silent pseudo fallback; no default switch
yet. See
[V130_DAE_DEFAULT_MIGRATION.md](V130_DAE_DEFAULT_MIGRATION.md).

## 1 — GetFEM 5.5 simplicial FEM contributions (SCRIMP-first)

**Scope decision (19 November 2026):** support only **triangles and
tetrahedra**. No quadrilateral, hexahedral or prism developments in v1.3.0.
SCRIMP users are the immediate audience. Our patches may later be handed
to the GetFEM maintainers **as a useful starting point**, not as a
community-ready, exhaustive and fully certified FEM contribution.
The independent source-level and native findings are in
[V130_GETFEM_FEM_AUDIT.md](V130_GETFEM_FEM_AUDIT.md).

- Reuse existing GetFEM 5.5 RT elements and lowest-order Nédélec.
- Fix demonstrably failing BDM tetrahedral constructors first; preserve
  existing working triangle cases. Build source patches live under
  `scripts/patches/getfem-5.5/` and are applied by
  `scripts/build-getfem-mpi.sh`. Do not modify the upstream source tree
  inside SCRIMP's repository or change point-0 production solvers.
- Add Nédélec first-kind and second-kind families on **simplices**;
  document order conventions and the distinction between local
  polynomial dimension and *globally conforming* DoFs. **Release-gated
  orders are k=1,2 only** (plus the pre-existing Nédélec-I k=0).
  Nédélec-I is coded parametrically but k>=3 remains unvalidated;
  **Nédélec-II is implemented only for k=1,2 and rejects k>=3**.
  Arbitrary-order support is still an enhancement, not a completed
  deliverable. The per-entity degree-of-freedom formulas and future
  tests are recorded in
  [V130_GETFEM_HIGHER_ORDER_TODO.md](V130_GETFEM_HIGHER_ORDER_TODO.md).
  Do not advertise unverified H(curl) elements.
- Keep checks **lean and purposeful**: constructors and expected local
  dimensions (initially k<=2, plus one higher-order smoke), a two-cell
  orientation/normal-or-tangent trace check, and a small assembly test.
  Preserve existing SCRIMP regression tests. Validate 1/2-rank MPI
  integration when the Maxwell workstream starts; broader convergence
  and energy/power tests belong to the actual Maxwell example (point 3),
  not an exhaustive independent GetFEM upstream certification campaign.
- Build the patched source in a **new test Docker image**, inspect the
  report before replacing the qualified point-0 image, and record the
  image ID and source patch identity. Packaging/infrastructure refinement
  remains point 2.

**Current evidence (19 November 2026):** the *unmodified*
Docker GetFEM 5.5 MPI-2 image constructed 12/14 descriptors; both
`FEM_BDMK(3,1)` and `FEM_BDMK(3,2)` failed at the upstream
`ipoint == nddl` assertion. **BDM tetrahedron corrections are now
validated for SCRIMP's lean k=1,2 scope** in disposable Docker:
the modified Python extension loads the new core, 12/30 local DoFs,
21/54 shared-face global DoFs, positive mass matrices and four
physical normal-trace comparisons pass (max relative mismatch
4.529e-6 at 1e-7 offsets). The approved `0001` BDM and `0002`
Python-link source patches remain in `scripts/patches/getfem-5.5/`.
The Nédélec-I `0003` patch is **now accepted for SCRIMP's
lean k=1,2 scope** and has been moved to the approved top-level patch
series. Its disposable native smoke passed on triangles/tetrahedra
(k=1/2): expected two-cell connectivity, SPD mass matrices and
tangential trace differences at most 2.385e-6 relative.
Nédélec-II `0004` has now **passed** the same native targeted
checks for k=1,2 on triangles/tetrahedra (6/12 local DoFs in 2D,
12/30 local DoFs in 3D; two-cell shared DoFs 10/21 in 2D and
18/48 in 3D; SPD mass matrices; relative tangential mismatch
at most 1.446e-6). It has been promoted to the approved top-level
patch sequence, **with the explicit k=1,2 constructor restriction**.
A general `k>=3` Nédélec-II implementation requires new face and
interior functionals and is not included. No claim of all-order
family support is made. 
The old-host-conda 5/14 run is
non-qualifying and must not be cited as a GetFEM 5.5 regression.

Exit for the lean GetFEM gate: SCRIMP's selected simplex families at
the explicitly accepted k=1,2 orders can be constructed and used in
small representative sequential tests with no known interface
misorientation, and the GetFEM source changes remain explicit,
version-pinned and available as self-contained patches. Upstream
acceptance is **not** an SCRIMP v1.3.0 release gate.

## 2 — Docker, native install and platform entry points

**Implementation in progress (19 November 2026).** The BuildKit
candidate Dockerfile, safe candidate build driver, and lean 1/2-rank
native smoke are now prepared. The **first build attempt** reached
GetFEM after a 1.60 MB Docker context transfer and successful
Conda/PETSc/SLEPc ABI checks, but stopped at a root-owned BuildKit
tarball cache. Cache ownership was repaired and the unintended CUDA
PETSc/SLEPc variant was excluded. A **second build confirmed**
writable caches, CPU complex PETSc/SLEPc 3.26.1, NumPy-2 Python
bindings and Gmsh imports, then stopped at a **GetFEM 5.5 release
archive SHA-256 discrepancy**. The Savannah 5.5 release archive was
subsequently compared to immutable Git tag `v5.5`:
1,170 common files, two explainable distribution differences
(`CMakeLists.txt` version metadata, `gf_util.cc` equivalent MUMPS
availability macro), six omitted files unused by our Autotools
build, and **all five approved source patches apply cleanly**.
The reviewed Savannah digest
`0e042770b86152b09be4bfd0ee160b8d10d49c3f8a0b9b06d55f1ffaa1d30fc4`
is now **strictly pinned**.

**Candidate build and first smoke PASS (19 Nov):** image
`scrimp-mpi:v1.3.0-rc1`, revision
`96961d647c2caf41dd77d2a7b5926e86e86c4b69`,
manifest digest
`sha256:ab08d9beaa60a1685b7e8df1e69f5a085f5571ca46c2b7d8554f53ec42ace3e8`.
The GetFEM 5.5 native build completed; MPI parallelism level 2,
complex PETSc/SLEPc 3.26.1, NumPy 2.5.5 and distributed MUMPS
solve passed on 1/2 ranks. RT, BDM and Nédélec I/II supported
simplex constructor checks passed. A newly scaffolded **user-owned
project**, mounted read-write, successfully ran 10 time steps of
a four-DoF port-Hamiltonian oscillator and wrote
`outputs/oscillator.csv` **to the host** (11 time points,
including Hamiltonian); its 2-rank MPI startup also passed.

The sample produced a nonfatal NumPy `ComplexWarning` in real
GetFEM state synchronization. A small source-only fix now selects
the real part explicitly, consistent with the existing TS monitor;
the candidate image above predates this change. The PETSc
`-info` unused-option warning requires separate triage.
**VS Code GUI attempt (19 Nov):** Docker now works without
`sudo` and the candidate image imports SCRIMP correctly as
`mambauser`. Dev Containers initially failed *before startup*
because the starter's colon-bearing `${localEnv:...:fallback}`
resolved to `scrimp-mpi` and attempted to pull nonexistent
`scrimp-mpi:latest`. The template now embeds
`"image": "scrimp-mpi:v1.3.0-rc1"` explicitly; generated-project
tests enforce a full image tag. Previously generated independent
projects require a one-line local update.

**VS Code Dev Containers Linux execution PASSED (19 Nov):**
The user opened the SCRIMP project in a real VS Code Dev Container
as `mambauser` with the `(scrimp)` Conda environment active
and mounted `/work`. In its integrated terminal,
`python simulations/first_simulation.py` completed 10 accepted
steps (0 rejected, 0 nonlinear failures), saved 11 rows to
`/work/outputs/oscillator.csv`, and
`mpiexec -n 2 python simulations/mpi_check.py` printed
rank 1/2, rank 2/2, and PASS.

A **separate editor discovery issue** was observed: on first
entry VS Code offered **"No Python found — Install uv and Python"**,
despite Python already existing and working at
`/opt/conda/envs/scrimp/bin/python`. The user should **dismiss
the offer**, not install a second Python; select the existing
interpreter with **Python: Select Interpreter** if necessary.
This is now covered in the [VS Code tutorial](source/vscode.rst)
and starter README. The template already declares the exact
Python path. The PETSc `-info` warning and the older-image
`ComplexWarning` are nonfatal and separately tracked.

**VS Code Python debugging PASS (20 Nov):** an Ubuntu Dev
Container session used the image's exact interpreter
`/opt/conda/envs/scrimp/bin/python` with the `debugpy`
launcher. F5 started the first SCRIMP simulation under the
debugger, an active breakpoint at `system = S.DPHS("real")`
was registered, repeated **F10 / Step Over** advanced execution
through the initialization, and VS Code stopped at line 31,
`p = S.State("p", "Momentum", "lumped")`. The screenshot
shows the live `q` and `system` objects in **Locals** and
the `main()` **Call Stack**. This qualifies interactive
serial Python breakpoints, stepping and variable inspection,
**not** multiprocess MPI debugging or completion after resuming.

**F5 resume-to-completion PASS (20 Nov):** the developer
resumed the paused and stepped-through `first_simulation.py`
with F5. The numerical run finished normally inside VS Code.

**VS Code interpreter and Pylance issues — resolved on Ubuntu
(20 Nov):** the independent project initially displayed
`No Python found / Install uv and Python`, followed by
`An Invalid Python interpreter is selected`, despite the
working `/opt/conda/envs/scrimp/bin/python`. The Python
extension log confirmed the path as active but reported
`No Python environment is set for this resource`; the
new `ms-python.vscode-python-envs` manager was in charge.

The **confirmed solution on the user's existing Dev Container**
is `"python.useEnvironmentsExtension": false` in
`customizations.vscode.settings` of its own
`.devcontainer/devcontainer.json`. The user reported **no
remaining interpreter warning** after adding this missing
line. The independent `python.analysis.extraPaths` setting
(`/workspace/scrimp`) had already eliminated the false Pylance
`import scrimp` diagnostic. The model template in the GitHub
branch contains both working settings and retains the explicit
`python.defaultInterpreterPath`. The earlier trial
`python-envs.globalSearchPaths` was insufficient and has
been removed from the template; it is unnecessary when the
new environment manager is disabled.

**Ubuntu VS Code first-user workflow: PASS as reported** for
opening the Dev Container, running a SCRIMP simulation with
host-persistent output, launching 2 MPI ranks, using F5/F10
breakpoints and variables, resuming to completion, resolving
editable imports in Pylance and selecting the preinstalled
Python interpreter without errors after applying the setting.
The full acceptance matrix, image provenance, exact Dev Container
settings and limitations are recorded in
[V130_UBUNTU_VSCODE_ACCEPTANCE.md](V130_UBUNTU_VSCODE_ACCEPTANCE.md).
A further fresh-project / fresh-user onboarding acceptance is
desirable, but this developer's observed blocking issue is
closed. This does not validate multi-process MPI debugging.

**20 November 2026 candidate audit:** the safe exact-SHA rc2 rebuild,
first/warm BuildKit evidence, fatal NumPy ComplexWarning regression,
numerically checked host CSV and PETSc `-info` diagnosis are specified
in [V130_DOCKER_BUILD.md](V130_DOCKER_BUILD.md), section 4a.
These are *new qualification gates*, not completed runtime results.
Source CI failure at `170f26d1` was traced to an obsolete GetFEM
timing helper path (not a native or physical failure); the workflow
and test were corrected and require confirmation by a new green run.
The Ubuntu/VS Code acceptance stays qualified only for the previously
tested image and editor workflow. Preserve `scrimp-mpi:local`.

**20 November Ubuntu Docker rc2 native/cache acceptance — PASS:**
two BuildKit builds at fixed commit `ae96cb7334d7c88cd2e0dbd4d9615af66f1f4a61`
staged 190 tracked files (1.54 MiB; transfer 1.62 MB), reused
all Conda, PETSc/SLEPc and GetFEM 5.5 native layers on the first
build and all Dockerfile layers on the warm rebuild. Native GetFEM
MPI level 2, complex PETSc/SLEPc, RT/BDM/Nédélec constrained
constructors, MUMPS at 1/2 ranks, 10-step independent oscillator,
strict ComplexWarning, 11-sample persistent Hamiltonian CSV and
2-rank project MPI startup **all pass**. Exact evidence and
attestation-vs-runtime digest distinction:
[V130_DOCKER_RC2_ACCEPTANCE.md](V130_DOCKER_RC2_ACCEPTANCE.md).
The two exports have identical main image manifests and configs
but different BuildKit attestation/manifests-list digests; this
does **not** certify clean-room reproducibility.

**20 November rc3 targeted native qualification — PASS:**
`scrimp-mpi:v1.3.0-rc3` built on Ubuntu from immutable
`1d6d0c90115213298047b28237a7f38087563a35`,
with only source-layer rebuilding and full native GetFEM/PETSc
cache hits. The strict `ComplexWarning` oscillator, 1/2-rank
distributed MUMPS, RT/BDM/Nédélec constructors, host-persistent
Hamiltonian CSV and independently scaffolded two-rank MPI startup
passed. The **PETSc `-info` source fix is confirmed at runtime**:
neither rank's import-only test contained that option. The
independent before/after `scrimp-mpi:local` image-ID `cmp`
completed without a difference. See
[V130_DOCKER_RC3_ACCEPTANCE.md](V130_DOCKER_RC3_ACCEPTANCE.md).

**Final rc3 PETSc warning-free diagnostic — PASS (20 November 2026).**
The corrected host-side smoke now runs a temporary mounted Python
file rather than `python -c`. The developer supplied both two-rank
results (`PASS: PETSc import options clean (rank 1/2)` and
`(rank 2/2)`), followed by `PASS: PETSc -info and -c absent,
no unused-option warning.` This closes the import-time PETSc
diagnostic without rebuilding rc3 or any native layers. The
**Ubuntu Docker native/BuildKit and independent starter smoke
milestone is closed**, combining the previously observed complete
rc3 smoke and this final diagnostic excerpt; see
[V130_DOCKER_RC3_ACCEPTANCE.md](V130_DOCKER_RC3_ACCEPTANCE.md).

**Still pending:** reproducible clean-room build with pinned
base-image digest/platform and tested explicit Conda lock;
fresh-user onboarding and cross-platform installation; scientific
full-trajectory/DAE consistency, discrete Hamiltonian/port-power
identities, Maxwell 3D and long-time memory/performance release
qualification. The current Python package metadata remains
`scrimp-1.2.0`; update it during the deliberate v1.3.0 packaging
freeze.
No GHCR image has been published or promoted.
The reproducible build procedure, cache boundaries, evidence paths
and deferred lock/platform work are documented in
[V130_DOCKER_BUILD.md](V130_DOCKER_BUILD.md). Following the
**>45 GB accidental Docker context transfer**, the supported build
now fetches GitHub source by default and creates a **physically isolated,
small, Git-tracked allowlisted context**. A developer may explicitly
select local source, but even then user outputs never enter the build.
No image has been promoted and `scrimp-mpi:local` remains untouched.

A **self-contained end-user project template** is now staged at
`project-template/` with a generator, future release ZIP packager,
VS Code Dev Containers configuration, Python debugging and MPI
startup task. A short port-Hamiltonian oscillator exports a persistent
host CSV. **Both host-mounted CLI and a real Ubuntu VS Code Dev
Container integrated-terminal run have passed**, including
the simulation and 2-rank MPI startup. **Serial Python F5/debugpy
breakpoints, F10 stepping, inspection of DPHS/State objects
and F5 resume-to-completion have now all passed**.
The interpreter-manager opt-out and Pylance paths have now
both been validated in the developer's Ubuntu VS Code project.
Clean-room onboarding, cross-platform behavior and multi-rank
debugging remain pending.
See the [VS Code onboarding guide](source/vscode.rst).

- Rework the Dockerfile into cache-friendly dependency layers:
  package download and native PETSc/SLEPc/GetFEM build layers
  **before** copying frequent SCRIMP changes.
- Provide BuildKit cache mounts for apt/micromamba and source
  tarballs, pinned hashes and a resumable/reproducible download
  strategy. Test fail/restart without unnecessary redownload.
- Keep Docker build context minimal; no historical artifacts,
  generated mesh/figures, large ignored outputs or local
  simulation archives should enter it.
- Ship an installer/launcher for **Ubuntu/Debian**, **Windows**
  (document Docker Desktop/WSL2 requirements) and **macOS**
  (including Apple Silicon image constraints). Separate
  supported/validated platforms from best-effort instructions.
- Ship the validated **independent VS Code project template** as a
  small GitHub release ZIP. Its Dev Container must select the actually
  published GHCR image (or an explicitly rebuilt local image);
  installed SCRIMP and user simulations must remain separate.
  Continue to offer native Conda + Python interpreter selection.
- Rebuild the image **locally on the user's computer**, prove
  SCRIMP, GetFEM, petsc4py/slepc4py and MPI smoke tests, then
  remove obsolete image(s) as requested. Avoid deleting the only
  functional image before its replacement passes.
- Provide version-pinned image metadata, licensing, examples and
  volume/permission handling.

Exit: a reproducible clean image, repeatable warm-cache build
and documented cross-platform entry points. This covers #25.

## 3 — New 3D Maxwell structure-preserving example

- Develop a scientifically meaningful 3D electromagnetic
  port-Hamiltonian/PFEM formulation with compatible **H(curl)**
  and **H(div)** spaces, appropriately oriented Nédélec and
  RT/BDM functions, and explicit boundary power/Poynting flux.
- Use a simple geometry with analytically or independently
  checkable initialization, boundary conditions and power identity.
- Establish convergence, time integration, material coefficients,
  energy evolution and discrete power-balance residual.
- Test at 1/2/4 (and, when useful, 6) MPI ranks and at multiple
  mesh and polynomial orders; qualify GetFEM element orientation.
- Generate ParaView-compatible results and explanatory visuals.

Exit: reproducible example, accurate port-power identity, working
MPI scaling and educational narrative. No unsupported claim that
spatial discretization is structure preserving solely because
its Hamiltonian is positive.

## 4 — Performance campaign and documentation figures

- Compare *equivalent* v1.2.0 serial and v1.3 serial wave cases
  only after gate 0; separately report current 1/2/4/6 MPI scaling.
- Repeat comparable 3D Maxwell workloads. Scale by mesh resolution,
  DOF count, number of timesteps and polynomial degree.
- Report setup, GetFEM assembly, GetFEM→PETSc bridge, initial
  consistency, main CN step times, symbolic/numeric LU,
  triangular solve, maximum RSS (per process **and total**),
  and parallel efficiency. Profile first step vs steady steps
  without summing inclusive timers.
- At least three repetitions in quiet host conditions.
- Publish accurate PNG and SVG plots under the documented output
  directory and selected *qualified* figures in Sphinx.

Exit: scientifically justified, reproducible performance results
with complete raw machine-readable data; known 4→6 rank scaling
limits displayed rather than concealed.

## 5 — Rebuild Sphinx documentation, preserving valid mathematics

- Edit/correct retained scientific derivations and tutorials; remove
  obsolete attempts and dead links, not useful mathematical content.
- Explain the complete default SCRIMP stack: PFEM/GetFEM,
  distributed COO bridge, PETSc TS/SNES/KSP, LU/MUMPS,
  parallel symbolic analysis, exact `on_change` numeric and
  symbolic reuse, initialization and memory/parallel limits.
- Document precisely which geometries, elements, nonlinear
  models, MPI operations, boundary conditions and power
  balances have actually passed validation.
- Provide a **VS Code Dev Containers** tutorial: host Docker
  installation, bind mounts, persistent outputs, MPI terminals,
  Python interpreter/debugging, UID/GID permission model
  and CPU/NUMA configuration.
  A first version is already drafted in `docs/source/vscode.rst`,
  with source-level project template tests. End-to-end Linux and
  cross-platform validation, screenshots and the released GHCR image
  remain pending.
- Ensure HTML, LaTeX/PDF, navigation, examples and screenshots
  render correctly. Maintain the user-level and developer-level
  reading paths.

Exit: reviewed Sphinx documentation, green HTML/PDF build and
reproducible container walkthrough. This covers #26.

## 6 — All examples / optional installation self-check

- Inventory **every** shipped example, then run each with
  its documented inputs in serial and at `mpiexec -n 4`
  (where the model's MPI contract makes sense).
- Check expected physical invariants, convergence/failure codes,
  power balance, short trajectory parity and controlled outputs.
- Fix/mark unsupported models explicitly; do not silently
  proclaim success when an example only imports.
- Supply a bounded opt-in CLI, e.g. `scrimp --check-install`
  or `python -m scrimp --check-install`, with quick serial
  and MPI smoke subsets, dependency information, machine-readable
  summary and meaningful exit code. Extended test suite is
  separate so installation checks remain practical.

Exit: complete example matrix and exercised installation check
on the validated platforms.

## 7 — GUI review / script prototyping role

- Audit which current API calls and options the GUI generates,
  then correct incompatibilities introduced by v1.3.0.
- Demonstrate GUI→Python SCRIMP script generation using at least
  one basic example and a clear round-trip reproducibility check.
- Document the GUI explicitly as a **prototyping aid for SCRIMP
  scripts**, not as a universal PDE editor or replacement for
  advanced scripted MPI/model development.
- Test startup, forms, file dialogs, generated scripts and
  failure messages; document supported desktops, Qt/display
  requirements and known restrictions.

Exit: working GUI smoke, explicit scope and updated screenshots/docs.

## 8 — Clean **tracked** sandbox, not local user outputs

- Preserve the documented `wave_1D` sample as the only
  **tracked sandbox example**.
- Move useful benchmark/research regressions to `tests/`,
  `benchmarks/` or `scripts/`, and update imports and CI.
  Do **not** delete an active test merely to make sandbox empty.
- Drop redundant experiments from the tracked tree after recording
  their findings in documentation and Git history.
- Never delete a user's local ignored `sandbox/outputs/`, mesh,
  animation, archives or benchmarks.

Exit: `git ls-files sandbox` contains only the intended documented
sample (and any narrowly required package/build metadata, if
explicitly justified), with all scientific tests functioning
elsewhere.

## 9 — Prepare and publish GitHub v1.3.0 release

- Run all unit, example, scientific, MPI and GUI gates; verify
  package headers/licenses, clean install, Python dependencies,
  wheels/source distribution, Docker provenance and Sphinx.
- Audit and close only genuinely resolved issues; reference
  fixes and regressions in release notes, retain open issues
  with documented limitations.
- Freeze packaging/docs versions only after full qualification.
- Merge the single reviewed v1.3.0 integration branch into
  `master`; tag `v1.3.0`, publish GitHub release
  (plus container artifact, when qualified), provide
  reproducible install + changelog instructions.

Exit: tagged and reproducible release, source and documentation
versions consistent, passing CI and scientifically validated demos.

## Git and work-stream strategy

- **Now:** finish (0) on
  `feature/getfem55-mpi-container-scalability`.
- **Once gate 0 passes:** create *one*
  `release/1.3.0-preparation` integration branch.
  Keep implementation commits logically separated; no branch
  per checklist item.
- **GetFEM upstream code:** use a separate branch/PR in the
  GetFEM fork because it is a different repository. Pin the
  exact patch SHA in the SCRIMP integration branch.
- **Before release:** review and merge only once into
  `master`. No premature release tag.
- Existing issues #6, #23, #25 and #26 are tracked individually.
  #23 is a distinct GetFEM MPI-2 root-cause investigation,
  not automatically solved by implementing new element families.
  Do not close #6 before its full physical/MPI tests pass.
- The tangential-velocity-controlled Saint-Venant 2D study for
  **LHMNC 2027 starts only after completing all gates 0–9.**

## Nonlinear shallow-water gyroscopic brick: true explicit RHS (20 November 2026)

**User decision:** keep `examples/shallow_water.py` gyroscopic
`G` **fully explicit**. The previous
`explicit=True` code only excluded `G` from
the implicit Jacobian while still reevaluating its
contribution at the current Newton trial state;
this broke the actual `IFunction`-`IJacobian`
derivative contract, producing a zero-`zd`
relative Jacobian directional-FD error of
`2.02366589e-5` on both MPI 1/2 ranks.
Frozen first full native output:

```text
sandbox/outputs/v130-shallow-mass-J-c9bd39e4-3B7wkKE3/
```

The **production** `explicit=True` semantics
have now been corrected:
`TS.setPreStep()` assembles the original
nonlinear weak brick at the **last accepted**
`(t_n,z_n)`; the resulting calibrated
distributed PETSc RHS is **frozen across all
Newton and line-search iterations** of the
following step. The term is excluded from
the current implicit RHS assembly and Jacobian.
The original (G) expression, sign and
`explicit=True` option remain unchanged.
This is a **lagged first-order explicit term**,
not a discrete-energy-preserving implicit
midpoint gyroscopic term. Unsupported
`explicit=True` forms are rejected.

**Acceptance status:** source-only isolated
tests and CI exercise the full accepted-state
freeze contract; **native FD/physical
qualification is pending**. The strict
`scripts/run-v130-shallow-mass-jacobian.sh`
does a corrected-default/legacy comparison
at **nonzero frozen explicit G** on MPI 1/2
with the **original 1e-5 threshold** and
requires zero cache drift across FD trials.
A separate short real PETSc time integration
and signed Hamiltonian/port-work audit must
verify the new split before closing
[#27](https://github.com/g-haine/scrimp/issues/27).
See [release sprint](V130_RELEASE_SPRINT_2026W48.md).
Do not confuse this with the already-qualified
linear algebraic DAE default or claim that
state-dependent mass Jacobian is universally
correct until the new native results pass.

## 20 November 2026 — shipped nonlinear shallow-water Jacobian native gate PASS

The issue #27 *directional Jacobian* acceptance sub-gate
is **CLOSED/PASS**, from the user's original frozen native
2D inviscid shallow-water model, MPI 1 and 2, at source
`fae066a64b649131351e686ac9383ea11c13c0a6`.
Full production `IFunction` versus `IJacobian`
directional FD with independent time derivative,
GetFEM state-dependent mass, explicit legacy A/B
and genuinely frozen **nonzero** gyroscopic forcing:

| MPI | Default maximum relative FD error | Legacy | Explicit G RHS norm |
| ---: | ---: | ---: | ---: |
| 1 | 1.06105781e-8 | 3.84378444e-6 | 0.00536124096 |
| 2 | 1.06105781e-8 | 3.84378444e-6 | 0.00536124096 |

Both meet the original `1e-5` limit; default
improves ~362×; the `explicit=True` cache
was unchanged throughout FD/Newton probes,
and the original gyroscopic weak form/sign
were preserved. Native output:
`sandbox/outputs/v130-shallow-mass-J-fae066a6-VvvbTRp4/`.

**IMPORTANT:** this experiment deliberately did
not execute a physical PETSc TS solve
(`main_TS_executed=False`). The remaining
original #27 acceptance gate is a separate
bounded genuine nonlinear BDF2 simulation
with accepted-state explicit G refresh,
nonvacuous manufactured initial momentum,
MPI full-state parity and measured
(not presumed zero) Hamiltonian drift.
`scripts/run-v130-shallow-explicit-physical-ts.sh`
is ready for native use, but as yet
**NOT NATIVE EXECUTED**. Keep #27 open
until the physical smoke passes.
General nonlinear DAE consistent
initialization remains separate #28.

## 20 November 2026 — shipped SWE physical TS/explicit G native PASS; issue #27 qualified

Following the previously recorded directional-FD PASS at
`fae066a64b649131351e686ac9383ea11c13c0a6`,
the independent **real PETSc TS** runner
`scripts/run-v130-shallow-explicit-physical-ts.sh`
passed on the same unchanged
`scrimp-mpi:v1.3.0-rc3` image (OCI source
`1d6d0c90115213298047b28237a7f38087563a35`),
using current source
`9ecd173983252795a3cba13b5879d29127b57638`.

Host-owned immutable evidence:
`sandbox/outputs/v130-shallow-explicit-TS-9ecd1739-LZjvda5C/`
(`manifest.json`, four per-case PETSc logs and JSON trajectory records,
`compare.log`, `shallow_explicit_TS_comparison.json`).
Original inviscid shallow-water weak bricks/signs are unchanged.
Only mesh, short horizon and, for one independent control, initial
momentum data are overridden. Physical main TS uses BDF2
with fixed `dt=1e-4`, four accepted steps, `t_f=4e-4`.
The historical pseudo initializer is requested explicitly;
**no general nonlinear DAE initial-data consistency qualification** is inferred.

| Initial momentum | MPI ranks | Accepted / rejected / SNES failures | Max explicit G RHS norm | Signed final Hamiltonian change |
| --- | ---: | ---: | ---: | ---: |
| Published `[0,0]` | 1 | 4 / 0 / 0 | 2.4381721e-6 | -0.00206566661 |
| Published `[0,0]` | 2 | 4 / 0 / 0 | 2.4381721e-6 | -0.00206566662 |
| Manufactured `[y,-x]` | 1 | 4 / 0 / 0 | 0.000688335575 | -0.00206572620 |
| Manufactured `[y,-x]` | 2 | 4 / 0 / 0 | 0.000688335575 | -0.00206572615 |

All main `TSConvergedReason=1`; the explicit RHS cache is rebuilt
at each accepted state and has **exactly zero recorded drift**
within Newton residual evaluations.
Maximum scaled full-state MPI 1/2 discrepancy:
`6.392155562495291e-16`; maximum relative Hamiltonian discrepancy:
`1.2677184090553238e-15`.
Relative absolute final Hamiltonian changes are
`6.5438023e-8` (published) and `6.54399075e-8` (manufactured).
They are **measured changes**, **not** evidence of exact pH power
preservation, nor do they isolate error from the lagged `G`.

**Release decision:** #27 is qualified in its tested native scope
and can be closed with this pinned evidence.
Keep broader MPI/port-power and 1/2/4-rank viscous SWE as
independent issue #6; upstream GetFEM assembly/partition diagnosis
as #23; nonlinear algebraic initial-data projection as #28.
A fixed-horizon `dt` refinement plus signed boundary/dissipation
power accounting is a separate scientific follow-up.

## 20 November — issue #6 native harness ready, **not yet physically qualified**

After completing and closing #27, a distinct issue #6 qualification
campaign was prepared on `release/1.3.0-preparation`:
source commits `81900e5012b293cfc679fed2f995f7630c0c0659`,
`ceeff736525a126feca563aaa59be5f0b17496d9`,
and `8a20d551760cb1e187f7e6911e485138ffb78fbd`.
The *Python source checks* workflow passed on `8a20d551`
(including negative controls for missing rank evidence, rank power
mismatch, state mismatch, nonpositive sampled h and invalid
Hamiltonian/physical power defects). **No native MPI scientific solve
has yet passed this new gate.**

`scripts/run-v130-issue6-viscous-mpi.sh` retains immutable
`scrimp-mpi:v1.3.0-rc3`, overlays only current Python
modules and independent diagnostics, enforces per-case timeouts,
and writes per (case, rank-count) work and output directories.
Two explicitly different modes:
`rest` (stationary no-feedback control, good for detecting
MPI deadlocks but scientifically vacuous as a dissipation gate);
`physical` (rotating original-form viscous Saint-Venant,
both open and closed loop, genuine nonzero-power requirement).
Both check MPI 1/2/4. Each MPI rank writes its own generic
algebraic-port-power JSON, while rank 0 saves full endpoint
state vectors and h-weighted boundary/viscous power histories.
`sandbox/compare_v130_issue6_viscous_mpi.py`
checks rank-local and cross-MPI parity, positive saved height
DOF coefficients, nonnegative measured dissipations,
true signed physical power and volume balances, and finite
PETSc completed time; it **does not** equate the generic
unweighted port powers with the model's h-weighted physics.

Native commands, **not run here**:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue6-viscous-mpi.sh rest
# Only once rest succeeds:
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue6-viscous-mpi.sh physical
```

These are bounded 0.02 s / `dt=0.002`, `h=0.2`
**diagnostic smokes**, not published 250 s reproduction.
The numerical thresholds are provisional short-smoke gates,
not permission to relabel failed physical tests as PASSED.
If a nonzero-rotation run fails its initialization, first
inspect native logs and the physical facet-tangential boundary
compatibility; report a real FAIL rather than replace it
with a stationary benchmark or weaken the criteria.
Preserve all host-owned `sandbox/outputs`.
Full issue #6 (including 4-rank powers/outputs,
no races, longer-time pH energy) remains OPEN.
GetFEM raw MPI assembly root cause remains separate #23.

## 20 November — rest MPI 1/2/4 actual native results: physics PASS, comparator false FAIL corrected

Original immutable host archive at source `10c5d22cd6e214b2ef9729cbc125a20769a4fab7`:
`sandbox/outputs/v130-issue6-viscous-rest-10c5d22c-6KsOMm1G/`.
All three native PETSc solves actually **completed** with 10 fixed
accepted steps to `t=0.02`, zero rejected and SNES failures, and
equal per-rank generic-port-power traces. Full-state end-point
MPI 1 vs MPI 2/4 scaled L2 differences
`3.58948e-17`/`3.64345e-17`.
Relative signed Hamiltonian-balance defects over `H0`:
MPI 1 `7.96408e-16`, MPI 2 `1.47904e-15`, MPI 4
`2.61677e-15`; volume balances `<1.6e-15` relative.
Height FE coefficients remain `~5`.
The original comparator **incorrectly** applied a 5%-of-actual
energy-exchange limit to a zero-activity rest state, where
work `~1e-26`, dissipation `~1e-33` and `~1e-14`
Hamiltonian variation are all roundoff-level.
The reported 31%/58%/100% exchange ratios were therefore
**numerically ill-conditioned, not physical energy failures**.

The rest gate now has an **explicit `--expect-rest`** branch with
near-zero physical activity, height equilibrium, strict absolute
Hamiltonian-balance and MPI state/port-power equality, while
the **nonzero-flow `--require-activity` branch keeps the
original 5% relative-to-exchange acceptance**. No original
physical equations or solver settings were changed.
The archived rest native evidence can be passed to the
new comparator without re-executing the physical simulations.
The rotating open/closed MPI 1/2/4 physical gate remains
unexecuted; issue #6 remains OPEN.

**Comparator detail (20 Nov):** `relative_balance_exchange` is also
excluded from MPI 1/2/4 **parity comparisons in explicit rest mode**,
because the archive shows `0.31/0.58/1.0` for this ill-conditioned
quantity. All physical arrays, full endpoint states, volume, absolute
Hamiltonian/port-power balance, and individual rank-port traces are
**still compared**. In nontrivial `physical` mode the relative-to-exchange
metric remains mandatory for both acceptance and MPI parity.

## 20 November — actual rotating SWE failed at nonlinear BDF2 initialization (MPI 1)

The first physical #6 campaign, source `50b8e41be0beae483df81ecf9fb2ef68ce9a80d5`,
immutable `scrimp-mpi:v1.3.0-rc3` image, failed **before main TS integration**,
on the first **open-loop MPI=1** case. Original archive:
`sandbox/outputs/v130-issue6-viscous-physical-50b8e41b-mLLBUGBV/`.
The full-state norm rose from `110.087965` to `59994902.3`
(amplification ~544,976) during the fixed BDF2 consistency step
`t=-1e-9` to `0`, and the unchanged guard (`max factor=100`)
correctly rejected this invalid initial state. No physical
energy/port-power or MPI 2/4 result was generated by this campaign.
This does **not** challenge the completed mass-Jacobian #27.

Independent geometric analysis of the **archived native `Disk.msh`**
with 225 triangles and 33 straight boundary segments confirms
the circular initial momentum `p0/rho=0.2[y,-x]/(r+1e-6)`
does **not** have exactly the prescribed polygon-facet trace
`U_wall=0.2 Rot*Normal`. At segment endpoints,
`max |p0/rho-U_wall|=0.0190327569` (9.52% of 0.2);
the maximum unwanted normal component is `0.0190111897`.
At mid-facets the discrepancy is only about `2e-7`.
This is confirmed **geometric incompatibility**, not yet a
demonstrated complete explanation for the enormous BDF2 norm.
Zero initial algebraic efforts/control variables are another
plausible source of an initial DAE impulse.

New **diagnostic-only** runner on the frozen Docker environment:
`scripts/run-v130-issue6-init-probe.sh`.
It enables existing `SCRIMP_INIT_DYNAMIC_STATE_AUDIT=1` and
`--newton-block-trace` and preserves the original physical weak
forms, nonlinear mass Jacobian, 100× norm guard and MPI calibration.
All probes use one rank to isolate initialization independently
of MPI; every probe has a distinct preserved host-owned directory.
Modes: `baseline` (original), `dt1e-7` (only initialization
step modified), `rigid-facet` (smooth rotating momentum, physical
facet boundary), `warm-ep` (diagnostic algebraic guess),
`unphysical-matched` (`Rot*X` polygon-wall trace, deliberately
**unphysical**), `from-rest` (different initial-value experiment,
not the original rotating p0). Neither a successful altered-data
probe nor a forcefully skipped initializer qualifies #6.

Begin with the original `baseline` to capture block norms and
the exact h/p changes, then use `dt1e-7` and `rigid-facet`
to discriminate inverse-step amplification from geometric
initial/boundary incompatibility. Do not increase the guard
or silently switch to `--no-init-step`. The nonlinear algebraic
DAE initializer (#28) remains separately unqualified.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue6-init-probe.sh baseline
```

## 20 November 2026 — priority change: #28 before further nontrivial #6 MPI work

The native baseline diagnostic at source `e5d4e6e32b5c912eccd15786987ab1ab79fcb7f2`
(`sandbox/outputs/v130-issue6-initial-baseline-e5d4e6e3-KXzesS0i/`)
confirms a defective rotating viscous-SWE **initial consistency** step,
not a main time-stepping or MPI deadlock: the original negative-time BDF2
initializer at `dt=1e-9` produces `||Y||=5.99949e7`, while the physical
height stays at 5 and the prescribed initial momentum changes by
`||p_after-p0||/||p0||=0.760984859`. SCRIMP correctly rejects the
state with its factor-100 amplification guard.

**Decision:** issue [#28](https://github.com/g-haine/scrimp/issues/28)
is now the prerequisite to resuming the **nontrivial rotating** viscous
Saint-Venant qualification under [#6](https://github.com/g-haine/scrimp/issues/6).
Do not launch additional MPI physical/energy campaigns until
fixed-differential initial feasibility and DAE index are understood.
The already validated stationary MPI 1/2/4 gate and nonlinear
mass-Jacobian #27 are unaffected.

**Mathematical scope:** for `z=(z_d,z_a)` with prescribed `z_d(t0)=z_d0`,
the candidate constrained solve uses
`G(z_a)=F_a(t0,z_d0,z_a,0)=0` with the *current* native Jacobian
`J_aa = dF_a/dz_a`. This projection is not a general-purpose cure:
for boundary velocity constraints with Lagrange multipliers the
semi-discrete system can be higher-index or rank deficient; hidden
differentiated constraints may prevent a suitable algebraic-only
projection. The faceted-wall prescribed velocity and the original
circular rotating initial momentum are not exactly compatible;
the solver must reject mathematically inconsistent initial data
instead of altering `p0`, boundary conditions or Hamiltonian.

**Ordered acceptance gates:**

1. Structural rank/consistency audit at **unaltered** original
   `h0,p0`: zero algebraic mass rows, index-1 eligibility, correct
   algebraic DOF partition, rank/conditioning of `J_aa`, residual
   after constitutive algebraic guesses, and boundary constraint
   compatibility. Explicitly report if original data are unsatisfiable.
2. Develop a **separate** nonlinear MPI-safe, safeguarded Newton
   initializer for a demonstrably regular supported DAE class,
   using production residuals/Jacobians and a fixed differential
   initial vector; validate directional FD, descent, rollback and
   collective rejection. **Do not** simply remove the linearity
   guards on the existing native-qualified linear initializer.
3. Native positive and negative nonlinear tests, MPI 1/2, exact
   bitwise preservation of differential coefficients, nontrivial
   physical port work and no silent fallback to pseudo/BDF.
4. Only then resume #6 with an explicitly compatible physical
   rotating initial condition, original weak forms, and MPI 1/2/4
   h-weighted boundary power and viscous dissipation audits.

**Release policy:** this dependency is for the problematic nontrivial
viscous-SWE qualification. It is **not** a claim that #28 already
works or that a full generic nonlinear initializer must enter v1.3.0.
If it cannot be natively qualified within release scope, retain
the linear-only algebraic default and explicitly reject unsupported
nonlinear initializations; do not weaken the scientific gate to
meet a date.

## 20 November 2026 — issue #28.0 structural audit: source ready, NATIVE PENDING

Source-only, independently reviewable diagnostic:
`sandbox/probe_v130_issue28_structure.py` and
`scripts/run-v130-issue28-structure.sh`.
It uses the exact published experimental viscous-SWE weak bricks,
physical `p0=tangent`, faceted wall and unchanged `h0=5`.
The runner pins `scrimp-mpi:v1.3.0-rc3` and overlays only Python
modules; no Docker rebuild and no pseudo/BDF or physical time step.

Mathematical audit: dynamic state intervals `h,p`, algebraic intervals
`e_h,e_p,Y_n,U_n,Y,U`, exact DOF coverage, nonlinear tangent mass
algebraic rows, initial native `F_a`, actual GetFEM stiffness
`K(z0)=dF(t0,z0,zdot=0)/dz`, its `K_aa` algebraic block and
`K_{a,Y}` / `K_{d,Y}` subblocks, plus a central directional
`IFunction` difference when perturbing `Y` alone. Original dynamic
coefficients are checked bitwise unchanged after the probe.

**Static prediction from the open-loop weak forms:** the boundary
traction multiplier `Y` appears in `h*Y.Test_p` (dynamic momentum
row), but in no algebraic weak-form row when `gamma=0`.
Thus `dF_a/dY` is an identically zero block and `J_aa` is
structurally singular. A naive Newton that solves only `F_a=0`
cannot determine `Y`, regardless of the initial algebraic guess.
This does **not** by itself prove differentiation index 2;
differentiated boundary constraints and the coupled mass system
must still be examined before designing a supported initializer.

Runner requires **native MPI 1 and 2** and compares pinned
diagnostics, but **has NOT yet been executed on the user host**.
Any failure is diagnostic evidence, not permission to change the
original PDE or relax the norm guard. The qualified linear default
initializer and the closed #27 Jacobian are unaffected.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-structure.sh
```

## 20 November 2026 — #28.0 true native MPI 1/2 structural audit PASSED

The original single- and two-rank GetFEM/PETSc assembly audit
`scripts/run-v130-issue28-structure.sh` has now passed natively at
`bc0b42846a1fd86c58952def9a0c9163c9bf7307`, unchanged image
`scrimp-mpi:v1.3.0-rc3` (OCI revision
`1d6d0c90115213298047b28237a7f38087563a35`).
Host-owned evidence:
`sandbox/outputs/v130-issue28-structure-bc0b4284-UkXgPXar/`.
The user-visible native output reports:

```text
AUDIT MPI=1: original h0,p0, mass rows, Jacobian Y and IFunction FD
ISSUE28 STRUCTURE PASS: Y absent from Fa, present in Fd
AUDIT MPI=2: original h0,p0, mass rows, Jacobian Y and IFunction FD
ISSUE28 STRUCTURE PASS: Y absent from Fa, present in Fd
PASS: original SWE structural MPI 1/2 parity; no time step/initialization performed
```

The actual production-assembled operator exhibits a traction `Y`
**null algebraic column** and **nonzero dynamic coupling**.
The gate checked zero mass in algebraic rows, independent
`IFunction` directional finite differences in the traction direction,
bitwise preservation of prescribed `h0,p0` and MPI parity.
It did **not** execute a nonlinear DAE initializer or claim an
index classification or compatible rotating initial data.
No numerical block norms or detailed matrix extrema were printed
in the shared console; these remain in the per-rank native
`result_1/structure.json` and `result_2/structure.json` records.

### Implication: the naive nonlinear algebraic Newton is insufficient

The weak forms of the **open-loop** model contain the following
algebraic conditions (signs follow the stored bricks and SCRIMP
`position=flow` convention):

- Interior velocity constitutive relation:
  `F_{e_p} = <h*(p/rho-e_p), Test_e_p> = 0`.
- Boundary traction-test relation:
  `F_Y = <h*(e_p-U), Test_Y>_Gamma = 0`.
- Prescribed wall velocity:
  `F_U = <-U + U_n*Normal + omega*f(t)*Rot*Normal, Test_U>_Gamma = 0`,
  with normal control `U_n=0` for `gamma=0`.

Eliminating `e_p,U,U_n` therefore yields a **weak boundary constraint
on the prescribed initial differential variable**
`p0/rho|Gamma = omega*Rot*Normal` for the faceted polygon boundary.
The initial circular rotating field agrees with the ideal circular
tangent but does not exactly satisfy this faceted trace in general.
This mismatch is independently plausible from the geometry
measurements, but the rank/size of the actual **discrete weak
compatibility defect has not yet been measured**.

The traction `Y` is a multiplier associated with this boundary
velocity constraint: it appears in the dynamic momentum equation
through `B_Y^T Y`, but not in any open-loop algebraic residual.
Hence `J_aa=dF_a/dz_a` has a zero `Y` column and cannot be inverted
as proposed for a regular index-1 nonlinear DAE. This alone
**does not establish index 2**, nor prove existence/uniqueness
of a differentiated-constraint solution.

**Next gate #28.1 (not yet implemented/qualified):** extract the
native algebraic-block norms and full intervals from both
`structure.json` files, assess algebraic constitutive projection
and the remaining boundary weak residual at fixed `h0,p0`.
Only if those constraints are feasible should the team attempt
a differentiated constraint (schematically
`C_q(q,t)*dot(q)+C_t(q,t)=0`) coupled to the dynamic mass
and multiplier traction. If prescribed data are incompatible,
fail explicitly rather than alter `p0` or prescribe an
nonphysical polygon-wall trace. Keep #28 and nontrivial #6 OPEN.

## 20 November 2026 — #28.1 fixed-h,p compatibility diagnostic: SOURCE ONLY / NATIVE PENDING

After the fully native #28.0 MPI 1/2 PASS the actual per-rank
`structure.json` files show the following identical structure:
3,300 total DOFs, 1,452 differential (`h=484`, `p=968`),
1,848 algebraic and 132 traction `Y` DOFs.
`max|M_a|=0`, `max|dF_a/dY|=0`,
`max|dF_d/dY|=0.5069655658773655`.
The independent directional FD for the full dynamic `Y` block
has norm `0.5009893370080248`, absolute Jacobian discrepancy
`5.003707553108401e-17`; algebraic FD sensitivity is exactly zero.
Original `||F_a||=0.2500610865473142` for MPI1
(`0.25006108654731424` on MPI2) with zero initial algebraic
guesses. Per-block initial `||F_a||`:
`e_h=0.0880456422545`, `e_p=0.1682909140783`,
`U=0.1626550955977`, and `Y,Y_n,U_n=0`.
**These initial-zero-guess residuals do not establish boundary
incompatibility.** In particular `F_Y=0` with `e_p=U=0` is vacuous.

New source-only, no-state-acceptance feasibility audit:
`sandbox/probe_v130_issue28_feasibility.py` and
`scripts/run-v130-issue28-feasibility.sh`.
At the unchanged original `h0,p0` it builds the true
`F_a(z0)` and `J=dF/dz` using production GetFEM/PETSc,
then removes ONLY the 132 *known null* traction `Y` columns
from the otherwise algebraic unknown set. The resulting
`1848 x 1716` rectangular native algebraic Jacobian
defines a **diagnostic** distributed PETSc LSQR least-squares
problem for the 1716 remaining algebraic unknowns.
The candidate is strictly read-only: `h0,p0` and original
`Y` coefficients are checked bitwise unchanged, and every
GetFEM unknown is restored even on error.

The audit recomputes the true `IFunction` at the projected
algebraic candidate, reports **all six** algebraic block
residual norms, and checks actual `F_a(z0+delta)` against
`F_a(z0)+J_{a,*}delta` (the original open-loop algebraic
forms are affine in algebraic unknowns for fixed `h,p`).
It reports LSQR convergence, normal-equation residual
`||J^T F_a||`, and distinguishes:
`FEASIBLE_AT_TOLERANCE`,
`INCOMPATIBLE_AT_TOLERANCE` (only with independently
certified LSQR and normal-equation stationarity), and
`INCONCLUSIVE_LSQR`. A mathematical incompatibility
claim must **not** be inferred from an unconverged
or rank-deficient numerical attempt alone.

The native harness compares MPI 1/2 on
**both** the original rotation `omega=0.2` and
a separate zero-rotation `omega=0` positive control.
Both use the same physical weak forms and polygon wall;
the zero-rotation case is **not** the original rotating
experiment. No BDF2 pseudo-step or physical TS is launched;
no nonlinear Newton, higher-index method or v1.3.0
release gate has been qualified by the source-only commit.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-feasibility.sh
```

## 20 November 2026 — proposed correction via *geometric* P2, not altered wall velocity

After stage #28.1 fixed-differential-data LSQR, the native dataset
`sandbox/outputs/v130-issue28-feasibility-9082bbee-1vi0TtA9/`
reported original omega=0.2 minimal algebraic residual
`||F_a||~0.02656138`, equally for MPI 1/2, with
the dominant e_p residual `~0.02646368`.
The separate stationary control had a far smaller
LSQR residual `5e-9--7e-9` but incorrectly exceeded the
very strict `~9.47e-10` threshold; its label
`INCOMPATIBLE_AT_TOLERANCE` must **NOT** be taken
as proof that a resting equilibrium does not exist.

**Key distinction:** the experiment had field FEM P2 but
GEOMETRIC order 1 (straight polygon edges). On the actual
polygon, a circular tangent initial `p0` has a small
normal trace; simply replacing `U=omega Rot*Normal`
with `U=omega Rot*X` would match some initial traces
but introduce **nonzero wall penetration**, so it is not
an acceptable physical correction.

An **opt-in** `geometry_order=2` parameter now exists on
`scrimp.utils.mesh.Disk` and in the viscous-SWE sandbox.
Only when selected, the Gmsh mesh is promoted with
`model.mesh.setOrder(2)` *after* mesh generation/refinement.
The existing `geometry_order=1` behavior is unchanged by
default; the PDE, boundary velocity law `U=omega Rot*Normal`,
FEM field orders, parameter values, pressure/viscosity
and initial rotating `h0,p0` are UNMODIFIED.
The quadratic isoparametric boundary approximates the
ideal circle, but is not an exact circle: **do not**
assume exact DAE consistency from geometry order alone.

New bounded native diagnostic:
`scripts/run-v130-issue28-P2-geometry.sh`.
It pins the same immutable rc3 Docker image, mounts the
changed mesh generator and sandbox modules read-only, and
runs **P2 geometric order**, separately for stationary and
original rotating p0 (both MPI 1 and MPI 2).
It reuses the native fixed-h0,p0 LSQR feasibility probe;
the reports now show the requested order, imported GetFEM
geometric point counts / element counts / transformations,
true original `IFunction` residuals, and an **informational**
comparison to P1 original `0.02656138`.
Full physics, historical BDF2 pseudo initialization and
main TS steps are **not** executed here.

The runner checks MPI equality and explicit state rollback.
It does **not** assert that P2 cures the initial-state
incompatibility; its primary qualification is whether
GetFEM really retains curved P2 geometry and what
the resulting native constraint residual is.
If the original tangential initial data remain incompatible,
we must diagnose and explain the specific discrete constraint
and advise a compatible p0 or a smooth startup from rest.
No automatic alteration of user-prescribed initial data.

Source commits:
`7875617c8e6ba6193f0cc14ffe9ec9e6594bcebc`
(opt-in curved mesh) and
`a9ebc1766a3814c2215681ae74af7385ece0aa19`
(native P2 feasibility runner and compatibility suggestions).
**Native P2 evidence remains PENDING.**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-P2-geometry.sh
```

## 20 November 2026 — corrected P2 geometry report NameError; follow-up queued

The first native order-2 attempt at source `40303b1` produced a
**genuine Gmsh P2 mesh** (225 type-9 six-node triangles and 33 type-8
three-node boundary edges). Native GetFEM reported 485 geometry points
for 225 elements (versus 131 P1 points); direct inspection of the
archived `Disk.msh` shows all 33 mid-edge boundary nodes lie on the
unit circle, ~0.004528 units away from their straight chord midpoints.
The test stopped on MPI 1 at zero rotation **after assembly**, when
constructing the JSON report: `_audit()` incorrectly referenced the
CLI-local name `args.geometry_order` and raised `NameError`.
This is a probe bookkeeping bug, NOT an observed GetFEM/PETSc physical
failure or a P2 compatibility conclusion. Source corrected to use
`system._issue28_geometry_order` and a dependency-free AST regression
suite now tests that the reporting function cannot reference the
out-of-scope CLI `args`. No residual from this aborted P2 run is valid.

**Sequencing explicitly agreed:** finish and interpret the P1/P2
geometric compatibility experiment **first**. **AFTER that**, return
to #28 with the original `examples/shallow_water.py` as a distinct
nonlinear-initialization test **without the additional curved-boundary
complication**. Do not launch that example or infer its suitability
before the present native geometry gate is resolved. Review whether
its algebraic subsystem is structurally regular before choosing a
Newton method; hold its prescribed dynamic state fixed. This deferred
task is recorded here and in issue #28 so it is not forgotten.

## 20 November 2026 — correct tangent sign on Gmsh P1/P2, GetFEM check PENDING

The user asked whether the initial rotation and imposed wall tangent
have accidentally opposite orientations. **Static source reading says
NO:** `Rot=[[0,1],[-1,0]]`, `Tangent=Rot*Normal`, and
`p0/rho=omega*(y,-x)/(sqrt(x*x+y*y)+1e-6)` all rotate clockwise.
For an outward unit radial normal at (1,0), both point down.

Independently computed from the native user archives
`v130-issue28-feasibility-9082bbee-1vi0TtA9.zip` and
`v130-issue28-P2-geometry-7286bec2-Rp90btp3.zip`:
3-point / 2-point P2/P1 curved/straight Gmsh edge map,
12-point Gaussian integration per boundary element,
outward radial-choice geometric normal from the boundary
parametric derivative, all **33 original boundary edges**.
For omega=.2 original circular `p0`, the relative boundary
`L2` norm `||p0/rho-omega*Rot*n||/||p0/rho||` is:

- **P1 geometry:** `0.0550176860263` (~5.50%).
- **P2 geometry:** `0.0000696599723462` (~0.00697%).
- For the **opposite tangent sign**, the corresponding relative
  mismatch is nearly **2.0 (200%)** on both meshes.

Thus curved Gmsh P2 reduces the **geometric** initial
tangential-trace mismatch by ~790x, whereas the previous
GetFEM native algebraic LSQR residual increased from
`0.0265613800268` (P1) to `0.0277336388365` (P2).
The latter is **not explained by a reversed analytic tangent**
or simply by the straight-sided polygon approximation.
Note: Gmsh-based geometry normals are **not** the same thing as
measuring the actual GetFEM GWFL `Normal` used in the weak bricks.

New read-only native sign diagnostic
`sandbox/probe_v130_issue28_normal_orientation.py` +
`scripts/run-v130-issue28-native-normal.sh`:
uses the **actual GetFEM generic boundary integral** on the
original model, with original rotating p0 and the production
`Tangent` macro. Directly measures both signs,
`(p/rho).Normal`, `X.Normal` and the imported P1/P2
geometry in MPI 1/2. No LSQR or TS time stepping.
It can isolate a GetFEM mesh-region normal sign/assembly issue
which is not visible from the exported Gmsh coordinates.
The native GetFEM normal test is **PENDING**; do not claim
that the Gmsh geometry-only numbers prove the sign seen
by the actual weak-form assembler.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-native-normal.sh
```

**Sequence retained:** complete this boundary/GWFL normal sign
and geometry investigation first; *afterwards* use original
`examples/shallow_water.py` for a separate nonlinear fixed-dynamic
initialization study, as requested. No changes to model physics,
initial p0, existing LSQR thresholds, or BDF2 safety guard.

## 20 November 2026 — ROOT CAUSE: GetFEM GWFL column-wise Rot reversed wall actuation

**Native evidence PASS, before correction:** user uploaded
`v130-issue28-native-normal-b903865b-91WeMgdE.zip`
(original `b903865bee829d18a329928fa04dd68e0b536eb7`,
unchanged image `scrimp-mpi:v1.3.0-rc3`).
Native GetFEM integration using **the actual production
`Normal`, `Tangent=Rot*Normal`, original p0**, MPI 1 and 2:

| Geometry | `||p/rho - omega*Tangent||/||p/rho||` | `||p/rho + omega*Tangent||/||p/rho||` | `||p/rho.Normal||/||p/rho||` |
|---|---:|---:|---:|
| P1 straight | 1.9992482297 | 0.0549180824 | 0.0548807074 |
| P2 curved | 2.0000023663 | 0.0000697293 | 0.0000696876 |

The GetFEM outward normal itself is **correct**:
`integral(X.Normal)` = `6.2452910639` (P1) and
`6.2831681228` (P2), near the expected `2*pi`.
MPI1/2 are equal to floating-point precision.
This **overturns** the previous assumption that the old
`Tangent` macro matched `p0` inside GetFEM. The geometric
P1/P2 edge-normal analysis was not sufficient to determine
the sign of the literal `Rot` as interpreted by GWFL.

**Exact root cause:** GetFEM GWFL writes a matrix
`[[a,b],[c,d]]` with the nested vectors as **COLUMNS**.
Thus the previous `Rot = [[0,1],[-1,0]]` represents the
**actual** linear operator `[[0,-1],[1,0]]`, i.e. CCW.
Original `p0/rho = omega*(y,-x)/(r+1e-6)` is CW.
The wall law therefore had the opposite tangential sign,
regardless of P1/P2 geometric accuracy. This explains why
curved P2 did **not** reduce the previous algebraic residual.
Reference: https://www.getfem.org/userdoc/gasm_high.html
(explicit matrices section).

**Corrected source commit:** `6db7a79cfe6f7b20d854c4b45913922f10fda0a9`.
The experimental viscous-SWE `Rot` macro is now written
`[[0,-1],[1,0]]` (GWFL columns), meaning the intended
row-wise clockwise matrix `[[0,1],[-1,0]]`.
Original p0, h0, PDE weak forms, viscous terms, normal
impermeability and all SCRIMP solver internals are unchanged.
The nonlinear gyroscopic form
`Gyro(p) = div(Rot*p) * Rot` is invariant under a **global
`Rot -> -Rot` sign flip**: it contains two factors of Rot.
The correction therefore flips the actual wall tangential
velocity but leaves the gyroscopic physics invariant.

A source-only test `tests/test_v130_issue28_rot_convention.py`
checks the actual column-wise matrix representation, clockwise
tangent at `(1,0)`, unchanged `Gyro`/`Curl2D` macro formulas
and formal sign invariance; CI runs it. The original
`scripts/run-v130-issue28-native-normal.sh` now rejects a
`REVERSED` or `AMBIGUOUS` native sign instead of printing
a warning and exiting 0.

**Old LSQR numbers are obsolete as physical evidence:**
the pre-fix `||F_a||_min` ~0.02656 (P1) and ~0.02773 (P2)
were assembled with a *wrong-sign wall velocity*.
Repeat fixed-state P2 feasibility with the correctly signed
wall before concluding anything about residual boundary
incompatibility, and retain the separate stationary control
and MPI 1/2 check. Even if fixed-state `F_a` becomes small,
`Y` remains a null column of `J_aa` in the open-loop system;
a complete index-aware DAE initializer is still unqualified.

**Ordered native experiments:**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-native-normal.sh
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-P2-geometry.sh
```

**Do not combine workstreams:** After the native
geometry/initial-boundary compatibility issue is interpreted,
and ONLY THEN, investigate nonlinear constrained
initialization using `examples/shallow_water.py`, as
explicitly requested. No native post-fix data yet.


## 20 November 2026 — post-fix GetFEM GWFL matrix convention: native P1/P2 evidence

**Important Python user pitfall:** the GetFEM Generic Weak-Form Language
(GWFL) writes nested explicit matrices as *column vectors*, whereas
Python/NumPy list-of-lists conventions use *rows*.

- **Authoritative upstream reference:**
  [GetFEM GWFL — Explicit matrices](https://getfem.org/userdoc/gasm_high.html#explicit-matrices).
- **User-facing Sphinx guide:** `docs/source/getfem_gwfl.rst`
  (reachable from the website index, getting started, initialization
  page, and repository README).
- **Incorrect for clockwise rotation:** GWFL `[[0,1],[-1,0]]`
  (counterclockwise as an actual matrix).
- **Correct clockwise:** GWFL `[[0,-1],[1,0]]`
  or the more readable row/semicolon form `[0,1;-1,0]`.
- Actual GWFL `Normal` is outward; in the unit disk
  `Rot*Normal` now agrees with the original
  prescribed clockwise velocity `p0/rho`.

**Native post-correction qualification of sign and reduced algebraic
feasibility (not of complete DAE initialization).** Branch source
`120c01bd8aa206ea0986ba876d65d35489a825b6`, immutable
Docker `scrimp-mpi:v1.3.0-rc3`, true GetFEM/PETSc MPI 1 and 2:

| Quantity | P1 geometric mesh | P2 curved geometric mesh |
|---|---:|---:|
| Actual GWFL tangent orientation | ALIGNED (MPI 1/2) | ALIGNED (MPI 1/2) |
| Relative aligned trace mismatch | `0.0549180824` | `6.97292855e-5` |
| Relative normal leakage | `0.0548807074` | `6.96876238e-5` |
| GWFL integral `X.Normal` | `6.24529106` | `6.28316812` |

The P2 geometric transformation was actually imported as
`GT_PK(2,2)`, not simply simulated with P2 state FEM on straight edges.

**Independent original P2 algebraic LSQR**, fixed `h0,p0`:
before the GWFL sign correction (source `7286bec`)
`||F_a||_min=0.0277336388`; after the correction (source
`120c01bd`) `4.1850871648724585e-8` (MPI 1) and
`4.183695287884114e-8` (MPI 2): a reduction of
**more than 660,000 times**. The post-fix rotating feasibility tolerance
is `2.60963115021321e-9`, so the strict label is still
`INCOMPATIBLE_AT_TOLERANCE` and **must not be advertised as PASS**.
The rotating per-block projected `e_p` norm is
`4.101738589874305e-8` on MPI 1.

At rest (P2), the original known equilibrium yielded a
diagnostic-only LSQR residual of `4.746567828624845e-9`
(MPI 1) and `4.721221465777258e-9` (MPI 2), with a stricter
acceptance tolerance `9.507859931505774e-10`.
**A fail against this numerical stopping threshold is not a proof
of physical incompatibility**, particularly for a resting
state already qualified as an equilibrium. Precision/conditioning
and subsequent high-index initialization need separate treatment.

Evidence: user-supplied archive `v130-issue28-.zip`,
the host native post-fix sign output directory
`sandbox/outputs/v130-issue28-native-normal-120c01bd-sq5HgESS/`,
and the native post-fix feasibility output
`sandbox/outputs/v130-issue28-P2-geometry-120c01bd-sKH7lgKN/`.

**Do not close #28 yet:** although the tangential sign error is fixed,
the open-loop traction `Y` remains absent from the algebraic equations
(`dF_a/dY=0`), so a complete consistent initial traction generally
needs a different, possibly constraint-differentiated formulation.
As agreed with the user, test `examples/shallow_water.py` as a
separate nonlinear-initializer experiment **only after** the current
boundary/geometric work is concluded.

## 20 November 2026 — #28 resting LSQR false-negative: exact-reference native gate PENDING

The corrected original curved-P2 fixed-state LSQR audit returned
`||F_a||=4.746567828624845e-9` (MPI 1) at `omega=0`, exceeding
its own strict numerical threshold `9.507859931505774e-10`,
even though a true resting equilibrium is known to exist.
An LSQR stop at relative tolerance does **not** certify a positive
best-residual lower bound in an ill-conditioned algebraic projection.
Do NOT interpret `INCOMPATIBLE_AT_TOLERANCE` as mathematical
infeasibility without an independent reference/certificate.

New independent native reference (source only, NOT YET RUN):
`sandbox/probe_v130_issue28_exact_rest.py` with
`scripts/run-v130-issue28-exact-rest.sh`.
This benchmark does **not use LSQR**. It holds the ORIGINAL
`h0=5,p0=0`, original physical wall controls and all weak forms.
It sets only algebraic unknowns to the exact expected stationary
values `e_h=rho*g*h0=0.5`, `e_p=0`, `U_n=U=Y=0`,
`Y_n=-e_h=-0.5`, reflecting the SCRIMP flow-brick sign
convention (`M_n[h]` registers `-h*Y_n*Test_Y_n`
and `-B_n[h]^T` registers `-h*e_h*Test_Y_n`).
Native `IFunction(t0,z,zdot=0)` checks **both**
`||F_a||` and `||F_d|| <= 1e-11` on P1 and P2 geometric
disks, MPI 1 and 2. The diagnostic includes a
wrong-sign `Y_n=+0.5` negative control, verifies exact
prescribed dynamic coefficients and complete GetFEM rollback.
No physical/initialization TS is run and no model/solver changes
are made. If the reference passes, the previous LSQR
`INCOMPATIBLE_AT_TOLERANCE` **rest label was a numerical
false negative** and numerical conditioning/termination must
be addressed independently from physical incompatibility.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-exact-rest.sh
```

**Next order:** close the geometry/sign audit on the evidence
once the known solution is calibrated (do not claim the
original rotating DAE is fully initialized), then and only
then start the separately planned `examples/shallow_water.py`
nonlinear DAE initialization qualification. This model has
state-dependent mass and a nonlinear pressure function;
first check whether its algebraic residual is actually
nonlinear in the *algebraic* variables at fixed h,p.

## 20 November 2026 — #28 exact stationary baseline NATIVE PASSED (MPI 1/2, geometry P1/P2)

**Host evidence:** `v130-issue28-exact-rest-b556f2ea-3E78L59H.zip`.
Source `b556f2ead84397b7158b1755d2dfa72ede8122df`,
immutable `scrimp-mpi:v1.3.0-rc3` OCI revision
`1d6d0c90115213298047b28237a7f38087563a35`.

The diagnostic inserted the independently known exact **algebraic**
stationary state WITHOUT LSQR, pseudo initialization or a physical TS:
`h0=5` and `p0=0` (prescribed original differential coefficients),
`e_h=0.5`, `Y_n=-0.5`, `e_p=U_n=U=Y=0`. It evaluated the
real production `IFunction(t0,z,zdot=0)` and rolled back all states.

| Geometry | MPI | `||F_a||_2` | `||F_d||_2` | Wrong `Y_n=+0.5`: `||F_a||_2` |
|---|---:|---:|---:|---:|
| P1 | 1 | `1.4633950546719204e-17` | `2.6084923555407635e-16` | `4.0700566141663` |
| P1 | 2 | `1.470126998300503e-17` | `2.6084923555407635e-16` | `4.0700566141663` |
| P2 | 1 | `1.4763422151112202e-17` | `2.699640468033266e-16` | `4.075223213705125` |
| P2 | 2 | `1.4744250640694402e-17` | `2.699640468033266e-16` | `4.075223213705125` |

**The stationary reference passes** (all < `1e-11`), with a
nonvacuous negative control. This proves that previous resting
`INCOMPATIBLE_AT_TOLERANCE` classifications (`||F_a||~4.7e-9`)
were **numerical LSQR false negatives**: identical fixed h0,p0
have an exact solution of the assembled FE constraints.
Do not loosen production tolerances merely to hide the LSQR failure.

**The #28 boundary geometry/sign *subgate* is closed as a
diagnosed and corrected GWFL literal error**, including true
P1/P2, MPI 1/2 sign and resting-equilibrium validation.
This does NOT close issue #28: for original nonzero rotation,
the corrected fixed-h,p algebraic LSQR residual remains
`~4.18e-8`, above the strict diagnostic threshold; the
traction `Y` is a null algebraic Jacobian column and consistent
initial `Y0` may require differentiated boundary constraints.
No complete rotating DAE initializer has been qualified.

**Next, as explicitly ordered by the user:** study the original
`examples/shallow_water.py` for nonlinear consistent
initialization on a straight-edged rectangular domain.
Verify the algebraic Jacobian and the existing
`explicit=True` frozen gyroscopic treatment. Preserve h0,p0.
Nonlinearity in h,p and state-dependent mass do not
automatically mean that the algebraic equations are
nonlinear in their algebraic unknowns; distinguish the
two levels before claiming a nonlinear algebraic Newton.

## 20 November 2026 — #28 second workstream BEGINS: shipped shallow-water nonlinear DAE, source only

The requested **AFTER geometry/normal sign** gate is now scheduled:
`examples/shallow_water.py` is the original published
**inviscid rectangular** PFEM shallow-water example.
The diagnosed sign convention in the separate viscous disk
sandbox is now documented and its exact P1/P2 resting
equilibrium has passed native MPI 1/2. This does not qualify
the full nonzero rotating disk initial state.

Source inspection of `examples/shallow_water.py`:

- Rectangular domain (original mesh `h=0.1`), no curved
  boundary geometry.
- Differential states `h,p`, with original prescribed
  `h0=3-(sign(x-.5)+1)/3`, `p0=0` and zero normal
  boundary inputs on all four sides.
- Nonlinear, h-dependent momentum mass (`M_p=h*p.Test_p`),
  nonlinear Dirac terms and quadratic pressure contribution
  `0.5*(p.p)/rho` in the effort constitutive relation.
- Original `G=(Gyro(p)*e_p).Test_p` is
  `explicit=True`; the *entire* term is frozen at
  the previously accepted state, and is absent from the
  implicit Jacobian, exactly as explicitly specified
  in the shipped model. The initial original `p0=0`
  makes its initial forcing zero; previous #27
  nonvacuous tests qualified its different-state freezing.
- Crucial nuance: for **fixed h,p**, the full original
  algebraic residual is **affine in e_h,e_p,U_i,Y_i**.
  This is a *nonlinear full DAE* with a linear
  algebraic subproblem, not a demonstration of
  nonlinear-in-algebraic-variables Newton.

New diagnostic-only harness:
`sandbox/probe_v130_issue28_shallow_initial.py`,
`scripts/run-v130-issue28-shallow-initial.sh`.
It reuses the source-proven original shallow-water model
builder from `sandbox/check_v130_shallow_water_mass_jacobian.py`,
which suppresses only the 5000-step physical run and plotting.
No changes to published weak forms, original h0,p0,
mesh or control law. It freezes the original
`explicit=True` gyro at the accepted initial state,
then assembles native GetFEM/PETSc `IFunction` and
`J_aa`. It checks algebraic mass rows zero,
a square algebraic block, then attempts **one**
distributed PETSc direct algebraic correction via MUMPS,
re-evaluating the true `F_a`. Independent algebraic
directional finite differences check the native tangent,
and the audit certifies `h0,p0` bitwise unchanged,
frozen gyro unchanged and complete GetFEM state rollback.

**Native run NOT YET EXECUTED; source-only staging.**
Exact finite residual acceptance remains strict `1e-9`.
The native harness compares MPI 1 and 2, pins
immutable `scrimp-mpi:v1.3.0-rc3`, records local logs
and JSON, rejects any missing or divergent original
model property, and does no initial/main TS steps.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-shallow-initial.sh
```

Even a successful diagnostic projection does **NOT** qualify
a general production nonlinear DAE initializer or
the nonlinear algebraic Newton class, and does
not replace accepted physical-time/energy-power
qualification. These must remain explicit future gates.
