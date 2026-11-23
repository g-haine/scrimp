# SCRIMP v1.3.0 — release sprint, week 48 (23–29 November 2026)

**Decision recorded 20 November 2026.** The user's delivery objective is
**during the week of 23 November 2026**; use **Friday 27 November**
as the planning checkpoint, with the weekend as contingency rather
than as an assumption that work can run unattended.
This is an **accelerated gate-based plan, not an unconditional release
promise**. The authoritative scope remains
[V130_RELEASE_PLAN.md](V130_RELEASE_PLAN.md), points 0–9.
Do not silently delete a committed gate to make the deadline.

## Binding reprioritization — Saturday 21 November 2026

**New P0 release blocker: [#29 accepted-state checkpoint/restart](https://github.com/g-haine/scrimp/issues/29)**,
necessary infrastructure for repeated LHMNC Saint-Venant experiments.
The v1.3.0 deliverable is a **generic, fail-closed cold restart
from a validated complete pH-DAE state**, with original absolute
physical time, all differential + algebraic FE coefficients,
exact mesh/global ordering, original boundary controls, MPI
and per-step PETSc/DAE/physical energy acceptance.
This is NOT a generic high-index initializer or a promise
to serialize PETSc TS/BDF histories. See new release
[gate 2C](V130_RELEASE_PLAN.md) and [technical roadmap](V130_DEVELOPMENT_ROADMAP.md).

**New Maxwell 3D demonstration moves to P2/optional pre-release.**
Only pursue after checkpoint + remaining compulsory v1.3.0
gates pass with spare time; otherwise perform Maxwell
immediately after the release. Keep its original
scientific criteria, existing published examples
and lean GetFEM RT/BDM/Nédélec checks; do not
claim the new 3D example was completed.

**Fresh scientific status** from the user's pinned original
12-case complete physical C1 spinup archive
(`c0b8f3d2`, `t_f=.02`,
`v130-issue28-full-spinup-c0b8f3d2-NSdC2dDW`):
12/12 GetFEM/PETSc native runs PASS,
MPI1/2 full FE states, true endpoint TSCN,
DAE and complete wall plateau.
Observed dynamic h,p temporal order
**1.84153 / 1.96182** (nominal h=.2/.1).
Endpoint weighted power-balance defect relative
to actual exchanged energy at N=40:
**2.80e-5 / 6.94e-5**. Those are
non-negligible compared with early-time ~1e-8;
**the flow need not be stationary merely because the wall
has reached its speed plateau**. Keep physical energy
and nonlinearity caveats open (#6/#28).

**Date rule:** retain Friday 27 November as the
go/no-go planning checkpoint, **not** a reason
to waive checkpoint/physical gates. If #29
cannot be properly qualified, the release
slips; Maxwell does not become a blocker.

## What is already qualified

- **Gate 0 (linear solver performance, native 2D acoustic reference): PASS**;
  bounded 1/2/4/6-rank physical strong scaling and independent v1.2.0
  comparison are archived in the release plan.
- **Gate 1 (selected RT, BDM, Nédélec I/II k=1,2, simplex GetFEM 5.5
  patch work): lean constructor/trace/assembly native tests PASS**.
  Broader all-orders or upstream GetFEM certification is not claimed.
- **Gate 2 (Docker/Ubuntu + VS Code): targeted rc3 native and
  debugging/onboarding smoke PASS**, including caching, MPI startup,
  persistent host output and imported PETSc options warning fix.
  Clean-room multi-platform reproducibility is not yet qualified.
- **DAE consistent initialization: full native implicit DEFAULT PASS**
  at source SHA `3a17a65b7a17fbdade82bd3ea2bd3bb8cd8a53f8`,
  qualified rc3 Docker source
  `1d6d0c90115213298047b28237a7f38087563a35`.
  Conservative and dissipative acoustic models: MPI 1/2, 10
  accepted fixed CN steps, exact q0/p0, ~1e-13 true algebraic
  residual, signed port work and nonzero dissipative power
  balances at ~1e-15, full CN state parity, explicit negative
  Gaussian rejection and independent pseudo/BDF2 regression.
  Output:
  `sandbox/outputs/v130-algebraic-default-h0p2-3a17a65b-H38miGDO/`.
  **This qualifies the two tested linear models, not nonlinear DAE
  initial data or all higher-index formulations.**

## Critical-path table

| Gate in master roadmap | Status as of 20 Nov | Before declaring v1.3.0 releasable |
| --- | --- | --- |
| 0 — Numerical solver / MPI performance | **PASS in tested wave scope** | Retain frozen benchmarks; no full repeat without relevant core change. |
| 1 — GetFEM 5.5 simplex FE | **Lean native PASS** | Integrate supported RT/BDM/Nédélec patches consistently into final image and retain their regressions; do not require a new Maxwell model. |
| 2 — Docker/native/VS Code | **rc3 native PASS** | Final source revision image/wheel smoke and provenance; do not rebuild on every Python edit. |
| DAE default (new) | **FULL native PASS, two linear pH-DAEs** | Preserve explicit pseudo/BDF2, documented incompatible Gaussian, and short memory/exceptional MPI cleanup gate. |
| **#27** state-dependent mass Jacobian | **FULL native gate PASS (MPI 1/2)**: directional FD plus genuine 4-step physical TS with accepted-state frozen G | Freeze the two evidenced campaigns; close #27 with explicit limitations. No automatic claims about exact energy conservation, order two, nonlinear DAE initialization or 4-rank MPI. |
| **#6** viscous SWE MPI/port powers | **Native qualification PENDING.** Instrumented experimental runner and source-only negative controls prepared; source-check CI PASS at `8a20d551`. | Execute short 1/2/4-rank rest then nontrivial physical open/closed smokes, investigate initialization or rank mismatch on FAIL; only later pursue long-time physics, output races and full issue closure. |
| **2C / #29** — Generic accepted-state checkpoint/restart (NEW P0) | **SOURCE/IMPLEMENTATION NOT STARTED**; native full-spinup states available for pilot | Before tag: atomically save full accepted FE state at absolute t*, exact mesh/DoF/control provenance, fail-closed GetFEM/PETSc cold restore, MPI 1→1/1→2/2→1/2→2 continuous-vs-split physical/DAE/power checks, and linear DAE regression. |
| 3 — New 3D Maxwell **(P2 optional)** | **NOT QUALIFIED, DEFERRED FROM CRITICAL PATH** | Only if all P0s are done before tag; else immediately post-release. Preserve 3D H(curl)/H(div), Poynting, orientation and MPI acceptance; published Maxwell examples/FE checks remain in regressions. |
| 4 — Performance figures | **Wave evidence PASS; publication figures incomplete** | Curate reproducible wave/MPI plots and metadata; defer new Maxwell scaling if its new model is postponed. |
| 5 / **#26** — Sphinx and research docs | **Partially complete** | HTML and LaTeX/PDF green on final source, accurately scoped MPI/jacobian/DAE/default, usable first-run path. |
| 6 — All shipped examples / install check | **Incomplete** | Inventory, bounded native smoke per supported model, explicitly report unsupported cases; implement or deliberately defer optional install-check CLI with a documented scope decision. |
| 7 — GUI | **Source migration of historical pseudo PASS; full GUI smoke PENDING** | Generate/run a small script, check interactive startup and error handling; label scope accurately. |
| 8 — Tracked sandbox | **Incomplete** | Preserve qualified scientific runners while moving/archiving as planned; never delete user-local `sandbox/outputs`. |
| 9 — Release package | **NOT STARTED** | Metadata version 1.3.0, wheel/sdist/import smoke, full native suite, release notes, merge/tag and issue decisions after gates. |

## Nonlinear mathematics: two separate workstreams

### P0 for v1.3.0: correct Jacobian of a state-dependent mass

For
`F(t,z,zd) = M(z)zd + R(t,z)`, the true implicit PETSc TS Jacobian is

```text
J_TS(z,zd,sigma) =
    sigma M(z) + dR/dz + D_z[M(z)](•) zd
```

The final term is **not** produced by merely assembling GetFEM's
state-dependent mass tangent and adding `sigma*M`. The current
`scrimp/dphs.py` already defaults to
`SCRIMP_NONLINEAR_MASS_JACOBIAN=finite_difference`; it computes
the contracted term by a centered finite difference in
the `zd` direction of the GetFEM tangent mass, with bounded
perturbations and PETSc storage reuse. **This symmetry-based
construction assumes the GetFEM nonlinear dt brick really
defines a smooth weak-form primitive `Q(z)` whose derivative
is `M(z)=DQ(z)`**. For an arbitrary non-variational
state-dependent mass operator there is no general
`D M(z)[v] zd = D M(z)[zd] v` identity; such models
would need their true contracted derivative assembled
directly or a matrix-free Jacobian action.
Therefore qualify the supported GetFEM representation,
not an unrestricted mathematical class of `M(z)`.
`SCRIMP_NONLINEAR_MASS_JACOBIAN=legacy` is available for
the negative-control A/B. The previous viscous SWE full-Jacobian
directional checks reported relative errors ~1e-11 on MPI 1/2
**with correct/GetFEM-calibrated RHS**; this does not certify
the original shipped `examples/shallow_water.py`.
**Issue #27 is the dedicated v1.3.0 release gate** for that
model, using native original weak bricks and independent `z`/`zd`
directional FD. A native PASS is required before closing it.

### P1 / explicitly conditional for v1.3.0: nonlinear DAE initialization

The qualified linear-only algebraic projection solves

```text
K_aa delta_a = -F_a(t0,z_d0,z_a,0)
```

with fixed prescribed dynamic variables `z_d0`.
Nonlinear extension (**issue #28**) replaces the frozen linear
block by an iterated real algebraic residual Jacobian

```text
J_aa(z_d0,z_a[k]) delta_a = -F_a(t0,z_d0,z_a[k],0)
z_a[k+1] = z_a[k] + lambda[k] delta_a
```

with distributed sparse PETSc, line search, true native residual
re-evaluation, dynamic-state exactness, rank-collective rejection,
rollback, explicit unsupported index/layout guards and memory
cleanup. It requires native positive/negative nonlinear
DAE examples. **Linearizing once is not enough for general
nonlinear consistency**; a damped Newton/continuation strategy
is required. No silent pseudo fallback.

**Deadline rule:** do not insert an unqualified nonlinear Newton
initializer into the v1.3.0 default purely to meet a date.
It may enter v1.3.0 only if the self-contained native physical
gates actually PASS before package freeze. Otherwise keep
the validated linear default and explicit historical pseudo/BDF2,
and advance #28 in the first subsequent feature release.
The shipped `shallow_water.py` requests pseudo explicitly
for historical reproducibility.

## GitHub issue triage and closure contracts

- **#6 (open bug):** parallel port powers/GetFEM/PETSc correctness.
  Much progress is real (MPI wave and dissipative power, bridge,
  Jacobian, state sync), but its original scope also covers
  **nonlinear Saint-Venant, rank 4, powers/outputs and MPI race/deadlock
  regressions**. Close only after those items actually PASS,
  or split/de-scope explicitly with a reviewed narrower issue.
  This is the main MPI correctness risk for a broad v1.3 claim.
- **#23 (open upstream/root-cause investigation):** GetFEM MPI-2
  duplicate full raw RHS and tangent. Keep **open** until the
  underlying assembly/partition root cause is located; the
  calibrated compatibility mode is not equivalent to an upstream fix.
  Document the exact build and residual calibration in the release.
- **#25 (open enhancement):** Docker cache/build-context and
  developer workflow. rc3 targeted native build/cache/VS Code
  gates passed; assess remaining clean-room/cache recovery
  acceptance separately. Close only when those explicit criteria
  are satisfied; otherwise document the residual follow-up.
- **#26 (open enhancement):** full Sphinx docs/examples/developer
  refresh. Close only after final HTML/LaTeX, examples, numerical
  MPI guide and first-run install walk-through are verified.
- **#27 (P0 gate passed 20 Nov):** state-dependent mass Jacobian
  on shipped shallow water: native MPI 1/2 directional FD and
  default/legacy control **PASS**; independent real BDF2 PETSc TS
  on MPI 1/2 with nonvacuous frozen explicit G **PASS**.
  Close with pinned native provenance and limitations.
- **#28 (new P1):** nonlinear constrained DAE initializer.
  Keep open until genuinely qualified; not a pretext
  to delay/relabel the passed linear default.

**No issue is closed by documentation alone.**
Every closure requires a commit, native test path,
environment provenance and explicit scope/remaining limits.

## Week 48 sequence — no unnecessary repeats

- **Monday 23 Nov — correctness and checkpoint design:** retain
  already native-PASSED #27 and new complete-ramp SWE results;
  keep #6/#28 remaining caveats explicit. Define #29 format,
  mesh/FE/absolute-time/control fingerprint and strict
  accepted-state-only import/export; first independent
  continuous-vs-split sandbox proof. Do not delete any
  existing #6/#23 correctness gate.
- **Tuesday 24 Nov — checkpoint implementation/native gates:** 
  qualify #29 actual full-state cold restart at fixed mesh
  first, then cross-MPI 1→2 and 2→1 plus negative controls
  and a second linear DAE. Integrate already-qualified
  GetFEM simplex patches and candidate Docker build/
  dependency provenance. Do **not** spend this day
  developing the new Maxwell 3D model while #29 is red.
- **Wednesday 25 Nov — checkpoint generalization and compatibility:**
  native generic accepted-state restore outside Saint-Venant,
  documentation/test of incomplete/incorrect checkpoint rejection,
  representative examples, small GUI round-trip, tests
  and code/header hygiene. Preserve supported modes and
  explicitly mark BDF/IMEX history restore as unsupported
  until independently qualified.
- **Thursday 26 Nov — reproducibility/docs and conditional Maxwell:**
  final checkpoint API/docs and a copyable restart demo,
  profiling figures, Sphinx HTML/LaTeX, installation
  and release notes/README. **Only if #29 and all other
  mandatory gates are green**, consider a separate,
  bounded Maxwell 3D attempt; otherwise queue it
  immediately after release without weakening gates.
- **Friday 27 Nov — go/no-go gate:** full source/CI
  and native smoke including **#29 cross-MPI restart**,
  wheel/sdist and candidate Docker reproducibility,
  resolved-issue closure evidence, final review/merge/tag
  **only if all P0 release-blocking gates are green**.
  Maxwell new 3D demo is optional; if deferred, list it
  in the next-release backlog. Otherwise publish blockers
  and move the release, not the physical thresholds.

## Immediate next native command

**Historical command — ALREADY PASSED (20 November):** this frozen
Jacobian FD run is retained for reproducibility; do not
rerun it just to close the remaining physical-step gate:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-shallow-mass-jacobian.sh
```

**Native FD result PASSED**, source
`fae066a64b649131351e686ac9383ea11c13c0a6`,
host output `sandbox/outputs/v130-shallow-mass-J-fae066a6-VvvbTRp4/`.
The corrected default vs historical legacy
max errors were `1.06105781e-8` vs `3.84378444e-6`
in BOTH MPI 1 and 2. The remaining gate is the
ACTUAL short physical PETSc TS test, not another
Jacobian FD campaign. No nonlinear DAE
consistent-initialization claim is implied.

**Commit discipline:** keep one integration branch
`release/1.3.0-preparation`; keep native qualification
outputs user-owned and Docker rc3 immutable. Do not
rebuild native dependencies for each Python-only edit.

## 20 November native #27 first trial — frozen rc3 source-overlay import failure

Native source under test
`1f2698f677200eb9a3874bd0d445c7f6513d7683`,
qualified unchanged `scrimp-mpi:v1.3.0-rc3`
(OCI source `1d6d0c90115213298047b28237a7f38087563a35`),
original host evidence:

```text
sandbox/outputs/v130-shallow-mass-J-1f2698f6-5OeuwCYG/
```

**Result: FAIL in Python import, before any GetFEM
assembly, Jacobian or numerical analysis.**
The runner mounted the updated `scrimp/dphs.py` into
the frozen rc3 image but omitted its new companion
`scrimp/utils/dae_initialization.py`, which the
old image predates. Therefore the first
`RUN SWE MASS J mode=default MPI=1` raised:

```text
ModuleNotFoundError: No module named 'scrimp.utils.dae_initialization'
```

**Fix:** mount the exact current
`scrimp/utils/dae_initialization.py` read-only beside
`dphs.py`, matching the previously **native-qualified**
`scripts/run-v130-dae-default-promotion.sh` design.
The runner now checks the helper exists and runs a
lightweight `import scrimp.dphs` preflight before
launching any MPI physics cases. A source-only test
guards both the overlay and preflight against this
exact regression.

No SCRIMP production solver or mathematical
model was modified. The same image and all FD,
MPI and legacy control acceptance criteria remain.
**Issue #27 stays OPEN; native shipped shallow-water
Jacobian qualification PENDING after rerun.**

## 20 November 2026 — fully explicit gyroscopic correction (native pending)

**Decision clarified by user:** the `G` gyroscopic brick
in `examples/shallow_water.py` must be treated
**completely explicitly**, NOT converted into a
fully implicit gyroscopic term merely to make the
Jacobian directional test pass. The previous
`explicit=True` was internally inconsistent:
`IFunction` reevaluated `G` at each *current*
Newton trial state while `IJacobian` omitted
its derivative. At `zd=0` the shipped model
exposed a common directional-FD discrepancy
`2.02366589e-5` in both corrected-mass and
legacy modes, while the corrected mass
Jacobian improved the nonzero-`zd` result
only slightly to that baseline. The
MPI1/2 results agreed to floating-point
roundoff. First full native output:

```text
sandbox/outputs/v130-shallow-mass-J-c9bd39e4-3B7wkKE3/
```

**Production change, not a tolerance relaxation:**

- The original published `S.Brick("G", "(Gyro(p) * e_p) . Test_p",
  ..., linear=False, explicit=True, position="effort")`
  expression/sign is restored unchanged.
- For each PETSc main or pseudo/BDF initialization step,
  `TS.setPreStep()` now assembles the **entire** nonlinear
  explicit GetFEM weak-form RHS at the **previously accepted**
  time and state `(t_n,z_n)` and stores it in a distributed
  PETSc vector. Its value stays *constant* for every
  `IFunction` Newton and line-search evaluation of the
  following step, including step rejection/retry, and
  is **not included in `IJacobian`**.
- Ordinary `assemble_rhs` excludes explicit bricks,
  `assemble_rhs(explicit_only=True)` computes only
  the cached lagged contribution, and the direct residual
  adds it once with its historical GetFEM sign.
  Generic-term MPI calibration is retained,
  including when source-RHS policy is `raw`.
- Unsupported `explicit=True` usages (linear, `dt=True`,
  or source bricks) now raise an error rather than
  silently receiving only half-explicit semantics.
  Explicit-source handling and full-IMEX splitting
  for those unsupported classes are **not** claimed.
- Source-only tests cover accepted-state freezing,
  non-recomputation at Newton trial states, update at
  new accepted times, original (G) form/sign,
  inclusion of both PETSc initial and physical TS
  PreStep hooks and strict comparator controls.

The strict native directional-FD gate
`scripts/run-v130-shallow-mass-jacobian.sh`
now freezes a **nonzero** explicit `G` at a
state deliberately different from the Jacobian
trial state. The diagnostic asserts the cache's
norm is nonzero and its coefficient drift is
**exactly zero** through all finite-difference
and Jacobian trials, and the original
`1e-5` maximum relative FD threshold
and default-vs-legacy negative control
are **unchanged**. Native physics remains
**PENDING** on this new production code.

**Important mathematical tradeoff:** fully lagging
the entire gyroscopic contribution is generally
a **first-order approximation for that term**,
even if the implicit integrator is BDF2.
A skew-symmetric continuous Dirac structure
does not automatically ensure exact discrete
Hamiltonian conservation when this complete
term is evaluated at the preceding state.
A suitable skew/energy-preserving implicit or
specialized split would be a *different* numerical
scheme. We must not claim unchanged old
shallow-water trajectories or exact Hamiltonian
balance from this modification without native
time-integration and power evidence.

**Issue #27** (and wider MPI issue #6) remains
open until both native directional FD and
actual short accepted-step TS/explicit-RHS
tests are satisfied. The previously passed
**linear DAE default** qualification is
unaffected by this narrow nonlinear-only change,
but the general release-wide regression
suite is still required.

## 20 November native #27 second trial — full mass-Jacobian directional FD PASS (historical: before TS run)

**Source:** `fae066a64b649131351e686ac9383ea11c13c0a6` on
`release/1.3.0-preparation`; unchanged qualified image
`scrimp-mpi:v1.3.0-rc3` (OCI source
`1d6d0c90115213298047b28237a7f38087563a35`).
The host-owned native output is:

```text
sandbox/outputs/v130-shallow-mass-J-fae066a6-VvvbTRp4/
```

**PASS:** native full `IFunction` versus `IJacobian` independent
directional difference on the **shipped original shallow-water weak
forms**, MPI 1/2, 1,246 DOFs, with nonzero independent `zd`,
shifts 0 and 1, prescribed perturbation `epsilon=1e-6`,
nonzero gyroscopic explicit RHS norm `0.00536124096` frozen at
a *different accepted state* and **exactly zero cache drift**
during all Jacobian and FD evaluations. No environment flag
was set for the production mass Hessian correction.
Both MPI-rank counts agree to printed precision.

| MPI | `finite_difference` default max relative directional error | Explicit `legacy` control |
| ---: | ---: | ---: |
| 1 | **1.06105781e-8** | **3.84378444e-6** |
| 2 | **1.06105781e-8** | **3.84378444e-6** |

The default is approximately **362 times more accurate** than
`legacy` at this test state. The unchanged strict
`1e-5` maximum error limit, negative-control separation and MPI
comparison all PASSED. The physical gyroscopic weak form and
sign were not modified; `explicit=True` means the complete term
is genuinely lagged at the last accepted state instead of
silently evaluated at Newton trial states without a tangent.

The **Jacobian directional accuracy sub-gate #27 is CLOSED/PASS**.
The **overall issue #27 remains OPEN** under its original
acceptance contract: still requires a short actual PETSc TS
integration with accepted-step explicit-RHS refresh, convergence
and (as feasible) signed energy/power audit. The Jacobian-only
runner explicitly reported `main_TS_executed=False`;
do not substitute this PASS for a production time-integration
or energy-preservation guarantee. In particular, the lagged
`G(t_n,z_n)` term has generally **first-order temporal
accuracy**, even alongside BDF2 or CN.

**Next action:** bounded native short physical TS on original
shallow-water weak forms in MPI 1/2, check positive convergence
reason, final time, accepted/rejected count, explicit cache
update from *accepted states only*, nonvacuous activity if
feasible, final-state MPI parity and a carefully scoped
Hamiltonian/port-work comparison. If the initial data prevent
a safe short solve, report that explicitly rather than
mark #27 fully closed or relaxing the Jacobian tolerances.

## Historical next native command — short *real* shallow-water PETSc TS (now PASSED)

The previously PASSED native Jacobian FD gate is frozen.
Run only the remaining physical-step test on the unchanged
qualified rc3 image. It uses original weak forms and
the published initial data, plus a separately labeled
manufactured rotating initial momentum to ensure
a **nonzero** explicit gyroscopic contribution.
It limits the physical solve to four BDF2 steps, logs
accepted-state explicit cache refresh and Newton
cache invariance, computes the actual Hamiltonian
and checks complete MPI 1/2 state/energy parity.
It **records**, rather than forces to zero, the
signed Hamiltonian drift of first-order lagged `G`.

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-shallow-explicit-physical-ts.sh
~~~

As with other source/native qualification gates, keep
issue #27 OPEN if native TS fails, diverges, the
explicit RHS is vacuous or not frozen, or MPI parity
fails. No release threshold is loosened to meet
the calendar deadline. Successful short-TS
evidence is **not** a proof of order 2 or exact
discrete Hamiltonian balance; a dedicated
time-refinement/energy study is separate.

## 20 November — native #27 physical-step gate PASS (supersedes historical PENDING sections above)

Native source: `9ecd173983252795a3cba13b5879d29127b57638`,
unchanged `scrimp-mpi:v1.3.0-rc3`, OCI source
`1d6d0c90115213298047b28237a7f38087563a35`.
Frozen host evidence:
`sandbox/outputs/v130-shallow-explicit-TS-9ecd1739-LZjvda5C/`.

The original inviscid shallow-water weak forms passed genuine
PETSc BDF2 TS integrations for both published momentum
`p0=[0,0]` and independent nonzero-G manufactured
`p0=[y,-x]`, at MPI 1 and 2. Each run reached
`t_f=4e-4` with four accepted fixed `dt=1e-4` steps,
zero rejected steps/SNES failures and `TSConvergedReason=1`.
The accepted-state gyroscopic RHS is nonzero (up to
`2.4381721e-6` published and `6.88335575e-4` manufactured)
and **never changes within Newton** (`cache drift = 0`).
Across MPI 1/2, maximum scaled full-state difference is
`6.392155562495291e-16`, Hamiltonian difference
`1.2677184090553238e-15`.

Signed `H_f - H_0` was `-0.00206566661` (published)
and `-0.00206572620` (manufactured); both are
`~6.54e-8` relative. These are **recorded defects**,
not a claim of exact structure-preserving discrete power
balance or global BDF2 order for first-order-lagged `G`.
The general nonlinear algebraic DAE initializer remains
issue #28, unqualified. The earlier PENDING/next-command
paragraphs in this sprint are an **historical audit trail**
and must not be used to schedule redundant native reruns.

**Decision:** original #27 release acceptance is satisfied
in its exact tested scope; close #27 with the two SHA/output
provenances. Next P0: original issue #6 on **genuine 1/2/4-rank
viscous Saint-Venant power diagnostics and MPI-safe outputs**,
while #23 remains a separate open GetFEM root-cause investigation.
For #6 use `sandbox/viscous_swe_2d_mpi.py`: preserve
its physical h-weighted boundary powers, volume balance,
both dissipations, non-negativity and true TS completion
checks. Compare trajectories/diagnostics to serial at
**equal physical parameters**. Never treat
`plot_Hamiltonian(with_powers=True)` alone as the physical
viscous SWE power identity.

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

## 20 November 2026 — #28 shallow-water first native FD FAIL; mesh-aware criteria and rotating DAE design

The original rectangular `examples/shallow_water.py` fixed-state
diagnostic at native MPI 1, original mesh 0.1, source
`63bb31b27ef5af720f9ba4d83be823271245a394`,
failed at its directional finite-difference Jacobian test:
`relative_error=5.53585698e-5` using only `epsilon=1e-6`.
Raw host evidence: `v130-issue28-shallow-initial-63bb31b2-28COqKkK.zip`.
This is a **numerical test failure, not proof** of a faulty
production Jacobian or a failed algebraic solve.
The initial published model has `rho=1000`, `g=10`
and fluid height about 3, so coenergy pressure
`e_h` is O(30000), and the tiny epsilon may suffer
severe subtraction cancellation.

New source-only native diagnostic commits `86faa546` and
`c7e3cd05`: scan central FD steps 0.1, 0.01, 0.001,
0.0001 and 0.000001, record both absolute and relative
errors and the native `IFunction` algebraic residual.
Require BOTH the 0.1 and 0.01 steps to satisfy the
original unchanged `1e-7` relative FD tolerance;
also require true post-correction algebraic residual
`<=1e-9` and the independently assembled affine identity.
The model, dynamic initial data, frozen explicit gyro,
FE families and control ports are **unmodified**.
Crucially, write diagnostic JSON even on a numerical
rejection so the next failure is interpretable.
New native run has NOT been executed.

### Mesh size changes how to interpret `||F_a||_2`

The present implementation uses
`||F_{a,h}||_2 <= atol + rtol ||F_{a,h}^{(0)}||_2`.
It does **not** explicitly mention spatial mesh step,
but the residual is a coefficient vector from weak
forms, so its Euclidean norm changes with FE basis,
DOF count, shape regularity, integration weights
and conditioning. On a quasiuniform 2D mesh, a
fixed smooth volumetric defect with nodal load
components O(h²) and O(h⁻²) unknowns can show
Euclidean residual O(h), *without any improvement*
in physical compatibility. Volume and boundary
residuals have different natural scalings and units.

Future release acceptance must combine original
GetFEM/PETSc residuals and numerical backward errors
with **blockwise dual/Riesz norms**
`sqrt(F_a^T R_a^{-1} F_a)` using suitable Gram
operators, physical relative boundary
`L²(Γ_h)` constraint errors, and h-refinement
studies. This has **not** been implemented as
a generic production criterion; do NOT loosen
current DAE thresholds just to force a PASS.
The user-facing treatment is in
`docs/source/initialization.rst`.

### Original rotating viscous SWE: still NOT initialized

Fixing the GWFL matrix sign and validating the
exact P1/P2 stationary equilibrium do not determine
the rotating initial traction `Y0`. Native
structure confirms `∂F_a/∂Y=0` and nonzero
`∂F_d/∂Y`; naive square Jaa Newton cannot determine
this boundary multiplier. The next intended
mathematical work, after the simpler rectangular
shallow-water nonlinear DAE gate:

1. Verify compatibility of the **reduced discrete**
   algebraic boundary constraints for fixed h0,p0.
   Reject rather than silently modify prescribed data.
2. Compute regular algebraic variables (efforts and
   boundary controls) at the fixed original h0,p0.
3. Differentiate the reduced discrete trace
   constraints `C(q,t)=0`: `C_q qdot + C_t=0`.
4. Assess rank/solvability of a bordered linearized
   system coupling dynamic mass equations to
   unknown `qdot0,Y0`, with gyroscopic, viscous
   and prescribed-control contributions intact.
5. Verify original boundary power, dissipation,
   MPI 1/2 and spatial refinement before claiming
   a physical rotating DAE initializer.

This is a **planned** possible index-aware
initialization method, not an established index
class or a native-qualified production algorithm.

**Native next command:**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-shallow-initial.sh
```


## 20 November 2026 — original rectangular shallow-water fixed-dynamics algebraic diagnostic NATIVE PASS, MPI 1/2

**Authoritative user-supplied native archive:**
`v130-issue28-shallow-initial-2098e68c-4vxWCDRh.zip`.
Source `2098e68c191662103c3a404e38b7bef53298a307`,
Docker `scrimp-mpi:v1.3.0-rc3` with OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Original `examples/shallow_water.py`; original `mesh=0.1`,
original dynamic initial coefficients `h0,p0` and wall controls,
original explicit gyro.

| Metric | MPI 1 | MPI 2 |
|---|---:|---:|
| Total DOFs | `4270` | `4270` |
| Differential DOFs | `2035` | `2035` |
| Algebraic DOFs | `2235` | `2235` |
| Initial algebraic `||F_a||_2` | `1281.7812562535455` | `1281.7812562535455` |
| After **one direct algebraic solve**, native `||F_a||_2` | `2.199213556173762e-12` | `2.3277177796294783e-12` |
| Affine relation `||F_a(z+delta)-F_a(z)-J_aa delta||_2` | `1.571067363774428e-13` | `1.590053793547698e-13` |
| Best FD/Jacobian relative error | `4.246216760258867e-10` | `4.607378233546647e-10` |
| Algebraic mass-row max | `0` | `0` |
| PETSc MUMPS direct reason | `4` | `4` |

**The previously suspicious FD error at epsilon=1e-6
was a numerical differencing artifact, not a GetFEM
Jacobian defect.** At `epsilon=0.1` the relative
difference is `4.25e-10` (MPI 1), at `0.01` it is
`4.93e-9`, at `0.001` `4.73e-8`,
at `0.0001` `5.59e-7` and at
`1e-6` `5.5359e-5`. MPI 2 shows the same
loss-of-precision trend. Both sufficiently sized
perturbations pass the original unchanged `1e-7`
FD requirement; the native algebraic residual
passes the unchanged `1e-9` threshold by more
than 400x. Diagnostic status `PASS_DIAGNOSTIC`.

Native JSON further confirms: `h0,p0` bitwise
preserved, entire GetFEM model bitwise restored,
explicit gyroscopic RHS frozen and unchanged
(its norm is **zero** for original `p0=0`);
both initial/physical TS remain unexecuted.

**Qualified scope:** one direct algebraic correction
in the original full **nonlinear-state-dependent DAE**
with **affine-in-algebraic-unknowns** constraints,
fixed prescribed differential data and nonlinear mass.
This is a **read-only diagnostic only**: no
production nonlinear DAE initializer, general
nonlinear algebraic Newton, physical integration
or power balance has yet been qualified.

**Next gate:** extend SCRIMP's opt-in consistent
initializer to nonlinear *state-dependent* DAE
systems with **regular algebraic blocks affine
in algebraic unknowns**, ensuring the same
real production `IFunction` and Jacobian,
no hidden changes to initial dynamic coefficients,
explicit-gyro accepted-state freezing, robust
MPI 1/2 error/rollback paths, numerical/physical
acceptance separation, and no fallthrough to
historical pseudo TS. Do NOT silently enable
this path for the rotating viscous disk, where
`Y` is an algebraic null-column and differentiated
boundary constraints are required. #28 stays open.


## 20 November 2026 — next gate staged: opt-in *production* nonlinear-state / affine-algebraic DAE initialization

After the **original native diagnostic** of
`examples/shallow_water.py` passed at MPI1/2
(`||F_a||=2.20e-12`/`2.33e-12` after a
single direct algebraic correction), production code
now exposes an **opt-in**, distinct initializer
`init_step_ts_type="nonlinear_affine"`.

- Existing default linear `algebraic` and historical
  `pseudo`/`bdf` behavior remain unchanged.
- Original differential FE coefficients `q0=(h0,p0)`
  remain bitwise unchanged. Only a **regular,
  affine-in-algebraic-variables** block is eligible.
- True native `IFunction` initial and final residuals
  and analytic `J_aa` are used, not a manufactured proxy;
  nonlinear mass and stiffness are allowed at fixed `q0`.
- The initializer freezes the original accepted-state
  explicit forcing, performs **one direct PETSc
  PREONLY+MUMPS algebraic solve**, checks the true
  affine identity and strict algebraic accuracy,
  and commits only accepted algebraic coefficients.
- Zero algebraic Jacobian columns are rejected before
  factorization, with full state rollback on failure:
  this **explicitly excludes** the original rotating
  viscous disk multiplier `Y` whose full `J_aa` is
  singular. It does NOT handle high-index initial
  differentiated boundary constraints.
- It does not implement or claim a general nonlinear
  algebraic Newton; it is a new narrow regular
  nonlinear-*state* DAE class.
- No automatic default migration or physical TS
  release gate is claimed.

New production implementation:
`scrimp/utils/nonlinear_affine_dae.py`,
with narrowly scoped dispatch in `scrimp/dphs.py`.
New native harness:
`sandbox/probe_v130_issue28_production_nonlinear_affine.py`
and `scripts/run-v130-issue28-production-nonlinear-affine.sh`.
The runner uses the original rectangular shallow-water
model, its actual production `DPHS.init_step()`,
and the immutable rc3 Docker image, MPI 1/2.
It explicitly verifies the **published** updated
algebraic GetFEM coefficients via native `IFunction`,
original differential coefficients, and caches.
No main physical TS is executed in this stage.

**SOURCE ONLY / NATIVE PRODUCTION GATE STILL PENDING.**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-production-nonlinear-affine.sh
```

After this native gate, separate negative controls
(singular/high-index disk and incompatible fixed
differential data), normal TS short-run with power
balance, MPI parity, and h-refinement remain release
dependencies. This production stage should not be
listed as fully qualified until user-supplied native
evidence is checked.

## 20 November 2026 — Production `nonlinear_affine` initial projection: NATIVE MPI1/2 PASS; host-only comparator typo fixed

The user ran the *actual* production `DPHS.init_step()` with
`init_step_ts_type="nonlinear_affine"` on the original
`examples/shallow_water.py` with original mesh `0.1`,
prescribed `h0,p0`, homogeneous boundary controls and
original explicit gyroscopic term, **without any physical TS**.
Source **`de76ff77b9ab33a0bded5b78096025511211589a`**
and pinned immutable image `scrimp-mpi:v1.3.0-rc3`
(OCI revision `1d6d0c90115213298047b28237a7f38087563a35`).

**Evidence archive personally supplied by the user:**
`v130-issue28-production-na-de76ff77-l7JjFNvy.zip`
(manifest, native MPI1/MPI2 logs and separately saved
`production_nonlinear_affine.json` for each rank).

| Metric | Native MPI 1 | Native MPI 2 |
|---|---:|---:|
| Total, differential, algebraic DOFs | 4270, 2035, 2235 | 4270, 2035, 2235 |
| Original `||F_a||` | `1281.7812562535455` | `1281.7812562535455` |
| **Published** true `||F_a||` | **`2.199213556173762e-12`** | **`2.3277177796294783e-12`** |
| Accuracy target | `1e-11` | `1e-11` |
| Affine identity error | `1.571067363774428e-13` | `1.590053793547698e-13` |
| PETSc MUMPS convergence reason | `4` | `4` |
| `max_abs_algebraic_mass_rows` | `0` | `0` |
| Algebraic Jacobian null columns | `0` | `0` |
| Prescribed differential `h0,p0` preserved bitwise | `true` | `true` |
| Algebraic corrections **published** to GetFEM | `true` | `true` |
| Accepted original-state explicit gyro cache | `ready` | `ready` |

Both reports have `status=PASS_NATIVE_PRODUCTION_INITIALIZER`,
with native final residual below the strict target.

**Host-only end-of-script mistake, not a native solver failure:**
after both runs printed `PASS`, the runner called
`math.isclose(a,b,rtol=...,abs_tol=...)`; the Python standard
library expects **`rel_tol`**, not `rtol`, so this final
comparison raised `TypeError` and made the shell exit nonzero.
Independent replay **directly against the supplied archive**
confirms that the two final residuals differ by only
`1.2850422345571612e-13` and satisfy the intended tolerance.

Commit **`e389cd8bbf419b182f58b420e45e496f30fd6a19`**
fixes this by extracting a pure-Python archived native
evidence comparator:
`sandbox/compare_v130_issue28_production_nonlinear_affine.py`,
adds positive and negative source-only regression tests,
and supports **offline revalidation without rerunning Docker**:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash scripts/run-v130-issue28-production-nonlinear-affine.sh --compare-only \
  sandbox/outputs/v130-issue28-production-na-de76ff77-l7JjFNvy
```

**Qualified scope:** the single opt-in real production initial
algebraic projection for an original **nonlinear-state**
but **affine-algebraic** DAE with rectangular geometry,
distributed MPI 1/2, strict residuals, bitwise fixed
differential data, published algebraic coefficients,
and native PETSc MUMPS convergence. No actual physical TS
step was taken, no numerical power balance was checked,
no nonlinear-in-algebraic Newton was exercised, and the
rotating viscous disk traction `Y0` remains **UNINITIALIZED**.

**Next qualification gates:** explicit native negative control
for original rotating open-loop `Y` null columns,
negative/rollback on genuinely incompatible original
initial data, short physical-time integration with
power/Hamiltonian budget, and mesh-refinement /
physically scaled residual norms. Keep issue #28 OPEN.

## 20 November 2026 — next native gate: reject high-index rotating disk (source-only staged)

The actual original rectangular shallow-water **production**
opt-in `nonlinear_affine` initialization passed natively in
MPI 1 and 2 (`||F_a||~2.2e-12`, original h0,p0
bitwise fixed). Its last host-only comparator typo was fixed
in `e389cd8` and the archived native JSON replayed
without rerunning GetFEM/PETSc.

The **next release-safety requirement** is a negative
control using the *original rotating viscous Saint-Venant*
open-loop disk with corrected GWFL wall velocity,
nonzero `omega=0.2`, and **curved geometric P2**.
Native #28 structural audit established that its
boundary traction `Y` has **zero columns in J_aa**
but nonzero dynamic coupling, so the new regular
`nonlinear_affine` algorithm MUST NOT accept or modify
its prescribed rotating h0,p0. A singular square
algebraic solve must be rejected *before* MUMPS,
without fallback to pseudo TS.

Staged (source-only; **native run pending**):
`sandbox/probe_v130_issue28_negative_high_index.py`
and `scripts/run-v130-issue28-negative-high-index.sh`.
The probe calls the real production `system.solve()`
with only `init_step_ts_type=nonlinear_affine` overridden,
catches **only** the expected
`UnsupportedNonlinearAlgebraicDAE` exception containing
`null columns`, checks at least all `Y` DOFs are detected
as null, verifies whole original GetFEM coefficient vector
is bitwise restored and no physical/pseudo TS ran.
It writes independent MPI1/MPI2 JSON and checks rank parity.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-negative-high-index.sh
```

A PASS here means **expected safe rejection**, not
that the rotating high-index `Y0` has been initialized.
The differentiable boundary-constraint/bordered
initialization of the original rotating problem,
negative incompatible-data rollback, short physical TS,
energy/power budget and mesh-aware residual norms
remain separate pending release gates.

## 20 November 2026 — actual high-index `Y` rejection MPI 1/2 PASS; post-solve rollback gate staged

**Native PASS of EXPECTED REJECTION, verified from uploaded archive:**
`v130-issue28-negative-high-index-bfe45f1f-oLsALoUF.zip`
(manifest, logs, P1 not involved; original *curved P2*
viscous rotating disk `omega=0.2`, MPI1/2).
Source `bfe45f1fdc1adf55776d6e05f7bf793ba26c374e`,
immutable rc3 OCI `1d6d0c90115213298047b28237a7f38087563a35`.

Actual `DPHS.solve()` calling production `nonlinear_affine`
raised exactly the expected `UnsupportedNonlinearAlgebraicDAE`
with **132 numerically null algebraic columns for 132
traction-Y DOFs**. The model is 3300 DOFs, 1452 dynamic,
1848 algebraic. GetFEM entire state was bitwise restored
under both MPI ranks, no pseudo/BDF or physical time step
executed; the report is `REJECTED`, and JSON confirms no
high-index rotating initialization was claimed.
Initial native `||F_a||_2 ~ 0.250963115`, tolerance
`2.60963115021321e-9`; the rejection was **structural
before factorization**, not a statement on numerical LSQR
feasibility or the differentiation index. The same
132 columns and the same exception were found in MPI1/2.

**Important limitation of that safety gate:** rejection
for singular `J_aa` occurs **before publishing an
algebraic Newton candidate**. It therefore does not
by itself validate rollback of PETSc `GlobalVec` when
an algebraic trial has been computed and fails the
final residual / precision acceptance check.

Code fix `ceb0da9e9fd3bfa7f4c657537123896297a099c8`
now snapshots the original `system.GlobalVec` before
the nonlinear-affine trial and restores it along with
the original GetFEM coefficients and explicit-RHS
metadata/cache on rejection. Successful accepted
initialization remains unchanged.

New dedicated **numerical**, NOT **physical-data**,
negative test staged:
`sandbox/probe_v130_issue28_negative_accuracy.py`
and `scripts/run-v130-issue28-negative-accuracy.sh`.
It uses unchanged shipped rectangular shallow-water
`h0,p0` and wall controls but deliberately configures
`init_algebraic_accuracy=1e-16`, well below the
previous native achievable `~2.2e-12` algebraic
residual. It demands a real computed MUMPS correction,
then an expected `InconsistentDAEInitialData` rejection;
verifies *both* GetFEM and PETSc `GlobalVec` bitwise
rollback, full explicit RHS metadata and previous
contents preserved, and no time steps, in MPI1/2.

**NATIVE NEGATIVE-ACCURACY TEST PENDING**, separate
from the previously passed high-index rejection.
This does not establish a test with genuinely
physically incompatible fixed differential data;
the latter remains an independent future gate.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-negative-accuracy.sh
```

After these safety gates, stage short actual time
integration on the original rectangular SWE and
check Hamiltonian, boundary power and explicit gyro
freeze properties before claiming physical runtime
qualification. The rotating high-index Y0 problem
still needs index-aware differentiated constraints.


## 21 November 2026 — #28 native numerical-rejection rollback gate PASS (MPI 1/2)

**User-supplied immutable native evidence**:
`v130-issue28-negative-accuracy-20b39eb2-d3D8JQf1.zip`
including manifest, rank-specific JSON and complete logs.
Source `20b39eb2e4a512b7536324331802c223dc924849`;
Docker image `scrimp-mpi:v1.3.0-rc3`,
OCI revision `1d6d0c90115213298047b28237a7f38087563a35`.

The original `examples/shallow_water.py` rectangle/mesh=0.1,
original differential h0,p0 and zero boundary controls
were **not modified**. On the real production
`nonlinear_affine` initializer, `init_algebraic_accuracy`
was deliberately set to `1e-16`, substantially below
the independently verified attainable ~`2.2e-12`
native floating-point residual.

| Metric | MPI 1 | MPI 2 |
|---|---:|---:|
| Native candidate `||F_a||_2` after successful MUMPS correction | `2.199213556173762e-12` | `2.3277177796294783e-12` |
| Strict imposed accuracy | `1e-16` | `1e-16` |
| PETSc linear solver converged reason | `4` | `4` |
| Production status | `REJECTED` | `REJECTED` |
| Exception | `InconsistentDAEInitialData` | `InconsistentDAEInitialData` |
| Bitwise GetFEM rollback | PASS | PASS |
| Bitwise PETSc `GlobalVec` rollback | PASS | PASS |
| Explicit forcing cache ready/time/content restoration | PASS | PASS |
| No physical/pseudo TS | PASS | PASS |

**PASS_EXPECTED_NUMERICAL_REJECTION** in MPI1 and MPI2.
This specifically validates rollback **after a candidate has
actually been computed**, complementing the separate native
high-index rotating `Y` rejection (132 null columns), which
occurs before attempting the algebraic solve.

The test demonstrates *numerical accuracy rejection*, **NOT**
physical incompatibility of prescribed differential initial
data, and does not determine the rotating viscous `Y0`.

Next release gate: a **short physical-time trajectory**
from the now qualified positive production initialization
on the **original rectangular** shallow-water model,
keeping original h0,p0 and controls. Independently measure
mass, Hamiltonian and zero-control boundary power;
track accepted TS steps, nonlinear failures and MPI
1/2 parity. BDF2 is **not** a discrete energy-preserving
scheme: any measured energy drift must be reported rather
than declared exactly power-preserving. Do not invent
a conservation certificate or automatically approve
the higher-index rotating DAE.

The original rotating viscous P2 case remains an open
high-index *differentiated-constraint* initialization
problem. #28 remains OPEN.


## 21 November 2026 — native MPI1/2 PASS: rejected computed algebraic trial fully rolled back

**Native user evidence archive:**
`v130-issue28-negative-accuracy-20b39eb2-d3D8JQf1.zip`,
source `20b39eb2e4a512b7536324331802c223dc924849`,
immutable image `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
This is the **unchanged original physical rectangular**
`examples/shallow_water.py` with an **artificially strict
numerical accuracy** `init_algebraic_accuracy=1e-16`
(no changes to physical `h0,p0` or wall conditions).

| Native metric | MPI 1 | MPI 2 |
|---|---:|---:|
| Original `||F_a||_2` | `1281.7812562535455` | `1281.7812562535455` |
| Computed algebraic candidate `||F_a||_2` | `2.199213556173762e-12` | `2.3277177796294783e-12` |
| Required accuracy | `1e-16` | `1e-16` |
| Native algebraic MUMPS convergence reason | `4` | `4` |
| Rejected computed correction norm | `643262.3257684928` | `643262.3257684928` |
| Algebraic affine-identity error | `1.571067363774428e-13` | `1.590053793547698e-13` |
| Initial/GetFEM rollback bitwise | PASS | PASS |
| PETSc `GlobalVec` rollback bitwise | PASS | PASS |
| Frozen explicit RHS previous ready/time/contents restored | PASS | PASS |

Both ranks returned `PASS_EXPECTED_NUMERICAL_REJECTION`
via the intended `InconsistentDAEInitialData` exception
**after** successfully computing an algebraic Newton
correction. The production report is `REJECTED` with
`accuracy_target_met=false`, no pseudo-initialization
TS and no physical time step.

**Qualified:** accepted production affine-algebraic projection
on the original nonlinear-state shallow water (earlier MPI 1/2);
safe structural rejection of rotating high-index `Y`
(earlier P2 MPI 1/2); **and now safe post-solve numerical
rejection with complete model, PETSc and explicit cache rollback**.

**Not qualified:** rejection due to genuinely physically
incompatible prescribed differential initial data;
general nonlinear-in-algebraic-variable Newton;
a physical-time run, Hamiltonian/power budget and mesh
independent physical acceptance; and the *original*
rotating viscous SWE initial multiplier `Y_0`.
Issue #28 remains OPEN.

**Next independent native gate:** a deliberately short
physical-time simulation of the original rectangular
`examples/shallow_water.py` using explicit opt-in
`nonlinear_affine` initialization. Compare full native
Hamiltonian and port-power balance with the original
inviscid homogeneous-boundary equations. Do not
modify source initial values, mesh, control laws or
the explicit-gyro accepted-state contract merely
to force an energy PASS. Keep rank-1/rank-2
comparison and explicit actual time-step counts.


## 21 November 2026 — staged first SHORT PHYSICAL TS qualification after nonlinear_affine (NATIVE PENDING)

The original published `examples/shallow_water.py` passed
native `nonlinear_affine` production algebraic initialization
on MPI1/2 with `||F_a||~2.2e-12` and bitwise fixed original h0,p0.
The independent **numerical post-solve rejection** and
full GetFEM/PETSc/explicit-cache rollback have now also
passed native MPI1/2, evidence archive
`v130-issue28-negative-accuracy-20b39eb2-d3D8JQf1.zip`.
See the preceding qualification record.

**NEXT physical gate** is a deliberately short 3-step run,
rather than another projection-only check:
`sandbox/probe_v130_issue28_shallow_short_physical.py`,
`scripts/run-v130-issue28-shallow-short-physical.sh`.

Test design:
- Load the original shipped model **without changing its weak
  forms**, rectangular mesh `h=0.1`, P2 h / discontinuous P1 p,
  original physical h0,p0, four prescribed **zero** wall controls,
  original `explicit=True` gyroscopic term, original energy terms.
- Suppress **only** the published 5000-step solve and
  plotting/export during model construction.
- Opt into `init_step_ts_type=nonlinear_affine` and run the
  **actual** production `system.solve()` and main PETSc BDF2
  integrator for `dt=1e-4` through `t_f=3e-4`.
  Use distributed PETSc/MUMPS, identical compiled pinned
  `scrimp-mpi:v1.3.0-rc3` image, MPI 1 and 2.
- Require accepted/converged TS, no rejected steps/SNES
  failures, correct final timestamp and initial saved h0,p0
  **bitwise**; original production initialization status PASS.
- Independently calculate each saved-sample total mass
  `int_Omega h`, kinetic `int h|p|²/(2rho)` and
  potential `int rho*g*h²/2`. Compare against SCRIMP's
  registered original Hamiltonian and measure drift.
- Independently calculate boundary powers `int_{Gamma_i}
  U_i*Y_i` on all four regions. With homogeneous imposed
  controls, supplied boundary power must be ~0.
- First *pilot* sanity gates: mass relative drift `<=1e-5`,
  Hamiltonian relative drift `<=5e-3`, absolute boundary
  power `<=1e-7`. Record full timeseries and defects.
  These values are **initial pilot thresholds**, not
  inferred spatial/temporal convergence certificates.
- **IMPORTANT:** BDF2 is **not symplectic/energy-preserving**
  in general. This run can check small short-time drift,
  but must never be labelled an exact discrete
  port-Hamiltonian power-preservation proof. A proper
  refinement/energy-exchange convergence study is a
  separate release gate.
- Retain rank-specific JSON/logs on failure and assert
  MPI1/2 trajectory parity. If the first full TS fails,
  diagnose actual TS reason and physical operator
  assembly; do not silently switch to pseudo TS.

**Source-only staged; native physical qualification PENDING.**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-shallow-short-physical.sh
```

This is the original **rectangular** inviscid shallow-water
problem. The original rotating **viscous, curved P2** case
has a different high-index constraint structure; its
initial traction `Y0` still requires an independent
index-aware differentiated-constraint algorithm,
and is NOT qualified by this test.


## 21 November 2026 — #28 opt-in production DAE initial projection + real short BDF2 TS staged (native PENDING)

**Previously native qualified (immutable `scrimp-mpi:v1.3.0-rc3`):**

- True `nonlinear_affine` initial projection on the
  original rectangular `examples/shallow_water.py`:
  `||F_a||=2.1992e-12` (MPI1) and `2.3277e-12`
  (MPI2), original differential h0,p0 unmodified.
- Native structural **expected rejection** of all 132
  algebraic null columns for rotating P2 traction Y.
- Native expected **post-computed-correction**
  numerical rejection at `accuracy=1e-16`,
  including exact rollback of all GetFEM fields,
  PETSc `GlobalVec`, and accepted-state explicit RHS
  metadata/content in true MPI1/2.

**Next source-only staged gate, NOT YET NATIVE PASSED:**
`sandbox/probe_v130_issue28_production_short_physical.py`
and `scripts/run-v130-issue28-production-short-physical.sh`,
with a separate archived JSON comparator
`sandbox/compare_v130_issue28_production_short_physical.py`.
This deliberately reuses the already tested original
weak-form, explicit-gyro and actual PETSc TS harness
`sandbox/probe_v130_shallow_explicit_ts.py`.

It retains the **original** `examples/shallow_water.py`
inviscid rectangular physical weak forms,
original `mesh=0.1`, original `h0=3-(sign(x-.5)+1)/3`,
original `p0=0`, and zero controls on all four
boundaries. Only the short-test time scheme
(BDF2 with four physical steps `dt=1e-4`)
and explicit selection of production
`init_step_ts_type="nonlinear_affine"` are test
overrides; there is **no manufactured nonzero
momentum**, no pseudo-initialization TS, and
no changed physical coefficients.

The probe calls the actual production `init_step()`
from within `solve()`, captures its original-state
publication report, then verifies:
- exactly four genuine accepted PETSc BDF2 physical
  steps, with no TS rejections or SNES failures;
- native `IFunction` actually evaluated throughout;
- accepted-time explicit gyro forcing frozen during
  all implicit Newton residual iterations;
- complete finite full-state trajectory and
  positive finite mechanical Hamiltonian;
- MPI1/MPI2 relative **full-state** and Hamiltonian
  parity `<=5e-8`, with identical sample times;
- the original h0,p0 unchanged by initialization and
  original published physical boundary controls.

**Important scientific boundary:** this test will
only **measure** the signed Hamiltonian drift.
BDF2 with a lagged explicit gyroscopic brick does
not guarantee an exact time-discrete
port-Hamiltonian power balance. The test must
not assert zero signed drift as a physics theorem.
A physically scaled energy, integrated boundary
power and time/mesh convergence study is a
**separate** subsequent release gate. No claim
of full rotating viscous DAE `Y0` consistency
or mesh-independent initialization follows.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-production-short-physical.sh
```

On host-side comparator issues, **replay saved
native results without Docker** using
`--compare-only <evidence-directory>` rather
than rerunning expensive physical steps.


## 21 November 2026 — TRUE shallow-water physical BDF2 after nonlinear_affine: MPI 1/2 NATIVE PASS

**User native evidence:** `v130-issue28-shallow-short-physical-d0cde0af-4DlgCMdU.zip`.
Source `d0cde0af5190290e0567bb06424d4ce8e1d2bfae`,
pinned Docker `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Original shipped `examples/shallow_water.py` and original
rectangular mesh step `0.1`, weak forms, FEM, initial
`h0=3-(sign(x-0.5)+1)/3`, `p0=0` and all four
homogeneous boundary controls. Production
`init_step_ts_type="nonlinear_affine"`; **no
pseudo-initialization**, followed by true BDF2 PETSc
TS with `dt=1e-4, t_f=3e-4`.

| Genuine native observable | MPI 1 | MPI 2 |
|---|---:|---:|
| TS convergence reason | `1` | `1` |
| Physical accepted / rejected steps / SNES failures | `3 / 0 / 0` | `3 / 0 / 0` |
| Original initial projection `||F_a||_2` | `2.199213556173762e-12` | `2.3277177796294783e-12` |
| Mass-relative drift | `2.310153069550858e-15` | `3.1986734809165835e-15` |
| Hamiltonian-relative drift / H0 | `1.437729942579708e-7` | `1.4377299161310965e-7` |
| Homogeneous boundary power (each edge and total) | `0` | `0` |
| Original H0 | `31636.2594029492` | `31636.259402949203` |
| Final H | `31636.254854509458` | `31636.254854509545` |
| Signed final H-H0 | `-0.0045484397414838895` | `-0.004548439657810377` |

Registered Hamiltonian agrees with native kinetic/potential
integrals. MPI absolute Hamiltonian sample discrepancy
`<=8.73115e-11`; relative `<=2.76e-15`.
The true mass is `~2.49903774955191` and stays at
floating-point precision. Both MPI implementations retain
original dynamic initial data exactly.

**Critical interpretation: the returned
`relative_energy_defect_exchanged=1` is an
UNINFORMATIVE metric for this closed problem, neither
a 100% failure nor a preservation certificate.**
Here the boundary work `W(t)=∫Σ_i U_iY_i dt≡0`,
hence the balance defect `D_H(t)=H(t)-H(0)-W(t)`
equals `H(t)-H(0)`. The chosen denominator is
`max(max|H-H0|,max|W|,1e-14)`, which equals
`max|D_H|` whenever the energy drifts and no
input work occurs. The ratio equals **1 by definition**.
The scientifically informative measurements are
`ΔH=-4.54844e-3` and `|ΔH|/H0≈1.43773e-7`,
plus their change under refinement in **time**.
For BDF2, exact discrete energy conservation
is NOT guaranteed; no exact pH power balance
can be claimed from this short smoke test.

**Sampling caveat:** native TS accepted `3` steps
but the stored trajectory has only timestamps
`[0,1e-4,3e-4]`, missing `2e-4`.
This does not invalidate the actual accepted
steps or final total energy and mass. However
subsequent power integration with nonzero
controls MUST use actual saved intervals,
not assume every accepted time was recorded.

**Qualified:** opt-in production nonlinear-affine
initialization followed by **real** short nonlinear
PETSc BDF2 physical-time integration,
mass to rounding, correct zero boundary power,
computed Hamiltonian and native MPI 1/2 parity.

**Not yet qualified:** exact time-discrete power
conservation, quantitative energy time/mesh
convergence, general nonlinear-in-algebraic
Newton, physically inconsistent fixed differential
data rejection, or the higher-index rotating
viscous disk initial traction `Y0`.
**Issue #28 remains OPEN.**

**Next independent pre-existing gate (already staged):**
four real physical BDF2 steps on original
`examples/shallow_water.py` using opt-in
`nonlinear_affine`, requiring complete full-state
MPI1/2 parity and checking freeze of original
explicit gyro throughout the implicit Newton
iterations. It measures Hamiltonian drift but
does not assert exact conservation.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-production-short-physical.sh
```

Afterward, time-refine BDF2 at fixed
physical `T` and track the signed
`H(T)-H(0)-∫P_boundary dt`, energy and
mass as well as full-state changes,
rather than this self-normalized
closed-system ratio.


## 21 November 2026 — Original shallow-water 4-step BDF2 after production nonlinear_affine: MPI 1/2 NATIVE PASS

**Actual user-uploaded native ZIP:**
`v130-issue28-production-short-ts-f36c6042-3Fpcl3DV.zip`,
source `f36c6042181470688f29fd4c1a75d1640f2d0c08`,
immutable `scrimp-mpi:v1.3.0-rc3` OCI revision
`1d6d0c90115213298047b28237a7f38087563a35`.

Original published `examples/shallow_water.py`,
rectangular physical mesh `h_mesh=0.1`, original
h0,p0 and zero boundary controls, unchanged
weak forms and original explicit gyro.
Production opt-in `DPHS.init_step()` with
`init_step_ts_type=nonlinear_affine` followed
by **four genuine physical BDF2 PETSc TS**
steps (`dt=1e-4`, `T=4e-4`).

| Actual native result | MPI 1 | MPI 2 |
|---|---:|---:|
| Post-projection original `||F_a||_2` | `2.199213556173762e-12` | `2.3277177796294783e-12` |
| Accepted/rejected steps | `4 / 0` | `4 / 0` |
| SNES failures / IFunction evaluations | `0 / 15` | `0 / 15` |
| Signed final Hamiltonian change | `-0.0046648074239783455` | `-0.004664807369408663` |
| Relative Hamiltonian drift / H0 | `1.474512951914752e-7` | `1.4745129346656575e-7` |
| Explicit gyro RHS drift within Newton | `0` | `0` |
| Max nonzero explicit gyro RHS norm | `2.9865184694522523e-5` | `2.9865184694484807e-5` |

The gyro is **not vacuous**: norms at accepted
step times 0, `1e-4`, `2e-4`, `3e-4` were
`0`, `4.866269214e-7`, `6.474541182e-6`,
`2.986518469e-5`, while its value remained
bitwise unchanged throughout Newton.
Native MPI1/2 full **4270 coefficient** trajectory
relative difference was `3.235700308035909e-16`;
Hamiltonian relative difference
`2.529867524588143e-15`.
All five sample timestamps `[0,1e-4,2e-4,3e-4,4e-4]`
were actually saved in this experiment, unlike
the earlier three-step pilot which missed
one intermediate saved time.

**Native qualified scope:** real production
nonlinear-state affine-algebraic initialization,
short physical BDF2 integration, nonvacuous
accepted-step explicit gyro freeze, full-state
MPI1/2 parity and measured Hamiltonian.
**Not qualified:** exact BDF2 power balance,
physical energy-defect temporal convergence,
general nonlinear algebraic Newton, or
original rotating viscous Y0.

### Next native gate: fixed T temporal refinement (SOURCE ONLY; PENDING)

To test the *actual* signed defect
`D_H(T)=H(T)-H(0)-integral_0^T P_boundary dt`,
NOT the self-normalized closed-system
ratio that equals 1 whenever
`P_boundary=0` but H drifts, stage
a fixed-final-time original physics campaign:

| dt | Required physical BDF2 steps | MPI ranks |
|---:|---:|---|
| `0.0001` | `4` | `1,2` |
| `0.00005` | `8` | `1,2` |
| `0.000025` | `16` | `1,2` |

All cases use original `h_mesh=0.1`,
exactly `T=0.0004`, original h0,p0
and homogeneous controls.
`sandbox/probe_v130_issue28_shallow_short_physical.py`
now supports `--dt` and `--final-time`
with **unchanged legacy default three-step
runtime**. It records actual accepted TS counts,
signed/relative final energy defects,
mass and all boundary powers,
and **actual** saved sample timestamps.
New runner
`scripts/run-v130-issue28-shallow-energy-refinement.sh`
retains native JSON/logs and MPI rank separation.
New host comparator
`sandbox/compare_v130_issue28_shallow_energy_refinement.py`
validates authentic 6-run evidence and outputs
adjacent refinement ratios and *observed*
log2 estimates, WITHOUT enforcing a
theoretical second-order convergence rate.
Unit tests guard missing or inconsistent
native data without expensive Docker runs.

The original h0 is discontinuous at x=0.5
and the gyro is lagged explicitly: do
not assume the formal smooth-solution BDF2
convergence rate or exact power conservation.
A monotonically smaller absolute energy
defect is a **measurement**, not itself
a complete proof of consistent physical
port-Hamiltonian integration.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-shallow-energy-refinement.sh
```

The runner provides `--compare-only <folder>` to
revalidate existing saved native evidence
without executing GetFEM/PETSc again.
Issue #28 remains OPEN for
differentiated-constraint Y0 on the original
rotating viscous disk and the remaining
mesh-dependent/physical power gates.


## 21 November 2026 — ORIGINAL shallow-water fixed-T energy-defect temporal refinement: NATIVE MPI 1/2 PASS

**Authoritative user-supplied archive inspected:** `v130-issue28-shallow-energy-refine-ae1146b8-wnHbWYHn.zip`.
Source `ae1146b89217e2ebaad5a876a9f755569aeae18e`,
immutable Docker `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.

Six independently logged **real** native GetFEM/PETSc
BDF2 time integrations, using the original published
rectangular `examples/shallow_water.py` (unchanged
weak forms, original `h0,p0`, physical boundary
controls all zero, explicitly lagged gyro, and
production opt-in `nonlinear_affine` initial
projection). The **spatial mesh is fixed**
`h_mesh=0.1`, the **physical final time is fixed**
`T=0.0004`; only temporal resolution changes.
All six reports independently confirm unchanged
prescribed initial data, true accepted step count,
no rejected steps or SNES failures, and original
algebraic initial projection `||F_a||=2.1992e-12`
(MPI 1) / `2.3277e-12` (MPI 2).

The closed-system port-power defect is computed from
the **actual signed** quantity
`D_H(T)=H(T)-H(0)-W_boundary(T)` with
`W_boundary(T)=integral_0^T Σ_i U_iY_i dt=0`
since all four original boundary controls vanish.

| `dt` | Actual steps | MPI1 signed `D_H(T)` | MPI2 signed `D_H(T)` | MPI1 final `|D_H|/H0` | MPI1 max relative mass drift |
|---:|---:|---:|---:|---:|---:|
| `1e-4` | `4` | `-0.0046648074239783455` | `-0.004664807369408663` | `1.474512951914752e-7` | `4.2648979745554295e-15` |
| `5e-5` | `8` | `-0.0011776470091717783` | `-0.0011776469655160327` | `3.7224597073002744e-8` | `3.198673480916572e-15` |
| `2.5e-5` | `16` | `-0.0002905913206632249` | `-0.0002905912660935428` | `9.185388100469784e-9` | `3.198673480916572e-15` |

At fixed `T`, the absolute signed-defect refinement
ratios are `3.9611253522046774` (4 to 8 steps)
and `4.0525883790472506` (8 to 16 steps),
giving **observed log2 orders `1.985910356921035`
and `2.018843648450496`**. Thus the evidence
strongly supports **approximately second-order
temporal decrease of this particular closed-system
energy defect over these dt values**. No BDF2
*exact* energy-conservation claim is justified:
the energy actually decreases for each finite step.

The native MPI rank-1/rank-2 maximum relative
Hamiltonian discrepancy is
`3.449818887630128e-15`, and the maximum
mass discrepancy `3.55408164546286e-15`.
Boundary power is exactly zero throughout
all reports; all original h0,p0 remain
unchanged at initialization.

**Sampling caveat:** 4 accepted steps
saved 5 states, 8 accepted steps saved 9,
16 accepted steps saved **16**, not 17.
The 16-step saved series includes both
`t=0` and `t=T` but is missing at least
one intermediate accepted step (`t=0.000375`).
The comparator validates actual saved timestamps,
not an imaginary uniform sample count. For
nonzero future boundary controls, numerically
integrating power requires the true saved time
intervals (or appropriate per-step accounting).

**Release qualification now passed:** nonlinear-state /
affine-algebraic *production* initial projection,
short true physical BDF2 integration,
nonvacuous accepted-step explicit gyro freeze
(previous run), MPI parity, mass conservation to
roundoff, and **fixed-spatial-mesh temporal
convergence of the closed-system signed Hamiltonian
defect with observed order close to two**.

**Scientific scope limitations:** the original
initial height is nonsmooth, so this short
empirical sequence does not prove a general BDF2
energy-error theorem. It does not qualify
nonzero boundary-work quadrature, mesh-independent
spatial acceptance/convergence, general nonlinear
algebraic Newton, or original viscous **rotating**
high-index traction `Y0`. The previously reported
self-normalized `relative_energy_defect_exchanged=1`
remains tautological in this closed system; use
the **signed `D_H(T)`** and the `H0`-scaled
quantity instead.

**Next substantive priority:** return to original
rotating viscous P2 disk issue #28: at fixed original
`h0,p0`, identify the *reduced discrete* boundary
compatibility constraints and determine whether their
time derivative coupled to dynamic FE equations has
sufficient rank to solve for initial `Y0`. Include
corrected GWFL rotation, full gyro/viscosity and
native MPI 1/2 parity; avoid changing original
physical data or introducing a pseudo TS.
Only then attempt a qualified consistent
high-index initialization.


## 21 November 2026 — #28 rotating viscous disk: first-differentiated discrete-constraint rank pilot (SOURCE ONLY)

The previous native six-run fixed-T shallow-water
`dt`-refinement has passed in MPI 1/2:
closed-system signed energy defect ratios
`3.9611253522` and `4.0525883790`,
observed orders `1.9859103569` and
`2.0188436485` at original h_mesh `0.1`
and final time `0.0004`.
This **empirical** result does not solve the
original rotating *viscous* P2 disk initial
traction. Focus returns to #28 high-index
initialization.

**Original native structural fact**:
`dF_a/dY=0` for the **132 traction-Y DOFs**,
while `dF_d/dY` is nonzero, so a square
ordinary `J_aa` Newton is invalid. The
rotating original `q0=(h0,p0)` is NOT
yet proven to satisfy the *exact discrete*
algebraic boundary constraints, even after
the GWFL orientation correction and curved
geometric P2 refinement. The post-correction
rotating reduced LSQR residual `~4.18e-8`
is not a mathematical nonexistence certificate:
the very same LSQR had a false-negative rest
classification, whereas the exact resting state
was subsequently verified at `~1e-17`.

At a **consistent** stationary or rotating
algebraic initial state, differentiate
`F_a(q,a_r,t)=0` once, with
`q=(h,p)`, `a_r=(e_h,e_p,Y_n,U_n,U)`
and `Y` the **unconstrained** traction.
The candidate first-differentiated block
operator is

```text
             columns: qdot      a_regular_dot      Y0
rows F_d :           M_dd           0             J_dY
rows dF_a/dt :       J_aq         J_a,a_reg          0
```

Its expected size on the original
P2 geometry `mesh=0.2` is `3300 x 3300`:
1452 original differential DOFs,
1716 regular-algebraic DOFs and
132 traction-Y DOFs.
The lower row represents
`J_aq*qdot+J_ar*a_r_dot+F_a,t=0`.
The dynamic right-hand side must
eventually include **all original**
gyro, viscosity, pressure and wall
terms; it is **NOT computed in this pilot**.

New source-only native audit:
`sandbox/probe_v130_issue28_differentiated_rank.py`
and `scripts/run-v130-issue28-differentiated-rank.sh`.
They build the *real* P2 GetFEM mass and
nonlinear tangent blocks and a PETSc sparse
two-row/three-column MatNest converted to
AIJ; perform two independent MUMPS
**synthetic** linear reconstruction tests,
reporting both x-recovery and residual
accuracy. Actual base cases:
(1) *exact resting equilibrium* `omega=0`
with only known algebraic values set to
`e_h=0.5` and `Y_n=-0.5`;
(2) original `omega=0.2` rotating physical
`h0,p0` with the original **zero**
algebraic guess and corrected wall
`Rot`. Both MPI1/2 native runs measure
`F_a` at the linearization point and
restore **all** GetFEM coefficients bitwise.
No pseudo TS, physical TS, artificial
physical forcing, algebraic projection or
actual initial Y0 solve is performed.

**Interpretation is intentionally narrow:**
successful inversion on two synthetic
right-hand sides is **a numerical
regularity diagnostic at that linearization**;
it does not rigorously prove full rank,
the differentiation index, or existence
of a compatible rotating initial state.
If the original rotating `F_a(q0,a0,t0)`
cannot be driven to the required
physically scaled tolerance **without
modifying the prescribed q0**, there
is no consistent DAE start at those exact
discrete data, regardless of the
bordered matrix rank.
Rank inference must therefore always be
reported alongside the native initial
algebraic residual.

**NATIVE RUN NOT YET EXECUTED.**
A numerically singular or inconclusive
bordered operator is a legitimate
diagnostic observation, **not a PASS of
rotating Y0 initialization**.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-differentiated-rank.sh
```

After inspecting native MPI1/2 rank/compatibility
data, design a mathematically defensible
consistent initial traction `Y0` (or
report infeasibility of exactly fixed
`h0,p0`), and only then attempt
a physical rotating TS with power
and mesh refinement. Issue #28 OPEN.


## 21 November 2026 — #28 differentiated bordered P2 viscous disk: NATIVE numerical regularity MPI 1/2 PASS (not Y0!)

**Authoritative native evidence uploaded and inspected:** `v130-issue28-differentiated-rank-f9e72ecc-pNH6VN8J.zip`.
Source `f9e72ecc0a999c76c89b3504461ad3c1c1214f94`, immutable
`scrimp-mpi:v1.3.0-rc3` OCI `1d6d0c90115213298047b28237a7f38087563a35`.
Four separate native tests: `rest` and `rotating`, each MPI 1 and 2,
on genuine P2 geometric disk, mesh step `0.2`, original FE weak forms
and corrected GWFL clockwise rotation. No pseudo/physical TS,
no physical initial traction `Y0`, no modification of initial `h0,p0`.

The matrix

```text
                  q_dot     a_regular_dot     Y0
dynamic rows        Mdd           0           JdY
d(Fa)/dt rows       Jad          Jar           0
```

has `3300 x 3300` coefficients with `1452` differential,
`1716` regular algebraic and `132` Y-traction unknowns.
The actual `J_aY` is identically zero; the actual max
absolute `J_dY` entry is `0.5072940586040593`. Algebraic
mass rows have maximum absolute value `0`.

| Native diagnostic | Rest MPI1 | Rest MPI2 | Rotating MPI1 | Rotating MPI2 |
|---|---:|---:|---:|---:|
| Original `||Fa||_2` at linearization | `1.4763422151112202e-17` | `1.4744250640694402e-17` | `0.250963115021321` | `0.2509631150213209` |
| Original `||Fd(q,zdot=0)||_2` | `2.699640468033266e-16` | `2.699640468033266e-16` | `0` | `0` |
| Maximum relative reconstruction error, two synthetic x | `4.151511441367386e-11` | `7.351843090655861e-12` | `1.7389642798601183e-10` | `2.8475843052613367e-12` |
| Maximum relative synthetic residual | `5.579034333565234e-12` | `2.4889542225783863e-13` | `1.587544485253625e-11` | `1.1058788590821502e-12` |

All eight synthetic MUMPS solves across 4 native tests converged
(PETSc reason `4`). Each result reports
`NUMERIC_REGULAR_AT_LINEARIZATION` and full initial GetFEM rollback.
The exact resting control retains Fa/Fd at numerical zero.

**Critical qualification boundary:** the rotating value
`||F_a||=0.250963115` is the residual of the **original zero
algebraic guess** at prescribed rotating `h0,p0`, not a minimized
compatibility residual. This pilot used *synthetic* right-hand
sides `A*x_ref` for rank diagnostics, **not the actual
dynamic/differentiated-constraint RHS**; the reported rotating
`||F_d||=0` is therefore not the physical force residual after
solving efforts and boundary controls. Numerical reconstruction
of two synthetic vectors is not proof of the differentiation index,
global invertibility or physical consistency of rotating
initial data. In particular, matrix invertibility at an
**incompatible** algebraic starting point cannot create a valid
DAE initial state.

The earlier original fixed-q rotating P2 reduced-LSQR algebraic
candidate reached a much smaller native `||F_a||≈4.18e-8`,
but missed the strict unscaled `2.60963e-9` coefficient
norm tolerance. This may involve numerical conditioning
and curved boundary trace incompatibility; **neither
physical compatibility nor nonexistence has been proven**.
The exact resting solution at machine precision previously
showed why a bare LSQR `INCOMPATIBLE_AT_TOLERANCE` label
is insufficient.

**NEXT:** audit the regular algebraic variables at fixed
original rotating q0 and the reduced boundary compatibility
in *proper physical norms*. Once and only once
`F_a(q0,a_r0)=0` is numerically justified, assemble
the real differentiated constraint RHS and solve for
`qdot0,a_reg_dot0,Y0` using the original dynamics
(gyro, viscosity, boundary law included). Keep q0
unchanged and require MPI 1/2 parity, rollback on
failure, no pseudo-TS, and subsequent physical energy
and boundary-power checks. #28 OPEN.


## 21 November 2026 — #28 ORIGINAL rotating disk: first genuine physical RHS / diagnostic Y0 pilot (NATIVE PENDING)

This follows the passed native P2 MPI 1/2 rank audit
(`3300 x 3300`, 132 traction Y variables, synthetic
MUMPS x-recovery at relative `~1e-10` or better).
**The rank audit did not compute the actual
physical forcing or a compatible rotating
algebraic initial state**. In particular
its `F_d=0` at rotation was evaluated at the
**original zero-algebraic guess**, not at
the physical constitutive efforts and
prescribed nonzero wall velocity.

Next source-only stage:
`sandbox/probe_v130_issue28_physical_rhs_y0.py`,
`scripts/run-v130-issue28-physical-rhs-y0.sh`,
`sandbox/compare_v130_issue28_physical_rhs_y0.py`
with host-only negative regression tests.
It uses the **unchanged original model**
`sandbox/viscous_swe_2d_mpi.py`,
real curved geometric P2 disk mesh `0.2`,
physical initial `q0=(h0,p0)`, corrected
clockwise GWFL convention, unmodified normal
and tangential boundary controls and full gyro
and viscous terms. There are **no accepted
physical or pseudo-initialization time steps**.

- **Known rest baseline (`omega=0`)**:
  initialize only the known analytical
  algebraic efforts `e_h=0.5`,
  `Y_n=-0.5`, with the exact original
  `h0=5,p0=0`, and require real
  GetFEM `F_a,F_d` at numerical zero.
- **Original rotating case (`omega=0.2`)**:
  hold original `h0,p0` and traction
  `Y=0` exactly fixed; compute only
  regular algebraic efforts/control fields
  `a_r0` with the true *rectangular*
  PETSc distributed LSQR, then evaluate
  the actual original `F_a` and
  per-block residuals. The previous
  comparable pilot got `||F_a||≈4.18e-8`:
  this does **not** establish exact
  discrete boundary compatibility.
- At that **projected** regular algebraic
  candidate, evaluate the **actual**
  differential residual `F_d(q0,a_r0,Y=0,zdot=0)`
  rather than taking the meaningless
  zero-algebraic-guess dynamic RHS.
  Independently finite-difference
  `partial_t F_a` within the actual
  wall ramp's flat initial interval.
- Assemble original GetFEM
  `Mdd, JdY, Jad, Jar` and solve
  the actual first-differentiated
  bordered system

  `Mdd qdot0 + JdY Y0 = -F_d`,

  `Jad qdot0 + Jar a_r_dot0 = -partial_t F_a`.

  Confirm the original full nonlinear
  `IFunction` at the resulting
  `qdot0,Y0` and the original
  `F_a(q0,a_r0,0)`. Never silently
  alter original `h0,p0` or the
  prescribed wall control, and restore
  all GetFEM/PETSc/cache state.
- Report actual `||Y0||_2` and
  `||qdot0||_2` as **DIAGNOSTIC
  ONLY**. Until independent,
  physically scaled rotating
  `F_a=0` compatibility is established,
  even a tiny dynamic residual does
  **NOT** qualify `Y0` as consistent.

The archived native evidence comparator
enforces this separation and supports
`--compare-only <saved-evidence-directory>`
to avoid rerunning MPI after a pure
host comparison failure.

**THIS PILOT IS SOURCE-ONLY. Native MPI1/2
execution and all Y0 candidate values PENDING.**
No claim of definite DAE index,
a physical high-index initializer,
a viable rotating trajectory,
or exact energy/power conservation.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-physical-rhs-y0.sh
```

The critical **next gate after seeing this
pilot** remains physical rotating
boundary-constraint compatibility,
including mesh-scaled norms of the
true discrete normal/tangential traces,
before *accepting* any computed Y0.

## 21 November 2026 — #28 true physical Y0 RHS: MPI2 MatNest ordering bug discovered and fixed in source (native PENDING)

**Failed original native diagnostic archive** `v130-issue28-physical-rhs-Y0-64eb6b2b-DB2mNSvT.zip`, source `64eb6b2b76007f4417c22bccd37df4f85e9d1ef5`, immutable rc3 OCI. The **exact rest reference passes**: native `||F_a||~1.47e-17`, `||Y0||~1e-15`, full dynamic residual `~1e-15`. With the **original rotating** P2 physical profile, both MPI ranks have identical original fixed-q initial states and nearly identical LSQR-projected `||F_a||~4.18e-8`, with actual physical `||F_d||~0.08654` at zero derivative and traction. But:

| Native rotating trial | MPI 1 | MPI 2 |
|---|---:|---:|
| Bordered linear KSP relative residual | `3.37e-13` | `1.68e-12` |
| Candidate `||qdot0||` | `6.57793` | `17.41334` |
| Candidate `||Y0||` | `0.0063778555` | `0.460843533` |
| **Original native GetFEM dynamic `||F_d(qdot0,Y0)||`** | **`1.65e-14`** | **`0.398585756`** |

The old runner stopped at the *correct* host safeguard: original GetFEM MPI2 `F_d` did not satisfy the actual PDE. **This is a real diagnostic MPI error**: small PETSc KSP residuals did not mean the correct physical RHS was assigned or the correct traction read.

**Source diagnosis and fix:** the old physical-RHS diagnostic creates a PETSc parallel `MatNest`, converts it to `AIJ`, then erroneously assumes the global converted row and column ordering is the serial block concatenation `[d; a]` and `[qdot; a_reg_dot; Y0]`. In MPI, MatNest block ranges follow its **native IS** permutations and can be interleaved across ranks; using `np.concatenate` for its RHS and slicing the globally gathered solution at fixed contiguous block offsets is invalid. The previous **synthetic** rank pilot was unaffected since it always inverted an internal `A*x_ref` without interpreting physics in GetFEM variable order.

The new native probe reads `matnest.getNestISs()` **before conversion** and uses its original parallel **row IS** to scatter the actual dynamic and differentiated-algebraic RHS. It checks a complete rank-local IS partition and a numerical RHS round-trip after assembly. It uses the **column IS**, plus the original compact subblock PETSc ownership ranges, to reconstruct `qdot0`, `a_regular_dot0`, and `Y0`, rather than naive global slices. It rejects trials if the injected candidate does not satisfy the **original full nonlinear native `IFunction` dynamics** at `1e-8`, even if PETSc's converted-matrix KSP reported convergence. The schema is incremented to **v2** so an old, incorrectly indexed archive cannot be silently relabelled as a success.

The comparator additionally checks **full physical qdot and Y0 coefficients**, not just the norm of Y0, with independent synthetic corruption unit tests for mapping and MPI parity. Original `h0,p0`, gyroscopic/viscous and boundary weak forms and frozen-gyro semantics are **unchanged**. No Y0 state is accepted, no pseudo/physical TS is performed.

**NATIVE RE-RUN PENDING:** `scripts/run-v130-issue28-physical-rhs-y0.sh`, MPI1/2 rest/rotating, from the current source. Only if the new *native* MPI2 original IFunction is below tolerance and true physical-vector parity passes can we qualify the **index-mapping diagnostic**. Even then rotating algebraic `||F_a||~4e-8` is still not a certified constraint-compatible initial DAE state; **no qualified rotating Y0 or time integration yet**.

## 21 November 2026 — #28 true physical-RHS/Y0 MPI 1/2 native gate PASSED (diagnostic, NOT initialized)

**Host evidence:** `sandbox/outputs/v130-issue28-physical-rhs-Y0-d937cdc6-pl0yI36l/`
(archived as `v130-issue28-physical-rhs-Y0-d937cdc6-pl0yI36l.zip`).
**Source:** `d937cdc6855b815fabc4108c4cad76d980271d3a`;
**immutable image:** `scrimp-mpi:v1.3.0-rc3`, OCI revision
`1d6d0c90115213298047b28237a7f38087563a35`.
The host comparator returns `COMPLETE_DIAGNOSTIC_NOT_INITIALIZER`.
All 7 host regression tests passed. All four genuine
PETSc/GetFEM MPI runs (rest/rotating × 1/2 ranks) completed.

Before this run, the probe incorrectly called
`matnest.convert("aij")` without `out=`, which performs a
petsc4py **in-place** conversion and destroys the native Nest IS;
this caused PETSc error 56 when `getNestISs()` was called.
The source fix at `e27c47f` uses
`matnest.convert("aij", out=PETSc.Mat())` and retains
the original Nest plus its row/column IS for parallel physical
RHS packaging and solution reconstruction. A host AST regression
test at `d937cdc` guards against repeating the in-place conversion.

| Real native diagnostic | Rest MPI1 | Rest MPI2 | Rotating MPI1 | Rotating MPI2 |
|---|---:|---:|---:|---:|
| Native projected `||Fa||_2` | 1.47634e-17 | 1.47443e-17 | 4.18508716e-8 | 4.18369529e-8 |
| Physical `||Fd(q0,a_r0,Y=0,zdot=0)||_2` | 2.69964e-16 | 2.69964e-16 | 0.0865426750 | 0.0865426757 |
| Solved `||qdot0||_2` | 1.25651e-14 | 1.25651e-14 | 6.57792793 | 6.57792791 |
| Diagnostic `||Y0||_2` | 2.25798e-16 | 2.25798e-16 | 0.00637785550 | 0.00637785651 |
| **Unmodified GetFEM native dynamic `||Fd(qdot0,Y0)||_2`** | **1.00534e-26** | **3.57982e-30** | **3.37799e-14** | **2.74326e-15** |
| Bordered relative linear residual | 1.56896e-26 | 2.52383e-28 | 1.46053e-13 | 2.85171e-14 |

The rotating LSQR projection used **5,149 / 5,154** iterations
(MPI1/2), with convergence reason 1, and reduced its unprojected
`||Fa||≈0.250963115` to `≈4.18e-8`. The residual is mainly
`e_p` (`4.10174e-8` MPI1, `4.10164e-8` MPI2), followed by
`U` (`≈6.72e-9`). This is **above** the previously imposed raw
fixed-mesh `2.60963e-9` coefficient-norm tolerance. It is not an
independent proof of exact discrete or physically scaled boundary
compatibility; it is also not sufficient to prove incompatibility,
given conditioning and projection accuracy.

The reconstructed full rotating `qdot0` and `Y0` coefficient
vectors agree between MPI1/2 to relative scales `6.00659e-7`
and `6.68769e-8`, respectively; the relative difference in
their `Y0` norms is `1.00769e-9`. Native MatNest RHS
round-trip errors are exactly zero for all four cases.
The prescribed physical `h0,p0` are unchanged bitwise,
and the GetFEM state is restored. No pseudo, BDF initialization
or true physical time step was executed. The previous bad MPI2
physical dynamic residual `0.398585756` has disappeared.

**QUALIFIED:** on this original P2 model and mesh, the true
post-projection physical RHS is packaged and solved consistently
with the native GetFEM/PETSc dynamics for MPI1/2, including
the exact rest control, correct MPI index mapping and MPI parity.
This **does not** prove first-differentiated full rank, generic
nonlinear DAE initialization, a consistent rotating `Y0`,
nor a physically admissible rotating trajectory.

**NEXT GATE (#28 remains OPEN, #6 nontrivial rotation blocked):**
at fixed original rotating `h0,p0`, independently measure actual
boundary normal/tangential velocity mismatch and interior
constitutive residuals using *physically scaled mesh-dependent
L2/dual norms* rather than only raw algebraic coefficient norms;
check residual sensitivity to LSQR stopping tolerance,
conditioning and mesh refinement. Separate numerical projection
error from any genuine incompatibility with the exact original
wall ramp and weak constraints. Preserve curved P2 geometry,
full viscous/gyro physics, fixed `h0,p0`, original controls,
rollback and MPI1/2 comparisons. If algebraic feasibility can be
established independently, only then promote the differentiated
solution toward an actual consistent initial state, followed by
a physical CN energy/port-power validation. Do not introduce
pseudo-time correction or alter the prescribed rotating state
merely to force acceptance.

The host logs still emit a petsc4py `ComplexWarning` when
reading the nominally real vectors into float arrays; a future
strict real/imaginary-value check should establish whether any
imaginary content exists instead of silently discarding it.
Unused PETSc BDF options arise because the probe intercepts the
initializer before TS; neither warning changes the gate result.

## 21 November 2026 — #28 physical constraint / LSQR-sensitivity gate staged (SOURCE ONLY)

**Next scientific prerequisite after the successful native physical-RHS MPI1/2
diagnostic**: distinguish an actual incompatibility of the *prescribed*
rotating momentum from inaccurate algebraic least squares and geometric
discretization. New independent, read-only probe
`sandbox/probe_v130_issue28_physical_compatibility.py`,
comparator `sandbox/compare_v130_issue28_physical_compatibility.py`,
negative controls `tests/test_v130_issue28_physical_compatibility.py`,
and pinned-Docker runner
`scripts/run-v130-issue28-physical-compatibility.sh`.

Using the unchanged `sandbox/viscous_swe_2d_mpi.py`,
`gamma=0`, `initial_profile=tangent`,
`wall_velocity_profile=facet_tangent`,
`geometry_order=height_order=2`, original clockwise GWFL
`Rot`, and `omega=0.2` (plus the `omega=0` rest control),
the probe evaluates **actual GetFEM integration-region scalar
integrals**, not the Euclidean norm of a PETSc RHS. At `t=0`
the original wall ramp is identically one. For
`v_0=p_0/rho` and `t=Rot*n`, it records:

- the domain area `|Omega|` and boundary length `|Gamma|`,
  for native MPI scalar assembly calibration;
- `||v_0.n||_{L2(Gamma)}` and
  `||v_0.t-omega||_{L2(Gamma)}`,
  plus the full `||v_0-omega*t||_{L2(Gamma)}`.
  These are measured at *prescribed* original `q0` **before
  any algebraic projection**; they are therefore independent
  of LSQR conditioning, `e_p` or `Y0`.
- at the analytic resting control, or projected rotating baseline,
  the physically scaled interior
  `||e_p-p/rho||_{L2(Omega)}`,
  `||sqrt(h)(e_p-p/rho)||_{L2(Omega)}`,
  `||e_h-rho*g*h-0.5*(p.p)/rho||_{L2(Omega)}`,
  and boundary `||e_p-U||_{L2(Gamma)}`,
  `||U-U_n*n-omega*t||_{L2(Gamma)}`, `||U_n||_{L2(Gamma)}`.
  Actual native algebraic `IFunction Fa` residual per block
  and `||A^T Fa||_2` are reported separately.

Each squared integral is archived with raw individual MPI-rank
GetFEM scalars and the raw sum. The RMS norm is normalized by
its **corresponding assembled domain/boundary measure**, so it
is invariant whether GetFEM scalar assembly is rank-local or
globally replicated. The host comparator independently identifies
which scalar MPI reduction convention applies by calibrating
area **and** perimeter against the MPI1 reference; it does
not silently sum a replicated global scalar. It checks the
vector/normal/tangent squared-norm decomposition and MPI1/2
physical RMS parity.

Physical mesh comparison: nominal `h=0.2` and `h=0.1`,
both curved P2, rest/rotating and MPI1/2. At `h=0.2`
the original rotating `h0,p0` and zero `Y` are fixed,
and regular algebraic variables are projected twice on the same
actual rectangular native PETSc system:
`LSQR(rtol=1e-12,atol=1e-13,max_it=10000)` and
`LSQR(rtol=1e-13,atol=1e-14,max_it=20000)`.
This is a **stopping-tolerance sensitivity test** that records
`||A^T Fa||`, individual native `Fa` blocks, correction
`||delta_a||` and the relative coefficient change between
the two LSQR tolerances; iteration counts alone are NOT a
spectral condition number estimate.
At refined `h=0.1` the rotating audit measures ONLY the
unprojected fixed-q boundary data, explicitly avoiding
conflation with a more expensive/different LSQR projection.

The harness rejects changed physical data, altered geometry,
nonzero exact-rest residual, bad MPI physical assembly, missing
independent trace diagnostics, and any claimed accepted
rotating `Y0` or time stepping. It intentionally does **not**
assert that mesh refinement must improve the mismatch or that
a small physical RMS implies exact weak DAE compatibility:
those are questions to answer from native evidence.
The physical FE trace diagnostic is not a proof of a unique
constraint-compatible initial state, nonlinear/high-index
initializer, KSP conditioning, or physical CN integration.

**SOURCE STAGED / NATIVE PENDING:** neither native GetFEM
`asm_generic` expressions nor MPI1/2 physical results have been
qualified by this source commit. No modification was made to
SCRIMP's production initializers, original SWE weak forms,
original `h0,p0`, the faceted wall law, or the image.

Run on the existing immutable image:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_physical_compatibility.py
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-physical-compatibility.sh
```

Use `bash scripts/run-v130-issue28-physical-compatibility.sh
--compare-only <saved-evidence-directory>` to replay archived evidence
without restarting expensive native runs. On completion, interpret
actual fixed-q normal/tangential trace errors, LSQR stopping
sensitivity and mesh trend *before* deciding whether rotating
constraints are numerically/physically feasible.
**#28 remains OPEN; nontrivial rotating #6 remains blocked.**

## 21 November 2026 — #28 original curved-P2 physical compatibility: NATIVE MPI 1/2 PASSED (TRACE DIAGNOSTIC ONLY)

Source `c1988939be5bd3836628295fe88ad6d59f385547`,
immutable `scrimp-mpi:v1.3.0-rc3` with OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Host evidence directory:
`sandbox/outputs/v130-issue28-physical-compatibility-c1988939-WqmfWd2x/`;
matching archive `v130-issue28-physical-compatibility-c1988939-WqmfWd2x.zip`.
**All 7 host-side negative/fixture tests passed** and all
8 genuine GetFEM runs (original rest/rotation, geometry P2,
nominal mesh 0.2/0.1, MPI1/2) completed. The host
comparator returned:
`COMPLETE_PHYSICAL_TRACE_DIAGNOSTIC_NOT_INITIALIZER`.

### Fixed original rotating physical momentum: native boundary RMS

| `omega=0.2`, independent GetFEM boundary measure | h=0.2, 3,300 DOFs | h=0.1, 10,272 DOFs | h-refinement factor |
|---|---:|---:|---:|
| `RMS((p0/rho).n)` | 1.3937491763e-5 | 2.0057694048e-6 | 0.1439117912 |
| `RMS((p0/rho).t-omega)` | 4.8200908963e-7 | 2.2071385559e-7 | 0.4579039282 |
| `RMS(p0/rho-omega*t)` | 1.3945824085e-5 | 2.0178764857e-6 | 0.1446939581 |
| `RMS((p0/rho).t)` | 0.1999995260 | 0.1999997794 | — |

Here `t=Rot*Normal` with the unchanged original clockwise GWFL
sign and wall law. MPI1 and MPI2 physically normalized RMS values
agree to rounding for **every** case; both rest meshes give exact
zero normal and tangential defects. Comparison against nominal mesh
halving gives descriptive two-point slopes `~2.797` (normal),
`~1.127` (tangential), `~2.789` (full vector), **not**
established convergence orders. These nonzero **strong**
boundary mismatch integrals are independent of the LSQR projection.
They are consistent with geometric/trace approximation effects
and show that the original rotating `q0` does not satisfy the
imposed trace **pointwise/L2-exactly** at these mesh sizes.
They do **NOT** by themselves prove failure of the *weak finite-
element* algebraic compatibility equations: an L2 nonzero field
can still have zero moments against the discrete test space.

**MPI scalar-integration convention independently checked:**
on all four case/mesh combinations GetFEM returns the same
**already global** scalar integral on both MPI processes.
The comparator correctly classifies MPI2 as
`replicated-global` after calibration against the actual
area and perimeter from MPI1. Never sum those two copies when
performing physical power, boundary or area audits.

### h=0.2 rotating algebraic LSQR sensitivity, actual native F_a

| Native quantity | MPI1 standard | MPI1 tight | MPI2 standard | MPI2 tight |
|---|---:|---:|---:|---:|
| LSQR iterations | 5171 | 5852 | 5101 | 5723 |
| `||Fa||_2` | 4.1794109732e-8 | 4.1573736753e-8 | 4.1929570788e-8 | 4.1576937032e-8 |
| `||A^T Fa||_2` | 5.0372e-11 | 5.4445e-12 | 4.3682e-11 | 5.3671e-12 |
| `||delta_a||_2` | 13.0655182281 | 13.0655182281 | 13.0655182283 | 13.0655182280 |
| relative `delta_a` change from standard | 0 | 3.01744e-7 | 0 | 3.90119e-7 |

Both KSP runs returned convergence reason 1. The smallest
rotating `Fa` remains `~4.157e-8`, **above**
the earlier unscaled `2.60963e-9` raw-coefficient threshold.
Tightening LSQR decreases the normal-equation gradient by
nearly an order of magnitude but only reduces `||Fa||`
by ~0.5%. The dominant tight residual block is `e_p`
(`~4.10055e-8`), followed by `U` (`~6.704e-9`).
This suggests a nonzero discrete least-squares floor,
but neither proves impossibility of exact weak algebraic
compatibility nor establishes the spectrum/conditioning
of the rectangular problem. Raw `Fa` Euclidean
coefficient norms must **not** be compared across
meshes without a consistent physical/dual scaling.

At the MPI1 tight projection, separately measured *strong*
physical RMS were `||e_p-p/rho||_Omega/√|Omega| =
1.36311987e-7`, weighted velocity constitutive RMS
`3.04802869e-7`,
`||e_p-U||_Gamma/√|Gamma|=6.10108e-10`,
`||U-U_n*n-omega*t||_Gamma/√|Gamma|=
1.39334856e-5`, `||U_n||_Gamma/√|Gamma|=
3.63294e-11`, and pressure constitutive RMS
`||e_h-rho*g*h-|p|²/(2rho)||_Omega/√|Omega|=
5.11671926e-4`.
The latter strong pressure error is **not**
itself evidence of an inconsistent *weak* constitutive
equation: the non-polynomial kinetic expression need
not belong to the P2 co-energy space even when its
Galerkin projection satisfies the test equations.

**QUALIFIED:** immutable native #28 diagnostic probe,
true original curved P2 strong boundary mismatch
and constitutive measures, calibrated MPI scalar
assembly, exact-rest positive control and MPI1/2
parity on both meshes, read-only original `h0,p0`.
There were no BDF/pseudo/physical time steps, no
accepted traction `Y0`, and no changed wall law.
The prior first-differentiated true-physical-RHS
MPI diagnostic remains separately qualified.

**NOT QUALIFIED:** weak discrete compatibility
`Fa(q0,a0)=0` for original rotation, unique consistent
`Y0`, differentiation index, physical rotating
time evolution or its Hamiltonian/port-power balance.
Original #28 and nontrivial rotating #6 remain OPEN.

**NEXT DISCRIMINATING GATE:** independently measure the
*moments of the prescribed wall mismatch* against the
actual GetFEM boundary test space, using a verified
MPI-consistent boundary Riesz mass matrix / dual norm
rather than inferring weak inconsistency from a
strong L2 mismatch. In parallel, check a *properly
scaled* h=0.1 regular-algebraic projection to separate
the stationary LSQR residual floor from the
geometry-dependent mismatch. Compare nonlinear
`IFunction Fa` blockwise after both projections,
including the interior co-energy weak projection.
Only if both algebraic feasibility and first-
differentiated consistency are independently verified
may a diagnostic `Y0` be promoted toward a DAE
initializer. Do not alter original `p0`, geometry
or imposed physical wall velocity to obtain a PASS.

## 21 November 2026 — #28 direct *weak* U-trace obstruction diagnostic staged (SOURCE ONLY)

The successful native **strong** physical-compatibility gate
`v130-issue28-physical-compatibility-c1988939-WqmfWd2x`
does **not** by itself establish weak finite-element inconsistency.
A nonzero `L2(Gamma)` strong mismatch can be orthogonal to the
actual discrete boundary test space. The tight LSQR residual
`~4.157e-8` provides circumstantial evidence only; it is
not a rigorously scaled or independently corroborated
minimum-distance certificate.

New read-only native diagnostic:
`sandbox/probe_v130_issue28_weak_wall.py`,
`sandbox/compare_v130_issue28_weak_wall.py`,
`scripts/run-v130-issue28-weak-wall.sh` and host
negative controls `tests/test_v130_issue28_weak_wall.py`.
It retains the original `h0=5,p0` tangent,
`omega=0.2` (with exact rest `omega=0` positive control),
P2 geometry/P2 unknowns, `gamma=0`, physical
`wall_velocity_profile=facet_tangent`, and pinned
immutable rc3 OCI. It runs on original nominal meshes
`h=0.2/0.1` and MPI1/2 with **no physical/initialization TS**,
LSQR or accepted traction.

At the *exact prescribed initial FE momentum*, assemble
with the actual GetFEM **order-one test** `Test_U`:

```text
b_full(v) = int_Gamma (p0/rho - omega*Tangent).v
b_n(v)    = int_Gamma ((p0/rho).Normal)*(Normal.v)
b_tau(v)  = int_Gamma (((p0/rho).Tangent)-omega)*(Tangent.v)
```

The probe checks `b_full=b_n+b_tau` in the actual vector
of U test-function moments. It also extracts the assembled
native `J_UU` block from the original `IFunction` tangent.
The physical control law contains `-U.Test_U`, so its
negative is the native boundary `U` mass matrix
`M_U=-J_UU`. Solve `M_U x=b` with MUMPS and report the
**discrete Riesz-projected physical boundary RMS**

```text
||Pi_U(p0/rho - omega*Tangent)||_{L2(Gamma)}/sqrt(|Gamma|)
  = sqrt(b^T M_U^{-1} b / |Gamma|).
```

For the exact rest state both strong and weak mismatch
should be numerically zero. In rotation, a nonzero
Riesz norm demonstrates nonvanishing moments against the
actual U test space and rules out the special explanation
that *all* of the original strong wall mismatch lies in
its discrete orthogonal complement. It still does NOT
alone establish full algebraic DAE incompatibility:
that conclusion additionally needs a verified statement
that exact interior `F_ep=0` enforces `e_p=p0/rho`,
exact `F_Y=0` enforces the matching U trace,
`F_Un=0` fixes `U_n=0`, and that the corresponding
FE/test/trace mass maps are nondegenerate, without any
unmodelled kernel.

The native comparator checks archived full first-order
U moment vectors, MPI rank assembly/reduction modes,
normal/tangential recomposition, serial-vs-parallel
coefficients, positive boundary mass Riesz solve,
Riesz RMS `<=` independently integrated strong RMS,
exact rest, bitwise read-only rollback, and strict
**NO Y0 INITIALIZATION**. It records diagnostic
moment size and mesh trend, never claims exact
DAE admissibility from a numerical threshold alone.
Host comparator fixtures are synthetic and must
not be confused with native quantitative results.

**SOURCE ONLY / NATIVE PENDING:** do not claim this
order-one GetFEM form, Riesz block, or MPI semantics have
yet been tested on the host. No change to production
SCRIMP solver, original weak bricks, initial state,
controls, or published release guarantees.

Run next native qualification:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_weak_wall.py
bash -n scripts/run-v130-issue28-weak-wall.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-weak-wall.sh
```

If the Python host comparator alone later needs a replay,
`bash scripts/run-v130-issue28-weak-wall.sh --compare-only
<saved-evidence-directory>` does not repeat GetFEM/MPI work.
Preserve original rotating `q0` and no physical time
integration: #28 and nontrivial rotating #6 remain OPEN.

## 21 November 2026 — #28 original weak U boundary moments: NATIVE MPI1/2 PASSED; exact elimination gate STAGED

**Native source:** `538f7e92bcda072594e6e13409a964c4bb4a415a`,
immutable `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Archived host evidence:
`sandbox/outputs/v130-issue28-weak-wall-538f7e92-86jIvI0n/`.
All six host regression tests passed and all eight native cases
(rest/rotation × P2 mesh nominal 0.2/0.1 × MPI 1/2) passed,
followed by the independently archived comparator:
`COMPLETE_WEAK_WALL_DIAGNOSTIC_NOT_INITIALIZER`.

The original `p0`, `h0=5`, `gamma=0` and prescribed
`U=U_n*n+omega*Rot*n` with `U_n=0` are unchanged.
For the actual GetFEM port-U functions `Test_U`, integrate
`b(v)=int_Gamma (p0/rho - omega*Tangent).v`.
Use the native `M_U=-J_UU` to measure
`sqrt(b^T M_U^{-1}b/|Gamma|)` (full and normal/tangential
components). No accepted Y0, no pseudo/BDF or physical time step.

| Native original rotating `omega=.2` | h=0.2, U=132 DOFs | h=0.1, U=252 DOFs |
|---|---:|---:|
| Strong wall mismatch RMS | 1.394582408498e-5 | 2.017876485655e-6 |
| **Native weak Riesz RMS** | **5.867041664403e-7** | **2.221552966076e-7** |
| Normal weak Riesz RMS | 2.598473916159e-7 | 1.969987243180e-8 |
| Tangential weak Riesz RMS | 4.819480885555e-7 | 2.207136474009e-7 |
| Full GetFEM wall moment coefficient `||b||_2` | 3.679724232681e-7 | 1.276350903720e-7 |
| `h0*||b||_2` expected Y-equation coefficient defect (conditional) | 1.839862116341e-6 | 6.381754518600e-7 |

Both MPI configurations have identical native vector moments to
printing precision; the parallel native order-one vectors and
scalar perimeter are **replicated globally**. The independent
full/normal/tangent moment decomposition errors are
`2.2213e-17` (coarse) and `1.3382e-17` (fine).
Riesz mass solves have relative residuals `<=9e-23`.
At rest, strong and weak residuals vanish exactly on both meshes.

**Qualified conclusion:** the **original prescribed rotating
momentum has a nonzero weak mismatch against the actual finite-
element wall U test space**, not only a strong L2 mismatch
orthogonal to that space. Its Riesz RMS falls only by a factor
`~2.641` between these nominal meshes; the weak tangential
component dominates on the fine mesh. The initial expression
`omega*(y,-x)/(sqrt(x*x+y*y)+1e-6)` includes an intentional
`1e-6` center regularization that also slightly changes the
trace near `r=1`; the contribution of this regularization
versus geometric interpolation has **not** been isolated.
Two mesh sizes do not establish asymptotic convergence.

**Important remaining logical step:** a nonzero `b` alone
does not prove the **full original algebraic system** inconsistent.
At constant `h0=5`, the original `p` and `e_p` FEM are
identical P2 spaces, so the interior weak constitutive
constraint `F_ep=0` should force `e_p=p0/rho` by the
nondegenerate h-weighted mass. Original `F_U=0` should
force the port U coefficients via its boundary mass
`M_U=-J_UU`. If `F_Y` in this unique joint solution is
`+/-h0*b != 0`, exact simultaneous satisfaction of
`F_ep,F_U,F_Y` is impossible for that fixed rotating q0
on the tested meshes. This is a **conditional finite-dimensional
block elimination argument**, not a generic DAE index theorem.

To verify this chain independently in the actual **unmodified
GetFEM IFunction**, new source-only diagnostic
`sandbox/probe_v130_issue28_exact_weak_elimination.py`
copies original `p0/rho` FE coefficients directly to
`e_p`, solves the native `F_U` with `-J_UU` (MUMPS),
checks native `F_ep=0`, `F_U=0`, and compares actual
`F_Y` coefficient by coefficient to the independently
assembled signed `h0*b`. Its companion
`sandbox/compare_v130_issue28_exact_weak_elimination.py`
checks MPI1/2 parity, original fixed q0, the previously
frozen full `||b||_2` values, positivity of the rotating
defect, and zero-rest positive controls. Paired host
negative tests:
`tests/test_v130_issue28_exact_weak_elimination.py`.
The pinned-Docker runner is
`scripts/run-v130-issue28-exact-weak-elimination.sh`.
Source changes do not touch the production initializers
or original PDE. **SOURCE STAGED / NATIVE NOT YET RUN**;
the block identity remains a hypothesis pending this gate.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_exact_weak_elimination.py
bash -n scripts/run-v130-issue28-exact-weak-elimination.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-exact-weak-elimination.sh
```

Even if the three native blocks prove inconsistent at original
rotating q0, do **not** change the original prescribed
initial momentum or physical wall control to force acceptance.
This evidence would justify explicitly rejecting that
original incompatible DAE initial condition and designing a
**separate, explicitly specified compatible physical state**
or a wall spin-up from validated rest for #6.
Never accept physical rotating Y0 or run a physical CN/energy
campaign until a consistent initial DAE candidate is actually
verified. #28 remains OPEN; nontrivial rotating #6 remains blocked.

### Analytical circular-wall limit of the original regularized tangent profile

There is an additional **exact continuous-data observation**, distinct
from the unqualified discrete FE elimination. The original prescribed
momentum is
`p0/rho = omega*(y,-x)/(sqrt(x*x+y*y)+epsilon)`, with
`epsilon=1e-6`. On the **ideal unit circle** (`r=1`,
`Tangent=(y,-x)`), the imposed wall is `omega*Tangent`,
so the original continuous initial trace differs by

```text
(p0/rho - U_wall)|_{r=1}
  = -omega*epsilon/(1+epsilon) * Tangent,
omega=.2, epsilon=1e-6
  => tangential deficit = 1.999998000002e-7.
```

Thus geometric P2 refinement alone cannot remove this exact
regularization-induced tangential mismatch **in the continuum
reference expression**: even a perfect circle retains it.
The measured fine-mesh native tangential weak RMS
`2.207136474e-7` is consistent with that scale, but
its decomposition into regularization versus FE interpolation
is not yet certified. The full original FE weak algebraic
incompatibility is still being checked independently by
the staged exact `F_ep/F_U/F_Y` block gate. Do **not**
silently alter `epsilon` or `p0` in the original
validation campaign; an explicitly new compatible-data
experiment must be separately labelled and qualified.

## 21 November 2026 — #28 exact F_ep/F_U elimination: native 6/8 passed; h=.1 floating-roundoff false rejection fixed

**Archived partial native campaign:**
`sandbox/outputs/v130-issue28-exact-weak-elimination-a28aee14-p8XVIVz0/`,
original source `a28aee142b2d1798e225535f9dab217a91938a0c`,
immutable `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
**Six of eight native results completed and are preserved**, all
six host-only comparator fixture tests passed. The original
analytic resting state passes MPI1/2 at both mesh 0.2 and 0.1.
The original rotating P2 mesh 0.2 passes native MPI1/2:

| Native diagnostic at fixed original rotating h0=5,p0, zero Y | h=.2 MPI1 | h=.2 MPI2 |
|---|---:|---:|
| `||F_ep||_2` after exact coefficient copy and U elimination | 5.91914e-18 | 5.91914e-18 |
| `||F_U||_2` after unmodified wall mass solve | 2.75512e-17 | 2.91411e-17 |
| **`||F_Y||_2`** | **1.839862116389e-6** | **1.839862116390e-6** |
| **Independent `h0*||b||_2`** | **1.839862116341e-6** | **1.839862116341e-6** |
| `||F_Y-h0*b||_2` | 1.58918e-16 | 1.70371e-16 |

The *seventh* native case, original rotating **h=0.1
MPI1**, computed the expected nonzero
`h0*||b||_2=6.38175452e-7` and
an independent `F_Y` agreeing **to
`1.09940161e-16` coefficient-norm error**, but
the old source incorrectly rejected the identity at
`error > 1e-10 * max(1e-6, ||h0*b||)`,
which is precisely `1e-16` on this small load.
The error exceeded the fixed arbitrary threshold by
`9.94e-18`; there is **no PETSc/GetFEM/MUMPS solver
failure** and no evidence of a different physical
solution. The seventh run did not write a JSON report;
the eighth, h=.1 rotating MPI2, did not start because
the shell runner correctly stopped on the false rejection.

**Source-only rounding fix:** the independent native
`F_Y=+/-h0*b` comparison now uses a
strict floating-roundoff-aware bound

```text
||F_Y - sign*h0*b||_2 <=
    5e-15 + 2e-9 * ||h0*b||_2
```

and archives the *actual* absolute and relative
coefficient errors and the tested tolerance. The
host comparator independently enforces this same
bound and remains backwards-compatible with the
six original native reports, which were verified
against the new bound *from their original saved
coefficient vectors*, not fabricated. Two new
host negative controls distinguish a valid
`~1.10e-16` rounding discrepancy from
an intentionally corrupted `2e-12` identity.

**Recover, do not repeat:** the native runner
`scripts/run-v130-issue28-exact-weak-elimination.sh`
now accepts `--resume <existing evidence directory>`.
It checks the immutable OCI revision, the exact
old/current source provenance, the existing original
manifest and its no-time-step/no-Y0 invariants;
each old report is independently validated before
it is **skipped**. Only the two missing h=.1 rotating
MPI1/2 native cases are rerun, with current source.
The original manifest `source` remains untouched;
`resume_metadata.json` separately records the
original versus resumed source commits. No
numerical result is silently relabelled as a
new native execution.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_exact_weak_elimination.py
bash -n scripts/run-v130-issue28-exact-weak-elimination.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-exact-weak-elimination.sh \
  --resume sandbox/outputs/v130-issue28-exact-weak-elimination-a28aee14-p8XVIVz0
```

**NATIVE RE-RUN PENDING** for the last two cases.
The pre-existing h=.2 rotating evidence provides
strong native support for **incompatibility at
fixed original q0** of the original
`F_ep,F_U,F_Y` block triple: exact P2
`e_p=p0/rho` and original control-port U
elimination leave a nonzero independently
verified weak `F_Y=h0*b` remainder.
The continuum regularized original tangent
profile has the independent nonzero
`~2e-7` tangential wall deficit even on
the ideal circle. These facts **do not**
qualify a general DAE index, a physical Y0
for this inconsistent initial state, nor a
time-integrated rotating trajectory.
Do not alter original h0,p0/epsilon,
wall law, geometric P2, or run any TS.
#28 and rotating #6 remain OPEN until the full
native archive and an explicitly compatible
initial-condition strategy are reviewed.

## 21 November 2026 — #28 native 8/8 qualified: exact weak block obstruction at the prescribed rotating q0

**FINAL NATIVE GATE PASSED:** original finite-element curved-P2 disk, original
`h0=5`, original regularized `p0/rho = omega*(y,-x)/(r+1e-6)`,
original imposed rotating wall `U=U_n*n+omega*Rot*n` and
`U_n=0`, `omega=.2` (plus exact-rest `omega=0` control),
nominal mesh `h=.2/.1`, MPI1/2, immutable
`scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Completed mixed-provenance (strictly preserved) native evidence
`sandbox/outputs/v130-issue28-exact-weak-elimination-a28aee14-p8XVIVz0/`.
Six validated native reports originate from source
`a28aee142b2d1798e225535f9dab217a91938a0c`;
the two `h=.1` rotating reports originate from resumed source
`9e4dd40ec36482484641ae35fc510df9bf3c6e7b`,
with **only** the absolute/relative roundoff guard corrected.
The original manifest was not relabelled, and
`resume_metadata.json` records the resumed source.
All **eight** native reports now pass the independent host comparator,
and all **eight** host regression/negative controls pass.

### Certified original weak block elimination (fixed q0)

At constant `h0=5`, `p` and `e_p` use the same P2
FE basis and the `h0`-weighted positive-definite interior
mass. Therefore the native interior constitutive block
`F_ep=0` uniquely implies `e_p=p0/rho`;
this was implemented by direct coefficient copy, **without
altering original p0**. Solving the original native
`F_U=0` equation using its boundary mass `M_U=-J_UU`
fixes U uniquely on the filtered U trace space. The
original `Y`-row residual after these two eliminations
is then independently identified coefficient-by-coefficient
as the unchanged weak wall-mismatch functional:

```text
b(v) = int_Gamma (p0/rho - omega*Rot*Normal).v,
F_Y = sign * h0 * b,  sign validated from actual IFunction.
```

| Original diagnostic | h=.2 rotating MPI1 | h=.2 rotating MPI2 | h=.1 rotating MPI1 | h=.1 rotating MPI2 |
|---|---:|---:|---:|---:|
| `||F_ep||_2` | 5.91914125e-18 | 5.91914125e-18 | 3.21638724e-18 | 3.21638724e-18 |
| `||F_U||_2` | 2.75511808e-17 | 2.91411181e-17 | 2.10639674e-17 | 2.14560291e-17 |
| **`||F_Y||_2`** | **1.83986212e-6** | **1.83986212e-6** | **6.38175452e-7** | **6.38175452e-7** |
| independent `h0||b||_2` | 1.83986212e-6 | 1.83986212e-6 | 6.38175452e-7 | 6.38175452e-7 |
| `||F_Y-sign*h0*b||_2` | 1.58918079e-16 | 1.70370660e-16 | 1.09940161e-16 | 1.15389289e-16 |

All four rest/MPI/mesh control cases have
`F_ep=F_U=F_Y=0` exactly. No numerical iteration
or nonlinear mass approximation was needed for
`e_p=p0/rho`; U was solved with the native
boundary mass. The independent weak wall
moment agrees with the unchanged native IFunction
to absolute `~1e-16`, while the rotating residual
is `~1e-6..1e-7`. The verified
roundoff criterion is
`||F_Y-sign*h0*b|| <= 5e-15+2e-9*||h0*b||`;
the former `1e-16` accidental false failure
at `h=.1` is resolved.

**VALIDATED MATHEMATICAL/NATIVE CONCLUSION (scope restricted
to this original discretization):** for the **fixed
original rotating h0,p0**, there is no simultaneous zero
of the three original algebraic constraint blocks
`F_ep,F_U,F_Y` on the two P2 meshes considered,
given the verified nondegenerate same-basis interior
and boundary mass operators. The contradiction
is **independent of LSQR stopping tolerances** and
is not a PETSc/MUMPS breakdown; it arises from
the prescribed initial data and wall actuation.
Even the ideal geometric unit circle retains
a `~1.999998e-7` tangential mismatch because the
original regularization uses `r+1e-6`.

**NOT qualified:** a global proof for arbitrary meshes
or initial states, general nonlinear DAE differentiation
index, physical rotating `Y0`, physical rotating
CN time evolution, or Hamiltonian/port-power balance.
No accepted initial state, no TS initialization or
physical time step was performed. Production SWE
weak bricks, geometry, physical initial values, and
SCRIMP production initializer were **not changed**.

### Implementation decision after this completed gate

The *scientific* obstruction for the **original** rotating
q0 is now resolved: do **not** continue iterating LSQR
or force a physical `Y0` on inconsistent original data.
The software issue #28 remains **OPEN** until SCRIMP
has a safe, explicit and verified numerical pathway:

1. **Failure semantics:** fail-fast informative diagnostics
   for incompatible prescribed differential data rather
   than accepting an ill-defined algebraic Y0 or an
   initialization transient. Develop a general interface
   and tests using the original rotating model as the
   negative control, the exact rest state as the
   positive control; do not unconditionally hard-code
   SWE-specific tolerances or turn the mass diagnostic
   into an undocumented solver default.
2. **Compatible experiment, deliberately separate:**
   construct explicitly specified physically compatible
   initial data and wall actuation (e.g. fully validated
   rest-to-rotation ramp, or a separately labelled
   constrained FE projection of rotating momentum),
   with exact all-block `Fa=0` **before** a differentiated
   physical RHS/Y0 consistency test. Never overwrite
   original `p0` or smooth wall law silently.
3. **Once accepted:** recheck true physical RHS
   (including `qdot0,Y0`), MPI parity, then validate
   actual physical Crank–Nicolson power/energy
   balances and eventually dissipative rotating
   Saint-Venant #6.

**Release status:** original rotating q0 algebraic
incompatibility is **qualified**; no generic
initialization or physical rotation is qualified.
#28 remains OPEN as an implementation item.
Nontrivial rotating #6 is still BLOCKED.

## 21 November 2026 — staged #28/#6 safe startup gate (SOURCE ONLY)

The original fixed-q0 P2 rotating incompatibility is now natively
qualified 8/8 with exact `F_ep,F_U,F_Y` elimination (see
`docs/V130_DAE_ALGEBRAIC_INITIALIZER_DESIGN.md`).
Next native experiment, **not yet run**, adds
`sandbox/probe_v130_issue28_startup_gate.py`,
`sandbox/compare_v130_issue28_startup_gate.py`,
`scripts/run-v130-issue28-startup-gate.sh`, and
`tests/test_v130_issue28_startup_gate.py`.
Pinned `scrimp-mpi:v1.3.0-rc3` OCI as before.
At original already-rotating q0, opt-in
`DPHS.init_step(... nonlinear_affine ...)` must
raise typed `UnsupportedNonlinearAlgebraicDAE`
for higher-index null algebraic columns with
full rollback and no TS. Separately, the explicitly
new physical rest-to-rotation experiment
(`start_from_rest_duration=.02`,
P2 curved disk, facet-tangential physical
wall, `omega=.2`, `p0=0`,
`e_h=.5,Y_n=-.5,e_p=U_n=U=Y=0`)
must satisfy the **full native IFunction**
at `t0` with `zdot=0` and the original
smoothstep `g(0)=g'(0)=0` must appear
as `O(t²)` in native boundary control.
Runs: two cases × mesh .2/.1 × MPI1/2,
read-only, no Y0/TS acceptance.
The original physical experiment and production
SCRIMP initializers are not modified.
Remaining release gates: new positive t0
native feasibility; differentiated
initial Y0/qdot and real CN
energy/port power. #28 OPEN, #6 blocked.

## 21 November 2026 — #28/#6 native dual startup gate: **8/8 PASSED** (compatible spin-up t0 + production typed rejection)

**Native qualification:** source
`fee8626b4634070e4e20686bc7033d9795351e20`,
image `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`,
archived
`sandbox/outputs/v130-issue28-startup-gate-fee8626b-0LcChumB/`.
All eight host-only negative/fixture tests PASSED; all
eight genuine GetFEM/PETSc probes PASSED
(two cases × nominal mesh 0.2/0.1 × MPI1/2),
and independent comparator returned
`NATIVE_STARTUP_SAFETY_AND_T0_FEASIBILITY_DIAGNOSTIC_ONLY`.
**DO NOT** confuse synthetic test displays
(`Fa=1e-17,Fd=1e-17,null_cols=8`)
with actual native measurements below.

The deliberately **new**, physically distinct open
Saint-Venant experiment (`h0=5,p0=0`,
P2 curved disk, original physical
facet-tangential wall, `omega=0.2`,
`start_from_rest_duration=.02`,
original viscous/gyro bricks) has
the analytic t0 algebraic candidate
`e_h=0.5, Y_n=-0.5, e_p=U=U_n=Y=0`
and `qdot0=0`.
Native full original GetFEM `IFunction`
with `zdot=0` evaluates:

| Positive spin-up t0 | h=.2 MPI1 | h=.2 MPI2 | h=.1 MPI1 | h=.1 MPI2 |
|---|---:|---:|---:|---:|
| `||Fa||_2` | 1.4763422151e-17 | 1.4744250641e-17 | 8.0033824547e-18 | 8.0057171923e-18 |
| `||Fd||_2` | 2.6996404680e-16 | 2.6996404680e-16 | 2.6280928039e-16 | 2.6280928039e-16 |
| `||F_U(.0002)|| / ||F_U(.0001)||` | 3.98662207358 | 3.98662207358 | 3.98662207358 | 3.98662207358 |

The physical rising smoothstep
`g(s)=3s²-2s³` with `s=t/.02`
satisfies `g(0)=g'(0)=0`,
and the **native wall form** ratio exactly
matches its two-time analytical ratio.
The full original `F(t0,z0,zdot0=0)`
is thus numerically zero in **both dynamic
and algebraic rows**, with MPI1/2 parity.
This qualifies a **t0-compatible FE state**,
not a computed high-index index analysis,
nor differentiated `Y0`, nor physical
time integration.

For the **unchanged incompatible original**
already-rotating tangent data (f(0)=1),
the existing opt-in production
`init_step_ts_type=nonlinear_affine`
path was called through the **real**
`DPHS.init_step` dispatch and threw
the expected typed
`UnsupportedNonlinearAlgebraicDAE`
on null algebraic columns:
132 (h=.2), 252 (h=.1),
identical MPI1/2. The native
`initialization_report.status` is
`REJECTED`, all original coefficients
and `GlobalVec` are restored bitwise,
and no pseudo/BDF/physical TS executed.
This is a **structural refusal of a
higher-index class**, NOT a generalized
runtime physical infeasibility diagnosis;
the earlier exact FE elimination is
the independent proof for the original
rotating q0. Importantly the same
structural limitation also applies to
some **physically compatible** states;
the high-index initializer remains
unimplemented.

**Next stringent gate:** verify the
first-differentiated DAE constraint system
at the **new compatible spin-up t0** with
correct true physical `-Fd` and correct
`-partial_t Fa`. Critical numerical
caveat: a naive finite difference of the
C1 smoothstep at its clipped start
`g(0)=g'(0)=0` produces a spurious
nonzero `O(dt)` derivative despite
the exact derivative being zero.
Use verified analytic `partial_t Fa(0)=0`
and an independent one-sided `dt`
convergence check rather than
reusing an unqualified fixed central
difference. Then compare MPI1/2
the genuine differentiated bordered
GetFEM/PETSc solve for `qdot0,Y0`,
and recheck native full `IFunction`
and differential-state preservation
before accepting an initial state.
**NO physical TS until this check passes.**
#28 remains OPEN as an actual high-index
production implementation; nontrivial
physical rotating #6 remains blocked.

## 21 November 2026 — #28/#6 true first-differentiated compatible C1 spinup t0 gate staged (SOURCE ONLY)

**Previously NATIVELY PASSED:** the separately specified
physical rest-to-rotation spinup state (h0=5,
p0=0, curved P2, original facet-tangent wall,
viscosity and gyroscopic/port-Hamiltonian
forms preserved, omega=.2 following
`start_from_rest_duration=.02`) satisfies
the original full `IFunction` at
`t0=0` to ~1e-16 on both mesh .2/.1,
MPI1/2 (8/8 previous dual-startup cases).
The existing opt-in `nonlinear_affine`
initializer separately rejects the old
already-rotating high-index negative control
with 132 / 252 null algebraic columns
and exact rollback. Those are **two
distinct results**. The old rotating fixed-q0
data were earlier proven algebraically
incompatible by exact `F_ep,F_U,F_Y`
elimination and remain untouched.

**New diagnostic now STAGED (native not yet run)**:
`sandbox/probe_v130_issue28_spinup_first_derivative.py`,
`sandbox/compare_v130_issue28_spinup_first_derivative.py`,
`scripts/run-v130-issue28-spinup-first-derivative.sh`,
`tests/test_v130_issue28_spinup_first_derivative.py`.

This is a deliberately **new** physically compatible initial
condition, NOT a hidden projection of the original rotating
momentum. The exact t0 candidate uses
`e_h=rho*g*h0=.5`,
`Y_n=-rho*g*h0=-.5`,
`e_p=U_n=U=Y=0` and `p0=0`.
It does NOT run the production
`nonlinear_affine` initializer,
which correctly declares this
high-index structure unsupported.
It reuses and adapts the already
natively qualified PETSc first-
differentiated physical-RHS MatNest
diagnostic, with the crucial
MPI native block IS remapping.

The first differentiated candidate is computed from

```text
  [ M_dd    0     J_dY ] [ qdot0  ]   [ -F_d(t0,z0,0)  ]
  [ J_ad   J_ar   0   ] [ a_dot0 ] = [ -partial_t F_a ]
                           [ Y0    ]
```

where `a_regular` excludes
the dynamic-only traction `Y`.
GetFEM's actual `M_dd,J_dY,J_ad,J_ar`,
original nonlinear IFunction
and MPI-aware PETSc MatNest/AIJ
are used; independent native
`IFunction(t0,z_trial,zdot_trial)`
must confirm algebraic and
dynamic residuals. No artificial
zero forcing is substituted for
the actual dynamic `-F_d`.

**VERY IMPORTANT C1 time derivative:** at the
new physical control `f(t)=f_decreasing(t)*g(t)`
with `g(t)=3s²-2s³`,
`s=clip(t/.02,0,1)`,
`g'(0)=0` **exactly** and
`f_decreasing(t)=1` for
`0<=t<.03`.
Consequently
`partial_t F_a(t0,z0,0)=0`
analytically for this model.
A naive finite difference with fixed
`dt` at the clipped initial time
would introduce spurious `O(dt)`
time forcing into the differentiated
constraints and a false `Y0`.
The new probe therefore uses
the mathematically exact value zero
**and independently verifies** the
actual native GetFEM `U`-block
time response at `delta=1e-4`
and `2delta` matches the
exact `g(2delta)/g(delta)`
to ~2e-6 relative; all other rows
must remain time invariant under
that frozen original t0 state.

Native checks (planned):
- nominal mesh .2 and .1, MPI1/2,
  pinned immutable `scrimp-mpi:v1.3.0-rc3`
  OCI `1d6d0c90115213298047b28237a7f38087563a35`;
- exact t0 algebraic `Fa<1e-11`,
  physical dynamic forcing below 1e-10,
  original nonlinear/P2 discrete mass
  and stiffness, true `J_dY` block
  and zero algebraic mass rows;
- MatNest global MPI row+column
  index partition and RHS roundtrip,
  direct MUMPS convergence and
  independently verified original
  `IFunction` with candidate
  `qdot0,Y0`;
- candidate norms `qdot0,Y0,a_dot0<1e-8`,
  full original dynamic residual
  <1e-10 and algebraic <1e-11;
- MPI1/2 physical and coefficient
  parity, bitwise preservation of
  differential q0 and rollback,
  no physical or pseudo/BDF TS,
  **no accepted physical Y0**
  and no generic index or
  high-index production claim.

The host fixture tests are synthetic
and no substitute for native GetFEM.
This diagnostic does **not** establish
a globally nonsingular first
differentiated operator for all
models, nor an accepted production
state or a physical Crank–Nicolson
energy/port-power trajectory.

Run:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_spinup_first_derivative.py
bash -n scripts/run-v130-issue28-spinup-first-derivative.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-spinup-first-derivative.sh
```

**Source-only pending native proof.**
Next if this pilot passes: choose
a clearly scoped high-index
consistent-initial-state
implementation strategy, distinguish
it from the existing
`nonlinear_affine` low-index
initializer; validate full rollback
and negative incompatible q0.
Only then start an actual physical
CN step, with Hamiltonian, boundary
port work and viscous dissipation
power-balance checks.
#28 OPEN, #6 nontrivial rotation BLOCKED.

## 21 November 2026 — #28/#6 compatible rest-spinup: physical differentiated t0 MatNest MPI 1/2 **NATIVE 4/4 PASSED**

**Native archive:** `sandbox/outputs/v130-issue28-spinup-first-differentiated-54bafa2c-y1es00lx/`,
original diagnostic source `54bafa2c10b50eb195d92c6a4b00b53970487d5a`,
immutable `scrimp-mpi:v1.3.0-rc3`, OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
All four native runs and the strict full-vector MPI comparator **PASSED**:
two curved-P2 meshes h=.2/.1 × MPI1/2.
The archived native comparator status is
`COMPLETE_DIFFERENTIATED_SPINUP_T0_DIAGNOSTIC_NO_TS`.

The *deliberately new* physical rest-to-rotation experiment is distinct
from the original incompatible rotating q0. Here `h0=5, p0=0`,
the original nonlinear conservative/viscous Saint-Venant weak bricks
are intact, `gamma=0`, and the **physical facet-tangent boundary
velocity** starts from rest with a C1 rising smoothstep
`start_from_rest_duration=.02`; `g(0)=g'(0)=0`.
At `t0` the exact initial algebraic candidate is
`e_h=.5, Y_n=-.5, e_p=U_n=U=Y=0`.

**Actual native values (NOT synthetic host fixtures):**

| Native quantities | h=.2 MPI1 | h=.2 MPI2 | h=.1 MPI1 | h=.1 MPI2 |
|---|---:|---:|---:|---:|
| `||F_a(t0,z0,0)||_2` | 1.476342215e-17 | 1.474425064e-17 | 8.003382455e-18 | 8.005717192e-18 |
| `||-F_d(t0,z0,0)||_2` | 2.699640468e-16 | 2.699640468e-16 | 2.628092804e-16 | 2.628092804e-16 |
| **`||candidate qdot0||_2`** | 1.256514249e-14 | 1.256514249e-14 | 4.150424431e-14 | 4.150424433e-14 |
| `||candidate regular a_dot0||_2` | 1.256514249e-14 | 1.256514249e-14 | 4.150424432e-14 | 4.150424433e-14 |
| **`||candidate Y0||_2`** | 2.257984145e-16 | 2.257984145e-16 | 3.044205831e-16 | 3.044205831e-16 |
| `||F_d(t0,z0,qdot0,Y0)||_2` | 1.005343326e-26 | 3.579807531e-30 | 3.335370174e-24 | 1.139478072e-29 |
| `||F_a(t0,z0,qdot0,Y0)||_2` | 1.476342215e-17 | 1.474425064e-17 | 8.003382455e-18 | 8.005717192e-18 |

Native direct PETSc/MUMPS bordered solve reports converged reason 4 for
all four runs, very small independent bordered residuals
(`1.57e-26`, `2.52e-28`, `3.76e-24`, `2.57e-27`);
full original `IFunction` verifies the resulting physical qdot0/Y0.
Native MatNest block index mappings/RHS roundtrips are exact, MPI1/2
coefficient parity passed, and original differential data and PETSc/
GetFEM state were restored after diagnostic. No accepted Y0 was
published and **no pseudo/BDF or physical TS step** occurred.

The true `partial_t F_a(t0,z0,0)=0` is analytically exact for
this new C1 wall onset: the original decreasing envelope is flat
for `t<.03`, and `g'(0)=0`.
Independently, the actual original GetFEM wall U control at
`delta=1e-4` and `2delta` has native ratio
`3.9866220735785958`, equal to the exact smoothstep quotient;
the other row temporal variation is **exactly zero**.
No spurious fixed-step derivative was substituted.

**Host fixture regression NOTE:** the original eight synthetic host
tests had 7 successes and one **test assertion mismatch only**:
`test_rejects_corrupt_MPI2_solution` perturbed the **scalar
norm** of MPI2 Y0, so the independent comparator correctly rejected
it at `Original native MPI first-differentiated mismatch in Y0`
*before* reaching the intended `vector MPI parity` exception.
This is NOT a native solver/MatNest failure.
A new source-only host fix at `f4c2303` separates
(1) scalar Y0 norm discrepancy from (2) a deliberate MPI1/MPI2
**coefficient permutation at matching norm**, ensuring the full-vector
gate is genuinely exercised. Now nine synthetic host controls
are expected; **host rerun pending**. No native output or numerical
acceptance criterion was modified.

**Qualified narrow conclusion:** the explicitly new physical
rest-to-spinup curved-P2 SWE has an exact t0 FE algebraic candidate
and a compatible first-differentiated *physical* qdot0/Y0 candidate
to near-machine precision for both meshes, MPI1/2.
The original already-rotating q0 remains proven algebraically
incompatible and was NOT altered.
**NOT qualified:** general high-index DAE initializer, generic
differentiation index, first **physical CN** step, sustained rotating
time evolution, port-Hamiltonian energy / boundary-work / viscous-
dissipation balance. The existing production nonlinear-affine
initializer still correctly rejects high-index layouts.
#28 OPEN, nontrivial rotating #6 physical time stepping BLOCKED
pending safe accepted-state publication and full physical CN audit.

**Host regression rerun (no native repeat required):**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_spinup_first_derivative.py
# Optional pure host recheck of UNMODIFIED archived native evidence:
bash scripts/run-v130-issue28-spinup-first-derivative.sh \
  --compare-only sandbox/outputs/v130-issue28-spinup-first-differentiated-54bafa2c-y1es00lx
```

## 21 November 2026 — #28/#6 first actual physical CN step: sandbox-only native gate STAGED

All previous source and native qualification gates are complete:
nine host negative/regression tests PASS and the complete archived
native 4/4 MPI1/2 compatible rest-spinup first-differentiated
physical RHS gate PASS at source `54bafa2` (archive
`sandbox/outputs/v130-issue28-spinup-first-differentiated-54bafa2c-y1es00lx/`).
The first-differentiated candidate for the **distinct compatible
physical spinup** is `qdot0≈0,Y0≈0`, with native residuals
~1e-24 or better; the unrelated original already rotating
`p0` remains proven algebraically incompatible and untouched.

**NEXT GATE — SOURCE STAGED / NATIVE NOT YET RUN**:
`sandbox/probe_v130_issue28_one_cn_spinup.py`,
`sandbox/compare_v130_issue28_one_cn_spinup.py`,
`tests/test_v130_issue28_one_cn_spinup.py`,
`scripts/run-v130-issue28-one-cn-spinup.sh`.

The new sandbox diagnostic uses the real original curved P2
Saint-Venant physics (unchanged SWE weak bricks, physical
facet-tangential wall, `gamma=0,h0=5,p0=0`), original
rotation target `omega=.2` driven by the separately
declared rising C1 smoothstep over `.02`.
Only within a single Python process, a temporary
`DPHS.init_step` monkeypatch publishes the analytically
native-verified t0 effort/traction coefficients
`e_h=.5,Y_n=-.5,e_p=U=U_n=Y=0`.
It evaluates the **true original IFunction** before
publication, demands zero original initial algebraic
and dynamic residual, keeps original q0 bitwise,
and writes exactly the same candidate to both GetFEM
and replicated PETSc `GlobalVec`. It never edits
production `scrimp/dphs.py`, its initializers or
the original incompatible experiment. This
patch is a **sandbox test fixture**, not a
new general high-index initializer.

The original SCRIMP PETSc TS then executes
`ts_type=cn,ts_adapt_type=none`
from `t0=0` to `tf=.002`
with exactly **one accepted step** and no
SNES failure/rejection. The probe fails closed
on wrong actual TS convergence reason, number of
steps, final physical time or missing exact
endpoint records.

**Physical audit before any energy assertion:**

- Re-evaluate the unchanged native GetFEM
  `IFunction` at actual CN midpoint
  `zmid=(z0+z1)/2`,
  `zdot=(z1-z0)/dt`,
  `tmid=dt/2`, independently of PETSc
  KSP/SNES convergence.
- Additionally archive original full and
  algebraic endpoint residuals at `tf`.
  **Subtle:** CN enforces the time-dependent
  algebraic constraint at the midpoint,
  not automatically at the endpoint.
  A nonzero endpoint `Fa(tf)` cannot
  silently be reclassified as compatibility;
  we will evaluate its magnitude and
  time-step sensitivity.
- Reuse the original SWE physical diagnostic
  `sandbox/viscous_swe_2d_mpi.py:diagnostics`,
  which correctly measures the **h-weighted**
  normal/viscous boundary powers and **separate**
  grad/div volume dissipation. The generic
  SCRIMP port power is *not* equivalent to
  these weighted physical SWE terms.
  Archive `H(tf)-H0`,
  `W=int(P_normal+P_viscous)dt`,
  `D=int(D_grad+D_div)dt`,
  `H(tf)-H0+D-W`,
  volume balance, and ratio of defect to
  **actual** exchanged energy rather
  than enormous background H0.
- No pass/fail *convergence* threshold
  is inferred from one dt. Physical energy
  convergence, sustained rotation and
  validated high-index initializer remain
  **explicitly UNQUALIFIED** until
  independent time refinement.

Two meshes nominal h=.2/.1 × MPI1/2,
immutable `scrimp-mpi:v1.3.0-rc3`
OCI `1d6d0c90115213298047b28237a7f38087563a35`.
Host-owned output stores per-rank raw logs,
two full physical endpoint vectors,
time traces of H and power terms and
JSON with original TS reason/steps/false
energy qualifiers. Comparator checks
both native scalar parity and **full**
initial/final FE coefficient MPI parity.
The host negative-test fixtures are
synthetic, not native evidence.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_one_cn_spinup.py
bash -n scripts/run-v130-issue28-one-cn-spinup.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-one-cn-spinup.sh
```

**STOP conditions:** if the real TS/SNES fails,
or if the independent CN midpoint equations
do not match the original IFunction,
do NOT claim any physical evolution.
If TS accepts one step, record the original
endpoint constraint and weighted power defects
without declaring energy conservation based
on a single step. Next gate will use temporal
refinement and possibly an endpoint
DAE-consistency policy that preserves
original physics.
#28 OPEN; physical rotating #6 remains
NOT YET QUALIFIED.

## 21 November 2026 — First physical CN pilot: PETSc TSCN midpoint false gate corrected, nonlinear-mass caution (SOURCE ONLY)

**Original native attempt at `d46c0302ba11a5cd4179998a37b71f5ba0577061`:**
host synthetic regression **9/9 PASS**; pinned
`scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
First native case h=.2 MPI1 published the exact
compatible physical `h0=5,p0=0`,
`e_h=.5,Y_n=-.5` state with
`Fa=1.47634222e-17,Fd=2.69964047e-16`,
then PETSc reported **1 ACCEPTED TSCN step**
`t=.002`, `0` rejected, `0` nonlinear solver
failures, 2 Newton and 2 linear iterations (elapsed ~.76s).
The sandbox then **incorrectly rejected** this actual
PETSc step because it evaluated an **implicit MIDPOINT
method** residual
`F(t=.001, (z0+z1)/2, (z1-z0)/.002)`
and found `||F_mid||_2=203.27337`.
This is a **wrong verification equation**, not
native evidence of failed PETSc TS or of physical correctness.
Only h=.2 MPI1 ran; h=.2 MPI2/h=.1 MPI1/2
were not started. The original incomplete evidence
archive is
`sandbox/outputs/v130-issue28-one-cn-spinup-d46c0302-q5YBlZXK/`
(log and manifest, **no accepted report or
power traces**).

**Authoritative PETSc source analysis**
([TSCN](https://petsc.org/release/manualpages/TS/TSCN/),
[TSTHETA](https://petsc.org/release/manualpages/TS/TSTHETA/),
`src/ts/impls/implicit/theta/theta.c`):
`ts_type=cn` is PETSc **ENDPOINT
trapezoidal theta=0.5**, *not* implicit midpoint.
In the source `shift=1/(theta*dt)=2/dt`,
`stage_time=t0+dt`, and the PETSc SNES
right-hand side is
`((theta-1)/theta)*F(t0,z0,0)=-F(t0,z0,0)`.
SNES thus solves the **actual** original
discrete equation

```text
F(t0+dt, z1, 2*(z1-z0)/dt) + F(t0, z0, 0) = 0.
```

**IMPORTANT: the sign is PLUS**; the earlier
verbal formula using a minus sign was corrected
after checking the actual PETSc source.

The updated native probe checks the **full
coefficient vector** from two independent calls
to the real SCRIMP `IFunction` for
both the dynamic and algebraic rows, with
`||F_TSCN||<1e-7`. It also measures
`||F_a(t1,z1,0)||<1e-7`, records
the old midpoint residual under
`OUT_OF_SCHEME` as a *nonacceptance diagnostic*,
and separately records original
`F(t0,z0,0)` and `F(t1,z1,2Δz/dt)`.
No production integrator or SWE weak form
has been modified.

**Independent methodological caution:**
PETSc's `theta.c` explicitly comments
`This formulation assumes linear time-independent mass matrix`
inside the **endpoint branch**. SCRIMP's
viscous Saint-Venant has actual nonlinear
`M(z)=M(h)`. Therefore even if the
correct PETSc TSCN equation passes,
that only qualifies what the native
*implementation* solved; it **does
NOT** qualify expected time order,
discrete pH power/energy balance,
general nonlinear-mass DAE suitability,
or sustained rotating time evolution.
An independent proper nonlinear-mass
and energy time-refinement study will
be mandatory. The PETSc implicit
midpoint method is a DIFFERENT method
(and per PETSc not stiffly accurate
for a general DAE), so must not be
silently substituted.

**Source-only corrections now staged:**
`sandbox/probe_v130_issue28_one_cn_spinup.py`,
`sandbox/compare_v130_issue28_one_cn_spinup.py`,
`tests/test_v130_issue28_one_cn_spinup.py`.
Tests now accept an arbitrary `203.27337`
*OUT_OF_SCHEME* midpoint residual but
fail on true TSCN endpoint inconsistency,
endpoint algebraic drift, incorrect
method metadata, incorrect boundary
power book-keeping, or MPI vector mismatch.
The previous h=.2 MPI1 failed gate needs
a **fresh native run**; it cannot be
relabelled PASS from PETSc's accepted-step
printout alone.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_one_cn_spinup.py
bash -n scripts/run-v130-issue28-one-cn-spinup.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-one-cn-spinup.sh
```

**Status:** original `t0` compatible state
and first differentiated DAE remain
native-qualified. **No** one-step physical
Hamiltonian/work/dissipation report has
yet passed. Both #28 and #6 remain OPEN.

## 21 November 2026 — native one physical PETSc TSCN step PASSED 4/4; fixed-time dt refinement pending

At source `a8422a9bef09cbd370b14eea8318f09bf2598d9e`
and pinned `scrimp-mpi:v1.3.0-rc3` OCI,
the **new compatible** curved-P2 Saint-Venant
rest-to-rotation C1 wall experiment `h0=5,p0=0`
has now completed **one real accepted** PETSc
TSCN endpoint-trapezoidal step `dt=.002`
on h=.2/.1 MPI1/2, **4/4 native PASS**.
The original IFunction endpoint TSCN residual
is ~4.81e-14 / 7.53e-14 and endpoint algebraic
~2.4e-17 / 1.6e-17 (h=.2/.1);
full original physical FEM coefficients change
(`||Δ(h,p)||=.047/.065`).
The **h-weighted physical** Hamiltonian/work/
grad+div viscosity balance has absolute defects
+6.293e-14 and −2.358e-14, respectively,
or +7.095e-9 and −5.038e-9 relative to
actual exchanged energy, with strict
MPI1/2 parity and volume defect <=1e-13.
Archived native evidence:
`sandbox/outputs/v130-issue28-one-cn-spinup-a8422a9b-2SJgTXLr/`.
The large implicit-midpoint residual is
out-of-scheme and not a PETSc TSCN criterion.

**Still not qualified:** nonlinear-mass TSCN
time accuracy (PETSc theta endpoint
assumes linear constant mass), long-time
spinup trajectory, general high-index
initialization or energy convergence.
Do NOT close #28/#6 on one step.

**SOURCE STAGED / NATIVE PENDING:** full
original fixed-time `t_f=.002`,
`N=1,2,4,8` TSCN steps,
h=.2/.1 and MPI1/2 (16 native cases).
New files:
`sandbox/probe_v130_issue28_cn_refinement.py`,
`sandbox/compare_v130_issue28_cn_refinement.py`,
`scripts/run-v130-issue28-cn-refinement.sh`,
`tests/test_v130_issue28_cn_refinement.py`.
Each step checks true original GetFEM
TSCN/DAE residual; endpoint states
and full h-weighted energy/port/volume
traces are archived. Host comparator
checks MPI full FE vectors and
h,p terminal coefficient differences
across dt; no premature order claim.
Original physical data and production
SCRIMP sources are unchanged.

## 21 November 2026 — #28/#6 NATIVE fixed-final-time PETSc TSCN refinement **16/16 PASS** (t_f=.002, N=1/2/4/8, curved P2 MPI1/2)

**COMPLETE original native evidence**: `sandbox/outputs/v130-issue28-cn-refinement-a9cc2d21-hhcttZc6/`, source `a9cc2d21ead982e5582ae77813a5f0558a1f694f`, pinned `scrimp-mpi:v1.3.0-rc3` OCI `1d6d0c90115213298047b28237a7f38087563a35`. Independently inspected the original uploaded 3.1 MB ZIP: all 16 native JSON and NPZ files, all 16 time/physical coefficient trajectories and original weighted power traces, original manifest and strict `cn_refinement_comparison.json` are present. All 11 host regression/negative tests PASS. Original comparator returns `NATIVE_TSCN_FIXED_FINAL_TIME_REFINEMENT_DIAGNOSTIC_ONLY`.

Experiment is the **NEW** explicitly compatible C1 rest-to-rotation physical SWE (`h0=5,p0=0`, original curved P2 weak physics, original physical facet-tangent wall, target `omega=.2`, rising `g(t)` over .02 seconds, current observation `t_f=.002`). The original incompatible *already rotating* initial p0 and production initializers are **not** modified. Temporary sandbox exact-effort/traction publication from earlier native qualified t0 is used for all runs. Each real PETSc TSCN (endpoint trapezoidal, `theta=.5`) step is independently evaluated using the TRUE original GetFEM IFunction equation `F(t_{k+1},z_{k+1},2(z_{k+1}-z_k)/dt)+F(t_k,z_k,0)=0`; the implicit-midpoint residual remains explicitly out-of-scheme.

**Native MPI1 raw data (MPI2 also PASS):**

| mesh | N | dt | max original TSCN residual | max endpoint Fa | ΔH | original h-weighted W | viscous D | ΔH+D−W | defect / actual exchange |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| .2 | 1 | .002 | 4.81e-14 | 2.41e-17 | 8.85816897e-6 | 8.86996911e-6 | 1.18002064e-8 | +6.29e-14 | +7.09e-9 |
| .2 | 2 | .001 | 9.99e-14 | 2.64e-17 | 8.85815305e-6 | 8.86790635e-6 | 9.75341625e-9 | +1.12e-13 | +1.26e-8 |
| .2 | 4 | .0005 | 1.95e-13 | 3.60e-17 | 8.85814976e-6 | 8.86781724e-6 | 9.66759114e-9 | +1.09e-13 | +1.23e-8 |
| .2 | 8 | .00025 | 3.96e-13 | 4.75e-17 | 8.85814897e-6 | 8.86782170e-6 | 9.67282062e-9 | +9.71e-14 | +1.10e-8 |
| .1 | 1 | .002 | 7.53e-14 | 1.59e-17 | 4.65876210e-6 | 4.68128340e-6 | 2.25212840e-8 | −2.36e-14 | −5.04e-9 |
| .1 | 2 | .001 | 5.26e-14 | 1.43e-17 | 4.65871765e-6 | 4.67735698e-6 | 1.86392535e-8 | −7.69e-14 | −1.64e-8 |
| .1 | 4 | .0005 | 1.07e-13 | 1.89e-17 | 4.65870844e-6 | 4.67719184e-6 | 1.84833457e-8 | −5.62e-14 | −1.20e-8 |
| .1 | 8 | .00025 | 2.11e-13 | 2.52e-17 | 4.65870624e-6 | 4.67720172e-6 | 1.84954129e-8 | −5.92e-14 | −1.27e-8 |

Max endpoint algebraic residual over all native MPI and temporal cases **<4.93e-17**, max TSCN residual **<3.97e-13**, max absolute mass balance defect **<1.60e-13**, and absolute relative power balance defect **<1.65e-8** of the *actually exchanged* physical energy, not of background H0. Independent raw trace reintegration reproduces all weighted wall/viscosity and Hamiltonian balances. MPI1/2 full FE vector parity and identical initial prescribed q0 PASS.

**Temporal original dynamic h,p FE coefficient differences** at fixed final physical time, same mesh:

| mesh | `||q1−q2||` | `||q2−q4||` | `||q4−q8||` | `p(1,2,4)` | `p(2,4,8)` |
|---|---:|---:|---:|---:|---:|
| .2 | 1.11668109484e-4 | 2.79155203799e-5 | 6.97878678223e-6 | **2.0000778803** | **2.0000192900** |
| .1 | 2.35596039991e-4 | 5.88720774443e-5 | 1.47163414212e-5 | **2.0006598472** | **2.0001644850** |

This is compelling evidence of **second-order temporal behavior for this specific early C1 spinup trajectory**, in the original physical dynamic FE coefficients on both tested spatial meshes. These are coefficient norms (not mesh-independent continuum norms) and numerical order estimates (not a general theorem). The power defect is near numerical roundoff, not systematically decreasing with dt, so **no energy-defect order** is inferred.

**Important still-open caveat:** PETSc endpoint TSCN implementation documents its assumption of linear time-independent mass; physical SWE `M(h)` is nonlinear. Empirical early-time order ≈2 does not remove the possible loss of general order/structure over other data, longer times, or generic DAE index. No general nonlinear-mass TSCN theorem/production initializer, nor full spinup/high-index long-time validation is claimed.

**Next physically meaningful gate:** keep original C1 ramp, original nonlinear weak bricks, h0=5,p0=0, but integrate THROUGH the full ramp `t_f=.02` (20 ms), with fixed dt .002/.001/.0005 = N=10/20/40, nominal meshes .2/.1 and MPI1/2. Check true original TSCN and Fa **at every step**, real TS acceptance, full FE coefficient parity, mass and h-weighted physical work+dissipation, and temporal endpoint differences/observed orders. Require positive h throughout; compare energy balance without arbitrary absolute tolerance masquerading as a proof. At t=.02 the rising C1 wall reaches full prescribed rotation `omega=.2` for the first time; the existing decreasing envelope starts only at t=.03, so the physical control remains comparable. **This longer gate is not yet native qualified.**

## 21 November 2026 — #28/#6 FULL C1 SPINUP to t_f=.02: separate 12-run native gate STAGED (source only)

The prior **16/16 native** fixed-time early-spinup TSCN
campaign at `t_f=.002` (source `a9cc2d21`) is fully
qualified as a *specific empirical* second-order
dynamic h,p coefficient convergence study, with
max original TSCN residual <3.97e-13,
endpoint Fa<4.93e-17,
physical h-weighted H+D-W defect absolute <1.12e-13,
and full MPI1/2 FE-vector parity. Original raw
archive `sandbox/outputs/v130-issue28-cn-refinement-a9cc2d21-hhcttZc6/`.

**New scientifically distinct gate (no native result yet):**
advance through the **entire rising C1 physical
wall rotation** by holding `t_final=.02`
instead of .002 and taking `N=10,20,40`
(`dt=.002,.001,.0005`), on each original
curved P2 nominal h=.2/.1 mesh with MPI1/2.
Total **12 planned native configurations**,
pinned original `scrimp-mpi:v1.3.0-rc3`
OCI `1d6d0c90115213298047b28237a7f38087563a35`.
This preserves original physical SWE viscous,
gyroscopic and boundary forms, `h0=5,p0=0`,
original `start_from_rest_duration=.02`
and the earlier native-qualified exact
sandbox-only initial GetFEM/GlobalVec publication.
It never fixes the incompatible already-rotating q0.

At `t=.02`, `g(t)=3s²-2s³`, clipped
`s=t/.02`, is exactly `g=1`;
`g'(t)=0` there. The original descending
envelope remains 1 until t=.03, so the
physical facet-tangent wall reaches full
target `omega=.2` exactly at the
end of the new simulation. New probe
independently re-evaluates original GetFEM
`U` port forcing with original frozen
stationary q0 at t=.02 and t=.02001,
requires `||F_U||>0` and matching
terminal plateau values to relative
<1e-10. This catches an accidentally
absent/mis-timed spinup control.

**Source files staged:**
- `sandbox/probe_v130_issue28_full_spinup.py`:
  original actual per-step TSCN endpoint
  IFunction, Fa(t_k), real accepted TS
  step counts/no SNES failure,
  strictly positive h coefficients,
  complete FE vector trajectory, original
  h-weighted wall normal/tangent work,
  grad/div dissipation, Hamiltonian
  and mass trace, and independent
  original physical wall plateau check.
- `sandbox/compare_v130_issue28_full_spinup.py`:
  strict 12-case provenance, original
  h,p FE intervals, fixed physical
  initial q0, same time horizon,
  weighted work/dissipation reintegration,
  full FE MPI1/2 parity and terminal
  dynamic h,p coefficient time differences
  `||q10-q20||`, `||q20-q40||`
  and **one empirical** observed order
  `log2(E10/E20)`.
- `tests/test_v130_issue28_full_spinup.py`:
  **12** synthetic host regression/negative
  cases; synthetic log numbers do
  not qualify native PETSc.
- `scripts/run-v130-issue28-full-spinup.sh`:
  12 native runs, pinned immutable
  original OCI, source code mounted
  read-only, host-owned evidence/logs,
  fail-closed native comparison.
  No production `scrimp/dphs.py` or
  original SWE weak physics changed.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_full_spinup.py
bash -n scripts/run-v130-issue28-full-spinup.sh

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-full-spinup.sh
```

**Scope:** This next gate tests one **completed
boundary acceleration ramp** but does NOT
cover the original descending envelope
after t=.03, a long-time rotating
Saint-Venant trajectory, general high-index
initialization, proof of temporal order
or discrete-gradient pH structure.
PETSc endpoint TSCN assumes linear
time-independent mass, unlike M(h);
continued empirical success does
not remove this documented caveat.
No native or host PASS has been
reported for the new full-spinup
source until the user runs it.

## 21 November 2026 — Full C1 spinup staging syntax repair; native gate still pending

User's first `scripts/run-v130-issue28-full-spinup.sh` run at source `e2ec401`
**did not execute PETSc at all**: both the host unittest importer
(`sandbox/compare_v130_issue28_full_spinup.py:346`) and the first native
probe (`sandbox/probe_v130_issue28_full_spinup.py:324`) failed Python
parse with an unterminated literal/f-string introduced while cloning the
0.002-second refinement scripts. Docker launched the first case but the
Python interpreter stopped before importing GetFEM or constructing a
simulation; all 12 native configurations remain **UNTESTED**. The
partial output directory
`sandbox/outputs/v130-issue28-full-spinup-e2ec401f-hCnrFHWt/`
is a *failed syntax-gate artifact*, not physical evidence.

The two literals have been repaired without modifying original SWE
physics, boundary control, PETSc TSCN solver, time-grid design, numeric
thresholds, or comparator mathematical tests. To prevent recurrence,
`.github/workflows/python-source-checks.yml` now compiles the newer
#28 native probes and runs the synthetic host comparison tests. A
separate fast workflow,
`.github/workflows/v130-issue28-full-spinup-host.yml`, executes
full-spinup Python compilation, `bash -n` and the 12 host tests
without requiring GetFEM/PETSc. These are source-only safety checks,
**not native qualification**.

Before any Docker rerun, perform the explicit local fail-fast checks:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python -m py_compile \
    sandbox/probe_v130_issue28_full_spinup.py \
    sandbox/compare_v130_issue28_full_spinup.py \
    tests/test_v130_issue28_full_spinup.py
python tests/test_v130_issue28_full_spinup.py
bash -n scripts/run-v130-issue28-full-spinup.sh
```

Only if all pass, run the same immutable
`scrimp-mpi:v1.3.0-rc3` 12-case full-spinup campaign.
The earlier **16/16 native PASSED** refinement at `t_f=.002` and
its empirical order-two result remain valid and unaffected.
#28 and #6 remain OPEN.

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
