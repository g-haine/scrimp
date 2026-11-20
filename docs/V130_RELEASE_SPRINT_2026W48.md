# SCRIMP v1.3.0 — release sprint, week 48 (23–29 November 2026)

**Decision recorded 20 November 2026.** The user's delivery objective is
**during the week of 23 November 2026**; use **Friday 27 November**
as the planning checkpoint, with the weekend as contingency rather
than as an assumption that work can run unattended.
This is an **accelerated gate-based plan, not an unconditional release
promise**. The authoritative scope remains
[V130_RELEASE_PLAN.md](V130_RELEASE_PLAN.md), points 0–9.
Do not silently delete a committed gate to make the deadline.

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
| 1 — GetFEM 5.5 simplex FE | **Lean native PASS** | Integrate supported patches consistently into final image; test Maxwell FE spaces. |
| 2 — Docker/native/VS Code | **rc3 native PASS** | Final source revision image/wheel smoke and provenance; do not rebuild on every Python edit. |
| DAE default (new) | **FULL native PASS, two linear pH-DAEs** | Preserve explicit pseudo/BDF2, documented incompatible Gaussian, and short memory/exceptional MPI cleanup gate. |
| **#27** state-dependent mass Jacobian | **FULL native gate PASS (MPI 1/2)**: directional FD plus genuine 4-step physical TS with accepted-state frozen G | Freeze the two evidenced campaigns; close #27 with explicit limitations. No automatic claims about exact energy conservation, order two, nonlinear DAE initialization or 4-rank MPI. |
| **#6** viscous SWE MPI/port powers | **Native qualification PENDING.** Instrumented experimental runner and source-only negative controls prepared; source-check CI PASS at `8a20d551`. | Execute short 1/2/4-rank rest then nontrivial physical open/closed smokes, investigate initialization or rank mismatch on FAIL; only later pursue long-time physics, output races and full issue closure. |
| 3 — New 3D Maxwell | **NOT QUALIFIED** | Minimal scientifically meaningful 3D H(curl)/H(div) example with boundary energy/Poynting power, element orientation and MPI correctness. **Largest new-model schedule risk.** |
| 4 — Performance figures | **Wave evidence PASS; publication figures incomplete** | Curate reproducible plots and metadata; Maxwell scaling only if Maxwell physics qualifies. |
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

- **Monday 23 Nov — correctness freeze:** retain the two
  already native-PASSED #27 gates, close the issue with pinned
  FD and short-TS evidence; prioritize #6 real MPI 1/2/4
  viscous SWE power, deadlock, output-integrity qualification
  and #23 upstream triage. Retain frozen DAE-default and wave
  physics results.
- **Tuesday 24 Nov — model & build gates:** GetFEM simplex
  patch integration and a narrowly scoped 3D Maxwell
  benchmark prototype, plus final candidate-image
  dependency/build context checks. Escalate early if
  Maxwell 3D cannot be physically qualified.
- **Wednesday 25 Nov — compatibility:** representative
  example matrix, small GUI round-trip, native test
  classification and code/header hygiene. Do not
  disguise nonlinear unsupported models as tested.
- **Thursday 26 Nov — reproducibility/docs:** final
  profiling figures, Sphinx HTML/LaTeX, installation
  and version metadata, release notes/README.
- **Friday 27 Nov — go/no-go gate:** full source/CI
  and native smoke, wheel/sdist and candidate Docker
  reproducibility, resolved-issue closure evidence,
  final review/merge/tag **only if all release-blocking
  gates are green**. Otherwise publish the blockers
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
