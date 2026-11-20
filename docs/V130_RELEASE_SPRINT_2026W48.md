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
