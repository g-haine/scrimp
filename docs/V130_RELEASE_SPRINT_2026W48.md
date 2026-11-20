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
| **#27** state-dependent mass Jacobian | **Corrected code already default**; native original shallow-water gate **PENDING** | Run `scripts/run-v130-shallow-mass-jacobian.sh`, verify real `examples/shallow_water.py` residual/Jacobian FD vs legacy, MPI 1/2, and archive evidence. |
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
- **#27 (new P0):** state-dependent mass Jacobian default on
  shipped shallow water. Close only with native 1/2-rank
  directional FD and legacy control logs.
- **#28 (new P1):** nonlinear constrained DAE initializer.
  Keep open until genuinely qualified; not a pretext
  to delay/relabel the passed linear default.

**No issue is closed by documentation alone.**
Every closure requires a commit, native test path,
environment provenance and explicit scope/remaining limits.

## Week 48 sequence — no unnecessary repeats

- **Monday 23 Nov — correctness freeze:** run #27 native
  shallow-water mass Jacobian check, pin FD and residual
  default/legacy evidence; triage #6/#23 and confirm no
  unaddressed P0 MPI deadlock or output race. Retain the
  already frozen DAE-default and wave physics results.
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

Use the **unchanged qualified rc3 image** to close the
nonlinear-mass Jacobian question without a long SWE run:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-shallow-mass-jacobian.sh
```

This runner will produce four *short, assembly/FD-only*
runs: production Jacobian default and explicit legacy,
each on MPI 1 and 2. It **must** compare the true
`IFunction` directional derivative with the production
`IJacobian` action at zero and nonzero independent
time derivative and shifts 0 and 1, then quantify
the legacy defect. Source-only CI is not a substitute.
It does not validate nonlinear DAE initialization.

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
