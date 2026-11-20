# SCRIMP v1.3.0 — Promotion of algebraic DAE initialization to default

**Status, 20 November 2026: ALGEBRAIC IMPLICIT DEFAULT
STAGED IN THE SOURCE CODE, NATIVE DEFAULT GATE PENDING.**
Two independent linear pH-DAE models have PASSED
their core opt-in MPI 1/2 physical tests. The
**library default is now algebraic on this
development branch**, but this change itself
is **not yet native-tested or release-approved**.
This document is the **migration acceptance
plan**, not a claim of general DAE compatibility. The user has
explicitly requested that the algebraic method become the default
**after qualification**, while preserving access to the legacy
pseudo-time and tiny-step BDF2 methods and documenting their effects.

See [mathematical design](V130_DAE_ALGEBRAIC_INITIALIZER_DESIGN.md),
[Sphinx initialization guide](source/initialization.rst) and
[chronological native science evidence](V130_DAE_INITIAL_VALIDATION.md).

## Explicit choice and non-negotiable behavior

**Phase A — COMPLETED native opt-in qualification:**
`init_step_ts_type="algebraic"` passes conservative
and dissipative native MPI 1/2 physical tests,
while leaving historical pseudo/BDF2 selectable.
It rejects unsupported layouts and incompatible
fixed differential data without silent fallback.

**Phase B — SOURCE CANDIDATE IMPLEMENTED; native PENDING:**
`algebraic` is now the implicit fallback in both
`DPHS.set_time_scheme` and `DPHS.init_step`
on the development branch. `pseudo` and
`bdf` are explicit alternatives. Historic
models depending on pseudo relaxation must
**explicitly opt in** to that behavior.
The new default is **not yet accepted for
release** until its separate native no-option
regressions and wider migration checks pass.

A *global* algebraic default is **not** safe without auditing
examples and unsupported model classes. A silent fallback to
pseudo is not an acceptable substitute for explicit user choice:
it would invalidate the promise that prescribed differential
variables remain unchanged.

## Preliminary example-impact inventory (source inspection)

| Example | Known relevant feature | Default-migration decision |
| --- | --- | --- |
| `examples/wave.py` | Linear wave but **published Gaussian fixed q0/p0 fails the discrete algebraic gate** at h=0.2, original left control | **Blocker**: make explicit pseudo for the historical published problem, or provide a separately documented compatible initial-data version. Do not silently change the Gaussian or alter its physical meaning. |
| `examples/wave_dissipative.py` | Extended wave with dissipative port; check algebraic classification and initial/control compatibility | Native model-specific projection tests needed; do not assume wave positive control covers dissipation. |
| `examples/wave_coenergy.py` | Substituted coenergy and dissipation ports | Currently explicitly **unsupported** by this initializer; retain explicit legacy strategy until a separate correct implementation exists. |
| `examples/heat.py` | Substituted costate, prescribed initial temperature | Structural/layout qualification needed; do not automatically label every heat formulation supported. |
| `examples/heat_wave.py` | Substituted coenergy variables and interface feedback controls | Currently outside validated linear, non-substituted scope; select explicit legacy or add a separately validated coupled algorithm. |
| `examples/mass_spring_damper.py` | `add_lumped_port` | **Unsupported** by initial algebraic module; historical mode needs explicit opt-in or future lumped support. |
| `examples/nonlinear_pendulum.py` | `add_lumped_port` and nonlinear sine constitutive relation | **Unsupported**; preserve explicit legacy choice and document why. |
| `examples/shallow_water.py` | State-dependent mass/stiffness | **Unsupported** by current linear LSQR method; a nonlinear constrained initializer is a separate scientific project. |

This is a **source-level preliminary classification**, not a
claim that any of these untested cases has passed a native
regression. Perform a source-wide scan of `examples/`,
`sandbox/` and tests that rely on unqualified old defaults
*before* changing `set_time_scheme`'s default.

## Acceptance gates BEFORE editing the default

1. **Core opt-in native correctness:** production
   `DPHS.init_step` (not a monkeypatched diagnostic)
   must reproduce the frozen MPI1/2 positive polynomial
   initial vectors, *all eleven* CN snapshots,
   nonzero signed port work, Hamiltonian balances
   and true algebraic constraint residuals.
   The published Gaussian must be explicitly
   rejected before the physical TS with its original
   state restored. The historical default pseudo
   must reproduce the independently measured
   original physical wave and the tiny negative-time
   BDF2 mode must run with its behavior *measured*
   rather than assumed to preserve initial data.
2. **Second independent compatible linear DAE:** do
   not claim generic linear support based solely
   on one manufactured polynomial wave. Test at
   least another nontrivial linear example with
   a distinct port or boundary arrangement.
3. **Unsupported-class behavior and errors:**
   missing variable intervals, nonzero algebraic
   mass rows, nonlinearity, lumped/substituted ports
   must fail explicitly **without MPI deadlocks**;
   user initial states must remain intact.
4. **Lifecycle, memory and performance:**
   demonstrate bounded PETSc temporary ownership,
   cleanup on success/rejection, no accidentally
   retained pseudo-TS, initialization CPU and
   per-rank peak RSS. The existing global-state
   `Scatter.toAll` remains a documented
   memory-scaling limitation.
5. **Migration of existing examples:** make the
   historical initialization method **explicit**
   where fixed-differential-data projection is not
   appropriate. Never rewrite physical initial
   conditions just to force a green test. Update
   the Sphinx guide and changelog simultaneously
   with the default switch and its compatibility note.
6. **New default native end-to-end test:**
   remove the explicit strategy option in a
   **compatible** new-model test and prove that
   algebraic projection is actually invoked and
   produces identical CN/power outcomes.
   Verify `pseudo` and tiny BDF2 still work
   when explicitly requested; if the mode is
   unsupported or incompatible, require a clear
   actionable exception rather than implicit TS.
7. **Release gate:** no Docker release or v1.3.0 tag
   until the migrated examples and the new
   default have passed their native regression.

## Reproducibility and frozen references

The immutable qualified *diagnostic* reference is from

`sandbox/outputs/v130-petsc-dae-mpi-h0p2-4d3b4d0c-gb0qNKB8/`

and the independent historical pseudo-state physical
comparison is from

`sandbox/outputs/v130-wave-init-power-h0p2-62e31d9b-zDa95p5J/`.

These directories are developer-owned native outputs,
not fabricated GitHub repository fixtures. The new
source-ready core launcher

`scripts/run-v130-core-dae-initial.sh`

requires them as inputs and must be executed on the
user's existing qualified native rc3 image. A passed
source-only CI workflow is **not** a substitute for
the native algebraic-default migration gate.


## 20 November 2026 — completed core gate, explicit legacy migration, second model pending

The production `DPHS.init_step(init_step_ts_type="algebraic")`
**PASSED** the full native reference qualification at
MPI 1 and 2 with the compatible polynomial wave:

- fixed original dynamic states `q0,p0` (both
  coefficient L2 drift **0**) and native true
  algebraic residuals `9.81418431e-14` /
  `9.61877785e-14`;
- complete 11-snapshot CN agreement with the
  immutable independent reference to scaled L2
  error **0** (MPI 1) and **5.30974918e-18** (MPI 2);
- driven signed Hamiltonian+port work defects
  `5.60853295e-15` / `6.05263063e-15`;
- published Gaussian **EXPECTED_REJECTION** at t0
  on both MPI ranks with `F_a=4.82244319e-6`
  versus `2.0689895e-10` threshold;
- original historical pseudo physical trajectory
  identical to its frozen reference;
- tiny negative-time BDF2 `init_step_dt=1e-6`
  **PASS** with measured *nonzero* initial-data drift
  `||Delta q0||=1.9595912757338825e-5`,
  `||Delta p0||=8.680599894760661e-10`.

Native host evidence:

```text
sandbox/outputs/v130-core-dae-h0p2-66fe02c2-WykmP41Y/
```

### Historical examples now request pseudo explicitly

The top-level examples
`heat.py`, `heat_wave.py`,
`mass_spring_damper.py`, `nonlinear_pendulum.py`,
`shallow_water.py`, `wave.py`,
`wave_coenergy.py`, and `wave_dissipative.py`
now specify **`init_step_ts_type="pseudo"`** in
their existing `set_time_scheme` call.
Their published physical equations, initial data
and time steps were **not changed**.
This is a **backward-compatibility preservation**
for the requested future default switch, not an
assertion that their historically relaxed initial
states are DAE-consistent while preserving
their prescribed differential coefficients.
The source-only regression
`tests/test_v130_default_migration.py`
verifies every explicit legacy choice and that
the published wave Gaussian remains unchanged.

**Important:** this source migration does not
establish native runtime equivalence of
all eight historical examples, nor does it
cover every `examples/published/`,
notebook, `sandbox/` or third-party script
using implicit defaults. That wider audit
remains a global default gate.

### Second independent linear port-Hamiltonian model

A native test is staged **but not yet executed**:
the *original* `examples/wave_dissipative.py`
constitutive and Dirac bricks, with a separately
manufactured compatible nonzero polynomial p0,
time-varying boundary control and a distinct
dissipative port. It evaluates the actual
production algorithm on 1/2 MPI ranks, and
tests the complete CN trajectories, actual
IFunction initial constraints, exact dynamic
coefficient preservation, and the **signed
Hamiltonian + boundary work + damping work**
identity. In particular, a dissipative wave
must not be assessed by a conservative-wave
boundary-only balance.

Next native command:

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \\
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \\
bash scripts/run-v130-second-linear-dae.sh
~~~

A 1/2-rank PASS for this second model plus a
remaining audit of default-dependent scripts,
memory/cleanup behavior and an independent
no-option **new-default** run are the next
acceptance conditions. If the dissipative
physics fails, archive the genuine diagnostics
and investigate; **do not weaken the
Hamiltonian/port/dissipation balance or replace
the example's constitutive equations to force
a success**.

**First second-model native attempt (20 November 2026):** the
production algebraic initializer **PASSED** the dissipative
wave's MPI-1 initial constraint gate in two passes, with
`||F_a||=9.99416399e-14`. All ten main CN steps
were accepted, with no rejections or nonlinear solver
failures, for 2,804 DOFs. The experiment then stopped
inside the *diagnostic* GetFEM power integration because
the interior damping port has `region=None`, while
`gf.asm` expects the full-domain integer
`region=-1`. This convention is already implemented
in SCRIMP's `Port.compute()`; the test harness is now
corrected to match it. **No dissipative energy balance
or MPI-2 result is yet qualified.** Original host
evidence:
`sandbox/outputs/v130-second-linear-dae-29aa01dd-ULacuhnd/`.
See section 22 of
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md).

**Default remains `pseudo` until the remaining
native migration gates are satisfied.**

## 20 November 2026 — second model fully qualified, default migration remaining

The dissipative-wave native retry, on the unchanged
qualified `scrimp-mpi:v1.3.0-rc3` image,
**PASSED all required scientific checks**. Host:

```text
sandbox/outputs/v130-second-linear-dae-7eec0172-Ostu9LKD/
```

It reused the original 2D `wave_dissipative.py`
weak/constitutive bricks with its independent
interior `Damping` port. With 2,804 DOFs,
fixed CN `dt=.002` for 10 steps, manufactured
FE-compatible `p0=x*(3-x)` and
`U_L(t)=0.08*sin(20*t)`, the production-owned
algebraic initializer reached true native
`||F_a||=9.99416399e-14` (MPI1) and
`9.93812525e-14` (MPI2), **preserving q0/p0
exactly** in both cases. Internal dissipated
work was **0.0132790026**, external signed
boundary work **-5.40249757e-5**, while
the complete cumulative Hamiltonian + signed
boundary work + internal dissipation defects
were **2.64805539e-15** and
**3.55097896e-15**. Full MPI1/2 CN
trajectory parity reached **6.3717319e-14**
(max scaled L2), well below **5e-8**.
See section 23 of
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md).

**Consequently, the second-independent-model
scientific acceptance gate is CLOSED / PASS.**
The global default should still be promoted
only together with explicit source migration,
a native no-option test, and actual memory/
exceptional-path checks. Changing this default
is intentionally *not* equivalent to treating
historical pseudo/BDF2 initialization as
physically identical: both remain supported
and documented as explicit relaxation choices.

**Source audit scope:** the repository also
contains `examples/published/` research
examples, `examples/notebooks/` and
`sandbox/` scripts which may depend on
the historical implicit default. The top-level
eight example protections cover only those
eight. A global change must explicitly
classify other production/published sources
before deciding which scripts need an
explicit legacy mode.

## Default source migration staged — native gate still required

The source branch has now been updated to use
`"algebraic"` as the implicit initializer
in **both** the initial time-scheme validation
and the actual `DPHS.init_step` dispatch.
The explicit `pseudo` and `bdf` paths
are preserved unchanged. The default does
**not** silently fall back to pseudo on
unsupported or numerically incompatible
prescribed initial data.

In addition to the eight top-level examples,
the source audit identified and protected:
the published `COMPEL_2025/Maxwell_2D.py`,
five historical notebooks (six time-scheme
calls), the legacy GUI-generated scripts
when the GUI has no explicit strategy
selection, the nonlinear/lumped
`sandbox/sloshing_tank.py`, and the historical
`sandbox/bench_acoustic_2d.py` timing baseline.
They now request `pseudo` explicitly to
preserve their previously accepted relaxed
initial states. The explicit list and source
tests are documented in
`tests/test_v130_default_migration.py`.
These migrations do **not** modify the
mathematical model, initial functions, physical
time integration scheme or other physical data.

**Special case:**
`set_time_scheme(init_step=False)` with no
initializer choice must remain legal. The
code tests the presence of an **explicit**
initializer option before rejecting the
contradictory combination
`init_step=False,init_step_ts_type="algebraic"`.
A separate lightweight native PETSc.Options
smoke also checks these semantics.

**Native default acceptance — NOT YET RUN:**

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-dae-default-promotion.sh
~~~

This new runner exercises **the actual absence of
`init_step_ts_type`** on:
- the manufactured conservative polynomial wave,
  MPI 1/2 and its frozen full-CN/power baselines;
- the published incompatible Gaussian, MPI 1/2,
  requiring atomic rejection before physical CN;
- the second dissipative pH wave, MPI 1/2, requiring
  complete `H-H0+W+D` balance, nonzero damping
  and full-state MPI parity;
- explicit historical `pseudo` and tiny-step
  negative-time `bdf` controls, plus
  `init_step=False` native options check.

The independent comparator refuses a passing
result if the initial-strategy keyword was
actually provided in either default-model
test. There is no tolerance relaxation.
All native outputs are saved separately
in a user-owned directory, while the
qualified rc3 image stays unchanged.

**Release still requires:** migration compatibility
tests beyond these smoke models, safe exceptional
MPI cleanup and memory/performance audit,
the source/header/Sphinx checks, and a
documented release compatibility notice.
No Docker release rebuild or tag has
been attempted.
