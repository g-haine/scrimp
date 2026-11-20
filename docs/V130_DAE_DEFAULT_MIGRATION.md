# SCRIMP v1.3.0 — Promotion of algebraic DAE initialization to default

**Status, 20 November 2026: NOT APPROVED YET.** This is the
**default-promotion acceptance plan**, not an assertion that the
new mode is qualified across the whole SCRIMP API. The user has
explicitly requested that the algebraic method become the default
**after qualification**, while preserving access to the legacy
pseudo-time and tiny-step BDF2 methods and documenting their effects.

See [mathematical design](V130_DAE_ALGEBRAIC_INITIALIZER_DESIGN.md),
[Sphinx initialization guide](source/initialization.rst) and
[chronological native science evidence](V130_DAE_INITIAL_VALIDATION.md).

## Explicit choice and non-negotiable behavior

**Phase A — current branch:** keep the historical default
`init_step_ts_type="pseudo"` while integrating and
qualifying an explicitly requested
`init_step_ts_type="algebraic"`. The new mode rejects unsupported
layouts and incompatible prescribed differential initial data.
It never silently runs pseudo-time TS or mutates differential
initial variables.

**Phase B — only after independent native qualification:**
make `algebraic` the default in `DPHS.set_time_scheme` and
in the no-argument fallback in `DPHS.init_step`.
Retain `pseudo` and `bdf` as explicit alternatives. In
particular, existing projects that deliberately rely on
pseudo-time regularization may need to specify
`init_step_ts_type="pseudo"`.

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
