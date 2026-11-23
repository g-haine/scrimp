# SCRIMP v1.3.0 — Preparing a compatible controlled rotating SWE state

**Status (23 November 2026):** native experiment **staged**, not yet
executed. This is a concrete initial-state generation workflow for
the *physical, dissipative* 2-D Saint-Venant PFEM example. It is not
a general stationary solver or a generic high-index DAE initializer.

## Why the prescribed originally rotating initial state cannot be used

On the original curved-P2 disk at nominal h=0.2/0.1, the prescribed
`h(0)=5`, pre-rotating tangential momentum and simultaneous moving-wall
constraints **do not** admit a common zero of the original algebraic
blocks `F_ep`, `F_U`, `F_Y`. This is established by the native 8/8
exact weak-elimination gate recorded in
[V130_RELEASE_SPRINT_2026W48.md](V130_RELEASE_SPRINT_2026W48.md).
The failure is not cured by changing PETSc tolerances or selecting a
different Newton solver. Never overwrite prescribed p0 silently.

Instead, begin from the **separately specified and algebraically
compatible** physical state

```text
h0 = 5, p0 = 0, e_h0 = rho*g*h0 = 0.5, Y_n0 = -rho*g*h0 = -0.5,
e_p0 = U0 = U_n0 = Y0 = 0.
```

The complete true original GetFEM `F(t0,z0,0)` is evaluated on all
MPI ranks BEFORE this analytic algebraic state is published. Publication
is a deliberately limited sandbox startup procedure; SCRIMP's
production high-index initialization is **not** modified.

## Separate the tangential spin-up from the normal feedback

The physical parameters in `sandbox/viscous_swe_2d_mpi.py` include
`rho=1`, `g=0.1`, `mu=1e-3`, `omega=0.2`, `h*=4`. The new
generation pilot uses `gamma=0.25`, curved-P2 geometry and the original
dissipative/gyroscopic GetFEM bricks.

Let

```text
r(t) = clamp(t/0.02, 0, 1)
f_up(t) = 3*r(t)^2 - 2*r(t)^3
s(t) = clamp((t-75)/125, 0, 1)
f_decay(t) = 1 - 3*s(t)^2 + 2*s(t)^3
```

The **declared preparation experiment** applies

```text
U_tangential(t) = omega * f_decay(t) * f_up(t) * Rot * n
U_normal(t) = -gamma * (1 - f_decay(t))
              * (Y_n + Y · n + rho*g*h*)
```

(with the original `U_n*n` contribution inside the vector wall law).

The crucial distinction is that **only** the tangential wall actuation
is multiplied by `f_up`. The existing closed-loop normal feedback
retains the original `f_decay`, rather than erroneously using
`f_up*f_decay`. Otherwise, when `gamma=.25`, the original h0=5 and
h*=4 produce nonzero normal forcing already at t=0, violating the
analytically compatible resting t0.

At t=0, `f_up=0` and `f_decay=1`: both wall controls are inactive;
at t>=0.02, `f_up=1`: the rotating wall law **exactly matches** the
original decaying-wall experiment. The normal feedback likewise matches
the original law **at every time**. In particular, its gain factor
`gamma*(1-f_decay)` remains zero for all t<75; this is a preparation
of a state FOR the later controlled experiment, **not yet a nontrivial
closed-loop feedback-response test**. The latter must be qualified at
times when the normal feedback is active.

This source change leaves the previously qualified **gamma=0**
open-loop wall/weak-form expressions unchanged. No extra general
high-index DAE initialization behavior is claimed.

## Reproducible native pilot (on the user's qualified host)

The original image `scrimp-mpi:v1.3.0-rc3` is pinned by its OCI
source revision, and SCRIMP production sources are mounted read-only.
No source rebuild, local output deletion or MPI oversubscription
performance claim is involved. The first pilot runs ten true
Crank–Nicolson steps from t=0 to t=0.02, on one rank:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 bash scripts/run-v130-rotating-controlled-ic.sh \
    seed 0.2 1 0.02
```

The runner prints a **unique host-owned** directory such as
`sandbox/outputs/v130-rotating-controlled-ic-<SHA>-<ID>/`, containing
`accepted_state.npz`, `initial_state_report.json`,
`initial_state_trace.npz`, `native.log` and the source manifest.
The checkpoint is written only after full physical TS acceptance,
finite positive height **coefficients**, per-step true `F_TSCN`,
original algebraic `F_a`, and h-weighted physical
Hamiltonian/work/dissipation and volume acceptance checks.

To extend an already **accepted** state in a fresh PETSc process
and a different MPI partition, preserving the **same** absolute
physical control law and original model fingerprint:

```bash
SCRIMP_DOCKER_NO_SUDO=1 bash scripts/run-v130-rotating-controlled-ic.sh \
    resume 0.2 2 0.024 \
    sandbox/outputs/<actual-seed-directory>/accepted_state.npz
```

Replace the placeholder directory by the one **actually printed** by
the seed run. The original checkpoint is never overwritten. The
continuation writes a NEW accepted checkpoint, at t=0.024 if the
native validation succeeds. Each command is capped at 100 CN steps
to avoid accidental unbounded in-memory trajectories; subsequent
continuations can use the latest checkpoint until the autonomous
wall plateau ends at t=75.

The native launch command tests actual GetFEM/PETSc code. The
`tests/test_v130_rotating_controlled_ic.py` suite and its CI job are
**host-only**, exercising the symbolic physical control schedule,
incompatible t0 regression, time-grid rejection and script structure.
Host-only PASS is **not** evidence of real physics convergence.

## Accepted state is not a stationary rotating state

An accepted full-state checkpoint satisfies the original DAE
constraints at its **absolute physical timestamp**, and can be
reused by `DPHS.load_checkpoint()` on the identical FE model and
future control law. It does **not** imply steady rotating motion.

The pilot records both necessary (but not sufficient) stationarity
diagnostics at the segment end:

- `||F_d(t, z, zdot=0)||_2`: true original dynamic FE residual for
  a hypothetical steady state.
- `||(z_d[n]-z_d[n-1])/dt||_2`: observed dynamic time derivative.

These are diagnostic **raw finite-element norms**, not a mesh-
independent stationarity threshold. They are reported together
with Hamiltonian, h-weighted boundary work and dissipation, volume
and original GetFEM CN/algebraic residuals. The report explicitly
writes `stationary_rotating_state_certified=false` even when an
accepted transient checkpoint exists.

For a defensible *steady rotating* initialization we must first
obtain native pilot results, then refine mesh/time and establish
small, stable residuals/derivatives over multiple windows while
normal feedback and external forcing remain constant. Because
`mu=0.001`, reaching an approximately stationary interior may
require much longer preparation; no unobserved convergence time
is asserted. A controlled response after t=75 requires a
separate native physical test with the original decaying wall law.

## Safety, provenance and release scope

- No change to `scrimp/dphs.py`, PETSc integrators, global FE numbering
  or scientific thresholds of the previously qualified #29 runs.
- No fresh checkpoint until **all** original physical validity gates
  pass; no silent pseudo/BDF projection or corrupted partial output.
- Resume checks the model/control fingerprint and full FE state,
  including algebraic variables; no BDF2/IMEX history is claimed.
- The experimental pilot belongs to #6/#28 scientific follow-up.
  Issue #29's previously completed **cold accepted-state API**
  qualification remains valid, not reinterpreted as steady state.
- Keep independent model/geometry controls: a gamma=0 checkpoint
  is *not* silently loaded into the gamma=0.25 model (different
  fingerprint, even when normal feedback is inactive).
