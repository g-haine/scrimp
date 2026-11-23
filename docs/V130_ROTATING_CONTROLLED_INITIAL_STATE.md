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

## 23 November 2026 — first controlled native run and rejected short-segment power gate

**Original source** `6df826a8`, pinned image
`scrimp-mpi:v1.3.0-rc3`, native user archives
`v130-rotating-controlled-ic-6df826a8-IDSkT9nu.zip` (seed)
and `v130-rotating-controlled-ic-6df826a8-nkyKGqxg.zip` (rejected
resume). Both ZIPs have valid CRC.

**Seed h=0.2, gamma=.25, MPI1, t=0..0.02: genuine native PASS.**
10 accepted TSCN steps, zero TS/SNES failures,
max original `||F_a||=3.242458e-14`,
max original `||F_TSCN||=8.1410743e-11`,
volume defect `3.907985e-14`. Actual port-Hamiltonian data:

```text
Delta_H  = +0.011303178934262359
W_wall   = +0.011521706709641357
D_visc   = +0.000218875523532252
Delta_H + D_visc - W_wall = +3.4774815e-7
relative to the seed's exchanged energy = 3.0182000e-5
```

The full accepted **3300-DOF** state was saved via production
`DPHS.save_checkpoint`. The final original hypothetical stationary
differential residual is `||F_d(z*,0)||=0.04617765493`, with
observed last-step dynamic derivative `27.8368`: clearly **not**
a stationary rotating equilibrium.

**First restart MPI1→2, t=.02..024: numerically integrated but
physically rejected.** The checkpoint passed production
`DPHS.load_checkpoint`; real PETSc accepted both TSCN steps with
zero rejections/SNES failures. Postprocessing reported a volume
defect `1.59872116e-14` and a **relative energy defect
`0.00309672446`**, exceeding the existing per-segment
`1e-3` acceptance threshold. Therefore **NO** new checkpoint,
accepted report or energy trace was written in source
`6df826a8`. Neither the absolute size of this segment's energy
defect nor its equality/difference from uninterrupted integration
can be established from that failed archive.

**Do not relax `1e-3` merely to force a PASS.** In a short interval
after the acceleration plateau, the segment's actual exchanged
energy may be extremely small, inflating `|Delta_H+D-W|/exchange`.
The physical nonlinear-mass CN discretization may also exhibit
truncation-induced discrete chain-rule defects. A proper test must
compare exact original **continuous vs split** native trajectories
and their per-step physical work and dissipation.

### New diagnosis: never lose a failed native physical trace

The follow-up code writes `physical_precheck_report.json` and
`physical_precheck_trace.npz` **before** the unchanged
power/volume acceptance gate. It keeps actual original FE endpoint
vectors, Hamiltonian, volume, integrated boundary work and viscous
dissipation, per-step `Delta_H`, `Delta_W`, `Delta_D`, energy
defects and original CN/Fa residual arrays, **even if the gate
later refuses checkpoint creation**. The report is explicitly a
*diagnostic*, NOT an accepted model state.

Native comparison procedure (on the original user host):

```bash
cd ~/Simulations/scrimp
git pull --ff-only

# Uninterrupted physical control: identical model, 0 -> 0.024.
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-rotating-controlled-ic.sh seed 0.2 1 0.024

# Retry the existing checkpoint as a diagnostic: MPI 1 -> 2.
# Nonzero exit remains possible/expected while energy gate is red.
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-rotating-controlled-ic.sh resume 0.2 2 0.024 \
  sandbox/outputs/v130-rotating-controlled-ic-6df826a8-IDSkT9nu/accepted_state.npz

# After copying the TWO NEW output directory names:
python sandbox/compare_v130_rotating_ic_continuous_split.py \
  --seed sandbox/outputs/v130-rotating-controlled-ic-6df826a8-IDSkT9nu \
  --continuous sandbox/outputs/<new-continuous-directory> \
  --resume sandbox/outputs/<new-restart-directory>
```

The comparator independently reads native source states/power traces,
checks model/mesh fingerprint, source checkpoint bitwise state, terminal
full-FE parity, Hamiltonian/volume continuity, and continuous-versus-
split work/dissipation/balance. Its **PASS is ONLY a restart parity
diagnostic**, not an override of the existing relative physical
energy acceptance threshold. If the same defect appears in both
uninterrupted and split integrations, investigate intrinsic
nonlinear-mass CN energy balance and normalization; if not, investigate
checkpoint/postprocessing clock/history. A subsequent independent
temporal-refinement comparison may be necessary.

The old checkpoint is immutable; failed diagnostic restarts must
not publish unqualified new accepted initial states.

## 23 November 2026 — qualification policy clarified: SCRIMP is responsible for replay equivalence

**Scientific division of responsibilities (user decision):** the question
of *why* a given CN/PFEM nonlinear-mass discretization has a nonzero
discrete Hamiltonian balance defect belongs to the mathematical model
and numerical-analysis research, **not** to this SCRIMP checkpoint
feature. The v1.3.0 restart test must establish that stopping and
resuming produces **the same discrete solution as running uninterrupted**,
including the same numerical balance defect; it must **not** require
that common defect to vanish. No changes to the physical bricks,
discrete model, Newton/MUMPS solver or TS integrator are warranted
by the observed local balance ratio alone.

### Already qualified native example (source `b660dbb2`)

For the *controlled rotating viscous Saint-Venant P2 disk* at
`gamma=0.25`, `mesh=0.2`, `dt=0.002`, pause at `t=0.02`:
the real original seed MPI1 and uninterrupted MPI1 up to `t=0.024`
were compared to a **genuine MPI1→MPI2 production API cold restart**
from `t=0.02` to `t=0.024`. The raw comparator reported:

```text
checkpoint bitwise restored           yes
final full FE state L2 difference     1.646398512702811e-13
initial Hamiltonian restart jump      -7.549516567451064e-15
initial volume restart jump           -7.105427357601002e-15
final Hamiltonian difference         -4.440892098500626e-15
integrated boundary work difference  -3.469446951953614e-18
integrated dissipation difference      3.7947076036992655e-19
continuous full H+D-W defect           7.234183384829951e-07
resumed full H+D-W defect              7.234183429273566e-07
difference of defects                  4.44436154545258e-15
original resumed TSCN residual max     5.5136210170999984e-14
original resumed algebraic max         2.9045221981055017e-15
```

This result demonstrates replay fidelity for **this** mesh,
time window and MPI permutation. The earlier `0.00309672446`
relative power gate over the two-step segment reflected the numerical
model's balance defect normalized by the *small work exchanged only
within that segment*. That number does **not** contradict restart
equivalence and is not the acceptance criterion for this feature.

### Generalizing to full-path equivalence (qualified gate)

`sandbox/probe_v130_rotating_controlled_ic.py` and
`scripts/run-v130-rotating-controlled-ic.sh` now support a
**separate explicitly requested** mode
`SCRIMP_V130_ROTATING_IC_QUALIFICATION=restart_parity`.

- `physical` remains the DEFAULT and preserves the historical
  relative energy-defect threshold, for those explicitly studying
  the time integrator.
- `restart_parity` records *exactly the same* Hamiltonian,
  original GetFEM port powers and dissipation, total and per-step
  defects, and mass/DAE residuals, but **never gates on the absolute
  model-energy defect**. A true PETSc failure, invalid native
  `F_a`/`F_TSCN`, wrong FE/control fingerprint, inconsistent
  physical clock, missing state or volume conservation still fails.
- Each accepted full FE `z(t_n)` is now archived, together with
  original instantaneous normal/viscous port powers and grad/div
  volume dissipation, not only endpoints. A production API checkpoint
  is saved only after numerical state/DAE acceptance.
- `sandbox/compare_v130_rotating_ic_parity.py` compares **every
  accepted full FE state** and physical time against the same
  uninterrupted trajectory, and checks pointwise work/dissipation,
  Hamiltonian and volume. Segment-integrated wall work, dissipation
  and balance are correctly *offset* across the pause, instead of
  erroneously comparing two separately zeroed integrals.
- The model's common `H+D-W` defect magnitude is a *reported
  diagnostic*; its **continuous-versus-restarted difference**
  is a **strict parity gate**.
- Unit tests have negative controls for an altered intermediate FE
  state, an altered boundary work record, a mismatched fingerprint,
  a corrupted checkpoint and an unauthorized pseudo-DAE
  reinitialization. These do not replace actual native runs.

**Native campaign:**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-rotating-ic-parity.sh
```

For **each** h=0.2 and h=0.1, the launcher independently runs
native seed trajectories from rest to `0.02` with MPI1 and MPI2,
uninterrupted `0→0.024` trajectories with MPI1/MPI2 and all
four cross-rank resumes `1→1, 1→2, 2→1, 2→2`. There are 16
native PETSc runs and **8 mandatory full-trajectory comparison gates**.
The campaign retains the host-owned native evidence from each run
and writes a unique `sandbox/outputs/v130-rotating-parity-.../`
directory with a manifest, per-case logs, `*.parity.json` verdicts
and an aggregate `summary.json`. No existing user-owned output
directory is deleted or overwritten.

### Full 8/8 NATIVE PASS — 23 November 2026

The user ran the actual complete native campaign from source
`2119bc551c6cc25a664e59fb2a6962af1a5202ed` under
`scrimp-mpi:v1.3.0-rc3`. Original host evidence directory:

```text
sandbox/outputs/v130-rotating-parity-2119bc55-OtGpBAwi
```

The uploaded ZIP `v130-rotating-.zip` has **valid CRC for all members**,
contains all 16 individual native run directories, their production
`accepted_state.npz` checkpoints and original `initial_state_trace.npz`
arrays, plus the campaign manifest, 8 comparison JSON reports and
`summary.json`. All of these results were additionally checked
directly from the original NPZ arrays independently of the runner.

Each case used original curved P2 disk, `gamma=.25`,
`dt=.002`, physical spinup from rest `0→.02`, and continuous
reference `0→.024` versus cold checkpoint continuation
`.02→.024`. Every native segment accepted all TSCN steps, with
zero rejected steps/SNES failures and native `F_a`/`F_TSCN`
within their independently qualified thresholds.

| Mesh | Full original FE DoFs | Checkpoint MPI | Max full-path FE L2 error | Max discrepancy in H+D−W |
|---|---:|---|---:|---:|
| 0.2 | 3300 | 1→1 | 1.65829742e-13 | 4.47558657e-16 |
| 0.2 | 3300 | 1→2 | 2.24551458e-13 | 1.02157866e-14 |
| 0.2 | 3300 | 2→1 | 2.32118847e-13 | 9.76822789e-15 |
| 0.2 | 3300 | 2→2 | 1.78483582e-13 | 8.88178420e-16 |
| 0.1 | 10272 | 1→1 | 1.45135103e-13 | 1.33226763e-15 |
| 0.1 | 10272 | 1→2 | 1.61193961e-13 | 5.46228644e-14 |
| 0.1 | 10272 | 2→1 | 1.59277967e-13 | 5.46228644e-14 |
| 0.1 | 10272 | 2→2 | 1.64218558e-13 | 4.45823933e-16 |

**Aggregate native verdict:
`PASS_COMPLETE_NATIVE_CHECKPOINT_PARITY_CAMPAIGN`; 8/8 PASS,
maximum trajectory error `2.321188465396335e-13`, maximum
`H+D−W` parity discrepancy `5.4622864391340453e-14`.**

This is SCRIMP's **pause/restart fidelity** gate, not a proof of
exact discrete energy conservation or a mathematical statement about
time-discretization consistency. The same nonzero numerical-model
energy-balance defect occurs in continuous and restarted trajectories
and is deliberately recorded, **not an independent gate** for
checkpoint equivalence. Original physical weak forms, numerical
time scheme, and production `DPHS` checkpoint implementation were
not altered for this campaign.

**Scope qualification:** these native results establish full-path
cold restart parity at h=.2/.1, MPI 1 and 2 and gamma=.25 in the
`0→.024` pre-feedback time window, for fixed TSCN dt=.002.
They do **not** claim MPI4, time-step-adaptive history replay,
long-time stationary rotating flows, t>75 normal-feedback-active
restarts, or model energy accuracy. Existing #29 production public
checkpoint evidence remains separately valid.

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
