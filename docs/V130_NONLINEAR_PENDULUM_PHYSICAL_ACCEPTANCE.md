# SCRIMP v1.3.0 — nonlinear lumped pendulum physical acceptance

## Current status — FULL PETSc CN 12-second MPI1/2 NATIVE PASS

**FULL NATIVE SCIENTIFIC AND VISUAL PASS**, 26 November 2026.
The original canonical nonlinear full-sine pendulum has been
qualified with PETSc TSCN (endpoint trapezoidal), exact
constitutive initial state, MUMPS MPI1/MPI2, two physical
configurations (free and torque amplitude 0.1), complete
12-second trajectories, temporal refinement, and real
canonical PNG/SVG/GIF output.

**User-host acceptance archive:**
`v130-pendulum-CN-12s-showcase-b93cfd27-Kd1WLVcv.zip`.
Source commit `b93cfd27dc3d69104da8d615d47fe6d1cc196302`,
image `scrimp-mpi:v1.3.0-rc3`;
SHA256 **`210b8fef9886c4966a42866c1d89bb786ead5f08f071741ddf16d62e581c47f0`**.
Independent ZIP verification: **74 members, CRC PASS**.
All eight planned native runs returned exit code 0.
All eight immutable published source examples remained unchanged.

| Full 12-second CN result | Free | Forced, torque amplitude 0.1 |
| --- | ---: | ---: |
| True accepted steps / states, dt=.01 | 1200 / 1201 | 1200 / 1201 |
| True accepted steps / states, dt=.005 | 2400 / 2401 | 2400 / 2401 |
| Max nonlinear sine constitutive error, dt=.01 | 3.464e-14 | 4.263e-14 |
| Max native H versus analytic H | 8.882e-16 | 8.882e-16 |
| Max CN co-energy power defect, dt=.01 | 4.836e-12 | 6.915e-12 |
| Max MPI1/MPI2 complete four-DOF discrepancy | 7.860e-14 | 1.015e-13 |
| Max theta/p RK4 error, dt=.01 | 5.735105e-3 | 5.732085e-3 |
| Max theta/p RK4 error, dt=.005 | 1.433858e-3 | 1.433163e-3 |
| RK4 error ratio coarse/refined | 3.999772 | 3.999606 |
| **Observed order log2(ratio)** | **1.999918** | **1.999858** |
| Native H(0) | 6.255270428583832 | 6.255270428583832 |
| Native H(12) | 2.4964883366297386 | 2.4909090620990706 |

Every run completed all its planned CN steps with the
original four true GetFEM lumped DOFs and physical gates
satisfied. MPI2 rank0/rank1 trajectories agree exactly;
MPI1/MPI2 deviations are numerical roundoff. The unforced
case displays strictly decreasing native H at every saved
coarse interval (maximum increment `-2.5674991e-9`).

The independent checker recomputed the true full-sine
constitutive error, CN endpoint equations, nonquadratic
co-energy power, accumulated Hamiltonian chain-rule
accounting and independent full-sine RK4 errors directly
from the raw NPZ records. The nonlinear Hamiltonian is
**not** asserted to be exactly conserved by CN.

**The canonical example itself** created the accepted
free and forced user graphics: three PNG/SVG pairs
(trajectory, phase portrait, physical energy and power),
plus a 36-frame schematic pendulum GIF for each complete
12-second simulation. These use actual saved GetFEM
states without inventing a spatial FE mesh or
interpolating physical solution states. The figure
manifests and native NPZ histories were checked against
the source trajectories, all PNGs opened successfully
and the plotted physical balances were visually inspected.
All graphical preview and full-run image gates PASS.

**Release decision: FULL SCIENTIFIC GATE CLOSED for this
configuration.** No additional native run is required
for the originally agreed free/forced CN MPI1/2 and
temporal-order acceptance. No claim is made for different
nonlinear parameters, other LU backends, MPI4+, or
a discrete-gradient exact-energy-preserving scheme.

## Earlier short PETSc CN acceptance

The **short original CN model also passed** in unforced
and forced configurations at `dt=.01`, `t_f=.08`,
MUMPS MPI1 and MPI2, prior to the full acceptance.

User-provided native evidence:
`v130-nonlinear-pendulum-physical-pilot-6621459f-PnkjIHci.zip`.
Original source revision `6621459f`,
image `scrimp-mpi:v1.3.0-rc3`.
**SHA256 `7c076d2c782b3702bcf046d7c849eb22860fe2c7199b57f68562c0d3aec63ec9`**,
**22 files**, archive **CRC PASS**. The original six NPZ trajectories
were read and their endpoint-CN differential equations, power pairing,
energy accounting, and complete MPI1/2 DOFs independently recomputed.

| Short native CN metric | Free | Forced, amplitude 0.1 |
| --- | ---: | ---: |
| Native exact initial state error | 0 | 0 |
| Max full-sine constitutive error | 1.776e-15 | 1.776e-15 |
| Max endpoint CN derivative residual | 1.533e-14 | 1.277e-14 |
| Max original co-energy CN power defect | 1.434e-13 | 1.114e-13 |
| Nonzero nonlinear H chain-rule defect | 2.385e-5 | 2.379e-5 |
| Max cumulative CN chain accounting defect | 3.345e-15 | 2.593e-15 |
| Max error to independent full-sine RK4 | 2.269e-5 | 2.262e-5 |
| All four original DAE DOFs MPI1/2 parity | 1.110e-16 | 1.110e-16 |
| Original H(0) | 6.255270428583832 | 6.255270428583832 |
| Original H(.08) | 6.254143694389639 | 6.253854348550849 |

All four original cases completed without numerical failures:
unforced/forced at MPI1 and MPI2. Both rank0 and rank1
MPI2 state histories are byte-identical. The initial state is
exact, including both nonlinear algebraic co-energies, and
no artificial initialization TS is used.

The short CN pilot remains part of the documented evidence;
the full twelve-second and visual acceptance above supersedes
its limited scope.

The **only canonical integration method** for
`examples/nonlinear_pendulum.py` is now
`ts_type="cn"` (PETSc 3.26 endpoint trapezoidal rule,
not implicit midpoint), with `ts_adapt_type="none"`,
fixed `dt=.01`, and exactly consistent initial
lumped co-energies with `init_step=False`.

The physical GetFEM nonlinear model is unchanged:
`m=1`, `length=1`, `I=1`, `g=9.81`,
`damping=.08`, `theta0=1.2`, `p0=0`.
There are four true lumped GetFEM scalar DOFs:
`theta,p,e_theta,e_p` and no spatial PDE.

```text
H(theta,p) = p²/(2*I) + m*g*length*(1-cos(theta))
e_theta = m*g*length*sin(theta)
e_p = p/I
theta_dot = e_p
p_dot = -e_theta - damping*e_p
        + torque_amplitude*sin(2*pi*0.3*t)
e_theta0 = m*g*length*sin(theta0)
e_p0 = 0
```

No old BDF2/pseudo choice is retained in the canonical
example. The current GetFEM/PETSc CN system has **passed the
eight-step MPI1/2 native gate** described above.

## Correct trapezoidal CN step identities

Let `delta a = (a_n-a_(n-1))/dt` and
`mean(a) = (a_n+a_(n-1))/2`.
PETSc `TSCN` uses the **endpoint** trapezoidal
variant of `TSTHETA`, with theta=1/2 and endpoint=true.

The independent checker must measure the true nonlinear
endpoint GetFEM co-energy variables and then test
**all physical time steps n>=1**:

```text
delta(theta) = mean(e_p)
delta(p) = -mean(e_theta) - damping*mean(e_p) + mean(torque)
mean(e_theta)*delta(theta) + mean(e_p)*delta(p)
      + damping*mean(e_p)^2 - mean(torque)*mean(e_p) = 0
```

This **discrete co-energy power identity** follows from
the actual CN differential residual and interconnection
skew-symmetry. It is **not** an exact equality between
Hamiltonian change and external port work minus damping:
`H(theta,p)` is not quadratic.

The actual nonquadratic chain-rule defect is

```text
Gamma_n = delta(H) - mean(e_theta)*delta(theta)
                      - mean(e_p)*delta(p)
delta(H) = mean(torque)*mean(e_p)
           - damping*mean(e_p)^2 + Gamma_n
```

Therefore the accumulated Hamiltonian difference from
the CN discrete work must equal the accumulated
`dt*Gamma_n` up to the measured CN residual.
**Gamma_n is measured, never assumed to vanish.**
Ordinary trapezoidal quadrature of nodal physical
`torque*e_p - damping*e_p²` is also different
from the genuine CN discrete work; the canonical
figures clearly label these quantities.

## Historical fixed-BDF2 qualification — archive only

**An earlier native user-host calculation was BDF2, NOT CN**:
source `2ae4bc7d`, original archive
`v130-nonlinear-pendulum-physical-pilot-2ae4bc7d-PVKCy6u4.zip`,
SHA256
`dcf6fca2ec56ff3f24fe39a8e55bfb790b29000a5fead5ec9cfba7b23bbd9d80`.
Archive opened: 22 files and CRC PASS. Four runs
(unforced/forced with torque amplitude 0.1; MPI1/MPI2),
eight BDF2 steps each, all four lumped DOFs present.
Physical initial data exact, full-sine relations accepted,
full-MPI state discrepancy 0, independent RK4 check PASS.
The old *BDF2* step co-energy residual was below
`5.26e-13`. These results **cannot** be used to
qualify or infer the behavior of the new CN scheme.

## Native short CN gate — completed

Run directly from the user-host repository:

```bash
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-nonlinear-pendulum-physical-pilot.sh
```

The script must run the *canonical example* in
four original-physics configurations:
unforced `torque_amplitude=0` and forced
`torque_amplitude=.1`, in each case MPI1/MPI2,
`t_f=.08`, `dt=.01`. Each rank stores the
nine true accepted four-variable GetFEM snapshots.

Scientific acceptance requires exact prescribed initial
states, true nonlinear end-point constitutive laws,
correct GetFEM H, genuine PETSc TSCN endpoint
residuals and their power pairing, nonzero nonlinear
chain-defect accounting, agreement with independent
full-sine RK4 ODE, complete MPI1/MPI2 parity,
all steps accepted, and immutable `examples/published`
sources. The pilot creates a CRC-checked evidence ZIP.

## First full native CN attempt — REJECTED, then corrected

The original 12-second campaign at source `9d7696ff`
was executed on the user host using Docker
`scrimp-mpi:v1.3.0-rc3`.
Evidence archive:
`v130-pendulum-CN-12s-showcase-9d7696ff-BaMf42Vz.zip`,
SHA256
`79faf60b7a9d3250b4f302e33a30c513ed06881a34012b392364b60485f9cbc1`,
**35 members, CRC PASS**. The two 8-step original
program-generated free and forced previews passed, including
PNG/SVG/GIF.

The **12-second free MPI1** simulation completed 1200 accepted
CN steps, with **0 rejected steps** and **0 nonlinear solver
failures** (2393 SNES iterations). Immediately afterward,
the original example's **physical constitutive gate rejected**
its genuine full GetFEM history:

```text
max |e_theta - 9.81 sin(theta)| = 1.3611093940824048e-7
acceptance limit             = 1e-8
max |GetFEM H - analytic H| = 8.881784197001252e-16
```

Because rendering is intentionally fail-closed, the full
program produced **no unqualified full figures**.
The full-host script correctly marked free MPI1 unqualified
and skipped subsequent MPI2, refinement, and forced runs.
**The 1200-step TS completion does not imply success
of the full scientific acceptance.**

The failure was compatible with loose nonlinear SNES
convergence thresholds, but the old complete state
trajectory was not archived, so that root cause was
not established conclusively. The canonical example
subsequently set `snes_atol=1e-12`,
`snes_rtol=1e-12`, `snes_stol=1e-14`,
without changing the physical GetFEM bricks, sine
constitutive law, CN integrator or `1e-8`
physical gate.

**The new native full campaign passed** all scientific
and visual gates, with maximum full-horizon sine
constitutive errors of only `4.263e-14`.
Thus the tightened-Newton configuration is validated
for the tested model, although this does not prove
which earlier SNES stopping condition caused the
original residual. Diagnostic archival of true
failed GetFEM histories remains enabled if any
future run rejects a constitutive law.

## Full twelve-second CN physical and visual gate — COMPLETED

The original canonical example now contains its own
`render_pendulum_visuals`, called by
`nonlinear_pendulum(plot=True)`. It generates
**directly from genuine accepted GetFEM lumped states**:
`pendulum_trajectory.{png,svg}`,
`pendulum_phase_portrait.{png,svg}`,
`pendulum_energy_and_power.{png,svg}` and
`pendulum_motion.gif` (schematic motion, no spatial FE mesh
or interpolated state). Both the free and forced
variants have their own `outputs/showcase/nonlinear_pendulum/`
subdirectories.

The canonical plots **do not** claim that nonlinear
CN conserves `H` exactly. They explicitly compare the
true Hamiltonian change, CN discrete work/dissipation
and the measured nonquadratic chain-rule defect.
These user figures **passed the complete native 12-second
campaign**, including bitwise equality of the original
GetFEM state histories, PNG/SVG presence and 36 true
accepted-time GIF frames in both physical variants.

To run the full acceptance after the short native PASS:

```bash
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-nonlinear-pendulum-CN-full-showcase.sh
```

The **successfully executed full acceptance runner** checked:
1. 8-step graphical **preflight for free and forced** original
   canonical examples (9 real GetFEM states and 9 GIF frames each).
2. Native `t_f=12, dt=.01`, 1200 accepted PETSc CN
   steps (**1201 true full DAE states**), MUMPS MPI1
   and MPI2, both physical torque cases.
3. Native `t_f=12, dt=.005`, 2400 accepted CN
   steps (**2401 true full DAE states**), MUMPS MPI1,
   again free and forced.
4. The source-independent verifier checks every recorded
   CN differential step, full-sine GetFEM constitutive
   law, exact co-energy power identity, Hamiltonian
   chain-rule accounting, complete MPI1/MPI2 DAE parity
   and both canonical figure inventories/actual GIFs.
5. A convergence gate compares full numerical
   trajectories to independent high-resolution
   full-sine RK4. The observed `max(q,p)` error ratio
   `error(dt=.01)/error(dt=.005)` must be within
   `[3.3,4.7]`, consistent with second order, for
   **both** forced and unforced cases. No fabricated
   time states or visualizations are allowed.
6. All native rank reports, raw original full histories,
   canonical figures, file integrity evidence and
   summary are packaged in a CRC-verified ZIP.

**NATIVE FULL PASS**: each original MPI run, both refined
runs, the exact CN co-energy power/chain-rule evidence,
second-order RK4 refinement and user-generated images
passed. The uploaded SHA256-pinned archive and independent
recomputation above close the agreed release gate.

No claims are made about other LU backends or MPI4+.
