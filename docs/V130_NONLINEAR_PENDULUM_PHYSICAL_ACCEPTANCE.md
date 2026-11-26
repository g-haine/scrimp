# SCRIMP v1.3.0 — nonlinear finite-dimensional pendulum acceptance

## Status — SHORT NATIVE MPI1/MPI2 PASS (26 November 2026)

**NATIVE PASS** for the original nonlinear sine pendulum in both
unforced and sinusoidally forced cases, with the single exact
constitutive initial state, MUMPS MPI1/MPI2, eight accepted
fixed BDF2 steps (`dt=.01`, `t_f=.08`), and **all four**
GetFEM lumped state/co-energy coefficients on every rank.

**User-host archive:** `v130-nonlinear-pendulum-physical-pilot-2ae4bc7d-PVKCy6u4.zip`,
source commit `2ae4bc7d`, Docker `scrimp-mpi:v1.3.0-rc3`.
SHA256 **`dcf6fca2ec56ff3f24fe39a8e55bfb790b29000a5fead5ec9cfba7b23bbd9d80`**.
Independently opened: **22 files**, ZIP CRC passed; native NPZ
trajectories and physical identities recomputed from raw coefficients.
All six MPI rank reports and all scientific gates passed.

| Native short-case result | Unforced | Forced (amplitude 0.1) |
| --- | ---: | ---: |
| Actual differential initial error | 0 | 0 |
| Actual nonlinear co-energy initial error | 0 | 0 |
| Max sine constitutive error | 3.193e-11 | 3.650e-11 |
| Max BDF2 differential residual | 5.440e-14 | 5.707e-14 |
| Max BDF2 co-energy power residual | 4.915e-13 | 5.257e-13 |
| Max full-sine RK4 theta/p absolute error | 2.232e-4 | 2.228e-4 |
| MPI1 vs MPI2 full DAE error (rank0/rank1) | 0 | 0 |
| Original GetFEM H(0) | 6.255270428583832 | 6.255270428583832 |
| Original GetFEM H(.08) | 6.252059671693924 | 6.251772872110134 |

The GetFEM Hamiltonian matches the analytic **nonquadratic**
Hamiltonian exactly at every saved state. In the unforced
damped pilot H is monotonically decreasing; with forcing,
the actual nonzero external power was recorded.

**Important scientific scope:** the discrete
co-energy-weighted BDF2 identity is not an exact
telescoping Hamiltonian difference. At the short
unforced horizon, the observed `H(.08)-H(0)`
is about `-0.00321076`, whereas trapezoidal integration
of the continuous `-d*e_p²` power at the saved times
is about `-0.00113893`. Such quadrature must not be
presented as an exact discrete energy ledger: the
time-discrete solution and PETSc BDF2 startup contribute
numerical error. This discrepancy does **not** contradict
the verified discrete equation, and does not justify a
claim of exact energy preservation or physical dissipation.
Temporal refinement remains a separate qualification.

**Decision:** SHORT NATIVE SCIENTIFIC PILOT **PASS**.
Complete `t_f=12` simulation, MPI parity over all 1201
states, and directly runnable original-example PNG/SVG/GIF
visualizations remain **PENDING**.

The canonical `examples/nonlinear_pendulum.py` models a genuine
**nonlinear** single pendulum with two dynamical scalar states
`theta,p`, two constitutive scalar co-energies
`e_theta,e_p`, and a full sine restoring torque.
Its reference interval is only the GetFEM weak-form
normalization support; there is no physical spatial FE mesh.

Original parameters and equations **unchanged**:

```text
m=1, length=1, I=m length²=1, g=9.81, damping=0.08
theta_dot = p/I
p_dot = -m*g*length*sin(theta) - damping*p/I
        + torque_amplitude*sin(2*pi*0.3*t)
H = p²/(2*I) + m*g*length*(1-cos(theta))
e_theta = m*g*length*sin(theta)
e_p = p/I
H_dot = -damping*e_p² + torque*e_p
```

The default prescribed physical IVP is
`theta(0)=1.2`, `p(0)=0`. In contrast to
an artificial pseudo-time relaxation, this
fully explicit constitutive model permits the exact,
reproducible GetFEM scalar co-energies

```text
e_theta(0)=9.81*sin(1.2), e_p(0)=0.
```

The canonical example now initializes both algebraic variables
explicitly and uses `init_step=False`, leaving its **physical
BDF2**, nonlinear full sine, all weak bricks, MUMPS and other
SCRIMP solver options unchanged. The example exposes **one**
exact constitutive initialization path, with no legacy mode.
This exact-initial-state change **passed** the native MPI1/MPI2
short pilot above; the full-horizon gate remains pending.

## Correct nonlinear BDF2 physics

Unlike the *quadratic* mass–spring Hamiltonian, the pendulum
Hamiltonian includes `1-cos(theta)`. The telescoping
quadratic BDF2 `H` formula from the linear oscillator
is **mathematically inapplicable**.

For original BDF2 steps `n>=2`, define

```text
D2 z[n] = (3*z[n]-4*z[n-1]+z[n-2])/(2*dt).
```

The correct, exactly checkable **nonquadratic co-energy
BDF2 power identity** follows by pairing the *actual
discrete differential equations* with the current
constitutive efforts:

```text
e_theta[n]*D2theta[n] + e_p[n]*D2p[n]
  + damping*e_p[n]² - torque(t[n])*e_p[n] = 0.
```

This is **not** an equality of a BDF2 finite
difference of `H(theta,p)` to power. The nonlinearity
precludes that inference. The **continuous** physical
Hamiltonian is independently compared to the original
GetFEM Hamiltonian from the saved accepted state.
For the unforced damped case, its overall trajectory
is required not to increase within solver tolerance.
For the driven case the nonzero forcing power is
retained with its correct sign and no monotonicity claim
is made.

The PETSc order-two **startup half-step is not
asserted to follow an unsupported BDF1 identity**.
No temporal refinement-order claim is made by the pilot.

## Native short release diagnostic — COMPLETED

`scripts/run-v130-nonlinear-pendulum-physical-pilot.sh`
runs four *original* GetFEM/PETSc/MUMPS simulations:

- No external torque, `t_f=.08, dt=.01`, MPI1 and MPI2.
- Genuine external torque `0.1*sin(2*pi*0.3*t)`,
  same eight physical steps, MPI1 and MPI2.

Each native rank saves all nine accepted full four-DOF
DAE states, independent RK4 reference to the **original
full sine ODE**, constitutive errors, actual nonlinear
Hamiltonian, viscous dissipation, forcing power, all
second-order BDF2 residuals and the exact
**co-energy-weighted nonlinear BDF2 step identity**.
The host checker enforces physical gates, complete
MPI1/MPI2 state and power parity, MPI2 replica agreement,
real source provenance, and an evidence ZIP with CRC.

The short native scientific qualification **PASSED** for
both torque cases, for all ranks, as reported above.
No full default `t_f=12` or original-example
PNG/SVG/GIF qualification is claimed until separate acceptance. No changes are
made to the eight immutable `examples/published`
sources.
