# SCRIMP v1.3.0 — nonlinear lumped pendulum physical acceptance

## Current status — SHORT PETSc CN MPI1/2 NATIVE PASS

The **short original CN model is qualified** in unforced and forced
configurations at `dt=.01`, `t_f=.08`, MUMPS MPI1 and MPI2.
This is a **scientific short pilot only**. Full 12-second, temporal
refinement and canonical user figures remain unqualified pending
the long native campaign.

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

**Only this short CN gate is closed.** The original canonical
12-second figures are not yet natively tested; source documentation
must not claim they are scientifically qualified.

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

**The last native user-host calculation was BDF2, NOT CN**:
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

**The default full `t_f=12` run, its complete MPI
history, independent temporal refinement and
PNG/SVG/GIF from the original example remain
unqualified** until the separate long native gate
has passed. The eight-step CN gate above is complete.

No claims are made about other LU backends or MPI4+.
