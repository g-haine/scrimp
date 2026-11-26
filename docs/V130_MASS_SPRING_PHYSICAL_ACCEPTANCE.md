# SCRIMP v1.3.0 — mass–spring–damper physical qualification

## Status, 26 November 2026

**Original canonical linear lumped pH model: SHORT MPI1/2 NATIVE PASS.**
The full original `t_f=8` BDF2 trajectory and canonical example-native
PNG/SVG/GIF user visualizations are **STAGED, NATIVE PENDING**.

The canonical original weak GetFEM/PETSc model has exactly four global
scalar unknowns `q,p,e_q,e_p`, with the actual prescribed
parameters `m=1.5, k=5, c=0.25, f(t)=0.2 sin(2 pi t)`,
`q(0)=0.4, p(0)=0`. Their exact initial co-energies are
`e_q(0)=2, e_p(0)=0`. The reference interval is a GetFEM
residual normalization device, NOT a finite-element discretization
of a spatial oscillator.

## Previous unqualified pseudo-time initialization

The first native short pilot, image
`scrimp-mpi:v1.3.0-rc3`, source `2cd517d8`,
`t_f=.08, dt=.01`, MPI1/2 MUMPS, completed all eight BDF2
physical time steps and preserved every original equation.
The maximum differential BDF2 residual was `2.28e-14`,
the exact step-energy defect `4.13e-14`, and all four
GetFEM DAE coefficients matched across MPI ranks to `8.88e-16`.

**But its original IVP initial-data gate FAILed.** The explicitly
selected legacy `init_step_ts_type="pseudo"` integrated
for a small `dt_init=dt**2=0.0001` and published the resulting
differential state as the physical initial state:
`q(0)-0.4=-1.33331e-8`, `p(0)-0=-0.00019999666`.
The latter is a finite, avoidable change of the physical
initial condition, NOT a small error of the BDF2 time discretization.
The legacy model remains reproducible explicitly by calling
`mass_spring_damper(initialization="legacy_pseudo")`, but
is **not** a qualified default for the prescribed IVP.

## Accepted exact constitutive initial state: NATIVE PASS

**Original user-host archive:**
`v130-mass-spring-physical-pilot-2568751b-mGLuqbrA.zip`;
SHA256
`e85c6b7ef97a8590e7b266719e6dd70760f85a081013b6fbb443af65074a636c`.
Checked original **12 members, ZIP CRC PASS** and directly
compared all native MPI rank0/1 NPZ trajectories.
Source `2568751b`, same Docker image, MUMPS MPI1 and MPI2.

Only the example's handling of **exact initial co-energy DOFs**
changed: `e_q=k*q0`, `e_p=p0/m` assigned directly with
GetFEM's `set_variable`, and `init_step=False` for the
default case. **No change** to constitutive, forcing,
Dirac, mass or damping bricks; no change to the subsequent fixed
PETSc BDF2 scheme. Generic SCRIMP DAE initialization is untouched.

- Native exact initial `(q0,p0,e_q0,e_p0)=(.4,0,2,0)`
  on all ranks, with zero error on either prescribed
  dynamic initial value. Pseudo-TS is not executed.
- Eight real accepted steps and nine GetFEM full DAE
  states; zero rejection or nonlinear-solver failures.
- True co-energy laws `e_q=kq`, `e_p=p/m` agree to
  `2.22e-16` and `1.39e-17` respectively.
- Original GetFEM Hamiltonian equals analytical
  `H=.5*k*q²+p²/(2*m)` to the last recorded bit.
- Exact time-discrete BDF2 differential state residuals
  at `k>=2`: `1.32e-14` (q) and `3.90e-15` (p).
- Exact quadratic **BDF2 step energy-power identity**
  at `k>=2`: max residual `4.2107377e-14`.
  The PETSc BDF2 startup half-step is NOT covered by
  a fictitious BDF1 identity.
- Independent RK4 integration of the **same forced ODE**:
  maximum original q/p absolute discrepancy
  `5.9975275e-5` over `t∈[0,.08]`.
- Complete four-variable MPI1 vs MPI2 history:
  max absolute coefficient discrepancy
  `3.46944695e-18`; MPI2 rank0 and rank1 have
  byte-identical original full-DAE histories.
- Original energy `H(0)=.4`, `H(.08)=.3995943677503393`.

**Decision:** the repaired **short** canonical original
mass-spring DAE/BDF2 initial-value problem has passed all
scientific native gates, in MPI1 and MPI2. This is not a
long-run or visual qualification.

## Full native example and user visualization release gate — PENDING

The source-only staged runner
`scripts/run-v130-mass-spring-full-showcase.sh` performs:

1. A true eight-step MUMPS MPI1 preview with the **canonical
   example's own PNG/SVG/GIF renderer** (no research renderer).
2. The complete `t_f=8,dt=.01` original system in
   MPI1 (800 steps, 801 accepted snapshots), with all
   canonical user figures saved under
   `outputs/showcase/mass_spring_damper/figures/`.
3. The identical full history in MPI2, with source
   independent host verification of all four original
   GetFEM DAE DOFs, co-energies, Hamiltonian,
   BDF2 derivative residuals and quadratic BDF2
   power identity for `k>=2`.
4. Independent RK4 forcing benchmark on the same ODE,
   full MPI rank0/1 parity, figure-file and GIF frame
   verification, and a ZIP containing the actual
   original FE/lumped traces, scientific evidence and
   **figures rendered by the canonical example**.

The default user command `python examples/mass_spring_damper.py`
runs `mass_spring_damper(plot=True)`; the new plot code
renders three PNG/SVG pairs for trajectories, phase space
and the BDF2 physical energy-power identity, plus
`mass_spring_motion.gif` using selected **actual
saved GetFEM lumped states**. The GIF is explicitly
labeled *schematic*: it does not invent a spatial
FE displacement mesh. No fabricated time interpolation.

This gate must **not** be called PASS before actual
user-host Docker execution. No claim about time-step
refinement, other LU backends, MPI4+, or nonlinear
pendulum qualification is made.
