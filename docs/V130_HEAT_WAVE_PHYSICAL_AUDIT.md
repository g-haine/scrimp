# SCRIMP v1.3.0 — Heat–Wave example: physical verification (OPEN)

**23 November 2026. The correct reference article was confirmed by the author;
its actual numerical experiments are NOT assumed to be the same as
`examples/heat_wave.py`. This document is an independent audit, not a claim
that the current example is faithful or physically qualified.**

## Definitive published reference (read from the actual PDF)

G. Haine, D. Matignon and F. Monteghetti,
*Long-Time Behavior of a Coupled Heat-Wave System Using a
Structure-Preserving Finite Element Method*,
Mathematical Reports 24(74), 1–2 (2022), pp. 187–215.
Original PDF: https://imar.ro/journals/Mathematical_Reports/Pdfs/2022/1-2/11.pdf

Section 2 eqs. (2.1)-(2.2): with unit material coefficients
and original interface normal orientations,
```text
d_t theta = Delta theta         in Omega_1 (heat)
d_tt w = Delta w                 in Omega_2 (wave)
theta = d_t w                    on Gamma_int
d_{n1} theta = -d_{n2} w         on Gamma_int
w = 0                            on Gamma_2 (external wave boundary)
```

Section 3 eqs. (3.6)-(3.15): the coenergy variables
are heat temperature `theta`, heat flux
`J_Q=-grad(theta)`, wave velocity `v=d_t(w)`, and
wave stress `sigma=grad(w)`. In the article's
**specific boundary-port convention**:
`J_Q.n1=-u1, y1=theta, u2=v, y2=sigma.n2`;
the gyrator interconnection (3.14) is
`u1=-y2, u2=y1`.
This yields
`H=0.5*int_heat theta^2+0.5*int_wave(v^2+sigma.sigma)`
and the closed-system identity (3.15):
```text
dH/dt = - int_Omega1 J_Q.J_Q = - int_Omega1 |grad(theta)|^2 <= 0
```
These are **physical** statements, separate from the
specific time integrator or LU factor solver.

The main example uses different boundary port **causality**
(temperature control on the heat side, traction on the wave
side). Consequently, the raw code
`U_T=Y_w`, `U_w=-Y_T` is not automatically inconsistent
with the article's `u1=-y2, u2=y1`.
They must be mapped by **physical trace definitions and
normal orientations**, not by matching names or changing signs.
The correct physical interface checks are
`theta-v=0` and `J_Q.n1-sigma.n2=0`
(equivalently `d_n1(theta)+d_n2(w)=0`) in appropriate
discrete trace dualities. An algebraic
`int(U_T*Y_T+U_w*Y_w)=0` may be tautologically true
for a gyrator and alone does NOT prove this interface
represents the physical PDE.

## Differences that must be analyzed, not erased

1. **Time discretization:** the 2022 paper's quantitative
   temporal energy-decay study uses **Crank–Nicolson**
   (Introduction, Section 5); current
   `examples/heat_wave.py` defaults to **adaptive BDF2**.
   For BDF2, `∫P dt` from saved nodal powers using
   trapezoids is an approximation, not the scheme's exact
   discrete power balance. A genuine CN physical check
   is scientifically desirable *after* inspecting the DAE.
2. **FEM spaces:** the paper's eqs. (4.3)-(4.5)
   use `DG_k` for heat temperature and wave velocity,
   `RT_k` for heat flux and wave stress, and
   `CG_(k+1)` on interface traces. Current
   `examples/heat_wave.py` chooses `DG_1` for
   `T`, **`CG_2` for `e_Q`**, **`CG_2`
   for `p`**, **`DG_1` for `q`**
   (interface ports `DG_1`).
   These are NOT the same mixed FEM spaces, so
   stability/inf-sup and discrete interface balance
   from the article cannot simply be assumed.
   These choices are hypotheses to assess, not
   automatically proof of an invalid physical model.
3. **Initialization:** the example chooses Gaussian
   `T(0),p(0)` centered near the interface and
   `q(0)=0`; an index-2 pHDAE imposes additional
   interface constraints. The original initial
   `pseudo` consistency step may change the
   accepted initial FE state. Preserve both prescribed
   and accepted `z(0)` for comparison and report it;
   do not treat a visually plausible Gaussian as
   compatible without the constraint check.
4. **Physical outputs:** inspect separate heat/wave
   Hamiltonian contributions, nonnegative heat
   dissipation, interface *physical* trace residual,
   boundary work, and the coupled balance. Check at
   multiple meshes and time steps. Solver factorization
   parity is necessary for implementation robustness
   but not a certificate of correctness of the PDE.

## What the first real native pilot actually established

At source `f294f0edde8311299d82b6eeec889756d7094d1a`,
the original `heat_wave_mumps_mpi1` produced:

| Original observed quantity | Value |
|---|---:|
| Full GetFEM/PETSc degrees of freedom | 26198 |
| Physical final time reached | 0.004 |
| Accepted adaptive BDF2 steps | 39 |
| Rejected steps | 4 |
| SNES nonlinear solver failures | 0 |
| Hamiltonian postprocessing | Completed |
| ParaView `p` and `T` export calls | Completed |

The original **test harness was wrong**: it required
`last_ts_steps==4`, `last_ts_rejections==0`, and
exactly five uniformly spaced saved times. The default
PETSc BDF2 adapter uses an initial dt=.001 then adjusts
its timestep; the monitor saves **accepted** physical
states rather than interpolating to exact uniform times.
The failure did **not** establish MUMPS divergence,
a wrong physical PDE or an invalid energy law.

## Diagnostic method — current native test PENDING

The new, separate non-invasive
`sandbox/measure_v130_heat_wave_physics.py`
uses the actual original SCRIMP `get_quantity()`
and accepted time/state history:

- `H_heat = 0.5*int_1(T*T)`;
- `H_wave = 0.5*int_2(p*p+q.q)`;
- `D_heat = int_1(e_Q.e_Q)`;
- signed `int_10 U_T*Y_T`, `int_10 U_w*Y_w`
  and outer `int_20 U_bnd*Y_bnd`;
- optional trace residuals `T-p`, `T+p`,
  `U_T-Y_w`, `U_w+Y_T`, each explicitly
  marked unavailable if GetFEM cannot evaluate the
  expression on that original FE trace.

It saves all original time histories **without changing
the equations** and records nonconservation and
dissipation without pretending a BDF2 trapezoid is an
exact energy identity. An unexpected failure in the
required original physical quantities must FAIL with
the full original log retained — it must not silently
replace missing physics with zeros.

The bounded backend test is therefore split into a
`heat_wave`-only first phase, then the inviscid
`shallow_water` model. Numerical MUMPS/SuperLU_DIST
parity can only be claimed on the original full FE
trajectories when the accepted time grids agree; if
adaptive grids differ, the test may compare endpoints,
and must explicitly disclose the reduced scope.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
SCRIMP_V130_2D_SCOPE=heat_wave \
bash scripts/run-v130-example-2d-lu.sh
```

The output ZIP includes the original full FE trajectories,
per-subsystem Hamiltonians, interface powers,
heat flux dissipation, optional trace residuals,
host logs and postprocessed figures/VTU files.
**No numeric or physical PASS claimed before
running these actual GetFEM native diagnostics.**
`examples/published/` is not modified, mounted,
or touched. The original paper's experiments are
not silently substituted for this example.

**Next scientific decision, after seeing the real traces:**
identify whether a demonstrated defect arises in
(1) port causality/signs, (2) mixed FE spaces,
(3) consistent initialization, (4) discretization/time
integration, or (5) postprocessing; change one element
at a time outside the immutable `published` subtree.
