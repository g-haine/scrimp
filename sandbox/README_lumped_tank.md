# Experimental coupled PDE–ODE sloshing tank

This branch adds the first **spatially global scalar ports** to SCRIMP,
with an experimental full-domain 1-D viscous shallow-water / translating-tank
example at [`sandbox/sloshing_tank.py`](sloshing_tank.py).

**Status: prototype.** No successful full GetFEM/PETSc run has yet been
performed in the development environment. Do not interpret the notebook's
previous DG(0) "total momentum" as a genuine lumped state. Please validate
the discretization, boundary conditions, convergence, and energy budget
before relying on this example for research results.

## How to use it

In an installed SCRIMP environment with GetFEM, PETSc, Gmsh, etc.:

```bash
python -m unittest tests.test_lumped_ports
python sandbox/sloshing_tank.py
```

The example prints the defects for

```text
Mt*v + integral(mom) dx - P = 0
P(t)-P(0) - integral_0^t force(s) ds = 0
```

and, when run as a script, plots the total Hamiltonian (including the
tank kinetic energy). For shorter simulations without plotting:

```python
from sandbox.sloshing_tank import sloshing_tank
sim = sloshing_tank(mesh_size=0.1, final_time=0.05, plot=False)
```

## New API

Use a dynamic global scalar state and effort as follows:

```python
P = S.State("P", "Total momentum", "lumped")
sim.add_state(P)
sim.add_costate(S.CoState("e_P", "Tank velocity", P))
sim.add_lumped_port("P")
sim.set_initial_value("P", 0.)
```

For one scalar algebraic unknown:

```python
sim.add_port(S.Port("velocity", "v", "unused", "lumped",
                    algebraic=True, substituted=True))
sim.add_lumped_port("velocity")
```

A lumped port has **one** GetFEM degree of freedom regardless of the
number of mesh elements. No FEM should be added to a lumped port. Such
ports have `region=None`; they can be used in GWFL expressions over a
distributed domain. Since a test function for a global scalar is constant,
an integrated ODE residual must be divided by the **domain measure**:

```text
(P - forcing)*Test_P/ell
```

when the domain is a one-dimensional interval of length `ell`. In
contrast, the algebraic constraint

```text
(mom + (Mt*v-P)/ell)*Test_v
```

enforces `integral(mom) dx + Mt*v - P = 0` in a **single global row**.

The backend also synchronizes the GetFEM state with each PETSc trial iterate
when global scalar ports are active. Without that update, nonlinear
PDE–ODE couplings are assembled using stale iterates. The original
nonlinear mass-matrix mechanism has **not** been altered in this branch.

## Why conservative fluid momentum?

The source PDF uses `alpha_h=b*h`, `alpha_2=rho*(u+v)`, and
`P=Mt*v+integral(alpha_h*alpha_2)`. Direct discretization of its second
equation gives `M(alpha_h)*d(alpha_2)/dt`, a state-dependent mass matrix.

The example uses instead

```text
a_h = alpha_h
mom = alpha_h * alpha_2
v = (P - integral(mom) dx)/Mt
```

with Hamiltonian

```text
H = integral(mom^2/(2*rho*a_h) + rho*g*a_h^2/(2*b)) dx
    + Mt*v^2/2
```

and efforts

```text
e_h = rho*g*a_h/b - mom^2/(2*rho*a_h^2)
e_m = mom/(rho*a_h) - v
e_P = v,  e_pos = 0
```

The associated noncanonical skew operator is

```text
J(h,m) = [[0,                          -d_x(a_h * .)],
          [-a_h*d_x, -d_x(mom * .) - mom*d_x]]
```

and `R` dissipates `integral(rho*nu*a_h*(d_x e_m)^2)`. The
finite-dimensional canonical pair is `(P,pos)`, and `P_dot=force(t)`.
The total-energy identity at the continuous level is

```text
dH/dt = force(t)*v - integral(rho*nu*a_h*(d_x e_m)^2) dx
```

for closed, impermeable walls with zero power flux.

## Important limitations and follow-up checks

1. The current finite-element weak forms omit boundary power terms. This
   models homogeneous zero-boundary-power / impermeable walls formally,
   but does **not** independently impose `e_m=0` pointwise. A follow-up
   should formulate explicit boundary ports/constraints and check the
   boundary traces. No-slip and viscosity-induced boundary stresses
   require particular care.
2. There is no positivity limiter for `a_h`; choose strictly positive
   initial fluid depth and monitor it. No wetting/drying.
3. Correct discrete skew symmetry and viscous dissipation must be
   checked numerically after running the coupled model, especially for
   different FE orders and mesh refinements.
4. Time integration with BDF2 is not energy conserving. In addition to
   the algebraic and momentum residuals, compare discrete energy changes
   with mechanical work minus viscous dissipation (and numerical losses).
5. The original file `swe_diss_tank.ipynb` treated `P` as a DG(0)
   field; this is *not* equivalent to one global DOF.
6. The nonstandard mass-brick issue for `a_h*d(alpha_2)/dt` remains
   separate. This sandbox example deliberately avoids that issue by
   using conservative momentum.
7. `DPHS.export_to_pv` is for FEM variables; do not call it on `P`,
   `pos`, or `v`. Their trajectories are available with
   `get_solution(...)`.

The aim is to keep the feature narrowly scoped until validated in a full
SCRIMP numerical environment. No change has been made to the repository's
default branch.
