# Experimental coupled PDE–ODE sloshing tank

This branch adds the first **spatially global scalar ports** to SCRIMP,
with an experimental full-domain 1-D viscous shallow-water / translating-tank
example at [`sandbox/sloshing_tank.py`](sloshing_tank.py).

**Status: experimental prototype.** The numerical workstation has run the
coupled model, including a 15-second Crank–Nicolson trajectory used for a
visualization. This does **not** settle the open physical questions: boundary
impermeability, spatial-field convergence and power-budget refinement.
The old notebook's DG(0) "total momentum" was not a genuine lumped state.

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

## Solver output and energy diagnostics

The GetFEM build used here may print `Assembly time ...` on every
assembly, even with `gf.util_trace_level(0)`. In GetFEM's C++ source
this line is written to `std::cout` under `GETFEM_PARA_LEVEL > 0`
without consulting the Python trace level. To filter **only this line**
from a terminal session without modifying SCRIMP or GetFEM:

```bash
set -o pipefail
python -u sandbox/sloshing_tank.py | sed -u -E 's/Assembly time [0-9.eE+-]+//g; /^[[:space:]]*$/d'
```

Use `2>&1` before the pipe if stderr should be filtered together with
stdout. The QT Wayland plugin warning is unrelated.

The script's `plot=True` mode now plots kinetic energy terms separately,
and compares `H(t)-H(0)` against accumulated mechanical work minus viscous
dissipation. The full absolute Hamiltonian can still be plotted via
`sim.plot_Hamiltonian(with_powers=False)`, but its gravitational offset
obscures the small variations in this example.

## Fixed-step convergence study

The first successful simulation had an energy-balance defect around
3.53e-5 at t about 0.19 (the plot stopped before the requested t_f=0.2
because the solver used PETSc's stepover mode and the extra final sample
was later discarded). The two energy curves nearly coincide, but this
is **not yet** a mesh/time convergence demonstration.

This branch now sets `ts_exact_final_time="matchstep"` for the sandbox
and records the final state explicitly when the exact endpoint is
reached. The example also accepts `time_step`, `save_step`, and
`adaptive` arguments. For reproducible fixed-step Crank–Nicolson tests (default):

```bash
python sandbox/tank_convergence.py
```

That script evaluates two meshes (h=0.05, 0.025) and three fixed time
steps (dt=0.008, 0.004, 0.002), prints mass, momentum and energy
defects, estimates refinement ratios and writes
`outputs/tank_convergence.csv`. Since the external force is known
analytically, the momentum diagnostic now compares to its **exact**
impulse as well as to the trapezoidal integration.

```bash
python sandbox/tank_convergence.py --mesh 0.05 0.025 --dt 0.008 0.004 0.002
```

To reproduce the earlier BDF2 power-defect studies, append
`--scheme bdf` (the new default is `--scheme cn`). Runtime simulation
has succeeded on the numerical workstation; the spatial-field and boundary
validation remains open.

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

## Why the first run stopped at t=0.0244 s

The first recorded run reached only 0.0243965 s and returned a huge
energy-balance defect, despite excellent mass conservation (9e-15) and
global constraint satisfaction (2e-16). Its total-momentum defect
(0.0003737) was almost exactly twice the **positive** applied impulse
(0.0001866), strongly suggesting the forcing sign was inverted.

GetFEM's `model.rhs()` from `build_rhs` is **minus the assembled
nonlinear weak-form residual**. SCRIMP's legacy `IFunction` computes
`K_tangent*z + rhs`, which is correct for certain homogeneous quadratic
weak forms but **incorrect for constant forcing and rational constitutive
terms**. In this branch, systems with global scalar ports and linear
time-derivative mass matrices now use the actual nonlinear residual

```text
IFunction(t,z,zdot) = M_linear*zdot + K_linear*z - GetFEM_rhs_nonlinear
IJacobian           = shift*M_linear + K_linear + K_nonlinear_tangent
```

The force-sign regression test in `tests/test_lumped_ports.py` checks the
simple ODE `P_dot = 2`. **This correction remains to be verified in the
target GetFEM/PETSc runtime**; it has not yet been rerun or validated.

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
4. Crank–Nicolson (and BDF2) is not in general exactly energy preserving
   for the present nonlinear system. In addition to
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

## Power-balance visualization

Run `python sandbox/animate_sloshing_tank.py --scheme cn` to generate
an MP4 with the true finite-element water surface, a tank-motion trace,
and synchronized plots of `H(t)-H(0)`, the trapezoidal time integral of
`force(t)*v(t)-D_visc(t)`, and their **numerical difference**. The video
also writes a CSV with the time-resolved power and energy diagnostics.
See `sandbox/README_sloshing_movie.md` for arguments and details.
