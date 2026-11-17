# Viscous rotational SWE: physical audit and continuation proposal

Status: **review / RFC, not a validated solver feature** (2026-11-17).
Reviewed source: `sandbox/viscous_swe_2d_mpi.py` on
`feature/getfem55-mpi-container-scalability`. Scientific comparison:
Cardoso-Ribeiro, Haine, Lefèvre, Matignon,
*Rotational shallow water equations with viscous damping and boundary
control: structure-preserving spatial discretization*,
Mathematics of Control, Signals, and Systems 37 (2025) 361–394,
DOI 10.1007/s00498-024-00404-6. Equations quoted below were checked
against the accessible Research Square v2 preprint (2024),
equations (32)–(43). Check changes in the version of record before
claiming exact numerical reproduction.

## Overall verdict

The **bulk/constitutive/port formulation** closely matches the
published weak PFEM model. No demonstrated error in the gyroscopic
term, signs of dissipation, or the physical Hamiltonian. The model
is **not yet validated as a numerical reproduction** of the rotating
tank case: one substantive FE-order difference, the two coupled
boundary constraints, differing IC/forcing regimes and integrators,
and insufficient independent invariants require tests.

## Audit of the actual Python bricks

- `h` is water height, `p = rho*v` is momentum per unit volume,
  with `e_h=rho*g*h+|p|^2/(2*rho)`, `e_p=p/rho`.
  `M_h` and `M_p[h]` yield `∂_t h` and `h∂_t p` respectively,
  not `∂_t(hp)` (paper Eq. 32, 35, 38).
- `Curl2D(p)=div(Rot*p)=∂_x p_y-∂_y p_x` and
  `Gyro(p)=Curl2D(p)*Rot`. The `G[h,p]` weak term
  `h*(Gyro(p)*e_p).Test_p` is pointwise skew in e_p
  and agrees with Eq. 32/38.
- The *symmetric* gradient `D(v)=Sym(Grad(v))`
  and the divergence give the two distinct viscous terms
  `-2mu*h*D(e_p):D(Test_p)` and
  `-2mu*h*div(e_p)*div(Test_p)`. These produce nonnegative
  volume dissipation for h>0 (paper Eq. 18, 35, 38).
- PFEM `D[h]`, `-D[h]^T` have matching h-weighted
  integration-by-parts structures, including the
  `-h*U_n*Test_h` boundary flux. `DPHS.add_brick`
  negates flow forms, a required sign convention;
  comparing raw string expressions without this convention
  leads to false sign errors.
- Boundary normal-pressure port `Y_n` enforces
  `Y_n=-e_h`; viscous port `Y` is the Lagrange multiplier
  associated with `h*(e_p-U).Test_Y`. The observable
  traction should be compared with `2mu*(D(e_p)n+
  div(e_p)n)` in a compatible weak sense; during
  inconsistent startup the multiplier may represent a large
  reaction and cannot simply be called physical stress.
- The paper requires `u^0 = u.n`; in this code the
  normal port `U_n` and vector port `U` are coupled by
  `U=U_n*Normal+omega*f(t)*Rot*Normal`. This is
  compatible for the *default* facet_tangent, open case.
  The `rigid_trace` diagnostic makes `U.n != U_n`
  on the faceted circle: useful as an isolated vector-law
  perturbation, but **NOT a fully consistent matched-boundary
  A/B test** (paper Remark 4 / Remark 5).
- For the open case `gamma=0` sets `U_n=0`.
  For the closed case the added feedback
  `U_n=-gamma*(1-f)*(Y_n+Y.Normal+rho*g*H_TARGET)`
  is structurally consistent with boundary passivity
  and shifted Hamiltonian
  `H_tilde = ∫ [h|p|²/(2rho) + rho*g*(h-h*)²/2]`,
  but this controller and its Lyapunov claim still require
  independent closed-loop numerical validation.
- The decreasing `ramp_expression` satisfies `f(0)=1`,
  `f(t>=200)=0`. The optional `start_from_rest_duration`
  adds a smooth rising envelope and changes the actual
  initial value problem; at t=0.01 with duration 0.02
  the wall has reached **half** the target angular speed.
- The code models a **flat-bottom**, constant density,
  constant viscosity, *positive-depth* SWE. It has no
  bathymetry, bed friction, or wet/dry handling. Testing
  `min(h)>0` is essential outside the current mild regime.

## Numerical discrepancies and outstanding gates

1. **Finite-element orders**: the published simulations use
   `P_{k+1}` for h/e_h and `P_k` for p/e_p/boundary
   variables, k=2, i.e. **P3/P2**. The sandbox historically
   used P2/P2 (all order 2). The opt-in
   `--height-order 3` now selects the published pair;
   the default stays 2 for reproducibility. This difference
   could affect inf-sup/trace compatibility and conditioning,
   but has not yet been shown to cause the initial impulse.
2. **Initial rotating state vs equilibrium**: uniform h=5
   and rigid rotating p= rho*omega*(y,-x) is generally
   **not** an exact steady SWE equilibrium, even on an exact
   disk. In an ideal unit disk with inviscid rigid rotation,
   balance of centripetal acceleration gives
   `h(r)=h_c+omega²*r²/(2*g)`; fixing initial mean h=5
   gives `h(r)=5+(omega²/(2*g))*(r²-1/2)`.
   For omega=0.02, g=0.1, this is only ±0.001 m,
   but it matters when asking for a *stationary* Newton
   initialization. A nonstationary uniform h at t=0
   remains legitimate if it satisfies *algebraic* constraints.
3. **Startup impulse**: at the original 0.02 rigid rotation,
   BDF2 with `init_dt=1e-9` gave
   `||Y||=5.99804e6`; at 1e-8 gave 599804 and
   at 1e-7 gave 59980.4. The multiplier scales
   with 1/init_dt. The zero-rotation equilibrium passes,
   with Y≈0. Matching only the vector velocity prescription
   (`rigid_trace`) and warm-starting e_p / U did not fix it.
   One pseudo step also leaves Y≈603189 and changes the
   rotating p state by 102.57% relative L2. Neither
   pseudo time nor the current BDF2 tiny-step procedure
   is demonstrated to be a DAE-consistency projection.
4. **Independent physical invariants**: the short transient
   spin-up from h=5,p=0 with default P2/P2 reached
   t=0.01 in 5 steps, no rejections or SNES failures,
   with maximum relative Hamiltonian balance defect
   6.0378585e-13 **normalized by H(0)**.
   This is encouraging but may be misleading if the
   actual energy exchange is orders of magnitude smaller
   than the large background potential energy.
   `diagnostics()` now reports an extra balance defect
   scaled by the observed changes/exchanges, plus
   an **independent volume conservation** audit:
   `V(t)-V(0)+∫_0^t∫_Γ h U_n =0`.
   These new diagnostics have not yet been executed.
5. **Temporal method**: the preprint reported BDF6 and
   dt<=0.01, whereas the script defaults to ARKIMEX-A2,
   with separate BDF2 startup. The PFEM power balance
   is semi-discrete; `get_quantity(CN=True)` samples the
   average of consecutive **states**, and
   `integrate_midpoint` is not an ARKIMEX-consistent
   exact quadrature. Test convergence under dt-halving;
   do not call the tiny defect exact time-preservation.
6. **No-flow, pressure/traction, and controller consistency**:
   check weak residuals separately at the final state.
   In the closed loop, test the shifted storage and
   integrated closed-loop supply, not only H.
   Validate on 1 and 4 ranks separately before trusting
   MPI scaling/performance measurements. The model
   explicitly declares `TS.EquationType.DAE_IMPLICIT_INDEX2`;
   that declaration is not a proof of the actual DAE index
   nor a substitute for constraint-residual diagnostics.
   `IJacobian` includes an opt-in/default finite-difference
   correction for the state-dependent mass-matrix Hessian:
   verify the fully assembled `IJacobian` against
   independent residual finite differences for both
   `zd=0` (stationary continuation) and nonzero `zd`.
   The startup logs also show a `ComplexWarning`
   when a PETSc vector is converted to real in
   `DPHS._sync_lumped_trial_state`; identify and log
   the imaginary component norm before assuming the
   conversion is harmless. None of these warnings alone
   demonstrates an error in the physical PDE.
7. **Parameters and target**: `RHO=1`, `GRAVITY=0.1`,
   `VISCOSITY=1e-3`, `H_INITIAL=5`,
   `H_TARGET=4`, `ROTATION=0.2` (script default),
   `gamma=0.25` (closed) are a specific model family,
   not an independently verified exact reproduction of
   the paper's particular fourth experiment. The
   experiments in this log use `--rotation=0.02`.
   Confirm exact values/forcing and initial field from the
   published source before claiming quantitative reproduction.

## Two different meanings of “wind-up”

### A. Transient physical spin-up (already implemented in sandbox)

Start with a physically compatible rest state h=5, p=0,
wall velocity 0; drive boundary rotation smoothly using
`a(t)=3(t/T)^2-2(t/T)^3` on [0,T], then a=1
(and multiply the long-time decreasing envelope).
Integrate the **real SWE** in real time. No artificial
equilibrium constraint. This is the appropriate mechanism
when one wants to model acceleration of the fluid.
A generic SCRIMP facility could generate smooth schedules
for boundary control data, e.g. `ControlSchedule` or
`Ramp`, but the actual dynamics stay under normal
PETSc TS. If physical t0 must have rotation already
established, complete the spin-up *before* the physical
t0 and restart at that endpoint. This gives a different
prepared state, not the arbitrary rotating, flat-h state.

### B. Stationary parameter continuation / homotopy (FSI-style)

Solve a *family of nonlinear stationary algebraic systems*
`F(z, zd=0; lambda)=0`, from lambda=0 to lambda=1,
with `omega(lambda)=lambda*omega_target` (or varying
viscosity or feedback gain). At each lambda:
- retain the last **accepted** complete state
  (including h,p,e_h,e_p,U_n,Y_n,U,Y) as the Newton guess;
- change the parameter as actual GetFEM data; reassemble
  the full residual/Jacobian at the new lambda;
- solve with PETSc SNES + line search, check converged
  reason and residual *and* `h_min>0`, port constraints,
  mass and maximum reasonable state norm;
- if the solve fails, roll back both lambda and state,
  halve the lambda step and retry; optionally increase
  the step on quick successes;
- use safe PC reset on changing sparsity with MPI-collective
  decisions; never silently substitute a stale symbolic LU;
- save checkpoints and per-stage residuals, Newton
  iterations, h minima, boundary mismatch, Hamiltonian
  and (for stationary states) dissipation versus supply.

**Crucial nullspace issue**: with `U_n=0` the stationary
system generally does not fix mean water height. The
continuation must enforce a prescribed total volume
`∫h=V0` (or an equivalent independent gauge), and remove
the associated redundant equation / introduce an appropriate
mass constraint. Otherwise stationary Newton may face a
singular Jacobian even at lambda=0.

This creates a **steady rotating equilibrium** at the
target parameter; it does NOT preserve a uniform initial
height 5 if the target rotation demands a radial pressure
profile. It is appropriate when the goal is a consistent
steady reference state, not reproduction of an arbitrary
nonstationary initial condition.

### C. DAE-consistent projection (a distinct feature)

Given target differential state (h0,p0), seek a nearby
`(h,p,e_h,e_p,U_n,Y_n,U,Y,hdot,pdot)` satisfying:
the algebraic constitutive equations, both compatible
boundary constraints and the full residual
`F(z,zd,t0)=0`. Optionally minimize a mass-weighted
distance of the *differential* components from (h0,p0),
with exact volume/no-flux constraints. This is a constrained
nonlinear solve (a KKT/saddle-point problem), not
`TSPSEUDO` and not ordinary stationary continuation.
Some IC requests may be infeasible and must fail
clearly rather than being silently projected far away.

## SCRIMP API proposal (not implemented)

Prefer two *separate* interfaces, avoiding the ambiguous
word “wind-up”:

- `S.Continuation(problem, parameter, start, stop,
  initial_step, min_step, max_step, validators, checkpoint)`:
  model-agnostic stage controller for nonlinear stationarity,
  accepts a callback for rebuilding/update of `F`, a
  callback for `SNES.solve`, stage acceptance, and MPI
  rollback. Start with natural-parameter continuation and
  PETSc SNES; add pseudo-arclength later if necessary.
- `S.Ramp(start, stop, t_start, t_end, shape="smoothstep")`:
  boundary/source-data generator for **physical time**
  as already prototyped by the sandbox spin-up.
- Subsequently `S.ConsistentInitializer(...)`: explicit
  constrained DAE initialization, with error norms and
  guaranteed mass/constraint checks.

Do not couple stationary continuation to `DPHS.init_step()`
or to the physical TS history. A converged stationary
state should be transferred in full to the initial
GetFEM/PETSc state at the *physical* start time, with
consistent algebraic port variables and verified residual.
Do not interpolate *energies* or raw DOFs without a
consistent state transfer when changing meshes.

**Existing numerical tools to consider**: GetFEM already
documents Moore-Penrose predictor-corrector continuation
in its C++ model API
(https://www.getfem.org/userdoc/model_continuation.html);
verify compatibility with Python bindings and SCRIMP's MPI
GetFEM/PETSc architecture rather than assuming it is
directly callable. PETSc SNES supports successive nonlinear
solves with warm-started vectors and therefore provides
a straightforward lower-level implementation route.

## Validated short spin-up: physical consistency gate (2026-11-17)

Both runs used one MPI rank, h=0.2 mesh, `rotation=0.02`,
`start_from_rest_duration=0.02`, BDF2 1e-8 s initial
consistency, ARKIMEX-A2 0.002 s main step, t_f=0.02,
and the conservative MUMPS PC reset. Both completed
successfully. The spin-up ramp reaches full angular
boundary speed at the requested endpoint.

| Metric | P2/P2 | P3/P2 (paper's mixed-order pair) |
| --- | ---: | ---: |
| Final time | 0.02 | 0.02 |
| Initial volume | 15.613227659735374 | 15.613227659735374 |
| Final volume | 15.613227659735394 | 15.613227659735387 |
| Maximum relative volume defect | 4.892219966223144e-15 | 3.754494392700675e-15 |
| Final Hamiltonian change | 0.00010770237425150597 | 0.00010777208511703051 |
| Integrated boundary work | 0.00010993744330331183 | 0.00011000439108787657 |
| Integrated viscous dissipation | 2.235019784423486e-6 | 2.2321542779310383e-6 |
| Max relative balance / H(0) | 1.2621959646808816e-11 | 3.8862666536424995e-11 |
| Max balance / measured energy exchange | 4.4814015033511895e-7 | 1.3789714530869638e-6 |

These figures establish a usable **nontrivial, physically
conservative smoke test** on 1 rank; they do not prove
time-step convergence, mesh convergence, or stability of the
original initially rotating profile. The relative defect
scaled by actual exchanges is the more meaningful metric
for this near-equilibrium flow.

## Explicit project priority: parallel SCRIMP performance

**Scope decision:** The immediate project goal is to optimize
parallel SCRIMP/GetFEM/PETSc/MUMPS behavior for PR #24,
**not** to complete the physical 2D Saint-Venant simulation
for LHMNC. The three proposed features
(A) stationary parameter continuation, (B) generic
time-dependent control ramp/spin-up, and (C) nonlinear
DAE-consistent initialization are **accepted as future ideas
but DEFERRED**. Do not implement these generic APIs or
pursue the unstable initially rotating SWE case as a
precondition to MPI optimization.

Use the now-working P2/P2 and especially P3/P2 short
spin-up as **correctness sentinels** for parallel work.
Suggested prioritized next gates:

1. Repeat the same deterministic full spin-up with 1, 2
   and 4 MPI ranks; compare final time, accepted/rejected
   steps, no SNES/KSP failure, final Hamiltonian and
   physical volume/work/dissipation against the one-rank
   baseline. Track differences with tolerances appropriate
   to floating-point MPI reductions; never infer accuracy
   from successful exit alone.
2. Profile `DPHS.solve()` with PETSc `-log_view` and
   SCRIMP's phase profiler, separating GetFEM assembly,
   GetFEM-to-PETSc bridge, Jacobian assembly, factor
   symbolic/numerical stages, linear solves and cleanup.
   Measure enough repeated samples and fixed CPU affinity
   to avoid attributing Python / logging / plotting overhead
   to parallel efficiency.
3. After correctness and phase profiling, increase the
   spatial problem size (mesh 0.1, then 0.05 where
   practical) for strong scaling; a 3300-DOF
   small-mesh case is NOT a credible speedup benchmark.
   Keep the 0.2 mesh as a fast correctness gate.
4. Compare conservative PC resets to an explicitly
   validated pattern-change-only policy in a targeted
   diagnostic. Existing saved 3300x3300 Jacobian replays
   on 1 and 4 ranks passed four persistent solves with
   only one reset on initial pattern growth (versus
   three unconditional resets), but this does NOT yet
   demonstrate a speedup or automatic safe integration
   in production `DPHS.IJacobian()`.
5. Leave `scrimp/dphs.py` default MUMPS reset semantics
   and production API unchanged until full multi-rank
   correctness, factor event counts, teardown and
   performance gates have passed.

**Avoid conflating a working transient physical spin-up
with a solution of the original rotating-initial-state
DAE. Neither the previously proposed continuation API nor
new LHMNC-specific physics are part of the current PR.**
