# 2D viscous Saint-Venant MPI experiment (work in progress)

This benchmark explores the circular-tank open-loop and closed-loop cases in
Cardoso-Ribeiro, Lefèvre, Haine, Matignon, *Passivity-Based Boundary Control
for the 2D Shallow Water Equations with Viscous Dissipation* (provided PDF).

It is **experimental, not validated**: no run should be considered a
reproduction of the figures until the spatial fields, boundary control,
time stepping and power balance have been checked on the actual GetFEM/PETSc
MPI environment.

## Code and implementation choices

Driver: `sandbox/viscous_swe_2d_mpi.py`.

- Disk of radius 1; CG2; density 1; gravity 0.1; viscosity 0.001;
  h(0) = 5; rotating tangential velocity approximately 0.2.
- Gamma = 0: normal wall velocity zero, tangential wall feedforward.
- Gamma = 0.25: implicit normal output feedback using the normal pressure
  port and viscous traction port, with the same tangential feedforward.
- PETSc `arkimex/a2` (TR-BDF2-type DIRK) with MUMPS for an MPI-capable LU
  solve and a maximum adaptive timestep of 0.01 s.
- The paper does not give the exact formula for f(t). The driver uses a C1
  cubic smoothstep from 1 to 0 between 75 and 200 s. The rotating field is
  regularized at the center of the disk, where a radial-constant tangential
  field is not continuous. These choices mean *quantitative reproduction is
  not established* even if the calculations converge.

## Smoke tests (start here)

Run from the **repository root**, in the SCRIMP conda environment with
GetFEM, PETSc, petsc4py and parallel MUMPS available. Each run should
finish without an MPI deadlock and save diagnostics for every case.

```bash
git switch release/1.3.0-preparation

# First validate in serial, then repeat using MPI
python sandbox/viscous_swe_2d_mpi.py --case open \
    --mesh 0.2 --final-time 0.5 --save-dt 0.1 --max-step 0.01

mpiexec -n 2 python sandbox/viscous_swe_2d_mpi.py --case open \
    --mesh 0.2 --final-time 0.5 --save-dt 0.1 --max-step 0.01

# Accelerate the ramp to exercise the feedback during a short smoke run
mpiexec -n 2 python sandbox/viscous_swe_2d_mpi.py --case closed \
    --mesh 0.2 --final-time 0.5 --save-dt 0.1 --max-step 0.01 \
    --ramp-start 0.1 --ramp-end 0.4
```

If the MPI launcher is named `mpirun`, use that instead of `mpiexec`.
To test 4 ranks, replace `-n 2` with `-n 4`. These are real MPI
computations, *not* separate independent processes running the same script.

## Full-length experiments

Once smoke tests pass and the energy defect is sensible, run the intended
250 s cases (these may be expensive at full resolution):

```bash
mpiexec -n 4 python sandbox/viscous_swe_2d_mpi.py --case open \
    --mesh 0.05 --final-time 250 --max-step 0.01 --save-dt 1.0

mpiexec -n 4 python sandbox/viscous_swe_2d_mpi.py --case closed \
    --mesh 0.05 --final-time 250 --max-step 0.01 --save-dt 1.0
```

The diagnostics are written to
`sandbox/outputs/viscous_swe_2d/{open,closed}/`. Each run writes
`energies.npz`, `metrics.json`, `power_balance.png` and
`lyapunov.png`. The SCRIMP built-in Hamiltonian plot is written separately
to `outputs/png/viscous_swe_{open,closed}_Hamiltonian.png` relative to
the working directory.

## Physics and acceptance criteria

In the notation of the paper, equation (13) states

```text
dH/dt = P_normal + P_viscous - D_grad - D_div
```

where D_grad and D_div are nonnegative. The driver computes the
time-integrated balance defect

```text
B(t) = H(t) - H(0)
       + integral_0^t (D_grad + D_div) dt
       - integral_0^t (P_normal + P_viscous) dt.
```

It uses midpoint *states* via `get_quantity(CN=True)` and the
rectangle rule on the resulting interval values: **do not apply a second
trapezoidal rule** to these midpoint-state quantities. Inspect
`metrics.json` and compare serial versus 2/4 ranks.

The paper reports maximum relative energy-balance defects around 0.3% for
the open-loop case and 1.2% for closed-loop, but those are **reference
results**, not assertions of identical errors with the approximations
described above.

Before v1.2.1, require (at minimum):

1. No deadlock when computing/plotting port powers in 1, 2 and 4 MPI ranks;
   `tests/test_parallel_powers.py` also checks rank-call symmetry.
2. A consistent, acceptably small energy-balance defect with no non-finite
   field values or negative depth.
3. Serial/MPI agreement (within solver tolerances) for Hamiltonian, power
   fluxes and final profiles.
4. Closed-loop stabilization near the target height (4), monotone
   Lyapunov decrease after the feedback starts, and qualitative agreement
   with the paper.
5. No tag/release until the real MPI tests have succeeded.

The generic `DPHS.plot_Hamiltonian(with_powers=True)` checks SCRIMP's
port-power postprocessing and exercises issue #6; it does **not**
automatically incorporate the height-weighted boundary fluxes for this
specific formulation. Use `power_balance.png` for the physical balance.

## v1.3.0 issue #6: bounded native qualification runner (source gate, NOT yet native-PASSED)

The maintained release-branch runner now separates a zero-rotation
**rest sanity check** from a **nontrivial physical smoke**. Both use
the existing experimental viscous 2D Saint-Venant weak forms, fixed
BDF2, mesh 0.2, ten nominal 0.002 s steps through 0.02 s,
and the previously qualified immutable rc3 image; no dependencies
are rebuilt. This is NOT a reproduction of the paper's long-time
250-second experiment or proof of pH discrete energy preservation.

```bash
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
    bash scripts/run-v130-issue6-viscous-mpi.sh rest

# Only after rest passes; NONZERO rotation, both open and closed loop.
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
    bash scripts/run-v130-issue6-viscous-mpi.sh physical
```

The run is intentionally fail-fast at the **first** failed native
case (including PETSc TS failure, timeout, nonpositive sampled height
coefficients, parity or power defect); retained host-owned logs remain
in `sandbox/outputs/v130-issue6-viscous-*/`. Each (case, MPI size)
uses a distinct working directory and artifact directory, preventing
cross-run mesh and output clobbering. Only rank 0 writes figures and
`metrics.json`; every MPI rank writes its **own** generic-port-power
record to `mpi_port_powers_rankN.json`. The comparator checks these
traces **within** the MPI run as well as **against serial**, compares
full endpoint state vectors and sampled physical h-weighted powers,
and verifies the signed energy and volume balances.

A check of `minimum_height_coefficient>0` is necessary but cannot
establish pointwise positivity for quadratic finite elements.
The short-smoke provisional gates are at most 5% error relative
to actual energy exchange, 1e-4 relative mass-balance error, and
rank agreement according to the comparator; thresholds are not
paper-reproduction claims. If the physical configuration does not
converge (historical initialization/geometry issues), keep #6 OPEN
and diagnose rather than switching silently to the stationary
case or relaxing thresholds. The CLI's generic
`plot_Hamiltonian(with_powers=True)` continues to exercise the
collective port-power code; its **unweighted port powers must not
be substituted** for the physical h-weighted energy identity.

## 20 November 2026 — rest gate: exchange-normalization false FAIL corrected

Native archive `sandbox/outputs/v130-issue6-viscous-rest-10c5d22c-6KsOMm1G/`
(source `10c5d22cd6e214b2ef9729cbc125a20769a4fab7`)
contains **successful** physical PETSc TS solves on MPI 1/2/4:
10/10 accepted steps, no SNES failure, no rejected steps,
rest equilibrium at `H0=3.90330691493384`, and a full-state
MPI agreement at end of `~3.6e-17` scaled L2 norm.
The previously unconditional 5% **exchange-scaled** defect gate
was mathematically inappropriate here: true work is `~1e-26`,
viscous dissipation `~1e-33`, and Hamiltonian differences
`~1e-14` are floating-point cancellation on a nonzero background.
Normalizing by `max(actual exchange,1e-14)` produced false
31%/58%/100% failures despite `max|B|/H0` being only
`7.96e-16`, `1.48e-15`, `2.62e-15` respectively.

The comparator now **explicitly** distinguishes `--expect-rest`
from `--require-activity`. In rest mode it **retains and reports**
the ill-conditioned observed ratio as a diagnostic, but does not
apply the 5% physical-exchange gate. Instead it requires relative
rest activity (Hamiltonian change, kinetic energy, boundary power,
dissipation, integrated work) `<=1e-10`, signed energy-balance
residual relative to `max(1,|H0|)` `<=1e-11`, constant
height coefficients to relative `1e-10`, true positive TS
completion, volume balance and strict MPI parity. Nonzero-flow
`physical` mode **retains the unchanged** 5%-of-actual-exchange
limit and activity requirement. Synthetic tests ensure rest mode
rejects both disguised kinetic flow and nonzero energy drift.

No GetFEM/SCRIMP physics, native solver, MPI semantics, or
previous numerical tolerance was altered. In particular a PASS
at rest is not a PASS for 2D viscous nontrivial physics or
long-time power balance. To reclassify the **existing** archived
results without rerunning native computations:

```bash
python sandbox/compare_v130_issue6_viscous_mpi.py \
    sandbox/outputs/v130-issue6-viscous-rest-10c5d22c-6KsOMm1G \
    --cases open --final-time 0.02 --expect-rest
```

**Comparator detail (20 Nov):** `relative_balance_exchange` is also
excluded from MPI 1/2/4 **parity comparisons in explicit rest mode**,
because the archive shows `0.31/0.58/1.0` for this ill-conditioned
quantity. All physical arrays, full endpoint states, volume, absolute
Hamiltonian/port-power balance, and individual rank-port traces are
**still compared**. In nontrivial `physical` mode the relative-to-exchange
metric remains mandatory for both acceptance and MPI parity.

## 20 November — actual rotating SWE failed at nonlinear BDF2 initialization (MPI 1)

The first physical #6 campaign, source `50b8e41be0beae483df81ecf9fb2ef68ce9a80d5`,
immutable `scrimp-mpi:v1.3.0-rc3` image, failed **before main TS integration**,
on the first **open-loop MPI=1** case. Original archive:
`sandbox/outputs/v130-issue6-viscous-physical-50b8e41b-mLLBUGBV/`.
The full-state norm rose from `110.087965` to `59994902.3`
(amplification ~544,976) during the fixed BDF2 consistency step
`t=-1e-9` to `0`, and the unchanged guard (`max factor=100`)
correctly rejected this invalid initial state. No physical
energy/port-power or MPI 2/4 result was generated by this campaign.
This does **not** challenge the completed mass-Jacobian #27.

Independent geometric analysis of the **archived native `Disk.msh`**
with 225 triangles and 33 straight boundary segments confirms
the circular initial momentum `p0/rho=0.2[y,-x]/(r+1e-6)`
does **not** have exactly the prescribed polygon-facet trace
`U_wall=0.2 Rot*Normal`. At segment endpoints,
`max |p0/rho-U_wall|=0.0190327569` (9.52% of 0.2);
the maximum unwanted normal component is `0.0190111897`.
At mid-facets the discrepancy is only about `2e-7`.
This is confirmed **geometric incompatibility**, not yet a
demonstrated complete explanation for the enormous BDF2 norm.
Zero initial algebraic efforts/control variables are another
plausible source of an initial DAE impulse.

New **diagnostic-only** runner on the frozen Docker environment:
`scripts/run-v130-issue6-init-probe.sh`.
It enables existing `SCRIMP_INIT_DYNAMIC_STATE_AUDIT=1` and
`--newton-block-trace` and preserves the original physical weak
forms, nonlinear mass Jacobian, 100× norm guard and MPI calibration.
All probes use one rank to isolate initialization independently
of MPI; every probe has a distinct preserved host-owned directory.
Modes: `baseline` (original), `dt1e-7` (only initialization
step modified), `rigid-facet` (smooth rotating momentum, physical
facet boundary), `warm-ep` (diagnostic algebraic guess),
`unphysical-matched` (`Rot*X` polygon-wall trace, deliberately
**unphysical**), `from-rest` (different initial-value experiment,
not the original rotating p0). Neither a successful altered-data
probe nor a forcefully skipped initializer qualifies #6.

Begin with the original `baseline` to capture block norms and
the exact h/p changes, then use `dt1e-7` and `rigid-facet`
to discriminate inverse-step amplification from geometric
initial/boundary incompatibility. Do not increase the guard
or silently switch to `--no-init-step`. The nonlinear algebraic
DAE initializer (#28) remains separately unqualified.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue6-init-probe.sh baseline
```

## 20 November 2026 — priority change: #28 before further nontrivial #6 MPI work

The native baseline diagnostic at source `e5d4e6e32b5c912eccd15786987ab1ab79fcb7f2`
(`sandbox/outputs/v130-issue6-initial-baseline-e5d4e6e3-KXzesS0i/`)
confirms a defective rotating viscous-SWE **initial consistency** step,
not a main time-stepping or MPI deadlock: the original negative-time BDF2
initializer at `dt=1e-9` produces `||Y||=5.99949e7`, while the physical
height stays at 5 and the prescribed initial momentum changes by
`||p_after-p0||/||p0||=0.760984859`. SCRIMP correctly rejects the
state with its factor-100 amplification guard.

**Decision:** issue [#28](https://github.com/g-haine/scrimp/issues/28)
is now the prerequisite to resuming the **nontrivial rotating** viscous
Saint-Venant qualification under [#6](https://github.com/g-haine/scrimp/issues/6).
Do not launch additional MPI physical/energy campaigns until
fixed-differential initial feasibility and DAE index are understood.
The already validated stationary MPI 1/2/4 gate and nonlinear
mass-Jacobian #27 are unaffected.

**Mathematical scope:** for `z=(z_d,z_a)` with prescribed `z_d(t0)=z_d0`,
the candidate constrained solve uses
`G(z_a)=F_a(t0,z_d0,z_a,0)=0` with the *current* native Jacobian
`J_aa = dF_a/dz_a`. This projection is not a general-purpose cure:
for boundary velocity constraints with Lagrange multipliers the
semi-discrete system can be higher-index or rank deficient; hidden
differentiated constraints may prevent a suitable algebraic-only
projection. The faceted-wall prescribed velocity and the original
circular rotating initial momentum are not exactly compatible;
the solver must reject mathematically inconsistent initial data
instead of altering `p0`, boundary conditions or Hamiltonian.

**Ordered acceptance gates:**

1. Structural rank/consistency audit at **unaltered** original
   `h0,p0`: zero algebraic mass rows, index-1 eligibility, correct
   algebraic DOF partition, rank/conditioning of `J_aa`, residual
   after constitutive algebraic guesses, and boundary constraint
   compatibility. Explicitly report if original data are unsatisfiable.
2. Develop a **separate** nonlinear MPI-safe, safeguarded Newton
   initializer for a demonstrably regular supported DAE class,
   using production residuals/Jacobians and a fixed differential
   initial vector; validate directional FD, descent, rollback and
   collective rejection. **Do not** simply remove the linearity
   guards on the existing native-qualified linear initializer.
3. Native positive and negative nonlinear tests, MPI 1/2, exact
   bitwise preservation of differential coefficients, nontrivial
   physical port work and no silent fallback to pseudo/BDF.
4. Only then resume #6 with an explicitly compatible physical
   rotating initial condition, original weak forms, and MPI 1/2/4
   h-weighted boundary power and viscous dissipation audits.

**Release policy:** this dependency is for the problematic nontrivial
viscous-SWE qualification. It is **not** a claim that #28 already
works or that a full generic nonlinear initializer must enter v1.3.0.
If it cannot be natively qualified within release scope, retain
the linear-only algebraic default and explicitly reject unsupported
nonlinear initializations; do not weaken the scientific gate to
meet a date.

## 20 November 2026 — proposed correction via *geometric* P2, not altered wall velocity

After stage #28.1 fixed-differential-data LSQR, the native dataset
`sandbox/outputs/v130-issue28-feasibility-9082bbee-1vi0TtA9/`
reported original omega=0.2 minimal algebraic residual
`||F_a||~0.02656138`, equally for MPI 1/2, with
the dominant e_p residual `~0.02646368`.
The separate stationary control had a far smaller
LSQR residual `5e-9--7e-9` but incorrectly exceeded the
very strict `~9.47e-10` threshold; its label
`INCOMPATIBLE_AT_TOLERANCE` must **NOT** be taken
as proof that a resting equilibrium does not exist.

**Key distinction:** the experiment had field FEM P2 but
GEOMETRIC order 1 (straight polygon edges). On the actual
polygon, a circular tangent initial `p0` has a small
normal trace; simply replacing `U=omega Rot*Normal`
with `U=omega Rot*X` would match some initial traces
but introduce **nonzero wall penetration**, so it is not
an acceptable physical correction.

An **opt-in** `geometry_order=2` parameter now exists on
`scrimp.utils.mesh.Disk` and in the viscous-SWE sandbox.
Only when selected, the Gmsh mesh is promoted with
`model.mesh.setOrder(2)` *after* mesh generation/refinement.
The existing `geometry_order=1` behavior is unchanged by
default; the PDE, boundary velocity law `U=omega Rot*Normal`,
FEM field orders, parameter values, pressure/viscosity
and initial rotating `h0,p0` are UNMODIFIED.
The quadratic isoparametric boundary approximates the
ideal circle, but is not an exact circle: **do not**
assume exact DAE consistency from geometry order alone.

New bounded native diagnostic:
`scripts/run-v130-issue28-P2-geometry.sh`.
It pins the same immutable rc3 Docker image, mounts the
changed mesh generator and sandbox modules read-only, and
runs **P2 geometric order**, separately for stationary and
original rotating p0 (both MPI 1 and MPI 2).
It reuses the native fixed-h0,p0 LSQR feasibility probe;
the reports now show the requested order, imported GetFEM
geometric point counts / element counts / transformations,
true original `IFunction` residuals, and an **informational**
comparison to P1 original `0.02656138`.
Full physics, historical BDF2 pseudo initialization and
main TS steps are **not** executed here.

The runner checks MPI equality and explicit state rollback.
It does **not** assert that P2 cures the initial-state
incompatibility; its primary qualification is whether
GetFEM really retains curved P2 geometry and what
the resulting native constraint residual is.
If the original tangential initial data remain incompatible,
we must diagnose and explain the specific discrete constraint
and advise a compatible p0 or a smooth startup from rest.
No automatic alteration of user-prescribed initial data.

Source commits:
`7875617c8e6ba6193f0cc14ffe9ec9e6594bcebc`
(opt-in curved mesh) and
`a9ebc1766a3814c2215681ae74af7385ece0aa19`
(native P2 feasibility runner and compatibility suggestions).
**Native P2 evidence remains PENDING.**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-P2-geometry.sh
```

## 20 November 2026 — correct tangent sign on Gmsh P1/P2, GetFEM check PENDING

The user asked whether the initial rotation and imposed wall tangent
have accidentally opposite orientations. **Static source reading says
NO:** `Rot=[[0,1],[-1,0]]`, `Tangent=Rot*Normal`, and
`p0/rho=omega*(y,-x)/(sqrt(x*x+y*y)+1e-6)` all rotate clockwise.
For an outward unit radial normal at (1,0), both point down.

Independently computed from the native user archives
`v130-issue28-feasibility-9082bbee-1vi0TtA9.zip` and
`v130-issue28-P2-geometry-7286bec2-Rp90btp3.zip`:
3-point / 2-point P2/P1 curved/straight Gmsh edge map,
12-point Gaussian integration per boundary element,
outward radial-choice geometric normal from the boundary
parametric derivative, all **33 original boundary edges**.
For omega=.2 original circular `p0`, the relative boundary
`L2` norm `||p0/rho-omega*Rot*n||/||p0/rho||` is:

- **P1 geometry:** `0.0550176860263` (~5.50%).
- **P2 geometry:** `0.0000696599723462` (~0.00697%).
- For the **opposite tangent sign**, the corresponding relative
  mismatch is nearly **2.0 (200%)** on both meshes.

Thus curved Gmsh P2 reduces the **geometric** initial
tangential-trace mismatch by ~790x, whereas the previous
GetFEM native algebraic LSQR residual increased from
`0.0265613800268` (P1) to `0.0277336388365` (P2).
The latter is **not explained by a reversed analytic tangent**
or simply by the straight-sided polygon approximation.
Note: Gmsh-based geometry normals are **not** the same thing as
measuring the actual GetFEM GWFL `Normal` used in the weak bricks.

New read-only native sign diagnostic
`sandbox/probe_v130_issue28_normal_orientation.py` +
`scripts/run-v130-issue28-native-normal.sh`:
uses the **actual GetFEM generic boundary integral** on the
original model, with original rotating p0 and the production
`Tangent` macro. Directly measures both signs,
`(p/rho).Normal`, `X.Normal` and the imported P1/P2
geometry in MPI 1/2. No LSQR or TS time stepping.
It can isolate a GetFEM mesh-region normal sign/assembly issue
which is not visible from the exported Gmsh coordinates.
The native GetFEM normal test is **PENDING**; do not claim
that the Gmsh geometry-only numbers prove the sign seen
by the actual weak-form assembler.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-native-normal.sh
```

**Sequence retained:** complete this boundary/GWFL normal sign
and geometry investigation first; *afterwards* use original
`examples/shallow_water.py` for a separate nonlinear fixed-dynamic
initialization study, as requested. No changes to model physics,
initial p0, existing LSQR thresholds, or BDF2 safety guard.

## 20 November 2026 — ROOT CAUSE: GetFEM GWFL column-wise Rot reversed wall actuation

**Native evidence PASS, before correction:** user uploaded
`v130-issue28-native-normal-b903865b-91WeMgdE.zip`
(original `b903865bee829d18a329928fa04dd68e0b536eb7`,
unchanged image `scrimp-mpi:v1.3.0-rc3`).
Native GetFEM integration using **the actual production
`Normal`, `Tangent=Rot*Normal`, original p0**, MPI 1 and 2:

| Geometry | `||p/rho - omega*Tangent||/||p/rho||` | `||p/rho + omega*Tangent||/||p/rho||` | `||p/rho.Normal||/||p/rho||` |
|---|---:|---:|---:|
| P1 straight | 1.9992482297 | 0.0549180824 | 0.0548807074 |
| P2 curved | 2.0000023663 | 0.0000697293 | 0.0000696876 |

The GetFEM outward normal itself is **correct**:
`integral(X.Normal)` = `6.2452910639` (P1) and
`6.2831681228` (P2), near the expected `2*pi`.
MPI1/2 are equal to floating-point precision.
This **overturns** the previous assumption that the old
`Tangent` macro matched `p0` inside GetFEM. The geometric
P1/P2 edge-normal analysis was not sufficient to determine
the sign of the literal `Rot` as interpreted by GWFL.

**Exact root cause:** GetFEM GWFL writes a matrix
`[[a,b],[c,d]]` with the nested vectors as **COLUMNS**.
Thus the previous `Rot = [[0,1],[-1,0]]` represents the
**actual** linear operator `[[0,-1],[1,0]]`, i.e. CCW.
Original `p0/rho = omega*(y,-x)/(r+1e-6)` is CW.
The wall law therefore had the opposite tangential sign,
regardless of P1/P2 geometric accuracy. This explains why
curved P2 did **not** reduce the previous algebraic residual.
Reference: https://www.getfem.org/userdoc/gasm_high.html
(explicit matrices section).

**Corrected source commit:** `6db7a79cfe6f7b20d854c4b45913922f10fda0a9`.
The experimental viscous-SWE `Rot` macro is now written
`[[0,-1],[1,0]]` (GWFL columns), meaning the intended
row-wise clockwise matrix `[[0,1],[-1,0]]`.
Original p0, h0, PDE weak forms, viscous terms, normal
impermeability and all SCRIMP solver internals are unchanged.
The nonlinear gyroscopic form
`Gyro(p) = div(Rot*p) * Rot` is invariant under a **global
`Rot -> -Rot` sign flip**: it contains two factors of Rot.
The correction therefore flips the actual wall tangential
velocity but leaves the gyroscopic physics invariant.

A source-only test `tests/test_v130_issue28_rot_convention.py`
checks the actual column-wise matrix representation, clockwise
tangent at `(1,0)`, unchanged `Gyro`/`Curl2D` macro formulas
and formal sign invariance; CI runs it. The original
`scripts/run-v130-issue28-native-normal.sh` now rejects a
`REVERSED` or `AMBIGUOUS` native sign instead of printing
a warning and exiting 0.

**Old LSQR numbers are obsolete as physical evidence:**
the pre-fix `||F_a||_min` ~0.02656 (P1) and ~0.02773 (P2)
were assembled with a *wrong-sign wall velocity*.
Repeat fixed-state P2 feasibility with the correctly signed
wall before concluding anything about residual boundary
incompatibility, and retain the separate stationary control
and MPI 1/2 check. Even if fixed-state `F_a` becomes small,
`Y` remains a null column of `J_aa` in the open-loop system;
a complete index-aware DAE initializer is still unqualified.

**Ordered native experiments:**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-native-normal.sh
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-P2-geometry.sh
```

**Do not combine workstreams:** After the native
geometry/initial-boundary compatibility issue is interpreted,
and ONLY THEN, investigate nonlinear constrained
initialization using `examples/shallow_water.py`, as
explicitly requested. No native post-fix data yet.


## 20 November 2026 — post-fix GetFEM GWFL matrix convention: native P1/P2 evidence

**Important Python user pitfall:** the GetFEM Generic Weak-Form Language
(GWFL) writes nested explicit matrices as *column vectors*, whereas
Python/NumPy list-of-lists conventions use *rows*.

- **Authoritative upstream reference:**
  [GetFEM GWFL — Explicit matrices](https://getfem.org/userdoc/gasm_high.html#explicit-matrices).
- **User-facing Sphinx guide:** `docs/source/getfem_gwfl.rst`
  (reachable from the website index, getting started, initialization
  page, and repository README).
- **Incorrect for clockwise rotation:** GWFL `[[0,1],[-1,0]]`
  (counterclockwise as an actual matrix).
- **Correct clockwise:** GWFL `[[0,-1],[1,0]]`
  or the more readable row/semicolon form `[0,1;-1,0]`.
- Actual GWFL `Normal` is outward; in the unit disk
  `Rot*Normal` now agrees with the original
  prescribed clockwise velocity `p0/rho`.

**Native post-correction qualification of sign and reduced algebraic
feasibility (not of complete DAE initialization).** Branch source
`120c01bd8aa206ea0986ba876d65d35489a825b6`, immutable
Docker `scrimp-mpi:v1.3.0-rc3`, true GetFEM/PETSc MPI 1 and 2:

| Quantity | P1 geometric mesh | P2 curved geometric mesh |
|---|---:|---:|
| Actual GWFL tangent orientation | ALIGNED (MPI 1/2) | ALIGNED (MPI 1/2) |
| Relative aligned trace mismatch | `0.0549180824` | `6.97292855e-5` |
| Relative normal leakage | `0.0548807074` | `6.96876238e-5` |
| GWFL integral `X.Normal` | `6.24529106` | `6.28316812` |

The P2 geometric transformation was actually imported as
`GT_PK(2,2)`, not simply simulated with P2 state FEM on straight edges.

**Independent original P2 algebraic LSQR**, fixed `h0,p0`:
before the GWFL sign correction (source `7286bec`)
`||F_a||_min=0.0277336388`; after the correction (source
`120c01bd`) `4.1850871648724585e-8` (MPI 1) and
`4.183695287884114e-8` (MPI 2): a reduction of
**more than 660,000 times**. The post-fix rotating feasibility tolerance
is `2.60963115021321e-9`, so the strict label is still
`INCOMPATIBLE_AT_TOLERANCE` and **must not be advertised as PASS**.
The rotating per-block projected `e_p` norm is
`4.101738589874305e-8` on MPI 1.

At rest (P2), the original known equilibrium yielded a
diagnostic-only LSQR residual of `4.746567828624845e-9`
(MPI 1) and `4.721221465777258e-9` (MPI 2), with a stricter
acceptance tolerance `9.507859931505774e-10`.
**A fail against this numerical stopping threshold is not a proof
of physical incompatibility**, particularly for a resting
state already qualified as an equilibrium. Precision/conditioning
and subsequent high-index initialization need separate treatment.

Evidence: user-supplied archive `v130-issue28-.zip`,
the host native post-fix sign output directory
`sandbox/outputs/v130-issue28-native-normal-120c01bd-sq5HgESS/`,
and the native post-fix feasibility output
`sandbox/outputs/v130-issue28-P2-geometry-120c01bd-sKH7lgKN/`.

**Do not close #28 yet:** although the tangential sign error is fixed,
the open-loop traction `Y` remains absent from the algebraic equations
(`dF_a/dY=0`), so a complete consistent initial traction generally
needs a different, possibly constraint-differentiated formulation.
As agreed with the user, test `examples/shallow_water.py` as a
separate nonlinear-initializer experiment **only after** the current
boundary/geometric work is concluded.

## 20 November 2026 — #28 exact stationary baseline NATIVE PASSED (MPI 1/2, geometry P1/P2)

**Host evidence:** `v130-issue28-exact-rest-b556f2ea-3E78L59H.zip`.
Source `b556f2ead84397b7158b1755d2dfa72ede8122df`,
immutable `scrimp-mpi:v1.3.0-rc3` OCI revision
`1d6d0c90115213298047b28237a7f38087563a35`.

The diagnostic inserted the independently known exact **algebraic**
stationary state WITHOUT LSQR, pseudo initialization or a physical TS:
`h0=5` and `p0=0` (prescribed original differential coefficients),
`e_h=0.5`, `Y_n=-0.5`, `e_p=U_n=U=Y=0`. It evaluated the
real production `IFunction(t0,z,zdot=0)` and rolled back all states.

| Geometry | MPI | `||F_a||_2` | `||F_d||_2` | Wrong `Y_n=+0.5`: `||F_a||_2` |
|---|---:|---:|---:|---:|
| P1 | 1 | `1.4633950546719204e-17` | `2.6084923555407635e-16` | `4.0700566141663` |
| P1 | 2 | `1.470126998300503e-17` | `2.6084923555407635e-16` | `4.0700566141663` |
| P2 | 1 | `1.4763422151112202e-17` | `2.699640468033266e-16` | `4.075223213705125` |
| P2 | 2 | `1.4744250640694402e-17` | `2.699640468033266e-16` | `4.075223213705125` |

**The stationary reference passes** (all < `1e-11`), with a
nonvacuous negative control. This proves that previous resting
`INCOMPATIBLE_AT_TOLERANCE` classifications (`||F_a||~4.7e-9`)
were **numerical LSQR false negatives**: identical fixed h0,p0
have an exact solution of the assembled FE constraints.
Do not loosen production tolerances merely to hide the LSQR failure.

**The #28 boundary geometry/sign *subgate* is closed as a
diagnosed and corrected GWFL literal error**, including true
P1/P2, MPI 1/2 sign and resting-equilibrium validation.
This does NOT close issue #28: for original nonzero rotation,
the corrected fixed-h,p algebraic LSQR residual remains
`~4.18e-8`, above the strict diagnostic threshold; the
traction `Y` is a null algebraic Jacobian column and consistent
initial `Y0` may require differentiated boundary constraints.
No complete rotating DAE initializer has been qualified.

**Next, as explicitly ordered by the user:** study the original
`examples/shallow_water.py` for nonlinear consistent
initialization on a straight-edged rectangular domain.
Verify the algebraic Jacobian and the existing
`explicit=True` frozen gyroscopic treatment. Preserve h0,p0.
Nonlinearity in h,p and state-dependent mass do not
automatically mean that the algebraic equations are
nonlinear in their algebraic unknowns; distinguish the
two levels before claiming a nonlinear algebraic Newton.
