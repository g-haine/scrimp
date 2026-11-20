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
