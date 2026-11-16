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
git switch fix/issue-6-parallel-powers-saint-venant-2d

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
