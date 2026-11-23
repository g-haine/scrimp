# SCRIMP v1.3.0 — Issue #6: active rotating SWE MPI 1/2/4 parity

**Status: SOURCE ONLY / NEW NATIVE CAMPAIGN PENDING.**
Prepared on 23 November 2026 on `release/1.3.0-preparation`.
The already qualified cold checkpoint 1/2-rank campaign from
`2119bc55` remains immutable and PASSED.

## Required question

SCRIMP must compute **the same discrete physics** and must not
deadlock/corrupt outputs when physical GetFEM and port-power operations
run on 1, 2 or 4 MPI ranks. The original end-to-end path
`DPHS.plot_Hamiltonian(with_powers=True)` is executed *on every rank*
so its collective `compute_powers()` cannot silently succeed on
rank 0 alone. The result is validated against the same original
physical Saint-Venant states and observables, not a fabricated
manufactured PDE or the inactive stationary zero-velocity test.

**Do not confuse two distinct power notions:**

- The original `Port.compute()` returns *unweighted* algebraic
  port powers, with `P[0]` instantaneous and `P[k>=1]` evaluated
  at the Crank–Nicolson mean state. This is useful to establish
  collective MPI execution and rank equality. It is NOT the
  h-weighted SWE energy flow.
- The **physical** original model computes
  `p_normal=h*U_n*Y_n`,
  `p_viscous=h*U.Y`,
  `d_grad=2*mu*h*D(e_p):D(e_p)` and
  `d_div=2*mu*h*div(e_p)*div(e_p)` in GetFEM at the CN
  mean state. These are checked per saved time point, with the
  original integrated external work, viscous dissipation,
  Hamiltonian and volume, against MPI1.
- Issue [#30](https://github.com/g-haine/scrimp/issues/30)
  independently describes `plot_powers()`'s erroneous trapezoid
  over **already midpoint-evaluated** step powers. No gate here
  relies on that plot's integrated numerical balance. The
  physical `H+D-W` defect magnitude is a discretization-model
  observation, not a checkpoint or MPI acceptance threshold;
  **its difference between MPI solutions** is a strict gate.

## Original model, precisely bounded native gate

`gamma=0.25`, curved quadratic disk, P2 fields,
`dt=0.002`, original weak forms, state-dependent mass, implicit
nonlinear gyroscopic term, original boundary sign and
original independent normal-feedback clock. The compatible
`t=0` DAE rest state is accepted exactly and the tangential
wall is smoothly accelerated over `0<t<.02`. Physical
solutions run until `t=.024`. No pseudo/BDF projection
is allowed for checkpoints, and PETSc CN must accept every
step with true `F_a/F_TSCN` tolerances.

**Important limitation:** on this time window the original
normal boundary feedback is still **inactive** (its change
begins at `t=75`). The wall is actively rotating and
producing nonzero viscous boundary power and dissipation,
but this is NOT qualification of a fully active closed-loop
normal-feedback control.

For each of `h=.2` and `h=.1`, run:

| Case | Native MPI | Physical interval |
|---|---|---|
| Compatible seed | 1, 4 | 0 → .02 |
| Uninterrupted physical reference | 1, 2, 4 | 0 → .024 |
| Cold checkpoint resume | 1 → 4 and 4 → 1 | .02 → .024 |

There are **14 independent physical GetFEM/PETSc jobs**, **4
cross-MPI native uninterrupted comparisons**, **4 additional
cold checkpoint-to-uninterrupted comparisons**, **34 independent
rank-specific generic-port files** and **14 root-only PNG physical
power plots**. This focused incremental MPI4 campaign is
**different from** the already completed MPI1/2 checkpoint
8/8 campaign and does not re-run all its permutations.

After every native job the underlying production pilot also
requires the original positive PETSc convergence reason, exact
step count/time, no rejection/SNES failure, full accepted
`z(t_n)`, `F_a < 1e-8`, `F_TSCN < 1e-7`, positive
FE height coefficients and conserved original volume
(<1e-8). The host rejects missing per-rank JSON or absent
root physical plot; the independent comparator checks every
MPI rank's original generic algebraic power against rank0,
then checks MPI2/MPI4 continuous physics against MPI1.
Full states use maximum L2 <=1e-8 and every recorded
power/energy/volume observable uses absolute parity <=1e-8.
True cross-rank cold restart must preserve the complete
source FE vector bitwise at the pause.

Root-only plot `mpi_physical_powers.png` is constructed
from original h-weighted physical GetFEM trace arrays, not
from the incorrect generic `plot_powers()` quadrature.
All output files are written under unique host-owned directories.
Failures retain all logs and native evidence. Source-only CI
checks the verifier's *negative controls*, **not** the MPI
physics itself.

## Reproduce on the already qualified Ubuntu Docker host

```bash
cd ~/Simulations/scrimp
git pull --ff-only

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue6-rotating-mpi4.sh
```

The launcher prints a unique
`sandbox/outputs/v130-issue6-active-mpi4-<source>-<suffix>`
directory. Within it, the host logs for all cases,
`native_cases.tsv`, `h0p2_native_mpi4.json`,
`h0p1_native_mpi4.json`, and aggregate `summary.json`
are immutable evidence. Each referenced native directory
has complete GetFEM full-state NPZ, accepted production
checkpoint, per-rank `mpi_port_powers_rank<N>.json` and
root-only `mpi_physical_powers.png`.

On complete native success, the runner also creates a
**self-contained, separately named ZIP** next to the campaign
folder, `v130-issue6-active-mpi4-<source>-<suffix>.zip`.
It includes the unchanged original native full-state NPZ traces,
production checkpoint arrays, native logs, every rank-local
port-power JSON, physical PNGs, source manifests and comparison
reports. It does **not** move, overwrite or delete any source
evidence. This ZIP may be shared for independent numerical review
without requiring access to the original host paths.

**Do not report this as PASS until the genuine host runs
complete with final `PASS_NATIVE_COMPATIBLE_ROTATING_SWE_MPI_1_2_4`.**
Do not close issue #6 automatically: other scope includes
longer-time, active normal feedback and parallel output
scenarios not tested here. No model, production solver,
GetFEM bridge, or numerical tolerances are changed.
