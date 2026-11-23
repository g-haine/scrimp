# SCRIMP v1.3.0 — Issue #6: active rotating SWE MPI 1/2/4 parity

**Status: NATIVE QUALIFIED / 14/14 PETSc runs PASS, MPI 1/2/4 parity PASS.**
Prepared and executed on 23 November 2026 on
`release/1.3.0-preparation` at source
`669540ef4805000612e48a1e81489e4590752bd8`.
The already qualified cold checkpoint 1/2-rank campaign from
`2119bc55` remains immutable and PASSED.

## Native results — 23 November 2026: COMPLETE PASS

Immutable original host evidence:
`sandbox/outputs/v130-issue6-active-mpi4-669540ef-JJl2w1Tm/`
and standalone ZIP with exactly that basename. All original
native logs, both mesh cases, 14 production DPHS checkpoints,
14 complete full-FE NPZ trajectories, 34 rank-local original
generic-port-power JSON files and root-only physical plots are
preserved. User-provided ZIP: **167 entries**, all pass CRC, original
data independently read/recomputed (not merely accepting the
campaign-comparator verdict).

Two meshes, original gamma=.25 physical rotating wall, curved
P2 domain, dt=.002: 14/14 native PETSc runs, zero rejected steps
or SNES failures, valid constraints, completed CN steps.

| Physical mesh | Complete FE DOFs | Continuous max whole-path L2 MPI1 vs MPI2 | Continuous max whole-path L2 MPI1 vs MPI4 | Maximum physical observables MPI1 vs MPI2/4 | Maximum generic port-power MPI1 vs MPI2/4 |
|---|---:|---:|---:|---:|---:|
| h=.2 | 3300 | 2.57507288e-13 | 2.1002e-13 | 4.79616347e-14 | 3.40005801e-16 |
| h=.1 | 10272 | 1.6614e-13 | 1.68086923e-13 | 2.32702746e-13 | 1.52655666e-16 |

Native cold restarts **1→4** and **4→1** match the full
uninterrupted trajectory, with bitwise source FE restoration:

| Mesh | Restart 1→4 max whole-path FE L2 | Restart 4→1 max whole-path FE L2 |
|---|---:|---:|
| h=.2 | 1.9078e-13 | 2.2930e-13 |
| h=.1 | 1.6358e-13 | 1.6878e-13 |

The **aggregate, independently reproduced** source-matched metrics:

- Whole-path MPI 1/2/4 FE discrepancy max:
  `2.5750728773395547e-13`.
- Original h-weighted physical observables discrepancy max:
  `2.32702746103467e-13`.
- Original unweighted generic port-power discrepancy max:
  `3.400058012914542e-16`; on each individual MPI run,
  all rank-local generic powers agree exactly to the saved precision.
- Four cross-rank cold restarts 1↔4 passed (two meshes), with
  original full FE state, time, physical powers and Hamiltonian.
- Actual nonzero generic port power activity (representative
  max magnitude h=.2: `0.2231`; h=.1: `0.1216`).
  The rank-local assembly test is not a vacuous resting solution.

Final campaign manifest reports
`PASS_NATIVE_COMPATIBLE_ROTATING_SWE_MPI_1_2_4`.

**Software verdict:** SCRIMP reproduces the original finite-element
trajectory, physical power observables and generic algebraic port
powers on MPI1, MPI2 and MPI4, with working collective plotting
and no concurrent rank-file overwrites in these outputs. The
production cold checkpoint also preserves the same simulation
for 1↔4 changes of MPI process count.

**Strict scope:** the experiment runs only t=0..0.024,
before normal feedback starts at t=75. The actual tangential wall
motion is nonzero, but this does not qualify normal-feedback-active
long-time dynamics, ParaView exports, adaptive multistep TS history
or a stationary rotating solution. Keep the broader issue #6 OPEN.
The model's absolute nonzero CN energy defect remains a
separate numerical-analysis question and has not been reclassified
as a SCRIMP MPI error.

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

**This criterion is now SATISFIED by the archived 14 native jobs and independent ZIP re-audit.**
Do not close issue #6 automatically: other scope includes
longer-time, active normal feedback and parallel output
scenarios not tested here. No model, production solver,
GetFEM bridge, or numerical tolerances are changed.
