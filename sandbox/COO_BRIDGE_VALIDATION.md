# Experimental PETSc COO bridge — validation log

Branch: `feature/getfem55-mpi-container-scalability` (PR #24).

The baseline GetFEM matrix conversion in `scrimp/utils/linalg.py` transfers
global rows with `Mat.setValuesIJV(..., ADD_VALUES)` and calls
`Mat.assemble()`. In GetFEM 5.5 MPI-2, contributions can target other
PETSc-owned rows and **must not** be discarded.

## Measured discovery (2026-11-17)

The existing `TimingProfiler` isolated `setup.bridge.mass` and
`setup.bridge.stiffness` on `examples/wave.py`, `h=0.08`, 12,682 DOF,
2 MPI ranks, fixed-step CN. Baseline setup bridge values were approximately
0.723 s (mass) and 0.705 s (stiffness). PETSc assembly, not GetFEM
copy/transposition, dominated.

The targeted sandbox prototype `probe_petsc_coo.py` used native PETSc
`Mat.setPreallocationCOO()/setValuesCOO()` with the *same* distributed
GetFEM contributions. It added explicit zero diagonal slots to match
SCRIMP's existing diagonal initialization. Each trial compared complete
PETSc **stored CSR** arrays, including indices, slots and values.

| Block | Trial | Historical IJV + assemble | COO incl. preparation | Speedup | CSR |
|---|---:|---:|---:|---:|---|
| mass | reference → COO | 0.610357 s | 0.003544 s | 172.2× | equal |
| mass | COO → reference | 0.283650 s | 0.003299 s | 86.0× | equal |
| stiffness | reference → COO | 0.689616 s | 0.010451 s | 66.0× | equal |
| stiffness | COO → reference | 0.501143 s | 0.010324 s | 48.5× | equal |

These are **isolated matrix assembly** timings, NOT complete `DPHS.solve()`
speedups. The reference is sensitive to ordering/warm-up. The COO
prototype started with newly created PETSc matrices, so it did **not**
establish safety or performance of COO on SCRIMP's already initialized
matrix handles. That is the next required gate.

## Opt-in production integration (unvalidated in native MPI)

- Default: `SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=ijv` (unchanged).
- Experimental: `SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo`.
- The COO path is called **only** from `DPHS.assemble_mass` and
  `DPHS.assemble_stiffness`, the initial linear setup operators.
- Nonlinear tangent matrices, potential structural changes during Newton,
  and all other matrices stay with the established IJV/ADD_VALUES
  transfer, regardless of this option.
- Replicated-matrix semantics are unchanged; their separate
  `SCRIMP_GETFEM_BRIDGE_REPLICATED` option is independent.
- COO preallocation includes the original owned diagonal zero slots.
- `setValuesCOO(..., INSERT_VALUES)` **sums repeated coordinates**
  internally, then replaces previous matrix coefficients and completes
  final assembly. No additional `Mat.assemble()` should be called.
- The original PETSc.Mat handle is reused; no assignment of new
  PETSc matrices or TS/SNES/KSP alterations.

## Native in-place COO test results (2026-11-17)

The user ran `test_real_getfem_distributed_coo_matches_existing_petsc_mat`
together with the distributed-contribution regression inside
`scrimp-mpi:local`, on the actual GitHub branch mounted read-only.
Both tests passed at 1 and 2 MPI ranks, on **every rank**. PETSc Mat
identity, repeated in-place transfers and stored-CSR equality therefore
have targeted native evidence for these ranks. This is not yet a
full-PDE or scaling result. Tests at 4 and 6 ranks remain outstanding.

Run the paired Docker-only acoustic benchmark and complete complex
trajectory parity on the **same commit** via:

```bash
bash sandbox/run_acoustic_bridge_pair.sh 0.08 2
```

The helper uses the actual `scrimp-mpi:local` image, the historical
`ijv` default and explicit `coo` opt-in, fixed 20-step CN parameters,
and `check_acoustic_state_parity.py` for every saved coefficient.
Host Python reads JSON reports only; it never simulates.

## Native end-to-end 2-rank acoustic validation (2026-11-17)

User-confirmed full Docker A/B with `bash sandbox/run_acoustic_bridge_pair.sh
0.08 2` on the *same* source commit `49bf3d59`, GetFEM 5.5,
`h=0.08`, 12,682 DoFs, 20 fixed-step CN steps, `t_f=0.04`,
`SCRIMP_GETFEM_MATRIX_MODE=distributed`, and automatically chosen MUMPS
reset policy `never`:

| Metric | Historical IJV | Experimental COO |
|---|---:|---:|
| `setup.bridge.mass` | 0.725561 s | 0.005618 s |
| `setup.bridge.stiffness` | 0.785340 s | 0.018183 s |
| Sum of those two setup bridge events | 1.510901 s | 0.023801 s |
| `DPHS.solve()` max rank wall | 2.159889 s | 0.658418 s |
| `main.petsc.ts_solve` | 0.509458 s | 0.494820 s |
| Max per-rank peak RSS | 409.9 MiB | 411.9 MiB |

Observed end-to-end speedup **3.2804×** (69.52% less wall time);
the two measured bridge transfers accelerated **~63.5×** in aggregate.
This is a **single paired trial**, not a statistically robust median;
the timing gain must be checked over more ranks/meshes before release.

The saved trajectories `state_ijv.npz` and `state_coo.npz` were independently
re-inspected from the user's archived run, not merely inferred from the
script's PASS status. Both have 20 timestamps and 20×12,682 coefficients
(dtype `float64`; the validator also supports complex arrays, but this
particular example is real). Timestamps are bitwise equal. Maximum
absolute coefficient difference: **6.372680161348399e-14**; relative
trajectory L2 difference: **7.476624355814846e-15**. Within each run,
the MPI rank-to-rank coefficient difference is zero. All 20 steps were
accepted; zero rejected steps and zero nonlinear solver failures; PETSc
and final saved time are both 0.04. Physical variable norms match within
round-off. No simulation errors occurred in the archived logs.

**No memory saving claim**: COO had +~2 MiB peak per-rank RSS in this
single run. No claim is made for nonconstant graph nonlinear Jacobians,
complex-valued states, other MPI sizes, dissipative Saint-Venant or
closed-loop heat-wave dynamics.

## Native end-to-end 2-rank acoustic, h=0.04 (2026-11-17)

User executed `bash sandbox/run_acoustic_bridge_pair.sh 0.04 2` on
source commit `c08b02c9` inside `scrimp-mpi:local`, GetFEM 5.5.
**47,942 DoFs**, 2 MPI ranks, 20 accepted fixed Crank–Nicolson steps,
`dt=0.002`, `t_f=0.04`, MUMPS policy `never`:

| Metric | Historical IJV | Experimental COO |
|---|---:|---:|
| `setup.bridge.mass` | 10.428504 s | 0.021403 s |
| `setup.bridge.stiffness` | 15.020396 s | 0.074983 s |
| Two setup bridges combined | 25.448900 s | 0.096386 s |
| `main.petsc.ts_solve` | 1.809227 s | 1.824937 s |
| `DPHS.solve()` max rank wall | 27.722780 s | 2.401666 s |
| Max rank peak RSS | 511.648 MiB | 515.570 MiB |

**End-to-end speedup 11.543×** (91.337% shorter solve),
setup-bridge speedup **264.03×**. One run per method, therefore no
median, confidence interval or broad scaling claim.

The user's archived run `bridge-pair-h0p04-2r-c08b02c9.zip` has
been independently inspected: both saved trajectories have
`times.shape == (20,)` and `states.shape == (20, 47942)`, real
`float64` coefficients. Timestamps are exactly equal.
**Max absolute coefficient difference = 0.0; relative L2 difference
= 0.0** across all 958,840 saved coefficients. Per-rank parity is
also exact; 20 accepted steps, zero rejected, zero Newton failures,
same final-time 0.04 and matching physical variable norms.

This is an unexpectedly strong result for this particular
problem/configuration, NOT a proof that all matrix patterns or
nonlinear models are identical. No reduction in memory is claimed:
peak RSS increased by ~3.92 MiB per rank in this run.

## Native in-place COO tests on 4 and 6 ranks (2026-11-17)

The user ran inside `scrimp-mpi:local` the in-place COO-vs-IJV
PETSc-Mat regression and the distributed contribution regression
at **4 and 6 ranks**. Each test passed on every rank. All algebraic
native-MPI tests have now passed at ranks **1, 2, 4 and 6**.

These are targeted tests, not complete PDE / MUMPS / energy-balance
validation at 4/6 ranks.

## Transparent joint strong-scaling campaign (new)

The scaling runner now crosses two **independent** axes:

- Distributed GetFEM bridge: `SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=csr`
  (historical CSR/IJV + `Mat.assemble()`, `ijv` remains a legacy
  alias), and `coo` (new native COO path on initial linear matrices).
  The **default remains `ijv`**, i.e. historical CSR transfer.
- MUMPS symbolic reset: `never`, `on_change`, `always` as the
  existing `SCRIMP_MUMPS_RESET_POLICY` values. Optional `auto` remains
  the benchmark alias for automatic policy selection (`never` for this
  linear wave case); `auto` must NEVER be assigned to the environment
  variable.

All methods use the same source commit, fixed mesh, fixed CN
parameters, and Docker GetFEM 5.5. New reports record the bridge,
effective MUMPS policy and commit explicitly; older results stay in
`scaling-mumps-h*`, unmodified. The analyzer **rejects** missing
`bridge_method`, mixed source commits, incomplete method/policy/rank
grids, unbalanced repetitions, nonconvergence or discrepant field
norms. Reports and plots identify `T_CSR(p)`, `T_COO(p)`,
`S(p)=T(1)/T(p)` independently per configuration, and
`T_CSR(p)/T_COO(p)` at equal rank and MUMPS policy.

```bash
# Fine-mesh 1/2/4/6 ranks × 3 repetitions × 2 bridges × 1 MUMPS policy
# = 24 Docker simulations, in a dedicated Git-SHA output folder.
SCRIMP_SCALING_TIMEOUT_S=1800 bash sandbox/run_acoustic_mumps_scaling.sh 0.02 sweep both never 3

# Only four runs: 1 and 6 ranks, one repetition, never policy
bash sandbox/run_acoustic_mumps_scaling.sh 0.04 pilot both never

# Explicit full comparison: 4 ranks × 5 reps × 2 bridges × 2 policies
# = 80 runs; pilot may be reused for the matching (never, 1/6) runs
bash sandbox/run_acoustic_mumps_scaling.sh 0.04 full both never,always

# Also include MUMPS 'on_change' (120 total configurations)
bash sandbox/run_acoustic_mumps_scaling.sh 0.04 full both never,on_change,always
```

The runner uses Docker for **both simulations and plotting**.
A plain Python 3 interpreter on the host only checks JSON metadata.
Run chart validation unit tests inside Docker via
`python -m unittest -v tests.test_acoustic_bridge_scaling`.
Chart outputs: `comparaison_bridge_speedup.png`,
`comparaison_bridge_times.png`, `comparaison_bridge_coo_gain.png`,
and `comparaison_bridge_summary.csv`.

## Validation gates before considering default activation

1. Native GetFEM 5.5 / complex PETSc 3.26 tests in Docker:
   `test_real_getfem_distributed_coo_matches_existing_petsc_mat`
   at 1, 2, 4 and 6 MPI ranks. It performs two transfers to the
   *same* matrix initialized by SCRIMP and checks full stored CSR.
   **If any test fails, do not enable the experimental option in a solve.**
2. One 2-rank, h=0.08, 20-step acoustic CN run with `ijv` and `coo`
   on the **same commit** using `sandbox/bench_acoustic_2d.py --profile`.
   Compare TS endpoint, failures, final variable norms, total wall,
   setup bridge phases and peak per-rank RSS.
3. Use existing `sandbox/check_acoustic_state_parity.py` to compare
   *complete complex trajectories* for IJV vs COO.
4. Repeat the 2-rank equivalence tests at h=0.04, then MPI 1/4/6,
   prior to the full repeated strong-scaling campaign.
5. Validate `examples/heat_wave.py` closed-loop control behavior
   separately (MUMPS policy change).
6. Validate dissipative, rotating Saint-Venant 2D with spin-up including
   nonlinear solves, power/Hamiltonian and volume balances before release.

No claim is made yet of end-to-end acceleration, memory improvement,
or physical parity for the opt-in production path.
