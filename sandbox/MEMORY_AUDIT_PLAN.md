# SCRIMP: first long-run MPI memory audit (2026-11-17)

## Scope

This is the next task after validating the native PETSc COO bridge.
**No changes to the DPHS physics, time integrator, solver policy or
core implementation** are made in this first diagnostic. We do not
claim a memory leak, or a memory reduction, until the data support one.

The relevant source code in `scrimp/dphs.py`:

- `DPHS.monitor()` gathers the complete PETSc state into a global
  vector and appends a full NumPy copy to `solution["z"]` at `dt_save`
  intervals, **on every MPI rank**.
- `_select_saved_solution_through_final_time()` later stacks the
  saved states with `np.asarray(states)` and filters them. A second
  temporary allocation can occur on the way to the final array.
- `DPHS.allocate_memory()` creates separate `mass`, `nl_mass`,
  `tangent_mass`, `stiffness`, `nl_stiffness`,
  `tangent_stiffness`, and vectors. The solver also duplicates
  the tangent Jacobian into `self.J` before the TS run.
- PETSc TS, KSP, PC and MUMPS may allocate native memory unrelated to
  Python's saved-state arrays.

A saved `float64` global state contains `8 * n_dofs` bytes, per
rank. For 47,942 DoFs and 250 snapshots, that is 91.4 MiB **per
rank** of array payload alone. This is retained output data, not
evidence of a leak.

## Paired long-run experiment

```bash
cd ~/Simulations/scrimp
git pull --ff-only

bash -n sandbox/run_acoustic_memory_audit.sh
bash sandbox/run_acoustic_memory_audit.sh 0.04 0.002 0.5 2
```

This launches **two simulations in `scrimp-mpi:local`**, with the
original `examples/wave.py` physics and these identical settings:

- GetFEM 5.5 MPI-2, 2 ranks pinned to physical cores `0,2`
- mesh `h=0.04` (previously 47,942 DoFs)
- PETSc CN with fixed `dt=0.002`, final time `0.5`, 250 steps
- `PC=LU/MUMPS`, reset policy `never`, distributed linear COO
- no extra collectives added by the probe inside PETSc callbacks

Only `dt_save` changes:

| Run | `dt_save` | Expected output storage |
|---|---:|---|
| `sparse` | `2*t_final` | Initial/final state only |
| `dense` | `dt` | Approximately every CN step |

All original `DPHS.monitor()` logic still runs. The probe merely
records local memory **after** selected callbacks every 10 accepted
steps, plus before model build, before solve and after the solve.

The analyzer requires identical source commit, DoFs, time-step
parameters, successful 250-step convergence, full rank participation,
and elementwise full initial/final state agreement between the two
runs. It fails rather than reporting an invalid physical comparison.

## Available measurements and limitations

Each rank independently writes `memory_<mode>_rank<i>.json` and two
full solution vectors in `memory_<mode>_rank<i>_state.npz`.

- Linux `/proc/self/status`: `VmRSS`, `VmHWM`, `RssAnon`,
  `RssFile`, `VmSize`, `VmData`, `VmSwap`, plus
  `resource.ru_maxrss`. These include **all** process allocations,
  including PETSc, MUMPS, GetFEM, Python, NumPy and shared mappings.
- Number of saved states and `numpy.ndarray.nbytes` (payload),
  distinguished from whole-process RSS. Neither is a measure of
  native PETSc/MUMPS RAM alone.
- Local PETSc `Mat.getInfo()` allocation metadata for the persistent
  matrices, collected after the solve. It does **not** include all
  MUMPS LU factorization workspace and is not a total process balance.
- The existing SCRIMP `TimingProfiler` for relevant costs during
  setup, initialization and main integration. Those times are
  inclusive/nested.

The analyzer produces:

- `acoustic_memory_rss.png`: maximum per-rank RSS versus CN step,
  dense and sparse saved-output policies
- `acoustic_memory_saved_payload.png`: saved NumPy payload size per
  rank versus CN step, for both policies
- `acoustic_memory_summary.json`: per-rank high-water marks, stored
  payload, PETSc matrix-local information, numerical parity
- `acoustic_memory_rank_summary.csv`: per-rank peak and ending RSS
- `acoustic_memory_timeseries.csv`: individual sampled time steps

Results live in a new git-commit-keyed folder under
`sandbox/outputs/acoustic-memory-h0p04-dt0p002-T0p5-2r-<sha>/`.
The legacy scaling benchmarks are not overwritten.

Use `SCRIMP_MEMORY_TIMEOUT_S=1800` for a longer explicit per-run
timeout. This prevents indefinite computation but is **not** a RAM
limit.

## First native result: 250-step baseline, archived 2026-11-17

**Original benchmark source commit:** `c6b8a06e`.
The user supplied the full run archive, including JSON and full initial/final
NPZ arrays from both ranks and both saving policies. Both solves accepted
**250 CN steps** without rejections or SNES failures, reaching `t_f=0.5`.

| Metric per rank | Sparse save | Dense save |
|---|---:|---:|
| Saved states | 2 | 250 |
| NumPy states payload | 0.732 MiB | 91.442 MiB |
| Max rank RSS at the final TS.monitor | 521.309 MiB | 610.172 MiB |
| Max rank resident-memory high-water mark (full solve) | 521.496 MiB | 792.992 MiB |
| Max rank RSS after full solve | 521.496 MiB | 701.867 MiB |

The **extra ~90.7 MiB** of retained samples corresponds closely
to the dense-minus-sparse RSS gap at the end of time integration.
But after the integration, the dense run also reaches a transient
HWM **~182.8 MiB above** its final monitor RSS, close to **twice**
the 91.4-MiB trajectory payload. Inspection of
`_select_saved_solution_through_final_time()` reveals a full
`np.asarray(states)` stack followed by *two* advanced-index copies
(sort and final-time filter). The hypothesis is therefore **avoidable
full-trajectory temporaries during final consolidation**, not an
established MUMPS memory leak.

Numerical comparison from the archived full NPZ coefficient vectors:

- Initial state is **bitwise equal** in both saving policies and ranks.
- Final state is **bitwise equal across ranks**, separately for each
  saving policy.
- Dense vs sparse final state: max absolute error **3.182885e-10**,
  relative L2 **3.178072e-11**, relative scaled Linf **2.604048e-11**.
- Both logs show exactly 250 CN steps and 250 linear iterations.
- The dense final snapshot was saved within the TS monitor, while the
  sparse endpoint was captured by the **post-TS forced monitor** at
  step `-1`. These different capture paths may account for tiny
  floating-point differences, but no causality has yet been proved.
- Separately, `set_time_scheme()` defaults `ts_adapt_dt_max` from
  `dt_save`, an unfortunate coupling that did not visibly affect
  this fixed-step `ts_adapt_type=none` experiment.

The initial analyzer's componentwise `np.allclose(rtol=1e-10,
atol=1e-11)` failed on a few components near zero despite the global
small errors. The revised analyzer requires **both** global relative
L2 and scaled Linf to be no more than `1e-10`; it records the exact
max absolute and relative differences. This is **numerical agreement**,
not bitwise equality and not proof of complete trajectory parity.

## First targeted memory improvement

Commit `0af04d7` rewrites the final-time trajectory selector in
`scrimp/dphs.py` to **select timestamps/indices first and stack
the retained state vectors only once**, rather than stacking and
making multiple full-array filtered copies. The same timestamps,
ordering, ULP correction and stored values are preserved by targeted
regression tests; it does not change GetFEM/PETSc, the time scheme,
MUMPS or the physical model.

The **prediction** is a lower post-solve peak with the dense storage
policy, possibly on the order of one or two 91.4-MiB full-array
temporaries. A quantitative saving claim requires a second native
Docker run at the same `h`, `dt`, `T` and MPI rank count.
Compare the original `c6b8a06e` directory with the updated output
directory; do not mix raw reports from distinct source commits within
one per-run analysis.

## Decision process after collecting the traces

1. If the dense/sparse RSS gap tracks retained saved-array bytes
   while the sparse run plateaus, focus on output storage policy and
   in-memory trajectory buffering. Preserve any required Hamiltonian,
   power, export and postprocessing APIs.
2. If the sparse run still grows with steps, investigate PETSc/MUMPS
   allocations, persistent Jacobians and reused vectors, with
   targeted measurements and a longer horizon.
3. If both runs have large setup/factorization peaks but little
   longitudinal growth, optimize **peak live memory** and repeated
   matrix allocations; do not call the behavior a leak.
4. Before each production change, compare the complete numerical
   trajectory, power/Hamiltonian and volume balances on appropriate
   physical models, and check time/memory tradeoffs.

We will use the same approach on dissipative 2D Saint-Venant spin-up
after the first acoustic characterization.
