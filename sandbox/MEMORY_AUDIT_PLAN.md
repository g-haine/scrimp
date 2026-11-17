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
