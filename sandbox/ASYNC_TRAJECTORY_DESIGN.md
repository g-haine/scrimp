# SCRIMP: asynchronous distributed binary trajectories

Status: **opt-in production integration implemented; short native PDE
validation still required**. The default remains `memory`.

Date: 2026-11-17. Branch: `feature/getfem55-mpi-container-scalability`.

## Goal and non-goals

Replace the O(P × global_DoFs × saved_snapshots) RAM requirement of
the current `self.solution["z"]` storage. A separate I/O process
writes binary data while physical PETSc TS integration proceeds.
Do **not** modify PETSc's solver memory, BDF/CN step history,
SNES/Newton iteration state, Hessians, mass/tangent matrices, GetFEM
nonlinear trial state or physical port-Hamiltonian equations.

The first experiment is a four-snapshot PETSc MPI **transport smoke
test**, not a long PDE run.

## Important distinction: solver state vs. archived state

In `scrimp/dphs.py`:

- PETSc TS manages its current `Vec`, internal time-scheme history
  and nonlinear solver workspaces; we must leave these untouched.
- `IFunction` and `IJacobian` call `_sync_lumped_trial_state(t,z)`.
  When the problem is nonlinear or contains lumped variables, GetFEM
  needs the **current Newton trial state**, not merely a previous
  saved snapshot. It uses a transient `VecScatter.toAll()` and
  `gf_model.to_variables()`. **This path must remain active** until a
  separate GetFEM-specific optimization is validated.
- `DPHS.monitor()` currently gathers a complete global Vec to
  `GlobalVec` even in **linear** problems, and saves a full NumPy
  copy to `self.solution["z"]` **on each MPI rank**. A streaming
  mode can use `z.getArray(readonly=True)` (its locally owned
  entries), avoiding a global gather *for output*.
- For nonlinear systems, transient global scatter may still be
  required to update GetFEM; that is **not** a reason to retain a
  complete trajectory in memory.
- `_select_saved_solution_through_final_time` must be adapted in
  streaming mode to filter the **time/index metadata**, without
  materializing the dense trajectory.
- `Hamiltonian.compute`, `Port.compute`, `get_solution`,
  `get_quantity` and `export_to_pv` index
  `self.solution["z"][k]`. A lazy sequence can reconstruct **one
  global state at a time** from partitioned binaries; existing API
  behavior must be tested for all operations, including midpoint
  states used in Crank–Nicolson power calculations.
- Any full-trajectory materialization via `np.asarray` must be
  explicitly opt-in and disclosed, not an automatic part of solve.

## Transport prototype (already committed)

`sandbox/async_binary_trajectory.py` provides:

1. **One completely separate Python I/O process**, started by the
   **container shell BEFORE `mpiexec`**. It imports no MPI, GetFEM or
   PETSc, and one I/O thread handles each rank's socket connection.
2. In each MPI rank, a background **transport thread** is the only
   component that performs blocking socket writes. A PETSc callback
   only makes one local NumPy copy and calls `queue.put_nowait`.
3. Only rank-owned values are transmitted, with global ownership
   `[lo, hi)` and an explicit scalar type (float64 or complex128).
   **There is no MPI all-gather for output.**
4. Per-rank binary files contain a fixed header and records
   `(step:int64, time:float64, owned_values[:])`. For the first
   local/single-node prototype, all ranks access one shared path.
5. A writer-generated **atomic JSON manifest** is published only
   after every rank finishes, global row partitions exactly cover
   `[0, N)`, all ranks report matching steps/timestamps, and scalar
   types/global sizes agree. Incomplete/failed runs never get a
   completed manifest.
6. `LazyTrajectory(manifest)[k]` memory-maps all local rank files
   and assembles just **one** global state (not the full time series).

The queue is bounded by `max_queued` (default 4) **per rank**. A
full queue raises `AsyncOutputBackpressure` immediately, rather
than silently dropping a state or growing memory without limit.
Writer/process errors propagate on closing/joining outside TS.
Thus callback-level disk waiting is avoided in normal operation.

**Physical limitation:** exact, lossless, indefinite streaming with
bounded RAM **cannot** guarantee that compute never waits or fails
when the disk remains slower than data production. The prototype
chooses **explicit failure** on saturation. An optional policy
permitting bounded waiting could be added after characterization;
automatic dropped snapshots are unacceptable for scientific output.

## Native fast smoke test (no PDE simulation)

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash -n sandbox/run_async_binary_smoke.sh
bash sandbox/run_async_binary_smoke.sh
```

The Docker runner first runs isolated subprocess tests for:

- exact rank-partitioned reconstruction (real and complex);
- correct step/time metadata and lazy indexing;
- shape/dtype enforcement and non-finite-state rejection;
- explicit queue saturation rather than data loss;
- rejection of mismatched rank timestamps.

It then launches the writer from the shell, followed by **2 native
MPI ranks** with one PETSc Vec of 17 global DOFs (uneven partition).
Each rank streams four small local snapshots; the readback must
reconstruct **all 17 coefficients exactly at all four times**.
No GetFEM/PDE integration or physical long run is involved.

Results are saved under
`sandbox/outputs/async-binary-smoke-<source-commit>/`.

## Integration status, safety gates and next short native check

1. Pass the native fast smoke test above. Inspect the writer log,
   file headers, all rank coverage, and exact full-vector readback.
2. **Implemented (opt-in):** `SCRIMP_TRAJECTORY_BACKEND=async_binary`
   (default remains `memory`). It requires a writer started **before**
   `mpiexec` and `SCRIMP_TRAJECTORY_DIR` and
   `SCRIMP_TRAJECTORY_SOCKET`. Optional bounded configuration:
   `SCRIMP_TRAJECTORY_MAX_QUEUED=4` and
   `SCRIMP_TRAJECTORY_COMMIT_TIMEOUT_S=30`.
   PETSc's matrix bridge setting is unchanged.
3. **Implemented:** at `DPHS.monitor()`'s original save decision,
   send just `z.getArray(readonly=True).real` with the PETSc
   ownership interval. The historical trajectory used real float64,
   so the opt-in backend deliberately preserves its output dtype.
   No full global gather occurs **for linear output**.
4. **Preserved:** nonlinear GetFEM trial-state synchronization
   (`_sync_lumped_trial_state`) and main monitor's nonlinear
   global scatter + `to_variables`. Initial-consistency pseudo TS
   still uses its original monitor behavior. PETSc TS/SNES/KSP/MUMPS
   memory and time-stepping history are untouched.
5. **Implemented:** `solution["t"]` accumulates timestamps only;
   `solution["z"]` remains an empty list during the integration.
   On successful completion, all sender threads drain, the external
   writer must atomically publish a verified manifest, and
   `solution["z"]` becomes a `LazyTrajectory` object. Index
   filtering preserves SCRIMP's ULP final-time logic.
   A writer timeout/failure prevents `solve_done=True`.
6. **Next gate:** validate the **short** fixed-step linear acoustic
   case (below), comparing every saved coefficient, Hamiltonian
   and port powers against the default memory mode. Then
   validate short nonlinear dissipative Saint-Venant 2D with
   spin-up, including power and mass balances. **No long run required.**
7. **Later:** measure bounded queue behavior, producer-copy cost,
   writer throughput and process peak RSS, and complete failure-path
   validation before changing the default.

## Short physical parity validation: acoustic PASS (2026-11-17)

The tiny PETSc native smoke test (17 DoFs, 4 frames, 2 ranks) has
**passed** on the user's Docker image (`29fbf9fd`); the user supplied
the binary archive. Manifest: `complex128` PETSc data, full
ownership `[0,9)` and `[9,17)`, exact reconstructed states.

Production integration is now in `scrimp/utils/async_binary_trajectory.py`
and `scrimp/dphs.py`. The sandbox launcher preserves the old
`sandbox/async_binary_trajectory.py` CLI via a thin proxy.

The user **completed** the six-step native physical test on source
commit `44b487be`, and supplied
`async-acoustic-parity-44b487be-he1O1w.zip`.
The archive was independently inspected:

- **12,682 DoFs, 2 MPI ranks, 6 accepted CN steps**, 7 saved
  states, no rejections or nonlinear failures.
- The two binary files carry exactly global ownership
  `[0, 6341)` and `[6341, 12682)`, 7 frames each.
- Direct decoding of raw binary frames reconstructs the streamed
  `async_binary_physical.npz` trajectory **bit for bit**.
- Both saved-time arrays are exactly equal, and full state arrays
  agree to **1.665334537e−16 max absolute error**,
  **1.524407066e−17 relative L2** between independent physical
  simulations (numerical roundoff only).
- The Hamiltonian and **all four boundary-power trajectories
  match exactly** between memory and asynchronous backends.
- The archive manifest and process-level frame checks validate
  that no global NumPy trajectory was written in the producer ranks.
  The test driver materialized the very short full trajectory only
  for parity comparison.

Reproduction command for the PASSED acoustic experiment:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash -n sandbox/run_async_acoustic_parity.sh
bash sandbox/run_async_acoustic_parity.sh
```

This starts one external writer before PETSc MPI, compares both modes
on the same `examples/wave.py` model (`h=0.08`, `dt=0.001`,
`T=0.006`, 2 ranks, COO/MUMPS `never`), checks the complete
trajectory, Hamiltonian and computed port powers, and requires a
committed `trajectory_manifest.json`. The short driver
**materializes the tiny trajectory only for validation**, not in
production. It writes a source-SHA-segregated report under
`sandbox/outputs/async-acoustic-parity-<sha>-<unique>/`.
The GitHub CI checks source syntax and the Python writer protocol but
**cannot** establish native GetFEM/PETSc physics parity on its own.

## Next: three-step nonlinear dissipative Saint-Venant spin-up

The next **short** native comparison is prepared, with no long run:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash -n sandbox/run_async_swe_parity.sh
bash sandbox/run_async_swe_parity.sh
```

This runs **two** physically identical very short (three-step)
tests using `sandbox/viscous_swe_2d_mpi.py:build_system` and its
unmodified `diagnostics()`:

- Open 2D viscous Saint-Venant, **nonlinear mass and stiffness**,
  `h_mesh=0.35`, two MPI ranks, Crank–Nicolson fixed
  `dt=0.001`, `t_f=0.003`, three steps and four archived states.
- Initial rigid profile at **rest**, with a smooth, nonzero
  facet-tangent wall-velocity spin-up over `0.03` time units;
  `gamma=0`, physical `facet_tangent` wall.
- Existing BDF2 initial consistency step and nonlinear GetFEM
  synchronization are preserved; `PC=LU/MUMPS`, reset
  policy `on_change`.
- Compare complete state snapshots, mechanical Hamiltonian,
  kinetic/potential energy, both external boundary-power
  contributions, both viscous dissipation terms, integrated
  input/dissipation, full energy and volume balance histories,
  plus Lyapunov samples, with explicit numeric gates.
- Reject incomplete manifests, mismatched simulation
  metadata, nonconvergence, nonfinite fields and negative
  water depth. The script reports balance defects; a PASS means
  **storage-backend equivalence**, not independent proof of
  published Saint-Venant physical accuracy.

This is a **new test definition**. No native Saint-Venant
results are claimed until the Docker script runs successfully.

The non-streaming default and existing model behavior are preserved.
Do not enable this backend automatically for all SCRIMP users until
the native physics tests, nonlinear Jacobian tests, writer failure
tests and physical power balances have passed.

### MPI/multinode caveat

The first writer uses Unix-domain sockets and a **shared local
filesystem**, suitable for the user's single-host Docker MPI setup.
On a multinode cluster, introduce **one writer per node** or MPI I/O
ranks, and use an appropriate shared-output layout. Do not silently
assume Unix sockets work across nodes.

### Why not fork a writer from PETSc?

Forking an initialized MPI/PETSc process can corrupt communicator
state or deadlock. The smoke-test shell launches a fresh writer
interpreter before `mpiexec`. Production launch and ownership must
be designed to preserve this separation.
