# SCRIMP: asynchronous distributed binary trajectories

Status: **standalone prototype**, NOT wired to `DPHS.solve()`.

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

## Necessary gates before changing production SCRIMP

1. Pass the native fast smoke test above. Inspect the writer log,
   file headers, all rank coverage, and exact full-vector readback.
2. Add an explicit **opt-in** backend, e.g.
   `SCRIMP_TRAJECTORY_BACKEND=async_binary`, default `memory`
   until validated. Do **not** change `SCRIMP_GETFEM_BRIDGE_DISTRIBUTED`.
3. Integrate at the **save decision in `DPHS.monitor`** only:
   copy `z.getArray(readonly=True)` into the bounded producer queue
   when the existing save condition is true. Preserve its initial,
   stop and final-lumped snapshot semantics, time updates and
   `gf_model.next_iter()`. In nonlinear cases retain the
   `_sync_lumped_trial_state` and nonlinear scatter paths.
4. A streaming run keeps `solution["t"]` and small metadata
   but never appends global vectors to `solution["z"]`.
   After TS completes, close transports, verify sidecar
   completion/manifest and substitute a lazy `solution["z"]`
   object. Only then mark `solve_done=True`.
5. Post-processing must read individual full snapshots lazily
   (Hamiltonian, powers and half-step CN interpolation) and match
   an equivalent memory-backend run. A failed I/O writer must not
   produce a solved-success marker or incomplete but plausible
   trajectory.
6. Validate on a **short** fixed-step linear acoustic case, then
   short nonlinear dissipative Saint-Venant 2D with spin-up, including
   exact/near-roundoff state parity and power, Hamiltonian, mass
   balances. No long benchmark is required to validate correctness.
7. Measure bounded transport queue, producer-copy cost, writer
   throughput, end-to-end peak RSS, and I/O error behavior.

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
