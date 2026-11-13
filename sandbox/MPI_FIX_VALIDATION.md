# MPI repair validation — issue #6

Worktree: `fix/issue-6-parallel-powers-saint-venant-2d`.

This is a **work-in-progress** correction for GetFEM MPI level 2 and
PETSc parallel layouts, not yet a validated numerical release.

## Reason for the change

Observed with `mpiexec -n 2`: last visible messages were

```text
Starting linear mass matrix assembly...
Linear mass matrix assembly done ...
Starting linear stiffness matrix assembly...
Assembly time ...
```

The failure is **before TS.solve**; the exact blocking call was not
identified from those logs alone. The instrumentation below distinguishes
GetFEM assembly, CSR transfer and PETSc Mat.assemble on every rank.

The user confirmed `gf.util_mpi_parallelism_level() == 2` for both
1-rank and 2-rank runs. GetFEM's documentation states that at level 2
model matrices carry **distributed contributions**, while assembled
model RHS vectors are **globally summed and replicated**:

https://getfem.org/userdoc/parallel.html

## Changes in this branch

- `scrimp/utils/linalg.py`: use global-row PETSc IJV insertion with
  ADD_VALUES for the *distributed* GetFEM matrix contributions; clear
  the target before each nonlinear reassembly.
- `scrimp/dphs.py`: insert globally reduced GetFEM RHS only on
  PETSc-owned rows; fix local/global initial vector indices in both
  TS initialization paths; allocate tolerances using the same
  distributed PETSc layout, not a full array per rank.
- `scrimp/dphs.py`: synchronize GetFEM time and the current Newton
  iterate for nonlinear spatially distributed PDEs, not just lumped
  systems. The historical state-dependent mass-matrix residual still
  requires a separate mathematical validation.
- `tests/test_mpi_algebra.py`: real PETSc MPI tests of the algebra
  bridge with synthetic distributed GetFEM contributions, RHS and
  tolerances. These do not replace physical GetFEM simulations.

## Local test sequence — run from repository root

Make sure the working tree is clean before switching branches.

```bash
git fetch origin
git switch fix/issue-6-parallel-powers-saint-venant-2d
git pull --ff-only

# Fast syntax checks; do not require the numerical solver.
python -m py_compile scrimp/utils/linalg.py scrimp/dphs.py \
    tests/test_parallel_powers.py tests/test_mpi_algebra.py

# Actual MPI tests (require getfem, PETSc, SLEPc and mpi4py).
# Use a timeout so a deadlock cannot run indefinitely.
timeout -k 5s 45s mpiexec -n 1 \
    python -u -m unittest -v tests.test_mpi_algebra
timeout -k 5s 45s mpiexec -n 2 \
    python -u -m unittest -v tests.test_mpi_algebra
timeout -k 5s 45s mpiexec -n 4 \
    python -u -m unittest -v tests.test_mpi_algebra

python -m unittest -v tests.test_parallel_powers
```

Do not proceed if algebra bridge tests fail. A failure means the
matrix transfer or vector layout needs adjustment before launching
costly Saint-Venant simulations.

## Diagnostic replay of original hang

With an MPI-safe parallel LU such as MUMPS, use the smallest smoke test:

```bash
SCRIMP_MPI_TRACE=1 timeout -k 5s 120s mpiexec -n 2 \
    python -u sandbox/viscous_swe_2d_mpi.py \
      --case open --mesh 0.2 --final-time 0.5 \
      --save-dt 0.1 --max-step 0.01 \
      > mpi_open_2ranks.log 2>&1

echo "exit code: $?"
tail -n 100 mpi_open_2ranks.log
```

`SCRIMP_MPI_TRACE=1` prints to standard error from **every rank**
before/after GetFEM build_matrix, before/after conversion, and
before/after Mat.assemble. It introduces **no additional MPI barrier**.

If the test fails, collect the complete file including each rank's last
stage, the command, GetFEM/PETSc versions and MPI rank count.

Only after the 2-rank run completes, repeat with 1 and 4 ranks,
then compare computed operator norms, final states, Hamiltonian and
power balance. For a credible release, the same discretization must
produce consistent results within numerical tolerances independently
of MPI size.

## 2026-11-13 — Real GetFEM test exposed MPI double-counting

The first user run of `test_real_getfem_mpi2_mass_matrix_is_rank_independent`
reported the **exact** real physical result on one rank
(`1^T M 1 = 1`) but **twice** the expected value with two ranks
(`1^T M 1 = 2`). The test never reached four ranks because the
loop stops on the first failed rank count.

This invalidates the prior assumption that every matrix exposed by
`gf.Model.tangent_matrix()` through Python necessarily carries
different rank-local partial contributions. The real mass test and the
equal coefficient counts on rank 0/1 are compatible with a **replicated
matrix** for the tested brick. GetFEM's general MPI-2 documentation
still describes partial *internal* matrix contributions; the Python
interface / specific brick may expose different semantics. We must
not conflate those contracts.

The converter now explicitly distinguishes two cases:

- **Replicated**: every rank holds an exactly identical matrix. PETSc
  inserts only rows locally owned by that rank with INSERT_VALUES.
- **Distributed**: ranks have different contributions; PETSc sums
  global-index entries with ADD_VALUES.

Historical (now superseded): the old `SCRIMP_GETFEM_MATRIX_MODE=auto`
selected between them by comparing exact CSR signatures. This approach
was subsequently found to **flap between modes during Newton iterations**
and has been removed. The default is now **replicated**. For a known
partial CSR matrix use `SCRIMP_GETFEM_MATRIX_MODE=distributed` explicitly. We have deliberately not divided entries by
the MPI size.

The synthetic matrix tests were also adjusted to exercise genuinely
distinct partial matrices (including a shared row) and replicated
matrices. Five tests were later validated at 1/2/4 ranks; a **sixth**
was added afterward for distinct PETSc operator/preconditioner matrices;
a **seventh** now exercises real nonlinear GetFEM tangents at two Newton
states. Both new tests require a real MPI run.

**Run this first; do not rerun Saint-Venant if this test fails:**

```bash
git fetch origin
git switch fix/issue-6-parallel-powers-saint-venant-2d
git pull --ff-only
for n in 1 2 4; do
  echo "MPI ranks = $n"
  timeout -k 5s 45s mpiexec -n "$n" \
    python -u -m unittest -v tests.test_mpi_algebra || break
done
```

## Observed SIGSEGV after synthetic tests passed (2026-11-13)

The user confirmed the **three original synthetic tests pass** with 1, 2
and 4 MPI ranks. A subsequent 2-rank Saint-Venant smoke still exits
**139 / SIGSEGV**, but it now completes multiple nonlinear GetFEM/PETSc
assemblies. The final printed trace was:

```text
[SCRIMP MPI 0/2] nl_stiffness: PETSc matrix assembled
[SCRIMP MPI 1/2] nl_stiffness: PETSc matrix assembled
```

This does **not** prove the assembly is mathematically correct or locate
the offending C/C++/Fortran call. A segfault is **not** a timeout.

The next patch adds traces around `tangent_stiffness.copy()`,
`Mat.axpy()`, `disable_bricks()`, MPI barrier, and the TS residual and
Jacobian callbacks. `linalg.convert_gmm_to_petsc` also prints the
number of genuinely nonzero sparse values and absolute-value sum,
not just the number of CSR slots. Equal CSR slot counts on both ranks
**do not** establish that the same matrix was assembled twice.

A **fourth test** `test_real_getfem_mpi2_mass_matrix_is_rank_independent`
uses an actual GetFEM mass brick and checks that the PETSc global
quadratic form for a constant function on the unit interval is exactly 1.
This is the first test of the real GetFEM model assembly; it has **not
yet run** in the target environment.

Test it **before** rerunning Saint-Venant:

```bash
git fetch origin
git switch fix/issue-6-parallel-powers-saint-venant-2d
git pull --ff-only

for n in 1 2 4; do
    echo "====== MPI ranks: $n ======"
    timeout -k 5s 45s mpiexec -n "$n" \
        python -u -m unittest -v tests.test_mpi_algebra || break
done
```

If all four tests pass, capture **a Python-level fault traceback** alongside
per-rank assembly traces (this does not replace a GDB native backtrace):

```bash
PYTHONFAULTHANDLER=1 SCRIMP_MPI_TRACE=1 \
timeout -k 5s 120s mpiexec -n 2 python -u \
    sandbox/viscous_swe_2d_mpi.py \
    --case open --mesh 0.2 --final-time 0.5 \
    --save-dt 0.1 --max-step 0.01 \
    > mpi_open_2ranks.log 2>&1

printf "exit code: %s\\n" "$?"
tail -n 160 mpi_open_2ranks.log
```

Interpret the new `solve:/init_step:` rank traces first: the
driver uses an **initial pseudo-TS consistency step**, before the main
ARKIMEX time integrator. If the crash is in the pseudo-TS initialization,
isolate this by running the same small problem on **one rank**, and
then compare with the diagnostic-only option `--no-init-step`:

```bash
# Same configuration, but only one MPI rank:
PYTHONFAULTHANDLER=1 SCRIMP_MPI_TRACE=1 \
timeout -k 5s 120s mpiexec -n 1 python -u \
    sandbox/viscous_swe_2d_mpi.py \
    --case open --mesh 0.2 --final-time 0.5 \
    --save-dt 0.1 --max-step 0.01 \
    > mpi_open_1rank.log 2>&1

# Diagnostic only: omission of consistency solve changes the DAE
# initialization. Results are NOT physical validation.
PYTHONFAULTHANDLER=1 SCRIMP_MPI_TRACE=1 \
timeout -k 5s 120s mpiexec -n 2 python -u \
    sandbox/viscous_swe_2d_mpi.py \
    --case open --mesh 0.2 --final-time 0.5 \
    --save-dt 0.1 --max-step 0.01 --no-init-step \
    > mpi_open_no_init_2ranks.log 2>&1
```

If faulthandler does not identify the native failing call, a GDB backtrace
from the process that crashed may be needed; do not try to fix a SIGSEGV by
blindly changing PETSc solver options or by inserting extra barriers.

## Pending work

- Confirm that the full GetFEM MPI-2 matrix CSR is correctly reconstructed
  under real model assembly, not only the synthetic distributed-matrix
  test; compare 1/2/4-rank operator-vector products.
- Investigate the nonlinear, state-dependent mass residual and Jacobian,
  including trial-state derivatives. No unsupported correctness claim.
- Eliminate all-rank ParaView writes and fixed-name Gmsh mesh races,
  after checking whether rank-0 GetFEM MPI-2 VTU contains a full domain.
- Check collective reductions for scalar GetFEM power/energy integrals,
  avoid redundant global gathers and consider MPI benchmarking.
- Run both published viscous Saint-Venant 2D cases, validate energy
  balances, Lyapunov evolution and serial/MPI agreement.
- Merge, tag or publish **only after** those gates are met.

## 2026-11-13 — Crash occurs AFTER IJacobian; matrix-mode flapping detected

The latest 2-rank smoke run still exited with `SIGSEGV` (139).
Both ranks printed `IJacobian: completed` immediately before termination.
The native crash is therefore **after** the Python Jacobian callback
returned; no Python exception was reported. This points to the native
PETSc TS/SNES/KSP/PC path, but a GDB backtrace is required to name the
faulting native call.

**Important additional observation:** repeated 3300-DOF nonlinear
stiffness assemblies toggled the previous SHA-based automatic mode
between `replicated` and `distributed`, even though both ranks
reported equal CSR counts and equal absolute-value sums. This can
spuriously double the matrix in some Newton iterations and destroy
numerical consistency. It is **not** proof of a native PETSc bug.

The converter no longer uses floating-point CSR hashes to decide
the insertion semantics. **Default: `replicated`**, which is supported
by the 1/2/4-rank real GetFEM mass integral tests. For an explicitly
known partial assembly, select `SCRIMP_GETFEM_MATRIX_MODE=distributed`.
The two synthetic matrix tests explicitly force the appropriate mode.

The `IJacobian` callback also now fills a distinct operator matrix
`A` from `P` rather than merely assembling an unchanged `A`.
With `SCRIMP_MPI_TRACE=1`, it reports PETSc-owned row dimensions,
nonzero counts and non-finite/zero diagonal diagnostics just before
the native factorization.

**Required checks after pulling:**

1. Re-run the **seven** algebra tests on 1/2/4 ranks.
2. Re-run the original 2-rank open-loop smoke with
   `PYTHONFAULTHANDLER=1 SCRIMP_MPI_TRACE=1`; confirm that every
   matrix transfer logs `mode=replicated` and inspect the last
   `IJacobian: operator ...` diagnostic.
3. If it still exits 139, compare the same case on **one rank**,
   and use the driver's `--no-init-step` switch on two ranks
   to isolate the pseudo-TS consistency initializer. This switch
   is *diagnostic only*: do not claim correct DAE initialization.

A native stack trace, if needed, can be obtained from the same
command launched through GDB on each rank (GDB must be installed):

```bash
SCRIMP_MPI_TRACE=1 mpiexec -n 2 gdb -q -batch \
    -ex run -ex 'thread apply all bt 20' --args \
    python -u sandbox/viscous_swe_2d_mpi.py \
    --case open --mesh 0.2 --final-time 0.5 \
    --save-dt 0.1 --max-step 0.01 \
    > mpi_gdb.log 2>&1
```

Use this only if the simplified tests still fail: GDB adds overhead
and running multiple debuggers in MPI can require local ptrace access.


## Seven-test gate (after fixed-mode transfer)

A seventh test now exercises an **actual nonlinear GetFEM model**:
`R(u)=integral(u^2 * Test_u)` on [0,1], with constant test states
u=1 and u=2. After every GetFEM assembly and conversion into the
**same** PETSc matrix, its tangent must satisfy
`1^T J(u) 1 = 2u` at 1, 2 and 4 ranks. This targets the nonlinear
stiffness mode flapping observed in the latest Saint-Venant log.

Both the new nonlinear test and the separate A/P Jacobian test remain
**unvalidated in the user's environment**. Syntax CI alone is not
sufficient to mark this gate complete.

## 2026-11-13 — Seven algebra tests pass; pseudo-TS aborts with allocator corruption

User supplied complete results for the revised suite:
- 1 rank: **7/7 OK**.
- 2 ranks: **7/7 OK** on *both* processes.
- 4 ranks: **7/7 OK** on *all four* processes.
- The real GetFEM MPI-2 mass test, the nonlinear tangent reassembly test,
  and the distinct PETSc A/P test all pass.

The subsequent open-loop Saint-Venant smoke exited **134 / SIGABRT**:
```text
IJacobian: completed                  # on BOTH ranks
malloc(): largebin double linked list corrupted (nextsize)
Fatal Python error: Aborted
  ... dphs.py, in init_step
```
At the last Jacobian callback, `P` was 3300×3300, with **zero
nonfinite matrix entries** on either rank; zero *diagonal* entries:
**132 on rank 0**, **0 on rank 1**. These diagonal zeros do **not**
prove singularity for a mixed saddle-point system. We must determine
whether any **entire numerical rows** are zero and identify their
variables. The fault is detected by the native memory allocator
**after** `IJacobian` returned, inside PETSc's pseudo-TS stage. The
native origin of corruption remains unknown.

### Additional targeted fix and diagnostic

The code review found `DPHS.init_matrix` initializing a PETSc vector
through `x.setArray(0)` — `Vec.setArray` expects a local-value
array, not a scalar. It now uses `x.set(0.0)`, creates the left
distributed vector (matching matrix rows), and explicitly destroys it
after setting/assembling the zero diagonal. This is a **plausible**
memory-corruption cause, **not yet demonstrated** as the culprit.

A new regression `test_init_matrix_preserves_distributed_vector_storage`
initializes matrices repeatedly and checks that they remain valid on
the relevant MPI layout. The MPI suite now has **8 tests**; the new
one has not yet been run.

With `SCRIMP_MPI_TRACE=1`, `IJacobian` additionally prints counts
and global index samples of **whole zero rows**, alongside zero
diagonal counts, to distinguish missing equations from expected
saddle-point structure.

**Next run (do not skip the tests):**

```bash
git pull --ff-only

for n in 1 2 4; do
  timeout -k 5s 45s mpiexec -n "$n" \
    python -u -m unittest -v tests.test_mpi_algebra || break
done
```

If all **8** tests pass:

```bash
PYTHONFAULTHANDLER=1 SCRIMP_MPI_TRACE=1 \
timeout -k 5s 120s mpiexec -n 2 \
  python -u sandbox/viscous_swe_2d_mpi.py \
    --case open --mesh 0.2 --final-time 0.5 \
    --save-dt 0.1 --max-step 0.01 \
  > mpi_open_2ranks.log 2>&1

echo "exit=$?"
tail -n 100 mpi_open_2ranks.log
```

If still exit 134/139, obtain a **native** GDB backtrace on the
crashing rank (see earlier instructions in this file). A Python
faulthandler frame at `TS.solve` does not identify the native
offending routine. Then run the **one-rank** smoke, followed by a
two-rank `--no-init-step` diagnostic run if needed; do not interpret
that latter run as physically valid.


## 2026-11-13 — Native GDB identifies a forced PETSc Mat destroy during TS monitor

The user provided `mpi_open_1rank.log` and `mpi_gdb.log`.

### One-rank smoke: clean process exit, **invalid numerical result**

```text
exit code = 0
pseudo-TS initial consistency step completed
main TS: 1 accepted step, 6 rejected
stopped at t = 1e-8 (requested t_f = 0.5)
nonlinear-mass abs_sum ~ 1.10e14
nonlinear-stiffness abs_sum ~ 5.26e38
max relative energy/power-balance defect = 48602.5983353432
```

Therefore *no* physical validation is claimed, even in serial. In
particular, the extreme Jacobian coefficients and automatic minimum-dt
stop are numerically unacceptable. The current experimental driver
now raises an explicit `RuntimeError` when it stops before `t_f`,
rather than exiting 0 and producing misleading metrics.

### Two-rank native stack: **specific destruction path**

On **both** ranks the last completed nonlinear callback was followed
by `free(): invalid pointer`. The GDB stack on the rank whose full
backtrace was captured is:

```text
libc free
 PetscFreeAlign -> PetscFreeA
 MatSeqXAIJFreeAIJ -> MatDestroy_SeqAIJ
 MatDestroy_MPIAIJ -> MatDestroy
 PetscObjectDestroy -> PetscGarbageCleanup
 petsc4py.PETSc.garbage_cleanup
 petsc4py.PETSc_TS_Monitor
 TSMonitor -> TSSolve
```

This identifies the **point of detection**: an MPIAIJ matrix is being
destroyed *inside the TS monitor*, not in a MUMPS factorization.
The invalid pointer may have been corrupted **before** destruction.
It does **not**, by itself, prove a PETSc, GetFEM or MUMPS library bug.
Gmsh had many parked worker threads, but they are not the crashing
stack in this GDB trace.

A code review then found:
- `DPHS.monitor()` manually called `PETSc.garbage_cleanup()` at
  **every** monitor event.
- `DPHS.assemble_nl_mass()` replaced `self.tangent_mass` with a *new*
  PETSc `Mat` via `self.mass.copy()` on **every** evaluation.
- `DPHS.assemble_nl_stiffness()` similarly replaced
  `self.tangent_stiffness` with `self.stiffness.copy()`.
- The old MPIAIJ matrices were thus queued for native destruction at
  the next forced `PetscGarbageCleanup`, inside `TS.solve`.

### Targeted experimental repair committed to dedicated branch

- **Persistent PETSc tangent matrices:** keep the original
  `self.tangent_mass` and `self.tangent_stiffness` Mat handles; zero
  and assemble them *in place* as linear + nonlinear contributions.
  `MatAXPY(..., DIFFERENT_NONZERO_PATTERN)` handles changing sparse
  patterns. Temporarily allow new nonzero allocation, then reinstate
  normal strict allocation checking.
- **No PETSc garbage collection inside `DPHS.monitor`:** remove the
  manually forced cleanup during a live TS callback.
- Avoid disposable `self.GlobalVec.copy()` allocations in the monitor;
  copy only the needed local NumPy values.
- New **ninth** MPI regression
  `test_nonlinear_tangents_reuse_mpi_matrices_without_accumulation`
  checks stable Mat identities, nonlinear coefficient and sparsity
  changes, repeatable matrix actions and no accumulation.
- The experimental SWE driver now explicitly errors if its main time
  integration does not reach the requested `--final-time`.

**These patches are hypotheses to test, NOT confirmed crash fixes.**
The local container has no GetFEM/PETSc/MPICH numerical runtime. The
GitHub source-check CI compiles syntax, not the distributed computation.

### Validation order

```bash
git pull --ff-only
for n in 1 2 4; do
  echo "MPI ranks: $n"
  timeout -k 5s 45s mpiexec -n "$n" \
    python -u -m unittest -v tests.test_mpi_algebra || break
done
```

Expected **9 tests on every MPI rank**; abort before the smoke test
if an algebra test fails.

```bash
PYTHONFAULTHANDLER=1 SCRIMP_MPI_TRACE=1 \
timeout -k 5s 120s mpiexec -n 2 \
  python -u sandbox/viscous_swe_2d_mpi.py \
    --case open --mesh 0.2 --final-time 0.5 \
    --save-dt 0.1 --max-step 0.01 \
  > mpi_open_2ranks.log 2>&1
echo "exit=$?"
tail -n 90 mpi_open_2ranks.log
```

If the native SIGABRT disappears but the integrator stops short or the
Hamiltonian/power defect is large, the **MPI memory problem might be
partially fixed, but the numerical Saint-Venant model is still wrong**.
Do not close #6, merge #22 or tag v1.2.1 until physical balances and
1/2/4-rank agreement are independently validated. The residual
`IFunction` in SCRIMP still uses a historical
`K_tangent(z) @ z + rhs(z)` reconstruction for state-dependent mass,
which is **not generically equal** to the nonlinear residual. This
requires a separate mathematical correction and finite-difference
Jacobian audit.

## 2026-11-13 — Nine tests pass; SIGSEGV persists after persistent matrices

After pulling the persistent-matrix/monitor-GC changes, the user ran
`tests.test_mpi_algebra`: **9/9 passing on every process** for 1,
2 and 4 MPI ranks. In particular,
`test_nonlinear_tangents_reuse_mpi_matrices_without_accumulation`
passes (changing coefficients and sparse patterns).

However the two-rank Saint-Venant smoke still exits with
**139 / SIGSEGV**, after both ranks print `IJacobian: completed`.
The matrix is 3300x3300; no NaNs/Infs and no completely zero rows;
rank 0 has 132 zero diagonal entries (expected possibility for mixed
algebraic equations). The new trace confirms `nl_mass:
tangent reassembled in place` and `nl_stiffness: tangent reassembled
in place` on both ranks, while `mode=replicated` stays stable.

**Therefore, removal of the monitor's forced `PETSc.garbage_cleanup()`
and retention of the tangent matrix handles have NOT eliminated the
crash.** The previous GDB trace identified an invalid `MatDestroy_MPIAIJ`
*at the forced cleanup*, but not necessarily the original corruption.
It is unsafe to infer a PETSc or MUMPS library defect without a
**new native trace of the CURRENT code**.

### Dynamic-pattern copy contract: candidate hardening

In the actual SWE log the nonlinear stiffness GetFEM CSR gains many
entries as Newton evolves (e.g. initial 86k to ~126k stored values).
Yet `self.J` is duplicated from an **initial** tangent matrix,
and `IJacobian()` copied a potentially different sparse structure
into that existing matrix. PETSc 3.21 documentation explicitly says
that `MatCopy` with `SAME_NONZERO_PATTERN` on different patterns
**may crash**:
https://web.cels.anl.gov/projects/petsc/vault/petsc-3.21/docs/manualpages/Mat/MatCopy.html

To make the contract unambiguous, `IJacobian` now specifies
`PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN` when copying the
tangent into P, adding the state-dependent mass with `MatAXPY`, and
copying P into a distinct A. This is a **defensive fix**; the previous
petsc4py default's exact behavior and whether it caused the crash
have not been established.

A new **tenth** MPI test
`test_jacobian_copy_with_changing_sparse_pattern` checks that A/P
initially have diagonal-only structure, then repeatedly accept
off-diagonal entries, including off-process columns.

### Validation commands

```bash
git pull --ff-only

for n in 1 2 4; do
  echo "MPI ranks: $n"
  timeout -k 5s 45s mpiexec -n "$n" \
    python -u -m unittest -v tests.test_mpi_algebra || break
done
```

Require **10/10 passing on all 1, 2 and 4 ranks**. If they pass,
rerun the two-rank SWE smoke:

```bash
PYTHONFAULTHANDLER=1 SCRIMP_MPI_TRACE=1 \
timeout -k 5s 120s mpiexec -n 2 \
  python -u sandbox/viscous_swe_2d_mpi.py \
    --case open --mesh 0.2 --final-time 0.5 \
    --save-dt 0.1 --max-step 0.01 \
  > mpi_open_2ranks.log 2>&1
echo "exit=$?"
tail -n 100 mpi_open_2ranks.log
```

If exit 139/134 persists, obtain a **new** native stack on the
corrected branch (show current HEAD with `git rev-parse --short HEAD`):

```bash
SCRIMP_MPI_TRACE=1 \
timeout -k 5s 120s mpiexec -n 2 \
  gdb -q -batch -ex 'set pagination off' -ex run -ex 'bt 40' \
    --args python -u sandbox/viscous_swe_2d_mpi.py \
    --case open --mesh 0.2 --final-time 0.5 \
    --save-dt 0.1 --max-step 0.01 \
  > mpi_gdb_after_pattern_fix.log 2>&1
echo "exit=$?"
```

For an independent **diagnostic-only** isolation of the initial
pseudo-TS (no claim of valid initial DAE consistency), also try
`--no-init-step` at 2 MPI ranks; if it does not crash but terminates
early with the new `RuntimeError`, that still distinguishes the
native failure stage from the invalid physical result.

## 2026-11-13 — Ten tests green; GDB pinpoints MUMPS factor destruction in TS.reset

**Latest user-verified results:** all **10 algebra tests pass** at MPI
sizes 1, 2 and 4. The original 2-rank SWE smoke still exits **139**
when not running under GDB, after both `IJacobian: completed`.

The user then captured the CURRENT native GDB stack with the hardened
MatCopy/MatAXPY code. **Crucially, the GDB run advances further than
the non-GDB run:** on both ranks:

```text
[SCRIMP MPI 0/2] init_step: after pseudo TS.solve
[SCRIMP MPI 1/2] init_step: after pseudo TS.solve
double free or corruption (!prev)
Thread 1 "python" received signal SIGABRT
```

GDB's native stack (same path on both ranks):

```text
__libc_free
  PetscFreeAlign
  MatDestroy_MUMPS
  MatDestroy
  PCReset_LU
  PCReset
  KSPReset
  SNESReset
  TSReset
  petsc4py.PETSc.TS.reset
```

This identifies **MUMPS factor destruction during explicit
`TS.reset()` AFTER `TS.solve()`** as the point where the allocator
detects a double free / heap corruption. This is NOT a claim that the
MUMPS library itself introduced the corruption. Since GDB affects
memory layout, the earlier non-GDB SIGSEGV after the Jacobian can be a
different manifestation of the same underlying native issue.

Also, the GDB run shows massive growth in GetFEM nonlinear term norms
during the pseudo-TS (nonlinear stiffness absolute sum about
`4.48e46`): the formulation/consistency solve is **numerically
invalid**, even when the native crash is isolated.

### First isolate the PETSc/MUMPS factor lifecycle without SCRIMP

A new standalone script, `sandbox/test_petsc_mumps_reset.py`, imports
**only PETSc**, assembles a well-conditioned distributed tridiagonal
matrix, solves with global LU/MUMPS and explicitly resets KSP
(`PCReset_LU -> MatDestroy_MUMPS`). The option `--solver bjacobi`
runs a control case without global MUMPS. Both run repeated cycles.

```bash
git pull --ff-only

timeout -k 5s 45s mpiexec -n 1 \
  python -u sandbox/test_petsc_mumps_reset.py \
  --solver mumps --cycles 3 \
  > mumps_minimal_1rank.log 2>&1
echo "minimal 1 rank exit=$?"
tail -n 35 mumps_minimal_1rank.log

timeout -k 5s 45s mpiexec -n 2 \
  python -u sandbox/test_petsc_mumps_reset.py \
  --solver mumps --cycles 3 \
  > mumps_minimal_2ranks.log 2>&1
echo "minimal 2 ranks exit=$?"
tail -n 35 mumps_minimal_2ranks.log
```

If global MUMPS crashes even on this small tridiagonal system, the
failure is independently reproducible in the PETSc/MUMPS stack and
should be treated as an environment/native library issue (preserving
the original environment). If it passes, that **does not rule out**
MUMPS on the large, indefinite 3300x3300 system.

### Then isolate the pseudo-TS preconditioner only

The code now accepts an **environment-variable diagnostic**:
`SCRIMP_INIT_PC_TYPE=bjacobi`, which overrides just the *initial*
TS's PC after its `setFromOptions()`, leaving the main TS unchanged.
It is NOT a production solver configuration and is NOT evidence of
valid numerical integration.

```bash
SCRIMP_INIT_PC_TYPE=bjacobi PYTHONFAULTHANDLER=1 SCRIMP_MPI_TRACE=1 \
timeout -k 5s 120s mpiexec -n 2 \
  python -u sandbox/viscous_swe_2d_mpi.py \
  --case open --mesh 0.2 --final-time 0.5 \
  --save-dt 0.1 --max-step 0.01 \
  > mpi_open_init_bjacobi_2ranks.log 2>&1

echo "pseudo-TS BJACOBI exit=$?"
grep -E 'init_step:|solve:|RuntimeError|Fatal Python error|SIGSEGV|Aborted|PC override' \
  mpi_open_init_bjacobi_2ranks.log | tail -n 45
tail -n 60 mpi_open_init_bjacobi_2ranks.log
```

Interpretation:
- **Minimal MUMPS fails**: strong signal of a native factor lifecycle
  issue independent of SCRIMP/GetFEM.
- **Minimal MUMPS passes; BJACOBI makes pseudo-TS reset succeed**:
  factor destruction is specific to MUMPS **in this nonlinear
  workload**, or earlier corruption is only revealed by MUMPS;
  gather state/stack diagnostics before claiming a MUMPS bug.
- **BJACOBI also crashes in TS.reset**: investigate other PETSc
  matrices, state copies and allocator damage. GDB on BJACOBI
  should show whether the factor-path changed.
- **BJACOBI diverges before reset**: inconclusive. Some DAE linear
  systems need global LU and block Jacobi may not converge.

The callback now prints explicit messages before and after `TS.reset`
and `TS.destroy`; logs can distinguish an actual reset corruption
from a later destroy. If further profiling is required, use
`-malloc_debug` in a dedicated PETSc test build; this is not enabled
by default in a Conda release build.

Do **not** remove `TS.reset()` merely to conceal the error: freeing
a corrupted MUMPS factor later would still be unsafe. Do not merge or
tag while the memory and numerical consistency problems persist.

## 2026-11-13 — Standalone MUMPS reset succeeds, preliminary BJACOBI test inconclusive

User ran `sandbox/test_petsc_mumps_reset.py --solver mumps --cycles 3`
on both **1 rank** and **2 ranks**. Every KSP solve converged
(`reason=4`), and all KSP reset/destroy cycles completed. Both exit
codes were **0**. Thus an ordinary MPIAIJ tridiagonal factorization and
MUMPS destruction **works** in the user's PETSc/MPI environment.
The native failure requires additional circumstances (the nonlinear
SWE Jacobian, TS/PC lifecycle, earlier heap corruption, etc.).

The first 2-rank SWE `SCRIMP_INIT_PC_TYPE=bjacobi` experiment ended
with PETSc **error 91 / DIVERGED_STEP_REJECTED**, not SIGSEGV. But
this experiment was **not a fair iterative BJACOBI control**:
the driver specifies `ksp_type=preonly` for LU/MUMPS, and the initial
diagnostic hook changed only `PC=bjacobi`. Thus there was no
GMRES iteration—only a one-shot block-Jacobi operation.

The branch diagnostic hook now changes **both** settings, and only
inside the initial pseudo-TS:
`SCRIMP_INIT_PC_TYPE=bjacobi` -> `KSP=gmres`,
`PC=bjacobi`, with `rtol=1e-8, atol=1e-12, max_it=500`.
Normal `LU/MUMPS + KSPPREONLY` operation is unchanged.

**Rerun, after pulling the dedicated branch:**

```bash
git pull --ff-only
SCRIMP_INIT_PC_TYPE=bjacobi PYTHONFAULTHANDLER=1 SCRIMP_MPI_TRACE=1 \
timeout -k 5s 120s mpiexec -n 2 \
  python -u sandbox/viscous_swe_2d_mpi.py \
  --case open --mesh 0.2 --final-time 0.5 \
  --save-dt 0.1 --max-step 0.01 \
  > mpi_open_init_gmres_bjacobi_2ranks.log 2>&1
echo "code=$?"
grep 'diagnostic linear solver override\|init_step: after pseudo TS.solve\|init_step: after TS.reset\|DIVERGED_' \
  mpi_open_init_gmres_bjacobi_2ranks.log | tail -n 25
tail -n 70 mpi_open_init_gmres_bjacobi_2ranks.log
```

If GMRES-BJACOBI still diverges it says little about the MUMPS
corruption: the coupled indefinite DAE may require a global direct
factorization or a physics-aware block preconditioner.

### Separate mathematical error in the nonlinear residual/Jacobian

The current historical `IFunction()` branch for state-dependent mass
builds `F = M_tangent(z) zd + K_tangent(z) z + rhs_nl(z)`.
For a GetFEM nonlinear term with residual `N(z)` and
`rhs_nl(z)=-N(z)`, `K_tangent(z)=N'(z)`, this gives
`N'(z)z-N(z)` instead of `N(z)`.
For a homogeneous cubic term `N(z)`, Euler's identity yields
`N'(z)z=3N(z)`, hence the implemented residual evaluates to
`2N(z)` rather than `N(z)`. The 2D SWE contains
`h*(Gyro(p)*e_p).Test_p`, which is cubic in model unknowns.

PETSc also requires `IJacobian=F_z+shift*F_zdot`.
For state-dependent `M(z)`, `F_z` contains
`(dM/dz)[delta_z] zd`, which is absent from the present
`shift*M_tangent + K_tangent` code. The high nonlinear coefficients
seen in the single-rank and GDB two-rank logs are consistent with an
inconsistent residual/Newton tangent, but do NOT prove all issues stem
from this alone. Before changing production semantics, create a
finite-difference directional Jacobian check for the full SWE system
at the actual initial state and nonlinear Newton iterates.

For PETSc's contract see:
https://petsc.org/main/manualpages/TS/TSSetIJacobian/

Do not merge PR #22 or tag release until the native memory issue,
time-integration convergence and port-Hamiltonian energy balance are
independently validated.

## 2026-11-13 — Correct GMRES/BJACOBI trial still diverges (PETSc 91)

The user reran the two-rank Saint-Venant initializer after the
diagnostic fix. **Both MPI ranks confirm
`KSP=gmres, PC=bjacobi`**; both then fail in `TS.solve()` with
PETSc **error 91 / `DIVERGED_STEP_REJECTED`**. No native SIGSEGV
or MUMPS factor teardown error is observed in the supplied output,
but the initial pseudo-TS did not complete, so this is not a
meaningful memory-safety acceptance test. It demonstrates only that
this GMRES+block-Jacobi configuration does not solve the nonlinear
DAE in this workload.

The independent PETSc/MUMPS tridiagonal test had already passed
three full factor/solve/reset cycles for 1 and 2 MPI ranks, so
there is no evidence of a generic MPI/MUMPS ABI incompatibility.

### Next gate: finite differences of the **full 2D SWE implicit residual**

A new diagnostic **`sandbox/check_swe_implicit_jacobian.py`**
reconstructs the full nonlinear 3300-DOF Saint-Venant model and
calls SCRIMP's actual `IFunction` and `IJacobian` without
creating a PETSc TS/SNES/KSP, using exactly the same GetFEM/MPI
transfer and nonlinear bricks as the simulation.

It checks, for a deterministic distributed direction `v`,
`J(z,zd,sigma) v` against the central difference

```text
[F(z+eps*v, zd+eps*sigma*v) - F(z-eps*v, zd-eps*sigma*v)] / (2*eps).
```

The tested points include `zd=0`, `zd!=0` and `sigma=0,1`,
which separately expose incorrect state derivatives of the
nonlinear force and missing `dM(z)/dz * zd` terms. The script
reports relative directional errors without failing by default:
**exit 0 does not mean the Jacobian is correct**.
Only `--fail-above=...` turns it into a pass/fail gate after a
mathematical correction.

```bash
git pull --ff-only

timeout -k 5s 120s mpiexec -n 1 \
  python -u sandbox/check_swe_implicit_jacobian.py \
  --mesh 0.2 --epsilon 1e-5 \
  > fd_swe_1rank.log 2>&1
echo "FD 1 rank exit=$?"
tail -n 75 fd_swe_1rank.log

timeout -k 5s 120s mpiexec -n 2 \
  python -u sandbox/check_swe_implicit_jacobian.py \
  --mesh 0.2 --epsilon 1e-5 \
  > fd_swe_2ranks.log 2>&1
echo "FD 2 ranks exit=$?"
tail -n 75 fd_swe_2ranks.log
```

If FD errors are large at 1 and 2 ranks, the residual/Jacobian
error predates MPI factorization and must be corrected before
interpreting a nonlinear convergence failure. Repeat the check
with `--epsilon 1e-4` and/or `--state-perturb 0.05` to check
robustness rather than over-interpreting floating-point roundoff.

**This tool is diagnostic only and has not been run in a native
GetFEM/PETSc environment**. If it fails to initialize or crashes
in PETSc/Mat operations, capture the log and treat that as its own
finding, not a numeric test result.

**Production formulas remain unchanged.** Avoid a blind rewrite of
nonlinear mass or RHS semantics until independent GetFEM residual
sign conventions and this directional finite-difference result
are checked. Even a consistent Jacobian is not sufficient to prove
that the model itself implements the intended pH equations; energy
and boundary-power balances must be separately validated.


## 2026-11-13 — Full nonlinear SWE FD results: 6.8% serial, 100% on 2 ranks

The user executed the full directional finite-difference audit
`sandbox/check_swe_implicit_jacobian.py` with `--mesh 0.2`,
`--epsilon 1e-5` at 1 and 2 ranks. Both runs produced structured
results without a native crash. For `zd=0`, `shift=0`:

| Quantity | 1 rank | 2 ranks |
| --- | ---: | ---: |
| N | 3300 | 3300 |
| `norm(J @ v)` | 0.3125827291 | 0.3125827291 |
| `norm(FD(F) @ v)` | 0.3122621337 | 0.0243832267 |
| `norm(J @ v - FD(F) @ v)` | 0.0213927492 | 0.3141442007 |
| Relative mismatch | **0.06843868** | **1.00499539** |

Other combinations `(zd=0 or nonzero, shift=0 or 1)` yield
approximately 6.8% mismatch at 1 rank and 100% at 2 ranks.
The PETSc Jacobian action is nearly rank-independent, but the
finite-difference response of the **actual implemented residual**
changes drastically with the MPI rank count.

This is strong evidence of BOTH a nonlinear residual-vs-Jacobian
inconsistency present already in serial **and** an additional MPI
assembly/state-transfer issue, though before concluding we must
validate the FD test and vary `epsilon`. The current source path
`K_tangent(z)z + gf_rhs(z)` is mathematically suspect, and its
`rhs(z)` may not be represented consistently across MPI ranks.

**An important source-contract check:** GetFEM's parallelization
documentation says its MPI level-2 RHS vectors are summed and
replicated on all ranks, but tangent matrices remain distributed:
https://getfem.org/userdoc/parallel.html .
Thus **do not assume** the RHS is partial merely because the
FD mismatch increases at 2 ranks. The difference might also be in
trial-state synchronization or how SCRIMP combines the nonlinear
GetFEM tangent, RHS and linear terms. The present
`SCRIMP_GETFEM_MATRIX_MODE=replicated` behavior passed targeted
mass/tangent tests but does not establish the full GetFEM operator
contract for every brick.

### Enhanced diagnostic, without production formula changes

`sandbox/check_swe_implicit_jacobian.py` now additionally records:

1. `petsc_to_getfem_state_scatter_error_by_rank`:
   reconstructs the full PETSc trial state from owned row arrays
   and compares it with GetFEM's current variables after sync.
2. `getfem_rhs_derivative_rank_difference` and
   `getfem_rhs_rank_signatures`: compare the **raw GetFEM full RHS**
   directional derivative between processes, before transfer to PETSc.
3. `getfem_rhs_derivative_vs_minus_nonlinear_tangent_relative`:
   separately check the expected GetFEM relation
   `d RHS(z) / dz = - K_nl(z)` by applying the actual nonlinear
   matrix to the same distributed direction.
4. `error_norm_by_equation`: identify the residual-equation blocks
   (`h,p,e_h,e_p,U_n,Y_n,U,Y`) responsible for the total FD mismatch.

The enhanced probe creates **no TS/SNES/KSP/MUMPS**, does not change
production code, and by default returns exit 0 even for large
numerical discrepancies.

Run on 1 and 2 ranks after `git pull --ff-only`:

```bash
for n in 1 2; do
    timeout -k 5s 120s mpiexec -n "$n" \
        python -u sandbox/check_swe_implicit_jacobian.py \
        --mesh 0.2 --epsilon 1e-5 \
        > "fd_swe_detailed_${n}ranks.log" 2>&1
    echo "FD detailed MPI $n: exit $?"
    tail -n 135 "fd_swe_detailed_${n}ranks.log"
done
```

Then vary `epsilon` on one rank with `--epsilon 1e-4`
and `--epsilon 1e-6` to distinguish truncation/round-off from
a genuine inconsistency. If MPI RHS signatures differ between ranks,
investigate actual GetFEM assembly; if they agree but FD differs,
check `K_tangent*z+rhs` and the PETSc state scatter with equation
blocks. **Do not change production sign conventions** based only on
these preliminary numbers.

## 2026-11-13 — Detailed 1/2-rank FD audit isolates **factor-two RHS**, not state scatter

The user completed the enhanced full Saint-Venant 2D FD audit at
`epsilon=1e-5`, on **one and two ranks**, exit code **0** in both cases.
The additional diagnostics show:

| Quantity | MPI 1 | MPI 2 |
| --- | ---: | ---: |
| `||Jv||` | 0.31258272908 | 0.31258272908 |
| `||FD(F)v||` | 0.31226213370 | 0.02438322670 |
| `||Jv - FD(F)v|| / max(norms)` | **0.06843868** | **1.00499539** |
| `||raw GetFEM RHS(z+)||` | 0.23407250883 | **0.46814501765** |
| `||D raw GetFEM RHS(z)v||` | 0.31236350142 | **0.62472700283** |
| `||D RHS(z)v + K_nl(z)v||/||K_nl(z)v||` | **6.14e-12** | **1.00000000** |
| PETSc→GetFEM trial-state scatter discrepancy | `[0]` | `[0, 0]` |
| RHS directional derivative discrepancy across ranks | 0 | 0 |

The *full* raw GetFEM RHS arrays have identical fingerprints across
the two ranks, but the norm is **twice the serial value**.
Moreover, the directional derivative equals `-K_nl v` to
`6e-12` in serial and is consistent with **`-2 K_nl v`** at
two ranks. Thus this is an **assembly-level factor-of-two
RHS/tangent inconsistency in this full GetFEM SWE model**, not a
PETSc state scatter failure or a double-insertion of the RHS
into PETSc (which uses owned-row INSERT_VALUES).

**Equation blocks:**
- Serial discrepancy `||Jv-FD(F)v||=0.021393` is almost entirely
  `U: 0.0200141` and `p: 0.0075554`.
- Two-rank discrepancy `0.314144` is dominated by
  `h: 0.19375`, `p: 0.19725`, `Y: 0.121025`,
  `Y_n: 0.0750314`, `U: 0.0400283`.
- Extra derivative in the PETSc Jacobian from state-dependent mass
  is still structurally missing; this check does not validate the
  physical equations.

### Diagnostic-only test of the factor, plus isolated GetFEM brick reproduction

The full probe now offers `--rhs-scale`; it applies a numeric
factor to SCRIMP's **PETSc copy** of the GetFEM RHS *after each
assembly* in this diagnostic instance ONLY, and leaves the raw GetFEM
RHS reporting unchanged. It does not alter the package or solver
behavior. In particular **`--rhs-scale=0.5` at MPI=2** tests
whether the FD error returns to the ~6.8% serial discrepancy.

A new `sandbox/check_getfem_mpi_rhs.py` independently builds
one-dimensional GetFEM `Model` objects with one `add_nonlinear_term`
brick each:
`u Test_u`, `u^2 Test_u`, `u^3 Test_u`.
It verifies that the integral of the derivative of negative
`gf.Model.rhs()` equals the integral of the GetFEM tangent
matrix after conversion to PETSc. The exact values at `u=1`
on [0,1] are 1, 2 and 3 respectively. It reports RHS/tangent
ratios without failing by default. Run with 1/2/4 MPI ranks to
distinguish a generic GetFEM MPI nonlinear RHS doubling
from an interaction specific to the multi-brick SWE.

```bash
git pull --ff-only

for n in 1 2 4; do
    timeout -k 5s 60s mpiexec -n "$n" \
        python -u sandbox/check_getfem_mpi_rhs.py \
        > "getfem_rhs_${n}rank.log" 2>&1
    echo "GetFEM RHS MPI ${n}: exit=$?"
    tail -n 65 "getfem_rhs_${n}rank.log"
done

timeout -k 5s 120s mpiexec -n 2 \
    python -u sandbox/check_swe_implicit_jacobian.py \
    --mesh 0.2 --epsilon 1e-5 --rhs-scale 0.5 \
    > fd_swe_2rank_half_rhs.log 2>&1
echo "2-rank half-RHS diagnostic: exit=$?"
tail -n 145 fd_swe_2rank_half_rhs.log
```

**Interpretation:** If the 1D models reproduce a factor `P`,
the issue is broadly in GetFEM RHS-vs-tangent assembly under this
installation; do not silently normalize every RHS without verifying
source, boundary, and externally imposed terms. If the 1D models
do NOT reproduce it, isolate the individual SWE GetFEM bricks and
assembly calls before modifying `scrimp/dphs.py`.
If the temporary `0.5` scaling returns the FD error to
~6.8%, the additional MPI discrepancy is explained by RHS
doubling; however the serial nonlinear residual/Jacobian mismatch
is still *not fixed*. No claims of physical validity or native
memory safety follow from this diagnostic.


## 2026-11-13 — New MUMPS 5.8.2 / PETSc 3.25.5 / MPICH 5.0.1 environment

The user created a coherent independent `scrimp-mumps-test` environment:

| Library | Original `scrimp` | `scrimp-mumps-test` |
| --- | --- | --- |
| MPICH | 4.3.1 | 5.0.1 |
| mpi4py | 4.1.0 | 4.1.2 |
| MUMPS-MPI | 5.7.3 | 5.8.2 |
| PETSc complex | 3.21.6 | 3.25.5 |
| petsc4py | 3.21.6 | 3.25.5 |

`sandbox/test_petsc_mumps_reset.py --solver mumps --cycles 10`
was run in the NEW environment at **1, 2 and 4 MPI ranks**.
Every process completed **all ten** solve/reset/destroy cycles
without error. Thus the new native numeric stack is functional on
the small tridiagonal benchmark. However, the ORIGINAL stack had
already passed the same sort of benchmark at 1/2 ranks (3 cycles).
So these results do **not** demonstrate that upgrading MUMPS fixes
the real Saint-Venant native crash.

A single-component `mumps-mpi=5.8.2` upgrade in the original env
was UNSATISFIABLE: PETSc 3.21.6 was built with the explicit conda
constraint `mumps-mpi >=5.7.3,<5.7.4`. The fresh environment
necessarily changes PETSc and MPICH too; it is not an isolated
MUMPS-only A/B experiment.

### Next: capture the *real* Newton matrix in the original environment

PETSc's `MatView` and `MatLoad` support collective binary
storage/reload of MPIAIJ matrices, independent of GetFEM:
https://petsc.org/main/manualpages/Mat/MatLoad/

We added an **opt-in** hook in `DPHS.IJacobian()` to export
assembled LU-preconditioner matrices **immediately before PETSc
factorization** when the env var
`SCRIMP_JACOBIAN_SNAPSHOT_PREFIX` is set. The hook automatically
writes a bounded set of power-of-two callback snapshots (numbers
1,2,4,8,16,32), stored as separate PETSc binary files. With the
env var **unset**, production execution is unchanged.

The snapshot captures the numerical operator, NOT the original
TS/SNES/KSP object history, the Newton state, or prior native heap
damage. A standalone replay cannot prove the full TS solver safe.

**Capture using ORIGINAL `scrimp` (locally built MPI GetFEM):**

```bash
git pull --ff-only
conda activate scrimp
mkdir -p sandbox/outputs/mumps_ab

SCRIMP_JACOBIAN_SNAPSHOT_PREFIX=sandbox/outputs/mumps_ab/jac \
SCRIMP_MPI_TRACE=1 PYTHONFAULTHANDLER=1 \
timeout -k 5s 120s mpiexec -n 2 \
  python -u sandbox/viscous_swe_2d_mpi.py \
    --case open --mesh 0.2 --final-time 0.5 \
    --save-dt 0.1 --max-step 0.01 \
  > mpi_original_mumps_snapshot.log 2>&1

echo "original SWE exit=$?"
ls -lh sandbox/outputs/mumps_ab/jac_*.petscbin
tail -n 35 mpi_original_mumps_snapshot.log
```

A 139/134 failure here is **expected** from the unresolved bug.
The important output is at least one complete binary matrix file
and corresponding `snapshot ... written` lines on both ranks.

**Replay one binary snapshot in BOTH environments using the
same MPI size (2).** This script imports *only* PETSc; it does
NOT need GetFEM, SCRIMP or Gmsh in the new environment.

```bash
MATRIX=sandbox/outputs/mumps_ab/jac_004.petscbin
# If 004 doesn't exist, choose one file printed by ls (001/002/008).
ls -lh "$MATRIX"

for env in scrimp scrimp-mumps-test; do
  conda activate "$env"
  timeout -k 5s 90s mpiexec -n 2 \
    python -u sandbox/replay_petsc_mumps_snapshot.py \
      --matrix "$MATRIX" --cycles 3 \
    > "mumps_replay_${env}.log" 2>&1
  echo "${env} replay exit=$?"
  tail -n 70 "mumps_replay_${env}.log"
done
```

The replay prints `size`, `Frobenius`, and the norm of
`b=A*ones`; these should agree across old/new environments.
It then performs a fresh LU/MUMPS KSP setup, solve, reset and
destroy in each cycle; compare native exit, PETSc errors and
convergence reasons.

If a snapshot fails to load because of an incompatible binary
format, index-size setting or scalar-type mismatch, STOP and report
the loader error; do not silently convert the file or infer a
MUMPS defect. The current envs both use complex PETSc, but still
verify the loaded norm before comparing solves.

**Interpretation**:
- Old replay crashes while new replay succeeds on the **same
  matrix**: strong evidence the problem is in the native
  PETSc/MUMPS/MPICH software stack, though not solely in MUMPS
  since the three versions changed together.
- Both replays succeed: the original crash likely depends on
  stateful TS/SNES lifecycle, past nonlinear assemblies,
  hidden heap damage, or a different Newton matrix snapshot.
  Repeat on later snapshots (008/016/032) before concluding.
- Both fail similarly: suspect a singular/unphysical Jacobian,
  native numerical edge case, or a common PETSc/MUMPS failure
  mode; compare PETSc diagnostics before drawing conclusions.

The independent GetFEM MPI RHS multiplication-by-rank-count
and the ~6.8% serial SWE residual/Jacobian mismatch are separate
blockers. A MUMPS upgrade by itself cannot fix either. No merge,
release or issue closure until physical balance, successful t_f,
MPI parity and native robustness are all demonstrated.


## 2026-11-13 — Both actual SWE Jacobians replay successfully on both MUMPS versions

The user replayed `jac_002.petscbin` at 2 ranks on:
- **Old:** MPICH 4.3.1, PETSc 3.21.6, MUMPS 5.7.3.
- **New:** MPICH 5.0.1, PETSc 3.25.5, MUMPS 5.8.2.

As for `jac_001`, **both** environments loaded `jac_002`
identically (`3300x3300`, Frobenius norm **1304526700.41**,
`||b=A*ones||=2.08510766e9`) and completed all **three**
fresh KSP_PREONLY + PCLU/MUMPS factor/solve/reset/destroy cycles
on **both ranks**, exit code **0**, converged reason=4.

**Strong outcome:** Neither the first nor second captured SWE
Jacobian independently causes a fault in old MUMPS 5.7.3, even when
factored and destroyed repeatedly. Upgrading MUMPS is **not**
experimentally demonstrated as a crash fix. The original full
2-rank SWE (with the GetFEM assembly, evolving Newton state,
TS/SNES/KSP object ownership) still SIGSEGVs after the second
`IJacobian`.

### Stronger lifecycle isolation: persistent KSP, mutable Mat, changing sparsity

The existing PETSc-only replay has a new `--next-matrix` mode.
It loads BOTH snapshots, creates **one persistent KSP/PCLU/MUMPS**
and **one mutable MPIAIJ Mat**, then alternately copies each
captured matrix *into the same Mat handle* using
`MatCopy(..., DIFFERENT_NONZERO_PATTERN)`. The KSP is **not reset
between Newton-like steps**. Each RHS is constructed consistently
as `b=A*ones`; the script checks convergence and a true linear
residual. Only after the sequence does it call `KSP.reset()`
and `KSP.destroy()`, the same native factor teardown implicated
in the original GDB stack. It still does NOT instantiate TS/SNES
or perform any GetFEM assembly.

Run from the repository directory, using the already captured files:

```bash
git pull --ff-only

FIRST=sandbox/outputs/mumps_ab/jac_001.petscbin
SECOND=sandbox/outputs/mumps_ab/jac_002.petscbin
for env in scrimp scrimp-mumps-test; do
  conda activate "$env"
  timeout -k 5s 90s mpiexec -n 2 \
    python -u sandbox/replay_petsc_mumps_snapshot.py \
      --matrix "$FIRST" --next-matrix "$SECOND" --cycles 6 \
    > "mumps_persistent_${env}.log" 2>&1
  echo "${env} persistent KSP exit=$?"
  tail -n 65 "mumps_persistent_${env}.log"
done
```

Interpret `--cycles 6` as SIX successive factorizations/solves
(`001,002,001,002,001,002`) **inside the same KSP/PC**.
The logs show whether failure occurs in `MatCopy`,
`KSP.solve` or final `KSP.reset`. If both pass, the
next suspect is the TS/SNES lifecycle or earlier corruption
during GetFEM/sCRIMP callbacks, NOT the intrinsic ability of
MUMPS 5.7.3 to factor the captured matrices. A new PETSc-only
TS/SNES lifecycle reproducer and a consistent nonlinear residual
check would then be a better next step than further MUMPS version
changes.

**Caution:** The driver uses a manufactured consistent RHS for
linear algebra testing; it cannot validate the physical SWE
residual, implicit Jacobian derivative or power balance.
The source-version smoke still exits 139; do not merge PR #22,
close #6 or tag v1.2.1 on the basis of these native tests.


## 2026-11-13 — Persistent KSP/PC/MUMPS with alternating SWE Jacobians also passes

The user ran the new `--next-matrix` PETSc-only replay on
`jac_001.petscbin` and `jac_002.petscbin` using TWO MPI ranks
in BOTH Conda environments:
- **Old** `scrimp`: PETSc 3.21.6, MUMPS 5.7.3, MPICH 4.3.1.
- **New** `scrimp-mumps-test`: PETSc 3.25.5, MUMPS 5.8.2,
  MPICH 5.0.1.

In **both** environments, 6/6 transitions
`001 → 002 → 001 → 002 → 001 → 002` succeeded with the **same**
persistent PETSc KSP/PC, and the **same mutable** MPIAIJ Mat handle.
Each update used `MatCopy(DIFFERENT_NONZERO_PATTERN)`, so KSP's
MUMPS factorization was reused/rebuilt appropriately.
Every KSP solve converged; maximum true relative linear residual:
- Old: **3.3047122e-16**.
- New: **3.10989255e-16**.

After all 6 transitions, both ranks completed the final
`KSP.reset()` and `KSP.destroy()` without native errors
(exit 0 in both environments). This greatly weakens the hypothesis
that **MUMPS 5.7.3 alone**, or a normal repeated LU factor
with changing sparsity, causes the original SWE SIGSEGV.
It does **not** eliminate TS/SNES callback lifetime errors or
GetFEM-native heap corruption in the full nonlinear calculation.

### Next priority: verify the mathematically correct nonlinear residual

Serial SWE directional finite-difference comparison previously
measured **~6.84%** mismatch between `IFunction` and
`IJacobian`, even before introducing MPI. The two-rank excess
mismatch was **fully explained** by the established GetFEM RHS
multiplication-by-MPI-size factor. The diagnostic already contains a
`--residual-form direct_getfem` alternative; **it overrides only
the isolated FD experiment's `IFunction`**, leaving SCRIMP
production equations and time solver untouched.

For fixed `zd=0`, one expects the GetFEM weak-form residual
`F(z,zd) = M(z) zd + K_linear z - rhs_nonlinear(z)`, hence
`F_z = K_linear + K_nonlinear` and
`F_zd = M(z)`. In contrast, SCRIMP's legacy nonlinear path
evaluates `K_tangent(z) z + rhs_nonlinear(z)`.
For `zd ≠ 0`, the direct-residual diagnostic still uses the
historical `IJacobian`, which is missing the
`(d M(z)/dz)[v] zd` term. Therefore the four test cases
`zd∈{0,nonzero}`, `shift∈{0,1}` *separately* test the
nonlinear force and state-dependent mass derivatives.

**Run in ORIGINAL `scrimp` (has locally compiled MPI GetFEM):**

```bash
git pull --ff-only
conda activate scrimp

timeout -k 5s 120s mpiexec -n 1 \
  python -u sandbox/check_swe_implicit_jacobian.py \
    --mesh 0.2 --epsilon 1e-5 \
    --residual-form direct_getfem \
  > fd_direct_1rank.log 2>&1
echo "direct FD 1 rank exit=$?"
tail -n 145 fd_direct_1rank.log

timeout -k 5s 120s mpiexec -n 2 \
  python -u sandbox/check_swe_implicit_jacobian.py \
    --mesh 0.2 --epsilon 1e-5 \
    --residual-form direct_getfem --rhs-scale 0.5 \
  > fd_direct_2ranks.log 2>&1
echo "direct FD 2 ranks exit=$?"
tail -n 145 fd_direct_2ranks.log
```

The `--rhs-scale 0.5` switch is an opt-in normalization applied to
PETSc's RHS **inside the test**, to compensate for the independently
measured GetFEM factor 2. In both reports, the
`getfem_rhs_derivative_vs_minus_nonlinear_tangent_relative`
diagnostic still compares the **unscaled raw GetFEM RHS** to
the tangent; it will therefore still report `1.0` at 2 ranks
even if the *PETSc residual's* derivative has been corrected.
Do **not** interpret that one raw diagnostic number as a failure of
the scaled FD experiment.

The key observations are `relative_error` and
`error_norm_by_equation`, especially for the `zd=0` rows.
If the direct weak-form residual's `zd=0` derivative agrees
with IJacobian in both MPI sizes, then a targeted change to
SCRIMP's production nonlinear `IFunction` is justified
(after verifying boundary/control source sign conventions).
`zd ≠ 0` will still require a mass-derivative Jacobian fix
and/or another consistent DAE treatment.

**No further MUMPS upgrades or merges on this evidence alone.**
Keep issue #6 open, PR #22 draft and SCRIMP master unchanged.


## 2026-11-13 — Direct GetFEM residual removes serial FD mismatch; mass Hessian isolated

**User-verified**, with the sandbox-only direct residual
`--residual-form direct_getfem` on the full nonlinear SWE model,
3300 DOFs, `epsilon=1e-5`:

| State derivative | MPI=1 | MPI=2 (sandbox RHS scale=0.5) |
| --- | ---: | ---: |
| `zd=0, shift=0` relative FD error | **6.48666598e-12** | **6.48666598e-12** |
| `zd=0, shift=1` relative FD error | **6.47095289e-12** | **6.47095289e-12** |
| `zd!=0, shift=0` relative FD error | **0.00199208730** | **0.00199208730** |
| `zd!=0, shift=1` relative FD error | **0.00198726742** | **0.00198726742** |

The existing production `IFunction` with state-dependent mass
uses `M(z)zd + K_tangent(z)z + rhs_GetFEM(z)`. The diagnostic
replaces it **only inside the probe** with
`M(z)zd + K_linear z - rhs_GetFEM(z)` (using the observed
GetFEM weak-form sign convention). The latter derivative matches
the existing `IJacobian` to ~`6.5e-12` for `zd=0` at 1 and
2 MPI ranks. The known factor-of-two raw GetFEM RHS at two ranks
is compensated by `--rhs-scale=0.5` inside the probe only.

For `zd!=0`, the discrepancy **0.0006227200101** is confined
to the `p` equation block, in the nonlinear mass brick
`h*p.Test_p` marked `dt=True`. These observations precisely
support a missing `(D_z M(z)[v])zd` contribution in `IJacobian`
and do NOT suggest a remaining error in the `U` boundary equation
after fixing residual evaluation.

### Experimental numerical test of the missing mass Hessian

For a smooth GetFEM flow primitive `Q(z)`, with
`M(z)=D Q(z)`, the time part of the SCRIMP residual is
`M(z)zd`. Its state derivative is
`D M(z)[v] zd = D M(z)[zd] v` by symmetry of the
second Fréchet derivative `D²Q`. Thus the missing **assembled
matrix** can be diagnosed with just two nonlinear mass tangent
assemblies at `z±eta*zd`:

```text
C(z,zd) = [M(z+eta*zd) - M(z-eta*zd)] / (2*eta)
J_corrected = J_legacy + C(z,zd).
```

This formula is exact up to floating-point error for the
bilinear `h*p.Test_p` brick, and is a second-order approximation
for smooth higher-order GetFEM flow primitives.

The FD audit now accepts **`--mass-jacobian finite_difference`**
and **`--mass-epsilon 1e-4`**. This is entirely sandbox-local:
it does NOT modify the production SCRIMP Jacobian. It also accepts
`--zd-direction independent` (default), using an independent
deterministic time derivative rather than `zd || v`, so the test
checks genuinely off-axis Hessian terms. `--zd-direction aligned`
reproduces the old `zd=11.0088*v` experiment.

Run with the original GetFEM-enabled `scrimp` env:

```bash
git pull --ff-only
conda activate scrimp
for n in 1 2; do
  scale=1.0
  if [ "$n" = 2 ]; then scale=0.5; fi
  timeout -k 5s 120s mpiexec -n "$n" \
    python -u sandbox/check_swe_implicit_jacobian.py \
      --mesh 0.2 --epsilon 1e-5 \
      --residual-form direct_getfem \
      --mass-jacobian finite_difference \
      --zd-direction independent \
      --rhs-scale "$scale" \
    > "fd_direct_mass_${n}ranks.log" 2>&1
  echo "MPI=$n; exit=$?"
  tail -n 140 "fd_direct_mass_${n}ranks.log"
done
```

Repeat with `--state-perturb 0.02` or `--epsilon 1e-4`
to avoid overfitting to a single state/direction. The key criterion
is a **small relative error in all four cases**, especially when
`zd!=0`. This test is not a replacement for a rigorously assembled
production Hessian: it is a diagnostic step to confirm the exact
missing term before changing any generic DPHS behavior.

**Warnings:** This check still does not validate the *intended*
Saint-Venant PDE semantics (whether its `dt=True` brick should
represent `d(hp)/dt` or `h dp/dt`), time integration success,
energy/power conservation, or native MPI memory safety. Those are
separate acceptance gates. No merge, release tag or issue closure.


## 2026-11-13 — Why a replicated GetFEM RHS must NOT automatically be divided by P

GetFEM's parallelization contract:
https://getfem.org/userdoc/parallel.html

GetFEM with MPI level 2 *normally* partitions mesh integration using
METIS. The local weak-form vectors `b_r` are summed in GetFEM, and
the global `b = sum_r b_r` is replicated **in full** on every
process. This is intentional and correct. The globally numbered
tangent matrix, in contrast, is documented as containing only
local contributions. PETSc must insert the globally reduced GetFEM
RHS in *owned rows* (INSERT_VALUES) and not perform a second SUM.

**The mathematical factor P arises only if contributions are already
complete on every rank before GetFEM's MPI reduction**, or if
some other stage performs an equivalent redundant reduction:

```text
Correct:   each rank computes b_r, SUM_MPI(b_r) = b
Observed:  SUM_MPI(b) = P*b
```

In the user's GetFEM 5.4.2 level-2 build, 1D
`u, u², u³` weak forms all had RHS integrated derivatives
multiplied by P=1,2,4 while SCRIMP's converted tangent had its
expected unscaled action. The 2D Saint-Venant FD mismatch on 2 ranks
vanished after the *sandbox* applied a 0.5 RHS multiplier. These
results establish the factor, but **do not yet establish the precise
GetFEM C++ path** that duplicates the contributions: it could be
an ineffective generic assembly partition, double reduction or a
feature of this local build. Do not claim the RHS *being replicated*
is a bug: it is expected and explicitly documented.

**Added raw GetFEM matrix diagnostic** to
`sandbox/check_getfem_mpi_rhs.py`:
`raw_getfem_tangent_by_rank` (pre-PETSc rankwise integral,
`nnz`, row-action digest), `sum_raw_getfem_tangent_integrals`
and `raw_rhs_minus_tangent_ratio_at_each_rank`. Compare at 1/2/4
ranks to see whether the *raw* GetFEM tangents are themselves
already complete replicas on every process instead of distributed
contributions.

**SCRIMP compatibility correction** (draft PR #22, not yet locally
tested): the default `SCRIMP_GETFEM_RHS_MODE=auto` now calibrates
the library behavior **once per process** by assembling an exactly
known scalar generic reference form on [0,1]:
`rhs_probe*Test_rhs_probe` with `rhs_probe=1`. The exact integral
is 1, so the GetFEM RHS entry sum must be -1. The observed ratio
may be **1** (properly reduced, scale=1) or **P** (over-reduced,
scale=1/P). An unexpected value or disagreement between MPI ranks
raises an error instead of silently applying a wrong correction.
The calibration is skipped in sequential mode and when MPI
level is not 2. Modes `raw` and `divide_by_mpi_size` remain
explicit overrides.

**Limitation**: this calibration applies to **GetFEM generic
nonlinear weak-form terms**. It does not prove that *all* models,
particularly explicit RHS, master-only constraints and atypical
source bricks, share exactly the same over-reduction. Check those
separately before generalizing the correction to arbitrary mixed
bricks. The proper long-term fix, if GetFEM generic assembly is
confirmed to skip MPI mesh partitioning, is upstream in GetFEM
and not a universal 1/P factor in SCRIMP.

### Reproduce and verify the mechanism with new diagnostics

```bash
git pull --ff-only
conda activate scrimp

for n in 1 2 4; do
  timeout -k 5s 120s mpiexec -n "$n" \
    python -u sandbox/check_getfem_mpi_rhs.py \
    > "getfem_rhs_partition_${n}rank.log" 2>&1
  echo "raw GetFEM MPI $n exit=$?"
  tail -n 90 "getfem_rhs_partition_${n}rank.log"
done

for n in 1 2 4; do
  timeout -k 5s 120s mpiexec -n "$n" \
    python -m unittest -v tests.test_mpi_algebra \
    > "algebra_default_${n}rank.log" 2>&1
  echo "SCRIMP tests MPI $n exit=$?"
  tail -n 25 "algebra_default_${n}rank.log"
done

# Test the real default production residual/Jacobian + auto RHS calibration
for n in 1 2; do
  timeout -k 5s 120s mpiexec -n "$n" \
    python -u sandbox/check_swe_implicit_jacobian.py \
    --mesh 0.2 --epsilon 1e-5 \
    --residual-form production --mass-jacobian production \
    --rhs-scale 1.0 --fail-above 1e-7 \
    > "fd_production_${n}rank.log" 2>&1
  echo "Production FD MPI $n exit=$?"
  tail -n 115 "fd_production_${n}rank.log"
done
```

For the last command, inspect `relative_error` for **all four**
`(zd, shift)` combinations, not simply the exit code. Then test
perturbed states (`--state-perturb 0.02`) and time integration
separately. The native MPI SIGSEGV and full power balance remain
unsolved; no merge, release or issue closure.


## 2026-11-13 — User-validated corrected DEFAULT production IFunction/IJacobian and MPI RHS calibration

The user supplied eight native log files from the current draft PR #22:
`fd_production_1rank.log`, `fd_production_2rank.log`,
`algebra_default_{1,2,4}rank.log`,
`getfem_rhs_partition_{1,2,4}rank.log`.

### Verified: production residual/Jacobian contract

**Original 3300-DOF nonlinear Saint-Venant model**, `epsilon=1e-5`,
**actual default** `DPHS.IFunction` and `DPHS.IJacobian` (no
sandbox residual override, no additional RHS factor), with
`--zd-direction independent`. All commands exited 0 and
`--fail-above 1e-7` was satisfied.

| (zd,shift) | 1 MPI rank | 2 MPI ranks |
| --- | ---: | ---: |
| (0,0) | 6.48666598165e-12 | 6.48666598053e-12 |
| (0,1) | 6.47095289044e-12 | 6.47095288755e-12 |
| (nonzero,0) | 8.91839504863e-12 | 8.91840080590e-12 |
| (nonzero,1) | 8.55862075486e-12 | 8.55862303144e-12 |

The GETFEM→PETSc Newton trial-state scatter discrepancy is 0
on both MPI ranks. The *raw*, unnormalized GetFEM RHS diagnostic
still has derivative/tangent discrepancy of ~1 for 2 ranks
(as expected: it reports GetFEM raw before SCRIMP's correction).

### Verified: native GetFEM/PETSc algebra tests

**All 11 tests per process PASS on MPI 1, 2 and 4**, including
real GetFEM MPI-2 mass/tangent assembly, the auto-calibrated
GetFEM nonlinear RHS, owned-row insertion, changing sparse patterns,
and persistent tangent matrix handles. PETSc prints harmless
warnings about unused `-info`; not test failures.

### Raw per-rank GetFEM tangent vs internally summed RHS

For the simple 1D GetFEM generic forms on [0,1],
the raw *pre-PETSc* tangent has integral and nnz:

| Test | MPI 1: tangent / raw RHS deriv. | MPI 2 | MPI 4 |
| --- | --- | --- | --- |
| `u*Test_u` | 1 / 1 | 1 / 2 | 1 / 4 |
| `u*u*Test_u` | 2 / 2 | 2 / 4 | 2 / 8 |
| `u*u*u*Test_u` | 3 / 3 | 3 / 6 | 3 / 12 |

In all multi-rank runs each process returned 49 sparse
tangent entries, the full serial tangent integral (1,2,3),
and exactly the SAME digest of the **tangent row-sum vector**
for a given form; RHS derivatives were **identical on all ranks**
and multiplied by P. This strongly suggests redundant
full-domain generic assembly before GetFEM's RHS reduction,
but a row-sum digest alone does **NOT** establish that every
matrix entry is identical. The diagnostic was strengthened to
also emit `raw_tangent_full_csr_digest` over row pointers,
column indices and all matrix coefficients; compare full
digests at 2 and 4 ranks before asserting literal full-matrix
replication.

GetFEM documentation (https://getfem.org/userdoc/parallel.html)
explicitly states that `GETFEM_PARA_LEVEL=2` partitions
integration using METIS, with distributed local tangent
contributions and a **globally summed / replicated RHS**.
A full RHS on all ranks is NORMAL. An extra factor P
indicates duplicated assembly contributions or another
equivalent extra reduction. We have not yet identified
the precise GetFEM C++ implementation path responsible
on this installation. Do not attribute the bug to
MUMPS (the simple test invokes no linear solver) or
to PETSc owned-row insertion.

### Next gates

1. Rerun the strengthened 1D probe at 2 and 4 ranks and
   compare `raw_tangent_full_csr_digest` across processes.
   Identical full hashes would establish that GetFEM returns
   equal sparse matrices before PETSc transfer, but still
   would not establish whether GetFEM already MPI-gathered
   them before exposing `tangent_matrix()`.
2. Test the actual production FD derivative at off-base
   states (e.g. `--state-perturb 0.02`) on 1 and 2 ranks.
3. Run the **full** nonlinear SWE temporal solver again
   (1/2 ranks); record time reached, native crash/safety,
   nonlinear convergence, actual Hamiltonian and integrated
   external/dissipation power. These are not proved by
   exact FD agreement or passing unit tests.

**No master merge, issue #6 closure or v1.2.1 release yet.**


## 2026-11-13 — Could SCRIMP simply be calling the wrong GetFEM assembly mode?

**Question from user:** SCRIMP calls
`gf_model.assembly(option="build_rhs")` and reads the complete
non-distributed `gf_model.rhs()`. Could that be the mistake?

**API/source check:** No wrong API call is established.
GetFEM explicitly documents `Model.assembly("build_rhs")`
as the RHS-only mode and `Model.assembly("build_all")` as
simultaneous RHS + tangent assembly:
https://getfem.org/python/cmdref_Model.html .
A global RHS vector replicated on every rank is ALSO the
documented intended GetFEM MPI-level-2 behavior, and is not
a cause for dividing by rank count in itself:
https://getfem.org/userdoc/parallel.html .

Nevertheless, **RHS-only and simultaneous RHS+matrix have
distinct internal code paths**. In the published GetFEM C++ source
(e.g. GetFEM 5.5 `src/getfem_models.cc`), `BUILD_RHS`
has its own assembly branch, while `BUILD_ALL` goes through
a combined branch. This is a legitimate possible source of
different MPI behavior in a specific GetFEM build. This
source inspection is NOT proof about the user's local 5.4.2
binary: we must run the comparison on that very installation.

**New controlled experiment:** `sandbox/check_getfem_mpi_rhs.py`
now runs both `build_rhs` (SCRIMP's actual option) and
`build_all` on independent identical 1D GetFEM models
for `u*Test_u`, `u^2*Test_u`, `u^3*Test_u` at
`u=1` on [0,1], without any TS/SNES/KSP/MUMPS. It reports
`assembly_option`, raw GetFEM RHS directional derivative
and pre-PETSc tangent, and their ratio. It also hashes
the entire raw GetFEM CSR sparse matrix on each rank,
not merely its action on the constant vector.

Run in the original `scrimp` environment:
```bash
git pull --ff-only
conda activate scrimp
for n in 1 2 4; do
  timeout -k 5s 120s mpiexec -n "$n" \
    python -u sandbox/check_getfem_mpi_rhs.py \
    > "getfem_build_modes_${n}rank.log" 2>&1
  echo "MPI $n: exit=$?"
  grep -E '"expression"|"assembly_option"|"rhs_to_tangent_integral_ratio"|"raw_tangent_full_csr_digest"' \
    "getfem_build_modes_${n}rank.log"
done
```

**Outcome interpretation:**
- `build_rhs` gives ratio P but `build_all` gives ratio 1:
  an actual RHS-only build-path defect. Investigate switching
  SCRIMP to `build_all` only after verifying source/constraint
  bricks, linear terms, and performance (unneeded tangents).
- BOTH modes yield ratio P and full raw CSR matrices are identical
  per rank: stronger evidence that complete generic contributions
  are being assembled per process (or previously globally
  accumulated) before RHS summation; no wrong user-level
  assembly option fixes it.
- BOTH modes yield ratio 1: initial diagnosis depends on
  context or another code path; re-check the previous test.
- Full CSR rank digests DIFFER: GetFEM returns distributed
  matrix contributions despite identical constant-vector actions;
  examine partitioned regions and reduction separately.

**Important:** The current SCRIMP `auto` RHS calibration
is a compatibility correction, not proof of an upstream bug
nor a fully generic fix for explicit RHS and constraints.
Do not remove it or replace it with a single unconditional
division without resolving these model-specific semantics.

## 2026-11-13 — Scope decision: defer GetFEM MPI root cause to #23

The user requested to **prioritize completion of 2D Saint-Venant now**.
A separate tracking issue **#23** contains the exact 1D MPI
reproducer (1/2/4 ranks), fully identical per-rank raw CSR tangents,
`build_rhs` versus `build_all` comparison, official GetFEM
MPI partitioning documentation and community research. The GetFEM
upstream partition/build root cause is **not a release gate for the
present Saint-Venant study** provided the calibrated normalization is
validated for its actual bricks and the full physical/temporal
acceptance checks pass.

**For this exact tested local GetFEM 5.4.2 build**, use the explicit
`SCRIMP_GETFEM_RHS_MODE=divide_by_mpi_size` override during the
Saint-Venant smoke and subsequent runs. This is an intentional
workaround for issue #23, *not a global assertion that all replicated
GetFEM RHS vectors must be divided*. The generic default remains
`auto` (calibrates against an exact integral).

### Next: run actual corrected Saint-Venant time integrator

```bash
git pull --ff-only
conda activate scrimp

for n in 1 2; do
    echo "====== SWE MPI ${n} ======"
    SCRIMP_GETFEM_RHS_MODE=divide_by_mpi_size \
    PYTHONFAULTHANDLER=1 \
    timeout -k 5s 120s mpiexec -n "$n" \
      python -u sandbox/viscous_swe_2d_mpi.py \
        --case open --mesh 0.2 \
        --final-time 0.05 --save-dt 0.01 --max-step 0.002 \
        --output-dir "sandbox/outputs/swe_corrected_${n}rank" \
      > "swe_corrected_${n}rank.log" 2>&1
    echo "MPI=${n} exit=$?"
    tail -n 60 "swe_corrected_${n}rank.log"
done
```

**Acceptance criteria**: no SIGSEGV or native heap failure; actually
reach `t_f=0.05` (the script explicitly raises otherwise);
consistent state/Hamiltonian and dissipation/port-work balance;
1-versus-2 MPI parity. Zero exit code **alone** is not sufficient,
and an early TS stop is a failure. Do not merge into master, tag a
release, or close #6 until these criteria and more demanding full
runs are met. If the first run still fails, use its log to
identify whether the failure is a Newton divergence, PETSc TS
failure, or native corruption. **Do not reopen #23** unless the
normalization is demonstrably wrong for these real SWE bricks.


## 2026-11-13 — First corrected SWE temporal run: still two independent failures

The user ran original `scrimp` + GetFEM with the explicitly verified
`SCRIMP_GETFEM_RHS_MODE=divide_by_mpi_size`, `--case open
--mesh 0.2 --final-time 0.05 --save-dt 0.01 --max-step 0.002`:

**MPI 1:** `exit=1`, but the immediate Python traceback was a
**test-driver mistake**, not the reason for integrator failure:
`if swe.solution["t"]` after `DPHS.solve()` was ambiguous because
`self.solution["t"]` has been converted from a list to an ndarray.
The benchmark now reads `np.asarray(...).size` and reports
early convergence as a clear `RuntimeError`. The physical TS
already reported an automatic stop at **t=1e-8**, dt minimum
**1e-8**, **six rejected steps**, **zero SNES failures**,
1 accepted step, seven Newton iterations; it never approached
the requested t_f=0.05. That is consistent with adaptive local
error control rejecting steps, rather than a proven nonlinear
Newton failure. No energy/balance diagnostics are valid.

**MPI 2:** `exit=139`, SIGSEGV after approximately 14 further
GetFEM assembly calls during the two-step `pseudo` TS initial
consistency procedure. Full native crash persists with the
mathematically consistent residual/Jacobian and P-normalized RHS;
it is not disproved by successful isolated MUMPS replays.

**Additional code inspection revealed two potentially relevant
problems in this branch:**
1. The experimental nonlinear mass Hessian FD method used
   `z±epsilon*zd` with a fixed epsilon=1e-4. During
   pseudo-initialization, dt=1e-9 and `zd` can become enormous,
   causing off-trajectory, possibly unphysical perturbations of
   the fluid height. The method now bounds the norm of the
   Hessian test displacement to
   `epsilon*max(1,||z||)`, while leaving the already validated
   moderate-`zd` FD cases unchanged; with
   `SCRIMP_MPI_TRACE=1` it reports `z_norm`,
   `zd_norm`, effective `eta` and `perturb_norm`.
   This is numerical safeguard, NOT an established crash fix.
2. In `DPHS.init_step()`, the temporary removal of the
   `ts_adapt_dt_min` option mistakenly deleted
   `ts_adapt_dt_max`; corrected the typo. The baseline
   pseudo time integration still uses a tiny dt and must be
   tested, not assumed to construct DAE-consistent data.

The benchmark now accepts **`--ts-type {arkimex,cn,bdf}`**,
**`--adapt-type {basic,none}`** and **`--petsc-monitor`**.
The original ARKIMEX-A2 / basic adaptive scheme is still the
default. TSAdapt **none** disables embedded-estimator-based
adaptivity, giving a diagnostic fixed-step comparison:
https://petsc.org/release/manualpages/TS/TSAdaptType/
`--no-init-step` was already supported as a diagnostic to
skip the crashing pseudo-TS. It does **not** guarantee a
consistent index-2 DAE initial state and is not an accepted
production simulation on its own.

### Minimal next experiment: separate adaptive error control and pseudo initialization

First, run **one rank**, with the usual initial pseudo-TS but
no main adaptive step rejection; keep the default ARKIMEX-A2:

```bash
git pull --ff-only
conda activate scrimp

SCRIMP_GETFEM_RHS_MODE=divide_by_mpi_size \
SCRIMP_MPI_TRACE=1 PYTHONFAULTHANDLER=1 \
timeout -k 5s 120s mpiexec -n 1 \
  python -u sandbox/viscous_swe_2d_mpi.py \
  --case open --mesh 0.2 --final-time 0.006 \
  --save-dt 0.002 --max-step 0.002 \
  --ts-type arkimex --adapt-type none --petsc-monitor \
  --output-dir sandbox/outputs/swe_fixed_1rank \
  > swe_fixed_1rank.log 2>&1
echo "serial fixed-dt exit=$?"
tail -n 80 swe_fixed_1rank.log
```

Second, try **two ranks** with the same integrator and fixed dt,
but skip the pseudo-TS initialization to isolate whether the
native fault is confined to that preliminary path:

```bash
SCRIMP_GETFEM_RHS_MODE=divide_by_mpi_size \
SCRIMP_MPI_TRACE=1 PYTHONFAULTHANDLER=1 \
timeout -k 5s 120s mpiexec -n 2 \
  python -u sandbox/viscous_swe_2d_mpi.py \
  --case open --mesh 0.2 --final-time 0.006 \
  --save-dt 0.002 --max-step 0.002 \
  --ts-type arkimex --adapt-type none --no-init-step \
  --petsc-monitor --output-dir sandbox/outputs/swe_fixed_2ranks \
  > swe_fixed_2ranks.log 2>&1
echo "two-rank bypass pseudo exit=$?"
tail -n 95 swe_fixed_2ranks.log
```

**Interpretation:** if one-rank fixed dt reaches 0.006
but the adaptive run rejected steps down to 1e-8, we have
isolated a problem involving adaptation/error estimator
or initial transient; it does NOT prove physical accuracy.
If two ranks avoid the SIGSEGV only with
`--no-init-step`, the native crash is strongly linked to
the pseudo-TS lifecycle; that still leaves the main TS
physics and MPI stability to be tested.
`SCRIMP_MPI_TRACE=1` reports the exact last callback
and TS convergence/termination reason.

Only when **correct DAE initial conditions, adaptive convergence,
MPI parity, native reliability and complete physical
Hamiltonian/power balance** have passed may the model be
considered validated. No merge or release yet.


## 2026-11-13 — PETSc error 56 in the first fixed-TS diagnostic is a test-driver bug

The user ran both planned `--adapt-type none --petsc-monitor`
experiments. **Both immediately failed with PETSc error 56** at
`TS.setFromOptions()`, **before either TS.solve()**:

```text
TSSetFromOptions -> TSMonitorSetFromOptions
 -> PetscOptionsGetViewer -> Unsupported viewer true
```

This is caused by **our diagnostic flag wiring**, not a Saint-Venant
Newton failure, MPI memory corruption, or a MUMPS factorization.
`--petsc-monitor` previously called
`swe.set_time_scheme(ts_monitor=True, ts_adapt_monitor=True,
snes_converged_reason=True, ksp_converged_reason=True)`.
SCRIMP's PETSc Options database stringifies `True` to
`"true"`, and PETSc interprets it as a VIEWER name rather than
as an argument-less monitor flag.

**Fixed** the script to call
`PETSc.Options().setValue("ts_monitor", None)`, and likewise for
the other three flags; PETSc receives options with no viewer
arguments. This correction needs a new native user run to verify
the user's PETSc 3.21.6 build, but matches the documented
options API. `--petsc-monitor` stays optional; omit it entirely
if any environment-specific monitor issue persists, since
`SCRIMP_MPI_TRACE=1` already traces critical solver milestones.

The previous runs with code 1 do **NOT** constitute valid tests
of fixed-step ARKIMEX, nor of `--no-init-step` on two ranks.
Retry the two commands documented directly above using the
new branch head, and do not infer anything about the native
SIGSEGV until a real nonlinear TS.solve() has occurred.


## 2026-11-13 — Replace Saint-Venant pseudo-TS initial solve with tiny BDF2 from negative time

**User request:** the nonlinear Saint-Venant mass matrix depends on
the state; test whether the disastrous pseudo-TS initial consistency
step is purely a numerical initialization artifact. Explicitly use
BDF2 with a *very small* step, starting at `-dt` instead of `0`.

**Reason for the test:** On a single MPI rank, a fixed-step
ARKIMEX run reached `t_f=0.006` in three steps without PETSc
step rejections or reported SNES failures, but the initial
pseudo-TS yielded `||z|| ~= 3.90477e24`; the power-balance defect
was `561036.152228463`, so that *successful exit code 0 was not
a valid physical solution*. On 2 ranks, bypassing pseudo-TS
reached the first main TS Jacobian, but then MUMPS reported
`Instance Error 1 in ZMUMPS_F77` / MPI_Abort(-99), exit 157.
That is a **different native failure** from the previous SIGSEGV,
and BDF2 is not proven to fix it.

**PETSc 3.21 technical detail:** BDF order 2 has an internal
history bootstrap, performing a first solve at half of the first
time step before the full step. Source:
https://web.cels.anl.gov/projects/petsc/vault/petsc-3.21/docs/src/ts/impls/bdf/bdf.c.html
(`TSBDF_Restart`). This permits *one accepted macro step*
from `t0-dt` to `t0` with `ts_bdf_order=2`, without falsely
calling a single ordinary backward-Euler step BDF2.

**Implementation in this draft branch:**
- `DPHS.init_step()` now accepts opt-in
  `init_step_ts_type="bdf"`,
  `init_step_bdf_order=2`,
  `init_step_from_negative_time=True`.
  It starts PETSc at `t0-init_step_dt`, sets target time `t0`,
  performs precisely **one fixed macro BDF2 step**, with
  `ts_adapt_type=none` (only for the initial TS).
- Restores `ts_bdf_order`, `ts_adapt_type`,
  `ts_adapt_dt_min` and `ts_adapt_dt_max` for the main TS
  after initial TS cleanup. No stale BDF settings leak into
  main ARKIMEX. Historical pseudo startup remains the generic
  SCRIMP default for unrelated models.
- `sandbox/viscous_swe_2d_mpi.py` now uses the above BDF2
  initialization by default, `dt=1e-9`, start `-1e-9`,
  target `0`, one step. Keep existing main ARKIMEX-A2
  and `--adapt-type` choice unchanged.
- Added explicit checks: positive TS convergence reason,
  final initial time numerically t0, and finite norm with
  `||z_final|| <= 100*max(1,||z_initial||)`. The norm guard
  is **diagnostic and model-specific**; it catches previous
  `1e24` blow-ups but is not proof of DAE consistency or
  physical height positivity.
- `SCRIMP_INIT_PC_TYPE` alternative remains available
  for isolating native MPI LU errors. The default uses MUMPS.

**Next two short tests (original `scrimp` environment only):**

```bash
git pull --ff-only
conda activate scrimp

for n in 1 2; do
  SCRIMP_GETFEM_RHS_MODE=divide_by_mpi_size \
  SCRIMP_MPI_TRACE=1 PYTHONFAULTHANDLER=1 \
  timeout -k 5s 120s mpiexec -n "$n" \
    python -u sandbox/viscous_swe_2d_mpi.py \
      --case open --mesh 0.2 \
      --final-time 0.006 --save-dt 0.002 --max-step 0.002 \
      --ts-type arkimex --adapt-type none \
      --output-dir "sandbox/outputs/swe_bdf2init_${n}rank" \
    > "swe_bdf2init_${n}rank.log" 2>&1
  echo "MPI=$n exit=$?"
  grep -E 'Initial consistency:|init_step: consistency TS completed|main TS completed|MUMPS|Instance Error|Initial BDF2|RuntimeError|maximum_relative_balance_error' \
    "swe_bdf2init_${n}rank.log" | tail -n 25
  tail -n 25 "swe_bdf2init_${n}rank.log"
done
```

**Interpretation gates:** If the initialization converges and produces
a bounded state on one rank, re-evaluate energy, dissipation and
boundary supply: **the relative power-balance defect should become
small on a truly valid solution**, and all quantities should remain
finite. Verify whether the 2-rank MUMPS Instance Error persists:
this failure occurred **even without any pseudo-TS**, so the new
initializer cannot be presumed to resolve it. If BDF2 fails to
initialize, do NOT simply increase the norm guard or force a
success. Analyse algebraic constraints, `h*p` temporal brick
semantics and physical consistency of the initial conditions.

This change remains only on feature branch PR #22. No master merge,
issue #6 closure, or release tag pending physical and MPI validation.


## 2026-11-13 — Negative-time BDF2 fails too; inspect FIRST runaway Newton block

The user's `swe_bdf2init_1rank.log` and
`swe_bdf2init_2rank.log` show the requested fixed
`dt=1e-9` BDF2 startup from `t=-1e-9` to `t=0`
does **not** reach t=0.

- **MPI 1:** `exit=1`; PETSc
  `TSBDF_Restart -> TSBDF_SNESSolve -> MatSolve_MUMPS`
  gives `INFOG(1)=-3` (PETSc error code 76).
  **MUMPS -3 denotes an invalid JOB/analysis-factor-solve
  sequence**, not the singular-matrix numerical error code.
  MUMPS manual: https://mumps-solver.org/doc/userguide_5.8.1.pdf .
  This happens only *after* an extremely large Newton step
  and multiple BDF startup retries.
- **MPI 2:** `exit=139`; native SIGSEGV after the second
  assembled 3300x3300 Jacobian in the BDF startup.
- On BOTH 1/2 ranks, after one Newton correction the trace
  reports **`||z||=1.31957e6`** and
  **`||zd||=2.63914e15`** (with BDF `shift=2e9`).
  The nonlinear stiffness magnitude jumps from
  `2058.69703` to **`533394.687`** at the
  same point. The two-rank result agrees exactly with
  serial up through this divergence; this phenomenon
  starts BEFORE an MPI-specific error.
- Even though TSAdapt was set to `none`, the failed
  initial nonlinear substep is internally retried with
  smaller intervals: the serial log shows BDF
  times `-5e-10`, then `-8.75e-10`, then
  `-9.6875e-10`; Jacobian shifts `2e9`, `8e9`,
  `3.2e10`. Therefore the startup is NOT successfully
  maintaining a constant step; an internal solver failure
  is triggering retries. Reducing the already-tiny
  startup timestep is not a demonstrated remedy.
- The initial matrices have 132 zero diagonal entries
  (but **no zero numerical rows**) on rank 0: these
  are not sufficient to diagnose singularity; an
  LU/MUMPS *solve* failure `-3` does not mean
  matrix singularity.

**Next observation, not another speculative numerical fix:**
`DPHS._trace_newton_blocks()` is enabled only by
`--newton-block-trace` (or
`SCRIMP_NEWTON_BLOCK_TRACE=1`). It logs full replicated
trial state norm, each GetFEM variable's norm/max absolute
value and, for fluid height `h`, min/max. It adds
NO MPI collectives or extra assembly. Use `--init-dt`
to vary negative-time BDF startup step only if needed.

Run first on **one rank only**:

```bash
git pull --ff-only
conda activate scrimp

SCRIMP_GETFEM_RHS_MODE=divide_by_mpi_size \
PYTHONFAULTHANDLER=1 \
timeout -k 5s 45s mpiexec -n 1 \
  python -u sandbox/viscous_swe_2d_mpi.py \
    --case open --mesh 0.2 \
    --final-time 0.006 --save-dt 0.002 --max-step 0.002 \
    --ts-type arkimex --adapt-type none \
    --newton-block-trace \
    --output-dir sandbox/outputs/swe_bdf2_blocks_1rank \
  > swe_bdf2_blocks_1rank.log 2>&1

echo "exit=$?"
grep -E 'Newton blocks|IJacobian: enter|MUMPS error|Initial consistency' \
  swe_bdf2_blocks_1rank.log | head -n 35
tail -n 22 swe_bdf2_blocks_1rank.log
```

The decisive datum is **which GetFEM block** jumps from
initial values (~`h=5`, tangential `p=0.2`) to
~1e6, and whether `h` becomes negative. This
distinguishes problems in the initial algebraic
constitutive/boundary constraints from dynamics and
incorrectly defined state-dependent mass. If Newton's
first step drives `h` negative or port variables explode,
solve consistent algebraic initial data explicitly
rather than forcing an ever smaller time step.

**Do not merge PR #22** or claim physically correct
energy balance. The native MPI MUMPS crash remains open;
GetFEM MPI RHS root cause remains deferred to #23.
