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
