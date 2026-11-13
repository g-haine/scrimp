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

