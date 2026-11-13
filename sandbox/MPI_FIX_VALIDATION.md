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
