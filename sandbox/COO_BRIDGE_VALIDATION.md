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
