# PETSc/MUMPS symbolic reuse validation (SCRIMP v1.3.0)

Status: **experimental tests; no production solver-policy change yet**.

## Reproduced changing-sparsity failure (2026-11-16)

In the standalone PETSc 3.26.1 complex-scalar MUMPS probe, cycles 0-2
preserve the CSR sparsity pattern and `native` succeeds for both 1 and 4 MPI
ranks, with 1 symbolic LU, 2 numeric LU and 2 solves in the initial
two-cycle run. With the second-neighbour couplings added at cycle 3,
the **unconditional native-reuse policy is unsafe in this test**:

- 1 rank, `native`, 4 cycles: KSP reports convergence, but the external
  relative residual is ~7.22 and the relative state error ~6.02.
- 4 ranks, `native`, 4 cycles: SIGSEGV in `KSPSolve`.
- 1 and 4 ranks, `always`, 4 cycles: all residual and state errors remain
  around 1e-16, with 4 symbolic and 4 numeric factorizations.
- The earlier `always` setup error (PETSc error 75) was a *probe defect*:
  after `PC.reset()` the standalone KSP must call
  `ksp.setOperators(operator)` to reattach the operator.

The exact-CSR diagnostic mode `--mode on_change` now resets PC and
rebinds KSP only on a collective change of *stored* sparsity. This is a
**hypothesis to validate**, not a production default or evidence that
the PETSc nonzero-state bookkeeping itself is necessarily defective.

The original no-reset vs always-reset acoustic benchmark used a
constant sparsity pattern, so its measured speedup remains valid for
that case. Do not extrapolate it to time-varying nonlinear sparsity.

**Confirmed four-cycle `on_change` run (1 and 4 MPI ranks):** both exit
successfully, `PC.reset()` executes once at cycle 3 (and on every rank in
MPI), exactly 2 symbolic LU / 4 numeric LU / 4 linear solves occur, and
all four full-state / residual errors are below 3e-16. The four-rank
`always` reference executes 4 symbolic LU / 4 numeric LU / 4 solves and
agrees to machine precision. These are tiny 240-by-240 probe matrices:
do not infer a meaningful wall-time speedup from these timings.

**Next check:** run all 6 cycles for `always` and `on_change` at 1, 2, 4,
and 6 ranks; specifically verify that the cycle-5 structural *removal*
(if it really changes the stored CSR pattern) is detected, the state
checks pass, and the solver shuts down normally. No production
`dphs.py` policy change until GetFEM nonlinear verification.


## Rationale

SCRIMP's existing `DPHS.IJacobian()` builds the preconditioner matrix by
`MatCopy(..., DIFFERENT_NONZERO_PATTERN)` followed by
`MatAXPY(..., DIFFERENT_NONZERO_PATTERN)`. These calls correctly refuse
to assume equality of the source and destination sparsity structures.

In PETSc 3.26.1, `PCSetUp()` compares the current nonzero state of its
preconditioning matrix with the previous state. `PCLU` redoes symbolic
factorization if that state changed; otherwise it can do a new numerical
factorization without repeating the symbolic phase. It is important not
to use `PC.setReusePreconditioner(True)`, because the **numeric**
factorization must track the current Jacobian.

The acoustic 4-rank single-run comparison gives 21 to 2 symbolic
factorizations and 2.0934 s to 0.8495 s for `DPHS.solve()` without
unconditional `PC.reset()`. That is **not** sufficient evidence to
remove the conservative workaround on nonlinear Saint-Venant.

Reference PETSc code:
- `PCSetUp()`: https://petsc.org/release/src/ksp/pc/interface/precon.c.html
- `PCSetUp_LU()`: https://petsc.org/release/src/ksp/pc/impls/factor/lu/lu.c.html
- `MatCopy()`: https://petsc.org/release/manualpages/Mat/MatCopy/
- `MatGetNonzeroState()`: https://petsc.org/release/manualpages/Mat/MatGetNonzeroState/

## 1. Standalone PETSc/MUMPS structural-change probe

Run from the repository root on the existing `scrimp-mpi:local` image.
These functional checks disable OpenMPI affinity to avoid the container
hwloc memory-binding warning; do not use their timings for scaling claims.
No native library needs to be rebuilt.

```bash
set -euo pipefail
out="$PWD/sandbox/outputs/symbolic-reuse"
mkdir -p "$out"

for n in 1 2 4 6; do
    case "$n" in
        1) cpus=0 ;;
        2) cpus=0,2 ;;
        4) cpus=0,2,4,6 ;;
        6) cpus=0,2,4,6,8,10 ;;
    esac
    for mode in always native; do
        sudo docker run --rm --cpuset-cpus="$cpus" \
            -v "$PWD/sandbox:/workspace/scrimp/sandbox:ro" \
            scrimp-mpi:local \
            mpiexec --bind-to none -n "$n" python -u \
            sandbox/check_petsc_mumps_symbolic_reuse.py \
            --mode "$mode" --petsc-log \
            > "$out/petsc_${n}_${mode}.log" 2>&1

        grep -aE '^MatLUFactorSym|^MatLUFactorNum|^MatSolve ' \
            "$out/petsc_${n}_${mode}.log"
    done
done
```

**Gates:** all six linear systems converge in each run, relative
state/residual errors < 1e-10, cycle 3 truly introduces a new stored
nonzero, and cycles 1/2/4 retain their patterns. Compare native and
always-reset PETSc event counts. Both execute six numerical LU solves;
the native branch should avoid symbolic reanalysis for fixed patterns.
Cycle 5 intentionally removes the extra coupling, but the test does
not require any particular symbolic event count for this transition.

## 2. Full acoustic state parity (reset versus native reuse)

`check_acoustic_state_parity.py` reconstructs exactly `examples/wave.py`
with the same mesh and time integrator as `bench_acoustic_2d.py`.
It stores every snapshot and *all* GetFEM state coefficients in NumPy
NPZ format (complex values preserved), compares full arrays rather than
just norms, and checks MPI rank-to-rank agreement after integration.

The container user needs write access to the host artifact directory
(no change to SCRIMP, conda, or Docker image).

```bash
set -euo pipefail
out="$PWD/sandbox/outputs/symbolic-reuse"
mkdir -p "$out"
chmod a+rwx "$out"

for n in 1 2 4 6; do
    case "$n" in
        1) cpus=0 ;;
        2) cpus=0,2 ;;
        4) cpus=0,2,4,6 ;;
        6) cpus=0,2,4,6,8,10 ;;
    esac
    for reset in 1 0; do
        reference=()
        label=reset
        if [ "$reset" = 0 ]; then
            label=native
            reference=(--reference "/artifacts/acoustic_${n}_reset.npz")
        fi
        sudo docker run --rm --cpuset-cpus="$cpus" \
            -e SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN="$reset" \
            -v "$PWD/scrimp:/workspace/scrimp/scrimp:ro" \
            -v "$PWD/sandbox:/workspace/scrimp/sandbox:ro" \
            -v "$PWD/examples:/workspace/scrimp/examples:ro" \
            -v "$out:/artifacts" \
            scrimp-mpi:local \
            mpiexec --bind-to none -n "$n" python -u \
            sandbox/check_acoustic_state_parity.py \
            --mesh 0.08 --dt 0.002 --final-time 0.04 \
            --output "/artifacts/acoustic_${n}_${label}.npz" \
            "${reference[@]}" \
            > "$out/acoustic_${n}_${label}.log" 2>&1
    done
done
```

**Gates:** complete state arrays agree within rtol=1e-10 and
atol=1e-11 at every saved time; 20 accepted CN steps, no rejections,
no SNES failures; final time exactly 0.04 within floating-point ULP.
The separately existing `check_acoustic_energy_mpi.py` must still be
run with **nonzero boundary power**, on 1/2/4/6 ranks and both reset
modes. Require the Hamiltonian + port-work identity, not merely
matching final solution norms.

## 3. Prerequisites for a SCRIMP production policy change

1. Verify the native PETSc structural test for all targeted MPI sizes.
2. Verify full acoustic trajectories and driven-port energy balances.
3. Reproduce a *real GetFEM* changing-sparsity nonlinear Jacobian and
   inspect PETSc symbolic/numeric counts and Newton residuals.
4. Validate nonlinear viscous Saint-Venant 2D and the factor lifecycle
   (repeated solve and teardown) with and without PC reset.
5. Only then integrate a documented default policy. Keep an explicit,
   safe fallback for configurations with unresolved native MUMPS defects.

Do **not** replace `DIFFERENT_NONZERO_PATTERN` with
`SAME_NONZERO_PATTERN` without a proof of structural equality.
Do **not** disable numeric factorization updates.
