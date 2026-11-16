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

**Completed six-cycle tests (2026-11-16):** all eight combinations
(1, 2, 4, 6 MPI ranks x `always` / `on_change`) return exit code 0.
Each performs six accurate solves and six numeric factorizations. `always`
performs six symbolic factorizations; `on_change` performs **two**, with
one reset at cycle 3 on all ranks. Maximum external algebraic residual
and full-state relative error are below 3.6e-16 for every configuration.

**Important coverage limit:** the cycle-5 source stiffness omits the
second-neighbour connections, but `MatCopy` into the persistent Pmat
*retains the already allocated CSR positions*. The observed stored CSR
change flag at cycle 5 is `false` on every rank: this case tests
numerically zeroed couplings, **not deletion of stored nonzeros**.
A genuine shrinking-pattern case remains untested. Do not interpret
cycle 5 as proof of safe symbolic reuse after structural removal.

The standalone test defaults to `on_change`, because `native` has already
crashed on the cycle-3 structural growth in this PETSc 3.26.1 build.
No production `dphs.py` policy change until GetFEM nonlinear verification.
Next prioritize the full acoustic state parity and power balance, then
real GetFEM nonlinear Jacobians with a deliberately changing pattern.

**Four-rank full acoustic trajectory parity confirmed (2026-11-16):**
the real `examples/wave.py` model with 12,682 DOFs, Crank–Nicolson
`dt=0.002` and `t_f=0.04` completed 20 accepted steps, 0 rejected
steps, 0 SNES failures, and 20 saved snapshots under both
`SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN=1` and `=0`. Both runs gave
zero rank-to-rank state differences. Full saved coefficient trajectories
agree: maximum absolute difference 2.8765878568e-13, relative L2
trajectory difference 6.0266261847e-14. Single-run reported solver
elapsed times were 1.816 s (reset) and 0.5967 s (reuse), an
approximately 3.04x indicative speedup; repeat runs before quoting a
statistically reliable performance improvement.

**Warning to investigate separately:** both runs issue
`ComplexWarning: Casting complex values to real discards the imaginary
part` in `scrimp/dphs.py:1046` during the initial consistency solve.
The trajectories agree in spite of this warning, but this test alone
does not establish that any imaginary data would be negligible in
general. Do not silently suppress the warning.

**Driven boundary-power audit completed at 4 MPI ranks (2026-11-16):**
the acoustic 2D problem with 2,230 DOFs and boundary control
`U_L=0.08*sin(20*t)` completed 10 Crank–Nicolson steps at `dt=0.002`
for both `SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN=1` and `=0` (exit 0, no
rejected steps, no SNES failures). The measured energy change was
`+3.352923209074335e-5` and integrated port work was
`-3.3529232090756814e-5`. The maximum absolute cumulative
power-balance defect was `1.8730439562634343e-17` (relative to H:
`1.3097473795212481e-15`); the maximum per-step defect was
`1.7550522667109103e-17`. In the uploaded run logs, **all 23
fields** of the JSON audit were *exactly equal*, including all saved
Hamiltonian samples, the separate energy terms, all four boundary
power series, every time-integrated port-work sample, and every
stepwise power-balance defect. The boundary work is appreciably
nonzero, so this is not a vacuous closed-boundary test.

Both runs emit the already documented real/complex cast warning
in initial consistency. They additionally print PETSc
'unused database options' messages for the script's argparse flags;
the JSON confirms the intended parameter values were applied and
this message is not a failed test.

**Nonlinear viscous Saint-Venant baseline (2026-11-16; not yet
converged):** 1 MPI rank, disk mesh=0.2, open case, rigid initial
rotation=0.02, `t_f=0.01`, fixed ARKIMEX step=0.005, and
`SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN=1`. The initial BDF2 consistency
step from `-1e-9` to `0` failed SCRIMP's state-amplification guard:
`||z_initial||=110.000484`, `||z_final||=5998040.66`, maximum
amplification 100. This is **not** evidence of a symbolic-reuse bug:
the conservative reset policy was enabled and the main ARKIMEX
integration never began. The run nonetheless generated two real
GetFEM/PETSc Jacobian snapshots,
`swe_jac_001.petscbin` and `swe_jac_002.petscbin`.

Use `replay_petsc_mumps_snapshot.py --inspect-only` to compare the
actual stored CSR patterns and matrix values of these two snapshots
*without factorization*. This new diagnostic avoids the old benchmark's
cycle-5 ambiguity from retained MatCopy sparsity and does not assume
that a failed nonlinear initialization produces well-conditioned LU
test matrices. Only once this inspection is understood should the
existing snapshot replay with explicit PC reset be considered.

**Independent nonlinear Jacobian snapshot inspection passed
(2026-11-16, 1 MPI rank):** snapshots `swe_jac_001` and `swe_jac_002`
are both 3300 x 3300, with 123,765 and 137,114 stored nonzeros,
respectively. Matrix 002 **adds exactly 13,349 distinct stored
locations**, removes none, and changes the stored pattern in 1,169
rows; this is genuine sparsity growth in *actual GetFEM* Jacobians,
not mere numerical zeros or lingering PETSc MatCopy slots. Both
snapshots have zero nonfinite entries; their stored numerical zeros
number 132 and 167. Their Frobenius norms are about 1.5133e9,
with a Frobenius difference of 1.082386e6 (relative difference
7.152434e-4). These matrices come from an **aborted, divergent BDF2
initialization** and are *not* evidence of a physically valid
Saint-Venant simulation; the extremely large matrix norm also
limits numerical interpretation.

**Real GetFEM MUMPS factorization/replay checks passed on one MPI rank
(2026-11-16):** snapshot 001, separately factored and solved, gave
relative residual `2.29575684e-16`; snapshot 002 gave
`1.99743508e-16`. With one persistent KSP/PC/Mat across 001 -> 002,
an explicit `PC.reset()` before the in-place MatCopy gave the same
two residuals, and `KSP.reset()` / `KSP.destroy()` completed normally.
All three runs exit 0. **Limit:** the original replay script checked
external residuals but did not yet compare the recovered state with
the manufactured exact solution `x=1`. These large-norm Jacobians
may be ill-conditioned, so tiny residuals alone do not prove a
small forward error. The replay script now reports and gates both
`relative_linear_residual` and `relative_state_error` at 1e-5.
Re-run the 001 -> 002 transition with this state-error check and
consider 4 cycles (001 -> 002 -> 001 -> 002) under explicit reset.
Do not run the unsafe `--transition-reset none`: synthetic changing-
pattern reuse already caused a SIGSEGV. Independently factoring
an already captured Jacobian does **not** validate the nonlinear
BDF2 step or the model's Hamiltonian balance.

**Four-solve real-GetFEM persistent-factor reference passed
(2026-11-16, 1 MPI rank):** the same KSP/PC/Mat was used for the
001 -> 002 -> 001 -> 002 sequence, with explicit `PC.reset()`
before each MatCopy transition. External relative residuals were
2.29575684e-16, 1.99743508e-16, 2.43279796e-16,
1.99743508e-16; the new manufactured-solution relative state
errors were **4.19515974e-7, 3.53538232e-7, 4.12679850e-7,
3.53538232e-7**, respectively. All were under the declared
1e-5 acceptance threshold, and full reset/destroy succeeded. This
is NOT 1e-16 accuracy in the unknowns: the state error indicates
nontrivial sensitivity, consistent with possible conditioning/scaling
issues in these matrices obtained during a divergent BDF2 solve.

**New targeted PETSc-only diagnostic:** `--transition-reset on_change`
compares the exact stored CSR on all ranks *after* in-place MatCopy and
before the next KSPSolve. It performs a collective `PC.reset()` and
`KSP.setOperators()` only if that stored pattern changed; it emits
`stored_pattern_changed` and `explicit_factor_resets`. This variant
has not yet been run on the captured nonlinear pair. Test it with
`--cycles 4` at 1 rank against the existing `pc` reference, then at
4 ranks (reloading/repartitioning the PETSc binary snapshots).
If the 4-rank load fails, first diagnose the binary MatLoad MPI
compatibility rather than altering SCRIMP.

**Real GetFEM conditional-reuse replay passed in serial and MPI
(2026-11-16):** all four cases (1 and 4 ranks, `pc` versus
`on_change`) completed four persistent-KSP solves and a clean
reset/destroy. Only the initial 001 -> 002 transition grows the
*stored* CSR; the later 002 -> 001 -> 002 transitions retain its
superset pattern. The `pc` reference resets 3 times, whereas
`on_change` resets just once **collectively on all ranks**.
In serial, both modes have maximum relative residual
2.43279796e-16 and maximum relative state error 4.19515974e-7.
With four ranks, the largest observed residual is 3.8289517e-16
(`on_change`) and the largest state error 1.90574948e-6
(both modes). The `pc` and `on_change` four-rank state errors
are close but not bitwise identical; this is consistent with
roundoff sensitivity and is not itself proof of a condition number.
All states passed the exploratory 1e-5 forward-error gate.

**What this does not prove:** these logs do not include PETSc
`MatLUFactorSym` event counts, so 3 versus 1 explicit PC resets
cannot yet be translated directly into symbolic-factorization
counts or an actual runtime speedup. The exact CSR comparison has
a per-Jacobian cost and an MPI collective. Furthermore, these
are replayed matrices from a BDF2 initialization that failed its
state-amplification safety gate, not a converged nonlinear
time-domain simulation. Keep the existing `dphs.py` behavior
and explicit opt-in flags unchanged.

**Next priority:** run the current 1-rank viscous Saint-Venant
baseline with `--newton-block-trace`, the same safe explicit
MUMPS reset, and the existing BDF2 amplification guard still
enabled. Review the per-GetFEM-block Newton norms and the
height minima/maxima to locate the source of initialization
divergence before retrying alternative time steps or solvers.

**Separately:** validate a converged nonlinear reference before
proposing any `dphs.py` policy change.


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
    for mode in always on_change; do
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
nonzero, and cycles 1/2/4/5 retain their stored patterns. Compare
`on_change` and `always` PETSc event counts. Both perform six numerical
LU factorizations/solves; `on_change` performs two symbolic factorizations
on the observed initial and expanded patterns. Cycle 5 omits numerical
second-neighbour couplings, but does **not** shrink stored CSR.

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
