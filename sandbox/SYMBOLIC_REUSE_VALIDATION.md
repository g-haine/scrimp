# PETSc/MUMPS symbolic reuse validation (SCRIMP v1.3.0)

Status: **experimental tests; no production solver-policy change yet**.

## End-to-end MPI on_change/always nonlinear parity (2026-11-17)

The user reran the already-validated P3/P2 physical Saint-Venant
smooth spin-up at `mesh=0.2`, `rotation=0.02`, `dt=0.002`,
`t_f=0.02`, with the new **automatic** nonlinear policy
(`on_change`) and `always`. Both policies passed on 2 and 4
MPI ranks, with 10 accepted steps, zero rejected steps, zero
nonlinear solver failures, 20 nonlinear and 20 linear iterations.

| MPI ranks | always elapsed (s) | automatic elapsed (s) | Indicative time reduction |
| ---: | ---: | ---: | ---: |
| 2 | 2.943 | 2.751 | 6.52% |
| 4 | 2.334 | 2.236 | 4.20% |

At 4 ranks, the maximum relative volume balance errors agree
within floating-point noise (1.365270688239657e-15 for always,
1.36527068823965e-15 for auto). The maximum energy residual
normalized by actual exchange is 1.3790925633682102e-6 for
always and 1.3790925632450105e-6 for auto. At 2 ranks the
corresponding energy-exchange residuals are
1.3790279710293058e-6 and 1.379023934139024e-6.
These are **single runs**, not statistical speedups.
PETSc symbolic factor counts were not collected in these
four runs. These successful tests validate this particular
nonlinear MPI physical trajectory, not arbitrary changing
CSR patterns in every solver or PC type.

## Repeating the five-repetition 2D acoustic scaling chart

The original acoustic-wave strong-scaling experiment
(`examples/wave.py`, fixed mesh `h=0.08`,
12,682 DOFs, Crank–Nicolson `dt=0.002`,
`t_f=0.04`, ranks 1/2/4/6 on distinct physical cores)
is used to visualize the impact of the **new all-linear
automatic `never` policy**. The existing
`sandbox/bench_acoustic_2d.py` accepts
`--series auto` or `--series always` for tagging
each JSON report and validates the actual selected policy
after the solve (auto must resolve to `never`; always
to `always`). The existing
`sandbox/analyze_acoustic_scaling.py` accepts
`--plot-comparison PREFIX` and all tagged JSON
reports. It checks identical physical numerical results,
time-stepping parameters and per-point repeat counts
before producing `PREFIX_speedup.png` and
`PREFIX_times.png`. Each curve uses **median**
end-to-end `wall_solve_seconds_max` values.
`--profile` invokes the preexisting SCRIMP profiler;
no new timing infrastructure was introduced.

For a legitimate comparison, run both policies with
the same P-core affinity at each MPI rank, five
repetitions per (policy, rank), preferably interleaved.
The previous image is a historical reference, not a
measured outcome for the newly integrated policy.

## Production DPHS reset policy: automatic linear fast path (2026-11-17)

**Implemented in `scrimp/dphs.py` on this PR; awaiting local runtime
validation.** The exact CSR `on_change` logic, already validated in
standalone and real-Jacobian PETSc tests, is now connected directly
to `DPHS.IJacobian()`. Each new `DPHS.solve()` resolves one reset policy
collectively and clears its initial/main TS-specific CSR histories.
After full implicit Jacobian assembly (including the optional
state-dependent mass-Hessian correction), `on_change` compares
detached copies of local PETSc stored CSR row pointers, columns and
ownership ranges, then uses a collective MPI logical OR. If the
global graph changes, it resets the MUMPS PC and rebinds the existing
KSP operators before the next factorization; numerical-value-only
changes do not reset MUMPS. This is an experimental integration,
not a demonstrated new end-to-end runtime speedup.

**Final default policy, including linearity optimization:**

| Registered SCRIMP bricks | No environment override | CSR check in IJacobian? |
| --- | --- | --- |
| All Jacobian bricks marked `linear=True`; controls certified prescribed | `never` | No |
| Any state-dependent feedback or other nonlinear brick | `on_change` | Yes, for MUMPS LU |
| Source brick not certified state-independent | `on_change` | Yes, for MUMPS LU |
| Empty model | `on_change` | Yes, if assembled with MUMPS |

Both the `initial` and `main` TS use their own CSR snapshots,
since their KSP/PC instances have independent factor lifecycles.

**Closed-loop clarification (2026-11-17):** an acoustic smoke run
failed its benchmark policy gate after *successfully completing*
20 CN steps because `set_control()` always created a source brick
with `linear=False`, misleading the former all-brick predicate
into `on_change`. Simply excluding **all** source bricks would
be unsafe: `examples/heat_wave.py` calls
`set_control("Interface Heat", "Y_w")` and
`set_control("Interface Wave", "-Y_T")` to interconnect
two subsystems. The GetFEM `add_source_term_brick` is intended
for prescribed RHS, not a differentiable feedback law:
https://getfem.org/userdoc/model_source_term.html .

The fixed `DPHS.set_control()` now distinguishes:
- **Certified prescribed** expression using known time, geometry
  and mathematical functions: unchanged GetFEM RHS source, tagged
  `_scrimp_prescribed_control_source=True`; does not force
  nonlinear stiffness, and qualifies for `never` if other
  Jacobian bricks are linear.
- **State-dependent expression**, including references to other
  port variables such as `Y_w`, `Y_T`, `e_p`:
  generated as a `linear=False` GetFEM *nonlinear variational
  constitutive term* `(expression).Test_U` (using appropriate
  scalar/vector/tensor contraction), so GetFEM contributes
  the correct feedback derivative to the Newton Jacobian.
  It selects `on_change` automatically.
- **Unknown custom macros/data or user source bricks**: never
  silently certified state-independent. Unknown expressions are
  conservatively routed through the variational nonlinear branch
  (and `on_change`), at the possible cost of extra assembly.
  This identifier-level detector is deliberately NOT a general
  GetFEM macro expansion or a proof of state-independence for
  arbitrary user-defined macro definitions.

The automatic fast path still assumes that any brick declared
`linear=True` has a stable stored sparsity graph. The feedback
route and its sign/Jacobian must be runtime-checked on the actual
interconnected heat/wave example before claiming full closed-loop
regression coverage.

The explicit environment override `SCRIMP_MUMPS_RESET_POLICY`
accepts `on_change`, `always`, or `never` and **always wins**
over automatic linearity detection. Legacy
`SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN=1` means `always` and
`=0` means `never`, but is only read when the new variable
is absent. For `always`, `PC.reset()` is performed on each
Jacobian as the conservative fallback. For `never`,
no repeated CSR snapshot/comparison or MPI change collective
is performed.

**Minimum local validation before acceptance of the optimization:**

```sh
conda activate scrimp
python -u sandbox/check_mumps_reset_policy_unit.py

# Linear system: should select never automatically, and agree with
# explicit always and on_change. Use the existing acoustic profiler.
# Unset the old flag to actually exercise the automatic choice:
env -u SCRIMP_MUMPS_RESET_POLICY \\
    -u SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN \\
    python -u sandbox/bench_acoustic_2d.py \\
    --mesh 0.08 --dt 0.002 --final-time 0.04 \\
    --save-dt 0.01 --profile --petsc-log

# Nonlinear PDE: the following is the already-working P3/P2 open
# physical spin-up. Repeat at MPI 1/2/4, compare against
# SCRIMP_MUMPS_RESET_POLICY=always for state invariants and
# MatLUFactorSym/Num counts. Use the existing power/volume reports.
env -u SCRIMP_MUMPS_RESET_POLICY \\
    -u SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN \\
    mpiexec -n 2 python -u sandbox/viscous_swe_2d_mpi.py \\
    --case open --mesh 0.2 --height-order 3 --rotation 0.02 \\
    --initial-profile rigid --wall-velocity-profile facet_tangent \\
    --start-from-rest-duration 0.02 --init-method bdf \\
    --init-dt 1e-8 --final-time 0.02 --max-step 0.002 \\
    --save-dt 0.002 --adapt-type none \\
    --output-dir sandbox/outputs/symbolic-reuse/swe_on_change_mpi2
```

Required gates: no nonlinear or linear solver failures, exact
physical endpoint, reference-close full state/invariants and power
balance, no factor corruption on 1/2/4 ranks, and fewer symbolic
factorizations than `always` when the stored CSR graph is
unchanged. Validate repeat-solve/teardown, and only then report
speedups from the already-existing SCRIMP/PETSc profiler. Avoid
claiming the new integration was executed in the user's container
until corresponding logs have been inspected.

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

**BDF2 state-block diagnosis (2026-11-16):** With the original
1-rank open-loop configuration (rotation 0.02, rigid initial
momentum, `init_dt=1e-9`), PETSc solves and Newton converge and
TS accepts the BDF2 startup step, but SCRIMP correctly rejects
the final state norm. The only wildly amplified block is the
boundary traction multiplier `Y`: 0 initially, `1.19961e7`
at the intermediate t=-5e-10, and `5.99804e6` at t=0.
Meanwhile h remains identically 5 (||h||=110), p changes
only from ||p||=0.326364 to 0.329979, e_p to 0.329979,
e_h to 11.0, Y_n to 4.05736 and prescribed boundary
velocity U to 0.161759; U_n stays near zero. Maximum local
Y magnitude is ~1.78704e6 at the intermediate step and
~8.93519e5 at t=0. This is **not** a SNES convergence
failure, a MUMPS factorization failure, or a physical height
blow-up: it is an unacceptably large algebraic boundary
traction that dominates the state norm.

**Rotation=0 control passed (2026-11-17):** repeat with the same
mesh=0.2, open case, BDF2 `init_dt=1e-9`, fixed ARKIMEX
`dt=0.005`, and rigid initial profile, but `rotation=0`,
completed initial consistency and the two main steps to
`t_f=0.01` (exit 0, no rejected steps/SNES failures).
At the final initial BDF2 state, `||h||=110`, `||e_h||=11`,
`||Y_n||=4.06202`, while `||p||` and `||e_p||` are
essentially zero and the wall traction `||Y||≈6.25e-13`.
The reported maximum relative power balance defect is
`2.1517886618086806e-27`; this is a **stationary**
equilibrium check with vanishing boundary power, not a
nontrivial driven nonlinear power-balance validation.

**Critical ramp clarification:** `ramp_expression()` is
**decreasing**, with `f(t)=1` up to 75 s and 0 after 200 s
by default. Thus `f(0)=1`, **not 0**. The finite-rotation
initial rigid-body field `p/rho=omega*(y,-x)` has the
correct orientation and amplitude for the *ideal* circular
boundary, but the discrete wall condition uses facetwise
`U=omega*Rot*Normal`, which does not coincide exactly
with that field on the polygonal mesh. Normal velocity
`U_n=0` is also imposed using the facet normals.
At rotation=0, the geometric mismatch disappears and
`Y` stays close to zero. This is consistent with, but does
not yet prove, a tiny-BDF2-step `O(1/init_dt)` constraint
reaction from an incompatible initial trace; possible
alternative causes include boundary-law sign/scaling,
mixed FEM inf-sup behavior and conditioning. The original
comment suggesting the initial p matches facet normals
exactly was misleading and has been corrected.

**BDF2 inverse-step-size scaling confirmed (2026-11-17):**
for `rotation=0.02` and the same rigid initial profile
and mesh=0.2, the final startup `||Y||` was
`5,998,040` at `init_dt=1e-9`, `599,804` at `1e-8`,
and `59,980.4` at `1e-7`. Each gives
`init_dt * ||Y|| = 0.00599804` to the precision shown.
Both new startups exit 1 under the existing safe
norm-amplification guard. This is strong evidence for a
discrete constraint reaction scaling as `1/init_dt`
rather than nonlinear Newton or MUMPS divergence.
It does **not** by itself establish that the mismatch
is exclusively due to faceted geometry.

**Important brick-sign inspection:** `DPHS.add_brick`
negates all `position="flow"` forms. The two boundary
bricks `M_b[h]: h*U.Test_Y` (flow) and
`B_mu[h]^T: h*e_p.Test_Y` (effort) therefore assemble
the weak constraint `h*(e_p-U).Test_Y`, not
`h*(e_p+U).Test_Y`. The sign is consistent with
`e_p=U` on the wall, so do **not** flip a brick
sign based only on the printed original expressions.

**Next discriminating check:** vary the disk mesh size,
holding `rotation=0.02`, rigid profile, startup
`dt=1e-8` and conservative MUMPS reset fixed. Compare
`init_dt * max|Y|` and `init_dt * ||Y||` across
mesh sizes, noting that an unweighted discrete 2-norm
is mesh-dependent. Directly evaluating the
initial boundary-trace discrepancy `p/rho-U`
(including the imposed `U_n=0`) would be even stronger.
Keep the state-amplification guard on, and do not
interpret failed BDF2 starts as a valid nonlinear
physical simulation.

**Mesh-refinement diagnostic (2026-11-17):** at
`init_dt=1e-8`, `rotation=0.02`, rigid initial profile,
the final startup traction block was:

| Nominal mesh h | `||Y||` (coefficient 2-norm) | `max|Y|` | `init_dt * max|Y|` |
| --- | ---: | ---: | ---: |
| 0.20 | ~599,804 | ~89,351.9 | ~8.93519e-4 |
| 0.15 | 556,404 | 81,939.2 | 8.19392e-4 |
| 0.10 | 450,722 | 54,893.9 | 5.48939e-4 |

The h=0.20 row is extrapolated from the prior h=0.20,
`init_dt=1e-9` value using the verified inverse-step scaling;
it was not a new h=0.20, 1e-8 solve. Both additional startups
are rejected by SCRIMP's existing state-amplification guard;
height stays exactly 5. From nominal h=0.20 to h=0.10,
the maximum traction coefficient is ~39% lower, consistent
with a reduced geometric mismatch but not yet proof of
its causality or an identifiable mesh-convergence order.
Coefficient norms are not mesh-independent physical traction
norms.

**Direct geometric-discrepancy probe added:**
`sandbox/check_swe_initial_boundary_mismatch.py` generates the
same SCRIMP/Gmsh Disk meshes without solving the PDE and uses
GetFEM boundary quadrature with `X` and `Normal`
to compute L2/RMS norms of the mismatch between
`omega*Rot*X` (rigid initial velocity) and
`omega*Rot*Normal` (imposed wall velocity at `t=0`),
plus the normal velocity violation
`omega*(X(2)*Normal(1)-X(1)*Normal(2))`.
The diagnostic verifies the normal/tangential orthogonal
decomposition. **The geometric test passed on the target Docker
environment (2026-11-17)** and returned the following independent
boundary integrals for `rotation=0.02`:

| Mesh h | Boundary length | Full trace mismatch L2 | Normal mismatch L2 | Tangential mismatch L2 | Full mismatch RMS |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.20 | 6.2736988581 | 0.002758570977 | 0.002749229120 | 0.000226832706 | 0.001101342322 |
| 0.15 | 6.2773278613 | 0.002166518726 | 0.002161983129 | 0.000140115465 | 0.000864719079 |
| 0.10 | 6.2805815932 | 0.001443787865 | 0.001442442853 | 0.000062305815 | 0.000576107403 |

From nominal h=0.20 to 0.10, the integrated L2 trace mismatch
drops ~47.7%; the mismatch is dominated by the **normal**
component (roughly 99.3% to 99.8% of squared L2 mismatch).
This independent geometric measurement supports the faceted-trace
hypothesis and is consistent with the drop in observed `max|Y|`,
but does NOT establish an equality between a velocity-trace norm
and the nonphysical-sized coefficient vector of the boundary
traction multiplier. The `-info` PETSc option warning in this
geometry-only test is nonfatal; the exit code is zero.

**Controlled causal test added:** the Saint-Venant sandbox now
accepts `--wall-velocity-profile rigid_trace`, which replaces the
default `Rot*Normal` wall actuation by `Rot*X` while keeping the
original rigid initial momentum. Because `ramp(0)=1`, this
makes the imposed vector U coincide **exactly** with the initial
`e_p=p/rho` on the mesh at t=0; no other equation is altered.
The new mode is restricted to the open case with rigid
initial momentum. It is **diagnostic only**: `Rot*X`
has a nonzero normal component on polygonal facets and does
not implement a physically admissible impermeable wall.
**Matched-vector-trace A/B test rejected the simple geometric
explanation (2026-11-17):** at mesh=0.2, rotation=0.02,
rigid initial profile and BDF2 `init_dt=1e-8`, the original
`facet_tangent` wall yielded final `||Y||≈599,804`.
The diagnostic `rigid_trace` wall, which exactly matches the
initial rigid velocity *vector* along the polygon boundary,
instead yielded `||Y||=603,188.64` with
`max|Y|=90,606.2`; the corresponding original
`max|Y|` was ~89,351.9. Both startups are rejected by the
unchanged safety guard. Changing the wall-vector profile
does **not** decrease Y; therefore the measured geometric
vector mismatch cannot be presented as its dominant
cause. In addition, `rigid_trace` violates the separate
normal no-flux condition on facets, so this A/B is not
a complete DAE-consistency test.

**New targeted DAE-guess diagnostic:** SCRIMP's standard
`set_initial_value()` only initializes dynamic states h,p;
GetFEM's algebraic e_h,e_p,U_n,Y_n,U,Y start at zero.
The initially nonzero p therefore disagrees with the
constitutive relation e_p=p/rho and the velocity
control U=omega*Rot*Normal (or omega*Rot*X in the
altered-wall experiment). The sandbox now accepts
`--algebraic-initial-guess ep` to prefill e_p=p/rho and
`--algebraic-initial-guess ep_wall` to prefill both
e_p and U with the exact rigid vector trace in the
open-case rigid_trace experiment. Both are *initial
guess* interventions and neither constructs guaranteed
consistent DAE initial data. The normal no-flux and
remaining algebraic constraints still require checking.
The script validates each GetFEM variable's length
and prints its initial norm, failing early on mismatch.

Compare three separate short runs at the same mesh=0.2,
rotation=0.02 and init_dt=1e-8:
`rigid_trace`+zero guess (already failed, Y≈603,189),
`rigid_trace`+ep, and `rigid_trace`+ep_wall.
Keep the strict state-amplification guard and safe MUMPS
reset enabled. An unchanged Y would suggest the impulse
is not merely caused by zero algebraic initial guesses,
and would motivate a proper DAE-consistent initialization
or a discrete boundary inf-sup/constraint audit.
Do not treat any one passed startup as physical validation.
```sh
python sandbox/check_swe_initial_boundary_mismatch.py \\
    --meshes 0.2 0.15 0.1 --rotation 0.02
```

**Algebraic-guess A/B tests rejected the initial-guess explanation
(2026-11-17):** with the open `rigid_trace` case,
rotation=0.02, mesh=0.2, and BDF2 `init_dt=1e-8`,
the baseline zero algebraic guess, the `ep` warm start
(prefilled ||e_p||=0.326363837), and the `ep_wall`
warm start (also prefilled ||U||=0.110561624) all
reached the same rejected BDF2 state
`||Y||≈603,189` (`max|Y|=90,606.2`).
The final total norm reported by the guard was
603188.64 in both warm-start runs. Initial guesses
alone do not remove the observed startup impulse.

**Pseudo-TS versus BDF2, controlled diagnostic (2026-11-17):**
The sandbox now exposes `--init-method {bdf,pseudo}`
(default bdf) and `--init-pseudo-steps N` (default 1),
using the existing `DPHS.init_step_ts_type` and
`init_step_nb_iter` API. The previous Saint-Venant
pseudo initialization had produced an enormous
(~1e24) state norm, so the opt-in
`init_step_max_norm_amplification=100` safety check is
now applied to **all** initial TS schemes rather than
only negative-start BDF2. Diagnostic environment
`SCRIMP_INIT_DYNAMIC_STATE_AUDIT=1` prints actual
changes in the differential states h and p after the
initial solve. This is essential because PETSc
`TSPSEUDO` targets `F(z,0)=0` (steady state):
it does **not** guarantee a projection preserving the
prescribed initial dynamic state or a physically valid
DAE-consistent state. Do not treat reaching the maximum
pseudo-step count as proof of steady-state convergence.

Suggested isolated comparison at one rank (guard enabled):
```sh
for guess in ep ep_wall; do
  SCRIMP_INIT_DYNAMIC_STATE_AUDIT=1 \\
  SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN=1 \\
    python sandbox/viscous_swe_2d_mpi.py \\
      --case open --mesh 0.2 --rotation 0.02 \\
      --initial-profile rigid --wall-velocity-profile rigid_trace \\
      --algebraic-initial-guess "$guess" \\
      --init-method pseudo --init-pseudo-steps 1 \\
      --init-dt 1e-8 --final-time 0.01 \\
      --max-step 0.005 --save-dt 0.005 \\
      --adapt-type none --newton-block-trace
done
```
**Warm-started pseudo test actually failed (2026-11-17):** the
open, rigid, `rigid_trace`, `ep_wall` configuration with
`init_method=pseudo`, one pseudo step of `1e-8`,
`mesh=0.2`, and `rotation=0.02` exited 1 under the
norm-amplification safeguard. Initial e_p and U guesses
were correctly installed (norms 0.326363837 and
0.110561624). The final `||Y||` remained ~603,189
and the total state norm was 603188.655 from an
initial 110.001024. The dynamic audit (using the *accepted
PETSc TS Vec*, not an intermediate GetFEM iterate) found:
`h`: absolute change 3.70929095e-8, relative
3.37208268e-10; `p`: absolute change 0.334753159,
relative **1.02570543** (norm from 0.326363837 to
0.329912882). Thus one very short pseudo step
substantially alters the rotating dynamic state while
still producing the same huge multiplier. The
warm-started pseudo strategy cannot be considered a
valid DAE-consistency projection. This does not prove a
bug in TSPSEUDO, nor that every possible choice of pseudo
parameters fails: only the stated configuration is tested.

**Recommended next controlled experiment:** test the previously
implemented **start-from-rest smooth spin-up** at
`--start-from-rest-duration 0.02`, `--mesh 0.2`,
`--rotation 0.02`, `--init-dt 1e-8`, `--max-step
0.002`, `--final-time 0.01`, open case, physical
`facet_tangent` wall, unchanged BDF initial TS and
guard. Start at p=0, so the zero-rotation control is
a meaningful reference initial state, and increase
the actual physical wall drive smoothly after t=0.
This is a *different initial-value problem*, not a fix
for rotating p at t=0. Inspect startup acceptance,
nonzero h/p evolution, power supply/dissipation, TS
endpoint and the energy balance. If it works, next
investigate DAE-constrained projection for the original
rotating field rather than ever disabling the guard.

Keep the standard BDF2 physical baseline unchanged.

**Separate alternate initial-value experiment:** The sandbox
also accepts `--start-from-rest-duration 0.01`, restricted
to the open case with rigid profile and physical
`facet_tangent` wall. This starts h=5, p=0 and
smoothly increases wall rotation from 0 at t=0 to
the target rotation at 0.01 s, preserving the original
long-time decaying ramp. It avoids prescribing an
initially rotating dynamic state and is useful for
testing transient dynamics and Hamiltonian/port balance;
it intentionally solves a *different physical initial-
value problem* and does not repair consistency of the
original rotating state. The guard remains enabled.

Never disable the state-amplification guard merely to
force the main simulation forward.

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
