# SCRIMP v1.3.0 — DAE consistent initialization: mathematical method, PETSc implementation, and integration contract

**Current state (20 November 2026):** the PETSc-distributed
algebraic DAE initializer is **implemented in production code,
scientifically qualified as an explicit option and NATIVE
QUALIFIED as the true implicit DEFAULT** on two independent
linear 2D port-Hamiltonian models: conservative and dissipative
acoustic waves, MPI 1/2. The frozen no-option acceptance
at source `3a17a65b7a17fbdade82bd3ea2bd3bb8cd8a53f8`
passes complete CN trajectories, signed boundary/damping
power, exact prescribed differential initial data and
incompatible Gaussian rejection; explicit
`pseudo` and negative-time tiny-step BDF2 remain available.
**This is not yet a tagged v1.3.0 release.**
Unsupported/nonlinear systems are rejected explicitly,
without silent pseudo fallback. Broader model, MPI/memory
and packaging gates remain open. See the release sprint.

See [the complete, chronological scientific evidence](V130_DAE_INITIAL_VALIDATION.md)
(section 19 for the final four-case PASS), and
[the release acceptance plan](V130_RELEASE_PLAN.md).

## 1. Problem: DAE consistency is not the same as the prescribed initial state

The PFEM discretization of a port-Hamiltonian model produces a DAE of
the generic form

\[
F(t,z,\dot z)=0,\qquad z=(z_d,z_a),
\]

where `z_d` contains **differential state coefficients** (in the wave
benchmark, `q` and `p`), and `z_a` contains **algebraic coefficients**
(coenergy/effort and boundary-port unknowns). The particular partition
is not inferred from variable names: it comes from the **actual**
SCRIMP state and port definitions and is checked against the
assembled matrix structure.

The user prescribes `z_d(t_0)=z_{d,0}`. Before starting physical-time
integration, the discrete constitutive and boundary relations must
also satisfy

\[
F_a(t_0,z_{d,0},z_a,\dot z=0)=0.
\]

Here `F_a` denotes the *actual rows of SCRIMP's PETSc IFunction*
classified as algebraic. Their rows in the **assembled mass matrix**
must vanish. There is **no claim** that the *entire*
`F(t_0,z,0)` should vanish: the differential equations may require
nonzero time derivatives at `t_0`.

The historical `init_step_ts_type=pseudo` solves a short artificial-time
PETSc TS problem to reduce the constraints. It can move the dynamic
state as well as the algebraic variables. In the manufactured
wave test it shifts the prescribed dynamic coefficient vectors by

| Coefficient-vector L2 drift after historical pseudo-TS | Value |
| --- | ---: |
| `q_0` | `1.95958812e-3` |
| `p_0` | `1.30196443e-5` |

The pseudo-TS solution still satisfies the signed energy/port-power
balance **after** it has chosen its modified starting state. This
does not restore the *user's specified initial conditions*.
Initial-data fidelity and port-Hamiltonian conservation are therefore
**independent acceptance criteria**.

## 2. What the tested algebraic projection actually solves

Let `z^0=(z_{d,0},z_a^0)` be the GetFEM initial coefficient vector
before any pseudo-time simulation. Set

\[
r_a^0=F_a(t_0,z^0,0).
\]

For the tested model, the algebraic residual is **linear in the
algebraic coefficients**. Denote the corresponding algebraic
block of the assembled stiffness by

\[
K_{aa} = K[I_a,I_a],
\]

where `I_a` lists the algebraic rows/columns. With all differential
coefficients **held fixed**, the correction is sought from

\[
K_{aa}\,\delta_a=-r_a^0,\qquad
z^{1}_a=z^0_a+\delta_a,\qquad
z^{1}_d=z^0_d.
\]

For a singular, redundant or incompatible algebraic block, a
least-squares linear solver can converge without satisfying the
original constraints. We **never** interpret an LSQR convergence
reason alone as a compatible DAE initial state.

After each correction the diagnostic re-evaluates the **original,
unaltered SCRIMP `IFunction(t_0,z,0)`** (including constitutive terms,
GetFEM RHS and prescribed control) and measures
`\|F_a\|_2` in the coefficient-space Euclidean norm. This
is *not* a physical energy norm or an a priori error bound.

The **unchanged scientific feasibility criterion** of the campaign is

\[
\|F_a(t_0,z,0)\|_2
 \le 10^{-10}+10^{-8}\|F_a(t_0,z^0,0)\|_2.
\]

To avoid differences in how MPI 1/2 stop once that comparatively loose
feasibility gate has passed, the **manufactured positive-control
diagnostic** separately asks for a tighter iterative refinement
target `min(feasibility_tolerance, 1e-11)`. This does **not**
relax the scientific acceptance threshold. It is currently a
*diagnostic-only target*, **not** a recommended universal default
for every future SCRIMP model.

At most **three** distributed LSQR correction passes are attempted.
Each pass archives the PETSc convergence reason, iterations,
reported KSP residual, correction norm, and true recomputed
`IFunction` algebraic residual. Refinement stops on target,
KSP failure, or stagnation. If a later step worsens the true
residual, the best evaluated state is restored and verified.
Unsupported or incompatible data produce an **explicit rejection**,
never a silent fallback to pseudo-TS.

### Structural assumptions, enforced in the prototype

The initial positive/negative tests cover only a model for which:

- GetFEM variable intervals form a **complete, disjoint** partition
  between differential states and algebraic efforts/ports;
- all algebraic rows have zero entries in the **actual mass matrix**
  up to the checked `1e-12` coefficient threshold;
- mass and stiffness are classified **linear** by SCRIMP,
  with the required assembled algebraic Jacobian valid at `t_0`;
- there are **no lumped variables** and **no substituted ports**;
- the real-valued distributed PETSc algebraic submatrix and
  the native residual use compatible variable indexing.

This deliberately does **not** qualify nonlinear state-dependent
mass matrices (even though SCRIMP v1.3.0 can support them),
higher-index DAEs, general feedback controls, or all possible
time-dependent stiffness structures. Before a production feature
is advertised for more models, those cases require explicit
structural checks and their own physical test problems.

## 3. PETSc/MPI implementation: what is distributed, what is not

The current diagnostic implementation is
`sandbox/probe_v130_petsc_mpi_initial.py`.
At initialization, `DPHS._solve_impl` has already assembled the
real GetFEM/PETSc mass, stiffness, RHS and PETSc residual work vectors.
The experimental script overrides **only the freshly constructed
`DPHS` instance's** `init_step` method, before SCRIMP creates its
unchanged main physical Crank–Nicolson TS. The source
`scrimp/dphs.py` has not been modified by this prototype.

The path is:

1. Read the actual GetFEM intervals through the diagnostic
   `_intervals` helper and classify algebraic/differential DOFs.
   Reject overlaps, uncovered DOFs and unsupported model features.
2. Verify each owned algebraic mass-matrix row is structurally zero.
   Use the owning PETSc row range to build a distributed `PETSc.IS`.
   Convert indices to **`PETSc.IntType`**, not assumed NumPy `int64`:
   the qualified rc3 PETSc build uses 32-bit `PetscInt`.
3. Assemble the actual original `IFunction` algebraic residual on
   the native PETSc-distributed vector. Extract
   `K_{aa}` by `Mat.createSubMatrix(algebraic_IS, algebraic_IS)`.
   The entire matrix is **never gathered to MPI rank 0**.
4. Solve the distributed correction with PETSc `KSP LSQR` and
   `PC NONE`. The trial uses the options prefix
   `scrimp_dae_trial_`, default internal KSP
   `rtol=1e-13`, `atol=1e-13`, and a bounded
   iteration count. The *actual final native residual*, not the
   internal KSP stopping criterion, decides compatibility.
5. Insert corrections **only into owned algebraic entries** of
   the distributed PETSc state. After assembly, SCRIMP's
   existing `Scatter.toAll` synchronizes an entire global state
   to the GetFEM model on every MPI rank.
   Verify each prescribed dynamic coefficient is **bitwise unchanged**.
6. Repeat if the strict diagnostic refinement target is not yet
   met, within the bounded iteration/stagnation safeguards.
   On rejection, restore the original prescribed GetFEM state
   and abort before any physical CN time step.
7. On success, publish the corrected algebraic state directly
   at `t_0`, **without constructing or executing any
   pseudo-time PETSc TS**, and let SCRIMP's original
   main CN integrator proceed unchanged.

**Important memory limitation:** although `K_{aa}` is PETSc
distributed and never assembled on rank 0, the prototype
**still replicates the full state vector on every rank**
using SCRIMP's existing `Scatter.toAll` and GetFEM variable
representation. The current validation proves distributed
linear algebra and MPI-2 agreement, **not end-to-end
constant-memory-per-rank scaling**. Production work must
measure memory and eventually address this replication where
GetFEM and the rest of SCRIMP permit it.

The implementation cleans up temporary `KSP`, submatrix, index set
and PETSc work vectors. The production path must be audited for
**all** early exits, including a collective MPI rejection that
cannot leave only one rank inside a PETSc collective.

## 4. Physical positive control and incompatible-data negative control

The benchmark is the published `examples/wave.py`: 2D acoustic
port-Hamiltonian equations, **2,230 total DOFs, 1,175 algebraic** at
mesh size `h=0.2`.

### Positive: exactly FE-representable initial data

The diagnostic changes **only** the wave initial momentum
inside model construction (not in `examples/wave.py`):

\[
q_0=(0,0),\quad p_0=x(3-x),\quad \rho(x)=3-x,\quad
e_p=\frac{p_0}{\rho}=x.
\]

The polynomial `p_0` (degree 2) and its coenergy `e_p`
(degree 1) are exactly representable in the tested CG2 space.
At the left boundary `x=0`, `e_p=0`, matching the
prescribed initial left control `U_L(0)=0`.
This establishes a concrete positive-control candidate for
the tested weak constitutive and boundary equations.

The diagnostic drives the left boundary after `t_0` with
`U_L(t)=0.08 sin(20t)`, and runs **10 fixed Crank–Nicolson
steps** of `dt=0.002` to `t_f=0.02`.
A nonzero boundary excitation is essential: a zero-work
test alone would not certify the boundary power bookkeeping.

### Negative: original published Gaussian

For the unmodified published control `U_L(0)=0` and

\[
p_0(x,y)=3^{-20((x-0.5)^2+(y-0.5)^2)},\quad q_0=0,
\]

the attempts to correct only algebraic variables leave a
nonzero native residual exceeding the campaign threshold.
The diagnostic must reject these data **before the physical TS
starts**, rather than silently adjusting `p_0` or `q_0`.
This is a tested **numerical rejection at a stated tolerance and
bounded correction policy**, **not** an exact analytic proof
that the Gaussian is infeasible in every possible sense or
that all reasonable discretizations must reject it.

## 5. Final native results (20 November 2026)

**Immutable run provenance:**

- Qualified image: `scrimp-mpi:v1.3.0-rc3`.
- OCI source revision:
  `1d6d0c90115213298047b28237a7f38087563a35`.
- Diagnostic source revision:
  `4d3b4d0cdee4e9b66d49170f5d128466c74e5ad8`.
- Developer-owned original output:
  `sandbox/outputs/v130-petsc-dae-mpi-h0p2-4d3b4d0c-gb0qNKB8/`.

| Measured quantity | MPI 1 | MPI 2 |
| --- | ---: | ---: |
| Positive-case initial algebraic residual `\|F_a^0\|_2` | 0.16250682 | 0.16250682 |
| Positive-case projected residual `\|F_a\|_2` | **9.8141843e-14** | **9.6187779e-14** |
| Feasibility tolerance (unchanged) | 1.7250682e-9 | 1.7250682e-9 |
| Strict positive-control accuracy target | 1e-11 | 1e-11 |
| Distributed LSQR passes | 2 | 2 |
| Last-pass KSP iterations | 410 | 399 |
| Prescribed `q_0,p_0` changed? | **No** | **No** |
| Absolute maximum signed `H-H_0+W` defect | **5.6085329e-15** | **6.0526306e-15** |
| Signed integrated external port work | **-5.3954589e-5** | **-5.3954589e-5** |
| Negative Gaussian residual after projection | 4.8224432e-6 | 4.8224432e-6 |
| Gaussian feasibility threshold | 2.0689895e-10 | 2.0689895e-10 |
| Gaussian rejected before CN? | **Yes** | **Yes** |

The **scaled L2 discrepancy of the fully projected states at
`t_0`** between 1 and 2 MPI ranks was
`9.21579389e-13`. The **maximum scaled L2 difference
across the eleven complete CN snapshots** (including `t_0`)
was `9.22182758e-13`, well below the unchanged
`5e-8` parity gate. The worst snapshot was the tenth step
at `t=0.02`.

The physical power identity in this convention is

\[
H^{n+1}-H^{n}+
\Delta t\sum_{\Gamma}P_{\Gamma}^{n+1/2}=0.
\]

The measured defect is close to floating-point precision,
with nonzero signed port work. Both MPI runs reproduce the
previous independent diagnostic algebraic-only result to
the prescribed comparator tolerance.

**Conclusion supported by native data:** the particular
linear acoustic model admits a distributed algebraic-only
initialization respecting the prescribed differential
coefficients, the main CN trajectory MPI parity, and the
port-Hamiltonian power balance with a driven port,
**without pseudo-TS**.

**Not established:** general nonlinear/higher-index DAE
initialization, all meshes and time schemes, large-scale
memory efficiency, absence of a nullspace/gauge issue on
other models, or production API/lifecycle behavior.

## 6. How to reproduce (frozen candidate image; diagnostic only)

On a host with the qualified image and branch checked out:

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-petsc-mpi-dae-initial.sh
~~~

The runner binds diagnostic files **read-only** into
the immutable image, executes the positive and
negative controls on MPI 1/2, saves per-run JSON
and (for the positive control) full `.npz` trajectories,
then runs
`sandbox/compare_v130_petsc_mpi_initial.py`. The
comparative report includes t0/trajectory MPI discrepancies,
true IFunction residuals, negative rejection, Hamiltonian
and signed nonzero port work. The provenance manifest
retains the Docker OCI revision and diagnostic Git commit.
The command does **not** rebuild the image or modify
`examples/wave.py` or `scrimp/dphs.py`.

## 7. Next engineering step: a guarded, explicit production option

**Recommended immediate work:** extract the diagnostic projection
into a small **production-owned helper** and add
an **explicit opt-in** dispatch in `DPHS.init_step`
through `init_step_ts_type=algebraic`.
Keep `pseudo` **the default** and do not alter its behavior.
The new mode must **not** construct a pseudo-TS. It must
write a properly synchronized accepted state at exactly `t_0`
before constructing the existing main TS.

The implementation should be staged and accepted through
independent gates:

1. **Module/API:** a production-owned, unit-testable PETSc
   projection helper; clear responsibilities for classification,
   mass-row checks, native residual evaluation, KSP LSQR,
   model synchronization and temporary object ownership.
   A linear-only model guard. Do not import from `sandbox`.
2. **Opt-in behavior:** accept `init_step_ts_type=algebraic`
   only with `init_step=true`, keep legacy `pseudo`/`bdf`
   unchanged, and give unsupported layouts/DAEs explicit
   errors. Distinguish *physical feasibility* from the
   separately user-configurable *accuracy target* and
   document that target's scope. No silent pseudo-TS fallback.
3. **Correctness and cleanup:** verify full initial q/p preservation
   and the true IFunction residual, record all refinement
   passes, restore the initial model state on rejection,
   destroy PETSc temporaries, and avoid rank-local
   exceptions that deadlock other ranks in collectives.
   In particular, carefully audit `set_time_scheme()` /
   `init_step` access to the `init_step_dt` setting:
   algebraic mode must **never depend on pseudo-time dt**.
4. **Native regression:** rerun the **same** polynomial/Gaussian
   1/2-rank experiment against the **immutable diagnostic**
   references above; compare full CN states, individual
   signed port powers, Hamiltonian defects, negative
   rejection and default pseudo unchanged. Qualify
   explicit t0 physical start (no hidden extra TS).
5. **Performance and memory:** compare initial-phase wall
   time, distributed operator size/nnz, PETSc iterations,
   peak per-rank RSS and lifetime of temporary matrices.
   Test at least a second linear pH-DAE, then more ranks
   and/or a larger wave. Treat replicated
   `Scatter.toAll`/GetFEM state as a real scalability
   limitation until measured/addressed.
6. **Release boundary:** only after the opt-in option
   passes its native regressions should it appear in
   Sphinx user-facing docs with a clear linear-only
   warning. Broader nonlinear, high-index, lumped,
   Saint-Venant, Maxwell and implicit-coupling support
   are **not** implied by this result.

**Next deliverable:** a narrowly scoped,
backward-compatible optional initializer, with focused
regression tests. **Not** a change to SCRIMP's default
DAE initialization, not a Docker rebuild of the qualified
rc3 image, and not a release until its own native gate passes.


## 8. First core implementation of the optional algebraic mode — native pending

**Implemented on `release/1.3.0-preparation`, 20 November 2026.**
The earlier sections describe the **qualified scientific diagnostic**.
The initial production-owned implementation is now present in:

- `scrimp/utils/dae_initialization.py`:
  `initialize_linear_algebraic(system)`, real linear-only
  structural checks, distributed PETSc `IS`/`K_aa`,
  KSP LSQR with bounded true-`IFunction` defect
  correction, explicit feasibility/rejection
  (`InconsistentDAEInitialData`), state restoration and
  PETSc temporary cleanup, plus `initialization_report`
  including per-pass diagnostics.
- `scrimp/dphs.py`:
  `DPHS.init_step()` routes only the explicitly requested
  `init_step_ts_type="algebraic"` to this helper,
  **before** the legacy `init_step_dt`/TS path.
  The current default in `set_time_scheme()` **remains
  `pseudo`** for backward compatibility until a separate
  native default-promotion gate passes.
- `tests/test_dae_algebraic_initialization_core.py`:
  lightweight source-level and isolated helper tests of
  options, complete variable layouts, known unsupported
  classes and the dispatch. Such checks **do not**
  substitute for native PETSc/GetFEM qualification.
- `sandbox/probe_v130_core_dae_initial.py` and
  `scripts/run-v130-core-dae-initial.sh`:
  run **the actual core implementation**, not a
  monkeypatched `DPHS.init_step`, in the unchanged
  qualified rc3 environment by mounting only the new
  core files at runtime.
- `sandbox/compare_v130_core_dae_initial.py`:
  independently compare the *complete* CN trajectory,
  signed per-boundary power, cumulative work, Hamiltonian
  and initial states against the frozen, previously
  qualified native prototype outputs; also check
  explicit Gaussian rejection and historical pseudo
  trajectory parity.

The core regression additionally runs the optional
**negative-time tiny-step BDF2** initializer with
`init_step_ts_type="bdf"`,
`init_step_bdf_order=2`,
`init_step_dt=1e-6`,
`init_step_from_negative_time=True`,
to document its observed behavior and detect breakage
of the legacy path. BDF2 is not required or expected
to preserve differential Cauchy data exactly.

**Native qualification command (not yet executed after
the core integration):**

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-core-dae-initial.sh
~~~

This command expects the original **developer-owned**
frozen reference output folders documented in
section 5. Do not rebuild the qualified image.
The new option is not a release feature until
**this core-specific native gate passes**.

### Explicit user request: change the default only AFTER qualification

The desired final v1.3.0 behavior is
`init_step_ts_type="algebraic"` when the user
does not specify an initialization strategy.
**That change has deliberately NOT been made.**
It would currently reject the published Gaussian
`examples/wave.py`, and the linear-only guard rejects
other historical examples with nonlinear mass/stiffness,
lumped or substituted ports. Changing the default
without a model-by-model migration plan would break
existing simulations.

The exact required default-promotion gates and the
preliminary impact inventory are recorded in
[V130_DAE_DEFAULT_MIGRATION.md](V130_DAE_DEFAULT_MIGRATION.md).
They require core native MPI parity, a second
independent linear DAE, explicit legacy choices for
unsupported examples, compatibility tests for BDF2
and pseudo, cleanup/memory measurements and an
independent native no-`init_step_ts_type` test
showing that the **new default** really does the
algebraic projection.

There is **no silent fallback to pseudo** on
unsupported or numerically incompatible data.
That is essential to preserve the meaning of
“the differential initial coefficients are fixed”.


## 9. Production-owned opt-in core successfully qualified on the wave

**The native core gate previously pending in section 8 has
now PASSED** (20 November 2026). With the qualified
`scrimp-mpi:v1.3.0-rc3` runtime, the
`DPHS.init_step(init_step_ts_type="algebraic")`
implementation reproduces all eleven *independent
frozen diagnostic* Crank–Nicolson states on MPI 1
and 2 to a maximum scaled error of **0** and
**5.30974918e-18**, respectively. The actual
native initial algebraic residuals are
**9.81418431e-14** and **9.61877785e-14**,
with exactly unchanged prescribed dynamic
coefficients, nonzero signed boundary work
`-5.39545892e-5` and Hamiltonian+work defects
`5.60853295e-15` / `6.05263063e-15`.

The Gaussian published initial data are rejected
before main CN at MPI 1 and 2 with measured
`||F_a||=4.82244319e-6` against
`2.0689895e-10` tolerance. In that case,
the GetFEM initial variables are restored.

**Native historical control measurements:**

- `pseudo`: coefficient norms
  `||Delta q0||=0.0019595881156396905`,
  `||Delta p0||=1.3019644341244006e-5`.
  The complete main CN trajectory is **bit-for-bit
  identical to the independent frozen pseudo
  physical reference** at the tested mesh/time/MPI.
- `bdf` (BDF2, negative start at `t0-1e-6`,
  fixed tiny step ending at `t0`):
  `||Delta q0||=1.9595912757338825e-5`,
  `||Delta p0||=8.680599894760661e-10`.
  This is a smaller, **nonzero** drift — not an
  exact initial data projection. Its maximum
  signed energy/work defect remains
  `4.74270349e-15`.

Exact native result:
`sandbox/outputs/v130-core-dae-h0p2-66fe02c2-WykmP41Y/`;
the exhaustive experiment and its independently
compared source baselines are archived in
[section 21 of the validation log](V130_DAE_INITIAL_VALIDATION.md).

**Scope:** the opt-in method is qualified
for *this* linear acoustic wave configuration,
including negative control and legacy regressions.
The requested global default promotion is a
**separate migration gate**:
[V130_DAE_DEFAULT_MIGRATION.md](V130_DAE_DEFAULT_MIGRATION.md).
Changing the default would otherwise break
currently unsupported examples and the
published incompatible Gaussian.

## 10. Second dissipative pH-DAE passed and default migration candidate staged

The native dissipative wave gate has now **PASSED**
(section 23 of
[V130_DAE_INITIAL_VALIDATION.md](V130_DAE_INITIAL_VALIDATION.md)).
It reuses the original `examples/wave_dissipative.py`
weak equations and resistive interior port, with a
manufactured FE-compatible `p0=x*(3-x)`,
driven boundary `U_L(t)=0.08*sin(20*t)`,
mesh `h=0.2`, 2,804 DOFs, ten fixed
Crank–Nicolson steps and MPI 1/2.

Measured native `||F_a(t0)||_2` after the
production algebraic projection:
`9.99416399e-14` and `9.93812525e-14`;
the prescribed differential state drift remains
exactly zero. The integrated dissipative work
was `D=0.0132790026>0`, and the signed boundary
work `W=-5.40249757e-5`. The maximum defects
of `H-H0+W+D` were `2.64805539e-15` and
`3.55097896e-15`. Complete CN MPI1/2
state parity reached max scaled L2
`6.3717319e-14`, well below the unchanged
`5e-8` threshold.

Frozen native output:

```text
sandbox/outputs/v130-second-linear-dae-7eec0172-Ostu9LKD/
```

With these **two independent models qualified**,
the development branch now selects `algebraic`
as the *default candidate* in both
`DPHS.set_time_scheme()` and `DPHS.init_step()`.
The published top-level examples, historical
notebooks, Maxwell paper code, legacy GUI script
generation and selected sandbox simulations
are pinned to `pseudo` explicitly where
their published behavior requires relaxation.
BDF2 with a tiny negative-time step is also
still available explicitly, and remains
documented as changing initial Cauchy data
by a small but nonzero amount.

**This is a candidate default in SOURCE ONLY.**
Before it is released, native tests must
call `set_time_scheme` *without*
`init_step_ts_type`, prove the compatible
models select the algebraic method, verify
the Gaussian is rejected without modifying
initial data, compare complete CN trajectories
and signed conservative/dissipative power
against independent frozen outputs, and run
legacy explicit pseudo and BDF2 controls.
Additionally, `init_step=False` must remain
valid without a strategy choice.

Run the staged native default acceptance gate:

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-dae-default-promotion.sh
~~~

The unchanged qualified rc3 image is used
with read-only source overlays; no new Docker
build or release tag. Broader example behavior,
memory lifetime and other DAE classes remain
separate release-quality gates.

## 11. Generalization to nonlinear algebraic constraints — proposed, NOT qualified

**Issue [#28](https://github.com/g-haine/scrimp/issues/28).**
The linear-only code is not automatically an index-independent
initializer for any nonlinear GetFEM model. The natural
extension is **Newton on the algebraic unknowns only**,
keeping user-prescribed dynamic coefficients exactly fixed.

For the full DAE residual
`F(t,z,zd)=0`, with `z=(z_d,z_a)`,
seek

```text
G(z_a) := F_a(t0, z_d0, z_a, zd=0) = 0.
```

At iteration `k`, assemble the true current-state
algebraic Jacobian (not the frozen linear
`stiffness[I_a,I_a]`):

```text
J_aa(z_a[k]) = (dF_a/dz_a)(t0,z_d0,z_a[k],0)
J_aa(z_a[k]) delta_a = -G(z_a[k])
z_a[k+1] = z_a[k] + lambda_k * delta_a
```

where `lambda_k` is chosen by a line search on a
true residual norm/merit function. The dynamic
entries must not move, including on rejected
Newton/line-search trials. On acceptance, the
**original unmodified IFunction** must satisfy
the algebraic feasibility norm, and its tolerance
must not be confused with KSP convergence.
The initial guess can use projected constitutive
values; lack of a compatible solution must be
reported rather than silently patched by pseudo-time
evolution. Singular/redundant systems need a
defined null-space or least-squares policy, and
higher-index constraints may require additional
consistency conditions. A single linearization
at the initial guess does **not** establish
nonlinear consistency.

**Production design:** factor out the already
tested distributed algebraic variable partition,
PETSc index sets, collective rank-safe failures,
rollback, exact dynamic-data checks, and
recomputed native residuals. Add a distinct
nonlinear constrained solve step only after
`J_aa` has been verified against directional
differences. Do not simply remove the
`linear_mass/linear_stiffness` guards
from the present linear LSQR code.

**Scheduling:** this is an eligible enhancement
for v1.3.0 only if it passes complete native
physical positive/negative tests before
package freeze. Otherwise preserve the
already-qualified linear default and
explicit legacy modes and continue #28
in the next feature cycle. No unsupported
nonlinear capability is advertised by
the current default.

## 12. Distinguish the separate nonlinear state-dependent mass Jacobian

**Release blocker [#27](https://github.com/g-haine/scrimp/issues/27).**

This is a DIFFERENT problem from t0 DAE initial consistency.
For a time-dependent implicit DAE residual with
state-dependent mass,

```text
F(t,z,zd)=M(z)zd+R(t,z),
```

PETSc SNES/TS requires

```text
J_TS = sigma*M(z) + dR/dz + [D_z M(z)(.)] zd.
```

Only calculating the GetFEM mass tangent `M(z)`
and adding `sigma M(z)` misses the last
term, which can prevent nonlinear convergence
and corrupt the Newton step. **SCRIMP's actual
`DPHS.IJacobian` already defaults to
`SCRIMP_NONLINEAR_MASS_JACOBIAN=finite_difference`**,
appending the contracted derivative via
`_add_nonlinear_mass_hessian_fd`.
The experimental `legacy` selection
omits that contribution. The existing
viscous SWE audit reported relative directional
FD errors around 1e-11 in MPI 1/2 on the
appropriately calibrated GetFEM MPI build.

The urgent v1.3.0 task is to **qualify the
already-implemented corrected default**
against the **actual shipped**
`examples/shallow_water.py` rather than
reimplement it or assume the prior sandbox
experiment proves the public example.
A new independent native gate uses:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-shallow-mass-jacobian.sh
```

The runner builds the unmodified shallow-water
weak forms and compares the real production
`IFunction` directional FD to the actual
`IJacobian` on MPI 1/2, with independent
nonzero `zd`, two shifts and explicit
historical `legacy` as a negative control.
It does **not** run a fake long physical
simulation or claim to initialize a nonlinear
DAE consistently. The result is **PENDING
NATIVE VALIDATION**, irrespective of a
source CI pass. See the dated critical-path
[V130_RELEASE_SPRINT_2026W48.md](V130_RELEASE_SPRINT_2026W48.md).

## 20 November 2026 — priority change: #28 before further nontrivial #6 MPI work

The native baseline diagnostic at source `e5d4e6e32b5c912eccd15786987ab1ab79fcb7f2`
(`sandbox/outputs/v130-issue6-initial-baseline-e5d4e6e3-KXzesS0i/`)
confirms a defective rotating viscous-SWE **initial consistency** step,
not a main time-stepping or MPI deadlock: the original negative-time BDF2
initializer at `dt=1e-9` produces `||Y||=5.99949e7`, while the physical
height stays at 5 and the prescribed initial momentum changes by
`||p_after-p0||/||p0||=0.760984859`. SCRIMP correctly rejects the
state with its factor-100 amplification guard.

**Decision:** issue [#28](https://github.com/g-haine/scrimp/issues/28)
is now the prerequisite to resuming the **nontrivial rotating** viscous
Saint-Venant qualification under [#6](https://github.com/g-haine/scrimp/issues/6).
Do not launch additional MPI physical/energy campaigns until
fixed-differential initial feasibility and DAE index are understood.
The already validated stationary MPI 1/2/4 gate and nonlinear
mass-Jacobian #27 are unaffected.

**Mathematical scope:** for `z=(z_d,z_a)` with prescribed `z_d(t0)=z_d0`,
the candidate constrained solve uses
`G(z_a)=F_a(t0,z_d0,z_a,0)=0` with the *current* native Jacobian
`J_aa = dF_a/dz_a`. This projection is not a general-purpose cure:
for boundary velocity constraints with Lagrange multipliers the
semi-discrete system can be higher-index or rank deficient; hidden
differentiated constraints may prevent a suitable algebraic-only
projection. The faceted-wall prescribed velocity and the original
circular rotating initial momentum are not exactly compatible;
the solver must reject mathematically inconsistent initial data
instead of altering `p0`, boundary conditions or Hamiltonian.

**Ordered acceptance gates:**

1. Structural rank/consistency audit at **unaltered** original
   `h0,p0`: zero algebraic mass rows, index-1 eligibility, correct
   algebraic DOF partition, rank/conditioning of `J_aa`, residual
   after constitutive algebraic guesses, and boundary constraint
   compatibility. Explicitly report if original data are unsatisfiable.
2. Develop a **separate** nonlinear MPI-safe, safeguarded Newton
   initializer for a demonstrably regular supported DAE class,
   using production residuals/Jacobians and a fixed differential
   initial vector; validate directional FD, descent, rollback and
   collective rejection. **Do not** simply remove the linearity
   guards on the existing native-qualified linear initializer.
3. Native positive and negative nonlinear tests, MPI 1/2, exact
   bitwise preservation of differential coefficients, nontrivial
   physical port work and no silent fallback to pseudo/BDF.
4. Only then resume #6 with an explicitly compatible physical
   rotating initial condition, original weak forms, and MPI 1/2/4
   h-weighted boundary power and viscous dissipation audits.

**Release policy:** this dependency is for the problematic nontrivial
viscous-SWE qualification. It is **not** a claim that #28 already
works or that a full generic nonlinear initializer must enter v1.3.0.
If it cannot be natively qualified within release scope, retain
the linear-only algebraic default and explicitly reject unsupported
nonlinear initializations; do not weaken the scientific gate to
meet a date.

## 20 November 2026 — issue #28.0 structural audit: source ready, NATIVE PENDING

Source-only, independently reviewable diagnostic:
`sandbox/probe_v130_issue28_structure.py` and
`scripts/run-v130-issue28-structure.sh`.
It uses the exact published experimental viscous-SWE weak bricks,
physical `p0=tangent`, faceted wall and unchanged `h0=5`.
The runner pins `scrimp-mpi:v1.3.0-rc3` and overlays only Python
modules; no Docker rebuild and no pseudo/BDF or physical time step.

Mathematical audit: dynamic state intervals `h,p`, algebraic intervals
`e_h,e_p,Y_n,U_n,Y,U`, exact DOF coverage, nonlinear tangent mass
algebraic rows, initial native `F_a`, actual GetFEM stiffness
`K(z0)=dF(t0,z0,zdot=0)/dz`, its `K_aa` algebraic block and
`K_{a,Y}` / `K_{d,Y}` subblocks, plus a central directional
`IFunction` difference when perturbing `Y` alone. Original dynamic
coefficients are checked bitwise unchanged after the probe.

**Static prediction from the open-loop weak forms:** the boundary
traction multiplier `Y` appears in `h*Y.Test_p` (dynamic momentum
row), but in no algebraic weak-form row when `gamma=0`.
Thus `dF_a/dY` is an identically zero block and `J_aa` is
structurally singular. A naive Newton that solves only `F_a=0`
cannot determine `Y`, regardless of the initial algebraic guess.
This does **not** by itself prove differentiation index 2;
differentiated boundary constraints and the coupled mass system
must still be examined before designing a supported initializer.

Runner requires **native MPI 1 and 2** and compares pinned
diagnostics, but **has NOT yet been executed on the user host**.
Any failure is diagnostic evidence, not permission to change the
original PDE or relax the norm guard. The qualified linear default
initializer and the closed #27 Jacobian are unaffected.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-structure.sh
```

## 20 November 2026 — #28.0 true native MPI 1/2 structural audit PASSED

The original single- and two-rank GetFEM/PETSc assembly audit
`scripts/run-v130-issue28-structure.sh` has now passed natively at
`bc0b42846a1fd86c58952def9a0c9163c9bf7307`, unchanged image
`scrimp-mpi:v1.3.0-rc3` (OCI revision
`1d6d0c90115213298047b28237a7f38087563a35`).
Host-owned evidence:
`sandbox/outputs/v130-issue28-structure-bc0b4284-UkXgPXar/`.
The user-visible native output reports:

```text
AUDIT MPI=1: original h0,p0, mass rows, Jacobian Y and IFunction FD
ISSUE28 STRUCTURE PASS: Y absent from Fa, present in Fd
AUDIT MPI=2: original h0,p0, mass rows, Jacobian Y and IFunction FD
ISSUE28 STRUCTURE PASS: Y absent from Fa, present in Fd
PASS: original SWE structural MPI 1/2 parity; no time step/initialization performed
```

The actual production-assembled operator exhibits a traction `Y`
**null algebraic column** and **nonzero dynamic coupling**.
The gate checked zero mass in algebraic rows, independent
`IFunction` directional finite differences in the traction direction,
bitwise preservation of prescribed `h0,p0` and MPI parity.
It did **not** execute a nonlinear DAE initializer or claim an
index classification or compatible rotating initial data.
No numerical block norms or detailed matrix extrema were printed
in the shared console; these remain in the per-rank native
`result_1/structure.json` and `result_2/structure.json` records.

### Implication: the naive nonlinear algebraic Newton is insufficient

The weak forms of the **open-loop** model contain the following
algebraic conditions (signs follow the stored bricks and SCRIMP
`position=flow` convention):

- Interior velocity constitutive relation:
  `F_{e_p} = <h*(p/rho-e_p), Test_e_p> = 0`.
- Boundary traction-test relation:
  `F_Y = <h*(e_p-U), Test_Y>_Gamma = 0`.
- Prescribed wall velocity:
  `F_U = <-U + U_n*Normal + omega*f(t)*Rot*Normal, Test_U>_Gamma = 0`,
  with normal control `U_n=0` for `gamma=0`.

Eliminating `e_p,U,U_n` therefore yields a **weak boundary constraint
on the prescribed initial differential variable**
`p0/rho|Gamma = omega*Rot*Normal` for the faceted polygon boundary.
The initial circular rotating field agrees with the ideal circular
tangent but does not exactly satisfy this faceted trace in general.
This mismatch is independently plausible from the geometry
measurements, but the rank/size of the actual **discrete weak
compatibility defect has not yet been measured**.

The traction `Y` is a multiplier associated with this boundary
velocity constraint: it appears in the dynamic momentum equation
through `B_Y^T Y`, but not in any open-loop algebraic residual.
Hence `J_aa=dF_a/dz_a` has a zero `Y` column and cannot be inverted
as proposed for a regular index-1 nonlinear DAE. This alone
**does not establish index 2**, nor prove existence/uniqueness
of a differentiated-constraint solution.

**Next gate #28.1 (not yet implemented/qualified):** extract the
native algebraic-block norms and full intervals from both
`structure.json` files, assess algebraic constitutive projection
and the remaining boundary weak residual at fixed `h0,p0`.
Only if those constraints are feasible should the team attempt
a differentiated constraint (schematically
`C_q(q,t)*dot(q)+C_t(q,t)=0`) coupled to the dynamic mass
and multiplier traction. If prescribed data are incompatible,
fail explicitly rather than alter `p0` or prescribe an
nonphysical polygon-wall trace. Keep #28 and nontrivial #6 OPEN.

## 20 November 2026 — #28.1 fixed-h,p compatibility diagnostic: SOURCE ONLY / NATIVE PENDING

After the fully native #28.0 MPI 1/2 PASS the actual per-rank
`structure.json` files show the following identical structure:
3,300 total DOFs, 1,452 differential (`h=484`, `p=968`),
1,848 algebraic and 132 traction `Y` DOFs.
`max|M_a|=0`, `max|dF_a/dY|=0`,
`max|dF_d/dY|=0.5069655658773655`.
The independent directional FD for the full dynamic `Y` block
has norm `0.5009893370080248`, absolute Jacobian discrepancy
`5.003707553108401e-17`; algebraic FD sensitivity is exactly zero.
Original `||F_a||=0.2500610865473142` for MPI1
(`0.25006108654731424` on MPI2) with zero initial algebraic
guesses. Per-block initial `||F_a||`:
`e_h=0.0880456422545`, `e_p=0.1682909140783`,
`U=0.1626550955977`, and `Y,Y_n,U_n=0`.
**These initial-zero-guess residuals do not establish boundary
incompatibility.** In particular `F_Y=0` with `e_p=U=0` is vacuous.

New source-only, no-state-acceptance feasibility audit:
`sandbox/probe_v130_issue28_feasibility.py` and
`scripts/run-v130-issue28-feasibility.sh`.
At the unchanged original `h0,p0` it builds the true
`F_a(z0)` and `J=dF/dz` using production GetFEM/PETSc,
then removes ONLY the 132 *known null* traction `Y` columns
from the otherwise algebraic unknown set. The resulting
`1848 x 1716` rectangular native algebraic Jacobian
defines a **diagnostic** distributed PETSc LSQR least-squares
problem for the 1716 remaining algebraic unknowns.
The candidate is strictly read-only: `h0,p0` and original
`Y` coefficients are checked bitwise unchanged, and every
GetFEM unknown is restored even on error.

The audit recomputes the true `IFunction` at the projected
algebraic candidate, reports **all six** algebraic block
residual norms, and checks actual `F_a(z0+delta)` against
`F_a(z0)+J_{a,*}delta` (the original open-loop algebraic
forms are affine in algebraic unknowns for fixed `h,p`).
It reports LSQR convergence, normal-equation residual
`||J^T F_a||`, and distinguishes:
`FEASIBLE_AT_TOLERANCE`,
`INCOMPATIBLE_AT_TOLERANCE` (only with independently
certified LSQR and normal-equation stationarity), and
`INCONCLUSIVE_LSQR`. A mathematical incompatibility
claim must **not** be inferred from an unconverged
or rank-deficient numerical attempt alone.

The native harness compares MPI 1/2 on
**both** the original rotation `omega=0.2` and
a separate zero-rotation `omega=0` positive control.
Both use the same physical weak forms and polygon wall;
the zero-rotation case is **not** the original rotating
experiment. No BDF2 pseudo-step or physical TS is launched;
no nonlinear Newton, higher-index method or v1.3.0
release gate has been qualified by the source-only commit.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-feasibility.sh
```

## 20 November 2026 — proposed correction via *geometric* P2, not altered wall velocity

After stage #28.1 fixed-differential-data LSQR, the native dataset
`sandbox/outputs/v130-issue28-feasibility-9082bbee-1vi0TtA9/`
reported original omega=0.2 minimal algebraic residual
`||F_a||~0.02656138`, equally for MPI 1/2, with
the dominant e_p residual `~0.02646368`.
The separate stationary control had a far smaller
LSQR residual `5e-9--7e-9` but incorrectly exceeded the
very strict `~9.47e-10` threshold; its label
`INCOMPATIBLE_AT_TOLERANCE` must **NOT** be taken
as proof that a resting equilibrium does not exist.

**Key distinction:** the experiment had field FEM P2 but
GEOMETRIC order 1 (straight polygon edges). On the actual
polygon, a circular tangent initial `p0` has a small
normal trace; simply replacing `U=omega Rot*Normal`
with `U=omega Rot*X` would match some initial traces
but introduce **nonzero wall penetration**, so it is not
an acceptable physical correction.

An **opt-in** `geometry_order=2` parameter now exists on
`scrimp.utils.mesh.Disk` and in the viscous-SWE sandbox.
Only when selected, the Gmsh mesh is promoted with
`model.mesh.setOrder(2)` *after* mesh generation/refinement.
The existing `geometry_order=1` behavior is unchanged by
default; the PDE, boundary velocity law `U=omega Rot*Normal`,
FEM field orders, parameter values, pressure/viscosity
and initial rotating `h0,p0` are UNMODIFIED.
The quadratic isoparametric boundary approximates the
ideal circle, but is not an exact circle: **do not**
assume exact DAE consistency from geometry order alone.

New bounded native diagnostic:
`scripts/run-v130-issue28-P2-geometry.sh`.
It pins the same immutable rc3 Docker image, mounts the
changed mesh generator and sandbox modules read-only, and
runs **P2 geometric order**, separately for stationary and
original rotating p0 (both MPI 1 and MPI 2).
It reuses the native fixed-h0,p0 LSQR feasibility probe;
the reports now show the requested order, imported GetFEM
geometric point counts / element counts / transformations,
true original `IFunction` residuals, and an **informational**
comparison to P1 original `0.02656138`.
Full physics, historical BDF2 pseudo initialization and
main TS steps are **not** executed here.

The runner checks MPI equality and explicit state rollback.
It does **not** assert that P2 cures the initial-state
incompatibility; its primary qualification is whether
GetFEM really retains curved P2 geometry and what
the resulting native constraint residual is.
If the original tangential initial data remain incompatible,
we must diagnose and explain the specific discrete constraint
and advise a compatible p0 or a smooth startup from rest.
No automatic alteration of user-prescribed initial data.

Source commits:
`7875617c8e6ba6193f0cc14ffe9ec9e6594bcebc`
(opt-in curved mesh) and
`a9ebc1766a3814c2215681ae74af7385ece0aa19`
(native P2 feasibility runner and compatibility suggestions).
**Native P2 evidence remains PENDING.**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-P2-geometry.sh
```

## 20 November 2026 — corrected P2 geometry report NameError; follow-up queued

The first native order-2 attempt at source `40303b1` produced a
**genuine Gmsh P2 mesh** (225 type-9 six-node triangles and 33 type-8
three-node boundary edges). Native GetFEM reported 485 geometry points
for 225 elements (versus 131 P1 points); direct inspection of the
archived `Disk.msh` shows all 33 mid-edge boundary nodes lie on the
unit circle, ~0.004528 units away from their straight chord midpoints.
The test stopped on MPI 1 at zero rotation **after assembly**, when
constructing the JSON report: `_audit()` incorrectly referenced the
CLI-local name `args.geometry_order` and raised `NameError`.
This is a probe bookkeeping bug, NOT an observed GetFEM/PETSc physical
failure or a P2 compatibility conclusion. Source corrected to use
`system._issue28_geometry_order` and a dependency-free AST regression
suite now tests that the reporting function cannot reference the
out-of-scope CLI `args`. No residual from this aborted P2 run is valid.

**Sequencing explicitly agreed:** finish and interpret the P1/P2
geometric compatibility experiment **first**. **AFTER that**, return
to #28 with the original `examples/shallow_water.py` as a distinct
nonlinear-initialization test **without the additional curved-boundary
complication**. Do not launch that example or infer its suitability
before the present native geometry gate is resolved. Review whether
its algebraic subsystem is structurally regular before choosing a
Newton method; hold its prescribed dynamic state fixed. This deferred
task is recorded here and in issue #28 so it is not forgotten.

## 20 November 2026 — correct tangent sign on Gmsh P1/P2, GetFEM check PENDING

The user asked whether the initial rotation and imposed wall tangent
have accidentally opposite orientations. **Static source reading says
NO:** `Rot=[[0,1],[-1,0]]`, `Tangent=Rot*Normal`, and
`p0/rho=omega*(y,-x)/(sqrt(x*x+y*y)+1e-6)` all rotate clockwise.
For an outward unit radial normal at (1,0), both point down.

Independently computed from the native user archives
`v130-issue28-feasibility-9082bbee-1vi0TtA9.zip` and
`v130-issue28-P2-geometry-7286bec2-Rp90btp3.zip`:
3-point / 2-point P2/P1 curved/straight Gmsh edge map,
12-point Gaussian integration per boundary element,
outward radial-choice geometric normal from the boundary
parametric derivative, all **33 original boundary edges**.
For omega=.2 original circular `p0`, the relative boundary
`L2` norm `||p0/rho-omega*Rot*n||/||p0/rho||` is:

- **P1 geometry:** `0.0550176860263` (~5.50%).
- **P2 geometry:** `0.0000696599723462` (~0.00697%).
- For the **opposite tangent sign**, the corresponding relative
  mismatch is nearly **2.0 (200%)** on both meshes.

Thus curved Gmsh P2 reduces the **geometric** initial
tangential-trace mismatch by ~790x, whereas the previous
GetFEM native algebraic LSQR residual increased from
`0.0265613800268` (P1) to `0.0277336388365` (P2).
The latter is **not explained by a reversed analytic tangent**
or simply by the straight-sided polygon approximation.
Note: Gmsh-based geometry normals are **not** the same thing as
measuring the actual GetFEM GWFL `Normal` used in the weak bricks.

New read-only native sign diagnostic
`sandbox/probe_v130_issue28_normal_orientation.py` +
`scripts/run-v130-issue28-native-normal.sh`:
uses the **actual GetFEM generic boundary integral** on the
original model, with original rotating p0 and the production
`Tangent` macro. Directly measures both signs,
`(p/rho).Normal`, `X.Normal` and the imported P1/P2
geometry in MPI 1/2. No LSQR or TS time stepping.
It can isolate a GetFEM mesh-region normal sign/assembly issue
which is not visible from the exported Gmsh coordinates.
The native GetFEM normal test is **PENDING**; do not claim
that the Gmsh geometry-only numbers prove the sign seen
by the actual weak-form assembler.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-native-normal.sh
```

**Sequence retained:** complete this boundary/GWFL normal sign
and geometry investigation first; *afterwards* use original
`examples/shallow_water.py` for a separate nonlinear fixed-dynamic
initialization study, as requested. No changes to model physics,
initial p0, existing LSQR thresholds, or BDF2 safety guard.

## 20 November 2026 — ROOT CAUSE: GetFEM GWFL column-wise Rot reversed wall actuation

**Native evidence PASS, before correction:** user uploaded
`v130-issue28-native-normal-b903865b-91WeMgdE.zip`
(original `b903865bee829d18a329928fa04dd68e0b536eb7`,
unchanged image `scrimp-mpi:v1.3.0-rc3`).
Native GetFEM integration using **the actual production
`Normal`, `Tangent=Rot*Normal`, original p0**, MPI 1 and 2:

| Geometry | `||p/rho - omega*Tangent||/||p/rho||` | `||p/rho + omega*Tangent||/||p/rho||` | `||p/rho.Normal||/||p/rho||` |
|---|---:|---:|---:|
| P1 straight | 1.9992482297 | 0.0549180824 | 0.0548807074 |
| P2 curved | 2.0000023663 | 0.0000697293 | 0.0000696876 |

The GetFEM outward normal itself is **correct**:
`integral(X.Normal)` = `6.2452910639` (P1) and
`6.2831681228` (P2), near the expected `2*pi`.
MPI1/2 are equal to floating-point precision.
This **overturns** the previous assumption that the old
`Tangent` macro matched `p0` inside GetFEM. The geometric
P1/P2 edge-normal analysis was not sufficient to determine
the sign of the literal `Rot` as interpreted by GWFL.

**Exact root cause:** GetFEM GWFL writes a matrix
`[[a,b],[c,d]]` with the nested vectors as **COLUMNS**.
Thus the previous `Rot = [[0,1],[-1,0]]` represents the
**actual** linear operator `[[0,-1],[1,0]]`, i.e. CCW.
Original `p0/rho = omega*(y,-x)/(r+1e-6)` is CW.
The wall law therefore had the opposite tangential sign,
regardless of P1/P2 geometric accuracy. This explains why
curved P2 did **not** reduce the previous algebraic residual.
Reference: https://www.getfem.org/userdoc/gasm_high.html
(explicit matrices section).

**Corrected source commit:** `6db7a79cfe6f7b20d854c4b45913922f10fda0a9`.
The experimental viscous-SWE `Rot` macro is now written
`[[0,-1],[1,0]]` (GWFL columns), meaning the intended
row-wise clockwise matrix `[[0,1],[-1,0]]`.
Original p0, h0, PDE weak forms, viscous terms, normal
impermeability and all SCRIMP solver internals are unchanged.
The nonlinear gyroscopic form
`Gyro(p) = div(Rot*p) * Rot` is invariant under a **global
`Rot -> -Rot` sign flip**: it contains two factors of Rot.
The correction therefore flips the actual wall tangential
velocity but leaves the gyroscopic physics invariant.

A source-only test `tests/test_v130_issue28_rot_convention.py`
checks the actual column-wise matrix representation, clockwise
tangent at `(1,0)`, unchanged `Gyro`/`Curl2D` macro formulas
and formal sign invariance; CI runs it. The original
`scripts/run-v130-issue28-native-normal.sh` now rejects a
`REVERSED` or `AMBIGUOUS` native sign instead of printing
a warning and exiting 0.

**Old LSQR numbers are obsolete as physical evidence:**
the pre-fix `||F_a||_min` ~0.02656 (P1) and ~0.02773 (P2)
were assembled with a *wrong-sign wall velocity*.
Repeat fixed-state P2 feasibility with the correctly signed
wall before concluding anything about residual boundary
incompatibility, and retain the separate stationary control
and MPI 1/2 check. Even if fixed-state `F_a` becomes small,
`Y` remains a null column of `J_aa` in the open-loop system;
a complete index-aware DAE initializer is still unqualified.

**Ordered native experiments:**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-native-normal.sh
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-P2-geometry.sh
```

**Do not combine workstreams:** After the native
geometry/initial-boundary compatibility issue is interpreted,
and ONLY THEN, investigate nonlinear constrained
initialization using `examples/shallow_water.py`, as
explicitly requested. No native post-fix data yet.


## 20 November 2026 — post-fix GetFEM GWFL matrix convention: native P1/P2 evidence

**Important Python user pitfall:** the GetFEM Generic Weak-Form Language
(GWFL) writes nested explicit matrices as *column vectors*, whereas
Python/NumPy list-of-lists conventions use *rows*.

- **Authoritative upstream reference:**
  [GetFEM GWFL — Explicit matrices](https://getfem.org/userdoc/gasm_high.html#explicit-matrices).
- **User-facing Sphinx guide:** `docs/source/getfem_gwfl.rst`
  (reachable from the website index, getting started, initialization
  page, and repository README).
- **Incorrect for clockwise rotation:** GWFL `[[0,1],[-1,0]]`
  (counterclockwise as an actual matrix).
- **Correct clockwise:** GWFL `[[0,-1],[1,0]]`
  or the more readable row/semicolon form `[0,1;-1,0]`.
- Actual GWFL `Normal` is outward; in the unit disk
  `Rot*Normal` now agrees with the original
  prescribed clockwise velocity `p0/rho`.

**Native post-correction qualification of sign and reduced algebraic
feasibility (not of complete DAE initialization).** Branch source
`120c01bd8aa206ea0986ba876d65d35489a825b6`, immutable
Docker `scrimp-mpi:v1.3.0-rc3`, true GetFEM/PETSc MPI 1 and 2:

| Quantity | P1 geometric mesh | P2 curved geometric mesh |
|---|---:|---:|
| Actual GWFL tangent orientation | ALIGNED (MPI 1/2) | ALIGNED (MPI 1/2) |
| Relative aligned trace mismatch | `0.0549180824` | `6.97292855e-5` |
| Relative normal leakage | `0.0548807074` | `6.96876238e-5` |
| GWFL integral `X.Normal` | `6.24529106` | `6.28316812` |

The P2 geometric transformation was actually imported as
`GT_PK(2,2)`, not simply simulated with P2 state FEM on straight edges.

**Independent original P2 algebraic LSQR**, fixed `h0,p0`:
before the GWFL sign correction (source `7286bec`)
`||F_a||_min=0.0277336388`; after the correction (source
`120c01bd`) `4.1850871648724585e-8` (MPI 1) and
`4.183695287884114e-8` (MPI 2): a reduction of
**more than 660,000 times**. The post-fix rotating feasibility tolerance
is `2.60963115021321e-9`, so the strict label is still
`INCOMPATIBLE_AT_TOLERANCE` and **must not be advertised as PASS**.
The rotating per-block projected `e_p` norm is
`4.101738589874305e-8` on MPI 1.

At rest (P2), the original known equilibrium yielded a
diagnostic-only LSQR residual of `4.746567828624845e-9`
(MPI 1) and `4.721221465777258e-9` (MPI 2), with a stricter
acceptance tolerance `9.507859931505774e-10`.
**A fail against this numerical stopping threshold is not a proof
of physical incompatibility**, particularly for a resting
state already qualified as an equilibrium. Precision/conditioning
and subsequent high-index initialization need separate treatment.

Evidence: user-supplied archive `v130-issue28-.zip`,
the host native post-fix sign output directory
`sandbox/outputs/v130-issue28-native-normal-120c01bd-sq5HgESS/`,
and the native post-fix feasibility output
`sandbox/outputs/v130-issue28-P2-geometry-120c01bd-sKH7lgKN/`.

**Do not close #28 yet:** although the tangential sign error is fixed,
the open-loop traction `Y` remains absent from the algebraic equations
(`dF_a/dY=0`), so a complete consistent initial traction generally
needs a different, possibly constraint-differentiated formulation.
As agreed with the user, test `examples/shallow_water.py` as a
separate nonlinear-initializer experiment **only after** the current
boundary/geometric work is concluded.

## 20 November 2026 — #28 resting LSQR false-negative: exact-reference native gate PENDING

The corrected original curved-P2 fixed-state LSQR audit returned
`||F_a||=4.746567828624845e-9` (MPI 1) at `omega=0`, exceeding
its own strict numerical threshold `9.507859931505774e-10`,
even though a true resting equilibrium is known to exist.
An LSQR stop at relative tolerance does **not** certify a positive
best-residual lower bound in an ill-conditioned algebraic projection.
Do NOT interpret `INCOMPATIBLE_AT_TOLERANCE` as mathematical
infeasibility without an independent reference/certificate.

New independent native reference (source only, NOT YET RUN):
`sandbox/probe_v130_issue28_exact_rest.py` with
`scripts/run-v130-issue28-exact-rest.sh`.
This benchmark does **not use LSQR**. It holds the ORIGINAL
`h0=5,p0=0`, original physical wall controls and all weak forms.
It sets only algebraic unknowns to the exact expected stationary
values `e_h=rho*g*h0=0.5`, `e_p=0`, `U_n=U=Y=0`,
`Y_n=-e_h=-0.5`, reflecting the SCRIMP flow-brick sign
convention (`M_n[h]` registers `-h*Y_n*Test_Y_n`
and `-B_n[h]^T` registers `-h*e_h*Test_Y_n`).
Native `IFunction(t0,z,zdot=0)` checks **both**
`||F_a||` and `||F_d|| <= 1e-11` on P1 and P2 geometric
disks, MPI 1 and 2. The diagnostic includes a
wrong-sign `Y_n=+0.5` negative control, verifies exact
prescribed dynamic coefficients and complete GetFEM rollback.
No physical/initialization TS is run and no model/solver changes
are made. If the reference passes, the previous LSQR
`INCOMPATIBLE_AT_TOLERANCE` **rest label was a numerical
false negative** and numerical conditioning/termination must
be addressed independently from physical incompatibility.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-exact-rest.sh
```

**Next order:** close the geometry/sign audit on the evidence
once the known solution is calibrated (do not claim the
original rotating DAE is fully initialized), then and only
then start the separately planned `examples/shallow_water.py`
nonlinear DAE initialization qualification. This model has
state-dependent mass and a nonlinear pressure function;
first check whether its algebraic residual is actually
nonlinear in the *algebraic* variables at fixed h,p.

## 20 November 2026 — #28 exact stationary baseline NATIVE PASSED (MPI 1/2, geometry P1/P2)

**Host evidence:** `v130-issue28-exact-rest-b556f2ea-3E78L59H.zip`.
Source `b556f2ead84397b7158b1755d2dfa72ede8122df`,
immutable `scrimp-mpi:v1.3.0-rc3` OCI revision
`1d6d0c90115213298047b28237a7f38087563a35`.

The diagnostic inserted the independently known exact **algebraic**
stationary state WITHOUT LSQR, pseudo initialization or a physical TS:
`h0=5` and `p0=0` (prescribed original differential coefficients),
`e_h=0.5`, `Y_n=-0.5`, `e_p=U_n=U=Y=0`. It evaluated the
real production `IFunction(t0,z,zdot=0)` and rolled back all states.

| Geometry | MPI | `||F_a||_2` | `||F_d||_2` | Wrong `Y_n=+0.5`: `||F_a||_2` |
|---|---:|---:|---:|---:|
| P1 | 1 | `1.4633950546719204e-17` | `2.6084923555407635e-16` | `4.0700566141663` |
| P1 | 2 | `1.470126998300503e-17` | `2.6084923555407635e-16` | `4.0700566141663` |
| P2 | 1 | `1.4763422151112202e-17` | `2.699640468033266e-16` | `4.075223213705125` |
| P2 | 2 | `1.4744250640694402e-17` | `2.699640468033266e-16` | `4.075223213705125` |

**The stationary reference passes** (all < `1e-11`), with a
nonvacuous negative control. This proves that previous resting
`INCOMPATIBLE_AT_TOLERANCE` classifications (`||F_a||~4.7e-9`)
were **numerical LSQR false negatives**: identical fixed h0,p0
have an exact solution of the assembled FE constraints.
Do not loosen production tolerances merely to hide the LSQR failure.

**The #28 boundary geometry/sign *subgate* is closed as a
diagnosed and corrected GWFL literal error**, including true
P1/P2, MPI 1/2 sign and resting-equilibrium validation.
This does NOT close issue #28: for original nonzero rotation,
the corrected fixed-h,p algebraic LSQR residual remains
`~4.18e-8`, above the strict diagnostic threshold; the
traction `Y` is a null algebraic Jacobian column and consistent
initial `Y0` may require differentiated boundary constraints.
No complete rotating DAE initializer has been qualified.

**Next, as explicitly ordered by the user:** study the original
`examples/shallow_water.py` for nonlinear consistent
initialization on a straight-edged rectangular domain.
Verify the algebraic Jacobian and the existing
`explicit=True` frozen gyroscopic treatment. Preserve h0,p0.
Nonlinearity in h,p and state-dependent mass do not
automatically mean that the algebraic equations are
nonlinear in their algebraic unknowns; distinguish the
two levels before claiming a nonlinear algebraic Newton.

## 20 November 2026 — #28 second workstream BEGINS: shipped shallow-water nonlinear DAE, source only

The requested **AFTER geometry/normal sign** gate is now scheduled:
`examples/shallow_water.py` is the original published
**inviscid rectangular** PFEM shallow-water example.
The diagnosed sign convention in the separate viscous disk
sandbox is now documented and its exact P1/P2 resting
equilibrium has passed native MPI 1/2. This does not qualify
the full nonzero rotating disk initial state.

Source inspection of `examples/shallow_water.py`:

- Rectangular domain (original mesh `h=0.1`), no curved
  boundary geometry.
- Differential states `h,p`, with original prescribed
  `h0=3-(sign(x-.5)+1)/3`, `p0=0` and zero normal
  boundary inputs on all four sides.
- Nonlinear, h-dependent momentum mass (`M_p=h*p.Test_p`),
  nonlinear Dirac terms and quadratic pressure contribution
  `0.5*(p.p)/rho` in the effort constitutive relation.
- Original `G=(Gyro(p)*e_p).Test_p` is
  `explicit=True`; the *entire* term is frozen at
  the previously accepted state, and is absent from the
  implicit Jacobian, exactly as explicitly specified
  in the shipped model. The initial original `p0=0`
  makes its initial forcing zero; previous #27
  nonvacuous tests qualified its different-state freezing.
- Crucial nuance: for **fixed h,p**, the full original
  algebraic residual is **affine in e_h,e_p,U_i,Y_i**.
  This is a *nonlinear full DAE* with a linear
  algebraic subproblem, not a demonstration of
  nonlinear-in-algebraic-variables Newton.

New diagnostic-only harness:
`sandbox/probe_v130_issue28_shallow_initial.py`,
`scripts/run-v130-issue28-shallow-initial.sh`.
It reuses the source-proven original shallow-water model
builder from `sandbox/check_v130_shallow_water_mass_jacobian.py`,
which suppresses only the 5000-step physical run and plotting.
No changes to published weak forms, original h0,p0,
mesh or control law. It freezes the original
`explicit=True` gyro at the accepted initial state,
then assembles native GetFEM/PETSc `IFunction` and
`J_aa`. It checks algebraic mass rows zero,
a square algebraic block, then attempts **one**
distributed PETSc direct algebraic correction via MUMPS,
re-evaluating the true `F_a`. Independent algebraic
directional finite differences check the native tangent,
and the audit certifies `h0,p0` bitwise unchanged,
frozen gyro unchanged and complete GetFEM state rollback.

**Native run NOT YET EXECUTED; source-only staging.**
Exact finite residual acceptance remains strict `1e-9`.
The native harness compares MPI 1 and 2, pins
immutable `scrimp-mpi:v1.3.0-rc3`, records local logs
and JSON, rejects any missing or divergent original
model property, and does no initial/main TS steps.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-shallow-initial.sh
```

Even a successful diagnostic projection does **NOT** qualify
a general production nonlinear DAE initializer or
the nonlinear algebraic Newton class, and does
not replace accepted physical-time/energy-power
qualification. These must remain explicit future gates.

## 20 November 2026 — #28 shallow-water first native FD FAIL; mesh-aware criteria and rotating DAE design

The original rectangular `examples/shallow_water.py` fixed-state
diagnostic at native MPI 1, original mesh 0.1, source
`63bb31b27ef5af720f9ba4d83be823271245a394`,
failed at its directional finite-difference Jacobian test:
`relative_error=5.53585698e-5` using only `epsilon=1e-6`.
Raw host evidence: `v130-issue28-shallow-initial-63bb31b2-28COqKkK.zip`.
This is a **numerical test failure, not proof** of a faulty
production Jacobian or a failed algebraic solve.
The initial published model has `rho=1000`, `g=10`
and fluid height about 3, so coenergy pressure
`e_h` is O(30000), and the tiny epsilon may suffer
severe subtraction cancellation.

New source-only native diagnostic commits `86faa546` and
`c7e3cd05`: scan central FD steps 0.1, 0.01, 0.001,
0.0001 and 0.000001, record both absolute and relative
errors and the native `IFunction` algebraic residual.
Require BOTH the 0.1 and 0.01 steps to satisfy the
original unchanged `1e-7` relative FD tolerance;
also require true post-correction algebraic residual
`<=1e-9` and the independently assembled affine identity.
The model, dynamic initial data, frozen explicit gyro,
FE families and control ports are **unmodified**.
Crucially, write diagnostic JSON even on a numerical
rejection so the next failure is interpretable.
New native run has NOT been executed.

### Mesh size changes how to interpret `||F_a||_2`

The present implementation uses
`||F_{a,h}||_2 <= atol + rtol ||F_{a,h}^{(0)}||_2`.
It does **not** explicitly mention spatial mesh step,
but the residual is a coefficient vector from weak
forms, so its Euclidean norm changes with FE basis,
DOF count, shape regularity, integration weights
and conditioning. On a quasiuniform 2D mesh, a
fixed smooth volumetric defect with nodal load
components O(h²) and O(h⁻²) unknowns can show
Euclidean residual O(h), *without any improvement*
in physical compatibility. Volume and boundary
residuals have different natural scalings and units.

Future release acceptance must combine original
GetFEM/PETSc residuals and numerical backward errors
with **blockwise dual/Riesz norms**
`sqrt(F_a^T R_a^{-1} F_a)` using suitable Gram
operators, physical relative boundary
`L²(Γ_h)` constraint errors, and h-refinement
studies. This has **not** been implemented as
a generic production criterion; do NOT loosen
current DAE thresholds just to force a PASS.
The user-facing treatment is in
`docs/source/initialization.rst`.

### Original rotating viscous SWE: still NOT initialized

Fixing the GWFL matrix sign and validating the
exact P1/P2 stationary equilibrium do not determine
the rotating initial traction `Y0`. Native
structure confirms `∂F_a/∂Y=0` and nonzero
`∂F_d/∂Y`; naive square Jaa Newton cannot determine
this boundary multiplier. The next intended
mathematical work, after the simpler rectangular
shallow-water nonlinear DAE gate:

1. Verify compatibility of the **reduced discrete**
   algebraic boundary constraints for fixed h0,p0.
   Reject rather than silently modify prescribed data.
2. Compute regular algebraic variables (efforts and
   boundary controls) at the fixed original h0,p0.
3. Differentiate the reduced discrete trace
   constraints `C(q,t)=0`: `C_q qdot + C_t=0`.
4. Assess rank/solvability of a bordered linearized
   system coupling dynamic mass equations to
   unknown `qdot0,Y0`, with gyroscopic, viscous
   and prescribed-control contributions intact.
5. Verify original boundary power, dissipation,
   MPI 1/2 and spatial refinement before claiming
   a physical rotating DAE initializer.

This is a **planned** possible index-aware
initialization method, not an established index
class or a native-qualified production algorithm.

**Native next command:**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-shallow-initial.sh
```


## 20 November 2026 — original rectangular shallow-water fixed-dynamics algebraic diagnostic NATIVE PASS, MPI 1/2

**Authoritative user-supplied native archive:**
`v130-issue28-shallow-initial-2098e68c-4vxWCDRh.zip`.
Source `2098e68c191662103c3a404e38b7bef53298a307`,
Docker `scrimp-mpi:v1.3.0-rc3` with OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Original `examples/shallow_water.py`; original `mesh=0.1`,
original dynamic initial coefficients `h0,p0` and wall controls,
original explicit gyro.

| Metric | MPI 1 | MPI 2 |
|---|---:|---:|
| Total DOFs | `4270` | `4270` |
| Differential DOFs | `2035` | `2035` |
| Algebraic DOFs | `2235` | `2235` |
| Initial algebraic `||F_a||_2` | `1281.7812562535455` | `1281.7812562535455` |
| After **one direct algebraic solve**, native `||F_a||_2` | `2.199213556173762e-12` | `2.3277177796294783e-12` |
| Affine relation `||F_a(z+delta)-F_a(z)-J_aa delta||_2` | `1.571067363774428e-13` | `1.590053793547698e-13` |
| Best FD/Jacobian relative error | `4.246216760258867e-10` | `4.607378233546647e-10` |
| Algebraic mass-row max | `0` | `0` |
| PETSc MUMPS direct reason | `4` | `4` |

**The previously suspicious FD error at epsilon=1e-6
was a numerical differencing artifact, not a GetFEM
Jacobian defect.** At `epsilon=0.1` the relative
difference is `4.25e-10` (MPI 1), at `0.01` it is
`4.93e-9`, at `0.001` `4.73e-8`,
at `0.0001` `5.59e-7` and at
`1e-6` `5.5359e-5`. MPI 2 shows the same
loss-of-precision trend. Both sufficiently sized
perturbations pass the original unchanged `1e-7`
FD requirement; the native algebraic residual
passes the unchanged `1e-9` threshold by more
than 400x. Diagnostic status `PASS_DIAGNOSTIC`.

Native JSON further confirms: `h0,p0` bitwise
preserved, entire GetFEM model bitwise restored,
explicit gyroscopic RHS frozen and unchanged
(its norm is **zero** for original `p0=0`);
both initial/physical TS remain unexecuted.

**Qualified scope:** one direct algebraic correction
in the original full **nonlinear-state-dependent DAE**
with **affine-in-algebraic-unknowns** constraints,
fixed prescribed differential data and nonlinear mass.
This is a **read-only diagnostic only**: no
production nonlinear DAE initializer, general
nonlinear algebraic Newton, physical integration
or power balance has yet been qualified.

**Next gate:** extend SCRIMP's opt-in consistent
initializer to nonlinear *state-dependent* DAE
systems with **regular algebraic blocks affine
in algebraic unknowns**, ensuring the same
real production `IFunction` and Jacobian,
no hidden changes to initial dynamic coefficients,
explicit-gyro accepted-state freezing, robust
MPI 1/2 error/rollback paths, numerical/physical
acceptance separation, and no fallthrough to
historical pseudo TS. Do NOT silently enable
this path for the rotating viscous disk, where
`Y` is an algebraic null-column and differentiated
boundary constraints are required. #28 stays open.


## 20 November 2026 — next gate staged: opt-in *production* nonlinear-state / affine-algebraic DAE initialization

After the **original native diagnostic** of
`examples/shallow_water.py` passed at MPI1/2
(`||F_a||=2.20e-12`/`2.33e-12` after a
single direct algebraic correction), production code
now exposes an **opt-in**, distinct initializer
`init_step_ts_type="nonlinear_affine"`.

- Existing default linear `algebraic` and historical
  `pseudo`/`bdf` behavior remain unchanged.
- Original differential FE coefficients `q0=(h0,p0)`
  remain bitwise unchanged. Only a **regular,
  affine-in-algebraic-variables** block is eligible.
- True native `IFunction` initial and final residuals
  and analytic `J_aa` are used, not a manufactured proxy;
  nonlinear mass and stiffness are allowed at fixed `q0`.
- The initializer freezes the original accepted-state
  explicit forcing, performs **one direct PETSc
  PREONLY+MUMPS algebraic solve**, checks the true
  affine identity and strict algebraic accuracy,
  and commits only accepted algebraic coefficients.
- Zero algebraic Jacobian columns are rejected before
  factorization, with full state rollback on failure:
  this **explicitly excludes** the original rotating
  viscous disk multiplier `Y` whose full `J_aa` is
  singular. It does NOT handle high-index initial
  differentiated boundary constraints.
- It does not implement or claim a general nonlinear
  algebraic Newton; it is a new narrow regular
  nonlinear-*state* DAE class.
- No automatic default migration or physical TS
  release gate is claimed.

New production implementation:
`scrimp/utils/nonlinear_affine_dae.py`,
with narrowly scoped dispatch in `scrimp/dphs.py`.
New native harness:
`sandbox/probe_v130_issue28_production_nonlinear_affine.py`
and `scripts/run-v130-issue28-production-nonlinear-affine.sh`.
The runner uses the original rectangular shallow-water
model, its actual production `DPHS.init_step()`,
and the immutable rc3 Docker image, MPI 1/2.
It explicitly verifies the **published** updated
algebraic GetFEM coefficients via native `IFunction`,
original differential coefficients, and caches.
No main physical TS is executed in this stage.

**SOURCE ONLY / NATIVE PRODUCTION GATE STILL PENDING.**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-production-nonlinear-affine.sh
```

After this native gate, separate negative controls
(singular/high-index disk and incompatible fixed
differential data), normal TS short-run with power
balance, MPI parity, and h-refinement remain release
dependencies. This production stage should not be
listed as fully qualified until user-supplied native
evidence is checked.
