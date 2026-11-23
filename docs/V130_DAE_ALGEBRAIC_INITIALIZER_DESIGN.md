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

## 20 November 2026 — Production `nonlinear_affine` initial projection: NATIVE MPI1/2 PASS; host-only comparator typo fixed

The user ran the *actual* production `DPHS.init_step()` with
`init_step_ts_type="nonlinear_affine"` on the original
`examples/shallow_water.py` with original mesh `0.1`,
prescribed `h0,p0`, homogeneous boundary controls and
original explicit gyroscopic term, **without any physical TS**.
Source **`de76ff77b9ab33a0bded5b78096025511211589a`**
and pinned immutable image `scrimp-mpi:v1.3.0-rc3`
(OCI revision `1d6d0c90115213298047b28237a7f38087563a35`).

**Evidence archive personally supplied by the user:**
`v130-issue28-production-na-de76ff77-l7JjFNvy.zip`
(manifest, native MPI1/MPI2 logs and separately saved
`production_nonlinear_affine.json` for each rank).

| Metric | Native MPI 1 | Native MPI 2 |
|---|---:|---:|
| Total, differential, algebraic DOFs | 4270, 2035, 2235 | 4270, 2035, 2235 |
| Original `||F_a||` | `1281.7812562535455` | `1281.7812562535455` |
| **Published** true `||F_a||` | **`2.199213556173762e-12`** | **`2.3277177796294783e-12`** |
| Accuracy target | `1e-11` | `1e-11` |
| Affine identity error | `1.571067363774428e-13` | `1.590053793547698e-13` |
| PETSc MUMPS convergence reason | `4` | `4` |
| `max_abs_algebraic_mass_rows` | `0` | `0` |
| Algebraic Jacobian null columns | `0` | `0` |
| Prescribed differential `h0,p0` preserved bitwise | `true` | `true` |
| Algebraic corrections **published** to GetFEM | `true` | `true` |
| Accepted original-state explicit gyro cache | `ready` | `ready` |

Both reports have `status=PASS_NATIVE_PRODUCTION_INITIALIZER`,
with native final residual below the strict target.

**Host-only end-of-script mistake, not a native solver failure:**
after both runs printed `PASS`, the runner called
`math.isclose(a,b,rtol=...,abs_tol=...)`; the Python standard
library expects **`rel_tol`**, not `rtol`, so this final
comparison raised `TypeError` and made the shell exit nonzero.
Independent replay **directly against the supplied archive**
confirms that the two final residuals differ by only
`1.2850422345571612e-13` and satisfy the intended tolerance.

Commit **`e389cd8bbf419b182f58b420e45e496f30fd6a19`**
fixes this by extracting a pure-Python archived native
evidence comparator:
`sandbox/compare_v130_issue28_production_nonlinear_affine.py`,
adds positive and negative source-only regression tests,
and supports **offline revalidation without rerunning Docker**:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash scripts/run-v130-issue28-production-nonlinear-affine.sh --compare-only \
  sandbox/outputs/v130-issue28-production-na-de76ff77-l7JjFNvy
```

**Qualified scope:** the single opt-in real production initial
algebraic projection for an original **nonlinear-state**
but **affine-algebraic** DAE with rectangular geometry,
distributed MPI 1/2, strict residuals, bitwise fixed
differential data, published algebraic coefficients,
and native PETSc MUMPS convergence. No actual physical TS
step was taken, no numerical power balance was checked,
no nonlinear-in-algebraic Newton was exercised, and the
rotating viscous disk traction `Y0` remains **UNINITIALIZED**.

**Next qualification gates:** explicit native negative control
for original rotating open-loop `Y` null columns,
negative/rollback on genuinely incompatible original
initial data, short physical-time integration with
power/Hamiltonian budget, and mesh-refinement /
physically scaled residual norms. Keep issue #28 OPEN.

## 20 November 2026 — next native gate: reject high-index rotating disk (source-only staged)

The actual original rectangular shallow-water **production**
opt-in `nonlinear_affine` initialization passed natively in
MPI 1 and 2 (`||F_a||~2.2e-12`, original h0,p0
bitwise fixed). Its last host-only comparator typo was fixed
in `e389cd8` and the archived native JSON replayed
without rerunning GetFEM/PETSc.

The **next release-safety requirement** is a negative
control using the *original rotating viscous Saint-Venant*
open-loop disk with corrected GWFL wall velocity,
nonzero `omega=0.2`, and **curved geometric P2**.
Native #28 structural audit established that its
boundary traction `Y` has **zero columns in J_aa**
but nonzero dynamic coupling, so the new regular
`nonlinear_affine` algorithm MUST NOT accept or modify
its prescribed rotating h0,p0. A singular square
algebraic solve must be rejected *before* MUMPS,
without fallback to pseudo TS.

Staged (source-only; **native run pending**):
`sandbox/probe_v130_issue28_negative_high_index.py`
and `scripts/run-v130-issue28-negative-high-index.sh`.
The probe calls the real production `system.solve()`
with only `init_step_ts_type=nonlinear_affine` overridden,
catches **only** the expected
`UnsupportedNonlinearAlgebraicDAE` exception containing
`null columns`, checks at least all `Y` DOFs are detected
as null, verifies whole original GetFEM coefficient vector
is bitwise restored and no physical/pseudo TS ran.
It writes independent MPI1/MPI2 JSON and checks rank parity.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-negative-high-index.sh
```

A PASS here means **expected safe rejection**, not
that the rotating high-index `Y0` has been initialized.
The differentiable boundary-constraint/bordered
initialization of the original rotating problem,
negative incompatible-data rollback, short physical TS,
energy/power budget and mesh-aware residual norms
remain separate pending release gates.

## 20 November 2026 — actual high-index `Y` rejection MPI 1/2 PASS; post-solve rollback gate staged

**Native PASS of EXPECTED REJECTION, verified from uploaded archive:**
`v130-issue28-negative-high-index-bfe45f1f-oLsALoUF.zip`
(manifest, logs, P1 not involved; original *curved P2*
viscous rotating disk `omega=0.2`, MPI1/2).
Source `bfe45f1fdc1adf55776d6e05f7bf793ba26c374e`,
immutable rc3 OCI `1d6d0c90115213298047b28237a7f38087563a35`.

Actual `DPHS.solve()` calling production `nonlinear_affine`
raised exactly the expected `UnsupportedNonlinearAlgebraicDAE`
with **132 numerically null algebraic columns for 132
traction-Y DOFs**. The model is 3300 DOFs, 1452 dynamic,
1848 algebraic. GetFEM entire state was bitwise restored
under both MPI ranks, no pseudo/BDF or physical time step
executed; the report is `REJECTED`, and JSON confirms no
high-index rotating initialization was claimed.
Initial native `||F_a||_2 ~ 0.250963115`, tolerance
`2.60963115021321e-9`; the rejection was **structural
before factorization**, not a statement on numerical LSQR
feasibility or the differentiation index. The same
132 columns and the same exception were found in MPI1/2.

**Important limitation of that safety gate:** rejection
for singular `J_aa` occurs **before publishing an
algebraic Newton candidate**. It therefore does not
by itself validate rollback of PETSc `GlobalVec` when
an algebraic trial has been computed and fails the
final residual / precision acceptance check.

Code fix `ceb0da9e9fd3bfa7f4c657537123896297a099c8`
now snapshots the original `system.GlobalVec` before
the nonlinear-affine trial and restores it along with
the original GetFEM coefficients and explicit-RHS
metadata/cache on rejection. Successful accepted
initialization remains unchanged.

New dedicated **numerical**, NOT **physical-data**,
negative test staged:
`sandbox/probe_v130_issue28_negative_accuracy.py`
and `scripts/run-v130-issue28-negative-accuracy.sh`.
It uses unchanged shipped rectangular shallow-water
`h0,p0` and wall controls but deliberately configures
`init_algebraic_accuracy=1e-16`, well below the
previous native achievable `~2.2e-12` algebraic
residual. It demands a real computed MUMPS correction,
then an expected `InconsistentDAEInitialData` rejection;
verifies *both* GetFEM and PETSc `GlobalVec` bitwise
rollback, full explicit RHS metadata and previous
contents preserved, and no time steps, in MPI1/2.

**NATIVE NEGATIVE-ACCURACY TEST PENDING**, separate
from the previously passed high-index rejection.
This does not establish a test with genuinely
physically incompatible fixed differential data;
the latter remains an independent future gate.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-negative-accuracy.sh
```

After these safety gates, stage short actual time
integration on the original rectangular SWE and
check Hamiltonian, boundary power and explicit gyro
freeze properties before claiming physical runtime
qualification. The rotating high-index Y0 problem
still needs index-aware differentiated constraints.


## 21 November 2026 — #28 native numerical-rejection rollback gate PASS (MPI 1/2)

**User-supplied immutable native evidence**:
`v130-issue28-negative-accuracy-20b39eb2-d3D8JQf1.zip`
including manifest, rank-specific JSON and complete logs.
Source `20b39eb2e4a512b7536324331802c223dc924849`;
Docker image `scrimp-mpi:v1.3.0-rc3`,
OCI revision `1d6d0c90115213298047b28237a7f38087563a35`.

The original `examples/shallow_water.py` rectangle/mesh=0.1,
original differential h0,p0 and zero boundary controls
were **not modified**. On the real production
`nonlinear_affine` initializer, `init_algebraic_accuracy`
was deliberately set to `1e-16`, substantially below
the independently verified attainable ~`2.2e-12`
native floating-point residual.

| Metric | MPI 1 | MPI 2 |
|---|---:|---:|
| Native candidate `||F_a||_2` after successful MUMPS correction | `2.199213556173762e-12` | `2.3277177796294783e-12` |
| Strict imposed accuracy | `1e-16` | `1e-16` |
| PETSc linear solver converged reason | `4` | `4` |
| Production status | `REJECTED` | `REJECTED` |
| Exception | `InconsistentDAEInitialData` | `InconsistentDAEInitialData` |
| Bitwise GetFEM rollback | PASS | PASS |
| Bitwise PETSc `GlobalVec` rollback | PASS | PASS |
| Explicit forcing cache ready/time/content restoration | PASS | PASS |
| No physical/pseudo TS | PASS | PASS |

**PASS_EXPECTED_NUMERICAL_REJECTION** in MPI1 and MPI2.
This specifically validates rollback **after a candidate has
actually been computed**, complementing the separate native
high-index rotating `Y` rejection (132 null columns), which
occurs before attempting the algebraic solve.

The test demonstrates *numerical accuracy rejection*, **NOT**
physical incompatibility of prescribed differential initial
data, and does not determine the rotating viscous `Y0`.

Next release gate: a **short physical-time trajectory**
from the now qualified positive production initialization
on the **original rectangular** shallow-water model,
keeping original h0,p0 and controls. Independently measure
mass, Hamiltonian and zero-control boundary power;
track accepted TS steps, nonlinear failures and MPI
1/2 parity. BDF2 is **not** a discrete energy-preserving
scheme: any measured energy drift must be reported rather
than declared exactly power-preserving. Do not invent
a conservation certificate or automatically approve
the higher-index rotating DAE.

The original rotating viscous P2 case remains an open
high-index *differentiated-constraint* initialization
problem. #28 remains OPEN.


## 21 November 2026 — native MPI1/2 PASS: rejected computed algebraic trial fully rolled back

**Native user evidence archive:**
`v130-issue28-negative-accuracy-20b39eb2-d3D8JQf1.zip`,
source `20b39eb2e4a512b7536324331802c223dc924849`,
immutable image `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
This is the **unchanged original physical rectangular**
`examples/shallow_water.py` with an **artificially strict
numerical accuracy** `init_algebraic_accuracy=1e-16`
(no changes to physical `h0,p0` or wall conditions).

| Native metric | MPI 1 | MPI 2 |
|---|---:|---:|
| Original `||F_a||_2` | `1281.7812562535455` | `1281.7812562535455` |
| Computed algebraic candidate `||F_a||_2` | `2.199213556173762e-12` | `2.3277177796294783e-12` |
| Required accuracy | `1e-16` | `1e-16` |
| Native algebraic MUMPS convergence reason | `4` | `4` |
| Rejected computed correction norm | `643262.3257684928` | `643262.3257684928` |
| Algebraic affine-identity error | `1.571067363774428e-13` | `1.590053793547698e-13` |
| Initial/GetFEM rollback bitwise | PASS | PASS |
| PETSc `GlobalVec` rollback bitwise | PASS | PASS |
| Frozen explicit RHS previous ready/time/contents restored | PASS | PASS |

Both ranks returned `PASS_EXPECTED_NUMERICAL_REJECTION`
via the intended `InconsistentDAEInitialData` exception
**after** successfully computing an algebraic Newton
correction. The production report is `REJECTED` with
`accuracy_target_met=false`, no pseudo-initialization
TS and no physical time step.

**Qualified:** accepted production affine-algebraic projection
on the original nonlinear-state shallow water (earlier MPI 1/2);
safe structural rejection of rotating high-index `Y`
(earlier P2 MPI 1/2); **and now safe post-solve numerical
rejection with complete model, PETSc and explicit cache rollback**.

**Not qualified:** rejection due to genuinely physically
incompatible prescribed differential initial data;
general nonlinear-in-algebraic-variable Newton;
a physical-time run, Hamiltonian/power budget and mesh
independent physical acceptance; and the *original*
rotating viscous SWE initial multiplier `Y_0`.
Issue #28 remains OPEN.

**Next independent native gate:** a deliberately short
physical-time simulation of the original rectangular
`examples/shallow_water.py` using explicit opt-in
`nonlinear_affine` initialization. Compare full native
Hamiltonian and port-power balance with the original
inviscid homogeneous-boundary equations. Do not
modify source initial values, mesh, control laws or
the explicit-gyro accepted-state contract merely
to force an energy PASS. Keep rank-1/rank-2
comparison and explicit actual time-step counts.


## 21 November 2026 — staged first SHORT PHYSICAL TS qualification after nonlinear_affine (NATIVE PENDING)

The original published `examples/shallow_water.py` passed
native `nonlinear_affine` production algebraic initialization
on MPI1/2 with `||F_a||~2.2e-12` and bitwise fixed original h0,p0.
The independent **numerical post-solve rejection** and
full GetFEM/PETSc/explicit-cache rollback have now also
passed native MPI1/2, evidence archive
`v130-issue28-negative-accuracy-20b39eb2-d3D8JQf1.zip`.
See the preceding qualification record.

**NEXT physical gate** is a deliberately short 3-step run,
rather than another projection-only check:
`sandbox/probe_v130_issue28_shallow_short_physical.py`,
`scripts/run-v130-issue28-shallow-short-physical.sh`.

Test design:
- Load the original shipped model **without changing its weak
  forms**, rectangular mesh `h=0.1`, P2 h / discontinuous P1 p,
  original physical h0,p0, four prescribed **zero** wall controls,
  original `explicit=True` gyroscopic term, original energy terms.
- Suppress **only** the published 5000-step solve and
  plotting/export during model construction.
- Opt into `init_step_ts_type=nonlinear_affine` and run the
  **actual** production `system.solve()` and main PETSc BDF2
  integrator for `dt=1e-4` through `t_f=3e-4`.
  Use distributed PETSc/MUMPS, identical compiled pinned
  `scrimp-mpi:v1.3.0-rc3` image, MPI 1 and 2.
- Require accepted/converged TS, no rejected steps/SNES
  failures, correct final timestamp and initial saved h0,p0
  **bitwise**; original production initialization status PASS.
- Independently calculate each saved-sample total mass
  `int_Omega h`, kinetic `int h|p|²/(2rho)` and
  potential `int rho*g*h²/2`. Compare against SCRIMP's
  registered original Hamiltonian and measure drift.
- Independently calculate boundary powers `int_{Gamma_i}
  U_i*Y_i` on all four regions. With homogeneous imposed
  controls, supplied boundary power must be ~0.
- First *pilot* sanity gates: mass relative drift `<=1e-5`,
  Hamiltonian relative drift `<=5e-3`, absolute boundary
  power `<=1e-7`. Record full timeseries and defects.
  These values are **initial pilot thresholds**, not
  inferred spatial/temporal convergence certificates.
- **IMPORTANT:** BDF2 is **not symplectic/energy-preserving**
  in general. This run can check small short-time drift,
  but must never be labelled an exact discrete
  port-Hamiltonian power-preservation proof. A proper
  refinement/energy-exchange convergence study is a
  separate release gate.
- Retain rank-specific JSON/logs on failure and assert
  MPI1/2 trajectory parity. If the first full TS fails,
  diagnose actual TS reason and physical operator
  assembly; do not silently switch to pseudo TS.

**Source-only staged; native physical qualification PENDING.**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-shallow-short-physical.sh
```

This is the original **rectangular** inviscid shallow-water
problem. The original rotating **viscous, curved P2** case
has a different high-index constraint structure; its
initial traction `Y0` still requires an independent
index-aware differentiated-constraint algorithm,
and is NOT qualified by this test.


## 21 November 2026 — #28 opt-in production DAE initial projection + real short BDF2 TS staged (native PENDING)

**Previously native qualified (immutable `scrimp-mpi:v1.3.0-rc3`):**

- True `nonlinear_affine` initial projection on the
  original rectangular `examples/shallow_water.py`:
  `||F_a||=2.1992e-12` (MPI1) and `2.3277e-12`
  (MPI2), original differential h0,p0 unmodified.
- Native structural **expected rejection** of all 132
  algebraic null columns for rotating P2 traction Y.
- Native expected **post-computed-correction**
  numerical rejection at `accuracy=1e-16`,
  including exact rollback of all GetFEM fields,
  PETSc `GlobalVec`, and accepted-state explicit RHS
  metadata/content in true MPI1/2.

**Next source-only staged gate, NOT YET NATIVE PASSED:**
`sandbox/probe_v130_issue28_production_short_physical.py`
and `scripts/run-v130-issue28-production-short-physical.sh`,
with a separate archived JSON comparator
`sandbox/compare_v130_issue28_production_short_physical.py`.
This deliberately reuses the already tested original
weak-form, explicit-gyro and actual PETSc TS harness
`sandbox/probe_v130_shallow_explicit_ts.py`.

It retains the **original** `examples/shallow_water.py`
inviscid rectangular physical weak forms,
original `mesh=0.1`, original `h0=3-(sign(x-.5)+1)/3`,
original `p0=0`, and zero controls on all four
boundaries. Only the short-test time scheme
(BDF2 with four physical steps `dt=1e-4`)
and explicit selection of production
`init_step_ts_type="nonlinear_affine"` are test
overrides; there is **no manufactured nonzero
momentum**, no pseudo-initialization TS, and
no changed physical coefficients.

The probe calls the actual production `init_step()`
from within `solve()`, captures its original-state
publication report, then verifies:
- exactly four genuine accepted PETSc BDF2 physical
  steps, with no TS rejections or SNES failures;
- native `IFunction` actually evaluated throughout;
- accepted-time explicit gyro forcing frozen during
  all implicit Newton residual iterations;
- complete finite full-state trajectory and
  positive finite mechanical Hamiltonian;
- MPI1/MPI2 relative **full-state** and Hamiltonian
  parity `<=5e-8`, with identical sample times;
- the original h0,p0 unchanged by initialization and
  original published physical boundary controls.

**Important scientific boundary:** this test will
only **measure** the signed Hamiltonian drift.
BDF2 with a lagged explicit gyroscopic brick does
not guarantee an exact time-discrete
port-Hamiltonian power balance. The test must
not assert zero signed drift as a physics theorem.
A physically scaled energy, integrated boundary
power and time/mesh convergence study is a
**separate** subsequent release gate. No claim
of full rotating viscous DAE `Y0` consistency
or mesh-independent initialization follows.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-production-short-physical.sh
```

On host-side comparator issues, **replay saved
native results without Docker** using
`--compare-only <evidence-directory>` rather
than rerunning expensive physical steps.


## 21 November 2026 — TRUE shallow-water physical BDF2 after nonlinear_affine: MPI 1/2 NATIVE PASS

**User native evidence:** `v130-issue28-shallow-short-physical-d0cde0af-4DlgCMdU.zip`.
Source `d0cde0af5190290e0567bb06424d4ce8e1d2bfae`,
pinned Docker `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Original shipped `examples/shallow_water.py` and original
rectangular mesh step `0.1`, weak forms, FEM, initial
`h0=3-(sign(x-0.5)+1)/3`, `p0=0` and all four
homogeneous boundary controls. Production
`init_step_ts_type="nonlinear_affine"`; **no
pseudo-initialization**, followed by true BDF2 PETSc
TS with `dt=1e-4, t_f=3e-4`.

| Genuine native observable | MPI 1 | MPI 2 |
|---|---:|---:|
| TS convergence reason | `1` | `1` |
| Physical accepted / rejected steps / SNES failures | `3 / 0 / 0` | `3 / 0 / 0` |
| Original initial projection `||F_a||_2` | `2.199213556173762e-12` | `2.3277177796294783e-12` |
| Mass-relative drift | `2.310153069550858e-15` | `3.1986734809165835e-15` |
| Hamiltonian-relative drift / H0 | `1.437729942579708e-7` | `1.4377299161310965e-7` |
| Homogeneous boundary power (each edge and total) | `0` | `0` |
| Original H0 | `31636.2594029492` | `31636.259402949203` |
| Final H | `31636.254854509458` | `31636.254854509545` |
| Signed final H-H0 | `-0.0045484397414838895` | `-0.004548439657810377` |

Registered Hamiltonian agrees with native kinetic/potential
integrals. MPI absolute Hamiltonian sample discrepancy
`<=8.73115e-11`; relative `<=2.76e-15`.
The true mass is `~2.49903774955191` and stays at
floating-point precision. Both MPI implementations retain
original dynamic initial data exactly.

**Critical interpretation: the returned
`relative_energy_defect_exchanged=1` is an
UNINFORMATIVE metric for this closed problem, neither
a 100% failure nor a preservation certificate.**
Here the boundary work `W(t)=∫Σ_i U_iY_i dt≡0`,
hence the balance defect `D_H(t)=H(t)-H(0)-W(t)`
equals `H(t)-H(0)`. The chosen denominator is
`max(max|H-H0|,max|W|,1e-14)`, which equals
`max|D_H|` whenever the energy drifts and no
input work occurs. The ratio equals **1 by definition**.
The scientifically informative measurements are
`ΔH=-4.54844e-3` and `|ΔH|/H0≈1.43773e-7`,
plus their change under refinement in **time**.
For BDF2, exact discrete energy conservation
is NOT guaranteed; no exact pH power balance
can be claimed from this short smoke test.

**Sampling caveat:** native TS accepted `3` steps
but the stored trajectory has only timestamps
`[0,1e-4,3e-4]`, missing `2e-4`.
This does not invalidate the actual accepted
steps or final total energy and mass. However
subsequent power integration with nonzero
controls MUST use actual saved intervals,
not assume every accepted time was recorded.

**Qualified:** opt-in production nonlinear-affine
initialization followed by **real** short nonlinear
PETSc BDF2 physical-time integration,
mass to rounding, correct zero boundary power,
computed Hamiltonian and native MPI 1/2 parity.

**Not yet qualified:** exact time-discrete power
conservation, quantitative energy time/mesh
convergence, general nonlinear-in-algebraic
Newton, physically inconsistent fixed differential
data rejection, or the higher-index rotating
viscous disk initial traction `Y0`.
**Issue #28 remains OPEN.**

**Next independent pre-existing gate (already staged):**
four real physical BDF2 steps on original
`examples/shallow_water.py` using opt-in
`nonlinear_affine`, requiring complete full-state
MPI1/2 parity and checking freeze of original
explicit gyro throughout the implicit Newton
iterations. It measures Hamiltonian drift but
does not assert exact conservation.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-production-short-physical.sh
```

Afterward, time-refine BDF2 at fixed
physical `T` and track the signed
`H(T)-H(0)-∫P_boundary dt`, energy and
mass as well as full-state changes,
rather than this self-normalized
closed-system ratio.


## 21 November 2026 — Original shallow-water 4-step BDF2 after production nonlinear_affine: MPI 1/2 NATIVE PASS

**Actual user-uploaded native ZIP:**
`v130-issue28-production-short-ts-f36c6042-3Fpcl3DV.zip`,
source `f36c6042181470688f29fd4c1a75d1640f2d0c08`,
immutable `scrimp-mpi:v1.3.0-rc3` OCI revision
`1d6d0c90115213298047b28237a7f38087563a35`.

Original published `examples/shallow_water.py`,
rectangular physical mesh `h_mesh=0.1`, original
h0,p0 and zero boundary controls, unchanged
weak forms and original explicit gyro.
Production opt-in `DPHS.init_step()` with
`init_step_ts_type=nonlinear_affine` followed
by **four genuine physical BDF2 PETSc TS**
steps (`dt=1e-4`, `T=4e-4`).

| Actual native result | MPI 1 | MPI 2 |
|---|---:|---:|
| Post-projection original `||F_a||_2` | `2.199213556173762e-12` | `2.3277177796294783e-12` |
| Accepted/rejected steps | `4 / 0` | `4 / 0` |
| SNES failures / IFunction evaluations | `0 / 15` | `0 / 15` |
| Signed final Hamiltonian change | `-0.0046648074239783455` | `-0.004664807369408663` |
| Relative Hamiltonian drift / H0 | `1.474512951914752e-7` | `1.4745129346656575e-7` |
| Explicit gyro RHS drift within Newton | `0` | `0` |
| Max nonzero explicit gyro RHS norm | `2.9865184694522523e-5` | `2.9865184694484807e-5` |

The gyro is **not vacuous**: norms at accepted
step times 0, `1e-4`, `2e-4`, `3e-4` were
`0`, `4.866269214e-7`, `6.474541182e-6`,
`2.986518469e-5`, while its value remained
bitwise unchanged throughout Newton.
Native MPI1/2 full **4270 coefficient** trajectory
relative difference was `3.235700308035909e-16`;
Hamiltonian relative difference
`2.529867524588143e-15`.
All five sample timestamps `[0,1e-4,2e-4,3e-4,4e-4]`
were actually saved in this experiment, unlike
the earlier three-step pilot which missed
one intermediate saved time.

**Native qualified scope:** real production
nonlinear-state affine-algebraic initialization,
short physical BDF2 integration, nonvacuous
accepted-step explicit gyro freeze, full-state
MPI1/2 parity and measured Hamiltonian.
**Not qualified:** exact BDF2 power balance,
physical energy-defect temporal convergence,
general nonlinear algebraic Newton, or
original rotating viscous Y0.

### Next native gate: fixed T temporal refinement (SOURCE ONLY; PENDING)

To test the *actual* signed defect
`D_H(T)=H(T)-H(0)-integral_0^T P_boundary dt`,
NOT the self-normalized closed-system
ratio that equals 1 whenever
`P_boundary=0` but H drifts, stage
a fixed-final-time original physics campaign:

| dt | Required physical BDF2 steps | MPI ranks |
|---:|---:|---|
| `0.0001` | `4` | `1,2` |
| `0.00005` | `8` | `1,2` |
| `0.000025` | `16` | `1,2` |

All cases use original `h_mesh=0.1`,
exactly `T=0.0004`, original h0,p0
and homogeneous controls.
`sandbox/probe_v130_issue28_shallow_short_physical.py`
now supports `--dt` and `--final-time`
with **unchanged legacy default three-step
runtime**. It records actual accepted TS counts,
signed/relative final energy defects,
mass and all boundary powers,
and **actual** saved sample timestamps.
New runner
`scripts/run-v130-issue28-shallow-energy-refinement.sh`
retains native JSON/logs and MPI rank separation.
New host comparator
`sandbox/compare_v130_issue28_shallow_energy_refinement.py`
validates authentic 6-run evidence and outputs
adjacent refinement ratios and *observed*
log2 estimates, WITHOUT enforcing a
theoretical second-order convergence rate.
Unit tests guard missing or inconsistent
native data without expensive Docker runs.

The original h0 is discontinuous at x=0.5
and the gyro is lagged explicitly: do
not assume the formal smooth-solution BDF2
convergence rate or exact power conservation.
A monotonically smaller absolute energy
defect is a **measurement**, not itself
a complete proof of consistent physical
port-Hamiltonian integration.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-shallow-energy-refinement.sh
```

The runner provides `--compare-only <folder>` to
revalidate existing saved native evidence
without executing GetFEM/PETSc again.
Issue #28 remains OPEN for
differentiated-constraint Y0 on the original
rotating viscous disk and the remaining
mesh-dependent/physical power gates.


## 21 November 2026 — ORIGINAL shallow-water fixed-T energy-defect temporal refinement: NATIVE MPI 1/2 PASS

**Authoritative user-supplied archive inspected:** `v130-issue28-shallow-energy-refine-ae1146b8-wnHbWYHn.zip`.
Source `ae1146b89217e2ebaad5a876a9f755569aeae18e`,
immutable Docker `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.

Six independently logged **real** native GetFEM/PETSc
BDF2 time integrations, using the original published
rectangular `examples/shallow_water.py` (unchanged
weak forms, original `h0,p0`, physical boundary
controls all zero, explicitly lagged gyro, and
production opt-in `nonlinear_affine` initial
projection). The **spatial mesh is fixed**
`h_mesh=0.1`, the **physical final time is fixed**
`T=0.0004`; only temporal resolution changes.
All six reports independently confirm unchanged
prescribed initial data, true accepted step count,
no rejected steps or SNES failures, and original
algebraic initial projection `||F_a||=2.1992e-12`
(MPI 1) / `2.3277e-12` (MPI 2).

The closed-system port-power defect is computed from
the **actual signed** quantity
`D_H(T)=H(T)-H(0)-W_boundary(T)` with
`W_boundary(T)=integral_0^T Σ_i U_iY_i dt=0`
since all four original boundary controls vanish.

| `dt` | Actual steps | MPI1 signed `D_H(T)` | MPI2 signed `D_H(T)` | MPI1 final `|D_H|/H0` | MPI1 max relative mass drift |
|---:|---:|---:|---:|---:|---:|
| `1e-4` | `4` | `-0.0046648074239783455` | `-0.004664807369408663` | `1.474512951914752e-7` | `4.2648979745554295e-15` |
| `5e-5` | `8` | `-0.0011776470091717783` | `-0.0011776469655160327` | `3.7224597073002744e-8` | `3.198673480916572e-15` |
| `2.5e-5` | `16` | `-0.0002905913206632249` | `-0.0002905912660935428` | `9.185388100469784e-9` | `3.198673480916572e-15` |

At fixed `T`, the absolute signed-defect refinement
ratios are `3.9611253522046774` (4 to 8 steps)
and `4.0525883790472506` (8 to 16 steps),
giving **observed log2 orders `1.985910356921035`
and `2.018843648450496`**. Thus the evidence
strongly supports **approximately second-order
temporal decrease of this particular closed-system
energy defect over these dt values**. No BDF2
*exact* energy-conservation claim is justified:
the energy actually decreases for each finite step.

The native MPI rank-1/rank-2 maximum relative
Hamiltonian discrepancy is
`3.449818887630128e-15`, and the maximum
mass discrepancy `3.55408164546286e-15`.
Boundary power is exactly zero throughout
all reports; all original h0,p0 remain
unchanged at initialization.

**Sampling caveat:** 4 accepted steps
saved 5 states, 8 accepted steps saved 9,
16 accepted steps saved **16**, not 17.
The 16-step saved series includes both
`t=0` and `t=T` but is missing at least
one intermediate accepted step (`t=0.000375`).
The comparator validates actual saved timestamps,
not an imaginary uniform sample count. For
nonzero future boundary controls, numerically
integrating power requires the true saved time
intervals (or appropriate per-step accounting).

**Release qualification now passed:** nonlinear-state /
affine-algebraic *production* initial projection,
short true physical BDF2 integration,
nonvacuous accepted-step explicit gyro freeze
(previous run), MPI parity, mass conservation to
roundoff, and **fixed-spatial-mesh temporal
convergence of the closed-system signed Hamiltonian
defect with observed order close to two**.

**Scientific scope limitations:** the original
initial height is nonsmooth, so this short
empirical sequence does not prove a general BDF2
energy-error theorem. It does not qualify
nonzero boundary-work quadrature, mesh-independent
spatial acceptance/convergence, general nonlinear
algebraic Newton, or original viscous **rotating**
high-index traction `Y0`. The previously reported
self-normalized `relative_energy_defect_exchanged=1`
remains tautological in this closed system; use
the **signed `D_H(T)`** and the `H0`-scaled
quantity instead.

**Next substantive priority:** return to original
rotating viscous P2 disk issue #28: at fixed original
`h0,p0`, identify the *reduced discrete* boundary
compatibility constraints and determine whether their
time derivative coupled to dynamic FE equations has
sufficient rank to solve for initial `Y0`. Include
corrected GWFL rotation, full gyro/viscosity and
native MPI 1/2 parity; avoid changing original
physical data or introducing a pseudo TS.
Only then attempt a qualified consistent
high-index initialization.


## 21 November 2026 — #28 rotating viscous disk: first-differentiated discrete-constraint rank pilot (SOURCE ONLY)

The previous native six-run fixed-T shallow-water
`dt`-refinement has passed in MPI 1/2:
closed-system signed energy defect ratios
`3.9611253522` and `4.0525883790`,
observed orders `1.9859103569` and
`2.0188436485` at original h_mesh `0.1`
and final time `0.0004`.
This **empirical** result does not solve the
original rotating *viscous* P2 disk initial
traction. Focus returns to #28 high-index
initialization.

**Original native structural fact**:
`dF_a/dY=0` for the **132 traction-Y DOFs**,
while `dF_d/dY` is nonzero, so a square
ordinary `J_aa` Newton is invalid. The
rotating original `q0=(h0,p0)` is NOT
yet proven to satisfy the *exact discrete*
algebraic boundary constraints, even after
the GWFL orientation correction and curved
geometric P2 refinement. The post-correction
rotating reduced LSQR residual `~4.18e-8`
is not a mathematical nonexistence certificate:
the very same LSQR had a false-negative rest
classification, whereas the exact resting state
was subsequently verified at `~1e-17`.

At a **consistent** stationary or rotating
algebraic initial state, differentiate
`F_a(q,a_r,t)=0` once, with
`q=(h,p)`, `a_r=(e_h,e_p,Y_n,U_n,U)`
and `Y` the **unconstrained** traction.
The candidate first-differentiated block
operator is

```text
             columns: qdot      a_regular_dot      Y0
rows F_d :           M_dd           0             J_dY
rows dF_a/dt :       J_aq         J_a,a_reg          0
```

Its expected size on the original
P2 geometry `mesh=0.2` is `3300 x 3300`:
1452 original differential DOFs,
1716 regular-algebraic DOFs and
132 traction-Y DOFs.
The lower row represents
`J_aq*qdot+J_ar*a_r_dot+F_a,t=0`.
The dynamic right-hand side must
eventually include **all original**
gyro, viscosity, pressure and wall
terms; it is **NOT computed in this pilot**.

New source-only native audit:
`sandbox/probe_v130_issue28_differentiated_rank.py`
and `scripts/run-v130-issue28-differentiated-rank.sh`.
They build the *real* P2 GetFEM mass and
nonlinear tangent blocks and a PETSc sparse
two-row/three-column MatNest converted to
AIJ; perform two independent MUMPS
**synthetic** linear reconstruction tests,
reporting both x-recovery and residual
accuracy. Actual base cases:
(1) *exact resting equilibrium* `omega=0`
with only known algebraic values set to
`e_h=0.5` and `Y_n=-0.5`;
(2) original `omega=0.2` rotating physical
`h0,p0` with the original **zero**
algebraic guess and corrected wall
`Rot`. Both MPI1/2 native runs measure
`F_a` at the linearization point and
restore **all** GetFEM coefficients bitwise.
No pseudo TS, physical TS, artificial
physical forcing, algebraic projection or
actual initial Y0 solve is performed.

**Interpretation is intentionally narrow:**
successful inversion on two synthetic
right-hand sides is **a numerical
regularity diagnostic at that linearization**;
it does not rigorously prove full rank,
the differentiation index, or existence
of a compatible rotating initial state.
If the original rotating `F_a(q0,a0,t0)`
cannot be driven to the required
physically scaled tolerance **without
modifying the prescribed q0**, there
is no consistent DAE start at those exact
discrete data, regardless of the
bordered matrix rank.
Rank inference must therefore always be
reported alongside the native initial
algebraic residual.

**NATIVE RUN NOT YET EXECUTED.**
A numerically singular or inconclusive
bordered operator is a legitimate
diagnostic observation, **not a PASS of
rotating Y0 initialization**.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-differentiated-rank.sh
```

After inspecting native MPI1/2 rank/compatibility
data, design a mathematically defensible
consistent initial traction `Y0` (or
report infeasibility of exactly fixed
`h0,p0`), and only then attempt
a physical rotating TS with power
and mesh refinement. Issue #28 OPEN.


## 21 November 2026 — #28 differentiated bordered P2 viscous disk: NATIVE numerical regularity MPI 1/2 PASS (not Y0!)

**Authoritative native evidence uploaded and inspected:** `v130-issue28-differentiated-rank-f9e72ecc-pNH6VN8J.zip`.
Source `f9e72ecc0a999c76c89b3504461ad3c1c1214f94`, immutable
`scrimp-mpi:v1.3.0-rc3` OCI `1d6d0c90115213298047b28237a7f38087563a35`.
Four separate native tests: `rest` and `rotating`, each MPI 1 and 2,
on genuine P2 geometric disk, mesh step `0.2`, original FE weak forms
and corrected GWFL clockwise rotation. No pseudo/physical TS,
no physical initial traction `Y0`, no modification of initial `h0,p0`.

The matrix

```text
                  q_dot     a_regular_dot     Y0
dynamic rows        Mdd           0           JdY
d(Fa)/dt rows       Jad          Jar           0
```

has `3300 x 3300` coefficients with `1452` differential,
`1716` regular algebraic and `132` Y-traction unknowns.
The actual `J_aY` is identically zero; the actual max
absolute `J_dY` entry is `0.5072940586040593`. Algebraic
mass rows have maximum absolute value `0`.

| Native diagnostic | Rest MPI1 | Rest MPI2 | Rotating MPI1 | Rotating MPI2 |
|---|---:|---:|---:|---:|
| Original `||Fa||_2` at linearization | `1.4763422151112202e-17` | `1.4744250640694402e-17` | `0.250963115021321` | `0.2509631150213209` |
| Original `||Fd(q,zdot=0)||_2` | `2.699640468033266e-16` | `2.699640468033266e-16` | `0` | `0` |
| Maximum relative reconstruction error, two synthetic x | `4.151511441367386e-11` | `7.351843090655861e-12` | `1.7389642798601183e-10` | `2.8475843052613367e-12` |
| Maximum relative synthetic residual | `5.579034333565234e-12` | `2.4889542225783863e-13` | `1.587544485253625e-11` | `1.1058788590821502e-12` |

All eight synthetic MUMPS solves across 4 native tests converged
(PETSc reason `4`). Each result reports
`NUMERIC_REGULAR_AT_LINEARIZATION` and full initial GetFEM rollback.
The exact resting control retains Fa/Fd at numerical zero.

**Critical qualification boundary:** the rotating value
`||F_a||=0.250963115` is the residual of the **original zero
algebraic guess** at prescribed rotating `h0,p0`, not a minimized
compatibility residual. This pilot used *synthetic* right-hand
sides `A*x_ref` for rank diagnostics, **not the actual
dynamic/differentiated-constraint RHS**; the reported rotating
`||F_d||=0` is therefore not the physical force residual after
solving efforts and boundary controls. Numerical reconstruction
of two synthetic vectors is not proof of the differentiation index,
global invertibility or physical consistency of rotating
initial data. In particular, matrix invertibility at an
**incompatible** algebraic starting point cannot create a valid
DAE initial state.

The earlier original fixed-q rotating P2 reduced-LSQR algebraic
candidate reached a much smaller native `||F_a||≈4.18e-8`,
but missed the strict unscaled `2.60963e-9` coefficient
norm tolerance. This may involve numerical conditioning
and curved boundary trace incompatibility; **neither
physical compatibility nor nonexistence has been proven**.
The exact resting solution at machine precision previously
showed why a bare LSQR `INCOMPATIBLE_AT_TOLERANCE` label
is insufficient.

**NEXT:** audit the regular algebraic variables at fixed
original rotating q0 and the reduced boundary compatibility
in *proper physical norms*. Once and only once
`F_a(q0,a_r0)=0` is numerically justified, assemble
the real differentiated constraint RHS and solve for
`qdot0,a_reg_dot0,Y0` using the original dynamics
(gyro, viscosity, boundary law included). Keep q0
unchanged and require MPI 1/2 parity, rollback on
failure, no pseudo-TS, and subsequent physical energy
and boundary-power checks. #28 OPEN.


## 21 November 2026 — #28 ORIGINAL rotating disk: first genuine physical RHS / diagnostic Y0 pilot (NATIVE PENDING)

This follows the passed native P2 MPI 1/2 rank audit
(`3300 x 3300`, 132 traction Y variables, synthetic
MUMPS x-recovery at relative `~1e-10` or better).
**The rank audit did not compute the actual
physical forcing or a compatible rotating
algebraic initial state**. In particular
its `F_d=0` at rotation was evaluated at the
**original zero-algebraic guess**, not at
the physical constitutive efforts and
prescribed nonzero wall velocity.

Next source-only stage:
`sandbox/probe_v130_issue28_physical_rhs_y0.py`,
`scripts/run-v130-issue28-physical-rhs-y0.sh`,
`sandbox/compare_v130_issue28_physical_rhs_y0.py`
with host-only negative regression tests.
It uses the **unchanged original model**
`sandbox/viscous_swe_2d_mpi.py`,
real curved geometric P2 disk mesh `0.2`,
physical initial `q0=(h0,p0)`, corrected
clockwise GWFL convention, unmodified normal
and tangential boundary controls and full gyro
and viscous terms. There are **no accepted
physical or pseudo-initialization time steps**.

- **Known rest baseline (`omega=0`)**:
  initialize only the known analytical
  algebraic efforts `e_h=0.5`,
  `Y_n=-0.5`, with the exact original
  `h0=5,p0=0`, and require real
  GetFEM `F_a,F_d` at numerical zero.
- **Original rotating case (`omega=0.2`)**:
  hold original `h0,p0` and traction
  `Y=0` exactly fixed; compute only
  regular algebraic efforts/control fields
  `a_r0` with the true *rectangular*
  PETSc distributed LSQR, then evaluate
  the actual original `F_a` and
  per-block residuals. The previous
  comparable pilot got `||F_a||≈4.18e-8`:
  this does **not** establish exact
  discrete boundary compatibility.
- At that **projected** regular algebraic
  candidate, evaluate the **actual**
  differential residual `F_d(q0,a_r0,Y=0,zdot=0)`
  rather than taking the meaningless
  zero-algebraic-guess dynamic RHS.
  Independently finite-difference
  `partial_t F_a` within the actual
  wall ramp's flat initial interval.
- Assemble original GetFEM
  `Mdd, JdY, Jad, Jar` and solve
  the actual first-differentiated
  bordered system

  `Mdd qdot0 + JdY Y0 = -F_d`,

  `Jad qdot0 + Jar a_r_dot0 = -partial_t F_a`.

  Confirm the original full nonlinear
  `IFunction` at the resulting
  `qdot0,Y0` and the original
  `F_a(q0,a_r0,0)`. Never silently
  alter original `h0,p0` or the
  prescribed wall control, and restore
  all GetFEM/PETSc/cache state.
- Report actual `||Y0||_2` and
  `||qdot0||_2` as **DIAGNOSTIC
  ONLY**. Until independent,
  physically scaled rotating
  `F_a=0` compatibility is established,
  even a tiny dynamic residual does
  **NOT** qualify `Y0` as consistent.

The archived native evidence comparator
enforces this separation and supports
`--compare-only <saved-evidence-directory>`
to avoid rerunning MPI after a pure
host comparison failure.

**THIS PILOT IS SOURCE-ONLY. Native MPI1/2
execution and all Y0 candidate values PENDING.**
No claim of definite DAE index,
a physical high-index initializer,
a viable rotating trajectory,
or exact energy/power conservation.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
  bash scripts/run-v130-issue28-physical-rhs-y0.sh
```

The critical **next gate after seeing this
pilot** remains physical rotating
boundary-constraint compatibility,
including mesh-scaled norms of the
true discrete normal/tangential traces,
before *accepting* any computed Y0.

## 21 November 2026 — #28 true physical Y0 RHS: MPI2 MatNest ordering bug discovered and fixed in source (native PENDING)

**Failed original native diagnostic archive** `v130-issue28-physical-rhs-Y0-64eb6b2b-DB2mNSvT.zip`, source `64eb6b2b76007f4417c22bccd37df4f85e9d1ef5`, immutable rc3 OCI. The **exact rest reference passes**: native `||F_a||~1.47e-17`, `||Y0||~1e-15`, full dynamic residual `~1e-15`. With the **original rotating** P2 physical profile, both MPI ranks have identical original fixed-q initial states and nearly identical LSQR-projected `||F_a||~4.18e-8`, with actual physical `||F_d||~0.08654` at zero derivative and traction. But:

| Native rotating trial | MPI 1 | MPI 2 |
|---|---:|---:|
| Bordered linear KSP relative residual | `3.37e-13` | `1.68e-12` |
| Candidate `||qdot0||` | `6.57793` | `17.41334` |
| Candidate `||Y0||` | `0.0063778555` | `0.460843533` |
| **Original native GetFEM dynamic `||F_d(qdot0,Y0)||`** | **`1.65e-14`** | **`0.398585756`** |

The old runner stopped at the *correct* host safeguard: original GetFEM MPI2 `F_d` did not satisfy the actual PDE. **This is a real diagnostic MPI error**: small PETSc KSP residuals did not mean the correct physical RHS was assigned or the correct traction read.

**Source diagnosis and fix:** the old physical-RHS diagnostic creates a PETSc parallel `MatNest`, converts it to `AIJ`, then erroneously assumes the global converted row and column ordering is the serial block concatenation `[d; a]` and `[qdot; a_reg_dot; Y0]`. In MPI, MatNest block ranges follow its **native IS** permutations and can be interleaved across ranks; using `np.concatenate` for its RHS and slicing the globally gathered solution at fixed contiguous block offsets is invalid. The previous **synthetic** rank pilot was unaffected since it always inverted an internal `A*x_ref` without interpreting physics in GetFEM variable order.

The new native probe reads `matnest.getNestISs()` **before conversion** and uses its original parallel **row IS** to scatter the actual dynamic and differentiated-algebraic RHS. It checks a complete rank-local IS partition and a numerical RHS round-trip after assembly. It uses the **column IS**, plus the original compact subblock PETSc ownership ranges, to reconstruct `qdot0`, `a_regular_dot0`, and `Y0`, rather than naive global slices. It rejects trials if the injected candidate does not satisfy the **original full nonlinear native `IFunction` dynamics** at `1e-8`, even if PETSc's converted-matrix KSP reported convergence. The schema is incremented to **v2** so an old, incorrectly indexed archive cannot be silently relabelled as a success.

The comparator additionally checks **full physical qdot and Y0 coefficients**, not just the norm of Y0, with independent synthetic corruption unit tests for mapping and MPI parity. Original `h0,p0`, gyroscopic/viscous and boundary weak forms and frozen-gyro semantics are **unchanged**. No Y0 state is accepted, no pseudo/physical TS is performed.

**NATIVE RE-RUN PENDING:** `scripts/run-v130-issue28-physical-rhs-y0.sh`, MPI1/2 rest/rotating, from the current source. Only if the new *native* MPI2 original IFunction is below tolerance and true physical-vector parity passes can we qualify the **index-mapping diagnostic**. Even then rotating algebraic `||F_a||~4e-8` is still not a certified constraint-compatible initial DAE state; **no qualified rotating Y0 or time integration yet**.

## 21 November 2026 — #28 true physical-RHS/Y0 MPI 1/2 native gate PASSED (diagnostic, NOT initialized)

**Host evidence:** `sandbox/outputs/v130-issue28-physical-rhs-Y0-d937cdc6-pl0yI36l/`
(archived as `v130-issue28-physical-rhs-Y0-d937cdc6-pl0yI36l.zip`).
**Source:** `d937cdc6855b815fabc4108c4cad76d980271d3a`;
**immutable image:** `scrimp-mpi:v1.3.0-rc3`, OCI revision
`1d6d0c90115213298047b28237a7f38087563a35`.
The host comparator returns `COMPLETE_DIAGNOSTIC_NOT_INITIALIZER`.
All 7 host regression tests passed. All four genuine
PETSc/GetFEM MPI runs (rest/rotating × 1/2 ranks) completed.

Before this run, the probe incorrectly called
`matnest.convert("aij")` without `out=`, which performs a
petsc4py **in-place** conversion and destroys the native Nest IS;
this caused PETSc error 56 when `getNestISs()` was called.
The source fix at `e27c47f` uses
`matnest.convert("aij", out=PETSc.Mat())` and retains
the original Nest plus its row/column IS for parallel physical
RHS packaging and solution reconstruction. A host AST regression
test at `d937cdc` guards against repeating the in-place conversion.

| Real native diagnostic | Rest MPI1 | Rest MPI2 | Rotating MPI1 | Rotating MPI2 |
|---|---:|---:|---:|---:|
| Native projected `||Fa||_2` | 1.47634e-17 | 1.47443e-17 | 4.18508716e-8 | 4.18369529e-8 |
| Physical `||Fd(q0,a_r0,Y=0,zdot=0)||_2` | 2.69964e-16 | 2.69964e-16 | 0.0865426750 | 0.0865426757 |
| Solved `||qdot0||_2` | 1.25651e-14 | 1.25651e-14 | 6.57792793 | 6.57792791 |
| Diagnostic `||Y0||_2` | 2.25798e-16 | 2.25798e-16 | 0.00637785550 | 0.00637785651 |
| **Unmodified GetFEM native dynamic `||Fd(qdot0,Y0)||_2`** | **1.00534e-26** | **3.57982e-30** | **3.37799e-14** | **2.74326e-15** |
| Bordered relative linear residual | 1.56896e-26 | 2.52383e-28 | 1.46053e-13 | 2.85171e-14 |

The rotating LSQR projection used **5,149 / 5,154** iterations
(MPI1/2), with convergence reason 1, and reduced its unprojected
`||Fa||≈0.250963115` to `≈4.18e-8`. The residual is mainly
`e_p` (`4.10174e-8` MPI1, `4.10164e-8` MPI2), followed by
`U` (`≈6.72e-9`). This is **above** the previously imposed raw
fixed-mesh `2.60963e-9` coefficient-norm tolerance. It is not an
independent proof of exact discrete or physically scaled boundary
compatibility; it is also not sufficient to prove incompatibility,
given conditioning and projection accuracy.

The reconstructed full rotating `qdot0` and `Y0` coefficient
vectors agree between MPI1/2 to relative scales `6.00659e-7`
and `6.68769e-8`, respectively; the relative difference in
their `Y0` norms is `1.00769e-9`. Native MatNest RHS
round-trip errors are exactly zero for all four cases.
The prescribed physical `h0,p0` are unchanged bitwise,
and the GetFEM state is restored. No pseudo, BDF initialization
or true physical time step was executed. The previous bad MPI2
physical dynamic residual `0.398585756` has disappeared.

**QUALIFIED:** on this original P2 model and mesh, the true
post-projection physical RHS is packaged and solved consistently
with the native GetFEM/PETSc dynamics for MPI1/2, including
the exact rest control, correct MPI index mapping and MPI parity.
This **does not** prove first-differentiated full rank, generic
nonlinear DAE initialization, a consistent rotating `Y0`,
nor a physically admissible rotating trajectory.

**NEXT GATE (#28 remains OPEN, #6 nontrivial rotation blocked):**
at fixed original rotating `h0,p0`, independently measure actual
boundary normal/tangential velocity mismatch and interior
constitutive residuals using *physically scaled mesh-dependent
L2/dual norms* rather than only raw algebraic coefficient norms;
check residual sensitivity to LSQR stopping tolerance,
conditioning and mesh refinement. Separate numerical projection
error from any genuine incompatibility with the exact original
wall ramp and weak constraints. Preserve curved P2 geometry,
full viscous/gyro physics, fixed `h0,p0`, original controls,
rollback and MPI1/2 comparisons. If algebraic feasibility can be
established independently, only then promote the differentiated
solution toward an actual consistent initial state, followed by
a physical CN energy/port-power validation. Do not introduce
pseudo-time correction or alter the prescribed rotating state
merely to force acceptance.

The host logs still emit a petsc4py `ComplexWarning` when
reading the nominally real vectors into float arrays; a future
strict real/imaginary-value check should establish whether any
imaginary content exists instead of silently discarding it.
Unused PETSc BDF options arise because the probe intercepts the
initializer before TS; neither warning changes the gate result.

## 21 November 2026 — #28 physical constraint / LSQR-sensitivity gate staged (SOURCE ONLY)

**Next scientific prerequisite after the successful native physical-RHS MPI1/2
diagnostic**: distinguish an actual incompatibility of the *prescribed*
rotating momentum from inaccurate algebraic least squares and geometric
discretization. New independent, read-only probe
`sandbox/probe_v130_issue28_physical_compatibility.py`,
comparator `sandbox/compare_v130_issue28_physical_compatibility.py`,
negative controls `tests/test_v130_issue28_physical_compatibility.py`,
and pinned-Docker runner
`scripts/run-v130-issue28-physical-compatibility.sh`.

Using the unchanged `sandbox/viscous_swe_2d_mpi.py`,
`gamma=0`, `initial_profile=tangent`,
`wall_velocity_profile=facet_tangent`,
`geometry_order=height_order=2`, original clockwise GWFL
`Rot`, and `omega=0.2` (plus the `omega=0` rest control),
the probe evaluates **actual GetFEM integration-region scalar
integrals**, not the Euclidean norm of a PETSc RHS. At `t=0`
the original wall ramp is identically one. For
`v_0=p_0/rho` and `t=Rot*n`, it records:

- the domain area `|Omega|` and boundary length `|Gamma|`,
  for native MPI scalar assembly calibration;
- `||v_0.n||_{L2(Gamma)}` and
  `||v_0.t-omega||_{L2(Gamma)}`,
  plus the full `||v_0-omega*t||_{L2(Gamma)}`.
  These are measured at *prescribed* original `q0` **before
  any algebraic projection**; they are therefore independent
  of LSQR conditioning, `e_p` or `Y0`.
- at the analytic resting control, or projected rotating baseline,
  the physically scaled interior
  `||e_p-p/rho||_{L2(Omega)}`,
  `||sqrt(h)(e_p-p/rho)||_{L2(Omega)}`,
  `||e_h-rho*g*h-0.5*(p.p)/rho||_{L2(Omega)}`,
  and boundary `||e_p-U||_{L2(Gamma)}`,
  `||U-U_n*n-omega*t||_{L2(Gamma)}`, `||U_n||_{L2(Gamma)}`.
  Actual native algebraic `IFunction Fa` residual per block
  and `||A^T Fa||_2` are reported separately.

Each squared integral is archived with raw individual MPI-rank
GetFEM scalars and the raw sum. The RMS norm is normalized by
its **corresponding assembled domain/boundary measure**, so it
is invariant whether GetFEM scalar assembly is rank-local or
globally replicated. The host comparator independently identifies
which scalar MPI reduction convention applies by calibrating
area **and** perimeter against the MPI1 reference; it does
not silently sum a replicated global scalar. It checks the
vector/normal/tangent squared-norm decomposition and MPI1/2
physical RMS parity.

Physical mesh comparison: nominal `h=0.2` and `h=0.1`,
both curved P2, rest/rotating and MPI1/2. At `h=0.2`
the original rotating `h0,p0` and zero `Y` are fixed,
and regular algebraic variables are projected twice on the same
actual rectangular native PETSc system:
`LSQR(rtol=1e-12,atol=1e-13,max_it=10000)` and
`LSQR(rtol=1e-13,atol=1e-14,max_it=20000)`.
This is a **stopping-tolerance sensitivity test** that records
`||A^T Fa||`, individual native `Fa` blocks, correction
`||delta_a||` and the relative coefficient change between
the two LSQR tolerances; iteration counts alone are NOT a
spectral condition number estimate.
At refined `h=0.1` the rotating audit measures ONLY the
unprojected fixed-q boundary data, explicitly avoiding
conflation with a more expensive/different LSQR projection.

The harness rejects changed physical data, altered geometry,
nonzero exact-rest residual, bad MPI physical assembly, missing
independent trace diagnostics, and any claimed accepted
rotating `Y0` or time stepping. It intentionally does **not**
assert that mesh refinement must improve the mismatch or that
a small physical RMS implies exact weak DAE compatibility:
those are questions to answer from native evidence.
The physical FE trace diagnostic is not a proof of a unique
constraint-compatible initial state, nonlinear/high-index
initializer, KSP conditioning, or physical CN integration.

**SOURCE STAGED / NATIVE PENDING:** neither native GetFEM
`asm_generic` expressions nor MPI1/2 physical results have been
qualified by this source commit. No modification was made to
SCRIMP's production initializers, original SWE weak forms,
original `h0,p0`, the faceted wall law, or the image.

Run on the existing immutable image:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_physical_compatibility.py
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-physical-compatibility.sh
```

Use `bash scripts/run-v130-issue28-physical-compatibility.sh
--compare-only <saved-evidence-directory>` to replay archived evidence
without restarting expensive native runs. On completion, interpret
actual fixed-q normal/tangential trace errors, LSQR stopping
sensitivity and mesh trend *before* deciding whether rotating
constraints are numerically/physically feasible.
**#28 remains OPEN; nontrivial rotating #6 remains blocked.**

## 21 November 2026 — #28 original curved-P2 physical compatibility: NATIVE MPI 1/2 PASSED (TRACE DIAGNOSTIC ONLY)

Source `c1988939be5bd3836628295fe88ad6d59f385547`,
immutable `scrimp-mpi:v1.3.0-rc3` with OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Host evidence directory:
`sandbox/outputs/v130-issue28-physical-compatibility-c1988939-WqmfWd2x/`;
matching archive `v130-issue28-physical-compatibility-c1988939-WqmfWd2x.zip`.
**All 7 host-side negative/fixture tests passed** and all
8 genuine GetFEM runs (original rest/rotation, geometry P2,
nominal mesh 0.2/0.1, MPI1/2) completed. The host
comparator returned:
`COMPLETE_PHYSICAL_TRACE_DIAGNOSTIC_NOT_INITIALIZER`.

### Fixed original rotating physical momentum: native boundary RMS

| `omega=0.2`, independent GetFEM boundary measure | h=0.2, 3,300 DOFs | h=0.1, 10,272 DOFs | h-refinement factor |
|---|---:|---:|---:|
| `RMS((p0/rho).n)` | 1.3937491763e-5 | 2.0057694048e-6 | 0.1439117912 |
| `RMS((p0/rho).t-omega)` | 4.8200908963e-7 | 2.2071385559e-7 | 0.4579039282 |
| `RMS(p0/rho-omega*t)` | 1.3945824085e-5 | 2.0178764857e-6 | 0.1446939581 |
| `RMS((p0/rho).t)` | 0.1999995260 | 0.1999997794 | — |

Here `t=Rot*Normal` with the unchanged original clockwise GWFL
sign and wall law. MPI1 and MPI2 physically normalized RMS values
agree to rounding for **every** case; both rest meshes give exact
zero normal and tangential defects. Comparison against nominal mesh
halving gives descriptive two-point slopes `~2.797` (normal),
`~1.127` (tangential), `~2.789` (full vector), **not**
established convergence orders. These nonzero **strong**
boundary mismatch integrals are independent of the LSQR projection.
They are consistent with geometric/trace approximation effects
and show that the original rotating `q0` does not satisfy the
imposed trace **pointwise/L2-exactly** at these mesh sizes.
They do **NOT** by themselves prove failure of the *weak finite-
element* algebraic compatibility equations: an L2 nonzero field
can still have zero moments against the discrete test space.

**MPI scalar-integration convention independently checked:**
on all four case/mesh combinations GetFEM returns the same
**already global** scalar integral on both MPI processes.
The comparator correctly classifies MPI2 as
`replicated-global` after calibration against the actual
area and perimeter from MPI1. Never sum those two copies when
performing physical power, boundary or area audits.

### h=0.2 rotating algebraic LSQR sensitivity, actual native F_a

| Native quantity | MPI1 standard | MPI1 tight | MPI2 standard | MPI2 tight |
|---|---:|---:|---:|---:|
| LSQR iterations | 5171 | 5852 | 5101 | 5723 |
| `||Fa||_2` | 4.1794109732e-8 | 4.1573736753e-8 | 4.1929570788e-8 | 4.1576937032e-8 |
| `||A^T Fa||_2` | 5.0372e-11 | 5.4445e-12 | 4.3682e-11 | 5.3671e-12 |
| `||delta_a||_2` | 13.0655182281 | 13.0655182281 | 13.0655182283 | 13.0655182280 |
| relative `delta_a` change from standard | 0 | 3.01744e-7 | 0 | 3.90119e-7 |

Both KSP runs returned convergence reason 1. The smallest
rotating `Fa` remains `~4.157e-8`, **above**
the earlier unscaled `2.60963e-9` raw-coefficient threshold.
Tightening LSQR decreases the normal-equation gradient by
nearly an order of magnitude but only reduces `||Fa||`
by ~0.5%. The dominant tight residual block is `e_p`
(`~4.10055e-8`), followed by `U` (`~6.704e-9`).
This suggests a nonzero discrete least-squares floor,
but neither proves impossibility of exact weak algebraic
compatibility nor establishes the spectrum/conditioning
of the rectangular problem. Raw `Fa` Euclidean
coefficient norms must **not** be compared across
meshes without a consistent physical/dual scaling.

At the MPI1 tight projection, separately measured *strong*
physical RMS were `||e_p-p/rho||_Omega/√|Omega| =
1.36311987e-7`, weighted velocity constitutive RMS
`3.04802869e-7`,
`||e_p-U||_Gamma/√|Gamma|=6.10108e-10`,
`||U-U_n*n-omega*t||_Gamma/√|Gamma|=
1.39334856e-5`, `||U_n||_Gamma/√|Gamma|=
3.63294e-11`, and pressure constitutive RMS
`||e_h-rho*g*h-|p|²/(2rho)||_Omega/√|Omega|=
5.11671926e-4`.
The latter strong pressure error is **not**
itself evidence of an inconsistent *weak* constitutive
equation: the non-polynomial kinetic expression need
not belong to the P2 co-energy space even when its
Galerkin projection satisfies the test equations.

**QUALIFIED:** immutable native #28 diagnostic probe,
true original curved P2 strong boundary mismatch
and constitutive measures, calibrated MPI scalar
assembly, exact-rest positive control and MPI1/2
parity on both meshes, read-only original `h0,p0`.
There were no BDF/pseudo/physical time steps, no
accepted traction `Y0`, and no changed wall law.
The prior first-differentiated true-physical-RHS
MPI diagnostic remains separately qualified.

**NOT QUALIFIED:** weak discrete compatibility
`Fa(q0,a0)=0` for original rotation, unique consistent
`Y0`, differentiation index, physical rotating
time evolution or its Hamiltonian/port-power balance.
Original #28 and nontrivial rotating #6 remain OPEN.

**NEXT DISCRIMINATING GATE:** independently measure the
*moments of the prescribed wall mismatch* against the
actual GetFEM boundary test space, using a verified
MPI-consistent boundary Riesz mass matrix / dual norm
rather than inferring weak inconsistency from a
strong L2 mismatch. In parallel, check a *properly
scaled* h=0.1 regular-algebraic projection to separate
the stationary LSQR residual floor from the
geometry-dependent mismatch. Compare nonlinear
`IFunction Fa` blockwise after both projections,
including the interior co-energy weak projection.
Only if both algebraic feasibility and first-
differentiated consistency are independently verified
may a diagnostic `Y0` be promoted toward a DAE
initializer. Do not alter original `p0`, geometry
or imposed physical wall velocity to obtain a PASS.

## 21 November 2026 — #28 direct *weak* U-trace obstruction diagnostic staged (SOURCE ONLY)

The successful native **strong** physical-compatibility gate
`v130-issue28-physical-compatibility-c1988939-WqmfWd2x`
does **not** by itself establish weak finite-element inconsistency.
A nonzero `L2(Gamma)` strong mismatch can be orthogonal to the
actual discrete boundary test space. The tight LSQR residual
`~4.157e-8` provides circumstantial evidence only; it is
not a rigorously scaled or independently corroborated
minimum-distance certificate.

New read-only native diagnostic:
`sandbox/probe_v130_issue28_weak_wall.py`,
`sandbox/compare_v130_issue28_weak_wall.py`,
`scripts/run-v130-issue28-weak-wall.sh` and host
negative controls `tests/test_v130_issue28_weak_wall.py`.
It retains the original `h0=5,p0` tangent,
`omega=0.2` (with exact rest `omega=0` positive control),
P2 geometry/P2 unknowns, `gamma=0`, physical
`wall_velocity_profile=facet_tangent`, and pinned
immutable rc3 OCI. It runs on original nominal meshes
`h=0.2/0.1` and MPI1/2 with **no physical/initialization TS**,
LSQR or accepted traction.

At the *exact prescribed initial FE momentum*, assemble
with the actual GetFEM **order-one test** `Test_U`:

```text
b_full(v) = int_Gamma (p0/rho - omega*Tangent).v
b_n(v)    = int_Gamma ((p0/rho).Normal)*(Normal.v)
b_tau(v)  = int_Gamma (((p0/rho).Tangent)-omega)*(Tangent.v)
```

The probe checks `b_full=b_n+b_tau` in the actual vector
of U test-function moments. It also extracts the assembled
native `J_UU` block from the original `IFunction` tangent.
The physical control law contains `-U.Test_U`, so its
negative is the native boundary `U` mass matrix
`M_U=-J_UU`. Solve `M_U x=b` with MUMPS and report the
**discrete Riesz-projected physical boundary RMS**

```text
||Pi_U(p0/rho - omega*Tangent)||_{L2(Gamma)}/sqrt(|Gamma|)
  = sqrt(b^T M_U^{-1} b / |Gamma|).
```

For the exact rest state both strong and weak mismatch
should be numerically zero. In rotation, a nonzero
Riesz norm demonstrates nonvanishing moments against the
actual U test space and rules out the special explanation
that *all* of the original strong wall mismatch lies in
its discrete orthogonal complement. It still does NOT
alone establish full algebraic DAE incompatibility:
that conclusion additionally needs a verified statement
that exact interior `F_ep=0` enforces `e_p=p0/rho`,
exact `F_Y=0` enforces the matching U trace,
`F_Un=0` fixes `U_n=0`, and that the corresponding
FE/test/trace mass maps are nondegenerate, without any
unmodelled kernel.

The native comparator checks archived full first-order
U moment vectors, MPI rank assembly/reduction modes,
normal/tangential recomposition, serial-vs-parallel
coefficients, positive boundary mass Riesz solve,
Riesz RMS `<=` independently integrated strong RMS,
exact rest, bitwise read-only rollback, and strict
**NO Y0 INITIALIZATION**. It records diagnostic
moment size and mesh trend, never claims exact
DAE admissibility from a numerical threshold alone.
Host comparator fixtures are synthetic and must
not be confused with native quantitative results.

**SOURCE ONLY / NATIVE PENDING:** do not claim this
order-one GetFEM form, Riesz block, or MPI semantics have
yet been tested on the host. No change to production
SCRIMP solver, original weak bricks, initial state,
controls, or published release guarantees.

Run next native qualification:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_weak_wall.py
bash -n scripts/run-v130-issue28-weak-wall.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-weak-wall.sh
```

If the Python host comparator alone later needs a replay,
`bash scripts/run-v130-issue28-weak-wall.sh --compare-only
<saved-evidence-directory>` does not repeat GetFEM/MPI work.
Preserve original rotating `q0` and no physical time
integration: #28 and nontrivial rotating #6 remain OPEN.

## 21 November 2026 — #28 original weak U boundary moments: NATIVE MPI1/2 PASSED; exact elimination gate STAGED

**Native source:** `538f7e92bcda072594e6e13409a964c4bb4a415a`,
immutable `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Archived host evidence:
`sandbox/outputs/v130-issue28-weak-wall-538f7e92-86jIvI0n/`.
All six host regression tests passed and all eight native cases
(rest/rotation × P2 mesh nominal 0.2/0.1 × MPI 1/2) passed,
followed by the independently archived comparator:
`COMPLETE_WEAK_WALL_DIAGNOSTIC_NOT_INITIALIZER`.

The original `p0`, `h0=5`, `gamma=0` and prescribed
`U=U_n*n+omega*Rot*n` with `U_n=0` are unchanged.
For the actual GetFEM port-U functions `Test_U`, integrate
`b(v)=int_Gamma (p0/rho - omega*Tangent).v`.
Use the native `M_U=-J_UU` to measure
`sqrt(b^T M_U^{-1}b/|Gamma|)` (full and normal/tangential
components). No accepted Y0, no pseudo/BDF or physical time step.

| Native original rotating `omega=.2` | h=0.2, U=132 DOFs | h=0.1, U=252 DOFs |
|---|---:|---:|
| Strong wall mismatch RMS | 1.394582408498e-5 | 2.017876485655e-6 |
| **Native weak Riesz RMS** | **5.867041664403e-7** | **2.221552966076e-7** |
| Normal weak Riesz RMS | 2.598473916159e-7 | 1.969987243180e-8 |
| Tangential weak Riesz RMS | 4.819480885555e-7 | 2.207136474009e-7 |
| Full GetFEM wall moment coefficient `||b||_2` | 3.679724232681e-7 | 1.276350903720e-7 |
| `h0*||b||_2` expected Y-equation coefficient defect (conditional) | 1.839862116341e-6 | 6.381754518600e-7 |

Both MPI configurations have identical native vector moments to
printing precision; the parallel native order-one vectors and
scalar perimeter are **replicated globally**. The independent
full/normal/tangent moment decomposition errors are
`2.2213e-17` (coarse) and `1.3382e-17` (fine).
Riesz mass solves have relative residuals `<=9e-23`.
At rest, strong and weak residuals vanish exactly on both meshes.

**Qualified conclusion:** the **original prescribed rotating
momentum has a nonzero weak mismatch against the actual finite-
element wall U test space**, not only a strong L2 mismatch
orthogonal to that space. Its Riesz RMS falls only by a factor
`~2.641` between these nominal meshes; the weak tangential
component dominates on the fine mesh. The initial expression
`omega*(y,-x)/(sqrt(x*x+y*y)+1e-6)` includes an intentional
`1e-6` center regularization that also slightly changes the
trace near `r=1`; the contribution of this regularization
versus geometric interpolation has **not** been isolated.
Two mesh sizes do not establish asymptotic convergence.

**Important remaining logical step:** a nonzero `b` alone
does not prove the **full original algebraic system** inconsistent.
At constant `h0=5`, the original `p` and `e_p` FEM are
identical P2 spaces, so the interior weak constitutive
constraint `F_ep=0` should force `e_p=p0/rho` by the
nondegenerate h-weighted mass. Original `F_U=0` should
force the port U coefficients via its boundary mass
`M_U=-J_UU`. If `F_Y` in this unique joint solution is
`+/-h0*b != 0`, exact simultaneous satisfaction of
`F_ep,F_U,F_Y` is impossible for that fixed rotating q0
on the tested meshes. This is a **conditional finite-dimensional
block elimination argument**, not a generic DAE index theorem.

To verify this chain independently in the actual **unmodified
GetFEM IFunction**, new source-only diagnostic
`sandbox/probe_v130_issue28_exact_weak_elimination.py`
copies original `p0/rho` FE coefficients directly to
`e_p`, solves the native `F_U` with `-J_UU` (MUMPS),
checks native `F_ep=0`, `F_U=0`, and compares actual
`F_Y` coefficient by coefficient to the independently
assembled signed `h0*b`. Its companion
`sandbox/compare_v130_issue28_exact_weak_elimination.py`
checks MPI1/2 parity, original fixed q0, the previously
frozen full `||b||_2` values, positivity of the rotating
defect, and zero-rest positive controls. Paired host
negative tests:
`tests/test_v130_issue28_exact_weak_elimination.py`.
The pinned-Docker runner is
`scripts/run-v130-issue28-exact-weak-elimination.sh`.
Source changes do not touch the production initializers
or original PDE. **SOURCE STAGED / NATIVE NOT YET RUN**;
the block identity remains a hypothesis pending this gate.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_exact_weak_elimination.py
bash -n scripts/run-v130-issue28-exact-weak-elimination.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-exact-weak-elimination.sh
```

Even if the three native blocks prove inconsistent at original
rotating q0, do **not** change the original prescribed
initial momentum or physical wall control to force acceptance.
This evidence would justify explicitly rejecting that
original incompatible DAE initial condition and designing a
**separate, explicitly specified compatible physical state**
or a wall spin-up from validated rest for #6.
Never accept physical rotating Y0 or run a physical CN/energy
campaign until a consistent initial DAE candidate is actually
verified. #28 remains OPEN; nontrivial rotating #6 remains blocked.

### Analytical circular-wall limit of the original regularized tangent profile

There is an additional **exact continuous-data observation**, distinct
from the unqualified discrete FE elimination. The original prescribed
momentum is
`p0/rho = omega*(y,-x)/(sqrt(x*x+y*y)+epsilon)`, with
`epsilon=1e-6`. On the **ideal unit circle** (`r=1`,
`Tangent=(y,-x)`), the imposed wall is `omega*Tangent`,
so the original continuous initial trace differs by

```text
(p0/rho - U_wall)|_{r=1}
  = -omega*epsilon/(1+epsilon) * Tangent,
omega=.2, epsilon=1e-6
  => tangential deficit = 1.999998000002e-7.
```

Thus geometric P2 refinement alone cannot remove this exact
regularization-induced tangential mismatch **in the continuum
reference expression**: even a perfect circle retains it.
The measured fine-mesh native tangential weak RMS
`2.207136474e-7` is consistent with that scale, but
its decomposition into regularization versus FE interpolation
is not yet certified. The full original FE weak algebraic
incompatibility is still being checked independently by
the staged exact `F_ep/F_U/F_Y` block gate. Do **not**
silently alter `epsilon` or `p0` in the original
validation campaign; an explicitly new compatible-data
experiment must be separately labelled and qualified.

## 21 November 2026 — #28 exact F_ep/F_U elimination: native 6/8 passed; h=.1 floating-roundoff false rejection fixed

**Archived partial native campaign:**
`sandbox/outputs/v130-issue28-exact-weak-elimination-a28aee14-p8XVIVz0/`,
original source `a28aee142b2d1798e225535f9dab217a91938a0c`,
immutable `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
**Six of eight native results completed and are preserved**, all
six host-only comparator fixture tests passed. The original
analytic resting state passes MPI1/2 at both mesh 0.2 and 0.1.
The original rotating P2 mesh 0.2 passes native MPI1/2:

| Native diagnostic at fixed original rotating h0=5,p0, zero Y | h=.2 MPI1 | h=.2 MPI2 |
|---|---:|---:|
| `||F_ep||_2` after exact coefficient copy and U elimination | 5.91914e-18 | 5.91914e-18 |
| `||F_U||_2` after unmodified wall mass solve | 2.75512e-17 | 2.91411e-17 |
| **`||F_Y||_2`** | **1.839862116389e-6** | **1.839862116390e-6** |
| **Independent `h0*||b||_2`** | **1.839862116341e-6** | **1.839862116341e-6** |
| `||F_Y-h0*b||_2` | 1.58918e-16 | 1.70371e-16 |

The *seventh* native case, original rotating **h=0.1
MPI1**, computed the expected nonzero
`h0*||b||_2=6.38175452e-7` and
an independent `F_Y` agreeing **to
`1.09940161e-16` coefficient-norm error**, but
the old source incorrectly rejected the identity at
`error > 1e-10 * max(1e-6, ||h0*b||)`,
which is precisely `1e-16` on this small load.
The error exceeded the fixed arbitrary threshold by
`9.94e-18`; there is **no PETSc/GetFEM/MUMPS solver
failure** and no evidence of a different physical
solution. The seventh run did not write a JSON report;
the eighth, h=.1 rotating MPI2, did not start because
the shell runner correctly stopped on the false rejection.

**Source-only rounding fix:** the independent native
`F_Y=+/-h0*b` comparison now uses a
strict floating-roundoff-aware bound

```text
||F_Y - sign*h0*b||_2 <=
    5e-15 + 2e-9 * ||h0*b||_2
```

and archives the *actual* absolute and relative
coefficient errors and the tested tolerance. The
host comparator independently enforces this same
bound and remains backwards-compatible with the
six original native reports, which were verified
against the new bound *from their original saved
coefficient vectors*, not fabricated. Two new
host negative controls distinguish a valid
`~1.10e-16` rounding discrepancy from
an intentionally corrupted `2e-12` identity.

**Recover, do not repeat:** the native runner
`scripts/run-v130-issue28-exact-weak-elimination.sh`
now accepts `--resume <existing evidence directory>`.
It checks the immutable OCI revision, the exact
old/current source provenance, the existing original
manifest and its no-time-step/no-Y0 invariants;
each old report is independently validated before
it is **skipped**. Only the two missing h=.1 rotating
MPI1/2 native cases are rerun, with current source.
The original manifest `source` remains untouched;
`resume_metadata.json` separately records the
original versus resumed source commits. No
numerical result is silently relabelled as a
new native execution.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_exact_weak_elimination.py
bash -n scripts/run-v130-issue28-exact-weak-elimination.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-exact-weak-elimination.sh \
  --resume sandbox/outputs/v130-issue28-exact-weak-elimination-a28aee14-p8XVIVz0
```

**NATIVE RE-RUN PENDING** for the last two cases.
The pre-existing h=.2 rotating evidence provides
strong native support for **incompatibility at
fixed original q0** of the original
`F_ep,F_U,F_Y` block triple: exact P2
`e_p=p0/rho` and original control-port U
elimination leave a nonzero independently
verified weak `F_Y=h0*b` remainder.
The continuum regularized original tangent
profile has the independent nonzero
`~2e-7` tangential wall deficit even on
the ideal circle. These facts **do not**
qualify a general DAE index, a physical Y0
for this inconsistent initial state, nor a
time-integrated rotating trajectory.
Do not alter original h0,p0/epsilon,
wall law, geometric P2, or run any TS.
#28 and rotating #6 remain OPEN until the full
native archive and an explicitly compatible
initial-condition strategy are reviewed.

## 21 November 2026 — #28 native 8/8 qualified: exact weak block obstruction at the prescribed rotating q0

**FINAL NATIVE GATE PASSED:** original finite-element curved-P2 disk, original
`h0=5`, original regularized `p0/rho = omega*(y,-x)/(r+1e-6)`,
original imposed rotating wall `U=U_n*n+omega*Rot*n` and
`U_n=0`, `omega=.2` (plus exact-rest `omega=0` control),
nominal mesh `h=.2/.1`, MPI1/2, immutable
`scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Completed mixed-provenance (strictly preserved) native evidence
`sandbox/outputs/v130-issue28-exact-weak-elimination-a28aee14-p8XVIVz0/`.
Six validated native reports originate from source
`a28aee142b2d1798e225535f9dab217a91938a0c`;
the two `h=.1` rotating reports originate from resumed source
`9e4dd40ec36482484641ae35fc510df9bf3c6e7b`,
with **only** the absolute/relative roundoff guard corrected.
The original manifest was not relabelled, and
`resume_metadata.json` records the resumed source.
All **eight** native reports now pass the independent host comparator,
and all **eight** host regression/negative controls pass.

### Certified original weak block elimination (fixed q0)

At constant `h0=5`, `p` and `e_p` use the same P2
FE basis and the `h0`-weighted positive-definite interior
mass. Therefore the native interior constitutive block
`F_ep=0` uniquely implies `e_p=p0/rho`;
this was implemented by direct coefficient copy, **without
altering original p0**. Solving the original native
`F_U=0` equation using its boundary mass `M_U=-J_UU`
fixes U uniquely on the filtered U trace space. The
original `Y`-row residual after these two eliminations
is then independently identified coefficient-by-coefficient
as the unchanged weak wall-mismatch functional:

```text
b(v) = int_Gamma (p0/rho - omega*Rot*Normal).v,
F_Y = sign * h0 * b,  sign validated from actual IFunction.
```

| Original diagnostic | h=.2 rotating MPI1 | h=.2 rotating MPI2 | h=.1 rotating MPI1 | h=.1 rotating MPI2 |
|---|---:|---:|---:|---:|
| `||F_ep||_2` | 5.91914125e-18 | 5.91914125e-18 | 3.21638724e-18 | 3.21638724e-18 |
| `||F_U||_2` | 2.75511808e-17 | 2.91411181e-17 | 2.10639674e-17 | 2.14560291e-17 |
| **`||F_Y||_2`** | **1.83986212e-6** | **1.83986212e-6** | **6.38175452e-7** | **6.38175452e-7** |
| independent `h0||b||_2` | 1.83986212e-6 | 1.83986212e-6 | 6.38175452e-7 | 6.38175452e-7 |
| `||F_Y-sign*h0*b||_2` | 1.58918079e-16 | 1.70370660e-16 | 1.09940161e-16 | 1.15389289e-16 |

All four rest/MPI/mesh control cases have
`F_ep=F_U=F_Y=0` exactly. No numerical iteration
or nonlinear mass approximation was needed for
`e_p=p0/rho`; U was solved with the native
boundary mass. The independent weak wall
moment agrees with the unchanged native IFunction
to absolute `~1e-16`, while the rotating residual
is `~1e-6..1e-7`. The verified
roundoff criterion is
`||F_Y-sign*h0*b|| <= 5e-15+2e-9*||h0*b||`;
the former `1e-16` accidental false failure
at `h=.1` is resolved.

**VALIDATED MATHEMATICAL/NATIVE CONCLUSION (scope restricted
to this original discretization):** for the **fixed
original rotating h0,p0**, there is no simultaneous zero
of the three original algebraic constraint blocks
`F_ep,F_U,F_Y` on the two P2 meshes considered,
given the verified nondegenerate same-basis interior
and boundary mass operators. The contradiction
is **independent of LSQR stopping tolerances** and
is not a PETSc/MUMPS breakdown; it arises from
the prescribed initial data and wall actuation.
Even the ideal geometric unit circle retains
a `~1.999998e-7` tangential mismatch because the
original regularization uses `r+1e-6`.

**NOT qualified:** a global proof for arbitrary meshes
or initial states, general nonlinear DAE differentiation
index, physical rotating `Y0`, physical rotating
CN time evolution, or Hamiltonian/port-power balance.
No accepted initial state, no TS initialization or
physical time step was performed. Production SWE
weak bricks, geometry, physical initial values, and
SCRIMP production initializer were **not changed**.

### Implementation decision after this completed gate

The *scientific* obstruction for the **original** rotating
q0 is now resolved: do **not** continue iterating LSQR
or force a physical `Y0` on inconsistent original data.
The software issue #28 remains **OPEN** until SCRIMP
has a safe, explicit and verified numerical pathway:

1. **Failure semantics:** fail-fast informative diagnostics
   for incompatible prescribed differential data rather
   than accepting an ill-defined algebraic Y0 or an
   initialization transient. Develop a general interface
   and tests using the original rotating model as the
   negative control, the exact rest state as the
   positive control; do not unconditionally hard-code
   SWE-specific tolerances or turn the mass diagnostic
   into an undocumented solver default.
2. **Compatible experiment, deliberately separate:**
   construct explicitly specified physically compatible
   initial data and wall actuation (e.g. fully validated
   rest-to-rotation ramp, or a separately labelled
   constrained FE projection of rotating momentum),
   with exact all-block `Fa=0` **before** a differentiated
   physical RHS/Y0 consistency test. Never overwrite
   original `p0` or smooth wall law silently.
3. **Once accepted:** recheck true physical RHS
   (including `qdot0,Y0`), MPI parity, then validate
   actual physical Crank–Nicolson power/energy
   balances and eventually dissipative rotating
   Saint-Venant #6.

**Release status:** original rotating q0 algebraic
incompatibility is **qualified**; no generic
initialization or physical rotation is qualified.
#28 remains OPEN as an implementation item.
Nontrivial rotating #6 is still BLOCKED.

## 21 November 2026 — #28/#6 dual native startup gate staged: typed high-index refusal and compatible rest-to-spinup t0 (SOURCE ONLY)

**Context:** the independent original P2 weak-block elimination native
campaign completed 8/8 tests (source `a28aee1` + roundoff-only
resumption `9e4dd40`). For prescribed rotating `h0=5,p0`,
the interior `F_ep=0` and original `F_U=0` are uniquely
eliminated, but `F_Y=+/-h0*b` remains nonzero; neither
LSQR convergence tweaks nor pseudo-time can make those
prescribed data algebraically feasible. This original
negative-control evidence remains unchanged.

**Two separate goals, neither a new production algorithm:**

1. **Existing software refusal.** Exercise the actual
   `DPHS.init_step` dispatch with explicit
   `init_step_ts_type=nonlinear_affine` on the *unchanged*
   original rotating `initial_profile=tangent` and
   `wall_velocity_profile=facet_tangent` (f(0)=1).
   The production `nonlinear_affine` initializer already
   has an opt-in `UnsupportedNonlinearAlgebraicDAE`
   guard for algebraic null columns (higher-index
   multipliers). The native test demands this exact
   **typed structural refusal**, a report with
   positive `algebraic_null_columns`,
   both GetFEM and PETSc GlobalVec rolled back
   bitwise, and **NO** pseudo/physical TS. This is
   a structural unsupported-class rejection, NOT
   a numerical certificate of incompatible q0 by itself:
   the earlier exact elimination supplies the
   independent mathematical incompatibility argument.
   The new test does **not** change SCRIMP's default
   initializer, nor does it assert that the nonlinear-affine
   guard solves the general #28 problem.

2. **Distinct compatible physical experiment at t0.**
   Keep the *same* curved P2 disk, original
   dissipative/gyro/viscous SWE bricks, `gamma=0`,
   `h0=5`, original physical
   `wall_velocity_profile=facet_tangent`,
   and `rotation=0.2` after ramp.
   Use the already implemented, opt-in
   `initial_profile=rigid` with
   `start_from_rest_duration=0.02` to prescribe
   **p0=0** and a genuinely stationary wall at t0.
   The physical wall actuation is the original
   decreasing envelope times a smoothstep rising
   `g(t)=3s²-2s³`, `s=clip(t/.02,0,1)`,
   with `g(0)=g'(0)=0`. This is an explicitly
   **new initial-value experiment**; it must never
   be conflated with or overwrite the incompatible
   original rotating experiment. The analytic
   t0 algebraic candidate is:
   ```text
   e_h = rho*g*h0 = 0.5
   e_p = 0
   Y_n = -rho*g*h0 = -0.5
   U_n = U = Y = 0
   z_dot(t0) = 0
   ```
   Evaluate **the original production IFunction** on
   this full candidate at t=0 and require **both**
   the differential and algebraic rows to have
   near-machine-zero residual, then check
   native `F_U(delta)`, `F_U(2delta)` on
   the same stationary frozen state for
   `delta=1e-4` so the ratio agrees with the
   exact `g(2delta)/g(delta)` of the wall
   smoothstep. This independently verifies the
   absence of an instantaneous boundary ramp impulse
   in the original GetFEM control law.
   This is **t0 compatibility only**:
   no differentiated Y0 solve, no nonlinear/high-index
   production initializer, and NO physical time step.

**Source staged now**, pending all native executions:

- `sandbox/probe_v130_issue28_startup_gate.py`;
- `sandbox/compare_v130_issue28_startup_gate.py`;
- `scripts/run-v130-issue28-startup-gate.sh`;
- `tests/test_v130_issue28_startup_gate.py`.

Harness keeps image `scrimp-mpi:v1.3.0-rc3`
and checks immutable OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
It runs **8** native cases:
rest-spinup and original-incompatible × nominal
P2 mesh `h=.2,.1` × MPI 1,2, and records native
reports, original source SHA, image revision,
host-owned working directory, and strict
independent host comparator. Host test fixtures
are **synthetic**, never native evidence.
Existing production `dphs.py`, default
linear-only `algebraic`, opt-in `nonlinear_affine`
and all SWE bricks remain UNCHANGED.

Run on the host:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_startup_gate.py
bash -n scripts/run-v130-issue28-startup-gate.sh

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-startup-gate.sh
```

**Qualification gates after this diagnostic:** if all native
checks pass, we may claim the existing opt-in
high-index structural fail-fast and the new spinup
full IFunction t0 positive-control candidate.
Neither finding is yet a general consistent-state
initializer. Before the first physical CN step,
separately verify the full differentiated DAE
initial equations including qdot0 and traction Y0,
and explicit energy / port-power balances.
#28 remains OPEN; rotating #6 remains BLOCKED.

## 21 November 2026 — #28/#6 native dual startup gate: **8/8 PASSED** (compatible spin-up t0 + production typed rejection)

**Native qualification:** source
`fee8626b4634070e4e20686bc7033d9795351e20`,
image `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`,
archived
`sandbox/outputs/v130-issue28-startup-gate-fee8626b-0LcChumB/`.
All eight host-only negative/fixture tests PASSED; all
eight genuine GetFEM/PETSc probes PASSED
(two cases × nominal mesh 0.2/0.1 × MPI1/2),
and independent comparator returned
`NATIVE_STARTUP_SAFETY_AND_T0_FEASIBILITY_DIAGNOSTIC_ONLY`.
**DO NOT** confuse synthetic test displays
(`Fa=1e-17,Fd=1e-17,null_cols=8`)
with actual native measurements below.

The deliberately **new**, physically distinct open
Saint-Venant experiment (`h0=5,p0=0`,
P2 curved disk, original physical
facet-tangential wall, `omega=0.2`,
`start_from_rest_duration=.02`,
original viscous/gyro bricks) has
the analytic t0 algebraic candidate
`e_h=0.5, Y_n=-0.5, e_p=U=U_n=Y=0`
and `qdot0=0`.
Native full original GetFEM `IFunction`
with `zdot=0` evaluates:

| Positive spin-up t0 | h=.2 MPI1 | h=.2 MPI2 | h=.1 MPI1 | h=.1 MPI2 |
|---|---:|---:|---:|---:|
| `||Fa||_2` | 1.4763422151e-17 | 1.4744250641e-17 | 8.0033824547e-18 | 8.0057171923e-18 |
| `||Fd||_2` | 2.6996404680e-16 | 2.6996404680e-16 | 2.6280928039e-16 | 2.6280928039e-16 |
| `||F_U(.0002)|| / ||F_U(.0001)||` | 3.98662207358 | 3.98662207358 | 3.98662207358 | 3.98662207358 |

The physical rising smoothstep
`g(s)=3s²-2s³` with `s=t/.02`
satisfies `g(0)=g'(0)=0`,
and the **native wall form** ratio exactly
matches its two-time analytical ratio.
The full original `F(t0,z0,zdot0=0)`
is thus numerically zero in **both dynamic
and algebraic rows**, with MPI1/2 parity.
This qualifies a **t0-compatible FE state**,
not a computed high-index index analysis,
nor differentiated `Y0`, nor physical
time integration.

For the **unchanged incompatible original**
already-rotating tangent data (f(0)=1),
the existing opt-in production
`init_step_ts_type=nonlinear_affine`
path was called through the **real**
`DPHS.init_step` dispatch and threw
the expected typed
`UnsupportedNonlinearAlgebraicDAE`
on null algebraic columns:
132 (h=.2), 252 (h=.1),
identical MPI1/2. The native
`initialization_report.status` is
`REJECTED`, all original coefficients
and `GlobalVec` are restored bitwise,
and no pseudo/BDF/physical TS executed.
This is a **structural refusal of a
higher-index class**, NOT a generalized
runtime physical infeasibility diagnosis;
the earlier exact FE elimination is
the independent proof for the original
rotating q0. Importantly the same
structural limitation also applies to
some **physically compatible** states;
the high-index initializer remains
unimplemented.

**Next stringent gate:** verify the
first-differentiated DAE constraint system
at the **new compatible spin-up t0** with
correct true physical `-Fd` and correct
`-partial_t Fa`. Critical numerical
caveat: a naive finite difference of the
C1 smoothstep at its clipped start
`g(0)=g'(0)=0` produces a spurious
nonzero `O(dt)` derivative despite
the exact derivative being zero.
Use verified analytic `partial_t Fa(0)=0`
and an independent one-sided `dt`
convergence check rather than
reusing an unqualified fixed central
difference. Then compare MPI1/2
the genuine differentiated bordered
GetFEM/PETSc solve for `qdot0,Y0`,
and recheck native full `IFunction`
and differential-state preservation
before accepting an initial state.
**NO physical TS until this check passes.**
#28 remains OPEN as an actual high-index
production implementation; nontrivial
physical rotating #6 remains blocked.

## 21 November 2026 — #28/#6 true first-differentiated compatible C1 spinup t0 gate staged (SOURCE ONLY)

**Previously NATIVELY PASSED:** the separately specified
physical rest-to-rotation spinup state (h0=5,
p0=0, curved P2, original facet-tangent wall,
viscosity and gyroscopic/port-Hamiltonian
forms preserved, omega=.2 following
`start_from_rest_duration=.02`) satisfies
the original full `IFunction` at
`t0=0` to ~1e-16 on both mesh .2/.1,
MPI1/2 (8/8 previous dual-startup cases).
The existing opt-in `nonlinear_affine`
initializer separately rejects the old
already-rotating high-index negative control
with 132 / 252 null algebraic columns
and exact rollback. Those are **two
distinct results**. The old rotating fixed-q0
data were earlier proven algebraically
incompatible by exact `F_ep,F_U,F_Y`
elimination and remain untouched.

**New diagnostic now STAGED (native not yet run)**:
`sandbox/probe_v130_issue28_spinup_first_derivative.py`,
`sandbox/compare_v130_issue28_spinup_first_derivative.py`,
`scripts/run-v130-issue28-spinup-first-derivative.sh`,
`tests/test_v130_issue28_spinup_first_derivative.py`.

This is a deliberately **new** physically compatible initial
condition, NOT a hidden projection of the original rotating
momentum. The exact t0 candidate uses
`e_h=rho*g*h0=.5`,
`Y_n=-rho*g*h0=-.5`,
`e_p=U_n=U=Y=0` and `p0=0`.
It does NOT run the production
`nonlinear_affine` initializer,
which correctly declares this
high-index structure unsupported.
It reuses and adapts the already
natively qualified PETSc first-
differentiated physical-RHS MatNest
diagnostic, with the crucial
MPI native block IS remapping.

The first differentiated candidate is computed from

```text
  [ M_dd    0     J_dY ] [ qdot0  ]   [ -F_d(t0,z0,0)  ]
  [ J_ad   J_ar   0   ] [ a_dot0 ] = [ -partial_t F_a ]
                           [ Y0    ]
```

where `a_regular` excludes
the dynamic-only traction `Y`.
GetFEM's actual `M_dd,J_dY,J_ad,J_ar`,
original nonlinear IFunction
and MPI-aware PETSc MatNest/AIJ
are used; independent native
`IFunction(t0,z_trial,zdot_trial)`
must confirm algebraic and
dynamic residuals. No artificial
zero forcing is substituted for
the actual dynamic `-F_d`.

**VERY IMPORTANT C1 time derivative:** at the
new physical control `f(t)=f_decreasing(t)*g(t)`
with `g(t)=3s²-2s³`,
`s=clip(t/.02,0,1)`,
`g'(0)=0` **exactly** and
`f_decreasing(t)=1` for
`0<=t<.03`.
Consequently
`partial_t F_a(t0,z0,0)=0`
analytically for this model.
A naive finite difference with fixed
`dt` at the clipped initial time
would introduce spurious `O(dt)`
time forcing into the differentiated
constraints and a false `Y0`.
The new probe therefore uses
the mathematically exact value zero
**and independently verifies** the
actual native GetFEM `U`-block
time response at `delta=1e-4`
and `2delta` matches the
exact `g(2delta)/g(delta)`
to ~2e-6 relative; all other rows
must remain time invariant under
that frozen original t0 state.

Native checks (planned):
- nominal mesh .2 and .1, MPI1/2,
  pinned immutable `scrimp-mpi:v1.3.0-rc3`
  OCI `1d6d0c90115213298047b28237a7f38087563a35`;
- exact t0 algebraic `Fa<1e-11`,
  physical dynamic forcing below 1e-10,
  original nonlinear/P2 discrete mass
  and stiffness, true `J_dY` block
  and zero algebraic mass rows;
- MatNest global MPI row+column
  index partition and RHS roundtrip,
  direct MUMPS convergence and
  independently verified original
  `IFunction` with candidate
  `qdot0,Y0`;
- candidate norms `qdot0,Y0,a_dot0<1e-8`,
  full original dynamic residual
  <1e-10 and algebraic <1e-11;
- MPI1/2 physical and coefficient
  parity, bitwise preservation of
  differential q0 and rollback,
  no physical or pseudo/BDF TS,
  **no accepted physical Y0**
  and no generic index or
  high-index production claim.

The host fixture tests are synthetic
and no substitute for native GetFEM.
This diagnostic does **not** establish
a globally nonsingular first
differentiated operator for all
models, nor an accepted production
state or a physical Crank–Nicolson
energy/port-power trajectory.

Run:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_spinup_first_derivative.py
bash -n scripts/run-v130-issue28-spinup-first-derivative.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-spinup-first-derivative.sh
```

**Source-only pending native proof.**
Next if this pilot passes: choose
a clearly scoped high-index
consistent-initial-state
implementation strategy, distinguish
it from the existing
`nonlinear_affine` low-index
initializer; validate full rollback
and negative incompatible q0.
Only then start an actual physical
CN step, with Hamiltonian, boundary
port work and viscous dissipation
power-balance checks.
#28 OPEN, #6 nontrivial rotation BLOCKED.

## 21 November 2026 — #28/#6 compatible rest-spinup: physical differentiated t0 MatNest MPI 1/2 **NATIVE 4/4 PASSED**

**Native archive:** `sandbox/outputs/v130-issue28-spinup-first-differentiated-54bafa2c-y1es00lx/`,
original diagnostic source `54bafa2c10b50eb195d92c6a4b00b53970487d5a`,
immutable `scrimp-mpi:v1.3.0-rc3`, OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
All four native runs and the strict full-vector MPI comparator **PASSED**:
two curved-P2 meshes h=.2/.1 × MPI1/2.
The archived native comparator status is
`COMPLETE_DIFFERENTIATED_SPINUP_T0_DIAGNOSTIC_NO_TS`.

The *deliberately new* physical rest-to-rotation experiment is distinct
from the original incompatible rotating q0. Here `h0=5, p0=0`,
the original nonlinear conservative/viscous Saint-Venant weak bricks
are intact, `gamma=0`, and the **physical facet-tangent boundary
velocity** starts from rest with a C1 rising smoothstep
`start_from_rest_duration=.02`; `g(0)=g'(0)=0`.
At `t0` the exact initial algebraic candidate is
`e_h=.5, Y_n=-.5, e_p=U_n=U=Y=0`.

**Actual native values (NOT synthetic host fixtures):**

| Native quantities | h=.2 MPI1 | h=.2 MPI2 | h=.1 MPI1 | h=.1 MPI2 |
|---|---:|---:|---:|---:|
| `||F_a(t0,z0,0)||_2` | 1.476342215e-17 | 1.474425064e-17 | 8.003382455e-18 | 8.005717192e-18 |
| `||-F_d(t0,z0,0)||_2` | 2.699640468e-16 | 2.699640468e-16 | 2.628092804e-16 | 2.628092804e-16 |
| **`||candidate qdot0||_2`** | 1.256514249e-14 | 1.256514249e-14 | 4.150424431e-14 | 4.150424433e-14 |
| `||candidate regular a_dot0||_2` | 1.256514249e-14 | 1.256514249e-14 | 4.150424432e-14 | 4.150424433e-14 |
| **`||candidate Y0||_2`** | 2.257984145e-16 | 2.257984145e-16 | 3.044205831e-16 | 3.044205831e-16 |
| `||F_d(t0,z0,qdot0,Y0)||_2` | 1.005343326e-26 | 3.579807531e-30 | 3.335370174e-24 | 1.139478072e-29 |
| `||F_a(t0,z0,qdot0,Y0)||_2` | 1.476342215e-17 | 1.474425064e-17 | 8.003382455e-18 | 8.005717192e-18 |

Native direct PETSc/MUMPS bordered solve reports converged reason 4 for
all four runs, very small independent bordered residuals
(`1.57e-26`, `2.52e-28`, `3.76e-24`, `2.57e-27`);
full original `IFunction` verifies the resulting physical qdot0/Y0.
Native MatNest block index mappings/RHS roundtrips are exact, MPI1/2
coefficient parity passed, and original differential data and PETSc/
GetFEM state were restored after diagnostic. No accepted Y0 was
published and **no pseudo/BDF or physical TS step** occurred.

The true `partial_t F_a(t0,z0,0)=0` is analytically exact for
this new C1 wall onset: the original decreasing envelope is flat
for `t<.03`, and `g'(0)=0`.
Independently, the actual original GetFEM wall U control at
`delta=1e-4` and `2delta` has native ratio
`3.9866220735785958`, equal to the exact smoothstep quotient;
the other row temporal variation is **exactly zero**.
No spurious fixed-step derivative was substituted.

**Host fixture regression NOTE:** the original eight synthetic host
tests had 7 successes and one **test assertion mismatch only**:
`test_rejects_corrupt_MPI2_solution` perturbed the **scalar
norm** of MPI2 Y0, so the independent comparator correctly rejected
it at `Original native MPI first-differentiated mismatch in Y0`
*before* reaching the intended `vector MPI parity` exception.
This is NOT a native solver/MatNest failure.
A new source-only host fix at `f4c2303` separates
(1) scalar Y0 norm discrepancy from (2) a deliberate MPI1/MPI2
**coefficient permutation at matching norm**, ensuring the full-vector
gate is genuinely exercised. Now nine synthetic host controls
are expected; **host rerun pending**. No native output or numerical
acceptance criterion was modified.

**Qualified narrow conclusion:** the explicitly new physical
rest-to-spinup curved-P2 SWE has an exact t0 FE algebraic candidate
and a compatible first-differentiated *physical* qdot0/Y0 candidate
to near-machine precision for both meshes, MPI1/2.
The original already-rotating q0 remains proven algebraically
incompatible and was NOT altered.
**NOT qualified:** general high-index DAE initializer, generic
differentiation index, first **physical CN** step, sustained rotating
time evolution, port-Hamiltonian energy / boundary-work / viscous-
dissipation balance. The existing production nonlinear-affine
initializer still correctly rejects high-index layouts.
#28 OPEN, nontrivial rotating #6 physical time stepping BLOCKED
pending safe accepted-state publication and full physical CN audit.

**Host regression rerun (no native repeat required):**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_spinup_first_derivative.py
# Optional pure host recheck of UNMODIFIED archived native evidence:
bash scripts/run-v130-issue28-spinup-first-derivative.sh \
  --compare-only sandbox/outputs/v130-issue28-spinup-first-differentiated-54bafa2c-y1es00lx
```

## 21 November 2026 — #28/#6 first actual physical CN step: sandbox-only native gate STAGED

All previous source and native qualification gates are complete:
nine host negative/regression tests PASS and the complete archived
native 4/4 MPI1/2 compatible rest-spinup first-differentiated
physical RHS gate PASS at source `54bafa2` (archive
`sandbox/outputs/v130-issue28-spinup-first-differentiated-54bafa2c-y1es00lx/`).
The first-differentiated candidate for the **distinct compatible
physical spinup** is `qdot0≈0,Y0≈0`, with native residuals
~1e-24 or better; the unrelated original already rotating
`p0` remains proven algebraically incompatible and untouched.

**NEXT GATE — SOURCE STAGED / NATIVE NOT YET RUN**:
`sandbox/probe_v130_issue28_one_cn_spinup.py`,
`sandbox/compare_v130_issue28_one_cn_spinup.py`,
`tests/test_v130_issue28_one_cn_spinup.py`,
`scripts/run-v130-issue28-one-cn-spinup.sh`.

The new sandbox diagnostic uses the real original curved P2
Saint-Venant physics (unchanged SWE weak bricks, physical
facet-tangential wall, `gamma=0,h0=5,p0=0`), original
rotation target `omega=.2` driven by the separately
declared rising C1 smoothstep over `.02`.
Only within a single Python process, a temporary
`DPHS.init_step` monkeypatch publishes the analytically
native-verified t0 effort/traction coefficients
`e_h=.5,Y_n=-.5,e_p=U=U_n=Y=0`.
It evaluates the **true original IFunction** before
publication, demands zero original initial algebraic
and dynamic residual, keeps original q0 bitwise,
and writes exactly the same candidate to both GetFEM
and replicated PETSc `GlobalVec`. It never edits
production `scrimp/dphs.py`, its initializers or
the original incompatible experiment. This
patch is a **sandbox test fixture**, not a
new general high-index initializer.

The original SCRIMP PETSc TS then executes
`ts_type=cn,ts_adapt_type=none`
from `t0=0` to `tf=.002`
with exactly **one accepted step** and no
SNES failure/rejection. The probe fails closed
on wrong actual TS convergence reason, number of
steps, final physical time or missing exact
endpoint records.

**Physical audit before any energy assertion:**

- Re-evaluate the unchanged native GetFEM
  `IFunction` at actual CN midpoint
  `zmid=(z0+z1)/2`,
  `zdot=(z1-z0)/dt`,
  `tmid=dt/2`, independently of PETSc
  KSP/SNES convergence.
- Additionally archive original full and
  algebraic endpoint residuals at `tf`.
  **Subtle:** CN enforces the time-dependent
  algebraic constraint at the midpoint,
  not automatically at the endpoint.
  A nonzero endpoint `Fa(tf)` cannot
  silently be reclassified as compatibility;
  we will evaluate its magnitude and
  time-step sensitivity.
- Reuse the original SWE physical diagnostic
  `sandbox/viscous_swe_2d_mpi.py:diagnostics`,
  which correctly measures the **h-weighted**
  normal/viscous boundary powers and **separate**
  grad/div volume dissipation. The generic
  SCRIMP port power is *not* equivalent to
  these weighted physical SWE terms.
  Archive `H(tf)-H0`,
  `W=int(P_normal+P_viscous)dt`,
  `D=int(D_grad+D_div)dt`,
  `H(tf)-H0+D-W`,
  volume balance, and ratio of defect to
  **actual** exchanged energy rather
  than enormous background H0.
- No pass/fail *convergence* threshold
  is inferred from one dt. Physical energy
  convergence, sustained rotation and
  validated high-index initializer remain
  **explicitly UNQUALIFIED** until
  independent time refinement.

Two meshes nominal h=.2/.1 × MPI1/2,
immutable `scrimp-mpi:v1.3.0-rc3`
OCI `1d6d0c90115213298047b28237a7f38087563a35`.
Host-owned output stores per-rank raw logs,
two full physical endpoint vectors,
time traces of H and power terms and
JSON with original TS reason/steps/false
energy qualifiers. Comparator checks
both native scalar parity and **full**
initial/final FE coefficient MPI parity.
The host negative-test fixtures are
synthetic, not native evidence.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_one_cn_spinup.py
bash -n scripts/run-v130-issue28-one-cn-spinup.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-one-cn-spinup.sh
```

**STOP conditions:** if the real TS/SNES fails,
or if the independent CN midpoint equations
do not match the original IFunction,
do NOT claim any physical evolution.
If TS accepts one step, record the original
endpoint constraint and weighted power defects
without declaring energy conservation based
on a single step. Next gate will use temporal
refinement and possibly an endpoint
DAE-consistency policy that preserves
original physics.
#28 OPEN; physical rotating #6 remains
NOT YET QUALIFIED.

## 21 November 2026 — First physical CN pilot: PETSc TSCN midpoint false gate corrected, nonlinear-mass caution (SOURCE ONLY)

**Original native attempt at `d46c0302ba11a5cd4179998a37b71f5ba0577061`:**
host synthetic regression **9/9 PASS**; pinned
`scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
First native case h=.2 MPI1 published the exact
compatible physical `h0=5,p0=0`,
`e_h=.5,Y_n=-.5` state with
`Fa=1.47634222e-17,Fd=2.69964047e-16`,
then PETSc reported **1 ACCEPTED TSCN step**
`t=.002`, `0` rejected, `0` nonlinear solver
failures, 2 Newton and 2 linear iterations (elapsed ~.76s).
The sandbox then **incorrectly rejected** this actual
PETSc step because it evaluated an **implicit MIDPOINT
method** residual
`F(t=.001, (z0+z1)/2, (z1-z0)/.002)`
and found `||F_mid||_2=203.27337`.
This is a **wrong verification equation**, not
native evidence of failed PETSc TS or of physical correctness.
Only h=.2 MPI1 ran; h=.2 MPI2/h=.1 MPI1/2
were not started. The original incomplete evidence
archive is
`sandbox/outputs/v130-issue28-one-cn-spinup-d46c0302-q5YBlZXK/`
(log and manifest, **no accepted report or
power traces**).

**Authoritative PETSc source analysis**
([TSCN](https://petsc.org/release/manualpages/TS/TSCN/),
[TSTHETA](https://petsc.org/release/manualpages/TS/TSTHETA/),
`src/ts/impls/implicit/theta/theta.c`):
`ts_type=cn` is PETSc **ENDPOINT
trapezoidal theta=0.5**, *not* implicit midpoint.
In the source `shift=1/(theta*dt)=2/dt`,
`stage_time=t0+dt`, and the PETSc SNES
right-hand side is
`((theta-1)/theta)*F(t0,z0,0)=-F(t0,z0,0)`.
SNES thus solves the **actual** original
discrete equation

```text
F(t0+dt, z1, 2*(z1-z0)/dt) + F(t0, z0, 0) = 0.
```

**IMPORTANT: the sign is PLUS**; the earlier
verbal formula using a minus sign was corrected
after checking the actual PETSc source.

The updated native probe checks the **full
coefficient vector** from two independent calls
to the real SCRIMP `IFunction` for
both the dynamic and algebraic rows, with
`||F_TSCN||<1e-7`. It also measures
`||F_a(t1,z1,0)||<1e-7`, records
the old midpoint residual under
`OUT_OF_SCHEME` as a *nonacceptance diagnostic*,
and separately records original
`F(t0,z0,0)` and `F(t1,z1,2Δz/dt)`.
No production integrator or SWE weak form
has been modified.

**Independent methodological caution:**
PETSc's `theta.c` explicitly comments
`This formulation assumes linear time-independent mass matrix`
inside the **endpoint branch**. SCRIMP's
viscous Saint-Venant has actual nonlinear
`M(z)=M(h)`. Therefore even if the
correct PETSc TSCN equation passes,
that only qualifies what the native
*implementation* solved; it **does
NOT** qualify expected time order,
discrete pH power/energy balance,
general nonlinear-mass DAE suitability,
or sustained rotating time evolution.
An independent proper nonlinear-mass
and energy time-refinement study will
be mandatory. The PETSc implicit
midpoint method is a DIFFERENT method
(and per PETSc not stiffly accurate
for a general DAE), so must not be
silently substituted.

**Source-only corrections now staged:**
`sandbox/probe_v130_issue28_one_cn_spinup.py`,
`sandbox/compare_v130_issue28_one_cn_spinup.py`,
`tests/test_v130_issue28_one_cn_spinup.py`.
Tests now accept an arbitrary `203.27337`
*OUT_OF_SCHEME* midpoint residual but
fail on true TSCN endpoint inconsistency,
endpoint algebraic drift, incorrect
method metadata, incorrect boundary
power book-keeping, or MPI vector mismatch.
The previous h=.2 MPI1 failed gate needs
a **fresh native run**; it cannot be
relabelled PASS from PETSc's accepted-step
printout alone.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_one_cn_spinup.py
bash -n scripts/run-v130-issue28-one-cn-spinup.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-one-cn-spinup.sh
```

**Status:** original `t0` compatible state
and first differentiated DAE remain
native-qualified. **No** one-step physical
Hamiltonian/work/dissipation report has
yet passed. Both #28 and #6 remain OPEN.

## 21 November 2026 — #28/#6 true first PETSc TSCN physical step **NATIVE 4/4 PASSED**, fixed-time temporal refinement source staged

**Complete immutable native evidence**: `sandbox/outputs/v130-issue28-one-cn-spinup-a8422a9b-2SJgTXLr/`,
source `a8422a9bef09cbd370b14eea8318f09bf2598d9e`,
`scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
Original host regression 11/11 passed, and native GetFEM/PETSc
h=.2/.1 × MPI1/2 **4/4 passed** with one TRUE accepted
`ts_type=cn` (PETSc TSCN **endpoint** trapezoidal) physical
step from `t=0` to `t=.002`; no rejections or SNES failures,
full original `IFunction` check of the actual TSCN equation
`F(t1,z1,2*(z1-z0)/dt)+F(t0,z0,0)=0`.
The separate stationary initial condition and its first
differentiation were previously qualified; this is the
**new** original-physics-compatible rest-to-rotation experiment
(`h0=5,p0=0`, curved P2 disk, original physical facet-tangent
wall, rising C1 ramp over 0.02), **NOT** a silent repair of the
original incompatible already-rotating `p0`.

| Native MPI1 at dt=.002 | h=.2 | h=.1 |
|---|---:|---:|
| `||F_TSCN||_2` | 4.8111341749e-14 | 7.5278201713e-14 |
| `||Fa(t1)||_2` | 2.4051939768e-17 | 1.5886709393e-17 |
| Out-of-scheme midpoint `||F_mid||_2` (not TSCN) | .00109999032 | .000796461984 |
| Actual `||Δ(h,p)||_2` | .04703883694 | .06471127850 |
| `H(t1)-H0` | +8.858168969450e-6 | +4.658762096543e-6 |
| Integrated original h-weighted boundary work `W` | +8.869969112874e-6 | +4.681283404140e-6 |
| Integrated grad+div viscous dissipation `D` | +1.180020635604e-8 | +2.252128401310e-8 |
| **`ΔH+D-W`** | **+6.293222880478e-14** | **−2.358338947305e-14** |
| **`(ΔH+D-W)/max(|ΔH|,|W|,|D|)`** | **+7.094977220771e-9** | **−5.037804259447e-9** |
| Volume-balance defect | +8.881784197001e-15 | −9.237055564881e-14 |

MPI2 TSCN residuals are `4.8116952862e-14` and
`7.5277495100e-14`, with algebraic residuals
`2.3471244434e-17` and `1.5854092554e-17`.
MPI2 energy errors are +6.1155871965e-14
and −1.5589783696e-14, consistent to rounding
with MPI1. Complete original physical endpoint
FE vectors (3300/10272 coefficients),
Hamiltonian and h-weighted boundary/dissipation
traces were inspected from all four native
archives. This is genuinely nonzero physical
evolution; no energy terms are synthetic.
The positive SWE height coefficients remain
positive. No production initializer or
original physics was changed.

**NARROW qualified conclusion:** a genuinely
accepted PETSc TSCN step reproduces its own
original discrete endpoint residual to ~1e-13,
keeps the full algebraic constraints to ~1e-17
and shows physical `ΔH+D-W` ~1e-14 (relative
to actual energy exchange ~1e-8) for the
*first* C1-spinup step on two P2 meshes,
MPI1/2. The out-of-scheme implicit midpoint
residual is NOT a failure criterion.
This is a very encouraging one-step
port-Hamiltonian diagnostic, but it
**DOES NOT yet prove temporal accuracy,
energy convergence or long rotating trajectories**.

**IMPORTANT scope:** PETSc `theta.c` comments
that its endpoint TSCN equation assumes a
*linear time-independent mass*, whereas actual
SWE has nonlinear `M(h)`. Passing the
native TSCN residual does NOT remove
this methodological caveat.
Neither #28 (general higher-index DAE)
nor #6 (physical long-time MPI SWE) is closed.

**NEXT: SOURCE-STAGED / NATIVE PENDING**
independent temporal refinement at *fixed physical*
`t_final=0.002` with `N=1,2,4,8`
steps (`dt=.002/N`), for
each nominal P2 mesh `h=.2/.1`, MPI1/2.
16 native runs are planned, retaining the
same original physical bricks, C1 wall,
`h0=5,p0=0` and exact sandbox-only
GetFEM/GlobalVec t0 publication.

- `sandbox/probe_v130_issue28_cn_refinement.py`
  evaluates **every** accepted endpoint original
  `F(t_k,z_k,0)+F(t_{k+1},z_{k+1},
  2Δz/Δt)` (actual TSCN), every
  `Fa(t_{k+1})` and the *out-of-scheme*
  midpoint for diagnostics only, while
  checking true PETSc step count, final
  time, MPI collective completion and
  positive h coefficients.
- It recomputes the original **h-weighted**
  physical wall powers, separate grad/div
  viscous dissipation, GetFEM mechanical H,
  and volume balance at every saved time;
  integrates each on the true accepted
  time grid and archives complete original
  physical FE vectors in `cn_refinement.npz`.
- `sandbox/compare_v130_issue28_cn_refinement.py`
  independently reconstructs the physical
  power integrals from raw saved traces,
  verifies MPI1/2 **full FE vector parity**,
  strict same-initial-state and same
  spatial layout across N, and compares
  only **original dynamic h,p** terminal
  coefficient differences:
  `||q_N-q_{2N}||_2` for N=1,2,4,
  with observed Richardson-style
  `log2(error_N/error_2N)` estimates.
  These are finite-dimensional coefficient
  diagnostics, NOT mesh-independent norms,
  and no order-2 gate is forced before
  inspecting native evidence.
- `scripts/run-v130-issue28-cn-refinement.sh`
  pins the original immutable Docker OCI
  and uses original sandbox-only
  sources read-only; `tests/test_v130_issue28_cn_refinement.py`
  contains eleven host synthetic positive/negative
  checks. **Host and native results
  are still pending.**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_cn_refinement.py
bash -n scripts/run-v130-issue28-cn-refinement.sh
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-cn-refinement.sh
```

**Remaining limitations:** the first 0.002 seconds
are only the rising start of a 0.02-second
wall spinup. Even if this time refinement
shows the expected slope and near-machine
power defect, it cannot by itself certify
production high-index initialization,
stability during full wall-rotation ramp,
or accurate/discrete-gradient pH balance
on long trajectories. Review actual
time-refinement evidence first.

## 21 November 2026 — #28/#6 NATIVE fixed-final-time PETSc TSCN refinement **16/16 PASS** (t_f=.002, N=1/2/4/8, curved P2 MPI1/2)

**COMPLETE original native evidence**: `sandbox/outputs/v130-issue28-cn-refinement-a9cc2d21-hhcttZc6/`, source `a9cc2d21ead982e5582ae77813a5f0558a1f694f`, pinned `scrimp-mpi:v1.3.0-rc3` OCI `1d6d0c90115213298047b28237a7f38087563a35`. Independently inspected the original uploaded 3.1 MB ZIP: all 16 native JSON and NPZ files, all 16 time/physical coefficient trajectories and original weighted power traces, original manifest and strict `cn_refinement_comparison.json` are present. All 11 host regression/negative tests PASS. Original comparator returns `NATIVE_TSCN_FIXED_FINAL_TIME_REFINEMENT_DIAGNOSTIC_ONLY`.

Experiment is the **NEW** explicitly compatible C1 rest-to-rotation physical SWE (`h0=5,p0=0`, original curved P2 weak physics, original physical facet-tangent wall, target `omega=.2`, rising `g(t)` over .02 seconds, current observation `t_f=.002`). The original incompatible *already rotating* initial p0 and production initializers are **not** modified. Temporary sandbox exact-effort/traction publication from earlier native qualified t0 is used for all runs. Each real PETSc TSCN (endpoint trapezoidal, `theta=.5`) step is independently evaluated using the TRUE original GetFEM IFunction equation `F(t_{k+1},z_{k+1},2(z_{k+1}-z_k)/dt)+F(t_k,z_k,0)=0`; the implicit-midpoint residual remains explicitly out-of-scheme.

**Native MPI1 raw data (MPI2 also PASS):**

| mesh | N | dt | max original TSCN residual | max endpoint Fa | ΔH | original h-weighted W | viscous D | ΔH+D−W | defect / actual exchange |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| .2 | 1 | .002 | 4.81e-14 | 2.41e-17 | 8.85816897e-6 | 8.86996911e-6 | 1.18002064e-8 | +6.29e-14 | +7.09e-9 |
| .2 | 2 | .001 | 9.99e-14 | 2.64e-17 | 8.85815305e-6 | 8.86790635e-6 | 9.75341625e-9 | +1.12e-13 | +1.26e-8 |
| .2 | 4 | .0005 | 1.95e-13 | 3.60e-17 | 8.85814976e-6 | 8.86781724e-6 | 9.66759114e-9 | +1.09e-13 | +1.23e-8 |
| .2 | 8 | .00025 | 3.96e-13 | 4.75e-17 | 8.85814897e-6 | 8.86782170e-6 | 9.67282062e-9 | +9.71e-14 | +1.10e-8 |
| .1 | 1 | .002 | 7.53e-14 | 1.59e-17 | 4.65876210e-6 | 4.68128340e-6 | 2.25212840e-8 | −2.36e-14 | −5.04e-9 |
| .1 | 2 | .001 | 5.26e-14 | 1.43e-17 | 4.65871765e-6 | 4.67735698e-6 | 1.86392535e-8 | −7.69e-14 | −1.64e-8 |
| .1 | 4 | .0005 | 1.07e-13 | 1.89e-17 | 4.65870844e-6 | 4.67719184e-6 | 1.84833457e-8 | −5.62e-14 | −1.20e-8 |
| .1 | 8 | .00025 | 2.11e-13 | 2.52e-17 | 4.65870624e-6 | 4.67720172e-6 | 1.84954129e-8 | −5.92e-14 | −1.27e-8 |

Max endpoint algebraic residual over all native MPI and temporal cases **<4.93e-17**, max TSCN residual **<3.97e-13**, max absolute mass balance defect **<1.60e-13**, and absolute relative power balance defect **<1.65e-8** of the *actually exchanged* physical energy, not of background H0. Independent raw trace reintegration reproduces all weighted wall/viscosity and Hamiltonian balances. MPI1/2 full FE vector parity and identical initial prescribed q0 PASS.

**Temporal original dynamic h,p FE coefficient differences** at fixed final physical time, same mesh:

| mesh | `||q1−q2||` | `||q2−q4||` | `||q4−q8||` | `p(1,2,4)` | `p(2,4,8)` |
|---|---:|---:|---:|---:|---:|
| .2 | 1.11668109484e-4 | 2.79155203799e-5 | 6.97878678223e-6 | **2.0000778803** | **2.0000192900** |
| .1 | 2.35596039991e-4 | 5.88720774443e-5 | 1.47163414212e-5 | **2.0006598472** | **2.0001644850** |

This is compelling evidence of **second-order temporal behavior for this specific early C1 spinup trajectory**, in the original physical dynamic FE coefficients on both tested spatial meshes. These are coefficient norms (not mesh-independent continuum norms) and numerical order estimates (not a general theorem). The power defect is near numerical roundoff, not systematically decreasing with dt, so **no energy-defect order** is inferred.

**Important still-open caveat:** PETSc endpoint TSCN implementation documents its assumption of linear time-independent mass; physical SWE `M(h)` is nonlinear. Empirical early-time order ≈2 does not remove the possible loss of general order/structure over other data, longer times, or generic DAE index. No general nonlinear-mass TSCN theorem/production initializer, nor full spinup/high-index long-time validation is claimed.

**Next physically meaningful gate:** keep original C1 ramp, original nonlinear weak bricks, h0=5,p0=0, but integrate THROUGH the full ramp `t_f=.02` (20 ms), with fixed dt .002/.001/.0005 = N=10/20/40, nominal meshes .2/.1 and MPI1/2. Check true original TSCN and Fa **at every step**, real TS acceptance, full FE coefficient parity, mass and h-weighted physical work+dissipation, and temporal endpoint differences/observed orders. Require positive h throughout; compare energy balance without arbitrary absolute tolerance masquerading as a proof. At t=.02 the rising C1 wall reaches full prescribed rotation `omega=.2` for the first time; the existing decreasing envelope starts only at t=.03, so the physical control remains comparable. **This longer gate is not yet native qualified.**

## 21 November 2026 — #28/#6 FULL C1 SPINUP to t_f=.02: separate 12-run native gate STAGED (source only)

The prior **16/16 native** fixed-time early-spinup TSCN
campaign at `t_f=.002` (source `a9cc2d21`) is fully
qualified as a *specific empirical* second-order
dynamic h,p coefficient convergence study, with
max original TSCN residual <3.97e-13,
endpoint Fa<4.93e-17,
physical h-weighted H+D-W defect absolute <1.12e-13,
and full MPI1/2 FE-vector parity. Original raw
archive `sandbox/outputs/v130-issue28-cn-refinement-a9cc2d21-hhcttZc6/`.

**New scientifically distinct gate (no native result yet):**
advance through the **entire rising C1 physical
wall rotation** by holding `t_final=.02`
instead of .002 and taking `N=10,20,40`
(`dt=.002,.001,.0005`), on each original
curved P2 nominal h=.2/.1 mesh with MPI1/2.
Total **12 planned native configurations**,
pinned original `scrimp-mpi:v1.3.0-rc3`
OCI `1d6d0c90115213298047b28237a7f38087563a35`.
This preserves original physical SWE viscous,
gyroscopic and boundary forms, `h0=5,p0=0`,
original `start_from_rest_duration=.02`
and the earlier native-qualified exact
sandbox-only initial GetFEM/GlobalVec publication.
It never fixes the incompatible already-rotating q0.

At `t=.02`, `g(t)=3s²-2s³`, clipped
`s=t/.02`, is exactly `g=1`;
`g'(t)=0` there. The original descending
envelope remains 1 until t=.03, so the
physical facet-tangent wall reaches full
target `omega=.2` exactly at the
end of the new simulation. New probe
independently re-evaluates original GetFEM
`U` port forcing with original frozen
stationary q0 at t=.02 and t=.02001,
requires `||F_U||>0` and matching
terminal plateau values to relative
<1e-10. This catches an accidentally
absent/mis-timed spinup control.

**Source files staged:**
- `sandbox/probe_v130_issue28_full_spinup.py`:
  original actual per-step TSCN endpoint
  IFunction, Fa(t_k), real accepted TS
  step counts/no SNES failure,
  strictly positive h coefficients,
  complete FE vector trajectory, original
  h-weighted wall normal/tangent work,
  grad/div dissipation, Hamiltonian
  and mass trace, and independent
  original physical wall plateau check.
- `sandbox/compare_v130_issue28_full_spinup.py`:
  strict 12-case provenance, original
  h,p FE intervals, fixed physical
  initial q0, same time horizon,
  weighted work/dissipation reintegration,
  full FE MPI1/2 parity and terminal
  dynamic h,p coefficient time differences
  `||q10-q20||`, `||q20-q40||`
  and **one empirical** observed order
  `log2(E10/E20)`.
- `tests/test_v130_issue28_full_spinup.py`:
  **12** synthetic host regression/negative
  cases; synthetic log numbers do
  not qualify native PETSc.
- `scripts/run-v130-issue28-full-spinup.sh`:
  12 native runs, pinned immutable
  original OCI, source code mounted
  read-only, host-owned evidence/logs,
  fail-closed native comparison.
  No production `scrimp/dphs.py` or
  original SWE weak physics changed.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_issue28_full_spinup.py
bash -n scripts/run-v130-issue28-full-spinup.sh

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue28-full-spinup.sh
```

**Scope:** This next gate tests one **completed
boundary acceleration ramp** but does NOT
cover the original descending envelope
after t=.03, a long-time rotating
Saint-Venant trajectory, general high-index
initialization, proof of temporal order
or discrete-gradient pH structure.
PETSc endpoint TSCN assumes linear
time-independent mass, unlike M(h);
continued empirical success does
not remove this documented caveat.
No native or host PASS has been
reported for the new full-spinup
source until the user runs it.

## 21 November 2026 — Full C1 spinup staging syntax repair; native gate still pending

User's first `scripts/run-v130-issue28-full-spinup.sh` run at source `e2ec401`
**did not execute PETSc at all**: both the host unittest importer
(`sandbox/compare_v130_issue28_full_spinup.py:346`) and the first native
probe (`sandbox/probe_v130_issue28_full_spinup.py:324`) failed Python
parse with an unterminated literal/f-string introduced while cloning the
0.002-second refinement scripts. Docker launched the first case but the
Python interpreter stopped before importing GetFEM or constructing a
simulation; all 12 native configurations remain **UNTESTED**. The
partial output directory
`sandbox/outputs/v130-issue28-full-spinup-e2ec401f-hCnrFHWt/`
is a *failed syntax-gate artifact*, not physical evidence.

The two literals have been repaired without modifying original SWE
physics, boundary control, PETSc TSCN solver, time-grid design, numeric
thresholds, or comparator mathematical tests. To prevent recurrence,
`.github/workflows/python-source-checks.yml` now compiles the newer
#28 native probes and runs the synthetic host comparison tests. A
separate fast workflow,
`.github/workflows/v130-issue28-full-spinup-host.yml`, executes
full-spinup Python compilation, `bash -n` and the 12 host tests
without requiring GetFEM/PETSc. These are source-only safety checks,
**not native qualification**.

Before any Docker rerun, perform the explicit local fail-fast checks:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python -m py_compile \
    sandbox/probe_v130_issue28_full_spinup.py \
    sandbox/compare_v130_issue28_full_spinup.py \
    tests/test_v130_issue28_full_spinup.py
python tests/test_v130_issue28_full_spinup.py
bash -n scripts/run-v130-issue28-full-spinup.sh
```

Only if all pass, run the same immutable
`scrimp-mpi:v1.3.0-rc3` 12-case full-spinup campaign.
The earlier **16/16 native PASSED** refinement at `t_f=.002` and
its empirical order-two result remain valid and unaffected.
#28 and #6 remain OPEN.

## 21 November 2026 — Full-ramp native 12/12 PASS; accepted-state checkpoint becomes v1.3.0 P0 (#29)

The **separate new compatible** rest-to-rotation
Saint-Venant curved-P2 experiment has now successfully
advanced through the complete physical **C1
wall acceleration to t=.02**. User archive
`sandbox/outputs/v130-issue28-full-spinup-c0b8f3d2-NSdC2dDW/`
at source `c0b8f3d2e4fe7f1ea8466020f1bdf163afa764da`
and pinned `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`
was inspected: **12 native JSON/NPZ cases**, h=.2/.1,
N=10/20/40 TSCN steps, MPI1/2, accepted physical
times, full original FE coefficients and original
h-weighted boundary work/viscous dissipation/Hamiltonian
traces. Strict comparator
`NATIVE_TSCN_FULL_C1_RAMP_REFINEMENT_DIAGNOSTIC_ONLY`
PASS, complete wall target control plateau verified.
Maximum `||F_TSCN||=1.598370025e-10`,
maximum endpoint `||F_a||=3.242346397e-14`,
maximum absolute volume defect
`9.414691249e-14`, full MPI1/2
FE parity PASS. Observed original dynamic
h,p coefficient order is **1.84153** (h=.2)
and **1.96182** (h=.1).

**Important:** h-weighted `ΔH+D-W` relative
to actual exchanged energy at t=.02 with N=40
is **+2.79935e-5** (h=.2) and
**+6.94182e-5** (h=.1), not the earlier
t=.002 ~1e-8 roundoff floor. The
defect reduces only modestly over this
temporal grid. No claim of general
nonlinear-mass TSCN pH power preservation,
stationary interior rotation, or general
index-two/high-index DAE initializer
follows from this native experiment.
The original already-rotating t0 differential
data are still infeasible on the tested meshes.

**New P0 #29:** preserve a previously
accepted complete DAE state (including all
efforts and boundary variables), its physical
time and verifiable actual FE mesh/DoF
ordering, and **restart without a new
pseudo-time initial consistency solve**.
The cold restarted state must satisfy the
original algebraic rows of `IFunction`
at `t=t_*` and unchanged physical control.
This is a safe workaround for repeatedly
preparing LHMNC initial data; it is NOT
a method to repair an inconsistent prescribed
initial differential state and does not
replace issue #28's independent nonlinear/
high-index feasibility and initializer work.
See [release gate 2C](V130_RELEASE_PLAN.md)
and [checkpoint issue #29](https://github.com/g-haine/scrimp/issues/29).

## 23 November 2026 — #29 original Saint-Venant native COLD checkpoint/restart pilot **16/16 PASSED**, production API staged separately

**Qualified native sandbox evidence**:
`sandbox/outputs/v130-issue29-cold-restart-d3c26930-nRhUY0qA/`,
source `d3c269304cf0fe7fe6b46d45844f50b1120a84c2`,
pinned immutable Docker `scrimp-mpi:v1.3.0-rc3`
OCI `1d6d0c90115213298047b28237a7f38087563a35`.
The user supplied actual final ZIP
`v130-issue29-cold-restart-d3c26930-nRhUY0qA-restart.zip`,
which was independently opened: **16 original
real PETSc/GetFEM phase reports, 16 full state+power
NPZ trajectories, 4 accepted full FE checkpoints,
one immutable source manifest and complete
`checkpoint_comparison.json`**. ZIP CRC intact,
comparison status
`NATIVE_COLD_CHECKPOINT_SPLIT_VS_CONTINUOUS_PILOT_PASSED`.
All 17/17 checkpoint-state + 13/13 continuous/restart
synthetic host regression controls PASS.

Experiment: original curved P2 dissipative
Saint-Venant with nonlinear mass, physical
facet-tangent wall C1 rest-to-rotation
from `h0=5,p0=0` (the previously
qualified compatible sandbox initial state),
`dt=.002`: real PETSc `0→.02` seed
(10 accepted), `0→.024` reference
(12 accepted), and `.02→.024`
fresh PETSc TSCN from the full saved state
(2 accepted). Both physical meshes
h=.2/.1 and **all 1→1,1→2,2→1,2→2**
source-to-target MPI ranks independently
tested (8 true cold restarts). No
PETSc TS/KSP/MUMPS history serialized;
no pseudo/BDF init/projection on restart.

**Bitwise checkpoint creation:** all 4
accepted-state full vectors match their
native seed's original `states[-1]`
exactly and SHA256 checksum passes:
3300 DOFs for nominal h=.2,
10272 DOFs for h=.1, at absolute
physical `t_*=0.02`. Each contains
full original dynamic and algebraic
GetFEM coefficients and immutable mesh/
FE/order/weak-form/control provenance.

**Observed genuine continuous-vs-split
FULL FE vector L2 differences at `t=.024`:**

| Nominal mesh | MPI 1→1 | 1→2 | 2→1 | 2→2 |
|---|---:|---:|---:|---:|
| h=.2 | 1.79515e-13 | 2.48915e-13 | 2.05195e-13 | 2.10372e-13 |
| h=.1 | 1.67651e-13 | 1.52415e-13 | 1.69267e-13 | 1.71498e-13 |

Hamiltonian **restart-time jump** max
`7.55e-15`. Independent original
h-weighted boundary+viscous+Hamiltonian
cumulative physical `ΔH+D−W`
difference between uninterrupted and
two-segment cold restart max
`3.02007e-14`; all algebraic rows
at resumed t* were checked to about
`1e-15`, and the two resumed
original endpoint TSCN steps and
volume balances were accepted
without failed SNES or rejected steps.
This establishes **native scientific
validity of the narrow sandbox
complete-accepted-state cold restart
for this P2 viscous SWE experiment**.

**Not established by this test:**
general arbitrary high-index DAE
initialization, long-term nonlinear-mass
TSCN power/accuracy, a stationary
interior rotating equilibrium, exact
BDF/IMEX history restart,
repartitioned/disparate FE numbering,
or the PRODUCTION `DPHS` API.
Do not close #29 on the sandbox pilot
alone or silently qualify a BDF restart.

**NEXT SOURCE ONLY: generic DPHS public
checkpoint API staged** (no native
results for this implementation yet):
- `scrimp/utils/dae_checkpoint.py`
  promotes the native-compatible
  versioned atomic full FE + model/
  mesh/FE DoF/weak-form fingerprint
  checkpoint reader/writer into the
  package.
- `scrimp/dphs.py` now exposes
  `save_checkpoint(path, extra=...)`
  after a fully accepted final TS
  and `load_checkpoint(path)`
  on a FRESH assembled model.
  The latter stages a complete state
  without mutating the physical initial
  differential variables until
  after the real PETSc/GetFEM
  assembly; before TS creation the
  original algebraic rows
  `F_a(t_*,z_*,0)` are rechecked.
  The solve then uses absolute t*
  and skips pseudo/BDF/algebraic
  consistency initialization.
  It is COLD state restart only,
  with no integrator history.
- Separate
  `sandbox/probe_v130_checkpoint_public_api.py`,
  `sandbox/compare_v130_checkpoint_public_api.py`,
  `sandbox/validate_v130_checkpoint_public_api_case.py`,
  `scripts/run-v130-checkpoint-public-api.sh`
  repeat true 16-case original
  physical SWE continuous/split and
  MPI cross-rank campaign using
  the actual **production DPHS
  save/load methods**. Native pending.
- `tests/test_v130_checkpoint_public_api.py`
  has host-only API-method AST,
  old/new file-schema interoperability
  and fail-closed physical-control
  regression checks. CI
  `v130-issue29-checkpoint-host.yml`
  now compiles the public API and
  runs the new suite.
- **Scientific gate still OPEN:**
  test new generic public API on a
  **second independent linear DAE
  model** and qualify its
  full native continuous/restart
  energy and algebraic constraints.

Preflight and **NEW API** native study:
```bash
cd ~/Simulations/scrimp
git pull --ff-only

python -m py_compile \
  scrimp/dphs.py scrimp/utils/dae_checkpoint.py \
  sandbox/probe_v130_checkpoint_public_api.py \
  sandbox/compare_v130_checkpoint_public_api.py \
  sandbox/validate_v130_checkpoint_public_api_case.py \
  tests/test_v130_checkpoint_public_api.py
python tests/test_v130_checkpoint_state.py
python tests/test_v130_checkpoint_restart.py
python tests/test_v130_checkpoint_public_api.py
bash -n scripts/run-v130-checkpoint-public-api.sh

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-checkpoint-public-api.sh
```

**Important:** public API is currently
source-only and may need further
corrections from first real
GETFEM/PETSc runs. The earlier
16/16 **sandbox** evidence is not
a substitute for native validation
of the new production entry points.
