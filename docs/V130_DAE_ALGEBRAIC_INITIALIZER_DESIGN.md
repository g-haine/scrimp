# SCRIMP v1.3.0 — DAE consistent initialization: mathematical method, PETSc implementation, and integration contract

**Document status (20 November 2026):** the **diagnostic**, opt-in-by-script,
MPI-distributed PETSc algebraic initialization method has passed its
native scientific acceptance tests for the specified **linear 2D acoustic
wave**. This document describes what is actually implemented, which
properties are measured, and what remains to turn it into a supported
SCRIMP option. **This is not yet a production feature.**

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
