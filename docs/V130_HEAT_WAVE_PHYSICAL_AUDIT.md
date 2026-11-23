# SCRIMP v1.3.0 — Heat–Wave 2D: physically meaningful example

**23 November 2026 — physical MODEL REVISION COMMITTED IN SOURCE;
native GetFEM/PETSc and CN power-balance qualification PENDING.**

The author confirmed the correct scientific reference:

G. Haine, D. Matignon and F. Monteghetti,
*Long-Time Behavior of a Coupled Heat-Wave System Using a
Structure-Preserving Finite Element Method*,
Mathematical Reports **24**(74), 1–2 (2022), 187–215.
https://imar.ro/journals/Mathematical_Reports/Pdfs/2022/1-2/11.pdf

The **original paper's numerical experiments are different** from
the SCRIMP `examples/heat_wave.py` example. No modifications,
even cosmetic, are allowed under `examples/published/`.

## Mathematical target (paper eqs. 2.1–2.2, 3.15)

With unit physical coefficients, heat inside `Omega_1`,
wave inside `Omega_2`, exterior wave displacement fixed:

```text
d_t T = Delta T                    in Omega_heat
d_tt w = Delta w                   in Omega_wave
T = d_t w = p                      at Gamma_interface
grad(T).n_heat = -grad(w).n_wave   at Gamma_interface
w = 0                              at Gamma_wave_exterior
```

In co-energy variables:
`e_Q=-grad T`, `p=d_t w`, `q=grad w`,
and the outward physical interface normals have opposite
directions. The energy is

```text
H = 1/2 int_heat T^2 + 1/2 int_wave (p^2+q.q)
dH/dt = - int_heat e_Q.e_Q = -D,   D>=0.
```

The outer Dirichlet `w=0` corresponds to fixed *velocity*
`p=0`, if this is compatible with its initial displacement.
It is **not** the Neumann condition `q.n=0`.

## What was incorrect or insufficient in the old example

The preceding `heat_wave.py` used mixed co-energy
variables but unsuitable/opaque combinations of FEMs
(`T:DG1`, `e_Q:CG2`, `p:CG2`, `q:DG1`,
interface/boundary `DG1`), an adaptive BDF2 time scheme,
and a purported Dirichlet outer boundary whose actual
wave boundary equations prescribed **zero traction**,
which is Neumann.

In addition, with the heat div-div weak form
`M_Q e_Q = int T div(test_e_Q) - int U_T test_e_Q.n`,
a **minus** sign is required for the temperature
Dirichlet input on the heat boundary. The former
`+U_T test_e_Q.Normal` instead prescribed
`T=-U_T`. The former `U_w=-Y_T` simultaneously
compensated this mistake in the algebraic energy
interconnection, so an apparently correct *zero
algebraic interface-power sum* did NOT by itself
establish the intended physical transmission law.
Both signs must be corrected together.

The previous native PETSc BDF2 Heat–Wave integration
already completed at t=.004 (39 accepted, 4 rejected,
0 SNES failures, 26,198 DoFs); it was the **old pilot's
false four-fixed-step acceptance test** that failed,
not MUMPS or the actual physical solver. That
historical evidence does not qualify the new model.

## New implementation in examples/heat_wave.py

The previously qualified SCRIMP/PETSc/GF **solvers and
PDE backend are NOT modified**. Only this **non-published
example** is revised.

### FEM choice and mathematical honesty

The paper's mixed FE spaces are (4.3)–(4.5):
```text
T: DG_k; e_Q: RT_k; p: DG_k; q: RT_k;
interface and outer boundary variables: CG_(k+1).
```

For its time-decay experiments the paper uses `k=0`.
The current SCRIMP example selects
```text
k=0
T      DG0
e_Q    RT0  (GetFEM FEM_RTK(2,0), intrinsic H(div))
p      CG1  (because the existing wave uses grad–grad weak form)
q      RT0  (intrinsic H(div), with proper normal traces)
U_T/Y_T, U_w/Y_w, U_bnd/Y_bnd: CG1 restricted to boundaries
```

**This is a mixed PFEM that remains different from the
paper's DG0/RT0 wave formulation.** To implement the
latter exactly, the wave weak forms and port causalities
must also be rebuilt to the paper's div-div system
(3.18–3.19) and its interface DAE (4.7–4.9).
An unverified partial replacement of `p:CG1` with
`p:DG0` while leaving `Grad(p)` or
`q.Grad(Test_p)` would remove the spatial wave
operator. Do **not** make that invalid shortcut.

Native `S.FEM("Heat flux",0,"RT")` and
`S.FEM("q",0,"RT")` require the current GetFEM build
to support `FEM_RTK(2,0)`. The example explicitly
checks that neither was silently replaced by CG.
GetFEM reference:
https://getfem.readthedocs.io/en/latest/userdoc/appendixA.html

### Boundary, coupling and model power

The heat div-div and wave grad-grad weak forms remain,
but their **physical signs** are now explicit:

- `U_T=Y_w=p`: temperature on the heat interface
  equals wave velocity.
- `Y_T=e_Q.Normal` uses the common GetFEM normal.
  To obtain physical outward heat flux,
  `orientation=+1` for heat in inner disk,
  `orientation=-1` for heat in outer annulus.
- `U_w=orientation*Y_T` imposes physical wave traction
  matching outward heat flux. The heat boundary
  brick is `-orientation*U_T*Test_e_Q.Normal`.
- The true *signed* interface power is
  `P_int=-orientation*U_T*Y_T+U_w*Y_w`, expected
  to vanish. An unweighted `U_T*Y_T+U_w*Y_w`
  sum is **not** the physical power balance here.
- For an outer **wave** boundary, a Lagrange-multiplier
  port enforces `p=U_bnd=0` (instead of setting
  zero traction); with heat on the exterior,
  `T=U_bnd=0` is imposed using the heat div-div
  temperature input. This follows the intended
  Dirichlet boundary conditions.

### Fixed Crank–Nicolson and DAE initialization

The example explicitly sets
`ts_type="cn"`, `ts_adapt_type="none"`,
`init_step_ts_type="algebraic"` and preserves
`dt=.001, dt_save=.05, t_f=2.` by default.
Users can still supply different `time_step`,
`save_step`, `final_time` and `factor_solver`
for bounded tests; if unspecified, LU=MUMPS as
in the qualified SCRIMP core. The initial Gaussian
profiles remain as before pending physical diagnostics
of the **accepted, algebraically consistent** initial
state; this does not claim that the raw prescribed
Gaussian initially satisfies every constraint.

No implicit adaptive BDF2 property is asserted.

## First native mixed-RT integration result: core initializer bug (OPEN)

Actual qualified Docker `scrimp-mpi:v1.3.0-rc3`, source
`8286f18dc87ba97b4bc690e5dbd85316f84e9f14`,
evidence ZIP `v130-main-2d-lu-8286f18d-zY4Twrpi.zip`
(11 entries, integrity PASS). The run was stopped at
`DPHS.set_initial_value("q","[0.,0.]")` **before**
mass/stiffness assembly, consistent initialization,
PETSc solve, Hamiltonian calculation or output validation.
The archive contains the source log and initial mesh,
but **no accepted numerical trajectory or physical
energy results**. Hence no new physical conclusion.

All intrinsic mixed FEMs instantiated successfully:
`T DG0`, heat flux `RT0`, wave velocity `CG1`,
wave stress `RT0`, interface and boundary `CG1`.
The wave-region RT0 restriction produced **3111
coefficients** (whole MeshFem `nbdof=4794`,
physical vector `qdim=2`). The legacy
`DPHS.set_initial_value` incorrectly forced
`reshape((qdim, int(3111/qdim)))` and raised
`ValueError: cannot reshape array of size 3111 into
shape (2,1555)`. RT/BDM intrinsic FE coefficients
are face flux moments, not component-wise nodal values.
A valid odd-length 3111-vector is completely normal.

SCRIMP core is now corrected to pass the exact
original **flat GetFEM coefficient vector** after
`basic_dof_on_region` filtering, with validation of
source evaluator coefficient count and region DOF
indices. No zero-padding, truncation, reordering,
FE element substitution, or weakening of native
scientific convergence checks. The change covers
any intrinsic-vector basis with a non-componentwise
DOF count, not merely this example.

Standalone regression
`tests/test_v130_rt_initial_values.py`
executes the **actual production method** using
controlled GetFEM-shaped fake arrays: RT0 odd
3111 DOFs, zero RT state, legacy vector CG,
scalar and unfiltered behavior, bad indices,
and mismatched evaluator size. Its source CI
is a separate guard: native GetFEM/PETSc must
still validate the actual RT0 initialization.

## Second native attempt: false nonlinear classification (OPEN)

Actual source `f97fcd1f866abf9e24591b068dc7a1651d271e69`,
original host archive
`v130-main-2d-lu-f97fcd1f-RpX0FgJi.zip` (11 ZIP
members, CRC PASS). The corrected general `DPHS.set_initial_value()`
successfully initializes **all 3111 RT0 wave coefficients**
without losing an entry. Original GetFEM linear mass/stiffness
assembly completes: model **7689** full DAE DoFs.
The failure occurs **before PETSc's physical CN integration**:

```text
ValueError: Algebraic initializer supports only linear mass and stiffness
```

Source diagnosis: `DPHS.set_control("Interface Heat","Y_w")`
and `DPHS.set_control("Interface Wave","Y_T")` created feedback
bricks with `linear=False` (historically essential for
nonlinear state-dependent feedback). Yet `U_T=Y_w` and
`U_w=Y_T` are **exact linear coupling relations**: the
weak terms `Y_w*Test_U_T` and `Y_T*Test_U_w` are
genuine GetFEM **bilinear** terms with constant tangent.
The global `linear_stiffness=False` flag was therefore
incorrect for this specific coupled model and correctly
caused the conservative *linear-only* algebraic initializer
to reject it.

**Source correction:** `DPHS.set_control()` now recognizes
only a **bare signed registered GetFEM model variable**
(`Y_w`, `Y_T`, `-Y_T`, etc.) as provably linear.
Such feedback uses `Brick(linear=True)`, so GetFEM
`add_linear_term()` includes the exact cross-variable
coupling in the constant stiffness matrix and in the
native `K_aa` projection. Arbitrary macros, sums,
products, coefficients or unknown expressions are **not**
guessed linear: their GetFEM nonlinear weak form and
Newton derivative remain unchanged. Prescribed open-loop
sources still use the existing independent RHS route.

Dedicated source-only tests execute actual production
`set_control()`, `add_brick()`, `set_linear_flags()`,
and classify source/linear/nonlinear GetFEM assembly
routes with negative controls. **No native PASS yet.**
A subsequent run must check the actual fixed-state
projection, its true algebraic `IFunction` residual,
the four fixed CN steps and the separate physical
heat/wave power balances. In particular, a structurally
singular coupled `K_aa` can still produce an explicit
initialization rejection and must not be masked by a
pseudo-time fallback.

### Intrinsic vector FE initial coefficients: generality and limits

The `DPHS.set_initial_value()` fix is **not RT0-specific**.
It flattens the exact `MeshFem.eval()` coefficient array
and retains `basic_dof_on_region()` indices in the
original GetFEM order; it does not reshape by `qdim`.
That structural rule also applies to intrinsic
`RT_k`, `BDM_k` and Nédélec edge basis functions.
However **coefficient-layout compatibility is not
complete family qualification**: element availability,
interpolation, Piola transformations, partial-region
index matching, normal/tangential traces, and actual
Maxwell/PDE weak-form solves require native tests.

The SCRIMP `FEM` class currently exposes native
`RT` and `BDM` string selectors; GetFEM documents
`FEM_RTK(P,K)` for simplices (K>=0),
`FEM_BDMK(P,K)` for simplices (K>=1), and
**`FEM_NEDELEC(P)` only for lowest-order edge
elements**. There is no native `Nedelec1`/`Nedelec2`
selector in the current SCRIMP `FEM` class.
GetFEM may accept an explicit descriptor through
SCRIMP's expert `order<0` escape hatch, but no
higher-order Nédélec first/second-family implementation
or support is thereby proven. Missing intrinsic
elements can currently fall back to CG in `FEM.set_fem()`;
production tests must require **strict element
identity** and never silently treat a fallback as PASS.

Therefore: `RT0` coefficient indexing observed in
the real example; higher `RT_k`, `BDM_k` and
Nédélec 1st/2nd families still **UNQUALIFIED**.
This is a pending separate FEM-family validation gate,
not a reason to regress to CG/DG or reshape/truncate
coefficient vectors.

An optional native, non-PDE **element-family coefficient
pilot** has now been prepared. It creates genuine
triangular GetFEM `MeshFem` objects, restricts their
DOFs to a subset of elements, constructs the same
`Model.add_filtered_fem_variable` that production SCRIMP
uses, calls the **actual**
`DPHS.set_initial_value()` on nonzero vector data,
then checks the resulting GetFEM model coefficients
against the independently evaluated ordered restriction.
Candidates: `FEM_RTK(2,0)`, `RTK(2,1)`,
`RTK(2,2)`, `FEM_BDMK(2,1)`,
`FEM_BDMK(2,2)`, and
`FEM_NEDELEC(2)`. A descriptor unavailable in this
particular native image is reported as **UNAVAILABLE**,
never passed. A descriptor that instantiates but fails
coefficient initialization makes the pilot FAIL.
Nédélec higher-order and second-kind elements are
explicitly **not implemented or qualified** by this
probe.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-intrinsic-fem-initialization.sh
```

This independent check produces host-owned JSON and a
small ZIP. It neither runs a physical PDE nor establishes
Maxwell/operator or MPI support. **Native pilot PENDING**.

## Genuine physical qualification — native PENDING

The revised `scripts/run-v130-example-2d-lu.sh`
by default runs only `heat_wave` at MPI1,
MUMPS and SuperLU_DIST (the unchanged
`shallow_water` stage is selected separately).
Each model uses its original geometry and
weak-form/constitutive parameters with explicit
short horizon t=0..0.004, fixed CN dt=.001.
For heat-wave, exactly 4 accepted steps are now
a correct **CN no-adaptation** criterion.
ParaView exporter is not proven MPI-safe: do
**not** run this top-level plotting example on
multiple ranks without a separate export audit.

The independent `sandbox/measure_v130_heat_wave_physics.py`
records the actual original GetFEM quantities:

- separate heat/wave energies and exact sum vs H;
- nonnegative heat flux dissipation;
- each raw interface port, true oriented net power,
  and outer boundary power;
- candidate physical temperature and flux trace errors,
  with unsupported evaluations clearly reported;
- **CN midpoint values** `get_quantity(..., CN=True)`
  for dissipation and port power on each accepted step;
- original per-step and cumulative identities,
  including `ΔH + dt*(D - P_signed)=0`
  and each subsystem's separate balance;
- all entire FE coefficient paths, actual field traces,
  Hamiltonians and public PNG/ParaView outputs;
- cross-LU full-path FE and power-observable
  agreement with an explicit `1e-8` comparator gate
  (MUMPS vs SuperLU_DIST, same fixed CN clock).

Diagnostics are **not automatically a physical-model
PASS**. Confirm energy and interface residuals
against explicit tolerances and physical meaning,
then qualify or correct the model. In particular,
the published article establishes a continuous
and semidiscrete identity, not an unconditional
PETSc time-integration proof for this different
specific weak form.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
SCRIMP_V130_2D_SCOPE=heat_wave \
bash scripts/run-v130-example-2d-lu.sh
```

**Status: source model changed; no genuine native
result yet.** `examples/published/` and every
historical local simulation output remain untouched.
