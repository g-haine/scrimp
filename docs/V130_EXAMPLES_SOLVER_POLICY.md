# SCRIMP v1.3.0 — solver policy and systematic example qualification

**23 November 2026 — actual MPI sparse LU preflight 7/7 PASS; all archived publication code protected by an immutable, byte-for-byte Git baseline.**

## 24 November — First tutorial: 1D vibrating string (NATIVE PASS)

The first scientific example in the current Sphinx tree is
`docs/source/started.rst`, not the 2D `examples/wave.py`.
The former exploratory 1D sandbox script did not match the
declared outward traction/positive observation port signs and
is replaced by the maintained `examples/wave_1d.py`.
The canonical initial strain is `q0=4*x*(1-x)`, exactly P2
and compatible with zero initial Neumann inputs.

**NATIVE PASS 24 November 2026:** MUMPS MPI1/MPI2,
SuperLU MPI1, SuperLU_DIST MPI1/MPI2: all **five actual
GetFEM/PETSc** cases, 128 DOFs, ten accepted fixed CN steps,
strict algebraic DAE projection, nontrivial signed port
work and maximal measured CN step energy defect
**1.791e-16**. Full accepted trajectories agree to
at most **2.22e-15**; all q/p PVD/VTU series are present.
Evidence `v130-first-tutorial-wave1d-2a170c82-cwkGUoXC.zip`,
source `2a170c8247f0f8aff09331747a3bb4812a04324f`.
Not a performance benchmark.

The Sphinx beginner chapter now contains actual native
energy-work and animated-strain/reconstructed-displacement
SVG illustrations. The optional rendering script generates
scientific PNG/SVG plots, user CSV/NPZ data and GIF from
the original accepted FE states and GetFEM outputs.
The separate **1.5 s long visual demonstration now
also passes natively** (MUMPS MPI1 only, source
`7e2e95a`, 150/150 accepted CN steps and
151 saved snapshots, maximum true per-step energy
defect `3.434752482434078e-16`). Its complete
`q`/`p` VTK output consists of **302 checked VTUs**;
the original 151-frame GIF and energy plots have
been inspected. Compact, data-derived SVG animations
and energy figures are included in the first Sphinx
chapter; full binary PNG/GIF/CSV/NPZ evidence remains
user-reproducible and is not bulk-committed.
The native proof archive is
`v130-wave1d-long-doc-7e2e95a4-75oCWG1T.zip`.
This long run is **not** an MPI2/SuperLU backend
qualification; those claims remain tied to the
short five-solver installation acceptance.

Scientific links: Haine–Matignon–Serhani (2023),
DOI 10.4208/ijnam2023-1005, and
Ferraro–Fournié–Haine (2024),
DOI 10.1016/j.ifacol.2024.08.267.
This short tutorial is an **illustration**, not a paper
figure or convergence-study reproduction.

Details: `docs/V130_EXAMPLE_01_WAVE_1D.md`.
Eight immutable `examples/published/` sources untouched.

## Pedagogical evidence contract for EACH subsequent example

The first vibrating-string example now establishes the standard
for the remaining release examples. Work through the **current Sphinx
order**, one example at a time; do not mechanically preserve an old
mesh, loading or source snippet when it obscures the mathematical
phenomenon. A model and its documentation are one jointly reviewed
scientific artifact. When a configuration changes, its mathematics,
source and published figure provenance must change together.

For every example, require:

1. **Scientific question and reproducible physics:** identify
   domain, state/co-state, physical parameters and dimension,
   boundary/control ports, expected mechanical/thermal balance
   and the numerical property being illustrated. Explicitly
   distinguish the genuine example from any cited ISAE-SUPAERO
   publication's original numerical experiment.
2. **Learnable SCRIMP API:** show real code that a first-time
   user can adapt; avoid burying the port-Hamiltonian mathematics
   under PETSc/VTK plumbing. Resolve conflicting signs, GWFL
   matrix conventions, incompatible initial data and FE spaces
   before reporting physical PASS.
3. **Actual install acceptance:** test MUMPS and SuperLU_DIST
   with 1/2 MPI ranks, and serial SuperLU with 1 rank where
   the model supports that backend. Compare original accepted
   FE states, converged PETSc time steps, nontrivial measured
   energies and port powers. Do not declare success solely
   from a LU factorization.
4. **Useful original outputs:** preserve native VTK/PVD as
   appropriate, and supply documented PNG/SVG figures, CSV/NPZ
   data, an optional animated GIF/SVG for online Sphinx and
   JSON provenance (image revision, source commit, solver,
   mesh, dt, accepted/rejected steps). Explain any
   postprocessed/reconstructed variable separately from
   the solved state, including its gauge and approximation.
5. **Evidence and maintenance:** retain one canonical source
   per example, remove obsolete duplicates once replacement
   is natively qualified, keep `examples/published/` byte-for-byte
   unchanged, validate both HTML/LaTeX Sphinx links, and
   record the exact native PASS/FAIL archive, including
   known scientific limits. Keep time-expensive demonstrations
   separate from the fast install smoke.

The first example satisfies the five-case native solver
and scientific physical gates and also has a **separately
validated 1.5-second MUMPS/MPI1 documentation simulation**.
The reviewed full GIF and native VTK data remain reproducible;
Sphinx ships compact energy and animated SVG figures from
that genuine long run. The extra long visualization is
not part of the five-backend installation smoke.

## Decision: stable release now, optional solver science later

The already-qualified `scrimp-mpi:v1.3.0-rc3` PETSc environment
has been queried with the **actual native** `PETSc.Sys.hasExternalPackage`:

| PETSc LU backend | Found in qualified rc3 | Intended use |
|---|---|---|
| `mumps` | True | Qualified default, serial and distributed |
| `superlu` | True | Sequential LU only |
| `superlu_dist` | True | MPI distributed LU: MPI1/2/4 tested PASS |

**No PETSc, GetFEM or Docker rebuild is necessary** to test
`superlu_dist` in the existing image. Package detection and real seven-case LU accuracy are now PASS;
this is still **not** a SCRIMP/GetFEM model-level qualification. The qualified solver policy remains `mumps` by default.
SCRIMP `DPHS.set_time_scheme()` already sets MUMPS only if no
user `pc_factor_mat_solver_type` was given, preserving user choice.

PETSc documents the parallel LU backend as
`pc_type=lu`, `pc_factor_mat_solver_type=superlu_dist` with
sparse AIJ matrices. The raw `superlu` backend is not
a parallel LU substitute for `superlu_dist`.

**Do not rebuild the native dependencies** merely because an old
example hardcodes `superlu`. Do not make a new mandatory
SuperLU_DIST vs MUMPS timing competition a blocker to v1.3.0.
After the example gates close, a separately versioned
reproducible performance campaign is worth doing on a curated
model selection (wave, heat, coupled heat-wave, genuinely nonlinear
Saint-Venant) with mesh/dt, memory, MPI rank placement and
symbolic/numeric factorization cost all controlled.
Such results are scientific evidence, not a blanket claim that
one backend is always faster.

## Stage B1 — update top-level examples (source STAGED; native pending)

The **four live, non-published** examples previously hardcoding
sequential `superlu` now take an optional **keyword-only**
`factor_solver=None` parameter:

```python
from examples.mass_spring_damper import mass_spring_damper

system = mass_spring_damper(final_time=0.04, time_step=0.01)
# Omitting factor_solver preserves SCRIMP's current MUMPS default.
system_dist = mass_spring_damper(
    final_time=0.04, time_step=0.01,
    factor_solver="superlu_dist",
)
# Serial SuperLU remains an explicit user choice on one rank.
```

Exactly the same optional parameter is available on
`nonlinear_pendulum()`, `heat_wave_eq()` and
`shallow_water_eq()`. No solver whitelist is hardcoded
into the example functions; the user's PETSc backend
selection passes through unchanged. For the two 2D
functions, explicit `final_time`, `time_step` and
`save_step` options permit bounded tests; **default
physical durations and initial timesteps remain configurable
without changing the original calls**. Heat–Wave now
**intentionally changes its FEMs, coupling signs,
external boundary condition, DAE initializer and
integrator** to make the example physically relevant;
see [V130_HEAT_WAVE_PHYSICAL_AUDIT.md](V130_HEAT_WAVE_PHYSICAL_AUDIT.md).
All four top-level entry points still
work without extra parameters and use MUMPS by default.
The other four top-level examples already inherit MUMPS
and have not been edited.

**Do not mistake this source change for a numerical PASS.**
Source-only tests validate the exact function signatures,
original defaults, preserved production MUMPS option
behavior and continued immutable `examples/published/`
guard. They cannot certify an actual nonlinear DAE
integration or solver performance.

### Bounded genuine physical pilot: original lumped systems

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-example-lumped-lu.sh
```

The driver overlays current Python source read-only onto
the previously qualified rc3 image, without rebuilding
GetFEM or PETSc. It executes **12 independent original
physical GetFEM/PETSc simulations**:

- mass–spring–damper and fully nonlinear pendulum,
  original weak forms, explicit pseudo initial relaxation,
  BDF2, `t=0..0.04`, `dt=0.01` (four accepted steps);
- for each example: default MUMPS at MPI1, explicitly
  selected MUMPS at MPI1/2, explicit SuperLU_DIST at
  MPI1/2 and sequential SuperLU at MPI1;
- strict PETSc convergence and accepted-time/step guards,
  no rejected steps or SNES failures; actual saved
  full-DAE state and model-specific Hamiltonian
  identities, nontrivial dynamics, rank-local state
  agreement and cross-solver whole-trajectory parity;
- completely separate host-owned working/output directories;
  output reports, original NPZ full states and complete
  host logs stored in a unique ZIP; **nothing is mounted
  for writing from `examples/published/`**.

This is an intentional **short pilot**: full default
8/12-second trajectories, the eight-model release matrix,
real 2D example output exports, extended MPI4 and
nonlinear Saint-Venant remain distinct subsequent gates.
**Native status: PENDING qualified-host execution.**

**New generic intrinsic-FEM coefficient preflight (source staged):**
`scripts/run-v130-intrinsic-fem-initialization.sh`
tests RT0/1/2, BDM1/2 and GetFEM's lowest-order
`FEM_NEDELEC(2)` on real filtered simplicial
MeshFem variables. Any truly unavailable element
is recorded as UNAVAILABLE, never passed; successful
native coefficient preservation is a much narrower
claim than full H(div)/H(curl) PDE support.
`examples/published/` remains immutable.
See [V130_HEAT_WAVE_PHYSICAL_AUDIT.md](V130_HEAT_WAVE_PHYSICAL_AUDIT.md).

### Stage B2 — Heat–Wave physics revised to RT/CN (NATIVE PENDING)

The old Heat–Wave example was not physically qualified:
it imposed zero *traction* on the outer wave boundary
rather than the specified Dirichlet condition, had
a sign convention enforcing `T=-p` at the interface
instead of `T=p`, used heat `CG2` flux and `DG1`
wave stress, and integrated in adaptive BDF2.
Its first native run did complete t=0.004 with
39 accepted, 4 rejected steps and zero nonlinear solver
failures; the first fixed-four-step smoke was an
incorrect **test protocol**, not evidence of divergent
MUMPS nor proof of correct physical transmission.

Because the purpose of this campaign is to provide
scientifically meaningful user examples, we have now
deliberately revised **only the live, non-published**
`examples/heat_wave.py` to:

- `T:DG0`, `e_Q:RT0`, `p:CG1`,
  `q:RT0`, interface/boundary traces `CG1`;
- fixed `CN` with the production constraint-consistent
  `algebraic` initializer and unchanged numerical
  default values `dt=.001,t_f=2.,dt_save=.05`;
- correct heat temperature input and common mesh-normal
  orientation, wave traction transmission and
  exterior zero-velocity Dirichlet enforced by a
  Lagrange multiplier;
- explicit MUMPS-default or optional PETSc
  `factor_solver`, with no changes to SCRIMP itself.

**Important scientific caveat:** the wave still uses
its established **grad–grad** partitioned form,
hence `p:CG1`, rather than the article's distinct
**div–div/DG0** wave discretization. This is not
claimed to reproduce the original article's numerical
method. Do not use DG0 wave velocity with the old
grad-grad bricks: its gradient would vanish elementwise.

The revised probe measures the actual original GetFEM
Hamiltonian contributions, CN midpoint heat dissipation,
raw interface port values with correct physical signs,
the two separate subsystem energy balances, outer
boundary work and optional traces. It uses
`get_quantity(...,CN=True)`, so no artificial
trapezoidal integration of already-stage-associated
CN power is introduced. It checks complete FE
trajectory and physical traces for MUMPS vs
SuperLU_DIST at **MPI=1 only** (parallel
ParaView I/O is not yet qualified).

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
SCRIMP_V130_2D_SCOPE=heat_wave \
bash scripts/run-v130-example-2d-lu.sh
```

The job must accept exactly four CN steps over 0..0.004,
with five stored times and no rejections.
A backend parity PASS **must not** be mislabeled
as physical correctness until its independent
energy and interface residuals are inspected.

The inviscid `shallow_water_eq()` is a separate,
**unchanged adaptive BDF2** model. It may be
run later using `SCRIMP_V130_2D_SCOPE=shallow_water`.
Do not simultaneously exercise unchanged public
ParaView writers on multiple MPI ranks; this
requires a dedicated safe export implementation.
`examples/published/` and all original historical
outputs remain strictly untouched.

Detailed formulation, comparison with the 2022
paper, and scientific acceptance criteria are in
[V130_HEAT_WAVE_PHYSICAL_AUDIT.md](V130_HEAT_WAVE_PHYSICAL_AUDIT.md).
**Native status: PENDING.** Second native attempt
at source `f97fcd1f`: original RT0 q 3111
coefficients initialize correctly, and full
7689-DOF GetFEM mass/stiffness assembly completes.
However the linear-only algebraic projector rejects
the constitutive feedback relations because
`set_control("...","Y_w")` and `"Y_T"`
were incorrectly marked `linear=False`.
They are bare signed model-variable **linear** feedback.
The core classification now selects the exact
constant linear GetFEM tangent for these controls,
with fail-closed nonlinear behavior for other
expressions and dedicated actual-method regression.
Numerical projection and physical CN balance remain
**unqualified**. First native MUMPS
attempt at source `8286f18d` established
DG0/RT0/CG1 MeshFem creation but was halted
by SCRIMP's *core* filtered RT0 initialization:
3111 wave RT DOFs with qdim=2 were erroneously
reshaped to (2,1555). This is now corrected
generically in `DPHS.set_initial_value()`,
using the raw 1D intrinsic-FE coefficient order,
with new host-only negative regressions.
No actual CN trajectory/power balance has
been obtained yet. Details in
[V130_HEAT_WAVE_PHYSICAL_AUDIT.md](V130_HEAT_WAVE_PHYSICAL_AUDIT.md).


## Inventory (source at current release/1.3.0-preparation; 23 November 2026)

Eight top-level executable examples:

| Script | Solver currently selected | Qualification action |
|---|---|---|
| `examples/heat.py` | no override; MUMPS default | Run physics + I/O |
| `examples/heat_wave.py` | now default MUMPS; optional backend | Native 2D qualification pending |
| `examples/mass_spring_damper.py` | now default MUMPS; optional backend | Short native pilot pending |
| `examples/nonlinear_pendulum.py` | now default MUMPS; optional backend | Short native pilot pending |
| `examples/shallow_water.py` | now default MUMPS; optional backend | Separate nonlinear DAE/Jacobian gate pending |
| `examples/wave.py` | no override; MUMPS default | Run physics + I/O |
| `examples/wave_coenergy.py` | no override; MUMPS default | Run physics + I/O |
| `examples/wave_dissipative.py` | no override; MUMPS default | Run physics + I/O |

### ABSOLUTE IMMUTABILITY RULE — `examples/published/`

All **eight Git-tracked files** inside `examples/published/`
are historical research snapshots exactly as they existed for
their associated articles. **NEVER edit, reformat, relocate,
delete, overwrite, add import shims, change solver options, or
commit new source files inside this directory.** Locally generated
outputs, figures, videos and ZIP evidence under historical
run directories are permitted and must be preserved exactly.
No automatic conversion,
`ruff/black` rewrite, license/header update, new entrypoint,
numerical hotfix or solver migration is permitted there.
The original `MCSS_2024/shallow_water.py` explicitly chooses
`superlu`; leave that exact choice intact even if it does
not execute with multiple ranks in the current Docker image.

Snapshot preservation is enforced by
`python scripts/check-v130-published-immutable.py` and
a dedicated pull-request/push CI. The guard checks the
**eight pinned Git blob IDs, the Git index, file permissions,
and the actual bytes of versioned sources**. It rejects
omissions, added/renamed **tracked sources**, and content changes
to tracked sources. It deliberately does **not** enumerate
or reject untracked local results. Its lock baseline is
`b6fe4933415fb19eb02a8448a9fc587432196227`.
The pinned contents are also verified by negative unit tests.
The guard is mandatory for release and must never be
bypassed to permit an archive SOURCE edit. The root `.gitignore`
also excludes `examples/published/**/outputs*/` and
`examples/published/**/outputs*.zip` from Git's untracked
reporting; this does not delete, move, stage, or rewrite them.
IMPORTANT: if any generated output has already been staged
by the user it must be inspected carefully; the guard rejects
all additions to the tracked source file inventory.

If a publication example needs modern API adaptations or
MUMPS/SuperLU_DIST experiments, construct a **new, clearly
named independent script outside `examples/published/`**
(e.g. under `scripts/`, `tests/`, or a dedicated modernized
examples directory). Link to the original source and record
the changes explicitly. Execute legacy scripts unchanged
where technically possible; otherwise record the historical
incompatibility without changing the snapshots.
The notebooks and helper scripts are inventoried separately;
do not blindly treat every helper as an executable model.

### Stage A — minimal actual PETSc backend preflight (7/7 NATIVE PASS)

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-petsc-lu-backends.sh
```

This does **seven** independent native factorization solves:
`mumps` and `superlu_dist` on MPI 1, 2, 4; sequential
`superlu` on MPI 1 only. Every case uses a true sparse
complex, nonsymmetric, diagonally dominant `MATAIJ` problem
of size 384, with actual distinct ownership on MPI ranks,
exact known complex vector, relative solution L2 error
and relative `Ax-b` residual gated at `1e-11`.
It validates the **actual** factor backend as reported
by PETSc, not just `hasExternalPackage()`. Missing factors,
non-convergence, wrong result and mismatched MPI matrix
format are hard failures. It stores the original host logs
and JSON for each case in a unique, immutable output
directory and makes a standalone ZIP.

**Stage A NATIVE PASS (user host; full original archive independently
read):** `scrimp-mpi:v1.3.0-rc3` at OCI revision
`1d6d0c90115213298047b28237a7f38087563a35`;
source `b6fe4933415fb19eb02a8448a9fc587432196227`.
Original ZIP `v130-petsc-lu-preflight-b6fe4933-0LlCfedE.zip`
contains 7 native PETSc JSON results, 7 corresponding
native logs, manifest and summary. All entries pass CRC.
Actual results:

| PETSc LU backend | MPI ranks accepted | Largest relative solution error | Largest relative Ax-b residual |
|---|---|---:|---:|
| `mumps` | 1, 2, 4 | 1.9386099094915744e-16 | 3.244033106455509e-16 |
| `superlu` (sequential) | 1 | 1.523376208552414e-16 | 2.449694660886362e-16 |
| `superlu_dist` | 1, 2, 4 | 1.7609301536061955e-16 | 2.947112437963227e-16 |

All 7 actual factorization paths returned PETSc KSP
convergence reason 4, correct requested factor backend,
and the two- and four-rank cases used genuine `mpiaij`
matrices. Largest relative error
`1.9386099094915744e-16`, largest relative
residual `3.244033106455509e-16`; threshold
`1e-11` was met with margin. This is **not** a
timing benchmark: the 384-dimensional test problem is
far too small to draw speed or scaling conclusions.

A successful Stage A means **PETSc native backend LU
compatibility on one small test matrix**. It does NOT mean
either backend is already qualified on any SCRIMP model,
that the solve is scalable, or that it preserves SCRIMP's
original DAE model behavior or power balance.
In particular, `DPHS` has MUMPS-specific factor reuse
logic; an actual SCRIMP run is required to validate its
non-MUMPS path.

### Stage B — example acceptance, P0 for release

1. **SOURCE STAGED:** four formerly hardcoded top-level
   scripts accept optional `factor_solver` with production
   MUMPS default. Native acceptance (first lumped pilot,
   then distributed PDE examples) is still mandatory.
   Do not require SuperLU_DIST on other installations.
2. Run each shipped top-level model with MUMPS on MPI 1;
   use MPI 4 only where the model and outputs have an
   explicitly tested MPI contract. Store PETSc status,
   rejected steps, physical time, algebraic residuals,
   model-appropriate power/energy, and actual output files.
3. Run a **limited** representative set with
   SuperLU_DIST MPI 2/4 if Stage A passes, compare
   *numerical trajectories and power observables* to MUMPS.
   Investigate/mark original model issues, not just
   factorization timings. `superlu` can remain an
   optional serial reference where supported.
4. The original `examples/published/` snapshots are **NEVER
   edited**. If a historical scenario must be compared,
   write an independent modernized reproduction OUTSIDE that
   directory with all departures recorded. Keep its large
   original time horizon separate from short acceptance tests.
5. Keep #30's `Port.compute` typed sampling policy and
   #6's genuine physical MPI and checkpoint fidelity
   tests unchanged; no exact BDF2 port-power identity
   is assumed in these examples.

### Stage C — research benchmark, explicitly outside v1.3.0 gate

With exact same physical models, mesh, MPI topology and
time step, record solver setup vs factorization vs KSP
solve, symbolic reuse behavior, memory peak and repeat
variability. The small Stage A matrix has almost no
predictive power for the large nonlinear pH-DAEs.
Report cases where SuperLU_DIST is faster/slower or
less/more robust than MUMPS; do not tune after seeing
a single timing. Keep any alternate Docker image
separately tagged and never overwrite the rc3 image
before fully qualifying a new candidate.
