# SCRIMP v1.3.0 — solver policy and systematic example qualification

**23 November 2026 — actual MPI sparse LU preflight 7/7 PASS; all archived publication code protected by an immutable, byte-for-byte Git baseline.**

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
physical durations, timesteps, meshes, weak forms,
controls and initialization all retain their original
values**. All four original top-level entry points still
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

### Stage B2 — real 2D Heat–Wave and inviscid Saint-Venant: original BDF2 adaptive behavior

The former hardcoded serial-SuperLU main scripts are now available
with MUMPS as the default, and `superlu_dist` as a pure
`factor_solver` opt-in. The first true model-scale test
is **not** an unrestricted solver benchmark:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
SCRIMP_V130_2D_SCOPE=heat_wave \
bash scripts/run-v130-example-2d-lu.sh
```

This launches the **actual shipped** `heat_wave_eq()` and
`shallow_water_eq()` models, with original weak forms,
original pseudo DAE initialization, meshes, constitutive
relations, boundary conditions and public figure/ParaView
export methods, with explicit `t=0..0.004`,
`dt=dt_save=0.001`. Other default values do not change.

The first native run `heat_wave_mumps_mpi1` reached
t=0.004 with **39 accepted adaptive BDF2 steps, 4 rejected
steps and no SNES failures**, 26,198 GetFEM DoFs,
and completed Hamiltonian and ParaView exports.
The old **test harness incorrectly required 4 fixed steps,
zero rejections and a uniform 5-point saved clock**.
This is a protocol false negative, *not* evidence of an
unsuccessful MUMPS solve or incorrect physical coupling.

The corrected gate checks actual accepted-time monotonicity,
physical final time, PETSc convergence and finite full DAE
states; it reports accepted and rejected steps explicitly.
For MUMPS vs SuperLU_DIST it gates the **full original
trajectory** only if both adaptive saved clocks match.
Otherwise it gates **original full FE and all physical
observables at identical start/end times**, explicitly
disclosing that this is endpoint-only and no
whole-trajectory parity has been claimed. All results
remain immutable and in separate workdirs.

`SCRIMP_V130_2D_SCOPE=heat_wave` is the first stage,
MUMPS and SuperLU_DIST MPI1; `=shallow_water`
is a later distinct stage, and `=all` runs both
only when the first-stage physics has been examined.
The heat-wave pilot records separate original
heat/wave energies, genuine GetFEM heat-flux
dissipation, raw pairings of both interface ports,
outer boundary power, candidate trace constraints
and BDF nodal trapezoidal balance **as diagnostics only**.
A numerical LU parity PASS is expressly NOT a physical
coupling model PASS. No physical qualification is
claimed at this point. See
[V130_HEAT_WAVE_PHYSICAL_AUDIT.md](V130_HEAT_WAVE_PHYSICAL_AUDIT.md).

For the later inviscid Saint-Venant stage, the probe
also saves the original GetFEM water-height integral
`∫h` and its time drift, without normalizing it
or interpreting approximate BDF energy conservation
as a factorization result. The separate scientific
audit and any model correction will follow the
native data, not precede them.

**MPI I/O scope warning:** current `DPHS.export_to_pv()`
opens the same `.pvd` on all ranks (and uses
root-only directory reset). The user previously
reported trouble locating generated outputs. It is
not yet proven safe to run these unchanged top-level
postprocessors simultaneously on MPI2/4. This
first stage therefore deliberately uses only MPI1
and independent host working directories. **No
simulation or document under `examples/published/`
is mounted, modified or touched.** After examining
both physical trajectories, qualify/fix generic
MPI export separately before extending actual
top-level 2D examples to MPI2/4.

**Scientific limitation:** numerical MUMPS/SuperLU_DIST
agreement does not itself validate correctness of the
original Heat–Wave coupling nor inviscid Saint-Venant
discretization. Physical Hamiltonian and volume
diagnostics are preserved for critical review; an
exact BDF2 power identity is **not** claimed.
On failure the launcher stops immediately and retains
the first original native failure logs and host
outputs; no forced PASS.

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
