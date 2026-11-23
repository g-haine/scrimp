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

## Inventory (source inspected at db6ad7a)

Eight top-level executable examples:

| Script | Solver currently selected | Qualification action |
|---|---|---|
| `examples/heat.py` | no override; MUMPS default | Run physics + I/O |
| `examples/heat_wave.py` | hardcoded `superlu` | Default MUMPS + optional solver choice |
| `examples/mass_spring_damper.py` | hardcoded `superlu` | Default MUMPS + optional solver choice |
| `examples/nonlinear_pendulum.py` | hardcoded `superlu` | Default MUMPS + optional solver choice |
| `examples/shallow_water.py` | hardcoded `superlu` | Default MUMPS + separate nonlinear DAE/Jacobian gate |
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

1. Add explicit opt-in solver choice to the four hardcoded
   top-level scripts while retaining MUMPS as the no-argument
   reference. Do not hardcode `superlu_dist` on all systems:
   users without that PETSc backend must still be able
   to run basic examples with their installed solver.
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
