# SCRIMP v1.3.0 — solver policy and systematic example qualification

**23 November 2026 — source preflight prepared; actual MPI direct-solver runs PENDING.**

## Decision: stable release now, optional solver science later

The already-qualified `scrimp-mpi:v1.3.0-rc3` PETSc environment
has been queried with the **actual native** `PETSc.Sys.hasExternalPackage`:

| PETSc LU backend | Found in qualified rc3 | Intended use |
|---|---|---|
| `mumps` | True | Qualified default, serial and distributed |
| `superlu` | True | Sequential LU only |
| `superlu_dist` | True | MPI distributed LU, runtime test pending |

**No PETSc, GetFEM or Docker rebuild is necessary** to test
`superlu_dist` in the existing image. These checks prove native
package registration, **not** successful factorization or runtime
stability. The qualified solver policy remains `mumps` by default.
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

**Published/research scripts:** `examples/published/MCSS_2024/shallow_water.py`
also hardcodes `superlu`, alongside historic BDF2, control
and integrator settings. Preserve published/archival defaults;
qualify an opt-in solver override separately so historical
science is not silently rewritten. Other published examples
and notebooks require a separate inventory of dependencies,
entry points and feasible MPI contracts; helpers such as
`u_tangent.py` are not blindly classified as executable
standalone solvers.

### Stage A — minimal actual PETSc backend preflight (PENDING)

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
4. Treat original published simulations separately:
   preserve their parameter provenance and do not
   confuse an unbounded 5-second, 50k-step historical
   scenario with a quick installation acceptance.
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
