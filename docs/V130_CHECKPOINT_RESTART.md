# SCRIMP v1.3.0 — Accepted-state checkpoint/restart (P0, issue #29)

## Scope and current evidence

The public methods `DPHS.save_checkpoint(path, *, extra=None)` and
`DPHS.load_checkpoint(path)` implement a *cold restart* from a **complete
accepted state** of a port-Hamiltonian DAE, not serialization of
PETSc's integrator history. The archive includes all differential
variables, algebraic efforts and boundary multipliers, their original
absolute physical time, exact FE mesh/convex/DoF numbering and
weak-form/parameter/control fingerprint and SHA-256 integrity hashes.

**Native qualified so far:** original curved-P2 viscous Saint-Venant,
2 FE meshes, MPI1/2 source/target cross-rank permutations (8
restarts; **16/16 real PETSc/GetFEM runs PASS**), continuous
`0→.024` vs seed `0→.02`+new PETSc `.02→.024`.
Full original FE state discrepancy ≤2.26e-13, energy-balance
difference ≤3.03e-14. Source
`260aeb558ffec8783044e7c93c2f02504d75eb41`,
immutable `scrimp-mpi:v1.3.0-rc3` OCI
`1d6d0c90115213298047b28237a7f38087563a35`.
This qualifies the *actual public API for this SWE model only*.

**Independent linear pH-DAE wave test:** source and host negative
controls are staged; native results **not yet available**.
See `scripts/run-v130-checkpoint-linear-dae.sh` and
`sandbox/probe_v130_checkpoint_linear_dae.py`.
A general release claim for #29 requires this second-model gate.

## A typical same-mesh restart

```python
# The user-defined build_model() reconstructs IDENTICAL GetFEM geometry,
# variable ordering, FEM families, bricks, parameters and control laws.
source = build_model(final_time=0.02)
source.solve()

# ALL MPI ranks call save_checkpoint() collectively.
# File must not already exist; publishing is atomic.
source.save_checkpoint(
    "checkpoints/accepted_spinup.npz",
    extra={"experiment": "LHMNC", "phase": "preparation"}
)

# A NEW process/DPHS with the original physical time-dependent control.
continuation = build_model(final_time=0.04)
# Configured TS type must be PETSc endpoint CN (ts_type="cn").
# load_checkpoint() reads/validates on all ranks, sets t_0=0.02,
# and disables the ordinary pseudo/BDF/algebraic initializer.
continuation.load_checkpoint("checkpoints/accepted_spinup.npz")
continuation.solve()
```

Call `load_checkpoint()` **after** the same FE model, physical
boundary law and future `set_time_scheme()` have been configured,
but **before** `solve()`. The checkpoint's time is preserved:
use `t-t_*` for plotting a relative phase, *not*
by resetting the PETSc physical clock. If the wall drive is
meant to remain rotating after `t_*`, encode that in the
original complete physical control law **before** the saved run
rather than silently changing the formula after loading:
the model fingerprint will refuse an incompatible law.

At load, SCRIMP checks physical model/FE identity, native MPI
agreement and source-state checksum; the actual full GetFEM
coefficient vector is installed **after physical operator assembly
and before PETSc TS**, and the original differential/algebraic
partition is used to evaluate `F_a(t_*,z_*,0)`. Rejection does
not silently project an inconsistent differential state.
A moving physical DAE state may have
`F_d(t_*,z_*,0) ≠ 0`; only algebraic compatibility is required.
All PETSc TS/KSP/MUMPS operators are rebuilt for the new run.

## Currently supported / not supported

- ✅ `ts_type="cn"` endpoint-trapezoidal PETSc
  **cold state** restart after a complete successful previous
  TS solve; original full native FE coefficient state preserved.
- ✅ Different MPI rank count, **provided actual FE geometry,
  basis and global DoF ordering remain identical**.
- ✅ Atomic versioned NPZ with SHA-256 full-state checksum and
  model fingerprint; existing accepted checkpoint never overwritten.
- ❌ Unknown, changed or independently regenerated mesh/global
  numbering, altered coefficients/weak form/boundary-control law,
  incomplete/unconverged source TS, missing algebraic state.
- ❌ Full PETSc integrator stage/history continuation for
  BDF2/multistep/IMEX, arbitrary index/high-index inconsistent
  q0 regularization and exact bitwise continuation of all
  subsequent floating-point steps.
- ❌ Automatically reconstructing meshes from the checkpoint:
  source model/mesh construction remains the user's responsibility.
- ⚠️ The MVP's fixed absolute algebraic tolerance
  (1e-9 in native FE coefficient space) is not yet a
  model-independent scale-invariant condition test.
  A second linear DAE and more negative gates are required
  before calling this generally release-qualified.
- ⚠️ A checkpoint at t*=.02 on the tested Saint-Venant
  model is the **end of the imposed wall acceleration**,
  not automatically a stationary interior rotating regime.

## Native test commands

```bash
cd ~/Simulations/scrimp
git pull --ff-only

python -m py_compile \
  sandbox/probe_v130_checkpoint_linear_dae.py \
  sandbox/compare_v130_checkpoint_linear_dae.py \
  tests/test_v130_checkpoint_linear_dae.py
python tests/test_v130_checkpoint_linear_dae.py
bash -n scripts/run-v130-checkpoint-linear-dae.sh

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-checkpoint-linear-dae.sh
```

The 8-run linear DAE campaign stores immutable host-owned evidence
under `sandbox/outputs/v130-issue29-linear-DAE-*/`.
If interrupted, use its `--resume <that directory>` option;
the runner independently checks real native provenance before
skipping completed cases. No old SWE run needs to be repeated.

See `docs/V130_DEVELOPMENT_ROADMAP.md`,
`docs/V130_RELEASE_PLAN.md`, and issue #29 for the
release-blocking and scientific acceptance status.

### 23 November 2026 — linear second-model host preflight: 6/7; diagnostic specificity corrected (native pending)

The user ran the staged second-independent-model
`tests/test_v130_checkpoint_linear_dae.py`
at source `5c9d924`:
Python `py_compile` for both linear DAE probes and
test file PASS, `bash -n` for the 8-case native
runner PASS; **7 synthetic host tests executed,
6 PASS, 1 FAIL**. The only failing test,
`test_rejects_wrong_model_fingerprint`,
correctly rejected a tampered resumed segment
report, but the comparator merged its model
SHA mismatch into a broad
`Actual DPHS linear DAE resume projected/changed
q0 or violated Fa` error, whereas the negative
test required the more accurate
`physically different model` reason.

**Fix staged on the source branch**:
`sandbox/compare_v130_checkpoint_linear_dae.py`
checks the accepted checkpoint's exact original
FE/weak-form model fingerprint against the
restored report *before* inspecting full state
byte equality or DAE constraints; seed source
fingerprint checks are likewise isolated.
This retains every previous strict
fail-closed safeguard. New explicit tests
separately cover (a) a **syntactically valid
but physically different 64-hex-digit model
hash** in restart, (b) wrong seed-model
fingerprint, and (c) a genuine false
`no_pseudo_BDF_or_projection` restart flag
with otherwise identical physical model.
There are now **9 synthetic host regression
tests** for the second model. No SWE/native
physics, PETSc solver, checkpoint format or
numerical acceptance threshold changed.

**Status**: new 9/9 host tests are NOT YET
confirmed by a rerun; **zero native second
linear DAE trials have yet been observed**.
Issue #29 remains OPEN. Actual production
DPHS checkpoint API's earlier original
Saint-Venant native **16/16 PASS** evidence
remains unchanged.

Revalidate before Docker:
```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_checkpoint_linear_dae.py
bash -n scripts/run-v130-checkpoint-linear-dae.sh

# Only if the host suite passes 9/9:
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-checkpoint-linear-dae.sh
```

### 23 November 2026 — First native independent linear wave DAE: 10 accepted steps; PETSc index ABI mismatch in diagnostic (fixed)

Host source `36aca5319a9770bc72157cf4978c18ab3d05f852`,
immutable `scrimp-mpi:v1.3.0-rc3`, actual partial
host evidence
`sandbox/outputs/v130-issue29-linear-DAE-36aca531-aJYhtX7N/`
and complete uploaded ZIP of that partial directory
were inspected directly. **All 9/9 linear DAE synthetic
host tests PASS**, alongside 17/17 checkpoint format and
5/5 public API host tests.

The first physical **second-model** run
`wave_seed_mpi1` constructed the genuine
2D linear-wave pH-DAE with 426 full
GetFEM DOFs (h=.5 Rectangle, 22
elements), and PETSc CN accepted
**10/10 physical steps from 0 to 0.02**,
0 rejections/0 SNES failures
(`Elapsed time: 0.01631s`).
It subsequently failed **after**
the accepted TS in the independent
post-solve audit `original_IFunction`,
while filling newly constructed PETSc
vectors via `Vec.setValues`. Python
`np.arange(lo,hi,dtype=int)` had
`dtype('int64')` on the user's host;
the pinned native PETSc build uses
`PetscInt` 32-bit. petsc4py correctly
rejected unsafe conversion:
`TypeError: Cannot cast array data
from dtype('int64') to dtype('int32')
according to the rule 'safe'`.
This is NOT a time-solver convergence,
FE physics or checkpoint data error,
but **no native complete second-model
residual/energy/port work nor accepted
checkpoint was written** in this
first failed run. Therefore *native*
second DAE checkpoint restart
remains **0/8 qualified**.

**Source fix:** fill the owned
`Vec.setValues` indices with
`np.arange(lo,hi,dtype=PETSc.IntType)`
for both state and time derivative
PETSc vectors. `PETSc.IntType`
selects the exact native `PetscInt`
width (32- or 64-bit), independent
of Python word size. A new host
regression inspects the true
`original_IFunction` source AST
for this index selection and
both `vec` and `vel` use sites;
**10 host second-model tests now
await rerun**. Solver/GetFEM
weak forms, DAE tolerances,
time stepping and runtime
public checkpoint API unchanged.

Because the old native phase
produced only an incomplete log
and no accepted-state file,
do **not** use its directory with
`--resume`: the strict
evidence-preserving runner must
reject that incomplete case.
Run a **fresh** native pilot:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python -m py_compile \
  sandbox/probe_v130_checkpoint_linear_dae.py \
  sandbox/compare_v130_checkpoint_linear_dae.py \
  tests/test_v130_checkpoint_linear_dae.py
python tests/test_v130_checkpoint_linear_dae.py
bash -n scripts/run-v130-checkpoint-linear-dae.sh

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-checkpoint-linear-dae.sh
```

**Evidence standard unchanged:** host
tests and 10 accepted forward
steps alone do not qualify checkpoint
restoration; the 8 full independent
real GetFEM/PETSc linear DAE phases
and 4 true cross-MPI restarts must
pass original algebraic/TSCN,
all FE state, Hamiltonian and
boundary-work comparisons before
closing #29.
