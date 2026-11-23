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
