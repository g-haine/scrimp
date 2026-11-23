# SCRIMP v1.3.0 — Issue #29: accepted-state cold-restart qualification

**Qualification date:** 23 November 2026. **Scientific acceptance:** PASS in the tested PETSc TSCN scope. **Release v1.3.0:** not yet implied by this gate.

## Scope and implementation

The production `DPHS.save_checkpoint(path)` / `DPHS.load_checkpoint(path)` API stores and restores a **complete accepted GetFEM finite-element state** (differential and algebraic coefficients), with its physical time and verified model/mesh/DoF/control provenance. The write is versioned and atomic; the load is fail-closed. After loading, SCRIMP recreates PETSc TS/operators, checks the original algebraic residual `F_a(t*, z*, 0)`, and resumes at the **absolute physical** time `t*`, without pseudo/BDF initialization or algebraic projection. Matching only the number of mesh cells is not sufficient: geometry and global coefficient ordering must match.

This is a **cold accepted-state restart**, **not** a saved PETSc TS/KSP/MUMPS process, factorization, BDF2 or IMEX history. Qualification covers CN with accepted endpoints and the exact original models below. It does not establish a universal high-index DAE initializer, exact multistep continuation, long-time SWE stationarity, or general nonlinear-mass energy conservation.

Relevant code: `scrimp/dphs.py`, `scrimp/utils/dae_checkpoint.py`; qualification tools: `sandbox/probe_v130_checkpoint_linear_dae.py`, `sandbox/compare_v130_checkpoint_linear_dae.py`, `scripts/run-v130-checkpoint-linear-dae.sh` and `tests/test_v130_checkpoint_linear_dae.py`.

## Native models and acceptance

1. **Original curved-P2 dissipative Saint-Venant spinup**, complete-state production API: previously archived **16/16 native checkpoints/restarts PASS** across both meshes and the required cross-MPI configurations. This is a real physical transient, not a proof of a stationary rotating interior; see the dated source evidence in `V130_RELEASE_SPRINT_2026W48.md`.
2. **Independent mixed-boundary 2-D linear wave pH-DAE**, original physical GetFEM/PETSc model, checkpoint at `t*=0.020`, `dt=0.002`, continuous reference `0→0.024`, seed `0→0.020`, cold continuation `0.020→0.024`. Two seed runs (MPI 1/2), two uninterrupted runs (MPI 1/2), four genuine cross-MPI restarts (1→1, 1→2, 2→1, 2→2): **8/8 native segments PASS**. Original accepted trajectories contain **426 FE DoFs** (189 dynamic, 237 algebraic), positive PETSc convergence reasons, no rejected steps or SNES failures; original `F_a` and true full-step CN residuals are checked at all applicable samples. The full initial restart state matches the saved source **bitwise**; this does not claim bitwise identity for the later numerical trajectory.

**Host evidence**, retained outside Git, immutable original run: `sandbox/outputs/v130-issue29-linear-DAE-e795ac44-Bo62JD0X/` (native source `e795ac44`, pinned `scrimp-mpi:v1.3.0-rc3`). The corrected host comparator supports original native report schema v1 and corrected future schema v2; the archived v1 trace/report are **not rewritten**.

User-confirmed final command on 23 November:

```bash
python tests/test_v130_checkpoint_linear_dae.py
bash scripts/run-v130-checkpoint-linear-dae.sh --compare-only \
  sandbox/outputs/v130-issue29-linear-DAE-e795ac44-Bo62JD0X
```

**Host tests:** `Ran 13 tests ... OK`. **Real native comparison:**

| Source → target MPI | `||z_cont(tf)-z_restart(tf)||_2` | Absolute difference in signed energy defect |
| --- | ---: | ---: |
| 1 → 1 | 1.1046626e-13 | 5.4210109e-20 |
| 1 → 2 | 5.8290712e-13 | 2.7105054e-19 |
| 2 → 1 | 6.4287550e-13 | 1.6263033e-19 |
| 2 → 2 | 4.1801798e-14 | 1.0842022e-19 |

All four comparator gates report **`LINEAR DAE CHECKPOINT PASS`**.

## Physical port-work correction and discrete energy identity

`Port.compute()` evaluates `P[k]` for `k>=1` on the midpoint of **step k**, while `P[0]` evaluates the initial state instantaneously. Therefore `P[0]` must **not** be treated as the initial sample of an endpoint-trapezoidal quadrature, especially across a restart. For this model's left-flow-port sign convention the correct CN work and energy identity are:

```text
W[n] = -sum_{k=1..n} (t[k]-t[k-1]) * P[k]
H[k]-H[k-1] + (t[k]-t[k-1])*P[k] = 0
```

The initial v1 native report used a trapezoidal rule and an incorrect sign; the corrected comparator independently recovers the proper per-step work from **unmodified native port-power arrays**, checks each original CN energy increment, and compares continuous versus split trajectories. The native linear-wave discrete power identity closes near machine precision (~1e-19 in the archived run); it is not merely a cancellation of two equally wrong energy defects. A new native v2 report writes the corrected quadrature directly. No change to the production `Port.compute()` semantics is made by this qualification.

## Negative controls and release interpretation

Host tests cover invalid/missing segments, corrupted source/mismatched model provenance, altered restored coefficients, algebraic constraint violation, incorrect boundary power and restart projection. The earlier production checkpoint gates additionally test accepted-state requirements, schema/checksums, mesh/FE/control/clock incompatibilities and unsupported integrator modes (see `tests/test_v130_checkpoint_state.py` and `tests/test_v130_checkpoint_public_api.py`).

**Decision:** Issue #29's **narrow scientific/API acceptance criteria are satisfied** by these native and host results. Other v1.3.0 workstreams and complete repository CI remain independent release gates. Any future broader integrator-history or nonlinear/high-index initialization feature requires a separate qualification.
