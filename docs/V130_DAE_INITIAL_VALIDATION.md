# SCRIMP v1.3.0 — Physical DAE initialization and full-wave-state qualification

**Status (20 November 2026): first native Ubuntu 1/2-rank experiment
**COMPLETED**. MPI full-trajectory parity and algebraic consistency
**PASS**; prescribed dynamic initial data **FAIL**. The historical
pseudo-TS is not physically qualified as a constrained projection.
No production solver patch is approved on this evidence alone.**

This is the first scientific acceptance investigation after closing
the Ubuntu Docker BuildKit/native-stack gate. The reference production
model is the unchanged **`examples/wave.py`** (2D acoustic, 4 boundary
control ports, implicit pH-DAE, fixed-step Crank–Nicolson).
The validated Docker candidate `scrimp-mpi:v1.3.0-rc3` remains
untouched; only diagnostic sources are mounted read-only at **runtime**.

## 1. Scientific question: what does pseudo-TS actually preserve?

The published initial data are:

```python
q_0 = "[0., 0.]"
p_0 = "3**(-20*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)))"
wave.set_initial_value("q", q_0)
wave.set_initial_value("p", p_0)
```

The default solver performs a tiny PETSc pseudo-TS step before
starting the physical CN integration. The accepted pseudo-state
is now **explicitly synchronized to GetFEM**, fixing a known
cross-version state-regression side effect. However, an earlier
native diagnosis of a larger wave reported a nonzero accepted
`q` norm (~`4.71761e-3`), even though `q(0)=0` was prescribed.
v1.2.0 equivalence and strong scaling **do not** establish
physical fidelity of these prescribed initial data.

Two independent properties must be distinguished:

1. **Algebraic DAE consistency:** effort and boundary/control
   constraints are equations with no time derivative. Their
   rows of `IFunction(t0,z,zd=0)` must have a small residual.
   The audit classifies these rows from the actual SCRIMP
   ports and verifies their **actual mass-matrix entries are
   zero** before interpreting the residual as algebraic.
2. **Differential initial-data fidelity:** the accepted pseudo-state
   must preserve the prescribed `q(0)` and `p(0)` within an
   explicitly chosen absolute/relative tolerance. We compare
   values **before** and **after** the pseudo-TS, not merely after
   its synchronization to GetFEM.

It is possible to improve (1) at the expense of (2).
The production implementation must not be changed just to recover
old trajectories or reduce one norm at the expense of physics.

## 2. Diagnostic implementation, deliberately non-invasive

- `sandbox/audit_v130_wave_dae_initial.py` uses the already
  production-qualified `bench_acoustic_2d._build_wave` to construct
  the actual model. It intercepts only `DPHS.init_step` to capture
  the actual pre/post GetFEM vectors; no solver/PETSc API changes,
  no fake initial solution and no mutation of the accepted state.
- The original `IFunction` measures **algebraic-only residuals**
  with zero derivative and real PETSc owned residual rows. A
  structural mass-row check rejects accidentally including a
  dynamical equation in the algebraic set. The GetFEM state is
  restored after each diagnostic evaluation.
- The monitor saves complete fixed-CN state vectors at each time
  (not just a final-state or norm comparison). The diagnostic
  checks that physical TS step zero matches the accepted
  pseudo-state and records norms by named variable.
- `sandbox/compare_v130_wave_dae_initial.py` independently
  compares `before`, `after`, and **all six saved time states**
  for the **same mesh** at MPI 1 and 2. It verifies both
  individual time states and the whole trajectory, plus
  the algebraic residuals and initial-fidelity verdict.
- The default mode is **diagnostic**. If the known pseudo-TS
  drift is reproduced, the report states
  `PHYSICAL INITIAL-DATA FIDELITY: FAIL` without preventing
  completion of the MPI parity and algebraic analyses.
  `--fail-on-dynamic-drift` on the native audit and
  `--require-initial-fidelity` on the comparator are available
  once a physically corrected initialization is being qualified.

### Scope and what this does NOT prove

These experiments inspect algebraic rows and the state of the
published wave. They do not test every nonlinear DAE index or
prove a rigorous general projection theorem. A small weak-form
residual is not evidence of unique DAE projection or conservation
of an arbitrary initial state's pH energy. The separately existing
`sandbox/check_acoustic_energy_mpi.py` and
`sandbox/analyze_acoustic_power_balance.py` are the next gate
for **actual midpoint control-port power and Hamiltonian identities**;
the present benchmark must not substitute for them.

## 3. Run the first small native experiment (no Docker rebuild)

From the Ubuntu workstation with the already accepted rc3 image:

~~~bash
cd ~/Simulations/scrimp
git switch release/1.3.0-preparation
git pull --ff-only

# This workstation already runs Docker without sudo:
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-dae-initial-audit.sh
~~~

The runner uses `h=0.2`, `dt=0.002`, `t_f=0.01`,
**five fixed CN steps**, original wave equations,
`SCRIMP_GETFEM_MATRIX_MODE=distributed`,
COO bridge, compact matrix memory and MUMPS
`on_change` for **1 and 2 MPI ranks**.

The source of the qualified solver remains the immutable rc3 OCI
revision `1d6d0c90115213298047b28237a7f38087563a35`;
the current Git commit identifies the **separately mounted
diagnostic**. Only two explicitly named Python files are mounted
read-only into the Docker runtime. **No working tree or
`sandbox/outputs` enters the Docker build context.**

Everything goes under a new unique host directory:

```text
sandbox/outputs/v130-dae-wave-h0p2-<audit_SHA>-XXXXXXXX/
  manifest.json
  wave_1r.log
  wave_2r.log
  wave_initial_1r.json
  wave_initial_1r.npz
  wave_initial_2r.json
  wave_initial_2r.npz
  compare.log
  wave_initial_dae_summary.json
```

There is **no need to run v1.2.0 or recompute the ~3-million-DoF
strong-scaling campaign**. Preserve the developer's unrelated local
edit to `sandbox/sloshing_tank.py`.

## 4. Acceptance and interpretation

| Criterion | Required result |
| --- | --- |
| Production state provenance | Immutable rc3 image label; current diagnostic Git SHA |
| Numerical integration | Exactly 5 accepted fixed CN steps, 0 rejections/failures |
| Start state | First saved `z(t0)` equals accepted pseudo-state (scaled L2 <= 5e-8) |
| MPI numerical parity | Initial/post state and every saved complete time state within 5e-8 scaled L2 plus 1e-9 absolute floor |
| Algebraic row classification | Actual PETSc mass entries on classified rows <= 1e-12 |
| DAE algebraic residual | Measure before and after pseudo-step; only decide a physical acceptance threshold after observing actual scales |
| Prescribed q/p initial data | Absolute drift <= 1e-8 + 1e-8 times the original field norm (**physical gate**, not relaxed merely to make it pass) |
| Hamiltonian / signed boundary power | **Separate subsequent** non-vacuous CN midpoint test; not certified here |

**Key interpretation:** a possible result
`MPI full-state parity = PASS`, `algebraic residual reduced`,
`physical initial-data fidelity = FAIL` is not contradictory.
It would precisely isolate the scientific deficiency that must be
corrected **before** certifying physical initial conditions and
reusing subsequent energy-balance results as meaningful physical
evidence.

**Next decision after real data:** review the weak-form algebraic
residual and q/p drift together. Then design a physically controlled
projection that adjusts only algebraic variables while honoring the
differential data, with explicit solvability checks, and requalify
the full trajectories and midpoint power balance. Do not silently
replace the accepted pseudo-TS state in production until this
projection is proven for the targeted model and MPI.

## 5. Native results — 20 November 2026 (measured, not predicted)

**Provenance:** exact native Ubuntu image
`scrimp-mpi:v1.3.0-rc3` (OCI revision
`1d6d0c90115213298047b28237a7f38087563a35`),
diagnostic commit `b56ceedbbba17ae3d79baf36b2e4e12a78040e1c`.
Published unchanged `wave.py`, `h=0.2`, 2,230 DOFs, fixed
CN `dt=0.002`, five accepted steps, 1/2 MPI ranks,
distributed GetFEM→PETSc COO, compact mode, MUMPS on_change.

**User-owned host output** (never staged into the Docker build):

~~~text
sandbox/outputs/v130-dae-wave-h0p2-b56ceedb-7TyrnhWS/
~~~

| Quantity (Euclidean coefficient norm) | 1 rank | 2 ranks |
| --- | ---: | ---: |
| Prescribed q initial L2 norm | 0 | 0 |
| Accepted pseudo q norm / q drift | 0.00114424447 | 0.00114424447 |
| Prescribed p norm | 3.02730717 | 3.02730717 |
| Accepted pseudo p norm | 3.02717936 | 3.02717936 |
| p initial-data L2 drift | 0.00751282216 | 0.00751282216 |
| Algebraic residual *before* pseudo-TS | 0.010689895 | 0.010689895 |
| Algebraic residual *after* pseudo-TS | 2.44609825e-18 | 2.31888379e-18 |

The maximum scaled full-state discrepancy at **every** saved
Crank–Nicolson time point is **1.78851051e-14** for MPI 1 vs 2.
Thus the MPI parity and post-pseudo algebraic residual tests
**PASS**, while physical fidelity of the independently prescribed
initial q/p **FAILS**. The printed `1.144e11` relative q change
divides by a `1e-14` floor because q originally vanishes: it
is **not** a physically meaningful relative percentage. The
relative coefficient-norm p drift is about 0.00248 (0.248%).

The original left port is an enforced flow input
`U_L(t,y)=0.1 sin(4t) sin(4 pi y) exp(-10(2.5-t)^2)`; hence
`U_L(0,y)=0`. The prescribed initial momentum is a Gaussian
centered at x=0.5, with a **nonzero mathematical trace at x=0**.
The constitutive constraint `e_p=p/rho` and boundary constraint
`e_p|Gamma_Left=U_L` could therefore be **incompatible**
when both are imposed at t=0. Discrete weak-form compatibility
is **not yet proven** (nor disproven) by this continuous trace
observation; a targeted algebraic-only projection solvability
experiment is required before proposing to preserve q/p exactly.

**Next scientific gate:** hold every prescribed dynamic coefficient
q and p fixed, solve *only* the actual discrete algebraic constraint
rows for all effort/control unknowns, evaluate the independent
production `IFunction(t0,z,zd=0)` residual after the solve, and
record any singularity, conditioning and nonzero minimum residual.
If that projection is unsolvable, redesign an **admissible initial
momentum**, not a silently altered pseudo-TS. If it is solvable,
test exact dynamic-data preservation and subsequent full
Crank–Nicolson power balances at 1/2 MPI ranks before any
production integration. Never equate mathematical field traces
with the actual tested discrete weak constraints without checking.

