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


## 6. Follow-up: algebraic-only feasibility with frozen differential initial data

**First native probe completed (20 November 2026); further A/B
diagnostics staged but not yet executed.** The first measured wave test
proves an important trade-off (algebraic residual ~1e-18 versus
changed q/p), but does not tell us whether the *discrete algebraic
system* can be solved with the originally prescribed q/p held fixed.
The nonzero pointwise left-boundary trace of the Gaussian gives a
**potential incompatibility**, not a conclusive discrete result.

The next experiment is **diagnostic only**, with identical published
`wave.py` equations, original `rc3` Docker image, h=0.2,
1/2 MPI ranks, and unchanged physical Crank–Nicolson solver. It does
not alter `DPHS.init_step`, the default initialization policy,
physical CN start state, sources or MUMPS configuration.

- `sandbox/probe_v130_wave_fixed_dynamics.py`: intercepts the
  initial step only to capture prescribed and accepted states;
  evaluates the actual distributed `IFunction(t0,z,zd=0)` on
  algebraic residual rows with prescribed dynamic coefficients
  held fixed; obtains the same rows/columns of the **assembled
  PETSc stiffness** as `K_aa`; gathers owned rows once; and
  forms `K_aa delta_a=-F_a`. All dynamic coefficients are left
  **bitwise unchanged**.
- `sandbox/solve_v130_dae_algebraic_block.py`: a pure-SciPy LSMR
  attempt on the (possibly rank-deficient) algebraic block, with
  explicit iteration, norm, condition and nonzero-minimum-residual
  diagnostics. A large residual with poor conditioning is **not**
  called a mathematical proof of incompatibility.
- The candidate is checked against the **original full production
  `IFunction`**, not just the extracted linear equation, to ensure
  the discrete Jacobian truly represents the algebraic rows.
  The source also measures the original accepted pseudo-TS residual
  and restores that state before the ordinary main physical solve.
- `sandbox/compare_v130_wave_fixed_dynamics.py` verifies
  1/2-rank agreement in the feasibility conclusion and actual
  residual, while requiring **zero change** in the originally
  prescribed differential coefficients and rejecting a failed
  Jacobian prediction.

Launch on the tested Ubuntu workstation, using its existing rc3
image and its existing no-sudo Docker permission:

~~~bash
cd ~/Simulations/scrimp
git switch release/1.3.0-preparation
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-fixed-dynamics-probe.sh
~~~

The runner verifies the immutable **rc3 image source label**;
stages no extra image layers; bind-mounts only the explicitly
named diagnostic Python files into the *runtime* container; writes
the reports under an isolated new directory in
`sandbox/outputs/v130-algebraic-only-h0p2-.../`.
The host's modified `sandbox/sloshing_tank.py` is not touched.

**Interpretation before changing SCRIMP:**

| Observation | Permitted inference | Next step |
| --- | --- | --- |
| `native F_a` near zero with q/p unchanged, validated `K_aa` | Algebraic-only projection is numerically feasible for **this** discrete wave | Verify uniqueness/conditioning, compare CN trajectories and port power before considering an opt-in initializer |
| `native F_a` remains nonzero with LSMR normal equations converged and acceptable conditioning | Strong evidence that **these** initial data do not satisfy all discretized algebraic constraints with q/p frozen | Inspect constitutive and boundary rows, identify admissible p0 and test with correct physical boundary traces |
| LSMR iteration/conditioning issues or mismatch between `K_aa` and native `IFunction` | Inconclusive; neither compatibility nor incompatibility established | Diagnose rank/nullspace or Jacobian/model assembly before deciding |

In all cases, an ordinary after-pseudo-TS CN solve must still
converge, and **no production code is patched on the basis of
these unexecuted diagnostics**. Source-only CI and synthetic tests
prove the probe's shape/decision semantics, not the actual
native wave result.

## 7. Measured fixed-dynamics feasibility and A/B follow-up (20 November 2026)

### 7.1 Completed native results: original published control

The actual native probe of the published wave at `h=0.2`,
2,230 DoFs (1,175 algebraic), 5 accepted CN steps, MPI 1/2,
using unchanged qualified `scrimp-mpi:v1.3.0-rc3` with
source revision `1d6d0c90115213298047b28237a7f38087563a35`,
and diagnostic commit
`de3b2c734d03f94bc11db3613a3c0daf259f2105`:

| Euclidean coefficient residual | 1 MPI rank | 2 MPI ranks |
| --- | ---: | ---: |
| Algebraic before correction | 0.010689895 | 0.010689895 |
| After solving K_aa delta_a=-F_a, q/p EXACTLY frozen | 4.82244319e-6 | 4.82244319e-6 |
| After standard pseudo-TS (q/p altered) | 2.44609825e-18 | 2.21270009e-18 |
| Difference actual residual vs K_aa prediction | 1.9331465e-18 | 1.8659101e-18 |
| Iterative LSMR condition estimate | 167.7 | 126.4 |

Actual host evidence:
`sandbox/outputs/v130-algebraic-only-h0p2-de3b2c73-RSTntEzH/`.

The algebraic-only correction reduced the residual by about
a factor 2,217, but the remaining residual was
`4.82e-6`, about `4.51e-4` of the original.
`FIXED-DYNAMICS ALGEBRAIC FEASIBILITY: UNRESOLVED`
at **both** MPI ranks. The reconstructed PETSc Jacobian accurately
predicts the actual GetFEM production `IFunction` to ~1e-18.
Thus, **this is a genuine remaining residual in the discrete
constraint equations**, not an assembly/Jacobian mismatch.
However, a *nonzero iterative least-squares residual does
not alone prove that no different, potentially larger algebraic
correction exists*; the estimator `condA` is not a formal rank
certificate.

### 7.2 Next experiment, source-ready; native execution pending

The enhanced diagnostic constructs a **computed left-nullspace
witness** `w = residual_after / ||residual_after||`, reports
`||K_aa^H w||` and `|w^H F_a|`, and evaluates the explicit bound

```text
||K_aa delta_a + F_a|| >=
  max(0, |w^H F_a| - ||K_aa^H w|| * R)
  for any ||delta_a|| <= R.
```

The bound is for a **specified finite correction radius**,
default `R=10*max(1, ||delta_a_LSMR||)`. It does *not*
establish infeasibility for arbitrarily large corrections,
nor an exact rank/nullspace theorem.

The diagnostic now also records remaining algebraic residuals
by **named field** (`e_q`, `e_p`, boundary input/output ports),
rather than relying on a global norm alone.

Finally, it compares **two different boundary-control problems**
on the exact same mesh and with **identical q0/p0**:

- **published:** original left velocity `U_L(0,y)=0`;
- **matched_left_control:** replace only the prescribed left
  velocity `U_L(t,y)` by `p0(0,y)/rho(0,y)
  = 3**(-20*(0.25+(y-0.5)^2))/3`, to ensure **pointwise
  initial left-boundary compatibility**.

The comparison does **not** modify `examples/wave.py`,
the FE spaces, mass/stiffness bricks, `q_0` or `p_0`;
the alternative is assembled with a **temporary diagnostic
interception of `DPHS.set_control`**. We keep the original
production image/source and its actual pseudo-TS/main CN
state untouched. Thus differences between scenarios test
the specific `p_0/rho` vs `U_L(0)` compatibility hypothesis;
they are not described as numerical equivalence of two
identical physical problems.

Run the **same command** after pulling the current branch:

~~~bash
cd ~/Simulations/scrimp
git switch release/1.3.0-preparation
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-fixed-dynamics-probe.sh
~~~

This updated runner performs four small native cases:
published/matched control × MPI 1/2, compares each scenario
across MPI and finally compares the two controls while
checking q0/p0 are unchanged. It archives independent JSON
and numerical reports in a new unique directory under
`sandbox/outputs/`, never in the Docker build context.

**Interpretation:**
a near-zero algebraic-only residual for `matched_left_control`
but a nonzero bounded dual obstruction for `published`
would strongly localize the defect to the initial boundary
compatibility of this wave example. If both problems remain
incompatible, investigate other boundary/constitutive
rows and possible singularity. In neither case should the
production pseudo-TS, initial state or solver be silently
changed from these diagnostic observations alone.


## 8. Follow-up A/B run and GetFEM parser fix (20 November 2026)

The developer ran the A/B diagnostic at immutable rc3 runtime revision
`1d6d0c90115213298047b28237a7f38087563a35`, diagnostic
source `7fd9c728495d98eb242778e214516519475a622c`.
The **published-control** 1/2-rank runs and their comparison **passed**.
The alternative `matched_left_control` run failed during the
GetFEM source-brick **parsing phase**, before any physical solver or
projection was executed. Its previous mathematical compatibility
status is therefore **NOT TESTED**, not FAIL.

Host evidence:

```text
sandbox/outputs/v130-algebraic-only-h0p2-7fd9c728-n2jTOuos/
```

### 8.1 What the completed published-control run proves

The residual after the fixed-q/p algebraic projection is
`4.82244319e-6` on both 1 and 2 ranks. The largest contribution
is specifically on the **`e_p` constitutive field**:

| Algebraic residual after algebraic-only projection | MPI 1 (coefficient L2) |
| --- | ---: |
| `e_p` | 4.82044425e-6 |
| `U_L` | 9.81721862e-8 |
| `Y_L` | 9.81721861e-8 |
| `e_q` and other controls | 0 to ~1e-17 |

LSMR stopped with code 2 on iterations 953/954 for MPI 1/2.
The recomputed normal-equation residuals are
`2.503872e-18` and `2.2877677e-18`. A computed
bounded-radius dual witness gives `4.82243693e-6` and
`4.82243747e-6` as lower bounds on algebraic residual norm
when `||delta_a|| <= 12.09304`. The witness adjoint norms are
`5.1921228e-13` and `4.744002e-13`. These are convincing
**numerical bounded-correction obstruction measurements**,
not exact rank-nullspace proofs and not infeasibility for all
unbounded corrections.

The pseudo-TS still reaches ~1e-18 algebraic residual **by
changing prescribed q/p**, so it is not a valid replacement for
the fixed-data compatibility investigation.

### 8.2 Harness-only GetFEM syntax error and correction

The diagnostic override had generated

```text
(3**(-20*(0.25+(y-0.5)*(y-0.5)))/3.)*Test_U_L
   ^
Unexpected token.
```

`3**...` was a **Python-style exponent**, whereas
`DPHS.set_control` passes the source to **GetFEM generic
assembly**; the parser rejected `**`. This failure occurred
during `add_source_term_brick`, not during initial consistency,
MUMPS, Jacobian extraction, CN, or physics. PETSc option warnings
following this exception are cleanup artifacts, not independently
qualified numerical regressions.

The diagnostic expression is now the equivalent
double-precision analytic trace, expressed in functions/operators
used by GetFEM generic assembly:

```text
exp(-21.972245773362197*(0.25+(y-0.5)*(y-0.5)))/3.
```

Here `21.972245773362197` is `20*log(3)` rounded to
double-precision. It represents `p0(0,y)/rho(0,y)`.
It does **not guarantee exact equality of the weak discrete
boundary trace** after FE interpolation, which remains what
the native A/B experiment must actually decide.

A regression test now checks that the diagnostic expression
contains no `**` and matches the intended Gaussian trace
at representative boundary points. **The test cannot replace
GetFEM parsing.** The native runner therefore first constructs
the **matched-control GetFEM model only**, validates its
source-brick grammar, and stops immediately with its
`control_preflight.log` if the parser rejects it. Only after
this preflight passes does it execute the four small A/B
MPI=1/2 physics/projection cases.

Repeat the **existing runner** after pulling the branch; no Docker
candidate rebuild, GitHub release, or production-code change:

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-fixed-dynamics-probe.sh
~~~

**A/B native rerun now completed:** see section 9.
The GetFEM preflight PASS and both MPI variants of each
scenario completed. The analytic trace-matched control
**increased**, rather than eliminated, the fixed-dynamic
algebraic residual. Neither configuration passed the
prescribed-data algebraic feasibility threshold.

## 9. Completed full native A/B experiment (20 November 2026)

**Observed native output:** immutable
`scrimp-mpi:v1.3.0-rc3` (source
`1d6d0c90115213298047b28237a7f38087563a35`),
diagnostic commit
`25f9bf27d260abca401899c0c8b2e1bd2c97300e`,
`h=0.2`, 2,230 DOFs, 1,175 algebraic, five accepted
Crank–Nicolson steps per case, identical prescribed Gaussian
`p0` and zero `q0`, MPI 1 and 2.

Original host output:

```text
sandbox/outputs/v130-algebraic-only-h0p2-25f9bf27-DbVTy7aw/
```

The **GetFEM matched-control source preflight passed**.
Both `published` and `matched_left_control` then completed
on 1/2 ranks; linearization verified against original PETSc
`IFunction` to approximately `1.9e-18`. All 1/2-rank
comparisons reported PASS with the same physical q/p
coefficients and equal residual norms.

| Euclidean coefficient residual | Published U_L(0)=0 | Analytic trace-matched U_L(0)=p0(0,y)/rho(0,y) |
| --- | ---: | ---: |
| Initial all algebraic rows | 0.0106898950 | 0.0106921438 |
| Fixed q/p, optimized algebraic-only | **4.82244319e-6** | **9.18141097e-6** |
| After original pseudo-TS (q/p changed) | ~2e-18 | ~2e-18 |
| Constitutive `e_p` residual after projection | 4.82044425e-6 | 9.17758251e-6 |
| Left input `U_L` residual after projection | 9.81721862e-8 | 1.87465454e-7 |
| Left output `Y_L` residual after projection | 9.81721861e-8 | 1.87465454e-7 |
| Computed dual lower bound (`||delta_a|| <= ~12.09`) | ~4.822437e-6 | ~9.181407e-6 |

The trace-matched control made the fixed-dynamic algebraic
residual **approximately 90.4% larger**, not smaller.
The numerical normal equation residuals are `~2e-18`
(published) and `~3e-18` (trace-matched), with bounded
adjoint-witness residuals near the measured least-squares
minimum. The original hypothesis that simply matching the
**analytic** trace `p0/rho` at the left boundary is a
sufficient discrete fix is **not supported by these data**.
This does **not** prove that the boundary constraint is
irrelevant or that no compatible control exists: the
constitutive relation is a **weak finite-element projection**,
and pointwise continuous traces do not generally commute
with that projection.

**Scientific next step:** manufacture a truly FE-representable
and left-boundary-compatible nonzero momentum for this
very same model:

```text
q0(x,y) = [0, 0]
p0(x,y) = x*(3-x)
rho(x,y) = 3-x
e_p(x,y) = x
U_L(0,y) = 0
```

Both `p0` (degree 2) and `e_p` (degree 1) are exactly
representable by the published CG2 space; `e_p=0` at
the left boundary. The weak constitutive relation
`-M_e_p e_p + CR_p p = 0` and left velocity constraint
have an explicit compatible candidate. The q/e_q
constitutive pair is initially zero. This is a
**different scientific initial-value problem** from the
published Gaussian, deliberately chosen as a positive
control for the algebraic-only projection algorithm.
In particular, even for this compatible initial state,
the historical pseudo-TS may still move differential
variables by integrating in pseudo-time; that movement
must be measured and reported separately.

A failure of the manufactured case, despite the exact
FE construction, would be evidence to investigate the
constraint extraction, the LSMR solve, or other algebraic
rows before proposing changes to the Gaussian data.
A PASS would isolate the obstruction to the original
initial data/control combination (not a universal proof
of the theoretical DAE constraint structure).

No change to `scrimp/dphs.py`, `examples/wave.py`,
the qualified Docker image, or user local files is
authorized by this A/B result alone.

## 10. Native exactly FE-compatible positive-control probe — ready, pending

The completed A/B evidence in section 9 rejects the simple
hypothesis that prescribing the **analytic** trace of the
Gaussian's p/rho at the left edge is sufficient to remove the
**discrete** algebraic obstruction. The smaller residual in the
published control is **not** proof that the original Gaussian
is an admissible DAE initial state. A different mathematical
experiment is needed to discriminate inconsistent data from a
fault in the extraction/projection algorithm.

The new source-only manufactured input to the otherwise
*unchanged published wave equations* is:

```text
q0 = [0., 0.]
p0 = x*(3-x)          (CG2 polynomial; nonzero)
rho = 3-x
e_p = p0/rho = x     (CG2 polynomial)
U_L(0,y) = 0         (unchanged published control)
```

This defines an explicit feasible **continuous and FE** candidate
for the constitutive and left-boundary algebraic constraints:
`p0` has degree two, `e_p=x` degree one, and
`e_p|Gamma_Left = 0`. The other three ports
are the unchanged force-driven ports. Unlike the
Gaussian A/B comparison, this deliberately changes
the initial data `p0`, but *only inside a runtime diagnostic
monkeypatch while constructing the same model*. It does
**not** change `examples/wave.py`, `dphs.py`, or
the qualified `scrimp-mpi:v1.3.0-rc3` image.

This test is intentionally a **positive control**, not
a prescription to replace the published Gaussian model.

The existing `sandbox/probe_v130_wave_fixed_dynamics.py`
now supports `--scenario polynomial_compatible` and
reports separately:

- actual algebraic-only feasibility and residual with
  prescribed q/p coefficients *exactly frozen*;
- named weak-form constitutive and boundary residuals;
- exact GetFEM `IFunction` vs extracted PETSc `K_aa`
  algebraic linearization;
- the **drift of prescribed dynamic variables in the
  historical pseudo-TS**, even when the original dynamic
  data are already algebraically consistent;
- MPI 1/2 parity and successful unchanged main CN solve.

The new runner executes **only 1 and 2 rank runs of this
polynomial positive control**, not the previously validated
four original/matched-control runs:

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-exact-fe-initial-probe.sh
~~~

It first builds the model with prescribed polynomial data
as a native GetFEM construction preflight, then runs two
complete small fixed-step CN simulations with the
unchanged production pseudo initializer. The feasibility
projection is **measured but never substituted** into the
production physical time integrator. It writes a unique
user-owned output directory under
`sandbox/outputs/v130-dae-fe-positive-h0p2-.../`.
Existing Docker tags, native layers, user workspace edits,
and the `sandbox/outputs` history remain unchanged.

**Interpretation depends on actual native evidence:**

| Positive control observation | Consequence |
| --- | --- |
| FE polynomial fixed-data algebraic residual <= 1e-10 + 1e-8 × initial residual | Strong evidence that the original Gaussian and boundary input, not general K_aa extraction/LSMR, cause its nonzero least-squares obstruction |
| FE polynomial residual above threshold | Investigate residuals by field, mesh interpolation and possible source/constraint implementation issues; do not assume incompatibility of Gaussian was proven |
| Pseudo-TS changes q/p even for exactly compatible data | The default pseudo-TS is not a pure algebraic projection; the accepted physical start must be reconsidered, independently of MPI parity |

No expected PASS is represented as measured until this
**new native experiment runs**.
