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

## 10. Native exactly FE-compatible positive-control probe — completed

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

**Native execution completed and positive-control PASS observed
at MPI 1 and 2.** See section 11 for measured evidence and interpretation.

## 11. Completed nonzero FE-compatible positive control — native evidence

**Run date:** 20 November 2026. Image:
`scrimp-mpi:v1.3.0-rc3`, immutable solver source
`1d6d0c90115213298047b28237a7f38087563a35`.
Diagnostic Git revision:
`1dfbc184bf82c3b859a64a95f9bbf950bb7f5af4`.
2,230 global DOFs, 1,175 algebraic, unchanged 2D-wave
bricks/FE/boundaries/MUMPS, `h=0.2`, CN `dt=0.002`,
five steps, MPI 1 and 2.
`p0=x*(3-x)`, `q0=0`, `U_L(0)=0`, `e_p=x`.

**Original host evidence (not committed/generated again):**

```text
sandbox/outputs/v130-dae-fe-positive-h0p2-1dfbc184-MEYxdCWf/
```

| Measured coefficient Euclidean L2 | 1 rank | 2 ranks |
| --- | ---: | ---: |
| Initial algebraic residual | 0.162506819 | 0.162506819 |
| After algebraic-only LSMR, q/p exactly frozen | **7.06518655e-12** | **7.13966896e-12** |
| After historical pseudo-TS | 1.88022801e-16 | 1.79232771e-16 |
| Original `IFunction` vs extracted linear predictor | 2.1793999e-17 | 2.2685362e-17 |
| q drift **caused by pseudo-TS**, despite FE-compatible data | **0.00195958812** | **0.00195958812** |
| p drift **caused by pseudo-TS**, despite FE-compatible data | **1.30196443e-5** | **1.30196443e-5** |

The fixed-dynamics projection **passes** the source-defined
algebraic tolerance at MPI 1 and MPI 2. The `e_p` row residual
after projection is ~`7.06e-12`; the other algebraic
fields remain at ~`1e-14` or less. The LSMR
correction is therefore not generically incapable of
solving SCRIMP's true weak constraints. This supports
diagnosing the Gaussian/left-boundary/control combination
**as the source of the earlier projection obstruction**,
but does not establish its exact mathematical nullspace or
rule out other numerical details for other DAE models.

The **separate pseudo-TS branch** still changes q and p
even when the same prescribed differential data are
known to admit an algebraically consistent state with
no such changes. Its role is time evolution, not
a projection that fixes differential initial data.
This is a reproducible physical initialization defect
relative to the stated requirement that explicitly
prescribed q/p are preserved. It does **not**
establish that PETSc or the actual full CN integrator
has a defect.

**Next independent gate:** from exactly the same
manufactured wave initial data, compare an unchanged
pseudo-TS initialization versus a **diagnostic-only
algebraic projection published as the main physical
start state**. Verify the full CN trajectories, 1/2-rank
parity and, independently, the **Hamiltonian plus signed
boundary power at every CN midpoint** with a driven
(non-vacuous) left port. Confirm the original physical
initial q/p are exactly preserved by the projected
branch. Do not enable a production option or change
`scrimp/dphs.py` until this gate passes; large
performance runs and Maxwell are not needed here.

## 12. Full CN initializer comparison with non-vacuous port power — source ready

**Status (20 November 2026): four-case native A/B COMPLETED.
Full CN trajectory MPI parity, signed non-vacuous Hamiltonian/port
power and projected prescribed-data fidelity all PASS. Historical
pseudo-TS fidelity FAIL remains reproducible.** See section 13
for complete measured native results.

### Controlled same-model A/B

The same published 2D-wave bricks, FE spaces and boundary geometry,
h=0.2, q0=[0,0], p0=x*(3-x), rho=3-x, are used in both modes,
with a genuinely driven left port U_L(t)=0.08*sin(20*t) (zero at t0).
Main TS: fixed Crank–Nicolson dt=0.002, t_f=0.02, ten steps,
MPI 1 and 2.

- Mode **pseudo** uses the unchanged production pseudo-TS accepted
  state as the first state of the main physical CN solve.
- Mode **algebraic** first runs that same original pseudo-TS, then
  uses a *diagnostic-only* real PETSc algebraic-block projection
  of the original prescribed dynamic state. The accepted projection
  replaces the initial GetFEM state **only before the main TS**
  is constructed. The actual original IFunction algebraic rows
  must pass tolerance, with all q/p coefficients bitwise preserved.
  No production initializer or global solver behavior is changed.

### Independent scientific checks

1. The actual first saved physical state must match the state
   published by each initializer. Algebraic mode preserves
   every prescribed dynamic coefficient (L2 <=1e-12).
2. The entire **11-step-including-t0 global CN trajectory**
   and time series of Hamiltonian, per-boundary midpoint power,
   accumulated port work and per-step defect must match
   between MPI 1 and MPI 2 **within each mode**.
   Different modes have different initial states, so their
   physical trajectories should not be required to coincide.
3. The independent existing native acoustic energy audit
   evaluates the same potential/kinetic terms as wave.py and
   the signed port-Hamiltonian midpoint identity:
   H[n]-H[n-1] + dt*sum_ports(flow*effort)[midpoint] = 0.
   The drive must yield nontrivial absolute boundary work
   exceeding 1e-9. The per-mode balance acceptance is
   max |H-H0+W| <= 1e-9 + 1e-5*balance_reference_scale.

A violation of one criterion cannot be hidden behind another
PASS: e.g., MPI parity does not establish correct initial data,
and conservation from a wrong t0 does not repair that data.

Run the small native experiment without rebuilding rc3:

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-wave-initializer-physics.sh
~~~

The runner executes only **pseudo/algebraic × MPI 1/2**,
10 CN steps per case. It stores unique host-side output in
sandbox/outputs/v130-wave-init-power-h0p2-.../ including
JSON summaries, complete .npz trajectories, four native logs,
the comparison and a manifest with the immutable source SHA
and current diagnostic Git revision.

Source-only regression tests cover bad later-time MPI parity,
an altered projected q0, and failure of an otherwise non-vacuous
signed power balance. The native program always preserves the
four per-case scientific records when a power check fails;
only the final strict comparison marks the overall experiment
as FAIL. There is no v1.2.0 rerun, large scaling benchmark,
Docker rebuild, GitHub release or change to the production
wave.py, dphs.py, scrimp-mpi:local or sandbox/sloshing_tank.py.

**Native evidence received and interpreted:** see section 13.
The previous pending status is superseded by measured results.

## 13. Full CN physical gate PASSED on two initializers — 20 November 2026

This **completed native experiment** validates (for the
specified manufactured, linear, FE-compatible acoustic model)
the separation between *physically preserving initial data*
and *preserving the port-Hamiltonian power balance*.

**Provenance:** unchanged qualified image
`scrimp-mpi:v1.3.0-rc3` with OCI source
`1d6d0c90115213298047b28237a7f38087563a35`,
diagnostic Git revision
`62e31d9b5a8f65637b48e791a708bd961567cc93`.
`p0=x*(3-x)`, `q0=0`,
`U_L(t)=0.08*sin(20*t)`,
`h=0.2`, `dt=0.002`, `t_f=0.02`,
**10 accepted fixed CN steps**, 11 complete saved states
including t0, 2,230 total DOFs, 1/2 MPI ranks.

Native host evidence, retained by the developer:

```text
sandbox/outputs/v130-wave-init-power-h0p2-62e31d9b-zDa95p5J/
```

| Measured criterion (Euclidean coefficient norms where relevant) | Historical pseudo-TS | Diagnostic algebraic-only initializer |
| --- | ---: | ---: |
| Prescribed q0 drift | **0.00195958812** | **0** |
| Prescribed p0 drift | **1.30196443e-5** | **0** |
| Max cumulative `|H-H0+W|`, MPI 1 | **4.19642484e-15** | **5.6085393e-15** |
| Max cumulative `|H-H0+W|`, MPI 2 | **3.75235596e-15** | **6.05263063e-15** |
| Port integrated work, MPI 1 (signed) | **-5.41131226e-5** | **-5.39545892e-5** |
| Absolute port work (non-vacuous), MPI 1 | **5.41131226e-5** | **5.39545892e-5** |
| Complete CN 1/2-rank trajectory **max scaled L2** | **1.75624262e-14** | **2.05308545e-12** |
| Original prescribed initial-data fidelity | **FAIL** | **PASS** |
| Non-vacuous signed midpoint Hamiltonian/port-power balance | **PASS** | **PASS** |
| Full trajectory MPI parity | **PASS** | **PASS** |

The power gate's previously specified threshold was
`1e-9+1e-5*balance_reference_scale`, about
`2.00015e-5`. Although this *tolerance* is intentionally
conservative, the actual native numerical defects were
approximately `4e-15` to `6e-15`, **ten orders of
magnitude below that threshold**. Port work is nonzero,
so this is not a vacuous zero-power test.

**Crucial scientific distinction:** the two main physical
Cauchy problems differ only through their initializer-
produced state at t0; consequently, their integrated
boundary works need not be equal. Correct port-Hamiltonian
power conservation from an unwanted initial state does
not make the prescribed input data faithful. The
algebraic-only branch verifies both relevant properties
on this **one** linear acoustic problem.

The tested algebraic branch currently uses a **diagnostic
centralized SciPy LSMR projection** of a gathered K_aa
block on rank 0 and, to preserve direct comparison, runs
the old pseudo-TS first before replacing its accepted
state. **This is NOT production-ready**: it does not
avoid the obsolete pseudo-TS or scale as distributed
MPI data grow. Do not promote the experimental code
unchanged to SCRIMP's core. The next engineering
gate is an independently validated PETSc-distributed
algebraic subproblem implementation, with no fake
PASS for incompatible Gaussian initial data.

## 14. Production DAE initializer design constraints (not yet implemented)

The next candidate must be explicit opt-in and preserve
the current default; at this stage the qualified API
behavior is only linear algebraic projection on models
with verified row/unknown classification.

**Algorithmic contract:**

1. Accept prescribed dynamic flow/state DOFs **bitwise
   unchanged**; determine algebraic unknowns using the
   actual GetFEM variable intervals / port substitution
   flags; check algebraic rows have zero derivative mass
   support in the assembled PETSc mass. Reject incomplete
   layouts, unsupported substituted/lumped/nonlinear
   models, rather than implicitly treating those as
   standard index-1 constraints.
2. At the exact physical t0, evaluate the **native**
   residual F(t0,z,zd=0) and its algebraic rows, using
   the existing assembled PETSc stiffness, source and
   GetFEM/PETSc data. The projected linearized system
   has the form `K_aa delta_a=-F_a`.
3. Avoid gathering global K_aa on one rank. Use a
   genuinely **MPI-distributed PETSc** algebraic
   submatrix and vector with ownership-defined row
   and column index sets; support LSQR or a justified
   rank-compatible linear solve, with a tested PETSc
   options prefix and explicit iteration/conditioning
   limits. Do not assume K_aa is nonsingular or
   that `KSPPREONLY+LU` solves least squares on a
   rank-deficient saddle-point block.
4. Re-evaluate the **original full IFunction** on
   the candidate. Require the actual algebraic
   residual <= an explicitly defined absolute/relative
   scale; verify every prescribed dynamic coefficient
   is unchanged. Reject incompatible initial data
   (as in the published Gaussian), **do not silently
   modify q0/p0 or fall back to pseudo-TS**.
5. Publish the candidate to GetFEM on all ranks
   with existing MPI scatter machinery, at t0.
   This occurs **instead of** running the pseudo-TS,
   so it really removes pseudo-time work, temporary
   TS/KSP/factorization and initial physical drift.
   Preserve the existing main CN/BDF TS solver, its
   options, factorization lifecycle and monitor,
   and the default `init_step_ts_type=pseudo`.
6. Qualify nonzero compatible wave initial data with
   full CN trajectories and signed power (1/2 MPI
   ranks) against the established diagnostic trace.
   Prove failure is explicit and atomic on the
   incompatible Gaussian. Record memory and elapsed
   cost separately from the old pseudo-TS; source-only
   tests never substitute for native GetFEM/PETSc.

**Beyond this qualified linear scope:** a nonlinear
DAE requires a different constrained nonlinear
Newton/projection strategy; higher-index or implicit
boundary couplings may require compatible derivatives.
Do not advertise this linear procedure as a generic
high-index/nonlinear initializer, and do not change
SCRIMP's default time scheme without separate evidence.

## 15. Next native gate: distributed PETSc algebraic projector (pending)

The full CN A/B power gate **PASSED** (section 13), but its
source-only diagnostic projection used rank-0 SciPy LSMR
on a gathered algebraic block and executed the old
pseudo-TS first. This validated physics, **not**
production-grade scalability or initialization lifecycle.

A new **diagnostic-only** candidate is staged:

- `sandbox/probe_v130_petsc_mpi_initial.py`: uses
  distributed PETSc `IS` index sets for algebraic
  DOFs and `Mat.createSubMatrix` on the **actual
  assembled PETSc stiffness**; solves
  `K_aa delta_a=-F_a` with PETSc `KSP LSQR`,
  `PC NONE`, and the `scrimp_dae_trial_`
  options prefix. No complete K_aa is gathered
  on rank 0. Linear mass and stiffness,
  complete disjoint flow/effort intervals, no
  lumped or substituted variables, and zero
  algebraic mass rows are required.
- The full production `IFunction(t0,z,0)` is
  recomputed after the correction, the
  exact prescribed differential coefficients
  are checked **bitwise**, and the model is
  synchronized via existing global PETSc Scatter.
  A KSP reason indicating least-squares stopping
  does **not** substitute for a true DAE residual PASS.
- **No initial pseudo-TS is executed at all**
  in this diagnostic: a temporary method
  intercepts `DPHS.init_step` before the
  production main CN TS. The actual immutable
  rc3 `dphs.py` and solver are unchanged.
- Positive control:
  `p0=x*(3-x)`, `q0=0`, driven left
  `U_L(t)=0.08*sin(20*t)`; require
  native algebraic feasibility, complete 10-step
  CN time series, MPI 1/2 parity, nonzero
  power and Hamiltonian+signed work balance.
  Cross-check the measured signed work against
  section 13's independent diagnostic.
- Negative control: the *published Gaussian*
  with original `U_L(0)=0`; the native
  algebraic residual must reject the
  incompatible-data test **before** any CN
  step. An expected residual rejection is not
  confused with a PETSc failure or claimed
  proof of global algebraic infeasibility.
- `sandbox/compare_v130_petsc_mpi_initial.py`
  checks both scenarios against MPI 1 and 2.
  Source-only tests inject power, late trajectory
  and incompatible-data acceptance failures.

**Next Ubuntu command** (same frozen rc3; runtime bind mounts
of diagnostic Python files only):

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-petsc-mpi-dae-initial.sh
~~~

The runner creates a new isolated host output under
`sandbox/outputs/v130-petsc-dae-mpi-h0p2-.../`
with four logs, a manifest, positive JSON/full
trajectories, negative JSON rejections and a
strict comparison. Previous full wave simulations
and the nearly 3M-DOF scaling tests are **not rerun**.

**STATUS: native PETSc submatrix and LSQR API paths
not yet exercised; do not claim numerical PASS.**
The next run may expose version/layout or convergence
issues, which must be resolved using real logs before
any production implementation. This is deliberately
a prototype, not a general linear/nonlinear DAE
initialization option in `scrimp/dphs.py`.

## 16. First native distributed PETSc probe: PetscInt width correction

**Actual native run received 20 November 2026, before patch:**
diagnostic Git commit
`62d92b6a53784623abd886f0b480974e66441e97`;
unchanged rc3 source
`1d6d0c90115213298047b28237a7f38087563a35`;
`h=0.2`, positive nonzero polynomial initial
`p0=x*(3-x)`, intended 10 fixed CN steps
and MPI 1/2.

Host evidence:

```text
sandbox/outputs/v130-petsc-dae-mpi-h0p2-62d92b6a-3urxt6Di/
```

**Observed stage:** the first positive MPI=1 run
constructs GetFEM, assembles the mass/stiffness
and advances past `PETSc.IS.createGeneral`,
`Mat.createSubMatrix`, and `KSP.solve`.
It fails **before evaluation of corrected residual or
main CN** when inserting the algebraic correction:

```text
z.setValues(local, correction.getArray(readonly=True).copy(),
            addv=PETSc.InsertMode.ADD_VALUES)
TypeError: Cannot cast array data from dtype('int64')
to dtype('int32') according to the rule 'safe'
```

The GetFEM index set preparation used NumPy
`int64`, whereas qualified rc3 PETSc was built
with 32-bit `PetscInt`. The `IS` constructor
already cast to `PETSc.IntType`, but the subsequent
`Vec.setValues` call did not.

**Targeted source correction:** check the original
owned global algebraic indices fit within
`np.iinfo(PETSc.IntType)`, cast to
`np.asarray(local,dtype=PETSc.IntType)`, and
use those in the additive PETSc vector insertion.
A source-only regression test exercises the 32-bit
conversion and guards the actual call site.

**Scientific status:** no qualified positive-state
projection, PETSc LSQR convergence, CN trajectory,
Hamiltonian/port power, or negative Gaussian
rejection was measured by this aborted first run.
The unused `-ts_adapt_type` and
`-ts_exact_final_time` options warning follows
the exception before the main TS was constructed;
it is not evidence of a main time-integrator error.

Next native command (same unchanged qualified rc3
runtime, no rebuild):

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-petsc-mpi-dae-initial.sh
~~~

A passing **source CI** alone cannot validate the
PETSc native correction or physical properties.

## 17. Second native distributed DAE projection — residual slightly above gate

**Native result received 20 November 2026:** diagnostic
`ce6e8edf99a03408390961a0906d6c1c6cd5c7b1`,
qualified immutable source
`1d6d0c90115213298047b28237a7f38087563a35`.
Run stopped on the first positive `polynomial_compatible`
MPI=1 scenario, before any physical CN step; other three
scenarios were **not executed**. Host evidence:

```text
sandbox/outputs/v130-petsc-dae-mpi-h0p2-ce6e8edf-yuonmmPp/
```

The previous `PETSc.IntType` correction now works, and
the distributed PETSc submatrix + LSQR path reaches
the actual GetFEM/PETSc algebraic acceptance check.
For `h=0.2`, `p0=x*(3-x)`, `q0=0`:

| Quantity | Native measurement |
| --- | ---: |
| Initial algebraic residual, Euclidean L2 | 0.162506819 |
| PETSc LSQR projected algebraic residual | **1.88801542e-9** |
| Unchanged acceptance threshold `1e-10+1e-8*initial` | **1.72506819e-9** |
| Excess above acceptance threshold | **9.45%** |
| Previously measured centralized SciPy positive residual | ~7e-12 |

Therefore, the **strict UNRESOLVED/FAIL is correct**;
the PETSc-distributed implementation has not yet
demonstrated an accurate enough projection, even
though this manufactured FE-compatible datum is
known numerically to admit one. Do not loosen the
physical residual tolerance or label PETSc LSQR
as validated on this evidence. The
`-ts_adapt_type`/`-ts_exact_final_time` warnings
occur because rejection prevents the main TS
construction, not because CN failed.

### Bounded distributed correction and failure evidence

The trial now performs **at most three total**
PETSc-distributed LSQR defect-correction passes,
using the same assembled distributed `K_aa`
operator and actual `IFunction` re-evaluation:

```text
r_k = F_a(t0, z_k, zdot=0)
K_aa * delta_k = -r_k                   (PETSc distributed LSQR)
z_(k+1)[algebraic] = z_k[algebraic] + delta_k
z_(k+1)[differential] = z_0[differential] (bitwise)
```

Each pass records PETSc's convergence reason,
iterations, reported residual, correction norm
and the **independent native algebraic residual**.
The procedure stops immediately on a true
residual PASS, PETSc divergence, or stagnation;
it never continues unboundedly or accepts
normal-equation convergence alone.
If a later pass worsens the true residual,
the best state is restored and re-evaluated.
The acceptance criterion **remains**
`1e-10+1e-8*||F_a(t0,z0,0)||`.

Importantly, a future positive-model rejection
now writes a structured JSON and prints its
entire LSQR history before exiting nonzero,
so the reason for failure can be inspected
without relying on a full exception trace.
The negative original Gaussian remains an
expected explicit t0 rejection only if the
actual residual exceeds threshold, without
starting physical CN.

**Status:** updated source and source-only regression
ready, but no post-refinement native results yet.
No `scrimp/dphs.py`, `examples/wave.py`,
Docker, MUMPS or main CN change.

Rerun the **same small native launcher** on rc3:

~~~bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-petsc-mpi-dae-initial.sh
~~~
