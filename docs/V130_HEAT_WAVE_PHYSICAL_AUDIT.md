# SCRIMP v1.3.0 — Heat–Wave 2D: physically meaningful example

**25 November 2026 — native Heat–Wave 2D qualification PASS:
MUMPS MPI1/MPI2 physical FE and CN power parity on four steps,
AND MUMPS MPI1 full 15-second energy-power identity on all
3,000 consecutive CN steps. All 3,001 original accepted FE states
were archived after the generic output-monitor correction.
Visual annular-gallery evidence is distinct: eigenmodes,
MPI2 long time, MPI-parallel VTK and SuperLU_DIST remain
unqualified.**

## Full 15-second MPI1/MPI2 first attempt — PHYSICS PASS, STRICT PARITY GATE FAIL (27 November 2026)

The native user-host full-length run on source `e6326198`
completed *both* MUMPS MPI2 and MPI1 original Gaussian
Heat–Wave solves, **3,000 accepted CN steps**, fixed
`dt=.005`, `t_f=15`, **3,001 original GetFEM
time states per rank**, **7,689 DOFs**. Source guards
passed and all eight published snapshots remained immutable.

**Original evidence (historical FAIL not rewritten):**
`v130-heat-wave-full-15s-MPI1-MPI2-e6326198-fhkvnmsR.zip`.
The native comparison status correctly remains
`FAIL_OR_UNQUALIFIED_FULL_HEAT_WAVE_MPI2`.
Both MPI2 ranks saved **bitwise identical**
full original FE coefficient and physical power traces.
Every native and independently recomputed
**per-rank physics gate passed** in all three runs:
positive thermal dissipation, nonnegative heat/wave
energy, canceling oriented interface power, homogeneous
outer work, individual CN subsystem and total energy
balance at every step, cumulative balance, and monotone
Hamiltonian.

| Original 15-second physical result | MPI2 (both ranks) | MPI1 |
| --- | ---: | ---: |
| H(0) | 0.7836421211876798 | 0.7836421211876778 |
| H(15) | 0.008268389052222184 | 0.008268389052222169 |
| Largest independent total CN step defect | 2.198e-15 | 4.018e-15 |
| Largest heat subsystem CN defect | 1.450e-15 | 1.381e-15 |
| Largest wave subsystem CN defect | 2.226e-15 | 3.979e-15 |
| Largest cumulative CN power defect | 2.670e-15 | 5.698e-15 |
| Largest signed interface power | 8.882e-16 | 3.109e-15 |

**Only the full space-time MPI parity checker rejected**
the evidence. It stopped on the final original
coefficient array `zf`:
`max|zf_MPI1-zf_MPI2|=2.758281700265819e-9`.
The original per-coefficient `np.allclose` gate uses
`atol=5e-10, rtol=2e-8`. Since near-zero FE
coefficients receive almost only the absolute
allowance, a small absolute difference can fail
even though full physical energies agree to
around `1.5e-17` at the endpoint.
**This is a possible explanation, not a
verified attribution.** The historical verifier
aborts on the first failing array, so its
summary does not establish the full space-time
maximum or any normalized FE error metric.

To avoid a scientifically unjustified tolerance
increase or re-running the already completed
native PDE, the source now includes the
**read-only offline diagnostic**

```bash
python tests/test_v130_heat_wave_full_mpi_diagnostic_source.py
python sandbox/diagnose_v130_heat_wave_full_mpi_parity.py \
  sandbox/outputs/v130-heat-wave-full-15s-MPI1-MPI2-e6326198-fhkvnmsR \
  --out /tmp/heat-wave-full-MPI-parity-diagnostic.json
```

It reads the **original** three native FE time histories,
verifies the exact MPI2 replica,
reports pointwise tolerance violations,
their count, time and DOF indices,
the maximum raw and scaled errors,
the normalized full-space-time and per-time
L2 differences, the twelve largest
final-state discrepancies and
independent mismatches for every
original GetFEM scalar/power trace.
The old evidence and original `FAIL`
are unchanged. The diagnostic cannot mark a
model **PASS**. First inspect these metrics,
then decide whether the existing *comparison
metric* is appropriate or whether a real
numerical MPI discrepancy needs investigation.

**Full MPI2 native parity remains PENDING.**
Do not conflate three separate correct
physical solves with a proven full FE parity.

## Next open gate — complete 15-second MPI1/MPI2 FE parity (27 November 2026, SOURCE-STAGED)

The reference Heat–Wave physics is **already qualified** for:
- Four real 0.001-second CN steps on 7,689 GetFEM DOFs,
  MPI1 versus MPI2 complete FE and heat/wave/interface
  physical parity.
- All 3,000 real 0.005-second CN steps through `t_f=15`,
  in MPI1, including 3,001 accepted FE states and
  the independent total and individual subsystem
  energy-power identities.
- Both 15-second authentic GetFEM visualization
  datasets (Gaussian and optional annular packet),
  without claiming modal purity.

**Unresolved independent question:** does the *entire*
15-second Gaussian coupled PDE retain every CN
energy/power identity, exact two-rank agreement, and
full FE trajectory parity when run with MUMPS MPI2?
The earlier four-step MPI check and 15-second MPI1
run are NOT a proof of this full-horizon comparison.

A **dedicated**, original-model-native experiment is
now source-staged, not yet run:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-heat-wave-full-mpi2-physics.sh
```

It deliberately runs the **unqualified MPI2 full
horizon first**. Only on MPI2 success does it
repeat the previously accepted MPI1 15-second case
as a fresh directly comparable baseline. This is
necessary to compare all 3,001 actual original
GetFEM state vectors (not merely their endpoints)
and **every** independently assembled heat, wave,
interface, boundary and cumulative power sample.

The probe is the already validated
`sandbox/probe_v130_heat_wave_mpi_extended.py`
with `--store-full-states`. The new checker
`scripts/verify-v130-heat-wave-full-mpi2.py`
retains the verified physical acceptance algorithm,
requires 3,001 consecutive GetFEM accepted-time
snapshots, all 7,689 coefficients at each time
and the original full native MUMPS physical gates
**on all three ranks**. MPI2 rank0/rank1 histories
must be bitwise identical, while MPI1 and MPI2
full histories must pass the established strict
`atol=5e-10, rtol=2e-8` FE/power comparison.
The complete original rank-local NPZ arrays,
logs and summary are preserved in a CRC-verified ZIP.
The large source histories (3 × 3,001 × 7,689
float64 FE coefficients) intentionally make this
an evidence-heavy experiment.

**PENDING NATIVE OUTPUT.** Do not claim full
15-second MPI2 parity or per-step power balance
before this new evidence has actually passed.
This scientific pilot deliberately suppresses
parallel ParaView and plotting output; no
parallel VTK, `SuperLU_DIST`, time/mesh-refinement
or whispering-gallery eigenmode qualification is
implied. The known MPI1 and short MPI1/2
qualification remain valid regardless of this
new experiment's result.

## Native qualification COMPLETE — MPI1/MPI2 parity and full 15-second CN power balance (25 November 2026)

**Canonical actual GetFEM/PETSc evidence:**
`v130-heat-wave-MPI2-longCN-97c7606b-7NBvjofW.zip`,
source commit
`97c7606bac14f7420140dd1fc2a1d88847301538`,
pinned image `scrimp-mpi:v1.3.0-rc3`.
Original ZIP SHA256:
`9aff9f1f6cf3708f8736f48a77731ad4d814c08fd0d54c915a32f8aa1ac516fe`.
The original `summary.json` contains
`PASS_HEAT_WAVE_MPI2_SHORT_AND_MPI1_LONG_CN_PHYSICAL`,
`full_long_time_power_balance_qualified=true`,
and `errors=[]`. All source tests passed, and
the eight tracked published examples remained immutable.

**Stage A: independently accepted MPI1/MPI2**
on the true 7,689-DOF canonical Heat–Wave
problem: four fixed CN steps at `dt=0.001`,
compatible initial Gaussian `p` after
physical projected-CG1 boundary clamp,
true default conforming interface compiler,
and GetFEM quadratures evaluated **on every rank**.
All original full state vectors and separate
heat/wave/interface power histories agree
between MPI1 and MPI2 within the prescribed
numerical tolerance.
The two MPI2 ranks saved **bitwise identical**
FE and power histories.

| MPI1/MPI2 check | Observed maximum |
|---|---:|
| Full FE coefficient difference | `2.2297719226571644e-12` |
| Total Hamiltonian difference | `5.662137425588298e-15` |
| Signed interface-power difference | `1.3322676295501878e-15` |
| Independent MPI2 CN total-balance defect | `9.610368056911511e-16` |

**Stage B: independently accepted full time-domain
energy-power identity** on the same canonical Gaussian
formulation, using MUMPS MPI1 and
`dt=0.005` to `t_f=15`:
**3,000 accepted CN steps**, with no rejected
steps or SNES failures, and **3,001 genuinely
saved consecutive states**. The `mpi1_long/preflight-rank0.json`
shows `accepted_PETSc_TS_steps=3000`,
`actual_saved_FE_states=3001`,
`missing_accepted_grid_indices=[]`
and `FE_state_shape=[3001,7689]`.
All data were obtained from real GetFEM/PETSc
assembly and quadratures, not interpolation
or a synthetic FE trajectory.

| Quantity, complete 15-second run | Measured native value |
|---|---:|
| `H(0)` | `0.7836421211876778` |
| `H(15)` | `0.008268389052222169` |
| Maximum independent total CN step energy-power residual | `4.017619570362285e-15` |
| Maximum heat-subsystem CN step residual | `1.3808398868775384e-15` |
| Maximum wave-subsystem CN step residual | `3.978506976989871e-15` |
| **Maximum cumulative energy-power residual** | **`5.697699100711716e-15`** |
| Maximum signed interface-power residual | `3.1086244689504383e-15` |
| Positive energy increases among accepted consecutive steps | **0** |
| Accepted original states absent from saved FE trajectory | **0** |

Independent offline inspection of the user-supplied
ZIP checked the native `rank0.npz` arrays:
the 3,001 monotone accepted timestamps,
all 7,689 unknowns in the first/final states,
finite and nonnegative heat/wave energies,
nontrivial positive thermal dissipation,
the signed interface sum,
and independently recomputed
`ΔH + Δt (D-P_interface-P_outer)`
**at every adjacent CN pair**.
The cumulative source `CN_cumulative_total_balance_error`
agrees with the independent acceptance.
The maximum discrepancy between adjacent
saved clock differences and prescribed `dt`
is `4.717294263545391e-13`,
consistent with accumulated floating-point
representation; **no accepted step is missing**.

The final Hamiltonian agrees **bitwise in printed
double precision** with the earlier separately
completed 15-second Gaussian visualization run:
`H(15)=0.008268389052222169`;
the independent source verifier reports
zero absolute difference at both endpoints.
The result establishes the discrete CN
**dissipation and power-preserving interconnection**
on the complete tested horizon; it does not
establish analytical polynomial decay exponents,
mesh/time convergence, general solver robustness,
or an annular whispering-gallery eigenmode.

**Fixed-CN sampling regression resolved:** the previous
source-`42116b58` native long calculation converged
3,000 TS steps but failed *scientific acceptance*
because accepted state `n=2999` was not written.
The source-`97c7606b` run, with the *output-only*
`fixed_cn_accepted_step_is_due` logic, archives
**all 3,001 real accepted FE states** and passes
the original uncompromised physical gates.
The old failed ZIP and diagnostic are preserved as
historical evidence; their status is not rewritten.

**Qualification scope:** reference MUMPS MPI1/2
physical parity is PASS on `t_f=0.004`;
MUMPS MPI1 full `t_f=15` true CN dissipation
and all subsystems power balances are PASS.
The `gallery_packet` 15-second run is a
**separate complete visualization dataset**, not a
full-horizon power-balance PASS or an eigenmode
certification. MPI2 15-second power,
parallel VTK export, SuperLU_DIST and
h/dt convergence remain explicitly unqualified.

## Historical: native MPI parity PASS but 15-second FE state 2999 missing (source 42116b58, superseded)

**Historical incident, resolved by source `97c7606b`:** the following section describes the earlier archive and its correct original FAIL status. It does not describe the latest qualified run.



**Original immutable proof:**
`v130-heat-wave-MPI2-longCN-42116b58-FCxHdEia.zip`,
source exact `42116b5803e5295225394582e6ddf930a4269afa`,
image `scrimp-mpi:v1.3.0-rc3` (pinned revision).
The source tests and the check of all eight
immutable published examples were PASS.

**Short stage — genuinely accepted:**
`short-verification.json` has
`PASS_HEAT_WAVE_MPI1_MPI2_SHORT_PHYSICAL`
for MUMPS MPI1 and MPI2, 7,689 FE unknowns,
four fixed CN steps at `dt=0.001`, with
strict compatible projected/clamped initial
wave data and production automatic
conforming-interface orientation.

All independent heat/wave/total balance,
true positive dissipation, interface
power cancellation, exterior zero power,
Hamiltonian monotonicity and cumulative
identity gates pass **on both ranks**.
The MPI2 ranks recorded **bitwise identical**
complete FE trajectories and physical
quadratures. Original MPI1/MPI2 maximum
absolute differences include:

| Independent native MPI1/MPI2 comparison | Absolute maximum |
|---|---:|
| Complete FE coefficients (`z`) | `2.2297719226571644e-12` |
| Hamiltonian `H` | `5.662137425588298e-15` |
| Interface signed power | `1.3322676295501878e-15` |
| Total CN step balance | `5.0584536559483695e-15` |
| MPI2 max **independent** total CN defect | `9.610368056911511e-16` |

Thus **canonical Heat–Wave physical MPI1/MPI2
parity is now qualified** for the actual
four-step model, not inferred from earlier
research probes or the PETSc success code.
Parallel VTK output remains separate.

**Long stage — scientifically NOT accepted yet:**
PETSc/MUMPS completed **3,000 actual CN
steps** at `dt=0.005` up to `t=15`,
in about 19.4 seconds of TS execution,
**zero rejected steps, zero nonlinear
solver failures**. The true FE
Hamiltonian integration completed in
another approximately 22.4 seconds.
However, the strict original-FE
time-grid gate rejected the run
**before** evaluating separate powers.

An independent audit of all **3,001
PETSc monitor records** in
`mpi1_long.host.log` identifies
**exactly one accepted state missing
from the saved FE trajectory**:

```text
step=2998 time=14.990 saved=YES
step=2999 time=14.995 saved=NO
step=3000 time=15.000 saved=YES
```

The monitor used the historical due-time
rule `accepted_t >= last_saved_t + dt_save`
(with a small 16-ULP guard). Rounding
accumulation made the genuine accepted
state 2999 fail that comparison.
The 3,000 saved FE states instead
of the required **3,001** cannot be
reinterpreted as 3,000 consecutive
CN intervals. No value has been
interpolated or manufactured, and the
original long ZIP retains its
historical `FAIL_OR_UNQUALIFIED_HEAT_WAVE_MPI2_LONG`.

**Generic output-only correction in SCRIMP:**
`scrimp/utils/trajectory_sampling.py` now
exposes `fixed_cn_accepted_step_is_due`.
For `ts_type="cn"`, `ts_adapt_type="none"`
and exactly `dt_save=dt`,
`DPHS.monitor` stores every real
**PETSc-accepted step** without relying
on a floating-time threshold. It
continues to use the original
`saved_time_is_due` schedule for
adaptive solvers, any other
time integrator or subsampled output.
Time integration, solution state,
interface equations, FE types,
DAE consistency and numerical power
acceptance **are unchanged**.

Regression tests reproduce the missing
accepted sample at step 2999 and
require retention of all 3,001 states,
while safeguarding adaptive,
non-CN and subsampled behavior.
The native probe now writes a
`preflight-rankN.json` before gate
evaluation recording accepted
PETSc steps, saved-state count,
missing grid indices, timestamps,
H and FE dimensions. An incomplete
trace **still fails closed**, with
its original evidence archived.

**The 3,000-step per-step physical
power identity remains PENDING** until
the corrected native run has archived
all 3,001 genuine adjacent FE states
and passed independent GetFEM
heat/wave/interface cumulative balances.

## Historical qualification plan (now completed for MPI2 short and MPI1 time long)

The original canonical Heat–Wave PDE has now been
extended with an **output-only** switch
`heat_wave_eq(..., physics_only=True)`,
which runs the same real GetFEM/PETSc equations,
strict algebraic DAE initializer and boundary
`clamp`, creates the same Hamiltonian terms, and
omits only the PNG/PVD writers. The existing
default tutorial output behavior is unchanged.
This resolves the known MPI file-writing conflict
without ever mocking the solve, the FEM or a
physical quadrature.

The two-stage fail-closed acceptance script is:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-heat-wave-mpi-long-physics.sh
```

**Stage A — distributed full FE equivalence:**
run the *same* gaussian canonical Heat–Wave model
with MUMPS MPI1 and MPI2, 4 fixed CN steps of
`0.001`. Require all oriented heat, wave, outer
power and dissipation gates **on every rank**;
archive the complete `5 x 7689` original GetFEM
state vectors from each rank; compare all
coefficient and FE quadrature traces across
MPI1 and MPI2 with explicit
`atol=5e-10, rtol=2e-8`, and require exact
agreement of the two independently saved MPI2
rank histories. An MPI2 physics or parity
failure **prevents the long run**.

**Stage B — the complete 15-second dissipative
energy-power identity:** MUMPS MPI1,
`dt=0.005`, `t_f=15`,
**3,000 genuine CN steps and 3,001 original
accepted-time observations**. The scientific
difference from the already completed long
**visual** simulations is crucial: their
`dt_save=0.05` stored only **301 states**.
They prove the long-time evolution of
`H(t)`, but the CN midpoint identity
`H^{n+1}-H^n=-dt D^{n+1/2}+dt P^{n+1/2}`
cannot be tested by treating their
*ten-step* gaps as one CN step.
The new stage sets **`dt_save=dt`**,
requiring every accepted FE state. GetFEM
reassembles physical observables on
the **real midpoint FE state of every
single CN step**, never a fabricated
subsampled midpoint. It measures:

- independent heat and wave energies
  and their sum, plus nonnegative and
  nontrivial true heat-flux dissipation;
- **each** heat and wave subsystem
  energy balance and independently
  recomputed total CN energy balance;
- signed cancellation of power at the
  automatically assembled two outward
  interface facets and zero outer work;
- **all-step cumulative**
  `H(t)-H(0)+sum dt*D-sum dt*P` and
  nonincreasing total Hamiltonian.

The native probe writes full original
FE histories for the **short MPI parity**
but keeps the **long** evidence compact
(the entire native time/H/energy/power
traces, original first/last FE vectors
and SHA256 of the complete original
state history). Rank-local original
data, run logs, image and source
manifest, independently recomputed
host-gate JSON and a ZIP are preserved
even when a criterion fails.

Default Stage B uses the historical
`gaussian` profile to allow an extra
independent end-state comparison with
the genuine pre-existing
`v130-heat-wave-long-gaussian-b9d56fcc-ySKGPg7G.zip`.
To *separately* audit the alternative
annular excitation, use
`SCRIMP_V130_HW_LONG_PROFILE=gallery_packet`.
For a cheap diagnostic dry run rather
than full-time qualification, set
`SCRIMP_V130_HW_LONG_STEPS=100`;
the report records the actual horizon
and makes no 15-second claim.

**Historical status at source `42116b58`:** the 15-second power proof was PENDING at that time because accepted state 2999 was absent. **This is resolved in source `97c7606b`**, whose 3,001 real original states pass the complete CN power audit. The four-step MPI1/MPI2 parity remains PASS. Successful
previous parallel-interface research
pilots and completed 15-second visual
runs are separate, still-valid
qualifications; none is substituted for
the requested new tests.
Parallel VTK writing, SuperLU_DIST,
modal whispering galleries and
spatiotemporal convergence remain out
of scope.

## Qualification after DAE initial correction — PASS (25 November 2026)

**Source/archived evidence:** canonical
`examples/heat_wave.py`,
commit `b9d56fccc18c1f9facacdd9b6c2645c4de25fc47`,
Docker `scrimp-mpi:v1.3.0-rc3`,
archive `v130-heat-wave-physical-b9d56fcc-3FQ8lFlC.zip`.
Its native `summary.json` reports
`PASS_HEAT_WAVE_PHYSICAL_MUMPS_MPI1`,
**ten strict independent gates PASS,
errors=[]**, and unchanged eight
`examples/published/` files.

The real `Concentric` mesh with physical
outward normals, `initial_wave_boundary="clamp"`
and strict algebraic initializer successfully
completed **four fixed Crank–Nicolson steps**
at `dt=0.001`, `t_f=0.004` using **MUMPS MPI1**.
The complete pH finite element system contains
**7 689 DOFs**. All positive heat/wave energies,
the independent `H_heat + H_wave = H` equality,
strictly positive thermal dissipation,
homogeneous outer work, signed interface cancellation,
two subsystem CN balances, total CN balance and
nonincreasing total Hamiltonian were verified
on actual GetFEM quadratures. Specifically:

| Observed real-GetFEM numerical quantity | MUMPS MPI1 |
|---|---:|
| `H(0)` | `0.7836421211876778` |
| `H(0.004)` | `0.6808984803105300` |
| Maximum CN midpoint heat dissipation rate | `32.04166757122644` |
| Maximum signed internal interface power | `8.881784197001252e-16` |
| Maximum outer-boundary power | `7.73365898273765e-24` |
| Maximum subsystem heat CN balance defect | `4.916206330918271e-15` |
| Maximum subsystem wave CN balance defect | `1.3851766955674805e-15` |
| Maximum independent total CN balance defect | `5.457440055423035e-15` |
| Maximum cumulative total energy-power defect | `3.8437031335739146e-15` |
| Separate physical acceptance gates | **10/10 PASS** |

**Interpretation:** the correction of the
Gaussian initial wave velocity is not
a tolerated inconsistency. It yields
an exactly compatible essential
boundary trace in the *dynamic FE data*
before the strict algebraic initializer.
The production conforming-interface
compiler correctly gives each side
its own outward normal. This native
result validates the combined correction
and true heat/wave interface energy
accounting in the specific short
four-step reference. It does **not**
establish independent MPI2 or
SuperLU_DIST acceptance, and does not
prove global long-time CN subsystem
power balances.

## 15-second true-GetFEM visual results — independently archived PASS

**Two authentic completed native runs,
same source `b9d56fcc` and pinned
Docker image.** The native process exit
code is 0 for each profile, and all
301 original `t=0,0.05,...,15`
Hamiltonian samples, 903 field VTU
(`p,q,T`, 301 each), four annular
PNG snapshots and an original
48-frame magnetic-free wave velocity
GIF are retained in each archive.

| Native long-run profile | Original ZIP | `H(0)` | `H(15)` | Relative decrease |
|---|---|---:|---:|---:|
| Historical `gaussian` | `v130-heat-wave-long-gaussian-b9d56fcc-ySKGPg7G.zip` | `0.7836421211876778` | `0.008268389052222169` | **98.9449%** |
| Exploratory `gallery_packet` | `v130-heat-wave-long-gallery_packet-b9d56fcc-jC3P6XEI.zip` | `0.46264406597912466` | `0.3779415175911232` | **18.3084%** |

**No FE Hamiltonian growth** was measured
between any successive saved times
for either initial profile. The *CG1
kinetic proxy* for localization within
the outer shell `r>=0.82` averages
approximately `0.3803` for the
gaussian versus `0.5871` for the
gallery packet on the 48 displayed
snapshots, with final proxy fractions
`0.3685` versus `0.6687`.
These are clear differences in
**transient kinetic localization**,
not a proof of persistent eigenmode
confinement, spectral modal purity or
theoretical decay rates. A
`cos(14 theta)` standing excitation
alone does not establish azimuthal
*direction of propagation*.
Neither long run has separately
passed all ten full-horizon
subsystem-power gates.

**Cause of original visual-driver FAIL:**
a Python/XML API mistake in the
post-simulation verifier:
`frame["file"]` was applied to
an `ElementTree.Element`, instead
of `frame.attrib["file"]`.
The `native_exit_code=0`, actual
`heat_wave_FE_H_long.npz`,
`provenance.json`, all original
VTK `p,q,T` frames, PNG/SVG and GIF
were **already written successfully**
for both profiles.
The final `summary.json` in each
original archive remains historically
marked `FAIL_OR_UNQUALIFIED_LONG_HEAT_WAVE_SHOWCASE`
because it documents that verifier bug,
**not a physical solver failure**.

**Corrections without resimulation:**
the current
`scripts/run-v130-heat-wave-long-showcase.sh`
uses `frame.attrib["file"]` and
its source AST tests reject the old error.
A new immutable-archive verifier
`scripts/verify-v130-heat-wave-long-archive.py`
checks the exact source and Docker
revision, original exit code,
complete 301-step saved FE clock,
nonincreasing positive true Hamiltonian,
PNG/SVG/GIF signatures and **48 frames**,
the original `p,q,T` PVD times
and **all 903 native VTU files**
(including sampled XML fields).
It emits a **separate `.reaudit.json`
sidecar**, preserving the source ZIPs
and their historical `summary.json`
unchanged. Independent examination of
the uploaded ZIPs confirms that all
these criteria are fulfilled. The
resulting classification is
`PASS_ARCHIVED_LONG_NATIVE_VISUALS_NOT_WG_MODE_PROOF`,
not a retrospective rewriting of the
original pilot's status.

Re-audit the existing ZIPs on the
original host, without recomputing
3 000 CN steps per case:

```bash
python scripts/verify-v130-heat-wave-long-archive.py \
  sandbox/outputs/v130-heat-wave-long-gallery_packet-b9d56fcc-jC3P6XEI.zip \
  sandbox/outputs/v130-heat-wave-long-gaussian-b9d56fcc-ySKGPg7G.zip
```

**Scope**: Heat–Wave MUMPS MPI1 physical
reference **PASS**, two completed
15-second native long simulations
with complete genuine visual
archives **PASS after independent
archive re-audit**. The
*whispering-gallery* modal
classification, decay-rate exponent,
full 15-second energy-power
identities, MPI2 and SuperLU_DIST
remain explicitly **unqualified**.



## 25 November 2026 — DAE initial data incompatibility and correction

**First real four-step MUMPS MPI1 pilot at source
`937718d2`: FAIL before CN.** The actual
Concentric geometry, automatic oriented interface
compiler, 7 689-DOF GetFEM mass and stiffness
matrices were constructed. The *strict* algebraic
initializer rejected the original L2 projection
of the Gaussian wave velocity because it does
not satisfy homogeneous Dirichlet data on
the outer annulus:

```text
InconsistentDAEInitialData:
||F_a|| = 0.00762918411 > 8.18210288e-09
Y_bnd = 0.005395, U_bnd = 0.005395
```

The native evidence is
`v130-heat-wave-physical-937718d2-9O1SrxEY.zip`.
This is **not** an error in the automatic
interface compiler, and no CN physics was
executed in that first run.

This exact difficulty was already resolved
and native-qualified in the prior strict
DAE compatibility studies:
`sandbox/probe_v130_heat_wave_boundary_compatibility.py`
and
`sandbox/probe_v130_auto_interface_mpi_cn.py`.
The model now integrates those two
*explicit physical initial-data* choices
instead of depending on a research
probe's monkey-patch of `DPHS.init_step`:

- `initial_wave_boundary="clamp"`
  (**default**): first project
  the analytic Gaussian into the real
  scalar CG1 wave space by GetFEM/PETSc
  L², then exactly set the **126**
  exterior Dirichlet DOFs to zero,
  without changing any free interior
  coefficient (known Concentric
  1/2 mesh); and
- `initial_wave_boundary="constrained_l2"`:
  keep the exact zero exterior trace
  but minimize the correction
  `||p_h-p_h^{L²}||_{M_p}` using the
  true sparse GetFEM mass `M_p`,
  the PETSc COO bridge and MUMPS;
  only **MPI1** is supported for this
  optional correction until separately
  validated in parallel.

The helper
`scrimp.utils.fe_initialization.enforce_zero_dirichlet_initial_values`
uses the **actual filtered GetFEM DOF
map** from `basic_dof_on_region`,
not assumed nodal ordering, and fails
on malformed scalar traces.
It records the number of fixed and
free DOFs, before/after trace norm,
correction amplitude and method.
Neither method changes the prescribed
heat `T` or stress `q` fields,
the weak forms, DAE tolerances or
the strict algebraic initializer.
`"unconstrained"` reproduces the
historical incompatibility **only when
explicitly requested** for diagnostics.
When the heat occupies the outer
annulus there is no exterior wave
Dirichlet condition to correct.

**Also corrected in the native pilot:**
the first original-example inspector
was evaluating interface power on
the old *doubled logical region* 10.
The updated inspector instead
collects independently on the two
**outward-oriented physical facet
regions** allocated by the production
compiler, using the already audited
`measure_v130_hw_dual_facet_physics.py`.
The physical gate now checks that
this collector and the actual
boundary-correction provenance were
used. The previous MPI interface
proofs remain valid; these edits
only repair the standalone main
example and its qualification.

**Status:** source corrected and
source tests prepared, but this
specific native physical pilot
must be **rerun**. No new PASS
Heat–Wave physics is claimed
from these changes alone.

## Long-time scientific visualizations (optional; native pending)

To show both physical time decay and the
**outer-annulus wave dynamics**, the *canonical*
`examples/heat_wave.py` now supports the
optional `--long-showcase` mode. A separate
driver executes it with `heat_region=1`,
`wave_region=2`, `t_f=15`,
`dt=0.005`, `dt_save=0.05` and
MUMPS MPI1:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-heat-wave-long-showcase.sh
```

The default `SCRIMP_V130_HW_SHOWCASE_PROFILE=gallery_packet`
is a high-azimuthal-wavenumber
`cos(14*theta)` wave-velocity packet,
radially concentrated around
`r≈0.86` and vanishing analytically
at both boundaries. It is *an excitation
designed to reveal a possible
whispering-gallery-like pattern*,
**NOT a certified eigenmode**.
The previous `gaussian` physical
initial condition remains the default
for ordinary `heat_wave_eq` and
the four-step acceptance pilot.
To render the latter over the
full 15-second horizon instead:

```bash
SCRIMP_V130_HW_SHOWCASE_PROFILE=gaussian \
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-heat-wave-long-showcase.sh
```

All figures are **generated in the
original Python example**, without
external rendering scripts:

- log-log of the total **original FE
  Hamiltonian**, excluding `t=0`,
  saved as PNG and SVG along with
  `heat_wave_FE_H_long.npz`;
- authentic original GetFEM
  `p`, `q` and `T` PVD/VTU
  sequences, one file per saved
  physical instant;
- four annular-wave PNGs and
  a GIF of at most 48 **actual
  saved FE states**, with fixed
  color range, no interpolation
  in time;
- an explicitly approximate
  concentration diagnostic: the
  fraction of **CG1 kinetic
  velocity proxy** `p²/2`
  represented by triangles with
  centroid at `r>=0.82`.
  It is *not* the full wave
  energy `(p²+|q|²)/2`
  nor a spectral localization test;
- `provenance.json` binding
  the figures to the source,
  original saved times, radial
  fraction sequence and native FE
  Hamiltonian values.

The native driver archives logs,
VTU/PVD, PNG/SVG/GIF, the real
FE trajectory and `summary.json`;
on failure it does not mislabel
the rendering as a PASS.
Only the *short pilot* verifies
the full signed CN physical
power-balance gates; neither
a long-time plot nor a
near-outer-boundary pattern
alone validates a whispering
gallery eigenmode or a
theoretical decay exponent.

**Status of new long showcase:
native execution still pending.**

## 25 November 2026 — next qualification stage prepared (native PENDING)

After the separate **Maxwell 3D Nédélec-II/BDM1 MPI1/2 PASS**,
Heat–Wave 2D is the next unqualified example in the current
Sphinx progression. A new dedicated, fail-closed driver is now
available:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
conda activate scrimp

SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-heat-wave-physical-pilot.sh
```

**Stage 1 only:** MUMPS MPI1, actual
`examples/heat_wave.py` with heat on the inner disk,
wave on the outer annulus, `dt=0.001`,
`t_f=0.004`, four genuine *fixed* Crank–Nicolson
steps. There are no modifications of the
physical equations, source FEM assignments,
algebraic initializer, solver stack,
or `examples/published/`.

Unlike the earlier diagnostic run,
**PETSc converged is not the acceptance condition**.
The new driver reads the original full FE trajectories
and independently sampled `get_quantity(...,CN=True)`
heat/wave observables, and requires simultaneously:

1. positive separate heat/wave Hamiltonians,
   independent equality to the original total `H`;
2. nonnegative **and nontrivial** thermal flux dissipation;
3. zero *signed* internal interface power
   `-U_T Y_T + U_w Y_w` and zero exterior power;
4. the **heat**, **wave** and **total** individual
   CN midpoint energy-power balances, both by accepted
   steps and cumulative total; and
5. nonincrease of the total energy and
   original `T`/`p` GetFEM VTK plus
   the public Hamiltonian figure.

The ten named gates use a bounded absolute
`1e-8 * max(1,max|H|)` threshold
(except positivity and a nontriviality
check). The original native
`sandbox/measure_v130_heat_wave_physics.py`
remains the source of physical observables;
no favorable values are manufactured.
Tests compile **the actual Python inside
the Bash heredoc before Docker**, check that
the run uses the true shipped physical model,
and enforce unchanged publication snapshots.
CI runs only these source tests, never claims
native integration success.

**On an unsuccessful native run** the driver still
generates `summary.json`, source/image manifest,
`native.host.log`, any GetFEM/PETSc files already
created, and
`sandbox/outputs/v130-heat-wave-physical-*.zip`.
The native summary is explicitly
`FAIL_OR_UNQUALIFIED_HEAT_WAVE_PHYSICAL`.
Only all ten physical gates on real GetFEM
data can produce
`PASS_HEAT_WAVE_PHYSICAL_MUMPS_MPI1`.

**Not yet qualified:** MPI2 ParaView output,
SuperLU_DIST, the full default 2 s horizon,
mesh/time convergence, and exact reproduction of
the cited paper's div-div wave formulation.
Further backend or MPI stages follow only
after the physical MUMPS MPI1 gate. The existing
published source snapshots remain byte-for-byte
untouched.

The author confirmed the correct scientific reference:

G. Haine, D. Matignon and F. Monteghetti,
*Long-Time Behavior of a Coupled Heat-Wave System Using a
Structure-Preserving Finite Element Method*,
Mathematical Reports **24**(74), 1–2 (2022), 187–215.
https://imar.ro/journals/Mathematical_Reports/Pdfs/2022/1-2/11.pdf

The **original paper's numerical experiments are different** from
the SCRIMP `examples/heat_wave.py` example. No modifications,
even cosmetic, are allowed under `examples/published/`.

## Mathematical target (paper eqs. 2.1–2.2, 3.15)

With unit physical coefficients, heat inside `Omega_1`,
wave inside `Omega_2`, exterior wave displacement fixed:

```text
d_t T = Delta T                    in Omega_heat
d_tt w = Delta w                   in Omega_wave
T = d_t w = p                      at Gamma_interface
grad(T).n_heat = -grad(w).n_wave   at Gamma_interface
w = 0                              at Gamma_wave_exterior
```

In co-energy variables:
`e_Q=-grad T`, `p=d_t w`, `q=grad w`,
and the outward physical interface normals have opposite
directions. The energy is

```text
H = 1/2 int_heat T^2 + 1/2 int_wave (p^2+q.q)
dH/dt = - int_heat e_Q.e_Q = -D,   D>=0.
```

The outer Dirichlet `w=0` corresponds to fixed *velocity*
`p=0`, if this is compatible with its initial displacement.
It is **not** the Neumann condition `q.n=0`.

## What was incorrect or insufficient in the old example

The preceding `heat_wave.py` used mixed co-energy
variables but unsuitable/opaque combinations of FEMs
(`T:DG1`, `e_Q:CG2`, `p:CG2`, `q:DG1`,
interface/boundary `DG1`), an adaptive BDF2 time scheme,
and a purported Dirichlet outer boundary whose actual
wave boundary equations prescribed **zero traction**,
which is Neumann.

In addition, with the heat div-div weak form
`M_Q e_Q = int T div(test_e_Q) - int U_T test_e_Q.n`,
a **minus** sign is required for the temperature
Dirichlet input on the heat boundary. The former
`+U_T test_e_Q.Normal` instead prescribed
`T=-U_T`. The former `U_w=-Y_T` simultaneously
compensated this mistake in the algebraic energy
interconnection, so an apparently correct *zero
algebraic interface-power sum* did NOT by itself
establish the intended physical transmission law.
Both signs must be corrected together.

The previous native PETSc BDF2 Heat–Wave integration
already completed at t=.004 (39 accepted, 4 rejected,
0 SNES failures, 26,198 DoFs); it was the **old pilot's
false four-fixed-step acceptance test** that failed,
not MUMPS or the actual physical solver. That
historical evidence does not qualify the new model.

## New implementation in examples/heat_wave.py

The previously qualified SCRIMP/PETSc/GF **solvers and
PDE backend are NOT modified**. Only this **non-published
example** is revised.

### FEM choice and mathematical honesty

The paper's mixed FE spaces are (4.3)–(4.5):
```text
T: DG_k; e_Q: RT_k; p: DG_k; q: RT_k;
interface and outer boundary variables: CG_(k+1).
```

For its time-decay experiments the paper uses `k=0`.
The current SCRIMP example selects
```text
k=0
T      DG0
e_Q    RT0  (GetFEM FEM_RTK(2,0), intrinsic H(div))
p      CG1  (because the existing wave uses grad–grad weak form)
q      RT0  (intrinsic H(div), with proper normal traces)
U_T/Y_T, U_w/Y_w, U_bnd/Y_bnd: CG1 restricted to boundaries
```

**This is a mixed PFEM that remains different from the
paper's DG0/RT0 wave formulation.** To implement the
latter exactly, the wave weak forms and port causalities
must also be rebuilt to the paper's div-div system
(3.18–3.19) and its interface DAE (4.7–4.9).
An unverified partial replacement of `p:CG1` with
`p:DG0` while leaving `Grad(p)` or
`q.Grad(Test_p)` would remove the spatial wave
operator. Do **not** make that invalid shortcut.

Native `S.FEM("Heat flux",0,"RT")` and
`S.FEM("q",0,"RT")` require the current GetFEM build
to support `FEM_RTK(2,0)`. The example explicitly
checks that neither was silently replaced by CG.
GetFEM reference:
https://getfem.readthedocs.io/en/latest/userdoc/appendixA.html

### Boundary, coupling and model power

The heat div-div and wave grad-grad weak forms remain,
but their **physical signs** are now explicit:

- `U_T=Y_w=p`: temperature on the heat interface
  equals wave velocity.
- `Y_T=e_Q.Normal` uses the common GetFEM normal.
  To obtain physical outward heat flux,
  `orientation=+1` for heat in inner disk,
  `orientation=-1` for heat in outer annulus.
- `U_w=orientation*Y_T` imposes physical wave traction
  matching outward heat flux. The heat boundary
  brick is `-orientation*U_T*Test_e_Q.Normal`.
- The true *signed* interface power is
  `P_int=-orientation*U_T*Y_T+U_w*Y_w`, expected
  to vanish. An unweighted `U_T*Y_T+U_w*Y_w`
  sum is **not** the physical power balance here.
- For an outer **wave** boundary, a Lagrange-multiplier
  port enforces `p=U_bnd=0` (instead of setting
  zero traction); with heat on the exterior,
  `T=U_bnd=0` is imposed using the heat div-div
  temperature input. This follows the intended
  Dirichlet boundary conditions.

### Fixed Crank–Nicolson and DAE initialization

The example explicitly sets
`ts_type="cn"`, `ts_adapt_type="none"`,
`init_step_ts_type="algebraic"` and preserves
`dt=.001, dt_save=.05, t_f=2.` by default.
Users can still supply different `time_step`,
`save_step`, `final_time` and `factor_solver`
for bounded tests; if unspecified, LU=MUMPS as
in the qualified SCRIMP core. The initial Gaussian
profiles remain as before pending physical diagnostics
of the **accepted, algebraically consistent** initial
state; this does not claim that the raw prescribed
Gaussian initially satisfies every constraint.

No implicit adaptive BDF2 property is asserted.

## First native mixed-RT integration result: core initializer bug (OPEN)

Actual qualified Docker `scrimp-mpi:v1.3.0-rc3`, source
`8286f18dc87ba97b4bc690e5dbd85316f84e9f14`,
evidence ZIP `v130-main-2d-lu-8286f18d-zY4Twrpi.zip`
(11 entries, integrity PASS). The run was stopped at
`DPHS.set_initial_value("q","[0.,0.]")` **before**
mass/stiffness assembly, consistent initialization,
PETSc solve, Hamiltonian calculation or output validation.
The archive contains the source log and initial mesh,
but **no accepted numerical trajectory or physical
energy results**. Hence no new physical conclusion.

All intrinsic mixed FEMs instantiated successfully:
`T DG0`, heat flux `RT0`, wave velocity `CG1`,
wave stress `RT0`, interface and boundary `CG1`.
The wave-region RT0 restriction produced **3111
coefficients** (whole MeshFem `nbdof=4794`,
physical vector `qdim=2`). The legacy
`DPHS.set_initial_value` incorrectly forced
`reshape((qdim, int(3111/qdim)))` and raised
`ValueError: cannot reshape array of size 3111 into
shape (2,1555)`. RT/BDM intrinsic FE coefficients
are face flux moments, not component-wise nodal values.
A valid odd-length 3111-vector is completely normal.

SCRIMP core is now corrected to pass the exact
original **flat GetFEM coefficient vector** after
`basic_dof_on_region` filtering, with validation of
source evaluator coefficient count and region DOF
indices. No zero-padding, truncation, reordering,
FE element substitution, or weakening of native
scientific convergence checks. The change covers
any intrinsic-vector basis with a non-componentwise
DOF count, not merely this example.

Standalone regression
`tests/test_v130_rt_initial_values.py`
executes the **actual production method** using
controlled GetFEM-shaped fake arrays: RT0 odd
3111 DOFs, zero RT state, legacy vector CG,
scalar and unfiltered behavior, bad indices,
and mismatched evaluator size. Its source CI
is a separate guard: native GetFEM/PETSc must
still validate the actual RT0 initialization.

## Second native attempt: false nonlinear classification (OPEN)

Actual source `f97fcd1f866abf9e24591b068dc7a1651d271e69`,
original host archive
`v130-main-2d-lu-f97fcd1f-RpX0FgJi.zip` (11 ZIP
members, CRC PASS). The corrected general `DPHS.set_initial_value()`
successfully initializes **all 3111 RT0 wave coefficients**
without losing an entry. Original GetFEM linear mass/stiffness
assembly completes: model **7689** full DAE DoFs.
The failure occurs **before PETSc's physical CN integration**:

```text
ValueError: Algebraic initializer supports only linear mass and stiffness
```

Source diagnosis: `DPHS.set_control("Interface Heat","Y_w")`
and `DPHS.set_control("Interface Wave","Y_T")` created feedback
bricks with `linear=False` (historically essential for
nonlinear state-dependent feedback). Yet `U_T=Y_w` and
`U_w=Y_T` are **exact linear coupling relations**: the
weak terms `Y_w*Test_U_T` and `Y_T*Test_U_w` are
genuine GetFEM **bilinear** terms with constant tangent.
The global `linear_stiffness=False` flag was therefore
incorrect for this specific coupled model and correctly
caused the conservative *linear-only* algebraic initializer
to reject it.

**Source correction:** `DPHS.set_control()` now recognizes
only a **bare signed registered GetFEM model variable**
(`Y_w`, `Y_T`, `-Y_T`, etc.) as provably linear.
Such feedback uses `Brick(linear=True)`, so GetFEM
`add_linear_term()` includes the exact cross-variable
coupling in the constant stiffness matrix and in the
native `K_aa` projection. Arbitrary macros, sums,
products, coefficients or unknown expressions are **not**
guessed linear: their GetFEM nonlinear weak form and
Newton derivative remain unchanged. Prescribed open-loop
sources still use the existing independent RHS route.

Dedicated source-only tests execute actual production
`set_control()`, `add_brick()`, `set_linear_flags()`,
and classify source/linear/nonlinear GetFEM assembly
routes with negative controls. **No native PASS yet.**
A subsequent run must check the actual fixed-state
projection, its true algebraic `IFunction` residual,
the four fixed CN steps and the separate physical
heat/wave power balances. In particular, a structurally
singular coupled `K_aa` can still produce an explicit
initialization rejection and must not be masked by a
pseudo-time fallback.

### Intrinsic vector FE initial coefficients: generality and limits

The `DPHS.set_initial_value()` fix is **not RT0-specific**.
It flattens the exact `MeshFem.eval()` coefficient array
and retains `basic_dof_on_region()` indices in the
original GetFEM order; it does not reshape by `qdim`.
That structural rule also applies to intrinsic
`RT_k`, `BDM_k` and Nédélec edge basis functions.
However **coefficient-layout compatibility is not
complete family qualification**: element availability,
interpolation, Piola transformations, partial-region
index matching, normal/tangential traces, and actual
Maxwell/PDE weak-form solves require native tests.

The SCRIMP `FEM` class currently exposes native
`RT` and `BDM` string selectors; GetFEM documents
`FEM_RTK(P,K)` for simplices (K>=0),
`FEM_BDMK(P,K)` for simplices (K>=1), and
**`FEM_NEDELEC(P)` only for lowest-order edge
elements**. There is no native `Nedelec1`/`Nedelec2`
selector in the current SCRIMP `FEM` class.
GetFEM may accept an explicit descriptor through
SCRIMP's expert `order<0` escape hatch, but no
higher-order Nédélec first/second-family implementation
or support is thereby proven. Missing intrinsic
elements can currently fall back to CG in `FEM.set_fem()`;
production tests must require **strict element
identity** and never silently treat a fallback as PASS.

Therefore: `RT0` coefficient indexing observed in
the real example; higher `RT_k`, `BDM_k` and
Nédélec 1st/2nd families still **UNQUALIFIED**.
This is a pending separate FEM-family validation gate,
not a reason to regress to CG/DG or reshape/truncate
coefficient vectors.

An optional native, non-PDE **element-family coefficient
pilot** has now been prepared. It creates genuine
triangular GetFEM `MeshFem` objects, restricts their
DOFs to a subset of elements, constructs the same
`Model.add_filtered_fem_variable` that production SCRIMP
uses, calls the **actual**
`DPHS.set_initial_value()` on nonzero vector data,
then checks the resulting GetFEM model coefficients
against the independently evaluated ordered restriction.
Candidates: `FEM_RTK(2,0)`, `RTK(2,1)`,
`RTK(2,2)`, `FEM_BDMK(2,1)`,
`FEM_BDMK(2,2)`, and
`FEM_NEDELEC(2)`. A descriptor unavailable in this
particular native image is reported as **UNAVAILABLE**,
never passed. A descriptor that instantiates but fails
coefficient initialization makes the pilot FAIL.
Nédélec higher-order and second-kind elements are
explicitly **not implemented or qualified** by this
probe.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-intrinsic-fem-initialization.sh
```

This independent check produces host-owned JSON and a
small ZIP. It neither runs a physical PDE nor establishes
Maxwell/operator or MPI support. **Native pilot PENDING**.

## Genuine physical qualification — native PENDING

The revised `scripts/run-v130-example-2d-lu.sh`
by default runs only `heat_wave` at MPI1,
MUMPS and SuperLU_DIST (the unchanged
`shallow_water` stage is selected separately).
Each model uses its original geometry and
weak-form/constitutive parameters with explicit
short horizon t=0..0.004, fixed CN dt=.001.
For heat-wave, exactly 4 accepted steps are now
a correct **CN no-adaptation** criterion.
ParaView exporter is not proven MPI-safe: do
**not** run this top-level plotting example on
multiple ranks without a separate export audit.

The independent `sandbox/measure_v130_heat_wave_physics.py`
records the actual original GetFEM quantities:

- separate heat/wave energies and exact sum vs H;
- nonnegative heat flux dissipation;
- each raw interface port, true oriented net power,
  and outer boundary power;
- candidate physical temperature and flux trace errors,
  with unsupported evaluations clearly reported;
- **CN midpoint values** `get_quantity(..., CN=True)`
  for dissipation and port power on each accepted step;
- original per-step and cumulative identities,
  including `ΔH + dt*(D - P_signed)=0`
  and each subsystem's separate balance;
- all entire FE coefficient paths, actual field traces,
  Hamiltonians and public PNG/ParaView outputs;
- cross-LU full-path FE and power-observable
  agreement with an explicit `1e-8` comparator gate
  (MUMPS vs SuperLU_DIST, same fixed CN clock).

Diagnostics are **not automatically a physical-model
PASS**. Confirm energy and interface residuals
against explicit tolerances and physical meaning,
then qualify or correct the model. In particular,
the published article establishes a continuous
and semidiscrete identity, not an unconditional
PETSc time-integration proof for this different
specific weak form.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
SCRIMP_V130_2D_SCOPE=heat_wave \
bash scripts/run-v130-example-2d-lu.sh
```

**Status: source model changed; no genuine native
result yet.** `examples/published/` and every
historical local simulation output remain untouched.
