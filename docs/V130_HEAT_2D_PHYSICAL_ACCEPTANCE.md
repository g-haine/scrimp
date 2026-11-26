# SCRIMP v1.3.0 — canonical 2-D heat equation: native qualification

**STATUS (26 November 2026): PASS — improved canonical CN/strict-algebraic
heat 2-D example.** Original GetFEM weak formulation, conductivity,
mass, FE spaces, mixed control ports and prescribed temperature data
retained. All 500 fixed CN steps, `h=.1`, `dt=.01`, `t_f=5`,
MPI1/2, strict all-DOF FE parity, dissipation, full quadratic CN
step-power identity, original user PNG/SVG/GIF and genuine GetFEM/VTK
history passed on the user host. Historical BDF2/pseudo simulation
remains separately qualified and reproducible.

## Full 500-step CN + strict algebraic canonical example: NATIVE PASS

**Native user-host evidence (original calculations at source revision
`92b9d870`):**
`sandbox/outputs/v130-heat-CN-full-showcase-92b9d870-gszcYbW2/`
and the companion
`v130-heat-CN-full-showcase-92b9d870-gszcYbW2.zip`.
Docker `scrimp-mpi:v1.3.0-rc3`, image revision
`1d6d0c90115213298047b28237a7f38087563a35`.
The **original, superseded ZIP** had 36 members and SHA256
`f7e42d747bb2b6059b73b50185cf31c3cfa489207894228e11916d83ced558e9`.
**Do not use that SHA256 for the newly repacked archive**: its bytes
changed when the selected raw VTUs were inserted. Its new SHA256
has not been independently supplied in the user-host transcript.

- Original `examples/heat.py`, not a research replacement, automatically
  chooses `ts_type=cn`, `init_step_ts_type=algebraic` with the **per-example**
  strict options `init_algebraic_accuracy=1e-13`,
  `init_algebraic_max_passes=8`, and both LSQR tolerances `1e-14`.
  No change to the generic SCRIMP DAE initializer.
- The original canonical initial state is preserved exactly in the
  differential temperature coefficients; initial algebraic
  residual `4.1672128285e-14` (MPI1) and
  `3.7133806177e-14` (MPI2), target `1e-13`: **PASS**.
  Only three of the eight permitted LSQR correction passes were
  required. No pseudo-time initialization.
- Original 5,830 FE DOFs, **501 authentic accepted snapshots**,
  zero rejected steps and zero nonlinear-solver failures.
  All **500** original GetFEM Hamiltonian values decrease
  monotonically from `7.02705436928772` to
  `5.764488338085319` (MPI1).
- All 501 Fourier dissipation samples are positive, minimum
  `0.13101852399680086`; the internal Fourier constitutive
  port agrees with the independently assembled quadratic form.
  The four *true weak control-port* powers are the relevant boundary
  terms: trace-based `T*(J_Q.Normal)` is not a substitute.
- Independently recomputed 500-step CN energy balance
  `H(5)-H(0) = -1.2625660312024012`,
  `sum(dt D_CN) = 0.9309314591751581`,
  `sum(dt P_ports,CN) = 0.33163457202721114`.
  The maximal **per-step** quadratic CN identity residual is
  `5.72745e-14` (MPI1) and `2.93489e-14` (MPI2);
  maximal cumulative residual `6.75016e-14` and
  `3.20e-14` respectively. **All 500 steps passed**, not
  merely the final aggregate.
- Full MPI1-to-MPI2 comparison across all 501 times and all
  **5,830 coefficients**: maximum absolute error
  `4.9181476669e-9`, below the **unchanged `1e-8`** gate;
  Hamiltonian difference `5.77316e-14`, Fourier dissipation
  difference `1.61537e-14`, and the largest original individual
  control-port power difference `1.25802e-14`.
  MPI2 rank 0 and rank 1 FE histories are **byte-identical**.
- The *canonical user example itself* generated
  `outputs/png/Hamiltonian_Heat_2D.png` and
  `outputs/showcase/heat/figures/` containing:
  `heat_temperature_snapshots.png/.svg` (four true DG1
  temperature maps with one common colormap),
  `heat_energy_and_ports.png/.svg` (quadratic Hamiltonian,
  Fourier dissipation, four-port sum and CN closure),
  and `heat_temperature_evolution.gif`
  (**36 genuine saved FE frames**, checked with Pillow).
  The PNGs were opened and visually inspected; no synthetic fields.
- `outputs/pv/T/T.pvd` lists exactly 501 genuine
  `T_00000.vtu` through `T_00500.vtu` files from GetFEM with
  exact original accepted timestamps `0.00,0.01,...,5.00`.
  The host run checked all 501 native file presences, sizes,
  SHA256 hashes and PVD timestamps, storing the complete hash
  inventory in `native_vtk_inventory.json`.
  Original VTK field files remain on the user host; the compact ZIP
  includes all real FE state histories and figures.

**Archival packaging correction: CLOSED — native user-host PASS
(26 November 2026).** The initially inspected 36-file ZIP inadvertently
contained **zero raw VTUs**, despite a manifest that claimed their
inclusion. The corrected
`sandbox/compare_v130_heat_cn_full_showcase.py` was subsequently
run **directly on the original evidence directory, without Docker or
another TS simulation**. The user supplied the full host output:

- `status=PASS_NATIVE_HEAT_CN_500_STEPS_MPI1_MPI2_AND_VISUALS`,
  `errors=[]`; all previously established MPI1/2 full-horizon
  scientific and visualization gates independently recomputed.
- `selected_actual_GetFEM_VTU_in_archive=8` (four genuine original
  preview VTUs and four genuine original long-run MPI1 VTUs).
- `ZIP_CRC_pass=true`; the script checks the actual raw VTU bytes
  extracted from the new ZIP against their unchanged user-host
  originals using SHA256, **before** atomically replacing the ZIP.
- The entire authentic 501-state `T.pvd`/`T_*.vtu` series
  remains in the original host evidence directory, indexed and
  hashed in the archive's `native_vtk_inventory.json`.
- The archive filename is unchanged, but the **pre-repack archive
  SHA256 is obsolete**. A new digest must be computed if the
  archive is moved or independently distributed.

No new simulation was required and none of the equations, time
integrator, initialization algorithm, or results was modified. The
remaining packaging concern is **resolved**. The rebuilt ZIP has
been validated by the user-host checker; no claim is made that
the rebuilt ZIP itself was independently uploaded and reopened
within this conversation.

**Acceptance decision: FULL PASS; evidence archive CLOSED.**
The canonical `examples/heat.py` v1.3.0 CN/strict-algebraic case is
scientifically qualified for the tested configuration, with
reproducible native PNG/SVG/GIF and complete user-host VTK series.
No remaining heat-example acceptance action is required for this
specified scope. The historical BDF2/pseudo case also remains
qualified under its own separately documented scope.

**Not claimed:** mesh or temporal convergence, other LU backends,
MPI4+, or generalization to unrelated physical models.

**26 November 2026 — NATIVE PASS (bounded: original long BDF2 model, MPI1/2,
all 499 order-two steps, port powers, full FE parity, and user figures).
PETSc BDF2 start-up half-step is not claimed as independently qualified.**

## 26 November 2026 — improved CN/strict algebraic short pilot: PASS

**Native archive:** `v130-heat-CN-strict-pilot-fcec569e-I3Fw63IG.zip`;
SHA256 `ea5aadd97dbfb622362b83b06c5335527ec85d31bb88e2f66e54951e024b4ad6`.
Checked ZIP CRC (16 authentic files), source `fcec569`,
`scrimp-mpi:v1.3.0-rc3` OCI revision
`1d6d0c90115213298047b28237a7f38087563a35`.

The original unmodified heat weak formulation `h=.1`, `dt=.01`,
`t_f=.04`, fixed CN, linear-algebraic DAE initializer, original
temperature initial value, all original four boundary controls, MUMPS
and COO were executed on MPI1 and MPI2 without rejected steps or
nonlinear failures. The stricter options were **actually transmitted
after DPHS construction**, confirmed by each original native JSON:

- `init_algebraic_max_passes=8` (upper limit; convergence required
  only three passes).
- `init_algebraic_accuracy=1e-13`.
- `init_algebraic_ksp_rtol=init_algebraic_ksp_atol=1e-14`.
- MPI1 algebraic residual: **4.16721e-14**; MPI2 **3.71338e-14**.
  Both have `accuracy_target_met=true`, and preserve every
  differential temperature DOF exactly without pseudo-time.
- Maximum stepwise exact quadratic CN power defect: **1.55891e-14**
  (MPI1) / **9.27557e-15** (MPI2); cumulative defects
  **2.33702e-14** / **1.04031e-14**.
- All **5,830 FE coefficients**, across all five accepted times,
  have maximum MPI1-to-MPI2 discrepancy **4.91798e-9**, below
  the **unchanged** `1e-8` threshold. Hamiltonian discrepancy:
  **2.04281e-14**. MPI2 ranks 0/1 are identical.
- The previously identified sign-alternating CN algebraic
  initialization mode is now at small, controlled magnitude:
  its MPI difference alternates with error at most
  **2.50525e-13**. Numerical parity passed, not merely
  coincident Hamiltonians.

**Decision:** the native four-step MPI1/MPI2 CN/strict pilot **PASS**.
The historical loose/default initializer had `1.99827e-6`
coefficient discrepancy and remains **UNQUALIFIED as the default for
this example**. To prevent shipping this known non-parity configuration,
`examples/heat.py` now uses the **strict initializer automatically
for CN only**; this does not alter the generic SCRIMP defaults.
The BDF2/pseudo configuration still selects its legacy initializer.
The canonical default without an explicit `--init-profile`
**has now passed** its own full-horizon native MPI1/2 verification,
as recorded above.

**Completed full-horizon release gate — NATIVE PASS:**
`scripts/run-v130-heat-cn-full-showcase.sh` validated the
automatically selected strict initializer in the canonical example,
500 steps, all original FE/power variables, all source-generated
figures and original GetFEM/VTK/PVD series in MPI1 and MPI2.
The raw-VTU ZIP packaging fix is also validated and closed as
documented above.

## Original source and authentic native evidence

- Canonical source: `examples/heat.py`, unchanged; source revision
  `9270b32b2863fdd678eb35e5b812412be92e6acc`.
- Docker: `scrimp-mpi:v1.3.0-rc3`, OCI revision
  `1d6d0c90115213298047b28237a7f38087563a35`.
- Genuine user-host archive:
  `v130-heat-physical-9270b32b-ti8W2pDK.zip`
  (SHA256 `e3e840bf9730109742fe37e18951dad475bd9201d810a9422c03269c15527d6f`).
  Independently opened and checked: **16 members; ZIP CRC PASS**.
  Original MPI1/MPI2 logs, full accepted FE trajectories, GetFEM
  quadrature histories for the Hamiltonian, Fourier dissipation and
  four boundary controls, original PNG figures, and JSON manifests.
- Reproduction of the historical BDF2 instrumentation from the current
  branch: `scripts/run-v130-heat-physical-diagnostic.sh`. The current
  canonical user entry point `python examples/heat.py` now defaults to
  CN/strict-algebraic; to reproduce the BDF2/pseudo baseline explicitly,
  run `python examples/heat.py --time-stepper bdf --initializer pseudo --no-showcase`.
- Eight tracked `examples/published/` sources remain immutable.

## Exact scope

The actual original model, **not** a replacement manufactured problem:
2-D rectangle `(0,2) x (0,1)`, `h=0.1`, DG1 temperature, CG2 heat flux,
Fourier conductivity `Lambda=0.01 I`, `rho=3`; three prescribed
temperature controls `U_B=U_R=U_T=1` and the original left mixed
Neumann control `U_L=0.2*T`. Original explicitly requested
`init_step_ts_type=pseudo`, original PETSc BDF2 `dt=0.01`,
`t_f=5`, MUMPS and COO, **5,830 total GetFEM degrees of freedom**.
The complete accepted FE time history `z[k]`, not merely the final
state or a post-processed synthetic trajectory, was compared.

| Evidence | MPI1 | MPI2 (both ranks) |
|---|---:|---:|
| Accepted steps | 500/500 | 500/500 |
| Rejected / SNES failures | 0 / 0 | 0 / 0 |
| Saved FE states | 501 | 501 on both ranks |
| Hamiltonian initial | 7.027001061737507 | 7.0270010617374865 |
| Hamiltonian final | 5.764470067122002 | 5.764470067122038 |
| Minimum Fourier dissipation | 0.13095032461114486 | 0.13095032461114492 |
| Max mismatch `D=f_Q.Lambda.f_Q` vs `f_Q.J_Q` | 1.11e-16 | 1.11e-16 |
| Max **exact quadratic BDF2** step power defect, `k>=2` | **3.7453e-12** | **2.0892e-12** |
| Max corrected **trapezoidal approximation**, entire trajectory | 3.2621e-5 | 3.2621e-5 |
| Original `Hamiltonian_Heat_2D.png` | PASS | PASS |

Numerically checked on original NPZ arrays:
`H(t)` is monotonically decreasing on **all 500 accepted time steps**,
the 501 Fourier dissipation samples are strictly positive, and
MPI2 ranks 0 and 1 have byte-identical full FE states. The maximal
absolute MPI1-to-MPI2 differences over **all 501 samples** are:

- Full FE coefficient: **2.31815e-13**.
- Hamiltonian: **5.77316e-14**.
- Fourier dissipation: **1.77636e-15**.
- Sum of four native boundary-control powers: **6.52256e-16**.
- Discrete BDF2 power defect: **4.97352e-12**.

The original PNG was generated by the canonical example itself. No
separate renderer is required by an end user.

## Power signs: use the true discrete control ports

A preliminary diagnostic incorrectly equated all four boundary
exchanges with pointwise `int_Gamma T (J_Q.n)`. That observable is **not**
equal to the port products in this mixed DG / control-multiplier FE
formulation: the imposed boundary effort and the trace of the DG field
are distinct discrete quantities.

The correct original semi-discrete power balance is

`dH/dt + D + P_B + P_R + P_T + P_L = 0`,

where `D = int_Omega f_Q.Lambda.f_Q` and the `P_*` are the original
signed GetFEM **control-port** flow-effort products (`U_B,Y_B`,
`U_R,Y_R`, `U_T,Y_T`, `U_L,Y_L`), retaining their original
orientations. In particular the fourth control is Neumann-like.

Using the pointwise traces gave a false accumulated residual
**0.0307890**. With the genuine port products, a **trapezoidal**
approximation to BDF2's continuous energy balance gives a maximum
cumulative residual **3.2621e-5** and final residual **1.7824e-5**.
This is an **approximate quadrature**, never an exact BDF2 claim.

## Stronger BDF2 quadratic energy test

Write `a=z_k`, `b=z_{k-1}`, `c=z_{k-2}`, and use the original
quadratic temperature Hamiltonian

`H(z) = (1/2) int_Omega rho T(z)^2`.

The verified GetFEM quadratic identity for each genuine fixed
step `k=2,...,500` is

```text
[ H(a) - H(b) + H(2a-b) - H(2b-c) + H(a-2b+c) ] / (2 dt)
    + D(a) + P_B(a) + P_R(a) + P_T(a) + P_L(a) = 0.
```

This follows from testing the *actual* original BDF2 semi-discrete
equation against the temperature at the accepted state. Every term
was assembled using the original GetFEM integration method and
actual FE coefficients; the defect is at solver roundoff scale,
not a favorable property of the plotting routine.

**Start-up limitation:** An auxiliary test previously labeled
`BDF1_initial_step_power_defect` gave about `-0.00452998` on
the first interval, but it assumed an unsupported single
backward-Euler start between `z[0]` and `z[1]`. PETSc
`TSBDF_Restart`, for BDF order > 1, creates an internal
intermediate implicit stage at `t_0 + dt/2`; its full startup
history was not recorded in this archive. Thus that **BDF1
surrogate is not a valid failure gate** for this BDF2 run, nor
does it constitute a proof of startup energy identity.
The original simulation still completes its first step and all
subsequent steps and its total Hamiltonian decreases.

Reference PETSc implementation:
<https://petsc.org/main/src/ts/impls/bdf/bdf.c.html>
(`TSBDF_Restart`).

## Acceptance decision and restrictions

**PASS for the original heat 2-D scientific example at the tested
configuration**, covering the original full trajectory, actual
physical Fourier dissipation, canonical port powers, all 499
order-two BDF2 quadratic power identities, MPI1/2 complete-state
parity, and original figures.

**NOT CLAIMED:** an exact total-horizon BDF2 quadrature identity,
an independently verified internal startup-stage power balance,
temporal/mesh convergence, a different solver backend, MPI4+, or
an independent partitioned-mesh scalability test.

This section qualifies specifically the historical BDF2/pseudo heat
model. The current canonical CN/strict-algebraic example was qualified
separately above. Neither configuration changes the original weak
equations, bricks, FE spaces or boundary controls.
