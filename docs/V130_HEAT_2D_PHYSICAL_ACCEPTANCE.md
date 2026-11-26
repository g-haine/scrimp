# SCRIMP v1.3.0 — canonical 2-D heat equation: native qualification

> **26 November 2026 — new v1.3 example candidate, SOURCE STAGED,
> native pending.** The original BDF2/pseudo demonstration documented
> below is the *historical reference*, fully qualified at source
> `9270b32`. The updated `examples/heat.py` now proposes CN with
> linear-algebraic consistency initialization; the user's original
> FEMs, weak forms, boundary controls and initial temperature are
> unchanged. Its new self-contained true-GetFEM visualizations will
> produce `outputs/showcase/heat/figures/heat_temperature_snapshots.png`,
> `heat_temperature_evolution.gif`, `heat_energy_and_ports.png`
> (plus SVG and genuine complete temperature VTK/PVD) when the
> canonical script is run. **None of these new CN numerical or
> visualization claims is native-qualified yet.** First run
> `scripts/run-v130-heat-cn-pilot.sh` (four steps, MPI1/2),
> then qualify the full canonical horizon and visual outputs.
> Preserve the BDF2 qualification independently using
> `scripts/run-v130-heat-physical-diagnostic.sh`.
>
> **Native CN short pilot 26 November — all physics PASS, global FE parity
> FAIL (not released):** original user archive
> `v130-heat-CN-pilot-1fc43520-LFg8YI8b.zip`; 16 ZIP members,
> CRC PASS. Four genuine fixed CN steps at `h=.1, dt=.01, t_f=.04`
> MPI1/2, zero TS rejections or SNES failures. Every per-rank genuine
> CN/port/Hamiltonian gate PASS, maximum CN step defects
> `1.5588e-14` (MPI1) and `9.2756e-15` (MPI2). Both MPI2 rank FE
> replicas are byte-identical. Nonetheless the global MPI1-vs-MPI2
> full-FE coefficient discrepancy is `1.99827153e-6` (criterion
> `1e-8`), so the pilot status is **UNQUALIFIED**, NOT PASS.
> The discrepancy already appears at the initialized `t=0` state
> and, remarkably, satisfies `delta_z[k] = (-1)^k delta_z[0]`
> on all 5830 DOFs to at most `2.3e-13`; it does not grow over time.
> This is consistent with the CN algebraic-constraint alternating
> parasitic mode excited by rank-dependent finite-precision LSQR
> consistency initialization. `init_algebraic_accuracy=1e-11`
> was not reached: MPI1 final residual `1.8196e-11`, MPI2
> `1.3157e-11`, after the default 3 passes; both are still
> within the looser physical feasibility criterion, without
> changing prescribed differential state data.
>
> **26 November, first “strict” run INVALID as a precision experiment:**
> archive `v130-heat-CN-strict-pilot-24229e1e-o5UV4ge2.zip`
> (ZIP CRC verified). Its three ranks faithfully reproduce the
> default initialization and CN result: 3 LSQR passes, accuracy target
> `1e-11`, `accuracy_target_met=false`, and global FE MPI discrepancy
> `1.99827157e-6`, despite excellent CN port-power balances.
> Root cause: the first probe set `PETSc.Options()` BEFORE calling
> `examples.heat.heat_eq()`, but the `DPHS` constructor invokes
> `get_cleared_TS_options()` and clears the PETSc options database.
> Hence **none of the requested stricter initializer settings was
> actually transmitted to SCRIMP**. This result does NOT establish
> that more accurate algebraic projection is ineffective.
> The obsolete pre-constructor probe has been removed.
>
> **Corrected STRICT diagnostic (SOURCE STAGED, native pending):**
> `SCRIMP_V130_HEAT_CN_INIT_PROFILE=strict bash scripts/run-v130-heat-cn-pilot.sh`
> now uses the **canonical** `examples.heat` main entry point with
> `--init-profile strict`, adds the four values to its
> `heat.set_time_scheme(**scheme)` call AFTER the `DPHS` constructor,
> and persists the actually selected values in `init_profile_applied`.
> The host checker verifies the four option values, the algebraic
> accuracy target, original energy balance and unchanged full-FE
> parities. A repeat native run is required before any numerical claim.
>
> **Next diagnostic (SOURCE STAGED, native pending):**
> `SCRIMP_V130_HEAT_CN_INIT_PROFILE=strict bash scripts/run-v130-heat-cn-pilot.sh`
> reuses the exact original CN model, FE spaces, boundary
> controls, plot and 4-step MPI1/2 checks but increases the
> already-supported algebraic initialization max passes to 8,
> targets accuracy `1e-13` and sets LSQR tolerances `1e-14`.
> The physical and strict all-DOF MPI parity gates are **not relaxed**;
> additionally the LSQR accuracy target must actually be met.
> This is an isolated acceptance experiment, **not** the newly
> qualified default. Preserve any failure archives and diagnose
> convergence before choosing final time-stepper/initialization.
>
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
The canonical default without an explicit `--init-profile` requires
its own independent full-horizon native verification.

**Next required release gate — NATIVE PENDING:**
`scripts/run-v130-heat-cn-full-showcase.sh`.
This first checks the true GetFEM VTU-to-PNG/GIF renderer with a
4-step original-model preview, then runs the canonical long
`t_f=5, dt=.01` CN example with the automatic strict initializer,
genuine 501-frame FE/VTK/ParaView output, and original-user PNG/SVG/GIF
figures in MPI1, followed by the full 501-state MPI2 simulation.
Its host checker tests each original quadrature-energy gate and
every MPI1/MPI2 FE coefficient with the unchanged `1e-8`
threshold, and verifies all VTK timestamps and genuine files.
All original VTK remain in the host evidence folder; only selected
hashed raw VTK are included in the more manageable evidence ZIP.
No PASS for the full horizon/showcase is claimed before user-host runs.

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
- Reproduction of this same instrumentation from the current branch:
  `scripts/run-v130-heat-physical-diagnostic.sh`. The original user
  entry point remains `python examples/heat.py` in the qualified environment.
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

This acceptance closes the historical BDF2/pseudo heat model's physical
gate. The upgraded canonical CN/algebraic user example has a separate
pending qualification and is NOT covered by this verdict. No change was made to
the original heat equations, bricks, FE spaces, boundary controls,
or integration settings.
