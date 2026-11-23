# SCRIMP v1.3.0 — Issue #30: physical CN port-power plotting quadrature

**23 November 2026: COMPLETE — original 8-segment CN data re-audit PASS; generic `Port.compute()` and public `plot_powers()` CN/BDF2 × MPI1/2 4/4 native PASS, independently re-audited.**

The separate checkpoint issue #29 has **already** passed. It exposed
a real SCRIMP postprocessing defect, not a problem with PETSc, the
GetFEM weak forms, nonlinearity, port-Hamiltonian mathematical
discretization, or any checkpoint state. The exact native source
is the original independent 2-D linear wave pH-DAE, eight actual
GetFEM/PETSc CN segments recorded at
`sandbox/outputs/v130-issue29-linear-DAE-e795ac44-Bo62JD0X/`.

## Technical explanation

SCRIMP's unchanged `Port.compute()` evaluates raw algebraic
`flow * effort` on each of the existing accepted physical
FE states:

- `P[0]` at `(t[0],z[0])`: **instantaneous power** at the initial
  state of that simulation segment.
- `P[k]` for `k>=1` at
  `(t[k],(z[k-1]+z[k])/2)`: existing **CN mean-state
  step power** indexed by the endpoint. In particular,
  time-dependent model expressions use `t[k]`, not
  `(t[k-1]+t[k])/2`; this fix does **not** pretend otherwise.

The *old* public `DPHS.plot_powers()` wrongly applied
`0.5*dt*(P[k-1]+P[k])`, as though these were endpoint powers.
The error is particularly visible when a new physical segment
starts at a checkpoint: its `P[0]` is instantaneous at the
restored `z(t*)`, and is **not** the preceding physical CN
step's `P[k]`.

For the existing CN values, the *signed raw port integral* must be:

```text
I[0] = 0
I[k] = I[k-1] + (t[k]-t[k-1]) * P[k]
```

This is NOT a generic physical boundary-work orientation.
For the archived independent wave model's specific left-flow
port, supplied boundary work is **`W = -I`**, and the
native discrete relation
`H[k]-H[k-1] + dt[k]*P[k] ≈ 0`
has already been independently validated per step at
roundoff. Other models have their own port signs, area/mass
weights and dissipation terms. The *plot's legacy Balance
curve* keeps the existing sign convention but is explicitly
labeled **model-dependent**; it must not be cited as physical
energy preservation without an independent model-specific
power audit.

## Generic `Port.compute()` sampling contract (new stage)

The first postprocessing fix correctly recovered the original
CN step quadrature in `DPHS.plot_powers()`, but retained the
**non-generic** legacy `Port.compute()`: it unconditionally averaged
`z[k-1]` and `z[k]`, regardless of PETSc's selected TS method.
For BDF, ARKIMEX, or any other TS this could silently label
a CN-averaged `P[k]` as a nodal value and then trapezoid it.

This is corrected structurally, preserving the production
CN physical results and original GetFEM weak forms:

```python
port.compute(solution, gf_model, domain)  # generic instantaneous default
port.compute(solution, gf_model, domain, sampling="instantaneous")
port.compute(solution, gf_model, domain, sampling="cn_step")

port.get_power()           # original signed raw series
port.get_power_sampling()  # "instantaneous" or "cn_step"
```

- **`instantaneous`**: evaluate original GetFEM form at each
  `(t[k], z[k])`, for every time scheme, including BDF and
  ARKIMEX. This is a mathematically well-defined *instantaneous
  port-power sample*, independent of a numerical integrator.
  Its nodal trapezoid is an **approximation**, not a claim
  about BDF's actual discrete quadrature, RK intermediate
  stages, or an exact Hamiltonian identity.
- **`cn_step`**: preserve historical native-qualified
  `Port.compute()` behavior bit-for-bit at each
  `(t[k],(z[k-1]+z[k])/2)`, k>=1. `P[0]` remains
  instantaneous; the original `t[k]` parameter evaluation
  remains **unchanged**, even for time-dependent controls.
  Correct discrete integral: `sum(dt[k]*P[k])`.
- **`DPHS.compute_powers()`** dispatches the sampling based
  on `self.time_scheme.ts_type` (CN vs other TS).
  Consequently, `get_power()` retains its previous CN
  meaning on an unmodified CN simulation, while the other
  time schemes no longer receive a hidden CN averaging
  operation. Direct `Port.compute()` now defaults to
  instantaneous rather than making an unadvertised CN
  assumption. This is an intentional correction of
  previously misleading non-CN outputs.
- **`DPHS.plot_powers()`** reads the stored sampling
  metadata from **each Port**, instead of presuming that
  all arrays share a convention merely because the
  `time_scheme` attribute says CN. It integrates CN-step
  samples by exact step sum and instantaneous samples
  by an explicit *approximate* nodal trapezoid.
  Unknown/external series refuse a misleading integral.
  Legacy external `port.set_power(values)` remains callable;
  callers must set `sampling=` explicitly before requesting
  a corresponding integrated plot.
- This does **not** expose PETSc internal BDF/IMEX stage
  states or build a generic stage quadrature. A future
  interface for exact multi-stage methods needs solver-stage
  output captured during time integration; it cannot be
  reconstructed from only accepted endpoints. Do not
  use energy-conservation assertions on those methods
  before native mathematical/numerical qualification.

The new host-only AST tests **execute the real method body**
with a fake GetFEM object, distinguishing endpoint values,
CN step means, explicit `t[k]`, sampling metadata,
invalid clocks and failed assembly. These tests are NOT
a native PETSc/GetFEM qualification.

### Original native data verification: PASS (user host)

On 23 November 2026 the user executed
`python tests/test_v130_power_integration.py`:
**7/7 PASS**, followed by
`sandbox/verify_v130_power_quadrature_native.py` on the
immutable original #29 eight-wave segments.
`PASS_PUBLIC_CN_STEP_INTEGRAL_ON_ORIGINAL_NATIVE_PORT_ARRAYS`:

- maximal true original per-step
  `abs(delta H + delta t * P[k]) = 2.439454888092385e-19`;
- maximal restarted vs continuous signed physical-work
  discrepancy `7.047314121155779e-19`;
- maximal incorrect historical trapezoid vs
  corrected stage-power work discrepancy
  `2.8536249430776475e-05`.

No physics rerun; original files remain immutable.
This proves the **CN raw port numerical quadrature**,
not the new generic BDF/IMEX sampling mode.

## Narrow code change, deliberately backwards-compatible

- `scrimp/utils/power_integration.py` adds a pure NumPy
  routine `integrate_port_power(times, samples, convention=...)`,
  which validates arrays and uses `cn_step` for true CN,
  intentionally discarding `P[0]` from step integration.
- `DPHS.plot_powers()` dispatches on **actual**
  `self.time_scheme.getString("ts_type")`.
  For CN it uses `cn_step`. For BDF, ARKIMEX and other
  methods it **retains the original trapezoid exactly**,
  labeled legacy/approximate: their time-stage semantics
  remain a separate research/implementation task.
- The generic `Port.compute()` **sampling API is now updated**,
  but the CN state/time evaluation is unchanged. Nonlinear FE
  bricks, physical signs, PETSc TS/SNES/MUMPS, the MPI bridge
  and the accepted-state checkpoint API are **unchanged**.
- Plotted legends explicitly identify the quadrature;
  no unwarranted universal “Balance must be constant” claim.
- Dedicated source-only tests include wrong positive
  trapezoid, initial `P[0]` discontinuity after checkpoint,
  physically signed wave identity, nonuniform accepted times,
  unqualified non-CN backward compatibility and deliberate
  invalid/nonfinite power inputs.

**Native reuse without repeating old PETSc integrations:**

```bash
cd ~/Simulations/scrimp
git pull --ff-only
python tests/test_v130_power_integration.py

python sandbox/verify_v130_power_quadrature_native.py \
  --evidence sandbox/outputs/v130-issue29-linear-DAE-e795ac44-Bo62JD0X
```

The standalone verification reads eight **immutable original**
`linear_states.npz` files and directly applies the
new public numerical helper to their original native CN
port samples. It compares per-step signed Hamiltonian
increments, original source+restart boundary work and
a negative control deliberately distinguishing the wrong
legacy trapezoid. It never edits or rewrites original
scientific data and requires no GetFEM/PETSc/Docker reruns.

### Reproducible original native results — 23 November 2026

Source `883b428febe5a317ab1832d57d02d2ad121765d4`; pinned
`scrimp-mpi:v1.3.0-rc3`, original host evidence:
`sandbox/outputs/v130-issue30-port-generic-native-883b428f-N4zJAtcX/`
and ZIP of identical name. The uploaded ZIP has 18 members,
all with valid CRC. Every original NPZ array and all four
PNG figures independently re-opened and checked; PNGs
measure 1350 × 750 pixels. Each physical wave simulation
used 426 complete FE/DAE DoFs and 12 accepted time steps.

| PETSc integrator | MPI ranks | Actual sampling | Max full-FE-trajectory MPI L2 difference | Max raw port power difference | Max cumulative port integral difference |
|---|---|---|---:|---:|---:|
| CN | 1, 2 | `cn_step` | 6.711684736112862e-13 | 3.191891195797325e-16 | 4.336808689942018e-19 |
| BDF2 | 1, 2 | `instantaneous` | 2.6174007281914428e-14 | 1.700029006457271e-16 | 1.0842021724855044e-19 |

All four original GetFEM/PETSc runs PASS, with zero
MPI intra-run generic port-power discrepancies, and
original DAE algebraic residuals <=3.02e-17. The physical
ports are nontrivial (max raw generic port-power magnitude
0.02836248 in CN and 0.02945672 in BDF2).
Independently recomputing every element of each signed
port integral from the original NPZ reproduces the archived
values exactly (maximum numerical difference **0.0**).

For CN, the genuine original wave step identity satisfies
`max |ΔH+Δt*P|=2.1006417091906648e-19`
(MPI1) and `9.147955830346444e-20` (MPI2);
the final model-specific supplied work balances the change
in Hamiltonian at roundoff.

For BDF2, the **approximate** signed nodal trapezoidal
integral differs from the final Hamiltonian change by
`-2.1638336660092855e-07` (MPI1). This is **not**
a checkpoint, MPI, GetFEM or plot coding failure, and is
not presented as an exact BDF stage-work discretization.
Such a property would require capture of the actual PETSc
time-stepping stages. The original physics, FE model,
and solver are unchanged.

The real native public plotting method was exercised in
each case; every Matplotlib port curve and the historical
model-dependent `Balance` curve were checked against
independent integrals in the original native process.
A separately encountered `Port.get_power()` malformed
logging call on the deliberately tested *failure* path
was fixed subsequently, without changing numerical outputs.

### New independent public plotting native smoke (PASS on qualified host)

The generic Port.compute AST tests alone did not constitute native evidence. The complementary real GetFEM/PETSc BDF2/CN solves and public Matplotlib artist checks have now been run and archived.

A **separate** qualified-host gate now invokes the existing
original 2-D linear wave DAE with the same domain, constitutive
parameters, time-dependent physical input and original rest
initial condition. For each of genuine PETSc `CN` and `BDF2`
at MPI 1 and MPI 2 (four actual solves, 12 steps each), it
verifies the original full FE state, algebraic constraints,
all rank-generic algebraic port samples, exact declared
`get_power_sampling()`, every public `DPHS.plot_powers()`
Matplotlib line versus an **independent** `integrate_port_power()`
numerical result, the legacy signed Balance plot, and
root-only PNG and full original NPZ publication. It
compares the complete MPI1/MPI2 paths and port powers for
each scheme. CN additionally checks the exact original
wave per-step Hamiltonian/port identity; BDF explicitly
**does not** assert a nonexistent exact stage work identity.

```bash
cd ~/Simulations/scrimp
git pull --ff-only
SCRIMP_DOCKER_NO_SUDO=1 \
SCRIMP_V130_DAE_IMAGE=scrimp-mpi:v1.3.0-rc3 \
bash scripts/run-v130-issue30-port-generic-native.sh
```

The driver pins the previously qualified `rc3` OCI
revision, overlays only the updated Python source, and
creates a unique host-owned
`sandbox/outputs/v130-issue30-port-generic-native-.../`
folder and a standalone ZIP of native logs, JSON, full FE
traces, port samples and public plots. It never overwrites
or deletes previous evidence. All four native runs completed successfully, and the original ZIP was independently checked on 23 November 2026.

**Issue #30's scoped acceptance conditions are now satisfied:** the generic sampling dispatch and genuine end-to-end native public `plot_powers()` checks passed on the qualified GetFEM/PETSc host. The native #29 re-audit is
already complete; preserve the honest "approximate" description
for other BDF/IMEX time integrators.
This targeted fix does not replace the separate #6
native-qualified physical h-weighted dissipative SWE
power audits, and does not imply a new proof of the
nonlinear Hamiltonian discretization.
