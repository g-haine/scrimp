# SCRIMP v1.3.0 — Issue #30: physical CN port-power plotting quadrature

**23 November 2026: SOURCE-LEVEL FIX STAGED, historical native
data verification PENDING, public plotting smoke PENDING.**

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
- `Port.compute()`, nonlinear FE bricks, original physical
  signs, PETSc TS/SNES/MUMPS, the MPI bridge and the
  accepted-state checkpoint API are **unchanged**.
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

**Do not mark issue #30 CLOSED until:** this re-audit and
at least one genuine end-to-end public `plot_powers()`
native plotting smoke demonstrate the corrected curve and
the unqualified BDF/IMEX modes remain backward compatible.
This targeted fix does not replace the separate #6
native-qualified physical h-weighted dissipative SWE
power audits, and does not imply a new proof of the
nonlinear Hamiltonian discretization.
