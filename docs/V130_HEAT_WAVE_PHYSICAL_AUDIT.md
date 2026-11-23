# SCRIMP v1.3.0 — Heat–Wave 2D: physically meaningful example

**23 November 2026 — physical MODEL REVISION COMMITTED IN SOURCE;
native GetFEM/PETSc and CN power-balance qualification PENDING.**

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
