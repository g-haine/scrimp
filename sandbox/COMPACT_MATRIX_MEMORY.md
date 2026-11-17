# SCRIMP: compact PETSc matrix memory (2026-11-17)

**Status:** implemented **opt-in**, native physical A/B validation
pending. Not enabled by default, and not yet a quantified memory gain.

## Context

After the validated rank-partitioned asynchronous trajectory storage,
the next opportunity concerns duplicate sparse PETSc Mat objects, not
the retained time samples.

Historical `DPHS.allocate_memory()` initializes six PETSc matrices:
`mass`, `nl_mass`, `tangent_mass`, `stiffness`,
`nl_stiffness`, `tangent_stiffness`; it then creates a
**separate** TS Jacobian `J`. For purely linear bricks, the tangent
operators coincide with their assembled linear components:

```text
tangent_mass      = mass
tangent_stiffness = stiffness
J                 = stiffness + sigma * mass
```

**Important:** `J` MUST remain a separate PETSc Mat for SNES/KSP/MUMPS,
because its values are updated in place when the TS shift changes.
PETSc TS owns its implicit integration and nonlinear histories; no
change to these objects is proposed.

For state-dependent mass, the FD correction required by the implicit
Jacobian is

```text
  (M(z+eta*zd) - M(z-eta*zd)) / (2*eta).
```

Historically, this creates two temporary PETSc state vectors and
**three full sparse copies** (plus, minus and correction) every
time the mass Hessian is evaluated. These are eventually destroyed,
so this is **allocation/peak-memory pressure**, not a leak.

## Opt-in policy

Run with `SCRIMP_MATRIX_MEMORY_MODE=compact`. Omit the variable
(or set `legacy`) for exactly the preexisting behavior.

In compact mode:

1. If `linear_mass=True`, alias `tangent_mass` to the **assembled**
   `mass` Mat; if `linear_stiffness=True`, similarly alias
   `tangent_stiffness` to `stiffness`. Destroy only the unused,
   **uninitialized** tangent Mat placeholders, before any TS/KSP
   references are created. Do not run empty nonlinear tangent
   assembly that would call `zeroEntries()` on an aliased base.
2. If either operator is nonlinear, its tangent stays independent
   and uses the original `zeroEntries()+AXPY` assembly. A mixed
   linear/nonlinear DPHS is supported.
3. In the nonlinear mass Hessian, use a **single persistent PETSc
   Mat** holding the + perturbation. Assemble the - perturbation
   directly into the existing `tangent_mass` and subtract it in
   place. Scale the result and add it to the TS Jacobian `P`.
   The Mat is **reused across IJacobian calls** using PETSc
   `Mat.copy(result=..., structure=DIFFERENT_NONZERO_PATTERN)`.
4. Use one persistent perturbation PETSc Vec. Copy the original
   Newton state into it **before each AXPY**, so both `z+eta*zd`
   and `z-eta*zd` have the same rounding behavior as before.
5. Restore GetFEM's original trial state and tangent mass in
   `finally`, including after exceptions. The nonlinear
   sparsity pattern may change: never assume `SAME_NONZERO_PATTERN`.

Two counters in DPHS expose the allocation pattern:
`_mass_hessian_fd_matrix_allocations` (target 1),
`_mass_hessian_fd_matrix_reuses` (growing after the first call),
`_mass_hessian_fd_probe_allocations` (target 1).

PETSc `Mat.copy(result=..., structure=...)` and
`Vec.copy(result=...)` are documented reuse APIs; these must still
be tested with the user's actual GetFEM 5.5/MPI-2 and PETSc build.
This change does not alter the nonlinear-mass Hessian formula, MUMPS
reset/reuse policy, COO bridge, or the time integrator.

The empty `nl_mass`/`nl_stiffness` placeholders still exist to
preserve the existing external DPHS matrix attributes. Optimizing
those small placeholders is lower priority than avoiding populated
tangent copies and repeated sparse factor-related allocations.

## Tests already implemented

`tests/test_matrix_memory.py` runs without PETSc/GetFEM. It checks:

- correct linear alias lifetime (base untouched; placeholder destroyed);
- correct central FD result against an analytically known affine
  mass tangent;
- exactly **one reusable matrix and Vec across five Jacobian
  evaluations**, rather than fresh copies;
- restoration of the original Newton/GetFEM state, including when a
  synthetic assembly raises.

`tests/test_compact_matrix_comparison.py` checks the comparison
tool's rejection of nonlinear balance drift, improper matrix sharing
and failure to reuse the Hessian workspace.

## Short native MPI physical parity (to run next)

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash -n sandbox/run_compact_matrix_parity.sh
bash sandbox/run_compact_matrix_parity.sh
```

The launcher uses the existing `scrimp-mpi:local` image and runs **two
physical cases**, both in `legacy` and `compact` mode (four short
solves total):

| Case | Mesh | MPI | TS | Steps | Storage |
|---|---|---:|---|---:|---|
| Linear acoustic waves 2D | h=0.04 | 2 | fixed CN, dt=0.001 | 3 | memory |
| Nonlinear viscous Saint-Venant open spin-up | h=0.35 | 2 | fixed CN, dt=0.001 | 3 | memory |

Both test variants use the **same in-memory output backend** so the
new asynchronous writer does not confound matrix-memory measurements.
The Saint-Venant test retains the original BDF2 consistency
initialization, the physical facet-tangent boundary, and the smooth
initial wall-velocity spin-up over `0.03` time units.

The A/B analyzer requires:

- convergence with **three** CN steps and unchanged time grid;
- full coefficient-wise parity of ALL saved states (not just norms);
- Hamiltonian and port powers (acoustics), and both physical power
  and volume balance histories (Saint-Venant);
- expected PETSc Mat aliases for linear operators, **no aliases**
  for nonlinear operators; one Hessian matrix+Vec across callbacks;
- per-rank `VmRSS` and `VmHWM` collected immediately after
  `solve()`, before any test-only full-trajectory materialization.

Output: `sandbox/outputs/compact-matrix-<sha>-<unique>/`,
including full per-case logs and
`matrix_memory_comparison.json`.

**Limitations:** The high-water RSS includes Python, GetFEM, PETSc,
MUMPS and the OS allocator, not merely CSR matrix payloads. The
matrix-memory mode changes only storage/assembly lifetime. Do not
claim any actual memory saving until a native A/B run passes, and do
not infer timing scaling from a single short run.

After validation, investigate whether the baseline allocation cost
is mainly PETSc AIJ matrices or MUMPS factorization workspace, before
making additional structural changes.
