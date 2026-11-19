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

## Native outcome: compact matrix A/B PASSED (2026-11-17)

The user ran the native four-case 3-step MPI experiment from
`17e38b4c` and supplied `compact-matrix-17e38b4c-lXSRAS.zip`.
The archive was read directly, including all per-rank JSON and
the complete NPZ trajectories.

| Model (2 MPI ranks) | Legacy max-rank VmHWM | Compact max-rank VmHWM | Change |
|---|---:|---:|---:|
| Acoustic 2D, 47,942 DDL | 508.949 MiB | 505.375 MiB | **−3.574 MiB** |
| Dissipative nonlinear SWE, 1,170 DDL | 386.238 MiB | 386.625 MiB | **+0.387 MiB** |

Both models passed the strict **full coefficient-by-coefficient**
parity test (not only endpoint state or norms). The maximum
state coefficient error was 1.110223e−16 for acoustics and
2.775558e−16 for SWE. The Hamiltonians were **exactly equal**;
SWE energy balance differed at most 4.235165e−22 and the
volume-balance histories matched exactly. Both SWE runs completed
3 CN steps with no rejects or nonlinear failures.

In the compact SWE run, **each rank created only one Hessian
matrix and one perturbation vector**, and reused the same matrix
**four times** across Jacobian evaluations. Legacy still allocated
temporary matrices repeatedly. The correct workspace behavior is
therefore independently verified on real PETSc/MPI.

**Interpretation:** On the 47k-DoF linear case, removing populated
tangent copies produces a measurable but relatively small
whole-process memory reduction (~0.7%). On the much smaller SWE
system, `+0.39 MiB` is NOT evidence of a meaningful memory
regression; the independent 1170-DoF copies are small relative to
library/native allocation overhead, and one-run VmHWM is sensitive
to allocator layout. Correctness and allocation-count goals PASS,
but nonlinear peak RSS savings are not demonstrated.

For acoustic compact, rank 1 contains approximately **327 MiB
RssAnon** and **190 MiB RssFile**. The original 47k-DoF
memory-mode process retains >500 MiB even after duplicate
operator removal. More small NumPy optimizations will not remove
the native solver/process baseline.

## Next diagnostic: isolate the dominant phase, no long run

Before tuning MUMPS or changing GetFEM matrix representation, identify
**where** memory actually increases. `scrimp/utils/phase_memory_audit.py`
and `DPHS._memory_phase` now provide *opt-in local-only* samples,
with no new MPI collectives, no forced GC and no solver changes.
Pass `SCRIMP_MEMORY_PHASE_AUDIT_DIR` per MPI-run output directory.

The probes cover the PETSc allocation stage, GetFEM mass/stiffness
and RHS assembly, creation of TS Jacobian/vectors, optional
BDF2 initialization, the main TS entry, the first monitor callbacks,
the final accepted solution, and TS cleanup. The **first accepted
main monitor step** occurs after the first linear factorization;
an RSS increment across this phase may implicate TS/KSP/MUMPS
workspace, but cannot alone prove ownership of individual bytes.

Run the **same short three-step** A/B command:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
bash sandbox/run_compact_matrix_parity.sh
```

In addition to the physical comparisons it now creates
`matrix_phase_memory_summary.json`,
`matrix_phase_memory_timeseries.csv` and
`{acoustic,swe}/{legacy,compact}/memory_phase_rank000{0,1}.json`.
Both rank-local `VmRSS` and `VmHWM` are recorded.

RSS measurements include Python, GetFEM, MPI, PETSc, MUMPS and
the OS allocator; memory deltas between phases do **not** precisely
assign bytes to any library. The previous native archive remains
the unchanged pre-instrumentation numerical reference.

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
