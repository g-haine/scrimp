# SCRIMP v1.3.0 — release plan and acceptance gates

**Scope agreed 18 November 2026.** This document is the authoritative
*v1.3.0 work plan*; [V130_DEVELOPMENT_ROADMAP.md](V130_DEVELOPMENT_ROADMAP.md)
retains the historical performance diagnostics, negative experiments
and outstanding technical limitations.

**Version principle:** v1.3.0 is a backward-compatible scientific
release. No general Saint-Venant 2D tangential-velocity boundary-control
development enters its critical path. That model is a separate,
**post-release** project for LHMNC 2027.

## 0 — Finish solver optimization and establish scientific equivalence

**Already achieved (production wave, h=0.04, 47,942 unknowns):**
three independent complete runs with 10 accepted fixed-CN steps at
1/2/4/6 MPI ranks. All current MPI final and initial states agree
within the prescribed tolerance, MUMPS parallel symbolic analysis
is effective, and native symbolic/numeric factor counts equal
exact collectively detected CSR pattern/value changes.
Current solution-time medians: **1.501 / 1.224 / 0.945 / 0.926 s**
at 1/2/4/6 ranks, respectively; efficiency at six ranks is
**27.0%**, a documented strong-scaling limitation.
The 1/2-rank native factor-reuse smoke passes with 1 symbolic,
2 numeric and 10 triangular solve operations over 10 CN steps.

**Resolved cross-version discrepancy (18 November 2026).** The
previous unqualified three-repeat campaign on the branch before the
initialization fix found scaled L2 state differences 0.376/0.397.
The cause was missing GetFEM synchronization of the accepted pseudo-TS
state after correctly excluding prescribed source bricks from
nonlinear Jacobian assembly (see diagnostic below).

**Native v1.2.0 → corrected current release comparison PASSED** on
commit `24e70e92f3db357121e6eb41b1e0e38dc6c94913`, with the
original tagged v1.2.0 commit
`88e2c1442d6aca3c7971c856e049c26d9d1caf0a`.
The user supplied the complete native Docker archive for direct
independent inspection. Fixed wave.py, h=0.04, **47,942 DOFs**,
CN dt=0.002, t_f=0.02, 10 accepted steps:

| Same-model sequential comparison | v1.2.0 native | Corrected candidate |
|---|---:|---:|
| Full wall solve including assembly | 73.267577 s | 1.512175 s |
| Setup (benchmark wrapper) | 0.095716 s | 0.096254 s |
| Mass assembly + PETSc transfer (logged) | 15.15 s | 0.0452 s |
| Stiffness assembly + PETSc transfer (logged) | 56.08 s | 0.1693 s |
| First accepted CN step | 0.39369 s | 0.38129 s |
| Mean CN steps 2–10 | 0.13113 s | 0.05357 s |
| Peak single-process RSS | 592.2 MiB | 582.6 MiB |

**One-run, state-qualified end-to-end ratio: 48.45x**. Most of
this gain is faster linear assembly/GetFEM-to-PETSc transfer, *not*
a 48x faster MUMPS solve; mean steady-CN step is ~2.45x faster.
Do not present the one-run ratio as a 3-repetition median or
extrapolate to nonlinear systems.

Independent full-vector comparison of archived NumPy snapshots:
the complete initial PETSc states are **bit-for-bit identical**;
final state L2 difference is **6.2150211e-30** (relative L2
**3.923971e-31**, maximum absolute component delta
**1.2809612e-30**). Current PETSc main CN records exactly
**1 symbolic LU, 2 numeric LU, 10 triangular solves**, matching
10 Jacobian callbacks with 2 exact numeric updates and 8 reuses.

Reference archive created by
`SCRIMP_V130_BENCH_RANKS=1 bash sandbox/run_v130_wave_benchmark.sh 0.04 1`,
in local `sandbox/outputs/v130-wave-h0p04-24e70e92-F8rBHU/`
(the output directory is deliberately ignored by Git). Both
sequential versions were rerun once from original source tags;
neither program was altered during the timed run.

**Scientific scope:** matching full endpoints validates the numerical
regression, **not yet** the continuous-in-time power/Hamiltonian
balance or fidelity of pseudo-TS initialization to the physically
prescribed `q(0)=0`. The corrected branch's MPI scaling at
1/2/4/6 ranks and full trajectories still need requalification.

**Focused source-level diagnosis (pending native confirmation).**
In v1.2.0, `set_control()` declares prescribed input source bricks
with `linear=False`. The old `DPHS.set_linear_flags()` therefore
marked the entire stiffness as nonlinear and `DPHS.monitor()`
copied every accepted pseudo-TS state back into GetFEM variables.
The current implementation correctly excludes RHS-only source bricks
when classifying Jacobian nonlinearity, so that *implicit* state
synchronization no longer occurs in the linear acoustic example.
This is a plausible explanation for the 17,640-component initial-state
difference, **not yet proof** that restoring the legacy synchronization
is physically correct.

A second change affects the residual: the current default is
`SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem` rather than v1.2.0's
legacy reconstruction from the tangent matrix. A quick A/B probe
is available, using the already archived v1.2.0 vectors; it runs only
two short **current-version** simulations and does not touch the
production solver:

~~~bash
bash sandbox/run_wave_initial_sync_probe.sh
~~~

It reports named GetFEM variable ranges, the pseudo-TS accepted
state, full initial/final vector differences against archived v1.2.0,
and `png/wave_initial_named_blocks.png`:

- current with accepted pseudo-step synchronization, direct residual;
- current with the same synchronization and legacy residual.

**Decision gate:** identify whether v1.2.0 agreement arises from
an undocumented pseudo-initialization side effect, changed residual
semantics, or both; then check initial conditions and DAE constraints
mathematically. **Do not change the production model to match an
unvalidated legacy numerical artifact.**

### 18 Nov diagnostic: root cause of cross-version state mismatch confirmed

The two current-version experimental runs with replayed pseudo-TS
monitor synchronization reproduce the archived **unmodified v1.2.0**
complete state exactly at initialization, and to ~6.22e-30 absolute
L2 at t_f, for **both** `direct_getfem` and `legacy` residual forms.
Without replay, the initial and final L2 differences were ~5.82
and ~6.29. The accepted pseudo-TS state is thus the actual source of
the state mismatch, **not** an effect of residual formulation.

The accepted pseudo-state has field norms
`q=0.00471761, p=14.36784, e_q=0.0202897, e_p=5.76279`,
while before the pseudo step `q=e_q=e_p=0`. The user-prescribed
`q(0)=0` is therefore perturbed by pseudo-time advancement; the
initialization's *physical* DAE-consistency and fidelity to prescribed
differential initial conditions still require a separate check. Do
not claim scientific validity simply from equality with v1.2.0.

**Production fix implemented:** after the *accepted* pseudo-TS solve,
unconditionally scatter its global PETSc state and transfer it to
GetFEM variables, before constructing the main TS initial vector.
This is independent of `linear_mass`/`linear_stiffness` and preserves
the correctly classified prescribed source bricks and MUMPS reuse.
A direct unit regression tests that both linear and nonlinear flag
configurations transfer the accepted state, without aliasing memory.

**Native production qualification PASSED:** the corrected
current source and untouched original v1.2.0 source were both run
once, with equal full initial/final states (numeric details above).
The optional standalone regression
`sandbox/run_wave_initial_state_regression.sh` remains useful
when both reference JSON and NPZ are available; the full clean-run
benchmark above is the authoritative regenerated reference.
The cross-release PNG/SVG comparison was generated and the 
original diagnostic mismatch is resolved.

Next: native current 1/2/4/6-rank full-vector parity after the fix,
then an independent scientific DAE/port-power check.

Exit gates:

- [x] Explain historical cross-version state mismatch:
  initial pseudo-TS state was not propagated to GetFEM after fixing
  the source/nonlinearity classification. Both experimental replays
  reproduce the native v1.2.0 complete trajectory endpoints.
- [x] Verify the production pseudo-state sync patch with a
  fresh native v1.2.0 and corrected-current sequential comparison;
  complete initial states bitwise identical and final states agree
  to ~4e-31 scaled L2; exact native MUMPS reuse audit PASS.
- [ ] Verify physical DAE consistency and whether pseudo-time
  advancement respects the intended prescribed initial q(0)=0.
- [ ] Correct any regression exposed; keep both versions' numerical
  outputs and the build/solve timing categories distinct.
- [ ] Compare full *time trajectories*, not just initial/final
  vectors, and check discrete Hamiltonian, control-port work and
  dissipative contribution when applicable.
- [x] Current 1/2/4/6 rank full-state and exact LU-reuse
  performance gate at h=0.04, three repetitions.
- [x] Re-run a fresh original-source cross-release sequential
  comparison after the initial-state synchronization fix (one
  repetition, fully matched initial/final states, ×48.45 measured).
- [ ] Verify corrected source at 1/2/4/6 MPI ranks and repeat
  qualified sequential performance before presenting a robust median;
  then qualify larger h=0.01 wave runs or document remaining limits.
- [ ] Document the additional independent candidate Jacobian's
  memory cost on a larger mesh and during long integrations.
- [ ] Inspect remaining native MUMPS/TS performance without
  introducing another solver optimization into the v1.3 scope.

**Git decision at this gate:** keep working on the existing
`feature/getfem55-mpi-container-scalability` branch until the
equivalence issue is understood. Then cut a single long-lived
`release/1.3.0-preparation` integration branch. Do not prematurely
merge an unqualified solver change to `master`.

### Final point-0 size-qualification campaign: h=0.005

As requested, run the **complete** current-version strong-scaling
campaign at **h=0.005** (approximately three million DOFs; actual
GetFEM `nbdof()` is authoritative). This is **not** a pilot and
never launches v1.2.0. The source is the corrected production
`examples/wave.py` with exact `on_change` Jacobian reuse,
fixed Crank–Nicolson (10 steps, dt=0.002), compact distributed COO
and parallel MUMPS symbolic ordering.

~~~bash
bash sandbox/run_v130_large_wave_scaling.sh 0.005 3
~~~

Runs **1, 2, 4, 6 MPI ranks**, **three repetitions** each, pinned
by default to physical CPUs 0,2,4,6,8,10. Real timing and
per-process memory outputs are recorded for *all* twelve runs,
with complete initial/final state parity across rank counts and
repetitions and a native PETSc event gate matching exact MUMPS
symbolic/numeric updates to collective Jacobian changes.
No release worktree, comparison or rebuilt v1.2.0 reference is
created. The script stops on an actual failed native test and
preserves all prior successful outputs.

Expected deliverables in the generated ignored
`sandbox/outputs/v130-large-wave-.../` folder:
`v130_large_wave_performance.md`, `v130_large_wave_summary.json`,
per-run CSV, scaling CSV, raw PETSc/SCRIMP logs, complete PETSc
reference states and **12 paired PNG/SVG figures**: run spread,
total/solve scaling, speedup, efficiency, CN startup/warm and
per-step evolution, phase/MUMPS costs, COO assembly/bridge, LU
reuse call counts and per-rank/summed RSS high-water marks.
**Summed per-rank `ru_maxrss` values are not a synchronized
whole-job memory peak**.

**Status:** prepared, synthetic parser/figure tests added; native
large-size host runs and performance measurements are not yet
present and must not be inferred from the earlier h=0.04 results.
The h=0.005 run is the requested final strong-scaling data gate;
it is not a substitute for the separate scientific DAE and
long-trajectory memory limitations in this roadmap.

## 1 — GetFEM 5.5 finite-element family audit and upstream patch

Do **not** reimplement an element simply because another project's
API uses a different name. The current GetFEM documentation describes
`FEM_RTK(P,K)` and `FEM_BDMK(P,K)` on simplices; it describes
`FEM_NEDELEC(P)` as a lowest-order edge element. Tensor-product
geometries have their own availability restrictions; distinguish
these from missing higher-order simplex elements.

1. Reproduce native GetFEM 5.5 element creation, dimensions,
   degrees of freedom, interpolation and assembly in 2D/3D.
2. Test Raviart–Thomas (start k=0..5), BDM (start k=1..5),
   and **both** first/second Nédélec families at their admissible
   orders. Distinguish triangular/tetrahedral from
   quadrilateral/hexahedral support explicitly.
3. Confirm the mathematical definitions and dimension formulas,
   orientation and facet matching, Piola pullbacks,
   H(div)/H(curl) continuity, convergence, exact sequence behavior,
   and serial vs 4-MPI rank results.
4. Implement only genuinely missing functionality as a minimal,
   independently reviewable GetFEM upstream patch, with native
   C++ and Python regressions. Do not carry an undocumented local
   binary patch.
5. Pin the tested GetFEM commit/patch in the Docker installation
   and disclose whether the patch was accepted upstream or
   vendored temporarily.

Exit: supported domain/order matrix with green native GetFEM tests
and reproducible source-level patch for any additions. “Arbitrary
order” must be a parameterized construction, not a hard-coded list
of only five degrees.

## 2 — Docker, native install and platform entry points

- Rework the Dockerfile into cache-friendly dependency layers:
  package download and native PETSc/SLEPc/GetFEM build layers
  **before** copying frequent SCRIMP changes.
- Provide BuildKit cache mounts for apt/micromamba and source
  tarballs, pinned hashes and a resumable/reproducible download
  strategy. Test fail/restart without unnecessary redownload.
- Keep Docker build context minimal; no historical artifacts,
  generated mesh/figures, large ignored outputs or local
  simulation archives should enter it.
- Ship an installer/launcher for **Ubuntu/Debian**, **Windows**
  (document Docker Desktop/WSL2 requirements) and **macOS**
  (including Apple Silicon image constraints). Separate
  supported/validated platforms from best-effort instructions.
- Rebuild the image **locally on the user's computer**, prove
  SCRIMP, GetFEM, petsc4py/slepc4py and MPI smoke tests, then
  remove obsolete image(s) as requested. Avoid deleting the only
  functional image before its replacement passes.
- Provide version-pinned image metadata, licensing, examples and
  volume/permission handling.

Exit: a reproducible clean image, repeatable warm-cache build
and documented cross-platform entry points. This covers #25.

## 3 — New 3D Maxwell structure-preserving example

- Develop a scientifically meaningful 3D electromagnetic
  port-Hamiltonian/PFEM formulation with compatible **H(curl)**
  and **H(div)** spaces, appropriately oriented Nédélec and
  RT/BDM functions, and explicit boundary power/Poynting flux.
- Use a simple geometry with analytically or independently
  checkable initialization, boundary conditions and power identity.
- Establish convergence, time integration, material coefficients,
  energy evolution and discrete power-balance residual.
- Test at 1/2/4 (and, when useful, 6) MPI ranks and at multiple
  mesh and polynomial orders; qualify GetFEM element orientation.
- Generate ParaView-compatible results and explanatory visuals.

Exit: reproducible example, accurate port-power identity, working
MPI scaling and educational narrative. No unsupported claim that
spatial discretization is structure preserving solely because
its Hamiltonian is positive.

## 4 — Performance campaign and documentation figures

- Compare *equivalent* v1.2.0 serial and v1.3 serial wave cases
  only after gate 0; separately report current 1/2/4/6 MPI scaling.
- Repeat comparable 3D Maxwell workloads. Scale by mesh resolution,
  DOF count, number of timesteps and polynomial degree.
- Report setup, GetFEM assembly, GetFEM→PETSc bridge, initial
  consistency, main CN step times, symbolic/numeric LU,
  triangular solve, maximum RSS (per process **and total**),
  and parallel efficiency. Profile first step vs steady steps
  without summing inclusive timers.
- At least three repetitions in quiet host conditions.
- Publish accurate PNG and SVG plots under the documented output
  directory and selected *qualified* figures in Sphinx.

Exit: scientifically justified, reproducible performance results
with complete raw machine-readable data; known 4→6 rank scaling
limits displayed rather than concealed.

## 5 — Rebuild Sphinx documentation, preserving valid mathematics

- Edit/correct retained scientific derivations and tutorials; remove
  obsolete attempts and dead links, not useful mathematical content.
- Explain the complete default SCRIMP stack: PFEM/GetFEM,
  distributed COO bridge, PETSc TS/SNES/KSP, LU/MUMPS,
  parallel symbolic analysis, exact `on_change` numeric and
  symbolic reuse, initialization and memory/parallel limits.
- Document precisely which geometries, elements, nonlinear
  models, MPI operations, boundary conditions and power
  balances have actually passed validation.
- Provide a **VS Code Dev Containers** tutorial: host Docker
  installation, bind mounts, persistent outputs, MPI terminals,
  Python interpreter/debugging, UID/GID permission model
  and CPU/NUMA configuration.
- Ensure HTML, LaTeX/PDF, navigation, examples and screenshots
  render correctly. Maintain the user-level and developer-level
  reading paths.

Exit: reviewed Sphinx documentation, green HTML/PDF build and
reproducible container walkthrough. This covers #26.

## 6 — All examples / optional installation self-check

- Inventory **every** shipped example, then run each with
  its documented inputs in serial and at `mpiexec -n 4`
  (where the model's MPI contract makes sense).
- Check expected physical invariants, convergence/failure codes,
  power balance, short trajectory parity and controlled outputs.
- Fix/mark unsupported models explicitly; do not silently
  proclaim success when an example only imports.
- Supply a bounded opt-in CLI, e.g. `scrimp --check-install`
  or `python -m scrimp --check-install`, with quick serial
  and MPI smoke subsets, dependency information, machine-readable
  summary and meaningful exit code. Extended test suite is
  separate so installation checks remain practical.

Exit: complete example matrix and exercised installation check
on the validated platforms.

## 7 — GUI review / script prototyping role

- Audit which current API calls and options the GUI generates,
  then correct incompatibilities introduced by v1.3.0.
- Demonstrate GUI→Python SCRIMP script generation using at least
  one basic example and a clear round-trip reproducibility check.
- Document the GUI explicitly as a **prototyping aid for SCRIMP
  scripts**, not as a universal PDE editor or replacement for
  advanced scripted MPI/model development.
- Test startup, forms, file dialogs, generated scripts and
  failure messages; document supported desktops, Qt/display
  requirements and known restrictions.

Exit: working GUI smoke, explicit scope and updated screenshots/docs.

## 8 — Clean **tracked** sandbox, not local user outputs

- Preserve the documented `wave_1D` sample as the only
  **tracked sandbox example**.
- Move useful benchmark/research regressions to `tests/`,
  `benchmarks/` or `scripts/`, and update imports and CI.
  Do **not** delete an active test merely to make sandbox empty.
- Drop redundant experiments from the tracked tree after recording
  their findings in documentation and Git history.
- Never delete a user's local ignored `sandbox/outputs/`, mesh,
  animation, archives or benchmarks.

Exit: `git ls-files sandbox` contains only the intended documented
sample (and any narrowly required package/build metadata, if
explicitly justified), with all scientific tests functioning
elsewhere.

## 9 — Prepare and publish GitHub v1.3.0 release

- Run all unit, example, scientific, MPI and GUI gates; verify
  package headers/licenses, clean install, Python dependencies,
  wheels/source distribution, Docker provenance and Sphinx.
- Audit and close only genuinely resolved issues; reference
  fixes and regressions in release notes, retain open issues
  with documented limitations.
- Freeze packaging/docs versions only after full qualification.
- Merge the single reviewed v1.3.0 integration branch into
  `master`; tag `v1.3.0`, publish GitHub release
  (plus container artifact, when qualified), provide
  reproducible install + changelog instructions.

Exit: tagged and reproducible release, source and documentation
versions consistent, passing CI and scientifically validated demos.

## Git and work-stream strategy

- **Now:** finish (0) on
  `feature/getfem55-mpi-container-scalability`.
- **Once gate 0 passes:** create *one*
  `release/1.3.0-preparation` integration branch.
  Keep implementation commits logically separated; no branch
  per checklist item.
- **GetFEM upstream code:** use a separate branch/PR in the
  GetFEM fork because it is a different repository. Pin the
  exact patch SHA in the SCRIMP integration branch.
- **Before release:** review and merge only once into
  `master`. No premature release tag.
- Existing issues #6, #23, #25 and #26 are tracked individually.
  #23 is a distinct GetFEM MPI-2 root-cause investigation,
  not automatically solved by implementing new element families.
  Do not close #6 before its full physical/MPI tests pass.
- The tangential-velocity-controlled Saint-Venant 2D study for
  **LHMNC 2027 starts only after completing all gates 0–9.**
