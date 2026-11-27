# What is [SCRIMP](https://g-haine.github.io/scrimp/)?

**[SCRIMP](https://g-haine.github.io/scrimp/)** (Simulation and ContRol of Interactions in Multi-Physics) is a python package of methods and classes for the structure-preserving discretization and simulation of multi-physics models, using the formalism of port-Hamiltonian systems ([van der Schaft and Maschke (2002)](https://doi.org/10.1016/S0393-0440(01)00083-3)). 

**[SCRIMP](https://g-haine.github.io/scrimp/)** aims at speeding the coding process of the **Partitioned Finite Element Method** on a wide range of (multi-)physical systems ([Ferraro *et al.* (2024)](https://doi.org/10.1016/j.ifacol.2024.08.267)), and scrimp and save time!

The main objective of a **structure-preserving discretization** in the port-Hamiltonian formalism is to obtain a discrete version of the power balance satisfied by the Hamiltonian functional.

**See the website for more information and documentation: https://g-haine.github.io/scrimp/**

## Important: GWFL matrices are not NumPy matrices

SCRIMP is a **Python** package, but the weak-form strings passed to
GetFEM use its **Generic Weak-Form Language (GWFL)**. In nested
GWFL matrices, inner lists are **columns**, not the rows one expects
from Python/NumPy. Confusing these conventions can reverse the
direction of a physical boundary control.

See the [SCRIMP matrix notation and sign-check guide](docs/source/getfem_gwfl.rst)
and the [official GetFEM GWFL documentation ("Explicit matrices")](https://getfem.org/userdoc/gasm_high.html#explicit-matrices).

# Installation — Docker first for v1.3.0

**The recommended v1.3.0 installation is the precompiled Docker image**, not a
host-side Conda environment. Docker supplies SCRIMP with its tested Linux
GetFEM 5.5, complex PETSc/SLEPc, MUMPS and MPI stack. User simulations and
outputs remain in a **separate local project**, never in the image or the
SCRIMP Git repository.

**Release status (27 November 2026):** a **pre-release** Docker image
`ghcr.io/g-haine/scrimp:1.3.0-dist-rc3` has been published and qualified
on a fresh GitHub-hosted Ubuntu 24.04/AMD64 runner **by immutable digest**:
`sha256:0cbbce6ad42b23552aff9d12555342f3864b50522054fe039c161a7145b2732a`.
The independent ZIP/oscillator/MPI-startup and genuine Heat--Wave MPI1/MPI2
Hamiltonian parity gates passed. This tag may have restricted GHCR visibility.
**The stable `ghcr.io/g-haine/scrimp:1.3.0` image, official project ZIP and
GitHub Release are NOT yet published.** Do not confuse the pre-release
candidate with the stable release.

**Once the public v1.3.0 release and image are qualified and published:**

1. Install Docker, plus VS Code and its **Dev Containers** extension if desired.
2. Download the **scrimp-project-template.zip** asset from the matching GitHub
   Release and extract it in your own simulations directory. The released ZIP
   will refer to the verified Docker image explicitly.
3. Pull the image and open the extracted project in VS Code:

   ```bash
   docker pull ghcr.io/g-haine/scrimp:1.3.0
   ```

4. Choose **Dev Containers: Reopen in Container** and execute:

   ```bash
   python simulations/first_simulation.py
   mpiexec -n 2 python simulations/mpi_check.py
   ```

The first command runs a real port-Hamiltonian oscillator and writes
`outputs/oscillator.csv` **on your computer**; the second checks MPI startup
but is not a PDE parity test. A terminal-only Docker workflow is also
supported, without VS Code. Linux/Ubuntu is the qualified baseline; Windows
WSL2 and macOS are best effort, and Apple Silicon needs an explicitly
compatible image/platform.

For a locally compiled image, release-candidate preparation, host volume
permissions and the terminal workflow, see the
[installation guide](docs/source/install.rst),
[first project guide](docs/source/vscode.rst) and
[Docker distribution gates](docs/V130_DOCKER_FIRST_RELEASE_GATE.md).

## Alternative: native Conda (existing installations and advanced users)

Conda remains supported as an **alternative**, especially for development and
users with a working existing native environment; it is not the v1.3.0 default.

1. Install [Conda](https://docs.conda.io/) or Micromamba.
2. Clone the repository: `git clone https://github.com/g-haine/scrimp`.
3. Enter it: `cd scrimp`.
4. Create the environment: `conda env create -f scrimp.yml`.
5. Activate: `conda activate scrimp`.
6. Install the package: `python -m pip install -e .`.

The historical `scrimp.yml` remains separate from the GetFEM 5.5/NumPy 2
Docker candidate; this change does not silently upgrade host dependencies.

## v1.3.0 release preparation: conforming interfaces and MPI ParaView

**On the `release/1.3.0-preparation` branch only**, SCRIMP now
automatically resolves conforming internal interfaces and, for MPI
with complete identical mesh/FE-state replicas, uses **rank 0 as the
single writer** of shared ParaView `.pvd`/`.vtu` outputs.
Both implementations are part of the default code paths rather
than opt-in pilots; **there is no published v1.3.0 release yet**.

An ambiguous/nonconforming interface or a ParaView model with
incomplete or inconsistent per-rank states **fails explicitly**,
rather than silently guessing orientation or exporting partial data.
The documented native qualification covers the 2D concentric
Heat-Wave model in both subdomain orientations at MPI 1 and 2;
other elements, geometries and truly distributed-only VTK data
are not implied.

To temporarily restore the previous behavior:

```bash
export SCRIMP_AUTO_INTERFACE=0
export SCRIMP_PV_EXPORT_MODE=legacy
```

**Warning:** the historical MPI `legacy` export may let
several ranks write the same file concurrently. For complete
per-rank replicas the new default is safer.

See the [dated production-path decision, successful native evidence,
corrected failures and explicit scope limitations](
docs/V130_INTERFACE_PV_PRODUCTION_DECISION_20261124.md).
The native **default-path smoke has PASSED**, without experimental
feature flags, at commit `241f58fb` in the pinned v1.3.0 RC3 image:
71 source tests, four MPI1/MPI2 genuine Heat-Wave DAE/CN runs,
energy/port-power balances and complete shared ParaView output.
The reproducible check is `scripts/run-v130-production-defaults.sh`;
evidence is `v130-auto-interface-MPI-CN-241f58fb-YEJXJhSK.zip`.
This proves only the documented Heat-Wave test scope; it is not
a tagged or published v1.3.0 release.

## v1.3.0 examples, energy/power and output files (pre-release)

The eight canonical Python examples now produce **their own real
GetFEM/PETSc plots and numerical data**, rather than requiring users
to invoke sandbox acceptance probes to render figures. The v1.3
user output convention is

```text
working_directory/
└── outputs/
    └── system_name/
        ├── png/     svg/     gif/
        ├── mesh/    log/
        └── npz/     json/
```

There is **one `outputs/` root per working directory**; examples
do not create generic `outputs/log`, `outputs/mesh`,
`outputs/png` or `outputs/pv` siblings. All eight canonical examples (eight MPI1 configurations and
Heat--Wave MPI2 as a ninth case) completed the **revised clean-layout**
native smoke on 27 November 2026, source `20d79d0d`:
`PASS_V130_ALL_NINE_REAL_USER_EXAMPLE_OUTPUT_CASES`, with
exactly one `outputs/` root, no legacy sibling folders,
genuine PVD/VTU time data and valid ZIP CRC. This
closes the canonical example output-layout gate.

`DPHS.plot_Hamiltonian()` exports actual total/partial energies
and raw port powers in PNG/SVG/NPZ. The helper
`compute_powers()` uses genuine CN step powers only when
**every accepted CN state is saved**; for decimated CN, BDF,
BDF2 and ARKIMEX it reports true saved-state nodal powers and an
**approximate trapezoidal work**, not hidden solver-stage work.
A generic plotted balance is not proof of energy conservation or
passivity for an arbitrary model.

See the [energy and time integration guide](docs/source/energy_and_power.rst),
[example output contract](docs/V130_EXAMPLES_OUTPUT_CONTRACT.md), and
[release-readiness status](docs/V130_RELEASE_READINESS_2026-11-27.md).
The documentation workflow builds HTML, LaTeX and PDF on pushes;
a maintained local reproducer is
`bash scripts/run-v130-docs-build.sh` (requires Sphinx and LaTeX).

## v1.3.0 preview: start your own SCRIMP project in VS Code

The upcoming v1.3.0 user workflow separates **installation** from
**simulations**. The goal is to choose a precompiled Docker image,
compile Docker locally, or use native Conda — and then work in the
**same independent VS Code project**. Your scripts, geometries and
results stay on your computer; they are never baked into Docker.

A reusable starter lives in [project-template](project-template/).
During development, create a project from a temporary SCRIMP checkout:

```bash
python3 scripts/create-scrimp-project.py "$HOME/MySCRIMPProject"
```

Open that new directory with VS Code and use **Dev Containers:
Reopen in Container** when using Docker. The template provides
serial Python debugging, an MPI startup task and a genuine short
mass–spring–damper simulation that saves
`outputs/oscillator.csv` locally. For native Conda, select the
native interpreter instead of reopening in a container.

See the [complete first-project / Dev Containers tutorial](docs/source/vscode.rst).
The future precompiled `ghcr.io/g-haine/scrimp:1.3.0` image
**is not published or qualified yet**, and the template is still
under pre-release validation. No SCRIMP source checkout should be
required of end users after the release project's ZIP asset ships.

## GetFEM 5.5 / NumPy 2 native MPI candidate (v1.3.0)

For a from-source **parallel GetFEM 5.5** build (Conda or Docker),
MPI dependency checks, actual wave scaling experiments and
an opt-in GHCR image publication workflow, see
[installation and MPI build instructions](docs/source/install.rst).

The local `scrimp-mpi:v1.3.0-rc3` image has passed the specified
native GetFEM/PETSc physical and MPI tests on Ubuntu. It is
**not** a publicly released image, nor evidence of a completed
clean-room installation or Windows/macOS qualification. The
conventional `scrimp.yml` environment remains unchanged.

## v1.3.0 development — PETSc/MUMPS performance

The development branch introduces an exact `on_change` MUMPS reuse policy
by default: the symbolic factor is retained while the sparse graph stays
unchanged, and the numerical factor is retained only when the actual PETSc
coefficients have not changed. On MPI runs, MUMPS parallel symbolic
analysis is enabled by default, subject to explicit solver options.

The v1.3.0 code is **not yet released**. The qualification compares the
original v1.2.0 wave example in serial with the current code in serial
and on 2/4/6 ranks, including actual full vectors and per-stage timings:

```bash
bash sandbox/run_v130_wave_benchmark.sh 0.04 3
```

See [solver reuse and benchmark methodology](docs/source/performance.rst)
and [v1.3.0 release plan and acceptance gates](docs/V130_RELEASE_PLAN.md),
with technical investigations and known limitations in the
[development history](docs/V130_DEVELOPMENT_ROADMAP.md).
