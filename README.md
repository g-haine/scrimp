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

# How to install
The easiest way to install SCRIMP is to use a conda environment.

1. Install <a href="https://docs.anaconda.com/free/anaconda/install/index.html"> Anaconda</a>
2. Clone the git repository: ```git clone https://github.com/g-haine/scrimp```
3. Enter the folder: ```cd scrimp```
4. Create the conda environment:  ```conda env create --file /path/to/scrimp/scrimp.yml```
5. Activate the environment:  ```conda activate scrimp```
6. Finish with pip: ```pip install -e .```

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
`outputs/png` or `outputs/pv` siblings. The eight original examples
completed a real nine-case native I/O smoke on 27 November 2026,
including Heat--Wave MPI2. The revised *clean* layout also passed a
targeted native Heat--Wave MPI2 run. Full nine-case requalification
**after** that directory fix remains a separate release gate.

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
