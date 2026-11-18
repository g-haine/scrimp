# What is [SCRIMP](https://g-haine.github.io/scrimp/)?

**[SCRIMP](https://g-haine.github.io/scrimp/)** (Simulation and ContRol of Interactions in Multi-Physics) is a python package of methods and classes for the structure-preserving discretization and simulation of multi-physics models, using the formalism of port-Hamiltonian systems ([van der Schaft and Maschke (2002)](https://doi.org/10.1016/S0393-0440(01)00083-3)). 

**[SCRIMP](https://g-haine.github.io/scrimp/)** aims at speeding the coding process of the **Partitioned Finite Element Method** on a wide range of (multi-)physical systems ([Ferraro *et al.* (2024)](https://doi.org/10.1016/j.ifacol.2024.08.267)), and scrimp and save time!

The main objective of a **structure-preserving discretization** in the port-Hamiltonian formalism is to obtain a discrete version of the power balance satisfied by the Hamiltonian functional.

**See the website for more information and documentation: https://g-haine.github.io/scrimp/**

# How to install
The easiest way to install SCRIMP is to use a conda environment.

1. Install <a href="https://docs.anaconda.com/free/anaconda/install/index.html"> Anaconda</a>
2. Clone the git repository: ```git clone https://github.com/g-haine/scrimp```
3. Enter the folder: ```cd scrimp```
4. Create the conda environment:  ```conda env create --file /path/to/scrimp/scrimp.yml```
5. Activate the environment:  ```conda activate scrimp```
6. Finish with pip: ```pip install -e .```

## Experimental GetFEM 5.5 / NumPy 2 experimental MPI stack

For a from-source **parallel GetFEM 5.5** build (Conda or Docker), MPI
dependency checks, 2D acoustic-wave scaling scripts, and an opt-in
GHCR image publication workflow, see
[installation and MPI build instructions](docs/source/install.rst).

The conventional `scrimp.yml` setup remains unchanged until
GetFEM 5.5 + NumPy 2 and physical parallel benchmarks pass in CI.
No experimental container is advertised as validated or released yet.

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
and [v1.3.0 release gates and open questions](docs/V130_DEVELOPMENT_ROADMAP.md).
