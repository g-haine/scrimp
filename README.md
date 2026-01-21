# What is [SCRIMP](https://g-haine.github.io/scrimp/)?

**[SCRIMP](https://g-haine.github.io/scrimp/)** (Simulation and ContRol of Interactions in Multi-Physics) is a python package of methods and classes for the structure-preserving discretization and simulation of multi-physics models, using the formalism of port-Hamiltonian systems ([van der Schaft and Maschke (2002)](https://doi.org/10.1016/S0393-0440(01)00083-3)). 

**[SCRIMP](https://g-haine.github.io/scrimp/)** aims at speeding the coding process of the **Partitioned Finite Element Method** on a wide range of (multi-)physical systems ([Ferraro *et al.* (2024)](https://doi.org/10.1016/j.ifacol.2024.08.267)), and scrimp and save time!

The main objective of a **structure-preserving discretization** in the port-Hamiltonian formalism is to obtain a discrete version of the power balance satisfied by the Hamiltonian functional.

**See the website for more information and documentation: https://g-haine.github.io/scrimp/**

# Schema-driven configuration

Starting from this release, SCRIMP can be configured using declarative **YAML** or **JSON** schema files. The new `scrimp.io.schema_loader` module provides Pydantic models describing meshes, integration rules, FEM spaces and material parameters. The helper CLI exposes a minimal template library so that you can bootstrap a project without writing a single line of Python code:

```bash
python -m scrimp.io.schema_loader list domain
python -m scrimp.io.schema_loader generate rectangle --output rectangle.yml
```

In your simulation script you can now load the schema and pass it directly to ``scrimp.Domain`` and ``scrimp.FEM``:

```python
import scrimp as S
from scrimp.io.schema_loader import load_schema

schema = load_schema("rectangle.yml")
domain = S.Domain(schema.domain)
dpHS = S.DPHS("real")
dpHS.set_domain(domain)

velocity = S.FEM(schema.fem.fields[0])
dpHS.add_FEM(velocity)
```

The imperative API showcased in the tutorials is still available, but the schema-based workflow enables reproducible studies and simpler collaboration.

# How to install
The easiest way to install SCRIMP is to use a conda environment.

1. Install <a href="https://docs.anaconda.com/free/anaconda/install/index.html"> Anaconda</a>
2. Clone the git repository: ```git clone https://github.com/g-haine/scrimp```
3. Enter the folder: ```cd scrimp```
4. Create the conda environment:  ```conda env create --file /path/to/scrimp/scrimp.yml```
5. Activate the environment:  ```conda activate scrimp```
6. Finish with pip: ```pip install -e .```
