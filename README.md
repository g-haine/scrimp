#  SCRIMP - Simulation and ContRol of Interactions in Multi-Physics

Git: https://github.com/g-haine/scrimp

Website: https://g-haine.github.io/scrimp/
___

## What is SCRIMP?

**SCRIMP** is a python collection of methods and classes for the structure-preserving 
discretization and simulation of multi-physics models, using the 
formalism of port-Hamiltonian systems. It intends to speed up the coding 
processus, and scrimp time!

**Notebooks** ([Jupyter](https://jupyter.org/)) are provided as examples.
___

## Development

Please report bug at: ghislain.haine@isae.fr

Current developers: Florian Monteghetti, Ghislain Haine, Antoine Bendhimerad-Hohl, Melvin Chopin

Past: Andrea Brugnoli, Anass Serhani, Xavier Vasseur

**Please read the LICENSE**

### Founding

ANR-DFG Project **INFIDHEM**: http://websites.isae.fr/infidhem

AID School Project **FAMAS**: https://www.defense.gouv.fr/aid

ANR Project **IMPACTS**: https://impacts.ens2m.fr/

### Third-party

[numpy](https://numpy.org/): A well-known package for scientific computing

[FEniCS](https://fenicsproject.org/): A finite elements library

[multiphenics](https://github.com/multiphenics/multiphenics): An easy way to prototype multiphysics problems using FEniCS

[PETSc](https://petsc.org/release/): The Portable, Extensible Toolkit for Scientific Computation

[GMSH](https://gmsh.info/): A three-dimensional finite element mesh generator

[matplotlib](https://matplotlib.org/): Visualization with Python

[Spyder](https://www.spyder-ide.org/): A scientific Python development environment

[Jupyter](https://jupyter.org/): A software for interactive computing across all programming languages

### How to cite SCRIMP?
    
Brugnoli, Andrea and Haine, Ghislain and Serhani, Anass and Vasseur, Xavier.

*Numerical Approximation of Port-Hamiltonian Systems for Hyperbolic or Parabolic PDEs with Boundary Control.*

(2021) **Journal of Applied Mathematics and Physics**, 09 (06). 1278-1321.
___

# How to install SCRIMP:

- Install **GMSH** client (not from the package manager! Outdated version!): https://gmsh.info/

- Install **conda**: https://docs.conda.io

- Install **mutliphenics**: https://github.com/multiphenics/multiphenics in your /path/to/scrimp/src/modules
```
mkdir /path/to/scrimp/src/modules/multiphenics
cd /path/to/scrimp/src/modules/multiphenics
git clone https://github.com/multiphenics/multiphenics.git
python setup.py install
```

- Create the *scrimp* environment:
```
conda install --file "/path/to/your/scrimp/directory/data/scrimp.yml"
```
___

# How to get started using SCRIMP:

- Activate the environment:

```
conda activate scrimp
```

- Import the module by starting your script.py with:

```
rootdir = "/path/to/your/scrimp/directory/"
import sys
sys.path.append(rootdir)
import scrimp
scrimp.setup_path(rootdir)
```
___

# How to test my installation:

- Run provided Notebooks.

** In your scrimp folder (do not forget to activate the scrimp environment before):

```
jupyter notebook &
```

** In [Jupyter](https://jupyter.org/), choose your Notebook.

- First notebook to test: 2D/wave.ipynb

