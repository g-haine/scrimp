#  SCRIMP - Simulation and ContRol of Interactions in Multi-Physics

Git: https://github.com/g-haine/scrimp
___

## What is SCRIMP?

**SCRIMP** is a python collection of methods and classes for the structure-preserving 
discretization and simulation of multi-physics models, using the 
formalism of port-Hamiltonian systems. It intends to speed up the coding 
processus.

**Notebooks** (Jupyter) are provided as examples.

## Development

Please report bug at: ghislain.haine@isae.fr

Current developers: Florian Monteghetti, Ghislain Haine, Antoine Bendhimerad-Hohl, Melvin Chopin

Past: Andrea Brugnoli, Anass Serhani, Xavier Vasseur

**Please read the LICENSE**

## Founding

ANR-DFG Project **INFIDHEM**: http://websites.isae.fr/infidhem

AID School Project **FAMAS**

ANR Project **IMPACTS**: https://impacts.ens2m.fr/

## How to cite SCRIMP?
    
Brugnoli, Andrea and Haine, Ghislain and Serhani, Anass and Vasseur, Xavier.

*Numerical Approximation of Port-Hamiltonian Systems for Hyperbolic or Parabolic PDEs with Boundary Control.*

(2021) **Journal of Applied Mathematics and Physics**, 09 (06). 1278-1321.
___

# How to install SCRIMP:

- Install **conda**: https://docs.conda.io

- Install **mutliphenics**: https://github.com/multiphenics/multiphenics

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

