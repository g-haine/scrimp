# What is SCRIMP?

**SCRIMP** (Simulation and ContRol of Interactions in Multi-Physics) is a python collection of methods and classes for the structure-preserving discretization and simulation of multi-physics models, using the formalism of port-Hamiltonian systems. 

**SCRIMP** aims at speeding the coding process of the **Partitioned Finite Element Method** on a wide range of (multi-)physical systems ([Brugnoli *et al.* (2021)](https://doi.org/10.4236/jamp.2021.96088)), and scrimp time!

The main objective of a **structure-preserving discretization** in the port-Hamiltonian formalism ([van der Schaft and Maschke (2002)](https://doi.org/10.1016/S0393-0440(01)00083-3)) is to obtain a discrete version of the power balance satisfied by the Hamiltonian functional.

# The Partitioned Finite Element Method

A recent scheme, known as the **Partitioned Finite Element Method** (PFEM) ([Cardoso-Ribeiro *et al.* (2021)](https://doi.org/10.1093/imamci/dnaa038)), achieves this goal. 

The strategy follows three steps, inspired by the Mixed Finite Element Method for steady-state problem with homogeneous boundary condition
* write the weak form of the system;
* integrate by parts a **partition** of the state (such that the control *appears*); and
* project on finite element spaces.

**Notebooks** ([Jupyter](https://jupyter.org/)) are provided as examples.

___

# Development

Please report bug at: ghislain.haine@isae.fr

Current developers: Florian Monteghetti, Ghislain Haine, Antoine Bendhimerad-Hohl, Melvin Chopin

Past: Andrea Brugnoli, Anass Serhani, Xavier Vasseur

**Please read the LICENSE**

# Founding

ANR-DFG Project **INFIDHEM**: http://websites.isae.fr/infidhem

AID School Project **FAMAS**: https://www.defense.gouv.fr/aid

ANR Project **IMPACTS**: https://impacts.ens2m.fr/

# Third-party

[numpy](https://numpy.org/): A well-known package for scientific computing

[FEniCS](https://fenicsproject.org/): A finite elements library

[multiphenics](https://github.com/multiphenics/multiphenics): An easy way to prototype multiphysics problems using FEniCS

[PETSc](https://petsc.org/release/): The Portable, Extensible Toolkit for Scientific Computation

[GMSH](https://gmsh.info/): A three-dimensional finite element mesh generator

[matplotlib](https://matplotlib.org/): Visualization with Python

[Spyder](https://www.spyder-ide.org/): A scientific Python development environment

[Jupyter](https://jupyter.org/): A software for interactive computing across all programming languages

___

# How to cite SCRIMP?
    
Brugnoli, Andrea and Haine, Ghislain and Serhani, Anass and Vasseur, Xavier.

*Numerical Approximation of Port-Hamiltonian Systems for Hyperbolic or Parabolic PDEs with Boundary Control.*

(2021) **Journal of Applied Mathematics and Physics**, 09 (06). 1278-1321.

___
