SCRIMP
======

# What is [SCRIMP](https://g-haine.github.io/scrimp/)?

**[SCRIMP](https://g-haine.github.io/scrimp/)** (Simulation and ContRol of Interactions in Multi-Physics) is a python package of methods and classes for the structure-preserving discretization and simulation of multi-physics models, using the formalism of port-Hamiltonian systems ([van der Schaft and Maschke (2002)](https://doi.org/10.1016/S0393-0440(01)00083-3)). 

**[SCRIMP](https://g-haine.github.io/scrimp/)** aims at speeding the coding process of the **Partitioned Finite Element Method** on a wide range of (multi-)physical systems ([Brugnoli *et al.* (2021)](https://doi.org/10.4236/jamp.2021.96088)), and scrimp and save time!

The main objective of a **structure-preserving discretization** in the port-Hamiltonian formalism is to obtain a discrete version of the power balance satisfied by the Hamiltonian functional.

**See the website for more information: https://g-haine.github.io/scrimp/**

# How to install
The easiest way to install SCRIMP is to use a conda environment.

1. Install <a href="https://docs.anaconda.com/free/anaconda/install/index.html"> Anaconda</a>
2. Clone the git repository: ```git clone https://github.com/g-haine/scrimp```
3. Pull the branch dev_v2: ```git fetch origin dev_v2:dev_v2``` 
4. Checkout dev_v2 branch: ```git checkout dev_v2```
5. Create the conda environment:  ```conda env create --file /path/to/scrimp/scrimp.yml```
6. Activate the environment:  ```conda activate scrimp```
