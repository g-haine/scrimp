# About SCRIMP

The main objective of a **structure-preserving discretization** in the port-Hamiltonian formalism ([van der Schaft and Maschke (2002)](https://doi.org/10.1016/S0393-0440(01)00083-3)) is to obtain a discrete version of the power balance satisfied by the Hamiltonian functional.

A recent scheme, known as the **Partitioned Finite Element Method** (PFEM) ([Cardoso-Ribeiro *et al.* (2021)](https://doi.org/10.1093/imamci/dnaa038)), achieves this goal.

**SCRIMP** (Simulation and ContRol of Interactions in Multi-Physics) aims at speeding the coding process of PFEM on a wide range of (multi-)physical systems ([Brugnoli *et al.* (2021)](https://doi.org/10.4236/jamp.2021.96088)).

## The Partitioned Finite Element Method

The strategy follows three steps, inspired by the Mixed Finite Element Method for steady-state problem with homogeneous boundary condition
* write the weak form of the system;
* integrate by parts a **partition** of the state (such that $u_\partial$ appears); and
* project on finite element spaces.
