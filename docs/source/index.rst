.. scrimp documentation master file, created by
   sphinx-quickstart on Sun Nov 27 14:16:43 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

##########
**SCRIMP**
##########

**Simulation and ContRol of Interactions in Multi-Physics**

.. image:: WorkFlow/workflow.png
    :width: 600px
    :height: 200px
    :alt: WorkFlow
    :align: center

***************
What is SCRIMP?
***************

**SCRIMP** (Simulation and ContRol of Interactions in Multi-Physics) is a python collection, *namely* a package, of *methods* and *classes* for the structure-preserving discretization and simulation of multi-physics models, using the formalism of port-Hamiltonian systems (`van der Schaft and Maschke (2002) <https://doi.org/10.1016/S0393-0440(01)00083-3>`_).

**SCRIMP** aims at speeding the coding process of the **Partitioned Finite Element Method** on a wide range of (multi-)physical systems (`Brugnoli *et al.* (2021) <https://doi.org/10.4236/jamp.2021.96088>`_), and scrimp and save time!

The documentation is `available in pdf <https://g-haine.github.io/scrimp/latex/scrimp.pdf>`_.

Table of Contents
=================

.. toctree::
    :maxdepth: 2
    
    install
    started
    examples
    scrimp
    References <biblio>
    ANR Impacts <https://impacts.ens2m.fr/>
    GitHub repo <https://github.com/g-haine/scrimp>
    PDF documentation <https://g-haine.github.io/scrimp/latex/scrimp.pdf>

************************
Port-Hamiltonian systems
************************

What are they?
==============

Let us sketch a rough portrait of port-Hamiltonian systems as they are considered in **SCRIMP**.

Port-Hamiltonian systems constitute a strongly structured class of control systems with collocated observation. It relies on a functional form :math:`\mathcal{H}` (the **Hamiltonian**), whose variables :math:`\alpha_i` are the **states** of the system. The **co-states** :math:`M_i e_i := \delta_{\alpha_i} \mathcal{H}` are defined as the variational derivative of the Hamiltonian with respect to the states, on the metric induced by the :math:`M_i` matrices.

The dynamics is provided *via* trajectories belonging in a **Dirac** structure, which can be represented by two matrices (of operators) :math:`M` symmetric and :math:`J` skew-symmetric as

.. math::

    M
    \begin{pmatrix} \frac{\rm d}{{\rm d}t} \alpha_1(t) \\ \vdots \\ \frac{\rm d}{{\rm d}t} \alpha_k(t) \\ f_R(t) \\ - y_{exp}(t) \\ u_{imp}(t) \end{pmatrix} 
    = J
    \begin{pmatrix} e_1(t) \\ \vdots \\ e_k(t) \\ e_R(t) \\ u_{exp}(t) \\ - y_{imp}(t) \end{pmatrix}

together with **constitutive relations**

.. math::

    M_i e_i(t) = \delta_{\alpha_i(t)} \mathcal{H}(\alpha_1(t), \cdots, \alpha_k(t))
    \qquad
    \mathcal{N}(t, f_R(t), e_R(t)) = 0
    
This structure allows to describe the evolution of the Hamiltonian along the trajectories

.. math::

    \frac{\rm d}{{\rm d}t} \mathcal{H}(\alpha_1(t), \cdots, \alpha_k(t)) = - e_R(t)^\top M_R f_R(t) + u_{exp}(t)^\top M_{exp} y_{exp}(t) + u_{imp}(t)^\top M_{imp} y_{imp}(t)

The first term of the right-hand side stands for a loss of *energy*, hence the name of *resistive port* for the couple :math:`(f_R,e_R)`. The other two terms stands for exchanges with the environment through the *control ports*. One is *explicit*, :math:`u_{exp}`, as a usual forcing term in the equations (its collocated output :math:`y_{exp}` plays no role in the dynamics). The other is *implicit*: :math:`u_{imp}` does not appear directly in the dynamics, and its collocated output :math:`y_{imp}` plays the role of the Lagrange multiplier imposing the value of :math:`u_{imp}`.

Each indexed matrix :math:`M_\ell` is the appropriate sub-matrix of :math:`M`.

A very important and useful fact is that the matrices :math:`M` and :math:`J` can depend on time and states!

The Partitioned Finite Element Method
=====================================

The main objective of a **structure-preserving discretization** in the port-Hamiltonian formalism is to obtain a discrete version of the power balance satisfied by the Hamiltonian functional.

A recent scheme, known as the **Partitioned Finite Element Method** (PFEM) (`Cardoso-Ribeiro *et al.* (2021) <https://doi.org/10.1093/imamci/dnaa038>`_), achieves this goal. 

The strategy follows three steps, inspired by the Mixed Finite Element Method for steady-state problem with homogeneous boundary condition

* write the weak form of the system;
* integrate by parts a **partition** of the state (such that *the control appears*); and
* project on finite element spaces.

*****************
Coding philosophy
*****************

**SCRIMP** assumes that the final user is not familiar with numerical simulations. The aim is to facilitate the first step from modelisation to simulation by sticking as much as possible to the port-Hamiltonian framework, getting rid of coding issues.

As such, these simplifications naturally imply a lack of optimization of the code. Nevertheless, the syntax of **SCRIMP** try to let confirmed users to reach finer tuning in order to perform more sophisticated simulations.

A basic usage of **SCRIMP** consists in a script with the following steps:

- Define a *domain*
- Define at least one *state*. And of course, its *co-state*, in order to get a *dynamical port*
- Define a Finite Element Method on this port: give at least an order, at first glance, default values are sufficient
- Define *algebraic ports* (not mandatory) and its FEM
- Define *control ports* (not mandatory) and its FEM
- Write down the forms on the *flow side* of the Dirac structure, *i.e.* the **brick** defining the matrix :math:`M`
- Write down the forms on the *effort side* of the Dirac structure, *i.e.* th **brick** defining the matrix :math:`J`
- Write down all the forms defining the *constitutive relations*, always with **bricks**
- Set up time scheme options: again, at first glance, default values are sufficient
- Solve
- Plot
- Export

We try to eliminate as much as possible the *computer-side* of the simulations, by following the port-Hamiltonian vocabulary, always by keeping the possibility of fine tuning available.

*******
Credits
*******

Development
===========

Please report bug at: ghislain.haine@isae.fr, Giuseppe.Ferrarro@isae-supaero.fr

Current developers: Antoine Bendhimerad-Hohl, Giuseppe Ferraro, Michel Fournié, Ghislain Haine

Past: Andrea Brugnoli, Melvin Chopin, Florian Monteghetti, Anass Serhani, Xavier Vasseur

**Please read the** `LICENSE <https://github.com/g-haine/scrimp/blob/master/LICENSE>`_

Funding
=======

- ANR Project `IMPACTS <https://impacts.ens2m.fr/>`_ -- IMplicit Port-hAmiltonian ConTrol Systems 

- AID School Project `FAMAS <https://www.defense.gouv.fr/aid>`_ -- Fast & Accurate MAxwell Solver 

- ANR-DFG Project `INFIDHEM <http://websites.isae.fr/infidhem>`_ -- INterconnected inFinite-Dimensional systems for HEterogeneous Media

Third-party
===========

The two **main** libraries used as core for SCRIMP are:

- `GetFEM <https://getfem.org/>`_ -- An open-source finite element library

- `PETSc <https://petsc.org/release/>`_ -- The Portable, Extensible Toolkit for Scientific Computation

Meshing is facilitated using (although not mandatory)
`GMSH <https://gmsh.info/>`_ -- A three-dimensional finite element mesh generator

Post-processing visualization is encouraged via
`ParaView <https://www.paraview.org/>`_ -- Post-processing visualization engine

and finally, SCRIMP also needs for some routines

- `matplotlib <https://matplotlib.org/>`_ -- Visualization with Python

- `numpy <https://numpy.org/>`_ -- A well-known package for scientific computing

One of our choice for IDE is
`Spyder <https://www.spyder-ide.org/>`_ -- A scientific Python development environment

How to cite SCRIMP?
===================

Brugnoli, Andrea and Haine, Ghislain and Serhani, Anass and Vasseur, Xavier.
*Numerical Approximation of Port-Hamiltonian Systems for Hyperbolic or Parabolic PDEs with Boundary Control.*
(2021) **Journal of Applied Mathematics and Physics**, 09 (06). 1278-1321.

.. code-block:: latex

    @article{Brugnoli2021,
    author = {Brugnoli, Andrea and Haine, Ghislain and Serhani, Anass and Vasseur, Xavier},
    title = { {Numerical Approximation of Port-Hamiltonian Systems for Hyperbolic or Parabolic PDEs with Boundary Control} },
    journal = {Journal of Applied Mathematics and Physics},
    volume = {09},
    issue = {06},
    pages = {1278--1321},
    year = {2021}
    }

