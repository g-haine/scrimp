Examples
========

.. _examples:

The examples illustrate both the **scientific structure** of the
models and the **practical SCRIMP API**. Start with the
:doc:`first 1D vibrating-string tutorial <started>`: it
introduces energy variables, co-energies, outward boundary
tractions, discrete power balance and a real installation check
using MUMPS, serial SuperLU and MPI SuperLU_DIST.

The examples follow a progression in **port-Hamiltonian structure**,
rather than a strict ordering by spatial dimension or mathematical
complexity. After the introductory 1D string, the wave and Maxwell
models illustrate distributed linear systems. The finite-dimensional
mass--spring--damper and nonlinear pendulum then isolate the roles of
interconnection, dissipation and nonlinear constitutive relations,
without requiring spatial finite elements. The heat equation introduces
the essential algebraic constraints of the mixed diffusion formulation,
before the coupled heat--wave and nonlinear shallow-water models.

The individual pages link to the relevant work in our
:doc:`publications <biblio>`. Their original paper figures are
**not automatically reproduced** by short installation checks;
each scientific qualification must state its physical model,
solver, mesh (where applicable), accepted time steps and
energy/port-power evidence.

.. toctree::
   :maxdepth: 1
   :glob:
   
   examples/wave
   examples/wave_coenergy
   examples/maxwell_3d
   examples/mass_spring_damper
   examples/nonlinear_pendulum
   examples/heat
   examples/heat_wave
   examples/shallow_water
   
