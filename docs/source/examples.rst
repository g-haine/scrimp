Examples
========

.. _examples:

The examples illustrate both the **scientific structure** of the
models and the **practical SCRIMP API**. Start with the
:doc:`first 1D vibrating-string tutorial <started>`: it
introduces energy variables, co-energies, outward boundary
tractions, discrete power balance and a real installation check
using MUMPS, serial SuperLU and MPI SuperLU_DIST.

The following models are presented in progressively richer spatial
and multiphysical settings, with a link to the relevant work in
our :doc:`publications <biblio>`. Their original paper figures
are **not automatically reproduced** by the short installation
smokes; each validated test must state its physical model, solver,
mesh, accepted time steps and energy/port-power evidence.

.. toctree::
   :maxdepth: 1
   :glob:
   
   examples/wave
   examples/heat
   examples/wave_coenergy
   examples/heat_wave
   examples/shallow_water
   examples/mass_spring_damper
   examples/nonlinear_pendulum
   
