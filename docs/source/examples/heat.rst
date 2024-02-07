The heat equation
=================

**!TO DO!**

Modelling
---------

The so-called *heat equation* is driven by the first law of thermodynamics.

Let :math:`\Omega \subset \mathbb{R}^n` be a bounded open connected set, with mass density :math:`\rho(\mathbf{x})`, for all :math:`\mathbf{x} \in \Omega`, and :math:`\mathbf{n}` be the outwards unit normal at the boundary :math:`\partial\Omega`. We assume that:

- The domain :math:`\Omega` does not change over time: *i.e.* we work at constant volume in a solid
- No chemical reaction is to be found in the domain
- Dulong-Petit's model: internal energy is proportional to temperature

Let us denotes:

- :math:`u` the internal energy density
- :math:`\mathbf{J}_Q` the heat flux
- :math:`T` the local temperature
- :math:`C_V := \left( \frac{d u}{d T} \right)_V` the isochoric heat capacity

The first law of thermodynamics, stating that in an isolated system, the energy is preserved, reads:

.. math::

    \rho(\mathbf{x}) \partial_t u(t, \mathbf{x}) = - {\rm div} \left( \mathbf{J}_Q(t, \mathbf{x}) \right), \qquad \forall t \ge 0, \mathbf{x} \in \Omega.

Under Dulong-Petit's model, one has :math:`u = C_V T`, which leads to

.. math::

    \rho(\mathbf{x}) C_V(\mathbf{x}) \partial_t T(t, \mathbf{x}) = - {\rm div} \left( \mathbf{J}_Q(t, \mathbf{x}) \right), \qquad \forall t \ge 0, \mathbf{x} \in \Omega.

As constitutive relation, the classical Fourier's law is considered:

.. math::

    \mathbf{J}_Q(t, \mathbf{x}) = - \lambda(\mathbf{x}) \cdot \textbf{grad}\left( T(t, \mathbf{x}) \right), \qquad \forall t \ge 0, \mathbf{x} \in \Omega,

where :math:`\lambda` is the **tensor-valued** heat conductivity of the medium.

Hamiltonian
~~~~~~~~~~~



State variable
~~~~~~~~~~~~~~



Causality
~~~~~~~~~



Discretization
--------------



Code explanation
----------------



.. automodule:: scrimp.examples.heat
   :members:
   :undoc-members:
   :show-inheritance:
