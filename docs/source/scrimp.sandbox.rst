scrimp.sandbox
==============

The ``sandbox`` folder is the **recommanded** folder to work in.

It already contains the example of the 1D wave equation presented `here <started>`_.


Experimental PDE--ODE coupling
-------------------------------

The sandbox also contains an experimental coupled shallow-water /
moving-tank system at ``sandbox/sloshing_tank.py`` and a convergence
driver at ``sandbox/tank_convergence.py``.

These models demonstrate global scalar ports coupled to distributed
states. The wall boundary conditions and spatial-field convergence
remain under investigation. The tank example is intentionally **not**
a release-validated example yet. See ``sandbox/README_lumped_tank.md``.
