Graphical user interface (experimental)
=======================================

The graphical interface is **not currently part of the qualified
v1.3.0 user workflow**. Older SCRIMP versions included a graphical
script generator, but its complete GUI-to-source-to-native-simulation
path has not yet passed the same release acceptance gates as the
Python examples. In particular, the published 1.3.0 interface must
not promise that every GUI-generated model and time scheme has
already been validated.

For a reproducible introduction, begin with
:doc:`started` and the directly executable
:doc:`examples`. The examples demonstrate the same public
port-Hamiltonian constructs using real PETSc/GetFEM solvers
and generate their own plots and data under
``outputs/<system_name>/``.

The GUI may be explored as an experimental feature only,
with its generated source reviewed before numerical use.
Independent validation of GUI startup, generated scripts,
error reporting and MPI output handling remains a
separate release-scope decision; it is **not** silently
treated as covered by the eight canonical example smokes.
