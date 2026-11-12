How to install
==============

Anaconda
--------

The easiest way to install SCRIMP is to use a conda environment.

1. Install `Anaconda <https://docs.anaconda.com/free/anaconda/install/index.html>`_
2. Clone the git repository: ``git clone https://github.com/g-haine/scrimp``
3. Enter the folder: ``cd scrimp``
4. Create the conda environment:  ``conda env create --file /path/to/scrimp/scrimp.yml``
5. Activate the environment:  ``conda activate scrimp``
6. Install the local package in editable mode: ``python -m pip install -e .``

Tests
-----

You can test the installation and the global scalar port extension with:

.. code-block:: bash

    python -m unittest -v tests.test_lumped_ports
    python examples/mass_spring_damper.py
    python examples/nonlinear_pendulum.py

The last two commands also require a graphical Matplotlib backend for plots.
The tests and examples need working GetFEM, PETSc, SLEPc and Gmsh installations.

Code structure
--------------

**SCRIMP** is developped as a *package*: the :code:`__init__.py` file of the */path/to/scrimp/* folder is the root file. Each subdirectory is a sub-package of **SCRIMP**. Files are called *module* in this framework and may be called *via* the command **import**. For instance the module *linalg* gathering linear algebra functions of the *subpackage* utils can be imported with :code:`import scrimp.utils.linalg`.

Documentation
-------------

You can find this documentation `here <https://g-haine.github.io/scrimp/latex/scrimp.pdf>`_.

It is automatically built upon the code comments using sphinx.

See `Sphinx <https://www.sphinx-doc.org/>`_ for further informations.

