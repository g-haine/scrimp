How to install
==============

Anaconda
--------

The easiest way to install SCRIMP is to use a conda environment.

1. Install `Anaconda <https://docs.anaconda.com/free/anaconda/install/index.html>`_
2. Clone the git repository: ```git clone https://github.com/g-haine/scrimp```
3. Enter the folder: ```cd scrimp```
4. Create the conda environment:  ```conda env create --file /path/to/scrimp/scrimp.yml```
5. Activate the environment:  ```conda activate scrimp```
6. Add scrimp to the PATH: ```conda develop /path/to/scrimp/```
7. Finish with pip: ```pip install -e .```

Tests
-----

You may test your installation by running avalaible examples in the ```examples``` folder.

Code structure
--------------

**SCRIMP** is developped as a *package*: the :code:`__init__.py` file of the */path/to/scrimp/* folder is the root file. Each subdirectory is a sub-package of **SCRIMP**. Files are called *module* in this framework and may be called *via* the command **import**. For instance the module *linalg* gathering linear algebra functions of the *subpackage* utils can be imported with :code:`import scrimp.utils.linalg`.

Documentation
-------------

You can build the documentation locally by running `sphinx-build` in the docs folder.

See `Sphinx <https://www.sphinx-doc.org/>`_ for further informations.

