How to install
==============

Anaconda
--------

The easiest way to install **SCRIMP** is to use *conda environment*.

#. Install `Anaconda <https://www.anaconda.com/>`_
#. Clone the git repository :code:`git clone https://github.com/g-haine/scrimp`
#. Create the conda environment :code:`conda create env --file /path/to/scrimp/data/scrimp.yml`
#. Activate the environment :code:`conda activate scrimp`
#. Add the *path/to/scrimp* folder to this environment :code:`conda develop /path/to/scrimp/`

Now you can :code:`import scrimp` in your scripts when they are launched from the conda environment *scrimp*.

Tests
-----

You may test your installation by running avalaible examples.

Importation of an example is possible with, *e.g.* for the wave system, :code:`from scrimp.examples.wave import wave`.

Alternatively, you may run all examples with

.. code-block:: python
    :linenos:
    
    from scrimp import test_install
    test_install()

Code structure
--------------

**SCRIMP** is developped as a *package*: the :code:`__init__.py` file of the */path/to/scrimp/* folder is the root file. Each subdirectory is a sub-package of **SCRIMP**. Files are called *module* in this framework and may be called *via* the command **import**. For instance the module *linalg* gathering linear algebra functions of the *subpackage* utils can be imported with :code:`import scrimp.utils.linalg`.

