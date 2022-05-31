#  SCRIMP - Simulation and ContRol of Interactions in Multi-Physics

__________________________________________________________________________________
__________________________________________________________________________________
   ____ _________________________ _______   |      ___ ___  ____   ____   ____
  |      ___  ___  __ __  _   _  |    |     |         |    |      /    \ |
  |___  |    |___|   |   | \ / | |____|     |         |    |____  |____| |___
      | |    |  \    |   |  V  | |          |         |         | |    | |
______|=|___=|   \=__|__=|     |=|          |      ___|___  ____| |    | |____
                               SCRIMP 0.5   |             == SUPAERO ==
__________________________________________________________________________________
__________________________________________________________________________________


            Simulation and ContRol of Interactions in Multi-Physics
            Git: https://github.com/g-haine/scrimp
__________________________________________________________________________________

           ANR-DFG Project INFIDHEM: http://websites.isae.fr/infidhem
                           AID School Project FAMAS
              ANR Project IMPACTS: https://impacts.ens2m.fr/
__________________________________________________________________________________

    Please report bug at: florian.monteghetti@isae.fr, ghislain.haine@isae.fr

    Current developers: Florian Monteghetti, Ghislain Haine,
                        Antoine Bendhimerad-Hohl, Melvin Chopin
    Past: Andrea Brugnoli, Anass Serhani, Xavier Vasseur

    Please read the LICENSE.txt
__________________________________________________________________________________

    SCRIMP is a python collection of methods and classes for the structure-
preserving discretization and simulation of multi-physics models, using the 
formalism of port-Hamiltonian systems. It intends to speed up the coding 
processus.
__________________________________________________________________________________

    How to cite SCRIMP?
    
    Brugnoli, Andrea and Haine, Ghislain and Serhani, Anass and Vasseur, Xavier.
    Numerical Approximation of Port-Hamiltonian Systems for Hyperbolic or 
        Parabolic PDEs with Boundary Control.
    (2021) Journal of Applied Mathematics and Physics, 09 (06). 1278-1321.

__________________________________________________________________________________

# How to install SCRIMP:

- After installing conda (https://docs.conda.io), install the scrimp environment:

    conda install --file "/path/to/your/scrimp/directory/data/scrimp.yml"

- Install mutliphenics: https://github.com/multiphenics/multiphenics

# How to get started using SCRIMP:

- Activate the environment:

    conda activate scrimp

- Import the module by starting your script.py by:

    rootdir = "/path/to/your/scrimp/directory/"
    import sys
    sys.path.append(rootdir)
    import scrimp
    scrimp.setup_path(rootdir)

