# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/published/MCSS_2023/launch_all_simulations.py
- authors:          Ghislain Haine
- date:             22 sep. 2023
- brief:            launch all simulations and post-processings of the paper about dissipative shallow water equations
"""

# Enter `mpirun -n 4 python launch_all_simulations.py` to run on 4 cores
import time
from petsc4py import PETSc
comm = PETSc.COMM_WORLD
rank = comm.getRank()
from shallow_water import shallow_water
from PV_trace_experiment import trace_experiment
start = time.perf_counter()
for formulation in ["grad","div"]:
    for k in [3,2,1,0]:
        try:
            shallow_water(k,formulation)
        except:
            print(f"Simulation of experiment {k} using formulation {formulation} fails.")
        try:
            if rank==0:
                trace_experiment(k,formulation)
        except:
            print(f"Post-processing of experiment {k} using formulation {formulation} fails.")
    comm.barrier()
if rank==0:
    print(f"All simulations have been done in {time.perf_counter() - start:1.4g}s")
