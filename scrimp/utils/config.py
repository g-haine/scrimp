# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             utils/config.py
- authors:          Ghislain Haine
- date:             23 jun. 2023
- brief:            functions to configure SCRIMP
"""

import os
import logging
import getfem as gf

def set_paths(path=None)
    global module_path, outputs_path
    
    if rank==0:
        if path is None:
            module_path = os.getcwd()
        else:
            module_path = path
        outputs_path = os.path.join(module_path, "outputs")
	for path in [os.path.join("outputs", "log"), os.path.join("outputs", "png"), os.path.join("outputs", "mesh"), os.path.join("outputs", "pv")]:
	    composed_path = os.path.join(module_path, path)
	    if not os.path.isdir(composed_path):
		os.makedirs(composed_path)

def set_verbose_gf(verbose):
    """Set the verbosity level of getfem
    
    Args:
        verbose (int): the level of verbosity
    """

    gf.util_trace_level(verbose)
    gf.util_warning_level(verbose)

def set_verbose(verbose=1):
    """Set the verbosity level of scrimp (0: quiet, 1: info, 2: debug)
    
    In `quiet` mode, debug are saved in a log file.
    
    Args:
        verbose (int): the level of verbosity, defaults to 1
    """
    
    try:
        assert verbose in [0, 1, 2]
    except AssertionError as err:
        logging.error(
            f"Verbosity levels are 0, 1 or 2, unknown {verbose} level."
        )
        raise err
    
    # Remove all previously setted handlers associated with the root logger object.
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    if verbose==0:
        # Set the log file
        logging.basicConfig(
            filename=os.path.join(module_path, "log", "scrimp.log"),
            encoding="utf-8",
            level=logging.DEBUG,
            filemode="w",
            format="",
        )
        set_verbose_gf(0)
    elif verbose==1:
        logging.basicConfig(
            level=logging.INFO,
            format="",
        )
        set_verbose_gf(0)
    elif verbose==2:
        logging.basicConfig(
            level=logging.DEBUG,
        )
        set_verbose_gf(2)
