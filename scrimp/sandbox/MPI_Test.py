#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 13:54:42 2023

@author: ghaine
"""

from IPython import get_ipython
ip = get_ipython()

# ip.run_cell("!python Allen-Cahn.py")

ip.run_cell("!mpiexec -n 2 python Allen-Cahn.py")