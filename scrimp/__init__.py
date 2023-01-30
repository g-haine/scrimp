# # SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# #
# # Copyright (C) 2015-2022 Ghislain Haine
# #
# # See the LICENSE file in the root directory for license information.
# #
# # github: https://github.com/g-haine/scrimp
#
# """
# - file:             __init__.py
# - author:           Ghislain Haine
# - date:             22 nov. 2022
# - last modified:    13 dec. 2022
# - brief:            main module init file
#
# Usage for installation:
#
# - conda create env --file /path/to/scrimp/data/scrimp.yml
# - conda activate scrimp
# - conda develop /path/to/scrimp/
# """
#
# SCRIMP_VERSION = '0.5.0'
# SCRIMP_VERSION_MAJOR = 0
# SCRIMP_VERSION_MINOR = 5
# SCRIMP_VERSION_PATCH = 0
#
# import os
# import sys
# import importlib.util
# import getfem as gf
#
# def get_version():
#     """
#     Get the current version of SCRIMP
#     """
#
#     print('SCRIMP version', SCRIMP_VERSION)
#
# def check_default_path():
#     """
#     Check if the default path exists, and create it otherwise
#     """
#
#     path = set_default_path()
#     if not os.path.exists(path):
#         os.makedirs(path)
#
# def set_default_path():
#     """
#     Set the default path folder for outputs to the path of this file + outputs
#     """
#
#     return os.path.join(os.path.dirname(os.path.realpath(__file__)), 'outputs')
#
# def set_verbose_gf(verbose):
#     """
#     Set the verbosity level of getfem
#
#     :param verbose: the level
#     :type verbose: int
#     """
#
#     gf.util_trace_level(verbose)
#     gf.util_warning_level(verbose)
#
# def test_install():
#     """
#     This function will try to run all scripts in the `examples` folder
#
#     This allows to test the installation. Be aware that this takes time!
#     """
#
#     path_to_examples = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'examples')
#     filenames = next(os.walk(path_to_examples), (None, None, []))[2]
#
#     for file in filenames:
#         if not file=='__init__.py':
#             spec = importlib.util.spec_from_file_location(os.path.splitext(file)[0], os.path.join(path_to_examples, file))
#             mod = importlib.util.module_from_spec(spec)
#             sys.modules[os.path.splitext(file)[0]] = mod
#             spec.loader.exec_module(mod)
#             function = getattr(mod, os.path.splitext(file)[0])
#             print('Start: '+os.path.splitext(file)[0])
#             try:
#                 function()
#                 print('\033[92m ========== \033[0m')
#                 print('\033[92m Example', os.path.splitext(file)[0], 'PASS \033[0m')
#                 print('\033[92m ========== \033[0m')
#             except:
#                 print('\033[91m ========== \033[0m')
#                 print('\033[91m Example', os.path.splitext(file)[0], 'FAILS \033[0m')
#                 print('\033[91m ========== \033[0m')
#