# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2025 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             setup.py
- authors:          Giuseppe Ferraro
- date:             31 may 2023
- brief:            setup file for install
"""

import setuptools

with open("scrimp/README.md", "r") as fh:

    long_description = fh.read()

setuptools.setup(

    name="scrimp", # Replace with your username

    version="1.1.1",

    author="Giuseppe Ferraro",

    author_email="giuseppe.ferraro@isae-supaero.fr",

    description="<Template Setup.py package>",

    long_description=long_description,

    long_description_content_type="text/markdown",

    url="https://github.com/g-haine/scrimp",

    packages=setuptools.find_packages(),

    classifiers=[

        "Programming Language :: Python :: 3",

        "License :: OSI Approved :: GNU General Public License version 3",

        "Operating System :: OS Independent",

    ],

    python_requires='>=3.6',

)
