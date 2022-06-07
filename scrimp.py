#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Elementary configuration script for SCRIMP.

Usage
------

- Activate the environment:

    conda activate scrimp

- Import the module by starting your script.py by:

    rootdir = "/path/to/your/scrimp/directory/"
    import sys
    sys.path.append(rootdir)
    import scrimp
    scrimp.setup_path(rootdir)
    
"""
import sys
import os

# Path management

def config_path_append(dir):
    """Append all subfolders of dir to sys.path."""
    n_dir = 0 # number of directory appended
    if (os.path.isdir(dir)):
        print("\t"+dir)
        sys.path.append(dir)
        n_dir += 1
    for root, dirs, files in os.walk(dir):
       for name in dirs:
          if (name != '__pycache__'):
              subdir = os.path.join(root, name)
              print("\t"+subdir)
              sys.path.append(subdir)
              n_dir += 1
    return n_dir

def setup_path(rootdir):
    """Setup the paths needed for SCRIMP."""
    print("Appending directories to sys.path using '"+rootdir+"' as root...")
    print("Data directories:")
    n_dir = config_path_append(os.path.join(rootdir,"data"))
    print(f"{n_dir} directories found.")    
    print("Source directories:")
    n_dir = config_path_append(os.path.join(rootdir,"src/modules"))
    print(f"{n_dir} directories found.")
   
    # If in rootdir, call setup_path
if os.path.isfile('scrimp.py'):

    print("--- Path configuration for SCRIMP")
    setup_path('.')
    print("--- Done.")
