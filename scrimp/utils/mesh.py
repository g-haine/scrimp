# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2022 Ghislain Haine
#
# See the LICENSE file in the root directory for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             utils/mesh.py
- author:           Ghislain Haine
- date:             22 nov. 2022
- last modified:    13 dec. 2022
- brief:            built-in geometries for direct use in SCRIMP
"""

import os
import math
import gmsh
import numpy as np
import getfem as gf
from scrimp import set_default_path, check_default_path

check_default_path()

def built_in_geometries():
    """
    A function to get all the infos about available built_in geometries
    
    :param: None
    :return: Print informations about all available built_in geometries
    """
    
    print('Available geometries with parameters:')
    print('=====================================')
    print('')
    print('* `Interval`, {"L": 1, "N": 20}, an segment of size L with N linearly spaced points')
    print('---> Domain `Omega`: 1, Left boundary `Gamma_Left`: 10, Right boundary `Gamma_Right`: 11')
    print('')
    print('* `Disk`, {"R": 1, "h": 0.1}, a disk of radius R with mesh size h')
    print('---> Domain `Omega`: 1, Boundary `Gamma`: 10')
    print('')
    print('* `Rectangle`, {"L": 2, "l": 1, "h": 0.1}, a rectangle of size L x l with mesh size h')
    print('---> Domain `Omega`: 1, Bottom boundary `Gamma_Bottom`: 10, Right boundary `Gamma_Right`: 11, Top boundary `Gamma_Top`: 12, Left boundary `Gamma_Left`: 13')
    
def Interval(parameters={'L': 1., 'h': 0.05}, terminal=1):
    """
    The geometry of a segment (0,L) with mesh size h
    
    :param parameters: The dictionary of parameters for the geometry
    :type parameters: dict
    :param terminal: An option to print meshing infos in the prompt, value `0` (quiet) or `1` (verbose, default)
    :type terminal: int
    
    :return: The mesh to use with getfem, the dimension, a dict of regions with getfem indices for dim n and a dict of regions with getfem indices for dim n-1
    :rtype: list[gf.Mesh, int, dict, dict]
    """
    
    L = parameters['L']
    h = parameters['h']
    
    mesh = gf.Mesh('regular simplices', np.arange(0,L+h,h,))
    cvids = mesh.cvid()
    all_faces = -1*np.ones(shape=cvids.size)
    mesh.set_region(1,[cvids,all_faces]) # Do we need to remove boundaries?
    mesh.set_region(10, mesh.faces_from_pid(mesh.pid_from_coords(0.)))
    mesh.set_region(11, mesh.faces_from_pid(mesh.pid_from_coords(L)))
    
    print('Interval (0,'+str(L)+') has been meshed')
    
    return [mesh, 1, {'Omega': 1}, {'Gamma_Left': 10, 'Gamma_Right': 11}]
    

def Disk(parameters={'R': 1., 'h': 0.1}, terminal=1):
    """
    The geometry of a Disk center in (0,0) with radius R and mesh size h
    
    :param parameters: The dictionary of parameters for the geometry
    :type parameters: dict
    :param terminal: An option to print meshing infos in the prompt, value `0` (quiet) or `1` (verbose, default)
    :type terminal: int
    
    :return: The mesh to use with getfem, the dimension, a dict of regions with getfem indices for dim n and a dict of regions with getfem indices for dim n-1
    :rtype: list[gf.Mesh, int, dict, dict]
    """
    
    h = parameters['h']
    R = parameters['R']
    
    # Init GMSH
    gmsh.initialize()
    # Ask GMSH to display information in the terminal
    gmsh.option.setNumber('General.Terminal', terminal)

    # Create a model and name it "MyCircle"
    model = gmsh.model
    model.add('Disk')

    # Create Point for the center of the circle
    center = model.geo.addPoint(0,0,0, h, 100)
    # Create 3 Points on the circle
    points = []
    for j in range(3):
      points.append(model.geo.addPoint(R*math.cos(2*math.pi*j/3), R*math.sin(2*math.pi*j/3), 0, h, j+101))
    # Create 3 circle arc
    lines = []
    for j in range(3):
      lines.append(model.geo.addCircleArc(points[j],center,points[(j+1)%3], j+10))
    # Curveloop and Surface
    curveloop = model.geo.addCurveLoop(lines, 13)
    model.geo.addPlaneSurface([curveloop], 1)

    # This command is mandatory and synchronize CAD with GMSH Model. The less you launch it, the better it is for performance purpose
    gmsh.model.geo.synchronize()
    
    # Mesh (2D)
    dim = 2
    model.mesh.generate(dim)
    # Write on disk
    path = set_default_path()
    gmsh.write(os.path.join(path, 'Disk.msh'))
    # Finalize GMSH
    gmsh.finalize()
    
    print('Disk with radius', R, 'has been meshed')
    
    # GetFEM mesh
    mesh = gf.Mesh('import', 'gmsh_with_lower_dim_elt', os.path.join(path, 'Disk.msh'))
    # Glue boundaries CircleArc
    mesh.region_merge(10,11)
    mesh.region_merge(10,12)
    # Delete those which have been glued
    mesh.delete_region([11,12])
    
    return [mesh, dim, {'Omega': 1}, {'Gamma': 10}]

def Rectangle(parameters={'L': 2., 'l': 1, 'h': 0.1}, terminal=1):
    """
    The geometry of a Rectangle (0,L)x(0,l) with mesh size h
    
    :param parameters: The dictionary of parameters for the geometry
    :type parameters: dict
    :param terminal: An option to print meshing infos in the prompt, value `0` (quiet) or `1` (verbose, default)
    :type terminal: int
    
    :return: The mesh to use with getfem, the dimension, a dict of regions with getfem indices for dim n and a dict of regions with getfem indices for dim n-1
    :rtype: list[gf.Mesh, int, dict, dict]
    """
    
    h = parameters['h']
    L = parameters['L']
    l = parameters['l']
    
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', terminal)

    model = gmsh.model
    model.add('Rectangle')
    
    A = model.geo.addPoint(0, 0, 0, h, 100)
    B = model.geo.addPoint(L, 0, 0, h, 101)
    C = model.geo.addPoint(L, l, 0, h, 102)
    D = model.geo.addPoint(0, l, 0, h, 103)
    
    Bottom = model.geo.addLine(A, B, 10)
    Right = model.geo.addLine(B, C, 11)
    Top = model.geo.addLine(C, D, 12)
    Left = model.geo.addLine(D, A, 13)
    
    curveloop = model.geo.addCurveLoop([Bottom, Right, Top, Left], 14)
    model.geo.addPlaneSurface([curveloop], 1)
    
    gmsh.model.geo.synchronize()
    
    dim = 2
    model.mesh.generate(dim)
    path = set_default_path()
    gmsh.write(os.path.join(path, 'Disk.msh'))
    gmsh.finalize()
    
    print('Rectangle (0,'+str(L)+')x(0,'+str(l)+') has been meshed')
    
    mesh = gf.Mesh('import', 'gmsh_with_lower_dim_elt', os.path.join(path, 'Disk.msh'))
    
    return [mesh, dim, {'Omega': 1}, {'Gamma_Bottom': 10, 'Gamma_Right': 11, 'Gamma_Top': 12, 'Gamma_Left': 13}]
    