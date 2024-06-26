# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2024 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             utils/mesh.py
- authors:          Ghislain Haine
- date:             22 nov. 2022
- brief:            built-in geometries for direct use in SCRIMP
"""

import petsc4py
import sys

petsc4py.init(sys.argv)
from petsc4py import PETSc

comm = PETSc.COMM_WORLD
rank = comm.getRank()

import os
import math
import gmsh
import numpy as np
import getfem as gf
import logging

import scrimp.utils.config
outputs_path = scrimp.utils.config.outputs_path

def built_in_geometries():
    """A function to get all the infos about available built_in geometries"""

    if rank == 0:
        logging.info("Available geometries with parameters:")
        logging.info("=====================================")
        logging.info("")
        logging.info(
            '* `Interval`, {"L": 1, "N": 20}, an segment of size L with N linearly spaced points'
        )
        logging.info(
            "---> Domain `Omega`: 1, Left boundary `Gamma_Left`: 10, Right boundary `Gamma_Right`: 11"
        )
        logging.info("")
        logging.info(
            '* `Disk`, {"R": 1, "h": 0.1}, a disk of radius R with mesh size h'
        )
        logging.info("---> Domain `Omega`: 1, Boundary `Gamma`: 10")
        logging.info("")
        logging.info(
            '* `Rectangle`, {"L": 2, "l": 1, "h": 0.1}, a rectangle of size L x l with mesh size h'
        )
        logging.info(
            "---> Domain `Omega`: 1, Bottom boundary `Gamma_Bottom`: 10, Right boundary `Gamma_Right`: 11, Top boundary `Gamma_Top`: 12, Left boundary `Gamma_Left`: 13"
        )
        logging.info("")
        logging.info(
            '* `Concentric`, {"R": 1, "r": 0.6, "h": 0.1}, a disk of radius r surrounded by an annulus of radii r and R with mesh size h'
        )
        logging.info(
            "---> Domain `Omega_Disk`: 1, `Omega_Annulus`: 2, Interface `Interface`: 10, Boundary `Gamma`: 20"
        )
        logging.info("")
        logging.info(
            '* `Ball`, {"R": 1, "h": 0.1}, a ball of radius R with mesh size h'
        )
        logging.info("---> Domain `Omega`: 1, Boundary `Gamma`: 10")
        logging.info("")


def Interval(parameters={"L": 1.0, "h": 0.05}, refine=0, terminal=1):
    """The geometry of a segment (0,L) with mesh size h

    - Domain `Omega`: 1,
    - Left boundary `Gamma_Left`: 10,
    - Right boundary `Gamma_Right`: 11

    Args:
        - parameters (dict): The dictionary of parameters for the geometry
        - refine (int): Ask for iterative refinements by splitting elements
        - terminal (int): An option to print meshing infos in the prompt, value `0` (quiet) or `1` (verbose, default)

    Returns:
        list[gf.Mesh, int, dict, dict]: The mesh to use with getfem, the dimension, a dict of regions with getfem indices for dim n and a dict of regions with getfem indices for dim n-1
    """

    L = parameters["L"]
    h = parameters["h"] / (refine + 1)

    mesh = gf.Mesh(
        "regular simplices",
        np.arange(
            0,
            L + h,
            h,
        ),
    )
    cvids = mesh.cvid()
    all_faces = -1 * np.ones(shape=cvids.size)
    mesh.set_region(1, [cvids, all_faces])  # Do we need to remove boundaries?
    mesh.set_region(10, mesh.faces_from_pid(mesh.pid_from_coords(0.0)))
    mesh.set_region(11, mesh.faces_from_pid(mesh.pid_from_coords(L)))

    if rank == 0:
        logging.info(f"Interval (0, {L}) has been meshed")

    return [
        mesh,
        1,
        {"Omega": 1},
        {"Gamma_Left": 10, "Gamma_Right": 11},
    ]


def Disk(parameters={"R": 1.0, "h": 0.1}, refine=0, terminal=1):
    """The geometry of a Disk center in (0,0) with radius R and mesh size h

    - Domain `Omega`: 1,
    - Boundary `Gamma`: 10

    Args:
        - parameters (dict): The dictionary of parameters for the geometry
        - refine (int): Ask for iterative refinements by splitting elements
        - terminal (int): An option to print meshing infos in the prompt, value `0` (quiet) or `1` (verbose, default)

    Returns:
        list[gf.Mesh, int, dict, dict]: The mesh to use with getfem, the dimension, a dict of regions with getfem indices for dim n and a dict of regions with getfem indices for dim n-1
    """

    h = parameters["h"]
    R = parameters["R"]

    # Init GMSH
    gmsh.initialize()
    # Ask GMSH to display information in the terminal
    gmsh.option.setNumber("General.Terminal", terminal)
    gmsh.option.setNumber("General.NumThreads", 0)  # Use system default

    # Create a model and name it "MyCircle"
    model = gmsh.model
    model.add("Disk")

    # Create Point for the center of the circle
    center = model.geo.addPoint(0, 0, 0, h, 100)
    # Create 3 Points on the circle
    points = []
    for j in range(3):
        points.append(
            model.geo.addPoint(
                R * math.cos(2 * math.pi * j / 3),
                R * math.sin(2 * math.pi * j / 3),
                0,
                h,
                j + 101,
            )
        )
    # Create 3 circle arc
    lines = []
    for j in range(3):
        lines.append(
            model.geo.addCircleArc(points[j], center, points[(j + 1) % 3], j + 10)
        )
    # Curveloop and Surface
    curveloop = model.geo.addCurveLoop(lines, 13)
    model.geo.addPlaneSurface([curveloop], 1)

    # This command is mandatory and synchronize CAD with GMSH Model. The less you launch it, the better it is for performance purpose
    gmsh.model.geo.synchronize()

    # Mesh (2D)
    dim = 2
    model.mesh.generate(dim)
    # Refine
    for i in range(refine):
        model.mesh.refine()
    # Write on disk
    if rank==0:
        gmsh.write(os.path.join(outputs_path, "mesh", "Disk.msh"))
    # Finalize GMSH
    comm.barrier()
    gmsh.finalize()

    if rank == 0:
        logging.info(f"Disk with radius {R} has been meshed")

    # GetFEM mesh
    mesh = gf.Mesh("import", "gmsh_with_lower_dim_elt", os.path.join(outputs_path, "mesh", "Disk.msh"))
    # Glue boundaries CircleArc
    mesh.region_merge(10, 11)
    mesh.region_merge(10, 12)
    # Delete those which have been glued
    mesh.delete_region([11, 12])

    return [
        mesh,
        dim,
        {"Omega": 1},
        {"Gamma": 10},
    ]


def Rectangle(parameters={"L": 2.0, "l": 1, "h": 0.1}, refine=0, terminal=1):
    """The geometry of a Rectangle (0,L)x(0,l) with mesh size h

    - Domain `Omega`: 1,
    - Bottom boundary `Gamma_Bottom`: 10,
    - Right boundary `Gamma_Right`: 11,
    - Top boundary `Gamma_Top`: 12,
    - Left boundary `Gamma_Left`: 13

    Args:
        - parameters (dict): The dictionary of parameters for the geometry
        - refine (int): Ask for iterative refinements by splitting elements
        - terminal (int): An option to print meshing infos in the prompt, value `0` (quiet) or `1` (verbose, default)

    Returns:
        list[gf.Mesh, int, dict, dict]: The mesh to use with getfem, the dimension, a dict of regions with getfem indices for dim n and a dict of regions with getfem indices for dim n-1
    """

    h = parameters["h"]
    L = parameters["L"]
    l = parameters["l"]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", terminal)
    gmsh.option.setNumber("General.NumThreads", 0)  # Use system default

    model = gmsh.model
    model.add("Rectangle")

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

    for i in range(refine):
        gmsh.model.mesh.refine()
    if rank==0:
        gmsh.write(os.path.join(outputs_path, "mesh", "Rectangle.msh"))
    # Finalize GMSH
    comm.barrier()
    gmsh.finalize()

    if rank == 0:
        logging.info("Rectangle (0," + str(L) + ")x(0," + str(l) + ") has been meshed")

    mesh = gf.Mesh("import", "gmsh_with_lower_dim_elt", os.path.join(outputs_path, "mesh", "Rectangle.msh"))

    return [
        mesh,
        dim,
        {"Omega": 1},
        {"Gamma_Bottom": 10, "Gamma_Right": 11, "Gamma_Top": 12, "Gamma_Left": 13},
    ]


def Concentric(parameters={"R": 1.0, "r": 0.6, "h": 0.1}, refine=0, terminal=1):
    """The geometry of a Disk of radius r surrounded by an annulus of radii r and R with mesh size h

    - Domain `Omega_Disk`: 1,
    - Domain `Omega_Annulus`: 2,
    - Interface `Interface`: 10,
    - Boundary `Gamma`: 20

    Args:
        - parameters (dict): The dictionary of parameters for the geometry
        - refine (int): Ask for iterative refinements by splitting elements
        - terminal (int): An option to print meshing infos in the prompt, value `0` (quiet) or `1` (verbose, default)

    Returns:
        list[gf.Mesh, int, dict, dict]: The mesh to use with getfem, the dimension, a dict of regions with getfem indices for dim n and a dict of regions with getfem indices for dim n-1
    """

    h = parameters["h"]
    R = parameters["R"]
    r = parameters["r"]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", terminal)
    gmsh.option.setNumber("General.NumThreads", 0)  # Use system default

    model = gmsh.model
    model.add("Concentric")

    center = model.geo.addPoint(0, 0, 0, h, 100)
    points_int = []
    points_ext = []
    for j in range(3):
        points_int.append(
            model.geo.addPoint(
                r * math.cos(2 * math.pi * j / 3),
                r * math.sin(2 * math.pi * j / 3),
                0,
                h,
                j + 101,
            )
        )
        points_ext.append(
            model.geo.addPoint(
                R * math.cos(2 * math.pi * j / 3),
                R * math.sin(2 * math.pi * j / 3),
                0,
                h,
                j + 201,
            )
        )

    lines_int = []
    lines_ext = []
    for j in range(3):
        lines_int.append(
            model.geo.addCircleArc(
                points_int[j], center, points_int[(j + 1) % 3], j + 10
            )
        )
        lines_ext.append(
            model.geo.addCircleArc(
                points_ext[j], center, points_ext[(j + 1) % 3], j + 20
            )
        )
    # Curveloop and Surface
    curveloop_int = model.geo.addCurveLoop(lines_int, 13)
    curveloop_ext = model.geo.addCurveLoop(lines_ext, 23)
    model.geo.addPlaneSurface([curveloop_int], 1)
    model.geo.addPlaneSurface([curveloop_ext, -curveloop_int], 2)

    gmsh.model.geo.synchronize()

    dim = 2
    model.mesh.generate(dim)
    for i in range(refine):
        gmsh.model.mesh.refine()
    if rank==0:
        gmsh.write(os.path.join(outputs_path, "mesh", "Concentric.msh"))
    # Finalize GMSH
    comm.barrier()
    gmsh.finalize()

    if rank == 0:
        logging.info(
            f"A disk of radius {r} surrounded by an annulus of radii {r} and {R} has been meshed"
        )

    mesh = gf.Mesh("import", "gmsh_with_lower_dim_elt", os.path.join(outputs_path, "mesh", "Concentric.msh"))
    mesh.region_merge(10, 11)
    mesh.region_merge(10, 12)
    mesh.delete_region([11, 12])

    mesh.region_merge(20, 21)
    mesh.region_merge(20, 22)
    mesh.delete_region([21, 22])
    
    # Ensure interface belongs to both domain 
    mesh.region_merge(1, 10)
    mesh.region_merge(2, 10)

    return [
        mesh,
        dim,
        {"Omega_Disk": 1, "Omega_Annulus": 2},
        {"Interface": 10, "Gamma": 20},
    ]


def Ball(parameters={"R": 1.0, "h": 0.1}, refine=0, terminal=1):
    """The geometry of a Ball of radius R centered in (0,0,0) with mesh size h

    - Domain `Omega`: 1,
    - Boundary `Gamma`: 10

    Args:
        - parameters (dict): The dictionary of parameters for the geometry
        - refine (int): Ask for iterative refinements by splitting elements
        - terminal (int): An option to print meshing infos in the prompt, value `0` (quiet) or `1` (verbose, default)

    Returns:
        list[gf.Mesh, int, dict, dict]: The mesh to use with getfem, the dimension, a dict of regions with getfem indices for dim n and a dict of regions with getfem indices for dim n-1
    """

    h = parameters["h"]
    R = parameters["R"]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", terminal)
    gmsh.option.setNumber("General.NumThreads", 0)  # Use system default

    model = gmsh.model
    model.add("Ball")

    model.geo.addPoint(0, 0, 0, h, 1001)

    model.geo.addPoint(1, 0, 0, h, 1002)
    model.geo.addPoint(0, 1, 0, h, 1003)
    model.geo.addPoint(0, 0, 1, h, 1004)
    model.geo.addPoint(-1, 0, 0, h, 1005)
    model.geo.addPoint(0, -1, 0, h, 1006)
    model.geo.addPoint(0, 0, -1, h, 1007)

    model.geo.addCircleArc(1002, 1001, 1003, 101)
    model.geo.addCircleArc(1003, 1001, 1005, 102)
    model.geo.addCircleArc(1005, 1001, 1006, 103)
    model.geo.addCircleArc(1006, 1001, 1002, 104)
    model.geo.addCircleArc(1002, 1001, 1007, 105)
    model.geo.addCircleArc(1007, 1001, 1005, 106)
    model.geo.addCircleArc(1005, 1001, 1004, 107)
    model.geo.addCircleArc(1004, 1001, 1002, 108)
    model.geo.addCircleArc(1006, 1001, 1007, 109)
    model.geo.addCircleArc(1007, 1001, 1003, 110)
    model.geo.addCircleArc(1003, 1001, 1004, 111)
    model.geo.addCircleArc(1004, 1001, 1006, 112)

    model.geo.addCurveLoop([101, 111, 108], 121)
    model.geo.addCurveLoop([102, 107, -111], 122)
    model.geo.addCurveLoop([103, -112, -107], 123)
    model.geo.addCurveLoop([104, -108, 112], 124)
    model.geo.addCurveLoop([105, 110, -101], 125)
    model.geo.addCurveLoop([-102, -110, 106], 126)
    model.geo.addCurveLoop([-103, -106, -109], 127)
    model.geo.addCurveLoop([-104, 109, -105], 128)

    model.geo.addSurfaceFilling([121], 10, 1001)
    model.geo.addSurfaceFilling([122], 11, 1001)
    model.geo.addSurfaceFilling([123], 12, 1001)
    model.geo.addSurfaceFilling([124], 13, 1001)
    model.geo.addSurfaceFilling([125], 14, 1001)
    model.geo.addSurfaceFilling([126], 15, 1001)
    model.geo.addSurfaceFilling([127], 16, 1001)
    model.geo.addSurfaceFilling([128], 17, 1001)

    model.geo.addSurfaceLoop([10, 11, 12, 13, 14, 15, 16, 17], 20)

    model.geo.addVolume([20], 1)

    gmsh.model.geo.synchronize()

    dim = 3
    model.mesh.generate(dim)
    for i in range(refine):
        gmsh.model.mesh.refine()
    if rank==0:
        gmsh.write(os.path.join(outputs_path, "mesh", "Ball.msh"))
    # Finalize GMSH
    comm.barrier()
    gmsh.finalize()

    if rank == 0:
        logging.info("A ball of radius " + str(R) + " has been meshed")

    mesh = gf.Mesh("import", "gmsh", os.path.join(outputs_path, "mesh", "Ball.msh"))
    mesh.region_merge(10, 11)
    mesh.region_merge(10, 12)
    mesh.region_merge(10, 13)
    mesh.region_merge(10, 14)
    mesh.region_merge(10, 15)
    mesh.region_merge(10, 16)
    mesh.region_merge(10, 17)
    mesh.delete_region([11, 12, 13, 14, 15, 16, 17])

    return [
        mesh,
        dim,
        {"Omega": 1},
        {"Gamma": 10},
    ]
