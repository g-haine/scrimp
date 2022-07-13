#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions related to mesh handling using meshio.

@author: Florian Monteghetti
"""

import meshio
import json
import os
import numpy as np


def get_cells_type_str(mshfile):
    """ Get string that lists mesh cell types. """
    mesh = meshio.read(mshfile)
    return '|'.join(list(mesh.cells_dict.keys()))


def get_dimension(mshfile):
    """ Get geometrical dimension. """    
    mesh = meshio.read(mshfile)
    dim = 0
    for i in range(mesh.points.shape[1]):
        if (np.any(mesh.points[:,i])):
            dim += 1
    return dim

############## GMSH meshes

def create_submesh_from_gmsh(mesh, cell_type, outfile, prune_z=False):
    """
    Create a submesh from a gmsh mesh by extracting one cell type. High-order
    meshes are supported.

    Parameters
    ----------
    mesh : meshio.Mesh
        Mesh read from a gmsh file.
    
    cell_type : str
        Cell type to extract, e.g. 'line', 'line3', 'triangle', triangle6'.

    outfile : str
        Mesh filename.
    
    prune_z : boolean
        Remove z coordinates.

    Returns
    -------
    mesh : meshio.Mesh
        Mesh that describes cell_type. For convenience, the names of all
        physical entities tags are included (provided meshio implements
        export of field_data).
    
    """
        # Vertices of all cell_type cells
    cells = mesh.get_cells_type(cell_type)
    try:
        label = "gmsh:physical"
            # Physical entities tags of all cell_type cells
        cell_data = mesh.get_cell_data(label, cell_type)
    except:
        raise ValueError(f"Cannot read cell_data '{label}' for '{cell_type}'")
        # dictionnary describing all physical entities
    field_data= mesh.field_data
        # Create submesh
    out_mesh = meshio.Mesh(points=mesh.points, # node coordinates
                            cells={cell_type: cells}, 
                            cell_data={"physical_entities_tag":[cell_data]},
                                # ! Not implemented in meshio 4.3.11,
                                # which silently ignores field_data
                            field_data={"physical_entities_name":field_data})    
    if prune_z:
        # out_mesh.prune_z_0()
        out_mesh = prune_z_0(out_mesh)
        
        # Write to file
    print(f"Exporting nodes and '{cell_type}' to {outfile}...")
    meshio.write(outfile, out_mesh,file_format="xdmf")
    print("Done.")    


def prune_z_0(mesh, tol: float = 1.0e-13):
    """
    Remove third (z) component of points if it is 0 everywhere (up to a
    tolerance). This function is deprecated in meshio 7.0.0.
    This implementation is copied from a previous version of meshio.
    https://github.com/nschloe/meshio/pull/1228 (29/11/2021)
    
    Parameters
    ----------
    mesh : meshio.Mesh
    tol : float
    """
    if mesh.points.shape[0] == 0:
        return
    if mesh.points.shape[1] == 3 and np.all(np.abs(mesh.points[:, 2]) < tol):
        mesh.points = mesh.points[:, :2]
    return mesh

def export_phys_gmsh(mesh, outfile):
    """
    Export physical entities names, tags, and dimensions to a JSON
    file.

    Parameters
    ----------
    mesh : meshio.Mesh
        Mesh read from a gmsh file.
    
    outfile : str
        Ouput file name

    """
        # Test field_data
    if len(mesh.field_data)==0:
        raise ValueError("Mesh has no physical entities (field_data is empty)")
    
        # Build dictionnary
    field_data = {}
    for key in mesh.field_data:
            # convert numpy.ndarray to list
        field_data[key] = mesh.field_data[key].tolist()
        # Export to json
    print(f"Exporting physical entities name/tag/dim to {outfile}...")
    with open(outfile,'w') as outfiled:
        json.dump(field_data, outfiled)
    print("Done.")
    return outfile
    
def read_phys_gmsh(jsonfile):
    """ Read gmsh physical entities from JSON file. Output format:
        gmsh_phys[dim-1]={name:tag} """
    d = dict()
    with open(jsonfile) as infile:
        d = json.load(infile)
        
        # Convert data[name] = [tag,dim] (dict(list))
        # to      data[dim-1] = {name:tag} (list(dict))
        
        # Find max. dimension
    dim_max = 0
    for name in d:
        if (len(d[name])>2):
            raise ValueError(f"JSON file contains more than two integers for '{name}'")
        dim = d[name][1]
        dim_max = max(dim_max,dim)
        # Build list of dict
    tags = [dict() for x in range(dim_max)]
    for name in d:
        tag= d[name][0]
        dim = d[name][1]
        tags[dim-1][name] = list([tag]) # important, must be a list

    return tags

def print_diagnostic_gmsh(meshfile):
    """
    Print information about a gmsh mesh file.
    
    Parameters
    ----------
    
    meshfile : str
        Path to gmsh mesh file (all formats, ASCII or binary).

    """
        # read gmsh file (format 2.2)
    msh = meshio.read(meshfile)
    print(f"---- Mesh {meshfile}")
    print(f"{msh.points.shape[0]} nodes")
    print(f"{len(msh.cells)} cell types:")
    for i in range(len(msh.cells)):
        print(f"\t#{i} type: '{msh.cells[i].type}' total: {msh.cells[i].data.shape[0]}")
    print(f"{len(msh.field_data)} physical entities:")
    for physName in msh.field_data:
        tag = msh.field_data[physName][0]
        print(f"\t'{physName}' tag: {tag} dim: {msh.field_data[physName][1]}")
        for i in range(len(msh.cells)):
            total = np.count_nonzero(msh.cell_data['gmsh:physical'][i]==tag)
            if (total>0):
                print(f"\t\t {total} '{msh.cells[i].type}'")
    print("----")

def convert_gmsh_to_XDMF(meshfile,prune_z=False):
    """
    Comvert gmsh mesh (all formats and orders) to XDMF using meshio.
    
    The following files are written: 
        meshfile_celltype.xdmf:
                nodes coordinates, celltype vertices and physical tags
        meshfile_phys.json:
                Names/tags/dimensions of gmsh physical entities.
    Example: 
        (1) For a 2D first-order mesh 'mesh.gmsh':
            mesh_triangle.xdmf, mesh_line.xdmf, mesh_phys.json
        (2) For a 3D second-order mesh 'mesh.gmsh':
            mesh_tetra10.xdmf, mesh_triangle6, mesh_phys.json
            
    Parameters
    ----------
    
    meshfile : str
        Path to gmsh mesh file (all formats, ASCII or binary).
        
    prune_z : boolean
        Remove z coordinates from mesh (only if all z coordinates are null).
        
    Returns
    -------
    xdmfiles : dict of str
        Dictionary of output filenames.

    """
    
        # Read the gmsh file using meshio
    mesh = meshio.read(meshfile)
        # Do not remove non-null z coordinates
    if (mesh.points.shape[1]==3) and (np.any(mesh.points[:,2])):
        prune_z = False
        # Export vertices of each cell type
    cell_types = mesh.cells_dict # dictionary of cell types
    print(f"Mesh '{meshfile}' contains cell types:\n\t{list(cell_types.keys())}")
    xdmfiles = dict()
    for cell_type in cell_types: # Loop over each cell type
        outfile = os.path.splitext(meshfile)[0]+"_"+cell_type+".xdmf"
        create_submesh_from_gmsh(mesh, cell_type, outfile, prune_z=prune_z)
        xdmfiles[cell_type]=outfile 
        # Export gmsh physical entities
    outfile = os.path.splitext(meshfile)[0]+"_phys.json"
    export_phys_gmsh(mesh, outfile)
    xdmfiles['gmsh_physical_entities'] = outfile
    return xdmfiles
