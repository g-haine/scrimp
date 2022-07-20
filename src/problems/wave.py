# -*- coding: utf-8 -*-
"""
Simplified methods for specific examples

@author: Ghislain Haine
"""

def set_forms(system,causality,parameters,geofile,h,refinement):
    """
    
    """

    forms = dict()
    info = ''
        
    mesh = mesh_from_geo(geofile,h,refinement)
    
    match [system,causality]:
        case ['wave','Neumann']:
            
        
        case ['wave','Dirichlet']:
            
        
        case ['heat','Neumann']:
            
        
        case ['wave','Dirichlet']:
            
        
        case _:
            
    
    return forms

def set_spaces():

    return W

def mesh_from_geo(geofile,h,refinement):
    
    assert geofile in ...
    
    
    return dmesh
