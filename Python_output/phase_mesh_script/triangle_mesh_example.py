#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 15:18:18 2024

@author: Erik Jacobsson
"""

import numpy as np
import gmsh
import sys
import triangle

def lines_to_coord(line_ex, line_ey):
    N             = np.shape(line_ex)[0]
    coord         = np.zeros((2 * N, 2))
    coord[:N, 0]  = line_ex[:, 0]
    coord[N:, 0]  = line_ex[:, 1]
    coord[:N, 1]  = line_ey[:, 0]
    coord[N:, 1]  = line_ey[:, 1]
    unique_coords = np.unique(coord, axis=0)
    return unique_coords

def xyz_from_coord(coord):
    N         = np.shape(coord)[0]
    xyz       = np.zeros(3 * N)
    xyz[::3]  = coord[:, 0]
    xyz[1::3] = coord[:, 1]
    xyz[2::3] = np.zeros_like(coord[:, 0])
    return xyz

def get_segments(coord,line_ex,line_ey):
    
    # No. line segments
    nseg = np.shape(line_ex)[0]
    
    # Define segement idex
    segment_idx = np.zeros([nseg,2],dtype=int)
    
    for iseg in range(nseg):
        P1 = np.array([line_ex[iseg,0],line_ey[iseg,0]])
        P2 = np.array([line_ex[iseg,1],line_ey[iseg,1]])
        P1idx = np.where((coord == P1).all(axis=1))[0][0]
        P2idx = np.where((coord == P2).all(axis=1))[0][0]
        segment_idx[iseg,:] = [P1idx,P2idx]
        
        
    s = 9
    
    return segment_idx
        

if __name__ == "__main__":

    # --- Nodes ---
    P1 = [0, 0]
    P2 = [1, 0]
    P3 = [1, 1]
    P4 = [0, 1]
    corner_coord = np.vstack((P1, P2, P3, P4))
    ncorners     = np.shape(corner_coord)[0]
    
    line_ex    = np.array([[0,0.3],[1.0,0.7],[0.3,0.7]])
    line_ey    = np.array([[0.45, 0.47], [0.47, 0.44], [0.47, 0.44]])
    line_coord = lines_to_coord(line_ex, line_ey)
    ncoord     = np.shape(line_coord)[0]
    
    # Total coord
    coord = np.vstack((line_coord, corner_coord))
    
    # Line segment connectivity of coords
    segment_idx = get_segments(coord,line_ex,line_ey)
    
    # --- Triangularize using the triangle software ---
    tri_input     = dict(vertices=coord, segment=segment_idx) # segment=segment_idx
    triangulation = triangle.triangulate(tri_input,opts='pca0.002YY') #pca0.02Y 
    tris          = triangulation['triangles']
    tris          = tris.flatten()
    tris          = tris + 1  # 1 based
    newcoord      = triangulation['vertices']
    N             = np.shape(newcoord)[0]

    # Extract xyz data
    xyz = xyz_from_coord(newcoord)
    
    # --- Define mesh ---
    
    # Initialize Gmsh
    gmsh.initialize()
    
    # Add triangles to mesh
    surf = gmsh.model.addDiscreteEntity(2)
    gmsh.model.mesh.addNodes(2, surf, range(1, N + 1), xyz)
    gmsh.model.mesh.addElementsByType(surf, 2, [], tris)
    
    # Generate mesh
    gmsh.model.mesh.generate(2)  # 2D mesh generation
    
    # --- Graphical illustration of mesh ---
    if 'close' not in sys.argv:
        gmsh.fltk.run()
    
    # Gmsh entities
    entities = gmsh.model.getEntities()
    
    # Finalize Gmsh
    gmsh.finalize()
