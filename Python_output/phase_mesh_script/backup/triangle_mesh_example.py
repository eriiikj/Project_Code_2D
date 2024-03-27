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

def get_line_conn(line_coord,line_ex,line_ey):
    
    # No. line segments
    nseg = np.shape(line_ex)[0]

    # Stack line_ex and line_ey into a single array for convenience
    lines_stacked = np.column_stack((line_ex.flatten(), line_ey.flatten()))
    
    # Find indices of line_segments in line_coord
    line_conn = np.where(np.all(line_coord == lines_stacked[:, None, :],axis=-1))[1]
    
    # Reshape to match line segments (nsegx2)
    line_conn = line_conn.reshape(nseg,2)
    
    return line_conn
        

if __name__ == "__main__":

    # --- Nodes ---
    P1 = [0, 0]
    P2 = [1, 0]
    P3 = [1, 1]
    P4 = [0, 1]
    corner_coord = np.vstack((P1, P2, P3, P4))
    ncorners     = np.shape(corner_coord)[0]


    # Many line segments
    ncoord          = 100
    xcoord          = np.linspace(0,1,ncoord)
    rand            = np.random.rand(ncoord)
    ylim1           = 0.51
    ylim2           = 0.5101
    ycoord          = rand*ylim1 + (1-rand)*ylim2
    line_coord      = np.zeros([ncoord,2])
    line_coord[:,0] = xcoord
    line_coord[:,1] = ycoord
    nseg            = ncoord - 1
    line_ex         = np.zeros([nseg,2])
    line_ex[:,0]    = xcoord[:-1]
    line_ex[:,1]    = xcoord[1:]
    line_ey         = np.zeros([nseg,2])
    line_ey[:,0]    = ycoord[:-1]
    line_ey[:,1]    = ycoord[1:]
    
    
    # A few line segments
    # line_ex    = np.array([[0,0.3],[1.0,0.7],[0.3,0.7]])
    # line_ey    = np.array([[0.45, 0.47], [0.47, 0.44], [0.47, 0.44]])
    # line_coord = lines_to_coord(line_ex, line_ey)
    ncoord     = np.shape(line_coord)[0]
    
    # Total coord (line_coord + domain corners)
    coord = np.vstack((line_coord, corner_coord))
    
    # Line segment connectivity of coords
    line_conn = get_line_conn(line_coord,line_ex,line_ey)
    
    # --- Constraint Delanuay triangularization using the triangle software ---
    tri_input     = dict(vertices=coord,segments=line_conn)
    triangulation = triangle.triangulate(tri_input,opts='pcq') #pca0.002YY 
    tris          = triangulation['triangles']
    tris          = tris.flatten()
    tris          = tris + 1
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
