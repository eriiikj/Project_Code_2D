#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 10:14:13 2024

@author: Erik Jacobsson
"""

import numpy as np
import gmsh
import sys
import triangle

# Initialize Gmsh
gmsh.initialize()

# Corner points
P1 = [0, 0, 0]
P2 = [1, 0, 0]
P3 = [1, 1, 0]
P4 = [0, 1, 0]

corner_coord = np.vstack((P1, P2, P3, P4))
ncorners     = np.shape(corner_coord)[0]

# Points along interpolated segment
xcoord = [0, 0.5, 1]
ycoord = [0.4, 0.5, 0.4]
zcoord = np.zeros_like(xcoord)
coord  = np.vstack((xcoord, ycoord, zcoord)).transpose()
nseg   = np.shape(coord)[0]
N      = nseg + ncorners

# xy points
xy = np.zeros(2 * N)
xy[:2 * ncorners:2] = corner_coord[:, 0]
xy[1:2 * ncorners:2] = corner_coord[:, 1]
xy[2 * ncorners::2] = coord[:, 0]
xy[2 * ncorners + 1::2] = coord[:, 1]

# # xyz points
# xyz = np.zeros(3 * N)
# xyz[:3 * ncorners:3]     = corner_coord[:, 0]
# xyz[1:3 * ncorners:3]    = corner_coord[:, 1]
# xyz[2:3 * ncorners:3]    = corner_coord[:, 2]
# xyz[3 * ncorners::3]     = coord[:, 0]
# xyz[3 * ncorners + 1::3] = coord[:, 1]
# xyz[3 * ncorners + 2::3] = coord[:, 2]

# Triangularize
tri_input     = dict(vertices=xy.reshape((N, 2)))
triangulation = triangle.triangulate(tri_input, 'a0.1')
tris1         = triangulation['triangles']
tris1         = tris1 + 1 # 1 based
tris1         = tris1.flatten()
newcoord      = triangulation['vertices']
N             = np.shape(newcoord)[0]

# xyz points
xyz = np.zeros(3 * N)
xyz[::3]  = newcoord[:, 0]
xyz[1::3] = newcoord[:, 1]
xyz[2::3] = np.zeros_like(newcoord[:, 0])


tris = gmsh.model.mesh.triangulate(xy)
# 
# v = [[0, 0], [0, 1], [1, 1], [1, 0]]
# tri_input2 = {'vertices': v}
# tri_input3 = dict(vertices=np.array(v).reshape((4, 2)))
# t = triangle.triangulate(tri_input3, 'a0.2')
# t['vertices'].tolist()
# [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0], [0.5, 0.5], [0.0, 0.5], [0.5, 0.0], [1.0, 0.5], [0.5, 1.0]]
# t['vertex_markers'].tolist()
# [[1], [1], [1], [1], [0], [1], [1], [1], [1]]
# t['triangles'].tolist()
# [[7, 2, 4], [5, 0, 4], [4, 8, 1], [4, 1, 5], [4, 0, 6], [6, 3, 4], [4, 3, 7], [4, 2, 8]]


# Add triangles to mesh
surf = gmsh.model.addDiscreteEntity(2)
gmsh.model.mesh.addNodes(2, surf, range(1, N + 1), xyz)
gmsh.model.mesh.addElementsByType(surf, 2, [], tris)


# Generate mesh
gmsh.model.mesh.generate(2)  # 2D mesh generation


# --- Graphical illustration of mesh ---
if 'close' not in sys.argv:
    gmsh.fltk.run()

entities = gmsh.model.getEntities()

# Finalize Gmsh
gmsh.finalize()
