# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 09:44:29 2020

@author: Erik Jacobsson
"""

import numpy as np
import json
import visvis as vv
import sys
import time
import os
import shutil
import scipy.io


import matplotlib.pyplot as plt
import itertools
import scipy.io

import pyvtk as vtk

import subprocess

import gmsh
import triangle

# Calfem
import calfem.core as cfc


class InputData(object):
    """Class for defining input data to our model."""
    def __init__(self):
        
        # --- Define private variabels, standrad values

        # Version
        self.version            = 1
        
        # Geometry [microns]
        self.grain              = 0
        self.a                  = 0
        self.ex                 = 0
        self.ey                 = 0 
        self.coord              = 0
        self.enod               = 0
        self.axisbc             = 0
        self.line_ex            = 0
        self.line_ey            = 0
        self.line_ex_all        = 0
        self.line_ey_all        = 0
        self.nseg               = 0
        self.line_seg_all       = 0
        self.sep_lines          = 0
        self.nsep_lines         = 0
        self.ngrains            = 0
        self.P1                 = 0
        self.P2                 = 0
        self.P3                 = 0
        self.P4                 = 0
        self.geoms              = 0
        self.interface_marker   = 0
        self.eq_comp            = 0
        self.material           = 0
        self.tot_imc_area_init  = 0
        self.tot_imc_area       = 0
        self.node_tags          = 0
        self.spline_tags        = 0
        self.line_coord         = 0
        self.line_conn          = 0
        self.all_line_conn      = 0
        self.boundary_coord     = 0
        self.boundary_conn      = 0
        self.mesh_points        = 0
        self.indNodBd           = 0
        self.indElemBd          = 0
        self.indLocalEdgBd      = 0

        # Location
        self.location           = None
        self.python_location    = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/Python_output'
        self.fortran_location   = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/build/src/main'
        
  
    
     # --- Define the geometry ---
    
    def geometry(self):
        """Define geometry in gmsh"""
        
        
        # --- Constraint Delanuay triangularization using triangle ---
        # tri_input     = dict(vertices=self.boundary_coord,segments=self.boundary_conn)
        # tri_input     = dict(vertices=self.line_coord,segments=self.line_conn) 
        tri_input     = dict(vertices=self.mesh_points,segments=self.all_line_conn) 
        # tri_input     = dict(vertices=self.boundary_coord,segments=self.boundary_conn) 
        print('Creating a constraint delanuay triangularization...')
        triangulation = triangle.triangulate(tri_input,opts='pq') #pqYY 
        print('Finished delanuay triangularization')
        tris          = triangulation['triangles']
        tris          = tris.flatten()
        tris          = tris + 1
        vertices      = triangulation['vertices']
        N             = np.shape(vertices)[0]
        
        
        
        # Extract xyz data
        xyz = self.xyz_from_coord(vertices)
            
        # Add triangles to gmsh
        surf = gmsh.model.addDiscreteEntity(2)
        gmsh.model.mesh.addNodes(2, surf, range(1, N + 1), xyz)
        gmsh.model.mesh.addElementsByType(surf, 2, [], tris)
        
        # Create the relevant Gmsh data structures from gmsh model
        gmsh.model.geo.synchronize()
        
        return triangulation
    
    def xyz_from_coord(self, coord):
        N         = np.shape(coord)[0]
        xyz       = np.zeros(3 * N)
        xyz[::3]  = coord[:, 0]
        xyz[1::3] = coord[:, 1]
        xyz[2::3] = np.zeros_like(coord[:, 0])
        return xyz
    
    
    
    def save(self, filename):
        """Save input data to dictionary."""
        input_data_file = {}
        
        # Version
        input_data_file["version"]                = self.version
        
        # Location
        input_data_file["location"]               = self.location  
        
        with open(filename, "w") as ofile:
            json.dump(input_data_file, ofile, sort_keys = True, indent = 4)
            
    

    def load(self, filename):
        """Read indata from file."""

        with open(filename, "r") as ifile:
            input_data_file = json.load(ifile)

        # Version
        self.version          = input_data_file["version"]
        
        # Location
        self.location         = input_data_file["location"]
        
        
        

        with open(filename, "r") as ifile:
            
            # Load the .mat file
            data = scipy.io.loadmat(filename)

            # Access the variable by its name
            self.ex     = np.transpose(data['ex'])*1e3
            self.ey     = np.transpose(data['ey'])*1e3
            self.coord  = np.transpose(data['coord'])*1e3
            

        axisbc = np.array([np.min(self.ex),np.max(self.ex),\
                            np.min(self.ey),np.max(self.ey)])
            
        # Identify domain corners:
        #  P4 ----- P3
        #  |         |
        #  |         |
        #  |         |
        #  P1 ----- P2
        
        P1   = np.array([axisbc[0],axisbc[2]])
        P2   = np.array([axisbc[1],axisbc[2]])
        P3   = np.array([axisbc[1],axisbc[3]])
        P4   = np.array([axisbc[0],axisbc[3]])
        
        # Mesh nodes
        self.P1nod = np.where(np.sqrt(((self.coord[:,0]-P1[0]))**2 + \
                                  ((self.coord[:,1]-P1[1]))**2)<1e-12)[0][0]
        self.P2nod = np.where(np.sqrt(((self.coord[:,0]-P2[0]))**2 + \
                                  ((self.coord[:,1]-P2[1]))**2)<1e-12)[0][0]
        self.P3nod = np.where(np.sqrt(((self.coord[:,0]-P3[0]))**2 + \
                                  ((self.coord[:,1]-P3[1]))**2)<1e-12)[0][0]
        self.P4nod = np.where(np.sqrt(((self.coord[:,0]-P4[0]))**2 + \
                                  ((self.coord[:,1]-P4[1]))**2)<1e-12)[0][0]

            
    def load_LsStep(self):
        """Read indata from file."""
        
        cwd = os.getcwd()
        path = self.python_location + '/single_study/mat_files'
        
        # Find file 'level_set_X.mat' with highest index X in mat_files
        files = []
        for i in os.listdir(path):
            if os.path.isfile(os.path.join(path,i)) and 'level_set_' in i:
                files.append(i)
        max_number = 0
        for file in files:
            file_number = int(file.split('level_set_')[1].split('.mat')[0])
            if (file_number>max_number):
                max_number = file_number
            
        filename = path + '/level_set_' + str(max_number) + '.mat'
        
        self.version = max_number
        

        with open(filename, "r") as ifile:
            
            # Access variables from level set file
            data = scipy.io.loadmat(filename)
            
            # Global variables
            self.a              = data['a']
            self.ngrains        = np.shape(self.a)[1]
            self.line_ex_all    = data['line_ex']*1e3
            self.line_ey_all    = data['line_ey']*1e3
            self.line_seg_all   = data['line_seg'].flatten()
            self.material       = data['material'].flatten()
            self.ex             = np.transpose(data['newex'])*1e3
            self.ey             = np.transpose(data['newey'])*1e3
            self.coord          = np.transpose(data['newcoord'])*1e3
            self.enod           = np.transpose(data['enod']) - 1
        
            
        # Round line coords to 8 decimals 
        self.line_ex_all = self.line_ex_all.round(decimals=6)
        self.line_ey_all = self.line_ey_all.round(decimals=6)
        self.coord       = self.coord.round(decimals=6)
            
        # # Domain corners
        # self.P1 = self.coord[self.P1nod]
        # self.P2 = self.coord[self.P2nod]
        # self.P3 = self.coord[self.P3nod]
        # self.P4 = self.coord[self.P4nod]
        
        # All domain boundary nodes
        self.indNodBd, self.indElemBd, self.indLocalEdgBd, edges = self.boundary_nodes(self.coord, self.enod)
        self.indLocalEdgBd = np.asarray(self.indLocalEdgBd)
        self.indLocalEdgBd = self.indLocalEdgBd - 1
            
        
        # --- Extract all line coords ---
        line_coord = np.empty((0, 2), dtype=float)
        for grain in range(self.ngrains):
            
            # Local variables for the grain
            g_cols        = [2*grain, 2*grain + 1]
            nseg          = self.line_seg_all[grain]
            line_ex       = self.line_ex_all[:nseg,g_cols]
            line_ey       = self.line_ey_all[:nseg,g_cols]

            # Stack line coords
            line_coord    = np.vstack((line_coord,self.lines_to_coord(line_ex, line_ey)))
            
            
        # Unique line coords
        line_coord      = np.unique(line_coord, axis=0)
        self.line_coord = line_coord
            
    
        # --- Extract line connectivity matrix line_conn connecting the 
        #     line_coords---
        line_conn = np.empty((0, 2), dtype=int)
        for grain in range(self.ngrains):
            
            # Grain boundary
            g_cols        = [2*grain, 2*grain + 1]
            nseg          = self.line_seg_all[grain]
            line_ex       = self.line_ex_all[:nseg,g_cols]
            line_ey       = self.line_ey_all[:nseg,g_cols]
 
            # Line segment connectivity of coords
            line_conn     = np.vstack((line_conn, self.get_line_conn(self.line_coord, line_ex, line_ey)))
            
            
            
        # Domain corners
        # corner_coords = np.vstack((self.P1, self.P2, self.P3, self.P4))
        
        # Domain nodes
        boundary_coord      = self.coord[self.indNodBd,:]
        
        # --- Extract boundary connectivity ---
        nboundary_segs = np.size(self.indElemBd)
        boundary_conn = np.zeros([nboundary_segs,2],dtype=int)
        for iseg in range(nboundary_segs):
            k1=self.indLocalEdgBd[iseg]
            k2=k1+1
            if (k1==3):
                k2=0
            nod1 = self.enod[self.indElemBd[iseg],k1]
            nod2 = self.enod[self.indElemBd[iseg],k2]
            boundary_conn[iseg,:] = [nod1,nod2]
            
        
        # Loop thorugh all boundary segs and find which line_coord intersect boundary
        iBcross = -1*np.ones(nboundary_segs,dtype=int)
        iBcrosseg = []
        for iBseg in range(nboundary_segs):
            inodA  = boundary_conn[iBseg,0]
            inodB  = boundary_conn[iBseg,1]
            A      = boundary_coord[inodA,:]
            B      = boundary_coord[inodB,:]
            tol = 1e-8
            xmin = np.min([A[0],B[0]]) - tol
            xmax = np.max([A[0],B[0]]) + tol
            ymin = np.min([A[1],B[1]]) - tol
            ymax = np.max([A[1],B[1]]) + tol
            c1 = xmin < line_coord[:,0]
            c2 = line_coord[:,0] < xmax
            c3 = ymin < line_coord[:,1]
            c4 = line_coord[:,1] < ymax
            a = np.where(np.logical_and(np.logical_and(c1,c2),np.logical_and(c3,c4)))
            if (np.size(a)>0):
                iBcross[iBseg] = a[0]
                iBcrosseg.append(iBseg)
                
        s= 9
        
        # Split line segments
        iBseg = iBcrosseg[0]
        # for iBseg in iBcrosseg:
        inodA  = boundary_conn[iBseg,0]
        inodB  = boundary_conn[iBseg,1]
        A      = boundary_coord[inodA,:]
        B      = boundary_coord[inodB,:]
        addLineCoordIdx = iBcross[iBseg]
        P      = line_coord[addLineCoordIdx,:]
        
        # Add point
        boundary_coord = np.vstack([boundary_coord,P])
        inodP          = nboundary_segs
        
        # Change boundary_conn[iBseg,:] from AB to AP
        boundary_conn[iBseg,1] = inodP
        
        # Add a new connection from P to B
        boundary_conn = np.vstack([boundary_conn,np.array([inodP,inodB],dtype=int)])
        
        # Update number of boundary segs
        nboundary_segs = nboundary_segs + 1
            
     
        # All mesh_points (line_coord + boundary_coords). OBS order important
        self.mesh_points = np.vstack((self.line_coord, boundary_coord))
        
        # All line connectivity
        self.line_conn      = np.unique(line_conn, axis=0)
        self.boundary_conn  = boundary_conn
        self.boundary_coord = boundary_coord
        self.all_line_conn  = np.vstack((self.line_conn, boundary_conn + np.shape(self.line_coord)[0]))
        
   
                
    def get_line_conn(self, line_coord, line_ex, line_ey):
        
        # No. line segments
        nseg = np.shape(line_ex)[0]

        # Stack line_ex and line_ey into a single array for convenience
        lines_stacked = np.column_stack((line_ex.flatten(), line_ey.flatten()))
        
        # Find indices of line_segments in line_coord
        line_conn = np.where(np.all(line_coord == lines_stacked[:, None, :],axis=-1))[1]
        
        # Reshape to match line segments (nsegx2)
        line_conn = line_conn.reshape(nseg,2)
        
        return line_conn
    
    def lines_to_coord(self, line_ex, line_ey):
        """ Function for extracting unique vertices in lines"""
        N             = np.shape(line_ex)[0]
        coord         = np.zeros((2 * N, 2))
        coord[:N, 0]  = line_ex[:, 0]
        coord[N:, 0]  = line_ex[:, 1]
        coord[:N, 1]  = line_ey[:, 0]
        coord[N:, 1]  = line_ey[:, 1]
        unique_coord  = np.unique(coord, axis=0)
        return coord
    
    
    
    def boundary_nodes(self, nodes, elem):
        numNod = nodes.shape[0]
        numElem, ndim = elem.shape
    
        # Triangle edges
        if (ndim==3):
            edges = np.unique(np.sort(np.concatenate((elem[:, [0, 1]], elem[:, [1, 2]], elem[:, [2, 0]]), axis=0), axis=1), axis=0)
        elif (ndim==4):
            edges = np.unique(np.sort(np.concatenate((elem[:, [0, 1]], elem[:, [1, 2]], elem[:, [2, 3]],elem[:, [3, 0]]), axis=0), axis=1), axis=0)
        indNodBd      = []
        indLocalEdgBd = []
        indElemBd     = []
    
        # Look for the edges belonging only to one element
        for i in range(edges.shape[0]):
            n1, n2 = edges[i, 0], edges[i, 1]
            indRow, indCol = np.where(elem == n1)  # Find elements owning the first node
            indElem, col = np.where(elem[indRow, :] == n2)  # Owning also the second one
    
            if len(indElem) == 1:  # Boundary edges
                indNodBd.extend([n1, n2])
                indElemBd.append(indRow[indElem])
                lloc1 = np.where((elem[indRow[indElem], :] == n1)[0])[0][0]
                lloc2 = np.where((elem[indRow[indElem], :] == n2)[0])[0][0]
    
                if (ndim==3):
                    aux = np.array([0, 0, 0])
                    aux[lloc1] = 1
                    aux[lloc2] = 1
                    number = aux[0] + 2 * aux[1] + 4 * aux[2]
    
                    if number == 3:
                        edgeBd = 1
                    elif number == 5:
                        edgeBd = 3
                    elif number == 6:
                        edgeBd = 2
                    else:
                        raise ValueError('Edge not allowed')
                        
                if (ndim==4):
                    aux = np.array([0, 0, 0, 0])
                    aux[lloc1] = 1
                    aux[lloc2] = 1
                    number = aux[0] + 2 * aux[1] + 4 * aux[2] + 8 * aux[3]
                    if number == 3:
                        edgeBd = 1
                    elif number == 6:
                        edgeBd = 2
                    elif number == 9:
                        edgeBd = 4
                    elif number == 12:
                        edgeBd = 3
                    else:
                        raise ValueError('Edge not allowed')
    
                indLocalEdgBd.append(edgeBd)
    
        indNodBd = np.unique(indNodBd)
        return indNodBd, np.concatenate(indElemBd), indLocalEdgBd, edges
            
        
class OutputData(object):
    """Class for storing output data from calculation."""
    def __init__(self):
        
        # Python output
        self.coord           = None
        self.enod            = None
        self.bcnod           = None
        self.bcval           = None
        self.bcval_idx       = None
        self.dofs_per_node   = None
        self.nelm            = None
        self.nnod            = None
        self.ndof            = None
        self.nodel           = None
        self.dofel           = None
        self.bot_c           = None
        self.top_c           = None
        
        # Locations - not saved
        self.single_location = None
        self.param_locations = None
        self.input_location  = None

        # Fortran loadstep output
        self.vm              = None
        self.biax            = None
        self.i_imc           = None
        
    def save(self, filename):
        """Saving output data to dictionary."""
        
        output_data_file = {}
        output_data_file["coord"]           = np.transpose(self.coord).tolist()
        output_data_file["enod"]            = np.transpose(self.enod).tolist()
        # output_data_file["bcnod"]           = self.bcnod.tolist()
        # output_data_file["bcval"]           = self.bcval.tolist()
        # output_data_file["bcval_idx"]       = self.bcval_idx.tolist()
        output_data_file["dofs_per_node"]   = self.dofs_per_node
        output_data_file["nelm"]            = self.nelm
        output_data_file["nnod"]            = self.nnod
        output_data_file["ndof"]            = self.ndof
        output_data_file["nodel"]           = self.nodel
        output_data_file["dofel"]           = self.dofel
        
        
        with open(filename, "w") as ofile:
            json.dump(output_data_file, ofile, sort_keys = True, indent = 4)

    def loadmesh(self, filename):
        """Read output data from file."""

        with open(filename, "r") as ifile:
            output_data_file = json.load(ifile)
            
        self.coord         = np.transpose(np.asarray(output_data_file["coord"]))
        self.enod          = np.transpose(np.asarray(output_data_file["enod"]))
        self.edof          = np.asarray(output_data_file["edof"])
        self.ex            = np.asarray(output_data_file["ex"])
        self.ey            = np.asarray(output_data_file["ey"])
        self.bcnod         = output_data_file["bcnod"]
        self.bcval         = np.asarray(output_data_file["bcval"])
        self.dofs          = np.asarray(output_data_file["dofs"])
        self.dofs_per_node = output_data_file["dofs_per_node"]
        self.el_type       = np.asarray(output_data_file["el_type"])
        self.nelm          = np.asarray(output_data_file["nelm"])
        self.nnod          = output_data_file["nnod"]
        self.ndof          = np.asarray(output_data_file["ndof"])
        self.nodel         = np.asarray(output_data_file["nodel"])
        self.dofel         = output_data_file["dofel"]

        
        
class Mesh(object):
    """Class for generating the mesh"""
    def __init__(self, input_data, output_data):
        self.input_data  = input_data
        self.output_data = output_data
        
    def generateMesh(self):
            
        # for grain in range(self.input_data.ngrains):
        #     g_cols   = [2*grain, 2*grain + 1]
        #     line_seg = self.input_data.line_seg_all[grain]
        #     line_ex  = self.input_data.line_ex_all[:line_seg,g_cols]
        #     line_ey  = self.input_data.line_ey_all[:line_seg,g_cols]
        
        
        # --- gmsh mesh generation ---
        
        # # Initialize gmsh
        gmsh.initialize()
   
        # Create geometry
        triangulation = self.input_data.geometry()
        
        # # Location
        # print('Generating mesh in folder: ', os.getcwd())
        
        # Generate mesh
        gmsh.model.mesh.generate(2)
        
        # Graphical illustration of mesh
        if 'close' not in sys.argv:
            gmsh.fltk.run()
        
        # Finalize Gmsh
        gmsh.finalize()
  
        # Extract coord and enod
        coord = triangulation['vertices']
        enod  = triangulation['triangles']
        
        # Make sure that enod is numbered counter-clockwise for 3-node elm
        self.reorder_enod_triangle(coord, enod)
        
        # --- Topology quantities ---
        dofs_per_node = 1
        nelm  = np.shape(enod)[0]
        nnod  = np.shape(coord)[0]
        ndof  = nnod*dofs_per_node
        nodel = np.size(enod[0])
        dofel = nodel*dofs_per_node


        # --- Transfer model variables to output data ---
        self.output_data.coord          = coord
        self.output_data.enod           = enod + 1
        self.output_data.dofs_per_node  = dofs_per_node
        self.output_data.nelm           = nelm
        self.output_data.nnod           = nnod
        self.output_data.ndof           = ndof
        self.output_data.nodel          = nodel
        self.output_data.dofel          = dofel
        
        
    def extract_mesh_data(self,inpmatrix):
        # Convert mesh data from .inp file to float
        
        # Use numpy.char.split to split each string into a list of substrings
        split_arrays = np.char.split(inpmatrix, sep=', ')
        
        # Use numpy.vstack to stack the resulting lists into a 2D array
        float_array = np.vstack([np.array(arr, dtype=float) for arr in split_arrays])
        
        return float_array    
    
    
    def get_boundary_nodes(self,coord,enod,ndof,nelm,line_ex,line_ey):
        """Finding boundary nodes along interphase"""
        
        # Nseg
        nseg = np.size(line_ex,0)
        
        # Node coordinates along interphase
        line_coord = np.zeros([nseg+1,2])
        line_coord[1:-1]
        line_coord[:-1,0] = line_ex[:,0]
        line_coord[-1,0]  = line_ex[-1,1]
        line_coord[:-1,1] = line_ey[:,0]
        line_coord[-1,1]  = line_ey[-1,1]
        
        # Extract boundary nodes
        boundary_nodes = np.zeros(np.size(line_ex,0)+1,dtype=int)
        k = 0
        for i in range(ndof):
            # For each coord, compute distance to all line_coords.
            # If distance smaller than tol, add node to upper b nodes
            d    = np.sqrt(np.sum((coord[i,:]-line_coord)**2,1));
            tol  = 1e-12
            dsum = np.sum(d<tol)
            if (dsum==1):
                boundary_nodes[k] = i+1
                k = k+1
            elif(dsum>1):
                sys.exit("Error in finding interphase nodes")
        

        # Sort boundary nodes
        cc = coord[boundary_nodes-1,0]
        boundary_nodes = boundary_nodes[np.argsort(cc)]
        
        # Return boundary nodes
        return boundary_nodes
        
    
    def reorder_enod_triangle(self, coord, enod):
        """Routine for reordering enod in a plane triangle s.t. the node
        numbering is counter-clockwise"""
        
        nelm = np.size(enod,0)
        
        
        for i in range(nelm):
            # Extract the node indices for the current element
            elm_nodes = enod[i]
    
            # Extract the coordinates for the current element
            elm_coords = coord[elm_nodes]
    
            if (i==2037):
                s = 9
            # Check if the nodes are in counter-clockwise order, if not, reverse the order
            if not self.is_ccw(*elm_coords):
                enod[i] = np.flip(elm_nodes)

        return enod
        

    def is_ccw(self,p1, p2, p3):
        """
        Check if the given three points (nodes) are in counter-clockwise order.
        """
        return (p2[0]-p1[0])*(p3[1]-p2[1])-(p2[1]-p1[1])*(p3[0]-p2[0])>0

    
    def generateSingleMesh(self):
        """Generate a single mesh"""
        
        # --- Make new directory 'phase_meshes' in single_study if not exist
        #     and enter ---     
        cwd             = os.getcwd() # Save path
        location        = 'single_study/phase_meshes'
        path            = os.path.join(self.input_data.python_location, location)
        isDir           = os.path.isdir(path)
        if not isDir: 
            os.mkdir(path)
            print("Directory '% s' created" % location) 
        os.chdir(path)
  
        # --- Save location ---
        self.input_data.location = os.getcwd()

        # --- Create mesh ---
        # gmsh.initialize()
        # gmsh.finalize()
        
        # Generate mesh
        self.generateMesh()
          
        # --- Save output data as json file ---
        self.output_data.save('triangle_mesh_' + str(self.input_data.version) + '.json')
        
