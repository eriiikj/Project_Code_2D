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
        self.ccsep_log          = 0
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
        self.P1nod              = 0
        self.P2nod              = 0
        self.P3nod              = 0
        self.P4nod              = 0
        self.node_tags          = 0
        self.spline_tags        = 0
        self.line_coord         = 0
        self.line_conn          = 0
        self.mesh_points        = 0

        # Location
        self.location           = None
        self.python_location    = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/Python_output'
        self.fortran_location   = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/build/src/main'
        
  
    
     # --- Define the geometry ---
    
    def geometry(self):
        """Define geometry in gmsh"""
        
        
        # --- Constraint Delanuay triangularization using triangle ---
        tri_input     = dict(vertices=self.mesh_points,segments=self.line_conn) # segments=self.line_conn
        print('Creating a constraint delanuay triangularization...')
        triangulation = triangle.triangulate(tri_input,opts='pcq') #pcqa0.03 
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
        
        # Mesh
        input_data_file["el_size_x"]              = self.el_size_x
        input_data_file["el_size_y"]              = self.el_size_y
        
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
        
        
    def load_LsInit(self):
        """Read indata from file."""
        
        cwd = os.getcwd()
        path = self.python_location + '/single_study/mat_files'
        
        # Filename initial level set data
        filename = path + '/level_setinit.mat'
        
        

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
            self.line_ex_all    = (data['line_ex']*1e3)
            self.line_ey_all    = data['line_ey']*1e3
            self.line_seg_all   = data['line_seg'].flatten()
            self.material       = data['material'].flatten()
            self.ex             = np.transpose(data['newex'])*1e3
            self.ey             = np.transpose(data['newey'])*1e3
            self.coord          = np.transpose(data['newcoord'])*1e3
            self.enod           = np.transpose(data['enod']) - 1
        
            
        # Round line coords to 6 decimals 
        self.line_ex_all = self.line_ex_all.round(decimals=6)
        self.line_ey_all = self.line_ey_all.round(decimals=6)
            
        # Domain corners
        self.P1   = self.coord[self.P1nod]
        self.P2   = self.coord[self.P2nod]
        self.P3   = self.coord[self.P3nod]
        self.P4   = self.coord[self.P4nod]
            
        
        
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
        self.line_coord = np.unique(line_coord, axis=0)
            
        # Domain corners
        corner_coords = np.vstack((self.P1, self.P2, self.P3, self.P4))
        
        # All mesh_points (line_coord + domain corners)
        self.mesh_points = np.unique(np.vstack((self.line_coord, corner_coords)), axis=0)
        
        
        
        
        # --- Extract line connectivity matrix line_conn ---
        line_conn = np.empty((0, 2), dtype=int)
        for grain in range(self.ngrains):
            
            # Grain boundary
            g_cols        = [2*grain, 2*grain + 1]
            nseg          = self.line_seg_all[grain]
            line_ex       = self.line_ex_all[:nseg,g_cols]
            line_ey       = self.line_ey_all[:nseg,g_cols]
 
            # Line segment connectivity of coords
            line_conn     = np.vstack((line_conn, self.get_line_conn(self.mesh_points, line_ex, line_ey)))
            
        # Unique line connectivity
        self.line_conn = np.unique(line_conn, axis=0)
                
                
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
        
        # Initialize gmsh
        gmsh.initialize()
   
        # Create gmsh geometry
        triangulation = self.input_data.geometry()
        
        # Location
        print('Generating mesh in folder: ', os.getcwd())
        
        # Generate mesh
        gmsh.model.mesh.generate(2)
        
        # Graphical illustration of mesh
        # if 'close' not in sys.argv:
        #     gmsh.fltk.run()
        
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
            d = np.sqrt(np.sum((coord[i,:]-line_coord)**2,1));
            tol = 1e-12
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
    
            # Check if the nodes are in counter-clockwise order, if not, reverse the order
            if self.is_ccw(*elm_coords):
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
        gmsh.initialize()
        gmsh.finalize()
        
        # Generate mesh
        self.generateMesh()
          
        # --- Save output data as json file ---
        self.output_data.save('triangle_mesh_' + str(self.input_data.version) + '.json')
        
        # # # --- Export location to Fortran ---
        # # self.exportLocationToFortran()
        
        # # --- Leave dir ---
        # os.chdir(cwd)
        
        
 
    def exportLocationToFortran(self):
        """Export location to Fortran program."""
        
        input_data_location_file                   = {}
        input_data_location_file["input_location"] = self.input_data.location

        exportFileName = 'python_phase_input_location.json'
        
        # --- Enter Fortran build ---
        current_dir      = os.getcwd()
        abs_fortran_path = self.input_data.fortran_location
        isDir            = os.path.isdir(abs_fortran_path) 
        if isDir: 
            os.chdir(abs_fortran_path)
        else:
            print("Fortran export input location not found.")

        print("Entering %s and specifying location %s." % (abs_fortran_path, self.input_data.location))
        with open(exportFileName, "w") as ofile:
            json.dump(input_data_location_file, ofile, sort_keys = True, indent = 4)
            
        # --- Go out of Fortran build ---
        os.chdir(current_dir)

    def makeDirinCurrent(self,dir_name):
        """Make dir and enter"""
        current_dir = os.getcwd()
        path        = os.path.join(current_dir, dir_name) 
        isDir       = os.path.isdir(path) 
        if isDir: shutil.rmtree(path)
        
        os.mkdir(path) 
        print("Directory '% s' created" % dir_name) 
        os.chdir(path)
