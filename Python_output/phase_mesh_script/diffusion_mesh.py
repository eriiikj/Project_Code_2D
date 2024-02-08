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

from contextlib import redirect_stdout
import io

class InputData(object):
    """Class for defining input data to our model."""
    def __init__(self):
        
        # --- Define private variabels, standrad values

        # Version
        self.version                = 1
        
        # Geometry
        self.line_ex                = 0
        self.line_ey                = 0
        self.P0                     = 0
        self.P1                     = 0
        self.P2                     = 0
        self.P3                     = 0
        
        # Mesh
        self.el_size_x              = 25e-3
        self.el_size_y              = 25e-3
        self.el_size_factor         = 0.2
        

        # Location
        self.location               = None
        
        
        # Markers
        self.lowerSupportMarkers    = None 
        self.upperSupportMarkers    = None 
  
    
     # --- Define the geometry ---
    
    def geometry(self):
        """Create geometry instance"""

    
        # Short form
        ngrains = self.ngrains
        
        for grain in range(ngrains):
            g_cols     = [2*grain, 2*grain + 1]
            nseg       = self.line_seg_all[grain]
            line_ex    = self.line_ex_all[:nseg,g_cols]
            line_ey    = self.line_ey_all[:nseg,g_cols]
            nsep_lines = self.nsep_lines_all[grain]
            sep_lines  = self.sep_lines_all[:nsep_lines,g_cols]
            
            
        
        
        
        
            # --- Add all line segement points ---
            npoint       = 0 
            line_counter = 0
            # Loop through all seperate lines
            for sep_idx in range(nsep_lines):
                
                # Cols of sep_idx
                sep_idx_cols = [2*sep_idx, 2*sep_idx + 1]
                
                # Segements to add
                nseg_sep     = sep_lines[sep_idx,1] - sep_lines[sep_idx,0] + 1
                # nseg_sep_add = line_seg_add[sep_idx]
      
                # Add points in primary interface
                for i in range(nseg_sep):
                    P = [line_ex[line_counter,0],line_ey[line_counter,0]]
                    line_counter = line_counter + 1
                    
                    # Add point
                    gmsh.model.geo.add_point(P[0], P[1], 0, tag=npoint)
                    # g.point(P, npoint)
                    npoint = npoint + 1
                    
                # Add last point if not closed curve
                start_point = np.array([line_ex[0,0], line_ey[0,0]])
                end_point   = np.array([line_ex[line_counter-1,1],line_ey[line_counter-1,1]])
                if (np.linalg.norm(start_point - end_point)>1e-15): # Closed curve
                    P = [line_ex[line_counter-1,1],line_ey[line_counter-1,1]]
                    gmsh.model.geo.add_point(P[0], P[1], 0, tag=npoint)
                    npoint = npoint + 1
                
                # # Add points in added interface
                # for i in range(nseg_sep_add-1):
                #     P = [line_ex_add[i,sep_idx_cols[1]],line_ey_add[i,sep_idx_cols[1]]]
                    
                #     # Add point
                #     gmsh.model.geo.add_point(P[0], P[1], 0, tag=npoint)
                #     # g.point(P, npoint)
                #     npoint = npoint + 1
        
    
    
    
            # --- Splines ---
            myid = 0
            
            # Markers
            self.interface_marker = 10
            
            
            # Loop through all seperate lines
            npoint_prev = 0
            interface_splines = np.zeros(nseg,dtype=int)
            marker_c = 0
            
            # for sep_idx in range(nsep_lines):
                
            #     # Cols of sep_idx
            #     sep_idx_cols = [2*sep_idx, 2*sep_idx + 1]
            #     nseg_sep     = sep_lines[sep_idx,1] - sep_lines[sep_idx,0] + 1
            #     # nseg_sep_add = line_seg_add[sep_idx]
                
            #     start_point = np.array([line_ex[sep_lines[sep_idx,0]-1,0], line_ey[sep_lines[sep_idx,0]-1,0]])
            #     end_point   = np.array([line_ex[sep_lines[sep_idx,1]-1,1], line_ey[sep_lines[sep_idx,1]-1,1]])
                    
                        
            #     # Interface splines
            #     for i in range(nseg_sep):
            #         left_point  = myid
            #         right_point = myid + 1
            #         # if (np.linalg.norm(start_point - end_point)<1e-15):
            #             if (sep_idx==nsep_lines-1 and nseg_sep_add==0 and i==nseg_sep-1):
            #                 left_point  = myid
            #                 right_point = 0
            #         # g.spline([left_point, right_point], ID=myid, el_on_curve=1)
            #         gmsh.model.geo.add_line(left_point, right_point,tag=myid)
            #         interface_splines[marker_c] = myid
            #         marker_c = marker_c + 1
            #         myid = myid + 1
                
                    
            #     # # Additional splines
            #     # for i in range(nseg_sep_add):
            #     #     left_point  = myid
            #     #     right_point = myid + 1
            #     #     if (sep_idx==nsep_lines-1 and i==nseg_sep_add-1):
            #     #         left_point  = myid
            #     #         right_point = 0
                        
        
            #         # g.spline([left_point, right_point], ID=myid)
            #         gmsh.model.geo.add_line(left_point, right_point,tag=myid)
            #         myid = myid + 1
         
            #     # npoint_prev = npoint_prev + nseg_sep + nseg_sep_add
    
            
            # --- Add unstructured surface ---
            splines = list(range(0,nseg_tot))
            splines = list(range(0,nseg_tot))
            face1 = gmsh.model.geo.add_curve_loop(splines)
            gmsh.model.geo.add_plane_surface([face1])
            
            gmsh.model.geo.addPhysicalGroup(dim=1, tags=interface_splines)
           
            
            # Create the relevant Gmsh data structures 
            # from Gmsh model.
            gmsh.model.geo.synchronize()
    
            
            return interface_splines
    
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
        
        # Geometry
        self.w                = input_data_file["w"]
        self.h                = input_data_file["h"]
        self.y2               = input_data_file["y2"]
        self.y3               = input_data_file["y3"]
        self.y4               = input_data_file["y4"]
        self.y5               = input_data_file["y5"]
        
        # Mesh
        self.el_size_x        = input_data_file["el_size_x"]
        self.el_size_y        = input_data_file["el_size_y"]
        
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
            

        self.axisbc = np.array([np.min(self.ex),np.max(self.ex),\
                                np.min(self.ey),np.max(self.ey)])
            
        axisbc = self.axisbc
        
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
            self.a              = data['a']*1e3
            self.line_seg       = data['line_seg'].flatten()
            self.line_seg_all   = self.line_seg 
            self.ngrains        = np.size(self.line_seg)
            self.line_ex_all    = data['line_ex']*1e3
            self.line_ey_all    = data['line_ey']*1e3
            self.sep_lines_all  = data['sep_lines']
            self.nsep_lines_all = np.sum(self.sep_lines_all[:,0:-1:2]>0,0)
            self.material       = data['material'].flatten()
            
            
            # IMC area - ?
            self.tot_imc_area_init = data['IMC_area_init'].flatten()[0]
            self.tot_imc_area      = data['IMC_area'].flatten()[0]
            
            if (self.tot_imc_area_init==0):
                self.sn_c0 = self.sn_c0_fix
            else:
                self.sn_c0 = (self.sn_c0_fix - self.sn_imc)*\
                    (self.tot_imc_area_init/self.tot_imc_area)\
                    +self.sn_imc
            
            left_to_right_prev = False
                    
                    
            # Local variables for the grain
            # Short form 
            grain     = self.grain
            g_cols    = [2*grain, 2*grain + 1]
            
            self.line_seg     = self.line_seg[grain]
            self.line_ex      = self.line_ex_all[:self.line_seg,g_cols]
            self.line_ey      = self.line_ey_all[:self.line_seg,g_cols]
            self.nsep_lines   = self.nsep_lines_all[self.grain]
            self.sep_lines    = self.sep_lines_all[:self.nsep_lines,g_cols]
                    
            
            self.el_size_x  = data['elsize_x'][0][0]*1e3
            self.el_size_y  = data['elsize_y'][0][0]*1e3
            
            # Short form 
            line_ex    = self.line_ex
            line_ey    = self.line_ey
            nseg       = self.line_seg
            sep_lines  = self.sep_lines
            nsep_lines = self.nsep_lines
            ngrains    = self.ngrains
            axisbc     = self.axisbc
            
            # Extra lines
            self.line_ex_add  = np.zeros([nseg*5,nsep_lines*2])
            self.line_ey_add  = np.zeros([nseg*5,nsep_lines*2])
            self.line_seg_add = np.zeros(nsep_lines,dtype=int)

            
            # Access the variable by its name
            self.ex     = np.transpose(data['newex'])*1e3
            self.ey     = np.transpose(data['newey'])*1e3
            self.coord  = np.transpose(data['newcoord'])*1e3
            
             
            # Identify points:
            #  P4 ----- P3
            #  |         |
            #  |         |
            #  |         |
            #  P1 ----- P2

            self.P1   = self.coord[self.P1nod]
            self.P2   = self.coord[self.P2nod]
            self.P3   = self.coord[self.P3nod]
            self.P4   = self.coord[self.P4nod]

        
            # Short form 
            P1        = self.P1
            P2        = self.P2
            P3        = self.P3
            P4        = self.P4
    
            for sep_idx in range(self.nsep_lines):
                
                # Cols of sep_idx
                sep_idx_cols = [2*sep_idx, 2*sep_idx + 1]
                
                # Start point
                iseg_start = sep_lines[sep_idx, 0] - 1 # 0 idx
                A = np.array([line_ex[iseg_start,0], \
                              line_ey[iseg_start,0]])
                
                # End point
                iseg_end = sep_lines[sep_idx, 1] - 1 # 0 idx
                B = np.array([line_ex[iseg_end,1], \
                              line_ey[iseg_end,1]])
                
                if (np.linalg.norm(A-B)>1e-13):
                    # A and B does not coincide. Both A and B located at 
                    # domain edge
                    
                    # Find point that should be added 
                    # Borders:
                    #  ---- 4 ----
                    #  |         |
                    #  1         3
                    #  |         |
                    #  ---- 2 ---
                    

                    # 1) Find domain border crossed by A
                    if (abs(A[0]-axisbc[0])<1e-13):
                        A_border = 1
                    elif (abs(A[1]-axisbc[2])<1e-13):
                        A_border = 2
                    elif (abs(A[0]-axisbc[1])<1e-13):
                        A_border = 3
                    elif (abs(A[1]-axisbc[3])<1e-13):
                        A_border = 4
                        
                    # 1) Find domain border crossed by B
                    if (abs(B[0]-axisbc[0])<1e-13):
                        B_border = 1
                    elif (abs(B[1]-axisbc[2])<1e-13):
                        B_border = 2
                    elif (abs(B[0]-axisbc[1])<1e-13):
                        B_border = 3
                    elif (abs(B[1]-axisbc[3])<1e-13):
                        B_border = 4
                    
                    if (A_border==B_border):
                        
                        # Add line segment from B to A (which goes along
                        # border) last in line list such that curve
                        # becomes closed
                        added_lines = self.line_seg_add[sep_idx]
                        self.line_ex_add[added_lines,sep_idx_cols] = [B[0], A[0]]
                        self.line_ey_add[added_lines,sep_idx_cols] = [B[1], A[1]]
                        self.line_seg_add[sep_idx] = self.line_seg_add[sep_idx] + 1
 
                        
                    else:
                            
                        
                        dl  = False
                        dl2 = False
                        
                        # Find closest lower node to A
                        closest_lower_node_A = np.argmin(np.sqrt(\
                        (self.coord[:,0]-A[0])**2 + \
                        (self.coord[:,1]-(A[1]-self.el_size_y))**2))
                            
                        # Find closest left node to A
                        closest_left_node_A = np.argmin(np.sqrt(\
                        (self.coord[:,0]-(A[0]-self.el_size_x))**2 + \
                        (self.coord[:,1]-A[1])**2))
                        
                        # A_border = 1
                        if (A_border==1):
                            if (B_border==2):
                                P = P1
                            elif (B_border==3):
                                
                                
                                if (self.nsep_lines==1):
                                    dl = True
                                    
                                    if (self.a[closest_lower_node_A,grain]<0):
                                        if (A[0]<B[0]): 
                                            P = [P2, P1]
                                        else:
                                            P = [P1, P2]
                                    else:
                                        if (A[0]<B[0]): 
                                            P = [P3, P4]
                                        else:
                                            P = [P4, P2]
                                elif(self.nsep_lines==2):
                                    dl2 = True
                                    if (sep_idx==0):
                                        next_sep = 1
                                    elif(sep_idx==1):
                                        next_sep = 0
                                    
                                    if (A[0]<B[0]):
                                        left_to_right=True
                                    else:
                                        left_to_right=False
                                    
                                    
                                    
                                    if (left_to_right==left_to_right_prev and sep_idx>0):
                                        # Reshuffle lines
                                        line_ex_copy = np.copy(line_ex[iseg_start:iseg_end+1,:])
                                        line_ex[iseg_start:iseg_end+1,0] = np.flip(line_ex_copy[:,1])
                                        line_ex[iseg_start:iseg_end+1,1] = np.flip(line_ex_copy[:,0])
                                        
                                        line_ey_copy = np.copy(line_ey[iseg_start:iseg_end+1,:])
                                        line_ey[iseg_start:iseg_end+1,0] = np.flip(line_ey_copy[:,1])
                                        line_ey[iseg_start:iseg_end+1,1] = np.flip(line_ey_copy[:,0])

                                        
                                    # Start point
                                    iseg_start = sep_lines[sep_idx, 0] - 1 # 0 idx
                                    A = np.array([line_ex[iseg_start,0], \
                                                  line_ey[iseg_start,0]])
                                    
                                    # End point
                                    iseg_end = sep_lines[sep_idx, 1] - 1 # 0 idx
                                    B = np.array([line_ex[iseg_end,1], \
                                                  line_ey[iseg_end,1]])  
                                        
                                        
                                    iseg_start2 = sep_lines[next_sep, 0] - 1
                                    C = np.array([line_ex[iseg_start2,0], \
                                                  line_ey[iseg_start2,0]])
                                    iseg_end2 = sep_lines[next_sep, 1] - 1 # 0 idx
                                    D = np.array([line_ex[iseg_end2,1], \
                                                  line_ey[iseg_end2,1]])
                                    
                                    
                                    if (D[0]==B[0]):
                                        P = D
                                    else:
                                        P = C
                                    left_to_right_prev = left_to_right    
                                    print('Assuming sep lines can be added.')
                                    
                                  
                                            
                                    
                                
                            elif (B_border==4):
                                P = P4
                                
          
                        # A_border = 2
                        if (A_border==2):
                            if (B_border==1):
                                P = P1
                            elif (B_border==3):
                                P = P3
                            elif (B_border==4):
                                dl = True
                                
                                if (self.a[closest_left_node_A,grain]<0):
                                    if (A[1]<B[1]): 
                                        P = [P4, P1]
                                    else:
                                        P = [P1, P4]
                                else:
                                    if (A[1]<B[1]): 
                                        P = [P3, P2]
                                    else:
                                        P = [P2, P3]
                        
                        # A_border = 3
                        if (A_border==3):
                            if (B_border==1):
                                
                                
                                if (self.nsep_lines==1):
                                    dl = True
                                    if (self.a[closest_lower_node_A,grain]<0):
                                        if (A[0]<B[0]): 
                                            P = [P2, P1]
                                        else:
                                            P = [P1, P2]
                                    else:
                                        if (A[0]<B[0]): 
                                            P = [P3, P4]
                                        else:
                                            P = [P4, P2]
                                else:
                                    print('Assuming sep lines can be added.')
                                
                                
                            elif (B_border==2):
                                P = P2
                            elif (B_border==4):
                                P = P3
                                
                                        
                        
                        # A_border = 4
                        if (A_border==4):
                            if (B_border==1):
                                P = P4
                            elif (B_border==2):
                                dl = True
                                if (self.a[closest_left_node_A,grain]<0):
                                    if (A[1]<B[1]): 
                                        P = [P4, P1]
                                    else:
                                        P = [P1, P4]
                                else:
                                    if (A[1]<B[1]): 
                                        P = [P3, P2]
                                    else:
                                        P = [P2, P3]
                            elif (B_border==3):
                                P = P3
                                
                                
                        if (dl):
                            PP1 = P[0]; PP2 = P[1]
                            
                            
                            
                            # Draw line from B to PP1
                            added_lines = self.line_seg_add[sep_idx]
                            self.line_ex_add[added_lines,sep_idx_cols] = [B[0], PP1[0]]
                            self.line_ey_add[added_lines,sep_idx_cols] = [B[1], PP1[1]]
                            self.line_seg_add[sep_idx] = self.line_seg_add[sep_idx] + 1
                            
                            # Draw line from PP1 to PP2
                            added_lines = self.line_seg_add[sep_idx]
                            self.line_ex_add[added_lines,sep_idx_cols] = [PP1[0], PP2[0]]
                            self.line_ey_add[added_lines,sep_idx_cols] = [PP1[1], PP2[1]]
                            self.line_seg_add[sep_idx] = self.line_seg_add[sep_idx] + 1
                            
                            
                            # Draw line from PP2 to A
                            added_lines = self.line_seg_add[sep_idx]
                            self.line_ex_add[added_lines,sep_idx_cols] = [PP2[0], A[0]]
                            self.line_ey_add[added_lines,sep_idx_cols] = [PP2[1], A[1]]
                            self.line_seg_add[sep_idx] = self.line_seg_add[sep_idx] + 1
 
                            
                        elif(dl2):
                            
                            # Draw line from B to P
                            added_lines = self.line_seg_add[sep_idx]
                            self.line_ex_add[added_lines,sep_idx_cols] = [B[0], P[0]]
                            self.line_ey_add[added_lines,sep_idx_cols] = [B[1], P[1]]
                            self.line_seg_add[sep_idx] = self.line_seg_add[sep_idx] + 1
                        else:
                                
                        
                            # Draw line from B to P
                            added_lines = self.line_seg_add[sep_idx]
                            self.line_ex_add[added_lines,sep_idx_cols] = [B[0], P[0]]
                            self.line_ey_add[added_lines,sep_idx_cols] = [B[1], P[1]]
                            self.line_seg_add[sep_idx] = self.line_seg_add[sep_idx] + 1
                            
                            # Draw line from P to A
                            added_lines = self.line_seg_add[sep_idx]
                            self.line_ex_add[added_lines,sep_idx_cols] = [P[0], A[0]]
                            self.line_ey_add[added_lines,sep_idx_cols] = [P[1], A[1]]
                            self.line_seg_add[sep_idx] = self.line_seg_add[sep_idx] + 1
 
    
 
    
        
class OutputData(object):
    """Class for storing output data from calculation."""
    def __init__(self):
        
        # Python output
        self.coord           = None
        self.enod            = None
        self.edof            = None
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
        output_data_file["edof"]            = self.edof.tolist()
        output_data_file["bcnod"]           = self.bcnod.tolist()
        output_data_file["bcval"]           = self.bcval.tolist()
        output_data_file["bcval_idx"]       = self.bcval_idx.tolist()
        output_data_file["dofs_per_node"]   = self.dofs_per_node
        output_data_file["nelm"]            = self.nelm
        output_data_file["nnod"]            = self.nnod
        output_data_file["ndof"]            = self.ndof
        output_data_file["nodel"]           = self.nodel
        output_data_file["dofel"]           = self.dofel
        output_data_file["input_location"]  = self.input_location
        
        
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

        
        
class Solver(object):
    """Class for generating the mesh and executing Fortran code"""
    def __init__(self, input_data, output_data):
        self.input_data  = input_data
        self.output_data = output_data
        
    def generateMesh(self):
            
        # Short form 
        grain        = self.input_data.grain
        line_ex      = self.input_data.line_ex
        line_ey      = self.input_data.line_ey
        nseg         = self.input_data.line_seg
        sep_lines    = self.input_data.sep_lines
        nsep_lines   = self.input_data.nsep_lines
        ngrains      = self.input_data.ngrains
        axisbc       = self.input_data.axisbc
        line_ex_add  = self.input_data.line_ex_add
        line_ey_add  = self.input_data.line_ey_add
        line_seg_add = self.input_data.line_seg_add
        nseg_tot     = np.sum(line_seg_add) + nseg
        
        
        # grain_idx
        grain_idx = np.linspace(0,ngrains-1,ngrains,dtype=int)
        
        
        for grain in range(self.input_data.ngrains):
            g_cols   = [2*grain, 2*grain + 1]
            line_seg = self.input_data.line_seg_all[grain]
            line_ex  = self.input_data.line_ex_all[:line_seg,g_cols]
            line_ey  = self.input_data.line_ey_all[:line_seg,g_cols]
    
        # # g_cols
        # g_cols         = [2*grain, 2*grain + 1]
        # grain_material = self.input_data.material[grain] - 1

        # Create geometry
        interface_splines = self.input_data.geometry()
        
        # Change location to phase_meshes
        # os.chdir(self.input_data.python_location + '/single_study/phase_meshes')
        print('where', os.getcwd())
        print('where', self.input_data.python_location)
        
        # Generate mesh:
        gmsh.model.mesh.generate(2) # 2D mesh generation
        
        # Write mesh data:
        filename = 'phase_mesh_' + str(self.input_data.version)
        # gmsh.write(filename s+ ".msh")
        gmsh.write(filename + ".inp")
        
        # --- Graphical illustration of mesh ---
        # if 'close' not in sys.argv:
        #     gmsh.fltk.run()
            
        
        # --- Extract coord and enod in the model ---
        entities = gmsh.model.getEntities()
        
        # --- Extract bcnod ---
        bcnods = []
        for line_tag in interface_splines:
            line_nods,_,_ = gmsh.model.mesh.getNodes(dim=1,tag=line_tag,includeBoundary=True)
            bcnods.extend(line_nods)
        bcnods = sorted(set(bcnods))
        bcnods = np.asarray(bcnods,dtype=int)
        bcval  = np.zeros(np.size(bcnods))
        
        # Number of bcnods
        n_interface_points = np.size(bcnods)
        
        with open(filename + '.inp') as infile:
            lines = [ln.strip() for ln in infile.readlines()]
            # Remove comments
            lines = [ln for ln in lines if not ln.startswith("**")]
            # Find section headers
            headings = [(ln[1:], n) for n, ln in enumerate(lines) if ln.startswith("*")]
            # Filter the headings so that every heading has a start-of-data and
            # end-of-data index.
            headings.append(("end", -1))
            ln = [h[1] for h in headings]
            headings = [
                (name, start + 1, end) for (name, start), end in zip(headings[:-1], ln[1:])
            ]
        
        for h in headings:
            name       = h[0]
            lname      = name.lower()
            start_line = h[1]
            end_line   = h[2]
            
            # Coord
            if lname.startswith("node"):
                coord = lines[start_line:end_line]
                
            # Connectivity
            if lname.startswith("element, type=cps3"):
                enod = lines[start_line:end_line]
                
                
        # Extract coord
        coord = self.extract_mesh_data(coord)
        coord = coord[:,1:3]
        
        # Extract enod
        enod  = self.extract_mesh_data(enod)
        enod  = enod[:,1:]
        enod  = enod.astype(int)
        edof  = enod
        
        # Make sure that enod is numbered counter-clockwise for 3-node elm
        self.reorder_enod_triangle(coord, enod)
        
        # Dofs per node
        dofs_per_node = 1

        # Line ex all
        line_ex_all  = self.input_data.line_ex_all
        line_ey_all  = self.input_data.line_ey_all
        line_seg_all = self.input_data.line_seg_all
        
        
        
        # --- Equilibrium concentrations of all materials ---
        cu_cu   = self.input_data.cu_cu
        cu_imc  = self.input_data.cu_imc
        cu_sn   = self.input_data.cu_sn
        imc_cu  = self.input_data.imc_cu
        imc_imc = self.input_data.imc_imc
        imc_sn  = self.input_data.imc_sn
        sn_cu   = self.input_data.sn_cu
        sn_imc  = self.input_data.sn_imc
        sn_sn   = self.input_data.sn_sn
        sn_c0   = self.input_data.sn_c0
        sn_c0_fix = self.input_data.sn_c0_fix
        
        
        # Short form
        cu_params     = self.input_data.cu_params
        imc_params    = self.input_data.imc_params
        sn_params     = self.input_data.sn_params
        molar_volumes = self.input_data.molar_volumes
        
        
        # Diffusion potentials
        mu_cu_cu   = cu_params[0]*(cu_cu  - cu_params[3]) + cu_params[1]
        mu_cu_imc  = cu_params[0]*(cu_imc - cu_params[3]) + cu_params[1]
        mu_cu_sn   = cu_params[0]*(cu_sn  - cu_params[3]) + cu_params[1]
        
        mu_imc_cu  = imc_params[0]*(imc_cu  - imc_params[3]) + imc_params[1]
        mu_imc_imc = imc_params[0]*(imc_imc - imc_params[3]) + imc_params[1]
        mu_imc_sn  = imc_params[0]*(imc_sn  - imc_params[3]) + imc_params[1]
        
        mu_sn_cu   = sn_params[0]*(sn_cu  - sn_params[3]) + sn_params[1]
        mu_sn_imc  = sn_params[0]*(sn_imc - sn_params[3]) + sn_params[1]
        mu_sn_sn   = sn_params[0]*(sn_sn  - sn_params[3]) + sn_params[1]
        
        # Equilibrium compositions
        # Row 1: Cu eq comps
        # Row 2: IMC eq comps
        # Row 3: Sn eq comps
        eq_comp = np.array([[mu_cu_cu , mu_cu_imc , mu_cu_sn],\
                            [mu_imc_cu, mu_imc_imc, mu_imc_sn],\
                            [mu_sn_cu , mu_sn_imc , mu_sn_sn]])
        
        
        grain_material = self.input_data.material[grain] - 1
        
        # Save triple junction points
        tp_points_x = []
        tp_points_y = []
        
        # Determine what other grain has smallest value at that node
        other_grains = grain_idx[grain_idx != grain]
        
        # Add bcval
        line_counter  = 0
        point_counter = 0
        for sep_idx in range(nsep_lines):
            
            # Cols of sep_idx
            sep_idx_cols = [2*sep_idx, 2*sep_idx + 1]
            
            # Number of segments lines in separate line
            nseg_sep     = sep_lines[sep_idx,1] - sep_lines[sep_idx,0] + 1
            line_counter_start = line_counter
            for k in range(nseg_sep):
                
                # Point P
                P = [line_ex[line_counter][0],line_ey[line_counter][0]]
                
                # See if point P is in other grain
                P_in_g = []
                for g in other_grains:
                    g_cols = [2*g, 2*g + 1]
                    line_ex_g = line_ex_all[:line_seg_all[g],[g_cols[0],g_cols[1]]]
                    line_ey_g = line_ey_all[:line_seg_all[g],[g_cols[0],g_cols[1]]]
                    if ((abs((line_ex_g - P[0]) + (line_ey_g - P[1]))<1e-15).any()):
                        P_in_g.append(g)
                
                P_in_g = np.array(P_in_g)
                    
                if (np.size(P_in_g)==2):
                    # triple junction
                    tp_points_x.append(P[0])
                    tp_points_y.append(P[1])
                    
                if (np.size(P_in_g)>1):
                    other_materials = self.input_data.material[P_in_g] - 1
                    if (any((other_materials-grain_material) == 0)):
                        P_in_g = P_in_g[other_materials != grain_material]
                    else:
                        P_in_g = P_in_g[0]
                            
                other_material = self.input_data.material[P_in_g] - 1

                        
                if (np.size(other_material)>1):
                    print('Error in bcval')

                # Find equilibrium composition grain, lowest_grain
                bcval[point_counter] = eq_comp[grain_material,other_material]
                
                # Update counter
                line_counter  = line_counter  + 1
                point_counter = point_counter + 1
     
        
     
            # Check if first and last point of current separate segment coincide
            A = np.array([line_ex[line_counter_start][0],line_ey[line_counter_start][0]])
            B = np.array([line_ex[line_counter-1][1],line_ey[line_counter-1][1]])
            
            if (np.linalg.norm(A-B)>1e-12):
                # If points dont coincide we must add one more bcval
                
                # Point P
                P = B
                
                
                # See if point P is in other grain
                P_in_g = []
                for g in other_grains:
                    g_cols = [2*g, 2*g + 1]
                    line_ex_g = line_ex_all[:line_seg_all[g],[g_cols[0],g_cols[1]]]
                    line_ey_g = line_ey_all[:line_seg_all[g],[g_cols[0],g_cols[1]]]
                    if ((abs((line_ex_g - P[0]) + (line_ey_g - P[1]))<1e-15).any()):
                        P_in_g.append(g)
                
                
                P_in_g = np.array(P_in_g)
                
                if (np.size(P_in_g)==2):
                    # triple junction
                    tp_points_x.append(P[0])
                    tp_points_y.append(P[1])
                    
                if (np.size(P_in_g)>1):
                    other_materials = self.input_data.material[P_in_g] - 1
                    if (any((other_materials-grain_material) == 0)):
                        P_in_g = P_in_g[other_materials != grain_material]
                    else:
                        P_in_g = P_in_g[0]
                            
                other_material = self.input_data.material[P_in_g] - 1
                
                
                if (np.size(other_material)>1):
                    print('Error in bcval')
                
                bcval[point_counter] = eq_comp[grain_material,other_material]
                point_counter = point_counter + 1
                
            
        # Rearrange bcnod
        nbc        = np.size(bcnods)
        bcnod      = np.zeros((nbc,2), dtype=np.int32)
        bcnod[:,0] = bcnods
        bcnod[:,1] = np.ones((nbc), dtype=np.int32)
        
        
        # Create a boolean mask and remove imc_imc values from bc
        # mask  = np.abs(bcval - self.input_data.imc_imc)>1e-12
        mask  = np.abs(bcval - mu_imc_imc) > 1e-12
        bcval = bcval[mask]
        bcnod = bcnod[mask]
        
        
        # Create a boolean mask and remove sn_sn values from bc
        # mask  = np.abs(bcval - self.input_data.imc_imc)>1e-12
        mask  = np.abs(bcval - mu_sn_sn)>1e-12
        bcval = bcval[mask]
        bcnod = bcnod[mask]
        
        
        # Remove tp points from bc
        tp_points_x = np.array(tp_points_x)
        tp_points_y = np.array(tp_points_y)
        ntp = np.size(tp_points_x)
        mask = np.full(np.size(bcnod[:,0]), True, dtype=bool)
        for (xtp,ytp) in zip(tp_points_x,tp_points_y):
            grain_nod = np.where(np.logical_and(abs(coord[bcnod[:,0]-1,0]-xtp) \
                      <1e-13,abs(coord[bcnod[:,0]-1,1]-ytp)<1e-13))[0][0]
            mask[grain_nod] = False 
                
            
        # Update bcnod and bcval
        bcval = bcval[mask]
        bcnod = bcnod[mask]
        
        bcval_idx = bcval.copy()

      
        # --- Topology quantities ---
            
        nelm  = np.shape(edof)[0]
        
        nnod  = np.shape(coord)[0]
        
        ndof  = nnod*dofs_per_node
        
        nodel = np.size(enod[0])
        
        dofel = nodel*dofs_per_node
    

        # --- Transfer model variables to output data ---
        self.output_data.coord          = coord
        self.output_data.enod           = enod
        self.output_data.edof           = edof
        self.output_data.bcnod          = bcnod
        self.output_data.bcval          = bcval
        self.output_data.bcval_idx      = bcval_idx
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
            elm_coords = coord[elm_nodes-1]
    
            # Check if the nodes are in counter-clockwise order, if not, reverse the order
            if self.is_ccw(*elm_coords):
                enod[i] = np.flip(elm_nodes)

        return enod
        

    def is_ccw(self,p1, p2, p3):
        """
        Check if the given three points (nodes) are in counter-clockwise order.
        """
        return (p2[1] - p1[1]) * (p3[0] - p2[0]) > (p2[0] - p1[0]) * (p3[1] - p2[1])

        
        
        
    def reorder_enod(self, coord,enod):
        """"Reordering enod. Start in lower left corner. Counter-clockwise."""
        
        nelm = np.size(enod,0)
        
        # Allocate enod_reshaped
        enod_reshaped = np.zeros(np.shape(enod)).astype(int)
        
        # Loop over all elms and sort
        for i in range(0, nelm):
            el_coord  = coord[enod[i]-1,:]
            el_enod   = enod[i,:]
            
            
            # 1) x and y idx sorted
            sort_idx_x = np.argsort(el_coord[:,0])
            sort_idx_y = np.argsort(el_coord[:,1])
            
            left_points  = el_enod[sort_idx_x[0:2]]
            right_points = el_enod[sort_idx_x[2:]]
            
            if (el_coord[sort_idx_x[0],1]<el_coord[sort_idx_x[1],1]):
                lower_left_point = left_points[0]
                upper_left_point = left_points[1]
            else:
                lower_left_point = left_points[1]
                upper_left_point = left_points[0]
                
            if (el_coord[sort_idx_x[2],1]<el_coord[sort_idx_x[3],1]):
                lower_right_point = right_points[0]
                upper_right_point = right_points[1]
            else:
                lower_right_point = right_points[1]
                upper_right_point = right_points[0]
                

            # Collect reshaped enod
            enod_reshaped[i,:] = np.array([lower_left_point, \
                                           lower_right_point,\
                                           upper_right_point,
                                           upper_left_point],dtype=int)
        
        # Return enod
        return enod_reshaped
    
    
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
        
        # Initialize gmsh:
        gmsh.initialize()
        
        # Set Mesh.SaveAll to 1
        gmsh.option.setNumber("Mesh.SaveAll", 1)
        
        # Generate mesh
        self.generateMesh()
        
        # Finalize the Gmsh API
        gmsh.finalize()
        
        # --- Save input data as json file ---
        # self.input_data.save('phase_input_' + str(self.input_data.version) + '.json')
        
        
        print('cwd: ', os.getcwd())
        # --- Save output data as json file ---
        self.output_data.save('phase_mesh_' + str(self.input_data.version) + '_' + str(self.input_data.grain+1) + '.json')
        
        # # --- Export location to Fortran ---
        # self.exportLocationToFortran()
        
        # --- Leave dir ---
        os.chdir(cwd)
        
        
 
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
