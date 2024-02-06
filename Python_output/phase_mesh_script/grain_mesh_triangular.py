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

import calfem.core as cfc
import calfem.geometry as cfg  
import calfem.mesh as cfm      
import calfem.vis as cfv
import calfem.utils as cfu

import matplotlib.pyplot as plt
import itertools
import scipy.io

import pyvtk as vtk

import subprocess

from contextlib import redirect_stdout
import io


cfu.enableLogging()

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
        self.cu_height              = 0
        
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
    
    def geometry(self, grain):
        """Create geometry instance"""


        # Define geoms
        self.geoms = np.zeros(self.ngrains)

    
        # --- Use local references for shorter form ---
        el_size_x = self.el_size_x
        el_size_y = self.el_size_y
        
            
        # Short form 
        grain        = self.grain
        line_ex      = self.line_ex
        line_ey      = self.line_ey
        nseg         = self.line_seg
        sep_lines    = self.sep_lines
        nsep_lines   = self.nsep_lines
        ngrains      = self.ngrains
        axisbc       = self.axisbc
        line_ex_add  = self.line_ex_add
        line_ey_add  = self.line_ey_add
        line_seg_add = self.line_seg_add
        nseg_tot     = np.sum(line_seg_add) + nseg
    
        # g_cols
        g_cols = [2*grain, 2*grain + 1]
    
        # Geometry object
        g = cfg.Geometry()
    
    

        # --- Add all line segement points ---
        npoint = 0


        line_counter     = 0
        line_add_counter = 0
        # Loop through all seperate lines
        for sep_idx in range(nsep_lines):
            
            # Cols of sep_idx
            sep_idx_cols = [2*sep_idx, 2*sep_idx + 1]
            
            # Segements to add
            nseg_sep     = sep_lines[sep_idx,1] - sep_lines[sep_idx,0] + 1
            nseg_sep_add = line_seg_add[sep_idx]
  
            # Add points in primary interface
            for i in range(nseg_sep):
                P = [line_ex[line_counter,0],line_ey[line_counter,0]]
                line_counter = line_counter + 1
                
                # Add point
                g.point(P, npoint)
                npoint = npoint + 1
                
            # Add last point if not closed curve
            start_point = np.array([line_ex[0,0], line_ey[0,0]])
            end_point   = np.array([line_ex[line_counter-1,1],line_ey[line_counter-1,1]])
            if (np.linalg.norm(start_point - end_point)>1e-15): # Closed curve
                P = [line_ex[line_counter-1,1],line_ey[line_counter-1,1]]
                g.point(P, npoint)
                npoint = npoint + 1
            
            # Add points in added interface
            for i in range(nseg_sep_add-1):
                P = [line_ex_add[i,sep_idx_cols[1]],line_ey_add[i,sep_idx_cols[1]]]
                
                # Add point
                g.point(P, npoint)
                npoint = npoint + 1
    

        # --- Splines ---
        myid = 0
        
        # Markers
        self.interface_marker = 10
        
        
        
        # Loop through all seperate lines
        npoint_prev = 0
        interface_splines = np.zeros(nseg,dtype=int)
        marker_c = 0
        
        for sep_idx in range(nsep_lines):
            
            # Cols of sep_idx
            sep_idx_cols = [2*sep_idx, 2*sep_idx + 1]
            nseg_sep     = sep_lines[sep_idx,1] - sep_lines[sep_idx,0] + 1
            nseg_sep_add = line_seg_add[sep_idx]
            
            start_point = np.array([line_ex[sep_lines[sep_idx,0]-1,0], line_ey[sep_lines[sep_idx,0]-1,0]])
            end_point   = np.array([line_ex[sep_lines[sep_idx,1]-1,1], line_ey[sep_lines[sep_idx,1]-1,1]])
                
                    
            # Interface splines
            for i in range(nseg_sep):
                left_point  = myid
                right_point = myid + 1
                if (np.linalg.norm(start_point - end_point)<1e-15):
                    if (sep_idx==nsep_lines-1 and nseg_sep_add==0 and i==nseg_sep-1):
                        left_point  = myid
                        right_point = 0
                g.spline([left_point, right_point], ID=myid, el_on_curve=1)
                interface_splines[marker_c] = myid
                marker_c = marker_c + 1
                myid = myid + 1
            
                
            # Additional splines
            for i in range(nseg_sep_add):
                left_point  = myid
                right_point = myid + 1
                if (sep_idx==nsep_lines-1 and i==nseg_sep_add-1):
                    left_point  = myid
                    right_point = 0
                    
    
                g.spline([left_point, right_point], ID=myid)
                myid = myid + 1
     
            npoint_prev = npoint_prev + nseg_sep + nseg_sep_add
            
        s = 9
            
        
        # for sep_idx in range(nsep_lines):
            
        #     # Cols of sep_idx
        #     sep_idx_cols = [2*sep_idx, 2*sep_idx + 1]
            
        #     nseg_sep     = sep_lines[sep_idx,1] - sep_lines[sep_idx,0] + 1
        #     nseg_sep_add = line_seg_add[sep_idx]
            
        #     # Segpoints of line separate segment (closed curve)
        #     loc_counter = 0
        #     segpoints = np.arange(npoint_prev,npoint_prev+nseg_sep+nseg_sep_add)
        #     segpoints = np.append(segpoints, npoint_prev)
        #     s = 9
  
            
        #     # Interface splines
        #     for i in range(nseg_sep):
        #         left_point  = segpoints[loc_counter]
        #         right_point = segpoints[loc_counter+1]
        #         g.spline([left_point, right_point], ID=myid, el_on_curve=1)
        #         interface_splines[marker_c] = myid
        #         marker_c = marker_c + 1
        #         myid = myid + 1
        #         loc_counter = loc_counter + 1
                

                
        #     # Additional splines
        #     for i in range(nseg_sep_add):
        #         left_point  = segpoints[loc_counter]
        #         right_point = segpoints[loc_counter+1]
        #         g.spline([left_point, right_point], ID=myid, el_on_curve=10)
        #         myid = myid + 1
        #         loc_counter = loc_counter + 1
     
        #     npoint_prev = npoint_prev + nseg_sep + nseg_sep_add
            
            
        

        
        
        # Set same marker (interface_marker) to all interface line segements
        for curveID in interface_splines:
            g.curveMarker(curveID, self.interface_marker)
            
        
        # --- Add unstructured surface ---
        splines = list(range(0,nseg_tot))
        g.addSurface(splines)
         
        
        return g
    
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
            
            
            # Short form 
            grain     = self.grain
            g_cols    = [2*grain, 2*grain + 1]
            
            
            # Load the .mat file
            data = scipy.io.loadmat(filename)
        

            left_to_right_prev = False
            
            # Access the variable by its name
            self.a            = data['a']*1e3
            self.line_seg     = data['line_seg'].flatten()
            self.line_seg_all = self.line_seg 
            self.ngrains      = np.size(self.line_seg)
            self.line_seg     = self.line_seg[grain]
            self.line_ex      = data['line_ex']*1e3
            self.line_ex      = self.line_ex[:self.line_seg,g_cols]
            self.line_ey      = data['line_ey']*1e3
            self.line_ey      = self.line_ey[:self.line_seg,g_cols]
            
            self.line_ex_all = data['line_ex']*1e3
            self.line_ey_all = data['line_ey']*1e3
            
            self.sep_lines  = data['sep_lines']
            self.nsep_lines = np.sum(self.sep_lines[:,0:-1:2]>0,0)
            self.nsep_lines = self.nsep_lines[self.grain]
            self.sep_lines  = self.sep_lines[:self.nsep_lines,g_cols]
            self.material   = data['material'].flatten()
            
            # IMC area
            self.tot_imc_area_init = data['IMC_area_init'].flatten()[0]
            self.tot_imc_area      = data['IMC_area'].flatten()[0]
            
            if (self.tot_imc_area_init==0):
                self.sn_c0 = self.sn_c0_fix
            else:
                self.sn_c0 = (self.sn_c0_fix - self.sn_imc)*\
                    (self.tot_imc_area_init/self.tot_imc_area)\
                    +self.sn_imc
                    
            
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
            
            # # Find min and max of lines
            # min_line_ex = np.where(self.line_ex == np.min(self.line_ex))
            # min_line_ex = np.array([min_line_ex[0][0],min_line_ex[1][0]])
            # max_line_ex = np.where(self.line_ex == np.max(self.line_ex))
            # max_line_ex = np.array([max_line_ex[0][0],max_line_ex[1][0]])
            
             
            
            # Identify points:
            #  P4 ----- P3
            #  |         |
            #  |         |
            #  |         |
            #  P1 ----- P2

            # self.P1   = np.array([axisbc[0],axisbc[2]])
            # self.P2   = np.array([axisbc[1],axisbc[2]])
            # self.P3   = np.array([axisbc[1],axisbc[3]])
            # self.P4   = np.array([axisbc[0],axisbc[3]])

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
        self.lower_bcnod     = None
        self.upper_bcnod     = None
        self.dofs            = None
        self.dofs_per_node   = None
        self.el_type         = None
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
        
        # Fortran main output
        self.ex              = None
        self.ey              = None

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
        output_data_file["ex"]              = self.ex.tolist()
        output_data_file["ey"]              = self.ey.tolist()
        output_data_file["bcnod"]           = self.bcnod.tolist()
        output_data_file["bcval"]           = self.bcval.tolist()
        output_data_file["bcval_idx"]       = self.bcval_idx.tolist()
        output_data_file["lower_bcnod"]     = self.lower_bcnod.tolist()
        output_data_file["upper_bcnod"]     = self.upper_bcnod.tolist()
        output_data_file["dofs"]            = self.dofs.tolist()
        output_data_file["dofs_per_node"]   = self.dofs_per_node
        output_data_file["el_type"]         = self.el_type
        output_data_file["nelm"]            = self.nelm
        output_data_file["nnod"]            = self.nnod
        output_data_file["ndof"]            = self.ndof
        output_data_file["nodel"]           = self.nodel
        output_data_file["dofel"]           = self.dofel
        output_data_file["input_location"]  = self.input_location
        output_data_file["bot_c"]           = self.bot_c
        output_data_file["top_c"]           = self.top_c
        
        
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
        
        self.ex, self.ey   = cfc.coordxtr(self.edof, self.coord, self.dofs)

        
        
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
        
        # g_cols
        g_cols = [2*grain, 2*grain + 1]
        
 
    
        # --- Call input_data for geometry ---
        g = self.input_data.geometry(grain)
        
        
        grain_material = self.input_data.material[grain] - 1
        if (grain_material==0 or grain_material==2):
            # Cu or Sn material - flux not important
            self.input_data.el_size_factor = 0.4
            

        # --- Mesh ---
        el_type        = 2        # <-- Triangle elements
        dofs_per_node  = 1        # <-- Concentration problem
        
        mesh                          = cfm.GmshMesh(g)   
        mesh.el_type                  = el_type
        mesh.dofs_per_node            = dofs_per_node
        mesh.return_boundary_elements = True
        mesh.el_size_factor           = self.input_data.el_size_factor # Unstructured mesh
        
        # set a trap and redirect stdout
        trap = io.StringIO()
        with redirect_stdout(trap):
            coord, edof, dofs, bdofs, elementmarkers, \
            boundaryElements  = mesh.create()
        
        
        # Obtain and rearrange enod
        enod = edof
        # enod   = self.reorder_enod(coord,enod)

        
        bcnod        = np.asarray([0])
        bcval        = np.asarray([0])
        lower_bcnods = np.asarray([0])
        upper_bcnods = np.asarray([0])
        bot_c        = np.asarray([0])
        top_c        = np.asarray([0])
        
        
        # Line ex all
        line_ex_all = self.input_data.line_ex_all
        line_ey_all = self.input_data.line_ey_all
        line_seg_all = self.input_data.line_seg_all
    

        # --- Extract bcnod ---
        bcnods = np.asarray(bdofs[self.input_data.interface_marker],\
                                  dtype=int)
        bcval = np.zeros(np.size(bcnods))
        
        
        # Number of bcnods
        n_interface_points = np.size(bcnods) 
        
     
        # Equilibrium concentrations of all materials
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
        
        
        # Diffusion potentials
        c_cu_cu   = cu_cu/molar_volumes[0]
        c_cu_imc  = cu_imc/molar_volumes[0]
        c_cu_sn   = cu_sn/molar_volumes[0]
        
        c_imc_cu  = imc_cu/molar_volumes[1]
        c_imc_imc = imc_imc/molar_volumes[1]
        c_imc_sn  = imc_sn/molar_volumes[1]
        
        c_sn_cu   = sn_cu/molar_volumes[2]
        c_sn_imc  = sn_imc/molar_volumes[2]
        c_sn_sn   = sn_sn/molar_volumes[2]
        
        
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
                
                s = 9
                
             
                # # Find closest node to interface point
                # closest_node = np.argmin(np.sqrt((self.input_data.coord[:,0]\
                #         -P[0])**2 + (self.input_data.coord[:,1]-P[1])**2 ))
                
                # # Find which of the other grains lowest value at node
                # other_grain = other_grains[np.argmin(self.input_data.a\
                #                                 [closest_node,other_grains])]
                
                # other_material = self.input_data.material[other_grain] - 1
                
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
                
                # # Determine what other grain has smallest value at that node
                # other_grains = grain_idx[grain_idx != grain]
                
                # # Find closest node to interface point
                # closest_node = np.argmin(np.sqrt(\
                # (self.input_data.coord[:,0]-P[0])**2 +\
                # (self.input_data.coord[:,1]-P[1])**2 ))
                
                # # Find which of the other grains lowest value at node
                # other_grain = other_grains[np.argmin(self.input_data.a\
                #                             [closest_node,other_grains])]
                    
                # other_material = self.input_data.material[other_grain] - 1
                
                # Find equilibrium composition grain, lowest_grain
                bcval[point_counter] = eq_comp[grain_material,other_material]
                point_counter = point_counter + 1
                
            
            
        
        
        # Rearrange bcnod
        nbc        = np.size(bcnods)
        bcnod      = np.zeros((nbc,2), dtype=np.int32)
        bcnod[:,0] = bcnods
        bcnod[:,1] = np.ones((nbc), dtype=np.int32)
        
    
        
        
        # Remove bc cond that is imc/imc interphase
        
        # Create a boolean mask and remove imc_imc values from bc
        # mask  = np.abs(bcval - self.input_data.imc_imc)>1e-12
        mask  = np.abs(bcval - mu_imc_imc)>1e-12
        bcval = bcval[mask]
        bcnod = bcnod[mask]
        
        
        # Remove bc cond that is sn/sn interphase
        
        # Create a boolean mask and remove imc_imc values from bc
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
        
        
        
        
        # If sn grain add top bc
        # if (grain_material==2):
        #     mu_sn_c0   = sn_params[0]*(sn_c0  - sn_params[3]) + sn_params[1]
        #     maxy      = np.max(coord[:,1])
        #     top_nodes = np.where(abs(coord[:,1] - maxy)<1e-12)
        #     top_nodes = top_nodes[0] + 1
        #     add_bcnod_top      = np.zeros((np.size(top_nodes),2), dtype=np.int32)
        #     add_bcnod_top[:,0] = top_nodes
        #     add_bcnod_top[:,1] = 1
        #     add_bcval_top      = np.zeros(np.size(top_nodes))
        #     add_bcval_top[:]   = mu_sn_c0
            
        #     # Define new bcnod and bcval
        #     bcnod = np.vstack((bcnod, add_bcnod_top))
        #     bcval = np.concatenate((bcval, add_bcval_top))
            
            
        bcval_idx = bcval.copy()
        # # If sn grain change conc at gb
        # if (grain_material==2):
        #     bccoord    = coord[bcnod[:,0]-1,:]
        #     elm_width  = self.input_data.el_size_factor
        #     mu_sn_c0   = sn_params[0]*(sn_c0_fix  - sn_params[3]) + sn_params[1]
        #     minx       = np.min(coord[:,0])
        #     maxx       = np.max(coord[:,0])
        #     miny       = np.min(coord[:,1])
        #     maxy       = np.max(coord[:,1])
            
        #     # Find lower left node in grain mesh
        #     min_y = coord[0][1]
        #     for x, y in coord:
        #         if x < minx + 1e-12:
        #             if y < min_y:
        #                 min_y = y
            
        #     Pll = np.array([minx, min_y])
            
        #     # Find nodes in bcnods
        #     line_bool      = line_ex[:,0]<Pll[0] + elm_width*7
        #     xcoord_include = line_ex[line_bool,0]
        #     ycoord_include = line_ey[line_bool,0]
        #     left_side_nodes = []
        #     for i in range(np.size(xcoord_include)):
        #         xcoord_i = xcoord_include[i]
        #         ycoord_i = ycoord_include[i]
        #         is_included = ((bccoord[:,0] - xcoord_i)**2 + (bccoord[:,1] - ycoord_i)**2)<1e-14
        #         if (np.any(is_included)):
        #             left_side_nodes.append(np.where(is_included)[0][0])
    
        #     left_side_nodes        = np.array(left_side_nodes)
        #     alpha = 5
        #     mu_sn_set = mu_sn_imc + (mu_sn_c0 - mu_sn_imc)*np.exp(-alpha*abs(xcoord_include - Pll[0]))

        #     bcval[left_side_nodes] = mu_sn_set
            
        #     print('as')
            
        
    
        
        # --- Topology quantities ---
            
        nelm  = np.shape(edof)[0]
        
        nnod  = np.shape(coord)[0]
        
        ndof  = nnod*dofs_per_node
        
        nodel = np.size(enod[0])
        
        dofel = nodel*dofs_per_node
        
        # ex and ey in mm
        ex, ey = cfc.coordxtr(enod, coord, dofs)
        ex     = ex*1e-3
        ey     = ey*1e-3
        
        # # Extract boundary nodes along line ex
        # line_ex = self.input_data.line_ex
        # line_ey = self.input_data.line_ey
        # boundary_nodes = self.get_boundary_nodes(coord,enod,ndof,nelm,\
                                                  # line_ex,line_ey)
            
        
        # --- Transfer model variables to output data ---
        self.output_data.coord          = coord
        self.output_data.enod           = enod
        self.output_data.edof           = edof
        self.output_data.ex             = ex
        self.output_data.ey             = ey
        self.output_data.bcnod          = bcnod
        self.output_data.bcval          = bcval
        self.output_data.bcval_idx      = bcval_idx
        self.output_data.lower_bcnod    = lower_bcnods
        self.output_data.upper_bcnod    = upper_bcnods
        self.output_data.dofs           = dofs
        self.output_data.dofs_per_node  = dofs_per_node
        self.output_data.el_type        = el_type
        self.output_data.nelm           = nelm
        self.output_data.nnod           = nnod
        self.output_data.ndof           = ndof
        self.output_data.nodel          = nodel
        self.output_data.dofel          = dofel
        
        
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
        self.generateMesh()
        
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
        
class Report(object):
    """Klass för presentation av indata och utdata i rapportform."""
    def __init__(self, input_data, output_data, save=False):
        self.input_data  = input_data
        self.output_data = output_data
        self.report      = ""

    def clear(self):
        self.report = ""

    def add_text(self, text=""):
        self.report += str(text)+"\n"

    def __str__(self):
        self.clear()
        
         # --- Wrirte report as table ---
        
        self.add_text()
        self.add_text("-------------- Model input ----------------------------------")
        self.add_text()
        self.add_text("Width [m]:")
        self.add_text(self.input_data.w)
        self.add_text()
        self.add_text("Height [m]:")
        self.add_text(self.input_data.h)
        self.add_text()
        self.add_text("y2 [m]:")
        self.add_text(self.input_data.y2)
        self.add_text()
        self.add_text("y3 [m]:")
        self.add_text(self.input_data.y3)
        self.add_text()
        self.add_text("y4 [m]:")
        self.add_text(self.input_data.y4)
        self.add_text()
        self.add_text("el_size_x [m]:")
        self.add_text(self.input_data.el_size_x)
        self.add_text()
        self.add_text("el_size_y_coarse [m]:")
        self.add_text(self.input_data.el_size_y_coarse)
        self.add_text()
        self.add_text("el_size_y_fine [m]:")
        self.add_text(self.input_data.el_size_y_fine)
        self.add_text("-------------- Model output ----------------------------------")
        self.add_text()
        self.add_text("Number of nodes:")
        self.add_text(self.output_data.nnod)
        self.add_text()
        self.add_text("Number of dofs:")
        self.add_text(self.output_data.ndof)
        self.add_text()
        self.add_text("Number of elements:")
        self.add_text(self.output_data.nelm)
        self.add_text()
        self.add_text("Number of nodes in an element:")
        self.add_text(self.output_data.nodel)
        self.add_text()
        self.add_text("Number of dofs in an element:")
        self.add_text(self.output_data.dofel)
        self.add_text()
        self.add_text("Number of dofs in a node:")
        self.add_text(self.output_data.dofs_per_node)
        self.add_text()
        self.add_text("Number of boundary conditions:")
        self.add_text(np.size(self.output_data.bcval))
        self.add_text("------------------------------------------------")
        
        return self.report
    
    def save(self, filename):
        with open(filename, "w", encoding = "UTF-8") as ofile:
            print(self, file = ofile)
    


class Visualisation(object):
    """Class for visualising data"""
    def __init__(self, input_data, output_data):
        self.input_data  = input_data
        self.output_data = output_data



    def showGeometry(self):
        """Show geometry, number of elms at segments etc"""
        # --- Get parameters needed for visualisation ---
        geometry      = self.input_data.geometry()
        
        # --- Draw Geometry ---
        
        cfv.figure()
        cfv.drawGeometry(geometry, draw_points=True, label_curves=True, title='Geometry')
        
        
    def showMatplotlibMesh(self):
        """Plot element and node labels"""
        # --- Get parameters needed for visualisation ---
        coord         = self.output_data.coord
        edof          = self.output_data.edof
        dofs          = self.output_data.dofs
        nnod          = self.output_data.nnod
        
        def drawNodes(nnod,coord):
            """Draw node labels"""
            offset = 0.02
            for i in range(0, nnod):
                cfv.add_text(str(i+1),coord[i]+offset)
                
        # --- Draw mesh with labels ---
        cfv.figure()
        ex, ey = cfc.coordxtr(edof, coord, dofs)
        elms   = np.asarray(range(np.shape(edof)[0])) + 1
        cfv.eldraw2(ex, ey, plotpar=[1, 2, 1], elnum=elms)
        drawNodes(nnod,coord)
        
        
        
    def showMesh(self):
        """Show mesh"""
        
        # --- Get parameters needed for visualisation ---
        coord         = self.output_data.coord
        edof          = self.output_data.edof
        dofs_per_node = self.output_data.dofs_per_node
        el_type       = self.output_data.el_type

        # --- Draw Mesh ---     
        f = vv.figure()
        cfv.drawMesh(
            coords        = coord,
            edof          = edof,
            dofs_per_node = dofs_per_node,
            el_type       = el_type,
            filled        = True,
            title         = "Mesh"
        )
        self.axes = f.currentAxes

    
    def saveMesh(self, filename='python_mesh'):
        """Save mesh as jpg file"""
        
        extension      = '.jpg'
        exportFileName = filename + extension 
    
        vv.callLater(1.0, vv.screenshot, exportFileName, self.axes, sf=2)
        

    def wait(self):
        """Denna metod ser till att fönstren hålls uppdaterade och kommer att returnera
        När sista fönstret stängs"""

        cfv.showAndWait()
