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

import subprocess


class InputData(object):
    """Class for defining input data to our model."""
    def __init__(self):
        
        # --- Define private variabels, standrad values

        # Version
        self.version            = 1
        
        # Geometry [microns]
        self.enod               = 0
        self.coord              = 0
        self.nelm               = 0
        self.nnod               = 0
        self.nodel              = 0
        self.ls                 = 0
        self.a                  = 0
        self.line_ex            = 0
        self.line_ey            = 0
        self.line_conn          = 0
        self.line_ex_all        = 0
        self.line_ey_all        = 0
        self.nseg               = 0
        self.line_seg_all       = 0
        self.ngrains            = 0
        self.eq_comp            = 0
        self.material           = 0
        self.tppoints           = 0
        self.ls                 = 0
        
        # Equilibrium compositions
        self.cu_cu              = 111
        self.cu_imc             = 0.1  # OBS not acc. to pd
        self.cu_sn              = 0.31117
        self.imc_cu             = 0.38214
        self.imc_imc            = 112
        self.imc_sn             = 0.45292
        self.sn_cu              = 0.999763
        self.sn_imc             = 0.99941
        self.sn_sn              = 113
        self.sn_c0_fix          = 0.999
        self.sn_c0              = 0
    
        
        self.cu_params          = [1.0133e5, -2.1146e4, -1.2842e4, 0.10569]
        self.imc_params         = [4e5     , -6.9892e3, -1.9185e4, 0.41753]
        self.sn_params          = [4.2059e6, 7.1680e3 , -1.5265e4, 0.99941]
        self.molar_volumes      = np.array([7.09, 10.7, 16.3])*1e3

        # Location
        self.location           = None
        self.python_location    = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/Python_output'
        self.fortran_location   = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/build/src/main'
        
    
    def save(self, filename):
        """Save input data to dictionary."""
        input_data_file = {}
        
        # Version
        input_data_file["version"]                = self.version
        
        # Location
        input_data_file["location"]               = self.location  
        
        with open(filename, "w") as ofile:
            json.dump(input_data_file, ofile, sort_keys = True, indent = 4)
            
    

            
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
            self.ngrains        = np.shape(self.a)[1]
            self.line_ex_all    = (data['line_ex']*1e3)
            self.line_ey_all    = data['line_ey']*1e3
            self.line_seg_all   = data['line_seg'].flatten()
            self.material       = data['material'].flatten()
            self.tppoints       = data['tppoints']*1e3
    
            
        # Round line coords to 6 decimals 
        self.line_ex_all = self.line_ex_all.round(decimals=6)
        self.line_ey_all = self.line_ey_all.round(decimals=6)
        
        self.tppoints = self.tppoints.round(decimals=6)
        
                
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
        return unique_coord
    
    def loadmesh(self, python_location):
        """Read output data from file."""

        cwd  = os.getcwd()
        path = python_location + '/single_study/phase_meshes'
        
        # Find file 'phase_mesh_X.mat' with highest index X in mat_files
        files = []
        for i in os.listdir(path):
            if os.path.isfile(os.path.join(path,i)) and 'triangle_mesh_' in i:
                files.append(i)
        max_number = 0
        for file in files:
            file_number = int(file.split('triangle_mesh_')[1].split('.json')[0])
            if (file_number>max_number):
                max_number = file_number
            
        filename = path + '/triangle_mesh_' + str(max_number) + '.json'
        
        
        with open(filename, "r") as ifile:
            output_data_file = json.load(ifile)
            
        self.coord         = np.transpose(np.asarray(output_data_file["coord"]))
        self.enod          = np.transpose(np.asarray(output_data_file["enod"]))-1
        self.nelm          = np.asarray(output_data_file["nelm"])
        self.nnod          = output_data_file["nnod"]
        self.nodel         = np.asarray(output_data_file["nodel"])

    def load_LSInterpolated(self,python_location):
        """Read indata from file."""
        
        cwd  = os.getcwd()
        path = python_location + '/single_study/mat_files'
        
        # Find file 'level_set_X.mat' with highest index X in mat_files
        files = []
        for i in os.listdir(path):
            if os.path.isfile(os.path.join(path,i)) and 'triangle_ls_' in i:
                files.append(i)
        max_number = 0
        for file in files:
            file_number = int(file.split('triangle_ls_')[1].split('.mat')[0])
            if (file_number>max_number):
                max_number = file_number
            
        filename = path + '/triangle_ls_' + str(max_number) + '.mat'
        
        self.version = max_number
        
    
        with open(filename, "r") as ifile:
            
            # Access variables from level set file
            data = scipy.io.loadmat(filename)
            
            # Global variables
            self.ls = data['ls']
        
class OutputData(object):
    """Class for storing output data from calculation."""
    def __init__(self):
    
        
        # Python output
        self.enod            = None
        self.coord           = None
        self.nelm            = None
        self.nnod            = None
        self.nodel           = None
        self.dofspernode     = None
        self.bcnod           = None
        self.bcval           = None
     
        # Locations - not saved
        self.single_location = None
        self.input_location  = None
        
    def reset(self):
        
        # Python output
        self.enod            = None
        self.coord           = None
        self.nelm            = None
        self.nnod            = None
        self.nodel           = None
        self.dofspernode     = None
        self.bcnod           = None
        self.bcval           = None
        
        
    def save(self, filename):
        """Saving output data to dictionary."""
        
        output_data_file = {}
        output_data_file["coord"]           = np.transpose(self.coord).tolist()
        output_data_file["enod"]            = np.transpose(self.enod).tolist()
        output_data_file["bcnod"]           = self.bcnod.tolist()
        output_data_file["bcval"]           = self.bcval.tolist()
        output_data_file["bcval_idx"]       = self.bcval_idx.tolist()
        output_data_file["nelm"]            = self.nelm
        output_data_file["nnod"]            = self.nnod
        output_data_file["nodel"]           = self.nodel
        output_data_file["dofspernode"]     = self.dofspernode
        
        
        with open(filename, "w") as ofile:
            json.dump(output_data_file, ofile, sort_keys = True, indent = 4)


class Mesh(object):
    """Class for generating the mesh"""
    def __init__(self, input_data, output_data):
        self.input_data  = input_data
        self.output_data = output_data
        
    def generateMesh(self, g):
            
     
        # Short form
        ngrains      = self.input_data.ngrains
        line_ex_all  = self.input_data.line_ex_all
        line_ey_all  = self.input_data.line_ey_all
        line_seg_all = self.input_data.line_seg_all
        tppoints     = self.input_data.tppoints
        enod         = self.input_data.enod
        coord        = self.input_data.coord
        nodel        = self.input_data.nodel
        ls           = self.input_data.ls
        
 
        # Extract partitioned mesh for each grain
        ls_ed      = ls[enod,g]
        line_elms  = np.sum((ls_ed<=0),1)==nodel
        enod_g     = np.asarray(enod[line_elms,:],dtype=np.int32)
        nelm_g     = np.shape(enod_g)[0]
        nods_g     = np.unique(enod_g)
        nnod_g     = np.size(nods_g)
        coord_g    = coord[nods_g,:]
        
        
        # New enod
        nods_g_new = np.linspace(0,nnod_g,nnod_g+1,dtype=int)
        for k in range(nnod_g):
            mask = (enod_g == nods_g[k])
            enod_g[mask] = nods_g_new[k]
        
        
        # --- Topology quantities ---
        nodel_g = np.size(enod_g,axis=1)
        
        # Boundary nodes - nodes shared with line_coord
        g_cols      = [2*g, 2*g + 1]
        line_seg    = line_seg_all[g]
        line_ex     = line_ex_all[:line_seg,g_cols]
        line_ey     = line_ey_all[:line_seg,g_cols]
        line_coord  = self.input_data.lines_to_coord(line_ex, line_ey)
        nline_coord = np.size(line_coord,0)
        
        

        bnods_logical = np.full(nnod_g,False)
        for k in range(np.shape(line_coord)[0]):
            P = line_coord[k,:]
            bnods_logical = np.logical_or(bnods_logical,np.linalg.norm(coord_g-P, axis=1)<1e-8)


        # Boundary nodes
        bcnods = np.where(bnods_logical)[0]
        nbc    = np.size(bcnods)
        bcval  = np.zeros_like(bcnods,dtype=float)
        
        # Equilibrium chemical potentials
        eq_comp,mu_imc_imc,mu_sn_sn = self.getChemicalPotentials()
        
        grain_material = self.input_data.material[g] - 1
        
        # Save triple junction points
        tp_points_x = []
        tp_points_y = []
        
        # Determine what other grain has smallest value at that node
        grain_idx    = np.linspace(0,ngrains-1,ngrains,dtype=int)
        other_grains = grain_idx[grain_idx != g]
        
        # Add bcval
        for i in range(nbc):
            
            # bcnod
            bcnod = bcnods[i]
            
            # Point P
            P = coord_g[bcnod,:]
                    
            # See if point P is in other grain
            P_in_g = []
            for gg in other_grains:
                gg_cols = [2*gg, 2*gg + 1]
                line_ex_g = line_ex_all[:line_seg_all[gg],[gg_cols[0],gg_cols[1]]]
                line_ey_g = line_ey_all[:line_seg_all[gg],[gg_cols[0],gg_cols[1]]]
                if ((abs((line_ex_g - P[0]) + (line_ey_g - P[1]))<1e-13).any()):
                    P_in_g.append(gg)
            
            P_in_g = np.array(P_in_g)
            
            if (np.size(P_in_g)>1):
                other_materials = self.input_data.material[P_in_g] - 1
                if (any((other_materials-grain_material) == 0)):
                    P_in_g = P_in_g[other_materials != grain_material]
                else:
                    P_in_g = P_in_g[0]
                        
            if (np.size(P_in_g)==0):
                print('Error in bcval')
                
            other_material = self.input_data.material[P_in_g] - 1
            
            
            if (np.size(other_material)>1):
                print('Error in bcval')

            # Find equilibrium composition grain, lowest_grain
            bcval[i] = eq_comp[grain_material,other_material]
      

            
        # Rearrange bcnod
        bcnod      = np.zeros((nbc,2), dtype=np.int32)
        bcnod[:,0] = bcnods + 1
        bcnod[:,1] = np.ones((nbc), dtype=np.int32)
  
        # Create a boolean mask and remove imc_imc values from bc
        mask  = np.abs(bcval - mu_imc_imc) > 1e-12
        bcval = bcval[mask]
        bcnod = bcnod[mask]
        
        
        # Create a boolean mask and remove sn_sn values from bc
        # mask  = np.abs(bcval - self.input_data.imc_imc)>1e-12
        mask  = np.abs(bcval - mu_sn_sn)>1e-12
        bcval = bcval[mask]
        bcnod = bcnod[mask]        
        
        # Remove tp points from bc
        mask = np.full(np.size(bcnod[:,0]), True, dtype=bool)
        for xtp, ytp in tppoints:        
            tp_loc = np.logical_and(np.abs(coord_g[bcnod[:,0]-1,0]-xtp)<1e-7, \
                                      np.abs(coord_g[bcnod[:,0]-1,1]-ytp)<1e-7)
            mask = np.logical_and(mask, tp_loc==False)    

        # Update bcnod and bcval
        bcval = bcval[mask]
        bcnod = bcnod[mask]
        
        bcval_idx = bcval.copy()
        bcnod_g   = bcnod
        bcval_g   = bcval

        # --- Transfer model variables to output data ---
        self.output_data.coord          = coord_g
        self.output_data.enod           = enod_g + 1
        self.output_data.bcnod          = bcnod_g
        self.output_data.bcval          = bcval_g
        self.output_data.bcval_idx      = bcval_idx
        self.output_data.nelm           = nelm_g
        self.output_data.nnod           = nnod_g
        self.output_data.nodel          = nodel_g
        self.output_data.dofspernode    = 1
        
    

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

        # Extract mesh
        for g in range(self.input_data.ngrains):
            # Reset output data
            self.output_data.reset()
            
            # Generate mesh
            self.generateMesh(g)
      
            # --- Save output data as json file ---
            self.output_data.save('phase_mesh_' + str(self.input_data.version) + '_g' +  str(g+1) + '.json')
        
        
 
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

    def getChemicalPotentials(self):
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
            
        return eq_comp, mu_imc_imc, mu_sn_sn
        