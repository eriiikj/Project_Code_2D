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


cfu.enableLogging()

class InputData(object):
    """Class for defining input data to our model."""
    def __init__(self):
        
        # --- Define private variabels, standrad values

        # Version
        self.version                = 1
        
        # Geometry
        self.geom_w                      = 1
        self.geom_h                      = 2.05
        self.geom_y2                     = 0.5   # lower ycoord mid region
        self.geom_y3                     = 0.6   # cu boundary
        self.geom_y4                     = 1.5   # upper ycoord mid region
        self.geom_y5                     = 1.7   # upper ycoord finer region
        
        
        # Mesh
        self.mesh_el_size_x              = 0.3
        self.mesh_el_size_y_coarse       = 0.5
        self.mesh_el_size_y_fine         = 0.3
        self.mesh_el_size_y_finest       = 0.1
        self.mesh_el_distrib_val         = 5
        
        # MFC
        self.mfc                    = True
        
        # IMC volume increase
        self.imc_vol_transf         = 1.1
        
        # Load loop
        self.imc_steps              = 4

        # Location
        self.location               = None
        
        # Param study
        self.param_steps            = 1
        self.param_mesh_el_size_y_fine   = True 
        self.mesh_el_size_y_fine_Start   = 0.02
        self.mesh_el_size_y_fine_End     = 0.01
        
        
        # Markers
        self.lowerSupportMarker    = 10
        self.sideSupportMarker     = 20
        self.upperSupportMarker    = 30
  
    
     # --- Define the geometry ---
    
    def geometry(self):
        """Create geometry instance"""

        g = cfg.Geometry()
        

        # --- Use local references for shorter form
        
        w  = self.geom_w
        h  = self.geom_h
        
        lowerSupportMarker = self.lowerSupportMarker
        sideSupportMarker  = self.sideSupportMarker
        upperSupportMarker = self.upperSupportMarker
        
        
        # Define coordinates 
        x1 = 0; x2 = w
        y1 = 0; y2 = self.geom_y2; y3 = self.geom_y3; y4 = self.geom_y4;
        y5 = self.geom_y5;  y6 = h
        mh = y4-y2

        # --- Points in model created with (...) method
        g.point([x1,y1], 1)
        g.point([x2,y1], 2)
        g.point([x2,y2], 3)
        g.point([x2,y3], 4)
        g.point([x2,y4], 5)
        g.point([x2,y5], 6)
        g.point([x2,y6], 7)
        g.point([x1,y6], 8)
        g.point([x1,y5], 9)
        g.point([x1,y4], 10)
        g.point([x1,y3], 11)
        g.point([x1,y2], 12)

        # --- Create curves with the spline(...) method ---
        mesh_el_size_x          = self.mesh_el_size_x
        mesh_el_size_y_coarse   = self.mesh_el_size_y_coarse
        mesh_el_size_y_fine     = self.mesh_el_size_y_fine
        mesh_el_size_y_finest   = self.mesh_el_size_y_finest

        nelmx              = np.round(w/mesh_el_size_x)
        nelmy_coarse_lower = np.round((y2-y1)/mesh_el_size_y_coarse) # coarse nelm on lower side
        nelmy_finest_lower = np.round((y3-y2)/mesh_el_size_y_finest)
        nelmy_finest_upper = np.round((y4-y3)/mesh_el_size_y_finest)
        nelmy_fine_upper   = np.round((y5-y4)/mesh_el_size_y_fine)
        nelmy_coarse_upper = np.round((y6-y5)/mesh_el_size_y_coarse) # coarse nelm on upper side
        
        el_distrb = self.mesh_el_distrib_val
        g.spline([1 , 2] , 1 , el_on_curve=nelmx, marker=lowerSupportMarker, \
                 el_distrib_type="bump", el_distrib_val=el_distrb)
        g.spline([2 , 3] , 2 , el_on_curve=nelmy_coarse_lower,marker=sideSupportMarker)
        g.spline([3 , 4] , 3 , el_on_curve=nelmy_finest_lower,marker=sideSupportMarker)
        g.spline([4 , 5] , 4 , el_on_curve=nelmy_finest_upper,marker=sideSupportMarker)
        g.spline([5 , 6] , 5 , el_on_curve=nelmy_fine_upper,marker=sideSupportMarker)
        g.spline([6 , 7] , 6 , el_on_curve=nelmy_coarse_upper,marker=sideSupportMarker)
        g.spline([7 , 8] , 7 , el_on_curve=nelmx,el_distrib_type="bump", \
                 el_distrib_val=el_distrb,marker=upperSupportMarker)
        g.spline([8 , 9] , 8 , el_on_curve=nelmy_coarse_upper,marker=sideSupportMarker)
        g.spline([9 , 10], 9 , el_on_curve=nelmy_fine_upper,marker=sideSupportMarker)
        g.spline([10, 11], 10, el_on_curve=nelmy_finest_upper,marker=sideSupportMarker)
        g.spline([11, 12], 11, el_on_curve=nelmy_finest_lower,marker=sideSupportMarker)
        g.spline([12, 1] , 12, el_on_curve=nelmy_coarse_lower,marker=sideSupportMarker)
        g.spline([12, 3] , 13, el_on_curve=nelmx,el_distrib_type="bump", \
                 el_distrib_val=el_distrb)
        g.spline([11, 4] , 14, el_on_curve=nelmx,el_distrib_type="bump", \
                 el_distrib_val=el_distrb)
        g.spline([10, 5] , 15, el_on_curve=nelmx,el_distrib_type="bump", \
                 el_distrib_val=el_distrb)
        g.spline([9 , 6] , 16, el_on_curve=nelmx,el_distrib_type="bump", \
                 el_distrib_val=el_distrb)

        
        # --- Define structured surfaces
        g.addStructuredSurface([1,2,13,12] , 1)
        g.addStructuredSurface([13,3,14,11], 2)
        g.addStructuredSurface([14,4,15,10], 3)
        g.addStructuredSurface([15,5,16,9] , 4)
        g.addStructuredSurface([16,6,7,8]  , 5)


        # --- Finally return the created geometry     
        return g

    
    def save(self, filename):
        """Save input data to dictionary."""
        input_data_file = {}
        
        # Version
        input_data_file["version"]                = self.version
        
        # Geometry
        input_data_file["geom_w"]                      = self.geom_w
        input_data_file["geom_h"]                      = self.geom_h
        input_data_file["geom_y2"]                     = self.geom_y2
        input_data_file["geom_y3"]                     = self.geom_y3
        input_data_file["geom_y4"]                     = self.geom_y4
        input_data_file["geom_y5"]                     = self.geom_y5
        
        # Mesh
        input_data_file["mesh_el_size_x"]              = self.mesh_el_size_x
        input_data_file["mesh_el_size_y_coarse"]       = self.mesh_el_size_y_coarse
        input_data_file["mesh_el_size_y_fine"]         = self.mesh_el_size_y_fine
        input_data_file["mesh_el_size_y_finest"]       = self.mesh_el_size_y_finest
        
        # MFC
        input_data_file["mfc"]                    = self.mfc
        
        # IMC volume increase
        input_data_file["imc_vol_transf"]         = self.imc_vol_transf
         
        # Load loop
        input_data_file["imc_steps"]              = self.imc_steps
        
        # Level set parameters
        input_data_file["ls_gamma"]               = self.ls_gamma
        input_data_file["ls_m"]                   = self.ls_m
        input_data_file["ls_zeta"]                = self.ls_zeta
        input_data_file["ls_mzetagb"]             = self.ls_mzetagb
        input_data_file["ls_phiw"]                = self.ls_phiw
        input_data_file["ls_malpha"]              = self.ls_malpha
        input_data_file["ls_DIMC"]                = self.ls_DIMC

        
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
        
        # Level set region
        self.lvlsety0         = input_data_file["lvlsety0"]
        self.lvlsety1         = input_data_file["lvlsety1"]
        
        # Mesh
        self.mesh_el_size_x        = input_data_file["mesh_el_size_x"]
        self.mesh_el_size_y_coarse = input_data_file["mesh_el_size_y_coarse"]
        self.mesh_el_size_y_fine   = input_data_file["mesh_el_size_y_fine"]
        self.mesh_el_size_y_finest = input_data_file["mesh_el_size_y_finest"]
        
        # MFC
        self.mfc              = input_data_file["mfc"]
        
        # IMC boundary
        self.imc_side_height_start  = input_data_file["imc_side_height_start"]
        self.imc_mid_height_start   = input_data_file["imc_mid_height_start"]
        self.imc_radius_inc         = input_data_file["imc_radius_inc"]
        self.interface_w            = input_data_file["interface_w"]
        self.imc_final_vol_per_area = input_data_file["imc_final_vol_per_area"]
        
        # IMC volume increase
        self.imc_vol_transf         = input_data_file["imc_vol_transf"]
        
        # Load loop
        self.imc_steps        = input_data_file["imc_steps"]
        self.load_steps       = input_data_file["load_steps"]
        self.load_steps_init  = input_data_file["load_steps_init"]
        
        # Location
        self.location         = input_data_file["location"]
        
class OutputData(object):
    """Class for storing output data from calculation."""
    def __init__(self):
        
        # Python output
        self.coord           = None
        self.enod            = None
        self.edof            = None
        self.bcnod           = None
        self.bcval           = None
        self.dofs            = None
        self.dofs_per_node   = None
        self.el_type         = None
        self.nelm            = None
        self.nnod            = None
        self.ndof            = None
        self.nodel           = None
        self.dofel           = None
        
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
        output_data_file["bcnod"]           = self.bcnod.tolist()
        output_data_file["bcval"]           = self.bcval.tolist()
        output_data_file["dofs"]            = self.dofs.tolist()
        output_data_file["dofs_per_node"]   = self.dofs_per_node
        output_data_file["el_type"]         = self.el_type
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

        
    # def loadFortranmain(self,filename):
    #     """Read fortran output mainout.json file."""

    #     with open(filename, "r") as ifile:
    #         output_data_file = json.load(ifile)
            
    #     self.ex             = np.transpose(np.asarray(output_data_file["ex"]))
    #     self.ey             = np.transpose(np.asarray(output_data_file["ey"]))
        
        
    # def loadFortranloadstep(self, filename):
    #     """Read fortran output from loadstep file."""

    #     with open(filename, "r") as ifile:
    #         output_data_file = json.load(ifile)
     
    #     self.vm              = np.asarray(output_data_file["vm"])
    #     self.biax_Sn         = np.asarray(output_data_file["biax_Sn"])
    #     self.i_imc           = output_data_file["i_IMC"]
        
    def loadMatLoadstep(self, filename):
        """Read fortran output from loadstep file."""

        matfile = scipy.io.loadmat(filename)
     
        self.vm      = np.asarray(matfile["vm"])
        self.biax_Sn = np.asarray(matfile["biax_Sn"])
        
        
class Solver(object):
    """Class for generating the mesh and executing Fortran code"""
    def __init__(self, input_data, output_data):
        self.input_data  = input_data
        self.output_data = output_data
        
    def generateMesh(self):
        # --- Call input_data for geometry ---
        geometry = self.input_data.geometry()

        # --- Mesh ---
        el_type        = 3        # <-- Quad elements
        dofs_per_node  = 2        # <-- 2D problem
        
        mesh                          = cfm.GmshMesh(geometry)   
        mesh.el_type                  = el_type
        mesh.dofs_per_node            = dofs_per_node
        mesh.return_boundary_elements = True
        coord, edof, dofs, bdofs, elementmarkers, \
        boundaryElements              = mesh.create()
        
        # Obtain and rearrange enod
        enod   = (edof[:,0:-1:2]+1)/2
        enod   = enod.astype(int)
        enod   = self.reorder_enod(coord,enod)
        
        print('Size enod:', np.shape(enod))
        print('nnods:', np.shape(coord))
        
        # --- Extract bcnod ---
        
        
        # -- Lower bc nods
        bcnods_lower = np.asarray(bdofs[self.input_data.lowerSupportMarker])\
            [0:-1:2]
        bcnods_lower = (bcnods_lower+1)/2
        bcnods_lower = bcnods_lower.astype('int')
        
        
        # -- Upper bc nods
        bcnods_upper = np.asarray(bdofs[self.input_data.upperSupportMarker])\
            [0:-1:2]
        bcnods_upper = (bcnods_upper+1)/2
        bcnods_upper = bcnods_upper.astype('int')
        
        # Remove nods of lower left and right corner from bcnods 
        # (included in MFC)
        if (self.input_data.mfc):
            x0 = min(coord[:,0])
            x1 = max(coord[:,0])
            
            remove_bcnods_bool = np.logical_or((abs(coord[bcnods_lower-1,0]-x0)<1e-6), \
                                (abs(coord[bcnods_lower-1,0]-x1)<1e-6))
            
            bcnods_lower = bcnods_lower[remove_bcnods_bool==False]
        
        
        nbc_upper     = np.size(bcnods_upper)
        bcdofsy_upper = 2*np.ones(nbc_upper)
        
        nbc_lower     = np.size(bcnods_lower)
        bcdofsy_lower = 2*np.ones(nbc_lower)
        bcdofsx_lower = 1*np.ones(nbc_lower)
        
        bcnods_lower = bcnods_lower.astype(int)
        
        
        
        # # -- Side bc nods
        # if (self.input_data.mfc==False):
        #     bcnods_sides = np.asarray(bdofs[self.input_data.sideSupportMarker])\
        #         [0:-1:2]
        #     bcnods_sides = bcnods_sides[2:] # Dont include corners again
        #     bcnods_sides = (bcnods_sides+1)/2
        #     bcnods_sides = bcnods_sides.astype('int')
     
        #     nbc_sides    = np.size(bcnods_sides)
        #     bcdofs_sides = 1*np.ones(nbc_sides)
            
        #     bcnods_sides = bcnods_sides.astype(int)
            
        #     # Concatenate bcnods and bcdofs
        #     bcnods = np.concatenate([bcnods_lower,bcnods_lower,bcnods_sides])
        #     bcdofs = np.concatenate([bcdofsx_lower,bcdofsy_lower,bcdofs_sides])
        # else:
        #     bcnods_sides = []
            
        #     # Concatenate bcnods and bcdofs
        #     bcnods = np.concatenate([bcnods_lower,bcnods_lower])
        #     bcdofs = np.concatenate([bcdofsx_lower,bcdofsy_lower])
            
            
        # -- Side bc nods
        if (self.input_data.mfc==False):
            bcnods_sides = np.asarray(bdofs[self.input_data.sideSupportMarker])\
                [0:-1:2]
            bcnods_sides = bcnods_sides[2:] # Dont include corners again
            bcnods_sides = (bcnods_sides+1)/2
            bcnods_sides = bcnods_sides.astype('int')
     
            nbc_sides    = np.size(bcnods_sides)
            bcdofs_sides = 1*np.ones(nbc_sides)
            
            bcnods_sides = bcnods_sides.astype(int)
            
            # Concatenate bcnods and bcdofs
            bcnods = np.concatenate([bcnods_lower,bcnods_sides])
            bcdofs = np.concatenate([bcdofsy_lower,bcdofs_sides])
        else:
            bcnods_sides = []
            
            # Concatenate bcnods and bcdofs
            bcnods = np.concatenate([bcnods_lower])
            bcdofs = np.concatenate([bcdofsy_lower])
                
    
            # Constraint x in middle node
            bcnod_xcoord = coord[bcnods-1,0]                # x coord of bcnod
            mw           = self.input_data.geom_w/2              # xmiddle
            bcnod_x_idx  = np.abs(bcnod_xcoord-mw).argmin() # idx of closest node
            bcnod_x      = bcnods[bcnod_x_idx]              # node
            bcdof_x      = 1                                # x-dir
            
            # Append x constraint
            bcnods = np.append(bcnods,bcnod_x)
            bcdofs = np.append(bcdofs,bcdof_x)
        
        
        
        
        
        nbc    = np.size(bcnods)
        
        bcnod      = np.zeros((nbc,2), dtype=np.int32)
        bcnod[:,0] = bcnods
        bcnod[:,1] = bcdofs
        
        bcval = np.zeros(np.shape(bcnod)[0])
        
        nelm  = np.shape(edof)[0]
        
        nnod  = np.shape(coord)[0]
        
        ndof  = nnod*dofs_per_node
        
        nodel = np.size(enod[0])
        
        dofel = nodel*dofs_per_node
        
        
        # --- Transfer model variables to output data ---
        self.output_data.coord         = coord
        self.output_data.enod          = enod
        self.output_data.edof          = edof
        self.output_data.bcnod         = bcnod
        self.output_data.bcval         = bcval
        self.output_data.dofs          = dofs
        self.output_data.dofs_per_node = dofs_per_node
        self.output_data.el_type       = el_type
        self.output_data.nelm          = nelm
        self.output_data.nnod          = nnod
        self.output_data.ndof          = ndof
        self.output_data.nodel         = nodel
        self.output_data.dofel         = dofel
        
  
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
        
        # --- Make new directory and enter ---     
        cwd = os.getcwd() # Save path
        single_dir = 'single_study'
        self.makeDirinCurrent(single_dir)
        
        # --- Save location ---
        self.input_data.location = os.getcwd()

        # --- Create mesh ---
        self.generateMesh()
        
        # --- Save input data as json file ---
        self.input_data.save('python_input.json')
        
        # --- Save output data as json file ---
        self.output_data.save('python_mesh.json')
        
        # --- Export location to Fortran ---
        self.exportLocationToFortran()
        
        # --- Leave dir ---
        os.chdir(cwd)
        
        
    def execute(self):
        abs_fortran_path = self.input_data.fortran_location
        isDir            = os.path.isdir(abs_fortran_path) 
        if isDir:
            
            # --- Enter Fortran directory ---
            os.chdir(abs_fortran_path)
            try:
                subprocess.run('./mainout')
            except OSError as e:
                e = str(e)
                print("File '%s' does not exist." % e[e.find('/')+1:-1])
            except KeyboardInterrupt:
                print("------------------------")
                print("Program code interrupted")
        else:
            print("Executable path does not exist")
        
    
    def executeParamStudy(self):
        """Generates several meshes. Even sized elms."""

        # --- Save initial input data ---
        old_param_mesh_el_size_y_finest = self.input_data.mesh_el_size_y_finest
        
        
        # --- Param study ---
        if self.input_data.param_mesh_el_size_y_finest:
            mesh_el_size_y_fine_Range = np.linspace(
                                self.input_data.mesh_el_size_y_finest_Start,
                                self.input_data.mesh_el_size_y_finest_End,
                                self.input_data.param_steps)
        
        
            # --- Make new directory and enter ---
            cwd1 = os.getcwd() # Save path
            self.makeDirinCurrent('param_study' )
            
            i = 1
            for mesh_el_size_y_fine in mesh_el_size_y_fine_Range:
                print("-------------------------------------------")    
                print("Executing for el_size = %g..." % mesh_el_size_y_fine)  
                
                # --- Make new directory and enter ---
                cwd2 = os.getcwd() # Save path
                dir_name = "param_study_mesh_{mesh_size:.3f}".format(mesh_size=mesh_el_size_y_fine)
                dir_name = dir_name.replace('.','')
                self.makeDirinCurrent(dir_name)
                
                # --- Assign parameter to InputData-instance
                self.input_data.mesh_el_size_y_finest = mesh_el_size_y_fine 
                
                # --- Save location ---
                self.input_data.location = os.getcwd()
    
                # --- Create mesh ---
                self.generateMesh()
                
                # --- Save input data as json file ---
                self.input_data.save("python_input.json")
                
                # --- Save mesh data as json file ---
                self.output_data.save("python_mesh.json")
                
                # --- Export location to Fortran ---
                self.exportLocationToFortran()
                
                # --- Execute Fortran program ---
                self.execute()
    
                # --- Leave dir ---
                os.chdir(cwd2)
    
                i += 1 
        
        # --- Reset input data ---
        self.input_data.mesh_el_size_y_finest = old_param_mesh_el_size_y_finest
        
        # --- Leave dir ---
        os.chdir(cwd1)
 
 
    def exportLocationToFortran(self):
        """Export location to Fortran program."""
        
        input_data_location_file                   = {}
        input_data_location_file["input_location"] = self.input_data.location

        exportFileName = 'python_input_location.json'
        
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
        self.add_text("mesh_el_size_x [m]:")
        self.add_text(self.input_data.mesh_el_size_x)
        self.add_text()
        self.add_text("mesh_el_size_y_coarse [m]:")
        self.add_text(self.input_data.mesh_el_size_y_coarse)
        self.add_text()
        self.add_text("mesh_el_size_y_fine [m]:")
        self.add_text(self.input_data.mesh_el_size_y_fine)
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
        
        
    def plot_mean_height_stress(self, filename, output_dir, plot_dir):
        """Routine for plotting mean stress"""
        
        cwd = os.getcwd() # Save path
        
        # Go to output dir
        
        # --- Go to py_files and loop over all load steps ---
        os.chdir(output_dir)
        os.chdir('mat_files')
        print("cwd: ", os.getcwd())
        
        # Total steps
        # total_steps = imc_steps*load_steps
        total_steps = np.size(os.listdir(os.getcwd()))-2
        # print('total_steps:', total_steps)
        
        
        # Short form
        imc_vol_transf   = self.input_data.imc_vol_transf
        mesh_el_size_y_finest = self.input_data.mesh_el_size_y_finest 
        imc_steps        = self.input_data.imc_steps
        load_steps       = self.input_data.load_steps
        
        
        # Short form
        coord = self.output_data.coord
        enod  = self.output_data.enod
        edof  = self.output_data.edof
        ex    = self.output_data.ex
        ey    = self.output_data.ey
        nelm  = self.output_data.nelm
        dofs  = self.output_data.dofs
        
        # --- Unit ---
        unit = 1e-3 # mm

        # --- Mean stress heights ---
        mean_stress_heights = np.asarray([1.2, 1.4, 1.6, 1.8])
       
        # --- Mean stress matrix. Row: loadstep. Col: height. ---
        vm_mean_stress_mat  = np.zeros([total_steps,
                                        np.size(mean_stress_heights)]);
    
        for i in range(total_steps): 
            loadstepfile = "".join(['stress_',str(i+1),'.mat'])   
            self.output_data.loadMatLoadstep(loadstepfile)
            
            # Short form vm
            vm = np.transpose(self.output_data.vm)
            
            for j in range(np.size(mean_stress_heights)):
                # Height
                h = mean_stress_heights[j]
                
                # Determine elms at mean stress height
                a = ey>h
                a = np.sum(a,1)
                c = np.logical_and(a>0,a<4)
                mean_stress_elms = np.asarray([i[0] for i in enumerate(c) if i[1]],dtype=int) # Obs 0 based
                
                # Mean stress at height=sum of stresses in all element gps/n.o gps
                vm_mean_stress_elms     = vm[mean_stress_elms,:]
                vm_mean_stress_mat[i,j] = np.sum(vm_mean_stress_elms)/ \
                                          np.size(vm_mean_stress_elms)
                                       
                                   
        # --- Plot mean stresses --- 
        #plt.style.use('_classic_test_patch')
        fig, ax = plt.subplots()
        marker  = itertools.cycle(('+', 'o', 11)) 
        xvec    = np.linspace(1,total_steps,total_steps)
        
        load_step_limits = load_steps*np.linspace(1,imc_steps,imc_steps,dtype=int)
        load_step_limits = np.insert(load_step_limits,0,1)
        
        
        for i in range(np.size(mean_stress_heights)): 
            
            
            yvec = vm_mean_stress_mat[:,i]
            h    = mean_stress_heights[i]
           
            # Plot mean stress
            plt.plot(xvec, yvec, '--',linewidth=1.0,marker = 'D',
                      markersize=4, label="h={:.1f}".format(h))
            ax.grid(color='k', linestyle='--', linewidth=0.2)
            xticks    = np.arange(0, total_steps+1,load_steps)
            xticks[0] = 1
            ax.set(xlim=(-0.5, total_steps+2), 
                    xticks=xticks, ylim=(0, 20),
                    yticks=np.arange(1, 20))
            
            plt.title(('Mean stress at different heights. Vol transform inc: {:.2f}. '\
                        'Mesh size: {:.3f}.').format(imc_vol_transf,mesh_el_size_y_finest))
                
            plt.xlabel('Load steps')
            plt.ylabel('Mean von Mises stress [MPa]')
            
            
        # --- Plot vertical lines ---
        for i in range(np.size(load_step_limits)-1):
            limit = load_step_limits[i]
            plt.vlines(limit,ymin=0,ymax=16,colors='k', 
                        linestyles='dashed',linewidth=0.4)
            plt.text(limit-1,0.2, "IMC step {:d}".format(i+1),
                      fontsize = 10, rotation="vertical")
 
        # Legend
        ax.legend(title='Heights:', loc='lower right')
            
        # Show
        plt.show()
        
        # Save
        os.chdir(plot_dir)
        print("plotting in: ", os.getcwd())
        plt.savefig(filename, bbox_inches='tight')

            
    
    def plot_mean_biax_stress(self, filename, output_dir, plot_dir):
        """Routine for plotting mean stress"""
        
        cwd = os.getcwd() # Save path

        # Short form
        imc_vol_transf   = self.input_data.imc_vol_transf
        mesh_el_size_y_finest = self.input_data.mesh_el_size_y_finest 
        imc_steps        = self.input_data.imc_steps
        load_steps       = self.input_data.load_steps
        
        # --- Go to py_files and loop over all load steps ---
        os.chdir(output_dir)
        os.chdir('mat_files')
        print("cwd: ", os.getcwd())
        
        # Total steps
        # total_steps = imc_steps*load_steps
        total_steps = np.size(os.listdir(os.getcwd()))-2
        print('total_steps:', total_steps)
        
        # Biax vector
        biax_Sn_mean_stress_steps = np.zeros(total_steps)
        
        # Short form
        coord = self.output_data.coord
        enod  = self.output_data.enod
        edof  = self.output_data.edof
        ex    = self.output_data.ex
        ey    = self.output_data.ey
        nelm  = self.output_data.nelm
        dofs  = self.output_data.dofs
       
        
        # --- Unit ---
        unit = 1e-3 # mm
        
        for i in range(total_steps): 
            loadstepfile = "".join(['stress_',str(i+1),'.mat'])   
            self.output_data.loadMatLoadstep(loadstepfile)
            
            # Short form biax_Sn
            biax_Sn = np.transpose(self.output_data.biax_Sn)
            
            biax_Sn_mean_stress_steps[i] = np.sum(biax_Sn)/np.size(biax_Sn)
                                      
                                   
        # --- Plot mean stresses --- 
        #plt.style.use('_classic_test_patch')
        fig, ax = plt.subplots()
        marker  = itertools.cycle(('+', 'o', 11)) 
        xvec    = np.linspace(1,total_steps,total_steps)
        yvec    = biax_Sn_mean_stress_steps
        
        load_step_limits = load_steps*np.linspace(1,imc_steps,imc_steps,dtype=int)
        load_step_limits = np.insert(load_step_limits,0,1)
        
           
        # Plot biaxial mean stress for all steps
        plt.plot(xvec, yvec, '--',linewidth=1.0,marker = 'D',
                  markersize=4)
        ax.grid(color='k', linestyle='--', linewidth=0.2)
        xticks    = np.arange(0, total_steps+1,load_steps)
        xticks[0] = 1
        ax.set(xlim=(-0.5, total_steps+2), 
                xticks=xticks, ylim=(-15, 7),
                yticks=np.arange(-15, 7,5))
        
        plt.title(('Mean biaxial stress. Vol transform inc: {:.2f}. '\
                    'Mesh size: {:.3f}.').format(imc_vol_transf,mesh_el_size_y_finest))
            
        plt.xlabel('Load steps')
        plt.ylabel('Mean biaxial stress [MPa]')
            
            
        # --- Plot vertical lines ---
        for i in range(np.size(load_step_limits)-1):
            limit = load_step_limits[i]
            plt.vlines(limit,ymin=-30,ymax=30,colors='k', 
                        linestyles='dashed',linewidth=0.4)
            plt.text(limit-1,0.2, "IMC step {:d}".format(i+1),
                      fontsize = 10, rotation="vertical")
 
        # Legend
        #ax.legend(title='Heights:', loc='lower right')
            
        # Show
        plt.show()
        
        # Save
        os.chdir(plot_dir)
        print("plotting in: ", os.getcwd())
        plt.savefig(filename, bbox_inches='tight')
        
    def plot_mean_stress_param(self, param_path, param_dirs, mean_height=False, mean_biax=False):
        """Plot mean stress for several files."""


        # --- Output dir ---
        

        # --- Create new fig dir, remove old
        if mean_height:
            fig_dir = "mean_height_figs"
        elif mean_biax:
            fig_dir = "mean_biax_figs"
            
        plot_dir = os.path.join(os.getcwd(), fig_dir) 
        isDir    = os.path.isdir(plot_dir) 
        if isDir: shutil.rmtree(plot_dir)
        os.mkdir(plot_dir)
        
        for param_dir in param_dirs:
            
            # --- Go into param dir ---
            os.chdir(param_dir)

            # --- Load python input data ---
            self.input_data.load('python_input.json')
            
            # --- Load python mesh data ---
            self.output_data.loadmesh('python_mesh.json')
            
            # --- Load fortran output data ---
            output_dir = 'py_files'
            os.chdir(output_dir)
            self.output_data.loadFortranmain('mainout.json')
            os.chdir(param_dir) # Go back
            
            # --- Plot param mean stress ---
            output_dir = param_path
            mesh_str = param_dir.split("mesh_")[1]
            if mean_height:
                filename = 'mean_height_stress_mesh_{}'.format(mesh_str)
                self.plot_mean_height_stress(filename, output_dir, plot_dir)
            elif mean_biax:
                filename = 'mean_biax__stress_mesh_{}'.format(mesh_str)
                self.plot_mean_biax_stress(filename, output_dir, plot_dir)
                     


    def wait(self):
        """Denna metod ser till att fönstren hålls uppdaterade och kommer att returnera
        När sista fönstret stängs"""

        cfv.showAndWait()

