# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 09:44:29 2022

@author: Erik Jacobsson
"""

import program_handle as ph
import json
import numpy as np
import math

import visvis as vv
from pathlib import Path
import shutil

import sys
import os
import subprocess

import logging
import threading
import time

import matplotlib.pyplot as plt
import itertools

from PyQt5.QtCore import pyqtSlot, pyqtSignal, QThread
from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QMainWindow, QFileDialog, QMessageBox
from PyQt5.QtGui import QFont
from PyQt5.uic import loadUi
from PyQt5.QtGui import *
from os import listdir
from os.path import isfile, join
import calfem.ui as cfui
import calfem.core as cfc
import calfem.vis as cfv
            

        
class MainClass(object):
    """Main class"""

    def __init__(self, app):
        """Constructor"""
        
        self.app = app
        
        # --- Create input data and output data instances ---
        self.input_data                        = ph.InputData()
        
        # Version
        self.input_data.version                = 1
        
        # Geometry [microns]
        el_size = 100e-3
        self.input_data.geom_w                 = 1
        self.input_data.geom_h                 = 2.05   
        self.input_data.geom_y2                = el_size*1    # lower ycoord mid region
        self.input_data.geom_y3                = el_size*4    # cu boundary
        self.input_data.geom_y4                = el_size*6    # upper ycoord mid region
        self.input_data.geom_y5                = el_size*8    # upper ycoord finer region
        
        
        # Mesh
        self.input_data.mesh_el_size_x         = 25e-3 # 25e-3
        self.input_data.mesh_el_size_y_coarse  = 25e-3
        self.input_data.mesh_el_size_y_fine    = 25e-3
        self.input_data.mesh_el_size_y_finest  = 25e-3
        self.input_data.mesh_el_distrib_val    = 1
        
        # MFC (remove lower right and left node from bc)
        self.input_data.mfc                    = False
        
        # IMC volume transformation
        self.input_data.imc_vol_transf         = 1.45
        
        # Load loop
        self.input_data.imc_steps              = 1
        
        # Level set function parameters
        self.input_data.ls_gamma               = 0.625*1e-6
        self.input_data.ls_m                   = 1e-18*(1e3)**3/(1/3600)
        self.input_data.ls_zeta                = 1/self.input_data.ls_m
        self.input_data.ls_mzetagb             = 40
        self.input_data.ls_phiw                = 25e-5
        self.input_data.ls_malpha              = 8e3
        self.input_data.ls_DIMC                = (4.575e-20)*(1e3)**2*3600
        
        # Location
        self.input_data.location               = None
        self.input_data.python_location        = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/Python_output'
        self.input_data.fortran_location       = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/build/src/main'
        
        # Param study
        self.input_data.param_steps            = 3
        self.input_data.param_el_size_y_finest = True
        self.input_data.el_size_y_finest_Start = 0.003
        self.input_data.el_size_y_finest_End   = 0.001
        
        # Output data
        self.output_data                       = ph.OutputData()
        
   
    def onGenerateSingleMesh(self):
        """Create a single mesh."""
        
        # --- Create solver instance ---
        self.solver = ph.Solver(self.input_data, self.output_data)
        
        # --- Create mesh ---
        self.solver.generateSingleMesh()


    def onGenerateReport(self):
        # --- Print report ---
        self.report = ph.Report(self.input_data, self.output_data)
        print(self.report)


    def onShowMesh(self,save=False):
        """Show mesh. Starts event loop"""

        # --- Plot visvis ---
        vis = ph.Visualisation(self.input_data, self.output_data)
        
        
        # --- Show geometry ---
        vis.showGeometry()
        
        # --- Show mesh ---
        # vis.showMesh()
        
        # --- Show node and elm labels ---
        # vis.showMatplotlibMesh()
        
        if save:
            vis.saveMesh()
            
        # --- Start event loop ---
        sys.exit(self.app.exec_())
        
        
    def onExecute(self):    
        """Execute single study"""
        
        # Start timer
        t0 = time.perf_counter()
        
        # --- Execute ---
        self.solver.execute()
        
        # Stop timer
        t1 = time.perf_counter()
        print(t1-t0, "seconds process time")
        
        
    def onExecuteParamStudy(self):
        """Execute Param study"""
        
        # --- Create a solver ---
        self.solver = ph.Solver(self.input_data, self.output_data)
        
        # --- Execute param study --
        self.solver.executeParamStudy()
            
            
    def onPlotMeanHeightStress(self):
        """Routine for plotting mean stress once"""
        
        single_path = self.input_data.python_location + '/single_study'
        os.chdir(single_path)
        
        # --- Load data ---
        self.onLoad(single_path)
        
        # --- Load visvis instance ---
        vis = ph.Visualisation(self.input_data, self.output_data)
        
        # --- Create new fig dir, remove old
        fig_dir  = "mean_height_figs"
        plot_dir = os.path.join(single_path, fig_dir) 
        isDir    = os.path.isdir(plot_dir) 
        if isDir: shutil.rmtree(plot_dir)
        os.mkdir(fig_dir)
        
        # --- Filename of plot ---
        filename = 'mean_height_stress'
        
        
        # --- Plot mean stresses ---
        output_dir = single_path
        vis.plot_mean_height_stress(filename, output_dir, plot_dir)
  
    def onPlotMeanBiaxStress(self):
        """Routine for plotting mean stress once"""
        
        single_path = self.input_data.python_location + '/single_study'
        os.chdir(single_path)
        
        # --- Load data ---
        self.onLoad(single_path)
        
        # --- Load visvis instance ---
        vis = ph.Visualisation(self.input_data, self.output_data)
        
        # --- Create new fig dir, remove old
        fig_dir  = "mean_biax_figs"
        plot_dir = os.path.join(single_path, fig_dir) 
        isDir    = os.path.isdir(plot_dir) 
        if isDir: shutil.rmtree(plot_dir)
        os.mkdir(fig_dir)
           
        # --- Filename of plot ---
        filename = 'mean_biax_stress'

        # --- Plot mean stresses ---
        output_dir = single_path
        vis.plot_mean_biax_stress(filename, output_dir, plot_dir)
        
        
    def onPlotStressParam(self, mean_height=False, mean_biax=False):
        """Routine for plotting param study mean stress for param study"""
        
        param_path = self.input_data.python_location + '/param_study'
        os.chdir(param_path)
        
        # --- Param dirs ---
        param_dirs = self.getParamDirs()
        
        # --- Load visvis instance ---
        vis = ph.Visualisation(self.input_data, self.output_data)

        # --- Plot mean stresses of all param studies---
        vis.plot_mean_stress_param(param_path, param_dirs, mean_height, mean_biax)
            

        
    def onLoad(self, path):
        """Loading output data"""
        
        cwd = os.getcwd()
        os.chdir(path)
        
        # --- Load python input data ---
        self.input_data.load('python_input.json')
        
        # --- Load python mesh data ---
        self.output_data.loadmesh('python_mesh.json')
        
        # # --- Load fortran output data ---
        # output_dir = 'py_files'
        # os.chdir(output_dir)
        # self.output_data.loadFortranmain('mainout.json')
        
        # os.chdir(cwd)
        
    def onLoadOutputData(self, filename):
        """Loading output data"""
        
        self.output_data.load(filename)
        
    
    def getParamDirs(self):
        """Routine for extracting param dirs"""
        
        # --- Go to param dir ---
        param_path  = self.input_data.python_location + '/param_study'
        isDir       = os.path.isdir(param_path) 
        if isDir: 
            # --- Find all param dirs ---
            param_dirs =  sorted(os.listdir(param_path))
        else:
            param_dirs = None
            print("param_dir does not exist.")
            
        # --- Only include dirs whose name include 'param' ---  
        param_dirs_out = []
        for p in param_dirs:
            if "param" in p:
                param_dirs_out.append(os.path.join(param_path, p))

        return param_dirs_out
    
 
    
if __name__ == "__main__":
    
    # --- Create application instance ---
    app = QApplication(sys.argv)
    
    # --- Create main instance ---
    main_instance = MainClass(app)

    
    # --- Single study ---
    
    # Create single mesh
    main_instance.onGenerateSingleMesh()
    
    # Visalize mesh
    # main_instance.onShowMesh(save=False)
    
    
    # Execute
    # main_instance.onExecute()
    
    # Plot mean stress at different heights of single mesh
    # main_instance.onPlotMeanHeightStress()
    
    # Plot mean biaxial stress of single mesh
    # main_instance.onPlotMeanBiaxStress()
    
    
    
    
    
    
    
    # --- Parameter study ---
    #main_instance.onExecuteParamStudy()
    
    # Plot mean stress of param mesh
    #main_instance.onPlotStressParam(mean_height=True)
    #main_instance.onPlotStressParam(mean_biax=True)
    
    
     