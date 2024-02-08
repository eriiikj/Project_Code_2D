# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 09:44:29 2022

@author: Erik Jacobsson
"""

import diffusion_mesh as ph
import json
import numpy as np
import math

import visvis as vv
from pathlib import Path
import shutil

import sys, getopt
import os
import subprocess

import logging
import threading
import time

import matplotlib.pyplot as plt
import itertools

#from PyQt5.QtCore import pyqtSlot, pyqtSignal, QThread
from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QMainWindow, QFileDialog, QMessageBox
from os.path import isfile, join

        
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
        self.input_data.grain                  = 0
        self.input_data.a                      = 0
        self.input_data.ex                     = 0
        self.input_data.ey                     = 0
        self.input_data.coord                  = 0
        self.input_data.axisbc                 = 0
        self.input_data.line_ex                = 0
        self.input_data.line_ey                = 0
        self.input_data.line_ex_all            = 0
        self.input_data.line_ey_all            = 0
        self.input_data.nseg                   = 0
        self.input_data.line_seg_all           = 0
        self.input_data.sep_lines              = 0
        self.input_data.nsep_lines             = 0
        self.input_data.ccsep_log              = 0
        self.input_data.ngrains                = 0
        self.input_data.P1                     = 0
        self.input_data.P2                     = 0
        self.input_data.P3                     = 0
        self.input_data.P4                     = 0
        self.input_data.geoms                  = 0
        self.input_data.line_ex_add            = 0
        self.input_data.line_ey_add            = 0
        self.input_data.line_seg_add           = 0
        self.interface_marker                  = 0
        self.eq_comp                           = 0
        self.material                          = 0
        self.input_data.tot_imc_area_init      = 0
        self.input_data.tot_imc_area           = 0
        self.input_data.P1nod                  = 0
        self.input_data.P2nod                  = 0
        self.input_data.P3nod                  = 0
        self.input_data.P4nod                  = 0
        
        # Equilibrium compositions
        self.input_data.cu_cu                  = 111
        self.input_data.cu_imc                 = 0.1  # OBS not acc. to pd
        self.input_data.cu_sn                  = 0.31117
        self.input_data.imc_cu                 = 0.38214
        self.input_data.imc_imc                = 112
        self.input_data.imc_sn                 = 0.45292
        self.input_data.sn_cu                  = 0.999763
        self.input_data.sn_imc                 = 0.99941
        self.input_data.sn_sn                  = 113
        self.input_data.sn_c0_fix              = 0.999
        self.input_data.sn_c0                  = 0
    
        
        self.input_data.cu_params  = [1.0133e5, -2.1146e4, -1.2842e4, 0.10569]
        self.input_data.imc_params = [4e5     , -6.9892e3, -1.9185e4, 0.41753]
        self.input_data.sn_params  = [4.2059e6, 7.1680e3 , -1.5265e4, 0.99941]
        self.input_data.molar_volumes = np.array([7.09, 10.7, 16.3])*1e3

        # Mesh
        self.input_data.el_size_x              = 15e-3
        self.input_data.el_size_y              = 15e-3
        self.input_data.el_size_factor         = 0.05

        # Location
        self.input_data.location               = None
        self.input_data.python_location        = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/Python_output'
        self.input_data.fortran_location       = '/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/build/src/main'
        
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
        # vis.showGeometry()
        
        # --- Show mesh ---
        vis.showMesh()
        
        # --- Show node and elm labels ---
        # vis.showMatplotlibMesh()
        
        if save:
            vis.saveMesh()
            
        # --- Start event loop ---
        sys.exit(self.app.exec_())
        
    def onLoad(self, path):
        """Loading output data"""
        
        cwd = os.getcwd()
        os.chdir(path)
        
        # --- Load python input data ---
        self.input_data.load('python_input.json')
        
        # --- Load python mesh data ---
        self.output_data.loadmesh('python_mesh.json')
        
        
    def onLoadLsInit(self):
        """Loading geometry data from .mat file"""
        

        # --- Load python input data ---
        self.input_data.load_LsInit()
        
    def onLoadLsStep(self):
        """Loading geometry data from .mat file"""
        
        # --- Load python input data ---
        self.input_data.load_LsStep()
    
        
    def onLoadOutputData(self, filename):
        """Loading output data"""
        
        self.output_data.load(filename)
        
    
if __name__ == "__main__":
    

    
    # --- Create application instance ---
    app = QApplication(sys.argv)
    
    # --- Create main instance ---
    main_instance = MainClass(app)
    

    # --- Load initial level set data
    main_instance.onLoadLsInit()

    # --- Load level set step
    main_instance.onLoadLsStep()
    
    # --- Single study ---
    
    # Create single mesh
    main_instance.onGenerateSingleMesh()
    
    # Visalize mesh
    # main_instance.onShowMesh(save=False)






    
     