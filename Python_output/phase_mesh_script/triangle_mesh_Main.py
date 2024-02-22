# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 09:44:29 2022

@author: Erik Jacobsson
"""

import triangle_mesh as ph
import json
import numpy as np
import math
from pathlib import Path
import shutil
import sys, getopt
import os
import subprocess
import time
import itertools
from os.path import isfile, join

        
class MainClass(object):
    """Main class"""

    def __init__(self):
        """Constructor"""
        
        # --- Create input data and output data instances ---
        self.input_data                        = ph.InputData()
        self.output_data                       = ph.OutputData()
         
    def onGenerateSingleMesh(self):
        """Create a single mesh."""
        
        # --- Create solver instance ---
        self.mesh = ph.Mesh(self.input_data, self.output_data)
        
        # --- Create mesh ---
        self.mesh.generateSingleMesh()


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

    
    # --- Create main instance ---
    main_instance = MainClass()
    
    # --- Load initial level set data
    main_instance.onLoadLsInit()


    # --- Load level set state ---
    main_instance.onLoadLsStep()
    

    # --- Create mesh ---
    main_instance.onGenerateSingleMesh()
