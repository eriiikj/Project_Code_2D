#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 08:34:18 2023

Script for converting .mat file to json format

@author: Erik Jacobsson
"""

import scipy.io
import json
import sys

class MatlabWorkspace(object):
    """Class for keeping info of matlab workspace."""
    
    def __init__(self, filename):
        
        # Filename
        self.filename = filename
        
        # Matlab workspace
        self.mat = scipy.io.loadmat(filename)
        
        
    def save_to_json(self, filename):
        """Save variable to json file."""
        matlab_save = {}
        
        # Variables to save
        matlab_save["a_ls"] = self.mat['a_ls'][:,0].tolist()
        
        with open(filename, "w") as ofile:
            json.dump(matlab_save, ofile, sort_keys = True, indent = 4)
        

if __name__ == "__main__":
    
    
    # Load .mat file provided from command line
    try:
        mat_filename = sys.argv[1]
    except IndexError as error:
        print("An exception occurred:", error)
        print('Provide the name of the .mat file on the command line')
        sys.exit(0)
    
    # Create matWorkspace object and load all data in .mat file
    lsinit = MatlabWorkspace(mat_filename)
    
    # Save in corresponding json file
    json_file_name = 'level_set_init.json'
    lsinit.save_to_json(json_file_name)