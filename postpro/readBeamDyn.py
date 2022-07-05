#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 10:10:54 2020

@author: ftrigaux

Read the BeamDyn driver output file
Requires numpy and pyplot
"""


import numpy as np
import matplotlib.pyplot as plt
#plt.close('all')

def readBD(filename):
    with open(filename) as fid:
        buffer = ""
        BD = {}
        
        # Read until header line and create dictionary keys
        while not(buffer.startswith("Time")):
            buffer = fid.readline();
            
        keys = buffer.split("\t")
        for i in range(len(keys)):
            keys[i] = keys[i].strip()
            BD[keys[i]] = [];
        
        # Read units line
        buffer = fid.readline();
        
        buffer = fid.readline()
        while not(buffer==""):
            # Read values from line and convert to float
            vals = buffer.split("\t")
            for i in range(len(vals)):
                vals[i] = vals[i].strip()
                vals[i] = float(vals[i])
            
            # Add the data to the dict
            for i in range(len(vals)):
                BD[keys[i]].append(vals[i]);
                
            buffer = fid.readline()
            
        for i in range(len(keys)):
            BD[keys[i]] = np.array(BD[keys[i]]);
        
        return BD
            

# Usage example:
#BD = readBD("bd_driver_dynamic_nrel_5mw.out")
#plt.plot(BD['Time'],BD['TipTDxr'],'k--');
#plt.plot(BD['Time'],BD['TipTDyr'],'k--');