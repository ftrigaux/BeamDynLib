#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 17:50:45 2020

@author: ftrigaux
"""
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp

import fsi_tools as fsi

x,bm,nt  = fsi.readStruct("/Users/ftrigaux/Documents/FSI/FSI_oneproc/struct_data.bin",17,61.5,v=1)

dr = x[1]-x[0];
Fy = bm['Fy'][nt-2]*dr;
Fz = bm['Fz'][nt-2]*dr;

with open("forces.dat","w") as fid:
    fid.write(str(Fy[i]));
    fid.write([str(Fz[i]) for i in range(len(Fz))]);
    



