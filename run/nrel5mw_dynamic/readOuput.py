#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 18:08:09 2022

@author: ftrigaux
"""

import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from beam3D_NU_v2 import createBeamFromCSV
from wrtDrvFile import writeDvrFile
from readBeamDyn import readBD

plt.close('all')

nx  = 17
L   = 61.5
eps = 0.00477465

dt = 4e-3
nt = 250
q1 = 5e3;
q2 = 1e4;
to = 0;


# --- Solve with BeamDynLib ---
nt2 = 25
nt_per_step = 10
u = [None]*nt2;
for i in range(nt2):
    u[i] = np.loadtxt("disp%d.dat"%(i+1))
u = np.array(u)
t1 = np.arange(1,nt2+1)*dt*nt_per_step

# --- Solve with PBeams from Fast ----
bm  = createBeamFromCSV("./nrel5mw_data.csv",nx=nx,L=L,eps=eps)
t2   = np.linspace(0,dt*nt,nt)

vtip = np.zeros(nt);
wtip = np.zeros(nt);
for i in range(nt):
    bm.solve_dyn(dt,1,q1=q1,q2=q2,t=to)
    vtip[i] = bm.uv['v'][-1];
    wtip[i] = bm.uv['w'][-1];

# --- Solve with BeamDyn from Fast ----
# write the driver file
dvrFileName="bd_dvr.inp"
outFileName="bd_dvr.out"
writeDvrFile(dvrFileName,dt=dt,t_final=dt*nt,distrLoad=[q1,q2,0.0,0.0,0.0,to],primaryFileName="bd_primary_nrel_5mw_dynamic.inp");

# Launch BeamDyn
sp.call(["/Users/ftrigaux/Documents/Beams/BeamDyn/BeamDyn",dvrFileName])
bd = readBD(outFileName);

#%%
plt.figure(1)
plt.plot(bm.x,bm.uv['u']);
plt.plot(bm.x,bm.uv['v']);
plt.plot(bm.x,bm.uv['w']);

for i in range(3,6):
    plt.plot(u[-1][2,:],u[-1][i,:],'.-')
    plt.grid(True);
    
#%%
plt.figure(2)
plt.plot(t2,vtip,'-')
plt.plot(t1,u[:,0+3,-1],'-')
plt.plot(bd['Time'],bd['TipTDxr'],'-')

plt.plot(t2,wtip,'--')
plt.plot(t1,u[:,1+3,-1],'--')
plt.plot(bd['Time'],bd['TipTDyr'],'--')
