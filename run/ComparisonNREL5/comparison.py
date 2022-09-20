import numpy  as np
import matplotlib.pyplot as plt

from beam3D_NU_v2 import createBeamFromCSV
from pyBeamDyn import PyBeamDyn


plt.close('all')



#%% ---- Definition of the dynamic problem  -----
dt = 0.05
nt = 100
q1 = 5e3
q2 = 5e3
to = 0.0

#%% ---- Solve the Dynamic Problem with PBeams -------
nx  = 24
L   = 61.5
eps = 0.00477465

bm  = createBeamFromCSV("./nrel5mw_data.csv",nx=nx,L=L,eps=eps)

vtip  = np.zeros(nt)
dvtip = np.zeros(nt)
wtip  = np.zeros(nt)
ptip  = np.zeros(nt)
t     = np.linspace(0,dt*nt,nt)

_q1 = q1 * bm.x/bm.L
_q2 = q2 * bm.x/bm.L
for i in range(nt):
    bm.solve_dyn(dt,1,q1=_q1,q2=_q2,t=to)
    vtip[i] = bm.uv['v'][-1];
    dvtip[i] = bm.duv['v'][-1];
    wtip[i] = bm.uv['w'][-1];
    ptip[i] = bm.uv['phi'][-1];
    
    
plt.figure(1)
plt.plot(t,vtip,c='r');
plt.plot(t,wtip,c='g');
plt.plot(t,ptip,c='b');
plt.legend(['Flap','Edge','Twist']);

print("Done solving using PBEAMS");

#%% ---- Solve the Dynamic Problem with BeamDyn -- python version -----
nBeam      = 1
inputFile  = "/Users/ftrigaux/Documents/Beams/BeamDynLib/run/nrel5mw_dynamic/bd_primary_nrel_5mw_dynamic.inp"

pbd = PyBeamDyn()
pbd.initBeamDyn(1,inputFile,1,nt=1,dt=dt,omega=-np.array([0.0,0,0],order='F'),DynamicSolve=1)

# Get the position of the distributed loads
xLoads, xDisp = pbd.getPositions()

# Set the loads
loads = np.zeros((6,pbd.nxLoads),order='F');
loads[0,:] = q1 * (xLoads[2,:]-1)/61.5;
loads[1,:] = q2 * (xLoads[2,:]-1)/61.5;
loads[5,:] = to * (xLoads[2,:]-1)/61.5;
pbd.setLoads(loads,1)

# Preallocate variables to extract displacement
vtip      = np.zeros(nt)
dvtip     = np.zeros(nt)
wtip      = np.zeros(nt)
ptip      = np.zeros(nt)
    
# Solve the dynamic
for i in range(nt):
    pbd.solve(1)
    
    # Preallocate variables to extract displacement
    x,u,du = pbd.getDisplacement(1)
    vtip[i]  =  u[0,-1];
    wtip[i]  =  u[1,-1];
    dvtip[i] = du[0,-1];
    ptip[i]  =  u[5,-1];

# Free the data
pbd.freeBeamDyn(1,1)

plt.figure(1)
plt.plot(t,vtip,'--',c='r');
plt.plot(t,wtip,'--',c='g');
plt.plot(t,ptip,'--',c='b');

#%% Read output from CBeamDyn

# data = np.fromfile("u.bin");
# u    = data.reshape((-1,6));
# plt.plot(u[:,0],'.-');