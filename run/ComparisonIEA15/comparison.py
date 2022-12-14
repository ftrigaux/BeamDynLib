import numpy  as np
import matplotlib.pyplot as plt

from beam3D_NU_v2 import Beam
from pyBeamDyn import PyBeamDyn
from scipy.interpolate import interp1d
from scipy.io import loadmat


plt.close('all')


#%% ---- Definition of the dynamic problem  -----
dt = 0.01
nt = 500
q1 = 0.0
q2 = 5e3
to = 0.0

#%% ---- Solve the Dynamic Problem with PBeams -------
data = loadmat('/Users/ftrigaux/Documents/IEA15model/Hawc_data.mat');

nx = 2;
ne = nx - 1;
dx = 1/ne;
xe  = np.arange(dx/2,1+dx/2,dx); # + is to include the last point
L  = 117

rhoA   = interp1d(data['r'][0]/L,data['m'][0]                     ,kind="cubic")(xe)
EIflap = interp1d(data['r'][0]/L,data['E'][0] * data['I_x'][0]    ,kind="cubic")(xe)
EIedg  = interp1d(data['r'][0]/L,data['E'][0] * data['I_y'][0]    ,kind="cubic")(xe)
EA     = interp1d(data['r'][0]/L,data['E'][0] * data['A'][0]      ,kind="cubic")(xe)
GJ     = interp1d(data['r'][0]/L,data['G'][0] * data['I_p'][0]    ,kind="cubic")(xe)
r2     = interp1d(data['r'][0]/L,data['ri_x'][0] + data['ri_y'][0],kind="cubic")(xe)
beta   = interp1d(data['r'][0]/L,data['pitch'][0] * np.pi / 180   ,kind="cubic")(xe)
eps    = 0.00477465

bm = Beam(nx=nx,L=L,rhoA=rhoA,EI1=EIflap,EI2=EIedg,EA=EA,GJ=GJ,r2=r2,eps=eps,twist=beta);

vtip  = np.zeros(nt)
dvtip = np.zeros(nt)
wtip  = np.zeros(nt)
ptip  = np.zeros(nt)
t     = np.linspace(0,dt*nt,nt)

_q1 = q1 * bm.x/bm.L
_q2 = q2 * bm.x/bm.L
for i in range(nt):
    bm.solve_dyn(dt,1,q1=_q1,q2=_q2,t=to)
    #bm.solve_static(q1=_q1,q2=_q2)
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
#inputFile  = "/Users/ftrigaux/Documents/Beams/BeamDynLib/run/IEA15MW_dynamic/IEA-15-240-RWT_BeamDyn.dat"
inputFile  = "/Users/ftrigaux/Documents/Beams/BeamDynLib/run/NREL5MW_dynamic_python/bd_primary_nrel_5mw_dynamic.inp"

pbd = PyBeamDyn()
pbd.initBeamDyn(1,inputFile,1,nt=1,dt=dt,omega=-np.array([1.0,0,0],order='F'),DynamicSolve=1)

# Get the position of the distributed loads
xLoads, xDisp = pbd.getPositions()

# Set the loads
loads = np.zeros((6,pbd.nxLoads),order='F');
loads[0,:] = q1 * (xLoads[2,:]-1)/L;
loads[1,:] = q2 * (xLoads[2,:]-1)/L;
loads[5,:] = to * (xLoads[2,:]-1)/L;
pbd.setLoads(loads,1)

# Preallocate variables to extract displacement
vtip      = np.zeros(nt)
dutip     = np.zeros(nt)
dvtip     = np.zeros(nt)
dwtip     = np.zeros(nt)
wtip      = np.zeros(nt)
ptip      = np.zeros(nt)
t         = np.linspace(dt,dt*nt,nt)
    
# Solve the dynamic
for i in range(nt):
    pbd.solve(1)
    
    # Preallocate variables to extract displacement
    x,u,du = pbd.getDisplacement(1)
    vtip[i]  =  u[0,-1];
    wtip[i]  =  u[1,-1];
    dutip[i] = du[2,-1];
    dvtip[i] = du[0,-1];
    dwtip[i] = du[1,-1];
    ptip[i]  =  u[5,-1];

# Free the data
pbd.freeBeamDyn(1,1)

plt.figure(1)
plt.plot(t,vtip,'--',c='r');
plt.plot(t,wtip,'--',c='g');
plt.plot(t,ptip,'--',c='b');


plt.figure()
plt.plot(dutip)
plt.plot(dvtip)
plt.plot(dwtip)
#%% Read output from CBeamDyn

# data = np.fromfile("u.bin");
# u    = data.reshape((-1,6));
# plt.plot(u[:,0],'.-');