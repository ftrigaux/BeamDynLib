import numpy  as np
import matplotlib.pyplot as plt

from pyBeamDyn import PyBeamDyn

#plt.close('all')

#%% ---- Definition of the dynamic problem  -----
dt = 0.05
nt = 100
q1 = 5e3
q2 = 5e3
to = 0.0

#%% ---- Solve the Dynamic Problem with BeamDyn -- python version -----
nBeam      = 1
L          = 120
inputFile  = "IEA-15-240-RWT_BeamDyn.dat"

pbd = PyBeamDyn()
pbd.initBeamDyn(1,inputFile,1,nt=1,dt=dt,omega=-np.array([0.0,0,0],order='F'),DynamicSolve=1)

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
t = np.linspace(dt,nt*dt,nt);
plt.plot(t,vtip,'--',c='r');
plt.plot(t,wtip,'--',c='g');
plt.xlabel('Time [s]');
plt.grid(True);

plt.figure(2)
plt.plot(t,ptip*180.0/np.pi,'--',c='b');
plt.xlabel('Time [s]');
plt.grid(True);

#%% Read output from CBeamDyn

# data = np.fromfile("u.bin");
# u    = data.reshape((-1,6));
# plt.plot(u[:,0],'.-');
