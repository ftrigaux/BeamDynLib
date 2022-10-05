import numpy  as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from time import time 

from pyBeamDyn import PyBeamDyn

plt.close('all')

#%% ---- Definition of the dynamic problem  -----
dt = 3e-2
nt = 200
q1 = 0.0
q2 = 0.0
to = 0.0

#%% Import the loads from file 

FAero = np.loadtxt("loads/loads_9ms_FAST.txt",skiprows=1);
Fx = interp1d(FAero[:,0],FAero[:,1],fill_value="extrapolate");
Fy = interp1d(FAero[:,0],FAero[:,2],fill_value="extrapolate");

#%% ---- Solve the Dynamic Problem with BeamDyn -- python version -----
nBeam      = 1
L          = 63
Rhub       = 1.5
omega      = 9*9/L;
inputFile  = "bd_primary_nrel_5mw_dynamic.inp"
DynamicSolve=1
nt_loc     = 5

if DynamicSolve==0:
    nt = 2

pbd = PyBeamDyn()
pbd.initBeamDyn(1,inputFile,1,nt=nt_loc,dt=dt,omega=-np.array([omega*1,0,0],order='F'),DynamicSolve=DynamicSolve)

# Get the position of the distributed loads
xLoads, xDisp = pbd.getPositions()

# Set the loads
loads = np.zeros((6,pbd.nxLoads));
loads[0,:] =  Fx((xLoads[2,:]+Rhub)/L);
loads[1,:] = -Fy((xLoads[2,:]+Rhub)/L);
loads[5,:] = 0.0;
pbd.setLoads(loads,1)

# Preallocate variables to extract displacement
vtip      = np.zeros(nt)
dvtip     = np.zeros(nt)
wtip      = np.zeros(nt)
ptip      = np.zeros(nt)
    
# Solve the dynamic
print("Solving...");
start = time()
for i in range(nt):
    print(i)
    pbd.solve(1)
    
    # Preallocate variables to extract displacement
    x,u,du = pbd.getDisplacement(1)
    vtip[i]  =  u[0,-1];
    wtip[i]  =  u[1,-1];
    dvtip[i] = du[0,-1];
    ptip[i]  =  u[5,-1];
end = time()
print("Done in %1.4e seconds -> %1.4e [s/iteration]"%(end-start,(end-start)/nt))

# Free the data
pbd.freeBeamDyn(1,1)
#%%

if DynamicSolve:
    plt.figure()
    t = np.linspace(dt,nt*dt*nt_loc,nt) * omega * 180/np.pi / 360;
    plt.plot(t,vtip,'--',c='r',label="Flapwise");
    plt.plot(t,wtip,'--',c='g',label="Edgewise");
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.legend()
    plt.ylabel("Tip displacement [m]")
    
    plt.figure()
    plt.plot(t,ptip*180.0/np.pi,'-',c='r');
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.ylabel('Tip torsion angle [deg]');
    
    print("flapwise tip displacement  = %f"%(vtip[-1]))
    print("edgewise tip displacement  = %f"%(wtip[-1]))
    print("torsion  tip displacement  = %f"%(ptip[-1]))
else:
    plt.figure()
    plt.plot(xDisp[2,:],u[0,:],'-',c='b',label='Flap');
    plt.plot(xDisp[2,:],u[1,:],'-',c='r',label='Edge');
    plt.xlabel('x');
    plt.grid(True);
    plt.ylabel('Displacement [m]');
    plt.legend()
    print("flapwise tip displacement  = %f"%(vtip[-1]))
    print("edgewise tip displacement  = %f"%(wtip[-1]))

#plt.show()

#%% Read output from CBeamDyn

# data = np.fromfile("u.bin");
# u    = data.reshape((-1,6));
# plt.plot(u[:,0],'.-');
