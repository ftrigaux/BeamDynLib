import numpy  as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interp2d

from time import time 

from pyBeamDyn import PyBeamDyn

#plt.close('all')
#%% ---- Definition of the dynamic problem  -----
dt = 2e-2
nt = 200
q1 = 1e4
q2 = 0.0
to = 0.0

#%% ---- Solve the Dynamic Problem with BeamDyn -- python version -----
nBeam      = 1
L          = 86.166
Rhub       = 2.8
omega      = 9*7.5/(L+Rhub);
inputFile  = "DTU10MW_BeamDyn_n048.dat"
DynamicSolve=1
nt_loc     = 1
GlbPos     = np.array([0.0,0.0,Rhub],order='F')
t          = np.linspace(dt,nt*dt*nt_loc,nt)

if DynamicSolve==0:
    nt = 2

pbd = PyBeamDyn()
pbd.initBeamDyn(1,inputFile,1,nt=nt_loc,dt=dt,omega=np.array([omega*1,0,0],order='F'),DynamicSolve=DynamicSolve,GlbPos=GlbPos,gravity=np.array([0,0,-9.81]))

# Get the position of the distributed loads
xLoads, xDisp = pbd.getPositions()

# Set the loads
loads = np.zeros((6,pbd.nxLoads));
loads[0,:] = q1;
loads[1,:] = q2;
loads[5,:] = to;
pbd.setLoads(loads,1)

# Preallocate variables to extract displacement
utip      = np.zeros((3,nt))
dvtip     = np.zeros(nt)
wtip      = np.zeros(nt)
ptip      = np.zeros(nt)
qtip      = np.zeros(nt)
rtip      = np.zeros(nt)
RootMr   = np.zeros((3,nt));

du   = np.zeros((6,xLoads[0].size));
    
# Solve the dynamic
print("Solving...");
start = time()
for i in range(nt):
    # Fake Damping
    # loads[0,:] = q1 - du[0,:]*10;
    # loads[1,:] = q2 - du[1,:]*10;
    # pbd.setLoads(loads,1)
    
    print(i)
    pbd.solve(1)
    
    # Preallocate variables to extract displacement
    x,u,du = pbd.getDisplacement(1)
    utip[:,i]  =  u[:3,-1];
    dvtip[i]   = du[0,-1];
    ptip[i]    =  u[5,-1];
    qtip[i]    =  u[3,-1];
    rtip[i]    =  u[4,-1];

    # Preallocate variables to extract displacement
    x,reactF = pbd.getReactionForce(1)
    RootMr[:,i]  =  reactF[3:,0];

end = time()
print("Done in %1.4e seconds -> %1.4e [s/iteration]"%(end-start,(end-start)/nt))

# Free the data
pbd.freeBeamDyn(1,1)
#%%

if DynamicSolve:
    plt.figure(1)
    t = np.linspace(dt,nt*dt*nt_loc,nt) * omega / 2.0 / np.pi;
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    plt.plot(t,utip.T);
    #plt.plot(t,wtip,'--',c='b',label="Edgewise");
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.legend()
    plt.ylabel("Tip displacement [m]")
    
    plt.figure(2)
    plt.plot(t,ptip*180.0/np.pi,'-',c='r');
    # plt.plot(t,qtip*180.0/np.pi,'-',c='g');
    # plt.plot(t,rtip*180.0/np.pi,'-',c='b');
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.ylabel('Tip angles [deg]');

    plt.figure(5)
    plt.plot(t,RootMr[:,:].T,'-');
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.ylabel('Root Moments [Nm]');

plt.figure(3)
plt.plot(xDisp[2,:],u[0,:],'-',c='b',label='Flap');
plt.plot(xDisp[2,:],u[1,:],'-',c='r',label='Edge');
plt.xlabel('x');
plt.grid(True);
plt.ylabel('Displacement [m]');
plt.legend()
print("flapwise tip displacement  = %f"%(utip[0,-1]))
print("edgewise tip displacement  = %f"%(utip[1,-1]))
print("torsion  tip displacement  = %f"%(ptip[-1]))

plt.figure(4)
plt.plot(xDisp[2,:],u[5,:]*180/np.pi,'-',c='b',label='twist');
plt.xlabel('x');
plt.grid(True);
plt.ylabel('twist [deg]');
plt.legend()

plt.figure(7)
plt.plot(xDisp[2,:],reactF.T,'-',c='b',label='twist');
plt.xlabel('x');
plt.grid(True);
plt.ylabel('twist [deg]');
plt.legend()

#%% Read output from CBeamDyn

data = np.fromfile("u.bin");
u    = data.reshape((-1,6));
plt.plot(u[:,0],'.-');
