import numpy  as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interp2d

from time import time 

from pyBeamDyn import PyBeamDyn

#plt.close('all')

#%% ---- Definition of the dynamic problem  -----
dt = 2e-2
nt = 400
q1 = 0.0
q2 = 0.0
to = 0.0

#%% Import the loads from file 

FAero = np.loadtxt("loads/loads_9ms_FAST.txt",skiprows=1);
Fx = interp1d(FAero[:,0],FAero[:,1],fill_value="extrapolate");
Fy = interp1d(FAero[:,0],FAero[:,2],fill_value="extrapolate");

Fx_t = np.fromfile("loads/loads_Fx_9ms_5e-4s.bin").reshape(50,10001);
Fy_t = np.fromfile("loads/loads_Fy_9ms_5e-4s.bin").reshape(50,10001);
t = np.linspace(0,1e4*5e-4,int(1e4+1));
Fx_t = interp2d(t,FAero[:,0],Fx_t);
Fy_t = interp2d(t,FAero[:,0],Fy_t);

#%% ---- Solve the Dynamic Problem with BeamDyn -- python version -----
nBeam      = 1
L          = 120.1
Rhub       = 3.0
omega      = 9*9/120;
inputFile  = "IEA-15-240-RWT_BeamDyn.dat"
DynamicSolve=1
nt_loc     = 1
GlbPos     = np.array([0.0,0.0,Rhub],order='F')
t          = np.linspace(dt,nt*dt*nt_loc,nt)

if DynamicSolve==0:
    nt = 2

pbd = PyBeamDyn()
pbd.initBeamDyn(1,inputFile,1,nt=nt_loc,dt=dt,omega=np.array([omega*1,0,0],order='F'),DynamicSolve=DynamicSolve,GlbPos=GlbPos)

# Get the position of the distributed loads
xLoads, xDisp = pbd.getPositions()

# Set the loads
loads = np.zeros((6,pbd.nxLoads));
#loads[0,:] = 5000;
#loads[1,:] = 800;
# loads[5,-1] = 5000;
loads[5,:] = -5000;
pbd.setLoads(loads,1)

# Preallocate variables to extract displacement
vtip      = np.zeros(nt)
utip      = np.zeros(nt)
dvtip     = np.zeros(nt)
wtip      = np.zeros(nt)
ptip      = np.zeros(nt)
qtip      = np.zeros(nt)
rtip      = np.zeros(nt)
    
# Solve the dynamic
print("Solving...");
start = time()
for i in range(nt):
    #loads[0,:] =  Fx_t(t[i],(xLoads[2,:]+Rhub)/L)[:,0];
    #loads[1,:] =  -Fy_t(t[i],(xLoads[2,:]+Rhub)/L)[:,0];
    #pbd.setLoads(loads,1)
    
    print(i)
    pbd.solve(1)
    
    # Preallocate variables to extract displacement
    x,u,du = pbd.getDisplacement(1)
    vtip[i]  =  u[0,-1];
    wtip[i]  =  u[1,-1];
    utip[i]  =  u[2,-1];
    dvtip[i] = du[0,-1];
    ptip[i]  =  u[5,-1];
    qtip[i]  =  u[3,-1];
    rtip[i]  =  u[4,-1];
end = time()
print("Done in %1.4e seconds -> %1.4e [s/iteration]"%(end-start,(end-start)/nt))

# Free the data
pbd.freeBeamDyn(1,1)
#%%

if DynamicSolve:
    plt.figure(1)
    t = np.linspace(dt,nt*dt*nt_loc,nt) * omega / 2.0 / np.pi;
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    plt.plot(t,vtip,'--',c='g',label="Flapwise");
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
    
    print("flapwise tip displacement  = %f"%(vtip[-1]))
    print("edgewise tip displacement  = %f"%(wtip[-1]))
    print("torsion  tip displacement  = %f"%(ptip[-1]))

plt.figure(3)
plt.plot(xDisp[2,:],u[0,:],'-',c='b',label='Flap');
plt.plot(xDisp[2,:],u[1,:],'-',c='r',label='Edge');
plt.xlabel('x');
plt.grid(True);
plt.ylabel('Displacement [m]');
plt.legend()
print("flapwise tip displacement  = %f"%(vtip[-1]))
print("edgewise tip displacement  = %f"%(wtip[-1]))

plt.figure(3)
plt.plot(xDisp[2,:],u[5,:]*180/np.pi,'-',c='b',label='twist');
plt.xlabel('x');
plt.grid(True);
plt.ylabel('twist [deg]');
plt.legend()


import scipy.fftpack
yf = scipy.fftpack.fft(ptip)
xf = scipy.fftpack.fftfreq(nt,d=dt)
plt.figure()
plt.plot(xf[:nt//2], 2.0/nt * np.abs(yf[:nt//2]))

#%% Read output from CBeamDyn

# data = np.fromfile("u.bin");
# u    = data.reshape((-1,6));
# plt.plot(u[:,0],'.-');
