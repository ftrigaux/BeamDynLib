import numpy  as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interp2d

plt.close('all')

from time import time 

from numpy import cos,sin

from pyBeamDyn import PyBeamDyn
from scipy.interpolate import PchipInterpolator

#plt.close('all')

from readBD import loadBladeProperties,inputNodes
from fast_io import readFast


eta,K,M = loadBladeProperties("IEA-15-240-RWT_BeamDyn_blade.dat")
xNodes  = inputNodes("IEA-15-240-RWT_BeamDyn_noPrebend.dat")
#%% ---- Definition of the dynamic problem  -----
dt = 1e-2
nt = 3000
q1 = 1000.0
q2 = 0.0
to = 0.0

#%% Import the loads from file 

# FAero = np.loadtxt("loads/loads_9ms_FAST.txt",skiprows=1);
# Fx = interp1d(FAero[:,0],FAero[:,1],fill_value="extrapolate");
# Fy = interp1d(FAero[:,0],FAero[:,2],fill_value="extrapolate");

# Fx_t = np.fromfile("loads/loads_Fx_9ms_5e-4s.bin").reshape(50,10001);
# Fy_t = np.fromfile("loads/loads_Fy_9ms_5e-4s.bin").reshape(50,10001);
# t = np.linspace(0,1e4*5e-4,int(1e4+1));
# Fx_t = interp2d(t,FAero[:,0],Fx_t);
# Fy_t = interp2d(t,FAero[:,0],Fy_t);

#%% ---- Solve the Dynamic Problem with BeamDyn -- python version -----
nBeam      = 1
L          = 120.0
Rhub       = 0*3.0
omega      = 9*9/120;
inputFile  = "IEA-15-240-RWT_BeamDyn_noPrebend.dat"
DynamicSolve=1
nt_loc     = 1
GlbPos     = np.array([0.0,0.0,Rhub],order='F')
t          = np.linspace(dt,nt*dt*nt_loc,nt)
gravity    = np.array([0.0,0.0,-0.0],order='F')

theta      = 0.0*np.pi/3.0;
RootOri    = np.array([[1.0,0.0,0.0],[0.0,cos(theta),sin(theta)],[0.0,-sin(theta),cos(theta)]],order='F');

tilt       = 0.0 * np.pi/180;
pitch      = 0.0;
#tiltMatrix = np.array([[cos(tilt),sin(tilt),0.0],[-sin(tilt),cos(tilt),0.0],[0.0,0.0,1.0]],order='F');
tiltMatrix = np.array([[cos(tilt),0.0,sin(tilt)],[0.0,1.0,0.0],[-sin(tilt),0.0,cos(tilt)]],order='F');
pitchMatrix = np.array([[cos(pitch),sin(pitch),0.0],[-sin(pitch),cos(pitch),0.0],[0.0,0.0,1.0]],order='F');
omega      = np.array([omega,0,0],order='F')

RootOri = tiltMatrix.T @ RootOri.T

# !! GlbPos must be expressed in the global domain (i.e. rotated with the reference orientation)
GlbPos  = RootOri @ GlbPos
omega   = RootOri @ omega

# Add Pitch to Root orientation
RootOri = pitchMatrix @ RootOri

GlbRotBladeT0=0

def pitch_command(t):
    X = [0.0,3.0+10,7.0+10,8.0+10,100.0]
    Y = [0.0,0.0,0.7,0.7,0.7]
    return PchipInterpolator(X,Y)(t)

# Print the simulation data 
print("\n----- Simulation Data ------\n");
print("Gravity = ", gravity);
print("GlbPos  = ", GlbPos);
print("omega   = ", omega);
print("Rotation Matrix:\n")
print('RootOri');
print(RootOri)

pbd = PyBeamDyn()
pbd.initBeamDyn(nBeam,inputFile,1,nt=nt_loc,dt=dt,omega=omega,DynamicSolve=DynamicSolve,GlbPos=GlbPos,gravity=gravity,RootOri=RootOri,GlbRotBladeT0=GlbRotBladeT0)

# Get the position of the distributed loads
xLoads, xDisp = pbd.getPositions()

# Set the loads
loads = np.zeros((6,pbd.nxLoads));
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
u_save    = np.zeros((nt,6,xDisp[0].size));

nt_slow = nt/2
# Solve the dynamic
print("Solving...");
start = time()
for i in range(nt):
    # Fake Damping
    loads[0,:] = q1 - du[0,:]*50;
    loads[1,:] = q2 - du[1,:]*50;
    pbd.setLoads(loads,1)
    
    pbd.setBC(1, omega, np.zeros(3),0,0.0,pitch_command(t[i]));
    
    print(i)
    pbd.solve(1)
    
    # Preallocate variables to extract displacement
    x,u,du = pbd.getDisplacement(1)
    utip[:,i]  =  u[:3,-1];
    dvtip[i]   = du[0,-1];
    ptip[i]    =  u[5,-1];
    qtip[i]    =  u[3,-1];
    rtip[i]    =  u[4,-1];
    u_save[i]  = np.copy(u)
    
    x,F = pbd.getReactionForce(1)
    RootMr[:,i] = F[3:,0];
    
    
end = time()
print("Done in %1.4e seconds -> %1.4e [s/iteration]"%(end-start,(end-start)/nt))

# Free the data
pbd.freeBeamDyn(1,nBeam)

# theta      = 4.0*np.pi/3.0;
# RootOri    = np.array([[1.0,0.0,0.0],[0.0,cos(theta),sin(theta)],[0.0,-sin(theta),cos(theta)]],order='F');
# pbd.initBeamDyn(nBeam,inputFile,2,nt=nt_loc,dt=dt,omega=omega,DynamicSolve=DynamicSolve,GlbPos=GlbPos,gravity=gravity,RootOri=RootOri,GlbRotBladeT0=1)

# # Get the position of the distributed loads
# xLoads2, xDisp2 = pbd.getPositions()

# # Free the data
# pbd.freeBeamDyn(1,nBeam)
# pbd.freeBeamDyn(2,nBeam)
#%%

# Compute the Root Moment due to gravity analytically
# Mg = integral( m * g * r * cos(twist) *dr )
mmmg = (M[1:,1,1] + M[:-1,1,1])/2.0 * np.linalg.norm(gravity)
Meta = (eta[1:] + eta[:-1])/2.0

ctwist = np.cos(interp1d(xNodes[:,2],xNodes[:,3]*np.pi/180.0)(Meta*117));
stwist = np.sin(interp1d(xNodes[:,2],xNodes[:,3]*np.pi/180.0)(Meta*117));

Mg_x = np.trapz(mmmg*Meta*117*stwist,x=Meta*117);
Mg_y = np.trapz(mmmg*Meta*117*ctwist,x=Meta*117);

if np.linalg.norm(omega) == 0:
    t = np.linspace(dt,nt*dt*nt_loc,nt) ;
    thb = theta * np.ones(t.shape);
    trot= 1
else:
    t = np.linspace(dt,nt*dt*nt_loc,nt) * np.linalg.norm(omega) / 2.0 / np.pi;
    thb  = t * 2.0 * np.pi;
    trot = 2.0 * np.pi / np.linalg.norm(omega)

plt.figure()
plt.plot(xNodes[:,2],xNodes[:,3]*np.pi/180.0)
plt.plot(Meta*117,interp1d(xNodes[:,2],xNodes[:,3]*np.pi/180.0)(Meta*117))

#%%
if DynamicSolve:
    plt.figure(1)
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    plt.plot(t,utip.T);
    #plt.plot(t,wtip,'--',c='b',label="Edgewise");
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.legend()
    plt.ylabel("Tip displacement [m]")
    
    plt.figure(2)
    plt.plot(t,ptip*180.0/np.pi,'-',c='r');
    #plt.plot(t,u_save[:,5,20]*180.0/np.pi,'-',c='r');
    plt.plot(t,pitch_command(t*trot),'k--')
    # plt.plot(t,qtip*180.0/np.pi,'-',c='g');
    # plt.plot(t,rtip*180.0/np.pi,'-',c='b');
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.ylabel('Tip angles [deg]');

    plt.figure(5)
    plt.plot(t,RootMr[:,:].T,'-');
    plt.plot(t,Mg_x*np.sin(thb),'k--');
    plt.plot(t,Mg_y*np.sin(thb),'k:');
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

# plt.figure(7)
# plt.plot(xDisp[2,:],reactF.T,'-',c='b',label='twist');
# plt.xlabel('x');
# plt.grid(True);
# plt.ylabel('twist [deg]');
# plt.legend()

#%% Read output from CBeamDyn

# data = np.fromfile("u.bin");
# u    = data.reshape((-1,6));
# plt.plot(u[:,0],'.-');


#%% Read Fast Output
if 0:
    FastData = readFast("IEA-15-240-RWT_BeamDyn_noPrebend_01.out")
    plt.figure()
    plt.plot(t,FastData['N001_MFbxr'],'-');
    plt.plot(t,FastData['N001_MFbyr'],'-');
    plt.plot(t,FastData['N001_MFbzr'],'-');
    plt.title("Inertial")
    plt.figure()
    plt.plot(t,FastData['N001_MFcxr'],'-');
    plt.plot(t,FastData['N001_MFcyr'],'-');
    plt.plot(t,FastData['N001_MFczr'],'-');
    plt.plot(t,Mg_x*np.sin(thb),'k--');
    plt.plot(t,Mg_y*np.sin(thb),'k:');
    plt.title("Elastic")
    plt.figure()
    plt.plot(t,FastData['N001_MFgxr'],'-');
    plt.plot(t,FastData['N001_MFgyr'],'-');
    plt.plot(t,FastData['N001_MFgzr'],'-');
    plt.plot(t,Mg_x*np.sin(thb),'k--');
    plt.plot(t,Mg_y*np.sin(thb),'k:');
    plt.title("Gravity")
    plt.figure()
    plt.plot(t,FastData['N001_MFdxr'],'-');
    plt.plot(t,FastData['N001_MFdyr'],'-');
    plt.plot(t,FastData['N001_MFdzr'],'-');
    plt.title("Fd")
    plt.figure()
    plt.plot(t,FastData['N001_MFixr'],'-');
    plt.plot(t,FastData['N001_MFiyr'],'-');
    plt.plot(t,FastData['N001_MFizr'],'-');
    plt.title("Inertial")
    plt.grid(True)