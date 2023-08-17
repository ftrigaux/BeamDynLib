import sys
import numpy  as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interp2d

plt.close('all')

from time import time 

from numpy import cos,sin

sys.path.append('/Users/ftrigaux/Documents/Beams/BeamDynLib')
sys.path.append('/Users/ftrigaux/Documents/BFtools')

from pyBeamDyn import PyBeamDyn
from scipy.interpolate import PchipInterpolator

#plt.close('all')

from readBD import loadBladeProperties,inputNodes
from fast_io import readFast


eta,K,M = loadBladeProperties("IEA-15-240-RWT_BeamDyn_blade.dat")
xNodes  = inputNodes("IEA-15-240-RWT_BeamDyn_noPrebend.dat")
#%% ---- Definition of the dynamic problem  -----
dt = 2e-2
nt = 400
q1 = 0.0
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
nBeam      = 2
L          = 120.0
Rhub       = 3.0
inputFile  = "IEA-15-240-RWT_BeamDyn_noPrebend.dat"
DynamicSolve=1
nt_loc     = 1
GlbPos     = np.array([0.0,0.0,Rhub],order='F')
t          = np.linspace(dt,nt*dt*nt_loc,nt)
gravity    = np.array([0.0,0.0,-0.0],order='F')

u_save = [0 for i in range(nBeam)]

def pitch_command(t):
    X = [0.0,3.0,7.0,8.0,100.0]
    Y = [0.0,0.0,0.7,0.7,0.7]
    return PchipInterpolator(X,Y)(t)

for iBeam in range(nBeam):
    theta      = iBeam*np.pi/3.0;
    RootOri    = np.array([[1.0,0.0,0.0],[0.0,cos(theta),sin(theta)],[0.0,-sin(theta),cos(theta)]],order='F');

    tilt       = 50.0 * np.pi/180;
    pitch      = 0.0;
    #tiltMatrix = np.array([[cos(tilt),sin(tilt),0.0],[-sin(tilt),cos(tilt),0.0],[0.0,0.0,1.0]],order='F');
    tiltMatrix = np.array([[cos(tilt),0.0,sin(tilt)],[0.0,1.0,0.0],[-sin(tilt),0.0,cos(tilt)]],order='F');
    pitchMatrix = np.array([[cos(pitch),sin(pitch),0.0],[-sin(pitch),cos(pitch),0.0],[0.0,0.0,1.0]],order='F');
    omega      = np.array([9*9/120.0,0,0],order='F')

    RootOri = tiltMatrix.T @ RootOri.T

    # !! GlbPos must be expressed in the global domain (i.e. rotated with the reference orientation)
    GlbPos  = RootOri @ GlbPos
    omega   = RootOri @ omega

    # Add Pitch to Root orientation
    RootOri = pitchMatrix @ RootOri

    GlbRotBladeT0=1

    # Print the simulation data 
    print("\n----- Simulation Data ------\n");
    print("Gravity = ", gravity);
    print("GlbPos  = ", GlbPos);
    print("omega   = ", omega);
    print("Rotation Matrix:\n")
    print('RootOri');
    print(RootOri)

    pbd = PyBeamDyn()
    pbd.initBeamDyn(nBeam,inputFile,iBeam+1,nt=nt_loc,dt=dt,omega=omega,DynamicSolve=DynamicSolve,GlbPos=GlbPos,gravity=gravity,RootOri=RootOri,GlbRotBladeT0=GlbRotBladeT0)

    # Get the position of the distributed loads
    xLoads, xDisp = pbd.getPositions()

    # Set the loads
    loads = np.zeros((6,pbd.nxLoads));
    loads[0,:] = q1;
    pbd.setLoads(loads,iBeam+1)

    # Preallocate variables to extract displacement
    utip      = np.zeros((3,nt))
    dvtip     = np.zeros(nt)
    wtip      = np.zeros(nt)
    ptip      = np.zeros(nt)
    qtip      = np.zeros(nt)
    rtip      = np.zeros(nt)
    RootMr   = np.zeros((3,nt));
    du   = np.zeros((6,xLoads[0].size));
    u_save[iBeam]    = np.zeros((nt,6,xDisp[0].size));

    nt_slow = nt/2
    # Solve the dynamic
    print("Solving...");
    start = time()
    for i in range(nt):
        # Fake Damping
        loads[0,:] = q1 - du[0,:]*50;
        loads[1,:] = q2 - du[1,:]*50;
        #pbd.setLoads(loads,iBeam+1)
        
        pbd.setBC(iBeam+1, omega, np.zeros(3),0,0.0,pitch_command(t[i]));
        
        print(i)
        pbd.solve(iBeam+1)
        
        # Preallocate variables to extract displacement
        x,u,du = pbd.getDisplacement(iBeam+1)
        utip[:,i]  =  u[:3,-1];
        dvtip[i]   = du[0,-1];
        ptip[i]    =  u[5,-1];
        qtip[i]    =  u[3,-1];
        rtip[i]    =  u[4,-1];
        u_save[iBeam][i]  = np.copy(u)
        
        x,F = pbd.getReactionForce(iBeam+1)
        RootMr[:,i] = F[3:,0];
        
        
    end = time()
    print("Done in %1.4e seconds -> %1.4e [s/iteration]"%(end-start,(end-start)/nt))

# Free the data
pbd.freeBeamDyn(iBeam+1,nBeam)

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

# if np.linalg.norm(omega) == 0:
#     t = np.linspace(dt,nt*dt*nt_loc,nt) ;
#     thb = theta * np.ones(t.shape);
#     trot= 1
# else:
#     t = np.linspace(dt,nt*dt*nt_loc,nt) * np.linalg.norm(omega) / 2.0 / np.pi;
#     thb  = t * 2.0 * np.pi;
#     trot = 2.0 * np.pi / np.linalg.norm(omega)

# t = np.linspace(dt,nt*dt*nt_loc,nt) ;

plt.figure()
plt.plot(xNodes[:,2],xNodes[:,3]*np.pi/180.0)
plt.plot(Meta*117,interp1d(xNodes[:,2],xNodes[:,3]*np.pi/180.0)(Meta*117))

#%%
if DynamicSolve:
    plt.figure()
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    plt.plot(t,u_save[0][:,0,50]);
    plt.plot(t,u_save[1][:,0,50]);
    #plt.plot(t,wtip,'--',c='b',label="Edgewise");
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.legend()
    plt.ylabel("Tip displacement [m]")

    plt.figure()
    plt.plot(t,u_save[0][:,5,50],'-',c='b');
    plt.plot(t,u_save[1][:,5,50],'-',c='r');
    #plt.plot(t,u_save[:,5,20]*180.0/np.pi,'-',c='r');
    #plt.plot(t,pitch_command(t),'k--')
    # plt.plot(t,qtip*180.0/np.pi,'-',c='g');
    # plt.plot(t,rtip*180.0/np.pi,'-',c='b');
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.ylabel('Tip angles [deg]');

if 0:
    plt.figure()
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