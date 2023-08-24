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
nt = 100
q1 = 2000.0
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
nBeam      = 3
L          = 120.0
Rhub       = 3.0
inputFile  = "IEA-15-240-RWT_BeamDyn_noPrebend.dat"
DynamicSolve=1
nt_loc     = 1
GlbPos     = np.array([0.0,0.0,Rhub],order='F')
t          = np.linspace(dt,nt*dt*nt_loc,nt)
gravity    = np.array([0.0,0.0,-0.0],order='F')

u_save = [0 for i in range(nBeam)]

leg = ["dt=%1.3f"%(dt/2),"dt=%1.3f then %1.3f"%(dt/2,dt),"dt=%1.3f"%(dt)]

for iBeam in range(nBeam):
    theta      = iBeam*np.pi/3.0;
    RootOri    = np.array([[1.0,0.0,0.0],[0.0,cos(theta),sin(theta)],[0.0,-sin(theta),cos(theta)]],order='F');

    tilt       = 0.0 * np.pi/180;
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
    if iBeam==0:
        pbd.initBeamDyn(nBeam,inputFile,iBeam+1,nt=nt_loc*2,dt=dt/2,omega=omega,DynamicSolve=DynamicSolve,GlbPos=GlbPos,gravity=gravity,RootOri=RootOri,GlbRotBladeT0=GlbRotBladeT0)
    else:
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

    # Solve the dynamic
    print("Solving...");
    start = time()
    for i in range(nt):
        # Fake Damping
        #loads[0,:] = q1 - du[0,:]*50;
        #loads[1,:] = q2 - du[1,:]*50;
        #pbd.setLoads(loads,iBeam+1)

        if i==0:
            pbd.refresh(2,dt=dt/2,nt=nt_loc*2)
        elif i==nt//2:
            pbd.refresh(2,dt=dt,nt=nt_loc)
        
        pbd.setBC(iBeam+1, omega, np.zeros(3),0,0.0,0.0);
        
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

#%%
if DynamicSolve:
    plt.figure()
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    for i in range(nBeam):
        plt.plot(t,u_save[i][:,0,50]);
    #plt.plot(t,wtip,'--',c='b',label="Edgewise");
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.legend(leg)
    plt.ylabel("Tip displacement [m]")

    plt.figure()
    for i in range(nBeam):
        plt.plot(t,u_save[i][:,5,50]);
    #plt.plot(t,u_save[:,5,20]*180.0/np.pi,'-',c='r');
    #plt.plot(t,pitch_command(t),'k--')
    # plt.plot(t,qtip*180.0/np.pi,'-',c='g');
    # plt.plot(t,rtip*180.0/np.pi,'-',c='b');
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.legend(leg)
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