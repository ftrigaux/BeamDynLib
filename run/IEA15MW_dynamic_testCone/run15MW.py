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
xNodes  = inputNodes("IEA-15-240-RWT_BeamDyn.dat")
#%% ---- Definition of the dynamic problem  -----
dt = 5e-2
nt = 400
q1 = 0.0
q2 = 1000.0
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

def pitch_command(t):
    #return np.pi/2.0 * t/10.0;
    return 0.0;

#%% ---- Solve the Dynamic Problem with BeamDyn -- python version -----
nBeam      = 1
L          = 120.0
Rhub       = 3.0
inputFile  = "IEA-15-240-RWT_BeamDyn.dat"
DynamicSolve=1
nt_loc     = 1
GlbPos     = np.array([0.0,0.0,Rhub],order='F')
t          = np.linspace(dt,nt*dt*nt_loc,nt)
gravity    = np.array([0.0,0.0,-0.0],order='F')

u_save = [0 for i in range(nBeam)]

theta0 = 0.0;
tiltA  = np.deg2rad(0.0);
yawA   = 0.0;
coneAngle = - np.deg2rad(30.0);

def computeRootOri(theta0,tiltA,yawA,coneAngle):
    RootOri       = np.zeros((3,3), order='F')

    RootOri[0][0] = -(sin(theta0) * sin(yawA) - sin(tiltA) * cos(theta0) * cos(yawA)) * sin(coneAngle) + cos(coneAngle) * cos(tiltA) * cos(yawA)
    RootOri[0][1] = -sin(theta0) * sin(tiltA) * cos(yawA) - sin(yawA) * cos(theta0)
    RootOri[0][2] = (sin(theta0) * sin(yawA) - sin(tiltA) * cos(theta0) * cos(yawA)) * cos(coneAngle) + sin(coneAngle) * cos(tiltA) * cos(yawA)

    RootOri[1][0] = -(-sin(theta0) * cos(yawA) - sin(tiltA) * sin(yawA) * cos(theta0)) * sin(coneAngle) + sin(yawA) * cos(coneAngle) * cos(tiltA)
    RootOri[1][1] = -sin(theta0) * sin(tiltA) * sin(yawA) + cos(theta0) * cos(yawA)
    RootOri[1][2] = (-sin(theta0) * cos(yawA) - sin(tiltA) * sin(yawA) * cos(theta0)) * cos(coneAngle) + sin(coneAngle) * sin(yawA) * cos(tiltA)

    RootOri[2][0] = -sin(coneAngle) * cos(theta0) * cos(tiltA) + sin(tiltA) * cos(coneAngle)
    RootOri[2][1] = sin(theta0) * cos(tiltA)
    RootOri[2][2] = sin(coneAngle) * sin(tiltA) + cos(coneAngle) * cos(theta0) * cos(tiltA)

    return RootOri;

for iBeam in range(nBeam):
    theta      = (iBeam+0)*np.pi/3.0;

    omega      = np.array([9*9/120.0,0,0],order='F')

    # !! GlbPos must be expressed in the global domain (i.e. rotated with the reference orientation)
    GlbPos  = computeRootOri(0.0,tiltA,yawA,coneAngle) @ GlbPos
    omega   = computeRootOri(0.0,tiltA,yawA,0.0) @ omega

    RootOri = computeRootOri(theta0+theta,tiltA,yawA,coneAngle);

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
    pbd.initBeamDyn(nBeam,inputFile,iBeam+1,nt=nt_loc,dt=dt,omega=omega,DynamicSolve=DynamicSolve,GlbPos=GlbPos,gravity=gravity,RootOri=RootOri,GlbRotBladeT0=GlbRotBladeT0,WrVTK=2,VTK_fps=10)

    # Get the position of the distributed loads
    xLoads, xDisp = pbd.getPositions()

    # Set the loads
    loads = np.zeros((6,pbd.nxLoads));
    loads[0,:] = q1;
    loads[1,:] = q2;
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
        loads[0,:] = q1 - du[0,:]*50;
        loads[1,:] = q2 - du[1,:]*50;
        #pbd.setLoads(loads,iBeam+1)
        
        pbd.setBC(iBeam+1, omega, np.zeros(3),0,0.0,pitch_command(t));
        
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
    plt.plot(t,u_save[0][:,0,50]);
    plt.plot(t,u_save[0][:,1,50]);
    #plt.plot(t,wtip,'--',c='b',label="Edgewise");
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.legend()
    plt.ylabel("Tip displacement [m]")

    plt.figure()
    plt.plot(t,u_save[0][:,5,50],'-',c='b');
    #plt.plot(t,u_save[:,5,20]*180.0/np.pi,'-',c='r');
    #plt.plot(t,pitch_command(t),'k--')
    # plt.plot(t,qtip*180.0/np.pi,'-',c='g');
    # plt.plot(t,rtip*180.0/np.pi,'-',c='b');
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.ylabel('Tip angles [deg]');

print(np.mean(u_save[0][:,0,-1]))
print(np.mean(u_save[0][:,1,-1]))
if 0:
    plt.figure()
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