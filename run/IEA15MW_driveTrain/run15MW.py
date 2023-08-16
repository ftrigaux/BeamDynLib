import numpy  as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interp2d

plt.close('all')

from time import time 

from numpy import cos,sin
from scipy.signal import welch

import sys
sys.path.append('/Users/ftrigaux/Documents/Beams/BeamDynLib');
sys.path.append('/Users/ftrigaux/Documents/BFtools');
from pyBeamDyn import PyBeamDyn

#plt.close('all')

from readBD import loadBladeProperties,inputNodes
from fast_io import readFast


eta,K,M = loadBladeProperties("IEA-15-240-RWT_BeamDyn_blade.dat")
xNodes  = inputNodes("IEA-15-240-RWT_BeamDyn.dat")
#%% ---- Definition of the dynamic problem  -----
dt = 1e-3
nt = 2000
q1 = 0.0
q2 = 300.0 
to = 0.0

#%% 
# Drive train inertia - computed from BeamDyn and OpenFAST data

Mblade = 68507.600 # elastodyn (kg-m^2)
Iblade = 101086645.135 # elastodyn (kg-m^2)
Ihub   = 973520
Igen   = 1836784

HubRad = 3 # m
nBlade = 3

IdriveTrain = nBlade * (Mblade * HubRad**2 + Iblade) + Ihub + Igen
IdriveTrain = 973520.000 + 1836784

Irotor = 340119198.835; # from elastodyn, without cone and tilt(kg-m^2)

#%% Import the loads from file 


#%% ---- Solve the Dynamic Problem with BeamDyn -- python version -----

omega_vec = np.zeros((3,nt+1),order='F');
dOmega    = np.zeros((3,nt+1),order='F');
theta     = np.zeros((3,nt+1),order='F');

nBeam      = 3

theta0 = np.zeros(nBeam)

for iBeam in range(nBeam):
    L          = 120.0
    Rhub       = 3.0
    omega      = 9*9/120;
    omega      = np.array([omega,0,0],order='F')
    inputFile  = "IEA-15-240-RWT_BeamDyn.dat"
    DynamicSolve=1
    nt_loc     = 1
    GlbPos     = np.array([0.0,0.0,Rhub],order='F')
    t          = np.linspace(dt,nt*dt*nt_loc,nt)
    gravity    = np.array([0.0,0.0,-0.0],order='F')

    theta0[iBeam]      = 2.0*iBeam*np.pi/3.0;
    theta[iBeam,0]     = 2.0*iBeam*np.pi/3.0;
    RootOri    = np.array([[1.0,0.0,0.0],[0.0,cos(theta0[iBeam]),sin(theta0[iBeam])],[0.0,-sin(theta0[iBeam]),cos(theta0[iBeam])]],order='F');

    tilt       = 0.0 * np.pi/180;
    pitch      = 0.0;
    #tiltMatrix = np.array([[cos(tilt),sin(tilt),0.0],[-sin(tilt),cos(tilt),0.0],[0.0,0.0,1.0]],order='F');
    tiltMatrix = np.array([[cos(tilt),0.0,sin(tilt)],[0.0,1.0,0.0],[-sin(tilt),0.0,cos(tilt)]],order='F');
    pitchMatrix = np.array([[cos(pitch),sin(pitch),0.0],[-sin(pitch),cos(pitch),0.0],[0.0,0.0,1.0]],order='F');

    RootOri = tiltMatrix.T @ RootOri.T

    # !! GlbPos must be expressed in the global domain (i.e. rotated with the reference orientation)
    GlbPos  = RootOri @ GlbPos
    omega   = RootOri @ omega
    omega_vec[:,0] = omega[:]

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
    print("init")
    pbd.initBeamDyn(nBeam,inputFile,iBeam+1,nt=nt_loc,dt=dt,omega=omega,DynamicSolve=DynamicSolve,GlbPos=GlbPos,gravity=gravity,RootOri=RootOri,GlbRotBladeT0=GlbRotBladeT0)

    # Get the position of the distributed loads
    xLoads, xDisp = pbd.getPositions()

    x,u,du = pbd.getDisplacement(iBeam+1)
    


# Preallocate variables to extract displacement
utip      = np.zeros((6,nt))
RootFM    = np.zeros((6,nt))
Mhub      = np.zeros(nt)
Fy        = np.zeros(nt)

# Solve the dynamic
print("Solving...");
start = time()
for i in range(nt):
    
    #pbd.setBC(1, omega * (nt_slow-min(i,nt_slow))/nt_slow, omega*(-1/nt_slow))
    Mgen = (q2 * (L-Rhub)**2.0/2.0 + q2*L*Rhub) * nBeam;
    for iBeam in range(nBeam):
        # Set the loads
        loads = np.zeros((6,pbd.nxLoads));
        #loads[0,:] = 1e4 * xLoads[2,:]/L; 
        loads[1,:] = -q2; 

        # Fake Damping
        #loads[0:3,:] = loads[0:3,:] - du[0:3,:]*100;

        pbd.setLoads(loads,iBeam+1)
        
        print(i)
        pbd.solve(iBeam+1)
        
        # Preallocate variables to extract displacement
        x,u,du = pbd.getDisplacement(iBeam+1)
        utip[:,i] = u[-6:,-1]
        
        x,F = pbd.getReactionForce(iBeam+1)

        # Send loads to driveTrain: 
        Mhub[i] = Mhub[i] + (F[3,0] - Rhub * F[1,0]);
    
        # Force in edgewise direction
        Fy[i] = F[1,0] * cos(theta[iBeam,i]);
        print(cos(theta[iBeam,i]))

    
    # Solve drive train dynamics
    dOmega[:,i+1] = dt * (Mhub[i] - Mgen) * np.array([1,0,0]) / IdriveTrain;

    ##dOmega[:,i+1] = np.zeros(3,order='F'); #TODO !! ! ! 

    omega_vec[:,i+1] = omega_vec[:,i] + dOmega[:,i+1];

    for iBeam in range(nBeam):
        theta[iBeam,i+1] =  theta[iBeam,i] + omega_vec[0,i+1] * dt;
        print('theta b%d : %1.15f'%(iBeam,theta[iBeam,i+1]))
    
        pbd.setBC(iBeam+1,omega=omega_vec[:,i+1],dOmega=dOmega[:,i+1], forceTheta=1, theta_rot=theta[iBeam,i+1]-theta0[iBeam]);

    print("Omega_new  = %f"%omega_vec[0,i+1]);
    print("dOmega     = %f"%dOmega[0,i+1]);

    RootFM[:,i] = F[:,0];
    
    
end = time()
print("Done in %1.4e seconds -> %1.4e [s/iteration]"%(end-start,(end-start)/nt))

    # Free the data
for iBeam in range(nBeam):
    print('Free Beam %d'%(iBeam))
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

if np.linalg.norm(omega) == 0:
    t = np.linspace(dt,nt*dt*nt_loc,nt) ;
else:
    t = np.linspace(dt,nt*dt*nt_loc,nt) #* np.linalg.norm(omega) / 2.0 / np.pi;
    thb  = t * 2.0 * np.pi;

plt.figure()
plt.plot(xNodes[:,2],xNodes[:,3]*np.pi/180.0)
plt.plot(Meta*117,interp1d(xNodes[:,2],xNodes[:,3]*np.pi/180.0)(Meta*117))

#%%

## np.savetxt("RootFM_Igen.dat",RootFM)
## np.savetxt("utip_Igen.dat",utip)
## np.savetxt("Mhub_Igen.dat",Mhub)
## np.savetxt("omega_vec_Igen.dat",omega_vec)
#%%
if DynamicSolve:
    plt.figure(1)
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    plt.plot(t,utip[:2,:].T);
    #plt.plot(t,wtip,'--',c='b',label="Edgewise");
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.legend()
    plt.ylabel("Tip displacement [m]")
    
    plt.figure(2)
    plt.plot(t,utip[5,:]*180.0/np.pi,'-',c='r');
    # plt.plot(t,qtip*180.0/np.pi,'-',c='g');
    # plt.plot(t,rtip*180.0/np.pi,'-',c='b');
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.ylabel('Tip angles [deg]');

    # plt.figure(5)
    # plt.plot(t,RootMr[:,:].T,'-');
    # plt.plot(t,Mg_x*np.sin(thb),'k--');
    # plt.plot(t,Mg_y*np.sin(thb),'k:');
    # plt.xlabel('Time [s]');
    # plt.grid(True);
    # plt.ylabel('Root Moments [Nm]');

plt.figure(3)
plt.plot(xDisp[2,:],u[0,:],'-',c='b',label='Flap');
plt.plot(xDisp[2,:],u[1,:],'-',c='r',label='Edge');
plt.xlabel('x');
plt.grid(True);
plt.ylabel('Displacement [m]');
plt.legend()
print("flapwise tip displacement  = %f"%(utip[0,-1]))
print("edgewise tip displacement  = %f"%(utip[1,-1]))
print("torsion  tip displacement  = %f"%(utip[5,-1]))

plt.figure(4)
plt.plot(xDisp[2,:],u[5,:]*180/np.pi,'-',c='b',label='twist');
plt.xlabel('x');
plt.grid(True);
plt.ylabel('twist [deg]');
plt.legend()

plt.figure()
plt.plot(t,omega_vec[0,:-1]);
plt.xlabel('t')
plt.ylabel('omega')
plt.grid(True);

plt.figure()
plt.plot(t,theta[0,:-1]);
plt.xlabel('t')
plt.ylabel('theta')
plt.grid(True);

plt.figure()
plt.plot(t,Mhub);
plt.xlabel('t')
plt.ylabel('M_hub')
plt.grid(True);

plt.figure()
plt.plot(t,RootFM[3,:].T);
plt.xlabel('Time [s]');
plt.grid(True);
plt.ylabel("RootFM");

plt.figure()
plt.plot(t,utip[0,:].T);
plt.plot(t,utip[1,:].T);
plt.xlabel('Time [s]');
plt.grid(True);
plt.ylabel("u_y");

plt.figure()
plt.plot(t,Fy);
plt.grid(True);

plt.figure()
edge1freq = 0.642;
freq, PQ = welch(Mhub,fs=1/dt,nperseg=Mhub.size/1);
freq_1rot = omega_vec[0,0]/2/np.pi
freq_rot = freq / (omega_vec[0,0]/2/np.pi);
plt.axvline(x=3,c='k')
plt.axvline(edge1freq/ (omega_vec[0,0]/2/np.pi),c='k',linestyle='--');
plt.loglog(freq_rot,PQ,c='black' ,label='Mhub');

freq, PQ = welch(Fy,fs=1/dt,nperseg=Fy.size/1);
plt.loglog(freq_rot,PQ,c='red' ,label='Fy');

freq, PQ = welch(utip[1,:],fs=1/dt,nperseg=utip[1,:].size/1);
plt.loglog(freq_rot,PQ,c='blue' ,label='utip');

plt.legend()
plt.xlim([1,100])
plt.grid(True)


#%%
import sys
sys.path.append('~/Documents/BFtools/')
from fast_io import readFast

FastData = readFast("/Users/ftrigaux/Documents/IEA15run/OpenFast/15MW_inertia/IEA15MW_AeroDyn.out")
t_rot = FastData['Time'] 
lastNode = 51;
dt_fst = FastData['Time'][1] - FastData['Time'][0]

#%%
# -- Time series of the tip displacement
plt.figure()
plt.plot(t,utip[0,:],label='flap BD');
plt.plot(t,utip[1,:],label='edge BD');
plt.plot(t_rot,FastData['B1N%03d_TDxr'%(lastNode)],'--',label='flap FAST');
plt.plot(t_rot,FastData['B1N%03d_TDyr'%(lastNode)],'--',label='edge FAST');
plt.xlabel('Time [s]');
plt.grid(True);
plt.legend()
plt.ylabel("Tip displacement [m]")


# -- Omega
fig= plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(t,omega_vec[0,:-1],label='omega BD');
ax1.plot(t_rot,FastData['RotSpeed']/30*np.pi,label='omega FAST');
ax1.set_xlabel('t')
ax1.set_ylabel('Rotation Speed (rad/s)');
ax1.grid(True);

# -- PSD of root moment
plt.figure()
edge1freq = 0.642;
freq, PQ = welch(Mhub,fs=1/dt,nperseg=Mhub.size/1);
freq_1rot = omega_vec[0,0]/2/np.pi
freq_rot = freq / (omega_vec[0,0]/2/np.pi);
plt.axvline(x=3,c='k')
plt.axvline(edge1freq/ (omega_vec[0,0]/2/np.pi),c='k',linestyle='--');
plt.loglog(freq_rot,PQ,c='black' ,label='Mhub BD');

freq, PQ = welch(FastData['RotTorq']*1e3,fs=1/dt_fst,nperseg=FastData['RotTorq'].size/1);
freq_rot = freq / freq_1rot;
plt.loglog(freq_rot,PQ,'r--',label='RotTorq FAST');



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
# FastData = readFast("IEA-15-240-RWT_BeamDyn_noPrebend_01.out")
# plt.figure()
# plt.plot(t,FastData['N001_MFbxr'],'-');
# plt.plot(t,FastData['N001_MFbyr'],'-');
# plt.plot(t,FastData['N001_MFbzr'],'-');
# plt.title("Inertial")
# plt.figure()
# plt.plot(t,FastData['N001_MFcxr'],'-');
# plt.plot(t,FastData['N001_MFcyr'],'-');
# plt.plot(t,FastData['N001_MFczr'],'-');
# plt.plot(t,Mg_x*np.sin(thb),'k--');
# plt.plot(t,Mg_y*np.sin(thb),'k:');
# plt.title("Elastic")
# plt.figure()
# plt.plot(t,FastData['N001_MFgxr'],'-');
# plt.plot(t,FastData['N001_MFgyr'],'-');
# plt.plot(t,FastData['N001_MFgzr'],'-');
# plt.plot(t,Mg_x*np.sin(thb),'k--');
# plt.plot(t,Mg_y*np.sin(thb),'k:');
# plt.title("Gravity")
# plt.figure()
# plt.plot(t,FastData['N001_MFdxr'],'-');
# plt.plot(t,FastData['N001_MFdyr'],'-');
# plt.plot(t,FastData['N001_MFdzr'],'-');
# plt.title("Fd")
# plt.figure()
# plt.plot(t,FastData['N001_MFixr'],'-');
# plt.plot(t,FastData['N001_MFiyr'],'-');
# plt.plot(t,FastData['N001_MFizr'],'-');
# plt.title("Inertial")
# plt.grid(True)