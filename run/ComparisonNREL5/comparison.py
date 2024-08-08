import numpy  as np
import matplotlib.pyplot as plt

from beam3D_NU_v2 import createBeamFromCSV
from pyBeamDyn import PyBeamDyn


plt.close('all')

##%matplotlib inline

#%% ---- Definition of the dynamic problem  -----
dt = 0.1
nt = 400
q1 = 10e3
q2 = 10e3
to = 10e2

Omega = 1

#%% ---- Solve the Dynamic Problem with PBeams -------
nx  = 48
L   = 61.5
eps = 0.00477465

bm  = createBeamFromCSV("./nrel5mw_data.csv",nx=nx,L=L,eps=eps)
bm.Omega = Omega
bm.dampK = 0.01

vtip0  = np.zeros(nt)
dvtip0 = np.zeros(nt)
wtip0  = np.zeros(nt)
ptip0  = np.zeros(nt)
RootMr0 = np.zeros((3,nt));
t     = np.linspace(0,dt*nt,nt)

_q1 = q1 * bm.x/bm.L
_q2 = q2 * bm.x/bm.L
for i in range(nt):
    print(i)
    bm.solve_dyn(dt,1,q1=_q1,q2=_q2,t=to)
    vtip0[i] = bm.uv['v'][-1];
    dvtip0[i] = bm.duv['v'][-1];
    wtip0[i] = bm.uv['w'][-1];
    ptip0[i] = bm.uv['phi'][-1]*np.rad2deg(1);
    F1,M1,F2,M2 = bm.getReactionLoads()
    RootMr0[1,i] =  M1
    RootMr0[0,i] =  M2
    
    
# plt.figure(1)
# plt.plot(t,vtip0,c='r');
# plt.plot(t,wtip0,c='g');
# plt.plot(t,ptip0,c='b');
# plt.legend(['Flap','Edge','Twist']);

# plt.figure(2);
# plt.plot(t,RootMr0[1,:].T,'-');

print("Done solving using PBEAMS");

#%% ---- Solve the Dynamic Problem with BeamDyn -- python version -----
nBeam      = 1
#inputFile  = "/Users/ftrigaux/Documents/Beams/BeamDynLib/run/nrel5mw_dynamic/bd_primary_nrel_5mw_dynamic.inp"
inputFile  = "./NRELOffshrBsline5MW_BeamDyn.dat"

pbd = PyBeamDyn()
pbd.initBeamDyn(1,inputFile,1,nt=1,dt=dt,omega=-np.array([Omega,0,0],order='F'),DynamicSolve=1)

# Get the position of the distributed loads
xLoads, xDisp = pbd.getPositions()

# Set the loads
loads = np.zeros((6,pbd.nxLoads),order='F');
loads[0,:] = q1 * (xLoads[2,:]-1)/61.5;
loads[1,:] = q2 * (xLoads[2,:]-1)/61.5;
loads[5,:] = to * (xLoads[2,:]-1)/61.5;
pbd.setLoads(loads,1)

# Preallocate variables to extract displacement
vtip      = np.zeros(nt)
dvtip     = np.zeros(nt)
wtip      = np.zeros(nt)
ptip      = np.zeros(nt)
RootMr   = np.zeros((3,nt));
    
# Solve the dynamic
for i in range(nt):
    pbd.solve(1)
    
    # Preallocate variables to extract displacement
    x,u,du = pbd.getDisplacement(1)
    vtip[i]  =  u[0,-1];
    wtip[i]  =  u[1,-1];
    dvtip[i] = du[0,-1];
    ptip[i]  =  u[5,-1]*np.rad2deg(1);

    x,F = pbd.getReactionForce(1)
    RootMr[:,i] = F[3:,0];

# Free the data
pbd.freeBeamDyn(1,1)

#%%
plt.figure(1)
plt.plot(t,vtip0,c='r');
plt.plot(t,wtip0,c='g');
plt.plot(t,ptip0,c='b');
plt.legend(['Flap','Edge','Twist']);
plt.plot(t,vtip,'--',c='r');
plt.plot(t,wtip,'--',c='g');
plt.plot(t,ptip,'--',c='b');
print("    %-5s\t   %-5s\t   %-5s"%('EB','BD','Diff'));
print("v   %4.2f\t %4.2f\t %4.2f %%"%(np.mean(vtip0),np.mean(vtip),(np.mean(vtip0)-np.mean(vtip))/np.mean(vtip0)*100))
print("w   %4.2f\t %4.2f\t %4.2f %%"%(np.mean(wtip0),np.mean(wtip),(np.mean(wtip0)-np.mean(wtip))/np.mean(wtip0)*100))
print("p   %4.2f\t %4.2f\t %4.2f %%"%(np.mean(ptip0),np.mean(ptip),(np.mean(ptip0)-np.mean(ptip))/np.mean(ptip0)*100))

plt.figure(2);
print("EB Mflap = %1.3f"%(np.mean(RootMr0[1,:])))
print("BD Mflap = %1.3f"%(np.mean(RootMr[1,:])))
print("Diff : %1.3f %%"%(np.mean(RootMr0[1,:]-RootMr[1,:])/np.mean(RootMr0[1,:])*100))

plt.plot(t,RootMr0[1,:].T,'b-');
plt.plot(t,RootMr[1,:].T,'b--');
plt.plot(t,-RootMr0[0,:].T,'g-');
plt.plot(t,RootMr[0,:].T,'g--');

#%% Read output from CBeamDyn

# data = np.fromfile("u.bin");
# u    = data.reshape((-1,6));
# plt.plot(u[:,0],'.-');