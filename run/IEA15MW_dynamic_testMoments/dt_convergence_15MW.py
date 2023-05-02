import numpy  as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interp2d

from time import time 

from pyBeamDyn import PyBeamDyn

plt.close('all')

#%% ---- Definition of the dynamic problem  -----
dt = 4e-2
nt = 100

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
inputFile  = "IEA-15-240-RWT_BeamDyn_testNoRefine.dat"
DynamicSolve=1
nt_loc     = 1
GlbPos     = np.array([0.0,0.0,Rhub],order='F')
t          = np.linspace(dt,nt*dt*nt_loc,nt)

if DynamicSolve==0:
    nt = 2

ndt = 5; # number of time-steps size

t          = [None]*ndt
lvtip      = [None]*ndt
lutip      = [None]*ndt
ldvtip     = [None]*ndt
lwtip      = [None]*ndt
lptip      = [None]*ndt
lqtip      = [None]*ndt
lrtip      = [None]*ndt

for it in range(ndt):
    nt = int(nt*2.0)
    dt =     dt/2.0

    t[it] = np.linspace(dt,nt*dt*nt_loc,nt);

    lvtip[it]      = np.zeros(nt)
    lutip[it]      = np.zeros(nt)
    ldvtip[it]     = np.zeros(nt)
    lwtip[it]      = np.zeros(nt)
    lptip[it]      = np.zeros(nt)
    lqtip[it]      = np.zeros(nt)
    lrtip[it]      = np.zeros(nt)
    
    pbd = PyBeamDyn()
    pbd.initBeamDyn(1,inputFile,1,nt=nt_loc,dt=dt,omega=np.array([omega*1,0,0],order='F'),DynamicSolve=DynamicSolve,GlbPos=GlbPos)
    
    # Get the position of the distributed loads
    xLoads, xDisp = pbd.getPositions()
    
    # Set the loads
    loads = np.zeros((6,pbd.nxLoads));
    loads[5,:] = -5000;
    pbd.setLoads(loads,1)
    
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
        lvtip[it][i]  =  u[0,-1];
        lwtip[it][i]  =  u[1,-1];
        lutip[it][i]  =  u[2,-1];
        ldvtip[it][i] = du[0,-1];
        lptip[it][i]  =  u[5,-1];
        lqtip[it][i]  =  u[3,-1];
        lrtip[it][i]  =  u[4,-1];
    end = time()
    print("Done in %1.4e seconds -> %1.4e [s/iteration]"%(end-start,(end-start)/nt))

    # Free the data
    pbd.freeBeamDyn(1,1)

#%%
for it in range(ndt):
    plt.figure(1)
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    plt.plot(t[it],lvtip[it],'--',c='g',label="Flapwise");
    #plt.plot(t,wtip,'--',c='b',label="Edgewise");
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.legend()
    plt.ylabel("Tip displacement [m]")
    
    plt.figure(2)
    plt.plot(t[it],lptip[it]*180.0/np.pi,'-',label='dt = %1.3e'%(t[it][0]));
    plt.xlabel('Time [s]');
    plt.grid(True);
    plt.ylabel('Tip angles [deg]');
    plt.legend()

plt.figure(2)
plt.savefig('./figure/BeamDyn_TimeConvergence_IEA15MW_rhoinf11.pdf',format='pdf');