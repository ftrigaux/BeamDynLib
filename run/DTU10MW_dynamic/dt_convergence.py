import numpy  as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interp2d

from time import time 

from pyBeamDyn import PyBeamDyn

plt.close('all')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Typewriter"],
    "font.size": 12})

#%% ---- Definition of the dynamic problem  -----
dt = 2e-2
nt = 50
loaddir = 5

#%% Import the loads from file 

#%% ---- Solve the Dynamic Problem with BeamDyn -- python version -----
nBeam      = 1
L          = 8.63655000e+01
Rhub       = 2.8
omega      = 7.5*9/(L+Rhub);
inputFile  = "DTU10MW_BeamDyn_n048.dat"#"../IEA15MW_dynamic/IEA-15-240-RWT_BeamDyn.dat"#
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
    nt_loc = int(nt*2.0**it)
    dt_loc = dt/(2.0**it)

    t[it] = np.linspace(dt_loc,nt_loc*dt_loc,nt_loc);

    lvtip[it]      = np.zeros(nt_loc)
    lutip[it]      = np.zeros(nt_loc)
    ldvtip[it]     = np.zeros(nt_loc)
    lwtip[it]      = np.zeros(nt_loc)
    lptip[it]      = np.zeros(nt_loc)
    lqtip[it]      = np.zeros(nt_loc)
    lrtip[it]      = np.zeros(nt_loc)
    
    pbd = PyBeamDyn()
    pbd.initBeamDyn(1,inputFile,1,nt=1,dt=dt_loc,omega=np.array([omega*0,0,0],order='F'),DynamicSolve=DynamicSolve,GlbPos=GlbPos)
    
    # Get the position of the distributed loads
    xLoads, xDisp = pbd.getPositions()
    
    # Set the loads
    loads = np.zeros((6,pbd.nxLoads));
    loads[loaddir,:] = 4000;
    pbd.setLoads(loads,1)
    
    # Solve the dynamic
    print("Solving...");
    start = time()
    for i in range(nt_loc):
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


fig, ax = plt.subplots(1,3,figsize=(12,3.5))
fig.add_subplot(111, frameon=False)

for it in range(ndt):
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    ax[0].plot(t[it],lvtip[it],'-',label='dt = %1.3e'%(t[it][0]));
    #ax[0].plot(t,wtip,'--',c='b',label="Edgewise");
ax[0].set_xlabel('Time [s]');
ax[0].set_ylabel(r'$x_{flap}$ [m]')
ax[0].grid(True);
ax[0].set_xlim([0,nt*dt])

for it in range(ndt):
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    ax[1].plot(t[it],lwtip[it],'-',label="Edgewise");
    #ax[1].plot(t,wtip,'--',c='b',label="Edgewise");
ax[1].set_xlabel('Time [s]');
ax[1].set_ylabel(r'$x_{edge}$ [m]')
ax[1].grid(True);
ax[1].set_xlim([0,nt*dt])
    
for it in range(ndt):
    ax[2].plot(t[it],lptip[it]*180.0/np.pi,'-',label='dt = %1.3e'%(t[it][0]));
ax[2].set_ylabel(r'$\phi_z$ [deg]');
ax[2].set_xlabel('Time [s]');
ax[2].grid(True);
ax[2].set_xlim([0,nt*dt])

plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

# Put a legend below current axis
plt.subplots_adjust(wspace=0.3,bottom=0.3)
handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center',ncol=6,bbox_to_anchor=(0.5,0.1),fontsize=11)
#plt.xlabel(r'$f/f_{rot}$ [-]')
#plt.ylabel('PSD')
plt.tight_layout();
plt.savefig('figure/dtConvergence_dir%d.pdf'%(loaddir));

print(np.mean(lptip[it]*180.0/np.pi))

# Get frequency
signal = lvtip;
plt.figure()
for it in range(ndt):
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    fft  = np.fft.fft(signal[it]);
    nt = signal[it].size;
    freq = np.fft.fftfreq(signal[it].size,d=dt/(2.0**it));
    plt.plot(freq[1:nt//2],np.abs(fft[1:nt//2]),label='dt = %1.3e'%(t[it][0]))
    plt.xlim([0,6]);
    plt.legend()
    plt.grid(True);
    iFmax = np.argmax(np.abs(fft[1:nt//2]))+1
    eigF = freq[iFmax];
    valF = np.abs(fft[iFmax]);
    plt.semilogy(eigF,valF,'.r');
    print('First eigen Frequency : %1.3f [Hz]\n'%(eigF));