import numpy  as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interp2d

from time import time 

from pyBeamDyn import PyBeamDyn
from scipy.signal import welch

plt.close('all')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Typewriter"],
    "font.size": 12})

#%% ---- Definition of the dynamic problem  -----
dt = 2e-2
nt = 100

#%% Import the loads from file 

#%% ---- Solve the Dynamic Problem with BeamDyn -- python version -----
nBeam      = 1
L          = 8.63655000e+01
Rhub       = 2.8
omega      = 7.5*9/(L+Rhub);
ncase      = [2,4,6,8,10] 
inputFile  = ["DTU10MW_BeamDyn_n048_order%d.dat"%(ix) for ix in ncase];
nx         = []
DynamicSolve=1
nt_loc     = 1
GlbPos     = np.array([0.0,0.0,Rhub],order='F')
t          = np.linspace(dt,nt*dt*nt_loc,nt)
nnx = len(inputFile)

load_dir = 5

if DynamicSolve==0:
    nt = 2

t          = [None]*nnx
lvtip      = [None]*nnx
lutip      = [None]*nnx
ldvtip     = [None]*nnx
lwtip      = [None]*nnx
lptip      = [None]*nnx
lqtip      = [None]*nnx
lrtip      = [None]*nnx

rtime      = [None]*nnx

for ix in range(nnx):

    t = np.linspace(dt,nt*dt*nt_loc,nt);

    lvtip[ix]      = np.zeros(nt)
    lutip[ix]      = np.zeros(nt)
    ldvtip[ix]     = np.zeros(nt)
    lwtip[ix]      = np.zeros(nt)
    lptip[ix]      = np.zeros(nt)
    lqtip[ix]      = np.zeros(nt)
    lrtip[ix]      = np.zeros(nt)
    
    pbd = PyBeamDyn()
    pbd.initBeamDyn(1,inputFile[ix],1,nt=nt_loc,dt=dt,omega=np.array([omega*1,0,0],order='F'),DynamicSolve=DynamicSolve,GlbPos=GlbPos)
    
    # Get the position of the distributed loads
    xLoads, xDisp = pbd.getPositions()
    
    # Set the loads
    loads = np.zeros((6,pbd.nxLoads));
    loads[load_dir,:] = 4000;
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
        lvtip[ix][i]  =  u[0,-1];
        lwtip[ix][i]  =  u[1,-1];
        lutip[ix][i]  =  u[2,-1];
        ldvtip[ix][i] = du[0,-1];
        lptip[ix][i]  =  u[5,-1];
        lqtip[ix][i]  =  u[3,-1];
        lrtip[ix][i]  =  u[4,-1];
    end = time()
    rtime[ix] = (end-start)/nt;
    print("Done in %1.4e seconds -> %1.4e [s/iteration]"%(end-start,(end-start)/nt))

    

    # Free the data
    pbd.freeBeamDyn(1,1)

#%%


fig, ax = plt.subplots(1,3,figsize=(12,3.5))
fig.add_subplot(111, frameon=False)

for ix in range(nnx):
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    ax[0].plot(t,lvtip[ix],'-',label='order = %d'%(ncase[ix]));
    #ax[0].plot(t,wtip,'--',c='b',label="Edgewise");
ax[0].set_xlabel('Time [s]');
ax[0].set_ylabel(r'$x_{flap}$ [m]')
ax[0].set_xlim([0,nt*dt])
ax[0].grid(True);
#ax[0].legend()

for ix in range(nnx):
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
        ax[1].plot(t,lwtip[ix],'-',label='order = %d'%(ncase[ix]));
    #ax[1].plot(t,wtip,'--',c='b',label="Edgewise");
ax[1].set_xlabel('Time [s]');
ax[1].set_ylabel(r'$x_{edge}$ [m]')
ax[1].set_xlim([0,nt*dt])
ax[1].grid(True);
#ax[1].legend()
    
for ix in range(nnx):
    ax[2].plot(t,lptip[ix]*np.rad2deg(1),'-',label='order = %d'%(ncase[ix]));
ax[2].set_xlabel('Time [s]');
ax[2].set_ylabel(r'$\phi_z$ [deg]');
ax[2].set_xlim([0,nt*dt])
ax[2].grid(True);
#ax[2].legend()

plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

# Put a legend below current axis
plt.subplots_adjust(wspace=0.3,bottom=0.3)
handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center',ncol=6,bbox_to_anchor=(0.5,0.1),fontsize=11)
#plt.xlabel(r'$f/f_{rot}$ [-]')
#plt.ylabel('PSD')
plt.tight_layout();
plt.savefig('figure/spatialConvergence_dir%d.pdf'%(load_dir));

# Get frequency - flap
signal = lvtip.copy();
plt.figure()
for ix in range(nnx):
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    fft  = np.fft.fft(signal[ix]);
    nt = signal[ix].size;
    freq = np.fft.fftfreq(signal[ix].size,d=dt);
    plt.plot(freq[1:nt//2],np.abs(fft[1:nt//2]),label='order = %d'%(ncase[ix]));
    plt.xlim([0,10]);
    plt.legend()
    plt.grid(True);
    iFmax = np.argmax(np.abs(fft[1:nt//2]))+1
    eigF = freq[iFmax];
    valF = np.abs(fft[iFmax]);
    plt.semilogy(eigF,valF,'.r');
    print('First eigen Frequency Flap : %1.3f [Hz]\n'%(eigF));

# Get frequency - flap
# signal = lwtip.copy();
# plt.figure()
# for ix in range(nnx):
#     #plt.plot(t,utip,'--',c='r',label="Spanwise");
#     fft  = np.fft.fft(signal[ix]);
#     nt = signal[ix].size;
#     freq = np.fft.fftfreq(signal[ix].size,d=dt);
#     plt.plot(freq[1:nt//2],np.abs(fft[1:nt//2]),label='order = %d'%(ncase[ix]));
#     plt.xlim([0,6]);
#     plt.legend()
#     plt.grid(True);
#     iFmax = np.argmax(np.abs(fft[1:nt//2]))+1
#     eigF = freq[iFmax];
#     valF = np.abs(fft[iFmax]);
#     plt.semilogy(eigF,valF,'.r');
#     print('First eigen Frequency Edge : %1.3f [Hz]\n'%(eigF));

# Get frequency - torsion
signal = lvtip.copy();
plt.figure();
fs = 1/dt;
for ix in range(nnx):
    #plt.plot(t,utip,'--',c='r',label="Spanwise");
    #fft  = np.fft.fft(signal[ix]);
    #nt = signal[ix].size;
    #freq = np.fft.fftfreq(signal[ix].size,d=dt);
    freq,Px = welch(signal[ix],fs=fs,nperseg=signal[ix].size/2);
    plt.semilogy(freq,np.abs(Px),label='order = %d'%(ncase[ix]))
    #plt.plot(freq[1:nt//2],np.abs(fft[1:nt//2]),label='nx = %d'%(nx[ix]));
    plt.xticks(np.arange(6)*2)
    plt.xlim([0,10]);
    plt.legend()
    plt.grid(True);
    #iFmax = np.argmax(np.abs(fft[1:nt//2]))+1
    #eigF = freq[iFmax];
    #valF = np.abs(fft[iFmax]);
    #plt.semilogy(eigF,valF,'.r');
    #print('First eigen Frequency Torsion: %1.3f [Hz]\n'%(eigF));

#%%
print('Runtime (s/iter) : ',rtime);
plt.figure(figsize=(4,3));
plt.plot(ncase,rtime,'.-');
plt.grid(True);
plt.xlabel('N nodes');
plt.ylabel('Time per iteration [s]');
plt.xlim([12,64]);
plt.tight_layout()
plt.savefig('figure/spatialConv_runtime_dir%d.pdf'%(load_dir));