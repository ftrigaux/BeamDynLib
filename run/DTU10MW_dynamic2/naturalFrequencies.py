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

# Use to suppress damping in input file
def replaceInFile(fname,str_old,str_new):
    # Read in the file
    with open(fname, 'r') as file :
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace(str_old, str_new)

    # Write the file out again
    with open(fname, 'w') as file:
        file.write(filedata)

#%% ---- Definition of the dynamic problem  -----
dt = 2.5e-3
nt = 20000

#%% Import the loads from file 

#%% ---- Solve the Dynamic Problem with BeamDyn -- python version -----
nBeam      = 1
L          = 8.63655000e+01
Rhub       = 2.8
omega      = 7.5*9/(L+Rhub);
inputFile  = "DTU10MW_BeamDyn_n048.dat"
DynamicSolve=1
nt_loc     = 1
GlbPos     = np.array([0.0,0.0,Rhub],order='F')
t          = np.linspace(dt,nt*dt*nt_loc,nt)

if DynamicSolve==0:
    nt = 2

loaddir = [5];
nl = len(loaddir)

# Remove damping
replaceInFile("DTU10MW_BeamDyn_n048_blade.dat","1  damp_type","0  damp_type")

t          = [None]*nl
lvtip      = [None]*nl
lutip      = [None]*nl
ldvtip     = [None]*nl
lwtip      = [None]*nl
lptip      = [None]*nl
lqtip      = [None]*nl
lrtip      = [None]*nl

rtime      = [None]*nl

for ix in range(len(loaddir)):

    t = np.linspace(dt,nt*dt*nt_loc,nt);

    lvtip[ix]      = np.zeros(nt)
    lutip[ix]      = np.zeros(nt)
    ldvtip[ix]     = np.zeros(nt)
    lwtip[ix]      = np.zeros(nt)
    lptip[ix]      = np.zeros(nt)
    lqtip[ix]      = np.zeros(nt)
    lrtip[ix]      = np.zeros(nt)
    
    pbd = PyBeamDyn()
    pbd.initBeamDyn(1,inputFile,1,nt=nt_loc,dt=dt,omega=np.array([omega*0,0,0],order='F'),DynamicSolve=DynamicSolve,GlbPos=GlbPos)
    
    # Get the position of the distributed loads
    xLoads, xDisp = pbd.getPositions()
    
    # Set the loads
    loads = np.zeros((6,pbd.nxLoads));
    loads[loaddir[ix],:] = 4000;
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


# Reset damping
replaceInFile("DTU10MW_BeamDyn_n048_blade.dat","0  damp_type","1  damp_type")
#%%
from scipy.signal import argrelextrema
# Get frequency - flap
signal = lvtip.copy();
ix = 0
plt.figure(figsize=(4,3))
#plt.plot(t,utip,'--',c='r',label="Spanwise");
fft  = np.fft.fft(signal[ix]);
nt = signal[ix].size;
freq = np.fft.fftfreq(signal[ix].size,d=dt);
plt.plot(freq[1:nt//2],np.abs(fft[1:nt//2]));
plt.xlim([0,10]);
plt.xlabel('Freq [Hz]')
plt.ylabel(r'FFT$(x_{tip,f})$')
plt.grid(True);
iFmax = np.argmax(np.abs(fft[1:nt//2]))+1
eigF = freq[iFmax];
valF = np.abs(fft[iFmax]);
plt.semilogy(eigF,valF,'.r');
print('First eigen Frequency Flap : %1.3f [Hz]\n'%(eigF));
plt.savefig('figure/natFreq_fft_flap.pdf')

#Get frequency - edge
signal = lwtip.copy();
ix = min(1,nl-1)
plt.figure(figsize=(4,3))
#plt.plot(t,utip,'--',c='r',label="Spanwise");
fft  = np.fft.fft(signal[ix]);
nt = signal[ix].size;
freq = np.fft.fftfreq(signal[ix].size,d=dt);
plt.plot(freq[1:nt//2],np.abs(fft[1:nt//2]));
plt.xlim([0,10]);
plt.xlabel('Freq [Hz]')
plt.grid(True);
iFmax = np.argmax(np.abs(fft[1:nt//2]))+1
eigF = freq[iFmax];
valF = np.abs(fft[iFmax]);
plt.ylabel(r'FFT$(x_{tip,e})$')
plt.semilogy(eigF,valF,'.r');
print('First eigen Frequency Edge : %1.3f [Hz]\n'%(eigF));
plt.savefig('figure/natFreq_fft_edge.pdf')


#Get frequency - Torsion
signal = lwtip.copy();
ix = min(2,nl-1)
plt.figure(figsize=(4,3))
#plt.plot(t,utip,'--',c='r',label="Spanwise");
fft  = np.fft.fft(signal[ix]);
nt = signal[ix].size;
freq = np.fft.fftfreq(signal[ix].size,d=dt);
plt.plot(freq[1:nt//2],np.abs(fft[1:nt//2]));
plt.ylabel(r'FFT$(x_{tip,\phi})$')
plt.xlabel('Freq [Hz]')
plt.xlim([0,10]);
plt.grid(True);
iFmax = np.argmax(np.abs(fft[1:nt//2]))+1
eigF = freq[argrelextrema(np.abs(fft),np.greater)];
valF = np.abs(fft[argrelextrema(np.abs(fft),np.greater)]);
plt.semilogy(eigF,valF,'.r');
print(eigF);
plt.savefig('figure/natFreq_fft_tors.pdf')

plt.figure();
plt.plot(t,lwtip[ix]*np.rad2deg(1));
plt.xlabel('T [s]');
plt.ylabel('phi_tip');