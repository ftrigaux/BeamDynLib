# 
#        ___        ___                         ___             
#       / _ \      / __\ ___  __ _ _ __ ___    /   \_   _ _ __  
#      / /_)/____ /__\/// _ \/ _` | '_ ` _ \  / /\ / | | | '_ \ 
#     / ___/_____/ \/  \  __/ (_| | | | | | |/ /_//| |_| | | | |
#     \/         \_____/\___|\__,_|_| |_| |_/___,'  \__, |_| |_|
#                                                   |___/       
#        A Python wrapper around the BeamDyn Code by NREL
#                       F.  Trigaux
#                        July 2022
# 
# 

import numpy  as np
import ctypes as ct
import matplotlib.pyplot as plt

# Define some custom types
ND_POINTER_1 = ct.POINTER(ct.c_double);
ND_POINTER_2 = ct.POINTER(ct.POINTER(ct.c_double));
ND_POINTER_3 = ct.POINTER(ct.POINTER(ct.POINTER(ct.c_double)));

ND_ARRAY_3     = ct.POINTER(ct.c_double*3)
ND_ARRAY_3x3   = ct.POINTER((ct.c_double*3)*3)

# Change a ND numpy array to a double pointer for fortran TYPE(C_PTR)
def NP_F_1D(array):
    return ct.byref(array.ctypes.data_as(ct.POINTER(ct.c_double)))
def NP_F_2D(array):
    return ct.byref(array.ctypes.data_as(ct.POINTER(ct.c_double)))

# Change a ND numpy array to a fortran array x(n)
def NP_F_3(array):
    A = ct.c_double * 3
    a = A()
    for i in range(3):
            a[i] = array[i]
    return a

def NP_F_3x3(array):
    A = (ct.c_double * 3) * 3
    a = A()
    for i in range(3):
        for j in range(3):
            a[i][j] = array[i][j]
    return a

# ---- Wrap the functions in Python ----
def opt(dt,ftype):
    if dt is not None:
        if np.isscalar(dt):
            return ct.byref(ftype(dt))
        elif dt.ndim==1:
            return NP_F_3(dt)
        else:
            return NP_F_3x3(dt)
    else:
        return None

# Specify args and result type for each function

class PyBeamDyn(ct.CDLL):

    def __init__(self,library_path='/Users/ftrigaux/Documents/Beams/BeamDynLib/libBeamDyn.so'):
        super().__init__(library_path)
        self.f_initBeamDyn.argtypes = [ct.c_int,ct.c_char_p,ct.c_int,
                                    ct.POINTER(ct.c_double),ct.POINTER(ct.c_int),ct.POINTER(ct.c_double),
                                    ct.POINTER(ct.c_int),
                                    ND_ARRAY_3,ND_ARRAY_3,ND_ARRAY_3,
                                    ND_ARRAY_3,ct.POINTER(ct.c_int),ND_ARRAY_3x3,
                                    ct.POINTER(ct.c_int), ct.POINTER(ct.c_int)]
        self.f_initBeamDyn.restype  = None

        self.f_refresh.argtypes = [ct.c_int, ct.POINTER(ct.c_double),ct.POINTER(ct.c_int),ct.POINTER(ct.c_double)]
        self.f_refresh.restype  = None

        self.f_getPositions.argtypes = [ND_POINTER_2,ND_POINTER_2]
        self.f_getPositions.restype  = None

        self.f_freeBeamDyn.argtypes = [ct.c_int,ct.c_int]
        self.f_freeBeamDyn.restype  = None

        self.f_setLoads.argtypes = [ND_POINTER_2, ct.c_int]
        self.f_setLoads.restype  = None

        self.f_solve.argtypes = [ct.c_int]
        self.f_solve.restype  = None

        self.f_getDisplacement.argtypes = [ND_POINTER_2,ND_POINTER_2,ND_POINTER_2, ct.c_int]
        self.f_getDisplacement.restype  = None

        self.f_getReactionForce.argtypes = [ND_POINTER_2,ND_POINTER_2, ct.c_int]
        self.f_getReactionForce.restype  = None

        self.f_setBC.argtypes = [ct.c_int,ND_ARRAY_3,ND_ARRAY_3,ct.c_int,ct.c_double,ct.c_double]
        self.f_setBC.restype  = None

        self.nxLoads = -1;
        self.nxDisp  = -1;
    
    def initBeamDyn(self,nBeam,inputFile,idxBeam,dt=None,nt=None,t=None,DynamicSolve=None,omega=None,domega=None,gravity=None,GlbPos=None, GlbRotBladeT0=None, GlbRot=None, RootOri=None):
        if (GlbRot is not None):
            print("!!! Warning : GlbRot is no more an input of BeamDyn !!! Use RootOri instead with GlbRotBladeT0=1 to test the global rotation matrix")
        nxLoads = ct.c_int(0); nxDisp = ct.c_int(0);
        self.f_initBeamDyn(nBeam,inputFile.encode('utf-8'),idxBeam,opt(dt,ct.c_double),opt(nt,ct.c_int),opt(t,ct.c_double),opt(DynamicSolve,ct.c_int),opt(omega,ct.c_double),opt(domega,ct.c_double),opt(gravity,ct.c_double),opt(GlbPos,ct.c_double),opt(GlbRotBladeT0,ct.c_int),opt(RootOri,ct.c_double),ct.byref(nxLoads),ct.byref(nxDisp))
        self.nxLoads = nxLoads.value;
        self.nxDisp  = nxDisp.value;
        return nxLoads.value, nxDisp.value

    def refresh(self,idxBeam,dt=None,nt=None,t=None):
        self.f_refresh(idxBeam,opt(dt,ct.c_double),opt(nt,ct.c_int),opt(t,ct.c_double));

    def getPositions(self):
        xLoads = np.zeros((3,self.nxLoads),order='F')
        xDisp  = np.zeros((3,self.nxDisp),order='F')
        self.f_getPositions(NP_F_2D(np.asfortranarray(xLoads)),NP_F_2D(np.asfortranarray(xDisp)))
        return xLoads, xDisp

    def setLoads(self,loads,idxBeam):
        self.f_setLoads(NP_F_2D(np.asfortranarray(loads)),idxBeam)

    def setBC(self,idxBeam,omega,dOmega,forceTheta=0,theta_rot=0.0,pitch_rad=0.0):
        self.f_setBC(idxBeam,NP_F_3(omega),NP_F_3(dOmega),forceTheta,theta_rot,pitch_rad)

    def getDisplacement(self,idxBeam):
        x  = np.zeros((3,self.nxDisp),order='F');
        u  = np.zeros((6,self.nxDisp),order='F');
        du = np.zeros((6,self.nxDisp),order='F'); 
        self.f_getDisplacement(NP_F_2D(x),NP_F_2D(u),NP_F_2D(du),idxBeam)
        return x,u,du
    
    def getReactionForce(self,idxBeam):
        x  = np.zeros((3,self.nxDisp),order='F');
        F  = np.zeros((6,self.nxDisp),order='F');
        self.f_getReactionForce(NP_F_2D(x),NP_F_2D(F),idxBeam)
        return x,F

    def solve(self,idxBeam):
        self.f_solve(idxBeam)

    def freeBeamDyn(self,idxBeam,nBeam):
        self.f_freeBeamDyn(idxBeam,idxBeam==nBeam)

    def getRotationMatrix(self,c):
        c1  = c[0]/4.0;
        c2  = c[1]/4.0;
        c3  = c[2]/4.0;
        c0  = 0.5 * (1.0-c1*c1-c2*c2-c3*c3);     # 1/4 the value of the the AIAA paper (after plugging in c1, c2, c3 conversions)


        tr0 = 1.0 - c0;                          # This is 1/4 the value of the AIAA paper, after converting c0.
        tr0 = 2.0/(tr0*tr0);                     # This is 32x the equivalent term from the AIAA paper.   This is well behaved and won't go to zero.

        # The following terms can be shown to match the transpose of the DCM given in the AIAA paper.
        Rot = np.zeros((3,3));
        Rot[0,0] = tr0*(c1*c1 + c0*c0) - 1.0;
        Rot[1,0] = tr0*(c1*c2 + c0*c3);
        Rot[2,0] = tr0*(c1*c3 - c0*c2);

        Rot[0,1] = tr0*(c1*c2 - c0*c3);
        Rot[1,1] = tr0*(c2*c2 + c0*c0) - 1.0;
        Rot[2,1] = tr0*(c2*c3 + c0*c1);

        Rot[0,2] = tr0*(c1*c3 + c0*c2);
        Rot[1,2] = tr0*(c2*c3 - c0*c1);
        Rot[2,2] = tr0*(c3*c3 + c0*c0) - 1.0;
        
        return Rot;

# Main function
if __name__ == "__main__":
    plt.close('all')
    # Initialize BeamDyn
    nBeam      = 1
    idxBeam    = 1
    inputFile  = "/Users/ftrigaux/Documents/Beams/BeamDynLib/run/nrel5mw_dynamic/bd_primary_nrel_5mw_dynamic.inp"
    nt_loc     = 1;
    nt         = 30;
    dt_loc     = 5e-2;
    
    DynamicSolve = 1
    
    omega=np.array([-1,0,0])
    GlbPos = np.array([0,0,1.5])
    
    theta = 0*np.pi/180;
    grav  = np.array([0,0,0]);
    GlbRotBladeT0 = 1;
    RootOri = np.array([[ 1, 0           , 0           ],
                         [ 0, np.cos(theta),-np.sin(theta)],
                         [ 0, np.sin(theta), np.cos(theta)]]);
    
    pbd = PyBeamDyn();

    pbd.initBeamDyn(nBeam, inputFile, idxBeam, nt=nt_loc, dt=dt_loc, GlbPos=GlbPos, omega=omega, DynamicSolve=DynamicSolve, gravity=grav, RootOri=RootOri, GlbRotBladeT0=GlbRotBladeT0)

    # Get the position of the distributed loads
    xLoads, xDisp = pbd.getPositions()

    # Set the loads
    loads = np.zeros((6,pbd.nxLoads));
    loads[0,:] = 1e4;
    pbd.setLoads(loads,idxBeam)
    
    # Preallocate variables to extract displacement
    t = np.linspace(dt_loc*nt_loc,nt*dt_loc*nt_loc,nt)

    vtip      = np.zeros((3,nt))
    vvtip     = np.zeros((3,nt))

    # Solve the dynamic
    plt.figure(1)
    for i in range(nt):
        pbd.setBC(idxBeam,omega,np.zeros(3))
        pbd.solve(idxBeam)
    
        # Preallocate variables to extract displacement
        x,u,du = pbd.getDisplacement(idxBeam)
        
        plt.plot(x[2,:],u[0,:],'.-r')
        plt.plot(x[2,:],u[1,:],'.-g')
        plt.plot(x[2,:],u[2,:],'.-b')
        
        vtip[:,i]  =  u[:3,-1];
        vvtip[:,i] = du[:3,-1];
    
    # Get Rotation matrix 
    c = u[3:,-1]
    Rot = pbd.getRotationMatrix(c);
    print(Rot);
    vdir = np.array([[0,0,0],[0,0.0,0.1]]).T * x[2,-1];
    vdirc = Rot @ vdir
    plt.plot(x[2,-1]+vdirc[2,:],u[0,-1]+vdirc[0,:],'k--');
    
    # Free the data
    pbd.freeBeamDyn(idxBeam,idxBeam==nBeam)
    
    plt.figure(2)
    plt.plot(t,vtip[0,:],'r');
    plt.plot(t,vtip[1,:],'g');
    plt.plot(t,vtip[2,:],'b');
    plt.grid(True);


    # Estimation of the velocity based on the blade positions
    vtest = (vtip[:,1:] - vtip[:,:-1])/(np.outer(np.ones(3),t[1:]-t[:-1]));
    tmean = (t[1:]+t[:-1])/2.0
    
    plt.figure(3)
    plt.plot(t,vvtip[0,:],'r');
    plt.plot(tmean,vtest[0,:],'r--');
    plt.plot(t,vvtip[1,:],'g');
    plt.plot(tmean,vtest[1,:],'g--');
    plt.plot(t,vvtip[2,:],'b');
    plt.plot(tmean,vtest[2,:],'b--');
    plt.grid(True);
