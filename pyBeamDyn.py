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

plt.close('all')

# import the shared library
bd = ct.CDLL('/Users/ftrigaux/Documents/Beams/BeamDynLib/libBeamDyn.so') 

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

# Specify args and result type for each function
bd.f_initBeamDyn.argtypes = [ct.c_int,ct.c_char_p,ct.c_int,
                             ct.POINTER(ct.c_double),ct.POINTER(ct.c_int),ct.POINTER(ct.c_double),
                             ct.POINTER(ct.c_int),
                             ND_ARRAY_3,ND_ARRAY_3,ND_ARRAY_3,
                             ND_ARRAY_3,ct.POINTER(ct.c_int),ND_ARRAY_3x3,ND_ARRAY_3x3,
                             ct.POINTER(ct.c_int), ct.POINTER(ct.c_int)]
bd.f_initBeamDyn.restype  = None

bd.f_refresh.argtypes = [ct.c_int, ct.POINTER(ct.c_double),ct.POINTER(ct.c_int),ct.POINTER(ct.c_double)]
bd.f_refresh.restype  = None

bd.f_getPositions.argtypes = [ND_POINTER_2,ND_POINTER_2]
bd.f_getPositions.restype  = None

bd.f_freeBeamDyn.argtypes = [ct.c_int,ct.c_int]
bd.f_freeBeamDyn.restype  = None

bd.f_setLoads.argtypes = [ND_POINTER_2, ct.c_int]
bd.f_setLoads.restype  = None

bd.f_solve.argtypes = [ct.c_int]
bd.f_solve.restype  = None

bd.f_getDisplacement.argtypes = [ND_POINTER_2,ND_POINTER_2,ND_POINTER_2, ct.c_int]
bd.f_getDisplacement.restype  = None

bd.f_setBC.argtypes = [ct.c_int,ND_ARRAY_3,ND_ARRAY_3]
bd.f_setBC.restype  = None

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
    
def py_initBeamDyn(nBeam,inputFile,idxBeam,dt=None,nt=None,t=None,DynamicSolve=None,omega=None,domega=None,gravity=None,GlbPos=None, GlbRotBladeT0=None, GlbRot=None, RootOri=None):
    nxLoads = ct.c_int(0); nxDisp = ct.c_int(0);
    bd.f_initBeamDyn(nBeam,inputFile.encode('utf-8'),idxBeam,opt(dt,ct.c_double),opt(nt,ct.c_int),opt(t,ct.c_double),opt(DynamicSolve,ct.c_int),opt(omega,ct.c_double),opt(domega,ct.c_double),opt(gravity,ct.c_double),opt(GlbPos,ct.c_double),opt(GlbRotBladeT0,ct.c_int),opt(GlbRot,ct.c_double),opt(RootOri,ct.c_double),ct.byref(nxLoads),ct.byref(nxDisp))
    return nxLoads.value, nxDisp.value

def py_refresh(idxBeam,dt=None,nt=None,t=None,omega=None,domega=None):
    nxLoads = ct.c_int(0); nxDisp = ct.c_int(0);
    bd.f_refresh(idxBeam,opt(dt,ct.c_double),opt(nt,ct.c_int),opt(t,ct.c_double));
    return nxLoads.value, nxDisp.value

def getRotationMatrix(c):
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
    # Initialize BeamDyn
    nBeam      = 1
    idxBeam    = 1
    inputFile  = "/Users/ftrigaux/Documents/Beams/BeamDynLib/run/nrel5mw_dynamic/bd_primary_nrel_5mw_dynamic.inp"
    nt_loc     = 1;
    nt         = 30;
    dt_loc     = 5e-2;
    
    DynamicSolve = 1
    
    omega=np.array([-1,0,0],order='F')
    GlbPos = np.array([0,0,1.5],order='F')
    
    theta = 0*np.pi/180;
    grav  = np.array([0,0,0],order='F');
    GlbRotBladeT0 = 1;
    RootOri = np.array([[ 1, 0           , 0           ],
                         [ 0, np.cos(theta),-np.sin(theta)],
                         [ 0, np.sin(theta), np.cos(theta)]]);
    nxL, nxD = py_initBeamDyn(nBeam, inputFile, idxBeam, nt=nt_loc, dt=dt_loc, GlbPos=GlbPos, omega=omega, DynamicSolve=DynamicSolve, gravity=grav, RootOri=RootOri, GlbRotBladeT0=GlbRotBladeT0)

    # Get the position of the distributed loads
    xLoads = np.zeros((3,nxL),order='F')
    xDisp  = np.zeros((3,nxD),order='F')
    bd.f_getPositions(NP_F_2D(xLoads),NP_F_2D(xDisp))

    # Set the loads
    loads = np.zeros((6,nxL),order='F');
    loads[0,:] = 1e4;
    bd.f_setLoads(NP_F_2D(loads),idxBeam)
    
    # Preallocate variables to extract displacement
    t = np.linspace(dt_loc*nt_loc,(nt+1)*dt_loc*nt_loc,nt)
    
    x  = np.zeros((3,nxD),order='F');
    u  = np.zeros((6,nxD),order='F');
    du = np.zeros((6,nxD),order='F'); 

    vtip      = np.zeros((3,nt))
    vvtip     = np.zeros((3,nt))

    # Solve the dynamic
    plt.figure(1)
    for i in range(nt):
        bd.f_setBC(idxBeam,NP_F_3(omega),NP_F_3(np.zeros(3,order='F')))
        bd.f_solve(idxBeam)
    
        # Preallocate variables to extract displacement
        bd.f_getDisplacement(NP_F_2D(x),NP_F_2D(u),NP_F_2D(du),idxBeam)
        
        plt.plot(x[2,:],u[0,:],'.-r')
        plt.plot(x[2,:],u[1,:],'.-g')
        plt.plot(x[2,:],u[2,:],'.-b')
        
        vtip[:,i]  =  u[:3,-1];
        vvtip[:,i] = du[:3,-1];
    
    # Get Rotation matrix 
    c = u[3:,-1]
    Rot = getRotationMatrix(c);
    print(Rot);
    vdir = np.array([[0,0,0],[0,0.0,0.1]]).T * x[2,-1];
    vdirc = Rot @ vdir
    plt.plot(x[2,-1]+vdirc[2,:],u[0,-1]+vdirc[0,:],'k--');
    
    # Free the data
    bd.f_freeBeamDyn(idxBeam,idxBeam==nBeam)
    
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
