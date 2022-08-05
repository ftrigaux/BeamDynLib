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

bd.f_refresh.argtypes = [ct.c_int, ct.POINTER(ct.c_double),ct.POINTER(ct.c_int),ct.POINTER(ct.c_double),ND_ARRAY_3,ND_ARRAY_3]
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
    bd.f_refresh(idxBeam,opt(dt,ct.c_double),opt(nt,ct.c_int),opt(t,ct.c_double),opt(omega,ct.c_double),opt(domega,ct.c_double));
    return nxLoads.value, nxDisp.value

# Main function
if __name__ == "__main__":
    # Initialize BeamDyn
    nBeam      = 1
    idxBeam    = 1
    inputFile  = "/Users/ftrigaux/Documents/Beams/BeamDynLib/run/nrel5mw_dynamic/bd_primary_nrel_5mw_dynamic.inp"
    nt_loc = 1;
    dt_loc = 1e-2;
    
    DynamicSolve = 1
    
    omega=np.array([-1,0,0],order='F')
    
    theta = 45*np.pi/180;
    grav  = np.array([0,0,-9.81],order='F');
    GlbRotBladeT0 = 1;
    RootOri = np.array([[ 1, 0           , 0           ],
                         [ 0, np.cos(theta),-np.sin(theta)],
                         [ 0, np.sin(theta), np.cos(theta)]]);
    nxL, nxD = py_initBeamDyn(nBeam, inputFile, idxBeam, nt=nt_loc, dt=dt_loc, omega=omega, DynamicSolve=DynamicSolve, gravity=grav, RootOri=RootOri, GlbRotBladeT0=GlbRotBladeT0)

    # Get the position of the distributed loads
    xLoads = np.zeros((nxL,3),order='F')
    xDisp  = np.zeros((nxD,3),order='F')
    bd.f_getPositions(NP_F_2D(xLoads),NP_F_2D(xDisp))

    # Set the loads
    loads = np.zeros((nxL,6),order='F');
    loads[:,0] = 1e3;
    bd.f_setLoads(NP_F_2D(loads),1)
    
    # Preallocate variables to extract displacement
    nt = 5;
    t = np.linspace(0,nt*dt_loc*nt_loc,nt)
    
    x  = np.zeros((nxD,3),order='F');
    u  = np.zeros((nxD,6),order='F');
    du = np.zeros((nxD,6),order='F'); 

    vtip      = np.zeros((3,nt))
    vvtip     = np.zeros((3,nt))

    # Solve the dynamic
    plt.figure(1)
    for i in range(nt):
        bd.f_solve(1)
    
        # Preallocate variables to extract displacement
        bd.f_getDisplacement(NP_F_2D(x),NP_F_2D(u),NP_F_2D(du),1)
        
        plt.plot(x[:,2],u[:,0],'.-r')
        plt.plot(x[:,2],u[:,1],'.-g')
        plt.plot(x[:,2],u[:,2],'.-b')
        
        vtip[:,i]  =  u[-1,:3];
        vvtip[:,i] = du[-1,:3];
    
    # Free the data
    bd.f_freeBeamDyn(1,1)
    
    plt.figure(2)
    plt.plot(t*180/np.pi,vtip[0,:],'r');
    plt.plot(t*180/np.pi,vtip[1,:],'g');
    plt.plot(t*180/np.pi,vtip[2,:],'b');
    plt.grid(True);
