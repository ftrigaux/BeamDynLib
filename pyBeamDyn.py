import numpy  as np
import ctypes as ct
import matplotlib.pyplot as plt

plt.close('all')

# import the shared library
bd = ct.CDLL('/Users/ftrigaux/Documents/Beams/BeamDynLib/libBeamDyn.so') 

ND_POINTER_1 = ct.POINTER(ct.c_double);
ND_POINTER_2 = ct.POINTER(ct.POINTER(ct.c_double));
ND_POINTER_3 = ct.POINTER(ct.POINTER(ct.POINTER(ct.c_double)));

# Change a ND numpy array to a double pointer for fortran
def NP_F_1D(array):
    return ct.byref(array.ctypes.data_as(ct.POINTER(ct.c_double)))
def NP_F_2D(array):
    return ct.byref(array.ctypes.data_as(ct.POINTER(ct.c_double)))


# Specify args and result type for each function
bd.f_initBeamDyn.argtypes = [ct.c_int,ct.c_char_p,
                           ct.POINTER(ct.c_double),ct.POINTER(ct.c_int),ct.POINTER(ct.c_double),
                           ct.POINTER(ct.c_int),
                           ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),
                           ct.POINTER(ct.c_int), ct.POINTER(ct.c_int)]
bd.f_initBeamDyn.restype  = None

bd.f_getPositions.argtypes = [ND_POINTER_2,ND_POINTER_2]
bd.f_getPositions.restype  = None

bd.f_freeBeamDyn.argtypes = []
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
        else: 
            return ct.byref(ftype(dt[0]))
    else:
        return None
    
def py_InitBeamDyn(nBeam,inputFile,dt=None,nt=None,t=None,DynamicSolve=None,omega=None,domega=None,gravity=None):
    nxLoads = ct.c_int(0); nxDisp = ct.c_int(0);
    bd.f_initBeamDyn(1,inputFile.encode('utf-8'),opt(dt,ct.c_double),opt(nt,ct.c_int),opt(t,ct.c_double),opt(DynamicSolve,ct.c_int),opt(omega,ct.c_double),opt(domega,ct.c_double),opt(gravity,ct.c_double),ct.byref(nxLoads),ct.byref(nxDisp))
    return nxLoads.value, nxDisp.value

# Main function
if __name__ == "__main__":
    # Initialize BeamDyn
    nBeam      = 1
    inputFile  = "/Users/ftrigaux/Documents/Beams/BeamDynLib/run/nrel5mw_dynamic/bd_primary_nrel_5mw_dynamic.inp"
    nt_loc = 1;
    dt_loc = 1e-2;
    nxL, nxD = py_InitBeamDyn(1,inputFile,nt=nt_loc,dt=dt_loc,omega=-np.array([0.0,0,0],order='F'),DynamicSolve=1)
    
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

    vtip      = np.zeros(nt)
    vtip_test = np.zeros(nt+1)
    vvtip     = np.zeros(nt)

    # Solve the dynamic
    plt.figure(1)
    for i in range(nt):
        bd.f_solve(1)
    
        # Preallocate variables to extract displacement
        bd.f_getDisplacement(NP_F_2D(x),NP_F_2D(u),NP_F_2D(du),1)
        
        plt.plot(x[:,2],u[:,0],'r')
        plt.plot(x[:,2],u[:,1],'g')
        plt.plot(x[:,2],u[:,2],'b')
        
        vtip[i]  =  u[-1,0];
        vvtip[i] = du[-1,0];
        vtip_test[i+1] = du[-1,0]*dt_loc + vtip_test[i] 
    
    # Free the data
    bd.f_freeBeamDyn()
    
    plt.figure(2)
    plt.plot(t,vtip);
    plt.plot(t,vvtip);
    
    plt.plot(t,vtip_test[1:]);