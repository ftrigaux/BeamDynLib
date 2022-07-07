import numpy  as np
import ctypes as ct
import matplotlib.pyplot as plt

plt.close('all')

# import the shared library
bd = ct.CDLL('BeamDyn_lib.so') 

ND_POINTER_1 = ct.POINTER(ct.c_double);
ND_POINTER_2 = ct.POINTER(ct.POINTER(ct.c_double));
ND_POINTER_3 = ct.POINTER(ct.POINTER(ct.POINTER(ct.c_double)));

# Change a ND numpy array to a double pointer for fortran
def NP_F_1D(array):
    return ct.byref(array.ctypes.data_as(ct.POINTER(ct.c_double)))
def NP_F_2D(array):
    return ct.byref(array.ctypes.data_as(ct.POINTER(ct.c_double)))


# Specify args and result type for each function
bd.initBeamDyn.argtypes = [ct.c_int,ct.c_char_p,
                           ct.POINTER(ct.c_double),ct.POINTER(ct.c_int),ct.POINTER(ct.c_double),
                           ct.POINTER(ct.c_int),
                           ct.POINTER(ct.c_double),ct.POINTER(ct.c_double)]
bd.initBeamDyn.restype  = None

bd.freeBeamDyn.argtypes = []
bd.freeBeamDyn.restype  = None

bd.setLoads.argtypes = [ND_POINTER_2, ct.c_int]
bd.setLoads.restype  = None

bd.solve.argtypes = [ct.c_int]
bd.solve.restype  = None

bd.getDisplacement.argtypes = [ND_POINTER_2,ND_POINTER_2,ND_POINTER_2, ct.c_int]
bd.getDisplacement.restype  = None

# ---- Wrap the functions in Python ----
def opt(dt,ftype):
    if dt is not None:
        if np.isscalar(dt):
            return ct.byref(ftype(dt))
        else: 
            return ct.byref(ftype(dt[0]))
    else:
        return None
    
def pyInitBeamDyn(nBeam,inputFile,dt=None,nt=None,t=None,DynamicSolve=None,omega=None,domega=None):
    bd.initBeamDyn(1,inputFile.encode('utf-8'),opt(dt,ct.c_double),opt(nt,ct.c_int),opt(t,ct.c_double),opt(DynamicSolve,ct.c_int),opt(omega,ct.c_double),opt(domega,ct.c_double))

# Main function
if __name__ == "__main__":
    # Initialize BeamDyn
    nBeam      = 1
    inputFile  = "/Users/ftrigaux/Documents/Beams/BeamDynLib/run/nrel5mw_dynamic/bd_primary_nrel_5mw_dynamic.inp"
    pyInitBeamDyn(1,inputFile,nt=25,omega=np.array([0.0,0,0],order='F'))
    

    # Set the loads
    loads = np.zeros((49,6),order='F');
    loads[:,1] = 1e3;
    bd.setLoads(NP_F_2D(loads),1)
    
    # Solve the dynamic
    plt.figure(1)
    for i in range(20):
        bd.solve(1)
    
        # Preallocate variables to extract displacement
        x  = np.zeros((49,3),order='F');
        u  = np.zeros((49,6),order='F');
        du = np.zeros((49,6),order='F'); 
        bd.getDisplacement(NP_F_2D(x),NP_F_2D(u),NP_F_2D(du),1)
        
        plt.plot(x[:,2],u[:,1],'b')
        plt.plot(x[:,2],u[:,0],'r')
        plt.plot(x[:,2],u[:,3],'g')
    
    bd.freeBeamDyn()