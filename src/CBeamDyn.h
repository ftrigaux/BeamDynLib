
/*
 *       ___        ___                         ___             
 *      / __\      / __\ ___  __ _ _ __ ___    /   \_   _ _ __  
 *     / /  _____ /__\/// _ \/ _` | '_ ` _ \  / /\ / | | | '_ \ 
 *    / /__|_____/ \/  \  __/ (_| | | | | | |/ /_//| |_| | | | |
 *    \____/     \_____/\___|\__,_|_| |_| |_/___,'  \__, |_| |_|
 *                                                  |___/       
 *        A C-wrapper around the BeamDyn Code by NREL
 *                       F.  Trigaux
 *                        July 2022
 */

#ifndef _CBEAMDYN_H_
#define _CBEAMDYN_H_

// Macro definition for Fortran array wrapping
#ifndef ALLOCATE1
    extern int i_,j_,k_;

    #define ALLOCATE1(ptr,type,n1) ptr = (type*)malloc((n1)*sizeof(type)) 
    #define ALLOCATE2(ptr,type,n1,n2) ptr = (type **)malloc((n1) * sizeof(type *)); \
                                    for (i_=0; i_<(n1); i_++){\
                                        ptr[i_] = (type *)malloc((n2) * sizeof(type));}
    #define ALLOCATE3(ptr,type,n1,n2,n3) ptr = (type ***)malloc((n1) * sizeof(type **)); \
                                    for (i_=0; i_<(n1); i_++){\
                                        ptr[i_] = (type **)malloc((n2) * sizeof(type *));\
                                        for (j_=0; j_<(n2); j_++){\
                                            ptr[i_][j_] = (type *)malloc((n3) * sizeof(type));}}

    #define DEALLOCATE1(ptr) free(ptr);
    #define DEALLOCATE2(ptr,n1) for (i_=0; i_<(n1); i_++){\
                                free(ptr[i_]);}\
                                free(ptr);
    #define DEALLOCATE3(ptr,n1,n2) for (i_=0; i_<(n1); i_++){\
                            for (j_=0; j_<(n2); j_++){\
                            free(ptr[i_][j_]);}\
                            free(ptr[i_]);}\
                            free(ptr);

    #define RAVEL2(obj,VEC,nx,ny,type) type*  rav_##VEC = (type*) malloc(sizeof(type)*(nx)*(ny)); \
    for (i_=0;i_<(nx);i_++){\
    for (j_=0;j_<(ny);j_++){\
    rav_##VEC[j_*(nx)+i_] = obj->VEC[i_][j_];\
    }}
    #define RAVEL3(obj,VEC,nx,ny,nz,type) type*  rav_##VEC = (type*) malloc(sizeof(type)*(nx)*(ny)*(nz)); \
    for (i_=0;i_<(nx);i_++){\
    for (j_=0;j_<(ny);j_++){\
    for (k_=0;k_<(nz);k_++){\
    rav_##VEC[k_*(nx)*(ny)+j_*(nx)+i_] = obj->VEC[i_][j_][k_];\
    }}}
    #define UNRAVEL2(obj,VEC,nx,ny,type)for (i_=0;i_<(nx);i_++){\
    for (j_=0;j_<(ny);j_++){\
    obj->VEC[i_][j_] = rav_##VEC[j_*(nx)+i_];\
    }}\
    free(rav_##VEC)
    #define UNRAVEL3(obj,VEC,nx,ny,nz,type) for (i_=0;i_<(nx);i_++){\
    for (j_=0;j_<(ny);j_++){\
    for (k_=0;k_<(nz);k_++){\
    obj->VEC[i_][j_][k_] = rav_##VEC[k_*(nx)*(ny)+j_*(nx)+i_];\
    }}}\
    free(rav_##VEC) 
#endif

#define FILLZERO2(ptr,n1,n2) for (i_=0; i_<(n1); i_++){\
                                    for (j_=0; j_<(n2); j_++){\
                                        ptr[i_][j_] = 0.0;\
                                    }}

typedef struct
{
    int     nBeam;                // Number of beams defined on this proc
    int     idx;                  // Index of this beam
    int     DynamicSolve;         // Dynmamic=1 or static=0
    char    *inputFile;           // Name of the primary input file containing the beam data

    int    GlbRotBladeT0;               // 1=The blade orientation at t=0 is the global frame
    double GlbPos[3];                   // Vector to the blade root, typically (0,0,Rtip)
    double GlbRot[3][3], RootOri[3][3]; // Global frame orientation and Blade root orientation

    int     nt;                   // Number of time step per solve
    double  dt;                   // Time increment [s]
    double  t;                    // Current structural time [s]

    int     nxL, nxD;             // Number of points for Load input and Displacement output
    double  **xLoads, **xDisp;    // Undeformed node position for loads and displacement [m]
    double  **loads;              // Load input (at xLoad) [N,Nm]
    double  **u, **du;            // Displacement and structural velocity output (at xDisp) [m]

    double  omega[3], domega[3], gravity[3];  // Angular velocity and acceleration [rad/s]
     
} BD_Data;

void BD_initBeamDyn(BD_Data *bd);

void BD_refresh(BD_Data *bd);

void BD_getPositions(BD_Data *bd);

void BD_freeBeamDyn(BD_Data *bd, int deallocAll);

void BD_setLoads(BD_Data *bd);

void BD_solve(BD_Data *bd);

void BD_getDisplacement(BD_Data *bd);

void BD_setBC(BD_Data *bd);

void BD_writeSolToBin(BD_Data *bd, char* fileName);

// Util functions
void BD_getRotationMatrix(double Rot[3][3], double c[3]);


// Fortran functions
void f_initBeamDyn(int,char*,int,double*,int*,double*,int*,double[3],double[3],double[3],double[3],int*,double[3][3],double[3][3],int*,int*);

void f_refresh(int,double*,int*,double*,double[3],double[3]);

void f_getPositions(double**,double**);

void f_freeBeamDyn(int,int);

void f_setLoads(double **, int);

void f_solve(int);

void f_getDisplacement(double**,double**,double**,int);

void f_setBC(int);

#endif