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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "CBeamDyn.h"

int i_,j_,k_; // iterator for macros

void BD_initBeamDyn(BD_Data *bd)
{
    f_initBeamDyn(bd->nBeam,bd->inputFile,bd->idx,
                &bd->dt, &bd->nt, &bd->t,          
                &bd->DynamicSolve,       
                bd->omega, bd->domega, bd->gravity,      
                &bd->nxL, &bd->nxD);

    ALLOCATE2(bd->xLoads,double,bd->nxL,3);
    ALLOCATE2(bd->xDisp, double,bd->nxD,3);
    ALLOCATE2(bd->loads, double,bd->nxL,6);
    ALLOCATE2(bd->u,     double,bd->nxD,6);
    ALLOCATE2(bd->du,    double,bd->nxD,6);

    FILLZERO2(bd->xLoads,bd->nxL,3);
    FILLZERO2(bd->xDisp, bd->nxD,3);
    FILLZERO2(bd->loads, bd->nxL,6);
    FILLZERO2(bd->u,     bd->nxD,6);
    FILLZERO2(bd->du,    bd->nxD,6);

}

void BD_refresh(BD_Data *bd)
{
    f_refresh(bd->idx,
                &bd->dt, &bd->nt, &bd->t,
                bd->omega, bd->domega);

}

void BD_getPositions(BD_Data *bd)
{
    RAVEL2(bd,xLoads,bd->nxL,3,double);
    RAVEL2(bd,xDisp ,bd->nxL,3,double);
    f_getPositions(&rav_xLoads, &rav_xDisp);

    UNRAVEL2(bd,xLoads,bd->nxL,3,double);
    UNRAVEL2(bd,xDisp ,bd->nxL,3,double);

}

void BD_freeBeamDyn(BD_Data *bd, int deallocAll)
{
    f_freeBeamDyn(bd->idx,deallocAll);

    DEALLOCATE2(bd->xLoads,bd->nxL);
    DEALLOCATE2(bd->xDisp ,bd->nxD);
    DEALLOCATE2(bd->loads ,bd->nxL);
    DEALLOCATE2(bd->u     ,bd->nxD);
    DEALLOCATE2(bd->du    ,bd->nxD);
}


void BD_setLoads(BD_Data *bd)
{
    int i,j;
    RAVEL2(bd,loads,bd->nxL,6,double);
    f_setLoads(&rav_loads,bd->idx);
    UNRAVEL2(bd,loads,bd->nxL,6,double);
}

void BD_solve(BD_Data *bd)
{
    f_solve(bd->idx);
}

void BD_getDisplacement(BD_Data *bd)
{
    RAVEL2(bd,xDisp ,bd->nxD,3,double);
    RAVEL2(bd,u,bd->nxD,6,double);
    RAVEL2(bd,du,bd->nxD,6,double);

    f_getDisplacement(&rav_xDisp,&rav_u,&rav_du,bd->idx);

    UNRAVEL2(bd,xDisp ,bd->nxD,3,double);
    UNRAVEL2(bd,u,bd->nxD,6,double);
    UNRAVEL2(bd,du,bd->nxD,6,double);
}

void BD_setBC(BD_Data *bd)
{
    f_setBC(bd->idx);
}

void BD_writeSolToBin(BD_Data *bd, char* fileName)
{
    int i;
    FILE *fid;
    char fname[128];
    sprintf(fname,"%s.bin",fileName);

    fid = fopen(fname,"ab");  // r for read, b for binary
    if (fid==NULL)
    {
        printf("Unable to open file %s to write the solution of the structural analysis\n",fileName);
        exit(EXIT_FAILURE);
    }

    for (i=0;i<bd->nxD;i++){
        fwrite(bd->u[i],sizeof(double),6,fid);
    }

    fclose(fid);
}

int main(int argc, char *argv[])
{

    BD_Data **bd;
    int i,j,k;
    clock_t start_t, end_t;

    int nBeam = 3;

    bd = (BD_Data**) malloc(nBeam*sizeof(BD_Data*));
    

    for (k=0;k<nBeam;k++)
    {
        bd[k] = (BD_Data*) malloc(sizeof(BD_Data));
        bd[k]->nBeam        = nBeam;
        bd[k]->idx          = k+1;
        bd[k]->DynamicSolve = 1;
        bd[k]->inputFile = "./run/nrel5mw_dynamic/bd_primary_nrel_5mw_dynamic.inp";

        bd[k]->dt        = 1e-3;
        bd[k]->nt        = 100;
        bd[k]->t         = 0;
        
        for (i=0;i<3;i++){
            bd[k]->omega[i]     = 0.0;
            bd[k]->domega[i]    = 0.0;
            bd[k]->gravity[i]   = 0.0;
        }
        BD_initBeamDyn(bd[k]);

        printf("nxD = %d \n",bd[k]->nxD);
        printf("nxL = %d \n",bd[k]->nxL);

        BD_getPositions(bd[k]);
    }

    for (k=0;k<nBeam;k++)
    {
        for (i=0;i<bd[k]->nxL;i++)
        {
            for (j=0; j<6;j++)
            {
                bd[k]->loads[i][j] = 0.0;
            }
            bd[k]->loads[i][0] = 1e3;
        }
    

        printf("Set the loads...\n");
        BD_setLoads(bd[k]);
        printf("Done!\n");
        
        start_t = clock();
        BD_solve(bd[k]);
        end_t = clock();
        printf("Done solving in %1.3e seconds!\n",(double)(end_t - start_t) / CLOCKS_PER_SEC);

        BD_getDisplacement(bd[k]);
        printf("Done!\n");
    }

    FILE *fid;

    fid = fopen("u.bin","wb");  // r for read, b for binary

    for (i=0;i<bd[0]->nxD;i++){
        
        fwrite(bd[0]->u[i],sizeof(double),6,fid);
    }

    fclose(fid);

    printf("Deallocating...\n");

    for (k=0;k<nBeam;k++)
    {
        BD_freeBeamDyn(bd[k],k==nBeam-1);
        free(bd[k]);
    }
    free(bd);

    printf("Done! :) \n\n");
}

/* Get the rotation matrix from the internal Wiener-Milenkovic parameters of GEBT. 
*/
void BD_getRotationMatrix(double out[3][3], double c[3])
{
    double c0,fac;
    c0 = 2.0 - (c[0]*c[0] + c[1]*c[1] + c[2]*c[2])/8.0;
    fac = 1.0/(4.0-c0)/(4.0-c0);

    out[0][0] =  fac * (c0*c0 + c[0]*c[0] - c[1]*c[1] - c[2]*c[2]);
    out[0][1] =  fac * (2*(c[0]*c[1] + c0*c[2]));
    out[0][2] =  fac * (2*(c[0]*c[2]-c0*c[1]));
    out[1][0] =  fac * (2*(c[0]*c[1] - c0*c[2]));
    out[1][1] =  fac * (c0*c0 - c[0]*c[0] + c[1]*c[1] - c[2]*c[2]);
    out[1][2] =  fac * (2*(c[1]*c[2] + c0*c[0]));
    out[2][0] =  fac * (2*(c[0]*c[2] + c0*c[1]));
    out[2][1] =  fac * (2*(c[1]*c[2] - c0*c[0]));
    out[2][2] =  fac * (c0*c0 - c[0]*c[0] - c[1]*c[1] + c[2]*c[2]);

}