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

int i_,j_; // iterator for macros

void BD_initBeamDyn(BD_Data *bd)
{
    f_initBeamDyn(bd->nBeam,bd->inputFile,
                &bd->dt, &bd->nt, &bd->t,          
                &bd->DynamicSolve,       
                bd->omega, bd->domega, bd->gravity,      
                &bd->nxL, &bd->nxD);

    ALLOCATE2(bd->xLoads,double,bd->nxL,3);
    ALLOCATE2(bd->xDisp, double,bd->nxD,3);
    ALLOCATE2(bd->loads, double,bd->nxL,6);
    ALLOCATE2(bd->u,     double,bd->nxD,6);
    ALLOCATE2(bd->du,    double,bd->nxD,6);

}

void BD_getPositions(BD_Data *bd)
{
    RAVEL2(bd,xLoads,bd->nxL,3,double);
    RAVEL2(bd,xDisp ,bd->nxL,3,double);
    f_getPositions(&rav_xLoads, &rav_xDisp);

    UNRAVEL2(bd,xLoads,bd->nxL,3,double);
    UNRAVEL2(bd,xDisp ,bd->nxL,3,double);

}

void BD_freeBeamDyn(BD_Data *bd)
{
    f_freeBeamDyn();

    DEALLOCATE2(bd->xLoads,bd->nxL);
    DEALLOCATE2(bd->xDisp ,bd->nxD);
    DEALLOCATE2(bd->loads ,bd->nxL);
    DEALLOCATE2(bd->u     ,bd->nxD);
    DEALLOCATE2(bd->du    ,bd->nxD);
}

void BD_setLoads(BD_Data *bd)
{
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

int main(int argc, char *argv[])
{
    BD_Data *bd;
    int i,j;
    double **loads, **x, **u, **du;
    clock_t start_t, end_t;

    bd = (BD_Data*) malloc(sizeof(BD_Data));

    bd->nBeam     = 1;
    bd->idx       = 1;
    bd->DynamicSolve = 1;
    bd->inputFile = "./run/nrel5mw_dynamic/bd_primary_nrel_5mw_dynamic.inp";

    bd->dt        = 1e-3;
    bd->nt        = 100;
    bd->t         = 0;
    
    for (i=0;i<3;i++){
        bd->omega[i]     = 0.0;
        bd->domega[i]    = 0.0;
        bd->gravity[i]   = 0.0;
    }

    BD_initBeamDyn(bd);

    printf("nxD = %d \n",bd->nxD);
    printf("nxL = %d \n",bd->nxL);

    BD_getPositions(bd);

    for (i=0;i<bd->nxL;i++)
    {
        for (j=0; j<6;j++)
        {
            bd->loads[i][j] = 0.0;
        }
        bd->loads[i][0] = 1e3;
    }

    printf("Set the loads...\n");
    BD_setLoads(bd);
    printf("Done!\n");
    
    start_t = clock();
    BD_solve(bd);
    end_t = clock();
    printf("Done solving in %1.3e seconds!\n",(double)(end_t - start_t) / CLOCKS_PER_SEC);

    BD_getDisplacement(bd);
    printf("Done!\n");

    FILE *fid;

    fid = fopen("u.bin","wb");  // r for read, b for binary

    for (i=0;i<bd->nxD;i++){
        fwrite(bd->u[i],sizeof(double),6,fid);
    }

    fclose(fid);

    printf("Deallocating...\n");

    BD_freeBeamDyn(bd);
    free(bd);

    printf("Done! :) \n\n");
}