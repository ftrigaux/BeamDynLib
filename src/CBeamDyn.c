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
                bd->GlbPos, &bd->GlbRotBladeT0, bd->GlbRot, bd->RootOri,
                &bd->nxL, &bd->nxD);

    ALLOCATE2(bd->xLoads,double,3,bd->nxL);
    ALLOCATE2(bd->xDisp, double,3,bd->nxD);
    ALLOCATE2(bd->loads, double,6,bd->nxL);
    ALLOCATE2(bd->u,     double,6,bd->nxD);
    ALLOCATE2(bd->du,    double,6,bd->nxD);

    FILLZERO2(bd->xLoads,3, bd->nxL);
    FILLZERO2(bd->xDisp, 3, bd->nxD);
    FILLZERO2(bd->loads, 6, bd->nxL);
    FILLZERO2(bd->u,     6, bd->nxD);
    FILLZERO2(bd->du,    6, bd->nxD);

}

void BD_refresh(BD_Data *bd)
{
    f_refresh(bd->idx,
                &bd->dt, &bd->nt, &bd->t);

}

void BD_getPositions(BD_Data *bd)
{
    RAVEL2(bd,xLoads,3,bd->nxL,double);
    RAVEL2(bd,xDisp ,3,bd->nxL,double);
    f_getPositions(&rav_xLoads, &rav_xDisp);

    UNRAVEL2(bd,xLoads,3,bd->nxL,double);
    UNRAVEL2(bd,xDisp ,3,bd->nxL,double);

}

void BD_freeBeamDyn(BD_Data *bd, int deallocAll)
{
    f_freeBeamDyn(bd->idx,deallocAll);

    DEALLOCATE2(bd->xLoads,3);
    DEALLOCATE2(bd->xDisp ,3);
    DEALLOCATE2(bd->loads ,6);
    DEALLOCATE2(bd->u     ,6);
    DEALLOCATE2(bd->du    ,6);
}


void BD_setLoads(BD_Data *bd)
{
    int i,j;
    RAVEL2(bd,loads,6,bd->nxL,double);
    f_setLoads(&rav_loads,bd->idx);
    UNRAVEL2(bd,loads,6,bd->nxL,double);
}

void BD_solve(BD_Data *bd)
{
    f_solve(bd->idx);

    bd->t += bd->nt * bd->dt;
}

void BD_getDisplacement(BD_Data *bd)
{
    RAVEL2(bd,xDisp ,3,bd->nxD,double);
    RAVEL2(bd,u     ,6,bd->nxD,double);
    RAVEL2(bd,du    ,6,bd->nxD,double);

    f_getDisplacement(&rav_xDisp,&rav_u,&rav_du,bd->idx);

    UNRAVEL2(bd,xDisp ,3,bd->nxD,double);
    UNRAVEL2(bd,u     ,6,bd->nxD,double);
    UNRAVEL2(bd,du    ,6,bd->nxD,double);
}

void BD_setBC(BD_Data *bd)
{
    f_setBC(bd->idx, bd->omega, bd->domega);
}

void BD_writeSolToBin(BD_Data *bd, char* fileName)
{
    int i;
    FILE *fid;
    char fname[128];
    sprintf(fname,"%s.bin",fileName);

    fid = fopen(fname,"ab");
    if (fid==NULL)
    {
        printf("Unable to open file %s to write the solution of the structural analysis\n",fileName);
        exit(EXIT_FAILURE);
    }

    for (i=0;i<6;i++){
        fwrite(bd->u[i],sizeof(double),bd->nxD,fid);
    }

    fclose(fid);
}

int main(int argc, char *argv[])
{

    BD_Data **bd;
    int i,j,k;
    int nt_glob;
    clock_t start_t, end_t;

    int nBeam = 1;
    nt_glob   = 30;

    bd = (BD_Data**) malloc(nBeam*sizeof(BD_Data*));
    

    for (k=0;k<nBeam;k++)
    {
        bd[k] = (BD_Data*) malloc(sizeof(BD_Data));

        // Set the number of beams, the index of each beam, dynamic solve and input file
        bd[k]->nBeam        = nBeam;
        bd[k]->idx          = k+1;
        bd[k]->DynamicSolve = 1;
        bd[k]->inputFile = "./run/nrel5mw_dynamic/bd_primary_nrel_5mw_dynamic.inp";

        // Set the time steps and the number of time steps
        bd[k]->dt        = 5e-2;
        bd[k]->nt        = 1;
        bd[k]->t         = 0;
        
        // Set the rotation speed and the gravity constrains
        for (i=0;i<3;i++){
            bd[k]->omega[i]     = 0.0;
            bd[k]->domega[i]    = 0.0;
            bd[k]->gravity[i]   = 0.0;
        }
        bd[k]->omega[0] = -1.0;

        // Set the position vector, the global orientation and the initial root orientation
        bd[k]->GlbRotBladeT0 = 1;
        for (i=0;i<3;i++){
            bd[k]->GlbPos[i]    = 0.0;
            for (j=0;j<3;j++){
                bd[k]->GlbRot[i][j]  = 0;
                bd[k]->RootOri[i][j] = 0;
                if (i==j)
                {
                    bd[k]->GlbRot[i][j]  = 1;
                    bd[k]->RootOri[i][j] = 1;
                }
            }
        }
        bd[k]->GlbPos[2]    = 1.5;

        BD_initBeamDyn(bd[k]);

        printf("nxD = %d \n",bd[k]->nxD);
        printf("nxL = %d \n",bd[k]->nxL);

        BD_getPositions(bd[k]);
    }

    for (k=0;k<nBeam;k++)
    {
        for (i=0;i<6;i++)
        {
            for (j=0; j<bd[k]->nxL;j++)
            {
                bd[k]->loads[i][j] = 0.0;
            }
        }
    

        printf("Set the loads...\n");
        BD_setLoads(bd[k]);
        printf("Done!\n");
        
        for (i=0;i<nt_glob;i++){
            start_t = clock();
            bd[k]->omega[0] = -1.0;
            BD_setBC(bd[k]);
            BD_solve(bd[k]);
            end_t = clock();
            printf("Done solving in %1.3e seconds!\n",(double)(end_t - start_t) / CLOCKS_PER_SEC);

            BD_getDisplacement(bd[k]);
        }
        printf("Done!\n");
    }

    FILE *fid;

    fid = fopen("u.bin","wb");  // r for read, b for binary

    for (i=0;i<6;i++){
        
        fwrite(bd[0]->u[i],sizeof(double),bd[0]->nxD,fid);
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

    return(EXIT_SUCCESS);
}

/* Get the rotation matrix from the internal Wiener-Milenkovic parameters of GEBT. 
*/
void BD_getRotationMatrix(double out[][3], double c[3])
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

    
    double c1,c2,c3,tr0;
    c1  = c[0]/4.0;
    c2  = c[1]/4.0;
    c3  = c[2]/4.0;
    c0  = 0.5 * (1.0-c1*c1-c2*c2-c3*c3);     // 1/4 the value of the the AIAA paper (after plugging in c1, c2, c3 conversions)


    tr0 = 1.0 - c0;                          // This is 1/4 the value of the AIAA paper, after converting c0.
    tr0 = 2.0/(tr0*tr0);                     // this is 32x the equivalent term from the AIAA paper.   This is well behaved and won't go to zero.

    // The following terms can be shown to match the transpose of the DCM given in the AIAA paper.
    out[0][0] = tr0*(c1*c1 + c0*c0) - 1.0;
    out[1][0] = tr0*(c1*c2 + c0*c3);
    out[2][0] = tr0*(c1*c3 - c0*c2);

    out[0][1] = tr0*(c1*c2 - c0*c3);
    out[1][1] = tr0*(c2*c2 + c0*c0) - 1.0;
    out[2][1] = tr0*(c2*c3 + c0*c1);

    out[0][2] = tr0*(c1*c3 + c0*c2);
    out[1][2] = tr0*(c2*c3 - c0*c1);
    out[2][2] = tr0*(c3*c3 + c0*c0) - 1.0;

}