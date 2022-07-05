#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct 
{
  int n;                           // Number of force point loads
  double t_initial, t_final, dt;   // Temporal values
  double DistrLoad[6];             // Values for distributed loads
  double *x;                       // Location of force points loads 
  double *F[6];                    // Values of force point loads 
} Param;


void wrtDrvFile(Param *pm)
{
  int i;
  FILE *fid = fopen("bd_driver_file.inp","w");

  fprintf(fid,"%s","------- BEAMDYN Driver with OpenFAST INPUT FILE -------\n");
  fprintf(fid,"%s","Dynamic analysis of rotating NREL 5MW blade under gravity force\n");

  fprintf(fid,"%s","---------------------- SIMULATION CONTROL ----------------------\n");
  fprintf(fid,"%-14s   %s\n","True","DynamicSolve");
  fprintf(fid,"%5e     %s\n",pm->t_initial,"t_initial");
  fprintf(fid,"%5e     %s\n",pm->t_final,"t_final");
  fprintf(fid,"%5e     %s\n",pm->dt,"dt");

  fprintf(fid,"%s","---------------------- GRAVITY PARAMETERS ----------------------\n");
  fprintf(fid,"%5e     %s\n",0.0,"Gx");
  fprintf(fid,"%5e     %s\n",0.0,"Gy");
  fprintf(fid,"%5e     %s\n",0.0,"Gz");

  fprintf(fid,"%s","----------------------  FRAME PARAMETERS  ----------------------\n");
  fprintf(fid,"%5e     %s\n",0.0,"GlbPos(1)");
  fprintf(fid,"%5e     %s\n",0.0,"GlbPos(2)");
  fprintf(fid,"%5e     %s\n",1.0,"GlbPos(3)");

  fprintf(fid,"%s","---The following 3 by 3 matrix is the direction cosine matirx ,GlbDCM(3,3),\n---relates global frame to reference blade frame\n");
  fprintf(fid,"%1.8e  %1.8e  %1.8e\n",1.0,0.0,0.0);
  fprintf(fid,"%1.8e  %1.8e  %1.8e\n",0.0,1.0,0.0);
  fprintf(fid,"%1.8e  %1.8e  %1.8e\n",0.0,0.0,1.0);
  fprintf(fid,"%-14s   %s\n","True","GlbRotBladeT0");

  fprintf(fid,"%s","----------------------  ROOT VELOCITY PARAMETER  ----------------------\n");
  fprintf(fid,"%5e     %s\n",0.0,"RootVel(4)");
  fprintf(fid,"%5e     %s\n",0.0,"RootVel(5)");
  fprintf(fid,"%5e     %s\n",0.0,"RootVel(6)");

  fprintf(fid,"%s","----------------------  APPLIED FORCE  ----------------------\n");
  fprintf(fid,"%5e     %s\n",0.0,"DistrLoad(1)");
  fprintf(fid,"%5e     %s\n",0.0,"DistrLoad(2)");
  fprintf(fid,"%5e     %s\n",0.0,"DistrLoad(3)");
  fprintf(fid,"%5e     %s\n",0.0,"DistrLoad(4)");
  fprintf(fid,"%5e     %s\n",0.0,"DistrLoad(5)");
  fprintf(fid,"%5e     %s\n",0.0,"DistrLoad(6)");
  fprintf(fid,"%5e     %s\n",0.0,"TipLoad(1)");
  fprintf(fid,"%5e     %s\n",0.0,"TipLoad(2)");
  fprintf(fid,"%5e     %s\n",0.0,"TipLoad(3)");
  fprintf(fid,"%5e     %s\n",0.0,"TipLoad(4)");
  fprintf(fid,"%5e     %s\n",0.0,"TipLoad(5)");
  fprintf(fid,"%5e     %s\n",0.0,"TipLoad(6)");
  fprintf(fid,"%-14d   %s\n",pm->n,"NumPointLoads");
  fprintf(fid,"%s\n","Non-dim blade-span eta   Fx          Fy            Fz           Mx           My           Mz\n(-)                      (N)         (N)           (N)          (N-m)        (N-m)        (N-m)");

  for (i=0;i<pm->n;i++)
  {
    fprintf(fid,"%1.8e   %1.8e   %1.8e   %1.8e   %1.8e   %1.8e   %1.8e\n",pm->x[i],pm->F[0][i],pm->F[1][i],pm->F[2][i],pm->F[3][i],pm->F[4][i],pm->F[5][i]);
  }
  fprintf(fid,"%s\n","---------------------- PRIMARY INPUT FILE --------------------------------------");
  fprintf(fid,"%s    %s","bd_primary_nrel_5mw_dynamic.inp","InputFile");


  fclose(fid);
}



int main(int argc, char *argv[])
{
  int i,j,n;
  double dr;
  Param *pm = (Param*) malloc(sizeof(Param));

  if (argc>1)
  {
    pm->n = atoi(argv[1]);
  }
  else
  {
    pm->n = 0;
  }

  if (argc>2)
  {
    pm->dt = atof(argv[1]);
  }
  else
  {
    pm->dt = 1.0;
  }
  
  pm->t_final   = 0.0;
  pm->t_initial = 0.0;

  dr = 1.0/((double) n);
  
  for (i=0;i<6;i++)
  {
    pm->F[i] = calloc(n,sizeof(double));
  }

  for (j=0;j<n;j++)
    {
      pm->F[0][j] = 1.0;
      pm->F[1][j] = 1.0;
      pm->F[2][j] = 1.0;
      pm->F[3][j] = 1.0;
      pm->F[4][j] = 1.0;
      pm->F[5][j] = 1.0;
    }

  pm->x = calloc(n,sizeof(double));
  for (i=0;i<n;i++)
  {
      pm->x[i] = dr/2+i*dr;
  }
  
  
  
  wrtDrvFile(pm);


  for (i=0;i<6;i++)
  {
    free(pm->F[i]);
  }
  free(pm->x);
  free(pm);
}



