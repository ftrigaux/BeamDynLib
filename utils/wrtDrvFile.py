
import numpy as np


def writeDvrFile(fileName,t_init=0.0,t_final=1.0,dt=0.1,distrLoad=[],x=np.zeros(0),F=np.zeros((0,0)),primaryFileName="bd_primary_nrel_5mw_dynamic.inp"):
    if distrLoad==[]:
        distrLoad = [0]*6;
        
    n = F.size
    
    fid = open(fileName,"w");
    fid.write("%s"%("------- BEAMDYN Driver with OpenFAST INPUT FILE -------\n"));
    fid.write("%s"%("Dynamic analysis of rotating NREL 5MW blade under gravity force\n"));
    
    fid.write("%s"%("---------------------- SIMULATION CONTROL ----------------------\n"));
    fid.write("%-14s   %s\n"%("True","DynamicSolve"));
    fid.write("%5e     %s\n"%(t_init,"t_initial"));
    fid.write("%5e     %s\n"%(t_final,"t_final"));
    fid.write("%5e     %s\n"%(dt,"dt"));
    
    fid.write("%s"%("---------------------- GRAVITY PARAMETERS ----------------------\n"));
    fid.write("%5e     %s\n"%(0.0,"Gx"));
    fid.write("%5e     %s\n"%(0.0,"Gy"));
    fid.write("%5e     %s\n"%(0.0,"Gz"));
    
    fid.write("%s"%("----------------------  FRAME PARAMETERS  ----------------------\n"));
    fid.write("%5e     %s\n"%(0.0,"GlbPos(1)"));
    fid.write("%5e     %s\n"%(0.0,"GlbPos(2)"));
    fid.write("%5e     %s\n"%(1.0,"GlbPos(3)"));
    
    fid.write("%s"%("---The following 3 by 3 matrix is the direction cosine matirx ,GlbDCM(3,3),\n---relates global frame to reference blade frame\n"));
    fid.write("%1.8e  %1.8e  %1.8e\n"%(1.0,0.0,0.0));
    fid.write("%1.8e  %1.8e  %1.8e\n"%(0.0,1.0,0.0));
    fid.write("%1.8e  %1.8e  %1.8e\n"%(0.0,0.0,1.0));
    fid.write("%-14s   %s\n"%("True","GlbRotBladeT0"));
    
    fid.write("%s"%("----------------------  ROOT VELOCITY PARAMETER  ----------------------\n"));
    fid.write("%5e     %s\n"%(0.0,"RootVel(4)"));
    fid.write("%5e     %s\n"%(0.0,"RootVel(5)"));
    fid.write("%5e     %s\n"%(0.0,"RootVel(6)"));
    
    fid.write("%s"%("----------------------  APPLIED FORCE  ----------------------\n"));
    fid.write("%5e     %s\n"%(distrLoad[0],"DistrLoad(1)"));
    fid.write("%5e     %s\n"%(distrLoad[1],"DistrLoad(2)"));
    fid.write("%5e     %s\n"%(distrLoad[2],"DistrLoad(3)"));
    fid.write("%5e     %s\n"%(distrLoad[3],"DistrLoad(4)"));
    fid.write("%5e     %s\n"%(distrLoad[4],"DistrLoad(5)"));
    fid.write("%5e     %s\n"%(distrLoad[5],"DistrLoad(6)"));
    fid.write("%5e     %s\n"%(0.0,"TipLoad(1)"));
    fid.write("%5e     %s\n"%(0.0,"TipLoad(2)"));
    fid.write("%5e     %s\n"%(0.0,"TipLoad(3)"));
    fid.write("%5e     %s\n"%(0.0,"TipLoad(4)"));
    fid.write("%5e     %s\n"%(0.0,"TipLoad(5)"));
    fid.write("%5e     %s\n"%(0.0,"TipLoad(6)"));
    fid.write("%-14d   %s\n"%(n,"NumPointLoads"));
    fid.write("%s\n"%("Non-dim blade-span eta   Fx          Fy            Fz           Mx           My           Mz\n(-)                      (N)         (N)           (N)          (N-m)        (N-m)        (N-m)"));
    
    for i in range(n):
      fid.write("%1.8e   %1.8e   %1.8e   %1.8e   %1.8e   %1.8e   %1.8e\n"%(x[i],F[i,0],F[i,1],F[i,2],F[i,3],F[i,4],F[i,5]));
    fid.write("%s\n"%("---------------------- PRIMARY INPUT FILE --------------------------------------"));
    fid.write("%s    %s"%(primaryFileName,"InputFile"));
    
    fid.close()