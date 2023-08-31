
MODULE BeamDynLib_CBind

    USE ISO_C_BINDING
    USE BeamDynLib

    IMPLICIT NONE

    TYPE(BD_UsrDataType), ALLOCATABLE :: BD_UsrData(:)
    
    CONTAINS  

    ! initBeamDyn : initiate the BeamDyn parameters.
    ! if GlbRotBladeT0 is True, then the GlbRot matrix will be set equal to RootOri. This mean that the initial Root Orientation is the global frame.
    ! if not, GlbRot will be set to [I], and the RootIni will be used as initial Orientation relative to the global frame.
    ! 
    SUBROUTINE initBeamDyn(nBeam,inputFile,idx,                    &
                           dt, nt, t,                              &
                           DynamicSolve,                           &
                           omega, domega, gravity,                 &
                           GlbPos, GlbRotBladeT0, RootOri,         &
                           WrVTK, VTK_fps,                         &
                           nxLoads, nxDisp) BIND(C,NAME="f_initBeamDyn")

        INTEGER(KIND=C_INT), INTENT(IN), VALUE           :: nBeam,idx
        CHARACTER(C_CHAR),   INTENT(IN)                  :: inputFile(*)
        REAL(KIND=C_DOUBLE), INTENT(IN), OPTIONAL        :: dt, t
        INTEGER(KIND=C_INT), INTENT(IN), OPTIONAL        :: nt, DynamicSolve, GlbRotBladeT0
        REAL(KIND=C_DOUBLE), INTENT(IN), OPTIONAL        :: omega(3), domega(3), gravity(3)
        REAL(KIND=C_DOUBLE), INTENT(IN), OPTIONAL        :: GlbPos(3), RootOri(3,3)
        INTEGER(KIND=C_INT), INTENT(IN), OPTIONAL        :: WrVTK, VTK_fps
        INTEGER(KIND=C_INT), INTENT(  OUT)               :: nxLoads,nxDisp        ! Returns the number of nodes for loads and displacement


        INTEGER(IntKi)  :: j
        CHARACTER(512) :: RootName

        IF (.NOT. ALLOCATED(BD_UsrData)) ALLOCATE(BD_UsrData(nBeam))

        
        IF(PRESENT(dt)) THEN; BD_UsrData(idx)%dt = dt; ELSE; BD_UsrData(idx)%dt = 4e-3;  ENDIF
        IF(PRESENT(nt)) THEN; BD_UsrData(idx)%nt = nt; ELSE; BD_UsrData(idx)%nt = 10  ;  ENDIF
        IF(PRESENT(t))  THEN; BD_UsrData(idx)%t  = t ; ELSE; BD_UsrData(idx)%t  = 0.0 ;  ENDIF

        IF(PRESENT(DynamicSolve))  THEN; BD_UsrData(idx)%DynamicSolve    = DynamicSolve==1;  ELSE; BD_UsrData(idx)%DynamicSolve = .TRUE.; ENDIF ! flag for dynamic or static solve (static:false)
        IF(PRESENT(GlbRotBladeT0)) THEN; BD_UsrData(idx)%GlbRotBladeT0   = GlbRotBladeT0==1; ELSE; BD_UsrData(idx)%GlbRotBladeT0 = .TRUE.; ENDIF ! Initial blade root orientation is also the GlbRot reference frame
            
        IF(PRESENT(GlbPos))      THEN; BD_UsrData(idx)%GlbPos(:)      = GlbPos(:)  ; ELSE; BD_UsrData(idx)%GlbPos      =  (/0.0, 0.0, 1.0/); ENDIF ! Initial vector position 
        IF(PRESENT(RootOri))     THEN; BD_UsrData(idx)%RootOri        = RootOri    ; ELSE; BD_UsrData(idx)%RootOri     = 0.0;                ENDIF ! DCM of the root orientation
        BD_UsrData(idx)%GlbRot      = 0.0; ! Will be set in BD_initFromUsrData()
        DO j=1,3
            BD_UsrData(idx)%GlbRot(j,j)      = 1.0
            IF( .NOT. PRESENT(RootOri))     BD_UsrData(idx)%RootOri(j,j)     = 1.0
        END DO 

        BD_UsrData(idx)%theta_rot = 0.0;
        BD_UsrData(idx)%PAngInp_rad = 0.0;
        IF(PRESENT(omega))   THEN; BD_UsrData(idx)%omega(:)    =  omega(:);   ELSE; BD_UsrData(idx)%omega    = 0.0; ENDIF    ! Angular velocity vector
        IF(PRESENT(domega))  THEN; BD_UsrData(idx)%domega(:)   = domega(:);   ELSE; BD_UsrData(idx)%domega   = 0.0; ENDIF    ! Angular acceleration vector
        IF(PRESENT(gravity)) THEN; BD_UsrData(idx)%grav(:)     = gravity(:);  ELSE; BD_UsrData(idx)%grav     = 0.0; ENDIF    ! Angular acceleration vector

        IF(PRESENT(gravity)) THEN; BD_UsrData(idx)%WrVTK     = WrVTK;    ELSE; BD_UsrData(idx)%WrVTK     = 0; ENDIF    ! Angular acceleration vector
        IF(PRESENT(gravity)) THEN; BD_UsrData(idx)%VTK_fps   = VTK_fps;  ELSE; BD_UsrData(idx)%VTK_fps   = 0.0; ENDIF    ! Angular acceleration vector

        ! copy string input file from C to Fortran
        j=1
        BD_UsrData(idx)%InputFile=" "
        DO WHILE (inputFile(j) /= C_NULL_CHAR .AND. j<LEN(BD_UsrData(idx)%InputFile))
            BD_UsrData(idx)%InputFile(j:j) = inputFile(j)
            j = j + 1
        ENDDO
        BD_UsrData(idx)%InputFile = TRIM(BD_UsrData(idx)%InputFile)

        ! Output file
        CALL GetRoot(BD_UsrData(idx)%InputFile,RootName)
        WRITE(BD_UsrData(idx)%OutputFile,"(A,A1,I0.2)") TRIM(RootName),"_",idx
        WRITE(*,*) "Output File Name:",BD_UsrData(idx)%OutputFile

        CALL BeamDyn_C_Init(BD_UsrData(idx))

        nxLoads = BD_UsrData(idx)%BD_Input(1)%DistrLoad%Nnodes
        nxDisp  = BD_UsrData(idx)%BD_Output%BldMotion%Nnodes
 
    END SUBROUTINE initBeamDyn

    SUBROUTINE BDrefresh(idx,                      &
                         dt, nt, t) BIND(C,NAME="f_refresh")

        INTEGER(KIND=C_INT), INTENT(IN), VALUE           :: idx
        REAL(KIND=C_DOUBLE), INTENT(IN), OPTIONAL        :: dt, t
        INTEGER(KIND=C_INT), INTENT(IN), OPTIONAL        :: nt
        
        IF(PRESENT(dt)) BD_UsrData(idx)%dt = dt
        IF(PRESENT(dt)) BD_UsrData(idx)%BD_Parameter%dt = dt
        IF(PRESENT(nt)) BD_UsrData(idx)%nt = nt
        IF(PRESENT(t))  BD_UsrData(idx)%t  = t

        CALL BD_TiSchmComputeCoefficients2(BD_UsrData(idx)%BD_Parameter) ! Need to recompute the generalized alpha time scheme coefficients
 
    END SUBROUTINE BDrefresh

    SUBROUTINE getPositions(xLoads, xDisp) BIND(C,NAME="f_getPositions")
        TYPE(c_ptr), INTENT(OUT) :: xLoads,xDisp

        REAL(R8Ki), POINTER :: f_xLoads(:,:), f_xDisp(:,:)
        INTEGER(IntKi) :: i

        CALL c_f_pointer(xLoads,f_xLoads,[3,BD_UsrData(1)%BD_Input(1)%DistrLoad%Nnodes])
        CALL c_f_pointer(xDisp, f_xDisp, [3,BD_UsrData(1)%BD_Output%BldMotion%Nnodes])

        DO i=1,BD_UsrData(1)%BD_Input(1)%DistrLoad%Nnodes
            f_xLoads(1:3,i) = MATMUL(BD_UsrData(1)%GlbRot, BD_UsrData(1)%BD_Input(1)%DistrLoad%Position(1:3,i) - BD_UsrData(1)%GlbPos(1:3))
        END DO

        DO i=1,BD_UsrData(1)%BD_Output%BldMotion%Nnodes
            f_xDisp(1:3,i)  = MATMUL(BD_UsrData(1)%GlbRot, BD_UsrData(1)%BD_Output%BldMotion%Position(1:3,i) - BD_UsrData(1)%GlbPos(1:3))
        END DO

        xLoads = C_LOC(f_xLoads(1,1))
        xDisp  = C_LOC(f_xDisp(1,1))

    END SUBROUTINE getPositions

    SUBROUTINE freeBeamDyn(idxBeam, deallocAll) BIND(C,NAME="f_freeBeamDyn")

        INTEGER(C_INT),INTENT(IN),VALUE:: idxBeam
        INTEGER(C_INT),INTENT(IN),VALUE:: deallocAll
        
        CALL BeamDyn_C_End(BD_UsrData(idxBeam))

        IF (deallocAll .eq. 1) THEN
            DEALLOCATE(BD_UsrData)
        ENDIF

    END SUBROUTINE freeBeamDyn
    

    SUBROUTINE setLoads(loads,idxBeam) BIND(C,NAME="f_setLoads")
        IMPLICIT NONE
        TYPE(C_PTR),INTENT(IN)         :: loads       ! loads in the blade root frame [r] (i.e. the loads rotate with the blade)
        INTEGER(C_INT),INTENT(IN),VALUE:: idxBeam

        REAL(ReKi), POINTER          :: f_loads(:,:)
        CALL c_f_pointer(loads,f_loads, [6,BD_UsrData%nxL])
        CALL BeamDyn_C_setLoads(BD_UsrData(idxBeam),f_loads)

    END SUBROUTINE setLoads


    SUBROUTINE solve(idxBeam) BIND(C,NAME="f_solve")

        INTEGER(C_INT),INTENT(IN),VALUE:: idxBeam

        REAL(DbKi) :: t_start, t_end

        CALL CPU_TIME(t_start)
        CALL BeamDyn_C_Solve(BD_UsrData(idxBeam))
        CALL CPU_TIME(t_end)
        !print '("Structural solver is done: CPU_Time = ",f6.3," [s]")',t_end-t_start

    END SUBROUTINE solve

    SUBROUTINE getDisplacement(x,u,du,idxBeam) BIND(C,NAME="f_getDisplacement")

        TYPE(C_PTR),INTENT(OUT)      :: x, u, du
        INTEGER(C_INT),INTENT(IN),VALUE :: idxBeam

        REAL(R8Ki), POINTER          :: f_x(:,:)
        REAL(R8Ki), POINTER          :: f_u(:,:)
        REAL(R8Ki), POINTER          :: f_du(:,:)

        INTEGER(IntKi) :: i

        CALL c_f_pointer(x ,f_x ,[3,BD_UsrData(idxBeam)%nxD])
        CALL c_f_pointer(u ,f_u ,[6,BD_UsrData(idxBeam)%nxD])
        CALL c_f_pointer(du,f_du,[6,BD_UsrData(idxBeam)%nxD])

        CALL BeamDyn_C_getDisp(BD_UsrData(idxBeam),f_u,f_du)

        DO i=1,BD_UsrData(idxBeam)%BD_Output%BldMotion%NNodes
            f_x(1:3,i)  = MATMUL(BD_UsrData(idxBeam)%GlbRot, BD_UsrData(idxBeam)%BD_Output%BldMotion%Position(1:3,i) - BD_UsrData(idxBeam)%GlbPos(1:3))
        ENDDO

        x  = C_LOC(f_x(1,1));
        u  = C_LOC(f_u(1,1));
        du = C_LOC(f_du(1,1));


    END SUBROUTINE getDisplacement

    SUBROUTINE getReactionForce(x,F,idxBeam) BIND(C,NAME="f_getReactionForce")

        TYPE(C_PTR),INTENT(OUT)      :: x, F
        INTEGER(C_INT),INTENT(IN),VALUE :: idxBeam

        REAL(R8Ki), POINTER          :: f_x(:,:)
        REAL(R8Ki), POINTER          :: f_F(:,:)

        INTEGER(IntKi) :: i

        CALL c_f_pointer(x ,f_x ,[3,BD_UsrData(idxBeam)%nxD])
        CALL c_f_pointer(F ,f_F ,[6,BD_UsrData(idxBeam)%nxD])

        CALL BeamDyn_C_getReactionForce(BD_UsrData(idxBeam),f_F)

        DO i=1,BD_UsrData(idxBeam)%BD_Output%ReactionForce%NNodes
            f_x(1:3,i)  = MATMUL(BD_UsrData(idxBeam)%GlbRot, BD_UsrData(idxBeam)%BD_Output%ReactionForce%Position(1:3,i) - BD_UsrData(idxBeam)%GlbPos(1:3))
        ENDDO

        x  = C_LOC(f_x(1,1));
        F  = C_LOC(f_F(1,1));

    END SUBROUTINE getReactionForce

    ! This function should be used before every "solve" to set the orientation and displacement of the root
    SUBROUTINE setBC(idxBeam,            &
                     omega, domega,      &
                     forceTheta, theta_rot, &
                     PAngInp_rad)  BIND(C,NAME="f_setBC")

        INTEGER(C_INT),INTENT(IN),VALUE :: idxBeam

        REAL(KIND=C_DOUBLE)          :: omega(3)
        REAL(KIND=C_DOUBLE)          :: dOmega(3)
        INTEGER(KIND=C_INT) , VALUE  :: forceTheta
        REAL(KIND=C_DOUBLE) , VALUE  :: theta_rot
        REAL(KIND=C_DOUBLE) , VALUE  :: PAngInp_rad

        BD_UsrData(idxBeam)%omega(1:3)  = omega(1:3);
        BD_UsrData(idxBeam)%dOmega(1:3) = dOmega(1:3);
        IF (forceTheta .eq. 1) THEN
            BD_UsrData(idxBeam)%theta_rot = REAL(theta_rot,ReKi);
        ENDIF
        BD_UsrData(idxBeam)%PAngInp_rad = REAL(PAngInp_rad,ReKi);

    END SUBROUTINE setBC


     ! This function writes a restart file
    SUBROUTINE writeRestartFile(idxBeam, restartFile)  BIND(C,NAME="f_writeRestartFile")

        INTEGER(KIND=C_INT), INTENT(IN), VALUE           :: idxBeam
        CHARACTER(C_CHAR),   INTENT(IN)                  :: restartFile(*)

        INTEGER(IntKi)         :: j
        CHARACTER(1024)        :: FRestartFile

        ! Read the restart file name from C 
        j=1
        FRestartFile=" "
        DO WHILE (restartFile(j) /= C_NULL_CHAR .AND. j<LEN(FRestartFile))
            FRestartFile(j:j) = restartFile(j)
            j = j + 1
        ENDDO
        FRestartFile = TRIM(FRestartFile)

        ! Create the file
        CALL BD_CreateCheckpoint_T(BD_UsrData(idxBeam)%t,BD_UsrData(idxBeam)%nt,1,BD_UsrData(idxBeam),FRestartFile,ErrStat,ErrMsg)
        CALL CheckError(BD_UsrData(idxBeam), ErrStat, ErrMsg, "during restart writing")
        

    END SUBROUTINE writeRestartFile


    ! This function reads a restart file
    SUBROUTINE readRestartFile(idxBeam, restartFile, t, nt)  BIND(C,NAME="f_readRestartFile")

        INTEGER(KIND=C_INT),  INTENT(IN), VALUE              :: idxBeam
        CHARACTER(C_CHAR),    INTENT(IN)                     :: restartFile(*)
        REAL(KIND=C_DOUBLE),  INTENT(OUT)                    :: t
        INTEGER(KIND=C_INT),  INTENT(OUT)                    :: nt

        INTEGER(IntKi)         :: j,nBeamLoc
        CHARACTER(1024)        :: FRestartFile

        ! Read the restart file name from C 
        j=1
        FRestartFile=" "
        DO WHILE (restartFile(j) /= C_NULL_CHAR .AND. j<LEN(FRestartFile))
            FRestartFile(j:j) = restartFile(j)
            j = j + 1
        ENDDO
        FRestartFile = TRIM(FRestartFile)

        ! Read the file
        CALL BD_RestoreFromCheckpoint_T(BD_UsrData(idxBeam)%t,BD_UsrData(idxBeam)%nt,nBeamLoc,BD_UsrData(idxBeam),FRestartFile,ErrStat,ErrMsg)
        CALL CheckError(BD_UsrData(idxBeam), ErrStat, ErrMsg, "during restart initialization")

        t  = REAL(BD_UsrData(idxBeam)%t, C_DOUBLE)
        nt = INT(BD_UsrData(idxBeam)%nt, C_INT)
        

    END SUBROUTINE readRestartFile

END MODULE BeamDynLib_CBind