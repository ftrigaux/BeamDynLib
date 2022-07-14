
MODULE CInterface

    USE ISO_C_BINDING
    USE BeamDyn_Usr

    IMPLICIT NONE

    TYPE(BD_UsrDataType), ALLOCATABLE :: BD_UsrData(:)
    INTEGER(IntKi)                    :: nBeam
    
    CONTAINS  

    SUBROUTINE initBeamDyn(nBeam_loc,inputFile,                 &
                           dt, nt, t,                           &
                           DynamicSolve,                        &
                           omega, domega,gravity,                &
                           nxLoads, nxDisp) BIND(C,NAME="f_initBeamDyn")

        INTEGER(KIND=C_INT), INTENT(IN), VALUE           :: nBeam_loc
        CHARACTER(C_CHAR),   INTENT(IN)                  :: inputFile(*)
        REAL(KIND=C_DOUBLE), INTENT(IN), OPTIONAL        :: dt, t
        INTEGER(KIND=C_INT), INTENT(IN), OPTIONAL        :: nt, DynamicSolve
        REAL(KIND=C_DOUBLE), INTENT(IN), OPTIONAL        :: omega(3), domega(3), gravity(3)
        INTEGER(KIND=C_INT), INTENT(  OUT)               :: nxLoads,nxDisp        ! Returns the number of nodes for loads and displacement


        INTEGER(IntKi) :: i,j

        ALLOCATE(BD_UsrData(nBeam_loc))

        DO i=1,nBeam_loc
            IF(PRESENT(dt)) THEN; BD_UsrData(i)%dt = dt; ELSE; BD_UsrData(i)%dt = 4e-3;  ENDIF
            IF(PRESENT(nt)) THEN; BD_UsrData(i)%nt = nt; ELSE; BD_UsrData(i)%nt = 10  ;  ENDIF
            IF(PRESENT(t))  THEN; BD_UsrData(i)%t  = t ; ELSE; BD_UsrData(i)%t  = 0.0 ;  ENDIF

            IF(PRESENT(DynamicSolve)) THEN; BD_UsrData(i)%DynamicSolve   = DynamicSolve==1; ELSE; BD_UsrData(i)%DynamicSolve = .TRUE.; ENDIF ! flag for dynamic or static solve (static:false)
            BD_UsrData(i)%GlbRotBladeT0  = .TRUE.! Initial blade root orientation is also the GlbRot reference frame
                
            BD_UsrData(i)%GlbPos      = (/0.0, 0.0, 1.0/)       ! Initial vector position 
            BD_UsrData(i)%RootOri     = 0.0 ! DCM of the initial root orientation
            BD_UsrData(i)%GlbRot      = 0.0 ! 
            BD_UsrData(i)%RootRelInit = 0.0 !
            DO j=1,3
                BD_UsrData(i)%RootOri(j,j)     = 1.0 ! DCM of the initial root orientation
                BD_UsrData(i)%GlbRot(j,j)      = 1.0 ! 
                BD_UsrData(i)%RootRelInit(j,j) = 1.0 !
            END DO 

            IF(PRESENT(omega))   THEN; BD_UsrData(i)%omega(:)    =  omega(:);   ELSE; BD_UsrData(i)%omega    = 0.0; ENDIF    ! Angular velocity vector
            IF(PRESENT(domega))  THEN; BD_UsrData(i)%omega(:)    = domega(:);   ELSE; BD_UsrData(i)%domega   = 0.0; ENDIF    ! Angular acceleration vector
            IF(PRESENT(gravity)) THEN; BD_UsrData(i)%grav(:)     = gravity(:);  ELSE; BD_UsrData(i)%grav     = 0.0; ENDIF    ! Angular acceleration vector

            ! copy string input file from C to Fortran
            j=1
            BD_UsrData(i)%InputFile=" "
            DO WHILE (inputFile(j) /= C_NULL_CHAR .AND. j<LEN(BD_UsrData(i)%InputFile))
                BD_UsrData(i)%InputFile(j:j) = inputFile(j)
                j = j + 1
            ENDDO
            BD_UsrData(i)%InputFile = TRIM(BD_UsrData(i)%InputFile)

            CALL BeamDyn_C_Init(BD_UsrData(i))

        ENDDO

        nBeam = nBeam_loc

        nxLoads = BD_UsrData(1)%BD_Input(1)%DistrLoad%Nnodes
        nxDisp  = BD_UsrData(1)%BD_Output%BldMotion%Nnodes
 
    END SUBROUTINE initBeamDyn

    SUBROUTINE getPositions(xLoads, xDisp) BIND(C,NAME="f_getPositions")
        TYPE(c_ptr), INTENT(OUT) :: xLoads,xDisp

        REAL(R8Ki), POINTER :: f_xLoads(:,:), f_xDisp(:,:)
        INTEGER(IntKi) :: i

        CALL c_f_pointer(xLoads,f_xLoads,[BD_UsrData(1)%BD_Input(1)%DistrLoad%Nnodes,3])
        CALL c_f_pointer(xDisp, f_xDisp, [BD_UsrData(1)%BD_Output%BldMotion%Nnodes,3])

        DO i=1,BD_UsrData(1)%BD_Input(1)%DistrLoad%Nnodes
            f_xLoads(i,1:3) = BD_UsrData(1)%BD_Input(1)%DistrLoad%Position(1:3,i)
        END DO

        DO i=1,BD_UsrData(1)%BD_Output%BldMotion%Nnodes
            f_xDisp(i,1:3) = BD_UsrData(1)%BD_Output%BldMotion%Position(1:3,i)
        END DO

        xLoads = C_LOC(f_xLoads(1,1))
        xDisp  = C_LOC(f_xDisp(1,1))

    END SUBROUTINE getPositions

    SUBROUTINE freeBeamDyn() BIND(C,NAME="f_freeBeamDyn")

        INTEGER(IntKi) :: i
        DO i=1,nBeam
            CALL BeamDyn_C_End(BD_UsrData(i))
            
        ENDDO

        DEALLOCATE(BD_UsrData)

    END SUBROUTINE freeBeamDyn
    

    SUBROUTINE setLoads(loads,idxBeam) BIND(C,NAME="f_setLoads")
        IMPLICIT NONE
        TYPE(C_PTR),INTENT(IN)      :: loads       ! coordinates for key points
        INTEGER(C_INT),INTENT(IN),VALUE:: idxBeam

        REAL(R8Ki), POINTER          :: f_loads(:,:)
        CALL c_f_pointer(loads,f_loads, [BD_UsrData%nxL,6]);
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

        INTEGER(IntKi) :: i,j

        CALL c_f_pointer(x ,f_x ,[BD_UsrData(idxBeam)%nxD,3])
        CALL c_f_pointer(u ,f_u ,[BD_UsrData(idxBeam)%nxD,6])
        CALL c_f_pointer(du,f_du,[BD_UsrData(idxBeam)%nxD,6])

        CALL BeamDyn_C_getDisp(BD_UsrData(1),f_u,f_du)

        !write(*,*) BD_UsrData(idxBeam)%nxD
        !write(*,*) BD_UsrData(idxBeam)%BD_Output%BldMotion%NNodes ! Should be equal!!

        DO i=1,BD_UsrData(idxBeam)%BD_Output%BldMotion%NNodes
            DO j=1,3
                f_x(i,j) = BD_UsrData(idxBeam)%BD_Output%BldMotion%Position(j,i)
            ENDDO
        ENDDO

        x  = C_LOC(f_x(1,1));
        u  = C_LOC(f_u(1,1));
        du = C_LOC(f_du(1,1));


    END SUBROUTINE getDisplacement

    ! This function should be used before every "solve" to set the orientation and displacement of the root
    SUBROUTINE setBC(idxBeam)  BIND(C,NAME="f_setBC")

        INTEGER(C_INT),INTENT(IN),VALUE :: idxBeam

        !sth = sin( theta )
        !cth = cos( theta )
        !
        !Orientation(1,1) =   cth
        !Orientation(2,1) =   -sth
        !Orientation(3,1) =   0.0_R8Ki
        !
        !Orientation(1,2) = sth
        !Orientation(2,2) = cth
        !Orientation(3,2) = 0.0_R8Ki
        !
        !Orientation(1,3) = 0.0_R8Ki
        !Orientation(2,3) = 0.0_R8Ki
        !Orientation(3,3) = 1.0_R8Ki
        !
        !BD_UsrData(idxBeam)%RotationCenter%TranslationDisp(:,1)   = 0.0_ReKi
        !BD_UsrData(idxBeam)%RotationCenter%Orientation(:,:,1)     = matmul(Orientation, matmul(BD_UsrData(idxBeam)%RotationCenter%RefOrientation(:,:,1),BD_UsrData(idxBeam)%RootRelInit))
        !BD_UsrData(idxBeam)%RotationCenter%TranslationVel(:,1)    = 0.0_ReKi
        !BD_UsrData(idxBeam)%RotationCenter%RotationVel(:,1)       = BD_UsrData(idxBeam)%BD_InitInput%RootVel(4:6)
        !BD_UsrData(idxBeam)%RotationCenter%TranslationAcc(:,1)    = 0.0_ReKi
        !BD_UsrData(idxBeam)%RotationCenter%RotationAcc(:,1)       = 0.0_ReKi

    END SUBROUTINE setBC

END MODULE CInterface