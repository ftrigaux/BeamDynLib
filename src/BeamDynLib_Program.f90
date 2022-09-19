PROGRAM BeamDynLib_Program

   USE BeamDynLib  ! all other modules inherited through this one

   IMPLICIT NONE

   INTEGER(IntKi) :: i,j,nBeam=1
   REAL :: t_start, t_end
   REAL(R8Ki),ALLOCATABLE :: loads(:,:),u(:,:),du(:,:)
   
   TYPE(BD_UsrDataType),ALLOCATABLE :: BD_UsrData(:)

   ALLOCATE(BD_UsrData(1))

   BD_UsrData(1)%dt = 4e-3
   BD_UsrData(1)%nt = 10
   BD_UsrData(1)%t  = 0.0

   BD_UsrData(1)%DynamicSolve   = .TRUE. ! flag for dynamic or static solve (static:false)
   BD_UsrData(1)%GlbRotBladeT0  = .TRUE.! Initial blade root orientation is also the GlbRot reference frame
      
   BD_UsrData(1)%GlbPos      = (/0.0, 0.0, 1.0/)       ! Initial vector position 
   BD_UsrData(1)%RootOri     = 0.0 ! DCM of the initial root orientation
   BD_UsrData(1)%GlbRot      = 0.0 ! 
   BD_UsrData(1)%RootRelInit = 0.0 !
   DO i=1,3
      BD_UsrData(1)%RootOri(i,i)     = 1.0 ! DCM of the initial root orientation
      BD_UsrData(1)%GlbRot(i,i)      = 1.0 ! 
      BD_UsrData(1)%RootRelInit(i,i) = 1.0 !
   END DO 

   BD_UsrData(1)%orientation  = 0.0      ! General beam orientation matrix
   BD_UsrData(1)%omega        = 0.0      ! Angular velocity vector
   BD_UsrData(1)%dOmega       = 0.0      ! Angular acceleration vector
   BD_UsrData(1)%grav         = 0.0      ! Angular acceleration vector

   BD_UsrData(1)%omega(1)        = 1.0      ! Angular velocity vector
   BD_UsrData(1)%grav(3)         = 9.0      ! Angular acceleration vector

   BD_UsrData(1)%InputFile = "./run/nrel5mw_dynamic/bd_primary_nrel_5mw_dynamic.inp"

   CALL BeamDyn_C_Init(BD_UsrData(1))

   ! Loads
   ALLOCATE(loads(6,BD_UsrData(1)%nxL))
   loads(:,:) = 0.0
   loads(1,:) = 5e3
   loads(2,:) = 1e4
   CALL BeamDyn_C_setLoads(BD_UsrData(1),loads)

   !write(*,*) BD_UsrData(1)%BD_Parameter%nqp
   !write(*,*) "node per elem: ", BD_UsrData(1)%BD_Parameter%nodes_per_elem

   ! To get displacement
   ALLOCATE( u(6,BD_UsrData(1)%nxD))
   ALLOCATE(du(6,BD_UsrData(1)%nxD))
   
   DO i=1,10
      CALL CPU_TIME(t_start)
      CALL BeamDyn_C_Solve(BD_UsrData(1))
      CALL CPU_TIME(t_end)
      !CALL writeOutputToFile(i)
      CALL BeamDyn_C_getDisp(BD_UsrData(1),u,du)
      print '("Run Time = ",f6.3," seconds. Tip displacement x-dir at t=",f6.3," [s] = ",f8.5," [m]")',t_end-t_start,BD_UsrData(1)%t,u(1,BD_UsrData(1)%nxD)
   END DO

   

   CALL BD_CreateCheckpoint_T(BD_UsrData(1)%t,BD_UsrData(1)%nt,nBeam,BD_UsrData,"Checkpoint",ErrStat,ErrMsg)

   CALL BeamDyn_C_End(BD_UsrData(1))
   
   DEALLOCATE(BD_UsrData)

   !---------------------------------------------------------------------------------------------------------------------------
   !---------------------------------------------------------------------------------------------------------------------------
   ! Restart the simulation from checkpoint
   !---------------------------------------------------------------------------------------------------------------------------
   !---------------------------------------------------------------------------------------------------------------------------

   ALLOCATE(BD_UsrData(1))

   BD_UsrData(1)%dt = 4e-3
   BD_UsrData(1)%nt = 10
   BD_UsrData(1)%t  = 0.0

   BD_UsrData(1)%DynamicSolve   = .TRUE. ! flag for dynamic or static solve (static:false)
   BD_UsrData(1)%GlbRotBladeT0  = .TRUE.! Initial blade root orientation is also the GlbRot reference frame
      
   BD_UsrData(1)%GlbPos      = (/0.0, 0.0, 1.0/)       ! Initial vector position 
   BD_UsrData(1)%RootOri     = 0.0 ! DCM of the initial root orientation
   BD_UsrData(1)%GlbRot      = 0.0 ! 
   BD_UsrData(1)%RootRelInit = 0.0 !
   DO i=1,3
      BD_UsrData(1)%RootOri(i,i)     = 1.0 ! DCM of the initial root orientation
      BD_UsrData(1)%GlbRot(i,i)      = 1.0 ! 
      BD_UsrData(1)%RootRelInit(i,i) = 1.0 !
   END DO 

   BD_UsrData(1)%orientation  = 0.0      ! Angular velocity vector
   BD_UsrData(1)%omega        = 0.0      ! Angular velocity vector
   BD_UsrData(1)%dOmega       = 0.0      ! Angular acceleration vector
   BD_UsrData(1)%grav         = 0.0      ! Angular acceleration vector

   BD_UsrData(1)%omega(1)        = 1.0      ! Angular velocity vector
   BD_UsrData(1)%grav(3)         = 9.0      ! Angular acceleration vector

   BD_UsrData(1)%InputFile = "./run/nrel5mw_dynamic/bd_primary_nrel_5mw_dynamic.inp"

   CALL BeamDyn_C_Init(BD_UsrData(1))

   CALL BD_RestoreFromCheckpoint_T(BD_UsrData(1)%t,BD_UsrData(1)%nt,nBeam,BD_UsrData,"Checkpoint",ErrStat,ErrMsg)

   ! Loads
   CALL BeamDyn_C_setLoads(BD_UsrData(1),loads)
   
   DO i=1,10
      CALL CPU_TIME(t_start)
      CALL BeamDyn_C_Solve(BD_UsrData(1))
      CALL CPU_TIME(t_end)
      !CALL writeOutputToFile(i)
      CALL BeamDyn_C_getDisp(BD_UsrData(1),u,du)
      print '("Run Time = ",f6.3," seconds. Tip displacement x-dir at t=",f6.3," [s] = ",f8.5," [m]")',t_end-t_start,BD_UsrData(1)%t,u(1,BD_UsrData(1)%nxD)
   END DO

   CALL BeamDyn_C_getDisp(BD_UsrData(1),u,du)

   print '("After restart: Tip displacement x-dir = ",f6.3," [m]")',u(1,BD_UsrData%nxD)

   CALL BeamDyn_C_End(BD_UsrData(1))

   DEALLOCATE(BD_UsrData)
   DEALLOCATE(loads)
   DEALLOCATE(u)
   DEALLOCATE(du)

   contains
   SUBROUTINE writeOutputToFile(no)

      INTEGER(intKi), INTENT(IN) :: no

      integer :: fu
      character(len=50) :: file_name
      CHARACTER(len=10) :: file_id

      ! Write the integer into a string:
      write(file_id, '(i0)') no

      ! Construct the filename:
      file_name = 'disp' // trim(adjustl(file_id)) // '.dat'

      ! Open the file with this name
      open(file = trim(file_name), unit = fu, status='new')
      do j=1,3
         write(fu,*) BD_UsrData(1)%BD_Output%BldMotion%Position(j,:)
      end do  
      do j=1,3
         write(fu,*) BD_UsrData(1)%BD_Output%BldMotion%TranslationDisp(j,:) 
      end do  
      
      close(fu) 
      END SUBROUTINE writeOutputToFile
   

END PROGRAM BeamDynLib_Program
