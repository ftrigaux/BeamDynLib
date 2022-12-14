PROGRAM BeamDynLib_Program

   USE BeamDynLib  ! all other modules inherited through this one

   IMPLICIT NONE

   INTEGER(IntKi) :: i,j,iBeam,nBeam=3
   REAL :: t_start, t_end
   REAL(R8Ki),ALLOCATABLE :: loads(:,:),u(:,:),du(:,:)
   character(len=1024) :: filename
   
   TYPE(BD_UsrDataType),ALLOCATABLE :: BD_UsrData(:)

   ALLOCATE(BD_UsrData(nBeam))

   DO iBeam=1,nBeam
   

   BD_UsrData(iBeam)%dt = 4e-3
   BD_UsrData(iBeam)%nt = 10
   BD_UsrData(iBeam)%t  = 0.0

   BD_UsrData(iBeam)%DynamicSolve   = .FALSE. ! flag for dynamic or static solve (static:false)
   BD_UsrData(iBeam)%GlbRotBladeT0  = .TRUE.! Initial blade root orientation is also the GlbRot reference frame
      
   BD_UsrData(iBeam)%GlbPos      = (/0.0, 0.0, 1.0/)       ! Initial vector position 
   BD_UsrData(iBeam)%RootOri     = 0.0 ! DCM of the initial root orientation
   BD_UsrData(iBeam)%GlbRot      = 0.0 ! 
   BD_UsrData(iBeam)%RootRelInit = 0.0 !
   DO i=1,3
      BD_UsrData(iBeam)%RootOri(i,i)     = 1.0 ! DCM of the initial root orientation
      BD_UsrData(iBeam)%GlbRot(i,i)      = 1.0 ! 
      BD_UsrData(iBeam)%RootRelInit(i,i) = 1.0 !
   END DO 

   BD_UsrData(iBeam)%orientation  = 0.0      ! General beam orientation matrix
   BD_UsrData(iBeam)%omega        = 0.0      ! Angular velocity vector
   BD_UsrData(iBeam)%dOmega       = 0.0      ! Angular acceleration vector
   BD_UsrData(iBeam)%grav         = 0.0      ! Angular acceleration vector

   BD_UsrData(iBeam)%omega(1)        = 1.0      ! Angular velocity vector
   BD_UsrData(iBeam)%grav(3)         = 9.0      ! Angular acceleration vector

   BD_UsrData(iBeam)%InputFile = "./run/IEA15MW_dynamic/IEA-15-240-RWT_BeamDyn.dat"

   CALL BeamDyn_C_Init(BD_UsrData(iBeam))

   ! Loads
   IF(.NOT. ALLOCATED(loads)) ALLOCATE(loads(6,BD_UsrData(iBeam)%nxL))
   loads(:,:) = 0.0
   loads(1,:) = 5e3
   loads(2,:) = 1e4
   CALL BeamDyn_C_setLoads(BD_UsrData(iBeam),loads)

   !write(*,*) BD_UsrData(iBeam)%BD_Parameter%nqp
   !write(*,*) "node per elem: ", BD_UsrData(iBeam)%BD_Parameter%nodes_per_elem

   ! To get displacement
   IF(.NOT. ALLOCATED(u))  ALLOCATE( u(6,BD_UsrData(iBeam)%nxD))
   IF(.NOT. ALLOCATED(du)) ALLOCATE(du(6,BD_UsrData(iBeam)%nxD))
   
   DO i=1,10
      CALL CPU_TIME(t_start)
      CALL BeamDyn_C_Solve(BD_UsrData(iBeam))
      CALL CPU_TIME(t_end)
      !CALL writeOutputToFile(i)
      CALL BeamDyn_C_getDisp(BD_UsrData(iBeam),u,du)
      print '("Run Time = ",f6.3," seconds. Tip displacement x-dir at t=",f6.3," [s] = ",f8.5," [m]")',t_end-t_start,BD_UsrData(iBeam)%t,u(1,BD_UsrData(iBeam)%nxD)
   END DO

   write (fileName, "(A12,I1)") "Checkpoint_b", iBeam

   CALL BD_CreateCheckpoint_T(BD_UsrData(iBeam)%t,BD_UsrData(iBeam)%nt,nBeam,BD_UsrData(iBeam),fileName,ErrStat,ErrMsg)

   CALL BeamDyn_C_End(BD_UsrData(iBeam))
   
   END DO

   DEALLOCATE(BD_UsrData)

   !---------------------------------------------------------------------------------------------------------------------------
   !---------------------------------------------------------------------------------------------------------------------------
   ! Restart the simulation from checkpoint
   !---------------------------------------------------------------------------------------------------------------------------
   !---------------------------------------------------------------------------------------------------------------------------

   ALLOCATE(BD_UsrData(nBeam))

   DO iBeam=1,nBeam

   BD_UsrData(iBeam)%dt = 4e-3
   BD_UsrData(iBeam)%nt = 10
   BD_UsrData(iBeam)%t  = 0.0

   BD_UsrData(iBeam)%DynamicSolve   = .TRUE. ! flag for dynamic or static solve (static:false)
   BD_UsrData(iBeam)%GlbRotBladeT0  = .TRUE.! Initial blade root orientation is also the GlbRot reference frame
      
   BD_UsrData(iBeam)%GlbPos      = (/0.0, 0.0, 1.0/)       ! Initial vector position 
   BD_UsrData(iBeam)%RootOri     = 0.0 ! DCM of the initial root orientation
   BD_UsrData(iBeam)%GlbRot      = 0.0 ! 
   BD_UsrData(iBeam)%RootRelInit = 0.0 !
   DO i=1,3
      BD_UsrData(iBeam)%RootOri(i,i)     = 1.0 ! DCM of the initial root orientation
      BD_UsrData(iBeam)%GlbRot(i,i)      = 1.0 ! 
      BD_UsrData(iBeam)%RootRelInit(i,i) = 1.0 !
   END DO 

   BD_UsrData(iBeam)%orientation  = 0.0      ! Angular velocity vector
   BD_UsrData(iBeam)%omega        = 0.0      ! Angular velocity vector
   BD_UsrData(iBeam)%dOmega       = 0.0      ! Angular acceleration vector
   BD_UsrData(iBeam)%grav         = 0.0      ! Angular acceleration vector

   BD_UsrData(iBeam)%omega(1)        = 1.0      ! Angular velocity vector
   BD_UsrData(iBeam)%grav(3)         = 9.0      ! Angular acceleration vector

   BD_UsrData(iBeam)%InputFile = "./run/IEA15MW_dynamic/IEA-15-240-RWT_BeamDyn.dat"

   CALL BeamDyn_C_Init(BD_UsrData(iBeam))

   write (fileName, "(A12,I1)") "Checkpoint_b", iBeam
   CALL BD_RestoreFromCheckpoint_T(BD_UsrData(iBeam)%t,BD_UsrData(iBeam)%nt,nBeam,BD_UsrData(iBeam),fileName,ErrStat,ErrMsg)

   ! Loads
   CALL BeamDyn_C_setLoads(BD_UsrData(iBeam),loads)
   
   DO i=1,10
      CALL CPU_TIME(t_start)
      CALL BeamDyn_C_Solve(BD_UsrData(iBeam))
      CALL CPU_TIME(t_end)
      !CALL writeOutputToFile(i)
      CALL BeamDyn_C_getDisp(BD_UsrData(iBeam),u,du)
      print '("Run Time = ",f6.3," seconds. Tip displacement x-dir at t=",f6.3," [s] = ",f8.5," [m]")',t_end-t_start,BD_UsrData(iBeam)%t,u(1,BD_UsrData(iBeam)%nxD)
   END DO

   CALL BeamDyn_C_getDisp(BD_UsrData(iBeam),u,du)

   print '("After restart: Tip displacement x-dir = ",f6.3," [m]")',u(1,BD_UsrData%nxD)

   CALL BeamDyn_C_End(BD_UsrData(iBeam))

   END DO

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
