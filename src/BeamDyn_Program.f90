PROGRAM BeamDyn_Program

   USE BeamDyn_Usr  ! all other modules inherited through this one

   IMPLICIT NONE

   INTEGER(IntKi) :: i,j
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

   BD_UsrData(1)%omega  = 0.0      ! Angular velocity vector
   BD_UsrData(1)%dOmega = 0.0      ! Angular acceleration vector

   BD_UsrData(1)%InputFile = "bd_primary_nrel_5mw_dynamic.inp"

   CALL BeamDyn_C_Init(BD_UsrData(1))

   ! Loads
   ALLOCATE(loads(BD_UsrData(1)%nxL,6))
   loads(:,:) = 0.0
   loads(:,1) = 5e3
   loads(:,2) = 1e4
   CALL BeamDyn_C_setLoads(BD_UsrData(1),loads)

   write(*,*) BD_UsrData(1)%BD_Parameter%nqp
   write(*,*) "node per elem: ", BD_UsrData(1)%BD_Parameter%nodes_per_elem
   
   DO i=1,25
      CALL CPU_TIME(t_start)
      CALL BeamDyn_C_Solve(BD_UsrData(1))
      CALL CPU_TIME(t_end)
      CALL writeOutputToFile(i)
      print '("Time = ",f6.3," seconds.")',t_end-t_start
   END DO
   
   ! Get displacement
   ALLOCATE(u(BD_UsrData(1)%nxD,6))
   ALLOCATE(du(BD_UsrData(1)%nxD,6))
   CALL BeamDyn_C_getDisp(BD_UsrData(1),u,du)

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
   

END PROGRAM BeamDyn_Program