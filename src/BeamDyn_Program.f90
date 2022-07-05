PROGRAM BeamDyn_Program

   USE BeamDyn_C  ! all other modules inherited through this one

   IMPLICIT NONE

   REAL :: t_start, t_end
   REAL :: loads(49,6)
   
   loads(:,:) = 0.0
   loads(:,1) = 1e2
   loads(:,2) = 0.0

   BD_UsrData%dt = 4e-3
   BD_UsrData%nt = 10
   BD_UsrData%t  = 0.0

   BD_UsrData%DynamicSolve   = .TRUE. ! flag for dynamic or static solve (static:false)
   BD_UsrData%GlbRotBladeT0  = .TRUE.! Initial blade root orientation is also the GlbRot reference frame
      
   BD_UsrData%GlbPos      = (/0.0, 0.0, 1.0/)       ! Initial vector position 
   BD_UsrData%RootOri     = 0.0 ! DCM of the initial root orientation
   BD_UsrData%GlbRot      = 0.0 ! 
   BD_UsrData%RootRelInit = 0.0 !
   DO i=1,3
      BD_UsrData%RootOri(i,i)     = 1.0 ! DCM of the initial root orientation
      BD_UsrData%GlbRot(i,i)      = 1.0 ! 
      BD_UsrData%RootRelInit(i,i) = 1.0 !
   END DO 

   BD_UsrData%omega  = 0.0      ! Angular velocity vector
   BD_UsrData%dOmega = 0.0     ! Angular acceleration vector

   BD_UsrData%InputFile = "bd_primary_nrel_5mw_dynamic.inp"

   CALL BeamDyn_C_Init()
   CALL BeamDyn_C_setLoads(loads)
   
   DO i=1,25
      CALL CPU_TIME(t_start)
      CALL BeamDyn_C_Solve()
      CALL CPU_TIME(t_end)
      CALL writeOutputToFile(i)
      print '("Time = ",f6.3," seconds.")',t_end-t_start
   END DO
   
   CALL BeamDyn_C_getDisp()

   CALL BeamDyn_C_End()

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
         write(fu,*) BD_Output%BldMotion%Position(j,:)
      end do  
      do j=1,3
         write(fu,*) BD_Output%BldMotion%TranslationDisp(j,:) 
      end do  
      
      close(fu) 
      END SUBROUTINE writeOutputToFile
   

END PROGRAM BeamDyn_Program