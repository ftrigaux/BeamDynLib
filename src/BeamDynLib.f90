!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2022  Francois Trigaux (UCLouvain, BELGIUM)
! 
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
MODULE BeamDynLib

   !USE BeamDyn_driver_subs  ! all other modules inherited through this one
   !USE BeamDyn
   !USE BeamDyn_Subs
   USE BeamDynLib_Types

   IMPLICIT NONE

   ! global glue-code-specific variables

   INTEGER(IntKi)                   :: ErrStat          ! Error status of the operation
   CHARACTER(1024)                  :: ErrMsg           ! Error message if ErrStat /= ErrID_None
   INTEGER(IntKi), parameter        :: BD_interp_order = 1  ! order of interpolation/extrapolation

   LOGICAL, PARAMETER               :: writeOutputFile = .FALSE. ! Whether to write the Fast output file for BeamDyn. Ok for debug but not ok for multiple turbine

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   !TYPE(BD_InitInputType)           :: BD_InitInput
   !TYPE(BD_ParameterType)           :: BD_Parameter
   !TYPE(BD_ContinuousStateType)     :: BD_ContinuousState
   !TYPE(BD_InitOutputType)          :: BD_InitOutput
   !TYPE(BD_DiscreteStateType)       :: BD_DiscreteState
   !TYPE(BD_ConstraintStateType)     :: BD_ConstraintState
   !TYPE(BD_OtherStateType)          :: BD_OtherState
   !TYPE(BD_MiscVarType)             :: BD_MiscVar
   !TYPE(BD_InputType) ,ALLOCATABLE  :: BD_Input(:)
   !REAL(DbKi),         ALLOCATABLE  :: BD_InputTimes(:)
   !TYPE(BD_OutputType)              :: BD_Output
   INTEGER(IntKi)                   :: DvrOut 
   
   !TYPE(BD_UsrDataType)             :: BD_UsrData
   
   !TYPE(BD_DriverInternalType)      :: DvrData

   ! local variables
   
   !CHARACTER(256)                   :: DvrInputFile
   !CHARACTER(256)                   :: RootName
   !INTEGER(IntKi)                   :: j             ! counter for various loops
   !INTEGER(IntKi)                   :: i             ! counter for various loops   
   INTEGER(IntKi)                   :: max_ld_step=8 ! maximum load steps for static runs.
   !REAL(DbKi)                       :: TiLstPrn      ! The simulation time of the last print (to file) [(s)]
   !REAL(ReKi)                       :: PrevClockTime ! Clock time at start of simulation in seconds [(s)]
   !REAL(ReKi)                       :: UsrTime1      ! User CPU time for simulation initialization [(s)]
   !REAL(ReKi)                       :: UsrTime2      ! User CPU time for simulation (without intialization) [(s)]
   !INTEGER(IntKi) , DIMENSION(1:8)  :: StrtTime      ! Start time of simulation (including intialization) [-]
   !INTEGER(IntKi) , DIMENSION(1:8)  :: SimStrtTime   ! Start time of simulation (after initialization) [-]
   !CHARACTER(200)                   :: git_commit    ! String containing the current git commit hash

   !TYPE(ProgDesc), PARAMETER        :: version   = ProgDesc( 'BeamDyn Driver', '', '' )  ! The version number of this program.
   
   CONTAINS

   ! Initialization
   SUBROUTINE BeamDyn_C_Init(usr)
      IMPLICIT NONE
      TYPE(BD_UsrDataType)   :: usr

      INTEGER(IntKi)                   :: i,j

      !CALL DATE_AND_TIME ( Values=StrtTime )                 ! Let's time the whole simulation
      !CALL CPU_TIME ( UsrTime1 )                             ! Initial time (this zeros the start time when used as a MATLAB function)
      !UsrTime1 = MAX( 0.0_ReKi, UsrTime1 )                   ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned

      CALL NWTC_Init()

      ! Read Driver input file is replaced by an initialization through the UsrData structure
      !CALL GET_COMMAND_ARGUMENT(1,DvrInputFile)
      !CALL GetRoot(DvrInputFile,RootName)
      !CALL BD_ReadDvrFile(DvrInputFile,dt_global,BD_InitInput,DvrData,ErrStat,ErrMsg)
      CALL BD_initFromUsrData(usr);
      CALL CheckError(usr,ErrStat,ErrMsg,"during BD_initFromUsrData() function of BeamDyn")

      !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
      ALLOCATE(usr%BD_Input(BD_interp_order + 1)) 
      ALLOCATE(usr%BD_InputTimes(BD_interp_order + 1)) 

      CALL BD_Init(usr%BD_InitInput            &
                     , usr%BD_Input(1)         &
                     , usr%BD_Parameter        &
                     , usr%BD_ContinuousState  &
                     , usr%BD_DiscreteState    &
                     , usr%BD_ConstraintState  &
                     , usr%BD_OtherState       &
                     , usr%BD_Output           &
                     , usr%BD_MiscVar          &
                     , usr%dt                  &
                     , usr%BD_InitOutput       &
                     , ErrStat             &
                     , ErrMsg )
         CALL CheckError(usr,ErrStat,ErrMsg,"during BD_Init() function of BeamDyn")
      
         ! If the Quasi-Static solve is in use, rerun the initialization with loads at t=0 
         ! (HACK: set in the driver only because computing Jacobians with this option [as in FAST glue code] is problematic)
      usr%BD_OtherState%RunQuasiStaticInit = usr%BD_Parameter%analysis_type == BD_DYN_SSS_ANALYSIS


      ! Set the Initial root orientation
      ! BD_Input(1)%RootMotion%Orientation(1:3,1:3,1) = DvrData%RootRelInit
      usr%BD_Input(1)%RootMotion%Orientation(1:3,1:3,1) = usr%RootRelInit
         
      CALL BDuser_InitRotationCenterMesh(usr)
         CALL CheckError(usr,ErrStat,ErrMsg,"during BDuser_InitRotationCenterMesh() function of BeamDyn")

      CALL initUsrData(usr)

      ! 
      ! CALL CreateMultiPointMeshes(DvrData,BD_InitInput,BD_InitOutput,BD_Parameter, BD_Output, BD_Input(1), ErrStat, ErrMsg)   
      ! CALL Transfer_MultipointLoads(DvrData, BD_Output, BD_Input(1), ErrStat, ErrMsg)   
      
      IF (writeOutputFile) THEN
         CALL Dvr_InitializeOutputFile(DvrOut,usr%BD_InitOutput,usr%OutputFile,ErrStat,ErrMsg)
         CALL CheckError(usr,ErrStat,ErrMsg,"during Dvr_InitializeOutputFile() function of BeamDyn")
      ENDIF

         
         ! initialize BD_Input and BD_InputTimes
      usr%BD_InputTimes(1) = usr%t
      !CALL BD_InputSolve( BD_InputTimes(1), BD_Input(1), DvrData, ErrStat, ErrMsg)
      CALL BDUsr_InputSolve(usr,1)
      
      DO j = 2,BD_interp_order+1
            ! create new meshes
         CALL BD_CopyInput (usr%BD_Input(1) , usr%BD_Input(j) , MESH_NEWCOPY, ErrStat, ErrMsg)
            CALL CheckError(usr,ErrStat,ErrMsg,"during BD_CopyInput() function of BeamDyn")
            
            ! solve for inputs at previous time steps
         usr%BD_InputTimes(j) = usr%t - (j - 1) * usr%dt
         ! CALL BD_InputSolve( BD_InputTimes(j), BD_Input(j), DvrData, ErrStat, ErrMsg)
         CALL BDUsr_InputSolve(usr,j)
            CALL CheckError(usr,ErrStat,ErrMsg,"during inputSolve() function of BeamDyn")
      END DO

      ! Write VTK reference if requested (ref is (0,0,0)
      if (usr%WrVTK > 0) then
         call SetVTKvars(usr)
         call MeshWrVTKreference( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), usr%BD_Output%BldMotion,   trim(usr%VTK_OutFileRoot)//'_BldMotion', ErrStat, ErrMsg );  CALL CheckError(usr,ErrStat,ErrMsg,"during MeshWrVTKreference() function of BeamDyn")
         call MeshWrVTKreference( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), usr%BD_Input(1)%PointLoad, trim(usr%VTK_OutFileRoot)//'_PointLoad', ErrStat, ErrMsg );  CALL CheckError(usr,ErrStat,ErrMsg,"during MeshWrVTKreference() function of BeamDyn")
         call MeshWrVTKreference( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), usr%BD_Input(1)%DistrLoad, trim(usr%VTK_OutFileRoot)//'_DistrLoad', ErrStat, ErrMsg );  CALL CheckError(usr,ErrStat,ErrMsg,"during MeshWrVTKreference() function of BeamDyn")
      endif
         ! Write VTK reference if requested (ref is (0,0,0)
      if (usr%WrVTK == 2) then
         usr%n_t_vtk = 0
         call MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), usr%BD_Output%BldMotion,   trim(usr%VTK_OutFileRoot)//'_BldMotion',  usr%n_t_vtk, .true., ErrStat, ErrMsg, usr%VTK_tWidth )
         call MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), usr%BD_Input(1)%PointLoad, trim(usr%VTK_OutFileRoot)//'_PointLoad',  usr%n_t_vtk, .true., ErrStat, ErrMsg, usr%VTK_tWidth )
         IF (usr%BD_Parameter%BldMotionNodeLoc == BD_MESH_QP) THEN
            CALL MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), usr%BD_Input(1)%DistrLoad, trim(usr%VTK_OutFileRoot)//'_DistrLoad',  usr%n_t_vtk, .true., ErrStat, ErrMsg, usr%VTK_tWidth, Sib=usr%BD_Output%BldMotion )
         ELSE
            CALL MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), usr%BD_Input(1)%DistrLoad, trim(usr%VTK_OutFileRoot)//'_DistrLoad',  usr%n_t_vtk, .true., ErrStat, ErrMsg, usr%VTK_tWidth)
         ENDIF
         CALL CheckError(usr,ErrStat,ErrMsg,"during MeshWrVTK() function of BeamDyn")
      endif


   END SUBROUTINE BeamDyn_C_Init

   SUBROUTINE BeamDyn_C_Solve(usr)
      IMPLICIT NONE
      TYPE(BD_UsrDataType) :: usr

      INTEGER(IntKi)  :: it
      INTEGER(IntKi)  :: i,j

      !.........................
      ! calculate outputs at t=0
      !.........................
   CALL BD_CalcOutput( usr%t, usr%BD_Input(1), usr%BD_Parameter, usr%BD_ContinuousState, usr%BD_DiscreteState, &
                           usr%BD_ConstraintState, usr%BD_OtherState,  usr%BD_Output, usr%BD_MiscVar, ErrStat, ErrMsg)
      CALL CheckError(usr,ErrStat,ErrMsg,"during BD_CalcOutput() function of BeamDyn")
   
   IF (writeOutputFile) THEN
      CALL Dvr_WriteOutputLine(usr%t,DvrOut,usr%BD_Parameter%OutFmt,usr%BD_Output)
   ENDIF 
   
      !.........................
      ! time marching
      !.........................
     
   
   DO it = 1, usr%nt  !Changed 0,nt to 1,nt

      ! Shift "window" of BD_Input 
      DO j = BD_interp_order, 1, -1
         CALL BD_CopyInput (usr%BD_Input(j),  usr%BD_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
            CALL CheckError(usr,ErrStat,ErrMsg,"during BD_CopyInput() function of BeamDyn")
         usr%BD_InputTimes(j+1) = usr%BD_InputTimes(j)
      END DO
      
      usr%BD_InputTimes(1)  = usr%t + usr%dt
      ! CALL BDUsr_UpdateOrientation(usr) ! This is not longer necessary as the orientation is updated in BD_inputSolve
      !CALL BD_InputSolve( BD_InputTimes(1), BD_Input(1), DvrData, ErrStat, ErrMsg)
      CALL BDUsr_InputSolve(usr,1)
         CALL CheckError(usr,ErrStat,ErrMsg,"during BD_InputSolve() function of BeamDyn")
      
                       
     IF(usr%BD_Parameter%analysis_type .EQ. BD_STATIC_ANALYSIS .AND. it > max_ld_step) EXIT

      ! update states from n_t_global to n_t_global + 1
     CALL BD_UpdateStates( usr%t, it, usr%BD_Input, usr%BD_InputTimes, usr%BD_Parameter, &
                               usr%BD_ContinuousState, &
                               usr%BD_DiscreteState, usr%BD_ConstraintState, &
                               usr%BD_OtherState, usr%BD_MiscVar, ErrStat, ErrMsg )
        CALL CheckError(usr,ErrStat,ErrMsg,"during BD_UpdateStates() function of BeamDyn")

        
      ! advance time
     usr%t = usr%t + usr%dt
           
      ! calculate outputs at n_t_global + 1
     CALL BD_CalcOutput( usr%t, usr%BD_Input(1), usr%BD_Parameter, usr%BD_ContinuousState, usr%BD_DiscreteState, &
                             usr%BD_ConstraintState, usr%BD_OtherState,  usr%BD_Output, usr%BD_MiscVar, ErrStat, ErrMsg)
        CALL CheckError(usr,ErrStat,ErrMsg,"during BD_CalcOutput() function of BeamDyn")

     !CALL Dvr_WriteOutputLine(t_global,DvrOut,BD_Parameter%OutFmt,BD_Output)

     if (usr%WrVTK == 2) then
         if ( MOD( NINT(usr%t / usr%dt), usr%n_VTKTime ) == 0 ) then
            usr%n_t_vtk = usr%n_t_vtk + 1
            call MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), usr%BD_Output%BldMotion,   trim(usr%VTK_OutFileRoot)//'_BldMotion', usr%n_t_vtk, .true., ErrStat, ErrMsg, usr%VTK_tWidth )
            call MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), usr%BD_Input(1)%PointLoad, trim(usr%VTK_OutFileRoot)//'_PointLoad', usr%n_t_vtk, .true., ErrStat, ErrMsg, usr%VTK_tWidth )
            IF (usr%BD_Parameter%BldMotionNodeLoc == BD_MESH_QP) THEN
               CALL MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), usr%BD_Input(1)%DistrLoad, trim(usr%VTK_OutFileRoot)//'_DistrLoad',  usr%n_t_vtk, .true., ErrStat, ErrMsg, usr%VTK_tWidth, Sib=usr%BD_Output%BldMotion )
            ELSE
               CALL MeshWrVTK( (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /), usr%BD_Input(1)%DistrLoad, trim(usr%VTK_OutFileRoot)//'_DistrLoad',  usr%n_t_vtk, .true., ErrStat, ErrMsg, usr%VTK_tWidth)
            ENDIF
            CALL CheckError(usr,ErrStat,ErrMsg,"during MeshWrVTK() function of BeamDyn")

         endif
      endif

   ENDDO

   !CALL RunTimes( StrtTime, UsrTime1, SimStrtTime, UsrTime2, BD_UsrData%t )

   END SUBROUTINE BeamDyn_C_Solve

   SUBROUTINE BeamDyn_C_End(usr)
      IMPLICIT NONE
      TYPE(BD_UsrDataType) :: usr

      INTEGER(IntKi)               :: i,j
      character(ErrMsgLen)         :: errMsg2                 ! temporary Error message if ErrStat /=
      integer(IntKi)               :: errStat2                ! temporary Error status of the operation
      character(*), parameter      :: RoutineName = 'Dvr_End'

      IF(DvrOut >0) CLOSE(DvrOut)

      IF ( ALLOCATED(usr%BD_Input) ) THEN
         CALL BD_End( usr%BD_Input(1), usr%BD_Parameter, usr%BD_ContinuousState, usr%BD_DiscreteState, &
               usr%BD_ConstraintState, usr%BD_OtherState, usr%BD_Output, usr%BD_MiscVar, ErrStat2, ErrMsg2 )
            CALL SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      
         DO i=2,BD_interp_order + 1
            CALL BD_DestroyInput( usr%BD_Input(i), ErrStat2, ErrMsg2 )
            CALL SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         ENDDO
         
         DEALLOCATE(usr%BD_Input)
      END IF

      IF(ALLOCATED(usr%BD_InputTimes )) DEALLOCATE(usr%BD_InputTimes )
      !if(allocated(DvrData%MultiPointLoad)) deallocate(DvrData%MultiPointLoad)

      CALL freeUsrData(usr)
      
      if (ErrStat >= AbortErrLev) then      
         CALL ProgAbort( 'BeamDyn Driver encountered simulation error level: '&
             //TRIM(GetErrStr(ErrStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      else
         ! We don't want to call NormStop() here because it calls ProgExit which ends the program. We just want to free data.
         !CALL NormStop()
         IF ( LEN_TRIM(ProgName) > 0 ) THEN
            CALL WrScr   ( NewLine//' '//TRIM( ProgName )//' terminated normally.' )
         ELSE
            CALL WrScr   ( NewLine//' Program terminated normally.' )
         END IF
         CALL WrScr    ( '' )
         !CALL ProgExit ( 0 )
      end if

   END SUBROUTINE BeamDyn_C_End

   SUBROUTINE CheckError(usr,ErrID,Msg,ErrLocMsg)
      TYPE(BD_UsrDataType)                 :: usr 
      INTEGER(IntKi), INTENT(IN)           :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN)           :: Msg         ! The error message (ErrMsg)
      CHARACTER(*),   INTENT(IN), OPTIONAL :: ErrLocMsg   ! an optional message describing the location of the error

      if (ErrStat /= ErrID_None) then
         CALL WrScr(TRIM(ErrMsg))
         
         if (ErrStat >= AbortErrLev) then
            CALL BeamDyn_C_End(usr)
            CALL ProgExit ( ErrStat )
         end if
      end if
         
   end SUBROUTINE CheckError

   SUBROUTINE SetVTKvars(usr)

      TYPE(BD_UsrDataType), INTENT(INOUT) :: usr
      real(R8Ki)  :: TmpTime
      real(R8Ki)  :: TmpRate

      usr%VTK_OutFileRoot = 'bd_vtk/'//trim(usr%OutputFile)
      usr%n_t_vtk = 0    ! first VTK output number

      ! now save the number of time steps between VTK file output:
      if (usr%VTK_fps == 0) then
         usr%n_VTKTime = HUGE(usr%n_VTKTime)
      else
         TmpTime = 1.0_R8Ki / usr%VTK_fps
         usr%n_VTKTime = NINT( TmpTime / usr%dt )
         ! I'll warn if p%n_VTKTime*p%DT is not TmpTime
         IF (usr%WrVTK == 2) THEN
            TmpRate = usr%n_VTKTime*usr%dt
            if (.not. EqualRealNos(TmpRate, TmpTime)) then
               call WrScr('1/VTK_fps is not an integer multiple of DT. FAST will output VTK information at '//&
                              trim(num2lstr(1.0_DbKi/TmpRate))//' fps, the closest rate possible.')
            end if
         end if
         CALL MKDIR('bd_vtk')
      end if

      !usr%VTK_tWidth = CEILING( log10( real(n_t_final, ReKi) / usr%n_VTKTime ) ) + 1
      usr%VTK_tWidth  = 6_IntKi  ! Hard-code this to 6 numbers because the total number of time-steps is unknown
   end SUBROUTINE SetVTKvars
 

   ! Creates a mesh for the user distributed loads
   ! SUBROUTINE initUsrData(Pos,NNodes)

   !    INTEGER(IntKi), INTENT(IN)      :: NNodes   
   !    REAL(r8Ki)    , INTENT(IN)      :: Pos(NNodes,3)

   !    real(r8Ki)                                 :: orientation(3,3)

   !    integer(intKi)                  :: ErrStat2          ! temporary Error status
   !    character(ErrMsgLen)            :: ErrMsg2           ! temporary Error message
   !    character(*), parameter         :: RoutineName = 'InitUsrLoads'

   !    ! Initialize mesh for the user loads
   !    CALL MeshCreate( BlankMesh    = usrLoadsMesh       &
   !                   , IOS          = COMPONENT_INPUT    &
   !                   , NNodes       = Nnodes   &
   !                   , Orientation  = .TRUE. &
   !                   , Force        = .TRUE. &
   !                   , Moment       = .True. &
   !                   , ErrStat      = ErrStat2 &
   !                   , ErrMess      = ErrMsg2)
   !    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !    if (ErrStat>=AbortErrLev) return

   !    ! Assign node position
   !    CALL eye(orientation,ErrStat2, ErrMsg2); ! Here, we assume no twist
   !    DO i=1,NNodes
   !       CALL MeshPositionNode ( Mesh    = usrLoadsMesh  &
   !                              ,INode   = i     &
   !                              ,Pos     = Pos(i,:)          &
   !                              ,ErrStat = ErrStat2     &
   !                              ,ErrMess = ErrMsg2      &
   !                              ,Orient  = orientation ) ! Orientation is set to global frame
   !          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !    ENDDO
      
   !    ! Construct elements
   !    DO i=1,NNodes-1
   !    CALL MeshConstructElement( Mesh      = usrLoadsMesh     &
   !                               ,Xelement = ELEMENT_LINE2    &
   !                               ,P1       = i                &
   !                               ,P2       = i+1              &
   !                               ,ErrStat  = ErrStat2         &
   !                               ,ErrMess  = ErrMsg2          )
   !       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !    ENDDO

   !    ! Commit the mesh
   !    CALL MeshCommit ( Mesh    = usrLoadsMesh      &
   !                      ,ErrStat = ErrStat2        &
   !                      ,ErrMess = ErrMsg2         )
   !    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !    ! initial guesses
   !    usrLoadsMesh%Force  = 0.0_ReKi
   !    usrLoadsMesh%Moment = 0.0_ReKi

   !    ! Initialize the mapping between user loads and blade mesh
   !    WRITE(*,*) "Create Laod Map"
   !    CALL MeshMapCreate( usrLoadsMesh, BD_Input(1)%DistrLoad, usrLoadsMap, ErrStat2, ErrMsg2 )
   !    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !    ! Copy the mesh for the blade displacement (same positions and orientation, but different fields)
   !    CALL MeshCopy( usrLoadsMesh, usrDispMesh, MESH_SIBLING, ErrStat2, ErrMsg2 &
   !                 , IOS              = COMPONENT_OUTPUT                        &
   !                 , TranslationDisp  = .TRUE.                                  &
   !                 , Orientation      = .TRUE.                                  &
   !                 , TranslationVel   = .TRUE.                                  &
   !                 , RotationVel      = .TRUE.)
   !    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !    ! Initialize the mapping between user displacement and blade displacement
   !    WRITE(*,*) "Create Disp Map"

   !    CALL MeshMapCreate( usrDispMesh, BD_Output%BldMotion, usrDispMap, ErrStat2, ErrMsg2 )
   !    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !    END SUBROUTINE initUsrData

   !    SUBROUTINE setUsrLoads(NNodes, DistrLoad)

   !       INTEGER(IntKi)    , INTENT(IN)      :: NNodes   
   !       REAL(DbKi)        , INTENT(IN)      :: DistrLoad(NNodes,6)

   !       integer(intKi)                  :: ErrStat2          ! temporary Error status
   !       character(ErrMsgLen)            :: ErrMsg2           ! temporary Error message
   !       character(*), parameter         :: RoutineName = 'setUsrLoads'

   !       DO i=1,usrLoadsMesh%NNodes
   !          usrLoadsMesh%Force(:,i)  =  DistrLoad(i,1:3)
   !          usrLoadsMesh%Moment(:,i) =  DistrLoad(i,4:6)
   !       ENDDO

   !       ! Transfer distributed loads from the user defined loads to the BD inputs
   !       CALL Transfer_Line2_to_Line2(usrLoadsMesh, BD_Input(1)%DistrLoad, usrLoadsMap, ErrStat2, ErrMsg2 &
   !                                    , SrcDisp = usrLoadsMesh &
   !                                    , DestDisp = BD_Output%BldMotion)
   !       CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !    END SUBROUTINE setUsrLoads

   !    SUBROUTINE getUsrDisp(Disp,Vel,NNodes)
   !       INTEGER(IntKi)    , INTENT(IN)      :: NNodes 
   !       Real(DbKi)        , INTENT(OUT)     :: Disp(NNodes,6)
   !       Real(DbKi)        , INTENT(OUT)     :: Vel(NNodes,6)

   !       Real(r8Ki) :: angles(3)

   !       integer(intKi)                  :: ErrStat2          ! temporary Error status
   !       character(ErrMsgLen)            :: ErrMsg2           ! temporary Error message
   !       character(*), parameter         :: RoutineName = 'getUsrDisp'

   !       CALL Transfer_Line2_to_Line2(BD_Output%BldMotion, usrDispMesh, usrDispMap, ErrStat2, ErrMsg2)
   !       CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !       DO i=1,NNodes
   !          Disp(i,1:3) = usrDispMesh%TranslationDisp(1:3,i)
   !          CALL BD_CrvExtractCrv(usrDispMesh%Orientation(1:3,1:3,i),angles,ErrStat,ErrMsg)
   !          Disp(i,4:6) = angles

   !          Vel(i,1:3) = usrDispMesh%TranslationVel(1:3,i)
   !          Vel(i,4:6) = usrDispMesh%RotationVel(1:3,i)

   !       ENDDO
   !    END SUBROUTINE getUsrDisp

   !    SUBROUTINE freeUsrData()
   !       CALL MeshDestroy(usrDispMesh,ErrStat,ErrMsg)
   !       CALL MeshDestroy(usrLoadsMesh,ErrStat,ErrMsg)

   !       CALL NWTC_Library_Destroymeshmaptype( usrDispMap, ErrStat, ErrMsg )
   !       CALL NWTC_Library_Destroymeshmaptype( usrLoadsMap, ErrStat, ErrMsg )
   !    END SUBROUTINE freeUsrData


   SUBROUTINE initUsrData(usr)
      TYPE(BD_UsrDataType) :: usr

      INTEGER(IntKi)                   :: i,j

      usr%nxD = usr%BD_Output%BldMotion%NNodes
      usr%nxL = usr%BD_Input(1)%DistrLoad%Nnodes

      ALLOCATE(usr%loads(6,usr%nxL))
      usr%loads(:,:) = 0.0_ReKi

   END SUBROUTINE initUsrData

   SUBROUTINE freeUsrData(usr)
      TYPE(BD_UsrDataType) :: usr
      
      CALL MeshDestroy(usr%RotationCenter, ErrStat, ErrMsg);
      IF (ALLOCATED(usr%loads)) THEN
         DEALLOCATE(usr%loads)
      ENDIF
   END SUBROUTINE freeUsrData


   SUBROUTINE BeamDyn_C_SetLoads(usr,loads)
      TYPE(BD_UsrDataType)             :: usr
      REAL(ReKi),INTENT(IN)            :: loads(:,:)
      INTEGER(IntKi)                   :: k

      DO k=1,usr%nxL
         usr%loads(1:6,k) = loads(1:6,k)
      ENDDO

   END SUBROUTINE BeamDyn_C_SetLoads

   SUBROUTINE BeamDyn_C_GetDisp(usr,u,du)
      
      TYPE(BD_UsrDataType),INTENT(IN)  :: usr
      REAL(R8Ki),          INTENT(OUT) :: u(:,:) 
      REAL(R8Ki),          INTENT(OUT) :: du(:,:)

      INTEGER(IntKi)                   :: i,j
      REAL(R8Ki)                       :: RootRelOrient(3,3),temp_vec(3),temp_vec2(3),d(3),d_ref(3),temp33(3,3),temp33_2(3,3)

      INTEGER(IntKi)               :: ErrStat2                     ! Temporary Error status
      CHARACTER(ErrMsgLen)         :: ErrMsg2                      ! Temporary Error message
      character(*), parameter      :: RoutineName = 'BeamDyn_C_GetDisp'

      DO i=1,usr%BD_Output%BldMotion%Nnodes
         ! Translation displacement
         RootRelOrient = matmul( transpose(usr%BD_Input(1)%RootMotion%Orientation(:,:,1)), usr%BD_Input(1)%RootMotion%RefOrientation(:,:,1))
         d     = usr%BD_Output%BldMotion%TranslationDisp(:, i) - usr%BD_Input(1)%RootMotion%TranslationDisp(:,1)
         d_ref = usr%BD_Output%BldMotion%Position(       :, i) - usr%BD_Input(1)%RootMotion%Position(       :,1)
         temp_vec2 = d + d_ref - matmul( RootRelOrient, d_ref ) ! tip displacement
         u(1:3,i) = MATMUL(usr%BD_Input(1)%RootMotion%Orientation(:,:,1),temp_vec2)

         ! Tip angular/rotational deflection Wiener-Milenkovic parameter (relative to the undeflected orientation) expressed in r
         CALL LAPACK_DGEMM('N', 'T', 1.0_BDKi, usr%BD_Output%BldMotion%RefOrientation(:,:,i), RootRelOrient,   0.0_BDKi, temp33_2, ErrStat2, ErrMsg2 )
         CALL LAPACK_DGEMM('T', 'N', 1.0_BDKi, usr%BD_Output%BldMotion%Orientation(   :,:,i), temp33_2,        0.0_BDKi, temp33,   ErrStat2, ErrMsg2 )
         CALL BD_CrvExtractCrv(temp33,temp_vec2, ErrStat2, ErrMsg2) ! temp_vec2 = the Wiener-Milenkovic parameters of the tip angular/rotational defelctions
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
         temp_vec = MATMUL(usr%BD_Input(1)%RootMotion%Orientation(:,:,1),temp_vec2) ! translate these parameters to the correct system for output

         u(4:6,i)  = temp_vec;

         ! Velocities are expressed in the blade root frame [r]
         ! du = (Rotation of the blade) * [TranslationVel - omega(r0+r)]
         d         = usr%BD_Output%BldMotion%Position(:, i) + usr%BD_Output%BldMotion%TranslationDisp(:, i)
         du(1:3,i) = MATMUL(usr%BD_Input(1)%RootMotion%Orientation(:,:,1), usr%BD_Output%BldMotion%TranslationVel(:,i) - cross_product(usr%BD_Input(1)%RootMotion%RotationVel(:,1),d))
         du(4:6,i) = MATMUL(usr%BD_Input(1)%RootMotion%Orientation(:,:,1), usr%BD_Output%BldMotion%RotationVel(:,i))
      ENDDO

   END SUBROUTINE BeamDyn_C_GetDisp

   SUBROUTINE BeamDyn_C_getReactionForce(usr, F)

   TYPE(BD_UsrDataType),INTENT(IN)  :: usr
   REAL(R8Ki),          INTENT(OUT) :: F(:,:)

   INTEGER(IntKi)                   :: i
   REAL(R8Ki)                       :: temp_vec(3)

   INTEGER(IntKi)               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                      ! Temporary Error message
   character(*), parameter      :: RoutineName = 'BeamDyn_C_GetDisp'

   DO i=1,usr%BD_Output%BldMotion%Nnodes
      SELECT CASE (usr%BD_Parameter%BldMotionNodeLoc)
      CASE (BD_MESH_FE)
         ! -- Expressed in local airfoil frame --
         !F(1:3,i) = MATMUL(usr%BD_Output%BldMotion%Orientation(:,:,i), usr%BD_MiscVar%BldInternalForceFE(1:3,i))
         !F(4:6,i) = MATMUL(usr%BD_Output%BldMotion%Orientation(:,:,i), usr%BD_MiscVar%BldInternalForceFE(4:6,i))
         ! -- Expressed in root frame --
         F(1:3,i) = MATMUL(usr%BD_MiscVar%u2%RootMotion%Orientation(:,:,1), usr%BD_MiscVar%BldInternalForceFE(1:3,i))
         F(4:6,i) = MATMUL(usr%BD_MiscVar%u2%RootMotion%Orientation(:,:,1), usr%BD_MiscVar%BldInternalForceFE(4:6,i))
      CASE (BD_MESH_QP)
         ! -- Expressed in local airfoil frame --
         !F(1:3,i) = MATMUL(usr%BD_Output%BldMotion%Orientation(:,:,i), usr%BD_MiscVar%BldInternalForceQP(1:3,i))
         !F(4:6,i) = MATMUL(usr%BD_Output%BldMotion%Orientation(:,:,i), usr%BD_MiscVar%BldInternalForceQP(4:6,i))
         ! -- Expressed in root frame --
         F(1:3,i) = MATMUL(usr%BD_MiscVar%u2%RootMotion%Orientation(:,:,1), usr%BD_MiscVar%BldInternalForceQP(1:3,i))
         F(4:6,i) = MATMUL(usr%BD_MiscVar%u2%RootMotion%Orientation(:,:,1), usr%BD_MiscVar%BldInternalForceQP(4:6,i))
      END SELECT
   END DO

   END SUBROUTINE BeamDyn_C_getReactionForce



   SUBROUTINE BD_initFromUsrData(usr)

   TYPE(BD_UsrDataType) :: usr

   INTEGER(IntKi)               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                      ! Temporary Error message
   character(*), parameter      :: RoutineName = 'BD_initFromUsrData'

   !---------------------- SIMULATION CONTROL --------------------------------------
   ! Set : DynamicSolve, t_initial, t_final, nt
      
   !---------------------- GRAVITY PARAMETER --------------------------------------
   usr%BD_InitInput%gravity(:) = 0.0_ReKi
   usr%BD_InitInput%gravity(1:3) = usr%grav(1:3)
   
   !---------------------- FRAME PARAMETER --------------------------------------
   ! Not set: GlbRotBladeT0
   usr%BD_InitInput%GlbPos(:)   = 0.0_ReKi
   usr%BD_InitInput%RootOri(:,:) = 0.0_R8Ki

   usr%BD_InitInput%GlbPos(1:3)      = usr%GlbPos(1:3)
   usr%BD_InitInput%RootOri(1:3,1:3) = usr%RootOri(1:3,1:3)

   usr%BD_InitInput%HubPos = 0.0_ReKi
   usr%BD_InitInput%HubRot = 0.0_R8Ki

      ! Use the initial blade root orientation as the GlbRot reference orientation for all calculations?
   IF ( usr%GlbRotBladeT0 ) THEN
      ! Set the GlbRot matrix
      usr%BD_InitInput%GlbRot = usr%BD_InitInput%RootOri
      CALL eye( usr%RootRelInit, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSE
      ! Initialize the GlbRot matrix as the identity.  Relative rotation for root to GlbRot
      usr%RootRelInit = usr%BD_InitInput%RootOri
      CALL eye( usr%BD_InitInput%GlbRot, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF

   usr%GlbRot(1:3,1:3)     = usr%BD_InitInput%GlbRot(1:3,1:3)
   usr%BD_InitInput%HubRot = usr%BD_InitInput%GlbRot(1:3,1:3)

   !---------------------- INITIAL VELOCITY PARAMETER --------------------------------
   usr%BD_InitInput%RootVel(4:6) = usr%omega(1:3)
   usr%BD_InitInput%RootVel(1:3) = cross_product(usr%BD_InitInput%RootVel(4:6),usr%BD_InitInput%GlbPos(:))

   ! Path to the structure input file
   usr%BD_InitInput%InputFile = usr%InputFile

   !---------------------- Outputs ---------------------------------------------------------
   !check that usr%WrVTK, usr%VTK_fps are set correctly
   if (usr%WrVTK < 0 .or. usr%WrVTK > 2) then
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   endif

   ! initialize the usr%BD_InitInput values not in the driver input file
   usr%BD_InitInput%RootName     = TRIM(usr%BD_Initinput%InputFile)
   usr%BD_InitInput%RootDisp     = MATMUL(usr%BD_InitInput%GlbPos(:),usr%RootRelInit) - usr%BD_InitInput%GlbPos(:)
   usr%BD_InitInput%RootVel(1:3) = Cross_Product(usr%BD_InitInput%RootVel(4:6),usr%BD_InitInput%GlbPos(:))
   usr%BD_InitInput%DynamicSolve = usr%DynamicSolve      ! QuasiStatic options handled within the BD code.
       
END SUBROUTINE BD_initFromUsrData


SUBROUTINE BDuser_InitRotationCenterMesh(usr)
   TYPE(BD_UsrDataType),        INTENT(INOUT) :: usr

   integer(intKi)                             :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                       :: ErrMsg2           ! temporary Error message
   character(*), parameter                    :: RoutineName = 'BDuser_InitRotationCenterMesh'
   
   real(r8Ki)                                 :: orientation(3,3)
   real(ReKi)                                 :: position(3)
   real(R8Ki)                                 :: z_hat(3)  ! unit-magnitude rotational velocity vector
   real(R8Ki)                                 :: Z_unit(3) ! unit vector in the Z direction
   real(R8Ki)                                 :: vec(3)    ! temporary vector
   real(R8Ki)                                 :: w         ! norm of omega

   ErrStat = ErrID_None
   ErrMsg = ''
   
   
   position = 0.0_ReKi  ! center of rotation
   
   w = TwoNorm( usr%BD_InitInput%RootVel(4:6) )
   
   if (EqualRealNos(w,0.0_R8Ki)) then
      w = 0.0_R8Ki
         ! the beam is not rotating, so pick an orientation
      call eye(orientation, ErrStat2, ErrMsg2)
   else
      z_hat = usr%BD_InitInput%RootVel(4:6) / w
      
      if ( EqualRealNos( z_hat(3), 1.0_R8Ki ) ) then
         call eye(orientation, ErrStat2, ErrMsg2)
      elseif ( EqualRealNos( z_hat(3), -1.0_R8Ki ) ) then
         orientation = 0.0_ReKi
         orientation(1,1) = -1.0_R8Ki
         orientation(2,2) =  1.0_R8Ki
         orientation(3,3) = -1.0_R8Ki
      else
         
         Z_unit = (/0.0_R8Ki, 0.0_R8Ki, 1.0_R8Ki/)
         
         vec = Z_unit - z_hat*z_hat(3) ! vec = matmul( eye(3) - outerproduct(z_hat,z_hat), (/ 0,0,1/) )
         vec = vec / TwoNorm(vec)      ! we've already checked that this is not zero
         orientation(1,:) = vec

         vec = cross_product(z_hat,Z_unit)
         vec = vec / TwoNorm(vec)
         orientation(2,:) = vec
                  
         orientation(3,:) = z_hat
         
      end if
      
   end if
   ! Overwrite with eye matrix
   call eye(orientation, ErrStat2, ErrMsg2)
   
   !.......................
   ! Mesh for center of rotation
   !.......................
   CALL MeshCreate( BlankMesh        = usr%RotationCenter &
                   ,IOS              = COMPONENT_OUTPUT          &
                   ,NNodes           = 1                         &
                   ,TranslationDisp  = .TRUE.                    &
                   ,Orientation      = .TRUE.                    &
                   ,TranslationVel   = .TRUE.                    &
                   ,RotationVel      = .TRUE.                    &
                   ,TranslationAcc   = .TRUE.                    &
                   ,RotationAcc      = .TRUE.                    &
                   ,ErrStat          = ErrStat2                  &
                   ,ErrMess          = ErrMsg2                    )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) return   
      
       ! set the reference position and orientation.
   CALL MeshPositionNode ( usr%RotationCenter, 1, position, ErrStat2, ErrMsg2, orientation )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL MeshConstructElement( usr%RotationCenter, ELEMENT_POINT, ErrStat2, ErrMsg2, 1 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL MeshCommit (usr%RotationCenter, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) return   
   
      
   ! note that the following fields do not change during the simulation;  orientation will change in BD_InputSolve.
   usr%RotationCenter%TranslationDisp  = 0.0_ReKi
   usr%RotationCenter%TranslationVel   = 0.0_ReKi
   usr%RotationCenter%RotationVel(:,1) = usr%BD_InitInput%RootVel(4:6)
   usr%RotationCenter%TranslationAcc   = 0.0_ReKi
   usr%RotationCenter%RotationAcc      = 0.0_ReKi
   
   
   !.......................
   ! initialize the mapping between the BD center of rotation and BD root motion input mesh:
   !.......................
      
   CALL MeshMapCreate( usr%RotationCenter, usr%BD_Input(1)%RootMotion, usr%Map_RotationCenter_to_RootMotion, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   usr%RotationCenter%remapFlag = .false.
   usr%BD_Input(1)%RootMotion%remapFlag = .false.
   
END SUBROUTINE BDuser_InitRotationCenterMesh

SUBROUTINE BDUsr_InputSolve(usr,i)
 
   TYPE(BD_UsrDataType) :: usr

   INTEGER,INTENT(IN) :: i
                                         
   ! local variables
   REAL(R8Ki)                                 :: Orientation(3,3)
   REAL(ReKi)                                 :: w,wt,theta,omega(3)    ! time from start start of simulation multiplied by magnitude of rotational velocity
   REAL(ReKi)                                 :: omegaP(3),alphaP(3)    ! pitch rotation velocity and acceleration
   REAL(ReKi)                                 :: temp_R(3,3), temp_cc(3), u_theta_pitch

   INTEGER(IntKi)                             :: j
   
   integer(intKi)                             :: ErrStat2         ! temporary Error status
   character(ErrMsgLen)                       :: ErrMsg2          ! temporary Error message
   character(*), parameter                    :: RoutineName = 'BDUsr_InputSolve'
   
   ErrStat = ErrID_None
   ErrMsg  = ''

   !.............................
   ! Set up u%RootMotion: 
   !.............................
   ! Compute the orientation at the center of rotation 
   w  = TwoNorm(usr%omega)
   wt = w * REAL((usr%BD_InputTimes(i) - usr%t),ReKi) ! Increment of rotation angle
   theta = usr%theta_rot + wt ! Current rotation angle (based on history + increment of rotation)
   IF (w .eq. 0.0_ReKi) THEN
      CALL eye(Orientation,ErrStat2,ErrMsg2)
   ELSE
      omega(1:3) = usr%omega(1:3) / w

      ! Rotation matrix obtained using Rodrigues formula: R = I + sin(theta) tilde(omega) + (1-cos(theta)) tilde(omega)^2
      CALL eye(Orientation,ErrStat2,ErrMsg);
      Orientation  = Orientation + sin(theta) * SkewSymMat(omega) + (1-cos(theta)) * MATMUL(SkewSymMat(omega),SkewSymMat(omega))
      Orientation  =  TRANSPOSE(Orientation)
   ENDIF

   if (i==1) then
      usr%theta_rot = usr%theta_rot + wt ! Update the rotation angle if we are not computing an estimate at another BD_inputTime
   endif
      
   usr%RotationCenter%Orientation(:,:,1) = matmul(Orientation, matmul(usr%RotationCenter%RefOrientation(:,:,1),usr%RootRelInit))
   usr%RotationCenter%RotationVel(:,1)   = usr%omega(1:3)
   usr%RotationCenter%RotationAcc(:,1)   = usr%dOmega(1:3)

   CALL Transfer_Point_to_Point( usr%RotationCenter, usr%BD_Input(i)%RootMotion, usr%Map_RotationCenter_to_RootMotion, ErrStat2, ErrMsg2)  
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   !.............................
   ! set up the point load input:
   ! @VA: if we want to apply these at different positions, we should call Transfer_MultipointLoads(); 
   ! putting the calculation here so we can have changing loads at some point...
   !.............................
   ! CALL Transfer_Point_to_Point( DvrData%mplLoads, u%PointLoad, DvrData%Map_mplLoads_to_PointLoad, ErrStat2, ErrMsg2, DvrData%mplMotion, DvrData%y_BldMotion_at_u_point)  
   !    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ! 
   ! u%PointLoad%Force(1:3,u%PointLoad%NNodes)  = u%PointLoad%Force(1:3,u%PointLoad%NNodes)  + DvrData%TipLoad(1:3)
   ! u%PointLoad%Moment(1:3,u%PointLoad%NNodes) = u%PointLoad%Moment(1:3,u%PointLoad%NNodes) + DvrData%TipLoad(4:6)

   !.............................
   ! Pitch angle input
   !.............................
   ! Rotation around pitch axis that must be added to the blade root motion.
   ! This must be done after the transfer from the rotation center mesh!

   ! If the rotationVel and acceleration due to pitch are important for the structural dynamics, they should be provided. Here, forced to zero
   omegaP = 0.0
   alphaP = 0.0
   usr%BD_Input(i)%RootMotion%RotationVel(:,1) = usr%BD_Input(i)%RootMotion%RotationVel(:,1) - omegaP * usr%BD_Input(i)%RootMotion%Orientation(3,:,1)
   usr%BD_Input(i)%RootMotion%RotationAcc(:,1) = usr%BD_Input(i)%RootMotion%RotationAcc(:,1) - alphaP * usr%BD_Input(i)%RootMotion%Orientation(3,:,1)

   ! Add the pitch angle
   temp_cc    = 0.0_ReKi
   temp_cc(3) = -usr%PAngInp_rad
   temp_R     = EulerConstruct(temp_cc)
   usr%BD_Input(i)%RootMotion%Orientation(:,:,1) = MATMUL(temp_R,usr%BD_Input(i)%RootMotion%Orientation(:,:,1)) ! Post-multiply around z

   !.............................
   ! LINE2 mesh: DistrLoad
   !.............................
   DO j=1,usr%BD_Input(i)%DistrLoad%NNodes
      usr%BD_Input(i)%DistrLoad%Force(:,j) =  MATMUL(TRANSPOSE(usr%BD_Input(i)%RootMotion%Orientation(:,:,1)), usr%loads(1:3,j))
      usr%BD_Input(i)%DistrLoad%Moment(:,j)=  MATMUL(TRANSPOSE(usr%BD_Input(i)%RootMotion%Orientation(:,:,1)), usr%loads(4:6,j))
   ENDDO

   ! ! This is just for check: 
   ! temp_R = MATMUL(usr%BD_Input(i)%RootMotion%Orientation(:,:,1),TRANSPOSE(usr%BD_Input(i)%HubMotion%Orientation(:,:,1)))
   ! !write(*,*) "Hub ori = ",usr%BD_Input(i)%HubMotion%Orientation(:,:,1)
   ! !write(*,*) temp_R
   ! write(*,*) "Root_orientation:",usr%BD_Input(i)%RootMotion%Orientation(:,:,1)
   ! temp_cc = EulerExtract(temp_R)   != Hub_theta_Root
   ! u_theta_pitch = temp_cc(3)

   ! WRITE(*,'(A,F10.5,F10.5)') "Effective Rotation angle VS demanded:",temp_cc(1), usr%theta_rot
   ! WRITE(*,'(A,F10.5)') "Effective y angle:",temp_cc(2)
   ! WRITE(*,'(A,F10.5)') "Effective pitch angle input:",u_theta_pitch
   ! WRITE(*,'(A,F10.5)') "Actual pitch               :",usr%BD_DiscreteState%thetaP


END SUBROUTINE BDUsr_InputSolve

SUBROUTINE BDUsr_UpdateOrientation(usr)

   TYPE(BD_UsrDataType) :: usr
   REAL(ReKi)                                 :: w,theta,omega(3),dOri(3,3)

   integer(intKi)                             :: ErrStat2         ! temporary Error status
   character(ErrMsgLen)                       :: ErrMsg2          ! temporary Error message

   w     = TwoNorm(usr%omega)
   IF (w .eq. 0.0_ReKi) THEN
      CALL eye(dOri,ErrStat2,ErrMsg2)
   ELSE
      omega(1:3) = usr%omega(1:3) / w
      theta = w * REAL(usr%dt,ReKi)

      ! Rotation matrix obtained using Rodrigues formula: R = I + sin(theta) tilde(omega) + (1-cos(theta)) tilde(omega)^2
      CALL eye(dOri,ErrStat2,ErrMsg);
      dOri  = dOri + sin(theta) * SkewSymMat(omega) + (1-cos(theta)) * MATMUL(SkewSymMat(omega),SkewSymMat(omega))
      dOri  =  TRANSPOSE(dOri)
   ENDIF

   !! Add the new orientation
   !usr%orientation = matmul(usr%orientation,dOri);

END SUBROUTINE

! Routine from the BeamDyn driver
SUBROUTINE Dvr_InitializeOutputFile(OutUnit,IntOutput,RootName,ErrStat,ErrMsg)


   INTEGER(IntKi),              INTENT(  OUT):: OutUnit
   TYPE(BD_InitOutputType),     INTENT(IN   ):: IntOutput     ! Output for initialization routine
   INTEGER(IntKi),              INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),                INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None
   CHARACTER(*),                INTENT(IN   ):: RootName

   integer(IntKi)                            :: i      
   integer(IntKi)                            :: numOuts
   INTEGER(IntKi)                            :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                      :: ErrMsg2                      ! Temporary Error message
   character(*), parameter                   :: RoutineName = 'Dvr_InitializeOutputFile'

   ErrStat = ErrID_none
   ErrMsg  = ""
   
   CALL GetNewUnit(OutUnit,ErrStat2,ErrMsg2)
   CALL OpenFOutFile ( OutUnit, trim(RootName)//'.out', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

   write (OutUnit,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim(GetNVD(IntOutput%Ver))
   write (OutUnit,'()' )    !print a blank line
   
   numOuts = size(IntOutput%WriteOutputHdr)
   !......................................................
   ! Write the names of the output parameters on one line:
   !......................................................

   write (OutUnit,'()')
   write (OutUnit,'()')
   write (OutUnit,'()')

   call WrFileNR ( OutUnit, 'Time' )

   do i=1,NumOuts
      call WrFileNR ( OutUnit, tab//IntOutput%WriteOutputHdr(i) )
   end do ! i

   write (OutUnit,'()')

      !......................................................
      ! Write the units of the output parameters on one line:
      !......................................................

   call WrFileNR ( OutUnit, '(s)' )

   do i=1,NumOuts
      call WrFileNR ( Outunit, tab//trim(IntOutput%WriteOutputUnt(i)) )
   end do ! i

   write (OutUnit,'()')  


END SUBROUTINE Dvr_InitializeOutputFile

! Routine from the beamdyn driver
SUBROUTINE Dvr_WriteOutputLine(t,OutUnit, OutFmt, Output)


   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   INTEGER(IntKi)         ,  intent(in   )   :: OutUnit              ! Status of error message
   CHARACTER(*)           ,  intent(in   )   :: OutFmt
!   real(ReKi)             ,  intent(in   )   :: output(:)            ! Rootname for the output file
   TYPE(BD_OutputType),      INTENT(IN   )   :: Output
      
   ! Local variables.

   integer(IntKi)                            :: errStat              ! Status of error message (we're going to ignore errors in writing to the file)
   character(ErrMsgLen)                      :: errMsg               ! Error message if ErrStat /= ErrID_None
   character(200)                            :: frmt                 ! A string to hold a format specifier
   character(15)                             :: tmpStr               ! temporary string to print the time output as text

   frmt = '"'//tab//'"'//trim(OutFmt)      ! format for array elements from individual modules
   
      ! time
   write( tmpStr, '(F15.6)' ) t
   call WrFileNR( OutUnit, tmpStr )
   call WrNumAryFileNR ( OutUnit, Output%WriteOutput,  frmt, errStat, errMsg )
   
     ! write a new line (advance to the next line)
   write (OutUnit,'()')
      
end SUBROUTINE Dvr_WriteOutputLine

! Recomputes the Time Stepping scheme coefficient when dt is updated in BD_refresh
! This is copy-pasted from BeamDyn.f90:BD_TiSchmComputeCoefficients. But the original function is PRIVATE, so I copied it here to avoid having to modify the BD sources
! If the routine in BeamDyn.90 is changed, you should changed this one too.
SUBROUTINE BD_TiSchmComputeCoefficients2(p)

      TYPE(BD_ParameterType), INTENT(INOUT) :: p 

      REAL(DbKi)                  :: tr0
      REAL(DbKi)                  :: tr1
      REAL(DbKi)                  :: tr2
      REAL(DbKi)                  :: alfam      ! \alpha_M
      REAL(DbKi)                  :: alfaf      ! \alpha_F
      REAL(DbKi)                  :: gama
      REAL(DbKi)                  :: beta
      REAL(DbKi)                  :: oalfaM     ! 1 - \alpha_M
      REAL(DbKi)                  :: deltat2    ! {\delta t}^2 = dt^2

      ! Bauchau equations 17.39
      tr0 = p%rhoinf + 1.0_BDKi
      alfam = (2.0_BDKi * p%rhoinf - 1.0_BDKi) / tr0
      alfaf = p%rhoinf / tr0

      ! Bauchau equations 17.40
      gama = 0.5_BDKi - alfam + alfaf
      beta = 0.25 * (1.0_BDKi - alfam + alfaf)**2

      ! The coefficents are then found using equations 17.41a - 17.41c
      deltat2 = p%dt**2
      oalfaM = 1.0_BDKi - alfam
      tr0 = alfaf / oalfaM
      tr1 = alfam / oalfaM
      tr2 = (1.0_BDKi - alfaf) / oalfaM

      p%coef(1) = beta * tr0 * deltat2
      p%coef(2) = (0.5_BDKi - beta/oalfaM) * deltat2
      p%coef(3) = gama * tr0 * p%dt
      p%coef(4) = (1.0_BDKi - gama / oalfaM) * p%dt
      p%coef(5) = tr0
      p%coef(6) = -tr1
      p%coef(7) = gama * tr2 * p%dt
      p%coef(8) = beta * tr2 * deltat2
      p%coef(9) = tr2

END SUBROUTINE

END MODULE BeamDynLib