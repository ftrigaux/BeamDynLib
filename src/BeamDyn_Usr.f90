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
MODULE BeamDyn_Usr

   !USE BeamDyn_driver_subs  ! all other modules inherited through this one
   USE BeamDyn
   USE BeamDyn_Subs

   IMPLICIT NONE

   !PUBLIC :: BeamDyn_C_Init
   !PUBLIC :: BeamDyn_C_SetLoads
   !PUBLIC :: BeamDyn_C_Solve
   !PUBLIC :: BeamDyn_C_GetDisp
   !PUBLIC :: BeamDyn_C_End

   TYPE , PUBLIC :: BD_UsrDataType
      LOGICAL        :: DynamicSolve  ! flag for dynamic or static solve (static:false)
      LOGICAL        :: GlbRotBladeT0 ! Initial blade root orientation is also the GlbRot reference frame

      INTEGER(IntKi)                          :: nxD,nxL;        ! Number of nodes for Displacement and Loads
      !REAL(ReKi),DIMENSION(:,:), ALLOCATABLE  :: x;         ! Nodes coordinates

      REAL(DbKi)     :: t             ! time
      REAL(DbKi)     :: dt            ! time increment
      INTEGER(IntKi) :: nt            ! number of substep
      
      REAL(DbKi)     :: GlbPos(3)        ! Initial vector position 
      REAL(DbKi)     :: RootOri(3,3)     ! DCM of the initial root orientation
      REAL(DbKi)     :: GlbRot(3,3)      ! 
      REAL(DbKi)     :: RootRelInit(3,3)

      REAL(DbKi)     :: orientation(3,3)   ! Orientation angles
      REAL(DbKi)     :: omega(3)           ! Angular velocity vector
      REAL(DbKi)     :: dOmega(3)          ! Angular acceleration vector

      REAL(ReKi),DIMENSION(:,:), ALLOCATABLE         :: loads    ! Forces and moment at the node positions; Size(NNodes,6)
      REAL(DbKi)                                     :: grav(3)  ! Gravity vector

      !REAL(DbKi),DIMENSION(:,:), ALLOCATABLE         :: u        ! Displacement variables (u,v,w,phi,th1,th2); Size(NNodes,6)
      !REAL(DbKi),DIMENSION(:,:), ALLOCATABLE         :: du       ! Velocity variables d(u,v,w,phi,th1,th2)/dt; Size(NNodes,6)

      CHARACTER(1024)  :: InputFile      !< Name of the input file; remove if there is no file [-]

      ! Contains the variables for the BeamDyn Analysis
      TYPE(BD_InitInputType)           :: BD_InitInput
      TYPE(BD_ParameterType)           :: BD_Parameter
      TYPE(BD_ContinuousStateType)     :: BD_ContinuousState
      TYPE(BD_InitOutputType)          :: BD_InitOutput
      TYPE(BD_DiscreteStateType)       :: BD_DiscreteState
      TYPE(BD_ConstraintStateType)     :: BD_ConstraintState
      TYPE(BD_OtherStateType)          :: BD_OtherState
      TYPE(BD_MiscVarType)             :: BD_MiscVar
      TYPE(BD_InputType) ,ALLOCATABLE  :: BD_Input(:)
      REAL(DbKi),         ALLOCATABLE  :: BD_InputTimes(:)
      TYPE(BD_OutputType)              :: BD_Output
      !INTEGER(IntKi)                   :: DvrOut 

      ! Private variables used inside BeamDyn code
      REAL(R8Ki)                                    :: w           ! magnitude of rotational velocity vector
      TYPE(MeshType)                                :: RotationCenter
      TYPE(MeshMapType)                             :: Map_RotationCenter_to_RootMotion

   END TYPE

   ! global glue-code-specific variables

   INTEGER(IntKi)                   :: ErrStat          ! Error status of the operation
   CHARACTER(1024)                  :: ErrMsg           ! Error message if ErrStat /= ErrID_None
   INTEGER(IntKi), parameter        :: BD_interp_order = 1  ! order of interpolation/extrapolation

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
   !!INTEGER(IntKi)                   :: DvrOut 
   
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
      !   CALL CheckError(usr)
      CALL BD_initFromUsrData(usr);
      CALL CheckError(usr)

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
         CALL CheckError(usr)
      
         ! If the Quasi-Static solve is in use, rerun the initialization with loads at t=0 
         ! (HACK: set in the driver only because computing Jacobians with this option [as in FAST glue code] is problematic)
      usr%BD_OtherState%RunQuasiStaticInit = usr%BD_Parameter%analysis_type == BD_DYN_SSS_ANALYSIS


      ! Set the Initial root orientation
      ! BD_Input(1)%RootMotion%Orientation(1:3,1:3,1) = DvrData%RootRelInit
      usr%BD_Input(1)%RootMotion%Orientation(1:3,1:3,1) = usr%RootRelInit
         
      CALL BDuser_InitRotationCenterMesh(usr)
         CALL CheckError(usr)

      CALL initUsrData(usr)

      ! 
      ! CALL CreateMultiPointMeshes(DvrData,BD_InitInput,BD_InitOutput,BD_Parameter, BD_Output, BD_Input(1), ErrStat, ErrMsg)   
      ! CALL Transfer_MultipointLoads(DvrData, BD_Output, BD_Input(1), ErrStat, ErrMsg)   
      
      ! CALL Dvr_InitializeOutputFile(DvrOut,BD_InitOutput,RootName,ErrStat,ErrMsg)
      !   CALL CheckError(usr)

         
         ! initialize BD_Input and BD_InputTimes
      usr%BD_InputTimes(1) = usr%t
      !CALL BD_InputSolve( BD_InputTimes(1), BD_Input(1), DvrData, ErrStat, ErrMsg)
      CALL BDUsr_InputSolve(usr,1)
      
      DO j = 2,BD_interp_order+1
            ! create new meshes
         CALL BD_CopyInput (usr%BD_Input(1) , usr%BD_Input(j) , MESH_NEWCOPY, ErrStat, ErrMsg)
            CALL CheckError(usr)
            
            ! solve for inputs at previous time steps
         usr%BD_InputTimes(j) = usr%t - (j - 1) * usr%dt
         ! CALL BD_InputSolve( BD_InputTimes(j), BD_Input(j), DvrData, ErrStat, ErrMsg)
         CALL BDUsr_InputSolve(usr,j)
            CALL CheckError(usr)
      END DO


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
      CALL CheckError(usr)
   
   ! CALL Dvr_WriteOutputLine(bd_usrdata%t,DvrOut,BD_Parameter%OutFmt,BD_Output)
   
      !.........................
      ! time marching
      !.........................
     
   
   DO it = 1, usr%nt  !Changed 0,nt to 1,nt

      ! Shift "window" of BD_Input 
      DO j = BD_interp_order, 1, -1
         CALL BD_CopyInput (usr%BD_Input(j),  usr%BD_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
            CALL CheckError(usr)
         usr%BD_InputTimes(j+1) = usr%BD_InputTimes(j)
      END DO
      
      usr%BD_InputTimes(1)  = usr%t + usr%dt
      CALL BDUsr_UpdateOrientation(usr)
      !CALL BD_InputSolve( BD_InputTimes(1), BD_Input(1), DvrData, ErrStat, ErrMsg)
      CALL BDUsr_InputSolve(usr,1)
         CALL CheckError(usr)
      
                       
     IF(usr%BD_Parameter%analysis_type .EQ. BD_STATIC_ANALYSIS .AND. it > max_ld_step) EXIT

      ! update states from n_t_global to n_t_global + 1
     CALL BD_UpdateStates( usr%t, it, usr%BD_Input, usr%BD_InputTimes, usr%BD_Parameter, &
                               usr%BD_ContinuousState, &
                               usr%BD_DiscreteState, usr%BD_ConstraintState, &
                               usr%BD_OtherState, usr%BD_MiscVar, ErrStat, ErrMsg )
        CALL CheckError(usr)

        
      ! advance time
     usr%t = usr%t + usr%dt
           
      ! calculate outputs at n_t_global + 1
     CALL BD_CalcOutput( usr%t, usr%BD_Input(1), usr%BD_Parameter, usr%BD_ContinuousState, usr%BD_DiscreteState, &
                             usr%BD_ConstraintState, usr%BD_OtherState,  usr%BD_Output, usr%BD_MiscVar, ErrStat, ErrMsg)
        CALL CheckError(usr)

     !CALL Dvr_WriteOutputLine(t_global,DvrOut,BD_Parameter%OutFmt,BD_Output)

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

      !IF(DvrOut >0) CLOSE(DvrOut)

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

   SUBROUTINE CheckError(usr)
      TYPE(BD_UsrDataType) :: usr

      if (ErrStat /= ErrID_None) then
         CALL WrScr(TRIM(ErrMsg))
         
         if (ErrStat >= AbortErrLev) then
            CALL BeamDyn_C_End(usr)
         end if
      end if
         
   end SUBROUTINE CheckError


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
   usr%BD_InitInput%GlbRot(:,:) = 0.0_R8Ki
   usr%BD_InitInput%RootOri(:,:) = 0.0_R8Ki

   usr%BD_InitInput%GlbPos(1:3)      = usr%GlbPos(1:3)
   usr%BD_InitInput%RootOri(1:3,1:3) = usr%RootOri(1:3,1:3)
   usr%BD_InitInput%GlbRot(1:3,1:3)  = usr%GlbRot(1:3,1:3)

   usr%orientation(1:3,1:3) = usr%RootOri(1:3,1:3)

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

   !---------------------- INITIAL VELOCITY PARAMETER --------------------------------
   usr%BD_InitInput%RootVel(4:6) = usr%omega(1:3)
   usr%BD_InitInput%RootVel(1:3) = cross_product(usr%BD_InitInput%RootVel(4:6),usr%BD_InitInput%GlbPos(:))

   ! Path to the structure input file
   usr%BD_InitInput%InputFile = usr%InputFile

   ! initialize the usr%BD_InitInput values not in the driver input file
   usr%BD_InitInput%RootName     = TRIM(usr%BD_Initinput%InputFile)
   usr%BD_InitInput%RootDisp     = matmul(transpose(usr%BD_InitInput%RootOri),usr%BD_InitInput%GlbPos) - usr%BD_InitInput%GlbPos
   usr%BD_InitInput%RootVel(1:3) = matmul(usr%BD_InitInput%RootOri, Cross_Product( usr%BD_InitInput%RootVel(4:6), usr%BD_InitInput%GlbPos ))  ! set translational velocities based on rotation and GlbPos.
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
   
   ErrStat = ErrID_None
   ErrMsg = ''
   
   
   position = 0.0_ReKi  ! center of rotation

   usr%w = TwoNorm( usr%BD_InitInput%RootVel(4:6) )
   
   call eye(orientation, ErrStat2, ErrMsg2) ! new version
   
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
   REAL(DbKi)                                 :: dt,w,theta,omega(3)          ! time from start start of simulation multiplied by magnitude of rotational velocity
   REAL(DbKi)                                 :: dOri(3,3)

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
   dt  = usr%BD_InputTimes(i) - usr%t
   IF (i == 1) THEN
      Orientation(1:3,1:3) = usr%orientation(1:3,1:3);
   ELSE
      w     = TwoNorm(usr%omega)
      IF (w .eq. 0.0_DbKi) THEN
         CALL eye(dOri,ErrStat2,ErrMsg2)
      ELSE
         theta = - w * (i-1) * usr%dt
         omega(1:3) = usr%omega(1:3) / w

         ! Rotation matrix obtained using Rodrigues formula: R = I + sin(theta) tilde(omega) + (1-cos(theta)) tilde(omega)^2
         CALL eye(dOri,ErrStat2,ErrMsg);
         dOri  = dOri + sin(theta) * SkewSymMat(omega) + (1-cos(theta)) * MATMUL(SkewSymMat(omega),SkewSymMat(omega))
         dOri = TRANSPOSE(dOri)
      ENDIF
      Orientation = matmul(usr%orientation,dOri);
   ENDIF
      
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
   ! LINE2 mesh: DistrLoad
   !.............................
   DO j=1,usr%BD_Input(i)%DistrLoad%NNodes
      usr%BD_Input(i)%DistrLoad%Force(:,j) =  MATMUL(TRANSPOSE(Orientation), usr%loads(1:3,j))
      usr%BD_Input(i)%DistrLoad%Moment(:,j)=  MATMUL(TRANSPOSE(Orientation), usr%loads(4:6,j))
   ENDDO

END SUBROUTINE BDUsr_InputSolve

SUBROUTINE BDUsr_UpdateOrientation(usr)

   TYPE(BD_UsrDataType) :: usr
   REAL(DbKi)                                 :: w,theta,omega(3),dOri(3,3)

   integer(intKi)                             :: ErrStat2         ! temporary Error status
   character(ErrMsgLen)                       :: ErrMsg2          ! temporary Error message

   w     = TwoNorm(usr%omega)
   IF (w .eq. 0.0_DbKi) THEN
      CALL eye(dOri,ErrStat2,ErrMsg2)
   ELSE
      omega(1:3) = usr%omega(1:3) / w
      theta = w * usr%dt

      ! Rotation matrix obtained using Rodrigues formula: R = I + sin(theta) tilde(omega) + (1-cos(theta)) tilde(omega)^2
      CALL eye(dOri,ErrStat2,ErrMsg);
      dOri  = dOri + sin(theta) * SkewSymMat(omega) + (1-cos(theta)) * MATMUL(SkewSymMat(omega),SkewSymMat(omega))
      dOri  =  TRANSPOSE(dOri)
   ENDIF

   ! Add the new orientation
   usr%orientation = matmul(usr%orientation,dOri);

END SUBROUTINE

END MODULE BeamDyn_Usr